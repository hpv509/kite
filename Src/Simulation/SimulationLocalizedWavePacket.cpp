#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include <cmath>
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
#include "Loop.hpp"
#include "Coefficients.hpp"
#include <stdexcept>
#include <algorithm>

// First entry: Row index of the probe in ProbeCoordinates
// Second entry: Coordinate index in local domain with ghosts
using LocalProbe = std::pair<std::size_t, std::size_t>;

#define INITIATE_PACKET                                                        \
  do {                                                                         \
    phi.initiate_phases();                                                     \
    phi.set_index(0);                                                          \
                                                                               \
    Coordinates<std::size_t, D + 1> global(r.Lt);                              \
    const bool gaussian = width > 0.;                                          \
                                                                               \
    if (!gaussian) {                                                           \
      auto x = static_cast<unsigned>(pos_[0]);                                 \
      auto y = static_cast<unsigned>(pos_[1]);                                 \
      auto o = static_cast<unsigned>(pos_[D]);                                 \
                                                                               \
      if constexpr (D == 2) {                                                  \
        global.set({x, y, o});                                                 \
      } else if constexpr (D == 3) {                                           \
        auto z = static_cast<unsigned>(pos_[2]);                               \
        global.set({x, y, z, o});                                              \
      }                                                                        \
                                                                               \
      phi.build_site(global.index);                                            \
    } else {                                                                   \
      Eigen::Matrix<T, -1, -1> psi0(r.Orb, 1);                                 \
      Eigen::Matrix<double, -1, -1> k0(D, 1);                                  \
      Eigen::Matrix<double, 1, D> r0;                                          \
                                                                               \
      psi0.setOnes();                                                          \
      for (unsigned i = 0; i < D; ++i) {                                       \
        k0(i, 0) = k0_[i];                                                     \
        r0(0, i) = pos_[i];                                                    \
      }                                                                        \
                                                                               \
      phi.build_wave_packet(k0, psi0, width, r0);                              \
    }                                                                          \
                                                                               \
    const bool filter = energy_window[0] != 0. || energy_window[1] != 0.;      \
    if (filter) {                                                              \
      const auto Emin = std::clamp<                                            \
        value_type>((energy_window[0] - energy_shift) / energy_scale, -1, 1);  \
      const auto Emax = std::clamp<                                            \
        value_type>((energy_window[1] - energy_shift) / energy_scale, -1, 1);  \
                                                                               \
      const Eigen::Array<value_type, -1, 1> coefs =                            \
        Coefficients::build_window<                                            \
          value_type>((Emin + Emax) * 0.5, std::abs(Emax - Emin), 256);        \
                                                                               \
      Eigen::Array<T, -1, 1> filtered(r.Sized);                                \
      filtered.setZero();                                                      \
                                                                               \
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {                     \
        phi.cheb_iteration(n);                                                 \
        filtered += coefs(n) * phi.v.col(phi.get_index()).array();             \
      }                                                                        \
      phi.v.col(0) = filtered;                                                 \
                                                                               \
      _Pragma("omp barrier") _Pragma("omp master") Global.soma = 0;            \
      _Pragma("omp barrier")                                                   \
                                                                               \
        _Pragma("omp critical")                                                \
          Global.soma += no_ghost_dot<T, D>(r, phi.v.col(0), phi.v.col(0));    \
                                                                               \
      _Pragma("omp barrier") phi.v.col(0) /= std::sqrt(Global.soma);           \
    }                                                                          \
                                                                               \
    phi.set_index(0);                                                          \
    phi.Exchange_Boundaries();                                                 \
  } while (0)

template <unsigned D, bool Global, typename F>
void for_each_orbital(LatticeStructure<D> &r, F &&func)
{
  Coordinates<std::size_t, D + 1> local(r.Ld);
  Coordinates<std::size_t, D + 1> global(r.Lt);
  std::array<unsigned, D> idx{};
  std::array<unsigned, D> start{};
  std::array<unsigned, D> final{};
  for (unsigned d = 0; d < D; ++d) {
    start[d] = NGHOSTS;
    final[d] = r.Ld[D - 1 - d] - NGHOSTS;
  }

  for (unsigned io = 0; io < r.Orb; ++io) {
    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2) {
        local.set({i[1], i[0], io});
      } else if constexpr (D == 3) {
        local.set({i[2], i[1], i[0], io});
      }

      if constexpr (Global) {
        r.convertCoordinates(global, local);
        func(local, global, io, i);
      } else {
        func(local, io, i);
      }
    };

    UnitCellLoop<D>::run(idx, start, final, body);
  }
}

template <typename T, unsigned D, typename V>
std::pair<Eigen::Array<V, D, 1>, Eigen::Array<V, D * D, 1>>
pos_moments(LatticeStructure<D> &r, const Eigen::Array<T, -1, 1> &psi)
{
  // Returns [< x >, < y >, < z >] and the second moments matrix in real space coordinates for a normalized |psi>
  // Does not consider the system's periodicity
  Eigen::Array<V, D, 1> moments1 = Eigen::Array<V, D, 1>::Zero();
  Eigen::Array<V, D, D> moments2 = Eigen::Array<V, D, D>::Zero();
  for_each_orbital<D, true>(
    r, [&](
         const Coordinates<std::size_t, D + 1> &local,
         const Coordinates<std::size_t, D + 1> &global, unsigned io,
         const std::array<unsigned, D> &
       ) {
      Eigen::Matrix<V, D, 1> R;
      for (unsigned j = 0; j < D; ++j) {
        R(j) = static_cast<V>(global.coord[j]);
      }

      Eigen::Array<V, D, 1> pos =
        (r.rLat.template cast<V>() * R + r.rOrb.col(io).template cast<V>())
          .array();

      const auto prob = static_cast<V>(std::norm(psi(local.index)));
      moments1 += prob * pos;
      moments2 += prob * (pos.matrix() * pos.matrix().transpose()).array();
    }
  );

  return {moments1, moments2.template reshaped<Eigen::ColMajor>(D * D, 1)};
}

template <typename T, unsigned D>
T no_ghost_dot(
  LatticeStructure<D> &r,
  const Eigen::Array<T, -1, 1> &left,
  const Eigen::Array<T, -1, 1> &right
)
{
  T result{0};
  for_each_orbital<D, false>(
    r, [&](
         const Coordinates<std::size_t, D + 1> &local, unsigned,
         const std::array<unsigned, D> &
       ) { result += std::conj(left(local.index)) * right(local.index); }
  );

  return result;
}

template <unsigned D>
std::vector<LocalProbe> build_local_probe_indices(
  LatticeStructure<D> &r,
  const Eigen::Array<std::size_t, -1, D + 1> &probe_coords
)
{
  const auto num_probes{static_cast<std::size_t>(probe_coords.rows())};
  std::vector<LocalProbe> local_probe_indices;

  for (std::size_t p = 0; p < num_probes; ++p) {
    Coordinates<std::size_t, D + 1> total_coords(r.Lt);
    Coordinates<std::size_t, D + 1> thread_coords(r.ld);
    Coordinates<std::size_t, D + 1> thread_coords_gh(r.Ld);
    Coordinates<std::size_t, D + 1> thread(r.nd);

    std::size_t T_thread[D + 1]{};
    std::size_t x_thread[D + 1]{};

    if constexpr (D == 2) {
      total_coords.set(
        {probe_coords(p, 0), probe_coords(p, 1), probe_coords(p, 2)}
      );
    } else if constexpr (D == 3) {
      total_coords.set(
        {probe_coords(p, 0), probe_coords(p, 1), probe_coords(p, 2),
         probe_coords(p, 3)}
      );
    }

    for (unsigned d = 0; d < D; ++d) {
    }

    // Orbital direction is not domain-decomposed
    T_thread[D] = 0;
    x_thread[D] = total_coords.coord[D];

    thread.set_index(T_thread);
    thread_coords.set_index(x_thread);

    if (thread.index == r.thread_id) {
      r.convertCoordinates(thread_coords_gh, thread_coords);
      local_probe_indices.push_back({p, thread_coords_gh.index});
    }
  }

  return local_probe_indices;
}

template <typename T, unsigned D, typename V>
Eigen::Array<V, -1, 1> spectral_function(
  LatticeStructure<D> &r,
  Simulation<T, D> &sim,
  unsigned num_moments,
  const Eigen::Array<T, -1, 1> &psi
)
{
  Eigen::Array<V, -1, 1> moments(num_moments);
  KPM_Vector<T, D> phi(2, sim);
  phi.initiate_phases();
  phi.set_index(0);
  phi.v.col(0) = psi;
  phi.Exchange_Boundaries();

  for (unsigned n = 0; n < num_moments; ++n) {
    phi.cheb_iteration(n);
    moments(n) = static_cast<V>(
      std::real(no_ghost_dot<T, D>(r, psi, phi.v.col(phi.get_index()).array()))
    );
  }

  return moments;
}

template <typename T, unsigned D>
void prepare_leads(
  LatticeStructure<D> &r,
  Hamiltonian<T, D> &h,
  const std::size_t sample_start,
  const std::size_t L
)
{
  for_each_orbital<D, true>(
    r, [&](
         const Coordinates<std::size_t, D + 1> &local,
         const Coordinates<std::size_t, D + 1> &global, unsigned io,
         const std::array<unsigned, D> &
       ) {
      const int address = h.Anderson_orb_address[io];
      if (address < 0)
        return;

      if (
        const auto x{global.coord[0]}; x < sample_start || x >= sample_start + L
      ) {
        const std::ptrdiff_t w_idx{
          static_cast<std::ptrdiff_t>(local.index) +
          (address - static_cast<std::ptrdiff_t>(io)) *
            static_cast<std::ptrdiff_t>(r.Nd)
        };
        h.U_Anderson[w_idx] = 0;
      }
    }
  );
}

template <typename T, unsigned D, typename V>
Eigen::Array<V, 3, 1> transmission_weight(
  LatticeStructure<D> &r,
  const std::size_t sample_start,
  const std::size_t L,
  const Eigen::Array<T, -1, 1> &psi
)
{
  Eigen::Array<V, 3, 1> weights = Eigen::Array<V, 3, 1>::Zero();

  for_each_orbital<D, true>(
    r, [&](
         const Coordinates<std::size_t, D + 1> &local,
         const Coordinates<std::size_t, D + 1> &global, unsigned,
         const std::array<unsigned, D> &
       ) {
      const auto x{global.coord[0]};
      const V prob = std::norm(psi(local.index));
      if (x < sample_start)
        weights(0) += prob;
      else if (sample_start <= x && x < sample_start + L)
        weights(1) += prob;
      else
        weights(2) += prob;
    }
  );

  return weights;
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_localized_wavepacket()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_localized_wavepacket\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_localized_wavepacket = false;
    try {
      int dummy_variable;
      get_hdf5<int>(
        &dummy_variable, file,
        (char *)"/Calculation/localized_wave_packet/Measurements"
      );
      Global.calculate_localized_wavepacket = true;
    } catch (H5::Exception &e) {
      debug_message("localized_wavepacket: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_wavepacket = false;
#pragma omp critical
  local_calculate_wavepacket = Global.calculate_localized_wavepacket;
#pragma omp barrier
  if (local_calculate_wavepacket) {
#pragma omp master
    std::cout << "Calculating time evolution of wave packet.\n";
#pragma omp barrier
    value_type time;
    unsigned num_measures;
    unsigned num_spectral_moments;
    std::array<value_type, D + 1> pos;
    std::array<value_type, D> k0;
    value_type width;
    std::array<value_type, 2> energy_window;
    unsigned long num_probes;
    Eigen::Array<unsigned long, -1, D + 1> probe_coords;
    unsigned long sample_start;
    unsigned long sample_L;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<value_type>(
        &time, file, (char *)"/Calculation/localized_wave_packet/Time"
      );
      get_hdf5<unsigned>(
        &num_measures, file,
        (char *)"/Calculation/localized_wave_packet/Measurements"
      );
      get_hdf5<unsigned>(
        &num_spectral_moments, file,
        (char *)"/Calculation/localized_wave_packet/NumSpectralMoments"
      );
      get_hdf5<value_type>(
        pos.data(), file,
        (char *)"/Calculation/localized_wave_packet/InitialPos"
      );
      get_hdf5<value_type>(
        k0.data(), file,
        (char *)"/Calculation/localized_wave_packet/InitialWaveVector"
      );
      get_hdf5<value_type>(
        &width, file, (char *)"/Calculation/localized_wave_packet/Width"
      );
      get_hdf5<value_type>(
        energy_window.data(), file,
        (char *)"/Calculation/localized_wave_packet/EnergyWindow"
      );
      get_hdf5<unsigned long>(
        &num_probes, file,
        (char *)"/Calculation/localized_wave_packet/NumProbes"
      );

      probe_coords.resize(num_probes, D + 1);
      if (num_probes > 0) {
        get_hdf5<unsigned long>(
          probe_coords.data(), file,
          (char *)"/Calculation/localized_wave_packet/ProbeCoordinates"
        );
      }
      get_hdf5<unsigned long>(
        &sample_start, file,
        (char *)"/Calculation/localized_wave_packet/SampleStart"
      );
      get_hdf5<unsigned long>(
        &sample_L, file,
        (char *)"/Calculation/localized_wave_packet/SampleLength"
      );

      file->close();
      delete file;
    }
    localized_wavepacket(
      time, num_measures, num_spectral_moments, pos, energy_window, k0, width,
      num_probes, probe_coords.template cast<std::size_t>(), sample_start,
      sample_L
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::localized_wavepacket(
  const value_type t,
  const unsigned measurements,
  const unsigned num_spectral_moments,
  const std::array<value_type, D + 1> &pos_,
  const std::array<value_type, 2> &energy_window,
  const std::array<value_type, D> &k0_,
  const value_type width,
  const std::size_t num_global_probes,
  const Eigen::Array<std::size_t, -1, D + 1> &probe_lattice_coords,
  const std::size_t sample_start,
  const std::size_t sample_L
)
  requires Complex<T>
{
  debug_message("Entered localized_wavepacket\n");
  value_type energy_scale;
  value_type energy_shift;
#pragma omp critical
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
    get_hdf5<value_type>(&energy_shift, file, (char *)"/EnergyShift");
    file->close();
    delete file;
  }
#pragma omp barrier
  const value_type measure_tau = energy_scale * t / measurements;
  const auto divisions = std::max<
    unsigned>(1, static_cast<unsigned>(std::ceil(measure_tau / value_type{32})));
  const value_type tau = measure_tau / divisions;
  Eigen::Array<T, -1, 1> coefs = Coefficients::build_cplx_exp<value_type>(tau);
  const T arg{0.0, -energy_shift * (tau / energy_scale)};
  coefs *= std::exp(arg);

  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(r.Sized);
  Eigen::Array<value_type, -1, -1> spectral_moments(num_spectral_moments, 2);
  Eigen::Array<value_type, -1, -1> moments1(D, measurements + 1);
  Eigen::Array<value_type, -1, -1> moments2(D * D, measurements + 1);
  Eigen::Array<T, -1, 1> return_amplitudes(measurements + 1);
  Eigen::Array<value_type, -1, -1> transmission_weights(3, measurements + 1);
  const std::vector<LocalProbe> local_probes =
    build_local_probe_indices<D>(r, probe_lattice_coords);
  Eigen::Array<T, -1, -1>
    local_propagators(local_probes.size(), measurements + 1);

  h.generate_twists();
  h.generate_disorder();
  std::
    fill(h.U_Anderson.begin(), h.U_Anderson.end(), static_cast<value_type>(0));
  INITIATE_PACKET;
  Eigen::Array<T, -1, 1> initial_state = phi.v.col(0);
  if (num_spectral_moments > 0)
    spectral_moments.col(0) = spectral_function<
      T, D, value_type>(r, *this, num_spectral_moments, initial_state);
  h.generate_disorder();
  prepare_leads(r, h, sample_start, sample_L);

  auto update_observables = [&](unsigned m) {
    const auto [m1, m2] = pos_moments<T, D, value_type>(r, phi.v.col(0));
    moments1.col(m) = m1;
    moments2.col(m) = m2;
    return_amplitudes(m) = no_ghost_dot<T, D>(r, initial_state, phi.v.col(0));
    transmission_weights.col(m) = transmission_weight<
      T, D, value_type>(r, sample_start, sample_L, phi.v.col(0));

    for (std::size_t i = 0; i < local_probes.size(); ++i)
      local_propagators(i, m) = phi.v(local_probes[i].second, 0);
  };

  update_observables(0);

  for (unsigned i = 0; i < measurements; ++i) {
    for (unsigned j = 0; j < divisions; ++j) {
      ket.setZero();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      phi.v.col(0) = ket;
      phi.set_index(0);
      phi.Exchange_Boundaries();
    }
    update_observables(i + 1);
  }

  std::
    fill(h.U_Anderson.begin(), h.U_Anderson.end(), static_cast<value_type>(0));
  if (num_spectral_moments > 0)
    spectral_moments.col(1) = spectral_function<
      T, D, value_type>(r, *this, num_spectral_moments, phi.v.col(0));

  store_localized_wavepacket(
    spectral_moments, moments1, moments2, return_amplitudes,
    transmission_weights, local_propagators, local_probes, num_global_probes
  );
}

template <typename T, unsigned D>
void Simulation<T, D>::store_localized_wavepacket(
  const Eigen::Array<value_type, -1, -1> &spectral_moments,
  const Eigen::Array<value_type, -1, -1> &moments1,
  const Eigen::Array<value_type, -1, -1> &moments2,
  const Eigen::Array<T, -1, 1> &return_amplitudes,
  const Eigen::Array<value_type, -1, -1> &transmission_weights,
  const Eigen::Array<T, -1, -1> &propagator_amplitudes,
  const std::vector<LocalProbe> &propagator_coords,
  const std::size_t num_global_probes
)
  requires Complex<T>
{
  debug_message("Entered store_localized_wavepacket\n");
  const std::size_t num_local_probes{propagator_coords.size()};
#pragma omp barrier
#pragma omp master
  {
    Global.results_1.resize(D, moments1.cols());
    Global.results_2.resize(D * D, moments2.cols());
    Global.results_3.resize(return_amplitudes.rows(), 1);
    Global.results_4.resize(num_global_probes, propagator_amplitudes.cols());
    Global.results_5.resize(spectral_moments.rows(), spectral_moments.cols());
    Global.results_6
      .resize(transmission_weights.rows(), transmission_weights.cols());
    Global.results_1.setZero();
    Global.results_2.setZero();
    Global.results_3.setZero();
    Global.results_4.setZero();
    Global.results_5.setZero();
    Global.results_6.setZero();
  }
#pragma omp barrier
#pragma omp critical
  {
    Global.results_1 += moments1.template cast<T>();
    Global.results_2 += moments2.template cast<T>();
    Global.results_3 += return_amplitudes;
    Global.results_5 += spectral_moments.template cast<T>();
    Global.results_6 += transmission_weights.template cast<T>();

    for (std::size_t p = 0; p < num_local_probes; ++p) {
      Global.results_4.row(propagator_coords[p].first) =
        propagator_amplitudes.row(p);
    }
  }
#pragma omp barrier
#pragma omp master
  {
    for (unsigned i = 0; i < moments1.cols(); ++i) {
      Eigen::Array<T, D, D> outer =
        (Global.results_1.col(i).matrix() *
         Global.results_1.col(i).matrix().transpose())
          .array();

      Global.results_2.col(i) -=
        outer.template reshaped<Eigen::ColMajor>(D * D, 1).template cast<T>();
    }
    Global.results_3 = Global.results_3.cwiseAbs2();
  }
#pragma omp barrier
#pragma omp master
  {
    Eigen::Array<value_type, -1, -1> results_1_real = Global.results_1.real();
    Eigen::Array<value_type, -1, -1> results_2_real = Global.results_2.real();
    Eigen::Array<value_type, -1, -1> results_3_real = Global.results_3.real();
    Eigen::Array<value_type, -1, -1> results_5_real = Global.results_5.real();
    Eigen::Array<value_type, -1, -1> results_6_real = Global.results_6.real();
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    const std::string name2("Calculation/localized_wave_packet/StatesMean");
    const std::string name3(
      "Calculation/localized_wave_packet/StatesCovariance"
    );
    const std::string name4(
      "Calculation/localized_wave_packet/ReturnProbability"
    );
    const std::string name5(
      "Calculation/localized_wave_packet/PropagatorAmplitude"
    );
    const std::string name6(
      "Calculation/localized_wave_packet/SpectralMoments"
    );
    const std::string name7(
      "Calculation/localized_wave_packet/TransmissionWeights"
    );
    write_hdf5(results_1_real, file, name2);
    write_hdf5(results_2_real, file, name3);
    write_hdf5(results_3_real, file, name4);
    if (num_global_probes > 0)
      write_hdf5(Global.results_4, file, name5);
    if (spectral_moments.rows() > 0)
      write_hdf5(results_5_real, file, name6);
    write_hdf5(results_6_real, file, name7);
    file->close();
    delete file;
  }
#pragma omp barrier
  debug_message("Left localized_wavepacket\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
