#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
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

template <typename T, unsigned D>
void Simulation<T, D>::calc_swave()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_swave\n");
  std::string base_grp = "/Calculation/s_wave/";
  std::string tmp = base_grp + "NumRandoms";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    Global.calculate_s_wave = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, &file, tmp);
      Global.calculate_s_wave = true;
    } catch (H5::Exception &e) {
      debug_message("s_wave: no need to calculate it.\n");
    }
    file.close();
  }
#pragma omp barrier
  bool local_calculate_s_wave = false;
#pragma omp critical
  local_calculate_s_wave = Global.calculate_s_wave;
#pragma omp barrier
  if (local_calculate_s_wave) {
#pragma omp master
    std::cout << "Calculating SWave.\n";
#pragma omp barrier
    int randoms;
    int num_itr;
    int prv_itr;
    value_type beta;
    value_type mu;
    value_type u;
    value_type weight_r;
    value_type weight_alpha;
    Eigen::Array<value_type, -1, 1> gamma(r.Orb);
    Eigen::Array<value_type, -1, 1> delta(r.Orb);
#pragma omp critical
    {
      H5::H5File file(name, H5F_ACC_RDONLY);
      std::string path = base_grp + "NumRandoms";
      get_hdf5<int>(&randoms, &file, path);
      path = base_grp + "NumIterations";
      get_hdf5<int>(&num_itr, &file, path);
      path = base_grp + "PrevIterations";
      get_hdf5<int>(&prv_itr, &file, path);
      path = base_grp + "Beta";
      get_hdf5<value_type>(&beta, &file, path);
      path = base_grp + "ChemicalPotential";
      get_hdf5<value_type>(&mu, &file, path);
      path = base_grp + "U";
      get_hdf5<value_type>(&u, &file, path);
      path = base_grp + "Gamma";
      get_hdf5<value_type>(gamma.data(), &file, path);
      path = base_grp + "Delta";
      get_hdf5<value_type>(delta.data(), &file, path);
      path = base_grp + "Wr";
      get_hdf5<value_type>(&weight_r, &file, path);
      path = base_grp + "Walpha";
      get_hdf5<value_type>(&weight_alpha, &file, path);
      file.close();
    }
    s_wave(
      randoms, num_itr, prv_itr, beta, mu, u, weight_r, weight_alpha, gamma,
      delta
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::s_wave(
  const int randoms_,
  const int num_itr_,
  const int prv_itr_,
  const value_type beta_,
  const value_type mu_,
  const value_type u_,
  const value_type weight_r_,
  const value_type weight_alpha_,
  const Eigen::Array<value_type, -1, 1> &gamma_,
  const Eigen::Array<value_type, -1, 1> &delta_
)
  requires Complex<T>
{
  debug_message("Entered SWave\n");
  value_type energy_scale;
#pragma omp critical
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    get_hdf5<value_type>(&energy_scale, &file, (char *)"/EnergyScale");
    file.close();
  }
#pragma omp barrier
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
  const value_type beta = beta_ * energy_scale;
  const value_type mu = mu_ / energy_scale;
  const value_type u = u_ / energy_scale;
  const Eigen::Array<value_type, -1, 1> gamma = gamma_ / energy_scale;
  const Eigen::Array<value_type, -1, 1> delta = delta_ / energy_scale;

  const value_type size = r.Sizet - r.SizetVacancies;

  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(beta, 0.0);
  // const Eigen::Array<value_type, -1, 1> cfree =
  //   Coefficients::build_free_energy<value_type>(beta, 0.0) * energy_scale;

  h.bdg.set_chemical_potential(mu);
  h.bdg.init_sw(delta, gamma);

  Eigen::Array<value_type, -1, 1> sum_gamma(r.Sized);
  sum_gamma.setZero();
  Eigen::Array<T, -1, 1> sum_delta(r.Sized);
  sum_delta.setZero();

  // h.bdg.free_energy.resize(num_itr_);
  // h.bdg.free_energy.setZero();

  value_type weight_avg = 1.0;
  value_type weight_sum = 0.0;
  for (unsigned k = 1; k <= prv_itr_; ++k) {
    weight_avg *= 1.0 + weight_r_ / std::pow(k, weight_alpha_);
    weight_sum += weight_avg;
  }
  if (prv_itr_ > 0) {
    sum_delta = weight_sum * h.bdg.s_delta;
    sum_gamma = weight_sum * h.bdg.hartree;
  }
  h.generate_disorder();
  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(2 * r.Sized);
  Eigen::Array<T, -1, 1> bra(2 * r.Sized);

  Eigen::Array<value_type, -1, 1> results_gamma(r.Sized);
  Eigen::Array<T, -1, 1> results_delta(r.Sized);

  for (unsigned itr = prv_itr_ + 1; itr <= prv_itr_ + num_itr_; ++itr) {
    results_gamma.setZero();
    results_delta.setZero();
    // value_type avr_free_tr = 0;
    for (int vec = 0; vec < randoms_; ++vec) {
      const value_type weight = 1.0 / (vec + 1);
      h.generate_twists();
      phi.initiate_phases();
      phi.set_index(0);
      phi.initiate_vector();
      phi.v.col(0) *= std::sqrt(size);

      // T free_en = 0;
      bra = phi.v.col(0);
      // Iteration for Free Energy - Will need if constexpr
      // phi.Exchange_Boundaries();
      // for (unsigned n = 0, N = cfree.size(); n < N; ++n) {
      //   phi.cheb_iteration(n);
      //   const T tmp =
      //     (bra.conjugate() * phi.v.col(phi.get_index()).array()).sum();
      //   free_en += cfree(n) * tmp;
      // }
      // avr_free_tr += weight * (free_en.real() - avr_free_tr);

      // Iteration for Diagonal Elements - Will need if constexpr
      ket.setZero();
      phi.v.setZero();
      phi.set_index(0);
      phi.v.col(0) = bra.matrix();
      phi.Exchange_Boundaries();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      const Eigen::Array<value_type, -1, 1> map_gamma =
        u * (bra.conjugate() * ket).abs2().head(r.Sized);
      results_gamma += weight * (map_gamma - results_gamma);

      // gamma = 1
      phi.v.setZero();
      phi.set_index(0);
      phi.v.col(0) = bra.matrix();
      phi.template pairing<-1>(1.0, phi.v.col(0));

      ket.setZero();
      phi.Exchange_Boundaries();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      phi.template pairing<1>(1.0, ket);

      const Eigen::Array<value_type, -1, 1> upsilon_1 =
        (bra.conjugate() * ket).abs2().head(r.Sized) -
        (bra.conjugate() * ket).abs2().tail(r.Sized);
      // gamma = i
      phi.v.setZero();
      phi.set_index(0);
      phi.v.col(0) = bra.matrix();
      phi.template pairing<-1>(T(0.0, 1.0), phi.v.col(0));

      ket.setZero();
      phi.Exchange_Boundaries();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      phi.template pairing<1>(T(0.0, 1.0), ket);

      const Eigen::Array<value_type, -1, 1> upsilon_i =
        (bra.conjugate() * ket).abs2().head(r.Sized) -
        (bra.conjugate() * ket).abs2().tail(r.Sized);

      Eigen::Array<T, -1, 1> map_delta(r.Sized);
      map_delta.real() = 0.5 * u * upsilon_1;
      map_delta.imag() = -0.5 * u * upsilon_i;
      results_delta += weight * (map_delta - results_delta);
    }
    weight_avg *= 1.0 + weight_r_ / std::pow(itr, weight_alpha_);
    weight_sum += weight_avg;
    sum_delta += weight_avg * results_delta;

    h.bdg.s_delta = sum_delta / weight_sum;
    sum_gamma += weight_avg * results_gamma;
    h.bdg.hartree = sum_gamma / weight_sum;
    h.bdg.update_onsite();

    // Needs if constexpr for grand canonical potential
    // value_type background = 0.5 * avr_free_tr;
    // std::array<unsigned, D> idx;
    // std::array<unsigned, D> start;
    // std::array<unsigned, D> final;
    // for (unsigned d = 0; d < D; ++d) {
    //   start[d] = NGHOSTS;
    //   final[d] = r.Ld[D - 1 - d] - NGHOSTS;
    // }
    // for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    //   auto body = [&](const std::array<unsigned, D> &i) {
    //     if constexpr (D == 2)
    //       local.set({i[1], i[0], io});
    //     else if constexpr (D == 3)
    //       local.set({i[2], i[1], i[0], io});
    //     r.convertCoordinates(global, local);
    //     const unsigned glob = global.index;
    //     const unsigned locl = local.index;
    //     const value_type dd =
    //       std::norm(h.bdg.s_delta(locl)) * energy_scale * energy_scale;
    //     const value_type tr_h = (h.bdg.onsite(locl) + h.bdg.mu) * energy_scale;
    //     const value_type e_c = (dd + tr_h * tr_h) / u_;
    //     background += tr_h - e_c;
    //   };
    //   UnitCellLoop<D>::run(idx, start, final, body);
    // }
    // h.bdg.free_energy(itr - prev_itr_ - 1) = background;
  }
  store_s_wave(energy_scale, u_, prv_itr_ + num_itr_);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_s_wave(
  const value_type energy_scale_,
  const value_type u_,
  const unsigned total_steps_
)
  requires Complex<T>
{
  debug_message("Entered store_swave\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  {
    Global.s_wave_map.resize(r.Sizet, 1);
    Global.free_energy.resize(h.bdg.free_energy.size(), 1);
    Global.free_energy.setZero();
  }
#pragma omp barrier
  std::array<unsigned, D> idx;
  std::array<unsigned, D> start;
  std::array<unsigned, D> final;
  for (unsigned d = 0; d < D; ++d) {
    start[d] = NGHOSTS;
    final[d] = r.Ld[D - 1 - d] - NGHOSTS;
  }
  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2)
        local.set({i[1], i[0], io});
      else if constexpr (D == 3)
        local.set({i[2], i[1], i[0], io});
      r.convertCoordinates(global, local);
      const unsigned glob = global.index;
      const unsigned locl = local.index;
      const value_type delta_i = std::real(h.bdg.s_delta(locl)) * energy_scale_;
      Global.s_wave_map(glob, 0) = delta_i;
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/s_wave/";
    std::string ng = base_grp + "Map";
    write_hdf5(Global.s_wave_map, &file, ng);
  }
#pragma omp barrier
  debug_message("Left store_swave\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
