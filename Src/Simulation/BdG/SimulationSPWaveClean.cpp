#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template <typename T, unsigned D> class Hamiltonian;
template <typename T, unsigned D> class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
#include "Loop.hpp"
#include "Coefficients.hpp"
#include "SigWeight.hpp"

template <typename T, unsigned D>
void Simulation<T, D>::calc_spwave_clean()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_spwave_clean\n");
  std::string base_grp = "/Calculation/sp_wave_c/";
  std::string tmp = base_grp + "NumRandoms";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    Global.calculate_sp_wave = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, &file, tmp);
      Global.calculate_sp_wave = true;
    } catch (H5::Exception &e) {
      debug_message("sp_wave: no need to calculate it.\n");
    }
    file.close();
  }
#pragma omp barrier
  bool local_calculate_sp_wave = false;
#pragma omp critical
  local_calculate_sp_wave = Global.calculate_sp_wave;
#pragma omp barrier
  if (local_calculate_sp_wave) {
#pragma omp master
    std::cout << "Calculating SPWave.\n";
#pragma omp barrier
    int randoms;
    int num_itr;
    int prv_itr;
    value_type N0;
    value_type tau;
    value_type u_init = 0.0;
    int avg_bonds = 0;
#pragma omp critical
    {
      H5::Exception::dontPrint();
      H5::H5File file(name, H5F_ACC_RDONLY);
      std::string path = base_grp + "NumRandoms";
      get_hdf5<int>(&randoms, &file, path);
      path = base_grp + "NumIterations";
      get_hdf5<int>(&num_itr, &file, path);
      path = base_grp + "PrevIterations";
      get_hdf5<int>(&prv_itr, &file, path);
      path = base_grp + "N0";
      get_hdf5<value_type>(&N0, &file, path);
      path = base_grp + "Tau";
      get_hdf5<value_type>(&tau, &file, path);
      path = base_grp + "InitDamping";
      try {
        get_hdf5<value_type>(&u_init, &file, path);
      } catch (H5::Exception &e) {
      }
      path = base_grp + "AverageBonds";
      try {
        get_hdf5<int>(&avg_bonds, &file, path);
      } catch (H5::Exception &e) {
      }
      file.close();
    }
#pragma omp master
    std::cout << "sp_wave_clean: bond average " << (avg_bonds ? "on" : "off")
              << "\n";
#pragma omp barrier
    sp_wave_clean(randoms, num_itr, prv_itr, N0, tau, u_init, avg_bonds != 0);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::sp_wave_clean(
  const int randoms_,
  const int num_itr_,
  const int prv_itr_,
  const value_type N0_,
  const value_type tau_,
  const value_type u_init_,
  const bool avg_bonds_
)
  requires Complex<T>
{
  debug_message("Entered SPWaveClean\n");

  if constexpr (!(pairing::is_s_wave && pairing::is_p_wave)) {
    (void)randoms_;
    (void)num_itr_;
    (void)prv_itr_;
    (void)N0_;
    (void)tau_;
    (void)u_init_;
    (void)avg_bonds_;
#pragma omp master
    std::cerr << "sp_wave_clean: this solver updates the on-site and the bond "
                 "channel together, so the build must carry both. Build with "
                 "PAIRING=3.\n";
#pragma omp barrier
    exit(1);
  } else {
    const value_type energy_scale = h.bdg.energy_scale;
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type n_cells = r.Nt;
    const unsigned max_bonds = h.pr.max_bonds;
    const unsigned n_slots = max_bonds * r.Orb;

    const Eigen::Array<value_type, -1, 1> coefs =
      Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);

    Eigen::Array<T, -1, 1> mean_s(r.Orb);
    Eigen::Array<T, -1, 1> map_s(r.Orb);
    Eigen::Array<T, -1, 1> local_s(r.Orb);
    Eigen::Array<T, -1, 1> per_orb(r.Orb);

    Eigen::Array<T, -1, -1> mean_b(max_bonds, r.Orb);
    Eigen::Array<T, -1, -1> map_b(max_bonds, r.Orb);
    Eigen::Array<T, -1, -1> local_b(max_bonds, r.Orb);
    Eigen::Array<T, -1, 1> per_orb_p(r.Orb);
    Eigen::Array<T, -1, 1> per_orb_h(r.Orb);

    value_type u_weight = u_init_;

    mean_s = h.pr.SDelta0;
    mean_b = h.pr.Delta0;
    h.pr.symmetrize_bonds(mean_b);
    if (avg_bonds_)
      h.pr.average_bonds(mean_b);
    h.pr.broadcast_s(mean_s, h.bdg.s_delta);
    h.pr.broadcast(mean_b, h.bdg.nn_delta);

    const auto record = [&](const unsigned row) {
      for (unsigned io = 0; io < r.Orb; ++io) {
        Global.s_delta_hist(row, io) = mean_s(io) * energy_scale;
        for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b)
          Global.nn_delta_hist(row, b + io * max_bonds) =
            mean_b(b, io) * energy_scale;
      }
    };
#pragma omp master
    {
      Global.orb_sum.resize(r.Orb);
      Global.nn_sum.resize(n_slots);
      Global.s_delta_hist.resize(num_itr_ + 1, r.Orb);
      Global.s_delta_hist.setZero();
      Global.nn_delta_hist.resize(num_itr_ + 1, n_slots);
      Global.nn_delta_hist.setZero();
      record(0);
    }
#pragma omp barrier
    h.generate_disorder();
    KPM_Vector<T, D> phi(2, *this);
    Eigen::Array<T, -1, 1> ket(2 * r.Sized);
    Eigen::Array<T, -1, 1> ket_ref(2 * r.Sized);

    for (unsigned itr = prv_itr_ + 1; itr <= prv_itr_ + num_itr_; ++itr) {
      local_s.setZero();
      local_b.setZero();

      for (int vec = 0; vec < randoms_; ++vec) {
        const value_type weight = 1.0 / (vec + 1);
        h.generate_twists();
        phi.initiate_phases();
        phi.set_index(0);
        phi.initiate_vector();
        phi.v.col(0) *= std::sqrt(size);
        ket.setZero();
        phi.Exchange_Boundaries();
        for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
          phi.cheb_iteration(n);
          ket += coefs(n) * phi.v.col(phi.get_index()).array();
        }
        ket_ref = ket;
        {
          phi.template pairing<1>(1.0, ket);
          const Eigen::Array<value_type, -1, 1> upsilon =
            ket.abs2().head(r.Sized) - ket.abs2().tail(r.Sized);
          const Eigen::Array<T, -1, 1> f_s = 0.5 * h.pr.U * upsilon;
          h.pr.orbital_sum(f_s, per_orb);
          local_s += weight * (per_orb - local_s);
        }
        for (unsigned io = 0; io < r.Orb; ++io)
          for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b) {
            ket = ket_ref;
            phi.template nn_pairing<1>(1.0, io, b, ket);
            const Eigen::Array<T, -1, 1> pp =
              ket.abs2().head(r.Sized).template cast<T>();
            const Eigen::Array<T, -1, 1> hh =
              ket.abs2().tail(r.Sized).template cast<T>();
            h.pr.orbital_sum(pp, per_orb_p);
            h.pr.orbital_sum(hh, per_orb_h);

            const unsigned jo = h.pr.target_orb(b, io);
            const value_type tmp_b = 0.5 * h.pr.V;
            const T f_b = tmp_b * (per_orb_p(io) - per_orb_h(jo));
            local_b(b, io) += weight * (f_b - local_b(b, io));
          }
      }
#pragma omp barrier
#pragma omp master
      {
        Global.orb_sum.setZero();
        Global.nn_sum.setZero();
      }
#pragma omp barrier
#pragma omp critical
      {
        for (unsigned io = 0; io < r.Orb; ++io) {
          Global.orb_sum(io) += local_s(io);
          for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b)
            Global.nn_sum(b + io * max_bonds) += local_b(b, io);
        }
      }
#pragma omp barrier
      map_b.setZero();
      for (unsigned io = 0; io < r.Orb; ++io) {
        map_s(io) = Global.orb_sum(io) / n_cells;
        for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b)
          map_b(b, io) = Global.nn_sum(b + io * max_bonds) / n_cells;
      }
      h.pr.symmetrize_bonds(map_b);
      if (avg_bonds_)
        h.pr.average_bonds(map_b);

      u_weight = 1.0 + u_weight / sig_weight_ratio<value_type>(itr, N0_, tau_);
      const value_type gamma_n = 1.0 / u_weight;
      mean_s += gamma_n * (map_s - mean_s);
      mean_b += gamma_n * (map_b - mean_b);
      h.pr.broadcast_s(mean_s, h.bdg.s_delta);
      h.pr.broadcast(mean_b, h.bdg.nn_delta);
#pragma omp barrier
#pragma omp master
      record(itr - prv_itr_);
#pragma omp barrier
    }
    store_sp_wave_clean(prv_itr_ + num_itr_, mean_s, mean_b, u_weight);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_sp_wave_clean(
  const unsigned total_steps_,
  const Eigen::Array<T, -1, 1> &mean_s_,
  const Eigen::Array<T, -1, -1> &mean_b_,
  const value_type u_weight_
)
  requires Complex<T>
{
  debug_message("Entered store_sp_wave_clean\n");
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/sp_wave_c/";
    std::string ng;

    Eigen::Array<value_type, -1, -1> hist_s = Global.s_delta_hist.real();
    ng = base_grp + "HistS";
    write_hdf5(hist_s, &file, ng);

    Eigen::Array<value_type, -1, -1> hist_b = Global.nn_delta_hist.real();
    ng = base_grp + "HistB";
    write_hdf5(hist_b, &file, ng);

    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);

    Eigen::Array<value_type, -1, -1> s_ri(2, r.Orb);
    for (unsigned io = 0; io < r.Orb; ++io) {
      s_ri(0, io) = mean_s_(io).real();
      s_ri(1, io) = mean_s_(io).imag();
    }
    ng = base_grp + "DeltaS";
    write_hdf5(s_ri, &file, ng);

    const unsigned max_bonds = h.pr.max_bonds;
    const unsigned n_slots = max_bonds * r.Orb;
    Eigen::Array<value_type, -1, -1> b_ri =
      Eigen::Array<value_type, -1, -1>::Zero(2, n_slots);
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b) {
        const std::size_t k = b + io * max_bonds;
        b_ri(0, k) = mean_b_(b, io).real();
        b_ri(1, k) = mean_b_(b, io).imag();
      }
    ng = base_grp + "DeltaB";
    write_hdf5(b_ri, &file, ng);

    Eigen::Array<value_type, -1, -1> damping(1, 1);
    damping(0, 0) = u_weight_;
    ng = base_grp + "Damping";
    write_hdf5(damping, &file, ng);
  }
#pragma omp barrier
  debug_message("Left store_sp_wave_clean\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
