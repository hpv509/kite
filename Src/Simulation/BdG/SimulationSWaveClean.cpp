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
void Simulation<T, D>::calc_swave_clean()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_swave_clean\n");
  std::string base_grp = "/Calculation/s_wave_c/";
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
    int randoms, num_itr, prv_itr;
    value_type N0, tau;
    value_type u_init = 0.0;
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
      file.close();
    }
    s_wave_clean(randoms, num_itr, prv_itr, N0, tau, u_init);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::s_wave_clean(
  const int randoms_,
  const int num_itr_,
  const int prv_itr_,
  const value_type N0_,
  const value_type tau_,
  const value_type u_init_
)
  requires Complex<T>
{
  debug_message("Entered SWaveClean\n");
  const value_type energy_scale = h.bdg.energy_scale;
  const value_type size = r.Sizet - r.SizetVacancies;
  const value_type n_cells = r.Nt;

  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);

  Eigen::Array<T, -1, 1> mean_delta(r.Orb);
  Eigen::Array<T, -1, 1> map_delta_orb(r.Orb);
  Eigen::Array<T, -1, 1> local_delta(r.Orb);
  Eigen::Array<T, -1, 1> per_orb(r.Orb);

  value_type u_weight = u_init_;
  mean_delta = h.pr.SDelta0;
  h.pr.broadcast_s(mean_delta, h.bdg.s_delta);
#pragma omp master
  {
    Global.orb_sum.resize(r.Orb);
    Global.s_delta_hist.resize(num_itr_ + 1, r.Orb);
    Global.s_delta_hist.setZero();
    for (unsigned io = 0; io < r.Orb; ++io)
      Global.s_delta_hist(0, io) = mean_delta(io) * energy_scale;
  }
#pragma omp barrier
  h.generate_disorder();
  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(2 * r.Sized);

  for (unsigned itr = prv_itr_ + 1; itr <= prv_itr_ + num_itr_; ++itr) {
    local_delta.setZero();
    for (int vec = 0; vec < randoms_; ++vec) {
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
      phi.template pairing<1>(1.0, ket);
      const Eigen::Array<value_type, -1, 1> upsilon =
        ket.abs2().head(r.Sized) - ket.abs2().tail(r.Sized);
      const Eigen::Array<T, -1, 1> map_delta = 0.5 * h.pr.U * upsilon;

      h.pr.orbital_sum(map_delta, per_orb);
      const value_type weight = 1.0 / (vec + 1);
      local_delta += weight * (per_orb - local_delta);
    }
#pragma omp barrier
#pragma omp master
    Global.orb_sum.setZero();
#pragma omp barrier
#pragma omp critical
    {
      for (unsigned io = 0; io < r.Orb; ++io)
        Global.orb_sum(io) += local_delta(io);
    }
#pragma omp barrier
    for (unsigned io = 0; io < r.Orb; ++io)
      map_delta_orb(io) = Global.orb_sum(io) / n_cells;

    u_weight = 1.0 + u_weight / sig_weight_ratio<value_type>(itr, N0_, tau_);
    const value_type gamma_n = 1.0 / u_weight;
    mean_delta += gamma_n * (map_delta_orb - mean_delta);
    h.pr.broadcast_s(mean_delta, h.bdg.s_delta);
#pragma omp barrier
#pragma omp master
    for (unsigned io = 0; io < r.Orb; ++io)
      Global.s_delta_hist(itr - prv_itr_, io) = mean_delta(io) * energy_scale;
#pragma omp barrier
  }
  store_s_wave_clean(prv_itr_ + num_itr_, mean_delta, u_weight);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_s_wave_clean(
  const unsigned total_steps_,
  const Eigen::Array<T, -1, 1> &mean_delta_,
  const value_type u_weight_
)
  requires Complex<T>
{
  debug_message("Entered store_s_wave_clean\n");
#pragma omp master
  {
    const Eigen::Array<value_type, -1, -1> s_wave_r =
      Global.s_delta_hist.real();
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/s_wave_c/";
    std::string ng = base_grp + "Hist";
    write_hdf5(s_wave_r, &file, ng);

    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);

    Eigen::Array<value_type, -1, -1> delta_ri(2, r.Orb);
    for (unsigned io = 0; io < r.Orb; ++io) {
      delta_ri(0, io) = mean_delta_(io).real();
      delta_ri(1, io) = mean_delta_(io).imag();
    }
    ng = base_grp + "Delta";
    write_hdf5(delta_ri, &file, ng);

    Eigen::Array<value_type, -1, -1> damping(1, 1);
    damping(0, 0) = u_weight_;
    ng = base_grp + "Damping";
    write_hdf5(damping, &file, ng);
  }
#pragma omp barrier
  debug_message("Left store_s_wave_clean\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
