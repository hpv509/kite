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
    s_wave_clean(
      randoms, num_itr, prv_itr, beta, mu, u, weight_r, weight_alpha, gamma,
      delta
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::s_wave_clean(
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
  Coordinates<std::size_t, D + 1> local(r.Ld);
  const value_type beta = beta_ * energy_scale;
  const value_type mu = mu_ / energy_scale;
  const value_type u = u_ / energy_scale;
  const Eigen::Array<value_type, -1, 1> gamma = gamma_ / energy_scale;
  const Eigen::Array<value_type, -1, 1> delta = delta_ / energy_scale;

  const value_type size = r.Sizet - r.SizetVacancies;

  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(beta, 0.0);

  h.bdg.set_chemical_potential(mu);
  h.bdg.init_fields(delta, gamma);

  Eigen::Array<T, -1, 1> sum_delta(r.Sized);
  sum_delta.setZero();
#pragma omp master
  {
    Global.orb_sum.resize(r.Orb);
    Global.s_delta_hist.resize(num_itr_ + 1, r.Orb);
    Global.s_delta_hist.setZero();
  }
#pragma omp barrier
  value_type weight_avg = 1.0;
  value_type weight_sum = 0.0;
  for (unsigned k = 1; k <= prv_itr_; ++k) {
    weight_avg *= 1.0 + weight_r_ / std::pow(k, weight_alpha_);
    weight_sum += weight_avg;
  }
  if (prv_itr_ > 0)
    sum_delta = weight_sum * h.bdg.s_delta;
  h.generate_disorder();
  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(2 * r.Sized);
  Eigen::Array<T, -1, 1> bra(2 * r.Sized);

  Eigen::Array<T, -1, 1> results_delta(r.Sized);
#pragma omp master
  for (unsigned io = 0; io < r.Orb; ++io) {
    local.set({r.Ld[0] / 2, r.Ld[1] / 2, io});
    Global.s_delta_hist(0, io) = h.bdg.s_delta(local.index) * energy_scale;
  }
#pragma omp barrier
  for (unsigned itr = prv_itr_ + 1; itr <= prv_itr_ + num_itr_; ++itr) {
    results_delta.setZero();
    for (int vec = 0; vec < randoms_; ++vec) {
      const value_type weight = 1.0 / (vec + 1);
      h.generate_twists();
      phi.initiate_phases();
      phi.set_index(0);
      phi.initiate_vector();
      phi.v.col(0) *= std::sqrt(size);

      phi.empty_ghosts(0); // use this for the averaging
      bra = phi.v.col(0);

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

      Eigen::Array<T, -1, 1> map_delta = 0.5 * u * upsilon_1;
#pragma omp barrier
#pragma omp master
      Global.orb_sum.setZero();
#pragma omp barrier
#pragma omp critical
      {
        for (unsigned io = 0; io < r.Orb; ++io)
          Global.orb_sum(io) +=
            map_delta.segment(io * r.Nd, r.Nd).sum() / static_cast<T>(r.Nt);
      }
#pragma omp barrier
      for (unsigned io = 0; io < r.Orb; ++io) {
        const T avr_orb = Global.orb_sum(io);
        map_delta.segment(io * r.Nd, r.Nd).setConstant(avr_orb);
      }
#pragma omp barrier
      results_delta += weight * (map_delta - results_delta);
    }
    weight_avg *= 1.0 + weight_r_ / std::pow(itr, weight_alpha_);
    weight_sum += weight_avg;
    sum_delta += weight_avg * results_delta;
    h.bdg.s_delta = sum_delta / weight_sum;
    h.bdg.update_onsite();
#pragma omp master
    for (unsigned io = 0; io < r.Orb; ++io) {
      local.set({r.Ld[0] / 2, r.Ld[1] / 2, io});
      Global.s_delta_hist(itr - prv_itr_, io) =
        h.bdg.s_delta(local.index) * energy_scale;
    }
#pragma omp barrier
  }
  store_s_wave_clean(prv_itr_ + num_itr_);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_s_wave_clean(const unsigned total_steps_)
  requires Complex<T>
{
  debug_message("Entered store_swave\n");
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
  }
#pragma omp barrier
  debug_message("Left store_swave\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
