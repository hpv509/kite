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

constexpr unsigned SWAVE_N_SYM = 16;

template <typename T> T sw_weight_ratio(const T n, const T N0, const T tau)
{
  if (std::isinf(N0))
    return std::exp(T(1.0) / tau);
  auto softplus = [](const T z) {
    return std::max(z, T(0.0)) + std::log1p(std::exp(-std::abs(z)));
  };
  return std::exp(softplus((N0 - n + T(1.0)) / tau) - softplus((N0 - n) / tau));
}

template <unsigned D, typename F>
void for_each_interior(LatticeStructure<D> &r_, F &&f_)
{
  Coordinates<std::size_t, D + 1> Ld(r_.Ld), Lt(r_.Lt);
  for (std::size_t i = 0; i < r_.Sized; ++i) {
    Ld.set_coord(i);
    bool ghost = false;
    for (unsigned k = 0; k < D; ++k)
      ghost |= Ld.coord[k] < NGHOSTS || Ld.coord[k] >= r_.Ld[k] - NGHOSTS;
    if (ghost)
      continue;
    r_.convertCoordinates(Lt, Ld);
    f_(i, Lt.index);
  }
}

template <typename T, typename G>
Eigen::Array<T, -1, 1>
reduce_orbitals(G &global_, const Eigen::Array<T, -1, 1> &per_orb_)
{
#pragma omp barrier
#pragma omp master
  {
    global_.orb_sum.resize(per_orb_.size());
    global_.orb_sum.setZero();
  }
#pragma omp barrier
#pragma omp critical
  global_.orb_sum += per_orb_;
#pragma omp barrier
  return global_.orb_sum;
}

template <typename T, unsigned D>
void erase_vacancies(
  const Vacancy_Operator<T, D> &hV_,
  Eigen::Array<T, -1, 1> &field_
)
{
  for (const auto &pos : hV_.position)
    for (const auto idx : pos)
      field_(idx) = 0.;
}

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
    std::cout << "Calculating SWave (spatially resolved).\n";
#pragma omp barrier
    int randoms, num_itr, prv_itr;
    value_type N0, tau;
    value_type u_init = 0.0;
    bool has_map = false;
    Eigen::Array<value_type, -1, 1> init_map(r.Sizet);
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
      path = base_grp + "InitMap";
      try {
        if (prv_itr > 0) {
          get_hdf5<value_type>(init_map.data(), &file, path);
          has_map = true;
        }
      } catch (H5::Exception &e) {
      }
      file.close();
    }
    h.pr.broadcast_s(h.pr.SDelta0, h.bdg.s_delta);
    if (has_map) {
      const value_type energy_scale = h.bdg.energy_scale;
      for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_) {
        h.bdg.s_delta(locl_) = T(init_map(glob_) / energy_scale);
      });
    }
    s_wave(randoms, num_itr, prv_itr, N0, tau, u_init);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::s_wave(
  const int randoms_,
  const int num_itr_,
  const int prv_itr_,
  const value_type N0_,
  const value_type tau_,
  const value_type u_init_
)
  requires Complex<T>
{
  debug_message("Entered SWave\n");
  const value_type energy_scale = h.bdg.energy_scale;
  const value_type size = r.Sizet - r.SizetVacancies;
  const unsigned Orb = r.Orb;

  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);

  Eigen::Array<T, -1, 1> per_orb(Orb);
  Eigen::Array<T, -1, 1> n_valid(Orb);
  Eigen::Array<T, -1, 1> mean(Orb);
  Eigen::Array<T, -1, 1> local_delta(r.Sized);

  value_type u_weight = u_init_;
#pragma omp master
  {
    Global.s_delta_hist.resize(num_itr_ + 1, Orb);
    Global.s_delta_hist.setZero();
  }
#pragma omp barrier
  h.generate_disorder();

  per_orb.setConstant(r.N);
  for (const auto &pos : h.hV.position)
    for (const auto idx : pos)
      per_orb(idx / r.Nd) -= 1.;
  n_valid = reduce_orbitals(Global, per_orb);

  erase_vacancies(h.hV, h.bdg.s_delta);
  h.pr.orbital_sum(h.bdg.s_delta, per_orb);
  mean = reduce_orbitals(Global, per_orb) / n_valid;
  if (prv_itr_ < SWAVE_N_SYM) {
    h.pr.broadcast_s(mean, h.bdg.s_delta);
    erase_vacancies(h.hV, h.bdg.s_delta);
  }
#pragma omp master
  Global.s_delta_hist.row(0) = (mean * energy_scale).transpose();
#pragma omp barrier
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

      const value_type weight = 1.0 / (vec + 1);
      local_delta += weight * (map_delta - local_delta);
    }
    erase_vacancies(h.hV, local_delta);
    if (itr <= SWAVE_N_SYM) {
      h.pr.orbital_sum(local_delta, per_orb);
      mean = reduce_orbitals(Global, per_orb) / n_valid;
      h.pr.broadcast_s(mean, local_delta);
      erase_vacancies(h.hV, local_delta);
    }
    u_weight = 1.0 + u_weight / sw_weight_ratio<value_type>(itr, N0_, tau_);
    const value_type gamma_n = 1.0 / u_weight;
    h.bdg.s_delta += gamma_n * (local_delta - h.bdg.s_delta);

    h.pr.orbital_sum(h.bdg.s_delta, per_orb);
    mean = reduce_orbitals(Global, per_orb) / n_valid;
#pragma omp master
    Global.s_delta_hist.row(itr - prv_itr_) = (mean * energy_scale).transpose();
#pragma omp barrier
  }
  store_s_wave(prv_itr_ + num_itr_, u_weight);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_s_wave(
  const unsigned total_steps_,
  const value_type u_weight_
)
  requires Complex<T>
{
  debug_message("Entered store_s_wave\n");
  const value_type energy_scale = h.bdg.energy_scale;
#pragma omp master
  {
    Global.s_delta_map.resize(r.Sizet);
    Global.s_delta_map.setZero();
  }
#pragma omp barrier
#pragma omp critical
  for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_) {
    Global.s_delta_map(glob_) = h.bdg.s_delta(locl_);
  });
#pragma omp barrier
#pragma omp master
  {
    const Eigen::Array<value_type, -1, -1> s_wave_r =
      Global.s_delta_hist.real();
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/s_wave/";
    std::string ng = base_grp + "Hist";
    write_hdf5(s_wave_r, &file, ng);

    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);

    Eigen::Array<value_type, -1, -1> delta_map(r.Sizet, 1);
    delta_map.col(0) = Global.s_delta_map.real() * energy_scale;
    ng = base_grp + "Map";
    write_hdf5(delta_map, &file, ng);

    Eigen::Array<value_type, -1, -1> damping(1, 1);
    damping(0, 0) = u_weight_;
    ng = base_grp + "Damping";
    write_hdf5(damping, &file, ng);
  }
#pragma omp barrier
  debug_message("Left store_s_wave\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
