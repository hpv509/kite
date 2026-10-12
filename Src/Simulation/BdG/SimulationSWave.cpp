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
#include "LocalUnitary.hpp"
#include "Loop.hpp"
#include "Coefficients.hpp"
#include "mpi_utils.hpp"

template <unsigned D, typename F>
void sw_for_each_interior(LatticeStructure<D> &r_, F &&f_)
{
  Coordinates<std::size_t, D + 1> Ld(r_.Ld), Lt(r_.Lt);
  for (std::size_t i = 0; i < r_.Sized; ++i)
    {
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
sw_reduce_orbitals_mpi(G &global_, const Eigen::Array<T, -1, 1> &per_orb_)
{
#pragma omp barrier
#pragma omp master
  {
    global_.orb_sum.resize(per_orb_.size());
    global_.orb_sum.setZero();
  }
#pragma omp barrier
#pragma omp critical
  {
    global_.orb_sum += per_orb_;
  }
#pragma omp barrier
#pragma omp master
  {
    kmpi::sum_all(global_.orb_sum);
  }
#pragma omp barrier
  return global_.orb_sum;
}

template <typename HT, unsigned D, typename FieldT>
void sw_erase_vacancies(const Vacancy_Operator<HT, D> &hV_, Eigen::Array<FieldT, -1, 1> &field_)
{
  for (const auto &pos : hV_.position)
    for (const auto idx : pos)
      field_(idx) = FieldT(0);
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_swave()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_swave\n"); 
  const std::string base_grp = "/Calculation/s_wave/";
  std::string tmp = base_grp + "NumRandomsTransient";
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
  if (!local_calculate_s_wave)
    return;
#pragma omp master
  std::cout << "Calculating spatially resolved SWave.\n";
#pragma omp barrier
  int randoms_transient, n_transient, n_stages, stage_itr, depth;
  value_type mixing;
  
  int prv_itr = 0;
  bool has_map = false;
  Eigen::Array<value_type, -1, 1> init_map(r.Sizet);
#pragma omp critical
  {
    H5::Exception::dontPrint();
    H5::H5File file(name, H5F_ACC_RDONLY);
    std::string path = base_grp + "NumRandomsTransient";
    get_hdf5<int>(&randoms_transient, &file, path);
    path = base_grp + "TransientIterations";
    get_hdf5<int>(&n_transient, &file, path);
    path = base_grp + "TailStages";
    get_hdf5<int>(&n_stages, &file, path);
    path = base_grp + "StageIterations";
    get_hdf5<int>(&stage_itr, &file, path);
    path = base_grp + "Depth";
    get_hdf5<int>(&depth, &file, path);
    path = base_grp + "Mixing";
    get_hdf5<value_type>(&mixing, &file, path);
    
    path = base_grp + "PrevIterations";
    try {
      get_hdf5<int>(&prv_itr, &file, path);
    } catch (H5::Exception &e) {
      prv_itr = 0;
    }
    path = base_grp + "InitMap";
    try {
      get_hdf5<value_type>(init_map.data(), &file, path);
      has_map = true;
    } catch (H5::Exception &e) {
      has_map = false;
    }
    file.close();
  }
  const value_type seed_imag = h.pr.SDelta0.imag().abs().maxCoeff();
  if (seed_imag > value_type(1e-12)) {
#pragma omp critical
    {
      std::cerr
	<< "s_wave: the real-space Anderson solver currently supports only "
	"real onsite Delta. Use a real onsite_delta seed.\n";
      exit(1);
    }
  }
  h.pr.SDelta0 = h.pr.SDelta0.real().template cast<T>().eval();
  h.pr.broadcast_s(h.pr.SDelta0, h.bdg.s_delta);
  if (has_map)
    sw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
    {
      h.bdg.s_delta(locl_) = T(init_map(glob_) / h.bdg.energy_scale);
    });
  s_wave(randoms_transient, n_transient, n_stages, stage_itr, depth, mixing, prv_itr);
}

template <typename T, unsigned D>
void Simulation<T, D>::s_wave(
  const int randoms_transient_,
  const int n_transient_,
  const int n_stages_,
  const int stage_itr_,
  const int depth_,
  const value_type mixing_,
  const int prv_itr_
)
  requires Complex<T>
{
  debug_message("Entered spatially resolved SWave\n");
  if constexpr (D != 2 || !pairing::is_s_wave) {
    (void)randoms_transient_;
    (void)n_transient_;
    (void)n_stages_;
    (void)stage_itr_;
    (void)depth_;
    (void)mixing_;
    (void)prv_itr_;
#pragma omp master
    std::cerr
      << (D != 2
	  ? "s_wave: only implemented for D = 2.\n"
	  : "s_wave: this build has no on-site pairing channel.\n");
#pragma omp barrier
    exit(1);
  } else {
    using RealVec = Eigen::Matrix<value_type, -1, 1>;
    using RealMat = Eigen::Matrix<value_type, -1, -1>;
    
    const value_type energy_scale = h.bdg.energy_scale;
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type n_ranks = static_cast<value_type>(kmpi::size());
    
    const int num_itr = n_transient_ + n_stages_ * stage_itr_;
    const unsigned m = depth_ > 0
      ? std::min<unsigned>(static_cast<unsigned>(depth_), static_cast<unsigned>(r.Sizet)) : 0;
    
    const Eigen::Array<value_type, -1, 1> coefs =
      Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);
    
    Eigen::Array<value_type, -1, 1> local_delta(r.Sized);
    Eigen::Array<T, -1, 1> per_orb(r.Orb);
    Eigen::Array<T, -1, 1> n_valid(r.Orb);
    Eigen::Array<T, -1, 1> mean_orb(r.Orb);
    
    struct SWaveState {
      Eigen::Array<value_type, -1, 1> mean;
      Eigen::Array<value_type, -1, 1> map;
      Eigen::Array<value_type, -1, 1> residual;
      Eigen::Array<value_type, -1, -1> hist;
      
      RealMat dX;
      RealMat dF;
      RealVec x_prev;
      RealVec f_prev;
      
      unsigned n_hist = 0;
      value_type res_prev = 0;
    };
    
    SWaveState *state = nullptr;
#pragma omp single copyprivate(state)
    {
      state = new SWaveState;
      
      state->mean = Eigen::Array<value_type, -1, 1>::Zero(r.Sizet);
      state->map = Eigen::Array<value_type, -1, 1>::Zero(r.Sizet);
      state->residual = Eigen::Array<value_type, -1, 1>::Zero(num_itr);
      state->hist = Eigen::Array<value_type, -1, -1>::Zero(num_itr + 1, r.Orb);
      
      state->dX = RealMat::Zero(r.Sizet, m);
      state->dF = RealMat::Zero(r.Sizet, m);
      state->x_prev = RealVec::Zero(r.Sizet);
      state->f_prev = RealVec::Zero(r.Sizet);
    }
    
    constexpr value_type max_log_step = value_type(0.405);
    constexpr value_type amp_min = value_type(1e-12);
    
    h.generate_disorder();
    sw_erase_vacancies(h.hV, h.bdg.s_delta);
    KPM_Vector<T, D> phi(2, *this);
    LocalUnitary<T, D> U(r, phi);
    Eigen::Array<T, -1, 1> ket(2 * r.Sized);
#pragma omp master
    state->mean.setZero();
#pragma omp barrier
#pragma omp critical
    {
      sw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
      {
        state->mean(glob_) = std::real(h.bdg.s_delta(locl_));
      });
    }
#pragma omp barrier
#pragma omp master
    {
      kmpi::sum_all(state->mean);
      state->mean /= n_ranks;
    }
#pragma omp barrier
    sw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
    {
      h.bdg.s_delta(locl_) = T(state->mean(glob_));
    });
    sw_erase_vacancies(h.hV, h.bdg.s_delta);
    
    per_orb.setConstant(r.N);
    for (const auto &pos : h.hV.position)
      for (const auto idx : pos)
        per_orb(idx / r.Nd) -= T(1);
    n_valid = sw_reduce_orbitals_mpi(Global, per_orb);
    h.pr.orbital_sum(h.bdg.s_delta, per_orb);
    mean_orb = sw_reduce_orbitals_mpi(Global, per_orb) / n_valid;
#pragma omp master
    state->hist.row(0) = mean_orb.real().transpose() * energy_scale;

    auto exchange_s_delta = [&]()
    {
      phi.v.setZero();
      phi.set_index(0);
      phi.v.col(0).head(r.Sized) = h.bdg.s_delta.matrix();
      phi.Exchange_Boundaries();
      h.bdg.s_delta = phi.v.col(0).head(r.Sized).array();
      h.bdg.s_delta = h.bdg.s_delta.real().template cast<T>().eval();
      sw_erase_vacancies(h.hV, h.bdg.s_delta);
    };
    
#pragma omp barrier
    exchange_s_delta();
#pragma omp barrier
    
    for (int itr = 0; itr < num_itr; ++itr) {
      rnd.init_random(seed_v);
      h.rnd.init_random(seed_h);
      const int stage = itr < n_transient_ ? 0 : 1 + (itr - n_transient_) / stage_itr_;
      const int n_vec = randoms_transient_ << stage;
      const bool fresh = itr == 0 || (itr >= n_transient_ &&
				      (itr - n_transient_) % stage_itr_ == 0);
      local_delta.setZero();
      for (int vec = 0; vec < n_vec; ++vec) {
        h.generate_twists();
        phi.initiate_phases();
        phi.set_index(0);
        phi.initiate_vector();
        phi.v.col(0) *= std::sqrt(size);
        phi.Exchange_Boundaries();
        phi.chebyshev_sum(coefs, ket.data());
        U.template onsite<1>(T(1), ket);
        const Eigen::Array<value_type, -1, 1> upsilon =
          ket.abs2().head(r.Sized) - ket.abs2().tail(r.Sized);
        const Eigen::Array<T, -1, 1> map_delta =
          value_type(0.5) * h.pr.U * upsilon;
        local_delta +=
          (map_delta.real() - local_delta) / value_type(vec + 1);
      }
      sw_erase_vacancies(h.hV, local_delta);
#pragma omp master
      state->map.setZero();
#pragma omp barrier
#pragma omp critical
      {
        sw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_) {
          state->map(glob_) = local_delta(locl_);
        });
      }
#pragma omp barrier
#pragma omp master
      {
        kmpi::sum_all(state->map);
        state->map /= n_ranks;
      }
#pragma omp barrier
#pragma omp master
      {
        const RealVec f_lin = (state->map - state->mean).matrix();
        state->residual(itr) = f_lin.cwiseAbs().maxCoeff() * energy_scale;
        const Eigen::Array<value_type, -1, 1> amp = state->mean.abs().max(amp_min);
        const Eigen::Array<value_type, -1, 1> amp_map = state->map.abs().max(amp_min);
        const RealVec x = amp.log().matrix();
        const RealVec f = (amp_map.log() - amp.log()).matrix();
        const value_type res = f.cwiseAbs().maxCoeff();
        if (fresh || res > value_type(1.5) * state->res_prev)
          state->n_hist = 0;
        else if (m > 0)
	  {
	    if (state->n_hist == m)
	      {
		for (unsigned j = 0; j + 1 < m; ++j)
		  {
		    state->dX.col(j) = state->dX.col(j + 1);
		    state->dF.col(j) = state->dF.col(j + 1);
		  }
		--state->n_hist;
	      }
	    state->dX.col(state->n_hist) = x - state->x_prev;
	    state->dF.col(state->n_hist) = f - state->f_prev;
	    ++state->n_hist;
	  }
        state->x_prev = x;
        state->f_prev = f;
        state->res_prev = res;
        RealVec step = mixing_ * f;
        if (state->n_hist > 0)
	  {
	    const auto dFh = state->dF.leftCols(state->n_hist);
	    const auto dXh = state->dX.leftCols(state->n_hist);
	    Eigen::CompleteOrthogonalDecomposition<RealMat> cod(dFh);
	    cod.setThreshold(1e-8);
	    const RealVec gamma = cod.solve(f);
	    const RealMat G = dXh + mixing_ * dFh;
	    step -= G * gamma;
	  }
        if (step.dot(f) <= value_type(0)) {
	  step = mixing_ * f;
	  state->n_hist = 0;
	}
        const value_type s_max = step.cwiseAbs().maxCoeff();
        if (s_max > max_log_step)
          step *= max_log_step / s_max;
        const Eigen::Array<value_type, -1, 1> phase = state->map / amp_map;
        state->mean = (x + step).array().exp() * phase;
      }
      
#pragma omp barrier
      sw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_) {
        h.bdg.s_delta(locl_) = T(state->mean(glob_));
      });
      sw_erase_vacancies(h.hV, h.bdg.s_delta);
#pragma omp barrier
      exchange_s_delta();
#pragma omp barrier
      
      h.pr.orbital_sum(h.bdg.s_delta, per_orb);
      mean_orb = sw_reduce_orbitals_mpi(Global, per_orb) / n_valid;
#pragma omp master
      state->hist.row(itr + 1) = mean_orb.real().transpose() * energy_scale;
#pragma omp barrier
    }
    
    store_s_wave(static_cast<unsigned>(prv_itr_ + num_itr), state->hist, state->residual, state->mean);
    
#pragma omp barrier
#pragma omp master
    delete state;
#pragma omp barrier
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_s_wave(
  const unsigned total_steps_,
  const Eigen::Array<value_type, -1, -1> &hist_,
  const Eigen::Array<value_type, -1, 1> &residual_,
  const Eigen::Array<value_type, -1, 1> &map_
)
  requires Complex<T>
{
  debug_message("Entered store_s_wave\n");
  
  const value_type energy_scale = h.bdg.energy_scale;
  
#pragma omp master
  if (kmpi::is_root()) {
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/s_wave/";
    
    std::string ng = base_grp + "Hist";
    write_hdf5(hist_, &file, ng);
    
    const Eigen::Array<value_type, -1, -1> res = residual_;
    ng = base_grp + "Residual";
    write_hdf5(res, &file, ng);
    
    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);
    
    Eigen::Array<value_type, -1, -1> delta_map(r.Sizet, 1);
    delta_map.col(0) = map_ * energy_scale;
    ng = base_grp + "Map";
    write_hdf5(delta_map, &file, ng);
    
    Eigen::Array<value_type, -1, -1> ranks(1, 1);
    ranks(0, 0) = static_cast<value_type>(kmpi::size());
    ng = base_grp + "NumRanks";
    write_hdf5(ranks, &file, ng);
    
    file.close();
  }
  
#pragma omp master
  kmpi::barrier();
#pragma omp barrier
  
  debug_message("Left store_s_wave\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
