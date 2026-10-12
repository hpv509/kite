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
void pw_for_each_interior(LatticeStructure<D> &r_, F &&f_)
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

template <unsigned D, typename F>
void pw_for_each_interior_orbital(LatticeStructure<D> &r_, const unsigned io_, F &&f_)
{
  Coordinates<std::size_t, D + 1> Ld(r_.Ld), Lt(r_.Lt);
  const std::size_t begin = std::size_t(io_) * r_.Nd;
  const std::size_t end = begin + r_.Nd;
  for (std::size_t i = begin; i < end; ++i)
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
pw_reduce_slots_mpi(G &global_, const Eigen::Array<T, -1, 1> &local_)
{
#pragma omp barrier
#pragma omp master
  {
    global_.orb_sum.resize(local_.size());
    global_.orb_sum.setZero();
  }
#pragma omp barrier
#pragma omp critical
  {
    global_.orb_sum += local_;
  }
#pragma omp barrier
#pragma omp master
  {
    kmpi::sum_all(global_.orb_sum);
  }
#pragma omp barrier
  return global_.orb_sum;
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_pwave()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_pwave\n");
  const std::string base_grp = "/Calculation/p_wave/";
  std::string tmp = base_grp + "NumRandomsTransient";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    Global.calculate_p_wave = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, &file, tmp);
      Global.calculate_p_wave = true;
    } catch (H5::Exception &e) {
      debug_message("p_wave: no need to calculate it.\n");
    }
    file.close();
  }
#pragma omp barrier
  bool local_calculate_p_wave = false;
#pragma omp critical
  local_calculate_p_wave = Global.calculate_p_wave;
#pragma omp barrier
  if (!local_calculate_p_wave)
    return;
#pragma omp master
  std::cout << "Calculating spatially resolved PWave.\n";
#pragma omp barrier
  int randoms_transient, n_transient, n_stages, stage_itr, depth;
  value_type mixing;
  int prv_itr = 0;
  const unsigned max_bonds = h.pr.bonds.max_bonds;
  const std::size_t n_slots = std::size_t(max_bonds) * r.Sizet;
  bool has_map = false;
  Eigen::Array<value_type, -1, 1> init_map(n_slots);
  init_map.setZero();
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
  
  const value_type seed_imag = h.pr.Delta0.imag().abs().maxCoeff();
  if (seed_imag > value_type(1e-12)) {
#pragma omp critical
    {
      std::cerr
	<< "p_wave: the real-space Anderson solver currently supports only "
	"real nearest-neighbour  Delta. Use a real p-wave seed.\n";
      exit(1);
    }
  }
  
  Eigen::Array<T, -1, -1> seed_delta = h.pr.Delta0;
  h.pr.symmetrize_bonds(seed_delta);
  h.pr.allocate(h.bdg.nn_delta, true);
  h.pr.broadcast(seed_delta, h.bdg.nn_delta);

  if (has_map)
    {
      pw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
      {
	const unsigned io = static_cast<unsigned>(locl_ / r.Nd);
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
	  const std::size_t k = b + std::size_t(max_bonds) * glob_;
	  h.bdg.nn_delta(b, locl_) = T(init_map(k) / h.bdg.energy_scale);
	}
      });
    }

  h.bdg.nn_delta = h.bdg.nn_delta.real().template cast<T>().eval();
  p_wave(randoms_transient, n_transient, n_stages, stage_itr, depth, mixing, prv_itr);
}

template <typename T, unsigned D>
void Simulation<T, D>::p_wave(
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
  debug_message("Entered spatially resolved PWave\n");
  if constexpr (D != 2 || !pairing::is_p_wave || pairing::is_s_wave) {
    (void)randoms_transient_;
    (void)n_transient_;
    (void)n_stages_;
    (void)stage_itr_;
    (void)depth_;
    (void)mixing_;
    (void)prv_itr_;
#pragma omp master
    {
      if constexpr (D != 2)
	std::cerr << "p_wave: only implemented for D = 2.\n";
      else if constexpr (!pairing::is_p_wave)
	std::cerr << "p_wave: nearest-neighbour pairing is not enabled.\n";
      else
	std::cerr << "p_wave: use PAIRING = 2. "
	  "PAIRING = 3 requires the coupled SPWaveClean solver.\n";
    }
#pragma omp barrier
    exit(1);
  }
  else {
    using RealVec = Eigen::Matrix<value_type, -1, 1>;
    using RealMat = Eigen::Matrix<value_type, -1, -1>;
    
    const value_type energy_scale = h.bdg.energy_scale;
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type n_ranks = static_cast<value_type>(kmpi::size());
    
    const unsigned max_bonds = h.pr.bonds.max_bonds;
    const std::size_t n_slots = std::size_t(max_bonds) * r.Sizet;
    const unsigned n_family_slots = max_bonds * r.Orb;

    const int num_itr = n_transient_ + n_stages_ * stage_itr_;
    const unsigned m = depth_ > 0
      ? std::min<unsigned>(static_cast<unsigned>(depth_), static_cast<unsigned>(n_slots)) : 0;
    
    const Eigen::Array<value_type, -1, 1> coefs =
      Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);
    
    Eigen::Array<T, -1, -1> local_delta(max_bonds, r.Sized);
    Eigen::Array<T, -1, -1> sym_delta(max_bonds, r.Sized);
    Eigen::Array<T, -1, 1> ket(2 * r.Sized);
    Eigen::Array<T, -1, 1> ket_ref(2 * r.Sized);

    // orbital average
    Eigen::Array<T, -1, 1> local_family_sum(n_family_slots);
    Eigen::Array<T, -1, 1> local_family_count(n_family_slots);
    Eigen::Array<T, -1, 1> family_count(n_family_slots);
    Eigen::Array<T, -1, 1> family_sum(n_family_slots);

    struct PWaveState {
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

    PWaveState *state = nullptr;
#pragma omp single copyprivate(state)
    {
      state = new PWaveState;

      state->mean = Eigen::Array<value_type, -1, 1>::Zero(n_slots);
      state->map = Eigen::Array<value_type, -1, 1>::Zero(n_slots);
      state->residual = Eigen::Array<value_type, -1, 1>::Zero(num_itr);
      state->hist =
        Eigen::Array<value_type, -1, -1>::Zero(num_itr + 1, n_family_slots);
      
      state->dX = RealMat::Zero(n_slots, m);
      state->dF = RealMat::Zero(n_slots, m);
      state->x_prev = RealVec::Zero(n_slots);
      state->f_prev = RealVec::Zero(n_slots);
    }

    constexpr value_type max_log_step = value_type(0.405);
    constexpr value_type amp_min = value_type(1e-12);

    h.generate_disorder();
    KPM_Vector<T, D> phi(2, *this);
    LocalUnitary<T, D> U(r, phi);
    
    Eigen::Array<T, -1, 1> alive_t = Eigen::Array<T, -1, 1>::Zero(r.Sized);
    pw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t)
    {
      alive_t(locl_) = T(1);
    });
    for (const auto &pos : h.hV.position)
      for (const auto idx : pos)
	alive_t(idx) = T(0);
    phi.v.setZero();
    phi.set_index(0);
    phi.v.col(0).head(r.Sized) = alive_t.matrix();
    phi.Exchange_Boundaries();
    alive_t = phi.v.col(0).head(r.Sized).array();
    const Eigen::Array<value_type, -1, 1> alive = alive_t.real();

    auto exchange_bond_field = [&](Eigen::Array<T, -1, -1> &field_)
    {
      for (unsigned b = 0; b < max_bonds; ++b) {
        phi.v.setZero();
        phi.set_index(0);
        phi.v.col(0).head(r.Sized) = field_.row(b).transpose().matrix();
        phi.Exchange_Boundaries();
        field_.row(b) = phi.v.col(0).head(r.Sized).array().transpose();
      }
    };
    auto valid_target = [&](const std::size_t locl_, const unsigned io_,
			    const unsigned b_, std::size_t &target_) -> bool
    {
      const std::ptrdiff_t shift = h.pr.bonds.dist_tile(b_, io_);
      const std::ptrdiff_t jt = static_cast<std::ptrdiff_t>(locl_) + shift;
      if (jt < 0 || jt >= static_cast<std::ptrdiff_t>(r.Sized))
        return false;
      target_ = static_cast<std::size_t>(jt);
      return alive(locl_) > value_type(0.5) && alive(target_) > value_type(0.5);
    };
    auto erase_invalid_bonds = [&](Eigen::Array<T, -1, -1> &field_)
    {
      pw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t)
      {
        const unsigned io = static_cast<unsigned>(locl_ / r.Nd);
        for (unsigned b = 0; b < max_bonds; ++b)
	  {
	    if (b >= h.pr.bonds.NBonds(io))
	      {
		field_(b, locl_) = T(0);
		continue;
	      }
	    std::size_t target = 0;
	    if (!valid_target(locl_, io, b, target))
	      field_(b, locl_) = T(0);
	  }
      });
    };
    auto exchange_nn_delta = [&]()
    {
      erase_invalid_bonds(h.bdg.nn_delta);
      exchange_bond_field(h.bdg.nn_delta);
      h.bdg.nn_delta = h.bdg.nn_delta.real().template cast<T>().eval();
      erase_invalid_bonds(h.bdg.nn_delta);
    };

#pragma omp barrier
    exchange_nn_delta();
#pragma omp barrier
    sym_delta.setZero();
    h.pr.symmetrize(h.bdg.nn_delta, sym_delta);
    h.bdg.nn_delta = sym_delta;
#pragma omp barrier
    exchange_nn_delta();
#pragma omp barrier

#pragma omp master
    state->mean.setZero();
#pragma omp barrier
#pragma omp critical
    {
      pw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
      {
        const unsigned io = static_cast<unsigned>(locl_ / r.Nd);
        for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
          const std::size_t k = b + std::size_t(max_bonds) * glob_;
          state->mean(k) = std::real(h.bdg.nn_delta(b, locl_));
        }
      });
    }
#pragma omp barrier
#pragma omp master
    {
      kmpi::sum_all(state->mean);
      state->mean /= n_ranks;
    }
#pragma omp barrier
    h.bdg.nn_delta.setZero();
    pw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
    {
      const unsigned io = static_cast<unsigned>(locl_ / r.Nd);
      for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
        const std::size_t k = b + std::size_t(max_bonds) * glob_;
        h.bdg.nn_delta(b, locl_) = T(state->mean(k));
      }
    });

#pragma omp barrier
    exchange_nn_delta();
#pragma omp barrier

    local_family_count.setZero();
    for (unsigned io = 0; io < r.Orb; ++io)
      pw_for_each_interior_orbital(r, io, [&](const std::size_t locl_, const std::size_t)
      {
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
	  std::size_t target = 0;
	  if (valid_target(locl_, io, b, target))
	    local_family_count(b + io * max_bonds) += T(1);
	}
      });
    family_count = pw_reduce_slots_mpi(Global, local_family_count);
    local_family_sum.setZero();
    for (unsigned io = 0; io < r.Orb; ++io)
      pw_for_each_interior_orbital(r, io, [&](const std::size_t locl_, const std::size_t)
      {
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
	  std::size_t target = 0;
	  if (valid_target(locl_, io, b, target))
	    local_family_sum(b + io * max_bonds) += h.bdg.nn_delta(b, locl_);
	}
      });
    family_sum = pw_reduce_slots_mpi(Global, local_family_sum);
#pragma omp master
    {
      for (unsigned k = 0; k < n_family_slots; ++k) {
        const value_type n = std::real(family_count(k));
        if (n > value_type(0))
          state->hist(0, k) = std::real(family_sum(k)) / n * energy_scale;
      }
    }
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
        phi.chebyshev_sum(coefs, ket_ref.data());
        for (unsigned io = 0; io < r.Orb; ++io)
          for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
            ket = ket_ref;
            U.template pairs<1>(h.pr.bonds, 1, io, b, -1, ket);
            const Eigen::Array<value_type, -1, 1> pp =
              ket.abs2().head(r.Sized);
            const Eigen::Array<value_type, -1, 1> hh =
              ket.abs2().tail(r.Sized);
	    
            pw_for_each_interior_orbital(r, io,
	    [&](const std::size_t locl_, const std::size_t)
	    {
	      std::size_t target = 0;
	      if (!valid_target(locl_, io, b, target)) {
		local_delta(b, locl_) = T(0);
		return;
	      }
	      const value_type f_b =
		value_type(0.5) * h.pr.V * (pp(locl_) - hh(target));
	      local_delta(b, locl_) +=
		(T(f_b) - local_delta(b, locl_)) / value_type(vec + 1);
	    });
          }
      }
      
#pragma omp barrier
      exchange_bond_field(local_delta);
#pragma omp barrier
      sym_delta.setZero();
      h.pr.symmetrize(local_delta, sym_delta);
      sym_delta = sym_delta.real().template cast<T>().eval();
      erase_invalid_bonds(sym_delta);
      
#pragma omp master
      state->map.setZero();
#pragma omp barrier
#pragma omp critical
      {
        pw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
	{
          const unsigned io = static_cast<unsigned>(locl_ / r.Nd);
          for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
            const std::size_t k = b + std::size_t(max_bonds) * glob_;
            state->map(k) = std::real(sym_delta(b, locl_));
          }
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
        const Eigen::Array<value_type, -1, 1> amp =
          state->mean.abs().max(amp_min);
        const Eigen::Array<value_type, -1, 1> amp_map =
          state->map.abs().max(amp_min);
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
        if (step.dot(f) <= value_type(0))
          {
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
      h.bdg.nn_delta.setZero();
      pw_for_each_interior(r, [&](const std::size_t locl_, const std::size_t glob_)
      {
        const unsigned io = static_cast<unsigned>(locl_ / r.Nd);
        for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
          const std::size_t k = b + std::size_t(max_bonds) * glob_;
          h.bdg.nn_delta(b, locl_) = T(state->mean(k));
        }
      }); 
#pragma omp barrier
      exchange_nn_delta();
#pragma omp barrier
      local_family_sum.setZero();
      for (unsigned io = 0; io < r.Orb; ++io)
        pw_for_each_interior_orbital(r, io,
	[&](const std::size_t locl_, const std::size_t)
	{
	  for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
	    std::size_t target = 0;
	    if (valid_target(locl_, io, b, target))
	      local_family_sum(b + io * max_bonds) += h.bdg.nn_delta(b, locl_);
	  }
	});
      family_sum = pw_reduce_slots_mpi(Global, local_family_sum);
#pragma omp master
      {
        for (unsigned k = 0; k < n_family_slots; ++k) {
          const value_type n = std::real(family_count(k));
          if (n > value_type(0))
            state->hist(itr + 1, k) =
              std::real(family_sum(k)) / n * energy_scale;
        }
      }
#pragma omp barrier
    }
    
    store_p_wave(
      static_cast<unsigned>(prv_itr_ + num_itr),
      state->hist,
      state->residual,
      state->mean
    );
    
#pragma omp barrier
#pragma omp master
    delete state;
#pragma omp barrier
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_p_wave(
  const unsigned total_steps_,
  const Eigen::Array<value_type, -1, -1> &hist_,
  const Eigen::Array<value_type, -1, 1> &residual_,
  const Eigen::Array<value_type, -1, 1> &map_
)
  requires Complex<T>
{
  debug_message("Entered store_p_wave\n");
  
  const value_type energy_scale = h.bdg.energy_scale;
  const unsigned max_bonds = h.pr.bonds.max_bonds;
  const std::size_t n_slots = std::size_t(max_bonds) * r.Sizet;
  
#pragma omp master
  if (kmpi::is_root()) {
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/p_wave/";
    
    std::string ng = base_grp + "Hist";
    write_hdf5(hist_, &file, ng);
    
    const Eigen::Array<value_type, -1, -1> res = residual_;
    ng = base_grp + "Residual";
    write_hdf5(res, &file, ng);
    
    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);

    Eigen::Array<value_type, -1, -1> delta_map(n_slots, 1);
    delta_map.col(0) = map_ * energy_scale;
    ng = base_grp + "Map";
    write_hdf5(delta_map, &file, ng);
    
    Eigen::Array<value_type, -1, -1> ranks(1, 1);
    ranks(0, 0) = static_cast<value_type>(kmpi::size());
    ng = base_grp + "NumRanks";
    write_hdf5(ranks, &file, ng);
    
    Eigen::Array<value_type, -1, -1> max_bonds_out(1, 1);
    max_bonds_out(0, 0) = static_cast<value_type>(max_bonds);
    ng = base_grp + "MaxBonds";
    write_hdf5(max_bonds_out, &file, ng);

    file.close();
  }
  
#pragma omp master
  kmpi::barrier();
#pragma omp barrier
  
  debug_message("Left store_p_wave\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
