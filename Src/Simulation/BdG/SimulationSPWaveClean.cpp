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
#include "LocalUnitary.hpp"
#include "Loop.hpp"
#include "Coefficients.hpp"
#include "mpi_utils.hpp"

template <typename T, unsigned D>
void Simulation<T, D>::calc_spwave_clean()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_spwave_clean\n");
  
  const std::string base_grp = "/Calculation/sp_wave_c/";
  std::string tmp = base_grp + "NumRandomsTransient";
  
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
    std::cout << "Calculating SPWaveClean.\n";
#pragma omp barrier
    int randoms_transient;
    int n_transient;
    int n_stages;
    int stage_itr;
    int depth;
    value_type mixing;
    int avg_bonds = 1;
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
      path = base_grp + "AverageBonds";
      get_hdf5<int>(&avg_bonds, &file, path);
      file.close();
    }
#pragma omp master
    std::cout << "sp_wave_clean: bond average " << (avg_bonds ? "on" : "off") << "\n";
#pragma omp barrier
    sp_wave_clean(
      randoms_transient,
      n_transient,
      n_stages,
      stage_itr,
      depth,
      mixing,
      avg_bonds != 0
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::sp_wave_clean(
  const int randoms_transient_,
  const int n_transient_,
  const int n_stages_,
  const int stage_itr_,
  const int depth_,
  const value_type mixing_,
  const bool avg_bonds_
)
  requires Complex<T>
{
  debug_message("Entered SPWaveClean\n");
  if constexpr (
    D != 2 || !(pairing::is_s_wave && pairing::is_p_wave)
  ) {
    (void)randoms_transient_;
    (void)n_transient_;
    (void)n_stages_;
    (void)stage_itr_;
    (void)depth_;
    (void)mixing_;
    (void)avg_bonds_;
#pragma omp master
    std::cerr
      << (D != 2
	  ? "sp_wave_clean: only implemented for D = 2.\n"
	  : "sp_wave_clean: this solver updates the on-site and bond "
	  "channels together. Build with PAIRING=3.\n");
    
#pragma omp barrier
    exit(1);
  } else {
    using Vec = Eigen::Matrix<T, -1, 1>;
    using Mat = Eigen::Matrix<T, -1, -1>;
    
    const value_type energy_scale = h.bdg.energy_scale;
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type n_cells = r.Nt;
    const value_type n_ranks = kmpi::size();
    const int num_itr = n_transient_ + n_stages_ * stage_itr_;
    
    const unsigned max_bonds = h.pr.bonds.max_bonds;
    const unsigned n_slots = max_bonds * r.Orb;
    unsigned n_active = 0;
    for (unsigned io = 0; io < r.Orb; ++io)
      n_active += h.pr.bonds.NBonds(io);
    
    const unsigned n_var = r.Orb + n_active;
    const unsigned m = std::min<unsigned>(depth_, n_var);
    const value_type max_step = value_type(0.5) / energy_scale; // Anderson
    
    const Eigen::Array<value_type, -1, 1> coefs =
      Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);
    
    Eigen::Array<T, -1, 1> mean_s(r.Orb);
    Eigen::Array<T, -1, 1> map_s(r.Orb);
    Eigen::Array<T, -1, 1> local_s(r.Orb);
    Eigen::Array<T, -1, 1> per_orb_s(r.Orb);

    Eigen::Array<T, -1, -1> mean_b(max_bonds, r.Orb);
    Eigen::Array<T, -1, -1> map_b(max_bonds, r.Orb);
    Eigen::Array<T, -1, -1> local_b(max_bonds, r.Orb);
    Eigen::Array<T, -1, 1> per_orb_p(r.Orb);
    Eigen::Array<T, -1, 1> per_orb_h(r.Orb);

    Eigen::Array<value_type, -1, 1> residual(num_itr);
    Eigen::Array<value_type, -1, 1> residual_s(num_itr);
    Eigen::Array<value_type, -1, 1> residual_b(num_itr);
    
    Eigen::Array<value_type, -1, -1> map_hist_s =
      Eigen::Array<value_type, -1, -1>::Zero(num_itr, r.Orb);
    Eigen::Array<value_type, -1, -1> map_hist_b =
      Eigen::Array<value_type, -1, -1>::Zero(num_itr, n_slots);
    Eigen::Array<value_type, -1, -1> hist_s =
      Eigen::Array<value_type, -1, -1>::Zero(num_itr + 1, r.Orb);
    Eigen::Array<value_type, -1, -1> hist_b =
      Eigen::Array<value_type, -1, -1>::Zero(num_itr + 1, n_slots);
    
    Mat dX = Mat::Zero(n_var, m);
    Mat dF = Mat::Zero(n_var, m);
    Vec x_prev(n_var);
    Vec f_prev(n_var);
    unsigned n_hist = 0;
    value_type res_prev = 0;
    
    const auto pack = [&](const Eigen::Array<T, -1, 1> &s_delta, const Eigen::Array<T, -1, -1> &bond_delta)
    {
      Vec v(n_var);
      for (unsigned io = 0; io < r.Orb; ++io)
	v(io) = s_delta(io);
      unsigned k = r.Orb;
      for (unsigned io = 0; io < r.Orb; ++io)
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
	  v(k++) = bond_delta(b, io);
      return v;
    };
    
    const auto unpack = [&](const Vec &v, Eigen::Array<T, -1, 1> &s_delta, Eigen::Array<T, -1, -1> &bond_delta)
    {
      for (unsigned io = 0; io < r.Orb; ++io)
	s_delta(io) = v(io);
      bond_delta.setZero();
      unsigned k = r.Orb;
      for (unsigned io = 0; io < r.Orb; ++io)
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
	  bond_delta(b, io) = v(k++);
    };
    
    mean_s = h.pr.SDelta0;
    mean_b = h.pr.Delta0;
    
    h.pr.symmetrize_bonds(mean_b);
    if (avg_bonds_)
      h.pr.average_bonds(mean_b);
    
    h.pr.broadcast_s(mean_s, h.bdg.s_delta);
    h.pr.broadcast(mean_b, h.bdg.nn_delta);
    
    hist_s.row(0) = mean_s.real().transpose() * energy_scale;
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
        hist_b(0, b + io * max_bonds) = mean_b(b, io).real() * energy_scale;
    
#pragma omp master
    {
      Global.orb_sum.resize(r.Orb);
      Global.nn_sum.resize(n_slots);
    }
#pragma omp barrier
    h.generate_disorder();
    KPM_Vector<T, D> phi(2, *this);
    LocalUnitary<T, D> U(r, phi);
    Eigen::Array<T, -1, 1> ket(2 * r.Sized);
    Eigen::Array<T, -1, 1> ket_ref(2 * r.Sized);
    for (int itr = 0; itr < num_itr; ++itr) {
      rnd.init_random(seed_v);
      h.rnd.init_random(seed_h);
      
      const int stage = itr < n_transient_ ? 0 : 1 + (itr - n_transient_) / stage_itr_;
      const int n_vec = randoms_transient_ << stage;
      const bool fresh = itr == 0 || (itr >= n_transient_ && (itr - n_transient_) % stage_itr_ == 0);
      
      local_s.setZero();
      local_b.setZero();
      
      for (int vec = 0; vec < n_vec; ++vec) {
        const value_type weight = value_type(1) / value_type(vec + 1);
        h.generate_twists();
        phi.initiate_phases();
        phi.set_index(0);
        phi.initiate_vector();
        phi.v.col(0) *= std::sqrt(size);
        phi.Exchange_Boundaries();
        phi.chebyshev_sum(coefs, ket_ref.data());

        // On-site
        ket = ket_ref;
        U.template onsite<1>(T(1), ket);
        const Eigen::Array<value_type, -1, 1> upsilon =
	  ket.abs2().head(r.Sized) - ket.abs2().tail(r.Sized);
        const Eigen::Array<T, -1, 1> f_s = 0.5 * h.pr.U * upsilon;
        h.pr.orbital_sum(f_s, per_orb_s);
        local_s += weight * (per_orb_s - local_s);
	
        // Bond
        for (unsigned io = 0; io < r.Orb; ++io)
          for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
            ket = ket_ref;
            U.template pairs<1>(h.pr.bonds, 1, io, b, -1, ket);
            const Eigen::Array<T, -1, 1> pp =
	      ket.abs2().head(r.Sized).template cast<T>();
            const Eigen::Array<T, -1, 1> hh =
	      ket.abs2().tail(r.Sized).template cast<T>();
            h.pr.orbital_sum(pp, per_orb_p);
            h.pr.orbital_sum(hh, per_orb_h);
            const unsigned jo = h.pr.bonds.target_orb(b, io);
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
          for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
            Global.nn_sum(b + io * max_bonds) += local_b(b, io);
        }
      }
#pragma omp barrier
#pragma omp master
      {
        kmpi::sum_all(Global.orb_sum);
        kmpi::sum_all(Global.nn_sum);
      }
#pragma omp barrier
      
      for (unsigned io = 0; io < r.Orb; ++io)
        map_s(io) = Global.orb_sum(io) / (n_cells * n_ranks);
      
      map_b.setZero();
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
          map_b(b, io) = Global.nn_sum(b + io * max_bonds) / (n_cells * n_ranks);
      h.pr.symmetrize_bonds(map_b);
      if (avg_bonds_)
        h.pr.average_bonds(map_b);
      
      const Vec current = pack(mean_s, mean_b);
      const Vec mapped = pack(map_s, map_b);
      
      const Vec f_lin = mapped - current;
      residual_s(itr) = (map_s - mean_s).matrix().cwiseAbs().maxCoeff() * energy_scale;
      residual_b(itr) = f_lin.tail(n_active).cwiseAbs().maxCoeff() * energy_scale;
      residual(itr) = std::max(residual_s(itr), residual_b(itr));
      
      map_hist_s.row(itr) = map_s.real().transpose() * energy_scale;
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
          map_hist_b(itr, b + io * max_bonds) = map_b(b, io).real() * energy_scale;
      
      const Vec x = current;
      const Vec f = f_lin;
      const value_type res = f.cwiseAbs().maxCoeff();
      if (fresh || res > 1.5 * res_prev)
        n_hist = 0;
      else if (m > 0)
	{
	  if (n_hist == m)
	    {
	      dX.leftCols(m - 1) = dX.rightCols(m - 1).eval();
	      dF.leftCols(m - 1) = dF.rightCols(m - 1).eval();
	      --n_hist;
	    }
	  dX.col(n_hist) = x - x_prev;
	  dF.col(n_hist) = f - f_prev;
	  ++n_hist;
	}
      x_prev = x;
      f_prev = f;
      res_prev = res;
      Vec step = mixing_ * f;
      if (n_hist > 0)
	{
	  const Mat dFh = dF.leftCols(n_hist);
	  const Mat dXh = dX.leftCols(n_hist);
	  Eigen::CompleteOrthogonalDecomposition<Mat> cod(dFh);
	  cod.setThreshold(1e-8);
	  const Vec gamma = cod.solve(f);
	  const Mat G = dXh + mixing_ * dFh;
	  step -= G * gamma;
	}
      if (!step.allFinite() ||
	  std::real(step.dot(f)) <= 0 ||
	  step.cwiseAbs().maxCoeff() > max_step) {
	step = mixing_ * f;
	n_hist = 0;
      }
      
      const value_type step_size = step.cwiseAbs().maxCoeff();
      if (step_size > max_step)
	step *= max_step / step_size;
      
      const Vec next = x + step;
      unpack(next, mean_s, mean_b);
      
      h.pr.symmetrize_bonds(mean_b);
      if (avg_bonds_)
        h.pr.average_bonds(mean_b);
      
      h.pr.broadcast_s(mean_s, h.bdg.s_delta);
      h.pr.broadcast(mean_b, h.bdg.nn_delta);
      
      hist_s.row(itr + 1) = mean_s.real().transpose() * energy_scale;
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0,  B = h.pr.bonds.NBonds(io);  b < B;  ++b)
          hist_b(itr + 1, b + io * max_bonds) = mean_b(b, io).real() * energy_scale;
      
#pragma omp barrier
    }
    
    store_sp_wave_clean(
      num_itr,
      mean_s,
      mean_b,
      residual,
      residual_s,
      residual_b,
      map_hist_s,
      map_hist_b,
      hist_s,
      hist_b
    );
  }
}


template <typename T, unsigned D>
void Simulation<T, D>::store_sp_wave_clean(
  const unsigned total_steps_,
  const Eigen::Array<T, -1, 1> &mean_s_,
  const Eigen::Array<T, -1, -1> &mean_b_,
  const Eigen::Array<value_type, -1, 1> &residual_,
  const Eigen::Array<value_type, -1, 1> &residual_s_,
  const Eigen::Array<value_type, -1, 1> &residual_b_,
  const Eigen::Array<value_type, -1, -1> &map_hist_s_,
  const Eigen::Array<value_type, -1, -1> &map_hist_b_,
  const Eigen::Array<value_type, -1, -1> &hist_s_,
  const Eigen::Array<value_type, -1, -1> &hist_b_
)
  requires Complex<T>
{
  debug_message("Entered store_sp_wave_clean\n");
  
#pragma omp master
  if (kmpi::is_root()) {
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/sp_wave_c/";
    std::string ng;
    
    ng = base_grp + "HistS";
    write_hdf5(hist_s_, &file, ng);
    ng = base_grp + "HistB";
    write_hdf5(hist_b_, &file, ng);
    ng = base_grp + "MapHistS";
    write_hdf5(map_hist_s_, &file, ng);
    ng = base_grp + "MapHistB";
    write_hdf5(map_hist_b_, &file, ng);
    
    const Eigen::Array<value_type, -1, -1> res = residual_;
    ng = base_grp + "Residual";
    write_hdf5(res, &file, ng);
    const Eigen::Array<value_type, -1, -1> res_s = residual_s_;
    ng = base_grp + "ResidualS";
    write_hdf5(res_s, &file, ng);
    const Eigen::Array<value_type, -1, -1> res_b = residual_b_;
    ng = base_grp + "ResidualB";
    write_hdf5(res_b, &file, ng);

    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);

    Eigen::Array<value_type, -1, -1> s_ri(2, r.Orb);
    for (unsigned io = 0; io < r.Orb; ++io)
      {
	s_ri(0, io) = mean_s_(io).real() * h.bdg.energy_scale;
	s_ri(1, io) = mean_s_(io).imag() * h.bdg.energy_scale;
      }
    ng = base_grp + "DeltaS";
    write_hdf5(s_ri, &file, ng);
    
    const unsigned max_bonds = h.pr.bonds.max_bonds;
    const unsigned n_slots = max_bonds * r.Orb;
    Eigen::Array<value_type, -1, -1> b_ri = Eigen::Array<value_type, -1, -1>::Zero(2, n_slots);
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
        const std::size_t k = b + io * max_bonds;
        b_ri(0, k) = mean_b_(b, io).real() * h.bdg.energy_scale;
        b_ri(1, k) = mean_b_(b, io).imag() * h.bdg.energy_scale;
      }
    ng = base_grp + "DeltaB";
    write_hdf5(b_ri, &file, ng);
    
    Eigen::Array<value_type, -1, -1> ranks(1, 1);
    ranks(0, 0) = static_cast<value_type>(kmpi::size());
    ng = base_grp + "NumRanks";
    write_hdf5(ranks, &file, ng);
  }
  
#pragma omp master
  kmpi::barrier();
  
#pragma omp barrier
  
  debug_message("Left store_sp_wave_clean\n");
}


#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
