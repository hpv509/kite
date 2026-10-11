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
void Simulation<T, D>::calc_pwave_clean()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_pwave_clean\n");
  std::string base_grp = "/Calculation/p_wave_c/";
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
  if (local_calculate_p_wave) {
#pragma omp master
    std::cout << "Calculating PWave.\n";
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
    std::cout << "p_wave_clean: bond average " << (avg_bonds ? "on" : "off")
              << "\n";
#pragma omp barrier
    p_wave_clean(
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
void Simulation<T, D>::p_wave_clean(
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
  debug_message("Entered PWaveClean\n");
  if constexpr (D != 2 || !pairing::is_p_wave || pairing::is_s_wave) {
    (void)randoms_transient_;
    (void)n_transient_;
    (void)n_stages_;
    (void)stage_itr_;
    (void)depth_;
    (void)mixing_;
    (void)avg_bonds_;
#pragma omp master
    {
      if constexpr (D != 2)
	std::cerr << "p_wave_clean: only implemented for D = 2.\n";
      else if constexpr (!pairing::is_p_wave)
	std::cerr << "p_wave_clean: nearest-neighbour pairing is not enabled.\n";
      else
	std::cerr << "p_wave_clean: use PAIRING = 2. "
	  "PAIRING = 3 requires the coupled SPWaveClean solver.\n";
    }
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
    
    const unsigned m = std::min<unsigned>(depth_, n_active);
    
    const Eigen::Array<value_type, -1, 1> coefs =
      Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);
    
    Eigen::Array<T, -1, -1> mean_delta(max_bonds, r.Orb);
    Eigen::Array<T, -1, -1> map_delta(max_bonds, r.Orb);
    Eigen::Array<T, -1, -1> local_delta(max_bonds, r.Orb);
    Eigen::Array<T, -1, 1> per_orb_p(r.Orb);
    Eigen::Array<T, -1, 1> per_orb_h(r.Orb);
    Eigen::Array<value_type, -1, 1> residual(num_itr);
    Eigen::Array<value_type, -1, -1> map_hist =
      Eigen::Array<value_type, -1, -1>::Zero(num_itr, n_slots);
    
    Mat dX = Mat::Zero(n_active, m);
    Mat dF = Mat::Zero(n_active, m);
    
    Vec x_prev(n_active);
    Vec f_prev(n_active);
    
    unsigned n_hist = 0;
    value_type res_prev = 0;
    
    constexpr value_type max_log_step = 0.405;
    constexpr value_type amp_min = 1e-12;

    const auto pack = [&](const Eigen::Array<T, -1, -1> &delta)
    {
      Vec v(n_active);
      unsigned k = 0;
      for (unsigned io = 0; io < r.Orb; ++io)
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
	  v(k++) = delta(b, io);
      return v;
    };
    
    const auto unpack = [&](const Vec &v, Eigen::Array<T, -1, -1> &delta)
    {
      delta.setZero();
      unsigned k = 0;
      for (unsigned io = 0; io < r.Orb; ++io)
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
	  delta(b, io) = v(k++);
    };
    
    mean_delta = h.pr.Delta0;
    h.pr.symmetrize_bonds(mean_delta);
    if (avg_bonds_)
      h.pr.average_bonds(mean_delta);
    h.pr.broadcast(mean_delta, h.bdg.nn_delta);

    const auto record = [&](const unsigned row) {
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
          Global.nn_delta_hist(row, b + io * max_bonds) =
            mean_delta(b, io) * energy_scale;
    };
#pragma omp master
    {
      Global.nn_sum.resize(n_slots);
      Global.nn_delta_hist.resize(num_itr + 1, n_slots);
      Global.nn_delta_hist.setZero();
      record(0);
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
      const bool fresh = itr == 0 ||
	(itr >= n_transient_ && (itr - n_transient_) % stage_itr_ == 0);
      
      local_delta.setZero();
      
      for (int vec = 0; vec < n_vec; ++vec) {
	const value_type weight = 1.0 / (vec + 1);
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
            const Eigen::Array<T, -1, 1> pp =
              ket.abs2().head(r.Sized).template cast<T>();
            const Eigen::Array<T, -1, 1> hh =
              ket.abs2().tail(r.Sized).template cast<T>();
            h.pr.orbital_sum(pp, per_orb_p);
            h.pr.orbital_sum(hh, per_orb_h);

            const unsigned jo = h.pr.bonds.target_orb(b, io);
            const value_type tmp_b = 0.5 * h.pr.V;
            const T f_b = tmp_b * (per_orb_p(io) - per_orb_h(jo));
            local_delta(b, io) += weight * (f_b - local_delta(b, io));
          }
      }
      
#pragma omp barrier
#pragma omp master
      Global.nn_sum.setZero();
#pragma omp barrier
      
#pragma omp critical
      {
        for (unsigned io = 0; io < r.Orb; ++io)
          for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
            Global.nn_sum(b + io * max_bonds) += local_delta(b, io);
      }
#pragma omp barrier

#pragma omp master
      kmpi::sum_all(Global.nn_sum);
#pragma omp barrier
      
      map_delta.setZero();
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
          map_delta(b, io) = Global.nn_sum(b + io * max_bonds) / (n_cells * n_ranks);
      h.pr.symmetrize_bonds(map_delta);
      if (avg_bonds_)
        h.pr.average_bonds(map_delta);

      const Vec current = pack(mean_delta);
      const Vec mapped = pack(map_delta);
      
      const Vec f_lin = mapped - current;
      residual(itr) = f_lin.cwiseAbs().maxCoeff() * energy_scale;
      
      for (unsigned io = 0; io < r.Orb; ++io)
	for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b)
	  map_hist(itr, b + io * max_bonds) = map_delta(b, io).real() * energy_scale;
      
      const Eigen::Array<value_type, -1, 1> amp = current.array().abs().max(amp_min);
      const Eigen::Array<value_type, -1, 1> amp_map = mapped.array().abs().max(amp_min);
      const Vec x = amp.log().template cast<T>().matrix();
      const Vec f = (amp_map.log() - amp.log()).template cast<T>().matrix();
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
      if (std::real(step.dot(f)) <= 0)
	{
	  step = mixing_ * f;
	  n_hist = 0;
	}
      const value_type s_max = step.cwiseAbs().maxCoeff();
      if (s_max > max_log_step)
	step *= max_log_step / s_max;
      const Eigen::Array<T, -1, 1> phase = mapped.array() / amp_map;
      
      const Vec next =
	((x + step).real().array().exp().template cast<T>() * phase).matrix();
      
      unpack(next, mean_delta);
      
      h.pr.symmetrize_bonds(mean_delta);
      if (avg_bonds_)
	h.pr.average_bonds(mean_delta);
      h.pr.broadcast(mean_delta, h.bdg.nn_delta);
      
#pragma omp barrier
#pragma omp master
      record(itr + 1);
#pragma omp barrier
    }
    store_p_wave_clean(num_itr, mean_delta, map_hist, residual);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_p_wave_clean(
  const unsigned total_steps_,
  const Eigen::Array<T, -1, -1> &mean_delta_,
  const Eigen::Array<value_type, -1, -1> &map_hist_,
  const Eigen::Array<value_type, -1, 1> &residual_
)
  requires Complex<T>
{
  debug_message("Entered store_p_wave_clean\n");
#pragma omp master
  if (kmpi::is_root()) {
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/p_wave_c/";
    std::string ng;
    
    Eigen::Array<value_type, -1, -1> hist_r = Global.nn_delta_hist.real();
    ng = base_grp + "Hist";
    write_hdf5(hist_r, &file, ng);
    
    ng = base_grp + "MapHist";
    write_hdf5(map_hist_, &file, ng);
    
    const Eigen::Array<value_type, -1, -1> res = residual_;
    ng = base_grp + "Residual";
    write_hdf5(res, &file, ng);
    
    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);

    const unsigned max_bonds = h.pr.bonds.max_bonds;
    const unsigned n_slots = max_bonds * r.Orb;
    Eigen::Array<value_type, -1, -1> delta_ri =
      Eigen::Array<value_type, -1, -1>::Zero(2, n_slots);
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0, B = h.pr.bonds.NBonds(io); b < B; ++b) {
        const std::size_t k = b + io * max_bonds;
        delta_ri(0, k) = mean_delta_(b, io).real() * h.bdg.energy_scale;
        delta_ri(1, k) = mean_delta_(b, io).imag() * h.bdg.energy_scale;
      }
    
    ng = base_grp + "Delta";
    write_hdf5(delta_ri, &file, ng);
    
    Eigen::Array<value_type, -1, -1> ranks(1, 1);
    ranks(0, 0) = static_cast<value_type>(kmpi::size());
    ng = base_grp + "NumRanks";
    write_hdf5(ranks, &file, ng);
  }

#pragma omp master
  kmpi::barrier();
#pragma omp barrier
  debug_message("Left store_p_wave_clean\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
