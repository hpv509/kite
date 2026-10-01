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
void Simulation<T, D>::calc_swave_anderson()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_swave_anderson\n");
  std::string base_grp = "/Calculation/s_wave_anderson/";
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
      debug_message("s_wave_anderson: no need to calculate it.\n");
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
    std::cout << "Calculating SWave (Anderson).\n";
#pragma omp barrier
    int randoms_transient, n_transient, n_stages, stage_itr, depth;
    value_type mixing;
#pragma omp critical
    {
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
      file.close();
    }
    s_wave_anderson(
      randoms_transient, n_transient, n_stages, stage_itr, depth, mixing
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::s_wave_anderson(
  const int randoms_transient_,
  const int n_transient_,
  const int n_stages_,
  const int stage_itr_,
  const int depth_,
  const value_type beta_
)
  requires Complex<T>
{
  debug_message("Entered SWaveAnderson\n");
  if constexpr (D != 2 || !pairing::is_s_wave) {
    (void)randoms_transient_;
    (void)n_transient_;
    (void)n_stages_;
    (void)stage_itr_;
    (void)depth_;
    (void)beta_;
#pragma omp master
    std::cerr
      << (D != 2 ? "s_wave_anderson: only implemented for D = 2.\n"
                 : "s_wave_anderson: this build has no on-site "
                   "channel. Build with PAIRING=1 or 3.\n");
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
    const unsigned m = std::min<unsigned>(depth_, r.Orb);

    const Eigen::Array<value_type, -1, 1> coefs =
      Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);

    Eigen::Array<T, -1, 1> mean_delta(r.Orb);
    Eigen::Array<T, -1, 1> map_delta_orb(r.Orb);
    Eigen::Array<T, -1, 1> local_delta(r.Orb);
    Eigen::Array<T, -1, 1> per_orb(r.Orb);
    Eigen::Array<value_type, -1, 1> residual(num_itr);
    Eigen::Array<value_type, -1, -1> map_hist(num_itr, r.Orb);
    Eigen::Array<value_type, -1, -1> hist(num_itr + 1, r.Orb);

    Mat dX = Mat::Zero(r.Orb, m), dF = Mat::Zero(r.Orb, m);
    Vec x_prev(r.Orb), f_prev(r.Orb);
    unsigned n_hist = 0;
    value_type res_prev = 0;
    constexpr value_type max_log_step = 0.405;
    constexpr value_type amp_min = 1e-12;

    mean_delta = h.pr.SDelta0;
    h.pr.broadcast_s(mean_delta, h.bdg.s_delta);
    hist.row(0) = mean_delta.real().transpose() * energy_scale;
#pragma omp master
    Global.orb_sum.resize(r.Orb);
#pragma omp barrier
    h.generate_disorder();
    KPM_Vector<T, D> phi(2, *this);
    LocalUnitary<T, D> U(r, phi);
    Eigen::Array<T, -1, 1> ket(2 * r.Sized);

    for (int itr = 0; itr < num_itr; ++itr) {
      rnd.init_random(seed_v);
      h.rnd.init_random(seed_h);
      const int stage =
        itr < n_transient_ ? 0 : 1 + (itr - n_transient_) / stage_itr_;
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
        const Eigen::Array<T, -1, 1> map_delta = 0.5 * h.pr.U * upsilon;

        h.pr.orbital_sum(map_delta, per_orb);
        local_delta += (per_orb - local_delta) / value_type(vec + 1);
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
#pragma omp master
      kmpi::sum_all(Global.orb_sum);
#pragma omp barrier
      for (unsigned io = 0; io < r.Orb; ++io)
        map_delta_orb(io) = Global.orb_sum(io) / (n_cells * n_ranks);

      const Vec f_lin = (map_delta_orb - mean_delta).matrix();
      residual(itr) = f_lin.cwiseAbs().maxCoeff() * energy_scale;
      map_hist.row(itr) = map_delta_orb.real().transpose() * energy_scale;

      const Eigen::Array<value_type, -1, 1> amp = mean_delta.abs().max(amp_min);
      const Eigen::Array<value_type, -1, 1> amp_map =
        map_delta_orb.abs().max(amp_min);
      const Vec x = amp.log().template cast<T>().matrix();
      const Vec f = (amp_map.log() - amp.log()).template cast<T>().matrix();
      const value_type res = f.cwiseAbs().maxCoeff();
      if (fresh || res > 1.5 * res_prev)
        n_hist = 0;
      else if (m > 0) {
        if (n_hist == m) {
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
      Vec step = beta_ * f;
      if (n_hist > 0) {
        const Mat dFh = dF.leftCols(n_hist);
        const Mat dXh = dX.leftCols(n_hist);
        Eigen::CompleteOrthogonalDecomposition<Mat> cod(dFh);
        cod.setThreshold(1e-8);
        const Vec gamma = cod.solve(f);
        const Mat G = dXh + beta_ * dFh;
        step -= G * gamma;
      }
      if (std::real(step.dot(f)) <= 0) {
        step = beta_ * f;
        n_hist = 0;
      }
      const value_type s_max = step.cwiseAbs().maxCoeff();
      if (s_max > max_log_step)
        step *= max_log_step / s_max;
      const Eigen::Array<T, -1, 1> phase = map_delta_orb / amp_map;
      mean_delta = (x + step).real().array().exp().template cast<T>() * phase;
      h.pr.broadcast_s(mean_delta, h.bdg.s_delta);
      hist.row(itr + 1) = mean_delta.real().transpose() * energy_scale;
#pragma omp barrier
    }
    store_s_wave_anderson(num_itr, mean_delta, residual, map_hist, hist);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_s_wave_anderson(
  const unsigned total_steps_,
  const Eigen::Array<T, -1, 1> &mean_delta_,
  const Eigen::Array<value_type, -1, 1> &residual_,
  const Eigen::Array<value_type, -1, -1> &map_hist_,
  const Eigen::Array<value_type, -1, -1> &hist_
)
  requires Complex<T>
{
  debug_message("Entered store_s_wave_anderson\n");
#pragma omp master
  if (kmpi::is_root()) {
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/s_wave_anderson/";
    std::string ng = base_grp + "Hist";
    write_hdf5(hist_, &file, ng);

    ng = base_grp + "MapHist";
    write_hdf5(map_hist_, &file, ng);

    const Eigen::Array<value_type, -1, -1> res = residual_;
    ng = base_grp + "Residual";
    write_hdf5(res, &file, ng);

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

    Eigen::Array<value_type, -1, -1> ranks(1, 1);
    ranks(0, 0) = static_cast<value_type>(kmpi::size());
    ng = base_grp + "NumRanks";
    write_hdf5(ranks, &file, ng);
  }
#pragma omp master
  kmpi::barrier();
#pragma omp barrier
  debug_message("Left store_s_wave_anderson\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
