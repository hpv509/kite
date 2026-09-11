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
void Simulation<T, D>::calc_pwave_clean()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_pwave_clean\n");
  std::string base_grp = "/Calculation/p_wave_c/";
  std::string tmp = base_grp + "NumRandoms";
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
    int randoms;
    int num_itr;
    int prv_itr;
    value_type weight_r;
    value_type weight_alpha;

    const unsigned max_bonds = h.pr.max_bonds;
    const unsigned n_bonds = max_bonds * r.Orb;
    Eigen::Array<T, -1, -1> sum_delta_init(max_bonds, r.Orb);
    value_type weight_sum_init = 0.0;
    value_type weight_avg_init = 1.0;
    sum_delta_init.setZero();

#pragma omp critical
    {
      H5::H5File file(name, H5F_ACC_RDONLY);
      std::string path = base_grp + "NumRandoms";
      get_hdf5<int>(&randoms, &file, path);
      path = base_grp + "NumIterations";
      get_hdf5<int>(&num_itr, &file, path);
      path = base_grp + "PrevIterations";
      get_hdf5<int>(&prv_itr, &file, path);
      path = base_grp + "Wr";
      get_hdf5<value_type>(&weight_r, &file, path);
      path = base_grp + "Walpha";
      get_hdf5<value_type>(&weight_alpha, &file, path);

      if (prv_itr > 0) {
        try {
          Eigen::Array<value_type, -1, -1> sum_delta_ri(2, n_bonds);
          path = base_grp + "SumDelta";
          get_hdf5<value_type>(sum_delta_ri.data(), &file, path);
          Eigen::Array<value_type, -1, -1> weights(2, 1);
          path = base_grp + "Weights";
          get_hdf5<value_type>(weights.data(), &file, path);

          for (unsigned io = 0; io < r.Orb; ++io)
            for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b) {
              const std::size_t k = b + io * max_bonds;
              sum_delta_init(b, io) = T(sum_delta_ri(0, k), sum_delta_ri(1, k));
            }
          weight_sum_init = weights(0, 0);
          weight_avg_init = weights(1, 0);
        } catch (H5::Exception &e) {
        }
      }
      file.close();
    }
    p_wave_clean(
      randoms, num_itr, prv_itr, weight_r, weight_alpha, weight_sum_init,
      weight_avg_init, sum_delta_init
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::p_wave_clean(
  const int randoms_,
  const int num_itr_,
  const int prv_itr_,
  const value_type weight_r_,
  const value_type weight_alpha_,
  const value_type weight_sum_init_,
  const value_type weight_avg_init_,
  const Eigen::Array<T, -1, -1> &sum_delta_init_
)
  requires Complex<T>
{
  debug_message("Entered PWaveClean\n");
  const value_type energy_scale = h.bdg.energy_scale;
  Coordinates<std::size_t, D + 1> local(r.Ld);
  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(h.bdg.beta, 0.0);

  const value_type size = r.Sizet - r.SizetVacancies;
  const value_type n_cells = r.Nt;
  const unsigned max_bonds = h.pr.max_bonds;
  const unsigned n_bonds = max_bonds * r.Orb;
  Eigen::Array<T, -1, -1> mean_delta(max_bonds, r.Orb);
  Eigen::Array<T, -1, -1> sum_delta(max_bonds, r.Orb);
  Eigen::Array<T, -1, -1> local_delta(max_bonds, r.Orb);

  value_type weight_avg = 1.0;
  value_type weight_sum = 0.0;
  sum_delta = sum_delta_init_;
  weight_sum = weight_sum_init_;
  weight_avg = weight_avg_init_;
  if (weight_sum > 0)
    mean_delta = sum_delta / weight_sum;
  else
    mean_delta = h.pr.Delta0;
  h.pr.broadcast(mean_delta, h.bdg.nn_delta);
#pragma omp master
  {
    Global.nn_sum.resize(n_bonds);
    Global.nn_delta_hist.resize(num_itr_ + 1, n_bonds);
    Global.nn_delta_hist.setZero();
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b)
        Global.nn_delta_hist(0, b + io * max_bonds) =
          mean_delta(b, io) * energy_scale;
  }
#pragma omp barrier
  h.generate_disorder();
  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(2 * r.Sized);
  Eigen::Array<T, -1, 1> bra(2 * r.Sized);

  for (unsigned itr = prv_itr_ + 1; itr <= prv_itr_ + num_itr_; ++itr) {
    local_delta.setZero();
    for (int vec = 0; vec < randoms_; ++vec) {
      h.generate_twists();
      phi.initiate_phases();

      phi.set_index(0);
      phi.initiate_vector();
      phi.v.col(0) *= std::sqrt(size);
      phi.empty_ghosts(0);
      bra = phi.v.col(0);
      const value_type weight = 1.0 / (vec + 1);

      for (unsigned io = 0; io < r.Orb; ++io) {
        for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b) {
          phi.v.setZero();
          phi.set_index(0);
          phi.v.col(0) = bra;
          // nn_pairing reads its partner out of the ghost cells
          phi.Exchange_Boundaries();
          phi.template nn_pairing<-1>(1.0, io, b, phi.v.col(0));

          ket.setZero();
          phi.Exchange_Boundaries();
          for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
            phi.cheb_iteration(n);
            ket += coefs(n) * phi.v.col(phi.get_index()).array();
          }
          phi.template nn_pairing<1>(1.0, io, b, ket);

          const Eigen::Array<value_type, -1, 1> ups =
            (bra.conjugate() * ket).abs2();
          const Eigen::Array<value_type, -1, 1> ups_p = ups.head(r.Sized);
          const Eigen::Array<value_type, -1, 1> ups_h = ups.tail(r.Sized);

          value_type acc = 0;
          std::size_t wb[D], we[D];

          phi.template pair_window<1>(io, b, wb, we);
          for (std::size_t i1 = wb[1]; i1 < we[1]; ++i1) {
            local.set({wb[0], i1, io});
            const std::size_t j0 = local.index;
            for (std::size_t i = j0, j1 = j0 + (we[0] - wb[0]); i < j1; ++i)
              acc += ups_p(i);
          }
          const std::size_t jo = h.pr.target_orb(b, io);
          phi.template pair_window<-1>(io, b, wb, we);
          for (std::size_t i1 = wb[1]; i1 < we[1]; ++i1) {
            local.set({wb[0], i1, jo});
            const std::size_t j0 = local.index;
            for (std::size_t j = j0, j1 = j0 + (we[0] - wb[0]); j < j1; ++j)
              acc -= ups_h(j);
          }
          acc *= 0.5;

          const T map_delta = h.pr.V * acc;
          local_delta(b, io) += weight * (map_delta - local_delta(b, io));
        }
      }
    }
#pragma omp master
    Global.nn_sum.setZero();
#pragma omp barrier
#pragma omp critical
    {
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b)
          Global.nn_sum(b + io * max_bonds) += local_delta(b, io);
    }
#pragma omp barrier
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b)
        mean_delta(b, io) = Global.nn_sum(b + io * max_bonds) / n_cells;

    h.pr.symmetrize_bonds(mean_delta);

    weight_avg *= 1.0 + weight_r_ / std::pow(itr, weight_alpha_);
    weight_sum += weight_avg;
    sum_delta += weight_avg * mean_delta;
    mean_delta = sum_delta / weight_sum;
    h.pr.broadcast(mean_delta, h.bdg.nn_delta);
#pragma omp barrier
#pragma omp master
    {
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b)
          Global.nn_delta_hist(itr - prv_itr_, b + io * max_bonds) =
            mean_delta(b, io) * energy_scale;
    }
#pragma omp barrier
  }
  store_p_wave_clean(prv_itr_ + num_itr_, sum_delta, weight_sum, weight_avg);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_p_wave_clean(
  const unsigned total_steps_,
  const Eigen::Array<T, -1, -1> &sum_delta_,
  const value_type weight_sum_,
  const value_type weight_avg_
)
  requires Complex<T>
{
  debug_message("Entered store_pwave_clean\n");
#pragma omp master
  {
    const Eigen::Array<value_type, -1, -1> p_wave_r =
      Global.nn_delta_hist.real();
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string base_grp = "/Calculation/p_wave_c/";
    std::string ng = base_grp + "Hist";
    write_hdf5(p_wave_r, &file, ng);

    Eigen::Array<value_type, -1, -1> total_steps(1, 1);
    total_steps(0, 0) = static_cast<value_type>(total_steps_);
    ng = base_grp + "TotalSteps";
    write_hdf5(total_steps, &file, ng);

    const unsigned max_bonds = h.pr.max_bonds;
    const unsigned n_bonds = max_bonds * r.Orb;
    Eigen::Array<value_type, -1, -1> sum_delta_ri(2, n_bonds);
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0, B = h.pr.NPairings(io); b < B; ++b) {
        const std::size_t k = b + io * max_bonds;
        sum_delta_ri(0, k) = sum_delta_(b, io).real();
        sum_delta_ri(1, k) = sum_delta_(b, io).imag();
      }
    ng = base_grp + "SumDelta";
    write_hdf5(sum_delta_ri, &file, ng);

    Eigen::Array<value_type, -1, -1> weights(2, 1);
    weights(0, 0) = weight_sum_;
    weights(1, 0) = weight_avg_;
    ng = base_grp + "Weights";
    write_hdf5(weights, &file, ng);
  }
#pragma omp barrier
  debug_message("Left store_pwave_clean\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
