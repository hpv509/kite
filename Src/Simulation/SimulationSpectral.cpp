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
#include "FFT/FFT.hpp"

template <typename T, unsigned D>
void Simulation<T, D>::calc_spectral(GlobalFFT<value_type> &global_fft_)
{
  debug_message("Entered Simulation::calc_spectral\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    Global.calculate_spectral_map = false;
    std::string tmp = "/Calculation/spectral_map/NumVectors";
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, &file, tmp);
      Global.calculate_spectral_map = true;
    } catch (H5::Exception &e) {
      debug_message("spectral map: no need to calculate it.\n");
    }
    file.close();
  }
#pragma omp barrier
  bool local_calculate_spectral_map = false;
#pragma omp critical
  local_calculate_spectral_map = Global.calculate_spectral_map;
#pragma omp barrier
  if (local_calculate_spectral_map) {
#pragma omp master
    std::cout << "Calculating Spectral Function.\n";
#pragma omp barrier
    int vectors;
    value_type energy;
    value_type sigma;
    int coef_id;
#pragma omp critical
    {
      H5::H5File file(name, H5F_ACC_RDONLY);
      const std::string base = "/Calculation/spectral_map/";
      std::string tmp;
      tmp = base + "NumVectors";
      get_hdf5<int>(&vectors, &file, tmp);
      tmp = base + "Energy";
      get_hdf5<value_type>(&energy, &file, tmp);
      tmp = base + "Sigma";
      get_hdf5<value_type>(&sigma, &file, tmp);
      tmp = base + "Coef_ID";
      get_hdf5<int>(&coef_id, &file, tmp);
      file.close();
    }
    spectral(vectors, energy, sigma, coef_id, global_fft_);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::spectral(
  const unsigned vectors_,
  const value_type energy_,
  const value_type sigma_,
  const unsigned coef_id_,
  GlobalFFT<value_type> &global_fft_
)
{
  debug_message("Entered ldos\n");
  if constexpr (is_tt<std::complex, T>::value) {
    value_type energy_scale;
#pragma omp critical
    {
      H5::H5File file(name, H5F_ACC_RDONLY);
      std::string tmp = "/EnergyScale";
      get_hdf5<value_type>(&energy_scale, &file, tmp);
      file.close();
    }
#pragma omp barrier
    const value_type target = energy_ / energy_scale;
    const value_type sigma = sigma_ / energy_scale;
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type factor =
      (coef_id_) ? 1.0 : std::sqrt(8 * M_PI) * sigma / energy_scale;
    const value_type fwhm = (coef_id_) ? sigma : std::sqrt(2) * sigma;
    const Eigen::Array<value_type, -1, 1> coefs =
      (coef_id_) ? Coefficients::build_window<value_type>(target, fwhm)
                 : Coefficients::build_gaussian<value_type>(target, fwhm);

    KPM_Vector<T, D> phi(2, *this);
    Eigen::Array<T, -1, 1> ket(r.Sized);
    Eigen::Array<T, -1, 1> bra(r.Sized);
    Eigen::Array<value_type, -1, -1> results(r.Sized, 2);
    results.setZero();
    Eigen::Array<value_type, -1, 1> prv(r.Sized);
    FFT<value_type, D> fft(r, global_fft_);

    h.generate_disorder();
    for (int vec = 0; vec < vectors_; ++vec) {
      h.generate_twists();
      phi.initiate_phases();
      phi.set_index(0);
      phi.initiate_vector();
      phi.v.col(0) *= std::sqrt(size);
      bra = phi.v.col(0);

      ket = bra;
      fft.inverse(ket);
      h.hV.erase_wavefunction(ket);
      phi.v.col(0) = ket;
      phi.Exchange_Boundaries();
      ket.setZero();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      fft.forward(ket);
      const value_type weight = 1.0 / (vec + 1);
      const Eigen::Array<value_type, -1, 1> map =
        factor * (bra.conjugate() * ket).abs2();

      prv = results.col(0);
      results.col(0) += weight * (map - results.col(0));
      results.col(1) +=
        weight * ((map - prv) * (map - results.col(0)) - results.col(1));
    }
    results.col(1) = results.col(1).sqrt() / std::sqrt(vectors_);
    store_spectral(results);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_spectral(
  const Eigen::Array<value_type, -1, -1> &results_
)
{
  debug_message("Entered store_spectral\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  Global.ldos_map.resize(r.Sizet, 2);
#pragma omp barrier
  std::array<unsigned, D> idx;
  std::array<unsigned, D> start;
  std::array<unsigned, D> final;
  for (unsigned d = 0; d < D; ++d) {
    start[d] = NGHOSTS;
    final[d] = r.Ld[D - 1 - d] - NGHOSTS;
  }
  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2)
        local.set({i[1], i[0], io});
      else if constexpr (D == 3)
        local.set({i[2], i[1], i[0], io});
      r.convertCoordinates(global, local);

      std::array<unsigned, D> shifted_coord;
      for (unsigned d = 0; d < D; ++d) {
        const unsigned N = r.Lt[d];
        shifted_coord[d] = (global.coord[d] + N / 2) % N;
      }
      unsigned global_idx = 0;
      unsigned stride = 1;
      for (int d = 0; d < D; ++d) {
        global_idx += shifted_coord[d] * stride;
        stride *= r.Lt[d];
      }
      global_idx += io * stride;
      Global.ldos_map.row(global_idx) = results_.row(local.index);
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  {
    const Eigen::Array<value_type, -1, -1> ldos_r = Global.ldos_map.real();
    H5::H5File file(name, H5F_ACC_RDWR);
    const std::string ng = "/Calculation/spectral_map/Map";
    write_hdf5(ldos_r, &file, ng);
  }
#pragma omp barrier
  debug_message("Left store_ldos\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
