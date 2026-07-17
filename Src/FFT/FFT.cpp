#include "FFT.hpp"
#include "Loop.hpp"

template <Real T, unsigned D>
FFT<T, D>::FFT(LatticeStructure<D> &r_, GlobalFFT<T> &global_fft_) :
  r(r_),
  global_fft(global_fft_),
  b_rcp((2.0 * M_PI * r_.rLat.inverse().transpose()).template cast<T>())
{
  int num_cells = 1;
  int rank_dims[D];
  for (unsigned i = 0; i < D; ++i) {
    num_cells *= r_.Lt[i];
    rank_dims[i] = r_.Lt[D - 1 - i];
  }
  global_in = reinterpret_cast<typename FFTWTraits<T>::cplx *>(global_fft.in);
  global_out = reinterpret_cast<typename FFTWTraits<T>::cplx *>(global_fft.out);
#pragma omp critical
  {
    fwd_plan = FFTWTraits<T>::plan_many_dft(
      static_cast<int>(D), rank_dims, static_cast<int>(r.Orb), global_in,
      nullptr, 1, num_cells, global_out, nullptr, 1, num_cells, FFTW_FORWARD,
      FFTW_ESTIMATE
    );
    bwd_plan = FFTWTraits<T>::plan_many_dft(
      static_cast<int>(D), rank_dims, static_cast<int>(r.Orb), global_in,
      nullptr, 1, num_cells, global_out, nullptr, 1, num_cells, FFTW_BACKWARD,
      FFTW_ESTIMATE
    );
  }
#pragma omp barrier
}

template <Real T, unsigned D>
FFT<T, D>::~FFT(void)
{
#pragma omp critical
  {
    if (fwd_plan)
      FFTWTraits<T>::destroy_plan(fwd_plan);
    if (bwd_plan)
      FFTWTraits<T>::destroy_plan(bwd_plan);
  }
#pragma omp barrier
}

template <Real T, unsigned D>
void FFT<T, D>::forward(Eigen::Array<cplx, -1, 1> &v_)
{
  Coordinates<std::size_t, D + 1> local(r.Ld);
  Coordinates<std::size_t, D + 1> global(r.Lt);
  cplx *vec_in = reinterpret_cast<cplx *>(global_in);
  cplx *vec_out = reinterpret_cast<cplx *>(global_out);

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
      vec_in[global.index] = v_(local.index);
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  FFTWTraits<T>::execute(fwd_plan);
#pragma omp barrier
  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    const Eigen::Matrix<T, D, 1> tau = r.rOrb.col(io).template cast<T>();

    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2)
        local.set({i[1], i[0], io});
      else if constexpr (D == 3)
        local.set({i[2], i[1], i[0], io});
      r.convertCoordinates(global, local);

      Eigen::Matrix<T, D, 1> k_fraction;
      for (unsigned i = 0; i < D; ++i) {
        const int g_idx = global.coord[i];
        const int k_idx = (g_idx < r.Lt[i] / 2) ? g_idx : g_idx - r.Lt[i];
        k_fraction(i) = k_idx / r.Lt[i];
      }
      const Eigen::Matrix<T, D, 1> k = b_rcp * k_fraction;
      const T arg = tau.dot(k);
      const cplx phase(std::cos(arg), -std::sin(arg));
      v_(local.index) = vec_out[global.index] * phase;
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
}

template <Real T, unsigned D>
void FFT<T, D>::inverse(Eigen::Array<cplx, -1, 1> &v_)
{
  Coordinates<std::size_t, D + 1> local(r.Ld);
  Coordinates<std::size_t, D + 1> global(r.Lt);
  std::complex<T> *vec_in = reinterpret_cast<std::complex<T> *>(global_in);
  std::complex<T> *vec_out = reinterpret_cast<std::complex<T> *>(global_out);

  std::array<unsigned, D> idx;
  std::array<unsigned, D> start;
  std::array<unsigned, D> final;
  T norm_factor = 1.0;
  for (unsigned d = 0; d < D; ++d) {
    start[d] = NGHOSTS;
    final[d] = r.Ld[D - 1 - d] - NGHOSTS;
    norm_factor *= r.Lt[d];
  }
  norm_factor = 1.0 / norm_factor;

  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    const Eigen::Matrix<T, D, 1> tau = r.rOrb.col(io).template cast<T>();

    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2)
        local.set({i[1], i[0], io});
      else if constexpr (D == 3)
        local.set({i[2], i[1], i[0], io});
      r.convertCoordinates(global, local);

      Eigen::Matrix<T, D, 1> k_fraction;
      for (unsigned d = 0; d < D; ++d) {
        const int g_idx = global.coord[d];
        const int k_idx = (g_idx < r.Lt[d] / 2) ? g_idx : g_idx - r.Lt[d];
        k_fraction(d) = k_idx / r.Lt[d];
      }
      const Eigen::Matrix<T, D, 1> k = b_rcp * k_fraction;
      const T arg = tau.dot(k);
      const std::complex<T> inv_phase(std::cos(arg), std::sin(arg));
      vec_in[global.index] = v_(local.index) * inv_phase;
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  FFTWTraits<T>::execute(bwd_plan);
#pragma omp barrier
  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2)
        local.set({i[1], i[0], io});
      else if constexpr (D == 3)
        local.set({i[2], i[1], i[0], io});
      r.convertCoordinates(global, local);
      v_(local.index) = vec_out[global.index] * norm_factor;
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
}

#define instantiate_real(type, dim) template class FFT<type, dim>;
#include "instantiate.hpp"
