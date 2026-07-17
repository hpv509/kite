#ifndef FFT_H_
#define FFT_H_
#include <fftw3.h>
#include <array>
#include "Generic.hpp"
#include <Eigen/Dense>
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "GlobalFFT.hpp"
#include "TraitsFFTW.hpp"

template <Real T, unsigned D>
struct FFT {
  LatticeStructure<D> &r;
  GlobalFFT<T> &global_fft;
  const Eigen::Matrix<T, D, D> b_rcp;
  using cplx = std::complex<T>;
  typename FFTWTraits<T>::cplx *global_in;
  typename FFTWTraits<T>::cplx *global_out;
  typename FFTWTraits<T>::plan fwd_plan;
  typename FFTWTraits<T>::plan bwd_plan;

  FFT(LatticeStructure<D> &, GlobalFFT<T> &);
  ~FFT();
  void forward(Eigen::Array<cplx, -1, 1> &);
  void inverse(Eigen::Array<cplx, -1, 1> &);
};

#endif
