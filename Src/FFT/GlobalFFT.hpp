#ifndef GLOBAL_FFT_H_
#define GLOBAL_FFT_H_
#include "ComplexTraits.hpp"

template <Real T>
struct GlobalFFT {
  using cplx = std::complex<T>;
  cplx *in = nullptr;
  cplx *out = nullptr;
  GlobalFFT() = default;
  ~GlobalFFT();
  void allocate(const std::size_t total_size);
};

#endif
