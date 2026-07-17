#include "GlobalFFT.hpp"
#include "TraitsFFTW.hpp"

template <Real T>
void GlobalFFT<T>::allocate(const std::size_t total_size)
{
  if (!in) {
    in = reinterpret_cast<cplx *>(fftw_malloc(sizeof(cplx) * total_size));
    out = reinterpret_cast<cplx *>(fftw_malloc(sizeof(cplx) * total_size));
  }
}

template <Real T>
GlobalFFT<T>::~GlobalFFT()
{
  if (in)
    fftw_free(in);
  if (out)
    fftw_free(out);
}

#define INSTANTIATE_GLOBAL_FFT(T) template class GlobalFFT<T>;
INSTANTIATE_GLOBAL_FFT(float)
INSTANTIATE_GLOBAL_FFT(double)
INSTANTIATE_GLOBAL_FFT(long double)
