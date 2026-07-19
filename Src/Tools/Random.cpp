#include <algorithm>
#include <array>
#include "Random.hpp"

template <Scalar T>
KPMRandom<T>::KPMRandom(const unsigned seed_)
{
  init_random(seed_);
}

template <Scalar T>
void KPMRandom<T>::init_random(const unsigned seed_)
{
  if (seed_ > 0)
    rng.seed(seed_);
  else {
    std::random_device r;
    std::array<int, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    rng.seed(seq);
  }
}

template <Scalar T>
auto KPMRandom<T>::get() -> real
{
  return dist(rng);
}

template <Scalar T>
auto KPMRandom<T>::uniform(const real mean, const real width) -> real
{
  return mean + std::sqrt(3.0) * width * (2 * dist(rng) - 1.0);
}

template <Scalar T>
auto KPMRandom<T>::gaussian(const real mean, const real width) -> real
{
  return mean + width * gauss(rng);
}

template <Scalar T>
T KPMRandom<T>::init()
{
  if constexpr (Complex<T>)
    return std::exp(T(0.0, 2 * M_PI * dist(rng)));
  else
    return (2 * dist(rng) - 1.0) * std::sqrt(3.0);
}

template class KPMRandom<float>;
template class KPMRandom<double>;
template class KPMRandom<long double>;
template class KPMRandom<std::complex<float>>;
template class KPMRandom<std::complex<double>>;
template class KPMRandom<std::complex<long double>>;
