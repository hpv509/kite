#include "ComplexTraits.hpp"

template <typename T>
T ComplexTraits<T>::assign_value(double x, double y)
{
  if constexpr (Complex<T>) {
    using value_type = typename T::value_type;
    return T(static_cast<value_type>(x), static_cast<value_type>(y));
  } else
    return static_cast<T>(x);
}

template <typename T>
T ComplexTraits<T>::myconj(const T x)
{
  if constexpr (Complex<T>)
    return std::conj(x);
  else
    return x;
}

template <typename T>
T ComplexTraits<T>::multEiphase(double phase)
{
  if constexpr (Complex<T>)
    return std::exp(T(0.0, phase));
  else
    return static_cast<T>(1.0);
}

template <typename T>
T ComplexTraits<T>::aux_wr(std::size_t x)
{
  if constexpr (Complex<T>) {
    using value_type = typename T::value_type;
    return T(static_cast<value_type>(x), static_cast<value_type>(2 * x));
  } else
    return static_cast<T>(x);
}

template class ComplexTraits<double>;
template class ComplexTraits<float>;
template class ComplexTraits<long double>;
template class ComplexTraits<std::complex<double>>;
template class ComplexTraits<std::complex<float>>;
template class ComplexTraits<std::complex<long double>>;
