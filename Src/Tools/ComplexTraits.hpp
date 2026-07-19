#ifndef TRAITS_H_
#define TRAITS_H_
#include <complex>

template <typename T>
struct extract_scalar {
  using type = T;
};

template <typename T>
struct extract_scalar<std::complex<T>> {
  using type = T;
};

template <typename T>
concept Real = std::is_floating_point_v<T>;

template <typename T>
concept Complex = requires {
  typename T::value_type;
  requires std::floating_point<typename T::value_type>;
};

template <typename T>
concept Scalar = Real<T> || Complex<T>;

template <typename T>
class ComplexTraits {
public:
  T myconj(T &x);
  T assign_value(double x, double y);
  T multEiphase(double phase);
  T aux_wr(std::size_t x);
};

#endif
