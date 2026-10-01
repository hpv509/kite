#ifndef LOCAL_CURRENT_H_
#define LOCAL_CURRENT_H_
#include "BondSweep.hpp"

template <typename T, unsigned D> class LocalCurrent;
template <typename T> class LocalCurrent<T, 2> : public BondSweep<T, 2> {
  using Base = BondSweep<T, 2>;
  using typename Base::value_type;

public:
  using Base::Base;
  template <typename Vec, typename Sink>
  void sample(
    const BondTable<2u> &bt,
    const T w,
    const unsigned io,
    const unsigned b,
    const int colour,
    const Vec &psi,
    Sink &&sink
  ) const
    requires Complex<T>
  {
    this->template for_pairs<1>(
      bt, io, b, colour,
      [&](std::size_t a, std::size_t c, T ph, std::size_t g0, std::size_t g1) {
        sink(
          g0, g1,
          value_type(2) *
            std::real(w * ph * std::conj(psi.coeff(a)) * psi.coeff(c))
        );
      }
    );
  }
};

#endif
