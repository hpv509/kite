#ifndef LOCAL_UNITARY_H_
#define LOCAL_UNITARY_H_
#include "BondSweep.hpp"

template <typename T, unsigned D> class LocalUnitary;
template <typename T> class LocalUnitary<T, 2> : public BondSweep<T, 2> {
  using Base = BondSweep<T, 2>;
  using Base::r;
  using typename Base::value_type;

public:
  template <typename Vec>
  LocalUnitary(LatticeStructure<2u> &r_, const Vec &v) : Base(r_, v)
  {}

  template <int S, typename Vec>
  void onsite(const T gamma_, Vec &&state_)
    requires Complex<T>
  {
    static_assert(S == -1 || S == 1);
    const value_type norm = 1 / std::sqrt(value_type(2));
    const T g = value_type(S) * gamma_;
    for (unsigned io = 0; io < r.Orb; ++io)
      for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
        const std::size_t p0 = io * r.Nd + i1 * r.Ld[0] + NGHOSTS;
        for (std::size_t p = p0; p < p0 + r.ld[0]; ++p) {
          const T pp = state_.coeff(p), hh = state_.coeff(p + r.offset);
          state_.coeffRef(p) = norm * (pp + std::conj(g) * hh);
          state_.coeffRef(p + r.offset) = norm * (-g * pp + hh);
        }
      }
#pragma omp barrier
  }

  template <int S, typename Vec>
  void pairs(
    const BondTable<2u> &bt,
    const T gamma_,
    const unsigned io,
    const unsigned b,
    const int colour,
    Vec &&state_
  )
    requires Complex<T>
  {
    static_assert(S == -1 || S == 1);
    constexpr value_type norm = 1 / std::sqrt(2.0);
    const T g = value_type(S) * gamma_;
    this->template for_pairs<
      0>(bt, io, b, colour, [&](std::size_t a, std::size_t c, T ph, auto, auto) {
      const T x = state_.coeff(a), y = state_.coeff(c);
      state_.coeffRef(a) = norm * (x + std::conj(g) * ph * y);
      state_.coeffRef(c) = norm * (-g * std::conj(ph) * x + y);
    });
#pragma omp barrier
  }
};

#endif
