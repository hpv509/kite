#ifndef COEFFICIENTS_H_
#define COEFFICIENTS_H_
#include "constraints.hpp"

struct Coefficients {
  arr<type, -1, 1>
  gauss_first(const type n, const arr<type, -1, 1> &, const type);
  arr<type, -1, 1>
  gauss_second(const type n, const arr<type, -1, 1> &, const type);

  mtx<type, -1, -1> build_gaussian(
    const mtx<type, -1, 1> &,
    const type,
    const type,
    const unsigned
  );
  mtx<cplx, -1, -1> build_dgreen(const mtx<cplx, -1, 1> &, const type, const unsigned);
};

namespace cst {
    constexpr cplx I{0.0, 1.0};
}

#endif
