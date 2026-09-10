#ifndef COND_H_
#define COND_H_

#include "coefficients.hpp"
#include <span>

template <typename DerivedY, typename DerivedX>
auto simpson_integrate(
  const Eigen::ArrayBase<DerivedY> &,
  const Eigen::ArrayBase<DerivedX> &
);

arr<type, -1, 1>
fermi_function(const arr<type, -1, 1> &, const type, const type);

mtx<type, -1, 1> calculate_conductivity(
  const std::string_view,
  const std::string_view,
  const std::span<const type>,
  const type,
  const std::span<const type>,
  const type,
  const type
);

#endif
