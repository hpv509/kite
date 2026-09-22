#ifndef SIGMOID_WEIGHT_H_
#define SIGMOID_WEIGHT_H_

template <typename T>
inline T sig_weight_ratio(const T n_, const T N0_, const T tau_)
{
  if (std::isinf(N0_))
    return std::exp(1.0 / tau_);
  auto softplus = [](const T z_) {
    return std::max(z_, T(0.0)) + std::log1p(std::exp(-std::abs(z_)));
  };
  return std::exp(
    softplus((N0_ - n_ + 1.0) / tau_) - softplus((N0_ - n_) / tau_)
  );
}
#endif
