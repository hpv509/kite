#ifndef COEFFICIENTS_H_
#define COEFFICIENTS_H_
#include <Eigen/Dense>

namespace Coefficients {
template <typename T>
T jackson(const int n_, const int polynomials_);

template <typename T>
T gauss_first(const int n_, const T mu_, const T sigma_);

template <typename T>
T gauss_second(const int n_, const T mu_, const T sigma_);

template <typename T>
Eigen::Array<T, -1, 1> build_gaussian(const T energy_, const T width_);

template <typename T>
Eigen::Array<T, -1, 1>
build_window(const T center_, const T width_, const T mult_ = 64);

template <typename T>
Eigen::Array<std::complex<T>, -1, 1> build_cplx_exp(const T t);

template <typename T>
Eigen::Array<T, -1, 1> build_fermi(const T beta_, const T mu_);

template <typename T>
Eigen::Array<T, -1, 1> build_fermi_sqrt(const T beta_, const T mu_);

template <typename T>
Eigen::Array<std::complex<T>, -1, 1> build_dgreen(const std::complex<T>);
} // namespace Coefficients

#endif
