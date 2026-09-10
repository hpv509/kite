#ifndef CONSTRAINTS_H_
#define CONSTRAINTS_H_
#include <concepts>
#if __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#elif __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#endif

template <typename T>
concept Real = std::is_floating_point_v<T>;

template <typename T>
concept Complex = requires {
  typename T::value_type;
  requires std::floating_point<typename T::value_type>;
};

template <typename T>
concept Scalar = Real<T> || Complex<T>;

template <Scalar T, int Rows, int Cols>
using arr = Eigen::Array<T, Rows, Cols>;
template <Scalar T, int Rows, int Cols>
using mtx = Eigen::Matrix<T, Rows, Cols>;

using type = double;
using cplx = std::complex<type>;

#endif
