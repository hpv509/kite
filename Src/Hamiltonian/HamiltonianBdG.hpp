#ifndef HAMILTONIAN_BDG_
#define HAMILTONIAN_BDG_
#include "ComplexTraits.hpp"
#include "Eigen/Dense"

template <Scalar T, unsigned D>
struct HamiltonianBdG {
  using real = typename extract_scalar<T>::type;
  Eigen::Array<real, -1, 1> hartree;
  Eigen::Array<T, -1, 1> s_delta;
  Eigen::Array<T, -1, -1> nn_delta;
  HamiltonianBdG(const std::size_t, const std::size_t max_hoppings_);
  void init_fields(const T s_delta_, const real ht_, const real mu_);
  void init_fields(const T s_delta_, const T nn_delta_, const real ht_, const real mu_);
};

#endif
