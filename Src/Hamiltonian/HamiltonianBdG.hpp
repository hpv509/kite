#ifndef HAMILTONIAN_BDG_
#define HAMILTONIAN_BDG_
#include "ComplexTraits.hpp"
#include "Eigen/Dense"

template <Scalar T, unsigned D>
struct HamiltonianBdG {
  using real = typename extract_scalar<T>::type;
  unsigned is_bdg;
  Eigen::Array<real, -1, 1> hartree;
  Eigen::Array<real, -1, 1> s_delta;
  HamiltonianBdG(const std::size_t);
  void init_fields(const real, const real);
};

#endif
