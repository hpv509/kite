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
  Eigen::Array<real, -1, 1> free_energy;
  Eigen::Array<real, -1, 1> onsite;
  real mu{0};

  HamiltonianBdG(const std::size_t, const std::size_t max_hoppings_);
  void update_onsite() { onsite = hartree - mu; }
  void set_chemical_potential(const real mu_)
  {
    mu = mu_;
    update_onsite();
  }
  void init_fields(
    const Eigen::Array<real, -1, 1> &s_delta_,
    const Eigen::Array<real, -1, 1> &ht_
  );
  void init_fields(const T s_delta_, const T nn_delta_, const real ht_);
};

#endif
