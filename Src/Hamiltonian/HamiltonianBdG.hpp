#ifndef HAMILTONIAN_BDG_
#define HAMILTONIAN_BDG_
#include "Generic.hpp"
#include "myHDF5.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "Pairing.hpp"

template <Scalar T, unsigned D>
struct HamiltonianBdG {
  using real = typename extract_scalar<T>::type;
  Eigen::Array<real, -1, 1> hartree;
  Eigen::Array<T, -1, 1> s_delta;
  Eigen::Array<T, -1, -1> nn_delta;
  Eigen::Array<real, -1, 1> free_energy;
  Eigen::Array<real, -1, 1> onsite;
  real energy_scale;
  real mu;
  real beta;

  HamiltonianBdG(char *, const LatticeStructure<D> &);
  void update_onsite() { onsite = hartree - mu; }
  void set_chemical_potential(real mu_rescaled)
  {
    mu = mu_rescaled;
    update_onsite();
  }

  void set_beta(real beta_rescaled) { beta = beta_rescaled; }

  void set_hartree(const Eigen::Array<real, -1, 1> &h_phys)
  {
    hartree = h_phys / energy_scale;
    update_onsite();
  }
  void
  init_sw(const Eigen::Array<real, -1, 1> &, const Eigen::Array<real, -1, 1> &);
};

#endif
