#include "HamiltonianBdG.hpp"

template <Scalar T, unsigned D>
HamiltonianBdG<T, D>::HamiltonianBdG(char *name, const LatticeStructure<D> &r)
{
  real mu_read = 0, beta_read = 0;
  const std::string base_dir = "/Hamiltonian/BdG/";
  std::string tmp;
#pragma omp critical
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    tmp = "/EnergyScale";
    get_hdf5<real>(&energy_scale, &file, tmp);
    try {
      H5::Exception::dontPrint();
      tmp = base_dir + "ChemicalPotential";
      get_hdf5<real>(&mu_read, &file, tmp);
      tmp = base_dir + "Beta";
      get_hdf5<real>(&beta_read, &file, (char *)"/Hamiltonian/BdG/Beta");
    } catch (H5::Exception &e) {
    }
    file.close();
  }

  hartree = Eigen::Array<real, -1, 1>::Zero(r.Sized);
  onsite = Eigen::Array<real, -1, 1>::Zero(r.Sized);
  s_delta = Eigen::Array<T, -1, 1>::Zero(r.Sized);

  set_beta(beta_read);
  set_chemical_potential(mu_read);
}

template <Scalar T, unsigned D>
void HamiltonianBdG<T, D>::init_sw(
  const Eigen::Array<real, -1, 1> &s_delta_,
  const Eigen::Array<real, -1, 1> &ht_
)
{
  const unsigned norb = s_delta_.size();
  const unsigned block = s_delta.size() / norb;
  for (unsigned i = 0; i < norb; ++i) {
    hartree.segment(i * block, block) = ht_(i);
    s_delta.segment(i * block, block) = static_cast<T>(s_delta_(i));
  }
  update_onsite();
}

#define instantiate(type, dim) template struct HamiltonianBdG<type, dim>;
#include "instantiate.hpp"
