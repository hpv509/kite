#include "HamiltonianBdG.hpp"

template <Scalar T, unsigned D>
HamiltonianBdG<T, D>::HamiltonianBdG(char *name, const LatticeStructure<D> &r)
{
  real mu_read = 0, beta_read = 0;
  const std::string base_dir = "/Hamiltonian/BdG/";
  std::string tmp;
  Eigen::Array<real, -1, 1> hartree_orb(r.Orb);
  hartree_orb.setZero();
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
      get_hdf5<real>(&beta_read, &file, tmp);
      tmp = base_dir + "Hartree";
      get_hdf5<real>(hartree_orb.data(), &file, tmp);
    } catch (H5::Exception &e) {
    }
    file.close();
  };
#pragma omp barrier
  hartree = Eigen::Array<real, -1, 1>::Zero(r.Sized);
  onsite = Eigen::Array<real, -1, 1>::Zero(r.Sized);
  s_delta = Eigen::Array<T, -1, 1>::Zero(r.Sized);

  const std::size_t block = r.Sized / r.Orb;
  for (unsigned io = 0; io < r.Orb; ++io)
    hartree.segment(io * block, block).setConstant(hartree_orb(io));

  set_beta(beta_read);
  set_chemical_potential(mu_read);
}

#define instantiate(type, dim) template struct HamiltonianBdG<type, dim>;
#include "instantiate.hpp"
