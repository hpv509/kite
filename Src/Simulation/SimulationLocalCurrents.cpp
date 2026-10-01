#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template <typename T, unsigned D> class Hamiltonian;
template <typename T, unsigned D> class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"
#include "BondTable.hpp"
#include "LocalCurrent.hpp"
#include "Coefficients.hpp"

namespace {
const std::string bond_grp = "/Calculation/bond_map/";

[[noreturn]] void bond_map_fail(const std::string &msg)
{
#pragma omp master
  std::cerr << "bond_map: " << msg << "\n";
#pragma omp barrier
  exit(1);
}
} // namespace

template <typename T, unsigned D>
void Simulation<T, D>::calc_bond_map()
  requires Complex<T>
{
  int randoms = 0;
  value_type beta = 0, mu = 0;
  bool present = false, complete = false;
#pragma omp critical
  {
    H5::Exception::dontPrint();
    H5::H5File file(name, H5F_ACC_RDONLY);
    std::string path = bond_grp + "NumRandoms";
    try {
      get_hdf5<int>(&randoms, &file, path);
      present = true;
      path = bond_grp + "Beta";
      get_hdf5<value_type>(&beta, &file, path);
      path = bond_grp + "ChemicalPotential";
      get_hdf5<value_type>(&mu, &file, path);
      complete = true;
    } catch (H5::Exception &) {
    }
    file.close();
  }
#pragma omp barrier
  if (!present)
    return;
  if (!complete || randoms < 1)
    bond_map_fail("needs NumRandoms >= 1, Beta and ChemicalPotential");
  bond_map(randoms, beta, mu);
}

template <typename T, unsigned D>
void Simulation<T, D>::bond_map(
  const int randoms_,
  const value_type beta_,
  const value_type mu_
)
  requires Complex<T>
{
  if constexpr (D != 2) {
    (void)randoms_, (void)beta_, (void)mu_;
    bond_map_fail("the bond kernels exist for D = 2 only");
  } else {
    if (r.MagneticField)
      bond_map_fail("magnetic field: fold the Peierls phase into W first");

    BondTable<D> bt(name, r, bond_grp, 0);
    if (bt.n_bonds == 0)
      bond_map_fail("empty bond table");
    bt.print("bond_map");

    Eigen::Array<T, -1, -1> W =
      Eigen::Array<T, -1, -1>::Zero(bt.max_bonds, r.Orb);
#pragma omp critical
    {
      H5::H5File file(name, H5F_ACC_RDONLY);
      std::string path = bond_grp + "Weights";
      get_hdf5<T>(W.data(), &file, path);
      file.close();
    }
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0; b < bt.NBonds(io); ++b)
        if (W(b, io) == T(0))
          bond_map_fail("every family needs a nonzero weight");

    // sqrt of the Fermi function, rescaled beta and mu; 3 beta moments,
    // truncation bias <= 2e-5 of the largest current (fermi-sqrt-truncation.md)
    const Eigen::Array<value_type, -1, 1> coefs =
      Coefficients::build_fermi_sqrt<value_type>(beta_, mu_);
    const value_type size = r.Sizet - r.SizetVacancies;
    const std::size_t L0 = r.Lt[0];
#pragma omp master
    {
      Global.bond_mean.setZero(bt.n_bonds, r.Nt);
      Global.bond_m2.setZero(bt.n_bonds, r.Nt);
      std::cout << "bond_map: " << coefs.size() << " Chebyshev moments, "
                << randoms_ << " random vectors\n";
    }
#pragma omp barrier

    h.generate_disorder();
    KPM_Vector<T, D> phi(2, *this);
    const LocalCurrent<T, D> currents(r, phi);
    Eigen::Array<T, -1, 1> ket(phi.v.rows());

    for (int vec = 0; vec < randoms_; ++vec) {
      h.generate_twists();
      phi.initiate_phases();
      phi.set_index(0);
      phi.initiate_vector();
      phi.v.col(0) *= std::sqrt(size);
      phi.Exchange_Boundaries();
      ket.setZero();
      for (unsigned n = 0; n < coefs.size(); ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      const value_type k = vec + 1;
      for (unsigned io = 0; io < r.Orb; ++io)
        for (unsigned b = 0; b < bt.NBonds(io); ++b) {
          const unsigned f = bt.family_index(b, io);
          for (int c = 0; c < (bt.coloured(b, io) ? 2 : 1); ++c)
            currents.sample(
              bt, W(b, io), io, b, c, ket,
              [&](std::size_t g0, std::size_t g1, value_type x) {
                value_type &mean = Global.bond_mean(f, g0 + L0 * g1);
                const value_type d = x - mean;
                mean += d / k;
                Global.bond_m2(f, g0 + L0 * g1) += d * (x - mean);
              }
            );
        }
    }
#pragma omp barrier
    store_bond_map(randoms_);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_bond_map(const int samples_)
  requires Complex<T>
{
#pragma omp master
  {
    H5::H5File file(name, H5F_ACC_RDWR);
    const value_type k = samples_;
    // standard error of the mean from the unbiased sample variance
    Eigen::Array<value_type, -1, -1> err = Global.bond_m2;
    if (samples_ > 1)
      err = (err / (k * (k - 1))).sqrt();
    else
      err.setZero();
    Eigen::Array<value_type, -1, -1> n_samples(1, 1);
    n_samples(0, 0) = k;
    std::string ng = bond_grp + "BondMean";
    write_hdf5(Global.bond_mean, &file, ng);
    ng = bond_grp + "BondStdErr";
    write_hdf5(err, &file, ng);
    ng = bond_grp + "Samples";
    write_hdf5(n_samples, &file, ng);
    file.close();
  }
#pragma omp barrier
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
