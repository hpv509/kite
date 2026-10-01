#include "Generic.hpp"
#include "myHDF5.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "Pairing.hpp"

template <typename T, unsigned D>
PairingStructure<T, D>::PairingStructure(char *name, LatticeStructure<D> &rr) :
  r(rr), bonds(name, rr, "/Pairing/", 1)
{
  has_bonds = bonds.n_bonds > 0;
  SDelta0 = Eigen::Array<T, -1, 1>::Zero(r.Orb);
  if (has_bonds)
    Delta0 = Eigen::Array<T, -1, -1>::Zero(bonds.max_bonds, r.Orb);
  bool bonds_complete = true;
#pragma omp critical
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    H5::Exception::dontPrint();
    std::string path = "/EnergyScale";
    get_hdf5<value_type>(&energy_scale, &file, path);
    try {
      path = "/Pairing/U";
      get_hdf5<value_type>(&U, &file, path);
      path = "/Pairing/SDelta0";
      get_hdf5<T>(SDelta0.data(), &file, path);
      has_onsite = true;
    } catch (H5::Exception &) {
      U = 0;
      SDelta0.setZero();
    }
    if (has_bonds) {
      try {
        path = "/Pairing/V";
        get_hdf5<value_type>(&V, &file, path);
        path = "/Pairing/Delta0";
        get_hdf5<T>(Delta0.data(), &file, path);
      } catch (H5::Exception &) {
        bonds_complete = false;
      }
    }
    file.close();
  }
}

template <typename T, unsigned D>
template <typename F>
void PairingStructure<T, D>::for_rows(unsigned io, F &&f) const
{
  Coordinates<std::size_t, D + 1> local(r.Ld);
  if constexpr (D == 1)
    f(local.set({std::size_t(NGHOSTS), std::size_t(io)}).index);
  else if constexpr (D == 2)
    for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1)
      f(local.set({NGHOSTS, i1, io}).index);
  else
    for (std::size_t i2 = NGHOSTS; i2 < r.Ld[2] - NGHOSTS; ++i2)
      for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1)
        f(local.set({NGHOSTS, i1, i2, io}).index);
}

template <typename T, unsigned D>
void PairingStructure<T, D>::broadcast_s(
  const Eigen::Array<T, -1, 1> &orb_values_,
  Eigen::Array<T, -1, 1> &field_
) const
{
  field_.setZero();
  for (unsigned io = 0; io < r.Orb; ++io)
    for_rows(io, [&](std::size_t j0) {
      field_.segment(j0, r.ld[0]).setConstant(orb_values_(io));
    });
#pragma omp barrier
}

template <typename T, unsigned D>
void PairingStructure<T, D>::allocate_s(Eigen::Array<T, -1, 1> &field) const
{
  field.setZero(r.Sized);
  if (has_onsite)
    broadcast_s(SDelta0, field);
}

template <typename T, unsigned D>
void PairingStructure<T, D>::orbital_sum(
  const Eigen::Array<T, -1, 1> &field_,
  Eigen::Array<T, -1, 1> &orb_values_
) const
{
  orb_values_.setZero(r.Orb);
  for (unsigned io = 0; io < r.Orb; ++io)
    for_rows(io, [&](std::size_t j0) {
      orb_values_(io) += field_.segment(j0, r.ld[0]).sum();
    });
}

template <typename T, unsigned D>
void PairingStructure<T, D>::broadcast(
  const Eigen::Array<T, -1, -1> &bond_values_,
  Eigen::Array<T, -1, -1> &field_
) const
{
  if (!has_bonds)
    return;
  if (field_.cols() == Eigen::Index(r.Orb))
    field_ = bond_values_; // one gap per family
  else {
    field_.setZero();
    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned b = 0; b < bonds.NBonds(io); ++b)
        for_rows(io, [&](std::size_t j0) {
          field_.block(b, j0, 1, r.ld[0]).setConstant(bond_values_(b, io));
        });
  }
#pragma omp barrier
}

template <typename T, unsigned D>
void PairingStructure<T, D>::allocate(
  Eigen::Array<T, -1, -1> &field,
  const bool per_site
) const
{
  if (!has_bonds)
    return;
  field.resize(bonds.max_bonds, per_site ? r.Sized : r.Orb);
  broadcast(Delta0, field);
}

template <typename T, unsigned D>
void PairingStructure<T, D>::symmetrize(
  const Eigen::Array<T, -1, -1> &raw,
  Eigen::Array<T, -1, -1> &result
) const
{
  if (!has_bonds)
    return;
  for (unsigned io = 0; io < r.Orb; ++io)
    for (unsigned b = 0; b < bonds.NBonds(io); ++b) {
      const unsigned jb = bonds.rev_bond(b, io);
      const std::ptrdiff_t s = bonds.dist_tile(b, io);
      for_rows(io, [&](std::size_t j0) {
        for (std::size_t i = j0; i < j0 + r.ld[0]; ++i)
          result(b, i) =
            value_type(0.5) *
            (raw(b, i) + raw(jb, std::size_t(std::ptrdiff_t(i) + s)));
      });
    }
#pragma omp barrier
}

template <typename T, unsigned D>
void PairingStructure<T, D>::symmetrize_bonds(
  Eigen::Array<T, -1, -1> &bond_values
) const
{
  if (!has_bonds)
    return;
  const Eigen::Array<T, -1, -1> map = bond_values;
  for (unsigned io = 0; io < r.Orb; ++io)
    for (unsigned b = 0; b < bonds.NBonds(io); ++b)
      bond_values(b, io) =
        value_type(0.5) *
        (map(b, io) + map(bonds.rev_bond(b, io), bonds.rev_orb(b, io)));
}

template <typename T, unsigned D>
void PairingStructure<T, D>::average_bonds(Eigen::Array<T, -1, -1> &bond_values
) const
{
  if (!has_bonds)
    return;
  for (unsigned io = 0; io < r.Orb; ++io) {
    const unsigned B = bonds.NBonds(io);
    bond_values.col(io).head(B).setConstant(
      bond_values.col(io).head(B).sum() / value_type(B)
    );
  }
}

#define instantiate(type, dim) template struct PairingStructure<type, dim>;
#include "instantiate.hpp"
