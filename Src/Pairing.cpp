#include "Generic.hpp"
#include "myHDF5.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "Pairing.hpp"

template <typename T, unsigned D>
PairingStructure<T, D>::PairingStructure(char *name, LatticeStructure<D> &rr) :
  r(rr)
{
  orb = r.Orb;
  max_bonds = 0;
  n_bonds = 0;
  U = 0;
  V = 0;
  NPairings = Eigen::Array<unsigned, -1, 1>::Zero(orb);
  SDelta0 = Eigen::Array<T, -1, 1>::Zero(orb);

  const std::string base_dir = "/Pairing/";
  std::string tmp;
#pragma omp critical
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    H5::Exception::dontPrint();
    tmp = "/EnergyScale";
    get_hdf5<value_type>(&energy_scale, &file, tmp);
    try {
      tmp = base_dir + "U";
      get_hdf5<value_type>(&U, &file, tmp);
      tmp = base_dir + "SDelta0";
      get_hdf5<T>(SDelta0.data(), &file, tmp);
      has_onsite = true;
    } catch (H5::Exception &e) {
      U = 0;
      SDelta0.setZero();
      has_onsite = false;
    }
    try {
      tmp = base_dir + "NPairings";
      get_hdf5<unsigned>(NPairings.data(), &file, tmp);
      has_bonds = (NPairings.sum() > 0);
    } catch (H5::Exception &e) {
      NPairings.setZero();
      has_bonds = false;
    }
    file.close();
  }

  if constexpr (pairing::is_s_wave)
    if (!has_onsite) {
#pragma omp master
      std::cerr
        << "PairingStructure: PAIRING has the onsite bit set but the "
           "configuration file carries no /Pairing/U or /Pairing/SDelta0.\n";
#pragma omp barrier
      exit(1);
    }
  if constexpr (pairing::is_p_wave)
    if (!has_bonds) {
#pragma omp master
      std::cerr << "PairingStructure: PAIRING has the nearest bit set but the "
                   "configuration file carries no bond table.\n";
#pragma omp barrier
      exit(1);
    }

  if (!has_bonds)
    return;

  max_bonds = NPairings.maxCoeff();
  n_bonds = NPairings.sum();

  dist.resize(max_bonds, orb);
  dist = Eigen::Array<int, -1, -1>::Zero(max_bonds, orb);
  rev_orb = Eigen::Array<int, -1, -1>::Constant(max_bonds, orb, -1);
  rev_bond = Eigen::Array<int, -1, -1>::Constant(max_bonds, orb, -1);
  Delta0 = Eigen::Array<T, -1, -1>::Zero(max_bonds, orb);
  SDelta0 = Eigen::Array<T, -1, 1>::Zero(orb);

#pragma omp critical
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    H5::Exception::dontPrint();
    tmp = "/EnergyScale";
    get_hdf5<value_type>(&energy_scale, &file, tmp);
    try {
      tmp = base_dir + "d";
      get_hdf5<int>(dist.data(), &file, tmp);
      tmp = base_dir + "ReverseOrbital";
      get_hdf5<int>(rev_orb.data(), &file, tmp);
      tmp = base_dir + "ReverseBond";
      get_hdf5<int>(rev_bond.data(), &file, tmp);
      tmp = base_dir + "Delta0";
      get_hdf5<T>(Delta0.data(), &file, tmp);
      tmp = base_dir + "V";
      get_hdf5<value_type>(&V, &file, tmp);
    } catch (H5::Exception &e) {
      std::cerr << "PairingStructure: incomplete bond table, cannot read "
                << tmp << "\n";
      exit(1);
    }
    try {
      tmp = base_dir + "U";
      get_hdf5<value_type>(&U, &file, tmp);
      tmp = base_dir + "SDelta0";
      get_hdf5<T>(SDelta0.data(), &file, tmp);
    } catch (H5::Exception &e) {
      U = 0;
      SDelta0.setZero();
    }
    file.close();
  }
  dist_tile.resize(max_bonds, orb);
  dist_tile.setZero();
  target_orb.resize(max_bonds, orb);
  target_orb.setZero();
  shift.resize(D * max_bonds, orb);
  shift.setZero();

  Coordinates<std::ptrdiff_t, D + 1> b3(r.lB3);
  Coordinates<std::size_t, D + 1> xd(r.Ld);

  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0, B = NPairings(io); b < B; ++b) {
      b3.set_coord(dist(b, io));
      const std::ptrdiff_t jo = b3.coord[D];
      target_orb(b, io) = jo;

      std::ptrdiff_t s = 0;
      for (unsigned k = 0; k < D; ++k) {
        const std::ptrdiff_t dr = b3.coord[k] - 1;
        shift(k * max_bonds + b, io) = dr;
        s += dr * xd.basis[k];
      }
      s += (jo - static_cast<std::ptrdiff_t>(io)) * xd.basis[D];
      dist_tile(b, io) = s;
    }
  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0, B = NPairings(io); b < B; ++b) {
      const int jo = rev_orb(b, io);
      const int jb = rev_bond(b, io);

      bool bad = (jo < 0) || (jb < 0) || (static_cast<unsigned>(jo) >= orb) ||
                 (static_cast<unsigned>(jb) >= NPairings(jo));
      if (!bad)
        bad = (rev_orb(jb, jo) != static_cast<int>(io)) ||
              (rev_bond(jb, jo) != static_cast<int>(b)) ||
              (dist_tile(jb, jo) != -dist_tile(b, io));

      if (bad) {
#pragma omp master
        std::cerr << "PairingStructure: inconsistent reverse bond for orbital "
                  << io << ", bond " << b << "\n";
#pragma omp barrier
        exit(1);
      }
    }
}

template <typename T, unsigned D>
void PairingStructure<T, D>::broadcast_s(
  const Eigen::Array<T, -1, 1> &orb_values_,
  Eigen::Array<T, -1, 1> &field_
) const
{
  field_.setZero();
  Coordinates<std::size_t, D + 1> local(r.Ld);

  for (unsigned io = 0; io < orb; ++io) {
    const T d0 = orb_values_(io);
    if constexpr (D == 2) {
      for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
        const std::size_t j0 = local.set({NGHOSTS, i1, io}).index;
        for (std::size_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
          field_(i) = d0;
      }
    } else if constexpr (D == 3) {
      for (std::size_t i2 = NGHOSTS; i2 < r.Ld[2] - NGHOSTS; ++i2)
        for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
          const std::size_t j0 = local.set({NGHOSTS, i1, i2, io}).index;
          for (std::size_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
            field_(i) = d0;
        }
    }
  }
#pragma omp barrier
}

template <typename T, unsigned D>
void PairingStructure<T, D>::allocate_s(Eigen::Array<T, -1, 1> &field) const
{
  field.resize(r.Sized);
  if (!has_onsite) {
    field.setZero();
    return;
  }
  broadcast_s(SDelta0, field);
}

template <typename T, unsigned D>
void PairingStructure<T, D>::orbital_sum(
  const Eigen::Array<T, -1, 1> &field_,
  Eigen::Array<T, -1, 1> &orb_values_
) const
{
  Coordinates<std::size_t, D + 1> local(r.Ld);
  orb_values_.resize(orb);
  orb_values_.setZero();

  for (unsigned io = 0; io < orb; ++io) {
    T acc = 0;
    if constexpr (D == 2) {
      for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
        const std::size_t j0 = local.set({NGHOSTS, i1, io}).index;
        for (std::size_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
          acc += field_(i);
      }
    } else if constexpr (D == 3) {
      for (std::size_t i2 = NGHOSTS; i2 < r.Ld[2] - NGHOSTS; ++i2)
        for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
          const std::size_t j0 = local.set({NGHOSTS, i1, i2, io}).index;
          for (std::size_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
            acc += field_(i);
        }
    }
    orb_values_(io) = acc;
  }
}

template <typename T, unsigned D>
void PairingStructure<T, D>::symmetrize(
  const Eigen::Array<T, -1, -1> &raw,
  Eigen::Array<T, -1, -1> &result
) const
{
  if (!has_bonds)
    return;

  constexpr value_type half = 0.5;
  Coordinates<std::size_t, D + 1> local(r.Ld);

  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0; b < NPairings(io); ++b) {
      const unsigned jb = static_cast<unsigned>(rev_bond(b, io));
      const std::ptrdiff_t s = dist_tile(b, io);

      if constexpr (D == 2) {
        for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
          local.set({NGHOSTS, i1, io});
          const std::size_t j0 = local.index;
          for (std::ptrdiff_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
            result(b, i) =
              half * (raw(b, i) + raw(jb, static_cast<std::size_t>(i + s)));
        }
      } else if constexpr (D == 3) {
        for (std::size_t i2 = NGHOSTS; i2 < r.Ld[2] - NGHOSTS; ++i2)
          for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
            local.set({NGHOSTS, i1, i2, io});
            const std::size_t j0 = local.index;
            for (std::ptrdiff_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
              result(b, i) =
                half * (raw(b, i) + raw(jb, static_cast<std::size_t>(i + s)));
          }
      }
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
  constexpr value_type half = 0.5;
  const Eigen::Array<T, -1, -1> map = bond_values;
  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0, B = NPairings(io); b < B; ++b)
      bond_values(b, io) =
        half * (map(b, io) + map(rev_bond(b, io), rev_orb(b, io)));
}

template <typename T, unsigned D>
void PairingStructure<T, D>::average_bonds(Eigen::Array<T, -1, -1> &bond_values
) const
{
  if (!has_bonds)
    return;
  for (unsigned io = 0; io < orb; ++io) {
    T mean = 0;
    const unsigned B = NPairings(io);
    for (unsigned b = 0; b < B; ++b)
      mean += bond_values(b, io);
    mean /= value_type(B);
    for (unsigned b = 0; b < B; ++b)
      bond_values(b, io) = mean;
  }
}

template <typename T, unsigned D>
void PairingStructure<T, D>::broadcast(
  const Eigen::Array<T, -1, -1> &bond_values_,
  Eigen::Array<T, -1, -1> &field_
) const
{
  if (!has_bonds)
    return;
  field_.setZero();
  Coordinates<std::size_t, D + 1> local(r.Ld);

  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0; b < NPairings(io); ++b) {
      const T d0 = bond_values_(b, io);
      if constexpr (D == 2) {
        for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
          const std::size_t j0 = local.set({NGHOSTS, i1, io}).index;
          for (std::size_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
            field_(b, i) = d0;
        }
      } else if constexpr (D == 3) {
        for (std::size_t i2 = NGHOSTS; i2 < r.Ld[2] - NGHOSTS; ++i2)
          for (std::size_t i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
            const std::size_t j0 = local.set({NGHOSTS, i1, i2, io}).index;
            for (std::size_t i = j0, j1 = j0 + r.ld[0]; i < j1; ++i)
              field_(b, i) = d0;
          }
      }
    }
#pragma omp barrier
}

template <typename T, unsigned D>
void PairingStructure<T, D>::allocate(Eigen::Array<T, -1, -1> &field) const
{
  if (!has_bonds)
    return;
  field.resize(max_bonds, r.Sized);
  broadcast(Delta0, field);
}

template <typename T, unsigned D> void PairingStructure<T, D>::print() const
{
#pragma omp master
  {
    if (!has_bonds) {
      std::cout << "Pairing table: absent.\n";
    } else {
      std::cout << "Pairing table: " << n_bonds << " directed bonds\n";
      for (unsigned io = 0; io < orb; ++io) {
        std::cout << "  orbital " << io << " (" << NPairings(io) << "):";
        for (unsigned b = 0; b < NPairings(io); ++b) {
          std::cout << "  [b=" << b << " -> o" << target_orb(b, io) << " dR=(";
          for (unsigned k = 0; k < D; ++k)
            std::cout << (k ? "," : "") << shift(k * max_bonds + b, io);
          std::cout << ")]";
        }
        std::cout << "\n";
      }
    }
  }
#pragma omp barrier
}

#define instantiate(type, dim) template struct PairingStructure<type, dim>;
#include "instantiate.hpp"
