#include "Generic.hpp"
#include "myHDF5.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "BondTable.hpp"
#include <set>

namespace {
[[noreturn]] void bond_table_fail(const std::string &msg)
{
#pragma omp master
  std::cerr << "BondTable: " << msg << "\n";
#pragma omp barrier
  exit(1);
}

template <typename T>
bool try_read(H5::H5File &file, const std::string &path, T *const dst)
{
  std::string p = path; // get_hdf5 takes a non-const reference
  try {
    get_hdf5<T>(dst, &file, p);
    return true;
  } catch (H5::Exception &) {
    return false;
  }
}
} // namespace

template <unsigned D>
BondTable<D>::BondTable(
  const char *name,
  const LatticeStructure<D> &rr,
  const std::string &group,
  const int default_sector
) :
  r(rr),
  orb(rr.Orb),
  NBonds(Eigen::Array<unsigned, -1, 1>::Zero(rr.Orb)),
  first(Eigen::Array<unsigned, -1, 1>::Zero(rr.Orb))
{
  using IArr = Eigen::Array<int, -1, -1>;
  IArr dist;
#pragma omp critical
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    H5::Exception::dontPrint();
    if (!try_read(file, group + "NBonds", NBonds.data()))
      try_read(file, group + "NPairings", NBonds.data());
    n_bonds = NBonds.sum();
    if (n_bonds) {
      max_bonds = NBonds.maxCoeff();
      dist = IArr::Constant(max_bonds, orb, -1);
      sector = IArr::Constant(max_bonds, orb, default_sector);
      try_read(file, group + "d", dist.data());
      try_read(file, group + "Sector", sector.data());
    }
    file.close();
  }
  if (n_bonds == 0)
    return;
  for (unsigned io = 1; io < orb; ++io)
    first(io) = first(io - 1) + NBonds(io - 1);

  target_orb = IArr::Zero(max_bonds, orb);
  shift = IArr::Zero(D * max_bonds, orb);
  parity_axis = IArr::Constant(max_bonds, orb, -1);
  rev_orb = IArr::Constant(max_bonds, orb, -1);
  rev_bond = IArr::Constant(max_bonds, orb, -1);
  dist_tile = Eigen::Array<std::ptrdiff_t, -1, -1>::Zero(max_bonds, orb);
  target_offset = Eigen::Array<std::size_t, -1, -1>::Zero(max_bonds, orb);

  unsigned lB3[D + 1], Ld[D + 1];
  std::copy_n(r.lB3, D + 1, lB3);
  std::copy_n(r.Ld, D + 1, Ld);
  Coordinates<std::ptrdiff_t, D + 1> b3(lB3);
  Coordinates<std::size_t, D + 1> xd(Ld);
  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0; b < NBonds(io); ++b) {
      const int sec = sector(b, io);
      if (dist(b, io) < 0)
        bond_table_fail(group + "d is missing or incomplete");
      b3.set_coord(dist(b, io));
      const std::ptrdiff_t jo = b3.coord[D];
      if (jo >= std::ptrdiff_t(orb) || (sec != 0 && sec != 1))
        bond_table_fail(group + ": bad target orbital or sector");
      target_orb(b, io) = jo;
      target_offset(b, io) = sec ? r.offset : 0;

      std::ptrdiff_t st = (jo - std::ptrdiff_t(io)) * xd.basis[D];
      for (unsigned k = 0; k < D; ++k) {
        shift(k * max_bonds + b, io) = b3.coord[k] - 1;
        st += (b3.coord[k] - 1) * std::ptrdiff_t(xd.basis[k]);
      }
      dist_tile(b, io) = st;
      if (st == 0 && sec == 0)
        bond_table_fail(
          group + ": a family from a site to itself is not a bond"
        );
      if (sec == 0 && jo == std::ptrdiff_t(io))
        for (unsigned k = 0; k < D; ++k)
          if (s(k, b, io)) {
            if (r.Lt[k] % 2)
              bond_table_fail(
                group +
                ": a same-orbital family needs an even length along its axis"
              );
            parity_axis(b, io) = k;
            break;
          }
    }
  has_reverse = true;
  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0; b < NBonds(io); ++b) {
      const unsigned jo = target_orb(b, io);
      for (unsigned jb = 0; jb < NBonds(jo); ++jb)
        if (target_orb(jb, jo) == int(io) && sector(jb, jo) == sector(b, io) &&
            dist_tile(jb, jo) == -dist_tile(b, io)) {
          if (rev_bond(b, io) >= 0)
            bond_table_fail(group + ": duplicated bond family");
          rev_orb(b, io) = jo;
          rev_bond(b, io) = jb;
        }
      has_reverse &= rev_bond(b, io) >= 0;
    }
  std::vector<std::set<std::pair<int, int>>> used;
  for (unsigned io = 0; io < orb; ++io)
    for (unsigned b = 0; b < NBonds(io); ++b) {
      const std::pair<int, int> src{int(io), 0},
        dst{target_orb(b, io), sector(b, io)};
      std::size_t t = 0;
      while (t < used.size() && (used[t].count(src) || used[t].count(dst)))
        ++t;
      if (t == used.size()) {
        used.emplace_back();
        tilings.emplace_back();
      }
      used[t].insert({src, dst});
      tilings[t].push_back({b, io});
    }
}

template <unsigned D> void BondTable<D>::print(const std::string &title) const
{
#pragma omp master
  {
    std::cout << title << ": " << n_bonds << " families in " << tilings.size()
              << " tilings\n";
    for (std::size_t t = 0; t < tilings.size(); ++t) {
      std::cout << "  tiling " << t << " (" << n_colours(tilings[t])
                << " colour(s)):";
      for (const auto &[b, io] : tilings[t]) {
        std::cout << "  o" << io << "->o" << target_orb(b, io)
                  << (sector(b, io) ? "h" : "") << "(";
        for (unsigned k = 0; k < D; ++k)
          std::cout << (k ? "," : "") << s(k, b, io);
        std::cout << ")";
      }
      std::cout << "\n";
    }
  }
#pragma omp barrier
}

template struct BondTable<1u>;
template struct BondTable<2u>;
template struct BondTable<3u>;
