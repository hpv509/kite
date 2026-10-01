#ifndef BOND_TABLE_H_
#define BOND_TABLE_H_

template <unsigned D> struct LatticeStructure;

template <unsigned D> struct BondTable {
  // family is one bond type: from orbital io, go to jo in the cell shifted
  // by (dx, dy, ...) in electron (0) or hole (1) sector
  using family = std::array<unsigned, 2>;

  const LatticeStructure<D> &r;
  unsigned orb = 0;
  unsigned max_bonds = 0;
  unsigned n_bonds = 0;
  bool has_reverse = false;

  // number of families starting at io
  Eigen::Array<unsigned, -1, 1> NBonds;
  // first(io) offset for linear indexing of families
  Eigen::Array<unsigned, -1, 1> first;
  // target_orb(b, io) = jo
  Eigen::Array<int, -1, -1> target_orb;
  // sector(b, io) 0: electron, 1: hole
  Eigen::Array<int, -1, -1> sector;
  // shift(k * max_bonds + b, io) family cell displacement along direction k {-1, 0, 1}
  Eigen::Array<int, -1, -1> shift;
  Eigen::Array<int, -1, -1> parity_axis;
  Eigen::Array<int, -1, -1> rev_orb;
  Eigen::Array<int, -1, -1> rev_bond;
  Eigen::Array<std::ptrdiff_t, -1, -1> dist_tile;
  // target_offset(b, io) start of the target block: 0 or r.offset (hole)
  Eigen::Array<std::size_t, -1, -1> target_offset;
  std::vector<std::vector<family>> tilings;

  BondTable(
    const char *name,
    const LatticeStructure<D> &r,
    const std::string &group,
    const int default_sector
  );

  unsigned family_index(const unsigned b, const unsigned io) const
  {
    return first(io) + b;
  }
  int s(const unsigned k, const unsigned b, const unsigned io) const
  {
    return shift(k * max_bonds + b, io);
  }
  bool coloured(const unsigned b, const unsigned io) const
  {
    return parity_axis(b, io) >= 0;
  }
  unsigned n_colours(const std::vector<family> &t) const
  {
    for (const auto &[b, io] : t)
      if (coloured(b, io))
        return 2;
    return 1;
  }

  template <int S>
  void window(
    const unsigned io,
    const unsigned b,
    std::size_t (&beg)[D],
    std::size_t (&end)[D]
  ) const
  {
    for (unsigned k = 0; k < D; ++k) {
      beg[k] = NGHOSTS;
      end[k] = r.Ld[k] - NGHOSTS;
      const int dr = S * s(k, b, io);
      if (dr < 0 && !r.boundary[k][0])
        beg[k] += 1;
      else if (dr > 0 && !r.boundary[k][1])
        end[k] -= 1;
    }
  }

  void print(const std::string &title) const;
};

#endif
