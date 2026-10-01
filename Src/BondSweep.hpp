#ifndef BOND_SWEEP_H_
#define BOND_SWEEP_H_
#include <array>
#include "BondTable.hpp"

template <typename T, unsigned D> class BondSweep;

template <typename T> class BondSweep<T, 2> {
protected:
  using value_type = typename extract_scalar<T>::type;
  LatticeStructure<2u> &r;
  T *const (*fb)[3]; // Fact_Bnd of the KPM_Vector
  std::array<std::size_t, 2> origin;

public:
  template <typename Vec>
  BondSweep(LatticeStructure<2u> &r_, const Vec &v) : r(r_), fb(v.Fact_Bnd)
  {
    Coordinates<std::size_t, 3> loc(r.Ld), glob(r.Lt);
    loc.set({std::size_t(NGHOSTS), std::size_t(NGHOSTS), std::size_t(0)});
    r.convertCoordinates(glob, loc);
    origin = {std::size_t(glob.coord[0]), std::size_t(glob.coord[1])};
  }

protected:
  template <int S, typename F>
  void
  for_pairs(const BondTable<2u> &bt, unsigned io, unsigned b, int colour, F &&f)
    const
  {
    static_assert(S == 0 || S == 1);
    const int k = bt.parity_axis(b, io);
    const std::ptrdiff_t s = bt.dist_tile(b, io);
    const std::size_t off = bt.target_offset(b, io);
    const std::size_t base = io * r.Nd;
    const int sh[2] = {bt.s(0, b, io), bt.s(1, b, io)};
    std::size_t beg[2], end[2], step[2] = {1, 1};
    bt.template window<1>(io, b, beg, end);
    if constexpr (S == 0) {
      std::size_t tb[2], te[2];
      bt.template window<-1>(io, b, tb, te);
      for (unsigned d = 0; d < 2; ++d) {
        beg[d] = std::min(beg[d], std::size_t(std::ptrdiff_t(tb[d]) - sh[d]));
        end[d] = std::max(end[d], std::size_t(std::ptrdiff_t(te[d]) - sh[d]));
      }
    }
    auto twist = [&](unsigned d, std::size_t x) -> T {
      if (S == 1 || (x >= NGHOSTS && x < r.Ld[d] - NGHOSTS))
        return fb[d][sh[d] + 1][x];
      return std::conj(fb[d][1 - sh[d]][x + sh[d]]);
    };
    if (k >= 0) {
      const std::size_t g = origin[k] + beg[k] - NGHOSTS;
      beg[k] += (g + unsigned(colour)) & 1u;
      step[k] = 2;
    }
    for (std::size_t i1 = beg[1]; i1 < end[1]; i1 += step[1]) {
      const T ty = twist(1, i1);
      for (std::size_t i0 = beg[0]; i0 < end[0]; i0 += step[0]) {
        const std::size_t a = base + i1 * r.Ld[0] + i0;
        f(a, std::size_t(std::ptrdiff_t(a) + s) + off, ty * twist(0, i0),
          origin[0] + i0 - NGHOSTS, origin[1] + i1 - NGHOSTS);
      }
    }
  }
};

#endif
