#ifndef PAIRING_H_
#define PAIRING_H_
#ifndef PAIRING
#define PAIRING 0
#endif
#include "ComplexTraits.hpp"
#include "BondTable.hpp"

namespace pairing {

enum : unsigned { none = 0, onsite = 1, nearest = 2, both = 3 };

static inline constexpr unsigned channel = PAIRING;

static inline constexpr bool is_bdg = (channel != none);
static inline constexpr bool is_s_wave = (channel & onsite) != 0;
static inline constexpr bool is_p_wave = (channel & nearest) != 0;

static_assert(channel <= both, "PAIRING must be 0, 1, 2 or 3");
}

template <unsigned D> struct LatticeStructure;

template <typename T, unsigned D> struct PairingStructure {
  using value_type = typename extract_scalar<T>::type;

  LatticeStructure<D> &r;
  BondTable<D> bonds;
  value_type energy_scale = 0;
  value_type U = 0, V = 0;
  Eigen::Array<T, -1, 1> SDelta0;
  Eigen::Array<T, -1, -1> Delta0;
  bool has_onsite = false;
  bool has_bonds = false;

  PairingStructure(char *, LatticeStructure<D> &);

  void allocate_s(Eigen::Array<T, -1, 1> &) const;
  void
  broadcast_s(const Eigen::Array<T, -1, 1> &, Eigen::Array<T, -1, 1> &) const;
  void
  orbital_sum(const Eigen::Array<T, -1, 1> &, Eigen::Array<T, -1, 1> &) const;

  void allocate(Eigen::Array<T, -1, -1> &, bool per_site = false) const;
  void
  broadcast(const Eigen::Array<T, -1, -1> &, Eigen::Array<T, -1, -1> &) const;
  void
  symmetrize(const Eigen::Array<T, -1, -1> &, Eigen::Array<T, -1, -1> &) const;
  void symmetrize_bonds(Eigen::Array<T, -1, -1> &) const;
  void average_bonds(Eigen::Array<T, -1, -1> &) const;
  void print() const { bonds.print("Pairing table"); }

  template <unsigned MULT>
  void diagonal(
    T *phi0,
    const T *phiM1,
    const std::size_t j0,
    const Eigen::Array<value_type, -1, 1> &onsite
  ) const
    requires(D == 2)
  {
    constexpr value_type order = MULT + 1;
    const std::size_t h = r.offset;
    for (std::size_t j = j0; j < j0 + TILE * r.Ld[0]; j += r.Ld[0])
      for (std::size_t i = j; i < j + TILE; ++i) {
        const T ht = order * onsite(i);
        phi0[i] += phiM1[i] * ht;
        phi0[i + h] -= phiM1[i + h] * ht;
      }
  }

  template <unsigned MULT>
  void s_wave(
    T *phi0,
    const T *phiM1,
    const std::size_t j0,
    const Eigen::Array<T, -1, 1> &delta
  ) const
    requires(D == 2)
  {
    constexpr value_type order = MULT + 1;
    const std::size_t h = r.offset;
    for (std::size_t j = j0; j < j0 + TILE * r.Ld[0]; j += r.Ld[0])
      for (std::size_t i = j; i < j + TILE; ++i) {
        const T sd = order * delta(i);
        phi0[i] += phiM1[i + h] * sd;
        phi0[i + h] += phiM1[i] * cj(sd);
      }
  }

  template <unsigned MULT>
  void p_wave(
    T *phi0,
    const T *phiM1,
    const std::size_t j0,
    const std::size_t io,
    T *const (*fb)[3],
    const Eigen::Array<T, -1, -1> &F
  ) const
    requires(D == 2 && Complex<T>)
  {
    constexpr value_type order = MULT + 1;
    const std::size_t h = r.offset;
    const bool per_family = F.cols() == r.Orb;
    const std::size_t gs = per_family ? 0 : F.rows();
    const std::size_t x0 = j0 % r.Ld[0];
    const std::size_t x1 = (j0 / r.Ld[0]) % r.Ld[1];
    for (unsigned b = 0, B = bonds.NBonds(io); b < B; ++b) {
      const std::ptrdiff_t s = bonds.dist_tile(b, io);
      const T *const fx = fb[0][bonds.s(0, b, io) + 1] + x0;
      const T *const fy = fb[1][bonds.s(1, b, io) + 1] + x1;
      const T *const g = &F(b, per_family ? io : 0);
      const bool bulk_x =
        std::all_of(fx, fx + TILE, [](const T &z) { return z == T(1); });
      for (std::size_t y = 0; y < TILE; ++y) {
        const std::size_t row = j0 + y * r.Ld[0];
        const T *const in = phiM1 + std::ptrdiff_t(row) + s;
        if (bulk_x && fy[y] == T(1))
          for (std::size_t x = 0, i = row; x < TILE; ++x, ++i) {
            const T f = order * g[i * gs];
            phi0[i] += f * in[x + h];
            phi0[i + h] += std::conj(f) * in[x];
          }
        else
          for (std::size_t x = 0, i = row; x < TILE; ++x, ++i) {
            const T bp = fy[y] * fx[x];
            const T f = order * g[i * gs];
            phi0[i] += f * bp * in[x + h];
            phi0[i + h] += std::conj(f) * bp * in[x];
          }
      }
    }
  }

private:
  template <typename F> void for_rows(unsigned io, F &&f) const;

  static T cj(const T &z)
  {
    if constexpr (Complex<T>)
      return std::conj(z);
    else
      return z;
  }
};

#endif
