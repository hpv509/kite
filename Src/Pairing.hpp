#ifndef PAIRING_H_
#define PAIRING_H_
#ifndef PAIRING
#define PAIRING 0
#endif
#include "ComplexTraits.hpp"

namespace pairing {

enum : unsigned { none = 0, onsite = 1, nearest = 2, both = 3 };

static inline constexpr unsigned channel = PAIRING;

static inline constexpr bool is_bdg = (channel != none);
static inline constexpr bool is_s_wave = (channel & onsite) != 0;
static inline constexpr bool is_p_wave = (channel & nearest) != 0;

static_assert(channel <= both, "PAIRING must be 0, 1, 2 or 3");
} // namespace pairing

template <unsigned D>
struct LatticeStructure;

template <typename T, unsigned D>
struct PairingStructure {
  using value_type = typename extract_scalar<T>::type;

  LatticeStructure<D> &r;

  unsigned orb;
  unsigned max_bonds;
  unsigned n_bonds;
  value_type energy_scale;

  value_type U;
  Eigen::Array<T, -1, 1> SDelta0;
  value_type V;
  Eigen::Array<T, -1, -1> Delta0;

  bool has_bonds = false;
  bool has_onsite = false;

  Eigen::Array<unsigned, -1, 1> NPairings;
  Eigen::Array<int, -1, -1> dist;
  Eigen::Array<int, -1, -1> target_orb;
  Eigen::Array<int, -1, -1> rev_orb;
  Eigen::Array<int, -1, -1> rev_bond;
  Eigen::Array<std::ptrdiff_t, -1, -1> dist_tile;
  Eigen::Array<int, -1, -1> shift; // (D * max_bonds, orb), each in {-1, 0, +1}

  PairingStructure(char *, LatticeStructure<D> &);

  void allocate(Eigen::Array<T, -1, -1> &) const;
  void
  broadcast(const Eigen::Array<T, -1, -1> &, Eigen::Array<T, -1, -1> &) const;
  void symmetrize_bonds(Eigen::Array<T, -1, -1> &) const;
  void
  symmetrize(const Eigen::Array<T, -1, -1> &, Eigen::Array<T, -1, -1> &) const;

  void allocate_s(Eigen::Array<T, -1, 1> &) const;
  void
  broadcast_s(const Eigen::Array<T, -1, 1> &, Eigen::Array<T, -1, 1> &) const;
  void
  orbital_sum(const Eigen::Array<T, -1, 1> &, Eigen::Array<T, -1, 1> &) const;
  void print() const;
};

#endif
