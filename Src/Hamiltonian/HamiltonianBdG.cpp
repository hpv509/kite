#include "HamiltonianBdG.hpp"

template <Scalar T, unsigned D>
HamiltonianBdG<T, D>::HamiltonianBdG(const std::size_t N_,
				     const std::size_t max_hoppings_) :
  hartree(N_), s_delta(N_), nn_delta(max_hoppings_, N_)
{
  hartree.setZero();
  s_delta.setZero();
  nn_delta.setZero();
}

template <Scalar T, unsigned D>
void HamiltonianBdG<T, D>::init_fields(
  const T s_delta_,
  const real ht_,
  const real mu_
)
{
  hartree += ht_ - mu_;
  s_delta += s_delta_;
}

template <Scalar T, unsigned D>
void HamiltonianBdG<T, D>::init_fields(
  const T s_delta_,
  const T nn_delta_,
  const real ht_,
  const real mu_
)
{
  init_fields(s_delta_, ht_, mu_);
  nn_delta += nn_delta_;
}

#define instantiate(type, dim) template struct HamiltonianBdG<type, dim>;
#include "instantiate.hpp"
