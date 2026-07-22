#include "HamiltonianBdG.hpp"

template <Scalar T, unsigned D>
HamiltonianBdG<T, D>::HamiltonianBdG(const std::size_t N_) :
  hartree(N_), s_delta(N_)
{
  hartree.setZero();
  s_delta.setZero();
}

template <Scalar T, unsigned D>
void HamiltonianBdG<T, D>::init_fields(
  const real delta_,
  const real ht_,
  const real mu_
)
{
  hartree += ht_ - mu_;
  s_delta += delta_;
}

#define instantiate(type, dim) template struct HamiltonianBdG<type, dim>;
#include "instantiate.hpp"
