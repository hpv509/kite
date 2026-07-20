#include "HamiltonianBdG.hpp"

template <Scalar T, unsigned D>
HamiltonianBdG<T, D>::HamiltonianBdG(const std::size_t N_) :
  hartree(N_), s_delta(N_)
{}

template <Scalar T, unsigned D>
void HamiltonianBdG<T, D>::init_fields(const real delta_, const real ht_)
{
  hartree.setOnes();
  hartree *= ht_;
  s_delta.setOnes();
  s_delta *= delta_;
}

#define instantiate(type, dim) template struct HamiltonianBdG<type, dim>;
#include "instantiate.hpp"
