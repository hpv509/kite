#include "HamiltonianBdG.hpp"

template <Scalar T, unsigned D>
HamiltonianBdG<T, D>::HamiltonianBdG(const std::size_t N_) :
  hartree(N_), s_delta(N_)
{
  hartree.setZero();
  s_delta.setZero();
}

#define instantiate(type, dim) template struct HamiltonianBdG<type, dim>;
#include "instantiate.hpp"
