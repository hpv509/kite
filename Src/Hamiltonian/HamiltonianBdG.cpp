#include "HamiltonianBdG.hpp"

template <Scalar T, unsigned D>
HamiltonianBdG<T, D>::HamiltonianBdG(
  const std::size_t N_,
  const std::size_t max_hoppings_
) :
  hartree(N_), s_delta(N_), nn_delta(max_hoppings_, N_)
{
  hartree.setZero();
  s_delta.setZero();
  nn_delta.setZero();
}

template <Scalar T, unsigned D>
void HamiltonianBdG<T, D>::init_fields(
  const Eigen::Array<real, -1, 1> &s_delta_,
  const Eigen::Array<real, -1, 1> &ht_
)
{
  const unsigned norb = s_delta_.size();
  const unsigned block = s_delta.size() / norb;
  for (unsigned i = 0; i < norb; ++i) {
    hartree.segment(i * block, block) = ht_(i);
    s_delta.segment(i * block, block) = static_cast<T>(s_delta_(i));
  }
  update_onsite();
}

template <Scalar T, unsigned D>
void HamiltonianBdG<T, D>::init_fields(
  const T s_delta_,
  const T nn_delta_,
  const real ht_
)
{
  hartree = ht_;
  s_delta = s_delta_;
  nn_delta = nn_delta_;
  update_onsite();
}

#define instantiate(type, dim) template struct HamiltonianBdG<type, dim>;
#include "instantiate.hpp"
