/***********************************************************/
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/
#ifndef RANDOM_H_
#define RANDOM_H_
#include <random>
#include "ComplexTraits.hpp"

template <scalar T>
class KPMRandom {
  std::mt19937 rng;
  std::uniform_real_distribution<double> dist;
  std::normal_distribution<> gauss;
public:
  using real = typename extract_scalar<T>::type;
  explicit KPMRandom(const unsigned);
  void init_random(const unsigned);
  real get();
  real uniform(const real, const real);
  real gaussian(const real, const real);
  T init();
};

#endif
