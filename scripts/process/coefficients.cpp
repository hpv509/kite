#include "coefficients.hpp"

arr<type, -1, 1>
Coefficients::gauss_first(type n, const arr<type, -1, 1> &mu, const type sigma)
{
  const arr<type, -1, 1> tmp = 1.0 - mu.square();
  const arr<type, -1, 1> numerator =
    n * mu * (sigma * sigma) * ((n * n * sigma * sigma) / tmp - 3.0);
  const arr<type, -1, 1> denominator = tmp.pow(-1.5);
  return numerator * denominator;
}

arr<type, -1, 1>
Coefficients::gauss_second(type n, const arr<type, -1, 1> &mu, const type sigma)
{
  const arr<type, -1, 1> term_1 = 7.0 * mu.square() - 4.0;
  const arr<type, -1, 1> tmp = 1.0 - mu.square();
  const arr<type, -1, 1> term_2 = 3.0 - 6.0 * n * n * sigma * sigma / tmp +
                                  std::pow(n * sigma, 4) / tmp.square();
  const arr<type, -1, 1> denominator = 1.0 / (24.0 * tmp);
  return sigma * sigma * term_1 * term_2 * denominator;
}

mtx<type, -1, -1> Coefficients::build_gaussian(
  const mtx<type, -1, 1> &energy_phys,
  const type sigma_phys,
  const type energy_scale,
  const unsigned Moments_G
)
{
  const type sigma = sigma_phys / energy_scale;
  const arr<type, -1, 1> E = energy_phys.array() / energy_scale;
  const int num_E = static_cast<int>(E.size());
  mtx<type, -1, -1> coefs(Moments_G, num_E);

  coefs.row(0) = 1.0 - gauss_second(0.0, E, sigma);
  for (int n = 1; n < Moments_G; ++n) {
    const type n_f = static_cast<type>(n);
    const arr<type, -1, 1> gaussian =
      (-0.5 * n_f * n_f * sigma * sigma / (1.0 - E.square())).exp();
    const arr<type, -1, 1> acos_E = E.acos();
    const arr<type, -1, 1> cossine = (n_f * acos_E).cos();
    const arr<type, -1, 1> sine = (n_f * acos_E).sin();
    const arr<type, -1, 1> g1 = gauss_first(n_f, E, sigma);
    const arr<type, -1, 1> g2 = gauss_second(n_f, E, sigma);
    coefs.row(n) = 2.0 * gaussian * (cossine * (1.0 - g2) - 0.5 * sine * g1);
  }
  const arr<type, -1, 1> prefactor =
    1.0 / (std::numbers::pi * (1.0 - E.square()).sqrt());
  coefs = coefs.array().rowwise() * prefactor.transpose();
  return coefs / energy_scale;
}

mtx<cplx, -1, -1> Coefficients::build_dgreen(
  const mtx<cplx, -1, 1> &z_phys,
  const type energy_scale,
  const unsigned Moments_D
)
{
  const arr<cplx, -1, 1> z = z_phys.array() / energy_scale;
  const int num_z = static_cast<int>(z.size());
  mtx<cplx, -1, -1> coefs(Moments_D, num_z);

  const arr<cplx, -1, 1> sq = 1.0 - z.square();
  const arr<cplx, -1, 1> sqr = sq.sqrt();
  const arr<cplx, -1, 1> diff = z - cst::I * sqr;
  arr<cplx, -1, 1> current_power = arr<cplx, -1, 1>::Ones(num_z);

  coefs.row(0) = -cst::I * z * current_power / (sq * sqr);
  current_power *= diff;
  for (int n = 1; n < Moments_D; ++n) {
    coefs.row(n) = 2.0 * (n - cst::I * z / sqr) * current_power / sq;
    current_power *= diff;
  }
  coefs /= (energy_scale * energy_scale);
  return coefs;
}
