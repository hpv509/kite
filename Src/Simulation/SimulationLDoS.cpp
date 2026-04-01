#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "Global.hpp"
#include "Random.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "queue.hpp"
#include "Simulation.hpp"
#include "Hamiltonian.hpp"
#include "KPM_VectorBasis.hpp"
#include "KPM_Vector.hpp"

template <typename T>
T gauss_first(const unsigned n_, const T mu_, const T sigma_)
{
  const T numerator = n_ * mu_ * sigma_ * sigma_ *
                      (n_ * n_ * sigma_ * sigma_ / (1 - mu_ * mu_) - 3);
  const T denominator = std::pow(1 - mu_ * mu_, -1.5);
  return numerator * denominator;
}

template <typename T>
T gauss_second(const unsigned n_, const T mu_, const T sigma_)
{
  const T term_1 = 7 * mu_ * mu_ - 4;
  const T tmp = 1 - mu_ * mu_;
  const T term_2 = 3 - 6 * n_ * n_ * sigma_ * sigma_ / tmp +
                   std::pow(n_ * sigma_, 4) / (tmp * tmp);
  const T denominator = 1 / (24 * tmp);
  return sigma_ * sigma_ * term_1 * term_2 * denominator;
}

template <typename T>
Eigen::Array<T, -1, 1> build_gaussian(const T energy_, const T width_)
{
  const unsigned number_polynomials = std::ceil(6.0 / width_);
  Eigen::Array<T, -1, 1> coefs(number_polynomials);
  coefs(0) = 1 - gauss_second(0, energy_, width_);

  for (int n = 1; n < number_polynomials; ++n) {
    const T gaussian =
      std::exp(-0.5 * n * n * width_ * width_ / (1 - energy_ * energy_));
    const T cossine = std::cos(n * std::acos(energy_));
    const T sine = std::sin(n * std::acos(energy_));
    coefs(n) = 2 * gaussian *
               (cossine * (1 - gauss_second(n, energy_, width_)) -
                0.5 * sine * gauss_first(n, energy_, width_));
  }
  const T prefactor = 1 / (M_PI * std::sqrt(1 - energy_ * energy_));
  coefs *= prefactor;
  return coefs;
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_ldos()
{
  debug_message("Entered Simulation::calc_ldos\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_ldos_map = false;
    try {
      int dummy_variable;
      get_hdf5<
        int>(&dummy_variable, file, (char *)"/Calculation/ldos_map/NumVectors");
      Global.calculate_ldos_map = true;
    } catch (H5::Exception &e) {
      debug_message("ldos: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_ldos_map = false;
#pragma omp critical
  local_calculate_ldos_map = Global.calculate_ldos_map;
#pragma omp barrier
  if (local_calculate_ldos_map) {
#pragma omp master
    std::cout << "Calculating LDoS.\n";
#pragma omp barrier
    int vectors;
    value_type energy;
    value_type sigma;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<int>(&vectors, file, (char *)"/Calculation/ldos_map/NumVectors");
      get_hdf5<
        value_type>(&energy, file, (char *)"/Calculation/ldos_map/Energy");
      get_hdf5<value_type>(&sigma, file, (char *)"/Calculation/ldos_map/Sigma");
      file->close();
      delete file;
    }
    ldos(vectors, energy, sigma);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::ldos(
  const int vectors_,
  const value_type energy_,
  const value_type sigma_
)
{
  debug_message("Entered ldos\n");
  if constexpr (is_tt<std::complex, T>::value) {
    value_type energy_scale;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
      file->close();
      delete file;
    }
#pragma omp barrier
    const value_type target = energy_ / energy_scale;
    const value_type sigma = sigma_ / energy_scale;
    const value_type fwhm = std::sqrt(2) * sigma;
    const value_type size = r.Sizet - r.SizetVacancies;
    const value_type factor = std::sqrt(8 * M_PI) * sigma;
    const Eigen::Array coefs = build_gaussian<value_type>(target, fwhm);

    KPM_Vector<T, D> phi(2, *this);
    Eigen::Array<T, -1, 1> ket(r.Sized);
    Eigen::Array<T, -1, 1> bra(r.Sized);
    Eigen::Array<T, -1, 1> results(r.Sized);

    std::array<KPM_Vector<T, D> *, 3> vectors;
    results.setZero();
    h.generate_disorder();
    for (int vec = 0; vec < vectors_; ++vec) {
      h.generate_twists();
      phi.initiate_phases();
      phi.set_index(0);
      phi.initiate_vector();
      phi.v.col(0) *= std::sqrt(size);
      bra = phi.v.col(0);
      ket.setZero();

      phi.Exchange_Boundaries();
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
        phi.cheb_iteration(n);
        ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      const value_type weight = 1.0 / (vec + 1);
      results += weight * (factor * (bra.conjugate() * ket).abs2() - results);
    }
    store_ldos(results);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::store_ldos(const Eigen::Array<T, -1, 1> &results_)
{
  debug_message("Entered store_ldos\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  Global.ldos_map.resize(r.Sizet, 1);
#pragma omp barrier
  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    if constexpr (D == 3) {
      for (unsigned i2 = NGHOSTS, I2 = r.Ld[2] - NGHOSTS; i2 < I2; ++i2)
        for (unsigned i1 = NGHOSTS, I1 = r.Ld[1] - NGHOSTS; i1 < I1; ++i1)
          for (unsigned i0 = NGHOSTS, I0 = r.Ld[0] - NGHOSTS; i0 < I0; ++i0) {
            local.set({i0, i1, i2, io});
            r.convertCoordinates(global, local);
            Global.ldos_map(global.index) = results_(local.index);
          }
    } else if constexpr (D == 2) {
      for (unsigned i1 = NGHOSTS, I1 = r.Ld[1] - NGHOSTS; i1 < I1; ++i1)
        for (unsigned i0 = NGHOSTS, I0 = r.Ld[0] - NGHOSTS; i0 < I0; ++i0) {
          local.set({i0, i1, io});
          r.convertCoordinates(global, local);
          Global.ldos_map(global.index) = results_(local.index);
        }
    }
  }
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    char buffer[200];
    std::sprintf(buffer, "/Calculation/ldos_map/Map");
    const std::string name(buffer);
    write_hdf5(Global.ldos_map, file, name);
    delete file;
  }
#pragma omp barrier
  debug_message("Left store_ldos\n");
}

#define INSTANTIATE_GAUSS(type)                                                \
  template type gauss_first<type>(const unsigned, const type, const type);     \
  template type gauss_second<type>(const unsigned, const type, const type);    \
  template Eigen::Array<type, -1, 1> build_gaussian(const type, const type);

INSTANTIATE_GAUSS(float)
INSTANTIATE_GAUSS(double)
INSTANTIATE_GAUSS(long double)
#undef INSTANTIATE_GAUSS

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
