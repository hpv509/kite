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
#include "Loop.hpp"
#include "Coefficients.hpp"

template <typename T, unsigned D>
void Simulation<T, D>::calc_swave()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_swave\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_s_wave = false;
    try {
      int dummy_variable;
      get_hdf5<
	int>(&dummy_variable, file, (char *)"/Calculation/s_wave/NumRandoms");
      Global.calculate_s_wave = true;
    } catch (H5::Exception &e) {
      debug_message("s_wave: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_s_wave = false;
#pragma omp critical
  local_calculate_s_wave = Global.calculate_s_wave;
#pragma omp barrier
  if (local_calculate_s_wave) {
#pragma omp master
    std::cout << "Calculating SWave.\n";
#pragma omp barrier
    int randoms;
    value_type beta;
    value_type mu;
    value_type u;
    value_type gamma;
    value_type delta;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<int>(&randoms, file, (char *)"/Calculation/s_wave/NumRandoms");
      get_hdf5<value_type>(&beta, file, (char *)"/Calculation/s_wave/Beta");
      get_hdf5<value_type>(&mu, file, (char *)"/Calculation/s_wave/ChemicalPotential");
      get_hdf5<value_type>(&u, file, (char *)"/Calculation/s_wave/U");
      get_hdf5<value_type>(&gamma, file, (char *)"/Calculation/s_wave/Gamma");
      get_hdf5<value_type>(&delta, file, (char *)"/Calculation/s_wave/Delta");

      file->close();
      delete file;
    }
    s_wave(randoms, beta, mu, u, gamma, delta);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::s_wave(
  const int randoms_,
  const value_type beta_,
  const value_type mu_,
  const value_type u_,
  const value_type gamma_,
  const value_type delta_
)
  requires Complex<T>
{
  debug_message("Entered SWave\n");
  value_type energy_scale;
#pragma omp critical
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    get_hdf5<value_type>(&energy_scale, file, (char *)"/EnergyScale");
    file->close();
    delete file;
  }
#pragma omp barrier
  const value_type beta = beta_ * energy_scale;
  const value_type mu = mu_ / energy_scale;
  const value_type u = u_ / energy_scale;
  const value_type gamma = gamma_ / energy_scale;
  const value_type delta = delta_ / energy_scale;

  const value_type size = r.Sizet - r.SizetVacancies;

  // falta escrever esta função
  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(beta, value_type(0.0));
  
  h.bdg.init_fields(delta, gamma, mu);

  int iteration = 1;
  const int max_iteration = 100000;
  Eigen::Array<value_type, -1, 1> sum_gamma(r.Sized / 2, 1);
  Eigen::Array<value_type, -1, 1> sum_delta(r.Sized / 2, 1);
  sum_gamma.setZero();
  sum_delta.setZero();

  value_type weight_avg = 1.0;
  value_type weight_sum = 0.0;
  value_type weight_r = 1.0;
  value_type weight_alpha = 1.0;

#pragma omp barrier

  h.generate_disorder();
  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(r.Sized);
  Eigen::Array<T, -1, 1> bra(r.Sized);

  // armazena no fim de cada ciclo de vetores aleatórios
  Eigen::Array<value_type, -1, 1> results_gamma(r.Sized / 2, 1);
  Eigen::Array<value_type, -1, 1> results_delta(r.Sized / 2, 1);

  while (iteration < max_iteration){

    results_gamma.setZero();
    results_delta.setZero();
    
    for (int vec = 0; vec < randoms_; ++vec) {
      h.generate_twists();
      
      // diagonal elements
      
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
      
      const Eigen::Array<value_type, -1, 1> map_gamma =
	u * (bra.conjugate() * ket).abs2().head(r.Sized/2);
      const value_type weight = 1.0 / (vec + 1);
      results_gamma += weight * (map_gamma - results_gamma);
      
      // off-diagonal elements
      
      phi.v.setZero();
      phi.set_index(0);
      phi.v.col(0) = bra.matrix();
      
      phi.template pairing<-1>(T(1.0), phi.v.col(0));
      
      ket.setZero();
      phi.Exchange_Boundaries();
      
      for (unsigned n = 0, N = coefs.size(); n < N; ++n) {
	phi.cheb_iteration(n);
	ket += coefs(n) * phi.v.col(phi.get_index()).array();
      }
      
      phi.template pairing<1>(T(1.0), ket);
      
      const Eigen::Array<value_type, -1, 1> map_delta =
	value_type(0.5) * u * ( (bra.conjugate() * ket).abs2().head(r.Sized/2)
				- (bra.conjugate() * ket).abs2().tail(r.Sized/2) );
      results_delta += weight * (map_delta - results_delta);
    }
    
    weight_avg *= 1.0 + weight_r / std::pow(value_type(iteration), weight_alpha);
    weight_sum += weight_avg;
    
    sum_gamma += weight_avg * results_gamma;
    sum_delta += weight_avg * results_delta;
    
    h.bdg.hartree = sum_gamma / weight_sum;
    h.bdg.s_delta = sum_delta / weight_sum;

    ++iteration;
  }
}
