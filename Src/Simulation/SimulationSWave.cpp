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
    Global.calculate_s_wave =false;
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
  
  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(beta, value_type(0.0));
  
  h.bdg.init_fields(delta, gamma, mu);

  int iteration = 1;
  const int max_iteration = 128;
  Eigen::Array<value_type, -1, 1> sum_gamma(r.Sized, 1);
  Eigen::Array<T, -1, 1> sum_delta(r.Sized, 1);
  sum_gamma.setZero();
  sum_delta.setZero();

  value_type weight_avg = 1.0;
  value_type weight_sum = 0.0;
  value_type weight_r = 1.0;
  value_type weight_alpha = 1.0;

#pragma omp barrier

  h.generate_disorder();
  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(2 * r.Sized);
  Eigen::Array<T, -1, 1> bra(2 * r.Sized);

  Eigen::Array<value_type, -1, 1> results_gamma(r.Sized, 1);
  Eigen::Array<T, -1, 1> results_delta(r.Sized, 1);
  
  while (iteration <= max_iteration)
    {

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
	for (unsigned n = 0, N = coefs.size(); n < N; ++n)
	  {
	    phi.cheb_iteration(n);
	    ket += coefs(n) * phi.v.col(phi.get_index()).array();
	  }
	
	const Eigen::Array<value_type, -1, 1> map_gamma =
	  u * (bra.conjugate() * ket).abs2().head(r.Sized) - mu;
	const value_type weight = 1.0 / (vec + 1);
	results_gamma += weight * (map_gamma - results_gamma);
	
	// off-diagonal elements

	// gamma = 1
	
	phi.v.setZero();
	phi.set_index(0);
	phi.v.col(0) = bra.matrix();
	
	phi.template pairing<-1>(T(1.0), phi.v.col(0));
	
	ket.setZero();
	phi.Exchange_Boundaries();
	
	for (unsigned n = 0, N = coefs.size(); n < N; ++n)
	  {
	    phi.cheb_iteration(n);
	    ket += coefs(n) * phi.v.col(phi.get_index()).array();
	  }
	
	phi.template pairing<1>(T(1.0), ket);
	
	const Eigen::Array<value_type, -1, 1> upsilon_1 =
	  (bra.conjugate() * ket).abs2().head(r.Sized)
	  -
	  (bra.conjugate() * ket).abs2().tail(r.Sized);

	// gamma = i

	phi.v.setZero();
	phi.set_index(0);
	phi.v.col(0) = bra.matrix();
	
	phi.template pairing<-1>(T(0.0, 1.0), phi.v.col(0));
	
	ket.setZero();
	phi.Exchange_Boundaries();
	
	for (unsigned n = 0, N = coefs.size(); n < N; ++n)
	  {
	    phi.cheb_iteration(n);
	    ket += coefs(n) * phi.v.col(phi.get_index()).array();
	  }
	
	phi.template pairing<1>(T(0.0, 1.0), ket);
	
	const Eigen::Array<value_type, -1, 1> upsilon_i =
	  (bra.conjugate() * ket).abs2().head(r.Sized)
	  -
	  (bra.conjugate() * ket).abs2().tail(r.Sized);

	Eigen::Array<T, -1, 1> map_delta(r.Sized);
	map_delta.real() = value_type(0.5) * u * upsilon_1;
	map_delta.imag() = - value_type(0.5) * u * upsilon_i;
	
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

  h.bdg.hartree += mu;
  store_s_wave(energy_scale);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_s_wave(const value_type energy_scale_)
  requires Complex<T>
{
  debug_message("Entered store_swave\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  Global.s_wave_map.resize(r.Sizet, 3);
#pragma omp barrier
  std::array<unsigned, D> idx;
  std::array<unsigned, D> start;
  std::array<unsigned, D> final;
  for (unsigned d = 0; d < D; ++d) {
    start[d] = NGHOSTS;
    final[d] = r.Ld[D - 1 - d] - NGHOSTS;
  }
  for (unsigned io = 0, Io = r.Orb; io < Io; ++io) {
    auto body = [&](const std::array<unsigned, D> &i) {
      if constexpr (D == 2)
        local.set({i[1], i[0], io});
      else if constexpr (D == 3)
        local.set({i[2], i[1], i[0], io});
      r.convertCoordinates(global, local);
      Global.s_wave_map(global.index,0) =
	h.bdg.hartree(local.index) * energy_scale_;
      Global.s_wave_map(global.index,1) =
	std::real(h.bdg.s_delta(local.index)) * energy_scale_;
      Global.s_wave_map(global.index,2) =
	std::imag(h.bdg.s_delta(local.index)) * energy_scale_;
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  {
    const Eigen::Array<value_type, -1, -1> s_wave_r = Global.s_wave_map.real();
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    char buffer[200];
    std::sprintf(buffer, "/Calculation/s_wave/Map");
    const std::string name(buffer);
    write_hdf5(s_wave_r, file, name);
    delete file;
  }
#pragma omp barrier
  debug_message("Left store_swave\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
