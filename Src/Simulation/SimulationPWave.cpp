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
void Simulation<T, D>::calc_pwave()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_pwave\n");
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_p_wave =false;
    try {
      int dummy_variable;
      get_hdf5<
	int>(&dummy_variable, file, (char *)"/Calculation/p_wave/NumRandoms");
      Global.calculate_p_wave = true;
    } catch (H5::Exception &e) {
      debug_message("p_wave: no need to calculate it.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_p_wave = false;
#pragma omp critical
  local_calculate_p_wave = Global.calculate_p_wave;
#pragma omp barrier
  if (local_calculate_p_wave) {
#pragma omp master
    std::cout << "Calculating PWave.\n";
#pragma omp barrier
    int randoms;
    value_type beta;
    value_type mu;
    value_type u;
    value_type v;
    value_type gamma;
    value_type s_delta;
    value_type nn_delta;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      get_hdf5<int>(&randoms, file, (char *)"/Calculation/p_wave/NumRandoms");
      get_hdf5<value_type>(&beta, file, (char *)"/Calculation/p_wave/Beta");
      get_hdf5<value_type>(&mu, file, (char *)"/Calculation/p_wave/ChemicalPotential");
      get_hdf5<value_type>(&u, file, (char *)"/Calculation/p_wave/U");
      get_hdf5<value_type>(&v, file, (char *)"/Calculation/p_wave/V");
      get_hdf5<value_type>(&gamma, file, (char *)"/Calculation/p_wave/Gamma");
      get_hdf5<value_type>(&s_delta, file, (char *)"/Calculation/p_wave/SDelta");
      get_hdf5<value_type>(&nn_delta, file, (char *)"/Calculation/p_wave/NNDelta");

      file->close();
      delete file;
    }
    p_wave(randoms, beta, mu, u, v, gamma, s_delta, nn_delta);
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::p_wave(
  const int randoms_,
  const value_type beta_,
  const value_type mu_,
  const value_type u_, // on-site
  const value_type v_, // nearest-neighbor
  const value_type gamma_,
  const value_type s_delta_,
  const value_type nn_delta_
)
  requires Complex<T>
{
  debug_message("Entered PWave\n");
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
  const value_type v = v_ / energy_scale;
  const value_type gamma = gamma_ / energy_scale;
  const value_type s_delta = s_delta_ / energy_scale;
  const value_type nn_delta = nn_delta_ / energy_scale;

  const value_type size = r.Sizet - r.SizetVacancies;

  const Eigen::Array<value_type, -1, 1> coefs =
    Coefficients::build_fermi_sqrt<value_type>(beta, value_type(0.0));

  h.bdg.init_fields(s_delta, nn_delta, gamma, mu);
  
  int iteration = 1;
  const int max_iteration = 1024;
  
  Eigen::Array<T, -1, 1> sum_s_delta(r.Sized, 1);
  Eigen::Array<T, -1, -1> sum_nn_delta(h.bdg.nn_delta.rows(), r.Sized);
  sum_s_delta.setZero();
  sum_nn_delta.setZero();

  value_type weight_avg = 1.0;
  value_type weight_sum = 0.0;
  value_type weight_r = 1.0;
  value_type weight_alpha = 1.0;

#pragma omp barrier

  h.generate_disorder();
  KPM_Vector<T, D> phi(2, *this);
  Eigen::Array<T, -1, 1> ket(2 * r.Sized);
  Eigen::Array<T, -1, 1> bra(2 * r.Sized);
  Eigen::Array<value_type, -1, 1> density(2 * r.Sized);

  Eigen::Array<T, -1, 1> results_s_delta(r.Sized, 1);
  Eigen::Array<T, -1, -1> results_nn_delta(h.bdg.nn_delta.rows(), r.Sized);
  Eigen::Array<T, -1, -1> raw_nn_delta(h.bdg.nn_delta.rows(), r.Sized);

  Eigen::Array<int, -1, -1> reverse_hopping(h.bdg.nn_delta.rows(), r.Orb);
  reverse_hopping.setConstant(-1);
  
  unsigned cell_basis = 1;
  for (unsigned d = 0; d < D; ++d)
    cell_basis *= 3;
  
  for (unsigned io = 0; io < r.Orb; ++io)
    for (unsigned ib = 0; ib < h.bdg.nn_delta.rows() && ib < h.hr.NHoppings(io); ++ib)
      {
	const unsigned d = h.hr.dist(ib, io);
	const unsigned final_io = d / cell_basis;
	const unsigned spatial_d = d % cell_basis;
	const unsigned reverse_spatial_d = cell_basis - 1 - spatial_d;
	const unsigned reverse_d = io * cell_basis + reverse_spatial_d;
	
	for (unsigned jb = 0; jb < h.bdg.nn_delta.rows() && jb < h.hr.NHoppings(final_io); ++jb)
	  if (h.hr.dist(jb, final_io) == reverse_d)
	    {
	      reverse_hopping(ib, io) = static_cast<int>(jb);
	      break;
	    }

	if (reverse_hopping(ib, io) < 0)
	  {
	    std::cerr
	      << "Could not find reverse hopping for "
	      << "orbital " << io
	      << ", hopping " << ib
	      << std::endl;
	    exit(1);
	  }
      }
  
  while (iteration <= max_iteration)
    {
      results_s_delta.setZero();
      results_nn_delta.setZero();

      for (int vec = 0; vec < randoms_; ++vec)
	{
	  h.generate_twists();
	  
	  // off-diagonal elements for on-site
	  
	  phi.initiate_phases();
	  phi.set_index(0);
	  phi.initiate_vector();
	  phi.v.col(0) *= std::sqrt(size);
	  bra = phi.v.col(0);
	  
	  // gamma = 1
	  phi.template pairing<-1>(T(1.0), phi.v.col(0));
	  
	  ket.setZero();
	  phi.Exchange_Boundaries();
	  
	  for (unsigned n = 0, N = coefs.size(); n < N; ++n)
	    {
	      phi.cheb_iteration(n);
	      ket += coefs(n) * phi.v.col(phi.get_index()).array();
	    }
	  
	  phi.template pairing<1>(T(1.0), ket);

	  density = (bra.conjugate() * ket).abs2();
	  const Eigen::Array<value_type, -1, 1> upsilon_1 =
	    density.head(r.Sized) - density.tail(r.Sized);
	  
	  phi.v.setZero();
	  phi.set_index(0);
	  phi.v.col(0) = bra.matrix();
	  
	  // gamma = i
	  phi.template pairing<-1>(T(0.0, 1.0), phi.v.col(0));
	  
	  ket.setZero();
	  phi.Exchange_Boundaries();
	  
	  for (unsigned n = 0, N = coefs.size(); n < N; ++n)
	    {
	      phi.cheb_iteration(n);
	      ket += coefs(n) * phi.v.col(phi.get_index()).array();
	    }
	  
	  phi.template pairing<1>(T(0.0, 1.0), ket);

	  density = (bra.conjugate() * ket).abs2();
	  const Eigen::Array<value_type, -1, 1> upsilon_i =
	    density.head(r.Sized) - density.tail(r.Sized);
	  
	  Eigen::Array<T, -1, 1> map_s_delta(r.Sized);
	  map_s_delta.real() = value_type(0.5) * u * upsilon_1;
	  map_s_delta.imag() = - value_type(0.5) * u * upsilon_i;
	  
	  const value_type weight = 1.0 / (vec + 1);
	  results_s_delta += weight * (map_s_delta - results_s_delta);
	  
	  // off-diagonal for nn
	  
	  for (unsigned hop = 0; hop < h.bdg.nn_delta.rows(); ++hop)
	    {
	      phi.v.setZero();
	      phi.set_index(0);
	      phi.v.col(0) = bra.matrix();
	      
	      // gamma = 1
	      phi.template pairing_nn<-1>(T(1.0), hop, phi.v.col(0));
	      
	      ket.setZero();
	      phi.Exchange_Boundaries();
	      
	      for (unsigned n = 0, N = coefs.size(); n < N; ++n)
		{
		  phi.cheb_iteration(n);
		  ket += coefs(n) * phi.v.col(phi.get_index()).array();
		}
	      
	      phi.template pairing_nn<1>(T(1.0), hop, ket);

	      density = (bra.conjugate() * ket).abs2();
	      const Eigen::Array<value_type, -1, 1> upsilon_1 =
		density.head(r.Sized) - density.tail(r.Sized);
	      
	      phi.v.setZero();
	      phi.set_index(0);
	      phi.v.col(0) = bra.matrix();
	      
	      //gamma = i
	      phi.template pairing_nn<-1>(T(0.0, 1.0), hop, phi.v.col(0));
	      
	      ket.setZero();
	      phi.Exchange_Boundaries();
	      
	      for (unsigned n = 0, N = coefs.size(); n < N; ++n)
		{
		  phi.cheb_iteration(n);
		  ket += coefs(n) * phi.v.col(phi.get_index()).array();
		}
	      
	      phi.template pairing_nn<1>(T(0.0, 1.0), hop, ket);

	      density = (bra.conjugate() * ket).abs2();
	      const Eigen::Array<value_type, -1, 1> upsilon_i =
		density.head(r.Sized) - density.tail(r.Sized);
	      
	      Eigen::Array<T, -1, 1> map_nn_delta(r.Sized);
	      map_nn_delta.real() = value_type(0.5) * v * upsilon_1;
	      map_nn_delta.imag() = - value_type(0.5) * v * upsilon_i;
	      
	      results_nn_delta.row(hop) += weight * (map_nn_delta.transpose() - results_nn_delta.row(hop));
	    } // end of cycle for hoppings
	} // end of cycle for random vectors

      // exchanging boundaries
      raw_nn_delta = results_nn_delta;
      for (unsigned ib = 0; ib < h.bdg.nn_delta.rows(); ++ib)
	{
	  phi.v.setZero();
	  phi.set_index(0);
	  phi.v.col(0).head(r.Sized) = raw_nn_delta.row(ib).transpose().matrix();
	  phi.Exchange_Boundaries();
	  raw_nn_delta.row(ib) = phi.v.col(0).head(r.Sized).array().transpose();
	}
      
      Coordinates<std::size_t, D + 1> local(r.Ld);
      for (unsigned io = 0; io < r.Orb; ++io)
	for (unsigned ib = 0; ib < h.bdg.nn_delta.rows() && ib < h.hr.NHoppings(io); ++ib)
	  {
	    const unsigned jb = static_cast<unsigned>(reverse_hopping(ib, io));
	    const std::ptrdiff_t d1 = h.hr.distance(ib, io);
	    
	    if constexpr (D == 2)
	      for (unsigned i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1)
		{
		  local.set({NGHOSTS, i1, io});
		  const std::size_t j0 = local.index;
		  const std::size_t j1 = j0 + r.ld[0];
		  
		  for (std::size_t i = j0; i < j1; ++i)
		    {
		      const std::size_t j = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(i) + d1);
		      results_nn_delta(ib, i) = value_type(0.5) * (raw_nn_delta(ib, i) + raw_nn_delta(jb, j));
		    }
		}

	    if constexpr (D == 3)
	      for (unsigned i2 = NGHOSTS; i2 < r.Ld[2] - NGHOSTS; ++i2)
		for (unsigned i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1)
		  {
		    local.set({NGHOSTS, i1, i2, io});
		    const std::size_t j0 = local.index;
		    const std::size_t j1 = j0 + r.ld[0];
		    
		    for (std::size_t i = j0; i < j1; ++i)
		      {
			const std::size_t j = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(i) + d1);
			results_nn_delta(ib, i) = value_type(0.5) * (raw_nn_delta(ib, i) + raw_nn_delta(jb, j));
		      }
		  }
	  }
      
      weight_avg *= 1.0 + weight_r / std::pow(value_type(iteration), weight_alpha);
      weight_sum += weight_avg;
      
      sum_s_delta += weight_avg * results_s_delta;
      sum_nn_delta += weight_avg * results_nn_delta;
      
      h.bdg.s_delta = sum_s_delta / weight_sum;
      h.bdg.nn_delta = sum_nn_delta / weight_sum;
      
      ++iteration;
    } // end of cycle of iterations

  store_p_wave(energy_scale);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_p_wave(const value_type energy_scale_)
  requires Complex<T>
{
  debug_message("Entered store_pwave\n");
  Coordinates<std::size_t, D + 1> global(r.Lt);
  Coordinates<std::size_t, D + 1> local(r.Ld);
#pragma omp master
  Global.p_wave_map.resize(r.Sizet, 2 + 2 * h.bdg.nn_delta.rows());
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
      Global.p_wave_map(global.index, 0) =
	std::real(h.bdg.s_delta(local.index)) * energy_scale_;
      Global.p_wave_map(global.index, 1) =
	std::imag(h.bdg.s_delta(local.index)) * energy_scale_;

      for (unsigned hop = 0; hop < h.bdg.nn_delta.rows(); ++hop)
	{
	  Global.p_wave_map(global.index, 2 + 2*hop) =
	    std::real(h.bdg.nn_delta(hop, local.index)) * energy_scale_;
	  Global.p_wave_map(global.index, 3 + 2*hop) =
	    std::imag(h.bdg.nn_delta(hop, local.index)) * energy_scale_;
	}
    };
    UnitCellLoop<D>::run(idx, start, final, body);
  }
#pragma omp barrier
#pragma omp master
  {
    const Eigen::Array<value_type, -1, -1> p_wave_r = Global.p_wave_map.real();
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDWR);
    char buffer[200];
    std::sprintf(buffer, "/Calculation/p_wave/Map");
    const std::string name(buffer);
    write_hdf5(p_wave_r, file, name);
    delete file;
  }
#pragma omp barrier
  debug_message("Left store_pwave\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
