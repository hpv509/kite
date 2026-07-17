/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/
#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "Global.hpp"
#include "Coordinates.hpp"
#include "LatticeStructure.hpp"
#include "myHDF5.hpp"
#include "Random.hpp"
template <typename T, unsigned D>
class Hamiltonian;
template <typename T, unsigned D>
class KPM_Vector;
#include "Simulation.hpp"
#include "SimulationGlobal.hpp"
#include "Hamiltonian.hpp"
template <typename T, unsigned D>
GlobalSimulation<T, D>::GlobalSimulation(char *name) : rglobal(name)
{
  debug_message("Entered global_simulation\n");

  // rglobal is an instance of Lattice Structure which contains all the information
  // about the periodic part of the lattice and the unit cell. It contains the
  // lattice vectors, the number of threads, the size of the lattice inside each
  // thread, thread id, orbital positions, etc

  // Global is an instance of GLOBAL_VARIABLES which contains global objects to be
  // shared among all threads

  Global.ghosts.resize(rglobal.get_BorderSize());
  std::fill(Global.ghosts.begin(), Global.ghosts.end(), 0);

  std::string path;
  H5::H5File file(name, H5F_ACC_RDONLY);
  path = "/EnergyScale";
  get_hdf5<double>(&EnergyScale, &file, path);
  unsigned s0;
  path = "/Seed0";
  get_hdf5<unsigned>(&s0, &file, path);
  unsigned s1;
  path = "/Seed1";
  get_hdf5<unsigned>(&s1, &file, path);
  file.close();
  omp_set_num_threads(rglobal.n_threads);

  using value_type = typename extract_scalar<T>::type;
  GlobalFFT<value_type> global_fft;
  global_fft.allocate(rglobal.Sizet);

  debug_message("Starting parallelization\n");
#pragma omp parallel default(shared)
  {
    const unsigned thread_id = omp_get_thread_num();
    const unsigned seed_h = (s0 != 0) * 2654435761 * (thread_id + 8191 * s0);
    const unsigned seed_v = (s1 != 0) * 3242174893 * (thread_id + 7919 * s1);
    Simulation<T, D> simul(name, Global, seed_v, seed_h);
    simul.calc_conddc();
    simul.calc_condopt();
    simul.calc_condopt2();
    simul.calc_singleshot();
    simul.calc_DOS();
    simul.calc_wavepacket();
    simul.calc_localized_wavepacket();
    simul.calc_LDOS();
    simul.calc_ARPES(); // fetches parameters from .h5 file and calculates ARPES
    simul.calc_ldos();
    simul.calc_spectral(global_fft);
    simul.calc_custom_one();
    simul.calc_custom_one_local();
    simul.calc_custom_two();
    simul.calc_custom_two_local();
    simul.calc_custom_ss_two();
    simul.calc_lcm();
    simul.calc_st_lcm();
  }
  debug_message("Left global_simulation\n");
}

#define instantiate(type, dim) template class GlobalSimulation<type, dim>;
#include "instantiate.hpp"
