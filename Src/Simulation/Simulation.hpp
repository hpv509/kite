/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/
#ifndef SIMULATION_H_
#define SIMULATION_H_
#include "FFT/GlobalFFT.hpp"

template <typename T, unsigned D>
class Simulation : public ComplexTraits<T> {
public:
  using ComplexTraits<T>::assign_value;
  using ComplexTraits<T>::myconj;
  using value_type = typename extract_scalar<T>::type;
  KPMRandom<T> rnd;
  std::vector<T> ghosts;
  LatticeStructure<D> r;
  GLOBAL_VARIABLES<T> &Global;
  char *name;
  Hamiltonian<T, D> h;
  std::unordered_map<char, unsigned> components_map;

  Simulation(char *, GLOBAL_VARIABLES<T> &, const unsigned, const unsigned);

  //void Measure_Gamma(measurement_queue);

  void Gamma1D(int, int, int, std::vector<std::vector<unsigned>>, std::string);
  void Gamma2D(
    int,
    int,
    std::vector<int>,
    std::vector<std::vector<unsigned>>,
    std::string
  );
  void Gamma3D(
    int,
    int,
    std::vector<int>,
    std::vector<std::vector<unsigned>>,
    std::string
  );
  void GammaGeneral(
    int,
    int,
    std::vector<int>,
    std::vector<std::vector<unsigned>>,
    std::string
  );
  void recursive_KPM(
    int,
    int,
    std::vector<int>,
    long *,
    long *,
    std::vector<std::vector<unsigned>>,
    std::vector<KPM_Vector<T, D> *> *,
    Eigen::Array<T, -1, -1> *
  );
  void store_gamma(
    Eigen::Array<T, -1, -1> *,
    std::vector<int>,
    std::vector<std::vector<unsigned>>,
    std::string
  );
  void store_gamma1D(Eigen::Array<T, -1, -1> *, std::string);
  void store_gamma3D(
    Eigen::Array<T, -1, -1> *,
    std::vector<int>,
    std::vector<std::vector<unsigned>>,
    std::string
  );
  std::vector<std::vector<unsigned>> process_string(std::string);
  double time_kpm(int);

  void calc_singleshot();
  void singleshot(
    Eigen::Array<double, -1, 1> energies,
    Eigen::Array<double, -1, 1> gammas,
    Eigen::Array<int, -1, 1> preserve_disorders,
    Eigen::Array<int, -1, 1> moments,
    int NDisorder,
    int NRandom,
    std::string direction_string
  );

  void calc_conddc();
  void CondDC(int, int, int, int);

  void calc_condopt();
  void CondOpt(int, int, int, int);

  void calc_condopt2();
  void CondOpt2(int, int, int, int, int);

  void calc_DOS();
  void DOS(int, int, int);
  void store_MU(Eigen::Array<T, -1, -1> *);

  void Gaussian_Wave_Packet();
  void calc_wavepacket();

  void LMU(int, int, Eigen::Array<unsigned long, -1, 1>);
  void calc_LDOS();
  void store_LMU(Eigen::Array<T, -1, -1> *);

  void calc_ARPES();
  void ARPES(
    int NDisorder,
    int NMoments,
    Eigen::Array<double, -1, -1> &k_vectors,
    Eigen::Matrix<T, -1, 1> &weight
  );
  void store_ARPES(Eigen::Array<T, -1, -1> *);

  // Custom Rank One
  static herr_t
  getMembers(H5::Group &group, const std::string &name, void *op_data);
  void calc_custom_one();
  void custom_one(
    const int,
    const int,
    const int,
    const Eigen::Array<T, -1, 1> &,
    const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &,
    const std::vector<std::string> &
  );
  void store_custom_one(const Eigen::Matrix<T, -1, 1> &, const unsigned);
  void multiply_orb_mtx(
    const Eigen::Matrix<std::complex<double>, -1, -1> &,
    const KPM_Vector<T, D> *,
    KPM_Vector<T, D> *
  );
  void act_with_stream(
    const std::vector<std::string> &,
    const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &,
    const Eigen::Array<T, -1, 1> &,
    std::array<KPM_Vector<T, D> *, 2> &,
    const unsigned
  );

  // Custom Rank One - Local
  void calc_custom_one_local();
  void custom_one_local(
    const int,
    const int,
    const int,
    const std::vector<int> &,
    const std::vector<int> &,
    const Eigen::Array<T, -1, 1> &,
    const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &,
    const std::vector<std::string> &
  );
  void store_custom_one_local(const Eigen::Matrix<T, -1, -1> &, const unsigned);

  // Custom Rank Two
  void calc_custom_two();
  void custom_two(
    const int,
    const int,
    const std::vector<int> &,
    const std::vector<std::vector<std::string>> &,
    const std::vector<Eigen::Array<T, -1, 1>> &,
    const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &
  );
  void store_custom_two(
    Eigen::Array<T, -1, -1> &,
    const unsigned,
    const std::vector<int> &
  );
  // Custom Rank Two - SS
  void calc_custom_ss_two()
    requires Complex<T>;
  void calc_custom_ss_two()
    requires Real<T>
  {};
  void custom_ss_two(
    const int,
    const int,
    const value_type,
    const value_type,
    const Eigen::Array<value_type, -1, 1> &energies_,
    const std::vector<std::vector<std::string>> &,
    const std::vector<Eigen::Array<T, -1, 1>> &,
    const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &
  )
    requires Complex<T>;
  void store_custom_ss_two(const Eigen::Array<T, -1, -1> &, const unsigned)
    requires Complex<T>;

  void calc_custom_two_local();
  void custom_two_local(
    const Eigen::Array<unsigned, -1, -1> &,
    const std::array<unsigned, 2> &,
    const std::vector<std::vector<std::string>> &,
    const std::vector<Eigen::Array<T, -1, 1>> &,
    const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &
  );
  void store_custom_two_local(
    const unsigned,
    const std::array<unsigned, 2> &,
    const Eigen::Array<unsigned, 1, D + 1> &,
    const Eigen::Array<T, -1, -1> &
  );
  // LDoS
  void calc_ldos()
    requires Complex<T>;
  void calc_ldos()
    requires Real<T>
  {};
  void ldos(const int, const value_type, const value_type, const int)
    requires Complex<T>;
  void store_ldos(const Eigen::Array<value_type, -1, -1> &)
    requires Complex<T>;
  // Spectral Function map
  void calc_spectral(GlobalFFT<value_type> &)
    requires Complex<T>;
  void calc_spectral(GlobalFFT<value_type> &)
    requires Real<T>
  {};
  void spectral(
    const unsigned,
    const value_type,
    const value_type,
    const unsigned,
    GlobalFFT<value_type> &
  )
    requires Complex<T>;
  void store_spectral(const Eigen::Array<value_type, -1, -1> &)
    requires Complex<T>;

  // LCM
  void calc_lcm();
  void lcm(
    const unsigned,
    const value_type,
    const value_type,
    const std::array<int, 2> &
  );
  void store_lcm(const Eigen::Array<T, -1, 1> &);

  // STLCM
  void calc_st_lcm();
  void st_lcm(const unsigned, const value_type, const value_type);
  void store_st_lcm(const Eigen::Array<T, -1, -1> &);

  // Wave-Packet Time Evol
  void calc_localized_wavepacket()
    requires Complex<T>;
  void calc_localized_wavepacket()
    requires Real<T>
  {};
  void localized_wavepacket(
    const value_type t,
    const unsigned measurements,
    const unsigned num_spectral_moments,
    const std::array<value_type, D + 1> &pos_,
    const std::array<value_type, 2> &energy_window,
    const std::array<value_type, D> &k0_,
    const value_type width,
    const std::size_t num_probes,
    const Eigen::Array<std::size_t, -1, D + 1> &prop_coords,
    const std::size_t sample_start,
    const std::size_t sample_L
  )
    requires Complex<T>;
  void store_localized_wavepacket(
    const Eigen::Array<value_type, -1, -1> &spectral_moments,
    const Eigen::Array<value_type, -1, -1> &moments1,
    const Eigen::Array<value_type, -1, -1> &moments2,
    const Eigen::Array<T, -1, 1> &return_amplitudes,
    const Eigen::Array<value_type, -1, -1> &transmission_weights,
    const Eigen::Array<T, -1, -1> &propagator_amplitudes,
    const std::vector<std::pair<std::size_t, std::size_t>> &propagator_coords,
    const std::size_t num_global_probes
  )
    requires Complex<T>;

  // S-Wave
  void calc_swave()
    requires Complex<T>;
  void calc_swave()
    requires Real<T>
  {};
  
  void s_wave(
    const int,
    const value_type,
    const value_type,
    const value_type,
    const value_type,
    const value_type
  )
    requires Complex<T>;
};

#endif
