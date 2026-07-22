/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

template <typename T>
struct GLOBAL_VARIABLES {
  std::vector<T> ghosts;
  std::vector<std::size_t>    element1;
  std::vector<std::ptrdiff_t> element2_diff;
  std::vector<T> hopping;
  std::vector<std::size_t>    element;
  std::vector<T> U;
  T soma; 

  // Averages
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> mu;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> lambda;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> singleshot_cond;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> general_gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> smaller_gamma;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_x;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_y;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_z;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_ident;
  Eigen::Array <T, Eigen::Dynamic, Eigen::Dynamic> avg_results;
  Eigen::Array <double,3,1> GlobBTwist; // Glob Boundary Twist Angles
  Eigen::Array<T, -1, -1> ldos_map;
  Eigen::Array<T, -1, -1> s_wave_map;
  Eigen::Array<T, -1, -1> localized_wavepacket;
  Eigen::Array<T, -1, -1> results_1;
  Eigen::Array<T, -1, -1> results_2;
  Eigen::Array<T, -1, -1> results_3;
  Eigen::Array<T, -1, -1> results_4;
  Eigen::Array<T, -1, -1> results_5;
  Eigen::Array<T, -1, -1> results_6;
  // Statistics
  std::vector<T> buffer;

  bool calculate_arpes;
  bool calculate_ldos;
  bool calculate_wavepacket;
  bool calculate_localized_wavepacket;
  bool calculate_dos;
  bool calculate_conddc;
  bool calculate_condopt;
  bool calculate_condopt2;
  bool calculate_singleshot;
  bool calculate_ldos_map;
  bool calculate_spectral_map;
  bool calculate_custom_one;
  bool calculate_custom_two;
  bool calculate_custom_two_local;
  bool calculate_custom_ss_two;
  bool calculate_lcm;
  bool calculate_st_lcm;
  bool calculate_s_wave;

  GLOBAL_VARIABLES();
  void addbond ( std::size_t, std::ptrdiff_t, T );
  void addlocal( std::size_t,  T);
};
