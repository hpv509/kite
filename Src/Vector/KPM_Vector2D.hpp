/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

#include <iomanip>
template <typename T>
class KPM_Vector <T, 2> : public KPM_VectorBasis <T,2> {
private:
  static const unsigned      D=2u;
  LatticeStructure<2u>        & r;
  std::size_t  	 *MemIndBeg[D][2];
  std::size_t  	 *MemIndEnd[D][2];
  std::size_t         block[D][2];
  std::size_t           tile[D];
  std::size_t    tile_ghosts[D];
  std::size_t       transf_max[D]; // [d][edged]
  std::size_t  transf_bound[D][2]; // [d][edged]
  Hamiltonian<T,2u>          & h;
  T               ***mult_t1_ghost_cor;
  Coordinates<std::size_t,3>   x;
  T                        *phi0;
  T                       *phiM1;
  T                       *phiM2;
  const std::size_t          std;
  const std::size_t          Io;
  const std::size_t          offset;
public:
  static inline constexpr unsigned is_bdg = LatticeStructure<D>::is_bdg;
  using value_type = typename extract_scalar<T>::type;
  T *Fact_Bnd[D][3]; //3 Modos [Salto Positivo, Não Salto, Salto Negativo]
  using KPM_VectorBasis<T,2>::simul;
  using KPM_VectorBasis<T,2>::index;
  using KPM_VectorBasis<T,2>::v;
  using KPM_VectorBasis<T,2>::memory;
  using KPM_VectorBasis<T,2>::aux_wr;
  using KPM_VectorBasis<T,2>::aux_test;
  using KPM_VectorBasis<T,2>::inc_index;
  using KPM_VectorBasis<T,2>::assign_value;
  using KPM_VectorBasis<T,2>::myconj;
  using KPM_VectorBasis<T,2>::multEiphase;
  
  KPM_Vector(int mem, Simulation<T,2> & sim);
  ~KPM_Vector(void);
  void initiate_vector();
  void initiate_phases();
  T get_point();
  void build_wave_packet(const Eigen::Matrix<double, -1, -1> &k, 
                         const Eigen::Matrix<T, -1, -1> &psi0, 
                         const double sigma,
                         const Eigen::Matrix<double, 1, 2> &vb);
  void build_planewave(Eigen::Matrix<double,-1,1> & k, Eigen::Matrix<T,-1,1> & weight);
  void build_site(unsigned long R);

  template < unsigned MULT,bool VELOCITY> 
  void build_regular_phases(int i1, unsigned axis);
  template < unsigned MULT> 
  void initiate_stride(std::size_t & istr);
  template < unsigned MULT> 
  void inline mult_local_disorder(const  std::size_t & j0, const  std::size_t & io);
  void inline mult_regular_hoppings(const  std::size_t & j0, const  std::size_t & io);
  void mult_position(const unsigned, KPM_Vector<T, 2> *);
  template <int S, typename Derived>
  void pairing(const T, Derived&& state_)
    requires Real<T>
  {};
  template <unsigned MULT>
  void mult_bdg_terms(const std::size_t);
  template <unsigned MULT, bool VELOCITY>
  void KPM_MOTOR(KPM_Vector<T, 2> *kpm_final, unsigned axis);

  template <unsigned MULT, bool VELOCITY>
  void multiply_defect(std::size_t, T *&, T *&, unsigned axis);

  void measure_wave_packet(T * bra, T * ket, T * results);
  void Exchange_Boundaries();
  void test_boundaries_system();
  void empty_ghosts(int mem_index);

  template <int S, typename Derived> // S = -1: partition -> lattice, S = 1: lattice -> partition
  void pairing(const T gamma_, Derived&& state_)
    requires Complex<T>
  {
    constexpr value_type norm = 1 / std::sqrt(2);
    const T gd = static_cast<T>(S) * gamma_;
    const T gc = std::conj(gd);
    Coordinates<std::size_t, 3> local(r.Ld);
    for (unsigned io = 0; io < r.Orb; ++io) {
      for (unsigned i1 = NGHOSTS, I1 = r.Ld[1] - NGHOSTS; i1 < I1; ++i1) {
	local.set({NGHOSTS, i1, io});
	unsigned pair_0 = local.index;
	unsigned pair_1 = pair_0 + r.offset;
	for (std::size_t i0 = 0, I0 = r.ld[0]; i0 < I0; ++i0) {
	  const T tmp_1 = state_.coeff(pair_0) + gc * state_.coeff(pair_1);
	  const T tmp_2 = -gd * state_.coeff(pair_0) + state_.coeff(pair_1);
	  state_.coeffRef(pair_0) = norm * tmp_1;
	  state_.coeffRef(pair_1) = norm * tmp_2;
	  ++pair_0;
	  ++pair_1;
	}
      }
    }
#pragma omp barrier
  }
  
  template <int S, typename Derived> // S = -1: partition -> lattice, S = 1: lattice -> partition
  void pairing_nn(const T gamma_, const unsigned hopping_, Derived&& state_)
    requires Complex<T>
  {
    static_assert(S == -1 || S == 1);
    
    constexpr value_type norm = 1 / std::sqrt(2);
    const T gd = static_cast<T>(S) * gamma_;
    const T gc = std::conj(gd);
    
    auto& global_state = simul.Global.nn_pairing_state;
    auto& global_result = simul.Global.nn_pairing_result;
    Coordinates<std::size_t, 3> local(r.Ld);
    Coordinates<std::size_t, 3> global(r.Lt);
    
    for (unsigned io = 0; io < r.Orb; ++io) {
      for (unsigned i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
	for (unsigned i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; ++i0) {
	  local.set({i0, i1, io});
	  r.convertCoordinates(global, local);
	  global_state[global.index] = state_.coeff(local.index);
	  global_state[global.index + r.Sizet] = state_.coeff(local.index + r.offset);
	}
      }
    }
    
#pragma omp barrier
    
#pragma omp master
    {
      std::copy(global_state.begin(), global_state.end(), global_result.begin());
      
      Coordinates<std::ptrdiff_t, 3> bond(r.lB3);
      Coordinates<std::size_t, 3> site_i(r.Lt);
      Coordinates<std::size_t, 3> site_j(r.Lt);
      
      for (unsigned io = 0; io < r.Orb; ++io) {
	if (hopping_ >= h.hr.NHoppings(io)) {
	  continue;
	}
	
	bond.set_coord(h.hr.dist(hopping_, io));
	
	const std::ptrdiff_t dx = bond.coord[0] - 1;
	const std::ptrdiff_t dy = bond.coord[1] - 1;
	const unsigned jo = static_cast<unsigned>(bond.coord[2]);
	
	for (unsigned i1 = 0; i1 < r.Lt[1]; ++i1) {
	  for (unsigned i0 = 0; i0 < r.Lt[0]; ++i0) {
	    
	    std::ptrdiff_t j0 = static_cast<std::ptrdiff_t>(i0) + dx;
	    std::ptrdiff_t j1 = static_cast<std::ptrdiff_t>(i1) + dy;
	    
	    // open-boundaries check
	    if (r.Bd[0] == 0 && (j0 < 0 || j0 >= static_cast<std::ptrdiff_t>(r.Lt[0])))
	      continue;
	    if (r.Bd[1] == 0 && (j1 < 0 || j1 >= static_cast<std::ptrdiff_t>(r.Lt[1])))
	      continue;

	    // periodic wrapping
	    j0 = (j0 + static_cast<std::ptrdiff_t>(r.Lt[0]))
	      % static_cast<std::ptrdiff_t>(r.Lt[0]);
	    j1 = (j1 + static_cast<std::ptrdiff_t>(r.Lt[1]))
	      % static_cast<std::ptrdiff_t>(r.Lt[1]);
	    
	    site_i.set({i0, i1, io});
	    site_j.set({static_cast<std::size_t>(j0), static_cast<std::size_t>(j1), jo});
	    
	    const std::size_t pair_0 = site_i.index;
	    const std::size_t pair_1 = site_j.index + r.Sizet;
	    const std::size_t partition_pair_1 = site_i.index + r.Sizet;

	    std::size_t input_0;
	    std::size_t input_1;
	    std::size_t output_0;
	    std::size_t output_1;
	    
	    if constexpr (S == -1)
	      {
		input_0 = pair_0;
		input_1 = partition_pair_1;
		output_0 = pair_0;
		output_1 = pair_1;
	      }
	    else
	      {
		input_0 = pair_0;
		input_1 = pair_1;
		output_0 = pair_0;
		output_1 = partition_pair_1;
	      }
	    
	    const T tmp_1 = global_state[input_0] + gc * global_state[input_1];
	    const T tmp_2 = -gd * global_state[input_0] + global_state[input_1];
	    
	    global_result[output_0] = norm * tmp_1;
	    global_result[output_1] = norm * tmp_2;
	  }
	}
      }
    }
    
#pragma omp barrier
    
    for (unsigned io = 0; io < r.Orb; ++io) {
      for (unsigned i1 = NGHOSTS; i1 < r.Ld[1] - NGHOSTS; ++i1) {
	for (unsigned i0 = NGHOSTS; i0 < r.Ld[0] - NGHOSTS; ++i0) {
	  local.set({i0, i1, io});
	  r.convertCoordinates(global, local);
	  state_.coeffRef(local.index) = global_result[global.index];
	  state_.coeffRef(local.index + r.offset) = global_result[global.index + r.Sizet];
	}
      }
    }
    
#pragma omp barrier
  }
    
};

