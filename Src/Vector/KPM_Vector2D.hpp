/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

#include <iomanip>
#include "Pairing.hpp"
template <typename T>
class KPM_Vector<T, 2> : public KPM_VectorBasis<T, 2> {
private:
  static const unsigned D = 2u;
  LatticeStructure<2u> &r;
  std::size_t *MemIndBeg[D][2];
  std::size_t *MemIndEnd[D][2];
  std::size_t block[D][2];
  std::size_t tile[D];
  std::size_t tile_ghosts[D];
  std::size_t transf_max[D];      // [d][edged]
  std::size_t transf_bound[D][2]; // [d][edged]
  Hamiltonian<T, 2u> &h;
  T ***mult_t1_ghost_cor;
  Coordinates<std::size_t, 3> x;
  T *phi0;
  T *phiM1;
  T *phiM2;
  T *pair_buf; // buffer for bond pairing rotations
  const std::size_t std;
  const std::size_t Io;
  const std::size_t offset;

public:
  static inline constexpr unsigned is_bdg = pairing::is_bdg;
  using value_type = typename extract_scalar<T>::type;
  T *Fact_Bnd[D][3]; //3 Modos [Salto Positivo, Não Salto, Salto Negativo]
  using KPM_VectorBasis<T, 2>::simul;
  using KPM_VectorBasis<T, 2>::index;
  using KPM_VectorBasis<T, 2>::v;
  using KPM_VectorBasis<T, 2>::memory;
  using KPM_VectorBasis<T, 2>::aux_wr;
  using KPM_VectorBasis<T, 2>::aux_test;
  using KPM_VectorBasis<T, 2>::inc_index;
  using KPM_VectorBasis<T, 2>::assign_value;
  using KPM_VectorBasis<T, 2>::myconj;
  using KPM_VectorBasis<T, 2>::multEiphase;

  KPM_Vector(int mem, Simulation<T, 2> &sim);
  ~KPM_Vector(void);
  void initiate_vector();
  void initiate_phases();
  T get_point();
  void build_wave_packet(
    const Eigen::Matrix<double, -1, -1> &k,
    const Eigen::Matrix<T, -1, -1> &psi0,
    const double sigma,
    const Eigen::Matrix<double, 1, 2> &vb
  );
  void build_planewave(
    Eigen::Matrix<double, -1, 1> &k,
    Eigen::Matrix<T, -1, 1> &weight
  );
  void build_site(unsigned long R);

  template <unsigned MULT, bool VELOCITY>
  void build_regular_phases(int i1, unsigned axis);
  template <unsigned MULT>
  void initiate_stride(std::size_t &istr);
  template <unsigned MULT>
  void inline mult_local_disorder(const std::size_t &j0, const std::size_t &io);
  template <unsigned MULT>
  void inline mult_regular_hoppings(
    const std::size_t &j0,
    const std::size_t &io
  );
  void mult_position(const unsigned, KPM_Vector<T, 2> *);
  template <int S, typename Derived>
  void pairing(const T, Derived &&state_)
    requires Real<T>
  {};
  template <unsigned MULT>
  void mult_diag_bdg_terms(const std::size_t);
  template <unsigned MULT>
  void mult_pairing_bonds(const std::size_t, const std::size_t)
    requires Complex<T>;
  template <unsigned MULT>
  void mult_pairing_bonds(const std::size_t, const std::size_t)
    requires Real<T>
  {};
  template <unsigned MULT, bool VELOCITY>
  void KPM_MOTOR(KPM_Vector<T, 2> *kpm_final, unsigned axis);

  template <unsigned MULT, bool VELOCITY>
  void multiply_defect(std::size_t, T *&, T *&, unsigned axis);

  void measure_wave_packet(T *bra, T *ket, T *results);
  void Exchange_Boundaries();
  void test_boundaries_system();
  void empty_ghosts(int mem_index);

  static inline void
  rotate_pair(T &p, T &hl, const T &gd, const T &gc, const value_type norm)
    requires Complex<T>
  {
    const T tmp_p = p + gc * hl;
    const T tmp_h = -gd * p + hl;
    p = norm * tmp_p;
    hl = norm * tmp_h;
  }

  template <int S, typename Derived>
  void pairing(const T gamma_, Derived &&state_)
    requires Complex<T>
  {
    static_assert(S == -1 || S == 1);
    // constexpr value_type norm = 1 / std::sqrt(2);
    const value_type norm = 1 / std::sqrt(2);
    const T gd = static_cast<T>(S) * gamma_;
    const T gc = std::conj(gd);
    Coordinates<std::size_t, 3> local(r.Ld);

    for (unsigned io = 0; io < r.Orb; ++io)
      for (unsigned i1 = NGHOSTS, I1 = r.Ld[1] - NGHOSTS; i1 < I1; ++i1) {
        local.set({std::size_t(NGHOSTS), std::size_t(i1), std::size_t(io)});
        std::size_t pair_0 = local.index;
        std::size_t pair_1 = pair_0 + offset;
        for (std::size_t i0 = 0, I0 = r.ld[0]; i0 < I0; ++i0) {
          T p = state_.coeff(pair_0);
          T hl = state_.coeff(pair_1);
          rotate_pair(p, hl, gd, gc, norm);
          state_.coeffRef(pair_0) = p;
          state_.coeffRef(pair_1) = hl;
          ++pair_0;
          ++pair_1;
        }
      }
#pragma omp barrier
  }

  template <int S>
  void pair_window(
    const unsigned io_,
    const unsigned b_,
    std::size_t (&beg)[D],
    std::size_t (&end)[D]
  ) const
  {
    static_assert(S == -1 || S == 1);
    for (unsigned k = 0; k < D; ++k) {
      beg[k] = NGHOSTS;
      end[k] = r.Ld[k] - NGHOSTS;
      int dr = h.pr.shift(k * h.pr.max_bonds + b_, io_);
      dr *= S;
      if (dr < 0 && !r.boundary[k][0])
        beg[k] += 1;
      else if (dr > 0 && !r.boundary[k][1])
        end[k] -= 1;
    }
  }

  // S = -1: partition -> lattice, S = 1: lattice -> partition
  template <int S, typename Derived>
  void nn_pairing(
    const T gamma_,
    const unsigned io_,
    const unsigned b_,
    Derived &&state_
  )
    requires Complex<T>
  {
    static_assert(S == -1 || S == 1);
    // constexpr value_type norm = 1 / std::sqrt(2);
    const value_type norm = 1 / std::sqrt(2);
    const T gd = static_cast<T>(S) * gamma_;
    const T gc = std::conj(gd);

    const std::size_t jo = h.pr.target_orb(b_, io_);
    const std::ptrdiff_t s = h.pr.dist_tile(b_, io_);

    const std::size_t ob = x.basis[2];
    const std::size_t io_base = io_ * ob;
    const std::size_t jo_base = jo * ob;

    for (std::size_t c = 0; c < r.Nd; ++c)
      pair_buf[c] = state_.coeff(jo_base + c + offset);

    // hole sector: j owned, partner i = j - s still holds its old value
    std::size_t p_beg[D];
    std::size_t p_end[D];

    pair_window<-1>(io_, b_, p_beg, p_end);
    for (std::size_t i1 = p_beg[1]; i1 < p_end[1]; ++i1) {
      const std::size_t row = i1 * std + jo_base;
      for (std::size_t i0 = p_beg[0]; i0 < p_end[0]; ++i0) {
        const std::size_t j = row + i0;
        const std::ptrdiff_t i = static_cast<std::size_t>(j - s);
        state_.coeffRef(j + offset) =
          norm * (-gd * state_.coeff(i) + pair_buf[j - jo_base]);
      }
    }
    // particle side: i owned, partner j = i + s read from the snapshot
    pair_window<1>(io_, b_, p_beg, p_end);
    for (std::size_t i1 = p_beg[1]; i1 < p_end[1]; ++i1) {
      const std::size_t row = i1 * std + io_base;
      for (std::size_t i0 = p_beg[0]; i0 < p_end[0]; ++i0) {
        const std::ptrdiff_t i = row + i0;
        const std::size_t jc = static_cast<std::size_t>(i + s) - jo_base;
        state_.coeffRef(i) = norm * (state_.coeff(i) + gc * pair_buf[jc]);
      }
    }
#pragma omp barrier
  }
};
