#ifndef LOOP_H
#define LOOP_H

template <unsigned D, unsigned I = 0>
struct UnitCellLoop {
  template <typename F>
  static void run(
    std::array<unsigned, D> &idx_,
    const std::array<unsigned, D> &start_,
    const std::array<unsigned, D> &final_,
    F &&f_
  )
  {
    for (idx_[I] = start_[I]; idx_[I] < final_[I]; ++idx_[I]) {
      UnitCellLoop<D, I + 1>::run(idx_, start_, final_, std::forward<F>(f_));
    }
  }
};

template <unsigned D>
struct UnitCellLoop<D, D> {
  template <typename F>
  static void run(
    std::array<unsigned, D> &idx_,
    const std::array<unsigned, D> &start_,
    const std::array<unsigned, D> &final_,
    F &&f_
  )
  {
    std::forward<F>(f_)(idx_);
  }
};

#endif
