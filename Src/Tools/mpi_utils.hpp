#ifndef KITE_MPI_UTILS_HPP
#define KITE_MPI_UTILS_HPP
// #include <climits>
// #include <complex>
// #include <concepts>
// #include <cstddef>
// #include <cstdio>
// #include <cstdlib>
// #include <type_traits>
#ifdef USE_MPI
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
#endif

namespace kmpi {
template <class T>
struct complex_traits : std::false_type {};
template <class R>
struct complex_traits<std::complex<R>> : std::true_type {
  using real = R;
};

template <class T>
concept Real = std::same_as<T, float> || std::same_as<T, double> ||
               std::same_as<T, long double>;

template <class T>
concept Complex =
  complex_traits<T>::value && Real<typename complex_traits<T>::real>;

template <class T>
concept Scalar = Real<T> || Complex<T>;

template <class C>
concept Contiguous = requires(C &c) {
  { c.data() } -> std::convertible_to<const volatile void *>;
  { c.size() } -> std::convertible_to<std::size_t>;
} && Scalar<std::remove_cvref_t<decltype(*std::declval<C &>().data())>>;

inline bool active() noexcept
{
#ifdef USE_MPI
  int init = 0, fin = 0;
  MPI_Initialized(&init);
  MPI_Finalized(&fin);
  return init && !fin;
#else
  return false;
#endif
}

inline int rank() noexcept
{
#ifdef USE_MPI
  if (active()) {
    int r = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    return r;
  }
#endif
  return 0;
}

inline int size() noexcept
{
#ifdef USE_MPI
  if (active()) {
    int s = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &s);
    return s;
  }
#endif
  return 1;
}

inline bool is_root() noexcept { return rank() == 0; }

inline void barrier() noexcept
{
#ifdef USE_MPI
  if (active())
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

[[noreturn]] inline void abort(int code = 1) noexcept
{
#ifdef USE_MPI
  if (active())
    MPI_Abort(MPI_COMM_WORLD, code);
#endif
  std::exit(code);
}

class Session {
public:
  Session(int &argc, char **&argv)
  {
#ifdef USE_MPI
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
      std::fprintf(stderr, "KITEx: MPI_THREAD_FUNNELED not supported\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
#else
    (void)argc;
    (void)argv;
#endif
  }
  ~Session()
  {
#ifdef USE_MPI
    if (active())
      MPI_Finalize();
#endif
  }
  Session(const Session &) = delete;
  Session &operator=(const Session &) = delete;
};

#ifdef USE_MPI
template <Real R>
MPI_Datatype mpi_type()
{
  if constexpr (std::same_as<R, float>)
    return MPI_FLOAT;
  else if constexpr (std::same_as<R, double>)
    return MPI_DOUBLE;
  else
    return MPI_LONG_DOUBLE;
}
#endif

template <Real R>
void sum_all(R *data, std::size_t n)
{
#ifdef USE_MPI
  if (!active() || size() == 1 || n == 0)
    return;
  if (n > static_cast<std::size_t>(INT_MAX))
    abort(2);
  const int cnt = static_cast<int>(n);
  const MPI_Datatype t = mpi_type<R>();
  if (is_root())
    MPI_Reduce(MPI_IN_PLACE, data, cnt, t, MPI_SUM, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce(data, nullptr, cnt, t, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(data, cnt, t, 0, MPI_COMM_WORLD);
#else
  (void)data;
  (void)n;
#endif
}

template <Complex C>
void sum_all(C *data, std::size_t n)
{
  using R = typename complex_traits<C>::real;
  sum_all(reinterpret_cast<R *>(data), 2 * n);
}

template <Contiguous C>
void sum_all(C &c)
{
  sum_all(c.data(), static_cast<std::size_t>(c.size()));
}

}
#endif
