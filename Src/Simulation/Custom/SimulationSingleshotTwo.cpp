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
#include "Coefficients.hpp"

template <typename T, unsigned D>
void Simulation<T, D>::calc_custom_ss_two()
  requires Complex<T>
{
  debug_message("Entered Simulation::calc_custom_single_two\n");
  std::string base_grp = "/Calculation/CustomSingleTwo/";
  std::string tmp = base_grp + "NumDisorder";
#pragma omp barrier
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_custom_ss_two = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, file, tmp);
      Global.calculate_custom_ss_two = true;
    } catch (H5::Exception &e) {
      debug_message(
        "CustomSingleTwo: no need to calculate Custom Singleshot Two.\n"
      );
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_custom_ss_two = false;
#pragma omp critical
  local_calculate_custom_ss_two = Global.calculate_custom_ss_two;
#pragma omp barrier
  if (local_calculate_custom_ss_two) {
#pragma omp master
    std::cout << "Calculating Custom Singleshot Two\n";
#pragma omp barrier
    unsigned number_samples;
    unsigned number_vectors;
    value_type sigma;
    value_type gamma;
    Eigen::Array<value_type, -1, 1> energies;
    std::vector<Eigen::Array<T, -1, 1>> coefs(2);
    std::vector<std::vector<std::string>> stream(2);
    std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> orb_operators;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      H5::H5File my_file = H5::H5File(name, H5F_ACC_RDONLY);
      debug_message(
        "Custom Single Two: checking if we need to calculate Custom Single "
        "Two.\n"
      );
      tmp = base_grp + "NumDisorder";
      get_hdf5<unsigned>(&number_samples, file, tmp);
      tmp = base_grp + "NumVectors";
      get_hdf5<unsigned>(&number_vectors, file, tmp);
      tmp = base_grp + "Gamma";
      get_hdf5<value_type>(&gamma, file, tmp);
      tmp = base_grp + "Sigma";
      get_hdf5<value_type>(&sigma, file, tmp);
      tmp = base_grp + "Energies";
      H5::DataSet dataset = file->openDataSet(tmp);
      H5::DataSpace dataspace = dataset.getSpace();
      hsize_t dims[1] = {0};
      dataspace.getSimpleExtentDims(dims, nullptr);
      energies.resize(dims[0]);
      get_hdf5<value_type>(energies.data(), file, tmp);

      std::string path;
      for (unsigned i = 0; i < 2; ++i) {
        const std::string base = base_grp + "Vertex" + std::to_string(i);
        path = base + "/NumCoefficients";
        int number_coefficients;
        get_hdf5<int>(&number_coefficients, file, path);
        coefs.at(i).resize(number_coefficients);
        path = base + "/Coefficients";
        get_hdf5<T>(coefs.at(i).data(), file, path);
        path = base + "/Operators";
        my_get_hdf5(stream.at(i), my_file, path);
      }
      tmp = base_grp + "CustomOperators/";
      H5::Group grp;
      grp = file->openGroup(tmp);
      hsize_t n = grp.getNumObjs();
      for (hsize_t i = 0; i < n; ++i) {
        H5std_string memberName = grp.getObjnameByIdx(i);
        if (
          auto err =
            this->getMembers(grp, std::string(memberName), &orb_operators)
        ) {
        }
      }
      file->close();
      delete file;
    }
#pragma omp barrier
    custom_ss_two(
      number_samples, number_vectors, gamma, sigma, energies, stream, coefs,
      orb_operators
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::custom_ss_two(
  const int samples_,
  const int vectors_,
  const value_type gamma_,
  const value_type sigma_,
  const Eigen::Array<value_type, -1, 1> &energies_,
  const std::vector<std::vector<std::string>> &stream_,
  const std::vector<Eigen::Array<T, -1, 1>> &coeffs_,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operators_
)
  requires Complex<T>
{
  constexpr T I{0.0, 1.0};
  value_type energy_scale;
#pragma omp critical
  {
    H5::H5File file(name, H5F_ACC_RDONLY);
    std::string path = "/EnergyScale";
    get_hdf5<value_type>(&energy_scale, &file, path);
    file.close();
  }
#pragma omp barrier
  const unsigned num_en = energies_.size();
  const value_type norm_0 = 1 / energy_scale;
  const value_type sigma = sigma_ * norm_0;
  const value_type gamma = gamma_ * norm_0;
  const value_type norm_1 = norm_0 * norm_0;

  std::array<Eigen::Array<T, -1, -1>, 2> coefs;
  unsigned count = 0;
  for (const auto &vv : energies_) {
    const value_type en = vv / energy_scale;
    const T z = en + I * gamma;
    const Eigen::Array<T, -1, 1> dgreen_tmp =
      Coefficients::build_dgreen<value_type>(z) * norm_1;
    const Eigen::Array<value_type, -1, 1> gauss_tmp =
      Coefficients::build_gaussian<value_type>(en, sigma) * norm_0;
    if (count == 0) {
      coefs[0].resize(dgreen_tmp.size(), num_en);
      coefs[1].resize(gauss_tmp.size(), num_en);
    }
    coefs[0].col(count) = dgreen_tmp;
    coefs[1].col(count) = gauss_tmp;
    ++count;
  }
  std::array<Eigen::Array<T, -1, 1>, 2> ket =
    {Eigen::Array<T, -1, 1>(r.Sized), Eigen::Array<T, -1, 1>(r.Sized)};
  Eigen::Array<T, -1, 1> bra(r.Sized);
  Eigen::Array<T, -1, 1> psi(r.Sized);
  Eigen::Array<T, -1, 1> result(num_en);
  result.setZero();
  std::array<KPM_Vector<T, D>, 2> vecs =
    {KPM_Vector<T, D>(1, *this), KPM_Vector<T, D>(2, *this)};
  std::array<KPM_Vector<T, D> *, 2> ptrs = {&vecs[0], &vecs[1]};

  unsigned average = 0;
  for (int disorder = 0; disorder < samples_; ++disorder) {
    h.generate_disorder();
    for (int vec = 0; vec < vectors_; ++vec) {
      h.generate_twists();
      for (auto &vec : vecs)
        vec.initiate_phases();
      vecs[0].initiate_vector();
      vecs[0].empty_ghosts(0);
      bra = vecs[0].v.col(0).conjugate();
      vecs[1].set_index(0);

      act_with_stream(stream_[0], operators_, coeffs_[0], ptrs, 0);
      psi = vecs[1].v.col(0);
      const value_type weight = 1.0 / (average + 1);

      for (unsigned i = 0; i < num_en; ++i) {
        for (unsigned p = 0; p < 2; ++p) {
          vecs[1].set_index(0);
          vecs[1].v.col(0) = psi;
          ket[p].setZero();
          const unsigned id_0 = p % 2;
          const unsigned id_1 = (p + 1) % 2;
          for (unsigned n = 0, N = coefs[id_0].rows(); n < N; ++n) {
            vecs[1].cheb_iteration(n);
            const unsigned idx_read = vecs[1].get_index();
            const T coef = coefs[id_0](n, i);
            ket[p] += coef * vecs[1].v.col(idx_read).array();
          }
          vecs[0].v.col(0) = ket[p];
          act_with_stream(stream_[1], operators_, coeffs_[1], ptrs, 0);
          ket[p].setZero();
          for (unsigned n = 0, N = coefs[id_1].rows(); n < N; ++n) {
            vecs[1].cheb_iteration(n);
            const unsigned idx_read = vecs[1].get_index();
            const T coef = std::conj(coefs[id_1](n, i));
            ket[p] += coef * vecs[1].v.col(idx_read).array();
          }
        }
        T res = (bra * ket[0] + (bra * ket[1]).conjugate()).sum();
        res *= 0.5;
        result(i) += (res - result(i)) * weight;
      }
      ++average;
    }
  }
  store_custom_ss_two(result, average);
}

template <typename T, unsigned D>
void Simulation<T, D>::store_custom_ss_two(
  const Eigen::Array<T, -1, -1> &result_,
  const unsigned avr_
)
  requires Complex<T>
{
  debug_message("Entered store_gamma\n");
#pragma omp master
  {
    Global.singleshot_cond.resize(result_.rows(), result_.cols());
    Global.singleshot_cond.setZero();
  }
#pragma omp barrier
#pragma omp critical
  Global.singleshot_cond += result_;
#pragma omp barrier
#pragma omp master
  {
    const std::array<std::string, 2> names{
      "/Calculation/CustomSingleTwo/Result", "/Calculation/CustomSingleTwo/Avr"
    };
    Eigen::Array<value_type, -1, -1> avr_arr;
    avr_arr.resize(1, 1);
    avr_arr(0) = avr_;
    unsigned count = 0;
    for (const auto &n : names) {
      H5::H5File file(name, H5F_ACC_RDWR);
      if (H5Lexists(file.getId(), n.c_str(), H5P_DEFAULT))
        file.unlink(n);
      if (count == 0)
        write_hdf5(Global.singleshot_cond, &file, n);
      else
        write_hdf5(avr_arr, &file, n);
      ++count;
    }
  }
#pragma omp barrier
  debug_message("Left store_gamma\n");
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
