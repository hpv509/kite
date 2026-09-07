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

template <typename T, unsigned D>
herr_t Simulation<T, D>::getMembers(
  H5::Group &group,
  const std::string &name,
  void *op_data
)
{
  using Complex = std::complex<double>;
  using MatrixRowMajor =
    Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  auto *operators =
    static_cast<std::vector<Eigen::Matrix<Complex, -1, -1>> *>(op_data);
  if (!operators) {
    std::cerr << "getMembers: op_data is null\n";
    return -1;
  }
  try {
    H5::DataSet dataset = group.openDataSet(name);
    H5::DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 2) {
      std::cerr << "getMembers: dataset '" << name << "' rank != 2 (got "
                << rank << ")\n";
      return -1;
    }
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims, nullptr);
    const size_t rows = static_cast<size_t>(dims[0]);
    const size_t cols = static_cast<size_t>(dims[1]);

    H5::CompType complexType(sizeof(Complex));
    complexType.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
    complexType.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);

    std::vector<Complex> buffer;
    buffer.resize(rows * cols);
    dataset.read(buffer.data(), complexType);

    Eigen::Map<MatrixRowMajor>
      mapped(buffer.data(), static_cast<int>(rows), static_cast<int>(cols));

    operators->emplace_back(mapped);

    return 0;
  } catch (const std::exception &e) {
    return -1;
  }
}

void split_string(
  const std::string &str_,
  const char delimiter_,
  std::vector<std::string> &operators_
)
{
  std::string op;
  std::istringstream ss(str_);
  while (std::getline(ss, op, delimiter_))
    operators_.push_back(op);
}

template <typename T, unsigned D>
void Simulation<T, D>::multiply_orb_mtx(
  const Eigen::Matrix<std::complex<double>, -1, -1> &orb_mtx_,
  const KPM_Vector<T, D> *v_1_,
  KPM_Vector<T, D> *v_2_
)
{
  v_2_->v.setZero();
  for (unsigned row = 0, Row = orb_mtx_.rows(); row < Row; ++row)
    for (unsigned col = 0, Col = orb_mtx_.cols(); col < Col; ++col) {
      if constexpr (std::is_same<T, value_type>::value) {
        const value_type coeff = orb_mtx_(row, col).real();
        v_2_->v.col(0).segment(row * r.Nd, r.Nd) +=
          coeff * v_1_->v.col(0).segment(col * r.Nd, r.Nd);
      } else {
        const T coeff = static_cast<T>(orb_mtx_(row, col));
        v_2_->v.col(0).segment(row * r.Nd, r.Nd) +=
          coeff * v_1_->v.col(0).segment(col * r.Nd, r.Nd);
      }
    }
  v_2_->Exchange_Boundaries();
}

template <typename T, unsigned D>
void Simulation<T, D>::act_with_stream(
  const std::vector<std::string> &stream_,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>>
    &orb_operators_,
  const Eigen::Array<T, -1, 1> &coeffs_,
  std::array<KPM_Vector<T, D> *, 2> &vectors_,
  const unsigned write_idx_
)
{
  vectors_[1]->set_index(write_idx_);
  vectors_[1]->v.col(write_idx_).setZero();
  std::array<KPM_Vector<T, D>, 2>
    aux{KPM_Vector<T, D>(1, *this), KPM_Vector<T, D>(1, *this)};
  for (auto &vv : aux) {
    vv.v.col(0).setZero();
    vv.initiate_phases();
  }

  const unsigned read_idx = vectors_[0]->get_index();
  vectors_[0]->Exchange_Boundaries();
  for (unsigned i = 0, I = stream_.size(); i < I; ++i) {
    const std::string stream = stream_.at(i);
    std::vector<std::string> operators;
    split_string(stream, '.', operators);
    aux[0].v.col(0) = vectors_[0]->v.col(read_idx);
    for (auto op = operators.rbegin(), Op = operators.rend(); op != Op; ++op) {
      if ((*op)[0] == 'v') {
        std::vector<unsigned> velocity_components;
        std::vector<std::vector<unsigned>> collection_velocity_components;
        for (unsigned j = 1, J = op->size(); j < J; ++j)
          velocity_components.push_back(components_map[(*op)[j]]);
        collection_velocity_components.push_back(velocity_components);

        h.build_velocity(velocity_components, 0);
        aux[0].Velocity(&aux[1], collection_velocity_components, 0);
      } else if ((*op)[0] == 'r') {
        aux[1].v.col(0) = aux[0].v.col(0);
        aux[0].Position(components_map[(*op)[1]], &aux[1]);
      } else if ((*op)[0] == 'l')
        multiply_orb_mtx(
          orb_operators_[std::stoi(std::string(1, (*op)[1]))], &aux[0], &aux[1]
        );
      aux[0].v.col(0) = aux[1].v.col(0);
      aux[1].v.col(0).setZero();
    }
    vectors_[1]->v.col(write_idx_) += coeffs_(i) * aux[0].v.col(0);
  }
  vectors_[1]->Exchange_Boundaries();
}

template <typename T, unsigned D>
void Simulation<T, D>::calc_custom_one()
{
  debug_message("Entered Simulation::calc_custom_one\n");
#pragma omp barrier
  const std::string base_grp = "/Calculation/CustomOne/";
  std::string tmp = base_grp + "NumCoefficients";
#pragma omp master
  {
    H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
    Global.calculate_custom_one = false;
    try {
      int dummy_variable;
      get_hdf5<int>(&dummy_variable, file, tmp);
      Global.calculate_custom_one = true;
    } catch (H5::Exception &e) {
      debug_message("Custom Density: no need to calculate Custom One.\n");
    }
    file->close();
    delete file;
  }
#pragma omp barrier
  bool local_calculate_custom_one = false;
#pragma omp critical
  local_calculate_custom_one = Global.calculate_custom_one;
#pragma omp barrier
  if (local_calculate_custom_one) {
#pragma omp master
    std::cout << "Calculating Custom One\n";
#pragma omp barrier
    int number_moments;
    int number_samples;
    int number_vectors;
    int number_coefficients;
    Eigen::Array<T, -1, 1> coefs;
    std::vector<std::string> stream;
    std::vector<Eigen::Matrix<std::complex<double>, -1, -1>>
      operator_collection;
#pragma omp critical
    {
      H5::H5File *file = new H5::H5File(name, H5F_ACC_RDONLY);
      H5::H5File my_file = H5::H5File(name, H5F_ACC_RDONLY);
      debug_message(
        "Custom Density: checking if we need to calculate Custom One.\n"
      );
      tmp = base_grp + "NumDisorder";
      get_hdf5<int>(&number_samples, file, tmp);
      tmp = base_grp + "NumVectors";
      get_hdf5<int>(&number_vectors, file, tmp);

      tmp = base_grp + "CustomOperators/";
      H5::Group grp;
      grp = file->openGroup(tmp);
      hsize_t n = grp.getNumObjs();
      for (hsize_t i = 0; i < n; ++i) {
        H5std_string memberName = grp.getObjnameByIdx(i);
        if (
          auto err =
            this->getMembers(grp, std::string(memberName), &operator_collection)
        ) {
          std::cerr << "getMembers failed for " << memberName << "\n";
        }
      }
      tmp = base_grp + "NumMoments";
      get_hdf5<int>(&number_moments, file, tmp);
      tmp = base_grp + "NumCoefficients";
      get_hdf5<int>(&number_coefficients, file, tmp);

      tmp = base_grp + "Coefficients";
      coefs.resize(number_coefficients);
      get_hdf5<T>(coefs.data(), file, tmp);
      tmp = base_grp + "Operators";
      my_get_hdf5(stream, my_file, tmp);
      file->close();
      delete file;
    }
#pragma omp barrier
    custom_one(
      number_moments, number_samples, number_vectors, coefs,
      operator_collection, stream
    );
  }
}

template <typename T, unsigned D>
void Simulation<T, D>::custom_one(
  const int number_moments_,
  const int number_samples_,
  const int number_random_vectors_,
  const Eigen::Array<T, -1, 1> &coeffs_,
  const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operators_,
  const std::vector<std::string> &stream_
)
{
  Eigen::Matrix<T, -1, 1> gamma =
    Eigen::Matrix<T, -1, 1>::Zero(number_moments_);
  Eigen::Matrix<T, 2, 1> tmp = Eigen::Matrix<T, 2, 1>::Zero();
  KPM_Vector<T, D> kpm_trc_0(1, *this);
  KPM_Vector<T, D> kpm_trc_1(2, *this);
  std::array<KPM_Vector<T, D> *, 2> vectors{&kpm_trc_0, &kpm_trc_1};
  unsigned average = 0;
  for (int rand_v = 0; rand_v < number_random_vectors_; ++rand_v) {
    h.generate_twists();
    kpm_trc_0.initiate_vector();
    kpm_trc_0.Exchange_Boundaries();
    for (int disorder = 0; disorder < number_samples_; ++disorder) {
      h.generate_disorder();
      for (unsigned i = 0, I = vectors.size(); i < I; ++i)
        vectors.at(i)->initiate_phases();
      kpm_trc_1.set_index(0);

      act_with_stream(stream_, operators_, coeffs_, vectors, 0);
      kpm_trc_0.empty_ghosts(0);
      for (int m = 0; m < number_moments_; m += 2) {
        kpm_trc_1.cheb_iteration(m);
        kpm_trc_1.cheb_iteration(m + 1);
        tmp.setZero();
        for (std::size_t ii = 0; ii < r.Sized; ii += r.Ld[0])
          tmp += kpm_trc_0.v.block(ii, 0, r.Ld[0], 1).adjoint() *
                 kpm_trc_1.v.block(ii, 0, r.Ld[0], 2);
        gamma.segment(m, 2) += (tmp - gamma.segment(m, 2)) / (average + 1);
      }
      ++average;
    }
  }
  store_custom_one(gamma, average);
}

// template <typename T, unsigned D>
// void Simulation<T,D>::custom_one(
//   const int number_moments_,
//   const int number_samples_,
//   const int number_random_vectors_,
//   const Eigen::Array<value_type, -1, 1> &coeffs_,
//   const std::vector<Eigen::Matrix<std::complex<double>, -1, -1>> &operators_,
//   const std::vector<std::string> &stream_)
// {
//   Eigen::Matrix<T, -1, 1> gamma = Eigen::Matrix<T, -1,
//   1>::Zero(number_moments_); Eigen::Matrix<T,  2, 1> tmp   = Eigen::Matrix<T,
//   2, 1>::Zero(); KPM_Vector<T, D> kpm_trc_0(1, *this); KPM_Vector<T, D>
//   kpm_trc_1(2, *this); KPM_Vector<T, D> kpm_aux_0(1, *this); KPM_Vector<T, D>
//   kpm_aux_1(1, *this); std::vector<KPM_Vector<T, D>*> vectors;
//   vectors.push_back(&kpm_trc_0);
//   vectors.push_back(&kpm_aux_0);
//   vectors.push_back(&kpm_aux_1);
//   vectors.push_back(&kpm_trc_1);
//   h.generate_disorder();
//   h.generate_twists();
//   for (unsigned i = 0, I = vectors.size(); i < I; ++i)
//     vectors.at(i) -> initiate_phases();
//   for (std::size_t pos = 0, Pos = r.Sizet; pos < Pos; ++pos) {
//     kpm_trc_0.build_site(pos);
//     kpm_trc_0.Exchange_Boundaries();
//     kpm_trc_1.set_index(0);
//     act_with_stream(stream_, operators_, coeffs_, vectors, 0);
//     kpm_trc_0.empty_ghosts(0);
//     for (int m = 0; m < number_moments_; m += 2) {
//       kpm_trc_1.cheb_iteration(m    );
//       kpm_trc_1.cheb_iteration(m + 1);
//       tmp.setZero();
//       for (std::size_t ii = 0; ii < r.Sized; ii += r.Ld[0])
//         tmp += kpm_trc_0.v.block(ii, 0, r.Ld[0], 1).adjoint() *
//                kpm_trc_1.v.block(ii, 0, r.Ld[0], 2);
//       gamma.segment(m, 2) += tmp;
//     }
//   }
//   store_custom_one(gamma, "/Calculation/CustomDensity/Gamma");
// }

template <typename T, unsigned D>
void Simulation<T, D>::store_custom_one(
  const Eigen::Matrix<T, -1, 1> &gamma_,
  const unsigned avr_
)
{
#pragma omp master
  Global.general_gamma = Eigen::Array<T, -1, 1>::Zero(gamma_.rows());
#pragma omp barrier
#pragma omp critical
  Global.general_gamma.matrix() +=
    0.5 * (gamma_.matrix() + gamma_.matrix().adjoint());
#pragma omp barrier
#pragma omp master
  {
    const std::array<std::string, 2>
      names{"/Calculation/CustomOne/Gamma", "/Calculation/CustomOne/Avr"};
    Eigen::Array<value_type, -1, -1> avr_arr;
    avr_arr.resize(1, 1);
    avr_arr(0) = avr_;
    unsigned count = 0;
    for (const auto &n : names) {
      H5::H5File file(name, H5F_ACC_RDWR);
      if (H5Lexists(file.getId(), n.c_str(), H5P_DEFAULT))
        file.unlink(n);
      if (count == 0)
        write_hdf5(Global.general_gamma, &file, n);
      else
        write_hdf5(avr_arr, &file, n);
      ++count;
    }
  }
#pragma omp barrier
}

#define instantiate(type, dim) template class Simulation<type, dim>;
#include "instantiate.hpp"
