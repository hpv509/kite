/***********************************************************/
/*                                                         */
/*   Copyright (C) 2018-2022, M. Andelkovic, L. Covaci,    */
/*  A. Ferreira, S. M. Joao, J. V. Lopes, T. G. Rappoport  */
/*                                                         */
/***********************************************************/

#include "Generic.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

template <>
H5::DataType DataTypeFor<int>::value = H5::PredType::NATIVE_INT;
template <>
H5::DataType DataTypeFor<unsigned int>::value = H5::PredType::NATIVE_UINT;
template <>
H5::DataType DataTypeFor<unsigned long>::value = H5::PredType::NATIVE_ULONG;
template <>
H5::DataType DataTypeFor<float>::value = H5::PredType::NATIVE_FLOAT;
template <>
H5::DataType DataTypeFor<double>::value = H5::PredType::NATIVE_DOUBLE;
template <>
H5::DataType DataTypeFor<long double>::value = H5::PredType::NATIVE_LDOUBLE;

template <typename T>
void get_hdf5(T *l, H5::H5File *file, char *name)
{
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));

  if constexpr (Complex<T>) {
    H5::CompType complex_data_type(sizeof(T));
    using value_type = typename T::value_type;

    complex_data_type.insertMember("r", 0, DataTypeFor<value_type>::value);
    complex_data_type
      .insertMember("i", sizeof(value_type), DataTypeFor<value_type>::value);
    dataset.read(l, complex_data_type);
  } else
    dataset.read(l, DataTypeFor<T>::value);
}

template <typename T>
void get_hdf5(T *l, H5::H5File *file, std::string &name)
{
  H5::DataSet dataset = H5::DataSet(file->openDataSet(name));

  if constexpr (Complex<T>) {
    H5::CompType complex_data_type(sizeof(T));
    using value_type = typename T::value_type;

    complex_data_type.insertMember("r", 0, DataTypeFor<value_type>::value);
    complex_data_type
      .insertMember("i", sizeof(value_type), DataTypeFor<value_type>::value);
    dataset.read(l, complex_data_type);
  } else
    dataset.read(l, DataTypeFor<T>::value);
}

void my_get_hdf5(
  std::vector<std::string> &v_,
  const H5::H5File &file_,
  const std::string &group_
)
{
  const H5::DataSet dataset = file_.openDataSet(group_);
  const H5::DataSpace dataspace = dataset.getSpace();
  const H5::StrType str_type = dataset.getStrType();

  hsize_t dims[1];
  dataspace.getSimpleExtentDims(dims, nullptr);
  const size_t num_elements = dims[0];
  v_.resize(num_elements);

  std::vector<char *> stream_data(num_elements, nullptr);
  dataset.read(stream_data.data(), str_type);
  for (size_t i = 0; i < num_elements; ++i) {
    v_[i] = stream_data[i];
    free(stream_data[i]);
  }
}

template <typename T>
void write_hdf5(
  const Eigen::Array<T, -1, -1> &mu,
  H5::H5File *file,
  const std::string name
)
{
  hsize_t dims[2], chunk_dims[2]; // dataset dimensions
  dims[0] = chunk_dims[0] = mu.cols();
  dims[1] = chunk_dims[1] = mu.rows();
  H5::DataSet dataset;
  H5::DataSpace dataspace = H5::DataSpace(2, dims);
  H5::DSetCreatPropList plist;

  plist.setChunk(2, chunk_dims);
  plist.setDeflate(6);

  if constexpr (Complex<T>) {
    using value_type = typename T::value_type;
    H5::CompType complex_datatype(sizeof(T));
    complex_datatype.insertMember("r", 0, DataTypeFor<value_type>::value);
    complex_datatype
      .insertMember("i", sizeof(value_type), DataTypeFor<value_type>::value);

    try {
      H5::Exception::dontPrint();
      dataset = file->createDataSet(name, complex_datatype, dataspace);
    } catch (H5::FileIException &E) {
      dataset = file->openDataSet(name);
    }
    dataset.write(mu.data(), complex_datatype);

  } else {
    try {
      H5::Exception::dontPrint();
      dataset = file->createDataSet(name, DataTypeFor<T>::value, dataspace);
    } catch (H5::FileIException &E) {
      dataset = file->openDataSet(name);
    }
    dataset.write(mu.data(), DataTypeFor<T>::value);
  }
}

#define instantiateTYPE(type)                                                  \
  template void get_hdf5<type>(type *, H5::H5File *, char *);                  \
  template void get_hdf5<type>(type *, H5::H5File *, std::string &);           \
  template void write_hdf5(                                                    \
    const Eigen::Array<type, -1, -1> &, H5::H5File *, const std::string        \
  );

instantiateTYPE(float);
instantiateTYPE(double);
instantiateTYPE(long double);
instantiateTYPE(std::complex<float>);
instantiateTYPE(std::complex<double>);
instantiateTYPE(std::complex<long double>);
instantiateTYPE(int);
instantiateTYPE(unsigned);
instantiateTYPE(unsigned long);
