#ifndef MYHDF5_H_
#define MYHDF5_H_

template <typename T>
class DataTypeFor {
public:
  static H5::DataType value;
};

template <typename T>
void get_hdf5(T *l, H5::H5File *file, char *name);

template <typename T>
void get_hdf5(T *l, H5::H5File *file, std::string &name);

void my_get_hdf5(
  std::vector<std::string> &v,
  const H5::H5File &file,
  const std::string &group
);

template <typename T>
void write_hdf5(
  const Eigen::Array<T, -1, -1> &mu,
  H5::H5File *file,
  const std::string name
);

#endif
