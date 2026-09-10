#include <H5Cpp.h>
#include <iostream>
#include "cond.hpp"

template <typename DerivedY, typename DerivedX>
auto simpson_integrate(
  const Eigen::ArrayBase<DerivedY> &y,
  const Eigen::ArrayBase<DerivedX> &x
)
{
  using ScalarType = typename DerivedY::Scalar;
  const int n = y.size();
  const type dx = x[1] - x[0];
  ScalarType sum = y[0] + y[n - 1];
  for (int i = 1; i < n - 1; i += 2)
    sum += 4.0 * y[i];
  for (int i = 2; i < n - 1; i += 2)
    sum += 2.0 * y[i];
  return sum * (dx / 3.0);
}

arr<type, -1, 1>
fermi_function(const arr<type, -1, 1> &E, const type mu, const type beta)
{
  arr<type, -1, 1> arg = beta * (E - mu);
  arg = arg.cwiseMax(-100.0).cwiseMin(100.0);
  return 1.0 / (1.0 + arg.exp());
}

mtx<type, -1, 1> calculate_conductivity(
  const std::string_view file_path,
  const std::string_view grp_path,
  const std::span<const type> mu_values,
  const type k_BT,
  const std::span<const type> E_grid,
  const type eta,
  const type sigma
)
{
  constexpr int pol_g = 512;
  constexpr int pol_d = 512;

  int num_orbitals = 0;
  type energy_scale = 0.0;
  mtx<type, -1, -1> latt_vecs;
  mtx<cplx, -1, -1> moments_matrix_full;

  Eigen::Map<const mtx<type, -1, 1>> mu_vec(mu_values.data(), mu_values.size());
  Eigen::Map<const mtx<type, -1, 1>> E_vec(E_grid.data(), E_grid.size());
  try {
    H5::Exception::dontPrint();
    H5::H5File file(std::string(file_path), H5F_ACC_RDONLY);
    H5::DataSet orb_ds = file.openDataSet("NOrbitals");
    orb_ds.read(&num_orbitals, H5::PredType::NATIVE_INT);
    H5::DataSet scale_ds = file.openDataSet("EnergyScale");
    scale_ds.read(&energy_scale, H5::PredType::NATIVE_DOUBLE);

    H5::DataSet latt_ds = file.openDataSet("LattVectors");
    H5::DataSpace latt_space = latt_ds.getSpace();
    hsize_t latt_dims[2];
    latt_space.getSimpleExtentDims(latt_dims, nullptr);

    std::vector<type> latt_data(latt_dims[0] * latt_dims[1]);
    latt_ds.read(latt_data.data(), H5::PredType::NATIVE_DOUBLE);
    latt_vecs = Eigen::Map<Eigen::Matrix<
      type, -1, -1,
      Eigen::RowMajor>>(latt_data.data(), latt_dims[0], latt_dims[1]);

    H5::DataSet mom_ds = file.openDataSet(std::string(grp_path));
    H5::DataSpace mom_space = mom_ds.getSpace();
    hsize_t mom_dims[2];
    mom_space.getSimpleExtentDims(mom_dims, nullptr);
    H5::CompType complex_datatype(sizeof(cplx));
    complex_datatype.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
    complex_datatype
      .insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);
    std::vector<cplx> mom_data(mom_dims[0] * mom_dims[1]);
    mom_ds.read(mom_data.data(), complex_datatype);

    mtx<cplx, -1, -1> temp_moments = Eigen::Map<Eigen::Matrix<
      cplx, -1, -1,
      Eigen::RowMajor>>(mom_data.data(), mom_dims[0], mom_dims[1]);
    moments_matrix_full = temp_moments.transpose();
  } catch (const H5::Exception &err) {
    std::cerr << "HDF5 Exception: " << err.getDetailMsg() << '\n';
    return mtx<type, -1, 1>::Zero(mu_values.size());
  }
  const type unit_cell_area = std::abs(
    latt_vecs(0, 0) * latt_vecs(1, 1) - latt_vecs(0, 1) * latt_vecs(1, 0)
  );
  constexpr type spin_degeneracy = 1.0;
  const mtx<cplx, -1, -1> moments_matrix =
    moments_matrix_full.topRows(pol_d).leftCols(pol_g);
  Coefficients coefs;
  const mtx<type, -1, -1> delta =
    coefs.build_gaussian(E_vec, sigma, energy_scale, pol_d);
  const mtx<cplx, -1, 1> E_grid_complex = E_vec.array() + cst::I * eta;
  const mtx<cplx, -1, -1> dgreenR =
    coefs.build_dgreen(E_grid_complex, energy_scale, pol_g);
  const mtx<cplx, -1, -1> tmp = moments_matrix * dgreenR;
  const arr<cplx, -1, 1> GammaE =
    (delta.cast<cplx>().array() * tmp.array()).colwise().sum();
  mtx<cplx, -1, 1> cond_dc = mtx<cplx, -1, 1>::Zero(mu_vec.size());
  const type beta = energy_scale / k_BT;
  const arr<type, -1, 1> E_grid_scaled = E_vec.array() / energy_scale;

  for (int i = 0; i < mu_vec.size(); ++i) {
    type mu_scaled = mu_vec[i] / energy_scale;
    const arr<cplx, -1, 1> integrand =
      GammaE * fermi_function(E_grid_scaled, mu_scaled, beta).cast<cplx>();
    cond_dc[i] = simpson_integrate(integrand, E_grid_scaled);
  }
  constexpr type units = 1.0 / (2.0 * std::numbers::pi);
  const type density_scale =
    (num_orbitals * spin_degeneracy) / (unit_cell_area * units);
  return 2.0 * cond_dc.imag() * density_scale;
}
