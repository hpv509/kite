#include <fstream>
#include <iomanip>
#include <iostream>
#include <format>
#include "cond.hpp"

int main(int argc, char *argv[])
{
  const int orb = std::stoi(argv[1]);
  constexpr int size = 512;
  constexpr int NEnergies = 500 + 1;
  constexpr double energy_scale = 2.7;
  std::vector<type> miu = {-0.44, -0.27, 0.0, 0.33};
  std::vector<type> E_grid(NEnergies);
  for (int i = 0; i < NEnergies; ++i) {
    E_grid[i] = (-0.99 + (2.0 * 0.99 * i) / (NEnergies - 1)) * energy_scale;
  }
  constexpr double kbT = 0.01;
  constexpr int pols = 512;

  const double eta = 8.0 * energy_scale / pols;
  const double sigma = 1.0 * eta;
  const double conc = 0.25;
  constexpr int window = 80;
  std::string base_dir = "Data/";

  std::vector<mtx<type, -1, -1>>
    map_matrices(miu.size(), mtx<type, -1, -1>::Zero(window, window));

  for (int y = 0; y < window; ++y) {
    for (int flag = 2; flag < 4; ++flag) {
    // for (int flag = 0; flag < 2; ++flag) {
      double sign = 1.0;
      if (flag < 2)
        sign = std::pow(-1.0, flag);
      else
        sign = std::pow(-1.0, flag + 1);

      std::string output_file = std::format(
        "{}LocalP_L{:03d}_F{:02d}_O{:02d}_Posy{:03d}.h5", base_dir, size, flag,
        orb, y
      );
      for (int x = 0; x < window; ++x) {
        std::string group =
          std::format("/Calculation/CustomTwoLocal/p_{:d}/Gamma", x);
        mtx<type, -1, 1> res = calculate_conductivity(
          output_file, group, miu, kbT, E_grid, eta, sigma
        );
        for (size_t i = 0; i < miu.size(); ++i)
          map_matrices[i](y, x) += res(i) * sign;
      }
    }
  }
  for (size_t i = 0; i < miu.size(); ++i) {
    std::string output_path = std::format(
      "SigmaXYY_pVacCentral_Window{:03d}_Eta{:.3f}_Fermi{:.3f}_Orb{:"
      "02d}."
      "dat",
      window, eta, miu[i], orb
    );
    std::ofstream out_file(output_path);
    out_file << std::scientific << std::setprecision(5);
    for (int y = 0; y < window; ++y) {
      for (int x = 0; x < window; ++x)
        out_file << map_matrices[i](y, x) << (x == window - 1 ? "" : " ");
      out_file << '\n';
    }
  }
  return 0;
}
