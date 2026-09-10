import process_cond as cond

import h5py
import numpy as np
import glob as gl
import sys

size = 512
NEnergies = 500 + 1
energy_scale = 2.7
# Integrand
# miu = np.array([-0.44, -0.4, 0.0, 0.33])
# E_grid = 1.0 * miu
# Integral
miu = np.array([-0.44, -0.27, 0.0, 0.33])
E_grid = np.linspace(-0.99, 0.99, NEnergies) * energy_scale

kbT = 0.01
pols = 512
orb = int(sys.argv[1])

eta = 8 * energy_scale / pols
sigma = 1.0 * eta
window = 2
conc = 0.5
base_dir = "Data/"

map = np.zeros((len(miu), window, window))
# for y in range(window):
for y in range(1):
    # for flag in range(4):
    for flag in range(1):
        if flag < 2:
            sign = (-1.0) ** flag
        else:
            sign = (-1.0) ** (flag + 1)
        output_file = (
            base_dir + f"Local_YZ_L{size:03d}_F{flag:02d}_O{orb:02d}_Posy{y:03d}.h5"
        )
        # for x in range(window):
        for x in range(1):
            group = f"/Calculation/CustomTwoLocal/p_{x:01d}/Gamma"
            res = cond.calculate_conductivity(
                output_file, group, miu, kbT, E_grid, eta, sigma
            )
            map[:, y, x] += res * sign
for i, fermi in enumerate(miu):
    output_path = base_dir + f"IntYZ_Clean_Window{window:03d}_Eta{eta:.3f}_Fermi{fermi:.3f}_Orb{orb:02d}.dat"
    np.savetxt(output_path, map[i, :, :], fmt="%.5e")
