import process_cond as cond
import h5py
import numpy as np
import glob as gl
import sys

values_seeds = np.arange(1, 400 + 1)
lx = 256
ly = 256
orb = 2
layer = 1
num_itr = 8
size = lx * ly * orb * layer

base_dir = "../slurm_tools/Data/"
# result = np.zeros(size)
# avr = np.zeros(size)
result = np.zeros((num_itr + 1, int(orb * layer)))
avr = np.zeros((num_itr + 1, int(orb * layer)))
var = np.zeros((num_itr + 1, int(orb * layer)))
prv = np.zeros((num_itr + 1, int(orb * layer)))

count = 0
for seed in values_seeds:
    path = base_dir + f"monolayer_twist_Id{seed:03d}.h5"
    try:
        with h5py.File(path, "r") as f:
            # data = f["/Calculation/s_wave/Map"][:].flatten()
            data = f["/Calculation/s_wave_c/Hist"][:].T
        new = data
        prv = avr.copy()
        avr += (new - avr) / (count + 1)
        var += ((new - prv) * (new - avr) - var) / (count + 1)

        count += 1
    except:
        continue
output = np.hstack((avr, np.sqrt(var / count)))
output_path = f"monolayer_graphene_checkg_twists_Itr{num_itr:04d}.dat"
np.savetxt(output_path, output, fmt="%.7e")
