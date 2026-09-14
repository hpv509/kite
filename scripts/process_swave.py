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
num_itr = 256
size = lx * ly * orb * layer
ud = 0.1
mu = 0.25

base_dir = "../slurm_tools/Data/"
result = np.zeros((num_itr + 1, int(orb * layer)))
avr = np.zeros((num_itr + 1, int(orb * layer)))
var = np.zeros((num_itr + 1, int(orb * layer)))
prv = np.zeros((num_itr + 1, int(orb * layer)))

count = 0
for seed in values_seeds:
    path = base_dir + f"monolayer_exp_twist_Id{seed:03d}_Ud{ud:.2f}_Fermi{mu:.2f}.h5"
    try:
        with h5py.File(path, "r") as f:
            data = f["/Calculation/s_wave_c/Hist"][:].T
        new = data
        prv = avr.copy()
        avr += (new - avr) / (count + 1)
        var += ((new - prv) * (new - avr) - var) / (count + 1)
        count += 1
    except:
        continue
output = np.hstack((avr, np.sqrt(var / count)))
output_path = f"monolayer_exp_twist_Ud{ud:.2f}_Fermi{mu:.2f}_1.dat"
np.savetxt(output_path, output, fmt="%.14e")
