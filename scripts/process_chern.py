import h5py as h5
import numpy as np
import glob as gl
import shutil as sh

base_dir = r"../slurm_tools/Data/"
save_dir = r""
group_name = "/Calculation/LCM/Marker"

length = 512
beta = 256
energy = 0.0
values_disorder = np.array([6.8, 6.9, 7.0, 7.1, 7.2, 7.25, 7.3, 7.4, 7.5, 7.6])

for i, w in enumerate(values_disorder):
    path = base_dir + f"LCM_L{length:04d}_Beta{beta:04d}*W{w:.1f}*"
    save = (
        save_dir + f"Single_Beta{beta:1.1f}_L{length:04d}_EF{energy:1.2f}_W{w:1.2f}.dat"
    )
    file_list = gl.glob(path)
    data = np.array([])
    for f in file_list:
        try:
            with h5.File(f, "r") as file:
                tmp = file[group_name][:]
                data = np.concatenate((data, tmp[0, :].imag))
        except:
            continue
        avg = 0.0
        var = 0.0
        prv = 0.0
        summed = data[::2] + data[1::2]
        for i, new in enumerate(summed):
            prv = avg
            avg += (new - avg) / (i + 1)
            var += ((new - avg) * (new - prv) - var) / (i + 1)
        with open(save, "w") as f:
            f.write(f"{avg:1.15e} {np.sqrt(var / len(summed)):1.15e}\n")
            for x in summed:
                f.write(f"{x:1.15e}\n")
