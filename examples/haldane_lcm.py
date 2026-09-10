import sys

sys.path.append("..")

import sys
from interfaces import kite
from interfaces import custom
from interfaces import lattice as latt
import numpy as np
import myaux

seed = int(sys.argv[1])
ll = int(sys.argv[2])
beta = int(sys.argv[3])
WW = float(sys.argv[4])

def haldane(onsite=(0, 0), t=1):
    a = 0.24595
    a_cc = 0.142
    t2 = t / 3
    phase = np.exp(1j * np.pi / 2)
    a1 = 1.0 * np.array([1, 0])
    a2 = 0.5 * np.array([1, np.sqrt(3)])
    a1 *= a
    a2 *= a
    lat = latt.Lattice(a1=a1, a2=a2)
    lat.add_sublattices(
        ("A", [0, -a_cc / 2], onsite[0]), ("B", [0, a_cc / 2], onsite[1])
    )
    lat.add_hoppings(
        ([0, 0], "A", "B", -t),
        ([1, -1], "A", "B", -t),
        ([0, -1], "A", "B", -t),
        ([1, 0], "A", "A", -t2 * phase),
        ([0, -1], "A", "A", -t2 * phase),
        ([-1, 1], "A", "A", -t2 * phase),
        ([1, 0], "B", "B", -t2 * np.conj(phase)),
        ([0, -1], "B", "B", -t2 * np.conj(phase)),
        ([-1, 1], "B", "B", -t2 * np.conj(phase)),
    )
    return lat


def main(onsite=(0, 0), t=1):
    lattice = haldane(onsite, t)
    nx = ny = 2
    lx = ly = ll
    std = WW / np.sqrt(12)
    energy_scale = (4.1 + 0.5 * WW) * t

    disorder = kite.Disorder(lattice)
    disorder.add_disorder("A", "Uniform", 0.0, std)
    disorder.add_disorder("B", "Uniform", 0.0, std)

    mode = "open"
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=0,
        spectrum_range=[-energy_scale, energy_scale],
    )
    calculation = kite.Calculation(configuration)
    calculation.local_chern(
        num_disorder_=8192, beta_=beta, miu_=0.0, pos_=np.array([0, 0])
    )
    output_file = f"Data/LCM_L{ll:04d}_Beta{beta:04d}_W{WW:2.1f}_Seed{seed:03d}.h5"
    kite.config_system(
        lattice, configuration, calculation, filename=output_file, disorder=disorder
    )
    return output_file


if __name__ == "__main__":
    output = main()
    print(output, file=sys.stderr)
