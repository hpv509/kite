import sys
sys.path.append("..")
from interfaces import kite
from interfaces import lattice as latt
import numpy as np

def weyl_semimetal(t=1.0, m=2.0):
    a1 = np.array([1, 0, 0])
    a2 = np.array([0, 1, 0])
    a3 = np.array([0, 0, 1])
    lat = latt.Lattice(a1 = a1, a2 = a2, a3 = a3)
    lat.add_sublattices(
        ('A', [0, 0, 0], m * t),
        ('B', [0, 0, 0], -m * t)
    )
    lat.add_hoppings(
        ([1, 0, 0], 'A', 'A', -0.5 * t),
        ([1, 0, 0], 'B', 'B',  0.5 * t),
        ([1, 0, 0], 'A', 'B', -0.5j * t),
        ([1, 0, 0], 'B', 'A', -0.5j * t),
        ([0, 1, 0], 'A', 'A', -0.5 * t),
        ([0, 1, 0], 'B', 'B',  0.5 * t),
        ([0, 1, 0], 'A', 'B', -0.5 * t),
        ([0, 1, 0], 'B', 'A',  0.5 * t),
        ([0, 0, 1], 'A', 'A', -0.5 * t),
        ([0, 0, 1], 'B', 'B',  0.5 * t)
    )
    return lat

def main():
    lattice = weyl_semimetal()
    # add scalar on-site Disorder (Box distribution in [W/2,W/2])
    # Number of Samples controlled by num_disorder
    # disorder = kite.Disorder(lattice)
    # disorder.add_disorder('A', 'Uniform', 0.0, anderson_w/np.sqrt(12))
    # disorder.add_disorder('B', 'Uniform', 0.0, anderson_w/np.sqrt(12))
    nx = 2
    ny = nz = 2
    lx = ly = lz = 32
    mode = "random"
    configuration = kite.Configuration(
        divisions=[nx, ny, nz],
        length=[lx, ly, lz],
        boundaries=["open", mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-8.0, 8.0]
    )
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=1000,
        num_moments=256,
        num_random=1,
        num_disorder=1
    )
    filename = "Data/Weyl.h5"
    kite.config_system(lattice, configuration, calculation, filename=filename)
    return filename

if __name__ == "__main__":
    output = main()
    print(output, file=sys.stderr)
