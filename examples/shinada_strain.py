import sys

sys.path.append("..")
from interfaces import kite
from interfaces import custom
from interfaces import lattice as latt
import numpy as np
import myaux

import warnings
from h5py.h5py_warnings import H5pyDeprecationWarning

warnings.filterwarnings("ignore", category=H5pyDeprecationWarning)

flag = int(sys.argv[1])
posy = int(sys.argv[2])
orbt = int(sys.argv[3])


def cuprate(t, r, t_prime):
    lat = latt.Lattice(a1=[1, 0], a2=[0, 1])
    lat.add_sublattices(
        ("d", [0.0, 0.0], 0.0),
        ("px", [0.5, 0.0], 0.0),
        ("py", [0.0, 0.5], 0.0),
    )
    hop_plus = (t + 1j * r) / 2
    hop_minus = (-t + 1j * r) / 2
    lat.add_hoppings(
        ([0, 0], "d", "px", hop_plus),
        ([-1, 0], "d", "px", hop_minus),
        ([0, 0], "d", "py", hop_plus),
        ([0, -1], "d", "py", hop_minus),
    )
    lat.add_hoppings(
        ([1, 0], "px", "py", -t_prime / 4),
        ([1, -1], "px", "py", t_prime / 4),
        ([0, 0], "px", "py", t_prime / 4),
        ([0, -1], "px", "py", -t_prime / 4),
    )
    return lat


def main():
    t = 1.0
    r = 1.5 * t
    t_nnn = 0.5 * t
    lattice = cuprate(t, r, t_nnn)
    nx = 2
    ny = 2
    size = 512
    lx = ly = size
    energy_scale = 2.7 * t

    node0 = [[+0, +0], "px"]
    node1 = [[+1, +0], "d"]
    node2 = [[+0, +0], "d"]

    struc_disorder_A = kite.StructuralDisorder(lattice, concentration=Nc / 100)
    struc_disorder_A.add_structural_disorder(
        (*node0, *node1, -0.5), (*node0, *node2, -0.5)
    )
    disorder_structural = [struc_disorder_A]
    mode = "open"
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-energy_scale, energy_scale],
        seed_h=2,
    )
    calculation = kite.Calculation(configuration)
    pols = 512
    if flag == 0:
        A = custom.Vertex(pols, [[1.0j, "vy.rx"], [1.0j, "rx.vy"]])
        B = custom.Vertex(pols, [[1.0j, "vy"]])
    elif flag == 1:
        A = custom.Vertex(pols, [[1.0j, "vy"]])
        B = custom.Vertex(pols, [[1.0j, "vy.rx"], [1.0j, "rx.vy"]])
    elif flag == 2:
        A = custom.Vertex(pols, [[1.0j, "vy.ry"], [1.0j, "ry.vy"]])
        B = custom.Vertex(pols, [[1.0j, "vx"]])
    else:
        A = custom.Vertex(pols, [[1.0j, "vy"]])
        B = custom.Vertex(pols, [[1.0j, "vx.ry"], [1.0j, "ry.vx"]])

    positions = []
    window = 32
    y = (ly - window) // 2 + posy
    for x in range((lx - window) // 2, (lx + window) // 2):
        positions.append([x, y, orbt])
    calculation.custom_two_local(
        stream_=[A, B],
        positions_=positions,
    )
    output_file = f"Data/Local_L{size:03d}_F{flag:02d}_O{orbt:02d}_Posy{posy:03d}.h5"
    kite.config_system(
        lattice,
        configuration,
        calculation,
        filename=output_file,
        disorder_structural=disorder_structural,
    )
    return output_file


if __name__ == "__main__":
    output = main()
    print(output, file=sys.stderr)
