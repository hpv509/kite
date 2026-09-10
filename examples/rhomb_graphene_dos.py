__all__ = ["main"]

import sys

sys.path.append("..")

from interfaces import kite
from interfaces import lattice as latt
import numpy as np


def graphene_lattice(N=3, onsite=None, t=1.0, t_perp=0.12):
    a = 0.24595  # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    c0 = 0.335  # [nm] interlayer spacing

    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)
    if onsite is None:
        onsite = [(0, 0)] * N

    v0 = np.array([0.0, a_cc])
    for l in range(N):
        shift = l * v0
        posA = np.array([0.0, -a_cc / 2 + shift[1], l * c0])
        posB = np.array([0.0, a_cc / 2 + shift[1], l * c0])
        lat.add_sublattices(
            (f"A{l+1}", posA, onsite[l][0]),
            (f"B{l+1}", posB, onsite[l][1]),
        )
    for l in range(1, N + 1):
        A, B = f"A{l}", f"B{l}"
        lat.add_hoppings(
            ([0, 0], A, B, -t),
            ([1, -1], A, B, -t),
            ([0, -1], A, B, -t),
        )
    for l in range(1, N):
        lat.add_hoppings(
            ([0, 0], f"B{l}", f"A{l+1}", t_perp),
        )
    return lat


def main(N=3, onsite=None, t=1.0, t_perp=0.12):
    lattice = graphene_lattice(N, onsite, t, t_perp)

    nx = ny = 2
    lx = ly = 512
    mode = "random"
    range = 3.1 * t + t_perp
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=False,
        precision=1,
        spectrum_range=[-range, range],
    )
    calculation = kite.Calculation(configuration)
    calculation.dos(
        num_points=4000,
        num_moments=2048,
        num_random=32,
        num_disorder=1,
    )
    output_file = f"Data/rhombo_graphene.h5"
    kite.config_system(lattice, configuration, calculation, filename=output_file)
    return output_file


if __name__ == "__main__":
    output = main(N=4, t_perp=0.6)
    print(output, file=sys.stderr)
