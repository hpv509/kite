import sys
sys.path.append("../../interfaces/")

import numpy as np

import lattice
import kite

t = 1.0
U = -3.0
mu = 0.5
beta = 20.0

Lx = 16
Ly = 16

num_random = 8
delta_initial = 0.5

filename = "ochoa_swave.h5"


def graphene_lattice():
    sqrt3 = np.sqrt(3.0)

    lat = lattice.Lattice(
        a1=[sqrt3, 0.0],
        a2=[sqrt3 / 2.0, 1.5],
    )

    lat.add_sublattices(
        ("A", [0.0, 0.0], 0.0),
        ("B", [0.0, 1.0], 0.0),
    )

    # Three A -> B nearest-neighbour bonds.
    # kite.config_system() adds the Hermitian conjugates.
    lat.add_hoppings(
        ([0, 0],  "A", "B", -t),
        ([0, -1], "A", "B", -t),
        ([1, -1], "A", "B", -t),
    )

    return lat


lat = graphene_lattice()

configuration = kite.Configuration(
    divisions=[2, 2],
    length=[Lx, Ly],
    boundaries=["periodic", "periodic"],
    is_complex=True,
    precision=1,
    spectrum_range=[-4.1, 4.1],
)

calculation = kite.Calculation(configuration)

calculation.s_wave(
    num_random=num_random,
    beta=beta,
    chemical_potential=mu,
    u=U,
    gamma=0.0,
    delta=delta_initial,
)

kite.config_system(
    lat,
    configuration,
    calculation,
    filename=filename,
)
