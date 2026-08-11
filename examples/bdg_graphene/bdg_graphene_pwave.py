import sys
sys.path.append("../../interfaces/")

import numpy as np

import lattice
import kite


t = 1.0
U = -2.0
V = -4.0

mu = 0.4
beta = 100.0

gamma_initial = 0.0
s_delta_initial = 0.2
nn_delta_initial = -0.3

num_random = 128

Lx = 32
Ly = 32

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
    divisions=[4, 4],
    length=[Lx, Ly],
    boundaries=["periodic", "periodic"],
    is_complex=False,  # kite.py should automatically change this
    precision=1,
    spectrum_range=[-10.0, 10.0],
)

calculation = kite.Calculation(configuration)

calculation.p_wave(
    num_random=num_random,
    beta=beta,
    chemical_potential=mu,
    u=U,
    v=V,
    gamma=gamma_initial,
    s_delta=s_delta_initial,
    nn_delta=nn_delta_initial,
)

filename = "graphene_pwave_cluster_Nr128.h5"

kite.config_system(
    lat,
    configuration,
    calculation,
    filename=filename,
)
