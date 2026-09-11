__all__ = ["main"]

import sys

sys.path.append("..")

from interfaces import kite
from interfaces import lattice as latt
import numpy as np

t = 1.0
U = -2.0
mu = 0.0
beta = 32.0
delta_initial = [0.042, 0.042]

id = int(sys.argv[1])

def graphene_lattice(N=1, onsite=None, t=1.0, t_perp=0.0):
    a = 0.24595  # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    c0 = 0.335  # [nm] interlayer spacing

    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)
    # onsite = [-1.0 * t, 0.0 * t, -0.2 * t]
    onsite = [0.0 * t]

    v0 = np.array([0.0, a_cc])

    for l in range(N):
        shift = l * v0
        posA = np.array([0.0, -a_cc / 2 + shift[1], l * c0])
        posB = np.array([0.0, a_cc / 2 + shift[1], l * c0])
        lat.add_sublattices(
            (f"A{l+1}", posA, onsite[l]),
            (f"B{l+1}", posB, onsite[l]),
        )
    for l in range(1, N + 1):
        A, B = f"A{l}", f"B{l}"
        lat.add_hoppings(
            ([0, 0], A, B, -t),
            ([1, -1], A, B, -t),
            ([0, -1], A, B, -t),
        )
    # for l in range(1, N):
    #     lat.add_hoppings(
    #         ([0, 0], f"B{l}", f"A{l+1}", 0.0),
    #     )
    return lat


def main(N=1, onsite=None, t=1.0, t_perp=0.0):
    lattice = graphene_lattice(N, onsite, t, t_perp)

    nx = ny = 2
    lx = 256
    ly = 256
    mode = "random"
    range = 4.6
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-range, range],
        seed_h=id,
        seed_v=id,
    )
    pairing = kite.Pairing(lattice, hubbard=U)
    pairing.add_onsite_pairing("A1", delta=delta_initial[0])
    pairing.add_onsite_pairing("B1", delta=delta_initial[1])

    bdg = kite.BdG(chemical_potential=mu, beta=beta, hartree=[0.0, 0.0])
    calculation = kite.Calculation(configuration)

    # calculation.dos(
    #     num_points=4000,
    #     num_moments=4096,
    #     num_random=64,
    #     num_disorder=1
    # )
    calculation.s_wave_clean(
            num_random=64,
            num_iterations=8,
        )
    output_file = f"Data/monolayer_twist_Id{id:03d}.h5"
    kite.config_system(
        lattice, configuration, calculation,
        pairing=pairing, bdg=bdg, filename=output_file,
    )
    return output_file

if __name__ == "__main__":
    output = main(N=1, t_perp=0.0)
    print(output, file=sys.stderr)
