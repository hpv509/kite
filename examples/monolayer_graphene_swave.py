__all__ = ["main"]

import sys

sys.path.append("..")

from interfaces import kite
from interfaces import lattice as latt
import numpy as np

t = 1.0
U = -2.0
beta = 32.0
delta_initial = [1.0893001e-01, 1.0896048e-01]

id = int(sys.argv[1])
ud = float(sys.argv[2])
mu = float(sys.argv[3])

def graphene_lattice(N=1, onsite=None, t=1.0, t_perp=0.5):
    a = 0.24595  # [nm] unit cell length
    a_cc = 0.142  # [nm] carbon-carbon distance
    c0 = 0.335  # [nm] interlayer spacing

    a1 = a * np.array([1, 0])
    a2 = a * np.array([1 / 2, 1 / 2 * np.sqrt(3)])

    lat = latt.Lattice(a1=a1, a2=a2)
    onsite = [-ud * t]

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
    return lat

def main(N=3, onsite=None, t=1.0, t_perp=0.0):
    lattice = graphene_lattice(N, onsite, t, t_perp)

    nx = ny = 2
    lx = 256
    ly = 256
    mode = "random"
    bound = np.abs(3 + t_perp + ud + np.abs(mu)) + 0.15 * t
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-bound, bound],
        seed_h=id,
        seed_v=id,
    )
    pairing = kite.Pairing(lattice, hubbard=U)
    for l in range(N):
        A, B = f"A{l + 1}", f"B{l + 1}"
        pairing.add_onsite_pairing(A, delta=delta_initial[2 * l])
        pairing.add_onsite_pairing(B, delta=delta_initial[2 * l + 1])

    bdg = kite.BdG(chemical_potential=mu, beta=beta, hartree=[0.0, 0.0])
    calculation = kite.Calculation(configuration)

    # calculation.dos(
    #     num_points=4000,
    #     num_moments=4096,
    #     num_random=64,
    #     num_disorder=1
    # )
    num_iterations = 256
    N0 = 128 * 4.0 / 5.0
    prev_iterations = 128
    calculation.s_wave_clean(
            num_random=16,
            num_iterations=num_iterations,
            prev_iterations=prev_iterations,
            N0=N0,
            tau=1.14,
        )
    output_file = f"Data/monolayer_exp_twist_Id{id:03d}_Ud{ud:.2f}_Fermi{mu:.2f}.h5"
    kite.config_system(
        lattice, configuration, calculation,
        pairing=pairing, bdg=bdg, filename=output_file,
    )
    return output_file

if __name__ == "__main__":
    output = main(N=1, t_perp=0.0)
    print(output, file=sys.stderr)
