"""Density of states of a Pi-Flux Model

##########################################################################
#                         Copyright 2022, KITE                           #
#                         Home page: quantum-kite.com                    #
##########################################################################

Units: Energy in units of hopping, |t| = 1
Lattice: Square lattice
Configuration: Periodic boundary conditions, double precision, automatic rescaling
Calculation type: Average DOS
Last updated: 28/07/2022
"""

__all__ = ["main"]

import sys
sys.path.append("..")

from interfaces import kite
from interfaces import lattice as latt
import myaux
import numpy as np

filename = "Data/Piflux.h5"

def square_lattice(onsite=[0], t=1):
    """Return lattice specification for a square lattice with nearest neighbor hoppings"""
    a1 = np.array([2, 0])
    a2 = np.array([0, 1])
    lat = latt.Lattice(a1 = a1, a2 = a2)
    lat.add_sublattices(
        ('A', [0, 0], onsite[0]),
        ('B', [1, 0], onsite[0])
    )
    lat.add_hoppings(
        ([0, 0], 'A', 'B', -t),
        ([0, 1], 'A', 'A', -t),
        ([0, 1], 'B', 'B',  t),
        ([1, 0], 'B', 'A', -t),
    )
    return lat

def main(onsite=[0, 0], t=1):
    """Prepare the input file for KITEx"""
    lattice = square_lattice(onsite)

    nx = 2
    ny = 2
    lx = 128
    ly = 256
    mode = "random"

    imp_posA = myaux.mycentral(lx, ly)
    struc_disorder_A = kite.StructuralDisorder(lattice, position=imp_posA)
    struc_disorder_A.add_vacancy("A")
    disorder_structural = [struc_disorder_A]

    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-3.2, 3.2],
    )
    calculation = kite.Calculation(configuration)
    # calculation.dos(
    #     num_points=4000,
    #     # num_moments=256,
    #     num_moments=512,
    #     num_random=4,
    #     num_disorder=1
    # )
    calculation.ldos_map(energy_=0.0, sigma_=0.01, vectors_=32)
    kite.config_system(
        lattice,
        configuration,
        calculation,
        filename=filename,
        disorder_structural=disorder_structural,
    )
    return filename


if __name__ == "__main__":
    output = main()
    print(output, file=sys.stderr)
