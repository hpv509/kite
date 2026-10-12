import sys
sys.path.insert(0, "../../interfaces")

import numpy as np
from lattice import Lattice
import kite

t = 1.0
a = 1.0
a1 = np.array([3*a/2,  np.sqrt(3)*a/2])
a2 = np.array([3*a/2, -np.sqrt(3)*a/2])
lattice = Lattice(a1=a1, a2=a2)
lattice.add_sublattices(
    ('A', [0.0, 0.0], 0.0),
    ('B', [a,   0.0], 0.0)
)
lattice.add_hoppings(
    ([ 0,  0], 'A', 'B', -t),
    ([-1,  0], 'A', 'B', -t),
    ([ 0, -1], 'A', 'B', -t)
)

mu = 0.4
beta = 100.0
bdg = kite.BdG(chemical_potential=mu, beta=beta)
delta_s0 = 0.3
delta_p0 = -0.3

pairing = kite.Pairing(
            lattice,
            interaction=-2.0,
            hubbard=-2.0,
            onsite_delta=delta_s0
            )

pairing.add_pairing([ 0,  0], 'A', 'B', delta_p0)
pairing.add_pairing([-1,  0], 'A', 'B', delta_p0)
pairing.add_pairing([ 0, -1], 'A', 'B', delta_p0)

configuration = kite.Configuration(
                    divisions=[2, 2],
                    length=[32, 32],
                    boundaries=['periodic', 'periodic'],
                    is_complex=True,
                    precision=1,
                    seed_h=1,
                    seed_v=2
                    )

calculation = kite.Calculation(configuration)

calculation.sp_wave_clean(
            num_random_transient=16,
            transient_iterations=20,
            tail_stages=3,
            stage_iterations=20,
            depth=2,
            mixing=0.3,
            average_bonds=True,
            )

filename = 'sp_wave_clean_graphene.h5'

kite.config_system(
    lattice,
    configuration,
    calculation,
    filename=filename,
    pairing=pairing,
    bdg=bdg
    )
