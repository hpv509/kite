""" Density of states of a Pi-Flux Model

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
import numpy as np
import myaux

# t = float(sys.argv[1])
# a = float(sys.argv[2])
# m = float(sys.argv[3])
# lambda_R = float(sys.argv[4])  #: [eV] Rashba SOC, Phys. Rev. B 93, 155104
# lambda_sv = float(sys.argv[5])
# Nd = int(sys.argv[6])
# sg = int(sys.argv[7])
# U0 = float(sys.argv[8])
# lx = int(sys.argv[9])
# ly = int(sys.argv[10])
# eta = float(sys.argv[11])
# vectors = int(sys.argv[12])
# energy = float(sys.argv[13])
# NThreads= int(sys.argv[14])
# mode = sys.argv[15]
# vec_id=int(sys.argv[16])

t = 1
a = 1
m = 1
lambda_R = 1   #: [eV] Rashba SOC, Phys. Rev. B 93, 155104
lambda_sv = 1
Nd = 4
sg = 1
U0 = 2
lx = 32
ly = 32
eta = 0.02
vectors = 32
energy = 0.0
NThreads= 4
mode = 1
vec_id=1

def graphene_lattice(onsite=[0, 0], t=1., a_cc=1., km_so_A=0.0, km_so_B=0.0, rashba_so=0.0):
    """Return lattice specification for a square lattice with nearest neighbor hoppings"""
    a = np.sqrt(3) * a_cc  # unit cell length

    # define basis vectors
    a1 = np.array([+a / 2, a * np.sqrt(3) / 2])  #: [nm] unit cell vectors graphene
    a2 = np.array([-a / 2, a * np.sqrt(3)/ 2])

    # define orbitals vectors
    posA = np.array([0, 0.0])
    posB = np.array([0, -a_cc])

    # create a lattice with 2 primitive vectors
    lat = latt.Lattice(a1 = a1, a2 = a2)

    # Add sublattices
    lat.add_sublattices(
        # name, position, and onsite potential
        ('Aup', posA, onsite[0]),
        ('Bup', posB, onsite[1]),
        ('Adown', posA, onsite[0]),
        ('Bdown', posB, onsite[1])
    )

    # Add hoppings
    lat.add_hoppings(
        # ([f - i ], i , f )
        # inside the main cell, between which atoms, and the value
        ([0, 0], 'Adown', 'Bdown', -t),
        ([0, 0], 'Aup', 'Bup', -t),
        # between neighboring cells, between which atoms, and the value
        ([0, +1], 'Aup', 'Bup', -t),
        ([0, +1], 'Adown', 'Bdown', -t),

        ([+1, 0], 'Aup', 'Bup', -t),
        ([+1, 0], 'Adown', 'Bdown', -t)
    )

    if np.abs(lambda_R) > 0:
        lat.add_hoppings(
            # Rashba nearest neighbor, spin flip
            # inside the main cell, between which atoms, and the value
            ([0, 0], 'Aup', 'Bdown', -1.0 * rashba_so),  # delta1
            ([0, +1], 'Aup', 'Bdown', (+0.5 - np.sqrt(3) / 2 * 1j) * rashba_so),  # delta2
            ([+1, 0], 'Aup', 'Bdown', (+0.5 + np.sqrt(3) / 2 * 1j) * rashba_so),  # delta3

            ([0, 0], 'Adown', 'Bup', -1.0 * rashba_so),  # delta1
            ([0, +1], 'Adown', 'Bup', (+0.5 + np.sqrt(3) / 2 * 1j) * rashba_so),  # delta2
            ([+1, 0], 'Adown', 'Bup', (+0.5 - np.sqrt(3) / 2 * 1j) * rashba_so)  # delta3
        )

    if np.abs(lambda_sv) > 0:
        # Kane-Mele SOC, same spin next-nearest
        # between neighboring cells, between which atoms, and the value
        lat.add_hoppings(
            ([0, +1], 'Aup', 'Aup', -km_so_A),
            ([0, +1], 'Adown', 'Adown', +km_so_A),

            ([+1, 0], 'Aup', 'Aup', +km_so_A),
            ([+1, 0], 'Adown', 'Adown', -km_so_A),

            ([1, -1], 'Aup', 'Aup', -km_so_A),
            ([1, -1], 'Adown', 'Adown', +km_so_A),

            ([0, +1], 'Bup', 'Bup', +km_so_B),
            ([0, +1], 'Bdown', 'Bdown', -km_so_B),

            ([+1, 0], 'Bup', 'Bup', -km_so_B),
            ([+1, 0], 'Bdown', 'Bdown', +km_so_B),

            ([1, -1], 'Bup', 'Bup', +km_so_B),
            ([1, -1], 'Bdown', 'Bdown', -km_so_B),
        )

    return lat

def main(onsite=[0, 0], t=1., a=1., lambda_sv=0.0, lambda_R=0.0, Nd=0, sg=0, U0=0., lx=1024, ly=1024, eta=0.002, vectors=16, energy=0.0, Nthreads=4, mode="random"):
    """Prepare the input file for KITEx"""

    lambda_I_A = -lambda_sv
    lambda_I_B = +lambda_sv

    rashba_so = lambda_R * 2.0 * 1.0j / 3.0  #: [eV] constant and geometrical factors that will define Rashba SOC

    km_so_A = lambda_I_A * 1j / (3 * np.sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic A SOC
    km_so_B = lambda_I_B * 1j / (3 * np.sqrt(3))  #: [eV] constant and geometrical factors that will define intrinsic B SOC

    lattice = graphene_lattice(onsite, t, a, km_so_A, km_so_B, rashba_so)

    nx = int(np.sqrt(Nthreads))
    ny = int(np.sqrt(Nthreads))

    if Nd == 1:
        if U0 == 0.:
            # --- vacancy ---
            imp_posA = myaux.mycentral(lx, ly)
            struc_disorder_A = kite.StructuralDisorder(lattice, position = imp_posA)
            struc_disorder_A.add_vacancy('Aup')
            struc_disorder_A.add_vacancy('Adown')
            disorder_structural = [struc_disorder_A]
        else:
            if sg == 0:
                # --- bridge impurity ---
                imp_posA = myaux.mycentral(lx, ly)
                struc_disorder_A = kite.StructuralDisorder(lattice, position = imp_posA)
                struc_disorder_A.add_structural_disorder(([0.,0.], "Aup", U0),
                                                         ([0.,0.], "Adown", U0),
                                                         ([0.,0.], "Bup", U0),
                                                         ([0.,0.], "Bdown", U0))
                disorder_structural = [struc_disorder_A]
            elif sg == 1:
                # --- onsite A impurity ---
                imp_posA = myaux.mycentral(lx, ly)
                struc_disorder_A = kite.StructuralDisorder(lattice, position = imp_posA)
                struc_disorder_A.add_structural_disorder(([0.,0.], "Aup", U0),
                                                         ([0.,0.], "Adown", U0))
                disorder_structural = [struc_disorder_A]
            elif sg == 3:
                # --- hollow impurity ---
                imp_posA = myaux.mycentral(lx, ly)
                struc_disorder_A = kite.StructuralDisorder(lattice, position = imp_posA)
                struc_disorder_A.add_structural_disorder(([0.,0.], "Aup", U0),
                                                         ([0.,0.], "Adown", U0),
                                                         ([0.,0.], "Bup", U0),
                                                         ([0.,0.], "Bdown", U0),
                                                         ([-1.,0.], "Aup", U0),
                                                         ([-1.,0.], "Adown", U0),
                                                         ([0.,1.], "Bup", U0),
                                                         ([0.,1.], "Bdown", U0),
                                                         ([-1.,1.], "Aup", U0),
                                                         ([-1.,1.], "Adown", U0),
                                                         ([-1.,1.], "Bup", U0),
                                                         ([-1.,1.], "Bdown", U0))
                disorder_structural = [struc_disorder_A]
            elif sg == 4:
                #hollow imbalanced impurity
                imp_posA = myaux.mycentral(lx, ly)
                struc_disorder_A = kite.StructuralDisorder(lattice, position = imp_posA)
                Ua = U0
                Ub = U0*0.01-U0*0.99
                struc_disorder_A.add_structural_disorder(([0.,0.], "Aup", Ua),
                                                         ([0.,0.], "Adown", Ua),
                                                         ([0.,0.], "Bup", Ub),
                                                         ([0.,0.], "Bdown", Ub),
                                                         ([-1.,0.], "Aup", Ua),
                                                         ([-1.,0.], "Adown", Ua),
                                                         ([0.,1.], "Bup", Ub),
                                                         ([0.,1.], "Bdown", Ub),
                                                         ([-1.,1.], "Aup", Ua),
                                                         ([-1.,1.], "Adown", Ua),
                                                         ([-1.,1.], "Bup", Ub),
                                                         ([-1.,1.], "Bdown", Ub))
                disorder_structural = [struc_disorder_A]
            elif sg == 5:
                #hollow staggered impurity
                imp_posA = myaux.mycentral(lx, ly)
                struc_disorder_A = kite.StructuralDisorder(lattice, position = imp_posA)
                Ua = U0
                Ub = -U0
                struc_disorder_A.add_structural_disorder(([0.,0.], "Aup", Ua),
                                                         ([0.,0.], "Adown", Ua),
                                                         ([0.,0.], "Bup", Ub),
                                                         ([0.,0.], "Bdown", Ub),
                                                         ([-1.,0.], "Aup", Ua),
                                                         ([-1.,0.], "Adown", Ua),
                                                         ([0.,1.], "Bup", Ub),
                                                         ([0.,1.], "Bdown", Ub),
                                                         ([-1.,1.], "Aup", Ua),
                                                         ([-1.,1.], "Adown", Ua),
                                                         ([-1.,1.], "Bup", Ub),
                                                         ([-1.,1.], "Bdown", Ub))
                disorder_structural = [struc_disorder_A]
    else:
        disorder_structural = []

    disorder = kite.Disorder(lattice)
    disorder.add_disorder('Aup', 'Uniform', 0.0, 0.0)
    disorder.add_disorder('Adown', 'Uniform', 0.0, 0.0)
    disorder.add_disorder('Bup', 'Uniform', 0.0, 0.0)
    disorder.add_disorder('Bdown', 'Uniform', 0.0, 0.0)
    # V0 = float(sys.argv[17])
    # Rc = int(sys.argv[18])
    V0 = 1.0
    Rc = 1
    # sigma_c = int(sys.argv[19])
    sigma_c = 1.0

    energy_scale = 1.1 * (3*np.abs(t) +
                          np.max(np.abs(np.array([U0, V0, m]))) +
                          2*np.abs(lambda_sv)/np.sqrt(3) +
                          2*np.abs(lambda_R))

    configuration = kite.Configuration(
        # divisions=[nx, ny],
        # length=[lx, ly],
        divisions=[2, 2],
        length=[512, 512],
        boundaries=["open", "open"],
        is_complex=True,
        precision=0,
        spectrum_range=[-energy_scale, energy_scale],
        custom_potential=1,
        seed_h=vec_id,
        seed_v=vec_id,
    )
    calculation = kite.Calculation(configuration)

    calculation.ldos_map(
        energy_ = energy,
        sigma_ = eta,
        vectors_ = vectors
    )

    # output_file = "/users/svf517/scratch/raw_data_impurity/data_%s_stag%.1f_Rc%d_sgc%d_t%.2f_a%.2f_m%.2e_sv%.2e_R%.2e_Nd%d_sg%d_U0%.2f_L%dx%d_eta%.2e_Av%d_E%.3e_id%d.h5" % (mode, V0, Rc, sigma_c, t, a, m, lambda_sv, lambda_R, Nd, sg, U0, lx, ly, eta, vectors, energy, vec_id)
    output_file = "data_%s_stag%.1f_Rc%d_sgc%d_t%.2f_a%.2f_m%.2e_sv%.2e_R%.2e_Nd%d_sg%d_U0%.2f_L%dx%d_eta%.2e_Av%d_E%.3e_id%d.h5" % (mode, V0, Rc, sigma_c, t, a, m, lambda_sv, lambda_R, Nd, sg, U0, lx, ly, eta, vectors, energy, vec_id)

    if Nd == 0:
        kite.config_system(
            lattice,
            configuration,
            calculation,
            filename=output_file
        )
    else:
        kite.config_system(
            lattice,
            configuration,
            calculation,
            filename=output_file,
            disorder_structural = disorder_structural
        )

    return output_file

if __name__ == "__main__":
    # print(sys.argv[10])
    output = main([m, -m], t, a, lambda_sv, lambda_R, Nd, sg, U0, lx, ly, eta, vectors, energy, NThreads, mode)
    print(output, file = sys.stderr)
