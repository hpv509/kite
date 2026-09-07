import sys

sys.path.append("..")

from interfaces import kite
from interfaces import custom
from interfaces import lattice as latt
import numpy as np

seed = int(sys.argv[1])
# flag = int(sys.argv[2])
size = int(sys.argv[2])
eta = float(sys.argv[3])

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
    lx = ly = size
    W = 1.0 * t
    sigma = eta / 16.0
    energy_scale = 2.7 * t + 0.5 * W
    mem = 2

    disorder = kite.Disorder(lattice)
    disorder.add_disorder("px", "Uniform", 0.0, W / np.sqrt(12.0))
    disorder.add_disorder("py", "Uniform", 0.0, W / np.sqrt(12.0))
    disorder.add_disorder("d", "Uniform", 0.0, W / np.sqrt(12.0))
    mode = "random"
    configuration = kite.Configuration(
        divisions=[nx, ny],
        length=[lx, ly],
        boundaries=[mode, mode],
        is_complex=True,
        precision=1,
        spectrum_range=[-energy_scale, energy_scale],
        seed=seed,
    )
    calculation = kite.Calculation(configuration)
    # if flag == 0:
    #     A = custom.Vertex(0, [[1.0j, "vy.rx"], [1.0j, "rx.vy"]])
    #     B = custom.Vertex(0, [[1.0j, "vy"]])
    # else:
    #     A = custom.Vertex(0, [[1.0j, "vy"]])
    #     B = custom.Vertex(0, [[1.0j, "vy.rx"], [1.0j, "rx.vy"]])
    A = custom.Vertex(0, [[1.0j, "vy"]])
    B = custom.Vertex(0, [[1.0j, "vy"]])
    calculation.custom_singleshot_two(
        stream_=[A, B],
        num_random_=2,
        num_disorder_=2,
        sigma_=sigma,
        energies_=[0.0],
        gamma_=eta,
    )
    # output_file = f"Data/SS_L{size:03d}_F{flag:02d}_Eta{eta:.3f}_Seed{seed:03d}.h5"
    output_file = f"Data/SS_L{size:03d}_Eta{eta:.3f}_Seed{seed:03d}.h5"
    kite.config_system(
        lattice, configuration, calculation, filename=output_file, disorder=disorder
    )
    return output_file


if __name__ == "__main__":
    output = main()
    print(output, file=sys.stderr)
