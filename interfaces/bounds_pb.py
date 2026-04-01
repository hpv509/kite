import numpy as np
import pybinding as pb

from scipy.spatial import cKDTree

def make_pybinding_model(lattice, disorder=None, disorder_structural=None, **kwargs):
    """Build a Pybinding model with disorder used in Kite. Bond disorder or magnetic field are not currently supported.

    Parameters
    ----------
    lattice : pb.Lattice
        Pybinding lattice object that carries the info about the unit cell vectors, unit cell cites, hopping terms and
        onsite energies.
    disorder : Disorder
        Class that introduces Disorder into the initially built lattice. For more info check the Disorder class.
    disorder_structural : StructuralDisorder
        Class that introduces StructuralDisorder into the initially built lattice. For more info check the
        StructuralDisorder class.
    **kwargs: Optional arguments like shape .

    """

    shape = kwargs.get('shape', None)
    if disorder_structural:
        # check if there's a bond disorder term
        # return an error if so
        disorder_struc_list = disorder_structural
        if not isinstance(disorder_structural, list):
            disorder_struc_list = [disorder_structural]

        for idx_struc, dis_struc in enumerate(disorder_struc_list):
            if len(dis_struc._sub_from):
                raise SystemExit(
                    'Automatic scaling is not supported when bond disorder is specified. Please select the scaling '
                    'bounds manually.')

    def gaussian_disorder(sub, mean_value, stdv):
        """Add gaussian disorder with selected mean value and standard deviation to the pybinding model.

        Parameters
        ----------
        sub : str
            Select a sublattice where disorder should be added.
        mean_value : float
            Select a mean value of the disorder.
        stdv : float
            Select standard deviation of the disorder.
        """

        @pb.onsite_energy_modifier
        def modify_energy(energy, sub_id):
            rand_onsite = np.random.normal(loc=mean_value, scale=stdv, size=len(energy[sub_id == sub]))
            energy[sub_id == sub] += rand_onsite
            return energy

        return modify_energy

    def deterministic_disorder(sub, mean_value):
        """Add deterministic disorder with selected mean value to the Pybinding model.

        Parameters
        ----------
        sub : str
            Select a sublattice where disorder should be added.
        mean_value : float
            Select a mean value of the disorder.

        """

        @pb.onsite_energy_modifier
        def modify_energy(energy, sub_id):
            onsite = mean_value * np.ones(len(energy[sub_id == sub]))
            energy[sub_id == sub] += onsite
            return energy

        return modify_energy

    def uniform_disorder(sub, mean_value, stdv):
        """Add uniform disorder with selected mean value and standard deviation to the Pybinding model.

        Parameters
        ----------
        sub : str
            Select a sublattice where disorder should be added.
        mean_value : float
            Select a mean value of the disorder.
        stdv : float
            Select standard deviation of the disorder.
        """

        @pb.onsite_energy_modifier
        def modify_energy(energy, sub_id):
            a = mean_value - stdv * np.sqrt(3)
            b = mean_value + stdv * np.sqrt(3)

            rand_onsite = np.random.uniform(low=a, high=b, size=len(energy[sub_id == sub]))
            energy[sub_id == sub] += rand_onsite

            return energy

        return modify_energy

    def vacancy_disorder(sub, concentration):
        """Add vacancy disorder with selected concentration to the Pybinding model.

        Parameters
        ----------
        sub : str
            Select a sublattice where disorder should be added.
        concentration : float
            Concentration of the vacancy disorder.
        """

        @pb.site_state_modifier(min_neighbors=2)
        def modifier(state, sub_id):
            rand_vec = np.random.rand(len(state))
            vacant_sublattice = np.logical_and(sub_id == sub, rand_vec < concentration)

            state[vacant_sublattice] = False
            return state

        return modifier

    def local_onsite_disorder(positions, value):
        """Add onsite disorder as a part of StructuralDisorder class to the Pybinding model.

        Parameters
        ----------
        positions : np.ndarray
            Select positions where disorder should appear
        value : np.ndarray
            Value of the disordered onsite term.
        """
        space_size = np.array(positions).shape[1]

        @pb.onsite_energy_modifier
        def modify_energy(x, y, z, energy):
            # all_positions = np.column_stack((x, y, z))[0:space_size, :]
            all_positions = np.stack([x, y, z], axis=1)[:, 0:space_size]

            kdtree1 = cKDTree(positions)
            kdtree2 = cKDTree(all_positions)

            d_max = 0.05
            # find the closest elements between two trees, with d < d_max. Finds the desired elements from the
            # parameters x, y, z being used inside the modifier function.
            coo = kdtree1.sparse_distance_matrix(kdtree2, d_max, output_type='coo_matrix')

            energy[coo.col] += value

            return energy

        return modify_energy

    if not shape:

        vectors = np.asarray(lattice.vectors)
        space_size = vectors.shape[0]

        # fix a size for 1D, 2D, or 3D
        referent_size = 1

        if space_size == 1:
            referent_size = 1000
        elif space_size == 2:
            referent_size = 200
        elif space_size == 3:
            referent_size = 50

        norm = np.sum(np.abs(vectors)**2, axis=-1)**(1./2)

        num_each_dir = (referent_size / norm).astype(int)

        shape_size = 2 * num_each_dir[0:space_size]
        symmetry = 1 * num_each_dir[0:space_size]

        shape = pb.primitive(*shape_size)
        trans_symm = pb.translational_symmetry(*symmetry)
        param = [shape, trans_symm]
    else:
        param = [shape]

    model = pb.Model(lattice, *param)

    if disorder:
        disorder_list = disorder
        if not isinstance(disorder, list):
            disorder_list = [disorder]
        for dis in disorder_list:
            for idx in range(len(dis._type)):
                if dis._type[idx].lower() == 'uniform':
                    model.add(
                        uniform_disorder(dis._sub_name[idx], dis._mean[idx], dis._stdv[idx]))
                if dis._type[idx].lower() == 'gaussian':
                    model.add(
                        gaussian_disorder(dis._sub_name[idx], dis._mean[idx], dis._stdv[idx]))
                if dis._type[idx].lower() == 'deterministic':
                    model.add(deterministic_disorder(dis._sub_name[idx], dis._mean[idx]))

    if disorder_structural:

        disorder_struc_list = disorder_structural
        if not isinstance(disorder_structural, list):
            disorder_struc_list = [disorder_structural]

        for idx_struc, dis_struc in enumerate(disorder_struc_list):
            num_sites = model.system.num_sites
            rand_vec = np.random.rand(num_sites)
            space_size = np.array(lattice.vectors).shape[0]

            for vac in dis_struc._vacancy_sub:
                model.add(vacancy_disorder(sub=vac, concentration=dis_struc._concentration))

            for idx in range(len(dis_struc._sub_onsite)):
                names, sublattices = zip(*model.lattice.sublattices.items())

                sublattice_alias = names.index(dis_struc._sub_onsite[idx])

                select_sublattice = model.system.sublattices == sublattice_alias
                sub_and_rand = np.logical_and(select_sublattice, rand_vec < dis_struc._concentration)

                # generates a set of random positions that will be added to the nodes in structural disorder, problem
                # because when only one sublattice is selected, effective concentration will be lower
                positions = np.stack([model.system.positions.x[sub_and_rand],
                                      model.system.positions.y[sub_and_rand],
                                      model.system.positions.z[sub_and_rand]], axis=1)[:, 0:space_size]

                # ref_pos = np.stack([model.system.positions.x[def_site_idx], model.system.positions.y[def_site_idx],
                #                     model.system.positions.z[def_site_idx]], axis=1)[:, 0:space_size]

                # get the position of onsite disordered sublattice
                vectors = np.array(lattice.vectors)[:, 0:space_size]
                pos_sub = lattice.sublattices[dis_struc._sub_onsite[idx]].position[0:space_size]

                # make an array of positions of sites where the onsite disorder will be added
                pos = pos_sub + np.dot(vectors.T, np.array(dis_struc._rel_idx_onsite[idx]))
                select_pos = positions + pos

                # add the onsite with value dis_struc._onsite[idx]
                model.add(local_onsite_disorder(select_pos, dis_struc._onsite[idx]))

    return model


def estimate_bounds(lattice, disorder=None, disorder_structural=None):
    model = make_pybinding_model(lattice, disorder, disorder_structural)
    kpm = pb.kpm(model)
    a, b = kpm.scaling_factors
    return -a/0.9 + b, a/0.9 + b
