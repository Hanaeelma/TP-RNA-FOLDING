import numpy as np

def calculate_gibbs_free_energy(distances):
    """
    Calculate the estimated Gibbs free energy for an RNA structure.

    Parameters:
        distances (dict): A dictionary containing the distances between residue pairs.

    Returns:
        float: The estimated Gibbs free energy.
    """
    energies = []

    for pair, dist_list in distances.items():
        # Interpolate distances between i and i+4
        interpolated_distances = np.interp(range(len(dist_list)), [0, len(dist_list) - 1], [dist_list[0], dist_list[-1]])

        # Calculate scores using linear interpolation
        scores = np.interp(interpolated_distances, [0, 20], [10, 0])

        # Sum all scores
        energy = np.sum(scores)
        energies.append(energy)

    # Calculate the estimated Gibbs free energy as the sum of all energies
    gibbs_energy = np.sum(energies)
    return gibbs_energy

