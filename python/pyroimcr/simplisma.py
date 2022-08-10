"""
Module implementing Simplisma algorithm (Winding et al. 1991).
Simplisma compute 'purest' spectra by beign spectra that explain the most data variance.

Functions
---------
simplisma
    Apply Simplsima algorithm
"""
import numpy as np
from numpy.typing import NDArray


def _compute_weights_matrix(
    coo_matrix: NDArray, max_purities_idx: NDArray, comp: int, var_idx: int
) -> NDArray:
    """Compute the weights matrix to find the next purest spectrum.

    Parameters
    ----------
    coo_matrix : NDArray
        The correlation around the origin matrix
    max_purities_idx : NDArray
        Array with the indices of purest spectra previously found.
    comp : int
        The researched purest spectrum number.
    var_idx : int
        The studied variable.

    Returns
    -------
    NDArray
        The weights matrix to take account of explained variance by other purest spectra.
    """
    weights_matrix = np.zeros((comp + 1, comp + 1), dtype=coo_matrix.dtype)
    weights_matrix[0, 0] = coo_matrix[var_idx, var_idx]
    for i in range(comp):
        purest_i_idx = np.int(max_purities_idx[i])
        weights_matrix[0, i + 1] = coo_matrix[var_idx, purest_i_idx]
        weights_matrix[i + 1, 0] = coo_matrix[purest_i_idx, var_idx]
        for j in range(comp):
            purest_j_idx = np.int(max_purities_idx[j])
            weights_matrix[i + 1, j + 1] = coo_matrix[purest_i_idx, purest_j_idx]
    return weights_matrix


def simplisma(data: NDArray, nb_purest: int, noise: float) -> NDArray:
    """Simplisma alogrithm [1]. Find 'purest' spectra based on the explained variance by spectra.

    Parameters
    ----------
    data : NDArray
        The input data.
    nb_purest : int
        The number of purest spectra to find.
    noise : float
        A small fraction of the maximmum mean in data to avoid division by zero.

    Returns
    -------
    NDArray
        The nb_purest purest spectra.

    References
    ----------
    [1] Windig, W., & Guilment, J. (1991). Interactive self-modeling mixture analysis. Analytical chemistry, 63(14), 1425-1432.
    """
    scaled_data = np.zeros(data.shape, dtype=data.dtype)
    max_purities_idx = np.zeros(nb_purest, dtype=np.int)
    w = np.zeros((nb_purest, data.shape[1]), dtype=data.dtype)

    mean = np.mean(data, axis=0)
    noise = noise / 100
    noise = np.max(mean) * noise
    std = np.std(data, axis=0)
    p_initial = std / (mean + noise)

    length = np.sqrt((std * std) + ((mean + noise) * (mean + noise)))
    for j in range(data.shape[1]):
        scaled_data[:, j] = data[:, j] / length[j]
    coo_matrix = np.dot(scaled_data.T, scaled_data) / data.shape[0]

    w[0] = (std * std) + (mean * mean)
    w[0] /= length * length
    p_initial *= w[0]

    max_purities_idx[0] = np.int(np.argmax(p_initial))

    for i in range(1, nb_purest):
        p_i = np.zeros(data.shape[1], dtype=data.dtype)
        for j in range(data.shape[1]):
            weights_matrix = _compute_weights_matrix(coo_matrix, max_purities_idx, i, j)
            w[i, j] = np.linalg.det(weights_matrix)
            p_i[j] = p_initial[j] * w[i, j]
        max_purities_idx[i] = np.int(np.argmax(p_i))

    purest = np.zeros((data.shape[0], nb_purest), dtype=data.dtype)
    for i in range(nb_purest):
        purest[:, i] = data[:, np.int(max_purities_idx[i])]
    return purest
