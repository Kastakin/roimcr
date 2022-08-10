from typing import Tuple

import numpy as np
from numpy.typing import NDArray


def mrdivide(a: NDArray, b: NDArray) -> NDArray:
    return np.linalg.solve(b.conj().T, a.conj().T).conj().T


def mldivide(a: NDArray, b: NDArray) -> NDArray:
    return np.linalg.lstsq(a, b)


def pcarep(
    og_x: NDArray, n_features: int
) -> Tuple[NDArray, NDArray, NDArray, NDArray, float]:
    u, s, vh = np.linalg.svd(og_x, full_matrices=False)
    u = u[:, :n_features]
    s = s[:n_features, :n_features]
    vh = vh[:, :n_features]

    x = u @ s @ vh
    res = og_x - x

    sst1 = (res**2).sum(dtype=float)
    sst2 = (og_x**2).sum(dtype=float)

    lof = (np.sqrt(sst1 / sst2)) * 100

    return u, s, vh, x, lof
