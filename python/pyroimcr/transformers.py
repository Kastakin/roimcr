"""
Module implementing the ROI extraction procedure (Winding et al. 1991).
This transformer tries to adhere to sklearn formalism for interoperability

Classes
---------
RoiExtractor
    Skleanr compliemnt transformer for roi extraction
"""

from typing import List, Tuple

from numpy.typing import NDArray

from .pyroimcr import get_roi


class RoiExtractor:
    def __init__(self) -> None:
        self.threshold = None
        self.t_factor = None
        self.tol = None

    def fit(self, *, X, y=None):
        return self

    def transform(
        self, *, X: List[List[NDArray]], y: List[NDArray]
    ) -> Tuple[NDArray, NDArray]:
        mzroi, msroi = get_roi(X, y, threshold=self.threshold)
        return mzroi, msroi

    def fit_transform(
        self, *, X: List[List[NDArray]], y: List[NDArray]
    ) -> Tuple[NDArray, NDArray]:
        return self.transform(X=X)
