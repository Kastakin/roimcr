from typing import Dict, List, Tuple

from numpy.typing import NDArray

def import_data(path: str, format: str) -> Tuple[List[NDArray], NDArray]:
    """
    Get peaks and relatives times from centroided LC-MS data in mzxml or netcdf format.

    Parameters
    ----------
    path : str
        the path to the LC-MS run file
    format : str
        file format of the LC-MS run file, only mzxml' and 'netcdf' files are currently supported

    Returns
    -------
    peaks : list[ndarray]
        list of of peaks present in  each scan
    times : ndarray
        elution times of each scan
    """

def get_roi(
    peaks_list: List[List[NDArray]],
    times_list: List[NDArray],
    *,
    threshold: float = 1000.0,
    t_factor: float = 1.0,
    tol: float = 0.5,
    tol_units: str = "da",
    update_method: str = "mean",
    min_occ: int = 1
) -> Tuple[NDArray, NDArray]:
    """Given the peaks lists of one or more LC-MS runs return their ROIs

    Args:
        peaks_list (List[List[ndarray]]): List of of the peaks lists of all the runs to consider
        times_list (List[ndarray]): List of the scan times for all the considered runs
        threshold (float, optional): Min intensity for a peak to be considered. Defaults to 1000.0.
        t_factor (float, optional): Minimum Signal to Noise ration for a ROI to be considered. Defaults to 1.0.
        tol (float, optional): Maximum difference allowed between to peaks mz value to consider them alike. Defaults to 0.5.
        tol_units (str, optional): Unit of `tol`; allowed values are `da` for Daltons and `ppm`. Defaults to "da".
        update_method (str, optional): Which metrics to use when computing the mz value of a ROI. Defaults to "mean".
        min_occ (int, optional): Minimum number of peaks to be present in a ROI for it to be considered. Defaults to 1.

    Returns:
    -------
    mzroi : ndarray
        list of mz value representative of each ROI
    msroi : ndarray
        matrix with dimensions (n_times x n_rois), element msroi[i,j] represent the intensity for ROI j at time i
    """
