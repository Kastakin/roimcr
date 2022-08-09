use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;
use roimcr::enums::{MzErrorType, MzRoiUpdater};
use roimcr::structs::options::RoiParams;
use roimcr::{compute_roi, load_data};

#[pyfunction]
/// Load a file .
fn import_data<'py>(
    py: Python<'py>,
    path: &str,
    format: &str,
) -> PyResult<(Vec<&'py PyArray2<f64>>, &'py PyArray1<f64>)> {
    let (peaks, times) = load_data(path, format);

    Ok((
        peaks
            .iter()
            .cloned()
            .map(|scan| scan.to_pyarray(py))
            .collect(),
        times.to_pyarray(py),
    ))
}

#[pyfunction(
    py_args = "*",
    threshold = "1000.0",
    t_factor = "1.0",
    tol = "0.5",
    tol_units = "\"da\"",
    update_method = "\"mean\"",
    min_occ = "1"
)]
fn get_roi<'py>(
    py: Python<'py>,
    peaks_list: Vec<Vec<PyReadonlyArray2<'_, f64>>>,
    times_list: Vec<PyReadonlyArray1<'_, f64>>,
    threshold: f64,
    t_factor: f64,
    tol: f64,
    tol_units: &str,
    update_method: &str,
    min_occ: u32,
) -> (&'py PyArray1<f64>, &'py PyArray2<f64>) {
    let mzroi_updater = match update_method.to_lowercase().as_str() {
        "mean" => MzRoiUpdater::Mean,
        "max" => MzRoiUpdater::Max,
        "median" => MzRoiUpdater::Median,
        "weighted" => MzRoiUpdater::Weighted,
        _ => panic!(
            "Unknown update method! Allowed methods are: 'mean', 'max', 'median' and 'weighted'"
        ),
    };

    let mz_error = match tol_units.to_lowercase().as_str() {
        "da" => MzErrorType::Dalton(tol),
        "ppm" => MzErrorType::Ppm(tol),
        _ => {
            panic!("Unknown update method! Allowed methods are: 'da' for Daltons and 'ppm' for ppm")
        }
    };

    let settings = RoiParams {
        threshold,
        t_factor,
        mz_error,
        mzroi_updater,
        min_occ,
    };

    let peaks_list: Vec<Vec<Array2<f64>>> = peaks_list
        .iter()
        .map(|run| {
            run.iter()
                .map(|peaks| {
                    peaks
                        .downcast::<PyArray2<f64>>()
                        .unwrap()
                        .readonly()
                        .as_array()
                        .to_owned()
                })
                .collect()
        })
        .collect();

    let times_list: Vec<Array1<f64>> = times_list
        .iter()
        .map(|times| {
            times
                .downcast::<PyArray1<f64>>()
                .unwrap()
                .readonly()
                .as_array()
                .to_owned()
        })
        .collect();

    let (mzroi, msroi, roicell) = compute_roi(&peaks_list, &times_list, settings);

    (mzroi.into_pyarray(py), msroi.into_pyarray(py))
}

/// A Python module implemented in Rust.
#[pymodule]
fn pyroimcr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(import_data, m)?)?;
    m.add_function(wrap_pyfunction!(get_roi, m)?)?;
    Ok(())
}
