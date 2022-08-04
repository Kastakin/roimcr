use numpy::{PyArray1, PyArray2, ToPyArray};
use pyo3::prelude::*;
use roimcr::load_data;

#[pyfunction]
/// Formats the sum of two numbers as string.
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

#[pyfunction]
/// Formats the sum of two numbers as string.
fn import_data<'py>(
    py: Python<'py>,
    path: &str,
) -> PyResult<(Vec<&'py PyArray2<f64>>, &'py PyArray1<f64>)> {
    let (peaks, times) = load_data(path);

    Ok((
        peaks
            .iter()
            .cloned()
            .map(|scan| scan.to_pyarray(py))
            .collect(),
        times.to_pyarray(py),
    ))
}

/// A Python module implemented in Rust.
#[pymodule]
fn pyroimcr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(import_data, m)?)?;
    Ok(())
}
