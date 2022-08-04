mod enums;
mod file;
mod structs;
mod utils;

use ndarray::prelude::*;
use std::path::Path;
use utils::ndarray_to_csv;

use crate::enums::{MzErrorType, MzRoiUpdater};
use crate::file::{load_mzxml, load_netcdf};
use crate::structs::data::Roicell;
use crate::structs::options::RoiParams;
use crate::utils::{argsort, sort_by_indices, subset};

pub fn load_data(path: &str) -> (Vec<Array2<f64>>, Array1<f64>) {
    let parsed = load_mzxml(Path::new(path)).expect("Error in parsing the file");
    let peaks = parsed.get_peaks();
    let times = parsed.get_times();
    (peaks, times)
}

pub fn compute_roi(
    peaks_list: &[Array2<f64>],
    times: &Array1<f64>,
    settings: RoiParams,
) -> (Vec<f64>, Roicell) {
    let mut mzroi: Vec<f64> = Vec::new();
    let mut roicell: Roicell = Roicell::new();

    for (scan, peaks) in peaks_list.iter().enumerate() {
        if peaks.is_empty() {
            continue;
        }
        let indices = peaks
            .axis_iter(Axis(0))
            .enumerate()
            .filter(|(_, r)| r[1] > settings.threshold)
            .map(|(index, _)| index)
            .collect::<Vec<_>>();

        if indices.is_empty() {
            continue;
        }

        let mz = subset(&peaks.index_axis(Axis(1), 0).to_vec(), &indices);
        let intensities = subset(&peaks.index_axis(Axis(1), 1).to_vec(), &indices);

        if scan == 1 {
            mzroi.push(mz[0]);
        }

        for i in 0..mz.len() {
            let compatible_rois: Vec<usize> = mzroi
                .iter()
                .enumerate()
                .filter(|(_, &roi_value)| {
                    (roi_value - mz[i]).abs() <= settings.mz_error.get_error(&mz[i])
                })
                .map(|(index, _)| index)
                .collect();

            if compatible_rois.is_empty() {
                roicell.new_roi();
                roicell.update_roi(
                    roicell.mzs.len() - 1,
                    &settings.mzroi_updater,
                    mz[i],
                    times[scan],
                    intensities[i],
                    scan,
                );
                mzroi.push(roicell.mz_roi[roicell.mzs.len() - 1])
            } else {
                for &roi in compatible_rois.iter() {
                    roicell.update_roi(
                        roi,
                        &settings.mzroi_updater,
                        mz[i],
                        times[scan],
                        intensities[i],
                        scan,
                    );
                    mzroi[roi] = roicell.mz_roi[roi];
                }
            }
        }
    }

    let sorted_roi_index = argsort(&mzroi);

    sort_by_indices(&mut mzroi, &sorted_roi_index);

    sort_by_indices(&mut roicell.mzs, &sorted_roi_index);
    sort_by_indices(&mut roicell.times, &sorted_roi_index);
    sort_by_indices(&mut roicell.intensities, &sorted_roi_index);
    sort_by_indices(&mut roicell.scans, &sorted_roi_index);
    sort_by_indices(&mut roicell.mz_roi, &sorted_roi_index);

    // Raggruppa roi con mz simili
    loop {
        let mz_error_roi: Vec<f64> = mzroi.windows(2).map(|w| (w[0] - w[1]).abs()).collect();

        let mergable_roi: Vec<usize> = mz_error_roi
            .iter()
            .enumerate()
            .filter(|(_, &mz_diff)| mz_diff < settings.mz_error.get_error(&mzroi[mzroi.len() - 2]))
            .map(|(index, _)| index)
            .collect();

        if mergable_roi.is_empty() {
            break;
        }

        mzroi[mergable_roi[0]] = roicell.mz_roi[mergable_roi[0] + 1];
        roicell.merge_roi(
            mergable_roi[0],
            mergable_roi[0] + 1,
            &settings.mzroi_updater,
        );
        mzroi.remove(mergable_roi[0] + 1);
    }

    // Filtrare per minimo numero di mz e minima intensità massima in roi
    let roi_num: Vec<u32> = roicell.mzs.iter().map(|i| i.len() as u32).collect();
    let max_int: Vec<f64> = roicell
        .intensities
        .iter()
        .map(|i| *i.iter().max_by(|a, b| a.total_cmp(b)).unwrap())
        .collect();

    let mut removable_roi: Vec<usize> = roi_num
        .iter()
        .zip(max_int.iter())
        .enumerate()
        .filter(|(_, (&n, &int))| {
            (n <= settings.min_occ) || (int <= settings.threshold * settings.t_factor)
        })
        .map(|(i, (_, _))| i)
        .collect::<Vec<usize>>();

    removable_roi.reverse();

    for roi in removable_roi {
        roicell.remove_roi(roi);
        mzroi.remove(roi);
    }

    let mut msroi: Array2<f64> = Array2::zeros((peaks_list.len(), mzroi.len()));

    for i in 0..mzroi.len() {
        let nval = roicell.scans[i].len();
        for j in 0..nval {
            let irow = roicell.scans[i][j] as usize;
            let msi = roicell.intensities[i][j];
            msroi[[irow, i]] += msi;
        }
    }
    return (mzroi, roicell);
}
