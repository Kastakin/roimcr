pub mod options {
    use crate::enums::*;
    // pub struct Settings {
    //     roi_params: RoiParams,
    // }

    pub struct RoiParams {
        pub threshold: f64,
        pub t_factor: f64,
        pub mz_error: MzErrorType,
        pub mzroi_updater: MzRoiUpdater,
        pub min_occ: u32,
    }
}

pub mod io {
    use ndarray::prelude::*;
    use serde::Deserialize;

    #[derive(Debug, Deserialize, PartialEq)]
    pub struct Peak {
        pub mz: f64,
        pub intensity: f64,
    }

    #[derive(Debug, Deserialize, PartialEq)]
    pub struct Scan {
        pub num: u64,
        #[serde(rename = "msLevel")]
        pub ms_level: u8,
        #[serde(rename = "peaksCount")]
        pub peaks_count: u64,
        #[serde(rename = "retentionTime", with = "read_rt")]
        pub retention_time: f64,
        #[serde(with = "read_peaks")]
        pub peaks: Vec<Peak>,
    }

    #[derive(Default, Debug, Deserialize, PartialEq)]
    pub struct MsRun {
        #[serde(skip)]
        pub name: String,
        #[serde(rename = "scan")]
        pub scans: Vec<Scan>,
    }

    impl MsRun {
        pub fn new() -> Self {
            Default::default()
        }

        pub fn get_peaks(&self) -> Vec<Array2<f64>> {
            let mut result: Vec<Array2<f64>> = Vec::new();

            for scan in &self.scans {
                if scan.ms_level != 1 {
                    continue;
                }
                let mut peaks: Vec<f64> = Vec::new();
                let n_peaks: usize = scan.peaks.len();

                for peak in &scan.peaks {
                    peaks.push(peak.mz);
                    peaks.push(peak.intensity);
                }

                result.push(
                    Array::from_shape_vec((n_peaks, 2), peaks)
                        .expect("Couldn't get correct number of peaks"),
                );
            }
            result
        }

        pub fn get_times(&self) -> Array1<f64> {
            let mut result: Vec<f64> = Vec::new();

            for scan in &self.scans {
                if scan.ms_level != 1 {
                    continue;
                }
                result.push(scan.retention_time);
            }

            Array::from_vec(result)
        }

        pub fn get_tic(&self) -> Array1<f64> {
            let mut result: Vec<f64> = Vec::new();

            for scan in &self.scans {
                if scan.ms_level != 1 {
                    continue;
                }

                let int_sum: f64 = scan.peaks.iter().map(|peak| peak.intensity).sum();

                result.push(int_sum)
            }

            Array::from_vec(result)
        }
    }

    #[derive(Debug, Deserialize, PartialEq)]
    pub struct MzXML {
        #[serde(rename = "msRun")]
        pub ms_run: MsRun,
    }

    mod read_rt {
        use serde::{self, Deserialize, Deserializer};

        pub fn deserialize<'de, D>(deserializer: D) -> Result<f64, D::Error>
        where
            D: Deserializer<'de>,
        {
            let input = String::deserialize(deserializer).expect("Can't deserialize RT");
            let s = input.trim_start_matches('P');
            let s = s.trim_start_matches('T');
            let s = s.trim_end_matches('S');

            let result: f64 = s.parse().expect("Couldn't convert RT {s} to f32");
            Ok(result)
        }
    }

    mod read_peaks {
        use crate::structs::io::Peak;
        use serde::{self, Deserialize, Deserializer};

        pub fn deserialize<'de, D>(deserializer: D) -> Result<Vec<Peak>, D::Error>
        where
            D: Deserializer<'de>,
        {
            let s = String::deserialize(deserializer)?;
            let mut buffer = [0u8; 4];
            let mut peaks = Vec::<Peak>::new();

            let decoded = base64::decode(s).expect("couldn't decode");

            for peak in 0..(decoded.len() / 8) {
                let mut new_peak = Peak {
                    mz: 0.0,
                    intensity: 0.0,
                };
                for i in 0..2 {
                    let start: usize = (peak * 8) + i * 4;
                    let finish: usize = (peak * 8) + 4 + (i * 4);
                    buffer.clone_from_slice(&decoded[start..finish]);
                    match i {
                        0 => {
                            new_peak.mz = f32::from_be_bytes(buffer) as f64;
                        }
                        1 => {
                            new_peak.intensity = f32::from_be_bytes(buffer) as f64;
                        }
                        _ => unreachable!(),
                    }
                }
                peaks.push(new_peak);
            }

            Ok(peaks)
        }
    }
}

pub mod data {
    use crate::enums::MzRoiUpdater;

    #[derive(Default, Debug)]
    pub struct Roicell {
        pub mzs: Vec<Vec<f64>>,
        pub times: Vec<Vec<f64>>,
        pub intensities: Vec<Vec<f64>>,
        pub scans: Vec<Vec<f64>>,
        pub mz_roi: Vec<f64>,
    }

    impl Roicell {
        pub fn new() -> Self {
            Default::default()
        }

        fn update_mzroi(&mut self, roi: usize, updater: &MzRoiUpdater) {
            let new_value = updater.calculate(&self.mzs[roi], &self.intensities[roi]);
            if self.mz_roi.get(roi).is_none() {
                self.mz_roi.push(new_value);
            } else {
                self.mz_roi[roi] = new_value;
            }
            // let avg = self.mzs[roi].iter().sum::<f64>() / self.mzs[roi].len() as f64;
            // if self.mz_roi.get(roi).is_none() {
            //     self.mz_roi.push(avg);
            // } else {
            //     self.mz_roi[roi] = avg;
            // }
        }

        pub fn update_roi(
            &mut self,
            roi: usize,
            updater: &MzRoiUpdater,
            mz: f64,
            time: f64,
            intensity: f64,
            scan: usize,
        ) {
            // Aggiungere mz[i] alle masse dell'mzroi
            self.mzs[roi].push(mz);
            // Aggiungere time[i] ai tempi dell'mzroi
            self.times[roi].push(time);
            // Aggiungere intensities[i] alle intensità dell'mzroi
            self.intensities[roi].push(intensity);
            // Aggiungere scan alle scan dell'mzroi
            self.scans[roi].push(scan as f64);

            // Aggiornare la media dell'mzroi
            self.update_mzroi(roi, updater);
        }

        pub fn new_roi(&mut self) {
            self.mzs.push(Vec::new());
            self.times.push(Vec::new());
            self.intensities.push(Vec::new());
            self.scans.push(Vec::new());
        }

        pub fn merge_roi(&mut self, roi_keep: usize, roi_remove: usize, updater: &MzRoiUpdater) {
            self.mzs[roi_keep] =
                [self.mzs[roi_keep].clone(), self.mzs[roi_remove].clone()].concat();
            self.times[roi_keep] =
                [self.times[roi_keep].clone(), self.times[roi_remove].clone()].concat();
            self.intensities[roi_keep] = [
                self.intensities[roi_keep].clone(),
                self.intensities[roi_remove].clone(),
            ]
            .concat();
            self.scans[roi_keep] =
                [self.scans[roi_keep].clone(), self.scans[roi_remove].clone()].concat();
            self.update_mzroi(roi_keep, updater);

            self.remove_roi(roi_keep);
        }

        pub fn remove_roi(&mut self, roi: usize) {
            self.mzs.remove(roi);
            self.times.remove(roi);
            self.intensities.remove(roi);
            self.scans.remove(roi);
            self.mz_roi.remove(roi);
        }
    }
}
