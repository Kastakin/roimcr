use crate::utils::argsort;

pub enum MzErrorType {
    Dalton(f64),
    Ppm(f64),
}

impl MzErrorType {
    pub fn get_error(&self, mz: &f64) -> f64 {
        match self {
            Self::Dalton(x) => *x,
            Self::Ppm(x) => (x * mz) / f64::powf(10.0, 6.0),
        }
    }
}

pub enum MzRoiUpdater {
    Mean,
    Median,
    Weighted,
    Max,
}

impl MzRoiUpdater {
    pub fn calculate(&self, values: &[f64], intensities: &[f64]) -> f64 {
        match self {
            Self::Mean => values.iter().sum::<f64>() / values.len() as f64,
            Self::Median => {
                let ord_indices = argsort(values);
                values[ord_indices[ord_indices.len() / 2]]
            }
            Self::Max => {
                let ord_indices = argsort(intensities);
                values[ord_indices[ord_indices.len() - 1]]
            }
            Self::Weighted => {
                values
                    .iter()
                    .zip(intensities.iter())
                    .map(|(mz, int)| mz * int)
                    .sum::<f64>()
                    / intensities.iter().sum::<f64>()
            }
        }
    }
}
