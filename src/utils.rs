use csv::WriterBuilder;
use ndarray::Array2;
use ndarray_csv::Array2Writer;

use std::fs::File;

pub fn subset(array: &[f64], indices: &[usize]) -> Vec<f64> {
    indices.iter().map(|&i| array[i]).collect::<Vec<_>>()
}

pub fn argsort(arr: &[f64]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..arr.len()).collect();
    indices.sort_unstable_by(move |&i, &j| arr[i].total_cmp(&arr[j]));
    indices
}

pub fn sort_by_indices<T>(data: &mut [T], indices: &[usize]) {
    let mut indices = indices.to_owned();
    for idx in 0..data.len() {
        if indices[idx] != idx {
            let mut current_idx = idx;
            loop {
                let target_idx = indices[current_idx];
                indices[current_idx] = current_idx;
                if indices[target_idx] == target_idx {
                    break;
                }
                data.swap(current_idx, target_idx);
                current_idx = target_idx;
            }
        }
    }
}

pub fn ndarray_to_csv(data: &Array2<f64>) {
    let file = File::create("examples/test.csv").expect("Couldn't create csv file!");
    let mut writer = WriterBuilder::new().has_headers(false).from_writer(file);
    writer
        .serialize_array2(data)
        .expect("Couldn't serialize array to csv file!");
}
