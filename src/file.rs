use crate::structs::io::{MsRun, MzXML, Peak, Scan};
use crate::utils::subset;
use netcdf::*;
use quick_xml::de::{from_str, DeError};
use std::{fs, path::Path};

pub fn load_mzxml(path: &Path) -> Result<MsRun, DeError> {
    let input = fs::read_to_string(path).expect("Error in reading the file");
    let mut result: MzXML = from_str(&input)?;
    result.ms_run.name = path
        .file_stem()
        .expect("No file name")
        .to_str()
        .expect("Can't convert filname to string")
        .to_string();

    Ok(result.ms_run)
}

pub fn load_netcdf(path: &Path) -> Result<MsRun, DeError> {
    let mut run: MsRun = MsRun::new();
    run.name = path
        .file_stem()
        .expect("No file name")
        .to_str()
        .expect("Can't convert filname to string")
        .to_string();
    let file = open(path).expect("Can't read netcdf file");
    // for var in file.variables(){
    //     println!("{}", var.name())
    // }

    let scan_index = file
        .variable("scan_index")
        .unwrap_or_else(|| panic!("Can't find scan_index in variables"))
        .values::<u32>(None, None)
        .expect("Can't read values")
        .into_raw_vec();

    let points_per_scan = file
        .variable("point_count")
        .unwrap_or_else(|| panic!("Can't find point_count in variables"))
        .values::<u32>(None, None)
        .expect("Can't read values")
        .into_raw_vec();

    let scans = scan_index
        .iter()
        .zip(points_per_scan.iter())
        .map(|(&a, &b)| ((a as usize)..((a + b) as usize)).collect())
        .collect::<Vec<Vec<usize>>>();

    let scan_time = file
        .variable("scan_acquisition_time")
        .unwrap_or_else(|| panic!("Can't find scan_acquisition_time in variables"))
        .values::<f64>(None, None)
        .expect("Can't read values")
        .into_raw_vec();

    let intensity_values = file
        .variable("intensity_values")
        .unwrap_or_else(|| panic!("Can't find intensity_values in variables"))
        .values::<f64>(None, None)
        .expect("Can't read values")
        .into_raw_vec();

    let mass_values = file
        .variable("mass_values")
        .unwrap_or_else(|| panic!("Can't find mass_values in variables"))
        .values::<f64>(None, None)
        .expect("Can't read values")
        .into_raw_vec();

    for (i, indices) in scans.iter().enumerate() {
        let mz = subset(&mass_values, indices);
        let intensities = subset(&intensity_values, indices);

        let peaks: Vec<Peak> = mz
            .iter()
            .zip(intensities.iter())
            .map(|(&a, &b)| Peak {
                mz: a,
                intensity: b,
            })
            .collect();

        run.scans.push(Scan {
            num: (i + 1) as u64,
            ms_level: 1,
            peaks_count: indices.len() as u64,
            retention_time: scan_time[i],
            peaks,
        });
    }
    Ok(run)
}
