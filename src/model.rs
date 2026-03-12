use std::path::PathBuf;

pub(crate) const SIZE: usize = 250_000;
pub(crate) const BXMX: usize = 200_000;
pub(crate) const CHAINDIF: i32 = 10_000;
pub(crate) const BOXSIZE: f64 = 4.0;
pub(crate) const RADIUS: f64 = 3.75;
pub(crate) const RADMIN: f64 = 3.25;
pub(crate) const MAXWIN: f64 = 100.694;
pub(crate) const LMT_95: f64 = 11.526_684_477_428_809;
pub(crate) const LMT_99: f64 = 17.190_823_041_860_433;

#[derive(Clone, Debug)]
pub(crate) struct AtomData {
    pub(crate) atmnum: usize,
    pub(crate) name: Vec<i32>,
    pub(crate) bnam: Vec<i32>,
    pub(crate) chain_id: Vec<u8>,
    pub(crate) res_seq: Vec<i32>,
    pub(crate) resnum: Vec<i32>,
    pub(crate) xyz_x: Vec<f64>,
    pub(crate) xyz_y: Vec<f64>,
    pub(crate) xyz_z: Vec<f64>,
    pub(crate) errat: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct FrameScoreRaw {
    pub(crate) center_residue: i32,
    pub(crate) error_value: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct ErratStats {
    pub(crate) stat: f64,
    pub(crate) pstat: f64,
    pub(crate) avg_probability: Option<f64>,
    pub(crate) overall_quality_factor: Option<f64>,
    pub(crate) errat: Vec<f64>,
    pub(crate) resnum: Vec<i32>,
    pub(crate) chain_id: Vec<u8>,
    pub(crate) atmnum: usize,
    pub(crate) warning_frames: Vec<i32>,
    pub(crate) scored_frames: Vec<FrameScoreRaw>,
}

#[derive(Clone, Debug)]
pub(crate) struct Paths {
    pub(crate) pdb: PathBuf,
    pub(crate) logf: PathBuf,
    pub(crate) plot: PathBuf,
}
