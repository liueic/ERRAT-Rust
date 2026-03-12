use std::path::{Path, PathBuf};

use crate::model::{ErratStats, LMT_95, LMT_99};

#[derive(Clone, Debug)]
pub struct RunOutput {
    pub logf: PathBuf,
    pub plot: PathBuf,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameStatus {
    Ok,
    Warning95,
    Warning99,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FrameScore {
    pub chain_id: String,
    pub center_residue: i32,
    pub error_value: f64,
    pub status: FrameStatus,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainSummary {
    pub chain_id: String,
    pub start_residue: i32,
    pub end_residue: i32,
}

#[derive(Clone, Debug)]
pub struct AnalysisResult {
    pub protein_id: String,
    pub input_path: PathBuf,
    pub scored_frame_count: usize,
    pub rejected_frame_count: usize,
    pub rejected_frame_ratio: Option<f64>,
    pub overall_quality_factor: Option<f64>,
    pub average_probability: Option<f64>,
    pub below_interaction_limit_frames: Vec<i32>,
    pub chain_summaries: Vec<ChainSummary>,
    pub frame_scores: Vec<FrameScore>,
    pub messages: Vec<String>,
    pub log_text: String,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ChainRange {
    pub(crate) chain_id: u8,
    pub(crate) start_residue: i32,
    pub(crate) end_residue: i32,
}

pub(crate) fn derive_file_string(input_pdb: &Path, protein_id: Option<&str>) -> String {
    protein_id
        .filter(|id| !id.trim().is_empty())
        .map(ToOwned::to_owned)
        .or_else(|| {
            input_pdb
                .file_stem()
                .and_then(|s| s.to_str())
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
        })
        .unwrap_or_else(|| "errat".to_string())
}

#[cfg(feature = "python")]
pub(crate) fn frame_status_name(status: FrameStatus) -> &'static str {
    match status {
        FrameStatus::Ok => "ok",
        FrameStatus::Warning95 => "warning95",
        FrameStatus::Warning99 => "warning99",
    }
}

fn chain_label(chain_id: u8) -> String {
    match chain_id {
        b' ' => String::new(),
        id => (id as char).to_string(),
    }
}

fn classify_frame(error_value: f64) -> FrameStatus {
    if error_value > LMT_99 {
        FrameStatus::Warning99
    } else if error_value > LMT_95 {
        FrameStatus::Warning95
    } else {
        FrameStatus::Ok
    }
}

pub(crate) fn compute_chain_ranges(stats: &ErratStats) -> Vec<ChainRange> {
    if stats.atmnum == 0 {
        return Vec::new();
    }

    let mut starts = [0i32; 100];
    let mut ends = [0i32; 100];
    let mut ids = [b' '; 100];
    let mut count = 1usize;

    starts[count] = stats.resnum[1] + 4;
    ids[count] = stats.chain_id[1];

    for idx in 1..stats.atmnum {
        if idx == stats.atmnum - 1 {
            ends[count] = stats.resnum[stats.atmnum] - 4;
        } else if stats.chain_id[idx] != stats.chain_id[idx + 1] && stats.resnum[idx] > 4 {
            ends[count] = stats.resnum[idx] - 4;
            count += 1;
            starts[count] = stats.resnum[idx + 1] + 4;
            ids[count] = stats.chain_id[idx + 1];
        }
    }

    (1..=count)
        .filter_map(|idx| {
            let start_residue = starts[idx];
            let end_residue = ends[idx];
            (end_residue >= start_residue).then_some(ChainRange {
                chain_id: ids[idx],
                start_residue,
                end_residue,
            })
        })
        .collect()
}

pub(crate) fn build_analysis_result(
    input_path: PathBuf,
    protein_id: String,
    stats: &ErratStats,
    log_text: String,
) -> AnalysisResult {
    let chain_ranges = compute_chain_ranges(stats);
    let chain_summaries = chain_ranges
        .iter()
        .map(|range| ChainSummary {
            chain_id: chain_label(range.chain_id),
            start_residue: range.start_residue,
            end_residue: range.end_residue,
        })
        .collect();

    let frame_scores = stats
        .scored_frames
        .iter()
        .map(|frame| {
            let chain_id = chain_ranges
                .iter()
                .find(|range| {
                    frame.center_residue >= range.start_residue
                        && frame.center_residue <= range.end_residue
                })
                .map(|range| chain_label(range.chain_id))
                .unwrap_or_default();
            FrameScore {
                chain_id,
                center_residue: frame.center_residue,
                error_value: frame.error_value,
                status: classify_frame(frame.error_value),
            }
        })
        .collect();

    let messages = log_text
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| line.to_string())
        .collect();

    AnalysisResult {
        protein_id,
        input_path,
        scored_frame_count: stats.stat as usize,
        rejected_frame_count: stats.pstat as usize,
        rejected_frame_ratio: if stats.stat > 0.0 {
            Some(stats.pstat / stats.stat)
        } else {
            None
        },
        overall_quality_factor: stats.overall_quality_factor,
        average_probability: stats.avg_probability,
        below_interaction_limit_frames: stats.warning_frames.clone(),
        chain_summaries,
        frame_scores,
        messages,
        log_text,
    }
}
