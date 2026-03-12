mod api;
mod model;
mod parser;
mod pipeline;
#[cfg(feature = "python")]
mod python;
mod render;
mod stats;

pub use api::{AnalysisResult, ChainSummary, FrameScore, FrameStatus, RunOutput};

use std::io;
use std::path::{Path, PathBuf};

use api::{build_analysis_result, derive_file_string};
use pipeline::{persist_outputs, process_structure_data, resolve_paths};

#[cfg(feature = "python")]
use pyo3::prelude::*;

#[derive(Clone, Debug)]
pub struct Config {
    pub file_string: String,
    pub job_id: String,
    pub base_path: PathBuf,
    pub input_pdb: Option<PathBuf>,
    pub output_dir: Option<PathBuf>,
    pub use_mmap: bool,
    pub output_pdf: bool,
}

pub fn default_base_path() -> PathBuf {
    if let Ok(val) = std::env::var("ERRAT_JOBS_PATH") {
        PathBuf::from(val)
    } else {
        std::env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join("outputs")
    }
}

pub fn analyze_file<P: AsRef<Path>>(
    input_pdb: P,
    protein_id: Option<&str>,
    use_mmap: bool,
) -> io::Result<AnalysisResult> {
    let input_path = input_pdb.as_ref().to_path_buf();
    let protein_id = derive_file_string(&input_path, protein_id);
    let (stats, log) = process_structure_data(&input_path, use_mmap)?;
    let log_text = String::from_utf8_lossy(&log).into_owned();
    Ok(build_analysis_result(
        input_path, protein_id, &stats, log_text,
    ))
}

pub fn analyze_and_write<P: AsRef<Path>, Q: AsRef<Path>>(
    input_pdb: P,
    output_dir: Q,
    protein_id: Option<&str>,
    use_mmap: bool,
    output_pdf: bool,
) -> io::Result<(AnalysisResult, RunOutput)> {
    let input_path = input_pdb.as_ref().to_path_buf();
    let protein_id = derive_file_string(&input_path, protein_id);
    let config = Config {
        file_string: protein_id.clone(),
        job_id: "cli".to_string(),
        base_path: default_base_path(),
        input_pdb: Some(input_path.clone()),
        output_dir: Some(output_dir.as_ref().to_path_buf()),
        use_mmap,
        output_pdf,
    };
    let paths = resolve_paths(&config);
    let (stats, mut log) = process_structure_data(&input_path, use_mmap)?;
    let log_text = String::from_utf8_lossy(&log).into_owned();
    let analysis = build_analysis_result(input_path, protein_id, &stats, log_text);
    persist_outputs(
        &paths,
        &config.file_string,
        &stats,
        &mut log,
        config.output_pdf,
    )?;
    Ok((
        analysis,
        RunOutput {
            logf: paths.logf,
            plot: paths.plot,
        },
    ))
}

pub fn write_report<P: AsRef<Path>, Q: AsRef<Path>>(
    input_pdb: P,
    output_dir: Q,
    protein_id: Option<&str>,
    use_mmap: bool,
    output_pdf: bool,
) -> io::Result<RunOutput> {
    let config = Config {
        file_string: derive_file_string(input_pdb.as_ref(), protein_id),
        job_id: "cli".to_string(),
        base_path: default_base_path(),
        input_pdb: Some(input_pdb.as_ref().to_path_buf()),
        output_dir: Some(output_dir.as_ref().to_path_buf()),
        use_mmap,
        output_pdf,
    };
    let paths = resolve_paths(&config);
    let (stats, mut log) = process_structure_data(&paths.pdb, config.use_mmap)?;
    persist_outputs(
        &paths,
        &config.file_string,
        &stats,
        &mut log,
        config.output_pdf,
    )?;
    Ok(RunOutput {
        logf: paths.logf,
        plot: paths.plot,
    })
}

pub fn run_file<P: AsRef<Path>, Q: AsRef<Path>>(
    input_pdb: P,
    output_dir: Q,
    protein_id: Option<&str>,
    use_mmap: bool,
    output_pdf: bool,
) -> io::Result<RunOutput> {
    write_report(input_pdb, output_dir, protein_id, use_mmap, output_pdf)
}

pub fn run(config: Config) -> io::Result<()> {
    let paths = resolve_paths(&config);
    let (stats, mut log) = process_structure_data(&paths.pdb, config.use_mmap)?;
    persist_outputs(
        &paths,
        &config.file_string,
        &stats,
        &mut log,
        config.output_pdf,
    )
}

#[cfg(feature = "python")]
#[pymodule]
#[pyo3(name = "_native")]
fn errat_rs_native(_py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    python::register(module)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn minimal_pdb() -> &'static str {
        concat!(
            "ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00 20.00           N\n",
            "ATOM      2  C   ALA A   1      11.504  13.607   2.500  1.00 20.00           C\n",
            "ATOM      3  O   ALA A   1      11.904  14.007   2.900  1.00 20.00           O\n",
            "ATOM      4  N   ALA A   2      12.304  14.407   3.300  1.00 20.00           N\n",
            "ATOM      5  C   ALA A   2      12.704  14.807   3.700  1.00 20.00           C\n",
            "ATOM      6  O   ALA A   2      13.104  15.207   4.100  1.00 20.00           O\n",
            "ATOM      7  N   ALA A   3      13.504  15.607   4.500  1.00 20.00           N\n",
            "ATOM      8  C   ALA A   3      13.904  16.007   4.900  1.00 20.00           C\n",
            "ATOM      9  O   ALA A   3      14.304  16.407   5.300  1.00 20.00           O\n",
        )
    }

    fn temp_test_dir(prefix: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!("{prefix}_{unique}"));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn analyze_file_returns_structured_result() {
        let temp_dir = temp_test_dir("errat_rs_analysis");
        let pdb_path = temp_dir.join("sample.pdb");
        fs::write(&pdb_path, minimal_pdb()).unwrap();

        let result = analyze_file(&pdb_path, Some("sample"), false).unwrap();

        assert_eq!(result.protein_id, "sample");
        assert_eq!(result.input_path, pdb_path);
        assert!(!result.messages.is_empty());

        let _ = fs::remove_dir_all(&temp_dir);
    }

    #[test]
    fn analyze_and_write_creates_pdf_report() {
        let temp_dir = temp_test_dir("errat_rs_report");
        let input_path = temp_dir.join("sample.pdb");
        let output_dir = temp_dir.join("out");
        fs::write(&input_path, minimal_pdb()).unwrap();

        let (analysis, report) =
            analyze_and_write(&input_path, &output_dir, Some("sample"), false, true).unwrap();

        assert_eq!(analysis.protein_id, "sample");
        assert!(report.logf.exists());
        assert!(report.plot.exists());
        assert_eq!(
            report.plot.extension().and_then(|ext| ext.to_str()),
            Some("pdf")
        );

        let _ = fs::remove_dir_all(&temp_dir);
    }
}
