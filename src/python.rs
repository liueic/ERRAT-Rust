use std::io;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3::wrap_pyfunction;

use crate::api::frame_status_name;
use crate::{AnalysisResult, analyze_and_write, analyze_file, write_report};

fn io_err_to_py(err: io::Error) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

fn analysis_to_pydict(py: Python<'_>, analysis: &AnalysisResult) -> PyResult<Py<PyDict>> {
    let result = PyDict::new(py);
    result.set_item("protein_id", &analysis.protein_id)?;
    result.set_item("input_path", analysis.input_path.to_string_lossy().as_ref())?;
    result.set_item("scored_frame_count", analysis.scored_frame_count)?;
    result.set_item("rejected_frame_count", analysis.rejected_frame_count)?;
    result.set_item("rejected_frame_ratio", analysis.rejected_frame_ratio)?;
    result.set_item("overall_quality_factor", analysis.overall_quality_factor)?;
    result.set_item("average_probability", analysis.average_probability)?;
    result.set_item(
        "below_interaction_limit_frames",
        &analysis.below_interaction_limit_frames,
    )?;
    result.set_item("messages", &analysis.messages)?;
    result.set_item("log_text", &analysis.log_text)?;

    let chain_summaries = PyList::empty(py);
    for chain in &analysis.chain_summaries {
        let item = PyDict::new(py);
        item.set_item("chain_id", &chain.chain_id)?;
        item.set_item("start_residue", chain.start_residue)?;
        item.set_item("end_residue", chain.end_residue)?;
        chain_summaries.append(item)?;
    }
    result.set_item("chain_summaries", chain_summaries)?;

    let frame_scores = PyList::empty(py);
    for frame in &analysis.frame_scores {
        let item = PyDict::new(py);
        item.set_item("chain_id", &frame.chain_id)?;
        item.set_item("center_residue", frame.center_residue)?;
        item.set_item("error_value", frame.error_value)?;
        item.set_item("status", frame_status_name(frame.status))?;
        frame_scores.append(item)?;
    }
    result.set_item("frame_scores", frame_scores)?;

    Ok(result.unbind())
}

fn report_to_pydict(
    py: Python<'_>,
    log_path: &str,
    plot_path: &str,
    output_format: &str,
) -> PyResult<Py<PyDict>> {
    let result = PyDict::new(py);
    result.set_item("log_path", log_path)?;
    result.set_item("plot_path", plot_path)?;
    result.set_item("output_format", output_format)?;
    Ok(result.unbind())
}

#[pyfunction(name = "analyze", signature = (input_path, protein_id=None, use_mmap=false))]
fn analyze_py(
    py: Python<'_>,
    input_path: &str,
    protein_id: Option<&str>,
    use_mmap: bool,
) -> PyResult<Py<PyDict>> {
    let analysis = analyze_file(input_path, protein_id, use_mmap).map_err(io_err_to_py)?;
    analysis_to_pydict(py, &analysis)
}

#[pyfunction(name = "write_report", signature = (input_path, output_dir, protein_id=None, output_format="ps", use_mmap=false))]
fn write_report_py(
    py: Python<'_>,
    input_path: &str,
    output_dir: &str,
    protein_id: Option<&str>,
    output_format: &str,
    use_mmap: bool,
) -> PyResult<Py<PyDict>> {
    let normalized = output_format.to_ascii_lowercase();
    let output_pdf = match normalized.as_str() {
        "ps" => false,
        "pdf" => true,
        _ => {
            return Err(PyValueError::new_err(
                "output_format must be either 'ps' or 'pdf'",
            ));
        }
    };

    let outputs = write_report(input_path, output_dir, protein_id, use_mmap, output_pdf)
        .map_err(io_err_to_py)?;

    report_to_pydict(
        py,
        outputs.logf.to_string_lossy().as_ref(),
        outputs.plot.to_string_lossy().as_ref(),
        &normalized,
    )
}

#[pyfunction(name = "analyze_and_write", signature = (input_path, output_dir, protein_id=None, output_format="ps", use_mmap=false))]
fn analyze_and_write_py(
    py: Python<'_>,
    input_path: &str,
    output_dir: &str,
    protein_id: Option<&str>,
    output_format: &str,
    use_mmap: bool,
) -> PyResult<Py<PyDict>> {
    let normalized = output_format.to_ascii_lowercase();
    let output_pdf = match normalized.as_str() {
        "ps" => false,
        "pdf" => true,
        _ => {
            return Err(PyValueError::new_err(
                "output_format must be either 'ps' or 'pdf'",
            ));
        }
    };

    let (analysis, outputs) =
        analyze_and_write(input_path, output_dir, protein_id, use_mmap, output_pdf)
            .map_err(io_err_to_py)?;
    let payload = analysis_to_pydict(py, &analysis)?;
    let report_paths = report_to_pydict(
        py,
        outputs.logf.to_string_lossy().as_ref(),
        outputs.plot.to_string_lossy().as_ref(),
        &normalized,
    )?;
    payload
        .bind(py)
        .set_item("report_paths", report_paths.bind(py))?;
    Ok(payload)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(analyze_py, module)?)?;
    module.add_function(wrap_pyfunction!(analyze_and_write_py, module)?)?;
    module.add_function(wrap_pyfunction!(write_report_py, module)?)?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
