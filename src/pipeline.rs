use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::Config;
use crate::model::{ErratStats, Paths};
use crate::parser::parse_structure;
use crate::render::{write_pdf, write_ps};
use crate::stats::compute_errat;

pub(crate) fn resolve_paths(config: &Config) -> Paths {
    if let (Some(input_pdb), Some(output_dir)) = (&config.input_pdb, &config.output_dir) {
        let base_name = input_pdb
            .file_stem()
            .and_then(|s| s.to_str())
            .filter(|s| !s.is_empty())
            .unwrap_or("errat");
        let mut logf = output_dir.clone();
        logf.push(format!("{base_name}.logf"));
        let mut plot = output_dir.clone();
        if config.output_pdf {
            plot.push(format!("{base_name}.pdf"));
        } else {
            plot.push(format!("{base_name}.ps"));
        }
        return Paths {
            pdb: input_pdb.clone(),
            logf,
            plot,
        };
    }

    let mut base = config.base_path.clone();
    base.push(&config.job_id);

    let mut pdb = base.clone();
    pdb.push("errat.pdb");

    let mut logf = base.clone();
    logf.push("errat.logf");

    let mut plot = base.clone();
    if config.output_pdf {
        plot.push("errat.pdf");
    } else {
        plot.push("errat.ps");
    }

    Paths { pdb, logf, plot }
}

pub(crate) fn process_structure_data(
    path: &std::path::PathBuf,
    use_mmap: bool,
) -> io::Result<(ErratStats, Vec<u8>)> {
    let mut log = Vec::new();
    let atom_data = parse_structure(path, &mut log, use_mmap)?;
    let stats = compute_errat(&atom_data, &mut log)?;
    Ok((stats, log))
}

pub(crate) fn persist_outputs(
    paths: &Paths,
    file_string: &str,
    stats: &ErratStats,
    log: &mut Vec<u8>,
    output_pdf: bool,
) -> io::Result<()> {
    if let Some(parent) = paths.logf.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let plotf = File::create(&paths.plot)?;
    let mut plotw = BufWriter::new(plotf);
    if stats.stat > 0.0 {
        if output_pdf {
            write_pdf(&mut plotw, log, file_string, stats)?;
        } else {
            write_ps(&mut plotw, log, file_string, stats)?;
        }
    }
    plotw.flush()?;

    let logf = File::create(&paths.logf)?;
    let mut logw = BufWriter::new(logf);
    logw.write_all(log)?;
    logw.flush()?;
    Ok(())
}
