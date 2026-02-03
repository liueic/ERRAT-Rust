use std::env;
use std::io;
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

fn print_usage() {
    eprintln!(
        "\nUsage:\n  errat <ProteinID> <JobID>\n  errat --input <pdb|cif> --out-dir <dir> [--protein-id <id>] [--mmap]\n  errat --input-dir <dir> --out-dir <dir> [--recursive] [--threads <n>] [--mmap]\n  errat --jobs-dir <dir> [--threads <n>] [--mmap]\n\nEnvironment:\n  ERRAT_JOBS_PATH   base directory for job folders (default: ./outputs)\n"
    );
}

struct BatchItem {
    label: String,
    config: errat::Config,
}

fn is_structure_file(path: &Path) -> bool {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    matches!(ext.as_str(), "pdb" | "cif" | "mmcif")
}

fn collect_inputs(dir: &Path, recursive: bool) -> io::Result<Vec<PathBuf>> {
    let mut inputs = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(next) = stack.pop() {
        for entry in std::fs::read_dir(&next)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                if recursive {
                    stack.push(path);
                }
            } else if is_structure_file(&path) {
                inputs.push(path);
            }
        }
    }
    inputs.sort();
    Ok(inputs)
}

fn run_batch(items: Vec<BatchItem>, threads: Option<usize>) -> io::Result<(usize, Vec<String>)> {
    let run_all = || {
        items
            .par_iter()
            .map(|item| (item.label.clone(), errat::run(item.config.clone())))
            .collect::<Vec<_>>()
    };

    let results = if let Some(threads) = threads {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
        pool.install(run_all)
    } else {
        run_all()
    };

    let mut success = 0usize;
    let mut errors = Vec::new();
    for (label, result) in results {
        match result {
            Ok(()) => success += 1,
            Err(err) => errors.push(format!("{label}: {err}")),
        }
    }
    Ok((success, errors))
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() == 1 || args.iter().any(|a| a == "-h" || a == "--help") {
        print_usage();
        return;
    }

    let mut input_pdb: Option<PathBuf> = None;
    let mut output_dir: Option<PathBuf> = None;
    let mut protein_id: Option<String> = None;
    let mut input_dir: Option<PathBuf> = None;
    let mut jobs_dir: Option<PathBuf> = None;
    let mut recursive = false;
    let mut threads: Option<usize> = None;
    let mut use_mmap = false;

    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "--input" => {
                i += 1;
                input_pdb = args.get(i).map(PathBuf::from);
            }
            "--out-dir" => {
                i += 1;
                output_dir = args.get(i).map(PathBuf::from);
            }
            "--protein-id" => {
                i += 1;
                protein_id = args.get(i).cloned();
            }
            "--input-dir" => {
                i += 1;
                input_dir = args.get(i).map(PathBuf::from);
            }
            "--jobs-dir" => {
                i += 1;
                jobs_dir = args.get(i).map(PathBuf::from);
            }
            "--recursive" => {
                recursive = true;
            }
            "--threads" => {
                i += 1;
                threads = args
                    .get(i)
                    .and_then(|v| v.parse::<usize>().ok())
                    .filter(|v| *v > 0);
            }
            "--mmap" => {
                use_mmap = true;
            }
            _ => {}
        }
        i += 1;
    }

    if protein_id.is_none() {
        if let Some(input_pdb) = input_pdb.as_ref() {
            if let Some(stem) = input_pdb.file_stem().and_then(|s| s.to_str()) {
                if !stem.is_empty() {
                    protein_id = Some(stem.to_string());
                }
            }
        }
    }

    if jobs_dir.is_some() && input_dir.is_some() {
        eprintln!("ERRAT failed: --jobs-dir and --input-dir cannot be used together.");
        std::process::exit(1);
    }

    if let Some(jobs_dir) = jobs_dir {
        let entries = match std::fs::read_dir(&jobs_dir) {
            Ok(entries) => entries,
            Err(err) => {
                eprintln!("ERRAT failed: {}", err);
                std::process::exit(1);
            }
        };

        let mut items = Vec::new();
        for entry in entries.flatten() {
            let path = entry.path();
            if !path.is_dir() {
                continue;
            }
            let job_id = entry.file_name().to_string_lossy().to_string();
            let mut pdb_path = path.clone();
            pdb_path.push("errat.pdb");
            if !pdb_path.exists() {
                continue;
            }
            items.push(BatchItem {
                label: job_id.clone(),
                config: errat::Config {
                    file_string: job_id.clone(),
                    job_id,
                    base_path: jobs_dir.clone(),
                    input_pdb: None,
                    output_dir: None,
                    use_mmap,
                },
            });
        }

        if items.is_empty() {
            eprintln!("ERRAT failed: no job folders with errat.pdb found.");
            std::process::exit(1);
        }

        match run_batch(items, threads) {
            Ok((success, errors)) => {
                for error in &errors {
                    eprintln!("ERRAT failed: {error}");
                }
                if !errors.is_empty() {
                    eprintln!(
                        "ERRAT batch completed with errors: {success} ok, {} failed.",
                        errors.len()
                    );
                    std::process::exit(1);
                }
            }
            Err(err) => {
                eprintln!("ERRAT failed: {err}");
                std::process::exit(1);
            }
        }
        return;
    }

    if let Some(input_dir) = input_dir {
        let output_dir = match output_dir {
            Some(dir) => dir,
            None => {
                eprintln!("ERRAT failed: --out-dir is required with --input-dir.");
                std::process::exit(1);
            }
        };

        if let Err(err) = std::fs::create_dir_all(&output_dir) {
            eprintln!("ERRAT failed: {}", err);
            std::process::exit(1);
        }

        let inputs = match collect_inputs(&input_dir, recursive) {
            Ok(inputs) => inputs,
            Err(err) => {
                eprintln!("ERRAT failed: {}", err);
                std::process::exit(1);
            }
        };

        if inputs.is_empty() {
            eprintln!("ERRAT failed: no input files found.");
            std::process::exit(1);
        }

        let items = inputs
            .into_iter()
            .filter_map(|input_pdb| {
                let stem = input_pdb
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .filter(|s| !s.is_empty())
                    .map(|s| s.to_string())?;
                Some(BatchItem {
                    label: stem.clone(),
                    config: errat::Config {
                        file_string: stem,
                        job_id: "cli".to_string(),
                        base_path: errat::default_base_path(),
                        input_pdb: Some(input_pdb),
                        output_dir: Some(output_dir.clone()),
                        use_mmap,
                    },
                })
            })
            .collect::<Vec<_>>();

        match run_batch(items, threads) {
            Ok((success, errors)) => {
                for error in &errors {
                    eprintln!("ERRAT failed: {error}");
                }
                if !errors.is_empty() {
                    eprintln!(
                        "ERRAT batch completed with errors: {success} ok, {} failed.",
                        errors.len()
                    );
                    std::process::exit(1);
                }
            }
            Err(err) => {
                eprintln!("ERRAT failed: {err}");
                std::process::exit(1);
            }
        }
        return;
    }

    let config = if let (Some(input_pdb), Some(output_dir), Some(protein_id)) =
        (input_pdb, output_dir, protein_id)
    {
        errat::Config {
            file_string: protein_id,
            job_id: "cli".to_string(),
            base_path: errat::default_base_path(),
            input_pdb: Some(input_pdb),
            output_dir: Some(output_dir),
            use_mmap,
        }
    } else if args.len() == 3 {
        let file_string = args[1].clone();
        let job_id = args[2].clone();
        errat::Config {
            file_string,
            job_id,
            base_path: errat::default_base_path(),
            input_pdb: None,
            output_dir: None,
            use_mmap,
        }
    } else {
        print_usage();
        std::process::exit(1);
    };

    if let Err(err) = errat::run(config) {
        eprintln!("ERRAT failed: {}", err);
        std::process::exit(1);
    }
}
