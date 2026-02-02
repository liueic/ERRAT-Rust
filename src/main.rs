use std::env;
use std::path::PathBuf;

fn print_usage() {
    eprintln!(
        "\nUsage:\n  errat <ProteinID> <JobID>\n  errat --input <pdb> --out-dir <dir> --protein-id <id>\n\nEnvironment:\n  ERRAT_JOBS_PATH   base directory for job folders (default: ./outputs)\n"
    );
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
            _ => {}
        }
        i += 1;
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
