use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!(
            "\n2 arguments required: er34 'string for file name phrase' JobID\npath is internally set: /var/www/Jobs/\n"
        );
        std::process::exit(1);
    }

    let file_string = args[1].clone();
    let job_id = args[2].clone();
    let base_path = errat::default_base_path();

    let config = errat::Config {
        file_string,
        job_id,
        base_path,
    };

    if let Err(err) = errat::run(config) {
        eprintln!("ERRAT failed: {}", err);
        std::process::exit(1);
    }
}
