use std::fs::{self, File};
use std::io::Write;
use std::process::Command;

fn write_minimal_pdb(path: &std::path::Path) {
    let mut file = File::create(path).unwrap();
    let pdb = concat!(
        "ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00 20.00           N\n",
        "ATOM      2  C   ALA A   1      11.504  13.607   2.500  1.00 20.00           C\n",
        "ATOM      3  O   ALA A   1      11.904  14.007   2.900  1.00 20.00           O\n",
        "ATOM      4  N   ALA A   2      12.304  14.407   3.300  1.00 20.00           N\n",
        "ATOM      5  C   ALA A   2      12.704  14.807   3.700  1.00 20.00           C\n",
        "ATOM      6  O   ALA A   2      13.104  15.207   4.100  1.00 20.00           O\n",
        "ATOM      7  N   ALA A   3      13.504  15.607   4.500  1.00 20.00           N\n",
        "ATOM      8  C   ALA A   3      13.904  16.007   4.900  1.00 20.00           C\n",
        "ATOM      9  O   ALA A   3      14.304  16.407   5.300  1.00 20.00           O\n",
    );
    file.write_all(pdb.as_bytes()).unwrap();
}

#[test]
fn cli_generates_outputs() {
    let temp_dir = std::env::temp_dir().join("errat_test_jobs");
    let _ = fs::remove_dir_all(&temp_dir);
    fs::create_dir_all(&temp_dir).unwrap();

    let job_id = "job1";
    let job_dir = temp_dir.join(job_id);
    fs::create_dir_all(&job_dir).unwrap();
    let pdb_path = job_dir.join("errat.pdb");
    write_minimal_pdb(&pdb_path);

    let exe = env!("CARGO_BIN_EXE_errat");
    let status = Command::new(exe)
        .arg("testfile")
        .arg(job_id)
        .env("ERRAT_JOBS_PATH", &temp_dir)
        .status()
        .expect("failed to run errat binary");
    assert!(status.success());

    let log_path = job_dir.join("errat.logf");
    let ps_path = job_dir.join("errat.ps");
    assert!(log_path.exists());
    assert!(ps_path.exists());

    let log_meta = fs::metadata(log_path).unwrap();
    assert!(log_meta.len() > 0);
}
