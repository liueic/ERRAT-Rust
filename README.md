# ERRAT-Rust

[中文说明](README.zh-CN.md)

## Overview
A Rust reimplementation of the classic ERRAT algorithm for protein structure validation. The Rust version reproduces the original logic and output formatting, and generates identical `.logf` and `.ps` outputs for the same input PDB.

## Project structure
```
ERRAT-Rust/
  Cargo.toml
  Cargo.lock
  README.md
  README.zh-CN.md
  LICENSE
  src/
  tests/
  examples/
  data/
    Hpyr004913.1-R65830.mRNA_relaxed_rank_001_alphafold2_ptm_model_5_seed_000.pdb
  outputs/
    job_run/
      errat.logf
      errat.ps
  scripts/
    run_sample.sh
    bench.sh
    compare_outputs.sh
```

## Build
```bash
cargo build --release
```

## Run
The program expects an input PDB named `errat.pdb` inside a job folder. Use `ERRAT_JOBS_PATH` to point to the base directory of job folders.

```bash
# prepare job folder
mkdir -p /path/to/jobs/my_job
cp /path/to/your.pdb /path/to/jobs/my_job/errat.pdb

# run
ERRAT_JOBS_PATH=/path/to/jobs \
  /path/to/ERRAT-Rust/target/release/errat <ProteinID> my_job
```

Outputs:
- `<job>/errat.logf`
- `<job>/errat.ps`

## Environment variable
- `ERRAT_JOBS_PATH`: base directory containing job folders. Default: `/var/www/Jobs/`.

## Tests
```bash
cargo test
```

## Scripts
- `scripts/run_sample.sh`: run the bundled sample PDB and write to `outputs/<job>`
- `scripts/bench.sh [runs] [protein_id] [job_cpp] [job_rs]`: run benchmarks (C++ if available, Rust always)
- `scripts/compare_outputs.sh [job_cpp] [job_rs]`: byte-wise compare `errat.logf` and `errat.ps`

## Reproducibility
This Rust version matches the original C++ output byte-for-byte for `errat.logf` and `errat.ps` when the same input PDB is used.

## References
- Colovos C, Yeates TO (1993). Verification of protein structures: patterns of nonbonded atomic interactions.
- PubMed: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8401235&dopt=Abstract

## Acknowledgements
- Original ERRAT C++ source: https://saves.mbi.ucla.edu/ERRAT/errat_src.tar.gz
