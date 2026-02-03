# ERRAT-Rust

[中文说明](README.zh-CN.md)

## Overview
A Rust reimplementation of the classic ERRAT algorithm for protein structure validation. The Rust version reproduces the original logic and output formatting, and generates identical `.logf` and `.ps` outputs for the same input PDB.

## Project structure
```
ERRAT-Rust/
  .github/
  .gitignore
  Cargo.toml
  Cargo.lock
  LICENSE
  README.md
  README.zh-CN.md
  src/
  tests/
```

## Build
```bash
cargo build --release
```

## Run
The program can run in job-folder mode or direct file mode.

### Job-folder mode (compatible with original layout)
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

### Batch job-folder mode
Process multiple job folders in parallel. Each subdirectory containing `errat.pdb` is treated as one job.

```bash
errat --jobs-dir /path/to/jobs --threads 8
```

### Direct file mode (CLI tool)
```bash
errat --input /path/to/input.pdb --out-dir /path/to/output --protein-id <ProteinID>
```

Minimal example (omit `--protein-id`):
```bash
errat --input /path/to/input.pdb --out-dir /path/to/output
```

Outputs:
- `<out-dir>/<input-stem>.logf`
- `<out-dir>/<input-stem>.ps`

Notes:
- `--input` supports `.pdb`, `.cif`, and `.mmcif`.
- If `--protein-id` is omitted, it defaults to the input filename without the extension.

### Batch direct file mode
Process all `.pdb`, `.cif`, and `.mmcif` files in a directory.

```bash
errat --input-dir /path/to/pdbs --out-dir /path/to/output --threads 8
```

Add `--recursive` to scan subdirectories.

### Optional memory mapping (PDB only)
Use `--mmap` to read PDB files via memory-mapped I/O.

```bash
errat --input /path/to/input.pdb --out-dir /path/to/output --mmap
```

## Environment variable
- `ERRAT_JOBS_PATH`: base directory containing job folders. Default: `./outputs`.

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
