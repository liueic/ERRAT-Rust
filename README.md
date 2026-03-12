# ERRAT-Rust

[中文说明](README.zh-CN.md)

## Overview
A Rust reimplementation of the classic ERRAT algorithm for protein structure validation. The Rust version reproduces the original logic and output formatting, and generates identical `.logf` and `.ps` outputs for the same input PDB. Use `--pdf` to generate PDF directly.

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

## Python wheel
This crate can also be built as a Python wheel named `errat-rs` and imported with `import errat_rs`.

Install from the repository root:

```bash
python3 -m pip install .
```

Analyze a structure and get structured results:

```python
import errat_rs

result = errat_rs.analyze(
    "/path/to/input.pdb",
    protein_id="MyProtein",
)
```

Write the classic ERRAT report files:

```python
report = errat_rs.write_report(
    "/path/to/input.pdb",
    "/path/to/output",
    protein_id="MyProtein",
    output_format="pdf",
)
```

Notes:
- `analyze()` returns structured Python dataclasses instead of only writing files.
- `write_report()` keeps the original `.logf` and `.ps` / `.pdf` outputs.
- `analyze_and_write()` is available when you want both in one call.

## Publish to PyPI
This repository includes a dedicated GitHub Actions workflow at `.github/workflows/pypi.yml`.
It builds an sdist plus platform wheels and publishes them with Trusted Publishing.

Before the first release:
- Create accounts on PyPI and TestPyPI.
- In GitHub repository settings, create `pypi` and `testpypi` environments.
- Configure required reviewers or manual approval for the `pypi` environment.
- In the PyPI project settings, add a pending Trusted Publisher for `errat-rs`.
- Use GitHub owner `liueic`, repository `ERRAT-Rust`, workflow filename `pypi.yml`, and environment `pypi`.
- Repeat the same setup on TestPyPI with environment `testpypi`.

Recommended release flow:
```bash
# 1. Bump versions in Cargo.toml and pyproject.toml

# 2. Verify locally
cargo test
cargo check --features python
python3 -m pip install . --force-reinstall

# 3. Push changes
git push origin main

# 4. Optional: publish to TestPyPI from GitHub Actions workflow_dispatch

# 5. Publish to PyPI by pushing a version tag
git tag v0.1.0
git push origin v0.1.0
```

The workflow enforces that:
- `Cargo.toml` and `pyproject.toml` use the same version.
- Tag `vX.Y.Z` matches the package version `X.Y.Z`.

After TestPyPI release, install from TestPyPI with:
```bash
python3 -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple errat-rs
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
- `<job>/errat.ps` (or `<job>/errat.pdf` with `--pdf`)

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
- `<out-dir>/<input-stem>.ps` (or `.pdf` with `--pdf`)

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

### Optional PDF output
Use `--pdf` to write PDF directly instead of PostScript.

```bash
errat --input /path/to/input.pdb --out-dir /path/to/output --pdf
```

## Environment variable
- `ERRAT_JOBS_PATH`: base directory containing job folders. Default: `./outputs`.

## Tests
```bash
cargo test
```

## Batch PDF vs PS+ps2pdf benchmark (example)
This is a reproducible way to compare direct PDF output (`--pdf`) against
PostScript plus `ps2pdf` in high-batch scenarios. It creates 100/1000 copies
of one PDB and measures both single-thread and multi-thread runs.

Notes:
- These timings are sensitive to OS cache state. For a “cold cache” test you
  need elevated privileges to drop caches, or reboot between runs.
- `ps2pdf` is run sequentially here; its cost dominates at high batch sizes.

```bash
# Build
cargo build --release

# Paths
PDB="/path/to/your.pdb"
BIN="/path/to/ERRAT-Rust/target/release/errat"
BASE="/tmp/errat_pdf_batch"
IN100="$BASE/in_100"
IN1000="$BASE/in_1000"
OUTD="$BASE/out_direct"
OUTP="$BASE/out_ps"

mkdir -p "$IN100" "$IN1000" "$OUTD" "$OUTP"

make_copies() {
  local count="$1"; local dir="$2"
  rm -f "$dir"/*.pdb
  for i in $(seq -w 1 "$count"); do
    cp "$PDB" "$dir"/sample_"$i".pdb
  done
}

make_copies 100 "$IN100"
make_copies 1000 "$IN1000"

run_direct() {
  local dir="$1"; local threads="$2"; local tag="$3"
  rm -f "$OUTD"/*
  /usr/bin/time -f "direct_${tag}_t${threads} real=%e user=%U sys=%S" \
    "$BIN" --input-dir "$dir" --out-dir "$OUTD" --threads "$threads" --pdf >/dev/null
}

run_ps_then_pdf() {
  local dir="$1"; local threads="$2"; local tag="$3"
  rm -f "$OUTP"/*
  /usr/bin/time -f "ps_${tag}_t${threads} real=%e user=%U sys=%S" \
    "$BIN" --input-dir "$dir" --out-dir "$OUTP" --threads "$threads" >/dev/null
  /usr/bin/time -f "ps2pdf_${tag}_t${threads} real=%e user=%U sys=%S" \
    bash -lc 'for f in "$0"/*.ps; do ps2pdf "$f" "${f%.ps}.pdf" >/dev/null; done' "$OUTP"
}

# Warmup (optional)
"$BIN" --input "$PDB" --out-dir "$OUTD" --pdf >/dev/null
"$BIN" --input "$PDB" --out-dir "$OUTP" >/dev/null
ps2pdf "$OUTP"/*.ps "$OUTP"/*.pdf >/dev/null 2>/dev/null || true

# 100 files
run_direct "$IN100" 1 "100"
run_ps_then_pdf "$IN100" 1 "100"
run_direct "$IN100" 8 "100"
run_ps_then_pdf "$IN100" 8 "100"

# 1000 files
run_direct "$IN1000" 1 "1000"
run_ps_then_pdf "$IN1000" 1 "1000"
run_direct "$IN1000" 8 "1000"
run_ps_then_pdf "$IN1000" 8 "1000"
```

## Scripts
- `scripts/run_sample.sh`: run the bundled sample PDB and write to `outputs/<job>`
- `scripts/bench.sh [runs] [protein_id] [job_cpp] [job_rs]`: run benchmarks (C++ if available, Rust always)
- `scripts/compare_outputs.sh [job_cpp] [job_rs]`: byte-wise compare `errat.logf` and `errat.ps`

## Reproducibility
This Rust version matches the original C++ output byte-for-byte for `errat.logf` and `errat.ps` when the same input PDB is used (PDF output is a separate code path).

## References
- Colovos C, Yeates TO (1993). Verification of protein structures: patterns of nonbonded atomic interactions.
- PubMed: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8401235&dopt=Abstract

## Acknowledgements
- Original ERRAT C++ source: https://saves.mbi.ucla.edu/ERRAT/errat_src.tar.gz
