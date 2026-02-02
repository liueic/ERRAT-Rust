# ERRAT-Rust

[English](#english) | [中文](#中文)

---

## English

### Overview
A Rust reimplementation of the classic ERRAT algorithm for protein structure validation. The Rust version reproduces the original logic and output formatting, and generates identical `.logf` and `.ps` outputs for the same input PDB.

### Project structure
```
ERRAT-Rust/
  Cargo.toml
  Cargo.lock
  README.md
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

### Build
```bash
cargo build --release
```

### Run
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

### Environment variable
- `ERRAT_JOBS_PATH`: base directory containing job folders. Default: `/var/www/Jobs/`.

### Tests
```bash
cargo test
```

### Scripts
- `scripts/run_sample.sh`: run the bundled sample PDB and write to `outputs/<job>`
- `scripts/bench.sh [runs] [protein_id] [job_cpp] [job_rs]`: run benchmarks (C++ if available, Rust always)
- `scripts/compare_outputs.sh [job_cpp] [job_rs]`: byte-wise compare `errat.logf` and `errat.ps`

### Reproducibility
This Rust version matches the original C++ output byte-for-byte for `errat.logf` and `errat.ps` when the same input PDB is used.

---

## 中文

### 简介
这是对经典 ERRAT 蛋白结构验证算法的 Rust 版本重写。该实现严格保持原始逻辑与输出格式，确保同一输入 PDB 下生成的 `.logf` 与 `.ps` 完全一致。

### 项目结构
```
ERRAT-Rust/
  Cargo.toml
  Cargo.lock
  README.md
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

### 编译
```bash
cargo build --release
```

### 运行
程序会在作业目录中读取名为 `errat.pdb` 的输入文件。通过 `ERRAT_JOBS_PATH` 指定作业目录的根路径。

```bash
# 准备作业目录
mkdir -p /path/to/jobs/my_job
cp /path/to/your.pdb /path/to/jobs/my_job/errat.pdb

# 运行
ERRAT_JOBS_PATH=/path/to/jobs \
  /path/to/ERRAT-Rust/target/release/errat <ProteinID> my_job
```

输出文件：
- `<job>/errat.logf`
- `<job>/errat.ps`

### 环境变量
- `ERRAT_JOBS_PATH`：作业目录根路径，默认 `/var/www/Jobs/`。

### 测试
```bash
cargo test
```

### 脚本
- `scripts/run_sample.sh`：使用内置示例 PDB 运行并写入 `outputs/<job>`
- `scripts/bench.sh [runs] [protein_id] [job_cpp] [job_rs]`：基准测试（如果存在 C++ 二进制则同时测试）
- `scripts/compare_outputs.sh [job_cpp] [job_rs]`：对 `errat.logf` 和 `errat.ps` 进行字节级比较

### 一致性
在相同输入 PDB 下，本 Rust 版本生成的 `errat.logf` 与 `errat.ps` 与原 C++ 版本保持字节级一致。
