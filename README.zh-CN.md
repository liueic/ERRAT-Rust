# ERRAT-Rust

[English](README.md)

## 简介
这是对经典 ERRAT 蛋白结构验证算法的 Rust 版本重写。该实现严格保持原始逻辑与输出格式，确保同一输入 PDB 下生成的 `.logf` 与 `.ps` 完全一致。使用 `--pdf` 可直接输出 PDF。

## 项目结构
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

## 编译
```bash
cargo build --release
```

## Python wheel
这个仓库现在也可以构建成 Python wheel，发行名为 `errat-rs`，导入名为 `errat_rs`。

在仓库根目录安装：

```bash
python3 -m pip install .
```

结构化分析示例：

```python
import errat_rs

result = errat_rs.analyze(
    "/path/to/input.pdb",
    protein_id="MyProtein",
)
```

导出经典 ERRAT 报告文件：

```python
report = errat_rs.write_report(
    "/path/to/input.pdb",
    "/path/to/output",
    protein_id="MyProtein",
    output_format="pdf",
)
```

说明：
- `analyze()` 返回结构化 Python dataclass，而不是只落地文件。
- `write_report()` 保留原来的 `.logf` 和 `.ps` / `.pdf` 输出能力。
- 如果两者都要，可以在一次调用里使用 `analyze_and_write()`。

## 发布到 PyPI
仓库里已经加入了专门的 GitHub Actions 工作流 `.github/workflows/pypi.yml`。
它会构建 sdist 和多平台 wheel，并通过 Trusted Publishing 发布到 PyPI / TestPyPI。

第一次发布前，你需要先做这几步：
- 注册 PyPI 和 TestPyPI 账号。
- 在 GitHub 仓库 Settings -> Environments 里创建 `pypi` 和 `testpypi` 两个 environment。
- 给 `pypi` environment 配置 required reviewers，或者至少开启逐次人工批准。
- 在 PyPI 项目设置里为 `errat-rs` 添加一个 pending Trusted Publisher。
- GitHub owner 填 `liueic`，repository 填 `ERRAT-Rust`，workflow filename 填 `pypi.yml`，environment 填 `pypi`。
- 在 TestPyPI 里重复一遍，environment 改成 `testpypi`。

推荐发布流程：
```bash
# 1. 同时更新 Cargo.toml 和 pyproject.toml 的版本号

# 2. 本地验证
cargo test
cargo check --features python
python3 -m pip install . --force-reinstall

# 3. 推送代码
git push origin main

# 4. 可选：在 GitHub Actions 里手动触发 workflow_dispatch，先发布到 TestPyPI

# 5. 正式发布到 PyPI
git tag v0.1.0
git push origin v0.1.0
```

这个工作流会强制检查：
- `Cargo.toml` 和 `pyproject.toml` 的版本号必须一致。
- Git tag `vX.Y.Z` 必须和包版本 `X.Y.Z` 一致。

TestPyPI 发布后，可以这样安装验证：
```bash
python3 -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple errat-rs
```

## 运行
支持作业目录模式与直接文件模式。

### 作业目录模式（兼容原始布局）
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
- `<job>/errat.ps`（或使用 `--pdf` 输出 `<job>/errat.pdf`）

### 作业目录批处理模式
并行处理多个作业目录。包含 `errat.pdb` 的子目录会被视为一个作业。

```bash
errat --jobs-dir /path/to/jobs --threads 8
```

### 直接文件模式（CLI 工具）
```bash
errat --input /path/to/input.pdb --out-dir /path/to/output --protein-id <ProteinID>
```

最简示例（省略 `--protein-id`）：
```bash
errat --input /path/to/input.pdb --out-dir /path/to/output
```

输出文件：
- `<out-dir>/<input-stem>.logf`
- `<out-dir>/<input-stem>.ps`（或使用 `--pdf` 输出 `.pdf`）

说明：
- `--input` 支持 `.pdb`、`.cif`、`.mmcif`。
- 若省略 `--protein-id`，默认使用输入文件名（去掉扩展名）。

### 直接文件批处理模式
批量处理目录下所有 `.pdb`、`.cif`、`.mmcif` 文件。

```bash
errat --input-dir /path/to/pdbs --out-dir /path/to/output --threads 8
```

如需扫描子目录，请加 `--recursive`。

### 可选内存映射（仅 PDB）
使用 `--mmap` 通过内存映射读取 PDB。

```bash
errat --input /path/to/input.pdb --out-dir /path/to/output --mmap
```

### 可选：直接输出 PDF
使用 `--pdf` 直接输出 PDF，避免额外的转换步骤。

```bash
errat --input /path/to/input.pdb --out-dir /path/to/output --pdf
```

## 环境变量
- `ERRAT_JOBS_PATH`：作业目录根路径，默认 `./outputs`。

## 测试
```bash
cargo test
```

## 批量 PDF vs PS+ps2pdf 基准测试（示例）
下面是可复现实验，用于比较“直接输出 PDF（`--pdf`）”与“先输出 PS 再用 `ps2pdf` 转换”的批量性能。
会复制同一个 PDB 为 100/1000 个文件，并测试单线程和多线程。

说明：
- 测试对操作系统缓存非常敏感。“冷缓存”需要使用管理员权限清缓存，或重启后测试。
- 这里 `ps2pdf` 是顺序执行，因此在高批量时它通常是瓶颈。

```bash
# 编译
cargo build --release

# 路径
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

# 预热（可选）
"$BIN" --input "$PDB" --out-dir "$OUTD" --pdf >/dev/null
"$BIN" --input "$PDB" --out-dir "$OUTP" >/dev/null
ps2pdf "$OUTP"/*.ps "$OUTP"/*.pdf >/dev/null 2>/dev/null || true

# 100 个文件
run_direct "$IN100" 1 "100"
run_ps_then_pdf "$IN100" 1 "100"
run_direct "$IN100" 8 "100"
run_ps_then_pdf "$IN100" 8 "100"

# 1000 个文件
run_direct "$IN1000" 1 "1000"
run_ps_then_pdf "$IN1000" 1 "1000"
run_direct "$IN1000" 8 "1000"
run_ps_then_pdf "$IN1000" 8 "1000"
```

## 脚本
- `scripts/run_sample.sh`：使用内置示例 PDB 运行并写入 `outputs/<job>`
- `scripts/bench.sh [runs] [protein_id] [job_cpp] [job_rs]`：基准测试（如果存在 C++ 二进制则同时测试）
- `scripts/compare_outputs.sh [job_cpp] [job_rs]`：对 `errat.logf` 和 `errat.ps` 进行字节级比较

## 一致性
在相同输入 PDB 下，本 Rust 版本生成的 `errat.logf` 与 `errat.ps` 与原 C++ 版本保持字节级一致（PDF 输出为独立路径）。

## 引用
- Colovos C, Yeates TO (1993). Verification of protein structures: patterns of nonbonded atomic interactions.
- PubMed: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=8401235&dopt=Abstract

## 致谢
- ERRAT 原始 C++ 源码: https://saves.mbi.ucla.edu/ERRAT/errat_src.tar.gz
