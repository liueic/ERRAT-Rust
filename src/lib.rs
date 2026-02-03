use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use rayon::prelude::*;

const SIZE: usize = 250_000;
const BXMX: usize = 200_000;
const CHAINDIF: i32 = 10_000;
const BOXSIZE: f64 = 4.0;
const RADIUS: f64 = 3.75;
const RADMIN: f64 = 3.25;
const MAXWIN: f64 = 100.694;
const LMT_95: f64 = 11.526_684_477_428_809;
const LMT_99: f64 = 17.190_823_041_860_433;

#[derive(Clone, Debug)]
pub struct Config {
    pub file_string: String,
    pub job_id: String,
    pub base_path: PathBuf,
    pub input_pdb: Option<PathBuf>,
    pub output_dir: Option<PathBuf>,
}

#[derive(Clone, Debug)]
struct AtomData {
    atmnum: usize,
    name: Vec<i32>,
    bnam: Vec<i32>,
    chain_id: Vec<u8>,
    res_seq: Vec<i32>,
    resnum: Vec<i32>,
    xyz_x: Vec<f64>,
    xyz_y: Vec<f64>,
    xyz_z: Vec<f64>,
    errat: Vec<f64>,
}

#[derive(Clone, Debug)]
struct ErratStats {
    stat: f64,
    pstat: f64,
    errat: Vec<f64>,
    resnum: Vec<i32>,
    chain_id: Vec<u8>,
    atmnum: usize,
}

#[derive(Clone, Debug)]
struct Paths {
    pdb: PathBuf,
    logf: PathBuf,
    ps: PathBuf,
}

pub fn default_base_path() -> PathBuf {
    if let Ok(val) = std::env::var("ERRAT_JOBS_PATH") {
        PathBuf::from(val)
    } else {
        std::env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join("outputs")
    }
}

pub fn run(config: Config) -> io::Result<()> {
    let paths = resolve_paths(&config);
    if let Some(parent) = paths.logf.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let logf = File::create(&paths.logf)?;
    let mut logw = BufWriter::new(logf);

    let psf = File::create(&paths.ps)?;
    let mut psw = BufWriter::new(psf);

    let pdbf = File::open(&paths.pdb)?;
    let mut pdb_reader = BufReader::new(pdbf);

    let atom_data = parse_structure(&paths.pdb, &mut pdb_reader, &mut logw)?;
    let stats = compute_errat(&atom_data, &mut logw)?;

    if stats.stat > 0.0 {
        write_ps(&mut psw, &mut logw, &config.file_string, &stats)?;
    }

    logw.flush()?;
    psw.flush()?;
    Ok(())
}

fn parse_structure<R: BufRead, W: Write>(
    path: &PathBuf,
    reader: &mut R,
    logw: &mut W,
) -> io::Result<AtomData> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if ext == "cif" || ext == "mmcif" {
        parse_mmcif(reader, logw)
    } else {
        parse_pdb(reader, logw)
    }
}

fn resolve_paths(config: &Config) -> Paths {
    if let (Some(input_pdb), Some(output_dir)) = (&config.input_pdb, &config.output_dir) {
        let base_name = input_pdb
            .file_stem()
            .and_then(|s| s.to_str())
            .filter(|s| !s.is_empty())
            .unwrap_or("errat");
        let mut logf = output_dir.clone();
        logf.push(format!("{base_name}.logf"));
        let mut ps = output_dir.clone();
        ps.push(format!("{base_name}.ps"));
        return Paths {
            pdb: input_pdb.clone(),
            logf,
            ps,
        };
    }

    let mut base = config.base_path.clone();
    base.push(&config.job_id);

    let mut pdb = base.clone();
    pdb.push("errat.pdb");

    let mut logf = base.clone();
    logf.push("errat.logf");

    let mut ps = base.clone();
    ps.push("errat.ps");

    Paths { pdb, logf, ps }
}

fn parse_pdb<R: BufRead, W: Write>(reader: &mut R, logw: &mut W) -> io::Result<AtomData> {
    let mut name = vec![0i32; SIZE + 2];
    let mut bnam = vec![0i32; SIZE + 2];
    let mut chain_id = vec![b' '; SIZE + 2];
    let mut res_seq = vec![0i32; SIZE + 2];
    let mut resnum = vec![0i32; SIZE + 2];
    let mut xyz_x = vec![0.0f64; SIZE + 2];
    let mut xyz_y = vec![0.0f64; SIZE + 2];
    let mut xyz_z = vec![0.0f64; SIZE + 2];
    let mut errat = vec![0.0f64; SIZE + 8];

    let mut i: usize = 0;
    let mut atmnum: usize = 0;
    let mut kadd: i32 = 0;
    let mut flag = false;
    let mut flag2 = false;
    let mut line = String::new();
    while !flag2 {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }
        let bytes = line.as_bytes();
        if bytes.len() < 6 {
            continue;
        }
        if &bytes[..6] != b"ATOM  " {
            continue;
        }
        if i + 1 > SIZE - 1 {
            writeln!(
                logw,
                "ERROR: PDB WITH TOO MANY ATOMS. CUT OFF FURTHER INPUT."
            )?;
            break;
        }
        i += 1;
        if bytes.len() < 54 {
            i -= 1;
            continue;
        }

        let name_temp = bytes[13];
        name[i] = match name_temp {
            b'C' => 1,
            b'N' => 2,
            b'O' => 3,
            _ => 0,
        };

        if bytes.len() < 16 {
            i -= 1;
            continue;
        }
        let name_temp2 = &bytes[13..16];
        bnam[i] = if name_temp2 == b"N  " || name_temp2 == b"C  " {
            1
        } else {
            0
        };

        let alt_loc = bytes[16] as char;
        let res_name = &bytes[17..20];
        chain_id[i] = bytes[21];

        let res_seq_temp = std::str::from_utf8(&bytes[22..26]).unwrap_or("");
        let res_seq_val = res_seq_temp.trim().parse::<f64>().unwrap_or(0.0);
        res_seq[i] = res_seq_val as i32;

        let x_temp = std::str::from_utf8(&bytes[30..38]).unwrap_or("");
        let y_temp = std::str::from_utf8(&bytes[38..46]).unwrap_or("");
        let z_temp = std::str::from_utf8(&bytes[46..54]).unwrap_or("");
        xyz_x[i] = x_temp.trim().parse::<f64>().unwrap_or(0.0);
        xyz_y[i] = y_temp.trim().parse::<f64>().unwrap_or(0.0);
        xyz_z[i] = z_temp.trim().parse::<f64>().unwrap_or(0.0);

        if !(alt_loc == ' ' || alt_loc == 'A' || alt_loc == 'a' || alt_loc == 'P') {
            writeln!(
                logw,
                "Reject 2' Conformation atom#\t{}\tchain\t{}",
                i,
                chain_id[i] as char
            )?;
            i -= 1;
            flag = true;
        }

        if !is_standard_residue(res_name) {
            i -= 1;
            flag = true;
            let res_name_str = std::str::from_utf8(res_name).unwrap_or("???");
            writeln!(
                logw,
                "***Warning: Reject Nonstardard Residue - {}",
                res_name_str
            )?;
        }

        if i >= 2 && !flag && chain_id[i] != chain_id[i - 1] {
            kadd += 1;
            writeln!(logw, "INCREMENTING CHAIN (kadd) {}", kadd)?;
        }

        if !flag {
            resnum[i] = res_seq[i] + (kadd * CHAINDIF);
            atmnum = i;
        }

        if i >= 2
            && !flag
            && chain_id[i] == chain_id[i - 1]
            && resnum[i] < resnum[i - 1]
        {
            writeln!(
                logw,
                "ERROR: RESNUM DECREASE. TERMINATE ANALYSIS{}\t{}",
                resnum[i], resnum[i - 1]
            )?;
            flag2 = true;
        }

        if i > 2
            && !flag
            && chain_id[i] == chain_id[i - 1]
            && resnum[i] != resnum[i - 1]
            && (resnum[i] - resnum[i - 1]) > 1
        {
            writeln!(
                logw,
                "WARNING: Missing Residues{}>>>{}",
                resnum[i - 1], resnum[i]
            )?;
        }

        if !flag {
            let idx = (resnum[i] + 4) as usize;
            if idx >= errat.len() {
                errat.resize(idx + 1, 0.0);
            }
            errat[idx] = 0.0;
        }

        flag = false;
    }

    Ok(AtomData {
        atmnum,
        name,
        bnam,
        chain_id,
        res_seq,
        resnum,
        xyz_x,
        xyz_y,
        xyz_z,
        errat,
    })
}

fn parse_mmcif<R: Read, W: Write>(reader: &mut R, logw: &mut W) -> io::Result<AtomData> {
    let mut name = vec![0i32; SIZE + 2];
    let mut bnam = vec![0i32; SIZE + 2];
    let mut chain_id = vec![b' '; SIZE + 2];
    let mut res_seq = vec![0i32; SIZE + 2];
    let mut resnum = vec![0i32; SIZE + 2];
    let mut xyz_x = vec![0.0f64; SIZE + 2];
    let mut xyz_y = vec![0.0f64; SIZE + 2];
    let mut xyz_z = vec![0.0f64; SIZE + 2];
    let mut errat = vec![0.0f64; SIZE + 8];

    let mut input = String::new();
    reader.read_to_string(&mut input)?;
    let tokens = tokenize_cif(&input);

    let mut i: usize = 0;
    let mut atmnum: usize = 0;
    let mut kadd: i32 = 0;
    let mut flag = false;
    let mut flag2 = false;

    let mut idx = 0;
    while idx < tokens.len() {
        if tokens[idx] != "loop_" {
            idx += 1;
            continue;
        }
        idx += 1;
        let mut cols = Vec::new();
        while idx < tokens.len() && tokens[idx].starts_with('_') {
            cols.push(tokens[idx].clone());
            idx += 1;
        }
        if cols.is_empty() {
            continue;
        }

        let is_atom_site = cols.iter().any(|c| c.starts_with("_atom_site."));
        let col_count = cols.len();

        if !is_atom_site {
            while idx + col_count <= tokens.len() {
                let t = &tokens[idx];
                if t == "loop_"
                    || t.starts_with('_')
                    || t.starts_with("data_")
                    || t.starts_with("save_")
                    || t == "stop_"
                {
                    break;
                }
                idx += col_count;
            }
            continue;
        }

        let col_index = |name: &str| -> Option<usize> {
            cols.iter().position(|c| {
                if c == name {
                    true
                } else if name.starts_with("_atom_site.") {
                    false
                } else {
                    c.ends_with(&format!(".{name}"))
                }
            })
        };

        let idx_group = col_index("group_PDB");
        let idx_atom = col_index("label_atom_id");
        let idx_type = col_index("type_symbol");
        let idx_alt = col_index("label_alt_id");
        let idx_res = col_index("label_comp_id");
        let idx_chain = col_index("auth_asym_id").or_else(|| col_index("label_asym_id"));
        let idx_seq = col_index("auth_seq_id").or_else(|| col_index("label_seq_id"));
        let idx_x = col_index("Cartn_x");
        let idx_y = col_index("Cartn_y");
        let idx_z = col_index("Cartn_z");

        if idx_atom.is_none() || idx_res.is_none() || idx_chain.is_none() || idx_seq.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "mmCIF missing required _atom_site columns",
            ));
        }
        if idx_x.is_none() || idx_y.is_none() || idx_z.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "mmCIF missing coordinate columns",
            ));
        }

        while idx + col_count <= tokens.len() {
            let t = &tokens[idx];
            if t == "loop_"
                || t.starts_with('_')
                || t.starts_with("data_")
                || t.starts_with("save_")
                || t == "stop_"
            {
                break;
            }

            let row = &tokens[idx..idx + col_count];
            idx += col_count;

            if let Some(g) = idx_group {
                let group = row[g].as_str();
                if group != "ATOM" {
                    continue;
                }
            }

            if i + 1 > SIZE - 1 {
                writeln!(
                    logw,
                    "ERROR: PDB WITH TOO MANY ATOMS. CUT OFF FURTHER INPUT."
                )?;
                break;
            }
            i += 1;

            let atom_name = row[idx_atom.unwrap()].as_str();
            let element = idx_type
                .and_then(|k| row.get(k))
                .map(|s| s.as_str())
                .unwrap_or(atom_name);
            let element_char = element.chars().next().unwrap_or(' ');
            name[i] = match element_char {
                'C' | 'c' => 1,
                'N' | 'n' => 2,
                'O' | 'o' => 3,
                _ => 0,
            };
            bnam[i] = if atom_name == "N" || atom_name == "C" { 1 } else { 0 };

            let alt_loc = idx_alt
                .and_then(|k| row.get(k))
                .map(|s| s.as_str())
                .unwrap_or(".");
            let alt_loc_char = alt_loc.chars().next().unwrap_or(' ');
            let alt_loc_char = match alt_loc_char {
                '.' | '?' => ' ',
                c => c,
            };

            let res_name_str = row[idx_res.unwrap()].as_str();
            let res_name_upper = res_name_str.to_ascii_uppercase();
            let res_name = res_name_upper.as_bytes();
            let chain = row[idx_chain.unwrap()].as_bytes();
            chain_id[i] = if chain.is_empty() { b' ' } else { chain[0] };

            let res_seq_val = row[idx_seq.unwrap()].parse::<f64>().unwrap_or(0.0);
            res_seq[i] = res_seq_val as i32;

            xyz_x[i] = row[idx_x.unwrap()].parse::<f64>().unwrap_or(0.0);
            xyz_y[i] = row[idx_y.unwrap()].parse::<f64>().unwrap_or(0.0);
            xyz_z[i] = row[idx_z.unwrap()].parse::<f64>().unwrap_or(0.0);

            if !(alt_loc_char == ' ' || alt_loc_char == 'A' || alt_loc_char == 'a' || alt_loc_char == 'P') {
                writeln!(
                    logw,
                    "Reject 2' Conformation atom#\t{}\tchain\t{}",
                    i,
                    chain_id[i] as char
                )?;
                i -= 1;
                flag = true;
            }

            if !is_standard_residue(res_name) {
                i -= 1;
                flag = true;
                let res_name_str = std::str::from_utf8(res_name).unwrap_or("???");
                writeln!(
                    logw,
                    "***Warning: Reject Nonstardard Residue - {}",
                    res_name_str
                )?;
            }

            if i >= 2 && !flag && chain_id[i] != chain_id[i - 1] {
                kadd += 1;
                writeln!(logw, "INCREMENTING CHAIN (kadd) {}", kadd)?;
            }

            if !flag {
                resnum[i] = res_seq[i] + (kadd * CHAINDIF);
                atmnum = i;
            }

            if i >= 2
                && !flag
                && chain_id[i] == chain_id[i - 1]
                && resnum[i] < resnum[i - 1]
            {
                writeln!(
                    logw,
                    "ERROR: RESNUM DECREASE. TERMINATE ANALYSIS{}\t{}",
                    resnum[i], resnum[i - 1]
                )?;
                flag2 = true;
            }

            if i > 2
                && !flag
                && chain_id[i] == chain_id[i - 1]
                && resnum[i] != resnum[i - 1]
                && (resnum[i] - resnum[i - 1]) > 1
            {
                writeln!(
                    logw,
                    "WARNING: Missing Residues{}>>>{}",
                    resnum[i - 1], resnum[i]
                )?;
            }

            if !flag {
                let idx = (resnum[i] + 4) as usize;
                if idx >= errat.len() {
                    errat.resize(idx + 1, 0.0);
                }
                errat[idx] = 0.0;
            }

            flag = false;
            if flag2 {
                break;
            }
        }

        if atmnum > 0 || flag2 {
            break;
        }
    }

    Ok(AtomData {
        atmnum,
        name,
        bnam,
        chain_id,
        res_seq,
        resnum,
        xyz_x,
        xyz_y,
        xyz_z,
        errat,
    })
}

fn tokenize_cif(input: &str) -> Vec<String> {
    let bytes = input.as_bytes();
    let mut tokens = Vec::new();
    let mut i = 0;
    let mut at_line_start = true;

    while i < bytes.len() {
        let c = bytes[i] as char;
        if c.is_whitespace() {
            if c == '\n' {
                at_line_start = true;
            }
            i += 1;
            continue;
        }

        if c == '#' {
            while i < bytes.len() && bytes[i] as char != '\n' {
                i += 1;
            }
            at_line_start = true;
            continue;
        }

        if c == ';' && at_line_start {
            i += 1;
            let start = i;
            let mut end = None;
            while i + 1 < bytes.len() {
                if bytes[i] as char == '\n' && bytes[i + 1] as char == ';' {
                    end = Some(i);
                    break;
                }
                i += 1;
            }
            if let Some(end_pos) = end {
                let val = &input[start..end_pos];
                tokens.push(val.to_string());
                i = end_pos + 2;
                while i < bytes.len() && bytes[i] as char != '\n' {
                    i += 1;
                }
                at_line_start = true;
                continue;
            }
        }

        if c == '\'' || c == '"' {
            let quote = c;
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i] as char != quote {
                i += 1;
            }
            let val = &input[start..i];
            tokens.push(val.to_string());
            i += 1;
            at_line_start = false;
            continue;
        }

        let start = i;
        while i < bytes.len() && !(bytes[i] as char).is_whitespace() {
            i += 1;
        }
        let val = &input[start..i];
        tokens.push(val.to_string());
        at_line_start = false;
    }

    tokens
}

#[derive(Clone, Copy)]
enum WindowOutcome {
    Warn(i32),
    Value { idx: usize, mtrx: f64 },
}

fn compute_window(
    i: usize,
    data: &AtomData,
    min: &[f64; 4],
    nbx: &[i32; 4],
    ibox_counts: &[i32],
    ibox_atoms: &[i32],
    box_slots: usize,
    rsq: f64,
    ssq: f64,
    ndelta: i32,
) -> Option<WindowOutcome> {
    let mut s = 1;
    let mut v = i;
    while s < 10 && v <= data.atmnum {
        let diff = data.resnum[v + 1] - data.resnum[v];
        if ((diff < 100) && (diff > 0)) || v == data.atmnum {
            s += 1;
        }
        v += 1;
    }
    if v > 0 {
        v -= 1;
    }

    if s != 10 || data.res_seq[v] <= data.res_seq[i] {
        return None;
    }

    let mut c = [[0.0f64; 4]; 4];
    for rer in i..=v {
        let jbx = ((data.xyz_x[rer] - (min[1] - 0.00001)) / BOXSIZE).floor() as i32;
        let jby = ((data.xyz_y[rer] - (min[2] - 0.00001)) / BOXSIZE).floor() as i32;
        let jbz = ((data.xyz_z[rer] - (min[3] - 0.00001)) / BOXSIZE).floor() as i32;

        let mut ibz1 = jbz - ndelta;
        if ibz1 < 0 {
            ibz1 = 0;
        }
        let mut ibz2 = jbz + ndelta;
        if ibz2 > nbx[3] - 1 {
            ibz2 = nbx[3] - 1;
        }

        let mut iby1 = jby - ndelta;
        if iby1 < 0 {
            iby1 = 0;
        }
        let mut iby2 = jby + ndelta;
        if iby2 > nbx[2] - 1 {
            iby2 = nbx[2] - 1;
        }

        let mut ibx1 = jbx - ndelta;
        if ibx1 < 0 {
            ibx1 = 0;
        }
        let mut ibx2 = jbx + ndelta;
        if ibx2 > nbx[1] - 1 {
            ibx2 = nbx[1] - 1;
        }

        let rer_x = data.xyz_x[rer];
        let rer_y = data.xyz_y[rer];
        let rer_z = data.xyz_z[rer];

        for j in ibz1..=ibz2 {
            for k in iby1..=iby2 {
                for l in ibx1..=ibx2 {
                    let ind = (1 + l + k * nbx[1] + j * nbx[1] * nbx[2]) as usize;
                    let count = ibox_counts[ind] as usize;
                    let limit = count.min(box_slots);
                    let base = ind * box_slots;
                    for m in 0..limit {
                        let n = ibox_atoms[base + m] as usize;

                        if data.resnum[rer] != data.resnum[n] {
                            let dx = data.xyz_x[n] - rer_x;
                            let dy = data.xyz_y[n] - rer_y;
                            let dz = data.xyz_z[n] - rer_z;
                            let dsq = dx * dx + dy * dy + dz * dz;
                            if dsq < rsq {
                                if data.bnam[rer] == 1
                                    && data.bnam[n] == 1
                                    && (((data.resnum[n] == data.resnum[rer] + 1)
                                        && (data.name[rer] == 1)
                                        && (data.name[n] == 2))
                                        || ((data.resnum[rer] == data.resnum[n] + 1)
                                            && (data.name[rer] == 2)
                                            && (data.name[n] == 1)))
                                {
                                    // skip backbone neighbor contact
                                } else if n >= i && n <= v {
                                    if data.resnum[rer] > data.resnum[n] {
                                        let temp1 = if dsq <= ssq {
                                            1.0
                                        } else {
                                            2.0 * (RADIUS - dsq.sqrt())
                                        };
                                        c[data.name[rer] as usize][data.name[n] as usize] +=
                                            temp1;
                                    }
                                } else {
                                    let temp1 = if dsq <= ssq {
                                        1.0
                                    } else {
                                        2.0 * (RADIUS - dsq.sqrt())
                                    };
                                    c[data.name[rer] as usize][data.name[n] as usize] += temp1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let mut temp2 = 0.0f64;
    for q in 1..=3 {
        for r in 1..=3 {
            temp2 += c[q][r];
        }
    }

    if temp2 > MAXWIN {
        let mut matrix = [0.0f64; 6];
        matrix[1] = c[1][1] / temp2;
        matrix[2] = (c[1][2] + c[2][1]) / temp2;
        matrix[3] = (c[1][3] + c[3][1]) / temp2;
        matrix[4] = c[2][2] / temp2;
        matrix[5] = (c[2][3] + c[3][2]) / temp2;

        let mtrx = matrixdb(&matrix);
        let idx = (data.resnum[i] + 4) as usize;
        Some(WindowOutcome::Value { idx, mtrx })
    } else {
        Some(WindowOutcome::Warn(data.resnum[i] + 4))
    }
}

fn compute_errat<W: Write>(data: &AtomData, logw: &mut W) -> io::Result<ErratStats> {
    let mut min = [0.0f64; 4];
    let mut max = [0.0f64; 4];
    for i in 1..=3 {
        min[i] = 999.0;
        max[i] = -999.0;
    }

    if data.atmnum == 0 {
        return Ok(ErratStats {
            stat: 0.0,
            pstat: 0.0,
            errat: data.errat.clone(),
            resnum: data.resnum.clone(),
            chain_id: data.chain_id.clone(),
            atmnum: data.atmnum,
        });
    }

    for i in 1..=data.atmnum {
        let vx = data.xyz_x[i];
        let vy = data.xyz_y[i];
        let vz = data.xyz_z[i];
        if vx < min[1] {
            min[1] = vx;
        }
        if vx > max[1] {
            max[1] = vx;
        }
        if vy < min[2] {
            min[2] = vy;
        }
        if vy > max[2] {
            max[2] = vy;
        }
        if vz < min[3] {
            min[3] = vz;
        }
        if vz > max[3] {
            max[3] = vz;
        }
    }

    for j in 1..=3 {
        write!(logw, "{}\t", fmt_sig6(min[j]))?;
    }
    for j in 1..=3 {
        write!(logw, "{}\t", fmt_sig6(max[j]))?;
    }
    writeln!(logw)?;

    let mut nbx = [0i32; 4];
    for i in 1..=3 {
        nbx[i] = ((max[i] - min[i]) / BOXSIZE) as i32 + 1;
    }

    let box_count = (nbx[1] * nbx[2] * nbx[3]) as i64;
    let mut flag2 = false;
    if box_count > (BXMX as i64 - 1) {
        writeln!(logw, "ERROR: TOO MANY BOXES")?;
        flag2 = true;
    }

    let box_slots = 15usize;
    let ibox_len = (box_count.max(0) as usize) + 1;
    let mut ibox_counts = vec![0i32; ibox_len];
    let mut ibox_atoms = vec![0i32; ibox_len * box_slots];

    if !flag2 {
        for i in 1..=data.atmnum {
            let ix = ((data.xyz_x[i] - (min[1] - 0.00001)) / BOXSIZE).floor() as i32;
            let iy = ((data.xyz_y[i] - (min[2] - 0.00001)) / BOXSIZE).floor() as i32;
            let iz = ((data.xyz_z[i] - (min[3] - 0.00001)) / BOXSIZE).floor() as i32;
            let ind = (1 + ix + iy * nbx[1] + iz * nbx[1] * nbx[2]) as usize;

            let temp = ibox_counts[ind] as usize;
            ibox_counts[ind] += 1;
            if temp < box_slots {
                let base = ind * box_slots;
                ibox_atoms[base + temp] = i as i32;
            }
        }

        for i in 1..ibox_counts.len() {
            if ibox_counts[i] > 15 {
                writeln!(logw, "TOO MANY ATOMS IN BOX #:\t{}", ibox_counts[i])?;
                flag2 = true;
            }
        }
    }

    let mut stat = 0.0f64;
    let mut pstat = 0.0f64;
    let mut mtrxstat = 0.0f64;
    let mut errat = data.errat.clone();

    if !flag2 {
        let rsq = RADIUS * RADIUS;
        let ssq = RADMIN * RADMIN;
        let ndelta = (RADIUS / BOXSIZE).ceil() as i32;
        let window_starts: Vec<usize> = (1..=data.atmnum)
            .filter(|&i| i == 1 || data.resnum[i] > data.resnum[i - 1])
            .collect();

        let results: Vec<Option<WindowOutcome>> = window_starts
            .par_iter()
            .map(|&i| {
                compute_window(
                    i,
                    data,
                    &min,
                    &nbx,
                    &ibox_counts,
                    &ibox_atoms,
                    box_slots,
                    rsq,
                    ssq,
                    ndelta,
                )
            })
            .collect();

        for outcome in results {
            if let Some(outcome) = outcome {
                match outcome {
                    WindowOutcome::Warn(frame) => {
                        writeln!(
                            logw,
                            "WARNING: Frame\t{}\tBelow Minimum Interaction Limit.",
                            frame
                        )?;
                    }
                    WindowOutcome::Value { idx, mtrx } => {
                        stat += 1.0;
                        mtrxstat += mtrx;

                        if mtrx > LMT_99 {
                            pstat += 1.0;
                        } else if mtrx > LMT_95 {
                            pstat += 1.0;
                        }

                        if idx >= errat.len() {
                            errat.resize(idx + 1, 0.0);
                        }
                        errat[idx] = mtrx;
                    }
                }
            }
        }
    }

    if stat > 0.0 {
        writeln!(
            logw,
            "Total frames: {}\tP frames {}\tNumber: {}",
            stat as i64,
            pstat as i64,
            fmt_sig6(pstat / stat)
        )?;
        writeln!(logw)?;
        writeln!(logw, "Avg Probability\t{}", fmt_sig6(mtrxstat / stat))?;
        writeln!(
            logw,
            "# Overall quality factor: {}",
            fmt_sig6(100.0 - (100.0 * pstat / stat))
        )?;
    }

    Ok(ErratStats {
        stat,
        pstat,
        errat,
        resnum: data.resnum.clone(),
        chain_id: data.chain_id.clone(),
        atmnum: data.atmnum,
    })
}

fn write_ps<P: Write, L: Write>(
    psw: &mut P,
    logw: &mut L,
    file_string: &str,
    stats: &ErratStats,
) -> io::Result<()> {
    let mut ir1 = [0i32; 100];
    let mut ir2 = [0i32; 100];
    let mut id_by_chain = [b' '; 100];

    let chainx = 1 + (stats.resnum[stats.atmnum] - 4) / CHAINDIF;

    let mut z2 = 1;
    ir1[z2] = stats.resnum[1] + 4;
    ir2[z2] = 0;
    id_by_chain[z2] = stats.chain_id[1];
    println!(
        "atn, chain#, chainID 1  {}  {}",
        z2,
        id_by_chain[z2] as char
    );

    for z1 in 1..stats.atmnum {
        if z1 == stats.atmnum - 1 {
            ir2[z2] = stats.resnum[stats.atmnum] - 4;
        } else if stats.chain_id[z1] != stats.chain_id[z1 + 1] && stats.resnum[z1] > 4 {
            ir2[z2] = stats.resnum[z1] - 4;
            z2 += 1;
            ir1[z2] = stats.resnum[z1 + 1] + 4;
            id_by_chain[z2] = stats.chain_id[z1 + 1];
        }
    }

    let mut mst = 0.0f64;
    for ich in 1..=chainx as usize {
        let mut ms = (ir2[ich] - ir1[ich] + 1) as f64 / (300.0 + 1.0);
        ms = (ir2[ich] - ir1[ich] + 1) as f64 / ms;
        if ms > mst {
            mst = ms;
        }
        if mst < 200.0 {
            mst = 200.0;
        }
    }

    let sz = 200.0 / mst;

    for ich in 1..=chainx as usize {
        let np = 1 + ((ir2[ich] - ir1[ich] + 1) as f64 / mst) as i32;
        for z1 in 1..=np {
            let ir0 = ir1[ich] + (mst as i32) * (z1 - 1);
            let mut ir = ir0 + (mst as i32) - 1;
            if ir > ir2[ich] {
                ir = ir2[ich];
            }

            let overall_quality = 100.0 - (100.0 * stats.pstat / stats.stat);

            writeln!(
                logw,
                "# Chain Label {}:    Residue range {} to {}",
                id_by_chain[ich] as char,
                ir0,
                ir
            )?;

            writeln!(psw, "%!PS")?;
            writeln!(psw, "%FIXED")?;
            writeln!(psw, "/sce {{8}} def /scr {{3}} def")?;
            writeln!(
                psw,
                "90 rotate 110 -380 translate /e95 {{11.527}} def /e99 {{17.191}} def"
            )?;
            writeln!(psw, "/Helvetica findfont 18 scalefont setfont 0.5 setlinewidth")?;
            writeln!(psw, "/bar1 {{/g {{1 1 1}} def bar}} def /bar2 {{/g {{1 1 0}} def bar}} def")?;
            writeln!(psw, "/bar3 {{/g {{1 0 0}} def bar}} def /bar {{sce mul /yval exch def")?;
            writeln!(psw, " scr mul /xval exch def")?;
            writeln!(psw, "newpath xval 0 moveto xval yval lineto scr -1 mul 0")?;
            writeln!(psw, " rlineto 0 yval -1 mul rlineto closepath gsave g setrgbcolor")?;
            writeln!(psw, " fill grestore stroke}} def")?;
            writeln!(psw, "/tick {{newpath 0.5 sub scr mul 0 moveto 0 -3 rlineto")?;
            writeln!(psw, " currentpoint stroke moveto -10 -12 rmoveto}} def")?;

            writeln!(psw, "% VARIABLE")?;
            writeln!(
                psw,
                "{:.3}   {:.3} scale /rlim {{{}}} def",
                sz,
                sz,
                ir - ir0 + 1
            )?;
            writeln!(psw, "gsave 0 30 sce mul 20 add translate ")?;
            writeln!(
                psw,
                "0 30 moveto (Chain#:{}) show ",
                id_by_chain[ich] as char
            )?;
            writeln!(psw, "0 50 moveto (File: {}) show ", file_string)?;
            writeln!(
                psw,
                "0 10 moveto (Overall quality factor**: {:.3})show",
                overall_quality
            )?;
            writeln!(psw, "0 70 moveto (Program: ERRAT2) show")?;
            writeln!(psw, "() show")?;

            writeln!(psw, "% FIXED")?;
            writeln!(psw, "grestore newpath 0 0 moveto 0 27 sce mul rlineto stroke")?;
            writeln!(psw, "newpath rlim scr mul 0 moveto 0 27 sce mul rlineto stroke")?;
            writeln!(psw, "newpath 0  0 moveto rlim scr mul 0 rlineto stroke")?;
            writeln!(psw, "newpath -3 e95 sce mul moveto rlim scr mul 3 add 0 rlineto")?;
            writeln!(psw, "stroke newpath -3 e99 sce mul moveto rlim scr mul 3 add 0")?;
            writeln!(psw, " rlineto stroke")?;
            writeln!(psw, "newpath 0  27  sce mul moveto rlim scr")?;
            writeln!(psw, " mul 0 rlineto stroke")?;
            writeln!(psw, "rlim scr mul 2 div 100 sub -34")?;
            writeln!(psw, " moveto (Residue # (window center)) show")?;
            writeln!(psw, "/Helvetica findfont 14 scalefont setfont 0.5 setlinewidth")?;
            writeln!(psw, "-34 e95 sce mul 4 sub moveto (95\\%) show")?;
            writeln!(psw, "-34 e99 sce mul 4 sub moveto (99\\%) show")?;
            writeln!(psw, "/Helvetica findfont 12 scalefont setfont 0.5 setlinewidth")?;
            writeln!(psw, "0 -70 moveto (*On the error axis, two lines are drawn to indicate the confidence with) show")?;
            writeln!(psw, "0 -82 moveto (which it is possible to reject regions that exceed that error value.) show")?;
            writeln!(psw, "0 -100 moveto (**Expressed as the percentage of the protein for which the calculated) show")?;
            writeln!(psw, "0 -112 moveto (error value falls below the 95\\% rejection limit.  Good high resolution) show")?;
            writeln!(psw, "0 -124 moveto (structures generally produce values around 95\\% or higher.  For lower) show")?;
            writeln!(psw, "0 -136 moveto (resolutions (2.5 to 3A) the average overall quality factor is around 91\\%. ) show")?;
            writeln!(psw, "/Helvetica findfont 18 scalefont setfont 0.5 setlinewidth")?;
            writeln!(psw, "gsave -40 -5 translate 90 rotate 80 0 moveto (Error value*)")?;
            writeln!(psw, "show grestore")?;
            writeln!(psw, "/Helvetica findfont 16 scalefont setfont 0.5 setlinewidth")?;

            for z2 in ir0..=ir {
                if z2 % 20 == 0 {
                    writeln!(psw, "{} tick        ", z2 - ir0 + 1)?;
                    writeln!(
                        psw,
                        "({}) show\t",
                        z2 - (CHAINDIF * (z2 / CHAINDIF))
                    )?;
                } else if z2 % 10 == 0 {
                    writeln!(psw, "{} tick\t", z2 - ir0 + 1)?;
                }
            }

            for z2 in ir0..=ir {
                let mut bar = "bar1";
                if stats.errat[z2 as usize] > LMT_95 {
                    bar = "bar2";
                }
                if stats.errat[z2 as usize] > LMT_99 {
                    bar = "bar3";
                }
                let mut val = stats.errat[z2 as usize];
                if val > 27.0 {
                    val = 27.0;
                }
                writeln!(psw, "{}\t{:.3} {}", z2 - ir0 + 1, val, bar)?;
            }
            writeln!(psw, "showpage")?;
        }
    }

    Ok(())
}

fn is_standard_residue(res_name: &[u8]) -> bool {
    matches!(
        res_name,
        b"GLY" | b"ALA" | b"VAL" | b"LEU" | b"ILE" | b"TYR" | b"CYS" | b"MET" | b"TRP"
            | b"PHE" | b"HIS" | b"PRO" | b"SER" | b"THR" | b"LYS" | b"ARG" | b"GLU"
            | b"ASP" | b"GLN" | b"ASN"
    )
}

fn matrixdb(matrix: &[f64; 6]) -> f64 {
    let b1: [[f64; 6]; 6] = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [
            0.0,
            5040.279078850848,
            3408.8051415836494,
            4152.904423767301,
            4236.20000417189,
            5054.7812102046255,
        ],
        [
            0.0,
            3408.805141583649,
            8491.906094010221,
            5958.88177787795,
            1521.3873527184862,
            4304.078200827222,
        ],
        [
            0.0,
            4152.9044237673015,
            5958.881777877952,
            7637.16708933505,
            6620.7157382230725,
            5287.691183798411,
        ],
        [
            0.0,
            4236.20000417189,
            1521.3873527184862,
            6620.7157382230725,
            18368.34377429841,
            4050.7978111188067,
        ],
        [
            0.0,
            5054.7812102046255,
            4304.078200827221,
            5287.69118379841,
            4050.7978111188067,
            6666.856740479165,
        ],
    ];

    let avg = [
        0.0,
        0.192765509919262,
        0.195575208778518,
        0.275322406824210,
        0.059102357035642,
        0.233154192767480,
    ];

    let mut v = [0.0f64; 6];
    for u in 1..6 {
        v[u] = matrix[u] - avg[u];
    }

    let mut c = [0.0f64; 6];
    for j in 1..6 {
        let mut x = 0.0;
        for k in 1..6 {
            x += v[k] * b1[k][j];
        }
        c[j] = x;
    }

    let mut total = 0.0;
    for k in 1..6 {
        total += c[k] * v[k];
    }

    total
}

fn fmt_sig6(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    let abs = value.abs();
    let exp = abs.log10().floor() as i32;
    let decimals = if exp >= 5 { 0 } else { (5 - exp) as usize };
    let mut s = format!("{:.*}", decimals, value);
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    if s == "-0" {
        s = "0".to_string();
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn matrixdb_zero_at_avg() {
        let mut matrix = [0.0f64; 6];
        matrix[1] = 0.192765509919262;
        matrix[2] = 0.195575208778518;
        matrix[3] = 0.275322406824210;
        matrix[4] = 0.059102357035642;
        matrix[5] = 0.233154192767480;
        let out = matrixdb(&matrix);
        assert!(out.abs() < 1e-9, "expected near-zero, got {}", out);
    }

    #[test]
    fn parse_pdb_rejects_nonstandard_and_altloc() {
        let pdb = b"\
ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00 20.00           N\n\
ATOM      2  CA  MSE A   2      12.000  13.000   2.000  1.00 20.00           C\n\
ATOM      3  CA BALA A   3      13.000  13.000   2.000  1.00 20.00           C\n";
        let mut reader = Cursor::new(pdb.as_ref());
        let mut log = Vec::new();
        let data = parse_pdb(&mut reader, &mut log).unwrap();
        assert_eq!(data.atmnum, 1);
    }
}
