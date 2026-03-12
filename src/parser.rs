use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::PathBuf;

use memmap2::MmapOptions;

use crate::model::{AtomData, CHAINDIF, SIZE};

pub(crate) fn parse_structure<W: Write>(
    path: &PathBuf,
    logw: &mut W,
    use_mmap: bool,
) -> io::Result<AtomData> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if ext == "cif" || ext == "mmcif" {
        let pdbf = File::open(path)?;
        let mut reader = BufReader::new(pdbf);
        parse_mmcif(&mut reader, logw)
    } else if use_mmap {
        parse_pdb_mmap(path, logw)
    } else {
        let pdbf = File::open(path)?;
        let mut reader = BufReader::new(pdbf);
        parse_pdb(&mut reader, logw)
    }
}

fn empty_atom_data() -> AtomData {
    AtomData {
        atmnum: 0,
        name: vec![0i32; SIZE + 2],
        bnam: vec![0i32; SIZE + 2],
        chain_id: vec![b' '; SIZE + 2],
        res_seq: vec![0i32; SIZE + 2],
        resnum: vec![0i32; SIZE + 2],
        xyz_x: vec![0.0f64; SIZE + 2],
        xyz_y: vec![0.0f64; SIZE + 2],
        xyz_z: vec![0.0f64; SIZE + 2],
        errat: vec![0.0f64; SIZE + 8],
    }
}

fn parse_pdb_mmap<W: Write>(path: &PathBuf, logw: &mut W) -> io::Result<AtomData> {
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    parse_pdb_bytes(&mmap, logw)
}

fn parse_pdb_bytes<W: Write>(bytes: &[u8], logw: &mut W) -> io::Result<AtomData> {
    let mut data = empty_atom_data();
    let mut i: usize = 0;
    let mut atmnum: usize = 0;
    let mut kadd: i32 = 0;
    let mut flag = false;
    let mut flag2 = false;

    let mut start = 0usize;
    while start < bytes.len() && !flag2 {
        let mut end = start;
        while end < bytes.len() && bytes[end] != b'\n' {
            end += 1;
        }
        let mut line = &bytes[start..end];
        if line.ends_with(b"\r") {
            line = &line[..line.len().saturating_sub(1)];
        }
        start = end + 1;

        if line.len() < 6 || &line[..6] != b"ATOM  " {
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
        if line.len() < 54 {
            i -= 1;
            continue;
        }

        let name_temp = line[13];
        data.name[i] = match name_temp {
            b'C' => 1,
            b'N' => 2,
            b'O' => 3,
            _ => 0,
        };

        if line.len() < 16 {
            i -= 1;
            continue;
        }
        let name_temp2 = &line[13..16];
        data.bnam[i] = if name_temp2 == b"N  " || name_temp2 == b"C  " {
            1
        } else {
            0
        };

        let alt_loc = line[16] as char;
        let res_name = &line[17..20];
        data.chain_id[i] = line[21];

        let res_seq_temp = std::str::from_utf8(&line[22..26]).unwrap_or("");
        let res_seq_val = res_seq_temp.trim().parse::<f64>().unwrap_or(0.0);
        data.res_seq[i] = res_seq_val as i32;

        let x_temp = std::str::from_utf8(&line[30..38]).unwrap_or("");
        let y_temp = std::str::from_utf8(&line[38..46]).unwrap_or("");
        let z_temp = std::str::from_utf8(&line[46..54]).unwrap_or("");
        data.xyz_x[i] = x_temp.trim().parse::<f64>().unwrap_or(0.0);
        data.xyz_y[i] = y_temp.trim().parse::<f64>().unwrap_or(0.0);
        data.xyz_z[i] = z_temp.trim().parse::<f64>().unwrap_or(0.0);

        if !(alt_loc == ' ' || alt_loc == 'A' || alt_loc == 'a' || alt_loc == 'P') {
            writeln!(
                logw,
                "Reject 2' Conformation atom#\t{}\tchain\t{}",
                i, data.chain_id[i] as char
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

        if i >= 2 && !flag && data.chain_id[i] != data.chain_id[i - 1] {
            kadd += 1;
            writeln!(logw, "INCREMENTING CHAIN (kadd) {}", kadd)?;
        }

        if !flag {
            data.resnum[i] = data.res_seq[i] + (kadd * CHAINDIF);
            atmnum = i;
        }

        if i >= 2
            && !flag
            && data.chain_id[i] == data.chain_id[i - 1]
            && data.resnum[i] < data.resnum[i - 1]
        {
            writeln!(
                logw,
                "ERROR: RESNUM DECREASE. TERMINATE ANALYSIS{}\t{}",
                data.resnum[i],
                data.resnum[i - 1]
            )?;
            flag2 = true;
        }

        if i > 2
            && !flag
            && data.chain_id[i] == data.chain_id[i - 1]
            && data.resnum[i] != data.resnum[i - 1]
            && (data.resnum[i] - data.resnum[i - 1]) > 1
        {
            writeln!(
                logw,
                "WARNING: Missing Residues{}>>>{}",
                data.resnum[i - 1],
                data.resnum[i]
            )?;
        }

        if !flag {
            let idx = (data.resnum[i] + 4) as usize;
            if idx >= data.errat.len() {
                data.errat.resize(idx + 1, 0.0);
            }
            data.errat[idx] = 0.0;
        }

        flag = false;
    }

    data.atmnum = atmnum;
    Ok(data)
}

fn parse_pdb<R: BufRead, W: Write>(reader: &mut R, logw: &mut W) -> io::Result<AtomData> {
    let mut data = empty_atom_data();
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
        if bytes.len() < 6 || &bytes[..6] != b"ATOM  " {
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
        data.name[i] = match name_temp {
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
        data.bnam[i] = if name_temp2 == b"N  " || name_temp2 == b"C  " {
            1
        } else {
            0
        };

        let alt_loc = bytes[16] as char;
        let res_name = &bytes[17..20];
        data.chain_id[i] = bytes[21];

        let res_seq_temp = std::str::from_utf8(&bytes[22..26]).unwrap_or("");
        let res_seq_val = res_seq_temp.trim().parse::<f64>().unwrap_or(0.0);
        data.res_seq[i] = res_seq_val as i32;

        let x_temp = std::str::from_utf8(&bytes[30..38]).unwrap_or("");
        let y_temp = std::str::from_utf8(&bytes[38..46]).unwrap_or("");
        let z_temp = std::str::from_utf8(&bytes[46..54]).unwrap_or("");
        data.xyz_x[i] = x_temp.trim().parse::<f64>().unwrap_or(0.0);
        data.xyz_y[i] = y_temp.trim().parse::<f64>().unwrap_or(0.0);
        data.xyz_z[i] = z_temp.trim().parse::<f64>().unwrap_or(0.0);

        if !(alt_loc == ' ' || alt_loc == 'A' || alt_loc == 'a' || alt_loc == 'P') {
            writeln!(
                logw,
                "Reject 2' Conformation atom#\t{}\tchain\t{}",
                i, data.chain_id[i] as char
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

        if i >= 2 && !flag && data.chain_id[i] != data.chain_id[i - 1] {
            kadd += 1;
            writeln!(logw, "INCREMENTING CHAIN (kadd) {}", kadd)?;
        }

        if !flag {
            data.resnum[i] = data.res_seq[i] + (kadd * CHAINDIF);
            atmnum = i;
        }

        if i >= 2
            && !flag
            && data.chain_id[i] == data.chain_id[i - 1]
            && data.resnum[i] < data.resnum[i - 1]
        {
            writeln!(
                logw,
                "ERROR: RESNUM DECREASE. TERMINATE ANALYSIS{}\t{}",
                data.resnum[i],
                data.resnum[i - 1]
            )?;
            flag2 = true;
        }

        if i > 2
            && !flag
            && data.chain_id[i] == data.chain_id[i - 1]
            && data.resnum[i] != data.resnum[i - 1]
            && (data.resnum[i] - data.resnum[i - 1]) > 1
        {
            writeln!(
                logw,
                "WARNING: Missing Residues{}>>>{}",
                data.resnum[i - 1],
                data.resnum[i]
            )?;
        }

        if !flag {
            let idx = (data.resnum[i] + 4) as usize;
            if idx >= data.errat.len() {
                data.errat.resize(idx + 1, 0.0);
            }
            data.errat[idx] = 0.0;
        }

        flag = false;
    }

    data.atmnum = atmnum;
    Ok(data)
}

fn parse_mmcif<R: Read, W: Write>(reader: &mut R, logw: &mut W) -> io::Result<AtomData> {
    let mut data = empty_atom_data();
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
                if row[g].as_str() != "ATOM" {
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

            let atom_name = row[idx_atom.expect("idx_atom checked above")].as_str();
            let element = idx_type
                .and_then(|k| row.get(k))
                .map(|s| s.as_str())
                .unwrap_or(atom_name);
            let element_char = element.chars().next().unwrap_or(' ');
            data.name[i] = match element_char {
                'C' | 'c' => 1,
                'N' | 'n' => 2,
                'O' | 'o' => 3,
                _ => 0,
            };
            data.bnam[i] = if atom_name == "N" || atom_name == "C" {
                1
            } else {
                0
            };

            let alt_loc = idx_alt
                .and_then(|k| row.get(k))
                .map(|s| s.as_str())
                .unwrap_or(".");
            let alt_loc_char = alt_loc.chars().next().unwrap_or(' ');
            let alt_loc_char = match alt_loc_char {
                '.' | '?' => ' ',
                c => c,
            };

            let res_name_str = row[idx_res.expect("idx_res checked above")].as_str();
            let res_name_upper = res_name_str.to_ascii_uppercase();
            let res_name = res_name_upper.as_bytes();
            let chain = row[idx_chain.expect("idx_chain checked above")].as_bytes();
            data.chain_id[i] = if chain.is_empty() { b' ' } else { chain[0] };

            let res_seq_val = row[idx_seq.expect("idx_seq checked above")]
                .parse::<f64>()
                .unwrap_or(0.0);
            data.res_seq[i] = res_seq_val as i32;

            data.xyz_x[i] = row[idx_x.expect("idx_x checked above")]
                .parse::<f64>()
                .unwrap_or(0.0);
            data.xyz_y[i] = row[idx_y.expect("idx_y checked above")]
                .parse::<f64>()
                .unwrap_or(0.0);
            data.xyz_z[i] = row[idx_z.expect("idx_z checked above")]
                .parse::<f64>()
                .unwrap_or(0.0);

            if !(alt_loc_char == ' '
                || alt_loc_char == 'A'
                || alt_loc_char == 'a'
                || alt_loc_char == 'P')
            {
                writeln!(
                    logw,
                    "Reject 2' Conformation atom#\t{}\tchain\t{}",
                    i, data.chain_id[i] as char
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

            if i >= 2 && !flag && data.chain_id[i] != data.chain_id[i - 1] {
                kadd += 1;
                writeln!(logw, "INCREMENTING CHAIN (kadd) {}", kadd)?;
            }

            if !flag {
                data.resnum[i] = data.res_seq[i] + (kadd * CHAINDIF);
                atmnum = i;
            }

            if i >= 2
                && !flag
                && data.chain_id[i] == data.chain_id[i - 1]
                && data.resnum[i] < data.resnum[i - 1]
            {
                writeln!(
                    logw,
                    "ERROR: RESNUM DECREASE. TERMINATE ANALYSIS{}\t{}",
                    data.resnum[i],
                    data.resnum[i - 1]
                )?;
                flag2 = true;
            }

            if i > 2
                && !flag
                && data.chain_id[i] == data.chain_id[i - 1]
                && data.resnum[i] != data.resnum[i - 1]
                && (data.resnum[i] - data.resnum[i - 1]) > 1
            {
                writeln!(
                    logw,
                    "WARNING: Missing Residues{}>>>{}",
                    data.resnum[i - 1],
                    data.resnum[i]
                )?;
            }

            if !flag {
                let idx = (data.resnum[i] + 4) as usize;
                if idx >= data.errat.len() {
                    data.errat.resize(idx + 1, 0.0);
                }
                data.errat[idx] = 0.0;
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

    data.atmnum = atmnum;
    Ok(data)
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

fn is_standard_residue(res_name: &[u8]) -> bool {
    matches!(
        res_name,
        b"GLY"
            | b"ALA"
            | b"VAL"
            | b"LEU"
            | b"ILE"
            | b"TYR"
            | b"CYS"
            | b"MET"
            | b"TRP"
            | b"PHE"
            | b"HIS"
            | b"PRO"
            | b"SER"
            | b"THR"
            | b"LYS"
            | b"ARG"
            | b"GLU"
            | b"ASP"
            | b"GLN"
            | b"ASN"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

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

    #[test]
    fn parse_mmcif_reads_basic_atom_loop() {
        let mmcif = "\
data_demo
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM N N . ALA A 1 11.104 13.207 2.100
ATOM C C . ALA A 1 11.504 13.607 2.500
ATOM O O . ALA A 1 11.904 14.007 2.900
";
        let mut reader = Cursor::new(mmcif.as_bytes());
        let mut log = Vec::new();
        let data = parse_mmcif(&mut reader, &mut log).unwrap();
        assert_eq!(data.atmnum, 3);
        assert_eq!(data.chain_id[1], b'A');
        assert_eq!(data.res_seq[1], 1);
    }

    #[test]
    fn tokenize_cif_preserves_semicolon_blocks() {
        let input = "data_demo\n;hello\nworld\n;\n";
        let tokens = tokenize_cif(input);
        assert_eq!(tokens[0], "data_demo");
        assert_eq!(tokens[1], "hello\nworld");
    }
}
