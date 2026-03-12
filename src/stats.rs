use std::io::{self, Write};

use rayon::prelude::*;

use crate::model::{
    AtomData, BOXSIZE, BXMX, ErratStats, FrameScoreRaw, LMT_95, LMT_99, MAXWIN, RADIUS, RADMIN,
};

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

                        if data.resnum[rer] == data.resnum[n] {
                            continue;
                        }

                        let dx = data.xyz_x[n] - rer_x;
                        let dy = data.xyz_y[n] - rer_y;
                        let dz = data.xyz_z[n] - rer_z;
                        let dsq = dx * dx + dy * dy + dz * dz;
                        if dsq >= rsq {
                            continue;
                        }

                        if data.bnam[rer] == 1
                            && data.bnam[n] == 1
                            && (((data.resnum[n] == data.resnum[rer] + 1)
                                && (data.name[rer] == 1)
                                && (data.name[n] == 2))
                                || ((data.resnum[rer] == data.resnum[n] + 1)
                                    && (data.name[rer] == 2)
                                    && (data.name[n] == 1)))
                        {
                            continue;
                        }

                        let temp1 = if dsq <= ssq {
                            1.0
                        } else {
                            2.0 * (RADIUS - dsq.sqrt())
                        };

                        if n >= i && n <= v {
                            if data.resnum[rer] > data.resnum[n] {
                                c[data.name[rer] as usize][data.name[n] as usize] += temp1;
                            }
                        } else {
                            c[data.name[rer] as usize][data.name[n] as usize] += temp1;
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

pub(crate) fn compute_errat<W: Write>(data: &AtomData, logw: &mut W) -> io::Result<ErratStats> {
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
            avg_probability: None,
            overall_quality_factor: None,
            errat: data.errat.clone(),
            resnum: data.resnum.clone(),
            chain_id: data.chain_id.clone(),
            atmnum: data.atmnum,
            warning_frames: Vec::new(),
            scored_frames: Vec::new(),
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
    let mut warning_frames = Vec::new();
    let mut scored_frames = Vec::new();

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
                        warning_frames.push(frame);
                        writeln!(
                            logw,
                            "WARNING: Frame\t{}\tBelow Minimum Interaction Limit.",
                            frame
                        )?;
                    }
                    WindowOutcome::Value { idx, mtrx } => {
                        stat += 1.0;
                        mtrxstat += mtrx;
                        scored_frames.push(FrameScoreRaw {
                            center_residue: idx as i32,
                            error_value: mtrx,
                        });

                        if mtrx > LMT_99 || mtrx > LMT_95 {
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

    let avg_probability = (stat > 0.0).then_some(mtrxstat / stat);
    let overall_quality_factor = (stat > 0.0).then_some(100.0 - (100.0 * pstat / stat));

    if stat > 0.0 {
        writeln!(
            logw,
            "Total frames: {}\tP frames {}\tNumber: {}",
            stat as i64,
            pstat as i64,
            fmt_sig6(pstat / stat)
        )?;
        writeln!(logw)?;
        writeln!(
            logw,
            "Avg Probability\t{}",
            fmt_sig6(avg_probability.unwrap_or(0.0))
        )?;
        writeln!(
            logw,
            "# Overall quality factor: {}",
            fmt_sig6(overall_quality_factor.unwrap_or(0.0))
        )?;
    }

    Ok(ErratStats {
        stat,
        pstat,
        avg_probability,
        overall_quality_factor,
        errat,
        resnum: data.resnum.clone(),
        chain_id: data.chain_id.clone(),
        atmnum: data.atmnum,
        warning_frames,
        scored_frames,
    })
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
    fn compute_errat_empty_structure_returns_zero_stats() {
        let data = AtomData {
            atmnum: 0,
            name: vec![0; 4],
            bnam: vec![0; 4],
            chain_id: vec![b' '; 4],
            res_seq: vec![0; 4],
            resnum: vec![0; 4],
            xyz_x: vec![0.0; 4],
            xyz_y: vec![0.0; 4],
            xyz_z: vec![0.0; 4],
            errat: vec![0.0; 8],
        };
        let mut log = Vec::new();
        let stats = compute_errat(&data, &mut log).unwrap();
        assert_eq!(stats.stat, 0.0);
        assert!(stats.avg_probability.is_none());
        assert!(stats.warning_frames.is_empty());
    }
}
