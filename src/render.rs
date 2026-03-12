use std::fmt::Write as FmtWrite;
use std::io::{self, Write};

use crate::api::compute_chain_ranges;
use crate::model::{CHAINDIF, ErratStats, LMT_95, LMT_99};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct PageSlice {
    chain_id: u8,
    start_residue: i32,
    end_residue: i32,
}

#[derive(Clone, Debug, PartialEq)]
struct PlotLayout {
    scale: f64,
    pages: Vec<PageSlice>,
}

fn build_plot_layout(stats: &ErratStats) -> PlotLayout {
    let chain_ranges = compute_chain_ranges(stats);
    if chain_ranges.is_empty() {
        return PlotLayout {
            scale: 1.0,
            pages: Vec::new(),
        };
    }

    let mut mst = 0.0f64;
    for chain in &chain_ranges {
        let span = (chain.end_residue - chain.start_residue + 1) as f64;
        let mut ms = span / 301.0;
        ms = span / ms;
        if ms > mst {
            mst = ms;
        }
        if mst < 200.0 {
            mst = 200.0;
        }
    }

    let mut pages = Vec::new();
    for chain in &chain_ranges {
        let np = 1 + ((chain.end_residue - chain.start_residue + 1) as f64 / mst) as i32;
        for page_index in 1..=np {
            let start_residue = chain.start_residue + (mst as i32) * (page_index - 1);
            let mut end_residue = start_residue + (mst as i32) - 1;
            if end_residue > chain.end_residue {
                end_residue = chain.end_residue;
            }
            pages.push(PageSlice {
                chain_id: chain.chain_id,
                start_residue,
                end_residue,
            });
        }
    }

    PlotLayout {
        scale: 200.0 / mst,
        pages,
    }
}

pub(crate) fn write_ps<P: Write, L: Write>(
    psw: &mut P,
    logw: &mut L,
    file_string: &str,
    stats: &ErratStats,
) -> io::Result<()> {
    let layout = build_plot_layout(stats);
    if layout.pages.is_empty() {
        return Ok(());
    }

    for page in &layout.pages {
        let overall_quality = stats.overall_quality_factor.unwrap_or(0.0);

        writeln!(
            logw,
            "# Chain Label {}:    Residue range {} to {}",
            page.chain_id as char, page.start_residue, page.end_residue
        )?;

        writeln!(psw, "%!PS")?;
        writeln!(psw, "%FIXED")?;
        writeln!(psw, "/sce {{8}} def /scr {{3}} def")?;
        writeln!(
            psw,
            "90 rotate 110 -380 translate /e95 {{11.527}} def /e99 {{17.191}} def"
        )?;
        writeln!(
            psw,
            "/Helvetica findfont 18 scalefont setfont 0.5 setlinewidth"
        )?;
        writeln!(
            psw,
            "/bar1 {{/g {{1 1 1}} def bar}} def /bar2 {{/g {{1 1 0}} def bar}} def"
        )?;
        writeln!(
            psw,
            "/bar3 {{/g {{1 0 0}} def bar}} def /bar {{sce mul /yval exch def"
        )?;
        writeln!(psw, " scr mul /xval exch def")?;
        writeln!(psw, "newpath xval 0 moveto xval yval lineto scr -1 mul 0")?;
        writeln!(
            psw,
            " rlineto 0 yval -1 mul rlineto closepath gsave g setrgbcolor"
        )?;
        writeln!(psw, " fill grestore stroke}} def")?;
        writeln!(psw, "/tick {{newpath 0.5 sub scr mul 0 moveto 0 -3 rlineto")?;
        writeln!(psw, " currentpoint stroke moveto -10 -12 rmoveto}} def")?;

        writeln!(psw, "% VARIABLE")?;
        writeln!(
            psw,
            "{:.3}   {:.3} scale /rlim {{{}}} def",
            layout.scale,
            layout.scale,
            page.end_residue - page.start_residue + 1
        )?;
        writeln!(psw, "gsave 0 30 sce mul 20 add translate ")?;
        writeln!(psw, "0 30 moveto (Chain#:{}) show ", page.chain_id as char)?;
        writeln!(psw, "0 50 moveto (File: {}) show ", file_string)?;
        writeln!(
            psw,
            "0 10 moveto (Overall quality factor**: {:.3})show",
            overall_quality
        )?;
        writeln!(psw, "0 70 moveto (Program: ERRAT2) show")?;
        writeln!(psw, "() show")?;

        writeln!(psw, "% FIXED")?;
        writeln!(
            psw,
            "grestore newpath 0 0 moveto 0 27 sce mul rlineto stroke"
        )?;
        writeln!(
            psw,
            "newpath rlim scr mul 0 moveto 0 27 sce mul rlineto stroke"
        )?;
        writeln!(psw, "newpath 0  0 moveto rlim scr mul 0 rlineto stroke")?;
        writeln!(
            psw,
            "newpath -3 e95 sce mul moveto rlim scr mul 3 add 0 rlineto"
        )?;
        writeln!(
            psw,
            "stroke newpath -3 e99 sce mul moveto rlim scr mul 3 add 0"
        )?;
        writeln!(psw, " rlineto stroke")?;
        writeln!(psw, "newpath 0  27  sce mul moveto rlim scr")?;
        writeln!(psw, " mul 0 rlineto stroke")?;
        writeln!(psw, "rlim scr mul 2 div 100 sub -34")?;
        writeln!(psw, " moveto (Residue # (window center)) show")?;
        writeln!(
            psw,
            "/Helvetica findfont 14 scalefont setfont 0.5 setlinewidth"
        )?;
        writeln!(psw, "-34 e95 sce mul 4 sub moveto (95\\%) show")?;
        writeln!(psw, "-34 e99 sce mul 4 sub moveto (99\\%) show")?;
        writeln!(
            psw,
            "/Helvetica findfont 12 scalefont setfont 0.5 setlinewidth"
        )?;
        writeln!(
            psw,
            "0 -70 moveto (*On the error axis, two lines are drawn to indicate the confidence with) show"
        )?;
        writeln!(
            psw,
            "0 -82 moveto (which it is possible to reject regions that exceed that error value.) show"
        )?;
        writeln!(
            psw,
            "0 -100 moveto (**Expressed as the percentage of the protein for which the calculated) show"
        )?;
        writeln!(
            psw,
            "0 -112 moveto (error value falls below the 95\\% rejection limit.  Good high resolution) show"
        )?;
        writeln!(
            psw,
            "0 -124 moveto (structures generally produce values around 95\\% or higher.  For lower) show"
        )?;
        writeln!(
            psw,
            "0 -136 moveto (resolutions (2.5 to 3A) the average overall quality factor is around 91\\%. ) show"
        )?;
        writeln!(
            psw,
            "/Helvetica findfont 18 scalefont setfont 0.5 setlinewidth"
        )?;
        writeln!(
            psw,
            "gsave -40 -5 translate 90 rotate 80 0 moveto (Error value*)"
        )?;
        writeln!(psw, "show grestore")?;
        writeln!(
            psw,
            "/Helvetica findfont 16 scalefont setfont 0.5 setlinewidth"
        )?;

        for residue in page.start_residue..=page.end_residue {
            if residue % 20 == 0 {
                writeln!(psw, "{} tick        ", residue - page.start_residue + 1)?;
                writeln!(
                    psw,
                    "({}) show\t",
                    residue - (CHAINDIF * (residue / CHAINDIF))
                )?;
            } else if residue % 10 == 0 {
                writeln!(psw, "{} tick\t", residue - page.start_residue + 1)?;
            }
        }

        for residue in page.start_residue..=page.end_residue {
            let mut bar = "bar1";
            if stats.errat[residue as usize] > LMT_95 {
                bar = "bar2";
            }
            if stats.errat[residue as usize] > LMT_99 {
                bar = "bar3";
            }
            let mut val = stats.errat[residue as usize];
            if val > 27.0 {
                val = 27.0;
            }
            writeln!(
                psw,
                "{}\t{:.3} {}",
                residue - page.start_residue + 1,
                val,
                bar
            )?;
        }
        writeln!(psw, "showpage")?;
    }

    Ok(())
}

pub(crate) fn write_pdf<P: Write, L: Write>(
    pdfw: &mut P,
    logw: &mut L,
    file_string: &str,
    stats: &ErratStats,
) -> io::Result<()> {
    let pages = build_pdf_pages(logw, file_string, stats)?;
    let pdf = build_pdf_document(&pages);
    pdfw.write_all(&pdf)?;
    Ok(())
}

fn build_pdf_pages<L: Write>(
    logw: &mut L,
    file_string: &str,
    stats: &ErratStats,
) -> io::Result<Vec<Vec<u8>>> {
    let layout = build_plot_layout(stats);
    if layout.pages.is_empty() {
        return Ok(Vec::new());
    }

    let mut pages = Vec::new();
    for page in &layout.pages {
        let overall_quality = stats.overall_quality_factor.unwrap_or(0.0);

        writeln!(
            logw,
            "# Chain Label {}:    Residue range {} to {}",
            page.chain_id as char, page.start_residue, page.end_residue
        )?;

        let mut page_buf = Vec::new();
        write_pdf_page(
            &mut page_buf,
            file_string,
            stats,
            page.start_residue,
            page.end_residue,
            page.chain_id,
            overall_quality,
            layout.scale,
        );
        pages.push(page_buf);
    }

    Ok(pages)
}

fn write_pdf_page(
    buf: &mut Vec<u8>,
    file_string: &str,
    stats: &ErratStats,
    ir0: i32,
    ir: i32,
    chain_id: u8,
    overall_quality: f64,
    sz: f64,
) {
    let scr = 3.0;
    let sce = 8.0;
    let e95 = 11.527;
    let e99 = 17.191;
    let rlim = (ir - ir0 + 1) as f64;

    let _ = write!(
        buf,
        "q\n0 1 -1 0 0 0 cm\n1 0 0 1 110 -380 cm\n{:.3} 0 0 {:.3} 0 0 cm\n0.5 w\n0 0 0 RG\n0 0 0 rg\n",
        sz, sz
    );

    let header_y = 30.0 * sce + 20.0;
    pdf_text(
        buf,
        0.0,
        header_y + 30.0,
        18.0,
        &format!("Chain#:{}", chain_id as char),
    );
    pdf_text(
        buf,
        0.0,
        header_y + 50.0,
        18.0,
        &format!("File: {}", file_string),
    );
    pdf_text(
        buf,
        0.0,
        header_y + 10.0,
        18.0,
        &format!("Overall quality factor**: {:.3}", overall_quality),
    );
    pdf_text(buf, 0.0, header_y + 70.0, 18.0, "Program: ERRAT2");

    pdf_line(buf, 0.0, 0.0, 0.0, 27.0 * sce);
    pdf_line(buf, rlim * scr, 0.0, rlim * scr, 27.0 * sce);
    pdf_line(buf, 0.0, 0.0, rlim * scr, 0.0);
    pdf_line(buf, -3.0, e95 * sce, rlim * scr + 3.0, e95 * sce);
    pdf_line(buf, -3.0, e99 * sce, rlim * scr + 3.0, e99 * sce);
    pdf_line(buf, 0.0, 27.0 * sce, rlim * scr, 27.0 * sce);

    pdf_text(
        buf,
        rlim * scr / 2.0 - 100.0,
        -34.0,
        18.0,
        "Residue # (window center)",
    );
    pdf_text(buf, -34.0, e95 * sce - 4.0, 14.0, "95%");
    pdf_text(buf, -34.0, e99 * sce - 4.0, 14.0, "99%");

    pdf_text(
        buf,
        0.0,
        -70.0,
        12.0,
        "*On the error axis, two lines are drawn to indicate the confidence with",
    );
    pdf_text(
        buf,
        0.0,
        -82.0,
        12.0,
        "which it is possible to reject regions that exceed that error value.",
    );
    pdf_text(
        buf,
        0.0,
        -100.0,
        12.0,
        "**Expressed as the percentage of the protein for which the calculated",
    );
    pdf_text(
        buf,
        0.0,
        -112.0,
        12.0,
        "error value falls below the 95% rejection limit.  Good high resolution",
    );
    pdf_text(
        buf,
        0.0,
        -124.0,
        12.0,
        "structures generally produce values around 95% or higher.  For lower",
    );
    pdf_text(
        buf,
        0.0,
        -136.0,
        12.0,
        "resolutions (2.5 to 3A) the average overall quality factor is around 91%. )",
    );

    let _ = write!(buf, "q 0 1 -1 0 -40 -5 cm\n");
    pdf_text(buf, 80.0, 0.0, 18.0, "Error value*");
    let _ = write!(buf, "Q\n");

    for residue in ir0..=ir {
        let x = (residue - ir0 + 1) as f64;
        if residue % 20 == 0 {
            let tick_x = (x - 0.5) * scr;
            pdf_line(buf, tick_x, 0.0, tick_x, -3.0);
            let label = residue - (CHAINDIF * (residue / CHAINDIF));
            pdf_text(buf, tick_x - 10.0, -15.0, 16.0, &label.to_string());
        } else if residue % 10 == 0 {
            let tick_x = (x - 0.5) * scr;
            pdf_line(buf, tick_x, 0.0, tick_x, -3.0);
        }
    }

    for residue in ir0..=ir {
        let mut bar = 1;
        if stats.errat[residue as usize] > LMT_95 {
            bar = 2;
        }
        if stats.errat[residue as usize] > LMT_99 {
            bar = 3;
        }
        let mut val = stats.errat[residue as usize];
        if val > 27.0 {
            val = 27.0;
        }
        let x = (residue - ir0 + 1) as f64 * scr;
        let y = val * sce;
        match bar {
            1 => pdf_set_fill_rgb(buf, 1.0, 1.0, 1.0),
            2 => pdf_set_fill_rgb(buf, 1.0, 1.0, 0.0),
            _ => pdf_set_fill_rgb(buf, 1.0, 0.0, 0.0),
        }
        pdf_rect_fill_stroke(buf, x - scr, 0.0, scr, y);
    }

    let _ = write!(buf, "Q\n");
}

fn build_pdf_document(pages: &[Vec<u8>]) -> Vec<u8> {
    let page_count = pages.len();
    let total_objects = 3 + page_count * 2;
    let mut buf = Vec::new();
    let mut offsets = Vec::with_capacity(total_objects);

    let _ = write!(buf, "%PDF-1.4\n%????\n");

    let catalog_id = 1;
    let pages_id = 2;
    let font_id = 3;
    let first_page_id = 4;
    let first_content_id = first_page_id + page_count;

    offsets.push(buf.len());
    let _ = write!(
        buf,
        "{} 0 obj\n<< /Type /Catalog /Pages {} 0 R >>\nendobj\n",
        catalog_id, pages_id
    );

    let mut kids = String::new();
    for i in 0..page_count {
        let _ = write!(kids, "{} 0 R ", first_page_id + i);
    }
    offsets.push(buf.len());
    let _ = write!(
        buf,
        "{} 0 obj\n<< /Type /Pages /Kids [{}] /Count {} >>\nendobj\n",
        pages_id, kids, page_count
    );

    offsets.push(buf.len());
    let _ = write!(
        buf,
        "{} 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n",
        font_id
    );

    for i in 0..page_count {
        let page_id = first_page_id + i;
        let content_id = first_content_id + i;
        offsets.push(buf.len());
        let _ = write!(
            buf,
            "{} 0 obj\n<< /Type /Page /Parent {} 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 {} 0 R >> >> /Contents {} 0 R >>\nendobj\n",
            page_id, pages_id, font_id, content_id
        );
    }

    for (i, content) in pages.iter().enumerate() {
        let content_id = first_content_id + i;
        let mut stream = Vec::new();
        let length = content.len() + 1;
        let _ = write!(stream, "<< /Length {} >>\nstream\n", length);
        stream.extend_from_slice(content);
        stream.push(b'\n');
        stream.extend_from_slice(b"endstream");
        offsets.push(buf.len());
        let _ = write!(buf, "{} 0 obj\n", content_id);
        buf.extend_from_slice(&stream);
        buf.extend_from_slice(b"\nendobj\n");
    }

    let xref_start = buf.len();
    let _ = write!(buf, "xref\n0 {}\n", total_objects + 1);
    buf.extend_from_slice(b"0000000000 65535 f \n");
    for offset in offsets.iter().take(total_objects) {
        let _ = write!(buf, "{:010} 00000 n \n", offset);
    }
    let _ = write!(
        buf,
        "trailer\n<< /Size {} /Root {} 0 R >>\nstartxref\n{}\n%%EOF\n",
        total_objects + 1,
        catalog_id,
        xref_start
    );

    buf
}

fn pdf_escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '(' => out.push_str("\\("),
            ')' => out.push_str("\\)"),
            _ => out.push(ch),
        }
    }
    out
}

fn pdf_text(buf: &mut Vec<u8>, x: f64, y: f64, size: f64, text: &str) {
    let escaped = pdf_escape(text);
    let _ = write!(
        buf,
        "BT /F1 {:.2} Tf 1 0 0 1 {:.3} {:.3} Tm ({}) Tj ET\n",
        size, x, y, escaped
    );
}

fn pdf_line(buf: &mut Vec<u8>, x1: f64, y1: f64, x2: f64, y2: f64) {
    let _ = write!(buf, "{:.3} {:.3} m {:.3} {:.3} l S\n", x1, y1, x2, y2);
}

fn pdf_rect_fill_stroke(buf: &mut Vec<u8>, x: f64, y: f64, w: f64, h: f64) {
    let _ = write!(buf, "{:.3} {:.3} {:.3} {:.3} re B\n", x, y, w, h);
}

fn pdf_set_fill_rgb(buf: &mut Vec<u8>, r: f64, g: f64, b: f64) {
    let _ = write!(buf, "{:.3} {:.3} {:.3} rg\n", r, g, b);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_pdf_document_writes_valid_header() {
        let pdf = build_pdf_document(&[b"BT ET\n".to_vec()]);
        assert!(pdf.starts_with(b"%PDF-1.4\n"));
        assert!(pdf.ends_with(b"%%EOF\n"));
    }

    #[test]
    fn build_plot_layout_splits_long_chains() {
        let stats = ErratStats {
            stat: 1.0,
            pstat: 0.0,
            avg_probability: Some(0.0),
            overall_quality_factor: Some(100.0),
            errat: vec![0.0; 600],
            resnum: {
                let mut v = vec![0; 4];
                v.push(1 - CHAINDIF);
                v.resize(10, 0);
                v[1] = -3;
                v[2] = 0;
                v[3] = 295;
                v
            },
            chain_id: vec![b' '; 10],
            atmnum: 3,
            warning_frames: Vec::new(),
            scored_frames: Vec::new(),
        };
        let layout = build_plot_layout(&stats);
        assert!(!layout.pages.is_empty());
        assert!(layout.scale > 0.0);
    }
}
