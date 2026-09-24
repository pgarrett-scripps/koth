//! Print koth's scan -> 1/K0 conversion for a Bruker `.d` run, for comparison
//! with Bruker's timsdata SDK (`tims_scannum_to_oneoverk0`).
//!
//! ```text
//! cargo run --release -p koth-ms --example tims_mobility_dump -- RUN.d [FRAME_ID ...]
//! ```
//!
//! Writes TSV `frame_id  scan  calibrated  linear` for every integer scan
//! `0..=max(NumScans)` of each requested frame (default: the first frame).

use std::path::PathBuf;

use koth_ms::io::tims_calibration::{load, MobilityScale};

fn main() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(
        args.next()
            .ok_or("usage: tims_mobility_dump RUN.d [FRAME_ID ...]")?,
    );
    let mut frames: Vec<usize> = args
        .map(|a| a.parse().map_err(|e| format!("frame id {a:?}: {e}")))
        .collect::<Result<_, _>>()?;
    if frames.is_empty() {
        frames.push(1);
    }
    let cal = load(&path, MobilityScale::Calibrated)?;
    let lin = load(&path, MobilityScale::Linear)?;
    let max_scans = max_num_scans(&path)?;
    println!("frame_id\tscan\tcalibrated\tlinear");
    for f in frames {
        for s in 0..=max_scans {
            println!(
                "{f}\t{s}\t{:.17e}\t{:.17e}",
                cal.convert(f, s),
                lin.convert(f, s)
            );
        }
    }
    Ok(())
}

fn max_num_scans(path: &std::path::Path) -> Result<u32, String> {
    let con = rusqlite::Connection::open_with_flags(
        path.join("analysis.tdf"),
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY,
    )
    .map_err(|e| e.to_string())?;
    con.query_row("SELECT MAX(NumScans) FROM Frames", [], |r| r.get(0))
        .map_err(|e| e.to_string())
}
