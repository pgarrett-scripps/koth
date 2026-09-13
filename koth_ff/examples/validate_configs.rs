//! Parse tuned feature-finder and alignment configs before a release.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    for path in std::env::args().skip(1) {
        let p = std::path::Path::new(&path);
        if p.file_name()
            .unwrap()
            .to_string_lossy()
            .starts_with("koth_align")
        {
            koth_ff::config::AlignConfig::from_toml(p)?;
        } else {
            koth_ff::config::KothConfig::from_toml(p)?;
        }
        println!("parsed {}", path);
    }
    Ok(())
}
