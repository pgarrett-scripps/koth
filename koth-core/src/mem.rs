/// Read current process RSS (resident set size) from /proc/self/status.
/// Returns MB on Linux, 0.0 on other platforms.
pub fn rss_mb() -> f64 {
    #[cfg(target_os = "linux")]
    {
        if let Ok(s) = std::fs::read_to_string("/proc/self/status") {
            for line in s.lines() {
                if let Some(rest) = line.strip_prefix("VmRSS:") {
                    if let Some(kb_str) = rest.split_whitespace().next() {
                        if let Ok(kb) = kb_str.parse::<u64>() {
                            return kb as f64 / 1024.0;
                        }
                    }
                }
            }
        }
        0.0
    }
    #[cfg(not(target_os = "linux"))]
    {
        0.0
    }
}

/// Log RSS at a named checkpoint.
pub fn log_mem(label: &str) {
    let mb = rss_mb();
    if mb > 0.0 {
        log::info!("[mem] {}: {:.0} MB RSS", label, mb);
    }
}
