// Add this as a debug function in builder.rs to trace the trim polyline

fn debug_print_trim_info(trim_polylines: &[Vec<DVec2>], uv_boundary: &[DVec2]) {
    eprintln!("=== DEBUG: TRIM POLYLINE INFO ===");
    eprintln!("UV Boundary rect: {:?}", uv_boundary);
    for (i, trim) in trim_polylines.iter().enumerate() {
        if trim.len() > 0 {
            let start = trim[0];
            let end = trim[trim.len() - 1];
            let is_closed = (start - end).length_squared() < 1e-6;
            let min_u = trim.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
            let max_u = trim.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
            let min_v = trim.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
            let max_v = trim.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
            
            eprintln!("Trim {}: len={} closed={} start={:.3?} end={:.3?}", 
                i, trim.len(), is_closed, start, end);
            eprintln!("  UV ranges: u=[{:.3}, {:.3}] v=[{:.3}, {:.3}]", 
                min_u, max_u, min_v, max_v);
            
            // Print first 5 and last 5 points
            eprintln!("  First 5 points:");
            for (j, p) in trim.iter().take(5).enumerate() {
                eprintln!("    [{}] u={:7.4} v={:7.4}", j, p.x, p.y);
            }
            if trim.len() > 10 {
                eprintln!("    ... {} more ...", trim.len() - 10);
            }
            eprintln!("  Last 5 points:");
            for (j, p) in trim.iter().skip(trim.len().saturating_sub(5)).enumerate() {
                let actual_idx = trim.len() - (5 - j);
                eprintln!("    [{}] u={:7.4} v={:7.4}", actual_idx, p.x, p.y);
            }
        }
    }
}
