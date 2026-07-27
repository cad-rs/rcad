// Debug utilities for boolean pipeline tracing.
//
// Usage:
//   dbg_seg!("outer wire", face_idx, seg);
//   dbg_verts!("result edges", result, edges);
//   dbg_block!("blocks", face_idx, &blocks);
//
// Enable at runtime: set DEBUG_BOOL=1 or DEBUG_BOOL=fi,fi,...

/// Print a WireSegment with its source name.
#[macro_export]
macro_rules! dbg_seg {
    ($label:expr, $face_idx:expr, $seg:expr) => {
        if $crate::bopalgo::builder::debug_utils::dbg_enabled_for($face_idx) {
            let __src = match &$seg.source {
                $crate::bopalgo::builder::WireEdgeSource::DsEdge(ei) => format!("DS{}", ei),
                $crate::bopalgo::builder::WireEdgeSource::IntersectionCurve(ci) => format!("IC{}", ci),
                $crate::bopalgo::builder::WireEdgeSource::SeamEdge => "Seam".into(),
            };
            eprintln!("[DBG] {} fi={} {} v{}->{} fwd={}",
                $label, $face_idx, __src, $seg.start_vertex, $seg.end_vertex, $seg.forward);
        }
    };
    ($label:expr, $face_idx:expr, $seg:expr, $extra:expr) => {
        if $crate::bopalgo::builder::debug_utils::dbg_enabled_for($face_idx) {
            let __src = match &$seg.source {
                $crate::bopalgo::builder::WireEdgeSource::DsEdge(ei) => format!("DS{}", ei),
                $crate::bopalgo::builder::WireEdgeSource::IntersectionCurve(ci) => format!("IC{}", ci),
                $crate::bopalgo::builder::WireEdgeSource::SeamEdge => "Seam".into(),
            };
            eprintln!("[DBG] {} fi={} {} v{}->{} fwd={} {}", 
                $label, $face_idx, __src, $seg.start_vertex, $seg.end_vertex, $seg.forward, $extra);
        }
    };
}

/// Print vertex positions from a ResultBuilder edge list.
#[macro_export]
macro_rules! dbg_result_verts {
    ($label:expr, $result:expr) => {{
        let __uniq: std::collections::BTreeSet<usize> = $result.edges.iter()
            .flat_map(|e| [e.0, e.1]).collect();
        eprintln!("[DBG] {} V={} E={}", $label, __uniq.len(), $result.edges.len());
        for &__vi in &__uniq {
            if __vi < $result.vertices.len() {
                let __p = $result.vertices[__vi];
                eprintln!("[DBG]   RV{}: ({:.6},{:.6},{:.6})", __vi, __p.x, __p.y, __p.z);
            }
        }
    }};
}

/// Print wire segment indices per wire.
#[macro_export]
macro_rules! dbg_wires {
    ($label:expr, $face_idx:expr, $wires:expr, $segments:expr) => {
        if $crate::bopalgo::builder::debug_utils::dbg_enabled_for($face_idx) {
            for (__wi, __w) in $wires.iter().enumerate() {
                let __srcs: Vec<String> = __w.iter().map(|&__si| {
                    match &$segments[__si].source {
                        $crate::bopalgo::builder::WireEdgeSource::DsEdge(ei) => format!("DS{}", ei),
                        $crate::bopalgo::builder::WireEdgeSource::IntersectionCurve(ci) => format!("IC{}", ci),
                        $crate::bopalgo::builder::WireEdgeSource::SeamEdge => "Seam".into(),
                    }
                }).collect();
                eprintln!("[DBG] {} fi={} wire[{}]: {} segs [{}]",
                    $label, $face_idx, __wi, __w.len(), __srcs.join(","));
            }
        }
    };
}

/// Print connexity block sizes.
#[macro_export]
macro_rules! dbg_blocks {
    ($label:expr, $face_idx:expr, $blocks:expr) => {
        if $crate::bopalgo::builder::debug_utils::dbg_enabled_for($face_idx) {
            for (__bi, __b) in $blocks.iter().enumerate() {
                eprintln!("[DBG] {} fi={} block[{}]: {} segs", $label, $face_idx, __bi, __b.len());
            }
        }
    };
}

/// Print SmartMap vertex degrees.
#[macro_export]
macro_rules! dbg_smartmap {
    ($label:expr, $face_idx:expr, $smart_map:expr) => {
        if $crate::bopalgo::builder::debug_utils::dbg_enabled_for($face_idx) {
            for (__v, __infos) in $smart_map.iter() {
                eprintln!("[DBG] {} fi={} v{}: in={} out={}",
                    $label, $face_idx, __v,
                    __infos.iter().filter(|ei| ei.in_flag).count(),
                    __infos.iter().filter(|ei| !ei.in_flag).count());
            }
        }
    };
}

/// Runtime filter: DEBUG_BOOL env var.
///   unset → disabled
///   "1"   → all faces
///   "0,1,3" → only faces 0, 1, 3
pub fn dbg_enabled_for(face_idx: usize) -> bool {
    use std::sync::OnceLock;
    static PARSED: OnceLock<DbgConfig> = OnceLock::new();
    let cfg = PARSED.get_or_init(|| {
        let raw = std::env::var("DEBUG_BOOL").unwrap_or_default();
        if raw == "1" { return DbgConfig::All; }
        let ids: Vec<usize> = raw.split(',').filter_map(|s| s.trim().parse().ok()).collect();
        if ids.is_empty() { DbgConfig::None } else { DbgConfig::Filter(ids) }
    });
    match cfg {
        DbgConfig::All => true,
        DbgConfig::Filter(ids) => ids.contains(&face_idx),
        DbgConfig::None => false,
    }
}

enum DbgConfig { All, Filter(Vec<usize>), None }
