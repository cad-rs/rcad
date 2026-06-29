//! Pipeline stage snapshot — dump intermediate DS + result state for
//! cross-checking against OCCT.
//!
//! Enabled by `RCAD_DUMP_PIPELINE=1`; output goes to `RCAD_DUMP_DIR`.
//! Grid/case labels from `RCAD_DUMP_GRID` / `RCAD_DUMP_CASE`.

use crate::bopds::ds::DS;
use serde_json::json;
use std::path::PathBuf;

pub fn is_enabled() -> bool {
    std::env::var("RCAD_DUMP_PIPELINE").is_ok()
}

pub fn dump_dir() -> PathBuf {
    match std::env::var("RCAD_DUMP_DIR") {
        Ok(d) => PathBuf::from(d),
        Err(_) => PathBuf::from("pipeline_dumps"),
    }
}

pub struct DumpCtx {
    grid: String,
    case: String,
    stage_counter: u32,
    enabled: bool,
}

impl DumpCtx {
    pub fn new(grid: &str, case: &str) -> Self {
        Self {
            grid: grid.to_string(),
            case: case.to_string(),
            stage_counter: 0,
            enabled: is_enabled(),
        }
    }

    pub fn snapshot(&mut self, stage: &str, ds: &DS, brep: Option<&rcad_kernel::topods::BRep>) {
        if !self.enabled { return; }
        self.stage_counter += 1;
        let seq = self.stage_counter;
        let prefix = format!("{}_{}", self.grid, self.case);
        let dir = dump_dir();
        let _ = std::fs::create_dir_all(&dir);
        let ds_v = serialize_ds(ds);
        let p = dir.join(format!("rcad_{}_s{:02}_{}_ds.json", prefix, seq, stage));
        if let Ok(s) = serde_json::to_string_pretty(&ds_v) {
            let _ = std::fs::write(&p, &s);
        }
        if let Some(b) = brep {
            let br_v = serialize_brep(b);
            let p2 = dir.join(format!("rcad_{}_s{:02}_{}_result.json", prefix, seq, stage));
            if let Ok(s) = serde_json::to_string_pretty(&br_v) {
                let _ = std::fs::write(&p2, &s);
            }
        }
        eprintln!("[DUMP] {} stage {}", prefix, stage);
    }
}

fn serialize_ds(ds: &DS) -> serde_json::Value {
    let nv = ds.vertices.len();
    let av = ds.a_vertex_count;
    let ee = cnt(&ds.interferences, |x| matches!(x, crate::bopds::ds::Interference::EdgeEdge{..}));
    let ef = cnt(&ds.interferences, |x| matches!(x, crate::bopds::ds::Interference::EdgeFace{..}));
    json!({"ds": {
        "nV": nv, "nA": av.min(nv), "nB": nv.saturating_sub(av),
        "nE": ds.edges.len(), "nF": ds.faces.len(),
        "nIC": ds.intersection_curves.len(),
        "nPB": ds.pave_blocks.len(), "nCB": ds.common_blocks.len(),
        "interf": { "EE": ee, "EF": ef, "total": ds.interferences.len() },
        "faces": ds.faces.iter().enumerate().map(|(fi, f)| {
            let st = format!("{:?}", f.surface);
            json!({"fi": fi, "surf": st, "nBE": f.boundary_edges.len(),
                "nIW": f.inner_boundary_edges.len(),
                "boundary_edges": f.boundary_edges,
                "nPBsIn": f.face_info.pave_blocks_in.len(),
                "nPBsSc": f.face_info.pave_blocks_sc.len(),
                "nCurvesSc": f.face_info.curves_sc.len(),
                "nVIn": f.face_info.vertices_in.len(),
                "curves_sc": f.face_info.curves_sc.iter().copied().collect::<Vec<_>>(),
                "vertices_in": f.face_info.vertices_in.iter().copied().collect::<Vec<_>>(),
            })
        }).collect::<Vec<_>>(),
        "edges": ds.edges.iter().enumerate().map(|(ei, e)| {
            let ct = format!("{:?}", e.curve);
            json!({"ei": ei, "sv": e.start_vertex, "ev": e.end_vertex,
                "curve": ct,
                "my_images": ds.my_images.get(ei).map(|v| v.clone()).unwrap_or_default(),
                "is_internal": e.is_internal,
                "nPBs": e.pave_blocks.len(),
            })
        }).collect::<Vec<_>>(),
        "intersection_curves": ds.intersection_curves.iter().enumerate().map(|(ci, ic)| {
            let ct = format!("{:?}", ic.curve);
            json!({"ci": ci, "curve": ct, "sv": ic.start_vertex, "ev": ic.end_vertex,
                "t_range": ic.t_range, "has_pca": ic.pcurve_on_a.is_some(),
                "has_pcb": ic.pcurve_on_b.is_some(),
                "n_pave_blocks": ic.pave_blocks.len(),
            })
        }).collect::<Vec<_>>(),
    }})
}

fn serialize_brep(b: &rcad_kernel::topods::BRep) -> serde_json::Value {
    let mut v=0u32;let mut e=0u32;let mut f=0u32;let mut sh=0u32;let mut so=0u32;let mut cp=0u32;let mut cs=0u32;
    for ts in &b.tshapes {
        match &**ts {
            rcad_kernel::topods::TShape::Vertex(_) => v+=1,
            rcad_kernel::topods::TShape::Edge(_) => e+=1,
            rcad_kernel::topods::TShape::Face(_) => f+=1,
            rcad_kernel::topods::TShape::Shell(_) => sh+=1,
            rcad_kernel::topods::TShape::Solid(_) => so+=1,
            rcad_kernel::topods::TShape::Compound(_) => cp+=1,
            rcad_kernel::topods::TShape::CompSolid(_) => cs+=1,
            _ => {}
        }
    }
    json!({"brep": {"V":v,"E":e,"F":f,"Shell":sh,"Solid":so,"Comp":cp,"CompSolid":cs,"total":v+e+f+sh+so+cp+cs}})
}

fn cnt<T>(v: &[T], f: fn(&T) -> bool) -> usize { v.iter().filter(|x| f(x)).count() }
