//! Pipeline stage snapshot — dump intermediate DS + result state for
//! cross-checking against OCCT.
//!
//! Enabled by `RCAD_DUMP_PIPELINE=1`; output goes to `RCAD_DUMP_DIR`.
//! Grid/case labels from `RCAD_DUMP_GRID` / `RCAD_DUMP_CASE`.

use crate::bopds::ds::DS;
use crate::bopds::face_info::FaceInfo;
use rcad_kernel::topods::ShapeType;
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
    module: &'static str,
    stage_counter: u32,
    enabled: bool,
}

impl DumpCtx {
    pub fn new(grid: &str, case: &str) -> Self {
        Self::new_with_module(grid, case, "")
    }

    /// Create a DumpCtx with a module suffix for the filename.
    /// E.g. module="pf" → "rcad_{grid}_{case}_pf_s01_..."
    pub fn new_with_module(grid: &str, case: &str, module: &'static str) -> Self {
        Self {
            grid: grid.to_string(),
            case: case.to_string(),
            module,
            stage_counter: 0,
            enabled: is_enabled(),
        }
    }

    pub fn snapshot(&mut self, stage: &str, ds: &DS, brep: Option<&rcad_kernel::topods::BRep>) {
        if !self.enabled { return; }
        self.stage_counter += 1;
        let seq = self.stage_counter;
        let mod_part = if self.module.is_empty() { String::new() } else { format!("{}_", self.module) };
        let prefix = format!("{}_{}", self.grid, self.case);
        let dir = dump_dir();
        let _ = std::fs::create_dir_all(&dir);
        let ds_v = serialize_ds(ds);
        let p = dir.join(format!("rcad_{}_{}s{:02}_{}_ds.json", prefix, mod_part, seq, stage));
        if let Ok(s) = serde_json::to_string_pretty(&ds_v) {
            let _ = std::fs::write(&p, &s);
        }
        if let Some(b) = brep {
            let br_v = serialize_brep(b);
            let p2 = dir.join(format!("rcad_{}_{}s{:02}_{}_result.json", prefix, mod_part, seq, stage));
            if let Ok(s) = serde_json::to_string_pretty(&br_v) {
                let _ = std::fs::write(&p2, &s);
            }
        }
        eprintln!("[DUMP] {} {}stage {}", prefix, mod_part, stage);
    }
}

/// Round a DVec3 to 6 decimal places for compact serialization.
fn round_pt(p: glam::DVec3) -> [f64; 3] {
    [format!("{:.6}", p.x).parse().unwrap_or(p.x),
     format!("{:.6}", p.y).parse().unwrap_or(p.y),
     format!("{:.6}", p.z).parse().unwrap_or(p.z)]
}

fn serialize_ds(ds: &DS) -> serde_json::Value {
    let nv = ds.vertices.len();
    let av = ds.a_vertex_count;
    let n_interf_total = ds.interf_vv.len() + ds.interf_ve.len() + ds.interf_vf.len()
        + ds.interf_ee.len() + ds.interf_ef.len() + ds.interf_ff.len();

    // Per-vertex details: coordinates, origin, is_internal
    let vertices: Vec<serde_json::Value> = ds.vertices.iter().enumerate().map(|(vi, v)| {
        json!({"vi": vi, "pt": round_pt(v.point), "tol": format!("{:.2e}", v.geom_tol),
            "origin": format!("{:?}", v.origin),
            "is_internal": v.is_internal,
        })
    }).collect();

    // Per-edge details: vertices, curve type, PBs, face reps, my_images
    let edges: Vec<serde_json::Value> = ds.edges.iter().enumerate().map(|(ei, e)| {
        let ct = format!("{:?}", e.curve);
        let my_img = ds.my_images.get(ei).cloned().unwrap_or_default();
        json!({"ei": ei, "sv": e.start_vertex, "ev": e.end_vertex,
            "curve": ct, "tol": format!("{:.2e}", e.geom_tol),
            "origin": format!("{:?}", e.origin),
            "is_internal": e.is_internal, "nPBs": e.pave_blocks.len(),
            "n_face_reps": e.face_reps.len(),
            "my_images": my_img,
        })
    }).collect();

    // Per-face details: surface type, edges, PBs, vertices in
    let faces: Vec<serde_json::Value> = ds.faces.iter().enumerate().map(|(fi, f)| {
        let st = format!("{:?}", f.surface);
        let fi_info = if fi < ds.faces.len() { &ds.faces[fi].face_info }
                      else { &FaceInfo::default() };
        json!({"fi": fi, "surf": st, "nBE": f.boundary_edges.len(),
            "nIW": f.inner_boundary_edges.len(),
            "boundary_edges": f.boundary_edges,
            "nPBsIn": fi_info.pave_blocks_in.len(),
            "nPBsSc": fi_info.pave_blocks_sc.len(),
            "nCurvesSc": fi_info.curves_sc.len(),
            "nVIn": fi_info.vertices_in.len(),
            "curves_sc": fi_info.curves_sc.iter().copied().collect::<Vec<_>>(),
            "vertices_in": fi_info.vertices_in.iter().copied().collect::<Vec<_>>(),
        })
    }).collect();

    // Per-intersection-curve details
    let ics: Vec<serde_json::Value> = ds.intersection_curves.iter().enumerate().map(|(ci, ic)| {
        let ct = format!("{:?}", ic.curve);
        json!({"ci": ci, "curve": ct, "sv": ic.start_vertex, "ev": ic.end_vertex,
            "t_range": [format!("{:.6}", ic.t_range[0]), format!("{:.6}", ic.t_range[1])],
            "has_pca": ic.pcurve_on_a.is_some(),
            "has_pcb": ic.pcurve_on_b.is_some(),
            "n_pave_blocks": ic.pave_blocks.len(),
        })
    }).collect();

    // PaveBlock pool detail: each PB's original_edge, indices, range, new_edge
    let pbs: Vec<serde_json::Value> = ds.pave_blocks.iter().enumerate().map(|(pi, pb)| {
        let r = pb.0.read().unwrap();
        let (v1, v2) = r.indices();
        let (t1, t2) = r.range();
        json!({"pbi": pi, "orig_edge": r.original_edge, "v1": v1, "v2": v2,
            "t_range": [format!("{:.6}", t1), format!("{:.6}", t2)],
            "new_edge": r.new_edge,
            "common_block": r.common_block_idx,
            "n_ext": r.ext_paves.len(),
        })
    }).collect();

    // Common block pool summary
    let cbs: Vec<serde_json::Value> = ds.common_blocks.iter().enumerate().map(|(cbi, cb)| {
        let faces: Vec<usize> = cb.faces().iter().copied().collect();
        let pbs: Vec<(usize, usize)> = cb.pave_blocks().iter().copied().collect();
        json!({"cbi": cbi, "faces": faces,
            "n_pbs": pbs.len(),
            "pbs": pbs.iter().map(|(pbi, _fi)| pbi).collect::<Vec<_>>(),
        })
    }).collect();

    // Vertex-Vertex interference details (VV): which vertex pairs are merged
    let interf_vv: Vec<serde_json::Value> = ds.interf_vv.iter().map(|ivv| {
        json!({"v1": ivv.v1, "v2": ivv.v2, "merged": ivv.merged_vertex})
    }).collect();

    // Edge interference details (EE): which edge pairs intersect
    let interf_ee: Vec<serde_json::Value> = ds.interf_ee.iter().map(|iee| {
        json!({"e1": iee.e1, "e2": iee.e2, "new_v": iee.new_vertex})
    }).collect();

    // Edge-face interference details (EF): which edge is on which face
    let interf_ef: Vec<serde_json::Value> = ds.interf_ef.iter().map(|ief| {
        json!({"edge": ief.edge, "face": ief.face, "new_v": ief.new_vertex})
    }).collect();

    // Face-face interference details (FF): which face pairs intersect
    let interf_ff: Vec<serde_json::Value> = ds.interf_ff.iter().map(|iff| {
        json!({"f1": iff.f1, "f2": iff.f2, "n_curves": iff.curves.len(),
            "n_points": iff.points.len(), "tangent": iff.tangent_faces,
            "curves": iff.curves,
        })
    }).collect();

    // my_images map: summarize how many images per original edge
    let my_images_summary: Vec<serde_json::Value> = ds.my_images.iter().enumerate()
        .filter(|(_, imgs)| !imgs.is_empty())
        .map(|(ei, imgs)| {
            json!({"orig_edge": ei, "n_images": imgs.len(), "images": imgs})
        }).collect();

    // Wire images: which wires have split edges
    let wire_images_summary: Vec<serde_json::Value> = ds.wire_images.iter().enumerate()
        .filter(|(_, wi)| wi.is_some())
        .map(|(wi, wire_img)| {
            json!({"wire_idx": wi, "n_edges": wire_img.as_ref().map(|v| v.len()).unwrap_or(0)})
        }).collect();

    json!({"ds": {
        // OCCT-aligned: count ShapeInfo entries by type (matches OCCT's
        // BOPDS_ShapeInfo::ShapeType counting in pipeline_dump.h).
        "nV": ds.shape_info.iter().filter(|si| si.shape_type == ShapeType::Vertex).count(),
        "nE": ds.shape_info.iter().filter(|si| si.shape_type == ShapeType::Edge).count(),
        "nF": ds.shape_info.iter().filter(|si| si.shape_type == ShapeType::Face).count(),
        "nPB": ds.pave_blocks.len(), "nCB": ds.common_blocks.len(),
        // rcad raw counts (total entities, all types)
        "nV_raw": nv, "nE_raw": ds.edges.len(), "nF_raw": ds.faces.len(),
        "nIC": ds.intersection_curves.len(),
        // OCCT-aligned: source & total shape count
        "nSource": ds.nb_source_shapes(),
        "nTotal": ds.shape_info.len(),
        "interf": { "VV": interf_vv.len(), "EE": interf_ee.len(), "EF": interf_ef.len(), "FF": interf_ff.len(), "total": n_interf_total },
        "vertices": vertices,
        "edges": edges,
        "faces": faces,
        "intersection_curves": ics,
        "pave_blocks": pbs,
        "common_blocks": cbs,
        "interf_vv": interf_vv,
        "interf_ee": interf_ee,
        "interf_ef": interf_ef,
        "interf_ff": interf_ff,
        "my_images": my_images_summary,
        "wire_images": wire_images_summary,
        "n_shells": ds.shells.len(),
        "n_solids": ds.solids.len(),
        "n_comp_solids": ds.comp_solids.len(),
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


