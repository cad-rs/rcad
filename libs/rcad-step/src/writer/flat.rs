//! Flat-array BRep representation for STEP export.
//!
//! Converts from topods::BRep (OCCT-aligned) to a flat array
//! representation that the writer uses internally. This is a
//! private submodule of the writer — no longer a public bridge.

use std::collections::HashMap;

use glam::DVec3;
use rcad_kernel::topods;

/// Flat vertex: just a 3D point.
#[derive(Debug, Clone)]
pub struct Vertex {
    pub point: DVec3,
}

/// Flat edge: indices into the vertices array.
#[derive(Debug, Clone)]
pub struct Edge {
    pub start: usize,
    pub end: usize,
}

/// An edge reference with explicit traversal direction inside a Wire.
#[derive(Debug, Clone)]
pub struct WireEdge {
    pub idx: usize,
    pub forward: bool,
}

/// A wire is an ordered list of (edge_idx, direction) pairs.
#[derive(Debug, Clone, Default)]
pub struct Wire {
    pub edges: Vec<WireEdge>,
}

impl Wire {
    pub fn new() -> Self { Self { edges: Vec::new() } }
}

/// A face: boundary wires, surface index, triangulation, sample point.
#[derive(Debug, Clone)]
pub struct Face {
    pub outer_wire: Wire,
    pub inner_wires: Vec<Wire>,
    pub normal: DVec3,
    pub triangles: Vec<[usize; 3]>,
    pub sample_point: Option<DVec3>,
    pub mesh_dirty: bool,
    pub surface_idx: Option<usize>,
}

impl Face {
    /// Shift all edge indices in wires by `ei_offset`.
    pub fn offset_wire_indices(&mut self, ei_offset: usize) {
        for we in &mut self.outer_wire.edges {
            we.idx += ei_offset;
        }
        for w in &mut self.inner_wires {
            for we in &mut w.edges {
                we.idx += ei_offset;
            }
        }
    }
}

/// A shell: ordered list of faces.
#[derive(Debug, Clone)]
pub struct Shell {
    pub faces: Vec<Face>,
}

/// A solid: ordered list of shells.
#[derive(Debug, Clone, Default)]
pub struct Solid {
    pub shells: Vec<Shell>,
    pub tag: Option<String>,
}

/// A compound holds labelled solids (and optionally sub-compounds / compsolids).
#[derive(Debug, Clone)]
pub struct Compound {
    pub solids: Vec<(Option<String>, Solid)>,
    pub comp_solids: Vec<(Option<String>, CompSolid)>,
    pub compounds: Vec<(Option<String>, Compound)>,
}

impl Compound {
    pub fn new() -> Self { Self { solids: Vec::new(), comp_solids: Vec::new(), compounds: Vec::new() } }
    pub fn add_solid(&mut self, tag: Option<String>, solid: Solid) {
        self.solids.push((tag, solid));
    }
    pub fn flatten_solids(&self) -> Vec<Solid> {
        self.solids.iter().map(|(_, s)| s.clone()).collect()
    }
}

/// A compsolid holds ordered solids.
#[derive(Debug, Clone)]
pub struct CompSolid {
    pub solids: Vec<Solid>,
}

/// Flat geometry store.
#[derive(Debug, Clone, Default)]
pub struct GeomStore {
    pub curves: Vec<rcad_kernel::Curve3>,
    pub surfaces: Vec<rcad_kernel::Surface3>,
    pub curve2ds: Vec<rcad_kernel::Curve2d>,
    pub edge_curve: Vec<Option<usize>>,
    pub edge_curve_range: Vec<Option<[f64; 2]>>,
    pub edge_pcurves: Vec<Vec<PCurve>>,
    pub face_surface: Vec<Option<usize>>,
    pub edge_same_parameter: Vec<bool>,
    pub edge_same_range: Vec<bool>,
    pub edge_degenerated: Vec<bool>,
    pub edge_tolerance: Vec<f64>,
    pub vertex_tolerance: Vec<f64>,
    pub face_tolerance: Vec<f64>,
    pub curve2d_range: Vec<Option<[f64; 2]>>,
    pub face_surface_range: Vec<Option<[f64; 4]>>,
    pub face_internal_vertices: Vec<Vec<usize>>,
    pub edge_vertex_params: Vec<Option<[f64; 2]>>,
}

/// A parameter-space curve binding.
#[derive(Debug, Clone)]
pub struct PCurve {
    pub surface_idx: usize,
    pub curve2d_idx: usize,
}

/// Flat BRep representation used internally by the STEP writer.
#[derive(Debug, Clone)]
pub struct FlatBRep {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub solids: Vec<Solid>,
    pub geom: GeomStore,
    pub compound: Option<Compound>,
    pub compsolid: Option<CompSolid>,
}

impl FlatBRep {
    /// Convert from topods::BRep to FlatBRep.
    pub fn from_topods(brep: &topods::BRep) -> Self {
        // Pass 1: collect vertices and edges
        let mut vertices = Vec::new();
        let mut edges = Vec::new();
        let mut v_map = HashMap::new();  // topods tshape index → flat vertex index
        let mut e_map = HashMap::new();  // topods tshape index → flat edge index

        for (ti, ts) in brep.tshapes.iter().enumerate() {
            match &**ts {
                topods::TShape::Vertex(vd) => {
                    v_map.insert(ti, vertices.len());
                    vertices.push(Vertex { point: vd.point });
                }
                topods::TShape::Edge(ed) => {
                    let start = v_map.get(&ed.first.index).copied().unwrap_or(0);
                    let end = v_map.get(&ed.last.index).copied().unwrap_or(0);
                    e_map.insert(ti, edges.len());
                    edges.push(Edge { start, end });
                }
                _ => {}
            }
        }

        // Pass 2: collect geometry stores
        let mut geom = GeomStore::default();
        geom.edge_curve.resize(edges.len(), None);
        geom.edge_curve_range.resize(edges.len(), None);
        geom.edge_pcurves.resize(edges.len(), Vec::new());

        for (ti, ts) in brep.tshapes.iter().enumerate() {
            if let topods::TShape::Edge(ed) = &**ts {
                if let Some(&fi) = e_map.get(&ti) {
                    let ci = ed.curve.as_ref().map(|crv| {
                        let idx = geom.curves.len();
                        geom.curves.push(crv.clone());
                        idx
                    });
                    geom.edge_curve[fi] = ci;
                    geom.edge_curve_range[fi] = Some(ed.range);
                }
            }
        }

        // Pass 2b: collect face surfaces and edge curve2ds
        let mut flat_fi = 0usize;
        for ts in &brep.tshapes {
            if let topods::TShape::Face(fd) = &**ts {
                let surf_idx = fd.surface.as_ref().map(|s| {
                    let idx = geom.surfaces.len();
                    geom.surfaces.push(s.clone());
                    idx
                });
                while geom.face_surface.len() <= flat_fi { geom.face_surface.push(None); }
                geom.face_surface[flat_fi] = surf_idx;
                flat_fi += 1;
            }
        }
        // Collect curve2ds from edge representations
        for ts in &brep.tshapes {
            if let topods::TShape::Edge(ed) = &**ts {
                for rep in &ed.representations {
                    match rep {
                        topods::CurveRepresentation::CurveOnSurface { pcurve, .. } |
                        topods::CurveRepresentation::CurveOnClosedSurface { pcurve1: pcurve, .. } => {
                            geom.curve2ds.push(pcurve.clone());
                        }
                        _ => {}
                    }
                }
            }
        }

        // Pass 3: collect solids/shells/faces/wires
        let mut solids = Vec::new();
        let compound = None;
        let compsolid = None;
        let mut flat_fi = 0usize;

        for ts in &brep.tshapes {
            if let topods::TShape::Solid(sd) = &**ts {
                let mut shells = Vec::new();
                for shell_sr in &sd.shells {
                    if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                        let mut faces = Vec::new();
                        for face_sr in &shd.faces {
                            if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                                let outer_wire = Self::wire_from_topods(brep, &e_map, &fd.outer_wire);
                                let inner_wires: Vec<Wire> = fd.inner_wires.iter()
                                    .map(|sr| Self::wire_from_topods(brep, &e_map, sr))
                                    .collect();
                                let sid = geom.face_surface.get(flat_fi).copied().flatten();
                                faces.push(Face {
                                    outer_wire,
                                    inner_wires,
                                    normal: DVec3::Z,
                                    triangles: Vec::new(),
                                    sample_point: fd.sample_point,
                                    mesh_dirty: true,
                                    surface_idx: sid,
                                });
                                flat_fi += 1;
                            }
                        }
                        shells.push(Shell { faces });
                    }
                }
                solids.push(Solid { shells, tag: None });
            }
        }

        FlatBRep {
            vertices,
            edges,
            solids,
            geom,
            compound,
            compsolid,
        }
    }

    fn wire_from_topods(
        brep: &topods::BRep,
        e_map: &HashMap<usize, usize>,
        wire_sr: &topods::ShapeRef,
    ) -> Wire {
        if let topods::TShape::Wire(wd) = &*brep.tshapes[wire_sr.index] {
            let edges: Vec<WireEdge> = wd.edges.iter().map(|sr| {
                let flat_ei = e_map.get(&sr.index).copied().unwrap_or(0);
                WireEdge {
                    idx: flat_ei,
                    forward: sr.orientation.is_forward(),
                }
            }).collect();
            Wire { edges }
        } else {
            Wire::new()
        }
    }

    pub fn is_compound(&self) -> bool {
        self.compound.is_some()
    }

    pub fn is_compsolid(&self) -> bool {
        self.compsolid.is_some()
    }

    pub fn as_compound(&self) -> Option<&Compound> {
        self.compound.as_ref()
    }

    pub fn as_compsolid(&self) -> Option<&CompSolid> {
        self.compsolid.as_ref()
    }
}

/// Type alias so STEP writer can use `BRep` directly.
pub type BRep = FlatBRep;

/// STEP export uncertainty for FlatBRep (analogous to `rcad_kernel::tolerance::step_export_uncertainty`).
pub fn step_export_uncertainty(brep: &FlatBRep) -> f64 {
    // Use max of vertex/edge tolerances if populated, else default 1e-6.
    let max_vert = brep.vertices.iter().enumerate()
        .filter_map(|(i, _)| brep.geom.vertex_tolerance.get(i))
        .copied()
        .fold(0.0_f64, f64::max);
    let max_edge = brep.edges.iter().enumerate()
        .filter_map(|(i, _)| brep.geom.edge_tolerance.get(i))
        .copied()
        .fold(0.0_f64, f64::max);
    let tol = max_vert.max(max_edge);
    if tol > 0.0 { tol * 10.0 } else { 1e-6 }
}
