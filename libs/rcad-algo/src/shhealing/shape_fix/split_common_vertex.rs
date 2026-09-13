//! 1:1 translation of OCCT `ShapeFix_SplitCommonVertex`
//! (`TKShHealing/ShapeFix/ShapeFix_SplitCommonVertex.hxx` L16-49 +
//! `ShapeFix_SplitCommonVertex.cxx` L1-167, docket row
//! `ShapeFix_SplitCommonVertex`).
//!
//! Function-count equation (OCCT ShapeFix_SplitCommonVertex.cxx = rcad
//! split_common_vertex.rs): `ShapeFix_SplitCommonVertex()` / `Init` /
//! `Perform` / `Shape` — 4 OCCT functions = 4 rcad functions.
//!
//! Architecture bridges:
//! 1. `ShapeFix_Root` inheritance — the rcad composition field `base`
//!    (Rust has no inheritance); the Context()/SetContext()/SendWarning
//!    root members are reached through it.
//! 2. `NCollection_Sequence<TopoDS_Shape>` — `Vec<Shape>`;
//!    `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape,
//!    TopTools_ShapeMapHasher>` — `IndexMap<(u64, u32), Shape>` keyed by
//!    (TShape ptr, location index) (the TopTools_ShapeMapHasher IsSame
//!    identity; the package-statics key convention).
//! 3. `TopoDS_Iterator(F, false)` — the no-cumulative-orientation child
//!    walk (`shape_build::brep_tool::iter_subshapes`).
//! 4. `Message_Msg("Fix.SplitCommonVertex.MSG0")` — the `MessageMsg`
//!    carrier (`shape_extend::msg`).

use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, ShapeType, TShape};

use crate::brep_algo::tool::brep_tool_tolerance;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_build::brep_tool::{iter_subshapes, topexp_explorer};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::root::ShapeFixRoot;

/// OCCT BRep_Tool::Pnt(V).
fn brep_tool_pnt(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT ShapeFix_SplitCommonVertex (hxx L29-47): two wires have a common
/// vertex — valid in the BRep model, invalid in STEP; before writing to
/// STEP the vertex must be split (each wire must have its own vertex).
pub struct ShapeFixSplitCommonVertex {
    /// OCCT the ShapeFix_Root base subobject (bridge #1).
    pub base: ShapeFixRoot,
    /// OCCT myShape (hxx L44).
    my_shape: Shape,
    /// OCCT myResult (hxx L45).
    my_result: Shape,
    /// OCCT myStatus (hxx L46) — set at construction; OCCT never reads it
    /// back in this class.
    #[allow(dead_code)]
    my_status: i32,
}

impl Default for ShapeFixSplitCommonVertex {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixSplitCommonVertex {
    /// OCCT ShapeFix_SplitCommonVertex::ShapeFix_SplitCommonVertex() (cxx
    /// L41-45).
    pub fn new() -> Self {
        let mut this = ShapeFixSplitCommonVertex {
            base: ShapeFixRoot::new(),
            my_shape: Shape::null(),
            my_result: Shape::null(),
            my_status: encode_status(ShapeExtendStatus::Ok),
        };
        // L44: SetPrecision(Precision::Confusion()).
        this.base.set_precision(CONFUSION);
        this
    }

    /// OCCT ShapeFix_SplitCommonVertex::Init (cxx L49-58).
    pub fn init(&mut self, brep: &mut BRep, s: &Shape) {
        // L51.
        self.my_shape = s.clone();
        // L52-55.
        if self.base.context().is_none() {
            self.base.set_context(ShapeBuildReShape::new());
        }
        // L56.
        self.my_result = self.my_shape.clone();
        // L57: Context()->Apply(myShape) — the result is discarded by OCCT.
        if let Some(ctx) = self.base.context_mut() {
            let _ = ctx.apply(brep, &self.my_shape, ShapeType::Shape);
        }
    }

    /// OCCT ShapeFix_SplitCommonVertex::Perform (cxx L62-160).
    pub fn perform(&mut self, brep: &mut BRep) {
        // L64-67.
        let st = self.my_shape.shape_type();
        if st as u32 > ShapeType::Face as u32 {
            return;
        }
        // L69-156: the per-face walk.
        for itf in topexp_explorer(brep, &self.my_shape, ShapeType::Face) {
            // L71-72: tmpFace = Context()->Apply(itf.Current()).
            let tmp_face = match self.base.context_mut() {
                Some(ctx) => ctx.apply(brep, &itf, ShapeType::Shape),
                None => itf.clone(),
            };
            let f = tmp_face;
            // L73-76.
            if f.is_null() {
                continue;
            }
            // L77-86: analys face and split if necessary — the child wires
            // (TopoDS_Iterator(F, false)).
            let mut wires: Vec<Shape> = Vec::new();
            for itw in iter_subshapes(brep, &f, false, true) {
                if itw.shape_type() != ShapeType::Wire {
                    continue;
                }
                wires.push(itw);
            }
            // L87-90.
            if wires.len() < 2 {
                continue;
            }
            // L91-92: MapVV — the split-vertex registry.
            let mut map_vv: IndexMap<(u64, u32), Shape> = IndexMap::new();
            // L93-152: the wire-pair walk.
            for nw1 in 1..wires.len() {
                // L95-96: (sewd1 is created by OCCT but never read after).
                let w1 = wires[nw1 - 1].clone();
                let sewd1 = WireData::new_from_wire(brep, &w1, true, true);
                let _ = &sewd1;
                // L97-151.
                for nw2 in (nw1 + 1)..=wires.len() {
                    // L99-100.
                    let w2 = wires[nw2 - 1].clone();
                    let sewd2 = WireData::new_from_wire(brep, &w2, true, true);

                    // L102-150.
                    for expv1 in topexp_explorer(brep, &w1, ShapeType::Vertex) {
                        let v1 = expv1;
                        for expv2 in topexp_explorer(brep, &w2, ShapeType::Vertex) {
                            let v2 = expv2;
                            // L108: common vertex exists (V1 == V2 —
                            // TopoDS operator== / IsEqual).
                            if !(v1.is_same(&v2) && v1.orientation == v2.orientation) {
                                continue;
                            }
                            // L110-123: the fresh split vertex.
                            let vnew =
                                if let Some(found) = map_vv.get(&(v2.ptr_id(), v2.location)) {
                                    found.clone()
                                } else {
                                    // L118-122.
                                    let p = brep_tool_pnt(&v2);
                                    let tol = brep_tool_tolerance(&v2);
                                    let mut b = BRepBuilder::new();
                                    let vnew = b.add_vertex(brep, p, tol);
                                    // L122: MapVV.Bind(V2, Vnew).
                                    map_vv.insert((v2.ptr_id(), v2.location), vnew.clone());
                                    vnew
                                };
                            // L124-147: replace V2 by Vnew on every edge of
                            // w2.
                            let sbe = ShapeBuildEdge;
                            let sae = ShapeAnalysisEdge::new();
                            let nb_edges = sewd2.nb_edges();
                            for ne2 in 1..=nb_edges {
                                // L128.
                                let e = sewd2.edge(ne2);
                                // L129-131.
                                let mut fv = sae.first_vertex(brep, &e);
                                let mut lv = sae.last_vertex(brep, &e);
                                let mut is_coinc = false;
                                // L132-136.
                                if fv.is_same(&v2) && fv.orientation == v2.orientation {
                                    fv = vnew.clone();
                                    is_coinc = true;
                                }
                                // L137-141.
                                if lv.is_same(&v2) && lv.orientation == v2.orientation {
                                    lv = vnew.clone();
                                    is_coinc = true;
                                }
                                // L142-146.
                                if is_coinc {
                                    // L144: NewE =
                                    // sbe.CopyReplaceVertices(E, FV, LV).
                                    let new_e = sbe.copy_replace_vertices(brep, &e, &fv, &lv);
                                    // L145.
                                    if let Some(ctx) = self.base.context_mut() {
                                        ctx.replace(brep, &e, &new_e);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // L153-156.
            if !map_vv.is_empty() {
                // L155.
                self.base
                    .send_warning_own(&MessageMsg::from_key("Fix.SplitCommonVertex.MSG0"));
            }
        }

        // L159.
        if let Some(ctx) = self.base.context_mut() {
            self.my_shape = ctx.apply(brep, &self.my_shape, ShapeType::Shape);
        }
    }

    /// OCCT ShapeFix_SplitCommonVertex::Shape (cxx L164-167).
    pub fn shape(&self) -> Shape {
        // L166.
        self.my_shape.clone()
    }
}
