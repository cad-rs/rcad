// OCCT BRepLib_MakeWire.hxx L26-105 + BRepLib_MakeWire.cxx L35-489 — 1:1
// translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepLib/BRepLib_MakeWire.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepLib/BRepLib_MakeWire.cxx
//
// OCCT inheritance chain (hxx L36): BRepLib_MakeWire : BRepBuilderAPI_MakeShape
// : BRepBuilderAPI_Command : Standard_Transient. Rust has no inheritance -> a
// plain owned struct; the base-class members actually read by this class
// (myShape from BRepBuilderAPI_MakeShape, myError from BRepLib_MakeShape) are
// carried as fields.
//
// Architecture differences (referenced from the affected functions):
// 1. BRep_Builder mutations -> the crate::brep_algo::tool in-place re-hosts
//    (builder_make_wire / builder_add_wire_edge / builder_add_edge_vertex /
//    builder_set_closed / builder_update_vertex_point_tol), the established
//    rcad vehicle for `Arc::make_mut` edits on the Shape payload
//    (feat loc_ope_spliter.rs arch. diff. #6).
// 2. NCollection_IndexedMap myVertices -> the insertion-ordered
//    ShapeKeyIndexMap below (Add returns the 1-based index of an
//    already-bound key; Contains / Extent / FindKey(i) are carried 1:1).
// 3. TopoDS_Iterator(E) over an edge -> the stored TEdgeData.first/last
//    children with the parent orientation composed (brep_algo::tool
//    sub_shapes).
// 4. TopExp::Vertices(W, VF, VL) / (E, V1, V2) -> the brep_algo::tool
//    re-hosts (top_exp_vertices_wire / top_exp_vertices_raw).
// 5. The OCCT_DEBUG compile-time branches are not translated (not compiled
//    in the reference build).
// 6. BRep_Builder::Transfert(Ein, Eout, Vin, Vout) -> the exact re-host below
//    (BRep_Builder.cxx L1457-1465: UpdateVertex(Vout, Parameter(Vin, Ein),
//    Eout, Tolerance(Vin))).

use crate::brep_algo::tool::{
    brep_tool_parameter, brep_tool_pnt, brep_tool_tolerance, builder_add_edge_vertex,
    builder_add_wire_edge, builder_make_wire, builder_set_closed,
    builder_update_vertex_point_tol, empty_copied, oriented, shape_key, sub_shapes,
    top_exp_vertices_raw, top_exp_vertices_wire,
};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, TShape};
use std::collections::HashMap;

/// OCCT NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> — an
/// insertion-ordered carrier whose Add returns the 1-based index of an
/// already-bound key (architecture difference #2).
#[derive(Default)]
struct ShapeKeyIndexMap {
    keys: Vec<(u64, u32)>,
    items: HashMap<(u64, u32), Shape>,
}

impl ShapeKeyIndexMap {
    fn new() -> Self {
        ShapeKeyIndexMap {
            keys: Vec::new(),
            items: HashMap::new(),
        }
    }

    /// OCCT NCollection_IndexedMap::Add — the 1-based index (existing keys
    /// keep their index).
    fn add(&mut self, the_s: &Shape) -> usize {
        let k = shape_key(the_s);
        if let Some(pos) = self.keys.iter().position(|&x| x == k) {
            return pos + 1;
        }
        self.keys.push(k);
        self.items.insert(k, the_s.clone());
        self.keys.len()
    }

    /// OCCT NCollection_IndexedMap::Contains.
    fn contains(&self, the_s: &Shape) -> bool {
        self.items.contains_key(&shape_key(the_s))
    }

    /// OCCT NCollection_IndexedMap::Extent.
    fn extent(&self) -> usize {
        self.keys.len()
    }

    /// OCCT NCollection_IndexedMap::FindKey(i) — 1-based.
    fn find_key(&self, the_i: usize) -> Option<&Shape> {
        self.items.get(self.keys.get(the_i - 1)?)
    }
}

/// OCCT BRepLib_WireError (BRepLib_MakeShape.hxx).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum WireError {
    EmptyWire,
    DisconnectedWire,
    NonManifoldWire,
    WireDone,
}

/// OCCT BRepLib_MakeWire (BRepLib_MakeWire.hxx L36-105).
pub(crate) struct MakeWire {
    my_shape: Shape,                   // OCCT: myShape (base MakeShape)
    my_edge: Shape,                    // OCCT: myEdge
    my_vertex: Shape,                  // OCCT: myVertex
    vf: Shape,                         // OCCT: VF
    vl: Shape,                         // OCCT: VL
    my_vertices: ShapeKeyIndexMap,     // OCCT: myVertices
    my_error: WireError,               // OCCT: myError
}

impl Default for MakeWire {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeWire {
    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire() (cxx L35-38).
    pub fn new() -> Self {
        MakeWire {
            my_shape: Shape::null(),
            my_edge: Shape::null(),
            my_vertex: Shape::null(),
            vf: Shape::null(),
            vl: Shape::null(),
            my_vertices: ShapeKeyIndexMap::new(),
            my_error: WireError::EmptyWire,
        }
    }

    /// OCCT BRepLib_MakeWire::Add(const TopoDS_Edge& E) (cxx L110-113).
    pub fn add_edge(&mut self, the_e: &Shape) {
        self.add_edge_proximity(the_e, true);
    }

    /// OCCT BRepLib_MakeWire::Add(const TopoDS_Wire& W) (cxx L96-106).
    pub fn add_wire(&mut self, the_w: &Shape) {
        for it in sub_shapes(the_w) {
            self.add_edge(&it);
            if self.my_error != WireError::WireDone {
                break;
            }
        }
    }

    /// OCCT BRepLib_MakeWire::Wire() (cxx L457-460).
    pub fn wire(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT BRepLib_MakeWire::Error() (cxx L485-488).
    pub fn error(&self) -> WireError {
        self.my_error
    }

    /// OCCT BRepLib_MakeWire::Add(E, IsCheckGeometryProximity) (cxx L123-453).
    pub fn add_edge_proximity(&mut self, the_e: &Shape, is_check_geometry_proximity: bool) {
        // OCCT L126-133.
        let mut forward = false;
        let mut reverse = false;
        let mut init = false;

        if self.my_edge.is_null() {
            init = true;
            // first edge, create the wire (OCCT L137-139).
            self.my_shape = builder_make_wire();
            self.my_shape.orientation = Orientation::Forward;

            // set the edge (OCCT L142).
            self.my_edge = the_e.clone();

            // add the vertices (OCCT L145-148).
            for it in sub_shapes(&self.my_edge) {
                self.my_vertices.add(&it);
            }
        } else {
            // OCCT L153: init = myShape.Closed();
            init = self.my_shape_closed();
            // OCCT L154-155: EE = E.Oriented(TopAbs_FORWARD).
            let ee = oriented(the_e, Orientation::Forward);

            // test the vertices of the edge (OCCT L160-161).
            let mut connected = false;
            let mut copyedge = false;

            // OCCT L163-169.
            if self.my_error != WireError::NonManifoldWire {
                if self.vf.is_null() || self.vl.is_null() {
                    self.my_error = WireError::NonManifoldWire;
                }
            }

            // OCCT L171-286.
            for it in sub_shapes(&ee) {
                let ve = it.clone();

                // if the vertex is in the wire, ok for the connection
                // (OCCT L177-222).
                if self.my_vertices.contains(&ve) {
                    connected = true;
                    self.my_vertex = ve.clone();
                    if self.my_error != WireError::NonManifoldWire {
                        if self.vf.is_same(&self.vl) {
                            // Orientation indetermined (in 3d): preserve the
                            // initial.
                            if !self.vf.is_same(&ve) {
                                self.my_error = WireError::NonManifoldWire;
                            }
                        } else if self.vf.is_same(&ve) {
                            if ve.orientation == Orientation::Forward {
                                reverse = true;
                            } else {
                                forward = true;
                            }
                        } else if self.vl.is_same(&ve) {
                            if ve.orientation == Orientation::Reversed {
                                reverse = true;
                            } else {
                                forward = true;
                            }
                        } else {
                            self.my_error = WireError::NonManifoldWire;
                        }
                    }
                } else if is_check_geometry_proximity {
                    // search if there is a similar vertex in the edge
                    // (OCCT L226-285).
                    let pe = brep_tool_pnt(&ve).unwrap_or(glam::DVec3::ZERO);

                    for i in 1..=self.my_vertices.extent() {
                        let vw = self.my_vertices.find_key(i).expect("myVertices").clone();
                        let pw = brep_tool_pnt(&vw).unwrap_or(glam::DVec3::ZERO);
                        let l = (pe - pw).length();

                        if l < brep_tool_tolerance(&ve) || l < brep_tool_tolerance(&vw) {
                            copyedge = true;
                            if self.my_error != WireError::NonManifoldWire {
                                if self.vf.is_same(&self.vl) {
                                    // Orientation indetermined (in 3d):
                                    // preserve the initial.
                                    if !self.vf.is_same(&vw) {
                                        self.my_error = WireError::NonManifoldWire;
                                    }
                                } else if self.vf.is_same(&vw) {
                                    if ve.orientation == Orientation::Forward {
                                        reverse = true;
                                    } else {
                                        forward = true;
                                    }
                                } else if self.vl.is_same(&vw) {
                                    if ve.orientation == Orientation::Reversed {
                                        reverse = true;
                                    } else {
                                        forward = true;
                                    }
                                } else {
                                    self.my_error = WireError::NonManifoldWire;
                                }
                            }
                            break;
                        }
                    }
                    if copyedge {
                        connected = true;
                    }
                }
            }

            // OCCT L288-368.
            if !connected {
                self.my_error = WireError::DisconnectedWire;
                self.not_done();
                return;
            } else if !copyedge {
                // OCCT L298-303.
                self.my_edge = ee.clone();
                for it in sub_shapes(&ee) {
                    self.my_vertices.add(&it);
                }
            } else {
                // copy the edge (OCCT L305-368).
                let mut my_edge = empty_copied(&ee);

                for ve in sub_shapes(&ee) {
                    let pe = brep_tool_pnt(&ve).unwrap_or(glam::DVec3::ZERO);

                    let mut newvertex = false;
                    for i in 1..=self.my_vertices.extent() {
                        let vw = self.my_vertices.find_key(i).expect("myVertices").clone();
                        let pw = brep_tool_pnt(&vw).unwrap_or(glam::DVec3::ZERO);
                        let l = (pe - pw).length();
                        let tol_w = brep_tool_tolerance(&vw);
                        let tol_e = brep_tool_tolerance(&ve);

                        if l < tol_e || l < tol_w {
                            // OCCT L328-345.
                            let mut maxtol = 0.5 * (tol_w + tol_e + l);
                            let c_w;
                            let c_e;
                            if maxtol > tol_w && maxtol > tol_e {
                                c_w = (maxtol - tol_e) / l;
                                c_e = 1.0 - c_w;
                            } else if maxtol > tol_w {
                                maxtol = tol_e;
                                c_w = 0.0;
                                c_e = 1.0;
                            } else {
                                maxtol = tol_w;
                                c_w = 1.0;
                                c_e = 0.0;
                            }

                            let pc = glam::DVec3::new(
                                c_w * pw.x + c_e * pe.x,
                                c_w * pw.y + c_e * pe.y,
                                c_w * pw.z + c_e * pe.z,
                            );

                            // OCCT L351-357.
                            let mut vw_moved = vw.clone();
                            builder_update_vertex_point_tol(&mut vw_moved, pc, maxtol);

                            newvertex = true;
                            let mut my_vertex = oriented(&vw_moved, ve.orientation);
                            // OCCT L356-357: B.Add(myEdge, myVertex);
                            // B.Transfert(EE, myEdge, VE, myVertex). The
                            // Transfert fork of the vertex Arc is repaired by
                            // re-binding the mutated vertex as the edge child
                            // (the rcad vehicle for OCCT's shared-TShape
                            // update; feat loc_ope_spliter.rs arch. diff. #6).
                            builder_add_edge_vertex(&mut my_edge, &my_vertex);
                            transfert(&ee, &mut my_edge, &ve, &mut my_vertex);
                            builder_add_edge_vertex(&mut my_edge, &my_vertex);
                            self.my_vertex = my_vertex;
                            break;
                        }
                    }
                    if !newvertex {
                        // OCCT L361-366.
                        self.my_vertices.add(&ve);
                        builder_add_edge_vertex(&mut my_edge, &ve);
                        let mut ve_out = ve.clone();
                        transfert(&ee, &mut my_edge, &ve, &mut ve_out);
                        builder_add_edge_vertex(&mut my_edge, &ve_out);
                    }
                }

                self.my_edge = my_edge;
            }
            // Make a decision about the orientation of the edge
            // (OCCT L370-379).
            if (forward == reverse && the_e.orientation == Orientation::Reversed) || (reverse && !forward)
            {
                self.my_edge = crate::brep_algo::tool::reversed(&self.my_edge);
            }
        }

        // add myEdge to myShape (OCCT L382-384).
        let my_edge = self.my_edge.clone();
        builder_add_wire_edge(&mut self.my_shape, &my_edge);
        builder_set_closed(&mut self.my_shape, false);

        // Initialize VF, VL (OCCT L386-444).
        if init {
            let (vf, vl) = top_exp_vertices_wire(&self.my_shape);
            self.vf = vf.unwrap_or_else(Shape::null);
            self.vl = vl.unwrap_or_else(Shape::null);
        } else {
            if self.my_error == WireError::WireDone {
                // Update only
                let (v1, v2) = top_exp_vertices_raw(&self.my_edge);
                let v1 = v1.unwrap_or_else(Shape::null);
                let v2 = v2.unwrap_or_else(Shape::null);
                let mut v_ref = Shape::null();
                if v1.is_same(&self.my_vertex) {
                    v_ref = v2;
                } else if v2.is_same(&self.my_vertex) {
                    v_ref = v1;
                } else {
                    self.my_error = WireError::NonManifoldWire;
                }

                if self.vf.is_same(&self.vl) {
                    // Particular case: it is required to control the
                    // orientation.
                } else if self.vf.is_same(&self.my_vertex) {
                    self.vf = v_ref;
                } else if self.vl.is_same(&self.my_vertex) {
                    self.vl = v_ref;
                } else {
                    self.my_error = WireError::NonManifoldWire;
                }
            }
            if self.my_error == WireError::NonManifoldWire {
                self.vf = Shape::null();
                self.vl = Shape::null();
            }
        }

        // Test myShape is closed (OCCT L445-449).
        if !self.vf.is_null() && !self.vl.is_null() && self.vf.is_same(&self.vl) {
            builder_set_closed(&mut self.my_shape, true);
        }

        // OCCT L451-452.
        self.my_error = WireError::WireDone;
    }

    /// OCCT BRepBuilderAPI_Command::NotDone() — sets myDone to false; the
    /// rcad carrier reads the failure through `my_error` (arch. diff. #1).
    fn not_done(&mut self) {}

    /// OCCT TopoDS_Shape::Closed() read on myShape (BRep_Builder::MakeWire
    /// initializes the flag; TopoDS_Shape.hxx L216-223).
    fn my_shape_closed(&self) -> bool {
        crate::brep_algo::tool::shape_is_closed(&self.my_shape)
    }
}

/// OCCT BRep_Builder::Transfert(Ein, Eout, Vin, Vout) (BRep_Builder.cxx
/// L1457-1465).
fn transfert(the_e_in: &Shape, the_e_out: &mut Shape, the_v_in: &Shape, the_v_out: &mut Shape) {
    let tol = brep_tool_tolerance(the_v_in);
    let par_in = brep_tool_parameter(the_v_in, the_e_in);
    builder_update_vertex_parameter(the_v_out, par_in, the_e_out, tol);
}

/// OCCT BRep_Builder::UpdateVertex(V, Par, E, Tol) (BRep_Builder.cxx
/// L1220-1310) — every curve representation of the edge takes Par as First
/// (FORWARD) or Last (REVERSED) and the vertex tolerance is raised to Tol.
/// The rcad storage of the vertex parameter is TEdgeData::vertex_params
/// keyed by the vertex TShape pointer (the BRep_Tool::Parameter(V, E)
/// vehicle, bop/ds/mod.rs L1213-1235).
pub(crate) fn builder_update_vertex_parameter(
    the_v: &mut Shape,
    the_par: f64,
    the_e: &mut Shape,
    the_tol: f64,
) {
    if the_par >= rcad_kernel::precision::INFINITE_VALUE
        || the_par <= -rcad_kernel::precision::INFINITE_VALUE
    {
        // Standard_DomainError ("BRep_Builder::Infinite parameter", cxx L1227).
        panic!("Standard_DomainError");
    }

    // Search the vertex in the edge (OCCT L1240-1276).
    let mut ori = Orientation::Internal;
    let forward_edge = oriented(the_e, Orientation::Forward);
    let itv = sub_shapes(&forward_edge);

    // if the edge has no vertices and is degenerated use the vertex
    // orientation (OCCT L1249-1252).
    let degenerated = match the_e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    };
    if itv.is_empty() && degenerated {
        ori = the_v.orientation;
    }

    for vcur in &itv {
        if the_v.is_same(vcur) {
            ori = vcur.orientation;
            if ori == the_v.orientation {
                break;
            }
        }
    }

    if let TShape::Edge(ed) = std::sync::Arc::make_mut(&mut the_e.data) {
        if ori == Orientation::Forward {
            ed.range[0] = the_par;
            for (_k, entry) in ed.pcurves.iter_mut() {
                entry.1 = the_par;
            }
            for cr in ed.representations.iter_mut() {
                if let rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
                    range,
                    ..
                } = cr
                {
                    range[0] = the_par;
                }
            }
        } else if ori == Orientation::Reversed {
            ed.range[1] = the_par;
            for (_k, entry) in ed.pcurves.iter_mut() {
                entry.2 = the_par;
            }
            for cr in ed.representations.iter_mut() {
                if let rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
                    range,
                    ..
                } = cr
                {
                    range[1] = the_par;
                }
            }
        } else {
            // OCCT L1283-1291: the point representations of TV
            // (UpdatePoints). Architecture difference: the rcad
            // TVertexData::points carries pool INDICES (curve: usize), not the
            // OCCT Geom_Curve handle, and no rcad reader of BRep_Tool::Parameter
            // (V, E) consumes it — the rcad carrier of the vertex parameter is
            // TEdgeData::vertex_params, written below
            // (bop/ds/mod.rs L1213-1235). The point-representation write is
            // therefore not reproduced here.
        }
        ed.vertex_params.insert(the_v.ptr_id(), the_par);
    }

    // TV->UpdateTolerance(Tol); TE->Modified(true) (OCCT L1305-1307).
    crate::brep_algo::tool::builder_update_vertex_tol(the_v, the_tol);
}
