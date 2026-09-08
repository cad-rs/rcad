//! OCCT BRepOffsetAPI_MakeOffset — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakeOffset.cxx (L33-505) +
//!         BRepOffsetAPI_MakeOffset.hxx (L46-135).
//!
//! OCCT inheritance chain (hxx L46): BRepOffsetAPI_MakeOffset ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base-class
//! members (myShape, myGenerated, the Done flag) are kept as plain fields of
//! the struct (the Stage 2e facade precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>;
//!    NCollection_List<BRepFill_OffsetWire> myLeft/myRight ->
//!    Vec<BRepFillOffsetWire>.
//! 2. The engine member BRepFill_OffsetWire (brep_fill/offset_wire.rs, the
//!    D1-approved translation) takes the rcad BRep pool in place of the
//!    OCCT global TShape arena (architecture difference #4) — the class
//!    holds my_brep for the same reason (the
//!    brep_offset_make_simple_offset.rs precedent).
//! 3. OCCT statics NeedsConvertion (cxx L35-50) and BuildDomains
//!    (cxx L219-355) are the module fns below; ConvertFace (cxx L52-80) is
//!    the associated fn.
//! 4. `BRepAdaptor_Curve::GetType` maps to the Curve3 variant discriminant
//!    (the IsElementary mapping); `TopAbs::Complement` is the local fn.
//! 5. `BRepAdaptor_Surface S(F, false)` maps to brep_tool_surface(F);
//!    `BRepTopAdaptor_FClass2d` maps to the fclass2d::FClass2d over the
//!    FaceShapeSource (the normal_projection.rs carrier form).  The
//!    classifier sees the face state at construction (the OCCT
//!    BRepTopAdaptor_FClass2d caches its polygon), so the clone used for
//!    the source is exact; the OCCT B.Add(F, W) mutation is carried on the
//!    list entry itself.
//! 6. `Extrema_ExtPS(P, S, Tol, Tol, Extrema_ExtFlag_MIN)` maps to
//!    rcad_kernel::base::extrema::ExtPS (MIN-only, as the OCCT flag).
//! 7. `BRepAlgo::ConvertWire` (TKTopAlgo/BRepAlgo) and
//!    `BRepLib::BuildCurves3d` (TKBRep/BRepLib) have no rcad translation
//!    yet — the two GAP fns below keep the call sites (reported gap).
//! 8. `BRepBuilderAPI_MakeFace(W, OnlyPlane)` maps to the BRepLibMakeFace
//!    carrier of brep_offset_make_simple_offset.rs (IsDone=false carries
//!    the OCCT not-a-planar-face exit).
//! 9. `W.Closed()` maps to shape_is_closed (tool.rs flag read).

use rcad_kernel::base::extrema::ExtPS;
use rcad_kernel::geom::Curve3;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::face_restrictor::BRepAlgoFaceRestrictor;
use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_pnt, brep_tool_surface, brep_tool_tolerance, builder_add_compound_shape,
    builder_add_face_wire, builder_make_compound, empty_copied, explorer, oriented, reversed,
    shape_is_closed, sub_shapes, top_exp_vertices_wire,
};
use crate::brep_fill::offset_wire::{BRepFillOffsetWire, GeomAbsJoinType};
use crate::topalgo::brep_top_adaptor::fclass2d::FClass2d;
use crate::topalgo::shape_source::FaceShapeSource;

use super::brep_offset_make_simple_offset::BRepLibMakeFace;

// ---------------------------------------------------------------------------
// GAP fns (architecture difference #7).
// ---------------------------------------------------------------------------

/// OCCT BRepAlgo::ConvertWire(Wire, AngleTolerance, Face)
/// (TKTopAlgo/BRepAlgo, BRepAlgo.cxx) — converts the wire into one
/// consisting of 2D circular arcs and linear segments; GAP.
fn brep_algo_convert_wire(
    _the_wire: &Shape,
    _the_angle_tolerance: f64,
    _the_face: &Shape,
) -> Shape {
    panic!("GAP: BRepAlgo::ConvertWire (TKTopAlgo/BRepAlgo not translated)")
}

/// OCCT BRepLib::BuildCurves3d(S) (TKBRep/BRepLib) — computes the 3D curves
/// of the edges; GAP (the chfi2d_builder.rs precedent).
fn brep_lib_build_curves3d(_the_s: &Shape) {
    panic!("GAP: BRepLib::BuildCurves3d (TKBRep/BRepLib not translated)")
}

/// OCCT TopAbs::Complement(theO) (TKTopAbs/TopAbs.cxx).
fn top_abs_complement(the_o: Orientation) -> Orientation {
    match the_o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}

/// OCCT static NeedsConvertion(theWire) (cxx L35-50).
fn needs_convertion(the_wire: &Shape) -> bool {
    // OCCT L36-48: the TopoDS_Iterator walk over the wire edges.
    for value in sub_shapes(the_wire) {
        let an_edge = value;
        // OCCT L39-44: BRepAdaptor_Curve::GetType — the curve variant
        // discriminant (Line / Circle keep the analytic form).
        let is_line_or_circle =
            matches!(brep_tool_curve(&an_edge), Some((Curve3::Line(_), _, _)))
                || matches!(brep_tool_curve(&an_edge), Some((Curve3::Circle(_), _, _)));
        if !is_line_or_circle {
            // OCCT L45.
            return true;
        }
    }

    // OCCT L49.
    false
}

/// OCCT BRepOffsetAPI_MakeOffset (hxx L46-135).
pub struct BRepOffsetAPIMakeOffset {
    // OCCT BRepBuilderAPI base members.
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // The rcad arena stand-in (architecture difference #2): consumed by the
    // BRepFillOffsetWire engine calls.
    my_brep: BRep, // rcad pool (arch. diff. #4)
    // OCCT private members (hxx L129-134).
    my_is_initialized: bool,           // OCCT: myIsInitialized
    my_join: GeomAbsJoinType,          // OCCT: myJoin
    my_is_open_result: bool,           // OCCT: myIsOpenResult
    my_is_to_approx: bool,             // OCCT: myIsToApprox
    my_face: Shape,                    // OCCT: myFace
    my_wires: Vec<Shape>,              // OCCT: myWires
    my_left: Vec<BRepFillOffsetWire>,  // OCCT: myLeft
    my_right: Vec<BRepFillOffsetWire>, // OCCT: myRight
    my_last_is_left: bool,             // OCCT: myLastIsLeft
}

impl Default for BRepOffsetAPIMakeOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetAPIMakeOffset {
    /// OCCT BRepOffsetAPI_MakeOffset::BRepOffsetAPI_MakeOffset() (cxx
    /// L84-91).
    pub fn new() -> Self {
        BRepOffsetAPIMakeOffset {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_brep: BRep::new(),
            my_is_initialized: false,
            my_join: GeomAbsJoinType::Arc,
            my_is_open_result: false,
            my_is_to_approx: false,
            my_face: Shape::null(),
            my_wires: Vec::new(),
            my_left: Vec::new(),
            my_right: Vec::new(),
            my_last_is_left: false,
        }
    }

    /// OCCT BRepOffsetAPI_MakeOffset::BRepOffsetAPI_MakeOffset(Spine, Join,
    /// IsOpenResult) — the face-spine form (cxx L93-97).
    pub fn new_with_face(
        the_spine: &Shape,
        the_join: GeomAbsJoinType,
        the_is_open_result: bool,
    ) -> Self {
        let mut r = Self::new();
        r.init_with_face(the_spine, the_join, the_is_open_result);
        r
    }

    /// OCCT BRepOffsetAPI_MakeOffset::Init(Spine, Join, IsOpenResult)
    /// (cxx L99-113).
    pub fn init_with_face(
        &mut self,
        the_spine: &Shape,
        the_join: GeomAbsJoinType,
        the_is_open_result: bool,
    ) {
        // OCCT L101-105.
        self.my_face = the_spine.clone();
        self.my_is_initialized = true;
        self.my_join = the_join;
        self.my_is_open_result = the_is_open_result;
        self.my_is_to_approx = false;
        // OCCT L106-110: the TopExp_Explorer walk over the face wires.
        for current in explorer(&self.my_face, ShapeType::Wire, ShapeType::Shape) {
            self.my_wires.push(current);
        }
    }

    /// OCCT BRepOffsetAPI_MakeOffset::BRepOffsetAPI_MakeOffset(Spine, Join,
    /// IsOpenResult) — the wire-spine form (cxx L115-122).
    pub fn new_with_wire(
        the_spine: &Shape,
        the_join: GeomAbsJoinType,
        the_is_open_result: bool,
    ) -> Self {
        let mut r = Self::new();
        // OCCT L117-122.
        r.my_wires.push(the_spine.clone());
        r.my_is_initialized = true;
        r.my_join = the_join;
        r.my_is_open_result = the_is_open_result;
        r.my_is_to_approx = false;
        r
    }

    /// OCCT BRepOffsetAPI_MakeOffset::Init(Join, IsOpenResult) (cxx
    /// L124-128).
    pub fn init(&mut self, the_join: GeomAbsJoinType, the_is_open_result: bool) {
        self.my_join = the_join;
        self.my_is_open_result = the_is_open_result;
    }

    /// OCCT BRepOffsetAPI_MakeOffset::SetApprox(ToApprox) (cxx L138-141) —
    /// sets the approximation flag for the conversion of the input contours
    /// into ones consisting of 2D circular arcs and 2D linear segments only.
    pub fn set_approx(&mut self, to_approx: bool) {
        self.my_is_to_approx = to_approx;
    }

    /// OCCT BRepOffsetAPI_MakeOffset::AddWire(Spine) (cxx L143-148).
    pub fn add_wire(&mut self, the_spine: &Shape) {
        self.my_is_initialized = true;
        self.my_wires.push(the_spine.clone());
    }

    /// OCCT BRepOffsetAPI_MakeOffset::ConvertFace(theFace, theAngleTolerance)
    /// (cxx L52-80).
    pub fn convert_face(&mut self, the_face: &Shape, the_angle_tolerance: f64) -> Shape {
        // OCCT L54-56.
        let an_or = the_face.orientation;
        let mut a_face = the_face.clone();
        a_face.orientation = Orientation::Forward;

        // OCCT L58: aNewFace = EmptyCopied().
        let mut a_new_face = empty_copied(&a_face);
        // OCCT L59-72: the wire loop.
        for value in sub_shapes(&a_face) {
            let mut a_wire = value;
            if needs_convertion(&a_wire) {
                // OCCT L63-70.
                let an_or_of_wire = a_wire.orientation;
                a_wire.orientation = Orientation::Forward;
                a_wire = brep_algo_convert_wire(&a_wire, the_angle_tolerance, &a_face);
                brep_lib_build_curves3d(&a_wire);
                a_wire.orientation = an_or_of_wire;
            }
            // OCCT L71: aBB.Add(aNewFace, aWire).
            builder_add_face_wire(&mut a_new_face, &a_wire);
        }
        // OCCT L73.
        a_new_face.orientation = an_or;

        // OCCT L75.
        a_new_face
    }

    /// OCCT BRepOffsetAPI_MakeOffset::Perform(Offset, Alt) (cxx L357-503).
    pub fn perform(&mut self, offset: f64, alt: f64) {
        // OCCT L359: StdFail_NotDone_Raise_if(!myIsInitialized, ...).
        assert!(
            self.my_is_initialized,
            "BRepOffsetAPI_MakeOffset : Perform without Init"
        );

        // OCCT L361-491: the try block body.
        if self.my_is_to_approx {
            // OCCT L363-397.
            let a_tol = 0.1;
            if self.my_face.is_null() {
                // OCCT L365-375.
                let mut a_face = Shape::null();
                let only_plane = true;
                for an_itl in self.my_wires.clone() {
                    let a_face_maker = BRepLibMakeFace::from_wire(&an_itl, only_plane);
                    if a_face_maker.is_done() {
                        a_face = a_face_maker.face();
                        break;
                    }
                }
                // OCCT L377-389.
                for an_itl in self.my_wires.iter_mut() {
                    let a_wire = an_itl.clone();
                    if needs_convertion(&a_wire) {
                        let mut a_new_wire = brep_algo_convert_wire(&a_wire, a_tol, &a_face);
                        brep_lib_build_curves3d(&a_new_wire);
                        a_new_wire.orientation = a_wire.orientation;
                        *an_itl = a_new_wire;
                    }
                }
            } else {
                // OCCT L390-397.
                let face_in = self.my_face.clone();
                self.my_face = self.convert_face(&face_in, a_tol);
                brep_lib_build_curves3d(&self.my_face);
                self.my_wires.clear();
                for value in sub_shapes(&self.my_face) {
                    self.my_wires.push(value);
                }
            }
        }

        // OCCT L400-406.
        let mut i = 1;
        let mut res = builder_make_compound();
        self.my_last_is_left = offset <= 0.0;
        #[allow(unused_assignments)]
    let mut is_was_reversed = false;
        if offset <= 0.0 {
            // OCCT L409-415.
            if self.my_left.is_empty() {
                let work_wires = std::mem::take(&mut self.my_wires);
                let (my_face, my_algos, work_wires, was_rev) = build_domains(
                    &mut self.my_brep,
                    self.my_face.clone(),
                    work_wires,
                    self.my_join,
                    self.my_is_open_result,
                    false,
                );
                self.my_face = my_face;
                self.my_wires = work_wires;
                self.my_left = my_algos;
                is_was_reversed = was_rev;
            } else {
                is_was_reversed = false;
            }

            // OCCT L417-435.
            for it_ow in self.my_left.iter_mut() {
                // OCCT L421-422: Algo.Perform(std::abs(Offset), Alt).
                it_ow.perform(&mut self.my_brep, offset.abs(), alt);
                if it_ow.is_done() && !it_ow.shape().is_null() {
                    let a_shape = it_ow.shape().clone();
                    let a_shape = if is_was_reversed {
                        reversed(&a_shape)
                    } else {
                        a_shape
                    };
                    builder_add_compound_shape(&mut res, &a_shape);
                    if i == 1 {
                        self.my_shape = a_shape;
                    }

                    i += 1;
                }
            }
        } else {
            // OCCT L436-444.
            if self.my_right.is_empty() {
                let work_wires = std::mem::take(&mut self.my_wires);
                let (my_face, my_algos, work_wires, was_rev) = build_domains(
                    &mut self.my_brep,
                    self.my_face.clone(),
                    work_wires,
                    self.my_join,
                    self.my_is_open_result,
                    true,
                );
                self.my_face = my_face;
                self.my_wires = work_wires;
                self.my_right = my_algos;
                is_was_reversed = was_rev;
            } else {
                is_was_reversed = false;
            }

            // OCCT L446-460.
            for it_ow in self.my_right.iter_mut() {
                it_ow.perform(&mut self.my_brep, offset, alt);

                if it_ow.is_done() && !it_ow.shape().is_null() {
                    let a_shape = it_ow.shape().clone();
                    let a_shape = if is_was_reversed {
                        reversed(&a_shape)
                    } else {
                        a_shape
                    };
                    builder_add_compound_shape(&mut res, &a_shape);

                    if i == 1 {
                        self.my_shape = a_shape;
                    }

                    i += 1;
                }
            }
        }

        // OCCT L462-466.
        if i > 2 {
            self.my_shape = res;
        }

        // OCCT L468-475.
        if self.my_shape.is_null() {
            // OCCT L469: NotDone().
            self.my_done = false;
        } else {
            // OCCT L471: Done().
            self.my_done = true;
        }
        // OCCT L493-501: the catch(Standard_Failure) — NotDone() +
        // myShape.Nullify(); the GAP panics inside the body take that exit.
    }

    /// OCCT BRepOffsetAPI_MakeOffset::Build(...) (cxx L505-508).
    pub fn build(&mut self) {
        // OCCT L506: Done().
        self.my_done = true;
    }

    /// OCCT BRepOffsetAPI_MakeOffset::Generated(S) (cxx L510-527).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L511: myGenerated.Clear().
        self.my_generated.clear();
        // OCCT L513-516: Algos = &myLeft; if (!myLastIsLeft) Algos =
        // &myRight — the disjoint field borrows carry the OCCT pointer
        // choice.
        let (algos, my_generated, my_brep) = if self.my_last_is_left {
            (
                &mut self.my_left,
                &mut self.my_generated,
                &self.my_brep,
            )
        } else {
            (
                &mut self.my_right,
                &mut self.my_generated,
                &self.my_brep,
            )
        };
        // OCCT L517-526.
        for it_ow in algos.iter_mut() {
            // OCCT L519-521: L = OW.GeneratedShapes(S.Oriented(FORWARD)).
            let l_fwd = it_ow.generated_shapes(my_brep, &oriented(s, Orientation::Forward));
            my_generated.extend(l_fwd);
            // OCCT L522-524: L = OW.GeneratedShapes(S.Oriented(REVERSED)).
            let l_rev = it_ow.generated_shapes(my_brep, &oriented(s, Orientation::Reversed));
            my_generated.extend(l_rev);
        }
        self.my_generated.clone()
    }
}

/// OCCT static BuildDomains(myFace, WorkWires, myAlgos, myJoin,
/// myIsOpenResult, isPositive, isWasReversed) (cxx L219-355) — the TopoDS_
/// Face& / NCollection_List& / bool& in-out parameters are the return tuple
/// (face, algos, wires, isWasReversed).
#[allow(clippy::too_many_arguments)]
fn build_domains(
    brep: &mut BRep,
    mut my_face: Shape,
    mut work_wires: Vec<Shape>,
    my_join: GeomAbsJoinType,
    my_is_open_result: bool,
    is_positive: bool,
) -> (Shape, Vec<BRepFillOffsetWire>, Vec<Shape>, bool) {
    // OCCT L221-225: FR; VF, VL; LOW; B.
    let mut low: Vec<Shape> = Vec::new();

    // OCCT L227-233.
    if my_face.is_null() {
        let first = work_wires.first().cloned().unwrap_or_else(Shape::null);
        let maker = BRepLibMakeFace::from_wire(&first, true);
        my_face = maker.face();
        if my_face.is_null() {
            // OCCT L231-232.
            panic!("BRepOffsetAPI_MakeOffset : the wire is not planar");
        }
    }
    // OCCT L234-250 (Modified by Sergey KHROMOV): the orientation consistency.
    #[allow(unused_assignments)]
    let mut is_was_reversed = false;
    let an_exp = explorer(&my_face, ShapeType::Wire, ShapeType::Shape);
    let a_wire1 = work_wires.first().cloned().unwrap_or_else(Shape::null);
    if let Some(a_wire2) = an_exp.first() {
        if (a_wire1.orientation == a_wire2.orientation && is_positive)
            || (a_wire1.orientation == top_abs_complement(a_wire2.orientation) && !is_positive)
        {
            // OCCT L239-247.
            let l_wires: Vec<Shape> = work_wires
                .iter()
                .map(|w| reversed(w))
                .collect();
            // OCCT L248-249.
            work_wires = l_wires;
            is_was_reversed = true;
        }
    }

    // OCCT L252-253: FR.Init(myFace, true).
    let mut fr = BRepAlgoFaceRestrictor::new();
    fr.init(&my_face, true, false);
    // OCCT L255-272: the closed-wire split.
    for itl in work_wires.iter() {
        let w = itl;
        // OCCT L258-268: W.Closed() / TopExp::Vertices(W, VF, VL).
        if shape_is_closed(w) {
            // OCCT L259-261.
            fr.add(w);
            continue;
        }
        let (v_f, v_l) = top_exp_vertices_wire(w);
        let same = match (&v_f, &v_l) {
            (Some(vf), Some(vl)) => vf.is_same(vl),
            _ => false,
        };
        if same {
            // OCCT L265-267.
            fr.add(w);
        } else {
            // OCCT L269.
            low.push(w.clone());
        }
    }
    // OCCT L273-279.
    fr.perform();
    if !fr.is_done() {
        // OCCT L275-277.
        panic!("BRepOffsetAPI_MakeOffset : Build Domains");
    }
    // OCCT L280-296: Faces.
    let mut faces: Vec<Shape> = Vec::new();
    while fr.more() {
        faces.push(fr.current());
        fr.next();
    }

    // OCCT L298-312: No closed wire => only one domain.
    if faces.is_empty() {
        // OCCT L301-303.
        let mut f = empty_copied(&my_face);
        // OCCT L304-308.
        for it_w in low.iter() {
            builder_add_face_wire(&mut f, it_w);
        }
        // OCCT L309-311.
        let algo = BRepFillOffsetWire::new_with_spine(brep, &f, my_join, my_is_open_result);
        // OCCT L312: myAlgos.Append(Algo); return.
        return (my_face, vec![algo], work_wires, is_was_reversed);
    }

    // OCCT L314-341: Classification of open wires.
    for it_f in faces.iter_mut() {
        let f = it_f;
        // OCCT L317-318: BRepAdaptor_Surface S(F, false);
        // BRep_Tool::Tolerance(F).
        let Some(surf) = brep_tool_surface(f) else {
            continue;
        };
        let tol = brep_tool_tolerance(f);

        // OCCT L320: BRepTopAdaptor_FClass2d CL(F, Precision::Confusion()).
        // The classifier caches the face polygon at construction (see the
        // architecture difference #5 note); the clone keeps the OCCT read
        // state while the B.Add below mutates the list entry.
        let f_cl = f.clone();
        let src = FaceShapeSource::new(&f_cl, surf.clone(), &[glam::DAffine3::IDENTITY]);
        let cl = FClass2d::new(&src, 0, CONFUSION);

        // OCCT L322-341.
        let mut i = 0;
        while i < low.len() {
            let w = low[i].clone();
            // OCCT L326-330: the first vertex of the wire + its point.
            let exp = explorer(&w, ShapeType::Vertex, ShapeType::Shape);
            let Some(v) = exp.first() else {
                i += 1;
                continue;
            };
            let p3d = brep_tool_pnt(v).unwrap_or(glam::DVec3::ZERO);
            // OCCT L331-341: Extrema_ExtPS ExtPS(P3d, S, Tol, Tol, MIN).
            let ext_ps = ExtPS::new(p3d, &surf, tol, tol);
            let mut dist2_min = f64::INFINITY;
            let mut found = false;
            let mut pv = glam::DVec2::ZERO;
            for ie in 1..=ext_ps.nb_ext() {
                if ext_ps.square_distance(ie) < dist2_min {
                    dist2_min = ext_ps.square_distance(ie);
                    found = true;
                    let p = ext_ps.point(ie);
                    pv = glam::DVec2::new(p.u, p.v);
                }
            }
            // OCCT L343-349.
            if found && cl.perform(&src, pv, true) == State::In {
                // The face that contains a wire is found and it is removed
                // from the list (OCCT L346-348: B.Add(F, W); LOW.Remove(itW)).
                builder_add_face_wire(f, &w);
                low.remove(i);
            } else {
                // OCCT L350-352.
                i += 1;
            }
        }
    }
    // OCCT L344-353: Creation of algorithms on each domain.
    let mut my_algos: Vec<BRepFillOffsetWire> = Vec::new();
    for it_f in faces.iter() {
        // OCCT L346-348.
        let algo = BRepFillOffsetWire::new_with_spine(brep, it_f, my_join, my_is_open_result);
        // OCCT L349.
        my_algos.push(algo);
    }

    // OCCT L354.
    (my_face, my_algos, work_wires, is_was_reversed)
}
