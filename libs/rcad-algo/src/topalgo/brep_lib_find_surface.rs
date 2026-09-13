//! OCCT BRepLib_FindSurface (TKTopAlgo/BRepLib — `BRepLib_FindSurface.hxx`
//! + `BRepLib_FindSurface.cxx` L51-621, 621 lines) — the planar-surface
//! finder over a shape: searches an existing shared surface of the edges,
//! then fits a least-squares plane over weighted sample points.
//!
//! Architecture bridges (the W1 shape_analysis/edge.rs numbering style):
//! 1. OCCT `BRep_Tool` / `TopExp` read the TShape graph through the shape
//!    document; the re-hosts below resolve TShape data through
//!    `rcad_kernel::topods::BRep` (the `topexp_explorer` re-host is the
//!    shape_build/brep_tool.rs precedent).
//! 2. `TopLoc_Location` -> `u32` (the `BRep.locations` table index,
//!    0 = identity; `IsEqual` -> `==`, `Identity()` -> 0).
//! 3. `Geom_Surface` / `Geom_Plane` handles -> `Surface3` / `Surface3::Plane`
//!    values (clone-on-handle); handle identity (`SS == mySurface`) -> the
//!    kernel `surface_same` value comparison; the
//!    `Geom_RectangularTrimmedSurface` basis walk -> `Surface3::Trimmed`.
//! 4. `gp_Pnt` / `gp_XYZ` / `gp_Vec` -> `DVec3`; `Geom_Plane::Coefficients`
//!    -> the local re-host below (the bop/int_tools/bean_face_intersector.rs
//!    precedent).
//! 5. `BRepAdaptor_Curve` -> the local re-host over the edge 3d curve +
//!    range; `NbIntervals(GeomAbs_C3)` -> the GeomAdaptor CN re-host (the
//!    brep_fill_sweep_c.rs precedent).
//! 6. `math_Matrix` / `math_Vector` / `math_Jacobi` -> the kernel 1-based
//!    `MatD` / `VecD` / `MathJacobi` (the same OCCT 1-based indexing).
//! 7. `BRepTools_WireExplorer` -> the brep_offset_inter2d.rs re-host (the
//!    edge walk follows the wire's own order; the OCCT explorer sorts edges
//!    by parametric connectivity on the face — for the natural-bounds tmp
//!    face of Is2DClosed the orders agree for non-seam wires).
//! 8. OCCT `OCC_CATCH_SIGNALS` + try/catch in Is2DClosed -> the rcad
//!    evaluation is panic-free; the OCCT `catch (Standard_Failure const&)
//!    { return false; }` exit is annotated at the early-return sites.

use rcad_kernel::math::math_jacobi::MathJacobi;
use rcad_kernel::math::{MatD, VecD};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, ShapeType, State, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::geom::{Curve2dEval, CurveEval, Surface3, SurfaceEval};

use crate::brep_algo::tool::brep_tool_tolerance;

// OCCT Standard_Real.hxx L146-151 (RealSmall) and gp.hxx L59-60:
// gp::Resolution() = RealSmall() = DBL_MIN.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;
// OCCT Standard_Real.hxx L182-185: RealLast() = DBL_MAX.
const REAL_LAST: f64 = f64::MAX;

// ---------------------------------------------------------------------------
// BRep_Tool / TopoDS re-hosts (architecture bridges #1-#3)
// ---------------------------------------------------------------------------

/// OCCT TopoDS_Shape::IsSame(S): same TShape and same Location.
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    std::sync::Arc::as_ptr(&a.data) == std::sync::Arc::as_ptr(&b.data)
        && a.location == b.location
}

/// The face surface registered in the pool under a TShape pointer (the
/// Geom_Surface handle identity stand-in; the shape_analysis/edge.rs
/// re-host).
fn face_surface_by_ptr(brep: &BRep, fptr: u64) -> Option<Surface3> {
    let ts = brep
        .tshapes
        .iter()
        .find(|ts| std::sync::Arc::as_ptr(ts) as u64 == fptr)?;
    match ts.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(E, C, S, L, f, l, I) — the index overload:
/// the I-th curve-on-surface representation of the edge; C is its pcurve and
/// (S, L) the paired surface + representation location.  Bridge #3: the
/// representation location is the pcurve-key hash (face-ptr, loc).
fn brep_tool_curve_on_surface_index(
    brep: &BRep,
    the_edge: &Shape,
    the_i: i32,
) -> (
    Option<rcad_kernel::geom::Curve2d>,
    Option<Surface3>,
    u32,
    f64,
    f64,
) {
    let reps = match the_edge.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => return (None, None, 0, 0.0, 0.0),
    };
    let mut k = 0i32;
    for r in reps {
        let (key, pc, range) = match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        k += 1;
        if k == the_i {
            let s = face_surface_by_ptr(brep, key.0);
            return (Some(pc.clone()), s, key.1, range[0], range[1]);
        }
    }
    (None, None, 0, 0.0, 0.0)
}

/// OCCT BRep_Tool::CurveOnSurface(E, S, L, f, l) — the stored (surface,
/// location) overload: matches the edge's curve representations by surface
/// value + location (the shape_analysis/edge.rs on-face matching model).
fn brep_tool_curve_on_surface_stored(
    brep: &BRep,
    the_edge: &Shape,
    the_surface: &Surface3,
    the_location: u32,
) -> Option<(rcad_kernel::geom::Curve2d, f64, f64)> {
    let reps = match the_edge.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => return None,
    };
    for r in reps {
        let (key, pc, range) = match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        let s = face_surface_by_ptr(brep, key.0);
        let matched = match s.as_ref() {
            Some(s) => rcad_kernel::topods::surface_same(s, the_surface) && key.1 == the_location,
            None => false,
        };
        if matched {
            return Some((pc.clone(), range[0], range[1]));
        }
    }
    None
}

/// OCCT TopExp::FirstVertex(E, CumOri = true) — the orientation-composed
/// first vertex (the loc_ope_wires_on_shape.rs re-host).
fn top_exp_first_vertex(the_edge: &Shape) -> Option<Shape> {
    let ed = match the_edge.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    if the_edge.orientation.compose(ed.first.orientation) == Orientation::Forward {
        return Some(ed.first.clone());
    }
    if the_edge.orientation.compose(ed.last.orientation) == Orientation::Forward {
        return Some(ed.last.clone());
    }
    None
}

/// OCCT TopExp::LastVertex(E, CumOri = true).
fn top_exp_last_vertex(the_edge: &Shape) -> Option<Shape> {
    let ed = match the_edge.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    if the_edge.orientation.compose(ed.last.orientation) == Orientation::Reversed {
        return Some(ed.last.clone());
    }
    if the_edge.orientation.compose(ed.first.orientation) == Orientation::Reversed {
        return Some(ed.first.clone());
    }
    None
}

// ---------------------------------------------------------------------------
// Statics of BRepLib_FindSurface.cxx
// ---------------------------------------------------------------------------

/// OCCT Controle(thePoints, thePlane) (cxx L51-69) — the maximal distance of
/// the points to the plane.
fn controle(the_points: &[glam::DVec3], the_plane: &Surface3) -> f64 {
    let mut df_max_dist = 0.0f64;
    let Surface3::Plane(the_pl) = the_plane else {
        return df_max_dist;
    };
    // Geom_Plane::Coefficients(a, b, c, d): n.x*X + n.y*Y + n.z*Z + d = 0
    // (the bean_face_intersector.rs re-host).
    let (a, b, c, d) = (
        the_pl.normal.x,
        the_pl.normal.y,
        the_pl.normal.z,
        -the_pl.normal.dot(the_pl.origin),
    );
    for xyz in the_points {
        let dist = (a * xyz.x + b * xyz.y + c * xyz.z + d).abs();
        if dist > df_max_dist {
            df_max_dist = dist;
        }
    }

    df_max_dist
}

/// OCCT Is2DConnected(theEdge1, theEdge2, theSurface, theLocation)
/// (cxx L76-99): true if the last vertex of theEdge1 coincides with the
/// first vertex of theEdge2 in the parametric space of theFace.
fn is_2d_connected(
    brep: &BRep,
    the_edge1: &Shape,
    the_edge2: &Shape,
    the_surface: &Surface3,
    the_location: u32,
) -> bool {
    // get 2D points
    let Some((a_curve1, f1, l1)) =
        brep_tool_curve_on_surface_stored(brep, the_edge1, the_surface, the_location)
    else {
        return false;
    };
    let p1 = Curve2dEval::point_at(
        &a_curve1,
        if the_edge1.orientation == Orientation::Forward {
            l1
        } else {
            f1
        },
    );
    let Some((a_curve2, f2, l2)) =
        brep_tool_curve_on_surface_stored(brep, the_edge2, the_surface, the_location)
    else {
        return false;
    };
    let p2 = Curve2dEval::point_at(
        &a_curve2,
        if the_edge2.orientation == Orientation::Forward {
            f2
        } else {
            l2
        },
    );

    // compare 2D points
    // GeomAdaptor_Surface aSurface(theSurface): UResolution + VResolution
    // (the kernel Surface3 re-hosts).
    let a_v = top_exp_first_vertex(the_edge2);
    let tol3d = a_v.as_ref().map(brep_tool_tolerance).unwrap_or(0.0);
    let tol2d = the_surface.u_resolution(tol3d) + the_surface.v_resolution(tol3d);
    let dist2 = p1.distance_squared(p2);
    dist2 < tol2d * tol2d
}

/// OCCT Is2DClosed(theShape, theSurface, theLocation) (cxx L107-155): true
/// if the edges of theShape form a closed wire in the parametric space of
/// theSurface.
fn is_2d_closed(
    brep: &mut BRep,
    the_shape: &Shape,
    the_surface: &Surface3,
    the_location: u32,
) -> bool {
    // OCC_CATCH_SIGNALS (bridge #8): the OCCT try/catch ->
    // catch (Standard_Failure const&) { return false; }.
    // get a wire theShape (OCCT L115-120).
    let a_wire_exp = crate::shhealing::shape_build::brep_tool::topexp_explorer(
        brep,
        the_shape,
        ShapeType::Wire,
    );
    let Some(a_wire) = a_wire_exp.into_iter().next() else {
        return false;
    };
    // a tmp face (OCCT L122): BRepLib_MakeFace(theSurface, PConfusion()).
    let a_tmp_face =
        BRepBuilder::new().make_face(brep, Some(the_surface.clone()), Shape::null());

    // check topological closeness using wire explorer, if the wire is not
    // closed the 1st and the last vertices of wire are different
    // (OCCT L126-132).
    let mut a_wire_explorer = crate::offset::brep_offset_inter2d::BRepToolsWireExplorer::new();
    a_wire_explorer.init(&a_wire, &a_tmp_face);
    if !a_wire_explorer.more() {
        return false;
    }
    // remember the 1st and the last edges of aWire.
    let a_first_edge = a_wire_explorer.current();
    let mut a_last_edge = a_first_edge.clone();
    // check if edges connected topologically (that is assured by
    // BRepTools_WireExplorer) are connected in 2D (OCCT L136-144).
    let mut a_prev_edge = a_first_edge.clone();
    a_wire_explorer.next();
    while a_wire_explorer.more() {
        a_last_edge = a_wire_explorer.current();
        if !is_2d_connected(brep, &a_prev_edge, &a_last_edge, the_surface, the_location) {
            return false;
        }
        a_prev_edge = a_last_edge.clone();
        a_wire_explorer.next();
    }
    // wire is closed if (1st vertex of aFisrtEdge) == (last vertex of
    // aLastEdge) in 2D (OCCT L147-149).
    let a_v1 = top_exp_first_vertex(&a_first_edge);
    let a_v2 = top_exp_last_vertex(&a_last_edge);
    let same_v = match (a_v1, a_v2) {
        (Some(v1), Some(v2)) => shape_is_same(&v1, &v2),
        _ => false,
    };
    same_v && is_2d_connected(brep, &a_last_edge, &a_first_edge, the_surface, the_location)
}

/// OCCT fillParams(theKnots, theDegree, theParMin, theParMax, theParams)
/// (cxx L178-214).
fn fill_params(
    the_knots: &[f64],
    the_degree: i32,
    the_par_min: f64,
    the_par_max: f64,
    the_params: &mut Vec<f64>,
) {
    let p_confusion = rcad_kernel::PCONFUSION;
    let mut a_prev_par = the_par_min;
    the_params.push(a_prev_par);

    let a_nb_p = std::cmp::max(the_degree, 1);

    // OCCT indexes NCollection_Array1 from 1: theKnots(i)/theKnots(i+1) ->
    // the_knots[i-1]/the_knots[i] (the 0-based slice).
    let mut i = 1usize;
    while i < the_knots.len() && the_knots[i] < the_par_max - p_confusion {
        if the_knots[i] < the_par_min + p_confusion {
            i += 1;
            continue;
        }

        let a_step = (the_knots[i] - the_knots[i - 1]) / a_nb_p as f64;
        for k in 1..=a_nb_p {
            let a_par = the_knots[i - 1] + k as f64 * a_step;
            if a_par > the_par_max - p_confusion {
                break;
            }

            if a_par > a_prev_par + p_confusion {
                the_params.push(a_par);
                a_prev_par = a_par;
            }
        }
        i += 1;
    }
    the_params.push(the_par_max);
}

/// OCCT fillPoints(theCurve, theParams, thePoints, theWeights)
/// (cxx L216-242).
fn fill_points(
    the_curve: &BRepAdaptorCurve,
    the_params: &[f64],
    the_points: &mut Vec<glam::DVec3>,
    the_weights: &mut Vec<f64>,
) {
    let mut a_dist_prev = 0.0f64;
    let mut a_dist_next;
    let mut a_p_prev = the_curve.value(the_params[0]);
    let mut a_p_next = a_p_prev;

    for i_p in 1..=the_params.len() {
        if i_p < the_params.len() {
            let a_param = the_params[i_p];
            a_p_next = the_curve.value(a_param);
            a_dist_next = a_p_prev.distance(a_p_next);
        } else {
            a_dist_next = 0.0;
        }

        the_points.push(a_p_prev);
        the_weights.push(a_dist_prev + a_dist_next);
        a_dist_prev = a_dist_next;
        a_p_prev = a_p_next;
    }
}

/// OCCT Standard_Real.hxx L240-247: Epsilon(theValue).
fn standard_real_epsilon(the_value: f64) -> f64 {
    if the_value >= 0.0 {
        let a_next = f64::from_bits(the_value.to_bits() + 1);
        a_next - the_value
    } else {
        let a_prev = f64::from_bits(the_value.to_bits() - 1);
        the_value - a_prev
    }
}

// ---------------------------------------------------------------------------
// BRepAdaptor_Curve re-host (architecture bridge #5)
// ---------------------------------------------------------------------------

/// OCCT BRepAdaptor_Curve(E) — the edge 3d curve adaptor with its range;
/// only the members Init consumes are carried.
#[derive(Debug, Clone)]
struct BRepAdaptorCurve {
    /// The underlying curve (BRepAdaptor_Curve::Curve).
    curve: rcad_kernel::geom::Curve3,
    /// OCCT myFirst.
    first: f64,
    /// OCCT myLast.
    last: f64,
}

impl BRepAdaptorCurve {
    /// OCCT BRepAdaptor_Curve(E): BRep_Tool::Curve(E, f, l).
    fn new(the_edge: &Shape) -> Option<Self> {
        match the_edge.data.as_ref() {
            TShape::Edge(ed) => ed.curve.as_ref().map(|c| BRepAdaptorCurve {
                curve: c.clone(),
                first: ed.range[0],
                last: ed.range[1],
            }),
            _ => None,
        }
    }

    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64 {
        self.first
    }

    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64 {
        self.last
    }

    /// OCCT Value(U).
    fn value(&self, the_u: f64) -> glam::DVec3 {
        CurveEval::point_at(&self.curve, the_u)
    }

    /// OCCT GetType() — the rcad Curve3 enum (the GeomAbs_CurveType
    /// dispatch).
    fn curve(&self) -> &rcad_kernel::geom::Curve3 {
        &self.curve
    }

    /// OCCT NbIntervals(GeomAbs_C3) — the GeomAdaptor CN re-host (the
    /// brep_fill_sweep_c.rs `geom_adaptor_nb_intervals_cn` precedent,
    /// GeomAdaptor_Curve.cxx L333-411): for the C3 (== CN) request the OCCT
    /// switch takes aCont = Degree and counts interior knot
    /// multiplicities > aCont.
    fn nb_intervals_c3(&self) -> i32 {
        match &self.curve {
            rcad_kernel::geom::Curve3::BSpline(bs) => {
                let a_cont = bs.degree as i32;
                let mut intervals = 0i32;
                let mut i = 0usize;
                let n = bs.knots.len();
                while i < n {
                    let mut j = i;
                    while j < n && (bs.knots[j] - bs.knots[i]).abs() < 1.0e-12 {
                        j += 1;
                    }
                    let mult = (j - i) as i32;
                    let k = bs.knots[i];
                    if k > self.first + 1.0e-12 && k < self.last - 1.0e-12 && mult > a_cont {
                        intervals += 1;
                    }
                    i = j;
                }
                intervals + 1
            }
            _ => 1,
        }
    }
}

// ---------------------------------------------------------------------------
// BRepLib_FindSurface
// ---------------------------------------------------------------------------

/// OCCT BRepLib_FindSurface (hxx L45-100) — finds a plane (or an existing
/// surface) on which the shape's edges lie.
pub struct BRepLibFindSurface {
    /// OCCT myTolerance.
    my_tolerance: f64,
    /// OCCT myTolReached.
    my_tol_reached: f64,
    /// OCCT isExisted.
    is_existed: bool,
    /// OCCT myLocation (bridge #2: the location table index, 0 = identity).
    my_location: u32,
    /// OCCT mySurface.
    my_surface: Option<Surface3>,
}

impl BRepLibFindSurface {
    /// OCCT BRepLib_FindSurface(S, Tol = -1, OnlyPlane = false) (cxx
    /// L168-174; the OnlyClosed argument takes its hxx default false).
    pub fn new(brep: &mut BRep, s: &Shape, tol: f64, only_plane: bool) -> Self {
        let mut this = BRepLibFindSurface {
            my_tolerance: 0.0,
            my_tol_reached: 0.0,
            is_existed: false,
            my_location: 0,
            my_surface: None,
        };
        this.init(brep, s, tol, only_plane, false);
        this
    }

    /// OCCT BRepLib_FindSurface(S, Tol, OnlyPlane, OnlyClosed) (cxx
    /// L168-174).
    pub fn new_closed(
        brep: &mut BRep,
        s: &Shape,
        tol: f64,
        only_plane: bool,
        only_closed: bool,
    ) -> Self {
        let mut this = BRepLibFindSurface {
            my_tolerance: 0.0,
            my_tol_reached: 0.0,
            is_existed: false,
            my_location: 0,
            my_surface: None,
        };
        this.init(brep, s, tol, only_plane, only_closed);
        this
    }

    /// OCCT Init(S, Tol, OnlyPlane, OnlyClosed) (cxx L248-579).
    pub fn init(
        &mut self,
        brep: &mut BRep,
        s: &Shape,
        tol: f64,
        only_plane: bool,
        only_closed: bool,
    ) {
        // OCCT L253-257.
        self.my_tolerance = tol;
        self.my_tol_reached = 0.0;
        self.is_existed = false;
        self.my_location = 0;
        self.my_surface = None;

        // compute the tolerance (OCCT L260-269)
        for e in
            crate::shhealing::shape_build::brep_tool::topexp_explorer(brep, s, ShapeType::Edge)
        {
            let t = brep_tool_tolerance(&e);
            if t > self.my_tolerance {
                self.my_tolerance = t;
            }
        }

        // search an existing surface (OCCT L272-276)
        let ex =
            crate::shhealing::shape_build::brep_tool::topexp_explorer(brep, s, ShapeType::Edge);
        let Some(e) = ex.first().cloned() else {
            return; // no edges ....
        };

        let mut i = 0i32;

        // iterate on the surfaces of the first edge (OCCT L286-342)
        loop {
            i += 1;
            let (a_pc, ss, _loc, _f, _l) = brep_tool_curve_on_surface_index(brep, &e, i);
            self.my_surface = ss;
            let _ = a_pc;
            if self.my_surface.is_none() {
                break;
            }
            // check the other edges (OCCT L295-321)
            for es in &ex {
                if !shape_is_same(&e, es) {
                    let mut ss: Option<Surface3> = None;
                    let mut j = 0i32;
                    loop {
                        j += 1;
                        let (_a_ppc, ss_j, loc_j, _ff, _ll) =
                            brep_tool_curve_on_surface_index(brep, es, j);
                        let Some(ss_val) = ss_j else {
                            break;
                        };
                        // OCCT L308: (SS == mySurface) && (L.IsEqual(myLocation))
                        let matched = match self.my_surface.as_ref() {
                            Some(ms) => {
                                rcad_kernel::topods::surface_same(&ss_val, ms)
                                    && loc_j == self.my_location
                            }
                            None => false,
                        };
                        if matched {
                            ss = Some(ss_val);
                            break;
                        }
                        // OCCT L312: SS.Nullify().
                        ss = None;
                    }

                    if ss.is_none() {
                        self.my_surface = None;
                        break;
                    }
                }
            }

            // if OnlyPlane, eval if mySurface is a plane (OCCT L324-331)
            if only_plane && self.my_surface.is_some() {
                if let Some(Surface3::Trimmed(ts)) = &self.my_surface {
                    // mySurface = down_cast<Geom_RectangularTrimmedSurface>(
                    //                mySurface)->BasisSurface();
                    self.my_surface = Some((*ts.basis).clone());
                }
                // mySurface = down_cast<Geom_Plane>(mySurface): the null
                // handle for a non-plane.
                if !matches!(self.my_surface.as_ref(), Some(Surface3::Plane(_))) {
                    self.my_surface = None;
                }
            }

            if self.my_surface.is_some() {
                // if S is e.g. the bottom face of a cylinder, mySurface can
                // be the lateral (cylindrical) face of the cylinder; reject
                // an improper mySurface (OCCT L335-341)
                if !only_closed
                    || is_2d_closed(
                        brep,
                        s,
                        self.my_surface.as_ref().unwrap(),
                        self.my_location,
                    )
                {
                    break;
                }
            }
        }

        if self.my_surface.is_some() {
            self.is_existed = true;
            return;
        }

        // no existing surface, search a plane
        // 07/02/02 akm vvv : (OCC157) changed algorithm
        //                    1. Collect the points along all edges of the shape
        //                       For each point calculate the WEIGHT = sum of
        //                       distances from neighboring points (_only_ same edge)
        //                    2. Minimizing the weighed sum of squared deviations
        //                       compute coefficients of the sought plane.

        let mut a_points: Vec<glam::DVec3> = Vec::new();
        let mut a_weight: Vec<f64> = Vec::new();

        // ======================= Step #1 (OCCT L362-423)
        for es in
            crate::shhealing::shape_build::brep_tool::topexp_explorer(brep, s, ShapeType::Edge)
        {
            // BRepAdaptor_Curve c(TopoDS::Edge(ex.Current())).
            let Some(c) = BRepAdaptorCurve::new(&es) else {
                continue;
            };

            let df_uf = c.first_parameter();
            let df_ul = c.last_parameter();
            // IsEqual(dfUf, dfUl) (Standard_Real.hxx L148-151).
            if (df_uf - df_ul).abs() < REAL_SMALL {
                // Degenerate
                continue;
            }
            let mut i_nb_points = 0i32;

            // Fill the parameters of the sampling points
            let mut a_params: Vec<f64> = Vec::new();
            match c.curve() {
                rcad_kernel::geom::Curve3::Bezier(g_c) => {
                    // GC->FirstParameter() / GC->LastParameter() - the
                    // natural Bezier range (the default domain).
                    let dom = CurveEval::default_domain(g_c);
                    let a_knots: [f64; 2] = [dom[0], dom[1]];
                    fill_params(
                        &a_knots,
                        g_c.control_points.len().saturating_sub(1) as i32,
                        df_uf,
                        df_ul,
                        &mut a_params,
                    );
                }
                rcad_kernel::geom::Curve3::BSpline(g_c) => {
                    // OCCT passes Geom_BSplineCurve::Knots() (the distinct
                    // knots); the rcad knots vector carries the expanded
                    // multiplicities — repeated entries produce a zero step
                    // and are skipped by the aPrevPar guard, so the emitted
                    // parameters are identical.
                    fill_params(&g_c.knots, g_c.degree as i32, df_uf, df_ul, &mut a_params);
                }
                rcad_kernel::geom::Curve3::Line(_) => {
                    // Two points on a straight segment
                    a_params.push(df_uf);
                    a_params.push(df_ul);
                }
                rcad_kernel::geom::Curve3::Circle(_)
                | rcad_kernel::geom::Curve3::Ellipse(_)
                | rcad_kernel::geom::Curve3::Hyperbola(_)
                | rcad_kernel::geom::Curve3::Parabola(_) => {
                    // Four points on other analytical curves
                    i_nb_points = 4;
                    // Put some points on other curves (the OCCT fallthrough)
                    if i_nb_points == 0 {
                        i_nb_points = 15 + c.nb_intervals_c3();
                    }

                    let a_bounds: [f64; 2] = [df_uf, df_ul];

                    fill_params(&a_bounds, i_nb_points - 1, df_uf, df_ul, &mut a_params);
                }
                _ => {
                    // Put some points on other curves
                    if i_nb_points == 0 {
                        i_nb_points = 15 + c.nb_intervals_c3();
                    }

                    let a_bounds: [f64; 2] = [df_uf, df_ul];

                    fill_params(&a_bounds, i_nb_points - 1, df_uf, df_ul, &mut a_params);
                }
            }

            // Add the points with weights to the sequences (OCCT L422)
            fill_points(&c, &a_params, &mut a_points, &mut a_weight);
        }

        if a_points.len() < 3 {
            return;
        }

        // ======================= Step #2 (OCCT L430-545)
        self.my_location = 0;
        let mut a_mat = MatD::new(3, 3);
        let mut a_vec = VecD::new(3);
        // Find the barycenter and normalize weights
        let mut df_max_weight = 0.0f64;
        let mut a_bary_center = glam::DVec3::ZERO;
        let mut df_sum_weight = 0.0f64;
        for (i_point, pnt) in a_points.iter().enumerate() {
            let df_w = a_weight[i_point];
            a_bary_center += df_w * *pnt;
            df_sum_weight += df_w;
            if df_w > df_max_weight {
                df_max_weight = df_w;
            }
        }
        a_bary_center /= df_sum_weight;

        // Fill the matrix and the right vector
        for (i_point, pnt) in a_points.iter().enumerate() {
            let p = *pnt - a_bary_center;
            let w = a_weight[i_point] / df_max_weight;
            a_mat.set(1, 1, a_mat.get(1, 1) + w * p.x * p.x);
            a_mat.set(1, 2, a_mat.get(1, 2) + w * p.x * p.y);
            a_mat.set(1, 3, a_mat.get(1, 3) + w * p.x * p.z);
            //
            a_mat.set(2, 2, a_mat.get(2, 2) + w * p.y * p.y);
            a_mat.set(2, 3, a_mat.get(2, 3) + w * p.y * p.z);
            //
            a_mat.set(3, 3, a_mat.get(3, 3) + w * p.z * p.z);
        }
        a_mat.set(2, 1, a_mat.get(1, 2));
        a_mat.set(3, 1, a_mat.get(1, 3));
        a_mat.set(3, 2, a_mat.get(2, 3));
        //
        let an_eignval = MathJacobi::new(&a_mat);
        let mut is_solved = an_eignval.is_done();
        let mut isol = 0usize;
        if is_solved {
            let an_e_vals = an_eignval.values();
            // We need vector with eigenvalue ~ 0.
            let mut an_e_min = REAL_LAST;
            let mut an_e_max = -REAL_LAST;
            for i in 1..=3usize {
                let an_e = an_e_vals.get(i).abs();
                if an_e_min > an_e {
                    an_e_min = an_e;
                    isol = i;
                }
                if an_e_max < an_e {
                    an_e_max = an_e;
                }
            }

            if isol == 0 {
                is_solved = false;
            } else {
                let eps = standard_real_epsilon(an_e_max);
                if an_e_min <= eps {
                    // anEignval.Vector(isol, aVec).
                    let v = an_eignval.vector(isol);
                    for i in 1..=3usize {
                        a_vec.set(i, v.get(i));
                    }
                } else {
                    // try using vector product of other axes
                    let mut ind = [0usize, 0usize];
                    for i in 1..=3usize {
                        if i == isol {
                            continue;
                        }
                        if ind[0] == 0 {
                            ind[0] = i;
                            continue;
                        }
                        if ind[1] == 0 {
                            ind[1] = i;
                        }
                    }
                    let a_vec1 = an_eignval.vector(ind[0]);
                    let a_vec2 = an_eignval.vector(ind[1]);
                    let a_v1 = glam::DVec3::new(a_vec1.get(1), a_vec1.get(2), a_vec1.get(3));
                    let a_v2 = glam::DVec3::new(a_vec2.get(1), a_vec2.get(2), a_vec2.get(3));
                    let a_n = a_v1.cross(a_v2);
                    a_vec.set(1, a_n.x);
                    a_vec.set(2, a_n.y);
                    a_vec.set(3, a_n.z);
                }
                // aVec.Norm2() < gp::Resolution().
                if glam::DVec3::new(a_vec.get(1), a_vec.get(2), a_vec.get(3)).length()
                    < REAL_SMALL
                {
                    is_solved = false;
                }
            }
        }

        if !is_solved {
            return;
        }
        // Removing very small values (OCCT L546-555)
        let a_max_v = a_vec
            .get(1)
            .abs()
            .max(a_vec.get(2).abs().max(a_vec.get(3).abs()));
        let eps = standard_real_epsilon(a_max_v);
        for i in 1..=3usize {
            if a_vec.get(i).abs() <= eps {
                a_vec.set(i, 0.0);
            }
        }
        let a_n = glam::DVec3::new(a_vec.get(1), a_vec.get(2), a_vec.get(3));
        let a_plane = Surface3::Plane(rcad_kernel::geom::Plane::new(a_bary_center, a_n));
        self.my_tol_reached = controle(&a_points, &a_plane);
        let a_weakness = 5.0f64;
        if self.my_tol_reached <= self.my_tolerance
            || (tol < 0.0 && self.my_tol_reached < self.my_tolerance * a_weakness)
        {
            self.my_surface = Some(a_plane.clone());
            // If S is wire, try to orient surface according to orientation
            // of wire. (OCCT L563-577)
            let wire_closed = match s.data.as_ref() {
                TShape::Wire(wd) => {
                    wd.flags & rcad_kernel::topods::tshape_flags::CLOSED != 0
                }
                _ => false,
            };
            if s.shape_type() == ShapeType::Wire && wire_closed {
                let a_w = s.clone();
                let a_tmp_face = BRepBuilder::new().make_face(
                    brep,
                    Some(a_plane.clone()),
                    Shape::null(),
                );
                BRepBuilder::new().add_to_face(brep, a_tmp_face.clone(), a_w);
                // BRepTopAdaptor_FClass2d FClass(aTmpFace, 0.) — the
                // FClass2dTopol re-host consumes an Arc<BRep> pool snapshot.
                let f_class = crate::topalgo::brep_top_adaptor::fclass2d_topol::FClass2dTopol::new(
                    std::sync::Arc::new(brep.clone()),
                    &a_tmp_face,
                    0.0,
                );
                if f_class.perform_infinite_point() == State::In {
                    let Surface3::Plane(the_pl) = &a_plane else {
                        return;
                    };
                    let a_norm = -the_pl.normal;
                    self.my_surface = Some(Surface3::Plane(rcad_kernel::geom::Plane::new(
                        the_pl.origin,
                        a_norm,
                    )));
                }
            }
        }
    }

    /// OCCT Found() (cxx L583-586).
    pub fn found(&self) -> bool {
        self.my_surface.is_some()
    }

    /// OCCT Surface() (cxx L590-593).
    pub fn surface(&self) -> Option<Surface3> {
        self.my_surface.clone()
    }

    /// OCCT Tolerance() (cxx L597-600).
    pub fn tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT ToleranceReached() (cxx L604-607).
    pub fn tolerance_reached(&self) -> f64 {
        self.my_tol_reached
    }

    /// OCCT Existed() (cxx L611-614).
    pub fn existed(&self) -> bool {
        self.is_existed
    }

    /// OCCT Location() (cxx L618-621).
    pub fn location(&self) -> u32 {
        self.my_location
    }
}
