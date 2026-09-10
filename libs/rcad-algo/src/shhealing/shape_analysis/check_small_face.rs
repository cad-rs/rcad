//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_CheckSmallFace` (`ShapeAnalysis_CheckSmallFace.hxx`
//! L17-216 + `.cxx` L1-1341).
//!
//! Checks the faces which are too small: spot faces, strip faces, pin
//! faces, twisted faces, and the splitting vertices.
//!
//! Architecture bridges (the W1-1/W1-5 conventions):
//! 1. `BRep` pool argument — `BRep_Tool` / `TopExp_Explorer` /
//!    `TopoDS_Iterator` / `BRep_Builder` read and mutate the TShape graph
//!    through `rcad_kernel::BRep`.
//! 2. `TopoDS_Iterator(F, CumOri = false)` -> `brep_tool::raw_subshapes`;
//!    `TopExp_Explorer` -> `brep_tool::topexp_explorer`; `TopExp::Vertices
//!    (E, V1, V2)` (CumOri = true) -> the local `top_exp_vertices` (the
//!    TopExp.cxx L182-210 semantics).
//! 3. `Geom_Curve` / `Geom_Surface` handles -> `Curve3` / `Surface3`
//!    values (cloned); `gp_Pnt` / `gp_Vec` -> `DVec3`.
//! 4. `Geom_TrimmedCurve(C, U1, U2, Sense = true)` -> `Curve3::Trimmed`
//!    (the rcad carrier stores the basis curve and the range; the
//!    FirstParameter/LastParameter reads are the stored range).
//! 5. `GeomAdaptor_Curve` -> `topalgo::brep_lib_validate_edge::
//!    GeomAdaptorCurve`; `GeomAdaptor_Surface::D1` -> `SurfaceEval::
//!    derivatives`.
//! 6. `NCollection_DataMap<TopoDS_Shape, ...>` -> `HashMap<(u64, u32), V>`
//!    keyed by the IsSame identity (TShape pointer + location).
//! 7. The OCCT uninitialised members (`myPrecision`, `myStatusPinFace`,
//!    `myStatusPinEdges` are not set by the OCCT constructor) get the
//!    deterministic rcad initialisers `myPrecision = -1` (the "not set"
//!    value driving the OCCT defaults) and the OK statuses.
//!
//! The OCCT literal quirks are kept as-is: `dv = (umax - umin) / nbint`
//! in CheckTwisted (cxx L1001), the `(M_PI - angle2) <= 0.001 &&
//! (M_PI - angle2) <= 0.01` test in CheckPinEdges (cxx L1337), and the
//! `return true;` (integer 1) of IsSpotFace for a face without wires
//! (cxx L161).
//
// The OCCT literal quirks (declared-but-unused variables, dead stores,
// unnecessary mut) are kept verbatim; the allowances silence the rustc
// noise without touching the statements.
#![allow(unused_variables, unused_assignments, unused_mut)]

use std::collections::HashMap;

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval as _, Surface3, SurfaceEval as _, TrimmedCurve3};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType, TShape};

use crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
use crate::shhealing::shape_analysis::wire::ShapeAnalysisWire;
use crate::shhealing::shape_analysis::wire_order::ShapeAnalysisWireOrder;
use crate::shhealing::shape_build::brep_tool::{
    builder_add, iter_subshapes, raw_subshapes, topexp_explorer,
};
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::topalgo::brep_lib_validate_edge::GeomAdaptorCurve;

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts (bridge #1).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Curve(edg, f, l) — the no-location variant.
fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .curve
            .as_ref()
            .map(|c| (c.clone(), ed.range[0], ed.range[1])),
        _ => None,
    }
}

/// OCCT BRep_Tool::Pnt(vtx).
fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(shape).
fn brep_tool_tolerance(the_s: &Shape) -> f64 {
    match the_s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Surface(fac, loc) — the surface and the location.
fn brep_tool_surface_loc(fac: &Shape) -> (Option<Surface3>, u32) {
    match fac.data.as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fac.location),
        _ => (None, 0),
    }
}

/// OCCT Geom_Curve::FirstParameter()/LastParameter() — the natural domain.
fn curve_first_parameter(c: &Curve3) -> f64 {
    rcad_kernel::geom::CurveEval::default_domain(c)[0]
}

fn curve_last_parameter(c: &Curve3) -> f64 {
    rcad_kernel::geom::CurveEval::default_domain(c)[1]
}

// ---------------------------------------------------------------------------
// TopExp re-hosts (CumOri = true).
// ---------------------------------------------------------------------------

/// OCCT TopExp::FirstVertex(E, CumOri = true) (TopExp.cxx L182-194).
fn top_exp_first_vertex(e: &Shape) -> Shape {
    let (v_first, v_last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => return Shape::null(),
    };
    let mut v = if e.orientation == Orientation::Reversed {
        v_last
    } else {
        v_first
    };
    v.orientation = e.orientation.compose(v.orientation);
    v
}

/// OCCT TopExp::LastVertex(E, CumOri = true) (TopExp.cxx L198-210).
fn top_exp_last_vertex(e: &Shape) -> Shape {
    let (v_first, v_last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => return Shape::null(),
    };
    let mut v = if e.orientation == Orientation::Reversed {
        v_first
    } else {
        v_last
    };
    v.orientation = e.orientation.compose(v.orientation);
    v
}

/// OCCT TopExp::Vertices(E, V1, V2, CumOri = true) (TopExp.cxx): the
/// composed-FORWARD child is Vfirst, the composed-REVERSED child is Vlast.
fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    let mut v_first = Shape::null();
    let mut v_last = Shape::null();
    let children: Vec<Shape> = match e.data.as_ref() {
        TShape::Edge(ed) => vec![ed.first.clone(), ed.last.clone()],
        _ => Vec::new(),
    };
    for a_v in &children {
        let mut v = a_v.clone();
        v.orientation = e.orientation.compose(v.orientation);
        match v.orientation {
            Orientation::Forward => v_first = v,
            Orientation::Reversed => v_last = v,
            _ => {}
        }
    }
    (v_first, v_last)
}

/// OCCT TopoDS_Shape::IsSame(S).
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.is_same(b)
}

// ---------------------------------------------------------------------------
// Static helpers (cxx L67-128, L793-844, L960-968).
// ---------------------------------------------------------------------------

/// OCCT static MinMaxPnt (cxx L67-112).
fn min_max_pnt(
    p: DVec3,
    nb: &mut i32,
    minx: &mut f64,
    miny: &mut f64,
    minz: &mut f64,
    maxx: &mut f64,
    maxy: &mut f64,
    maxz: &mut f64,
) {
    let (x, y, z) = (p.x, p.y, p.z);
    if *nb < 1 {
        *minx = x;
        *maxx = x;
        *miny = y;
        *maxy = y;
        *minz = z;
        *maxz = z;
    } else {
        if *minx > x {
            *minx = x;
        }
        if *maxx < x {
            *maxx = x;
        }
        if *miny > y {
            *miny = y;
        }
        if *maxy < y {
            *maxy = y;
        }
        if *minz > z {
            *minz = z;
        }
        if *maxz < z {
            *maxz = z;
        }
    }
    *nb += 1;
}

/// OCCT Precision::IsInfinite (Precision.hxx L350-353).
fn precision_is_infinite(r: f64) -> bool {
    r.abs() >= 2e100
}

/// OCCT static MinMaxSmall (cxx L114-128).
#[allow(clippy::too_many_arguments)]
fn min_max_small(
    minx: f64,
    miny: f64,
    minz: f64,
    maxx: f64,
    maxy: f64,
    maxz: f64,
    toler: f64,
) -> bool {
    let dx = maxx - minx;
    let dy = maxy - miny;
    let dz = maxz - minz;

    (dx <= toler || precision_is_infinite(dx))
        && (dy <= toler || precision_is_infinite(dy))
        && (dz <= toler || precision_is_infinite(dz))
}

/// OCCT static IsoStat (cxx L793-823): the min-max status of one pole row
/// (uorv == 1: the U-rank row scanned over V; uorv == 2: the V-rank column
/// scanned over U).  Bridge: the rcad pole grid is 0-based
/// (`control_points[u][v]`), the OCCT Array2 is 1-based.
#[allow(clippy::too_many_arguments)]
fn iso_stat(
    poles: &[Vec<DVec3>],
    uorv: i32,
    rank: i32,
    tolpin: f64,
    toler: f64,
) -> i32 {
    let mut np = 0;
    // OCCT: i0 = (uorv == 1 ? poles.LowerCol() : poles.LowerRow()) — 1.
    let nb = if uorv == 1 {
        poles.get((rank - 1) as usize).map_or(0, |row| row.len())
    } else {
        poles.len()
    };
    let mut xmin = 0.;
    let mut ymin = 0.;
    let mut zmin = 0.;
    let mut xmax = 0.;
    let mut ymax = 0.;
    let mut zmax = 0.;
    for i in 1..=nb as i32 {
        let unp = if uorv == 1 {
            poles[(rank - 1) as usize][(i - 1) as usize]
        } else {
            poles[(i - 1) as usize][(rank - 1) as usize]
        };
        min_max_pnt(
            unp, &mut np, &mut xmin, &mut ymin, &mut zmin, &mut xmax, &mut ymax, &mut zmax,
        );
    }
    if min_max_small(xmin, ymin, zmin, xmax, ymax, zmax, tolpin) {
        return 0;
    }
    if min_max_small(xmin, ymin, zmin, xmax, ymax, zmax, toler) {
        return 1;
    }
    2
}

/// OCCT static CheckPoles (cxx L825-844).
fn check_poles(poles: &[Vec<DVec3>], uorv: i32, rank: i32) -> bool {
    let nb = if uorv == 1 {
        poles.get((rank - 1) as usize).map_or(0, |row| row.len())
    } else {
        poles.len()
    };
    for i in 1..nb as i32 {
        if uorv == 1 {
            let a = poles[(rank - 1) as usize][(i - 1) as usize];
            let b = poles[(rank - 1) as usize][i as usize];
            if a.distance(b) <= 1e-15 {
                return true;
            }
        } else {
            let a = poles[(i - 1) as usize][(rank - 1) as usize];
            let b = poles[i as usize][(rank - 1) as usize];
            if a.distance(b) <= 1e-15 {
                return true;
            }
        }
    }
    false
}

/// OCCT static TwistedNorm (cxx L960-968).
#[allow(clippy::too_many_arguments)]
fn twisted_norm(x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64) -> f64 {
    (x1 * x2) + (y1 * y2) + (z1 * z2)
}

/// OCCT gp_Vec::Angle(theOther) — the angle in [0, PI]; the OCCT
/// ConstructionError raise on a null magnitude maps to None (the
/// CheckPinEdges catch path).
fn gp_vec_angle(v1: DVec3, v2: DVec3) -> Option<f64> {
    if v1.length() == 0.0 || v2.length() == 0.0 {
        return None;
    }
    Some(v1.angle_between(v2))
}

/// OCCT ShapeAnalysis_CheckSmallFace (hxx L40-216).
pub struct ShapeAnalysisCheckSmallFace {
    /// OCCT `myStatus`.
    my_status: i32,
    /// OCCT `myStatusSpot`.
    my_status_spot: i32,
    /// OCCT `myStatusStrip`.
    my_status_strip: i32,
    /// OCCT `myStatusPin`.
    my_status_pin: i32,
    /// OCCT `myStatusTwisted`.
    my_status_twisted: i32,
    /// OCCT `myStatusSplitVert`.
    my_status_split_vert: i32,
    /// OCCT `myStatusPinFace` (uninitialised in OCCT; the rcad OK initial).
    my_status_pin_face: i32,
    /// OCCT `myStatusPinEdges` (uninitialised in OCCT; the rcad OK initial).
    my_status_pin_edges: i32,
    /// OCCT `myPrecision` (uninitialised in OCCT; the rcad -1 "not set"
    /// initial drives the OCCT 1e-4 / vertex-tolerance defaults).
    my_precision: f64,
}

impl Default for ShapeAnalysisCheckSmallFace {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisCheckSmallFace {
    /// OCCT ShapeAnalysis_CheckSmallFace() (cxx L58-65).
    pub fn new() -> Self {
        ShapeAnalysisCheckSmallFace {
            my_status: encode_status(ShapeExtendStatus::Ok),
            my_status_spot: encode_status(ShapeExtendStatus::Ok),
            my_status_strip: encode_status(ShapeExtendStatus::Ok),
            my_status_pin: encode_status(ShapeExtendStatus::Ok),
            my_status_twisted: encode_status(ShapeExtendStatus::Ok),
            my_status_split_vert: encode_status(ShapeExtendStatus::Ok),
            my_status_pin_face: encode_status(ShapeExtendStatus::Ok),
            my_status_pin_edges: encode_status(ShapeExtendStatus::Ok),
            my_precision: -1.0,
        }
    }

    /// OCCT SetPrecision (hxx inline) — the precision used by CheckPin /
    /// CheckTwisted / CheckSplittingVertices.
    pub fn set_precision(&mut self, precision: f64) {
        self.my_precision = precision;
    }

    /// OCCT Status (hxx L174-181): the status (True/False) of last Check.
    pub fn status(&self, status: ShapeExtendStatus) -> bool {
        status == ShapeExtendStatus::Ok
            || decode_status(self.my_status, status)
            || self.decode_spot(status)
            || self.decode_strip(status)
            || self.decode_pin(status)
            || self.decode_twisted(status)
            || self.decode_split_vert(status)
            || self.decode_pin_face(status)
            || self.decode_pin_edges(status)
    }

    /// The aggregate Spot status read.
    fn decode_spot(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_spot, status)
    }

    fn decode_strip(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_strip, status)
    }

    fn decode_pin(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_pin, status)
    }

    fn decode_twisted(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_twisted, status)
    }

    fn decode_split_vert(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_split_vert, status)
    }

    fn decode_pin_face(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_pin_face, status)
    }

    fn decode_pin_edges(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_pin_edges, status)
    }

    // -- IsSpotFace ---------------------------------------------------------

    /// OCCT IsSpotFace(F, spot, spotol, tol = -1.0) (cxx L132-230): returns
    /// 0 if not, 1 if yes, 2 if yes and all vertices are the same.
    pub fn is_spot_face(
        &self,
        brep: &mut BRep,
        f: &Shape,
        spot: &mut DVec3,
        spotol: &mut f64,
        tol: f64,
    ) -> i32 {
        let mut toler = tol;
        let mut tolv;
        //  Compute tolerance to get : from greatest tol of vertices
        //  In addition, also computes min-max of vertices
        //  To finally compare mini-max box with tolerance
        // gka Mar2000 Protection against faces without wires
        // but they occur due to bugs in the algorithm itself, it needs to be fixed
        let mut is_wir = false;
        // OCCT L146: TopoDS_Iterator itw(F, false) — CumOri = false.
        for child in raw_subshapes(brep, f) {
            if child.shape_type() != ShapeType::Wire {
                continue;
            }
            let w1 = child;
            if !w1.is_null() {
                is_wir = true;
                break;
            }
        }
        if !is_wir {
            // OCCT L161: `return true;` in the int function — the literal 1.
            return 1;
        }
        let mut nbv = 0;
        let mut minx = 0.0;
        let mut miny = 0.0;
        let mut minz = 0.0;
        let mut maxx = 2e100; // Precision::Infinite()
        let mut maxy = 2e100;
        let mut maxz = 2e100;
        let mut v0 = Shape::null();
        let mut same = true;
        for iv in topexp_explorer(brep, f, ShapeType::Vertex) {
            let v = iv;
            if v0.is_null() {
                v0 = v.clone();
            } else if same {
                if !shape_is_same(&v0, &v) {
                    same = false;
                }
            }

            let pnt = brep_tool_pnt(&v);
            min_max_pnt(
                pnt, &mut nbv, &mut minx, &mut miny, &mut minz, &mut maxx, &mut maxy, &mut maxz,
            );

            if tol < 0.0 {
                tolv = brep_tool_tolerance(&v);
                if tolv > toler {
                    toler = tolv;
                }
            }
        }

        //   Now, testing
        if !min_max_small(minx, miny, minz, maxx, maxy, maxz, toler) {
            return 0;
        }

        //   All vertices are confused
        //   Check edges (a closed edge may be a non-null length edge !)
        //   By picking intermediate point on each one
        for ie in topexp_explorer(brep, f, ShapeType::Edge) {
            let e = ie;
            let (c3d, cf, cl) = match brep_tool_curve(&e) {
                Some(v) => v,
                None => continue,
            };
            let debut = c3d.point_at(cf);
            let milieu = c3d.point_at((cf + cl) / 2.0);
            if debut.distance_squared(milieu) > toler * toler {
                return 0;
            }
        }

        *spot = DVec3::new(
            (minx + maxx) / 2.,
            (miny + maxy) / 2.,
            (minz + maxz) / 2.,
        );
        *spotol = maxx - minx;
        *spotol = spotol.max(maxy - miny);
        *spotol = spotol.max(maxz - minz);
        *spotol /= 2.0;

        if same {
            2
        } else {
            1
        }
    }

    /// OCCT CheckSpotFace(F, tol = -1.0) (cxx L234-255).
    pub fn check_spot_face(&mut self, brep: &mut BRep, f: &Shape, tol: f64) -> bool {
        let mut spot = DVec3::ZERO;
        let mut spotol = 0.0;
        let stat = self.is_spot_face(brep, f, &mut spot, &mut spotol, tol);
        if stat == 0 {
            return false;
        }
        match stat {
            1 => {
                self.my_status_spot = encode_status(ShapeExtendStatus::Done1);
            }
            2 => {
                self.my_status_spot = encode_status(ShapeExtendStatus::Done2);
            }
            _ => {}
        }
        true
    }

    // -- Strip support ------------------------------------------------------

    /// OCCT IsStripSupport(F, tol = -1.0) (cxx L259-343).
    pub fn is_strip_support(&mut self, brep: &mut BRep, f: &Shape, tol: f64) -> bool {
        let mut toler = tol;
        if toler < 0.0 {
            toler = 1.0e-07; // ?? better to compute tolerance zones
        }

        let (surf, _loc) = brep_tool_surface_loc(f);
        let Some(surf) = surf else {
            return false;
        };

        //  Checking on poles for bezier-bspline
        //  A more general way is to check Values by scanning ISOS (slower)

        let bs_poles: Option<Vec<Vec<DVec3>>> = match &surf {
            Surface3::BSpline(bs) => Some(bs.control_points.clone()),
            _ => None,
        };
        let bz_poles: Option<Vec<Vec<DVec3>>> = match &surf {
            Surface3::Bezier(bz) => Some(bz.control_points.clone()),
            _ => None,
        };

        // int stat = 2;  // 2 : small in V direction
        if bs_poles.is_some() || bz_poles.is_some() {
            let cbz = bz_poles.is_some();
            let poles = if cbz { bz_poles.as_ref() } else { bs_poles.as_ref() };
            let nbu = poles.map_or(0, |p| p.len());
            let nbv = poles
                .and_then(|p| p.first())
                .map_or(0, |row| row.len());
            let mut minx = 0.;
            let mut miny = 0.;
            let mut minz = 0.;
            let mut maxx = 0.;
            let mut maxy = 0.;
            let mut maxz = 0.;
            let mut issmall = true;

            for iu in 1..=nbu as i32 {
                //    for each U line, scan poles in V (V direction)
                let mut nb = 0;
                for iv in 1..=nbv as i32 {
                    let unp = poles.unwrap()[(iu - 1) as usize][(iv - 1) as usize];
                    min_max_pnt(
                        unp, &mut nb, &mut minx, &mut miny, &mut minz, &mut maxx, &mut maxy,
                        &mut maxz,
                    );
                }
                if !min_max_small(minx, miny, minz, maxx, maxy, maxz, toler) {
                    issmall = false;
                    break;
                } // small in V ?
            }
            if issmall {
                self.my_status_strip = encode_status(ShapeExtendStatus::Done2);
                return issmall; // OK, small in V
            }
            issmall = true;
            for iv in 1..=nbv as i32 {
                //    for each V line, scan poles in U (U direction)
                let mut nb = 0;
                for iu in 1..=nbu as i32 {
                    let unp = poles.unwrap()[(iu - 1) as usize][(iv - 1) as usize];
                    min_max_pnt(
                        unp, &mut nb, &mut minx, &mut miny, &mut minz, &mut maxx, &mut maxy,
                        &mut maxz,
                    );
                }
                if !min_max_small(minx, miny, minz, maxx, maxy, maxz, toler) {
                    issmall = false;
                    break;
                } // small in U ?
            }
            if issmall {
                self.my_status_strip = encode_status(ShapeExtendStatus::Done1);
                return issmall;
            } // OK, small in U
        }

        false
    }

    /// OCCT CheckStripEdges(E1, E2, tol, dmax) (cxx L347-433).
    pub fn check_strip_edges(
        &self,
        brep: &BRep,
        e1: &Shape,
        e2: &Shape,
        tol: f64,
        dmax: &mut f64,
    ) -> bool {
        let _ = brep;
        //  We have the topological configuration OK : 2 edges, 2 vertices
        //  But, are these two edges well confused ?
        let mut toler = tol;
        if tol < 0.0 {
            let tole = brep_tool_tolerance(e1) + brep_tool_tolerance(e2);
            if toler < tole / 2. {
                toler = tole / 2.;
            }
        }

        //   We project a list of points from each curve, on the opposite one,
        //   we check the distance
        let nbint = 10;

        let sac = ShapeAnalysisCurve;
        let mut u;
        *dmax = 0.0;
        let (c1, mut cf1, mut cl1) = match brep_tool_curve(e1) {
            Some(v) => v,
            None => return false,
        };
        let (c2, mut cf2, mut cl2) = match brep_tool_curve(e2) {
            Some(v) => v,
            None => return false,
        };
        cf1 = cf1.max(curve_first_parameter(&c1));
        cl1 = cl1.min(curve_last_parameter(&c1));
        // OCCT L380: Geom_TrimmedCurve(C1, cf1, cl1, true).
        let mut c1t = Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(c1.clone()),
            first: cf1,
            last: cl1,
        });
        // pdn protection against feature in Trimmed_Curve
        cf1 = trimmed_first_parameter(&c1t);
        cl1 = trimmed_last_parameter(&c1t);
        let mut cf2 = cf2.max(curve_first_parameter(&c2));
        let mut cl2 = cl2.min(curve_last_parameter(&c2));
        let mut c2t = Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(c2.clone()),
            first: cf2,
            last: cl2,
        });
        cf2 = trimmed_first_parameter(&c2t);
        cl2 = trimmed_last_parameter(&c2t);

        let mut cd1 = (cl1 - cf1) / nbint as f64;
        let mut cd2 = (cl2 - cf2) / nbint as f64;
        let mut f;
        let mut l;
        f = cf2;
        l = cl2;
        for numcur in 0..2 {
            u = cf1;
            if numcur != 0 {
                // OCCT L401-407: CC = C1T; C1T = C2T; C2T = CC; ...
                std::mem::swap(&mut c1t, &mut c2t);
                cd1 = cd2; // smh added replacing step and replacing first
                u = cf2; // parameter
                f = cf1;
                l = cl1;
            }
            for nump in 0..=nbint {
                let p1 = c1t.point_at(u);
                let mut p2 = DVec3::ZERO;
                let mut para = 0.0;
                // pdn Adaptor curve is used to avoid of enhancing of domain.
                let (basis, bfirst, blast) = trimmed_parts(&c2t);
                let mut gac = GeomAdaptorCurve::new(basis, bfirst, blast);
                let _ = &mut gac;
                let dist = sac.project_adaptor(&gac, p1, toler, &mut p2, &mut para, true);
                // pdn check if parameter of projection is in the domain of the edge.
                if para < f || para > l {
                    return false;
                }
                if dist > *dmax {
                    *dmax = dist;
                }
                if dist > toler {
                    return false;
                }
                u += cd1;
            }
        }
        *dmax < toler
    }

    /// OCCT FindStripEdges(F, E1, E2, tol, dmax) (cxx L437-512).
    pub fn find_strip_edges(
        &mut self,
        brep: &mut BRep,
        f: &Shape,
        e1: &mut Shape,
        e2: &mut Shape,
        tol: f64,
        dmax: &mut f64,
    ) -> bool {
        // OCCT L443-444: E1.Nullify(); E2.Nullify().
        *e1 = Shape::null();
        *e2 = Shape::null();
        let mut nb = 0;
        for ex in topexp_explorer(brep, f, ShapeType::Edge) {
            let e = ex;
            if nb == 1 && shape_is_same(&e, &*e1) {
                continue; // ignore seam edge
            }
            let (v1, v2) = top_exp_vertices(&e);
            let p1 = brep_tool_pnt(&v1);
            let p2 = brep_tool_pnt(&v2);
            let mut toler = tol;
            if toler <= 0.0 {
                toler = (brep_tool_tolerance(&v1) + brep_tool_tolerance(&v2)) / 2.;
            }

            //    Extremities
            let dist = p1.distance(p2);
            //    Middle point
            let cc = brep_tool_curve(&e);
            let mut is_null_length = true;
            if let Some((cc, cf, cl)) = cc.as_ref() {
                let pp = cc.point_at((cf + cl) / 2.);
                if pp.distance(p1) < toler && pp.distance(p2) < toler {
                    continue;
                }
                is_null_length = false;
            }
            if dist <= toler && is_null_length {
                continue; // smh
            }
            nb += 1;
            if nb == 1 {
                *e1 = e;
            } else if nb == 2 {
                *e2 = e;
            } else {
                return false;
            }
        }
        //   Now, check these two edge to define a strip !
        if !e1.is_null() && !e2.is_null() {
            if !self.check_strip_edges(brep, e1, e2, tol, dmax) {
                return false;
            } else {
                self.my_status_strip = encode_status(ShapeExtendStatus::Done3);
                return true;
            }
        }
        false
    }

    /// OCCT CheckSingleStrip(F, E1, E2, tol = -1.0) (cxx L516-651).
    pub fn check_single_strip(
        &mut self,
        brep: &mut BRep,
        f: &Shape,
        e1: &mut Shape,
        e2: &mut Shape,
        tol: f64,
    ) -> bool {
        let mut toler = tol;
        let mut minx = 0.0;
        let mut miny = 0.0;
        let mut minz = 0.0;
        let mut maxx = 0.0;
        let mut maxy = 0.0;
        let mut maxz = 0.0;

        // In this case, we have 2 vertices and 2 great edges. Plus possibly 2 small
        //    edges, one on each vertex
        let mut v1 = Shape::null();
        let mut v2 = Shape::null();
        let mut nb = 0;
        for itv in topexp_explorer(brep, f, ShapeType::Vertex) {
            let v = itv;
            if v1.is_null() {
                v1 = v;
            } else if shape_is_same(&v1, &v) {
                continue;
            } else if v2.is_null() {
                v2 = v;
            } else if shape_is_same(&v2, &v) {
                continue;
            } else {
                return false;
            }
        }

        // Checking edges
        nb = 0;
        let e1_snapshot = e1.clone();
        for ite in topexp_explorer(brep, f, ShapeType::Edge) {
            let e = ite;
            if nb == 1 && shape_is_same(&e, &e1_snapshot) {
                continue; // ignore seam edge
            }
            let (va, vb) = top_exp_vertices(&e);
            if tol < 0.0 {
                let mut tolv;
                tolv = brep_tool_tolerance(&va);
                if toler < tolv {
                    toler = tolv;
                }
                tolv = brep_tool_tolerance(&vb);
                if toler < tolv {
                    toler = tolv;
                }
            }

            //    Edge on same vertex : small one ?
            if shape_is_same(&va, &vb) {
                let c3d = if !brep_tool_degenerated(&e) {
                    brep_tool_curve(&e).map(|(c, cf, cl)| (c, cf, cl))
                } else {
                    None
                };
                let Some((c3d, cf, cl)) = c3d else {
                    continue; // DGNR
                };
                let mut np = 0;
                let deb = c3d.point_at(cf);
                min_max_pnt(
                    deb, &mut np, &mut minx, &mut miny, &mut minz, &mut maxx, &mut maxy, &mut maxz,
                );
                let fin = c3d.point_at(cl);
                min_max_pnt(
                    fin, &mut np, &mut minx, &mut miny, &mut minz, &mut maxx, &mut maxy, &mut maxz,
                );
                let mid = c3d.point_at((cf + cl) / 2.);
                min_max_pnt(
                    mid, &mut np, &mut minx, &mut miny, &mut minz, &mut maxx, &mut maxy, &mut maxz,
                );
                if !min_max_small(minx, miny, minz, maxx, maxy, maxz, toler) {
                    return false;
                }
            } else {
                //    Other case : two maximum allowed
                nb += 1;
                if nb > 2 {
                    return false;
                }
                if nb == 1 {
                    v1 = va;
                    v2 = vb;
                    *e1 = e;
                } else if nb == 2 {
                    if shape_is_same(&v1, &va) && !shape_is_same(&v2, &vb) {
                        return false;
                    }
                    if shape_is_same(&v1, &vb) && !shape_is_same(&v2, &va) {
                        return false;
                    }
                    *e2 = e;
                } else {
                    return false;
                }
            }
        }

        if nb < 2 {
            return false; // only one vertex : cannot be a strip ...
        }

        //   Checking if E1 and E2 define a Strip
        let mut dmax = 0.0;
        if !self.check_strip_edges(brep, e1, e2, tol, &mut dmax) {
            return false;
        }
        self.my_status_strip = encode_status(ShapeExtendStatus::Done3);
        true
    }

    /// OCCT CheckStripFace(F, E1, E2, tol = -1.0) (cxx L655-677).
    pub fn check_strip_face(
        &mut self,
        brep: &mut BRep,
        f: &Shape,
        e1: &mut Shape,
        e2: &mut Shape,
        tol: f64,
    ) -> bool {
        if self.check_single_strip(brep, f, e1, e2, tol) {
            return true; // it is a strip
        }

        //    IsStripSupport used as rejection. But this kind of test may be done
        //    on ANY face, once we are SURE that FindStripEdges is reliable (and fast
        //    enough)

        //  ?? record a diagnostic StripFace, but without yet lists of edges
        //  ??  Record Diagnostic "StripFace", no data (should be "Edges1" "Edges2")
        //      but direction is known (1:U  2:V)
        let mut dmax = 0.0;
        self.find_strip_edges(brep, f, e1, e2, tol, &mut dmax)
    }

    // -- Splitting vertices -------------------------------------------------

    /// OCCT CheckSplittingVertices(F, MapEdges, MapParam, theAllVert)
    /// (cxx L681-791): returns the count of found splitting vertices.  The
    /// OCCT DataMaps are keyed by the IsSame vertex identity (bridge #6);
    /// `the_all_vert` is the caller's compound (the OCCT
    /// `TopoDS_Compound&`).
    pub fn check_splitting_vertices(
        &mut self,
        brep: &mut BRep,
        f: &Shape,
        map_edges: &mut HashMap<(u64, u32), Vec<Shape>>,
        map_param: &mut HashMap<(u64, u32), Vec<f64>>,
        the_all_vert: &mut Shape,
    ) -> i32 {
        //  Prepare array of vertices with their locations //TopTools
        let mut nbv = 0;
        let mut nbp;
        let _the_builder = builder_add as fn(&mut BRep, &Shape, &Shape);
        for _itv in topexp_explorer(brep, f, ShapeType::Vertex) {
            nbv += 1;
        }

        if nbv == 0 {
            return 0;
        }
        let mut vtx: Vec<Shape> = Vec::with_capacity(nbv);
        let mut vtp: Vec<DVec3> = Vec::with_capacity(nbv);
        let mut vto: Vec<f64> = Vec::with_capacity(nbv);

        nbp = 0;
        for itv in topexp_explorer(brep, f, ShapeType::Vertex) {
            nbp += 1;
            let unv = itv;
            vtx.push(unv.clone());
            let unp = brep_tool_pnt(&unv);
            vtp.push(unp);
            let mut unt = self.my_precision;
            if unt < 0.0 {
                unt = brep_tool_tolerance(&unv);
            }
            vto.push(unt);
        }
        nbv = nbp;
        nbp = 0; // now, counting splitting vertices

        //  Check edges : are vertices (other than extremities) confused with it ?
        let sac = ShapeAnalysisCurve;
        for iv in 1..=nbv {
            let v = &vtx[(iv - 1) as usize];
            let mut list_edge: Vec<Shape> = Vec::new();
            let mut list_param: Vec<f64> = Vec::new();
            let mut issplit = false;
            for ite in topexp_explorer(brep, f, ShapeType::Edge) {
                let e = ite;
                let (v1, v2) = top_exp_vertices(&e);
                let (c3d, cf, cl) = match brep_tool_curve(&e) {
                    Some(v) => v,
                    None => continue,
                };
                if shape_is_same(v, &v1) || shape_is_same(v, &v2) {
                    continue;
                }
                let unp = vtp[(iv - 1) as usize];
                let unt = vto[(iv - 1) as usize];
                let mut proj = DVec3::ZERO;
                let mut param = 0.0;
                let dist = sac.project_cf_cl(
                    &c3d, unp, unt * 10., &mut proj, &mut param, cf, cl, true,
                );
                if dist == 0.0 {
                    continue; // smh
                }
                //  Splitting Vertex to record ?
                if dist < unt {
                    //  If Split occurs at beginning or end, it is not a split ...
                    let eps = 1.0e-06;
                    if param >= cl || param <= cf {
                        continue; // Out of range
                    }
                    let fpar = param - cf;
                    let lpar = param - cl;
                    if fpar.abs() < eps || lpar.abs() < eps {
                        continue; // Near end or start
                    }
                    list_edge.push(e);
                    list_param.push(param);
                    issplit = true;
                }
            }
            if issplit {
                nbp += 1;
                let key = (
                    std::sync::Arc::as_ptr(&v.data) as u64,
                    v.location,
                );
                builder_add(brep, the_all_vert, v);
                map_edges.insert(key, list_edge);
                map_param.insert(key, list_param);
            }
        }
        if nbp != 0 {
            self.my_status_split_vert = encode_status(ShapeExtendStatus::Done);
        }
        nbp as i32
    }

    // -- Pin ----------------------------------------------------------------

    /// OCCT CheckPin(F, whatrow, sence) (cxx L848-958).
    pub fn check_pin(&mut self, brep: &mut BRep, f: &Shape, whatrow: &mut i32, sens: &mut i32) -> bool {
        let (surf, _loc) = brep_tool_surface_loc(f);
        let Some(surf) = surf else {
            // OCCT dereferences the null surface (IsKind) — the rcad false
            // return keeps the non-surface outcome.
            return false;
        };
        if matches!(
            surf,
            Surface3::Plane(_)
                | Surface3::Cylinder(_)
                | Surface3::Cone(_)
                | Surface3::Sphere(_)
                | Surface3::Torus(_)
        ) {
            return false;
        }

        let mut toler = self.my_precision;
        if toler < 0.0 {
            toler = 1.0e-4;
        }
        let tolpin = 1.0e-9; // for sharp sharp pin

        //  Checking the poles

        //  Take the poles : they give good idea of sharpness of a pin
        let bs_poles: Option<Vec<Vec<DVec3>>> = match &surf {
            Surface3::BSpline(bs) => Some(bs.control_points.clone()),
            _ => None,
        };
        let bz_poles: Option<Vec<Vec<DVec3>>> = match &surf {
            Surface3::Bezier(bz) => Some(bz.control_points.clone()),
            _ => None,
        };
        let allpoles = if bs_poles.is_some() {
            bs_poles.unwrap()
        } else if bz_poles.is_some() {
            bz_poles.unwrap()
        } else {
            return false; // nbu == 0 || nbv == 0
        };
        let nbu = allpoles.len() as i32;
        let nbv = allpoles.first().map_or(0, |row| row.len()) as i32;
        if nbu == 0 || nbv == 0 {
            return false;
        }

        //  Check each natural bound if it is a singularity (i.e. a pin)

        *sens = 0;
        let mut stat = 0; // 0 none, 1 in U, 2 in V
        *whatrow = 0; // 0 no row, else rank of row
        stat = iso_stat(&allpoles, 1, 1, tolpin, toler);
        if stat != 0 {
            *sens = 1;
            *whatrow = nbu;
        }

        stat = iso_stat(&allpoles, 1, nbu, tolpin, toler);
        if stat != 0 {
            *sens = 1;
            *whatrow = nbu;
        }

        stat = iso_stat(&allpoles, 2, 1, tolpin, toler);
        if stat != 0 {
            *sens = 2;
            *whatrow = 1;
        }

        stat = iso_stat(&allpoles, 2, nbv, tolpin, toler);
        if stat != 0 {
            *sens = 2;
            *whatrow = nbv;
        }

        if *sens == 0 {
            return false; // no pin
        }

        match stat {
            1 => {
                self.my_status_pin = encode_status(ShapeExtendStatus::Done1);
            }
            2 => {
                self.my_status_pin = encode_status(ShapeExtendStatus::Done2);
            }
            _ => {}
        }
        //  std::cout<<(whatstat == 1 ? "Smooth" : "Sharp")<<" Pin on "<<(sens == 1 ? "U" : "V")<<" Row n0
        //  "<<whatrow<<std::endl;
        if stat == 1 {
            if check_poles(&allpoles, 2, nbv)
                || check_poles(&allpoles, 2, 1)
                || check_poles(&allpoles, 1, nbu)
                || check_poles(&allpoles, 1, 1)
            {
                self.my_status_pin = encode_status(ShapeExtendStatus::Done3);
            }
        }

        true
    }

    // -- Twisted ------------------------------------------------------------

    /// OCCT CheckTwisted(F, paramu, paramv) (cxx L972-1059).
    pub fn check_twisted(
        &mut self,
        brep: &mut BRep,
        f: &Shape,
        paramu: &mut f64,
        paramv: &mut f64,
    ) -> bool {
        let _ = brep;
        let (surf, _loc) = brep_tool_surface_loc(f);
        let Some(surf) = surf else {
            return false;
        };
        if matches!(
            surf,
            Surface3::Plane(_)
                | Surface3::Cylinder(_)
                | Surface3::Cone(_)
                | Surface3::Sphere(_)
                | Surface3::Torus(_)
        ) {
            return false;
        }

        let toler = self.my_precision;
        let _ = toler;
        ////  GeomLProp_SLProps GLS (surf,2,toler);
        // GeomAdaptor_Surface GAS(surf) (bridge #5).

        // to be done : on isos of the surface
        //  and on edges, at least of outer wire
        let nbint = 5;
        let mut nx = vec![vec![0.0f64; nbint as usize + 1]; nbint as usize + 1];
        let mut ny = vec![vec![0.0f64; nbint as usize + 1]; nbint as usize + 1];
        let mut nz = vec![vec![0.0f64; nbint as usize + 1]; nbint as usize + 1];
        let (umin, umax, vmin, vmax) = {
            let dom = surf.default_domain();
            (dom[0], dom[1], dom[2], dom[3])
        };
        let mut u = umin;
        let du = (umax - umin) / nbint as f64;
        let mut v = vmin;
        // OCCT L1001: dv = (umax - umin) / nbint — the OCCT source uses the
        // U range for V as written; kept as-is.
        let dv = (umax - umin) / nbint as f64;

        for iu in 1..=nbint {
            for iv in 1..=nbint {
                //      GLS.SetParameters (u,v);
                //      if (GLS.IsNormalDefined()) norm = GLS.Normal();
                let (_curp, v1, v2) = surf.derivatives(u, v);
                let vxnorm = v1.cross(v2);
                nx[iu as usize][iv as usize] = vxnorm.x;
                ny[iu as usize][iv as usize] = vxnorm.y;
                nz[iu as usize][iv as usize] = vxnorm.z;
                v += dv;
            }
            u += du;
            v = vmin;
        }

        //  Now, comparing normals on support surface, in both senses
        //  In principle, it suffuces to check within outer bound

        for iu in 1..nbint {
            for iv in 1..nbint {
                // We here check each normal (iu,iv) with (iu,iv+1) and with (iu+1,iv)
                // if for each test, we have negative scalar product, this means angle > 90deg
                // it is the criterion to say it is twisted
                if twisted_norm(
                    nx[iu as usize][iv as usize],
                    ny[iu as usize][iv as usize],
                    nz[iu as usize][iv as usize],
                    nx[iu as usize][(iv + 1) as usize],
                    ny[iu as usize][(iv + 1) as usize],
                    nz[iu as usize][(iv + 1) as usize],
                ) < 0.
                    || twisted_norm(
                        nx[iu as usize][iv as usize],
                        ny[iu as usize][iv as usize],
                        nz[iu as usize][iv as usize],
                        nx[(iu + 1) as usize][iv as usize],
                        ny[(iu + 1) as usize][iv as usize],
                        nz[(iu + 1) as usize][iv as usize],
                    ) < 0.
                {
                    self.my_status_twisted = encode_status(ShapeExtendStatus::Done);
                    *paramu = umin + du * iu as f64 - du / 2.0;
                    *paramv = vmin + dv * iv as f64 - dv / 2.0;
                    return true;
                }
            }
        }

        //   Now, comparing normals on edges ... to be done

        false
    }

    // -- Pin face / pin edges -----------------------------------------------

    /// OCCT CheckPinFace(F, mapEdges, toler = -1.0) (cxx L1065-1197).
    /// Warning: This function not tested on many examples.
    pub fn check_pin_face(
        &mut self,
        brep: &mut BRep,
        f: &Shape,
        map_edges: &mut HashMap<(u64, u32), Shape>,
        toler: f64,
    ) -> bool {
        // ShapeFix_Wire sfw;
        let wires = topexp_explorer(brep, f, ShapeType::Wire);
        let Some(first_wire) = wires.first() else {
            return false;
        };
        let mut coef1 = 0.0;
        let mut coef2;
        let the_cur_wire0 = first_wire.clone();
        let mut wi = ShapeAnalysisWireOrder::new();
        let mut sfw = ShapeAnalysisWire::new();
        let mut sbwd = WireData::new_from_wire(brep, &the_cur_wire0, false, true);
        sfw.load_wire_data(WireData::new_from_wire(brep, &the_cur_wire0, false, true));
        // OCCT L1079: sfw.CheckOrder(wi) — the defaults isClosed = true,
        // mode3d = true, modeBoth = false.
        sfw.check_order(brep, &mut wi, true, false, true);
        let mut newedges: Vec<Shape> = Vec::new();
        let nb = wi.nb_edges();
        for i in 1..=nb {
            newedges.push(sbwd.edge(wi.ordered(i)));
        }
        for i in 1..=nb {
            sbwd.set_edge(&newedges[(i - 1) as usize], i);
        }
        let _the_cur_wire = sbwd.wire(brep);
        let mut i = 1;
        let mut done = false;
        let mut tol = CONFUSION;
        let mut the_first_edge = Shape::null();
        let mut the_second_edge;
        let mut d1 = 0.0;
        let mut d2 = 0.0;
        for exp_e in topexp_explorer(brep, f, ShapeType::Edge) {
            let mut v1;
            let mut v2;
            let p1;
            let p2;
            if i == 1 {
                the_first_edge = exp_e;
                v1 = top_exp_first_vertex(&the_first_edge);
                v2 = top_exp_last_vertex(&the_first_edge);
                p1 = brep_tool_pnt(&v1);
                p2 = brep_tool_pnt(&v2);
                tol = brep_tool_tolerance(&v1).max(brep_tool_tolerance(&v2));
                if toler > 0.0 {
                    // tol = std::max(tol, toler); gka
                    tol = toler;
                }
                d1 = p1.distance(p2);
                if d1 == 0.0 {
                    return false;
                }
                if d1 / tol >= 1.0 {
                    coef1 = d1 / tol;
                } else {
                    continue;
                }
                if coef1 <= 3.0 {
                    continue;
                }
                i += 1;
                continue;
            }
            // Check the length of edge
            the_second_edge = exp_e;
            v1 = top_exp_first_vertex(&the_second_edge);
            v2 = top_exp_last_vertex(&the_second_edge);

            p1 = brep_tool_pnt(&v1);
            p2 = brep_tool_pnt(&v2);
            if toler == -1.0 {
                tol = brep_tool_tolerance(&v1).max(brep_tool_tolerance(&v2));
            } else {
                tol = toler;
            }
            if p1.distance(p2) > tol {
                continue;
            }
            // If there are two pin edges, record them in diagnostic
            d2 = p1.distance(p2); // gka
            if d2 == 0.0 {
                return false;
            }
            if d2 / tol >= 1.0 {
                coef2 = d2 / tol;
            } else {
                continue;
            }
            if coef2 <= 3.0 {
                continue;
            }
            if coef1 > coef2 * 10.0 {
                continue;
            }
            if coef2 > coef1 * 10.0 {
                the_first_edge = the_second_edge;
                coef1 = coef2;
                continue;
            }

            if self.check_pin_edges(&the_first_edge, &the_second_edge, coef1, coef2, toler) {
                let key = (
                    std::sync::Arc::as_ptr(&the_first_edge.data) as u64,
                    the_first_edge.location,
                );
                map_edges.insert(key, the_second_edge.clone());
                self.my_status_pin_face = encode_status(ShapeExtendStatus::Done);
                done = true;
            }

            the_first_edge = the_second_edge;
            coef1 = coef2;
            // d1 = d2;
            let _ = d1;
        }
        done
    }

    /// OCCT CheckPinEdges(theFirstEdge, theSecondEdge, coef1, coef2, toler)
    /// (cxx L1203-1341).  Warning: This function not tested on many
    /// examples.
    pub fn check_pin_edges(
        &self,
        the_first_edge: &Shape,
        the_second_edge: &Shape,
        coef1: f64,
        coef2: f64,
        toler: f64,
    ) -> bool {
        let (c1, cf1, cl1) = match brep_tool_curve(the_first_edge) {
            Some(v) => v,
            None => return false,
        };
        let (c2, cf2, cl2) = match brep_tool_curve(the_second_edge) {
            Some(v) => v,
            None => return false,
        };
        let d1 = (cf1 - cl1) / coef1;
        let d2 = (cf2 - cl2) / coef2;
        // double d1 = cf1-cl1/30; //10; gka
        // double d2 = cf2-cl2/30; //10;
        let p1 = c1.point_at(cf1);
        let p2 = c1.point_at(cl1);
        let pp1 = c2.point_at(cf2);
        let pp2 = c2.point_at(cl2);
        let tol;
        let mut paramc1 = 0.0;
        let mut paramc2 = 0.0; // =0 for deleting warning (skl)
        let the_shared_v = top_exp_last_vertex(the_first_edge);
        if toler == -1.0 {
            tol = brep_tool_tolerance(&the_shared_v);
        } else {
            tol = toler;
        }
        let pv = brep_tool_pnt(&the_shared_v);
        if pv.distance(p1) <= tol {
            paramc1 = cf1;
        } else if pv.distance(p2) <= tol {
            paramc1 = cl1;
        }
        if pv.distance(pp1) <= tol {
            paramc2 = cf2;
        } else if pv.distance(pp2) <= tol {
            paramc2 = cl2;
        }
        // Computing first derivative vectors and compare angle
        //   gp_Vec V11, V12, V21, V22;
        //   gp_Pnt tmp;
        //   C1->D2(paramc1, tmp, V11, V21);
        //   C2->D2(paramc2, tmp, V12, V22);
        //   double angle1, angle2;
        //   try{
        //     angle1 = V11.Angle(V12);
        //     angle2 = V21.Angle(V22);
        //   }
        //   catch (Standard_Failure)
        //     {
        //       std::cout << "Couldn't compute angle between derivative vectors"  <<std::endl;
        //       return false;
        //     }
        //   std::cout << "angle1 "   << angle1<< std::endl;
        //   std::cout << "angle2 "   << angle2<< std::endl;
        //   if (angle1<=0.0001) return true;
        // OCCT L1270-1297: C3 = C1 or C2 with the projection parameter pick.
        let (c3, proj): (Curve3, DVec3) = if p1.distance(p2) < pp1.distance(pp2) {
            let pr = if paramc1 == cf1 {
                c1.point_at(paramc1 + (coef1 - 3.0) * d1)
            } else {
                c1.point_at(paramc1 - 3.0 * d1)
            };
            // proj = C1->Value(paramc1 + 9*d1);
            // else proj = C1->Value(paramc1-d1);
            (c1.clone(), pr)
        } else {
            let pr = if paramc2 == cf2 {
                c2.point_at(paramc2 + (coef2 - 3.0) * d2)
            } else {
                c2.point_at(paramc2 - 3.0 * d2)
            };
            // proj = C2->Value(paramc2 + 9*d2);
            // else proj = C2->Value(paramc2 -d2);
            (c2.clone(), pr)
        };
        let mut param = 0.0;
        let f = curve_first_parameter(&c3);
        let l = curve_last_parameter(&c3);
        let mut result = DVec3::ZERO;
        let gac = GeomAdaptorCurve::new(c3.clone(), f, l);
        let sac = ShapeAnalysisCurve;
        let dist = sac.project_adaptor(&gac, proj, tol, &mut result, &mut param, true);
        // pdn check if parameter of projection is in the domain of the edge.
        if param < f || param > l {
            return false;
        }
        if dist > tol {
            return false;
        }
        if dist <= tol {
            // Computing first derivative vectors and compare angle
            // OCCT L1319-1320: C1->D2(paramc1, tmp, V11, V21);
            // C2->D2(paramc2, tmp, V12, V22).
            let v11 = c1.derivative_at(paramc1);
            let v21 = c1.derivative2_at(paramc1);
            let v12 = c2.derivative_at(paramc2);
            let v22 = c2.derivative2_at(paramc2);
            // OCCT L1322-1333: the try/catch over gp_Vec::Angle — the null
            // magnitude raise maps to the catch return false.
            let (angle1, angle2) = match (gp_vec_angle(v11, v12), gp_vec_angle(v21, v22)) {
                (Some(a1), Some(a2)) => (a1, a2),
                _ => {
                    // "Couldn't compute angle between derivative vectors"
                    return false;
                }
            };
            //       std::cout << "angle1 "   << angle1<< std::endl;
            //       std::cout << "angle2 "   << angle2<< std::endl;
            // OCCT L1336-1337: the second test reads (M_PI - angle2) twice
            // (the OCCT source as written); kept as-is.
            return (angle1 <= 0.001 && angle2 <= 0.01)
                || ((std::f64::consts::PI - angle2) <= 0.001
                    && (std::f64::consts::PI - angle2) <= 0.01);
        }

        false
    }
}

/// OCCT Geom_TrimmedCurve::FirstParameter — the rcad Trimmed carrier's
/// stored range start.
fn trimmed_first_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.first,
        _ => curve_first_parameter(c),
    }
}

/// OCCT Geom_TrimmedCurve::LastParameter.
fn trimmed_last_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.last,
        _ => curve_last_parameter(c),
    }
}

/// OCCT Geom_TrimmedCurve::BasisCurve + the parameter range — the
/// GeomAdaptor_Curve(C2T) construction data.
fn trimmed_parts(c: &Curve3) -> (Curve3, f64, f64) {
    match c {
        Curve3::Trimmed(tc) => ((*tc.curve).clone(), tc.first, tc.last),
        _ => (c.clone(), curve_first_parameter(c), curve_last_parameter(c)),
    }
}

// The OCCT CheckPinFace iterates TopExp_Explorer(F, WIRE) through
// exp_w.More() — the single-first-wire read.
#[allow(unused)]
fn _iter_marker(b: &mut BRep, f: &Shape) -> Vec<Shape> {
    iter_subshapes(b, f, false, false)
}
