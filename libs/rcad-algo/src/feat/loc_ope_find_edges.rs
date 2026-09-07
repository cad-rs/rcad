// OCCT LocOpe_FindEdges.hxx L28-56 + LocOpe_FindEdges.cxx L37-409 +
// LocOpe_FindEdges.lxx L430-474 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_FindEdges.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_FindEdges.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_FindEdges.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. BRep_Tool::Curve(edg, Loc, f, l) is re-hosted below as brep_tool_curve:
//    rcad standalone feat shapes carry their geometry on the TShape
//    (TEdgeData.curve/range) and the locations pool of the producing BRep is
//    not reachable from a bare Shape, so the world transform of the OCCT
//    TopLoc_Location branch (cxx L54-58 / L75-79) applies only when the
//    edge carries an identity location — the feat pipeline shapes do.
// 2. NCollection_List<TopoDS_Shape> is Vec<Shape>; the list iterator
//    (myItFrom/myItTo) is an Option<usize> index — None carries the OCCT
//    "iterator not initialized" state (More() == false).
// 3. Geom curve types map to geom::Curve3 variants (Geom_Line -> Line3,
//    Geom_Circle -> Circle3, Geom_Ellipse -> Ellipse3, Geom_BSplineCurve ->
//    BSplineCurve3, Geom_BezierCurve -> BezierCurve3, Geom_TrimmedCurve ->
//    TrimmedCurve3). gp_Lin/gp_Circ/gp_Elips are the rcad Line3/Circle3/
//    Ellipse3 carriers (same convention as fillet/chfi_ds.rs).
// 4. Geom_BSplineCurve knots: rcad BSplineCurve3 stores the full knot vector
//    with multiplicities expanded; NbKnots/Knots(k)/Multiplicities(k) are
//    derived by run-length splitting (bspline_knots_and_multis below).
// 5. ElCLib statics (Value/Parameter/D1/InPeriod) are re-hosted below as
//    pure-math helpers with the standard gp formulas (ElCLib.cxx); they have
//    no independent rcad module yet.
//
// first consumer: BRepFeat_Form family (3b) — BRepFeat_Form::Propagate /
// BRepFeat_MakeDPrism use LocOpe_FindEdges to match profile edges with
// support edges.

use crate::feat::brep_feat_builder::explorer;
use glam::DVec3;
use rcad_kernel::geom::{BSplineCurve3, Curve3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ ShapeType, TShape };

/// OCCT Geom curve kind handle comparison (STANDARD_TYPE checks) restricted
/// to the types LocOpe_FindEdges::Set distinguishes (cxx L65-67 / L92 / L116 /
/// L178 / L235 / L336).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LocOpeCurveKind {
    Line,
    Circle,
    Ellipse,
    BSpline,
    Bezier,
}

/// OCCT Geom_Curve::DynamicType discrimination for the kinds above; None for
/// every other curve type (cxx L65-70: "continue").
pub(crate) fn curve_kind(c: &Curve3) -> Option<LocOpeCurveKind> {
    match c {
        Curve3::Line(_) => Some(LocOpeCurveKind::Line),
        Curve3::Circle(_) => Some(LocOpeCurveKind::Circle),
        Curve3::Ellipse(_) => Some(LocOpeCurveKind::Ellipse),
        Curve3::BSpline(_) => Some(LocOpeCurveKind::BSpline),
        Curve3::Bezier(_) => Some(LocOpeCurveKind::Bezier),
        _ => None,
    }
}

/// OCCT Geom_TrimmedCurve::BasisCurve() — strip the trimmed wrapper
/// (cxx L60-64 / L81-85).
pub(crate) fn basis_curve(c: &Curve3) -> Curve3 {
    match c {
        Curve3::Trimmed(tc) => (*tc.curve).clone(),
        _ => c.clone(),
    }
}

/// OCCT BRep_Tool::Curve(edg, Loc, f, l) (BRep_Tool.cxx) — the 3D curve and
/// parameter range of the edge (architecture difference #1).
pub(crate) fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::Tolerance(edg) (BRep_Tool.cxx).
pub(crate) fn brep_tool_tolerance(edg: &Shape) -> f64 {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.tolerance,
        _ => 0.0,
    }
}

/// OCCT ElCLib::InPeriod(U, Uf, Ul) — map U into [Uf, Ul) with the period
/// p = Ul - Uf (ElCLib.cxx ElCLib::InPeriod; pure math helper, architecture
/// difference #5).
pub(crate) fn elclib_in_period(u: f64, u_f: f64, u_l: f64) -> f64 {
    let p = u_l - u_f;
    let mut r = (u - u_f) % p;
    if r < 0.0 {
        r += p;
    }
    u_f + r
}

/// OCCT ElCLib::Value(U, gp_Lin) — the point at parameter U on the line
/// (ElCLib.cxx; P = P0 + U * V).
pub(crate) fn elclib_value_lin(lin: &rcad_kernel::geom::Line3, u: f64) -> DVec3 {
    lin.origin + lin.direction * u
}

/// OCCT ElCLib::Parameter(gp_Lin, P) — the parameter of the projection of P
/// on the line (ElCLib.cxx; U = (P - P0) . V).
pub(crate) fn elclib_parameter_lin(lin: &rcad_kernel::geom::Line3, p: DVec3) -> f64 {
    (p - lin.origin).dot(lin.direction)
}

/// OCCT ElCLib::Value(U, gp_Circ) — P = Loc + R*(cos(U)*XDir + sin(U)*YDir)
/// (ElCLib.cxx).
pub(crate) fn elclib_value_circle(c: &rcad_kernel::geom::Circle3, u: f64) -> DVec3 {
    c.center + (c.x_dir * u.cos() + c.y_dir * u.sin()) * c.radius
}

/// OCCT ElCLib::Parameter(gp_Circ, P) — U = atan2((P-Loc).YDir,
/// (P-Loc).XDir) (ElCLib.cxx).
pub(crate) fn elclib_parameter_circle(c: &rcad_kernel::geom::Circle3, p: DVec3) -> f64 {
    let d = p - c.center;
    d.dot(c.y_dir).atan2(d.dot(c.x_dir))
}

/// OCCT ElCLib::D1(U, gp_Circ, P, V1) — P as Value, V1 = R*(-sin(U)*XDir +
/// cos(U)*YDir) (ElCLib.cxx).
pub(crate) fn elclib_d1_circle(c: &rcad_kernel::geom::Circle3, u: f64) -> (DVec3, DVec3) {
    let p = elclib_value_circle(c, u);
    let v = (c.y_dir * u.cos() - c.x_dir * u.sin()) * c.radius;
    (p, v)
}

/// OCCT ElCLib::Value(U, gp_Elips) — P = Loc + MajR*cos(U)*XDir +
/// MinR*sin(U)*YDir (ElCLib.cxx; YDir = Normal ^ XDir).
pub(crate) fn elclib_value_ellipse(e: &rcad_kernel::geom::Ellipse3, u: f64) -> DVec3 {
    let y_dir = e.normal.cross(e.major_dir);
    e.center + (e.major_dir * u.cos() * e.major_radius + y_dir * u.sin() * e.minor_radius)
}

/// OCCT ElCLib::Parameter(gp_Elips, P) — U = atan2((P-Loc).YDir/MinR,
/// (P-Loc).XDir/MajR) (ElCLib.cxx).
pub(crate) fn elclib_parameter_ellipse(e: &rcad_kernel::geom::Ellipse3, p: DVec3) -> f64 {
    let d = p - e.center;
    let y_dir = e.normal.cross(e.major_dir);
    (d.dot(y_dir) / e.minor_radius).atan2(d.dot(e.major_dir) / e.major_radius)
}

/// OCCT ElCLib::D1(U, gp_Elips, P, V1) — V1 = MinR*cos(U)*YDir -
/// MajR*sin(U)*XDir (ElCLib.cxx).
pub(crate) fn elclib_d1_ellipse(e: &rcad_kernel::geom::Ellipse3, u: f64) -> (DVec3, DVec3) {
    let p = elclib_value_ellipse(e, u);
    let y_dir = e.normal.cross(e.major_dir);
    let v = y_dir * (u.cos() * e.minor_radius) - e.major_dir * (u.sin() * e.major_radius);
    (p, v)
}

/// OCCT Geom_BSplineCurve::Knots()/Multiplicities() derived from the rcad
/// full knot vector (architecture difference #4) — distinct knot values with
/// their run-length multiplicities.
pub(crate) fn bspline_knots_and_multis(b: &BSplineCurve3) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for &k in &b.knots {
        if let Some(last) = knots.last().copied() {
            if k == last {
                *mults.last_mut().expect("non-empty mults") += 1;
                continue;
            }
        }
        knots.push(k);
        mults.push(1);
    }
    (knots, mults)
}

/// OCCT Geom_BezierCurve::IsRational() — weights not all 1 (architecture
/// mapping: rcad stores 1.0 weights for the polynomial case).
pub(crate) fn bezier_is_rational(weights: &[f64]) -> bool {
    weights.iter().any(|&w| w != 1.0)
}

/// OCCT LocOpe_FindEdges (LocOpe_FindEdges.hxx L28-56).
pub struct LocOpeFindEdges {
    my_f_from: Shape,   // OCCT: myFFrom (TopoDS_Shape)
    my_f_to: Shape,     // OCCT: myFTo (TopoDS_Shape)
    my_l_from: Vec<Shape>, // OCCT: myLFrom (NCollection_List<TopoDS_Shape>)
    my_l_to: Vec<Shape>,   // OCCT: myLTo (NCollection_List<TopoDS_Shape>)
    my_it_from: Option<usize>, // OCCT: myItFrom (list iterator, arch. diff. #2)
    my_it_to: Option<usize>,   // OCCT: myItTo (list iterator, arch. diff. #2)
}

impl Default for LocOpeFindEdges {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeFindEdges {
    /// OCCT LocOpe_FindEdges::LocOpe_FindEdges() (lxx L430).
    pub fn new() -> Self {
        LocOpeFindEdges {
            my_f_from: Shape::null(),
            my_f_to: Shape::null(),
            my_l_from: Vec::new(),
            my_l_to: Vec::new(),
            my_it_from: None,
            my_it_to: None,
        }
    }

    /// OCCT LocOpe_FindEdges::LocOpe_FindEdges(FFrom, FTo) (lxx L434-437) —
    /// Rust has no overloading: the value constructor carries the `_from_to`
    /// suffix.
    pub fn new_from_to(the_f_from: &Shape, the_f_to: &Shape) -> Self {
        let mut res = LocOpeFindEdges::new();
        res.set(the_f_from, the_f_to);
        res
    }

    /// OCCT LocOpe_FindEdges::Set(FFrom, FTo) (cxx L37-409).
    pub fn set(&mut self, the_f_from: &Shape, the_f_to: &Shape) {
        self.my_f_from = the_f_from.clone();
        self.my_f_to = the_f_to.clone();
        self.my_l_from.clear();
        self.my_l_to.clear();

        // OCCT cxx L44-48: explorers, curves, ranges, type handles.
        // OCCT cxx L50:
        for edgf in explorer(&self.my_f_from, ShapeType::Edge, ShapeType::Shape) {
            // OCCT cxx L52-53: TopoDS::Edge(expf.Current()); BRep_Tool::Curve.
            let Some((mut cf, ff, lf)) = brep_tool_curve(&edgf) else {
                continue;
            };
            // OCCT cxx L54-58: location transform (architecture difference #1 —
            // standalone feat shapes carry identity locations).
            // OCCT cxx L59-64: Tf = Cf->DynamicType(); trimmed -> basis.
            if matches!(cf, Curve3::Trimmed(_)) {
                cf = basis_curve(&cf);
            }
            // OCCT cxx L65-70: only Line/Circle/Ellipse/BSpline/Bezier survive.
            let Some(tfk) = curve_kind(&cf) else {
                continue;
            };
            // OCCT cxx L71:
            for edgt in explorer(&self.my_f_to, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L73-74.
                let Some((mut ct, ft, lt)) = brep_tool_curve(&edgt) else {
                    continue;
                };
                // OCCT cxx L75-79: location transform (arch. diff. #1).
                // OCCT cxx L80-85: Tt = Ct->DynamicType(); trimmed -> basis.
                if matches!(ct, Curve3::Trimmed(_)) {
                    ct = basis_curve(&ct);
                }
                let ttk = curve_kind(&ct);
                // OCCT cxx L86-89: if (Tt != Tf) continue.
                if ttk != Some(tfk) {
                    continue;
                }
                // OCCT cxx L90-91: "On a presomption de confusion".
                let mut tol = rcad_kernel::precision::CONFUSION;
                match tfk {
                    LocOpeCurveKind::Line => {
                        // OCCT cxx L92-93.
                        let rcad_kernel::geom::Curve3::Line(lif) = &cf else {
                            continue;
                        };
                        let rcad_kernel::geom::Curve3::Line(lit) = &ct else {
                            continue;
                        };
                        // OCCT cxx L96-99.
                        let p1 = elclib_value_lin(lif, ff);
                        let p2 = elclib_value_lin(lif, lf);
                        let prm1 = elclib_parameter_lin(lit, p1);
                        let prm2 = elclib_parameter_lin(lit, p2);
                        if prm1 >= ft - tol && prm1 <= lt + tol && prm2 >= ft - tol && prm2 <= lt + tol {
                            // OCCT cxx L102-107.
                            tol *= tol;
                            let mut pt = elclib_value_lin(lit, prm1);
                            if pt.distance_squared(p1) <= tol {
                                pt = elclib_value_lin(lit, prm2);
                                if pt.distance_squared(p2) <= tol {
                                    self.my_l_from.push(edgf.clone());
                                    self.my_l_to.push(edgt.clone());
                                    break;
                                }
                            }
                        }
                    }
                    LocOpeCurveKind::Circle => {
                        // OCCT cxx L116-119.
                        let rcad_kernel::geom::Curve3::Circle(cif) = &cf else {
                            continue;
                        };
                        let rcad_kernel::geom::Curve3::Circle(cit) = &ct else {
                            continue;
                        };
                        if (cif.radius - cit.radius).abs() <= tol
                            && cif.center.distance_squared(cit.center) <= tol * tol
                        {
                            // OCCT cxx L123-128: start point, period calibration,
                            // same-direction detection.
                            let (p1, tgf) = elclib_d1_circle(cif, ff);
                            let p2 = elclib_value_circle(cif, lf);

                            let mut prm1 = elclib_parameter_circle(cit, p1);
                            let tol2d = rcad_kernel::precision::PCONFUSION;
                            if (prm1 - ft).abs() <= tol2d {
                                prm1 = ft;
                            }
                            prm1 = elclib_in_period(prm1, ft, ft + 2.0 * std::f64::consts::PI);
                            let (_, tgt) = elclib_d1_circle(cit, prm1);

                            let mut prm2 = elclib_parameter_circle(cit, p2);
                            if tgt.dot(tgf) > 0.0 {
                                // OCCT cxx L140-146: "meme sens".
                                while prm2 <= prm1 {
                                    prm2 += 2.0 * std::f64::consts::PI;
                                }
                            } else {
                                // OCCT cxx L147-157.
                                if (prm1 - ft).abs() <= rcad_kernel::precision::ANGULAR {
                                    prm1 += 2.0 * std::f64::consts::PI;
                                }
                                while prm2 >= prm1 {
                                    prm2 -= 2.0 * std::f64::consts::PI;
                                }
                            }

                            // OCCT cxx L159-175.
                            if prm1 >= ft - tol && prm1 <= lt + tol && prm2 >= ft - tol && prm2 <= lt + tol {
                                self.my_l_from.push(edgf.clone());
                                self.my_l_to.push(edgt.clone());
                                break;
                            }
                            // OCCT cxx L165-175: "Cas non traite : on est a
                            // cheval" — nothing appended.
                        }
                    }
                    LocOpeCurveKind::Ellipse => {
                        // OCCT cxx L178-181.
                        let rcad_kernel::geom::Curve3::Ellipse(cif) = &cf else {
                            continue;
                        };
                        let rcad_kernel::geom::Curve3::Ellipse(cit) = &ct else {
                            continue;
                        };
                        if (cif.major_radius - cit.major_radius).abs() <= tol
                            && (cif.minor_radius - cit.minor_radius).abs() <= tol
                            && cif.center.distance_squared(cit.center) <= tol * tol
                        {
                            // OCCT cxx L187-196.
                            let (p1, tgf) = elclib_d1_ellipse(cif, ff);
                            let p2 = elclib_value_ellipse(cif, lf);

                            let mut prm1 = elclib_parameter_ellipse(cit, p1);
                            prm1 = elclib_in_period(prm1, ft, ft + 2.0 * std::f64::consts::PI);
                            let (_, tgt) = elclib_d1_ellipse(cit, prm1);

                            let mut prm2 = elclib_parameter_ellipse(cit, p2);
                            if tgt.dot(tgf) > 0.0 {
                                // OCCT cxx L199-205: "meme sens".
                                while prm2 <= prm1 {
                                    prm2 += 2.0 * std::f64::consts::PI;
                                }
                            } else {
                                // OCCT cxx L206-216.
                                if (prm1 - ft).abs() <= rcad_kernel::precision::ANGULAR {
                                    prm1 += 2.0 * std::f64::consts::PI;
                                }
                                while prm2 >= prm1 {
                                    prm2 -= 2.0 * std::f64::consts::PI;
                                }
                            }

                            // OCCT cxx L218-233.
                            if prm1 >= ft - tol && prm1 <= lt + tol && prm2 >= ft - tol && prm2 <= lt + tol {
                                self.my_l_from.push(edgf.clone());
                                self.my_l_to.push(edgt.clone());
                                break;
                            }
                            // OCCT cxx L224-232: "Cas non traite" — nothing.
                        }
                    }
                    LocOpeCurveKind::BSpline => {
                        // OCCT cxx L235-238.
                        let rcad_kernel::geom::Curve3::BSpline(bf) = &cf else {
                            continue;
                        };
                        let rcad_kernel::geom::Curve3::BSpline(bt) = &ct else {
                            continue;
                        };
                        // OCCT cxx L240-246: NbPoles.
                        let mut is_same = true;
                        let nbpoles = bf.control_points.len();
                        if nbpoles != bt.control_points.len() {
                            is_same = false;
                        }

                        if is_same {
                            // OCCT cxx L248-254: NbKnots.
                            let nbknots = bspline_knots_and_multis(bf).0.len();
                            if nbknots != bspline_knots_and_multis(bt).0.len() {
                                is_same = false;
                            }

                            if is_same {
                                // OCCT cxx L256-269: poles within BRep_Tool::Tolerance.
                                let tol3d = brep_tool_tolerance(&edgt);
                                for p in 1..=nbpoles {
                                    if bf.control_points[p - 1].distance(bt.control_points[p - 1]) > tol3d {
                                        is_same = false;
                                        break;
                                    }
                                }

                                if is_same {
                                    // OCCT cxx L271-291: knots and multiplicities.
                                    let (kf, mf) = bspline_knots_and_multis(bf);
                                    let (kt, mt) = bspline_knots_and_multis(bt);
                                    for k in 1..=nbknots {
                                        if kf[k - 1] - kt[k - 1] > tol {
                                            is_same = false;
                                            break;
                                        }
                                        if (mf[k - 1] - mt[k - 1]) as f64 > tol {
                                            is_same = false;
                                            break;
                                        }
                                    }

                                    // OCCT cxx L293-306: rationality agreement.
                                    let rf = bf.weights.iter().any(|&w| w != 1.0);
                                    let rt = bt.weights.iter().any(|&w| w != 1.0);
                                    if !rf && rt {
                                        is_same = false;
                                    } else if rf && !rt {
                                        is_same = false;
                                    }

                                    // OCCT cxx L308-321: weights comparison.
                                    if is_same && rf {
                                        for w in 1..=nbpoles {
                                            if (bf.weights[w - 1] - bt.weights[w - 1]).abs() > tol {
                                                is_same = false;
                                                break;
                                            }
                                        }
                                    }

                                    // OCCT cxx L323-331.
                                    if is_same {
                                        self.my_l_from.push(edgf.clone());
                                        self.my_l_to.push(edgt.clone());
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    LocOpeCurveKind::Bezier => {
                        // OCCT cxx L336-339.
                        let rcad_kernel::geom::Curve3::Bezier(bf) = &cf else {
                            continue;
                        };
                        let rcad_kernel::geom::Curve3::Bezier(bt) = &ct else {
                            continue;
                        };
                        // OCCT cxx L341-347: NbPoles.
                        let mut is_same = true;
                        let nbpoles = bf.control_points.len();
                        if nbpoles != bt.control_points.len() {
                            is_same = false;
                        }

                        if is_same {
                            // OCCT cxx L349-361: poles within Tol.
                            for p in 1..=nbpoles {
                                if bf.control_points[p - 1].distance(bt.control_points[p - 1]) > tol {
                                    is_same = false;
                                    break;
                                }
                            }

                            if is_same {
                                // OCCT cxx L363-378: rationality agreement.
                                let rf = bezier_is_rational(&bf.weights);
                                let rt = bezier_is_rational(&bt.weights);
                                if !rf && rt {
                                    is_same = false;
                                } else if rf && !rt {
                                    is_same = false;
                                }

                                // OCCT cxx L380-393: weights comparison.
                                if is_same && rf {
                                    for w in 1..=nbpoles {
                                        if (bf.weights[w - 1] - bt.weights[w - 1]).abs() > tol {
                                            is_same = false;
                                            break;
                                        }
                                    }
                                }

                                // OCCT cxx L395-403.
                                if is_same {
                                    self.my_l_from.push(edgf.clone());
                                    self.my_l_to.push(edgt.clone());
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// OCCT LocOpe_FindEdges::InitIterator() (lxx L441-445).
    pub fn init_iterator(&mut self) {
        self.my_it_from = if self.my_l_from.is_empty() { None } else { Some(0) };
        self.my_it_to = if self.my_l_to.is_empty() { None } else { Some(0) };
    }

    /// OCCT LocOpe_FindEdges::More() (lxx L449-452).
    pub fn more(&self) -> bool {
        self.my_it_from.map_or(false, |i| i < self.my_l_from.len())
    }

    /// OCCT LocOpe_FindEdges::EdgeFrom() (lxx L456-459) — TopoDS::Edge cast;
    /// the list only ever receives edges (cxx L109/L161/...).
    pub fn edge_from(&self) -> Shape {
        self.my_l_from[self.my_it_from.expect("More() == true")].clone()
    }

    /// OCCT LocOpe_FindEdges::EdgeTo() (lxx L463-466).
    pub fn edge_to(&self) -> Shape {
        self.my_l_to[self.my_it_to.expect("More() == true")].clone()
    }

    /// OCCT LocOpe_FindEdges::Next() (lxx L470-474).
    pub fn next(&mut self) {
        self.my_it_from = self.my_it_from.map(|i| i + 1);
        self.my_it_to = self.my_it_to.map(|i| i + 1);
    }
}
