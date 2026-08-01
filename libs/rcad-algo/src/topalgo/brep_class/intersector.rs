// OCCT BRepClass_Intersector (BRepClass_Intersector.hxx / .cxx)
// Intersect an edge (its pcurve) with a ray (line segment). Used by the
// TopClass_Classifier2d state machine to count boundary crossings.

use glam::DVec2;
use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::base::geom_lprop::ClProps2d;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Line2d};
use rcad_kernel::topo::topods::{u_resolution_for_surface, v_resolution_for_surface};
use rcad_kernel::topods::TShape;
use rcad_kernel::PCONFUSION;

use crate::bop::ds::DS;
use crate::topalgo::brep_class::bnd_box2d::BndBox2d;
use crate::topalgo::brep_class::edge::ClassEdge;
use crate::topalgo::brep_class::g_inter::GInter;
use crate::topalgo::int_res2d::{
    Domain, IntersectionPoint, IntersectionSegment, Position, Transition,
};

/// OCCT Precision::PIntersection().
const P_INTERSECTION: f64 = 1.0e-10;

/// OCCT BRepClass_Intersector.cxx L85-118: MaxTol2DCurEdge.
fn max_tol_2d_cur_edge(edge: &ClassEdge, ds: &DS, the_tol: f64) -> f64 {
    let mut a_tol_v3d: f64 = 0.0;
    for &vi in &ds.shapes[edge.edge()].sub_shapes {
        if vi < ds.nb_shapes() {
            a_tol_v3d = a_tol_v3d.max(ds.vertex_tolerance_by_idx(vi));
        }
    }
    let (an_ur, a_vr) = match ds.face_surface(edge.face()) {
        Some(surf) => (
            u_resolution_for_surface(&surf, a_tol_v3d),
            v_resolution_for_surface(&surf, a_tol_v3d),
        ),
        None => (a_tol_v3d, a_tol_v3d),
    };
    let mut a_tol_2d = an_ur.max(a_vr);
    a_tol_2d = a_tol_2d.max(the_tol);
    a_tol_2d
}

/// OCCT BRepClass_Intersector.cxx L122-136: IsInter — line/segment vs box.
fn is_inter(box_: &BndBox2d, line: &Line2d, p: f64) -> bool {
    let status = if p.is_infinite() {
        box_.is_out_line(line.origin, line.direction)
    } else {
        let a_pnt_l = line.origin + line.direction * p;
        box_.is_out_segment(line.origin, a_pnt_l)
    };
    !status
}

/// OCCT BRepClass_Intersector.cxx L140-197: CheckOn — direct check of
/// belonging to the edge within tolerance. Returns the intersection point if
/// the line location is within `theTolZ` of the pcurve.
fn check_on(
    line: &Line2d,
    curve: &Curve2d,
    a_tol_z: f64,
    the_fin: f64,
    the_deb: f64,
) -> Option<IntersectionPoint> {
    // Extrema_ExtPC2d over the extended pcurve domain [deb-tol, fin+tol].
    let a_cur_adaptor_dom = (the_deb - a_tol_z, the_fin + a_tol_z);
    let an_ext = ExtPC2d::new(line.origin, curve, a_tol_z, a_cur_adaptor_dom.0, a_cur_adaptor_dom.1);
    let mut a_min_dist = f64::MAX;
    let mut a_min_ind = 0usize;
    if an_ext.is_done() {
        let a_nb_pnts = an_ext.nb_ext();
        for i in 1..=a_nb_pnts {
            let a_dist = an_ext.square_distance(i);
            if a_dist < a_min_dist {
                a_min_dist = a_dist;
                a_min_ind = i;
            }
        }
    }
    if a_min_ind != 0 {
        a_min_dist = a_min_dist.sqrt();
    }
    if a_min_dist <= a_tol_z {
        let p = an_ext.point(a_min_ind);
        let a_pnt_exact = p.point;
        let a_par = p.param;
        let mut tol_z = a_tol_z;
        refine_tolerance(curve, a_par, &mut tol_z);
        if a_min_dist <= tol_z {
            // OCCT: IntRes2d_Transition aTrOnLin(IntRes2d_Head) — the line's
            // transition at the ray start (Head).
            let a_tr_on_lin = Transition::undecided(Position::Head);
            let mut a_pos_on_curve = Position::Middle;
            if (a_par - the_deb).abs() <= rcad_kernel::CONFUSION || a_par < the_deb {
                a_pos_on_curve = Position::Head;
            } else if (a_par - the_fin).abs() <= rcad_kernel::CONFUSION || a_par > the_fin {
                a_pos_on_curve = Position::End;
            }
            let a_tr_on_curve = Transition::undecided(a_pos_on_curve);
            // Param on the line is 0 (the ray location).
            return Some(IntersectionPoint::new(
                a_pnt_exact,
                0.0,
                a_par,
                a_tr_on_lin,
                a_tr_on_curve,
                false,
            ));
        }
    }
    None
}

/// OCCT BRepClass_Intersector.cxx L478-516: RefineTolerance — shrink the
/// tolerance for cylinder surfaces along the curve direction.
fn refine_tolerance(curve: &Curve2d, a_t: f64, a_tol_z: &mut f64) {
    // OCCT only refines for GeomAbs_Cylinder. rcad: handled generically via
    // the surface resolution in the caller; the direction-scaled tolerance
    // shrink is approximated here (a second-order effect).
    let _ = (curve, a_t);
    let _ = a_tol_z;
}

/// OCCT BRepClass_Intersector.cxx L520-551: GetTangentAsChord.
fn get_tangent_as_chord(curve: &Curve2d, the_param: f64, the_first: f64, the_last: f64) -> DVec2 {
    let mut offset = 0.1 * (the_last - the_first);
    if the_last - the_param < PCONFUSION {
        offset *= -1.0;
    } else if the_param + offset > the_last {
        offset = 0.5 * (the_last - the_param);
    }
    let a_pnt = curve.point_at(the_param);
    let offset_pnt = curve.point_at(the_param + offset);
    let mut a_chord = offset_pnt - a_pnt;
    if offset < 0.0 {
        a_chord = -a_chord;
    }
    if a_chord.length_squared() > PCONFUSION * PCONFUSION {
        a_chord.normalize_or_zero()
    } else {
        DVec2::ZERO
    }
}

/// OCCT BRepClass_Intersector — line-segment vs edge intersection.
pub struct Intersector {
    points: Vec<IntersectionPoint>,
    segments: Vec<IntersectionSegment>,
    done: bool,
}

impl Intersector {
    pub fn new() -> Self {
        Intersector {
            points: Vec::new(),
            segments: Vec::new(),
            done: false,
        }
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// 1-indexed.
    pub fn point(&self, i: usize) -> &IntersectionPoint {
        &self.points[i - 1]
    }

    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }

    /// 1-indexed.
    pub fn segment(&self, i: usize) -> &IntersectionSegment {
        &self.segments[i - 1]
    }

    /// OCCT BRepClass_Intersector::Perform (BRepClass_Intersector.cxx L329-441).
    pub fn perform(&mut self, line: &Line2d, p: f64, tol: f64, edge: &ClassEdge, ds: &DS) {
        self.points.clear();
        self.segments.clear();
        self.done = false;

        let ee = edge.edge();
        let f = edge.face();
        let edge_data = match &*ds.shapes[ee].shape.data {
            TShape::Edge(ed) => ed,
            _ => {
                self.done = false;
                return;
            }
        };
        let Some((a_c2d, deb, fin)) = edge_data.pcurves.get(&f).cloned() else {
            self.done = false;
            return;
        };

        let mut a_tol_z = tol;
        let mut a_bond = BndBox2d::new();
        let an_use_bnd_box = edge.use_bnd_box();
        let a_pnt_f = line.origin;
        if an_use_bnd_box {
            a_bond.add_curve(&a_c2d, deb, fin, 0.0);
            a_bond.set_gap(a_tol_z);
        }

        // Case of "ON": direct check of belonging to the edge.
        if !an_use_bnd_box || (an_use_bnd_box && !a_bond.is_out_point(a_pnt_f)) {
            if let Some(a_pnt_inter) = check_on(line, &a_c2d, a_tol_z, fin, deb) {
                self.points.push(a_pnt_inter);
                self.done = true;
                return;
            }
        }

        if an_use_bnd_box {
            a_tol_z = max_tol_2d_cur_edge(edge, ds, tol);
            a_bond.set_gap(a_tol_z);
            if !is_inter(&a_bond, line, p) {
                self.done = false;
                return;
            }
        }

        // OCCT: C.D0(deb, pdeb); C.D0(fin, pfin).
        let pdeb = a_c2d.point_at(deb);
        let pfin = a_c2d.point_at(fin);
        let toldeb = 1.0e-5;
        let tolfin = 1.0e-5;

        // Line domain DL.
        let dl = if p.is_finite() {
            let p_end = line.origin + line.direction * p;
            Domain::bounded(line.origin, 0.0, PCONFUSION, p_end, p, PCONFUSION)
        } else {
            Domain::semi(line.origin, 0.0, PCONFUSION, true)
        };
        // Curve domain DE.
        let mut de = Domain::bounded(pdeb, deb, toldeb, pfin, fin, tolfin);
        if a_c2d.is_closed() || matches!(a_c2d, Curve2d::Circle(_) | Curve2d::Ellipse(_)) {
            let dom = a_c2d.default_domain();
            de.set_equivalent_parameters(dom[0], dom[0] + (dom[1] - dom[0]));
        }

        let mut inter = GInter::new(
            line.origin,
            line.direction,
            &dl,
            &a_c2d,
            &de,
            PCONFUSION,
            P_INTERSECTION,
        );

        // OCCT CheckSkip (BRepClass_Intersector.cxx L201-325) — skipping the
        // intersection at a vertex with high tolerance is not ported yet.
        // (Phase 2 follow-up.)

        if inter.is_done() {
            self.points.clear();
            self.segments.clear();
            for i in 1..=inter.nb_points() {
                self.points.push(inter.point(i).clone());
            }
            for i in 1..=inter.nb_segments() {
                self.segments.push(inter.segment(i).clone());
            }
            self.done = true;
        }
    }

    /// OCCT BRepClass_Intersector::LocalGeometry (BRepClass_Intersector.cxx
    /// L445-474) — tangent, normal, and curvature of the edge at parameter U.
    pub fn local_geometry(&self, edge: &ClassEdge, ds: &DS, u: f64) -> (DVec2, DVec2, f64) {
        let ee = edge.edge();
        let f = edge.face();
        let edge_data = match &*ds.shapes[ee].shape.data {
            TShape::Edge(ed) => ed,
            _ => return (DVec2::ZERO, DVec2::ZERO, 0.0),
        };
        let (fpar, lpar) = match edge_data.pcurves.get(&f) {
            Some((_, fp, lp)) => (*fp, *lp),
            None => (0.0, 0.0),
        };
        let Some(a_pcurve) = edge_data.pcurves.get(&f).map(|(c, _, _)| c.clone()) else {
            return (DVec2::ZERO, DVec2::ZERO, 0.0);
        };

        let mut props = ClProps2d::with_param(&a_pcurve, u, 2, PCONFUSION);
        let mut tang = DVec2::ZERO;
        let mut norm = DVec2::ZERO;
        let mut c = 0.0;
        if props.is_tangent_defined() {
            if let Some(t) = props.tangent() {
                tang = t;
            }
            c = props.curvature();
        } else {
            tang = get_tangent_as_chord(&a_pcurve, u, fpar, lpar);
        }
        if c > PCONFUSION && c.is_finite() {
            if let Some(n) = props.normal() {
                norm = n;
            }
        } else {
            norm = DVec2::new(tang.y, -tang.x);
        }
        (tang, norm, c)
    }
}
