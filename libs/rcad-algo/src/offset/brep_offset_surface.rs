// OCCT BRepOffset.cxx L1-365 + BRepOffset.hxx L31-65 + BRepOffset_Status.hxx
// L25-31 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset.cxx / .hxx / BRepOffset_Status.hxx
//
// OCCT inheritance chain (hxx L31): none — BRepOffset is a static-only
// auxiliary class.  The task name "BRepOffset_Surface" is the static
// BRepOffset::Surface member (cxx L44-202); the file also carries
// CollapseSingularities (cxx L206-365).
//
// Architecture differences:
// 1. Handle(Geom_Surface) polymorphism -> the rcad Surface3 enum; the OCCT
//    `TheType == STANDARD_TYPE(Geom_X)` DynamicType equality maps to an
//    exact enum-variant match (NOT IsKind — the OCCT chains check the exact
//    dynamic type, so a RectangularTrimmedSurface does not hit the plane
//    branch).
// 2. The rcad analytic surfaces carry the gp_Ax3 frame flattened
//    (location/axis/x_dir/y_dir); gp_Ax3::Rotate(Ax1, PI) on the X/Y
//    directions maps to the exact PI-rotation closed form v' = 2(v.n)n - v
//    (the OCCT vehicle is the gp_Trsf quaternion matrix — same exact result
//    for the angle-PI fixed form up to floating-point order).
// 3. gp_Ax3::Direct() — the rcad cylinder frame derives the sense as
//    (X ^ Y) . Z > 0 (gp_Ax3.cxx); the plane/sphere/torus frames are
//    right-handed by rcad construction.
// 4. `allowC0` — the OCCT flag feeds Geom_OffsetSurface's C0 handling; the
//    rcad OffsetSurface carries no C0 gate (annotated at the construction).
// 5. BRepOffset_Status lives in brep_offset_offset.rs (the batch-1 module);
//    re-exported here under its OCCT name mapping.

use glam::DVec3;
use rcad_kernel::geom::{
    ConicalSurface, CylindricalSurface, OffsetSurface, SphericalSurface, Surface3,
    ToroidalSurface, TrimmedSurface,
};
use rcad_kernel::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool::{brep_tool_pnt, brep_tool_tolerance, explorer, top_exp_vertices_raw};
// The BRep_Tool::Degenerated carrier home (the feat package precedent).
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_degenerated;

pub use super::brep_offset_offset::BRepOffsetStatus;

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Surface (BRepOffset.cxx L44-202).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset::Surface(Surface, Offset, theStatus, allowC0 = false)
/// (BRepOffset.cxx L44-202) — returns the Offset surface computed from
/// `the_surface` at `the_offset`.  If possible, the real surface type is
/// returned (e.g. an Offset of a plane is a plane); otherwise a
/// Geom_OffsetSurface.
pub fn brep_offset_surface(
    the_surface: &Surface3,
    the_offset: f64,
    the_status: &mut BRepOffsetStatus,
    allow_c0: bool,
) -> Surface3 {
    // OCCT L45: constexpr double Tol = Precision::Confusion();
    let tol = rcad_kernel::precision::CONFUSION;

    // OCCT L46-48: theStatus = BRepOffset_Good; Result; TheType.
    *the_status = BRepOffsetStatus::Good;
    let mut result: Option<Surface3> = None;

    // OCCT L50: occ::handle<Standard_Type> TheType = Surface->DynamicType();
    // OCCT L52-62: Geom_Plane branch.
    if let Surface3::Plane(p) = the_surface {
        // OCCT L57: T = XDirection() ^ YDirection().
        let mut t = p.u_dir.cross(p.v_dir);
        // OCCT L58: T *= Offset.
        t *= the_offset;
        // OCCT L59: Result = down_cast<Geom_Plane>(P->Translated(T)).
        result = Some(Surface3::Plane(rcad_kernel::geom::Plane {
            origin: p.origin + t,
            ..*p
        }));
    }
    // OCCT L63-90: Geom_CylindricalSurface branch.
    else if let Surface3::Cylinder(c) = the_surface {
        // OCCT L65-67: Radius = C->Radius(); Axis = C->Position().
        let mut radius = c.radius;
        let mut axis = cylinder_ax3(c);
        // OCCT L68-75: the Direct() sense selects the offset sign.
        if cylinder_ax3_direct(c) {
            radius += the_offset;
        } else {
            radius -= the_offset;
        }
        if radius >= tol {
            // OCCT L78: Result = new Geom_CylindricalSurface(Axis, Radius).
            result = Some(Surface3::Cylinder(CylindricalSurface {
                origin: axis.0,
                axis: axis.1,
                radius,
                ref_dir: axis.2,
                y_dir: c.y_dir,
            }));
        } else if radius <= -tol {
            // OCCT L82-85: rotate the frame by PI about the axis; the radius
            // is negated and the status reversed.
            let ax1_dir = axis.1;
            axis.2 = rotate_pi_about(axis.2, ax1_dir);
            if let Some(ref mut y) = axis.3 {
                *y = rotate_pi_about(*y, ax1_dir);
            }
            result = Some(Surface3::Cylinder(CylindricalSurface {
                origin: axis.0,
                axis: axis.1,
                radius: radius.abs(),
                ref_dir: axis.2,
                y_dir: axis.3,
            }));
            *the_status = BRepOffsetStatus::Reversed;
        } else {
            // OCCT L88: theStatus = BRepOffset_Degenerated.
            *the_status = BRepOffsetStatus::Degenerated;
        }
    }
    // OCCT L91-113: Geom_ConicalSurface branch.
    else if let Surface3::Cone(c) = the_surface {
        // OCCT L93-96: Alpha = SemiAngle(); Radius = RefRadius()
        //               + Offset * cos(Alpha); Axis = Position().
        let mut alpha = c.half_angle_rad;
        let mut radius = c.radius + the_offset * alpha.cos();
        let mut axis = cone_ax3(c);
        if radius >= 0. {
            // OCCT L97-102: translate along Z by -Offset*sin(Alpha).
            let mut z = axis.1;
            z *= -the_offset * alpha.sin();
            axis.0 += z;
        } else {
            // OCCT L103-111: negate the radius, translate, rotate by PI,
            // negate the semi-angle.
            radius = -radius;
            let mut z = axis.1;
            z *= -the_offset * alpha.sin();
            axis.0 += z;
            let ax1_dir = axis.1;
            axis.2 = rotate_pi_about(axis.2, ax1_dir);
            alpha = -alpha;
        }
        // OCCT L112: Result = new Geom_ConicalSurface(Axis, Alpha, Radius).
        result = Some(Surface3::Cone(ConicalSurface {
            apex: axis.0,
            axis: axis.1,
            radius,
            half_angle_rad: alpha,
            ref_dir: axis.2,
        }));
    }
    // OCCT L114-142: Geom_SphericalSurface branch.
    else if let Surface3::Sphere(s) = the_surface {
        // OCCT L116-118: Radius = S->Radius(); Axis = S->Position().
        let mut radius = s.radius;
        let mut axis = sphere_ax3(s);
        if sphere_ax3_direct(s) {
            radius += the_offset;
        } else {
            radius -= the_offset;
        }
        if radius >= tol {
            // OCCT L129: Result = new Geom_SphericalSurface(Axis, Radius).
            result = Some(Surface3::Sphere(SphericalSurface {
                center: axis.0,
                axis: axis.1,
                radius,
                ref_dir: axis.2,
            }));
        } else if radius <= -tol {
            // OCCT L133-137: rotate by PI, ZReverse the axis, the radius is
            // passed NEGATED (positive) and the status reversed.
            let ax1_dir = axis.1;
            axis.2 = rotate_pi_about(axis.2, ax1_dir);
            // OCCT L134: Axis.ZReverse() == gp_Ax3.hxx L131: axis.Reverse().
            axis.1 = -axis.1;
            result = Some(Surface3::Sphere(SphericalSurface {
                center: axis.0,
                axis: axis.1,
                radius: -radius,
                ref_dir: axis.2,
            }));
            *the_status = BRepOffsetStatus::Reversed;
        } else {
            // OCCT L140: theStatus = BRepOffset_Degenerated.
            *the_status = BRepOffsetStatus::Degenerated;
        }
    }
    // OCCT L143-172: Geom_ToroidalSurface branch.
    else if let Surface3::Torus(s) = the_surface {
        // OCCT L145-148: MajorRadius; MinorRadius; Axis = Position().
        let major_radius = s.major_radius;
        let mut minor_radius = s.minor_radius;
        let axis = torus_ax3(s);
        if minor_radius < major_radius {
            // OCCT L149-158: the Direct() sense selects the offset sign.
            if torus_ax3_direct(s) {
                minor_radius += the_offset;
            } else {
                minor_radius -= the_offset;
            }
            if minor_radius >= tol {
                // OCCT L161: Result = new Geom_ToroidalSurface(Axis,
                //             MajorRadius, MinorRadius).
                result = Some(Surface3::Torus(ToroidalSurface {
                    center: axis.0,
                    axis: axis.1,
                    ref_dir: axis.2,
                    major_radius,
                    minor_radius,
                }));
            } else if minor_radius <= -tol {
                // OCCT L163-166: theStatus = BRepOffset_Reversed (the Result
                // stays null and falls through to the OffsetSurface form).
                *the_status = BRepOffsetStatus::Reversed;
            } else {
                // OCCT L169: theStatus = BRepOffset_Degenerated.
                *the_status = BRepOffsetStatus::Degenerated;
            }
        }
    }
    // OCCT L173-178: Geom_SurfaceOfRevolution / Geom_SurfaceOfLinearExtrusion
    // branches — empty in OCCT (fall through to the OffsetSurface form).
    // OCCT L179-181: Geom_BSplineSurface branch — empty in OCCT.
    // OCCT L182-191: Geom_RectangularTrimmedSurface branch.
    else if let Surface3::Trimmed(s) = the_surface {
        // OCCT L186-187: S->Bounds(U1, U2, V1, V2).
        let (u1, u2, v1, v2) = (s.trim[0], s.trim[1], s.trim[2], s.trim[3]);
        // OCCT L188-189: recursive Surface(BasisSurface, Offset, ...).
        let off = brep_offset_surface(
            &s.basis,
            the_offset,
            the_status,
            allow_c0,
        );
        // OCCT L190: Result = new Geom_RectangularTrimmedSurface(Off, U1,
        //             U2, V1, V2).
        result = Some(Surface3::Trimmed(TrimmedSurface::new(off, u1, u2, v1, v2)));
    }
    // OCCT L192-194: Geom_OffsetSurface branch — empty in OCCT.

    // OCCT L196-199: if (Result.IsNull()) Result = new
    //                Geom_OffsetSurface(Surface, Offset, allowC0).
    let result = match result {
        Some(r) => r,
        None => {
            // Architecture difference #4: the OCCT allowC0 flag gates the
            // Geom_OffsetSurface C0 handling; the rcad OffsetSurface carries
            // no C0 gate.
            let _ = allow_c0;
            Surface3::Offset(OffsetSurface {
                basis: Box::new(the_surface.clone()),
                offset_distance: the_offset,
            })
        }
    };

    // OCCT L201: return Result.
    result
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset::CollapseSingularities (BRepOffset.cxx L206-365).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset::CollapseSingularities(theSurface, theFace, thePrecision)
/// (BRepOffset.cxx L206-365) — preprocess a bspline (or bezier/revolution —
/// future work, as in OCCT) surface to be offset by collapsing each singular
/// side to a single point.  Returns either the original surface or its
/// modified copy (if some poles have been moved).
pub fn collapse_singularities(
    the_surface: &Surface3,
    the_face: &Shape,
    the_precision: f64,
) -> Surface3 {
    // OCCT L212-218: check the surface type — for the moment only bspline
    // surfaces are treated.
    let a_bspline = match the_surface {
        Surface3::BSpline(b) => b,
        _ => return the_surface.clone(),
    };

    // OCCT L221-222: find singularities (vertices of degenerated edges).
    let mut a_degen_pnt: Vec<DVec3> = Vec::new();
    let mut a_degen_tol: Vec<f64> = Vec::new();
    // OCCT L223-239: for (TopExp_Explorer anExp(theFace, TopAbs_EDGE); ...).
    for an_edge in explorer(the_face, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L225-229: skip non-degenerated edges.
        if !brep_tool_degenerated(&an_edge) {
            continue;
        }
        // OCCT L230-231: TopExp::Vertices(anEdge, aV1, aV2) (CumOri=false).
        let (a_v1, a_v2) = top_exp_vertices_raw(&an_edge);
        let a_v1 = match a_v1 {
            Some(v) => v,
            None => continue,
        };
        let a_v2 = match a_v2 {
            Some(v) => v,
            None => continue,
        };
        // OCCT L232-235: skip edges whose extremities are not the same.
        if !a_v1.is_same(&a_v2) {
            continue;
        }

        // OCCT L237-238.
        a_degen_pnt.push(brep_tool_pnt(&a_v1).unwrap_or(DVec3::ZERO));
        a_degen_tol.push(brep_tool_tolerance(&a_v1));
    }

    // OCCT L242-245: iterate by sides of the surface (the BSpline guard is
    // the enum match above); the poles Array2 -> the rcad control-point grid
    // [u][v] (LowerRow/UpperRow = 0/nu-1, LowerCol/UpperCol = 0/nv-1;
    // RowLength = nv, ColLength = nu).  The OCCT Array2 always carries poles;
    // the rcad guard avoids the 0-1 underflow on an empty grid.
    let a_poles = &a_bspline.control_points;
    let nb_u = a_poles.len();
    let nb_v = if nb_u > 0 { a_poles[0].len() } else { 0 };
    if nb_u == 0 || nb_v == 0 {
        return the_surface.clone();
    }

    // OCCT L247: occ::handle<Geom_BSplineSurface> aCopy; (null until moved).
    let mut a_copy: Option<rcad_kernel::geom::BSplineSurface> = None;

    // OCCT L249-257: iterate by sides: {U=0; V=0; U=1; V=1}.
    let row_start: [usize; 4] = [0, 0, nb_u - 1, 0];
    let col_start: [usize; 4] = [0, 0, 0, nb_v - 1];
    let row_step: [usize; 4] = [0, 1, 0, 1];
    let col_step: [usize; 4] = [1, 0, 1, 0];
    let nb_steps: [usize; 4] = [nb_v, nb_u, nb_v, nb_u];

    for i_side in 0..4 {
        // OCCT L260-268: compute the center of gravity of the side poles.
        let mut a_sum = DVec3::ZERO;
        for i_pole in 0..nb_steps[i_side] {
            a_sum += a_poles[row_start[i_side] + i_pole * row_step[i_side]]
                [col_start[i_side] + i_pole * col_step[i_side]];
        }
        let a_center = a_sum / nb_steps[i_side] as f64;

        // OCCT L270-272: determine if all poles of the side fit into:
        let mut is_collapsed = true; // aCenter precisely (with gp::Resolution())
        let mut is_singular = true; // aCenter with thePrecision
        // OCCT L274-279: NCollection_LocalArray<bool,4>
        // isDegenerated(aDegenPnt.Extent()) — all entries initialized true.
        let mut is_degenerated: Vec<bool> = vec![true; a_degen_pnt.len()];
        // OCCT L280-307: the per-pole walk.
        for i_pole in 0..nb_steps[i_side] {
            let a_pole = &a_poles[row_start[i_side] + i_pole * row_step[i_side]]
                [col_start[i_side] + i_pole * col_step[i_side]];

            // OCCT L285-290: distance from CG.
            let a_dist_cg = a_center.distance(*a_pole);
            if a_dist_cg > f64::MIN_POSITIVE {
                // OCCT gp::Resolution() == RealSmall() == DBL_MIN.
                is_collapsed = false;
            }
            if a_dist_cg > the_precision {
                is_singular = false;
            }

            // OCCT L296-306: distances from the degenerated points (the two
            // NCollection_List iterators advance in lockstep with iDegen).
            for (i_degen, a_deg_pnt) in a_degen_pnt.iter().enumerate() {
                let a_deg_tol = a_degen_tol[i_degen];
                if is_degenerated[i_degen] && a_deg_pnt.distance(*a_pole) >= a_deg_tol {
                    is_degenerated[i_degen] = false;
                }
            }
        }
        if is_collapsed {
            // OCCT L308-311: already Ok, nothing to be done.
            continue;
        }

        // OCCT L313-315: decide to collapse the side: either if it is
        // singular with thePrecision, or if it fits into one (and only one)
        // degenerated point.
        if !is_singular {
            let mut a_nb_fit = 0;
            // OCCT L317-335: the lockstep list walk with in-flight removal
            // (Remove(Iterator&) drops the current item and advances; the
            // rcad vehicle keeps a read index over the mutating lists).
            let mut a_it_idx = 0usize;
            for i_degen in 0..is_degenerated.len() {
                if is_degenerated[i_degen] {
                    // OCCT L324-327: remove the degenerated point as soon as
                    // it fits at least one side, to prevent total collapse.
                    a_degen_pnt.remove(a_it_idx);
                    a_degen_tol.remove(a_it_idx);
                    a_nb_fit += 1;
                } else {
                    // OCCT L330-334: aDegPntIt.Next(); aDegTolIt.Next().
                    a_it_idx += 1;
                }
            }

            // OCCT L337-339: if the side fits more than one degenerated
            // vertex, do not collapse it — to be on the safe side.
            is_singular = a_nb_fit == 1;
        }

        // OCCT L342-355: do collapse.
        if is_singular {
            if a_copy.is_none() {
                // OCCT L345-348: aCopy = down_cast(theSurface->Copy()).
                a_copy = Some(a_bspline.clone());
            }
            let the_copy = a_copy.as_mut().expect("aCopy");
            for i_pole in 0..nb_steps[i_side] {
                // OCCT L351-353: aCopy->SetPole(RowStart + i*RowStep,
                //               ColStart + i*ColStep, aCenter).
                let r = row_start[i_side] + i_pole * row_step[i_side];
                let c = col_start[i_side] + i_pole * col_step[i_side];
                the_copy.control_points[r][c] = a_center;
            }
        }
    }

    // OCCT L358-364.
    match a_copy {
        Some(a_copy) => Surface3::BSpline(a_copy),
        None => the_surface.clone(),
    }
}

// ---------------------------------------------------------------------------
// Frame vehicles (architecture differences #2/#3).
// ---------------------------------------------------------------------------

/// OCCT Geom_CylindricalSurface::Position() — the flattened gp_Ax3:
/// (location, main direction, X direction, optional Y direction).
fn cylinder_ax3(c: &CylindricalSurface) -> (DVec3, DVec3, DVec3, Option<DVec3>) {
    (c.origin, c.axis, c.ref_dir, c.y_dir)
}

/// OCCT gp_Ax3::Direct() on the cylinder frame — (X ^ Y) . Z > 0
/// (gp_Ax3.cxx); the explicit y_dir (the left-handed swept-lateral frame)
/// participates in the sense.
fn cylinder_ax3_direct(c: &CylindricalSurface) -> bool {
    let y = c.y_axis();
    c.ref_dir.cross(y).dot(c.axis) > 0.0
}

/// OCCT Geom_SphericalSurface::Position() — the flattened gp_Ax3.
fn sphere_ax3(s: &SphericalSurface) -> (DVec3, DVec3, DVec3) {
    (s.center, s.axis, s.ref_dir)
}

/// OCCT gp_Ax3::Direct() on the sphere frame — (X ^ Y) . Z > 0.
fn sphere_ax3_direct(s: &SphericalSurface) -> bool {
    let y = s.ref_dir_perp();
    s.ref_dir.cross(y).dot(s.axis) > 0.0
}

/// OCCT Geom_ConicalSurface::Position() — the flattened gp_Ax3 (the rcad
/// `apex` is the gp_Cone Location point: the axis point at RefRadius).
fn cone_ax3(c: &ConicalSurface) -> (DVec3, DVec3, DVec3) {
    (c.apex, c.axis, c.ref_dir)
}

/// OCCT Geom_ToroidalSurface::Position() — the flattened gp_Ax3.
fn torus_ax3(s: &ToroidalSurface) -> (DVec3, DVec3, DVec3) {
    (s.center, s.axis, s.ref_dir)
}

/// OCCT gp_Ax3::Direct() on the torus frame — (X ^ Y) . Z > 0
/// (the Y direction is axis x ref_dir by the rcad convention).
fn torus_ax3_direct(s: &ToroidalSurface) -> bool {
    let y = s.axis.cross(s.ref_dir);
    s.ref_dir.cross(y).dot(s.axis) > 0.0
}

/// OCCT gp_Ax3::Rotate(gp_Ax1(Location, Direction), M_PI) applied to a frame
/// direction — the exact angle-PI rotation v' = 2(v.n)n - v (architecture
/// difference #2: the OCCT vehicle is the gp_Trsf quaternion matrix).
fn rotate_pi_about(v: DVec3, n: DVec3) -> DVec3 {
    2.0 * v.dot(n) * n - v
}
