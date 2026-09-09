//! OCCT GeomLib_IsPlanarSurface (TKGeomBase/GeomLib —
//! `GeomLib_IsPlanarSurface.hxx` + `GeomLib_IsPlanarSurface.cxx` L33-286)
//! — finds if a surface is a planar surface (the planarity probe of
//! BRepOffset_MakeOffset::LinearizeFaces).
//!
//! Architecture bridges (the brep_offset_make_offset.rs numbering style):
//! 1. `occ::handle<Geom_Surface>` -> `&Surface3` (clone-on-handle; the
//!    geom_lib_check_curve_on_surface.rs precedent).
//! 2. `gp_Pnt` / `gp_Vec` / `gp_Dir` -> `DVec3`; `gp_Dir` construction
//!    normalizes (`gp_Dir(const gp_Vec&)` raises for a null vector — the
//!    OCCT failure path is the panic of the normalize form).
//! 3. `gp_Pln` -> `rcad_kernel::geom::Plane`; `Plan.SetPosition(axe)` +
//!    `Plan.SetLocation(Bary)` fold into the Plane value construction
//!    (arch. diff. #53: the Ax3 position form folds into the surface struct
//!    fields; `Plane::with_axes` is the `gp_Pln(gp_Ax3(P, N, Vx))` form).
//! 4. `GeomAdaptor_Surface` -> the local type dispatch + the Bounds / D1 /
//!    D0 / NbUIntervals / NbVIntervals re-hosts (the SurfaceEval /
//!    GeomAdaptor translation precedents in rcad_kernel::geom::eval).
//! 5. `GeomAdaptor_Curve` -> the `Curve3` adaptor-layer re-hosts
//!    (`nb_intervals` — GeomAdaptor_Curve.cxx L371-460 translation; the
//!    `GetType` / `NbPoles` / parameter reads fold into the Curve3 match).
//! 6. `NCollection_Array1<gp_Pnt> Poles(1, NbU*NbV)` -> `Vec<DVec3>` (the
//!    OCCT 1-based indexing maps to the Vec order).
//! 7. `GeomLib::Inertia` -> `rcad_kernel::base::geom_lib::inertia` (the
//!    GeomLib.cxx L1976-2093 translation).

use glam::DVec3;
use rcad_kernel::base::geom_lib::inertia;
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve3, CurveEval, Plane, Surface3, SurfaceEval};
use rcad_kernel::math::GeomAbsShape;

// OCCT gp.hxx L59-60 / Standard_Real.hxx L146-151:
// gp::Resolution() = RealSmall() = DBL_MIN.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT gp_Dir::Angle(theOther) (gp_Dir.hxx/gp_Dir.lxx): the angle in
/// radians between the two unit directions — acos of the clamped dot of the
/// normalized carriers.
fn dir_angle(a: &DVec3, b: &DVec3) -> f64 {
    a.normalize_or_zero()
        .dot(b.normalize_or_zero())
        .clamp(-1.0, 1.0)
        .acos()
}

/// OCCT gp_Pln::Distance(P) (gp_Pln.cxx): the absolute value of the signed
/// distance of the point to the plane.
fn gp_pln_distance(plan: &Plane, p: &DVec3) -> f64 {
    (p - plan.origin).dot(plan.normal).abs()
}

/// OCCT Geom_Surface::Bounds — the rcad SurfaceEval stand-in (the
/// brep_offset_offset_b.rs surface_bounds4 form).
fn surface_bounds4(the_s: &Surface3) -> (f64, f64, f64, f64) {
    let [u1, u2, v1, v2] = the_s.default_domain();
    (u1, u2, v1, v2)
}

/// OCCT GeomAdaptor_Surface::NbUIntervals(S) (GeomAdaptor_Surface.cxx
/// L643-697).
///
/// GAP leaf (arch. diff. #52): the GeomAbs_BSplineSurface case builds the
/// V-iso basis curve at the first V knot and feeds
/// GeomAdaptor_Curve::NbIntervals — the Geom_BSplineSurface::VIso
/// construction is not translated (the brep_offset_offset_b.rs
/// surface_u_iso_gap precedent), so the case keeps the OCCT failure path.
fn nb_u_intervals(s: &Surface3, cont: GeomAbsShape) -> usize {
    match s {
        Surface3::BSpline(_) => {
            panic!(
                "GAP: GeomAdaptor_Surface::NbUIntervals (BSpline case needs \
                 Geom_BSplineSurface::VIso — iso-curve construction not translated)"
            );
        }
        Surface3::LinearExtrusion(ext) => {
            // OCCT L654-663: the basis curve; BSpline only, else 1.
            if matches!(&*ext.profile, Curve3::BSpline(_)) {
                ext.profile.nb_intervals(cont)
            } else {
                1
            }
        }
        Surface3::Offset(os) => {
            // OCCT L664-685: the BaseS mapping (GeomAbsShape carries no
            // G1/G2 members — the OCCT throw cases do not exist here).
            let base_s = match cont {
                GeomAbsShape::C0 => GeomAbsShape::C1,
                GeomAbsShape::C1 => GeomAbsShape::C2,
                GeomAbsShape::C2 => GeomAbsShape::C3,
                GeomAbsShape::C3 | GeomAbsShape::CN => GeomAbsShape::CN,
            };
            nb_u_intervals(&os.basis, base_s)
        }
        _ => 1,
    }
}

/// OCCT GeomAdaptor_Surface::NbVIntervals(S) (GeomAdaptor_Surface.cxx
/// L701-755).
///
/// GAP leaf (arch. diff. #52): the GeomAbs_BSplineSurface case needs the
/// Geom_BSplineSurface::UIso construction (see nb_u_intervals).
fn nb_v_intervals(s: &Surface3, cont: GeomAbsShape) -> usize {
    match s {
        Surface3::BSpline(_) => {
            panic!(
                "GAP: GeomAdaptor_Surface::NbVIntervals (BSpline case needs \
                 Geom_BSplineSurface::UIso — iso-curve construction not translated)"
            );
        }
        Surface3::Revolution(rev) => {
            // OCCT L710-721: the basis curve; BSpline only, else 1.
            if matches!(&*rev.profile, Curve3::BSpline(_)) {
                rev.profile.nb_intervals(cont)
            } else {
                1
            }
        }
        Surface3::Offset(os) => {
            // OCCT L722-753: the BaseS mapping.
            let base_s = match cont {
                GeomAbsShape::C0 => GeomAbsShape::C1,
                GeomAbsShape::C1 => GeomAbsShape::C2,
                GeomAbsShape::C2 => GeomAbsShape::C3,
                GeomAbsShape::C3 | GeomAbsShape::CN => GeomAbsShape::CN,
            };
            nb_v_intervals(&os.basis, base_s)
        }
        _ => 1,
    }
}

/// OCCT Controle(Poles, Tol, S, Plan) (GeomLib_IsPlanarSurface.cxx L33-89)
/// — the poles-inertia plane fit.
fn controle_poles(poles: &[DVec3], tol: f64, s: &Surface3, plan: &mut Plane) -> bool {
    let mut is_plan = false;
    // OCCT L42: constexpr double aTolSingular = Precision::Confusion();
    let a_tol_singular = CONFUSION;
    // OCCT L44: GeomLib::Inertia(Poles, Bary, DX, DY, gx, gy, gz);
    let mut bary = DVec3::ZERO;
    let mut dx = DVec3::ZERO;
    let mut dy = DVec3::ZERO;
    let mut gx = 0.0;
    let mut gy = 0.0;
    let mut gz = 0.0;
    inertia(poles, &mut bary, &mut dx, &mut dy, &mut gx, &mut gy, &mut gz);
    if gz < tol && gy > a_tol_singular {
        // OCCT L47-51: S->Bounds(umin, umax, vmin, vmax);
        // S->D1((umin + umax) / 2, (vmin + vmax) / 2, P, DU, DV);
        let (umin, umax, vmin, vmax) = surface_bounds4(s);
        let (_p, du, dv) = s.derivatives((umin + umax) / 2.0, (vmin + vmax) / 2.0);

        if du.length_squared() > GP_RESOLUTION && dv.length_squared() > GP_RESOLUTION {
            // Choose DX as close as possible to DU (OCCT L55-80)
            // OCCT L56: gp_Dir du(DU);
            let mut du_dir = du.normalize_or_zero();
            let mut dxx = dx;
            let mut dyy = dy;
            let mut angle1 = dir_angle(&du_dir, &dxx);
            let mut angle2 = dir_angle(&du_dir, &dyy);
            if angle1 > std::f64::consts::FRAC_PI_2 {
                angle1 = std::f64::consts::PI - angle1;
            }
            if angle2 > std::f64::consts::FRAC_PI_2 {
                angle2 = std::f64::consts::PI - angle2;
            }
            if angle2 < angle1 {
                // OCCT L68-72: du = DY; DY = DX; DX = du;
                let tmp = du_dir;
                du_dir = dyy;
                dyy = dxx;
                dxx = tmp;
            }
            // OCCT L73-76: if (DX.Angle(DU) > M_PI / 2) DX.Reverse();
            if dir_angle(&dxx, &du_dir) > std::f64::consts::FRAC_PI_2 {
                dxx = -dxx;
            }
            // OCCT L77-80: if (DY.Angle(DV) > M_PI / 2) DY.Reverse();
            if dir_angle(&dyy, &dv) > std::f64::consts::FRAC_PI_2 {
                dyy = -dyy;
            }

            // OCCT L82-85: gp_Ax3 axe(Bary, DX ^ DY, DX);
            // Plan.SetPosition(axe); Plan.SetLocation(Bary); IsPlan = true;
            *plan = Plane::with_axes(bary, dxx.cross(dyy), dxx);
            is_plan = true;
        }
    }
    is_plan
}

/// OCCT Controle(C, Plan, Tol) (GeomLib_IsPlanarSurface.cxx L91-141) — the
/// curve-on-plane sampling probe.
fn controle_curve(c: &Curve3, plan: &Plane, tol: f64) -> bool {
    let mut b = true;
    // OCCT L96-97: GeomAdaptor_Curve AC(C); Type = AC.GetType();
    // OCCT L99-127: the Nb switch (NbPoles = the poles count; the default
    // branch is 8 + 3 * AC.NbIntervals(GeomAbs_CN)).
    let nb: i32 = match c {
        Curve3::Line(_) => 2,
        Curve3::Circle(_) => 3,
        Curve3::Ellipse(_) | Curve3::Hyperbola(_) | Curve3::Parabola(_) => 5,
        Curve3::Bezier(bc) => bc.control_points.len() as i32,
        Curve3::BSpline(bs) => bs.control_points.len() as i32,
        _ => (8 + 3 * c.nb_intervals(GeomAbsShape::CN)) as i32,
    };

    // OCCT L129-132: f = AC.FirstParameter(); l = AC.LastParameter();
    // du = (l - f) / (Nb - 1);
    let domain = c.default_domain();
    let f = domain[0];
    let l = domain[1];
    let du = (l - f) / (nb - 1) as f64;
    // OCCT L133-138: for (ii = 1; ii <= Nb && B; ii++)
    //   u = (ii - 1) * du + f; d = Plan.Distance(C->Value(u)); B = d < Tol;
    let mut ii = 1;
    while ii <= nb && b {
        let u = (ii - 1) as f64 * du + f;
        let d = gp_pln_distance(plan, &c.point_at(u));
        b = d < tol;
        ii += 1;
    }

    b
}

/// OCCT GeomLib_IsPlanarSurface (GeomLib_IsPlanarSurface.hxx L28-45).
pub struct GeomLibIsPlanarSurface {
    /// OCCT: myPlan.
    my_plan: Plane,
    /// OCCT: IsPlan.
    is_plan: bool,
}

impl GeomLibIsPlanarSurface {
    /// OCCT GeomLib_IsPlanarSurface::GeomLib_IsPlanarSurface(S, Tol = 1.0e-7)
    /// (GeomLib_IsPlanarSurface.cxx L143-272).
    pub fn new(s: &Surface3, tol: f64) -> Self {
        // OCCT L147: GeomAdaptor_Surface AS(S); — the default gp_Pln of the
        // myPlan member is the OXYZ frame (gp_Ax3 default).
        let mut my_plan = Plane::new(DVec3::ZERO, DVec3::Z);
        // OCCT L148-150: Type = AS.GetType();
        let is_plan: bool;

        match s {
            Surface3::Plane(pl) => {
                // OCCT L154-157: case GeomAbs_Plane: IsPlan = true; myPlan = AS.Plane();
                my_plan = *pl;
                is_plan = true;
            }
            Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_)
            | Surface3::Torus(_) => {
                // OCCT L159-165
                is_plan = false;
            }
            Surface3::Revolution(rev) => {
                // OCCT L167-209: case GeomAbs_SurfaceOfRevolution
                // OCCT L168: bool Essai = true;
                let mut essai = true;
                // OCCT L171: gp_Dir Dir = AS.AxeOfRevolution().Direction();
                let mut dir = rev.axis_dir;
                // OCCT L172-174: S->Bounds(...); S->D1(mid, mid, P, DU, DV);
                let (umin, umax, vmin, vmax) = surface_bounds4(s);
                let (mut p, mut du, mut dv) =
                    s.derivatives((umin + umax) / 2.0, (vmin + vmax) / 2.0);
                // OCCT L175-180: the degenerate-D1 retry shifted by 10%.
                if du.length() <= GP_RESOLUTION || dv.length() <= GP_RESOLUTION {
                    let new_u = (umin + umax) / 2.0 + (umax - umin) * 0.1;
                    let new_v = (vmin + vmax) / 2.0 + (vmax - vmin) * 0.1;
                    let r = s.derivatives(new_u, new_v);
                    p = r.0;
                    du = r.1;
                    dv = r.2;
                }
                // OCCT L181: Dn = DU ^ DV;
                let dn = du.cross(dv);
                // OCCT L182-191
                if dn.length() > 1.0e-7 {
                    let mut angle = dir_angle(&dir, &dn);
                    if angle > std::f64::consts::FRAC_PI_2 {
                        angle = std::f64::consts::PI - angle;
                        // OCCT L188: Dir.Reverse();
                        dir = -dir;
                    }
                    // OCCT L190: Essai = (angle < 0.1);
                    essai = angle < 0.1;
                }

                if essai {
                    // OCCT L195-198: gp_Ax3 axe(P, Dir); axe.SetXDirection(DU);
                    // myPlan.SetPosition(axe); myPlan.SetLocation(P);
                    my_plan = Plane::with_axes(p, dir, du);
                    // OCCT L199-201: C = S->UIso(Umin);
                    // IsPlan = Controle(C, myPlan, Tol);
                    let c = surface_of_revolution_u_iso(s, umin);
                    is_plan = controle_curve(&c, &my_plan, tol);
                } else {
                    // OCCT L203-206
                    is_plan = false;
                }
            }
            Surface3::LinearExtrusion(ext) => {
                // OCCT L210-248: case GeomAbs_SurfaceOfExtrusion
                // OCCT L211: bool Essai = false;
                let mut essai = false;
                // OCCT L212-218
                let (umin, umax, vmin, vmax) = surface_bounds4(s);
                let (mut p, mut du, mut dv) =
                    s.derivatives((umin + umax) / 2.0, (vmin + vmax) / 2.0);
                // OCCT L219-224
                if du.length() <= GP_RESOLUTION || dv.length() <= GP_RESOLUTION {
                    let new_u = (umin + umax) / 2.0 + (umax - umin) * 0.1;
                    let new_v = (vmin + vmax) / 2.0 + (vmax - vmin) * 0.1;
                    let r = s.derivatives(new_u, new_v);
                    p = r.0;
                    du = r.1;
                    dv = r.2;
                }
                // OCCT L225-233: Dn = Du ^ Dv; norm = Dn.Magnitude();
                let mut dn = du.cross(dv);
                let norm = dn.length();
                if norm > 1.0e-15 {
                    // OCCT L229: Dn /= norm;
                    dn /= norm;
                    // OCCT L230: double angmax = Tol / (Vmax - Vmin);
                    let angmax = tol / (vmax - vmin);
                    // OCCT L231-232: gp_Dir D(Dn);
                    // Essai = (D.IsNormal(AS.Direction(), angmax));
                    // gp_Dir::IsNormal (gp_Dir.hxx L171-179):
                    // |M_PI/2 - Angle| <= theAngularTolerance.
                    let d = dn.normalize_or_zero();
                    essai = (std::f64::consts::FRAC_PI_2 - dir_angle(&d, &ext.direction)).abs()
                        <= angmax;
                }
                if essai {
                    // OCCT L235-238: gp_Ax3 axe(P, Dn, Du);
                    // myPlan.SetPosition(axe); myPlan.SetLocation(P);
                    my_plan = Plane::with_axes(p, dn, du);
                    // OCCT L239-241: C = S->VIso((Vmin + Vmax) / 2);
                    // IsPlan = Controle(C, myPlan, Tol);
                    let c = surface_of_extrusion_v_iso(s, (vmin + vmax) / 2.0);
                    is_plan = controle_curve(&c, &my_plan, tol);
                } else {
                    // OCCT L243-246
                    is_plan = false;
                }
            }
            _ => {
                // OCCT L250-270: the sampled-poles default (BSpline, Bezier,
                // Offset, ...).
                // OCCT L251-253: NbU = 8 + 3 * AS.NbUIntervals(GeomAbs_CN);
                // NbV = 8 + 3 * AS.NbVIntervals(GeomAbs_CN);
                let nb_u = (8 + 3 * nb_u_intervals(s, GeomAbsShape::CN)) as usize;
                let nb_v = (8 + 3 * nb_v_intervals(s, GeomAbsShape::CN)) as usize;
                // OCCT L254-257
                let (umin, umax, vmin, vmax) = surface_bounds4(s);
                let du = (umax - umin) / (nb_u - 1) as f64;
                let dv = (vmax - vmin) / (nb_v - 1) as f64;
                // OCCT L258-267: NCollection_Array1<gp_Pnt> Pnts(1, NbU * NbV);
                // for (ii = 0, kk = 1; ii < NbU; ii++) { U = Umin + du * ii;
                //   for (jj = 0; jj < NbV; jj++, kk++) { V = Vmin + dv * jj;
                //     S->D0(U, V, Pnts(kk)); } }
                let mut pnts: Vec<DVec3> = Vec::with_capacity(nb_u * nb_v);
                for ii in 0..nb_u {
                    let u = umin + du * ii as f64;
                    for jj in 0..nb_v {
                        let v = vmin + dv * jj as f64;
                        pnts.push(s.point_at(u, v));
                    }
                }

                // OCCT L269: IsPlan = Controle(Pnts, Tol, S, myPlan);
                is_plan = controle_poles(&pnts, tol, s, &mut my_plan);
            }
        }

        GeomLibIsPlanarSurface { my_plan, is_plan }
    }

    /// OCCT GeomLib_IsPlanarSurface::IsPlanar()
    /// (GeomLib_IsPlanarSurface.cxx L274-277) — return if the Surface is a plan.
    pub fn is_planar(&self) -> bool {
        self.is_plan
    }

    /// OCCT GeomLib_IsPlanarSurface::Plan()
    /// (GeomLib_IsPlanarSurface.cxx L279-286) — return the plan definition.
    pub fn plan(&self) -> Plane {
        if !self.is_plan {
            // OCCT L283: throw StdFail_NotDone(" GeomLib_IsPlanarSurface");
            panic!("StdFail_NotDone:  GeomLib_IsPlanarSurface");
        }
        self.my_plan
    }
}

/// OCCT Geom_SurfaceOfRevolution::UIso(U) (Geom_SurfaceOfRevolution.cxx) —
/// GAP leaf (arch. diff. #52; the brep_offset_offset_b.rs surface_u_iso_gap
/// precedent): the iso-curve construction is not translated; the OCCT
/// failure path is kept.
fn surface_of_revolution_u_iso(the_s: &Surface3, the_u: f64) -> Curve3 {
    let _ = (the_s, the_u);
    panic!("GAP: Geom_SurfaceOfRevolution::UIso (iso-curve construction not translated)");
}

/// OCCT Geom_SurfaceOfLinearExtrusion::VIso(V)
/// (Geom_SurfaceOfLinearExtrusion.cxx) — GAP leaf (arch. diff. #52; the
/// brep_offset_offset_b.rs surface_v_iso_gap precedent).
fn surface_of_extrusion_v_iso(the_s: &Surface3, the_v: f64) -> Curve3 {
    let _ = (the_s, the_v);
    panic!("GAP: Geom_SurfaceOfLinearExtrusion::VIso (iso-curve construction not translated)");
}
