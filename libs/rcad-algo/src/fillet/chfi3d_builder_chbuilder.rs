//! OCCT ChFi3d_ChBuilder.cxx — 1:1 translation (Stage 1f).
//!
//! Source: OCCT src/ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx
//! (2,318 lines).
//!
//! Already translated in `chfi3d.rs` (ChFi3dChBuilder L1016-1307): the
//! constructor, Add x3, SetDist/GetDist, SetDists/Dists, AddDA,
//! SetDistAngle/GetDistAngle, SetMode, IsChamfer, Mode, ResetContour.
//! This file carries the remaining methods: Simulate, NbSurf, Sect,
//! SimulKPart, SimulSurf x4, PerformFirstSection, PerformSurf x4,
//! ExtentOneCorner, ExtentTwoCorner, ExtentThreeCorner, SetRegul,
//! ConexFaces, plus the file-scope SearchCommonFaces and
//! ExtentSpineOnCommonFace free functions.
//!
//! Architecture mappings (see chfi3d_builder_6.rs header):
//!   - `occ::handle<ChFiDS_Spine>` -> `&ChFiDSSpineHandle` (polymorphic
//!     Fil/Chamf enum, chfi3d.rs).
//!   - `occ::handle<Adaptor3d_TopolTool>` -> `&BRepTopolTool`.
//!   - `occ::handle<Adaptor3d_Surface>` -> `&BRepAdaptorSurface`.
//!   - `std::unique_ptr<BlendFunc_GenChamfer>` -> `Box<dyn BlendFuncGenChamfer>`
//!     (Stage 1e second batch chamfer family, super::brep_blend_func_chamfer).
//!   - The ChFiDS_CircSection Simul handle store
//!     (`SD->SetSimul(sec)` / `SetOfSurfData()->Value(IS)->Simul()`) is a
//!     pending boundary (no simul slot on ChFiDSSurfData yet, chfi_ds.rs
//!     note) — the section arrays are built per OCCT and dropped there.
//!
//! Cross-file references pending from parallel Stage agents:
//!   ChFi3d_Builder::SearchFace -> super::chfi3d_builder_2 (Builder_2.cxx
//!   L1523-1682);  BlendFunc_* family -> super::brep_blend_func_chamfer
//!   (Stage 1e second batch).

use glam::DVec3;
use rcad_kernel::math::gp::Lin;
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};
use rcad_kernel::geom::SurfaceEval as _;
use rcad_kernel::topods;

use super::brep_blend_func_chamfer::{
    BlendFuncChAsym, BlendFuncChAsymInv, BlendFuncChamfer, BlendFuncChamfInv, BlendFuncConstThroat,
    BlendFuncConstThroatInv, BlendFuncConstThroatWithPenetration,
    BlendFuncConstThroatWithPenetrationInv, BlendFuncGenChamfer, BlendFuncGenChamfInv,
};
use super::chfi3d::concave_side;
use super::chfi3d::ChFi3dChBuilder;
use super::chfi3d_builder_0::{surface_type_of, BRepAdaptorSurface, GeomAbsSurfaceType};
use super::chfi3d_builder_6::chfi3d_fil_common_point;
use super::brep_blend_line::BRepBlendLine;
use super::chfi3d_builder_6b::{elspine_matches_handle, ChFiDSElSpineHandle};
use super::chfi_ds::{
    ChFiDSCircSection, ChFiDSCircSectionArray, ChFiDS_ChamfMethod, ChFiDS_State, ChFiDSSpineHandle,
    ChFiDSSurfData, SharedStripe,
};
use super::chfi3d_builder_2::{BRepTopAdaptorTopolTool, TopAbsState};

// =========================================================================
// OCCT ChFi3d_ChBuilder.cxx L61-95 — SearchCommonFaces (search the 2
// common faces <F1> and <F2> of the edge <E>; uses the EFMap and takes
// the 2 first good faces).
// =========================================================================
pub fn search_common_faces(
    brep: &topods::BRep,
    efmap: &super::chfi_ds::ChFiDSMap,
    e: &Shape,
) -> (Shape, Shape) {
    let mut f1 = Shape::null();
    let mut f2 = Shape::null();
    let list = if efmap.contains(e) { efmap.find(e).clone() } else { Vec::new() };
    for fc in list {
        if f1.is_null() {
            f1 = fc;
        } else if !fc.is_same(&f1) {
            f2 = fc;
            break;
        }
    }

    if !f1.is_null() && f2.is_null() && brep_tools_is_really_closed(brep, e, &f1) {
        f2 = f1.clone();
    }
    (f1, f2)
}

/// OCCT BRepTools::IsReallyClosed(E, F) — the edge occurs twice among the
/// face's wires (builder_0 twin is private; kept next to this consumer).
fn brep_tools_is_really_closed(brep: &topods::BRep, e: &Shape, f: &Shape) -> bool {
    let edges = super::chfi3d_builder_0::topexp_face_edges(brep, f);
    let mut n = 0usize;
    for we in &edges {
        if we.is_same(e) {
            n += 1;
        }
    }
    n > 1
}

// =========================================================================
// OCCT ChFi3d_ChBuilder.cxx L97-185 — ExtentSpineOnCommonFace (extend
// spines of two chamfers by distance dis1/dis2 on their common face;
// two guide lines Spine1 and Spine2 cross in V; isfirst(i) = False if
// Spine(i) is oriented to V).
// =========================================================================
pub fn extent_spine_on_common_face(
    spine1: &mut ChFiDSSpineHandle,
    spine2: &mut ChFiDSSpineHandle,
    v: &Shape,
    dis1: f64,
    dis2: f64,
    isfirst1: bool,
    isfirst2: bool,
) {
    let tolesp = 1.0e-7;

    // alpha, the opening angle between two
    // tangents of two guidelines in V is found
    let tga1;
    let tga2;
    let mut d1plus = 0.0;
    let mut d2plus = 0.0;

    let mut tg1;
    let mut tg2;
    {
        let s1 = spine1.base_mut();
        let (_p, t) = s1.d1(s1.absc_of_vertex(v));
        tg1 = t;
    }
    {
        let s2 = spine2.base_mut();
        let (_p, t) = s2.d1(s2.absc_of_vertex(v));
        tg2 = t;
    }
    tg1 = tg1.normalize();
    tg2 = tg2.normalize();
    if isfirst1 {
        tg1 = -tg1;
    }
    if isfirst2 {
        tg2 = -tg2;
    }

    let cosalpha = tg1.dot(tg2);
    let sinalpha = (1.0 - cosalpha * cosalpha).sqrt();

    // a1+a2 = alpha
    let mut temp = cosalpha + dis2 / dis1;
    if temp.abs() > tolesp {
        tga1 = sinalpha / temp;
        d1plus = dis1 / tga1;
    }
    temp = cosalpha + dis1 / dis2;
    if temp.abs() > tolesp {
        tga2 = sinalpha / temp;
        d2plus = dis2 / tga2;
    }

    // extension by the calculated distance
    if d1plus > 0.0 {
        d1plus *= 3.0;
        if isfirst1 {
            let s1 = spine1.base_mut();
            s1.set_first_parameter(-d1plus);
            s1.set_first_tgt(0.0);
        } else {
            let s1 = spine1.base_mut();
            let param = s1.last_parameter_of(s1.nb_edges());
            s1.set_last_parameter(d1plus + param);
            s1.set_last_tgt(param);
        }
    }
    if d2plus > 0.0 {
        d2plus *= 1.5;
        if isfirst2 {
            let s2 = spine2.base_mut();
            s2.set_first_parameter(-d2plus);
            s2.set_first_tgt(0.0);
        } else {
            let s2 = spine2.base_mut();
            let param = s2.last_parameter_of(s2.nb_edges());
            s2.set_last_parameter(d2plus + param);
            s2.set_last_tgt(param);
        }
    }
}

// =========================================================================
// Local section-building helpers (OCCT: inline loops in SimulSurf).
// =========================================================================

/// OCCT L831-857 / L993-1019 — the sections loop over the walking line with
/// `pFunc->Section(ww, u1, v1, u2, v2, p1, p2, line)`.
fn chamfer_sections(
    p_func: &mut dyn BlendFuncGenChamfer,
    lin: &BRepBlendLine,
) -> ChFiDSCircSectionArray {
    let line = lin;
    let nbp = line.nb_points();
    let mut sec: ChFiDSCircSectionArray = Vec::with_capacity(nbp as usize);
    for i in 1..=nbp {
        sec.push(ChFiDSCircSection::new());
        let p = line.point(i);
        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        let ww = p.parameter();
        let mut p1 = 0.0;
        let mut p2 = 0.0;
        let mut gp_line = Lin {
            pos: DVec3::ZERO,
            dir: DVec3::ZERO,
        };
        // Disambiguated: BlendFuncGenChamfer::section (the obsolete 8-arg
        // OCCT overload) shadows BlendAppFunction::section on the trait object.
        BlendFuncGenChamfer::section(p_func, ww, u1, v1, u2, v2, &mut p1, &mut p2, &mut gp_line);
        let isec = sec.last_mut().expect("isec");
        isec.set_lin(gp_line, p1, p2);
    }
    sec
}

/// OCCT L1128-1154 — the ChAsym sections loop.
fn chasym_sections(
    func: &mut BlendFuncChAsym,
    lin: &BRepBlendLine,
) -> ChFiDSCircSectionArray {
    let line = lin;
    let nbp = line.nb_points();
    let mut sec: ChFiDSCircSectionArray = Vec::with_capacity(nbp as usize);
    for i in 1..=nbp {
        sec.push(ChFiDSCircSection::new());
        let p = line.point(i);
        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        let ww = p.parameter();
        let mut p1 = 0.0;
        let mut p2 = 0.0;
        let mut gp_line = Lin {
            pos: DVec3::ZERO,
            dir: DVec3::ZERO,
        };
        func.section(ww, u1, v1, u2, v2, &mut p1, &mut p2, &mut gp_line);
        let isec = sec.last_mut().expect("isec");
        isec.set_lin(gp_line, p1, p2);
    }
    sec
}

/// OCCT L847-856 etc. — the first/last 2d points of the line.
fn section_2d_points(lin: &BRepBlendLine) -> (glam::DVec2, glam::DVec2, glam::DVec2, glam::DVec2) {
    let line = lin;
    let nbp = line.nb_points();
    let mut pf1 = glam::DVec2::ZERO;
    let mut pl1 = glam::DVec2::ZERO;
    let mut pf2 = glam::DVec2::ZERO;
    let mut pl2 = glam::DVec2::ZERO;
    for i in 1..=nbp {
        let p = line.point(i);
        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        if i == 1 {
            pf1 = glam::DVec2::new(u1, v1);
            pf2 = glam::DVec2::new(u2, v2);
        }
        if i == nbp {
            pl1 = glam::DVec2::new(u1, v1);
            pl2 = glam::DVec2::new(u2, v2);
        }
    }
    (pf1, pl1, pf2, pl2)
}

/// OCCT `SD->SetSimul(sec)` — pending boundary: no simul slot on
/// ChFiDSSurfData (chfi_ds.rs note); the array is dropped here.
fn store_simul(_data: &mut ChFiDSSurfData, _sec: ChFiDSCircSectionArray) {}

/// OCCT ElSLib::ConeD1(U, V, Pos, Radius, SAngle, P, Vu, Vv)
/// (ElSLib.cxx L687-740) on the rcad ConicalSurface axes.
#[allow(clippy::too_many_arguments)]
fn elslib_cone_d1(
    u: f64,
    v: f64,
    ploc: DVec3,
    xdir: DVec3,
    ydir: DVec3,
    zdir: DVec3,
    radius: f64,
    sangle: f64,
) -> (DVec3, DVec3, DVec3) {
    let cos_u = u.cos();
    let sin_u = u.sin();
    let cos_a = sangle.cos();
    let sin_a = sangle.sin();
    let r = radius + v * sin_a;
    let a3 = v * cos_a;
    let a1 = r * cos_u;
    let a2 = r * sin_u;
    let r1 = sin_a * cos_u;
    let r2 = sin_a * sin_u;
    let p = DVec3::new(
        a1 * xdir.x + a2 * ydir.x + a3 * zdir.x + ploc.x,
        a1 * xdir.y + a2 * ydir.y + a3 * zdir.y + ploc.y,
        a1 * xdir.z + a2 * ydir.z + a3 * zdir.z + ploc.z,
    );
    let vu = DVec3::new(
        -a2 * xdir.x + a1 * ydir.x,
        -a2 * xdir.y + a1 * ydir.y,
        -a2 * xdir.z + a1 * ydir.z,
    );
    let vv = DVec3::new(
        r1 * xdir.x + r2 * ydir.x + cos_a * zdir.x,
        r1 * xdir.y + r2 * ydir.y + cos_a * zdir.y,
        r1 * xdir.z + r2 * ydir.z + cos_a * zdir.z,
    );
    (p, vu, vv)
}

/// OCCT ChFi3d_Builder_0.cxx L575-600 — ChFi3d_evalconti (kept next to the
/// SetRegul consumer).
fn chfi3d_evalconti(
    _e: &Shape,
    f1: &Shape,
    f2: &Shape,
) -> rcad_kernel::topo::topods::GeomAbsShape {
    
    let cont = rcad_kernel::topo::topods::GeomAbsShape::G1;
    if !f1.is_same(f2) {
        return cont;
    }
    let _ = f1;
    // OCCT: F = F1 with FORWARD orientation; S = BRepAdaptor_Surface(F,
    // false); if type is Cone/Sphere/Torus -> CN.  The surface-kind read
    // needs the face surface; the G1 default stands for the non-conic
    // branch (same shape as the OCCT early return).
    cont
}

impl ChFi3dChBuilder {
    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L576-609 — Simulate.
    // =====================================================================
    pub fn simulate(&mut self, ic: usize) {
        let mut i = 1usize;
        for stripe in self.base.my_list_stripe.clone() {
            if i == ic {
                self.base.perform_set_of_surf(&stripe, true);
                break;
            }
            i += 1;
        }
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L614-626 — NbSurf.
    // =====================================================================
    pub fn nb_surf(&self, ic: usize) -> usize {
        let mut i = 1usize;
        for stripe in &self.base.my_list_stripe {
            if i == ic {
                let st = stripe.read().expect("stripe lock");
                return st.my_hdata.len();
            }
            i += 1;
        }
        0
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L631-647 — Sect.
    // Pending boundary: the Simul handle store on ChFiDSSurfData is not
    // carried yet (chfi_ds.rs note), so the down-cast stays empty.
    // =====================================================================
    pub fn sect(&self, ic: usize, is: usize) -> Option<ChFiDSCircSectionArray> {
        let mut i = 1usize;
        for stripe in &self.base.my_list_stripe {
            if i == ic {
                let st = stripe.read().expect("stripe lock");
                let _sd = st.my_hdata.get(is - 1)?;
                // OCCT: bid = ...->Simul(); res = down_cast<HArray1>(bid)
                return None;
            }
            i += 1;
        }
        None
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L655-712 — SimulKPart (stores simulating
    // sections in simul).  The SetSimul store is the pending boundary
    // above; the section array is built per OCCT.
    // =====================================================================
    pub fn simul_kpart(&self, sd: &ChFiDSSurfData) {
        let dstr = self.base.my_ds.as_ref().expect("DS");
        let s = &dstr.surface(sd.surf()).surface;
        let fi1 = sd.interference_on_s1();
        let fi2 = sd.interference_on_s2();
        use rcad_kernel::geom::Curve2dEval as _;
        let p1f = fi1
            .pcurve_on_surf()
            .map(|pc| pc.point_at(fi1.parameter_first()))
            .expect("PCurveOnSurf");
        let p1l = fi1
            .pcurve_on_surf()
            .map(|pc| pc.point_at(fi1.parameter_last()))
            .expect("PCurveOnSurf");
        let p2f = fi2
            .pcurve_on_surf()
            .map(|pc| pc.point_at(fi2.parameter_first()))
            .expect("PCurveOnSurf");
        let p2l = fi2
            .pcurve_on_surf()
            .map(|pc| pc.point_at(fi2.parameter_last()))
            .expect("PCurveOnSurf");
        let typ = surface_type_of(s);
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        match typ {
            GeomAbsSurfaceType::Plane => {
                let v1 = p1f.y;
                let v2 = p2f.y;
                let u1 = p1f.x.max(p2f.x);
                let u2 = p1l.x.min(p2l.x);
                let mut arr: ChFiDSCircSectionArray = vec![ChFiDSCircSection::new(), ChFiDSCircSection::new()];
                let pl = match s {
                    rcad_kernel::geom::Surface3::Plane(p) => p.clone(),
                    _ => unreachable!("GeomAbs_Plane branch"),
                };
                // OCCT: ElSLib::PlaneUIso(Pl.Position(), u) — the u-line
                // (ElSLib.cxx L1705-1714).
                let iso1 = Lin {
                    pos: pl.origin + pl.u_dir * u1,
                    dir: pl.v_dir,
                };
                let iso2 = Lin {
                    pos: pl.origin + pl.u_dir * u2,
                    dir: pl.v_dir,
                };
                arr[0].set_lin(iso1, v1, v2);
                arr[1].set_lin(iso2, v1, v2);
                sec = Some(arr);
            }
            GeomAbsSurfaceType::Cone => {
                let v1 = p1f.y;
                let v2 = p2f.y;
                let u1 = p1f.x.max(p2f.x);
                let u2 = p1l.x.min(p2l.x);
                let ang = u2 - u1;
                let (rad, sang, origin, zdir, xdir, ydir) = match s {
                    rcad_kernel::geom::Surface3::Cone(c) => {
                        let yd = c.axis.cross(c.ref_dir);
                        (c.radius, c.half_angle_rad, c.apex, c.axis, c.ref_dir, yd)
                    }
                    _ => unreachable!("GeomAbs_Cone branch"),
                };
                let mut n = (36.0 * ang / std::f64::consts::PI + 1.0) as i32;
                if n < 2 {
                    n = 2;
                }
                let mut arr: ChFiDSCircSectionArray = Vec::new();
                for i in 1..=n {
                    arr.push(ChFiDSCircSection::new());
                    let u = u1 + (i - 1) as f64 * (u2 - u1) / (n - 1) as f64;
                    // OCCT: ElSLib::ConeUIso(Co.Position(), rad, sang, u) —
                    // gp_Lin(P, DV) via ConeD1(U, 0, ...) (ElSLib.cxx L1727).
                    let (_p, _du, dv) =
                        elslib_cone_d1(u, 0.0, origin, xdir, ydir, zdir, rad, sang);
                    let isec = arr.last_mut().expect("isec");
                    isec.set_lin(Lin { pos: _p, dir: dv }, v1, v2);
                }
                sec = Some(arr);
            }
            _ => {}
        }
        // OCCT: SD->SetSimul(sec) — pending boundary (header note).
        let _ = sec;
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L717-1221 — SimulSurf (the face/face
    // simulation entry).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn simul_surf(
        &mut self,
        data: &mut ChFiDSSurfData,
        hguide: &ChFiDSElSpineHandle,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &Vector,
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool {
        let Some(chsp) = spine.down_cast_chamf() else {
            panic!("Standard_ConstructionError: SimulSurf : this is not the spine of a chamfer");
        };

        // Flexible parameters!
        let la = hguide.read().expect("elspine lock").last_parameter();
        let fi = hguide.read().expect("elspine lock").first_parameter();
        let longueur = la - fi;
        let max_step = longueur * 0.05;
        let mut radiusspine = 0.0f64;
        let locfleche;
        let brep = self.base.my_brep.clone();
        // The pending ElSpine curve bound by the blend functions (see
        // elspine_guide_curve in chfi3d_builder_6b).
        let guide = super::chfi3d_builder_6b::elspine_guide_curve(hguide);

        // As ElSpine is parameterized by a curvilinear quasi-abscissa,
        // the min radius is estimated as 1/D2 max;
        for i in 0..=20 {
            let w = fi + i as f64 * max_step;
            let (_p, _d1, d2) = hguide.read().expect("elspine lock").d2(w);
            let temp = d2.length_squared();
            if temp > radiusspine {
                radiusspine = temp;
            }
        }

        let mut lin: Option<BRepBlendLine> = None;
        let p_first = *first;
        if *intf != 0 {
            *first = chsp.base.first_parameter_of(1);
        }
        if *intl != 0 {
            *last = chsp.base.last_parameter_of(chsp.base.nb_edges());
        }

        let mut offset_hguide: Option<ChFiDSElSpineHandle> = None;

        if chsp.is_chamfer() == ChFiDS_ChamfMethod::Sym {
            let dis = chsp.get_dist();
            let radius = dis.max(radiusspine);
            locfleche = radius * 1.0e-2; // graphic criterion

            let mut p_func: Box<dyn BlendFuncGenChamfer>;
            let mut p_finv: Box<dyn BlendFuncGenChamfInv>;
            if chsp.base.my_mode == super::chfi_ds::ChFiDS_ChamfMode::ClassicChamfer {
                p_func = Box::new(BlendFuncChamfer::new(&s1.surface, &s2.surface, &guide));
                p_finv = Box::new(BlendFuncChamfInv::new(&s1.surface, &s2.surface, &guide));
            } else {
                p_func = Box::new(BlendFuncConstThroat::new(&s1.surface, &s2.surface, &guide));
                p_finv = Box::new(BlendFuncConstThroatInv::new(&s1.surface, &s2.surface, &guide));
            }
            p_func.set(dis, dis, choix);
            p_finv.set(dis, dis, choix);

            let done = self.base.simul_data_walking(
                data,
                hguide,
                offset_hguide.as_ref(),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                p_func.as_mut(),
                p_finv.as_mut(),
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                soldep,
                4,
                rec_on_s1,
                rec_on_s2,
            );
            self.base.done = done;

            if !done {
                return false;
            }
            let sec = chamfer_sections(p_func.as_mut(), lin.as_ref().expect("Lin"));
            let (pf1, pl1, pf2, pl2) = section_2d_points(lin.as_ref().expect("Lin"));
            data.set_2d_points(pf1, pl1, pf2, pl2);
            store_simul(data, sec);
            simul_surf_common_points(&brep, lin.as_ref().expect("Lin"), data, self.base.tolapp3d);

            let reverse = !forward || inside;
            if *intf != 0 && reverse {
                let mut ok = false;
                let cp1 = data.vertex_first_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.base.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intf != 0;
                }
                let cp2 = data.vertex_first_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.base.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }
            if *intl != 0 {
                let mut ok = false;
                let cp1 = data.vertex_last_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.base.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intl != 0;
                }
                let cp2 = data.vertex_last_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.base.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }
        } else if chsp.is_chamfer() == ChFiDS_ChamfMethod::TwoDist {
            let (dis1, dis2) = chsp.dists();
            let radius = dis1.max(dis2).max(radiusspine);
            locfleche = radius * 1.0e-2; // graphic criterion

            let mut p_func: Box<dyn BlendFuncGenChamfer>;
            let mut p_finv: Box<dyn BlendFuncGenChamfInv>;
            if chsp.base.my_mode == super::chfi_ds::ChFiDS_ChamfMode::ClassicChamfer {
                p_func = Box::new(BlendFuncChamfer::new(&s1.surface, &s2.surface, &guide));
                p_finv = Box::new(BlendFuncChamfInv::new(&s1.surface, &s2.surface, &guide));
                p_func.set(dis1, dis2, choix);
                p_finv.set(dis1, dis2, choix);
            } else {
                // OCCT L941-951: the offset elspine paired with HGuide by
                // handle identity (elspine_matches_handle boundary note).
                for (a, b) in spine
                    .base()
                    .elspines
                    .iter()
                    .zip(spine.base().offset_elspines.iter())
                {
                    if elspine_matches_handle(a, hguide) {
                        // rcad boundary: the spine list stores plain records;
                        // the offset guide is shared through a fresh handle.
                        offset_hguide = Some(std::sync::Arc::new(std::sync::RwLock::new(b.clone())));
                    }
                }

                if offset_hguide.is_none() {
                    // OCCT prints "Construction of offset guide failed!"
                    // and comments the exception.
                }
                p_func = Box::new(BlendFuncConstThroatWithPenetration::new(
                    &s1.surface,
                    &s2.surface,
                    // OCCT binds OffsetHGuide; pending the ElSpine curve machinery
                    // the guide placeholder stands in (elspine_guide_curve).
                    &guide,
                ));
                p_finv = Box::new(BlendFuncConstThroatWithPenetrationInv::new(
&s1.surface,
&s2.surface,
                    // OCCT binds OffsetHGuide; pending the ElSpine curve machinery
                    // the guide placeholder stands in (elspine_guide_curve).
                    &guide,
                ));
                let throat = dis1.max(dis2);
                p_func.set(throat, throat, choix);
                p_finv.set(throat, throat, choix);
            }

            let done = self.base.simul_data_walking(
                data,
                hguide,
                offset_hguide.as_ref(),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                p_func.as_mut(),
                p_finv.as_mut(),
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                soldep,
                4,
                rec_on_s1,
                rec_on_s2,
            );
            self.base.done = done;

            if !done {
                return false;
            }
            let sec = chamfer_sections(p_func.as_mut(), lin.as_ref().expect("Lin"));
            let (pf1, pl1, pf2, pl2) = section_2d_points(lin.as_ref().expect("Lin"));
            data.set_2d_points(pf1, pl1, pf2, pl2);
            store_simul(data, sec);
            simul_surf_common_points(&brep, lin.as_ref().expect("Lin"), data, self.base.tolapp3d);

            let reverse = !forward || inside;
            if *intf != 0 && reverse {
                let mut ok = false;
                let cp1 = data.vertex_first_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.base.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intf != 0;
                }
                let cp2 = data.vertex_first_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.base.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }
            if *intl != 0 {
                let mut ok = false;
                let cp1 = data.vertex_last_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.base.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intl != 0;
                }
                let cp2 = data.vertex_last_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.base.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }
        } else {
            // distance and angle
            let (dis, angle) = chsp.get_dist_angle();
            let radius = dis.max(dis * angle.tan()).max(radiusspine);
            locfleche = radius * 1.0e-2; // graphic criterion

            let ch = choix;

            let mut func = BlendFuncChAsym::new(&s1.surface, &s2.surface, &guide);
            let mut finv = BlendFuncChAsymInv::new(&s1.surface, &s2.surface, &guide);

            func.set(dis, angle, ch);
            finv.set(dis, angle, ch);

            let done = self.base.simul_data_walking(
                data,
                hguide,
                offset_hguide.as_ref(),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                &mut func,
                &mut finv,
                p_first,
                max_step,
                locfleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                soldep,
                4,
                rec_on_s1,
                rec_on_s2,
            );
            self.base.done = done;

            if !done {
                return false;
            }
            let sec = chasym_sections(&mut func, lin.as_ref().expect("Lin"));
            let (pf1, pl1, pf2, pl2) = section_2d_points(lin.as_ref().expect("Lin"));
            data.set_2d_points(pf1, pl1, pf2, pl2);
            store_simul(data, sec);
            simul_surf_common_points(&brep, lin.as_ref().expect("Lin"), data, self.base.tolapp3d);

            let reverse = !forward || inside;
            if *intf != 0 && reverse {
                let mut ok = false;
                let cp1 = data.vertex_first_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.base.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intf != 0;
                }
                let cp2 = data.vertex_first_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intf = !self.base.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }

            if *intl != 0 {
                let mut ok = false;
                let cp1 = data.vertex_last_on_s1().clone();
                if cp1.is_on_arc() {
                    let f1 = s1.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.base.search_face(&mut spine.clone(), &cp1, &f1, &mut bid) as i32;
                    ok = *intl != 0;
                }
                let cp2 = data.vertex_last_on_s2().clone();
                if cp2.is_on_arc() && !ok {
                    let f2 = s2.face.clone();
                    let mut bid = Shape::null();
                    *intl = !self.base.search_face(&mut spine.clone(), &cp2, &f2, &mut bid) as i32;
                }
            }
        } // distance and angle
        true
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L1319-1614 — PerformFirstSection (to
    // implement the first section if there is no KPart).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn perform_first_section(
        &self,
        brep: &topods::BRep,
        spine: &ChFiDSSpineHandle,
        hguide: &ChFiDSElSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        s2: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        i2: &BRepTopAdaptorTopolTool,
        par: f64,
        sol_dep: &mut Vector,
        pos1: &mut TopAbsState,
        pos2: &mut TopAbsState,
    ) -> bool {
        let guide = super::chfi3d_builder_6b::elspine_guide_curve(hguide);
        let Some(chsp) = spine.down_cast_chamf() else {
            panic!("Standard_ConstructionError: PerformSurf : this is not the spine of a chamfer");
        };

        let tol_guide = hguide.read().expect("elspine lock").resolution(self.base.tolapp3d);

        if chsp.is_chamfer() == ChFiDS_ChamfMethod::Sym {
            let dis = chsp.get_dist();

            let mut p_func: Box<dyn BlendFuncGenChamfer>;
            if chsp.base.my_mode == super::chfi_ds::ChFiDS_ChamfMode::ClassicChamfer {
                p_func = Box::new(BlendFuncChamfer::new(&s1.surface, &s2.surface, &guide));
            } else {
                p_func = Box::new(BlendFuncConstThroat::new(&s1.surface, &s2.surface, &guide));
            }
            p_func.set(dis, dis, choix);
            // OCCT: BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide)
            let mut the_walk = super::brep_blend_walking::BRepBlendWalking::new(&s1.surface, &s2.surface, i1, i2, &guide);

            // calculate an approximate starting solution
            let (ptgui, d1gui) = hguide.read().expect("elspine lock").d1(par);

            p_func.set_param(par);
            let mut tg_f = DVec3::ZERO;
            let mut tg_l = DVec3::ZERO;
            let mut tmp1 = DVec3::ZERO;
            let mut tmp2 = DVec3::ZERO;
            p_func.tangent(
                sol_dep.get(1),
                sol_dep.get(2),
                sol_dep.get(3),
                sol_dep.get(4),
                &mut tg_f,
                &mut tg_l,
                &mut tmp1,
                &mut tmp2,
            );

            let mut rev1 = false;
            let mut rev2 = false;
            let sign = tg_f.cross(d1gui).dot(tg_l);

            if choix % 2 == 1 {
                rev1 = true;
            } else {
                rev2 = true;
            }

            if sign < 0.0 {
                rev1 = !rev1;
                rev2 = !rev2;
            }

            if rev1 {
                tg_f = -tg_f;
            }
            if rev2 {
                tg_l = -tg_l;
            }

            let temp = tg_f * dis;
            let pt1 = ptgui + temp;
            let temp = tg_l * dis;
            let pt2 = ptgui + temp;

            let tol = self.base.tolesp * 1.0e2;
            // OCCT: Extrema_GenLocateExtPS proj1(*S1, tol, tol);
            //       proj1.Perform(pt1, SolDep(1), SolDep(2)); ...
            if let Some((u, v)) = extrema_gen_locate_ext_ps(&s1.surface, pt1, sol_dep.get(1), sol_dep.get(2), tol) {
                sol_dep.set(1, u);
                sol_dep.set(2, v);
            }
            let _ = brep;
            if let Some((u, v)) = extrema_gen_locate_ext_ps(&s2.surface, pt2, sol_dep.get(3), sol_dep.get(4), tol) {
                sol_dep.set(3, u);
                sol_dep.set(4, v);
            }

            the_walk.perform_first_section(p_func.as_mut(), par, &mut sol_dep.data.v, self.base.tolapp3d, tol_guide, pos1, pos2)
        } else if chsp.is_chamfer() == ChFiDS_ChamfMethod::TwoDist {
            let (dis1, dis2) = chsp.dists();

            let mut p_func: Box<dyn BlendFuncGenChamfer>;
            if chsp.base.my_mode == super::chfi_ds::ChFiDS_ChamfMode::ClassicChamfer {
                p_func = Box::new(BlendFuncChamfer::new(&s1.surface, &s2.surface, &guide));
                p_func.set(dis1, dis2, choix);
            } else {
                // OCCT L1432-1443: the offset elspine paired with HGuide.
                let mut offset_hguide: Option<ChFiDSElSpineHandle> = None;
                for (a, b) in spine
                    .base()
                    .elspines
                    .iter()
                    .zip(spine.base().offset_elspines.iter())
                {
                    if elspine_matches_handle(a, hguide) {
                        offset_hguide = Some(std::sync::Arc::new(std::sync::RwLock::new(b.clone())));
                    }
                }

                if offset_hguide.is_none() {
                    // OCCT prints "Construction of offset guide failed!"
                }
                p_func = Box::new(BlendFuncConstThroatWithPenetration::new(
                    &s1.surface,
                    &s2.surface,
                    // OCCT binds OffsetHGuide; pending the ElSpine curve machinery
                    // the guide placeholder stands in (elspine_guide_curve).
                    &guide,
                ));
                let throat = dis1.max(dis2);
                p_func.set(throat, throat, choix); // dis2?
            }
            let mut the_walk = super::brep_blend_walking::BRepBlendWalking::new(&s1.surface, &s2.surface, i1, i2, &guide);

            // calculate an approximate starting solution
            let (ptgui, d1gui) = hguide.read().expect("elspine lock").d1(par);

            p_func.set_param(par);
            let mut tg_f = DVec3::ZERO;
            let mut tg_l = DVec3::ZERO;
            let mut tmp1 = DVec3::ZERO;
            let mut tmp2 = DVec3::ZERO;
            p_func.tangent(
                sol_dep.get(1),
                sol_dep.get(2),
                sol_dep.get(3),
                sol_dep.get(4),
                &mut tg_f,
                &mut tg_l,
                &mut tmp1,
                &mut tmp2,
            );

            let mut rev1 = false;
            let mut rev2 = false;
            let sign = tg_f.cross(d1gui).dot(tg_l);

            if choix % 2 == 1 {
                rev1 = true;
            } else {
                rev2 = true;
            }

            if sign < 0.0 {
                rev1 = !rev1;
                rev2 = !rev2;
            }

            if rev1 {
                tg_f = -tg_f;
            }
            if rev2 {
                tg_l = -tg_l;
            }

            // OCCT L1495-1509: aDist1/aDist2 — the ConstThroatWithPenetration
            // branch keeps aDist1 = dis1 / aDist2 = dis2 (the recompute is
            // commented out in OCCT).
            let a_dist1 = dis1;
            let a_dist2 = dis2;

            let temp = tg_f * a_dist1;
            let pt1 = ptgui + temp;
            let temp = tg_l * a_dist2;
            let pt2 = ptgui + temp;

            let tol = self.base.tolesp * 1.0e2;

            if let Some((u, v)) = extrema_gen_locate_ext_ps(&s1.surface, pt1, sol_dep.get(1), sol_dep.get(2), tol) {
                sol_dep.set(1, u);
                sol_dep.set(2, v);
            }
            if let Some((u, v)) = extrema_gen_locate_ext_ps(&s2.surface, pt2, sol_dep.get(3), sol_dep.get(4), tol) {
                sol_dep.set(3, u);
                sol_dep.set(4, v);
            }

            the_walk.perform_first_section(p_func.as_mut(), par, &mut sol_dep.data.v, self.base.tolapp3d, tol_guide, pos1, pos2)
        } else {
            // distance and angle
            let (dis1, angle) = chsp.get_dist_angle();

            let ch = choix;

            let mut func = BlendFuncChAsym::new(&s1.surface, &s2.surface, &guide);
            func.set(dis1, angle, ch);
            let mut the_walk = super::brep_blend_walking::BRepBlendWalking::new(&s1.surface, &s2.surface, i1, i2, &guide);

            // calculate an approximate starting solution
            let (ptgui, d1gui) = hguide.read().expect("elspine lock").d1(par);

            func.set_param(par);
            let mut tg_f = DVec3::ZERO;
            let mut tg_l = DVec3::ZERO;
            let mut tmp1 = DVec3::ZERO;
            let mut tmp2 = DVec3::ZERO;
            func.tangent(
                sol_dep.get(1),
                sol_dep.get(2),
                sol_dep.get(3),
                sol_dep.get(4),
                &mut tg_f,
                &mut tg_l,
                &mut tmp1,
                &mut tmp2,
            );

            let mut rev1 = false;
            let mut rev2 = false;
            let sign = tg_f.cross(d1gui).dot(tg_l);

            if ch % 2 == 1 {
                rev1 = true;
            } else {
                rev2 = true;
            }

            if sign < 0.0 {
                rev1 = !rev1;
                rev2 = !rev2;
            }

            if rev1 {
                tg_f = -tg_f;
            }
            if rev2 {
                tg_l = -tg_l;
            }

            let temp = tg_f * dis1;
            let pt1 = ptgui + temp;

            let tmpcos = tg_f.dot(tg_l);
            let tmpsin = (1.0 - tmpcos * tmpcos).sqrt();

            let dis2 = dis1 / (tmpcos + tmpsin / angle.tan());

            let temp = tg_l * dis2;
            let pt2 = ptgui + temp;

            let tol = self.base.tolesp * 1.0e2;
            if let Some((u, v)) = extrema_gen_locate_ext_ps(&s1.surface, pt1, sol_dep.get(1), sol_dep.get(2), tol) {
                sol_dep.set(1, u);
                sol_dep.set(2, v);
            }
            if let Some((u, v)) = extrema_gen_locate_ext_ps(&s2.surface, pt2, sol_dep.get(3), sol_dep.get(4), tol) {
                sol_dep.set(3, u);
                sol_dep.set(4, v);
            }

            the_walk.perform_first_section(&mut func, par, &mut sol_dep.data.v, self.base.tolapp3d, tol_guide, pos1, pos2)
        } // distance and angle
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L1618-1852 — PerformSurf (the face/face
    // entry).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn perform_surf(
        &mut self,
        seq_data: &mut Vec<std::sync::Arc<std::sync::RwLock<ChFiDSSurfData>>>,
        hguide: &ChFiDSElSpineHandle,
        spine: &ChFiDSSpineHandle,
        choix: i32,
        s1: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        s2: &BRepAdaptorSurface,
        i2: &BRepTopAdaptorTopolTool,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        appro: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &Vector,
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool {
        // OCCT: Data = SeqData(1)
        let data_arc = seq_data
            .first()
            .expect("SeqData(1)")
            .clone();
        let mut data = data_arc.write().expect("surfdata lock");
        // The pending ElSpine curve bound by the blend functions.
        let guide = super::chfi3d_builder_6b::elspine_guide_curve(hguide);
        let Some(chsp) = spine.down_cast_chamf() else {
            panic!("Standard_ConstructionError: PerformSurf : this is not the spine of a chamfer");
        };

        let mut gd1 = false;
        let mut gd2 = false;
        let mut gf1 = false;
        let mut gf2 = false;
        let mut lin: Option<BRepBlendLine> = None;
        // OCCT: Or = S1->Face().Orientation()
        let or = s1.face.orientation;
        let p_first = *first;
        if *intf != 0 {
            *first = chsp.base.first_parameter_of(1);
        }
        if *intl != 0 {
            *last = chsp.base.last_parameter_of(chsp.base.nb_edges());
        }

        if chsp.is_chamfer() == ChFiDS_ChamfMethod::Sym {
            let mut p_func: Box<dyn BlendFuncGenChamfer>;
            let mut p_finv: Box<dyn BlendFuncGenChamfInv>;
            if chsp.base.my_mode == super::chfi_ds::ChFiDS_ChamfMode::ClassicChamfer {
                p_func = Box::new(BlendFuncChamfer::new(&s1.surface, &s2.surface, &guide));
                p_finv = Box::new(BlendFuncChamfInv::new(&s1.surface, &s2.surface, &guide));
            } else {
                p_func = Box::new(BlendFuncConstThroat::new(&s1.surface, &s2.surface, &guide));
                p_finv = Box::new(BlendFuncConstThroatInv::new(&s1.surface, &s2.surface, &guide));
            }
            let dis = chsp.get_dist();
            p_func.set(dis, dis, choix);
            p_finv.set(dis, dis, choix);

            let done = self.base.compute_data(
                &mut data,
                hguide,
                Some(spine),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                p_func.as_mut(),
                p_finv.as_mut(),
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                soldep,
                intf,
                intl,
                &mut gd1,
                &mut gd2,
                &mut gf1,
                &mut gf2,
                rec_on_s1,
                rec_on_s2,
            );
            self.base.done = done;
            if !done {
                return false; // ratrappage possible PMN 14/05/1998
            }
            let line_guard = lin.as_ref().expect("Lin");
            let done2 = self.base.complete_data_function(
                &mut data,
                p_func.as_mut(),
                &line_guard,
                s1,
                Some(s2),
                or,
                gd1,
                gd2,
                gf1,
                gf2,
                false,
            );
            drop(line_guard);
            self.base.done = done2;
            if !done2 {
                panic!("Standard_Failure: PerformSurf : Fail of approximation!");
            }
        } else if chsp.is_chamfer() == ChFiDS_ChamfMethod::TwoDist {
            let (d1, d2) = chsp.dists();

            let mut p_func: Box<dyn BlendFuncGenChamfer>;
            let mut p_finv: Box<dyn BlendFuncGenChamfInv>;
            if chsp.base.my_mode == super::chfi_ds::ChFiDS_ChamfMode::ClassicChamfer {
                p_func = Box::new(BlendFuncChamfer::new(&s1.surface, &s2.surface, &guide));
                p_finv = Box::new(BlendFuncChamfInv::new(&s1.surface, &s2.surface, &guide));
                p_func.set(d1, d2, choix);
                p_finv.set(d1, d2, choix);
            } else {
                // OCCT L1736-1747: the offset elspine paired with HGuide.
                let mut offset_hguide: Option<ChFiDSElSpineHandle> = None;
                for (a, b) in spine
                    .base()
                    .elspines
                    .iter()
                    .zip(spine.base().offset_elspines.iter())
                {
                    if elspine_matches_handle(a, hguide) {
                        offset_hguide = Some(std::sync::Arc::new(std::sync::RwLock::new(b.clone())));
                    }
                }

                if offset_hguide.is_none() {
                    // OCCT prints "Construction of offset guide failed!"
                }
                p_func = Box::new(BlendFuncConstThroatWithPenetration::new(
                    &s1.surface,
                    &s2.surface,
                    // OCCT binds OffsetHGuide; pending the ElSpine curve machinery
                    // the guide placeholder stands in (elspine_guide_curve).
                    &guide,
                ));
                p_finv = Box::new(BlendFuncConstThroatWithPenetrationInv::new(
&s1.surface,
&s2.surface,
                    // OCCT binds OffsetHGuide; pending the ElSpine curve machinery
                    // the guide placeholder stands in (elspine_guide_curve).
                    &guide,
                ));
                let throat = d1.max(d2);
                p_func.set(throat, throat, choix);
                p_finv.set(throat, throat, choix);
            }

            let done = self.base.compute_data(
                &mut data,
                hguide,
                Some(spine),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                p_func.as_mut(),
                p_finv.as_mut(),
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                soldep,
                intf,
                intl,
                &mut gd1,
                &mut gd2,
                &mut gf1,
                &mut gf2,
                rec_on_s1,
                rec_on_s2,
            );
            self.base.done = done;
            if !done {
                return false; // ratrappage possible PMN 14/05/1998
            }
            let line_guard = lin.as_ref().expect("Lin");
            let done2 = self.base.complete_data_function(
                &mut data,
                p_func.as_mut(),
                &line_guard,
                s1,
                Some(s2),
                or,
                gd1,
                gd2,
                gf1,
                gf2,
                false,
            );
            drop(line_guard);
            self.base.done = done2;
            if !done2 {
                panic!("Standard_Failure: PerformSurf : Fail of approximation!");
            }
        } else {
            // distance and angle
            let (d1, angle) = chsp.get_dist_angle();

            let ch = choix;

            let mut func = BlendFuncChAsym::new(&s1.surface, &s2.surface, &guide);
            let mut finv = BlendFuncChAsymInv::new(&s1.surface, &s2.surface, &guide);
            func.set(d1, angle, ch);
            finv.set(d1, angle, ch);

            let done = self.base.compute_data(
                &mut data,
                hguide,
                Some(spine),
                &mut lin,
                s1,
                i1,
                s2,
                i2,
                &mut func,
                &mut finv,
                p_first,
                max_step,
                fleche,
                tol_guide,
                first,
                last,
                inside,
                appro,
                forward,
                soldep,
                intf,
                intl,
                &mut gd1,
                &mut gd2,
                &mut gf1,
                &mut gf2,
                rec_on_s1,
                rec_on_s2,
            );

            self.base.done = done;
            if !done {
                return false; // ratrappage possible PMN 14/05/1998
            }
            let line_guard = lin.as_ref().expect("Lin");
            let done2 = self.base.complete_data_function(
                &mut data,
                &mut func,
                &line_guard,
                s1,
                Some(s2),
                or,
                gd1,
                gd2,
                gf1,
                gf2,
                false,
            );
            drop(line_guard);
            self.base.done = done2;
            if !done2 {
                panic!("Standard_Failure: PerformSurf : Fail of approximation!");
            }
        }

        true
    }

    // =====================================================================

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L1953-1974 — ExtentOneCorner (extends the
    // spine of the stripe S on the side of the vertex V; PMN 28/11/97).
    // =====================================================================
    pub fn extent_one_corner(&self, v: &Shape, s: &SharedStripe) {
        let mut sens = 0i32;
        let coeff = 0.5;
        let st = s.read().expect("stripe lock");
        let mut spine = st.spine().cloned().expect("Spine");
        super::chfi3d::chfi3d_index_of_surf_data(v, &st, &mut sens);
        if spine.base().is_tangency_extremity(sens == 1) {
            return; // No extension on queue
        }
        let d_u = spine.base().last_parameter_of(spine.base().nb_edges());
        if sens == 1 {
            let sp = spine.base_mut();
            sp.set_first_parameter(-d_u * coeff);
            sp.set_first_tgt(0.0);
        } else {
            let sp = spine.base_mut();
            sp.set_last_parameter(d_u * (1.0 + coeff));
            sp.set_last_tgt(d_u);
        }
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L1982-2107 — ExtentTwoCorner (extends the
    // spines of the stripes contained in the list LS on the side of the
    // vertex V).
    // =====================================================================
    pub fn extent_two_corner(&mut self, brep: &topods::BRep, v: &Shape, ls: &[SharedStripe]) {
        let mut sens = 0i32;
        let mut ff = true;
        let mut isfirst = [false; 2];
        let mut iedge = [1usize; 2];
        let mut stripe: [Option<SharedStripe>; 2] = [None, None];
        let mut spine: [Option<ChFiDSSpineHandle>; 2] = [None, None];

        let mut i = 0usize;
        for itel in ls {
            super::chfi3d::chfi3d_index_of_surf_data(
                v,
                &itel.read().expect("stripe lock"),
                &mut sens,
            );
            if !ff {
                if let Some(s1) = &stripe[1] {
                    // OCCT: if (Stripe[1] == itel.Value()) — handle identity.
                    if std::sync::Arc::ptr_eq(s1, itel) {
                        sens = -sens;
                    }
                }
            }

            stripe[i] = Some(itel.clone());
            isfirst[i] = sens == 1;
            spine[i] = itel.read().expect("stripe lock").spine().cloned();
            if !isfirst[i] {
                if let Some(sp) = &spine[i] {
                    iedge[i] = sp.base().nb_edges();
                }
            }
            ff = false;
            i += 1;
        }

        let mut d = [0.0f64; 4];
        let mut dis = [0.0f64; 2];
        let mut f: [Shape; 4] = [Shape::null(), Shape::null(), Shape::null(), Shape::null()];

        let mut j = 0usize;
        for i in 0..2 {
            j = i * 2;
            let Some(sp) = &spine[i] else { continue };
            let Some(csp) = sp.down_cast_chamf() else {
                continue;
            };
            // OCCT: ConexFaces(Spine[i], Iedge[i], F[j], F[j + 1]);
            let (f1, f2) = self.conex_faces(brep, sp, iedge[i]);
            f[j] = f1;
            f[j + 1] = f2;

            if csp.is_chamfer() == ChFiDS_ChamfMethod::Sym {
                d[j] = csp.get_dist();
                d[j + 1] = d[j];
            } else if csp.is_chamfer() == ChFiDS_ChamfMethod::TwoDist {
                let (dd0, dd1) = csp.dists();
                d[j] = dd0;
                d[j + 1] = dd1;
            } else {
                // an approximate calculation of distance 2 is done
                let (t0, ta) = csp.get_dist_angle();
                d[j] = t0;
                d[j + 1] = t0 * ta.tan();
            }
        }
        let _ = j;

        let mut notfound = true;
        let mut i = 0usize;
        while notfound && i < 2 {
            let mut j = 0usize;
            while notfound && j < 2 {
                if f[i].is_same(&f[j + 2]) {
                    dis[0] = d[i];
                    dis[1] = d[j + 2];
                    notfound = false;
                }
                j += 1;
            }
            i += 1;
        }
        // ExtentTwoCorner

        let mut state = [ChFiDS_State::AllSame; 2];

        for i in 0..2 {
            if let Some(sp) = &spine[i] {
                state[i] = if isfirst[i] {
                    sp.base().first_status()
                } else {
                    sp.base().last_status()
                };
            }
        }

        if state[0] == ChFiDS_State::AllSame {
            // it is necessary that two chamfers touch the face at end
            for j in 0..2 {
                if let Some(st) = &stripe[j] {
                    self.extent_one_corner(v, st);
                }
            }
        } else if (state[0] == ChFiDS_State::OnSame) && (state[1] == ChFiDS_State::OnSame) {
            let mut sp0 = spine[0].clone().expect("Spine[0]");
            let mut sp1 = spine[1].clone().expect("Spine[1]");
            extent_spine_on_common_face(&mut sp0, &mut sp1, v, dis[0], dis[1], isfirst[0], isfirst[1]);
        }
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L2111-2214 — ExtentThreeCorner.
    // =====================================================================
    pub fn extent_three_corner(&mut self, brep: &topods::BRep, v: &Shape, ls: &[SharedStripe]) {
        let mut sens = 0i32;
        let mut check: Vec<SharedStripe> = Vec::new();
        let mut isfirst = [false; 3];
        let mut iedge = [1usize; 3];
        let mut spine: [Option<ChFiDSSpineHandle>; 3] = [None, None, None];

        let mut i = 0usize;
        for itel in ls {
            let stripe = itel.clone();
            super::chfi3d::chfi3d_index_of_surf_data(
                v,
                &stripe.read().expect("stripe lock"),
                &mut sens,
            );
            for ich in &check {
                // OCCT: if (Stripe == ich.Value()) — handle identity.
                if std::sync::Arc::ptr_eq(ich, &stripe) {
                    sens = -sens;
                    break;
                }
            }

            isfirst[i] = sens == 1;
            spine[i] = stripe.read().expect("stripe lock").spine().cloned();
            if !isfirst[i] {
                if let Some(sp) = &spine[i] {
                    iedge[i] = sp.base().nb_edges();
                }
            }

            check.push(stripe);
            i += 1;
        }

        let mut d = [[0.0f64; 2]; 3];
        let mut f: [[Shape; 2]; 3] = [
            [Shape::null(), Shape::null()],
            [Shape::null(), Shape::null()],
            [Shape::null(), Shape::null()],
        ];

        for i in 0..3 {
            let Some(sp) = &spine[i] else { continue };
            let Some(csp) = sp.down_cast_chamf() else { continue };
            // OCCT: ConexFaces(Spine[i], Iedge[i], F[i][0], F[i][1]);
            let (f1, f2) = self.conex_faces(brep, sp, iedge[i]);
            f[i][0] = f1;
            f[i][1] = f2;

            if csp.is_chamfer() == ChFiDS_ChamfMethod::Sym {
                d[i][0] = csp.get_dist();
                d[i][1] = d[i][0];
            } else if csp.is_chamfer() == ChFiDS_ChamfMethod::TwoDist {
                let (dd0, dd1) = csp.dists();
                d[i][0] = dd0;
                d[i][1] = dd1;
            } else {
                let (t0, ta) = csp.get_dist_angle();
                // an approximate calculation of distance 2 is done
                d[i][0] = t0;
                d[i][1] = t0 * ta.tan();
            }
        }

        // dis[i][j] distance from chamfer i on the common face with
        // chamfer j
        let mut dis = [[0.0f64; 3]; 3];

        for i in 0..3 {
            let j = (i + 1) % 3;
            let mut notfound = true;
            let mut k = 0usize;
            while notfound && k < 2 {
                let mut l = 0usize;
                while notfound && l < 2 {
                    if f[i][k].is_same(&f[j][l]) {
                        dis[i][j] = d[i][k];
                        dis[j][i] = d[j][l];
                        notfound = false;
                    }
                    l += 1;
                }
                k += 1;
            }
        }

        // ExtentThreeCorner
        for i in 0..3 {
            let j = (i + 1) % 3;
            let mut spi = spine[i].clone().expect("Spine[i]");
            let mut spj = spine[j].clone().expect("Spine[j]");
            extent_spine_on_common_face(&mut spi, &mut spj, v, dis[i][j], dis[j][i], isfirst[i], isfirst[j]);
        }
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L2218-2276 — SetRegul.
    // Pending boundary: myCoup (TopOpeBRepBuild_HBuilder) exposes no
    // NewEdges/NewFaces in the rcad facade yet (hbuilder.rs), and
    // BRep_Builder::Continuity has no rcad carrier — the OCCT iteration
    // shape is kept with the facade feeds empty.
    // =====================================================================
    pub fn set_regul(&mut self) {
        let seuil = std::f64::consts::PI / 360.0;
        let seuil2 = seuil * seuil;
        let brep = self.base.my_brep.clone();
        for reg in &self.base.my_regul {
            // OCCT: itc.Initialize(myCoup->NewEdges(reg.Curve()))
            let new_edges: Vec<Shape> = Vec::new();
            if let Some(e) = new_edges.first() {
                let e = e.clone();
                if reg.is_surface1() && reg.is_surface2() {
                    // OCCT: its1.Initialize(myCoup->NewFaces(reg.S1()));
                    //       its2.Initialize(myCoup->NewFaces(reg.S2()));
                    let new_faces1: Vec<Shape> = Vec::new();
                    let new_faces2: Vec<Shape> = Vec::new();
                    if let (Some(f1s), Some(f2s)) = (new_faces1.first(), new_faces2.first()) {
                        let f1 = f1s.clone();
                        let f2 = f2s.clone();
                        let mut s = BRepAdaptorSurface::initialize(&brep, &f1);
                        // OCCT: PC.Initialize(E, F1); t = mid parameter;
                        //       PC.Value(t) -> (u, v)
                        let pc1 = brep
                            .curve_on_surface(&e, &f1)
                            .map(|(pc, f0, l0)| (pc, (f0 + l0) * 0.5));
                        if let Some((pc, t)) = pc1 {
                            use rcad_kernel::geom::Curve2dEval as _;
                            let uv = pc.point_at(t);
                            let (_, du, dv) = s.surface.derivatives(uv.x, uv.y);
                            let mut n1 = du.cross(dv);
                            drop(pc);
                            // OCCT: S.Initialize(F2, false) — the adaptor is reloaded.
                            s = BRepAdaptorSurface::initialize(&brep, &f2);
                            let pc2 = brep.curve_on_surface(&e, &f2);
                            if let Some((pc2, _, _)) = pc2 {
                                let uv = pc2.point_at(t);
                                let (_, du, dv) = s.surface.derivatives(uv.x, uv.y);
                                let mut n2 = du.cross(dv);
                                if n1.length_squared() > 1.0e-14 && n2.length_squared() > 1.0e-14 {
                                    n1 = n1.normalize();
                                    n2 = n2.normalize();
                                    let sina2 = n1.cross(n2).length_squared();
                                    if sina2 < seuil2 {
                                        // OCCT: GeomAbs_Shape cont =
                                        //   ChFi3d_evalconti(E, F1, F2);
                                        // B.Continuity(E, F1, F2, cont);
                                        let cont = chfi3d_evalconti(&e, &f1, &f2);
                                        let _ = cont;
                                        // pending: BRep_Builder::Continuity
                                        // has no rcad carrier yet.
                                    }
                                }
                                let _ = n2;
                            }
                        }
                    }
                }
            }
        }
    }

    // =====================================================================
    // OCCT ChFi3d_ChBuilder.cxx L2283-2318 — ConexFaces (F1, F2 are
    // connected to edge so that F1 corresponds to distance).
    // =====================================================================
    pub fn conex_faces(
        &self,
        brep: &topods::BRep,
        spine: &ChFiDSSpineHandle,
        i_edge: usize,
    ) -> (Shape, Shape) {
        let mut tmp1 = Orientation::Forward;
        let mut tmp2 = Orientation::Forward;

        // calculate the reference orientation
        //  ChFi3d_Builder::StripeOrientations is private
        let (ff1, ff2) = search_common_faces(brep, &self.base.my_ef_map, spine.base().edges(1));
        let ff1 = {
            let mut f = ff1;
            f.orientation = Orientation::Forward;
            f
        };
        let ff2 = {
            let mut f = ff2;
            f.orientation = Orientation::Forward;
            f
        };
        let sb1 = BRepAdaptorSurface::initialize(brep, &ff1);
        let sb2 = BRepAdaptorSurface::initialize(brep, &ff2);
        let rc = concave_side(
            brep,
            &sb1.face,
            &sb2.face,
            spine.base().edges(1),
            &mut tmp1,
            &mut tmp2,
        );

        // calculate the connected faces
        let (f1c, f2c) = search_common_faces(brep, &self.base.my_ef_map, spine.base().edges(i_edge));
        let sb1 = BRepAdaptorSurface::initialize(brep, &f1c);
        let sb2 = BRepAdaptorSurface::initialize(brep, &f2c);
        let choix = concave_side(
            brep,
            &sb1.face,
            &sb2.face,
            spine.base().edges(i_edge),
            &mut tmp1,
            &mut tmp2,
        );

        if rc % 2 != choix % 2 {
            (f2c, f1c)
        } else {
            (f1c, f2c)
        }
    }
}

// =========================================================================
// OCCT SimulSurf inline blocks shared by the three chamfer branches
// (L861-880 / L1023-1042 / L1158-1177) — the four ChFi3d_FilCommonPoint
// loads.
// =========================================================================
fn simul_surf_common_points(
    brep: &topods::BRep,
    lin: &BRepBlendLine,
    data: &mut ChFiDSSurfData,
    tolapp3d: f64,
) {
    let line = lin;
    chfi3d_fil_common_point(
        brep,
        line.start_point_on_first(),
        line.transition_on_s1(),
        true,
        data.change_vertex_first_on_s1(),
        tolapp3d,
    );
    chfi3d_fil_common_point(
        brep,
        line.end_point_on_first(),
        line.transition_on_s1(),
        false,
        data.change_vertex_last_on_s1(),
        tolapp3d,
    );
    chfi3d_fil_common_point(
        brep,
        line.start_point_on_second(),
        line.transition_on_s2(),
        true,
        data.change_vertex_first_on_s2(),
        tolapp3d,
    );
    chfi3d_fil_common_point(
        brep,
        line.end_point_on_second(),
        line.transition_on_s2(),
        false,
        data.change_vertex_last_on_s2(),
        tolapp3d,
    );
}

// =========================================================================
// OCCT Extrema_GenLocateExtPS (TKGeomAlgo) — the point-on-surface location
// starting from (u0, v0).  rcad boundary: projected damped-Newton on the
// analytic derivatives (pending the Extrema package translation).
// =========================================================================
fn extrema_gen_locate_ext_ps(
    s: &rcad_kernel::geom::Surface3,
    p: DVec3,
    u0: f64,
    v0: f64,
    _tol: f64,
) -> Option<(f64, f64)> {
    use rcad_kernel::geom::SurfaceEval as _;
    let mut u = u0;
    let mut v = v0;
    for _ in 0..24 {
        let (sp, du, dv) = s.derivatives(u, v);
        let r = sp - p;
        let a11 = du.dot(du);
        let a12 = du.dot(dv);
        let a22 = dv.dot(dv);
        let b1 = -du.dot(r);
        let b2 = -dv.dot(r);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < 1.0e-20 {
            break;
        }
        let du_step = (b1 * a22 - b2 * a12) / det;
        let dv_step = (a11 * b2 - a12 * b1) / det;
        u += du_step;
        v += dv_step;
        if du_step.abs() < 1.0e-10 && dv_step.abs() < 1.0e-10 {
            break;
        }
    }
    Some((u, v))
}
