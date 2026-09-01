//! OCCT Intf package (Intf_PIType / Intf_SectionPoint) plus the
//! IntCurveSurface_InterUtils::SectionPointToParameters template.
//!
//! 1:1 translations:
//!   - Intf_PIType.hxx — the intersection point dimension enum (numeric
//!     order matters: EXTERNAL=0 < FACE=1 < EDGE=2 < VERTEX=3, used by
//!     Intf_SectionPoint::Merge).
//!   - Intf_SectionPoint.cxx L24-256 — the section point data class.
//!   - IntCurveSurface_InterUtils.pxx L740-848 — SectionPointToParameters:
//!     barycentric parameter interpolation on the polyhedron face with the
//!     degenerate-face fallback to the longest edge, plus the polygon's
//!     curve parameter.

use glam::{DVec2, DVec3};
use rcad_kernel::base::int_ana2d::{AnaIntersection2d, Conic2d};
use rcad_kernel::core::precision::{ANGULAR, COMPUTATIONAL, INFINITE_VALUE, SQUARE_CONFUSION};
use rcad_kernel::geom::{
    Curve2dEval, CurveEval, Hyperbola2d, Hyperbola3, Line2d, Line3, Parabola2d, Parabola3,
};
use rcad_kernel::math::bnd::{BndBox, BndBox2d};

/// OCCT Intf_PIType (Intf_PIType.hxx) — the dimension of the intersection
/// support on an object.  The numeric order (EXTERNAL < FACE < EDGE < VERTEX)
/// is significant (Intf_SectionPoint::Merge).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum IntfPIType {
    External,
    Face,
    Edge,
    Vertex,
}

/// OCCT Intf_SectionPoint (Intf_SectionPoint.hxx / .cxx) — an intersection
/// point with its two supports (on the object and on the tool).
#[derive(Debug, Clone, Copy)]
pub struct IntfSectionPoint {
    my_pnt: DVec3,
    pub dimen_obje: IntfPIType,
    pub index_o1: i32,
    pub index_o2: i32,
    pub param_obje: f64,
    pub dimen_tool: IntfPIType,
    pub index_t1: i32,
    pub index_t2: i32,
    pub param_tool: f64,
    pub incide: f64,
}

impl Default for IntfSectionPoint {
    /// OCCT Intf_SectionPoint() (Intf_SectionPoint.cxx L137-149).
    fn default() -> Self {
        IntfSectionPoint {
            my_pnt: DVec3::ZERO,
            dimen_obje: IntfPIType::External,
            index_o1: 0,
            index_o2: 0,
            param_obje: 0.0,
            dimen_tool: IntfPIType::External,
            index_t1: 0,
            index_t2: 0,
            param_tool: 0.0,
            incide: 0.0,
        }
    }
}

impl IntfSectionPoint {
    /// OCCT Intf_SectionPoint(Where, Dim1, Addr1, Addr2, Param1, Dim2,
    /// Addr3, Addr4, Param2, Incid) (L153-174).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        where_: DVec3,
        dim1: IntfPIType,
        addr1: i32,
        addr2: i32,
        param1: f64,
        dim2: IntfPIType,
        addr3: i32,
        addr4: i32,
        param2: f64,
        incid: f64,
    ) -> Self {
        IntfSectionPoint {
            my_pnt: where_,
            dimen_obje: dim1,
            index_o1: addr1,
            index_o2: addr2,
            param_obje: param1,
            dimen_tool: dim2,
            index_t1: addr3,
            index_t2: addr4,
            param_tool: param2,
            incide: incid,
        }
    }

    /// OCCT Intf_SectionPoint::Pnt() (L24-27).
    pub fn pnt(&self) -> DVec3 {
        self.my_pnt
    }

    /// OCCT Intf_SectionPoint::InfoFirst(Dim, Add1, Add2, Param) (L31-37).
    pub fn info_first_full(&self) -> (IntfPIType, i32, i32, f64) {
        (self.dimen_obje, self.index_o1, self.index_o2, self.param_obje)
    }

    /// OCCT Intf_SectionPoint::InfoFirst(Dim, Add, Param) (L41-46).
    pub fn info_first(&self) -> (IntfPIType, i32, f64) {
        (self.dimen_obje, self.index_o2, self.param_obje)
    }

    /// OCCT Intf_SectionPoint::InfoSecond(Dim, Add1, Add2, Param) (L50-56).
    pub fn info_second_full(&self) -> (IntfPIType, i32, i32, f64) {
        (self.dimen_tool, self.index_t1, self.index_t2, self.param_tool)
    }

    /// OCCT Intf_SectionPoint::InfoSecond(Dim, Add, Param) (L60-65).
    pub fn info_second(&self) -> (IntfPIType, i32, f64) {
        (self.dimen_tool, self.index_t2, self.param_tool)
    }

    /// OCCT Intf_SectionPoint::Incidence() (L69-72).
    pub fn incidence(&self) -> f64 {
        self.incide
    }

    /// OCCT Intf_SectionPoint::IsOnSameEdge (L76-133).
    pub fn is_on_same_edge(&self, other: &IntfSectionPoint) -> bool {
        let mut is_on = false;
        if self.dimen_obje == IntfPIType::Edge {
            if other.dimen_obje == IntfPIType::Edge {
                is_on = self.index_o1 == other.index_o1 && self.index_o2 == other.index_o2;
            } else if other.dimen_obje == IntfPIType::Vertex {
                is_on = self.index_o1 == other.index_o1 || self.index_o2 == other.index_o1;
            }
        } else if self.dimen_obje == IntfPIType::Vertex {
            if other.dimen_obje == IntfPIType::Edge {
                is_on = self.index_o1 == other.index_o1 || self.index_o1 == other.index_o2;
            } else if other.dimen_obje == IntfPIType::Vertex {
                is_on = self.index_t1 == other.index_t1;
            }
        }
        if !is_on {
            if self.dimen_tool == IntfPIType::Edge {
                if other.dimen_tool == IntfPIType::Edge {
                    is_on = self.index_t1 == other.index_t1 && self.index_t2 == other.index_t2;
                } else if other.dimen_tool == IntfPIType::Vertex {
                    is_on = self.index_t1 == other.index_t1 || self.index_t2 == other.index_t1;
                }
            } else if self.dimen_tool == IntfPIType::Vertex {
                if other.dimen_tool == IntfPIType::Edge {
                    is_on = self.index_t1 == other.index_t1 || self.index_t1 == other.index_t2;
                } else if other.dimen_tool == IntfPIType::Vertex {
                    is_on = self.index_t1 == other.index_t1;
                }
            }
        }
        is_on
    }

    /// OCCT Intf_SectionPoint::Merge (L201-232) — merge the higher-dimension
    /// supports into `other`.
    pub fn merge(&mut self, other: &mut IntfSectionPoint) {
        other.my_pnt = self.my_pnt;
        if self.dimen_obje >= other.dimen_obje {
            other.dimen_obje = self.dimen_obje;
            other.index_o1 = self.index_o1;
            other.index_o2 = self.index_o2;
            other.param_obje = self.param_obje;
        } else {
            self.dimen_obje = other.dimen_obje;
            self.index_o1 = other.index_o1;
            self.index_o2 = other.index_o2;
            self.param_obje = other.param_obje;
        }
        if self.dimen_tool >= other.dimen_tool {
            other.dimen_tool = self.dimen_tool;
            other.index_t1 = self.index_t1;
            other.index_t2 = self.index_t2;
            other.param_tool = self.param_tool;
        } else {
            self.dimen_tool = other.dimen_tool;
            self.index_t1 = other.index_t1;
            self.index_t2 = other.index_t2;
            self.param_tool = other.param_tool;
        }
    }
}

/// The polyhedron interface of IntCurveSurface_InterUtils::SectionPointToParameters
/// (the OCCT template parameter PolyhedronType): a triangulated sampling of a
/// surface.  Indices are 1-based (the OCCT polyhedron convention).
pub trait PolyhedronLike {
    fn triangle(&self, t: usize) -> (usize, usize, usize);
    fn point(&self, index: usize) -> DVec3;
    fn parameters(&self, index: usize) -> (f64, f64);
}

/// The polygon interface (the OCCT template parameter PolygonType): the
/// sampled curve with the approximate parameter mapping.
pub trait PolygonLike {
    fn approx_param_on_curve(&self, index: usize, param_on_line: f64) -> f64;
}

/// OCCT IntCurveSurface_InterUtils::SectionPointToParameters
/// (IntCurveSurface_InterUtils.pxx L740-848) — converts a section point's
/// supports into the approximate surface (U, V) and curve (W) parameters.
/// Returns (U, V, W).
pub fn section_point_to_parameters<P: PolyhedronLike, G: PolygonLike>(
    sp: &IntfSectionPoint,
    polyhedron: &P,
    polygon: &G,
) -> (f64, f64, f64) {
    let (_typ, adr1, adr2, param) = sp.info_second_full();
    let p = sp.pnt();
    let mut u1 = 0.0;
    let mut v1 = 0.0;
    //----------------------------------------------------------------------
    //--          Approximate parameter calculation on surface            --
    //----------------------------------------------------------------------
    match _typ {
        //-- Adr1 is the vertex number
        IntfPIType::Vertex => {
            let (a, b) = polyhedron.parameters(adr1 as usize);
            u1 = a;
            v1 = b;
        }
        IntfPIType::Edge => {
            let (a, b) = polyhedron.parameters(adr1 as usize);
            u1 = a;
            v1 = b;
            let (u, v) = polyhedron.parameters(adr2 as usize);
            u1 += param * (u - u1);
            v1 += param * (v - v1);
        }
        IntfPIType::Face => {
            let (pt1, pt2, pt3) = polyhedron.triangle(adr1 as usize);
            let pa = polyhedron.point(pt1);
            let pb = polyhedron.point(pt2);
            let pc = polyhedron.point(pt3);
            let (ua, va) = polyhedron.parameters(pt1);
            let (ub, vb) = polyhedron.parameters(pt2);
            let (uc, vc) = polyhedron.parameters(pt3);
            let normale = (pb - pa).cross(pc - pa);
            let cc = (pb - pa).cross(p - pa).dot(normale);
            let ca = (pc - pb).cross(p - pb).dot(normale);
            let cb = (pc - pa).cross(p - pc).dot(normale);
            let cabc = ca + cb + cc;

            if cabc.abs() > COMPUTATIONAL {
                let ca = ca / cabc;
                let cb = cb / cabc;
                let cc = cc / cabc;
                u1 = ca * ua + cb * ub + cc * uc;
                v1 = ca * va + cb * vb + cc * vc;
            } else {
                let a_ab = pb - pa;
                let a_bc = pc - pb;
                let a_ca = pa - pc;
                let sq_ab = a_ab.length_squared();
                let sq_bc = a_bc.length_squared();
                let sq_ca = a_ca.length_squared();

                if sq_ab >= sq_bc && sq_ab >= sq_ca && sq_ab > SQUARE_CONFUSION {
                    let mut a_t = (p - pa).dot(a_ab) / sq_ab;
                    a_t = a_t.min(1.0).max(0.0);
                    u1 = ua + a_t * (ub - ua);
                    v1 = va + a_t * (vb - va);
                } else if sq_bc >= sq_ca && sq_bc > SQUARE_CONFUSION {
                    let mut a_t = (p - pb).dot(a_bc) / sq_bc;
                    a_t = a_t.min(1.0).max(0.0);
                    u1 = ub + a_t * (uc - ub);
                    v1 = vb + a_t * (vc - vb);
                } else if sq_ca > SQUARE_CONFUSION {
                    let mut a_t = (p - pc).dot(a_ca) / sq_ca;
                    a_t = a_t.min(1.0).max(0.0);
                    u1 = uc + a_t * (ua - uc);
                    v1 = vc + a_t * (va - vc);
                } else {
                    u1 = ua;
                    v1 = va;
                }
            }
        }
        IntfPIType::External => {}
    }
    //----------------------------------------------------------------------
    //--              Approximate point calculation on Curve              --
    //----------------------------------------------------------------------
    let (_typ, seg_index, param) = sp.info_first();
    let w = polygon.approx_param_on_curve(seg_index as usize, param);
    (u1, v1, w)
}

// =============================================================================
// Intf_Tool — box computation for infinite conics (Intf_Tool.hxx/.cxx)
// =============================================================================

/// 2D cross product (OCCT gp_XY::Crossed).
fn cross2d(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// OCCT ElCLib::Value(para, gp_Hypr2d) — hyperbola point (Curve2dEval).
fn hypr2d_value(h: &Hyperbola2d, t: f64) -> DVec2 {
    h.point_at(t)
}

/// OCCT ElCLib::D1(para, gp_Hypr2d, P, V).
fn hypr2d_d1(h: &Hyperbola2d, t: f64) -> (DVec2, DVec2) {
    (h.point_at(t), h.derivative_at(t))
}

/// OCCT ElCLib::Value(para, gp_Parab2d).
fn parab2d_value(p: &Parabola2d, t: f64) -> DVec2 {
    p.point_at(t)
}

/// OCCT ElCLib::D1(para, gp_Parab2d, P, V).
fn parab2d_d1(p: &Parabola2d, t: f64) -> (DVec2, DVec2) {
    (p.point_at(t), p.derivative_at(t))
}

/// OCCT gp_Pln(A, B, C, D) — the plane A·x + B·y + C·z + D = 0.
fn plane_from_coeffs(a: f64, b: f64, c: f64, d: f64) -> rcad_kernel::geom::Plane {
    let n = DVec3::new(a, b, c);
    let n2 = n.length_squared();
    let origin = if n2 > 0.0 { -d * n / n2 } else { DVec3::ZERO };
    rcad_kernel::geom::Plane::new(origin, n)
}

/// OCCT Standard_Real.hxx IsEqual — |a - b| < RealSmall() (DBL_MIN).
fn is_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < f64::MIN_POSITIVE
}

/// OCCT Intf_Tool — creates boxes for infinite lines/conics in a given
/// bounding domain and exposes the parameter segments of the curve portions
/// inside the domain (Intf_Tool.hxx L35-85, Intf_Tool.cxx whole).
pub struct IntfTool {
    nb_seg: usize,
    begin_on_curve: [f64; 6],
    end_on_curve: [f64; 6],
    bord: [i32; 12],
    xint: [f64; 12],
    yint: [f64; 12],
    zint: [f64; 12],
    parint: [f64; 12],
}

impl IntfTool {
    /// OCCT Intf_Tool() (Intf_Tool.cxx L38-48).
    pub fn new() -> Self {
        IntfTool {
            nb_seg: 0,
            begin_on_curve: [0.0; 6],
            end_on_curve: [0.0; 6],
            bord: [0; 12],
            xint: [0.0; 12],
            yint: [0.0; 12],
            zint: [0.0; 12],
            parint: [0.0; 12],
        }
    }

    /// OCCT Intf_Tool::Lin2dBox (Intf_Tool.cxx L52-206).
    pub fn lin_2d_box(&mut self, l2d: &Line2d, domain: &BndBox2d, box_lin: &mut BndBox2d) {
        self.nb_seg = 0;
        box_lin.set_void();
        if domain.is_whole() {
            let loc = l2d.origin;
            let dir = l2d.direction.normalize_or_zero();
            box_lin.add_point(loc);
            box_lin.add_point(loc + dir);
            box_lin.add_point(loc - dir);
            self.nb_seg = 1;
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            return;
        } else if domain.is_void() {
            return;
        }

        let (xmin, ymin, xmax, ymax) = domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0));
        let mut xmin_b = 0.0;
        let mut xmax_b = 0.0;
        let mut ymin_b = 0.0;
        let mut ymax_b = 0.0;
        let mut parmin = -INFINITE_VALUE;
        let mut parmax = INFINITE_VALUE;
        let (mut x_to_set, mut y_to_set) = (false, false);

        if l2d.direction.x > 0.0 {
            parmin = if domain.is_open_xmin() {
                -INFINITE_VALUE
            } else {
                (xmin - l2d.origin.x) / l2d.direction.x
            };
            parmax = if domain.is_open_xmax() {
                INFINITE_VALUE
            } else {
                (xmax - l2d.origin.x) / l2d.direction.x
            };
            x_to_set = true;
        } else if l2d.direction.x < 0.0 {
            parmin = if domain.is_open_xmax() {
                -INFINITE_VALUE
            } else {
                (xmax - l2d.origin.x) / l2d.direction.x
            };
            parmax = if domain.is_open_xmin() {
                INFINITE_VALUE
            } else {
                (xmin - l2d.origin.x) / l2d.direction.x
            };
            x_to_set = true;
        } else if l2d.origin.x < xmin || xmax < l2d.origin.x {
            return;
        } else {
            // Parallel to the X axis.
            xmin_b = l2d.origin.x;
            xmax_b = l2d.origin.x;
            x_to_set = false;
        }

        if l2d.direction.y > 0.0 {
            let parcur = if domain.is_open_ymin() {
                -INFINITE_VALUE
            } else {
                (ymin - l2d.origin.y) / l2d.direction.y
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymax() {
                INFINITE_VALUE
            } else {
                (ymax - l2d.origin.y) / l2d.direction.y
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else if l2d.direction.y < 0.0 {
            let parcur = if domain.is_open_ymax() {
                -INFINITE_VALUE
            } else {
                (ymax - l2d.origin.y) / l2d.direction.y
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymin() {
                INFINITE_VALUE
            } else {
                (ymin - l2d.origin.y) / l2d.direction.y
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else if l2d.origin.y < ymin || ymax < l2d.origin.y {
            return;
        } else {
            // Parallel to the Y axis.
            ymin_b = l2d.origin.y;
            ymax_b = l2d.origin.y;
            y_to_set = false;
        }

        self.nb_seg += 1;
        self.begin_on_curve[0] = parmin;
        self.end_on_curve[0] = parmax;

        if x_to_set {
            let par1 = l2d.origin.x + parmin * l2d.direction.x;
            let par2 = l2d.origin.x + parmax * l2d.direction.x;
            xmin_b = par1.min(par2);
            xmax_b = par1.max(par2);
        }
        if y_to_set {
            let par1 = l2d.origin.y + parmin * l2d.direction.y;
            let par2 = l2d.origin.y + parmax * l2d.direction.y;
            ymin_b = par1.min(par2);
            ymax_b = par1.max(par2);
        }

        box_lin.update(xmin_b, ymin_b, xmax_b, ymax_b);
    }

    /// OCCT Intf_Tool::Hypr2dBox (Intf_Tool.cxx L210-369).
    pub fn hypr_2d_box(&mut self, hypr: &Hyperbola2d, domain: &BndBox2d, box_hypr: &mut BndBox2d) {
        self.nb_seg = 0;
        box_hypr.set_void();
        if domain.is_whole() {
            box_hypr.set_whole();
            self.nb_seg = 1;
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            return;
        } else if domain.is_void() {
            return;
        }

        let nb_pi = self.inters_2d_hypr(hypr, domain);

        if nb_pi > 0 {
            let (mut xmin_b, mut ymin_b, mut xmax_b, mut ymax_b) =
                domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0));

            for npi in 0..nb_pi {
                xmin_b = xmin_b.min(self.xint[npi]);
                xmax_b = xmax_b.max(self.xint[npi]);
                ymin_b = ymin_b.min(self.yint[npi]);
                ymax_b = ymax_b.max(self.yint[npi]);
            }
            box_hypr.update(xmin_b, ymin_b, xmax_b, ymax_b);

            // Selection sort of the intersection parameters.
            for npi in 0..nb_pi {
                let mut npk = npi;
                for npj in npi + 1..nb_pi {
                    if self.parint[npj] < self.parint[npk] {
                        npk = npj;
                    }
                }
                if npk != npi {
                    let parmin = self.parint[npk];
                    self.parint[npk] = self.parint[npi];
                    self.parint[npi] = parmin;
                    let npj = self.bord[npk];
                    self.bord[npk] = self.bord[npi];
                    self.bord[npi] = npj;
                }
            }

            let mut sinan = 0.0;
            let mut out = true;

            for npi in 0..nb_pi {
                let (_, tan) = hypr2d_d1(hypr, self.parint[npi]);
                sinan = match self.bord[npi] {
                    1 => cross2d(DVec2::new(-1.0, 0.0), tan),
                    2 => cross2d(DVec2::new(0.0, -1.0), tan),
                    3 => cross2d(DVec2::new(1.0, 0.0), tan),
                    4 => cross2d(DVec2::new(0.0, 1.0), tan),
                    _ => 0.0,
                };
                if sinan.abs() > ANGULAR {
                    if sinan > 0.0 {
                        if self.nb_seg < 6 {
                            out = false;
                            self.begin_on_curve[self.nb_seg] = self.parint[npi];
                            self.nb_seg += 1;
                        }
                    } else {
                        if out && self.nb_seg < 6 {
                            self.begin_on_curve[self.nb_seg] = -INFINITE_VALUE;
                            self.nb_seg += 1;
                        }
                        if self.nb_seg > 0 {
                            self.end_on_curve[self.nb_seg - 1] = self.parint[npi];
                        }
                        out = true;

                        let mut ipmin = if self.begin_on_curve[self.nb_seg - 1] < -10.0 {
                            -10
                        } else {
                            self.begin_on_curve[self.nb_seg - 1] as i32
                        };
                        let mut ipmax = if self.end_on_curve[self.nb_seg - 1] > 10.0 {
                            10
                        } else {
                            self.end_on_curve[self.nb_seg - 1] as i32
                        };
                        ipmin = ipmin * 10 + 1;
                        ipmax = ipmax * 10 - 1;
                        let mut ip = ipmin;
                        let mut pas = 1;
                        while ip <= ipmax {
                            box_hypr.add_point(hypr2d_value(hypr, ip as f64 / 10.0));
                            if ip.abs() <= 10 {
                                pas = 1;
                            } else {
                                pas = 10;
                            }
                            ip += pas;
                        }
                    }
                }
            }
            if !out && self.nb_seg > 0 {
                self.end_on_curve[self.nb_seg - 1] = INFINITE_VALUE;
            }
        } else if !domain.is_out_point(hypr2d_value(hypr, 0.0)) {
            *box_hypr = domain.clone();
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            self.nb_seg = 1;
        }
    }

    /// OCCT Intf_Tool::Parab2dBox (Intf_Tool.cxx L477-638).
    pub fn parab_2d_box(&mut self, parab: &Parabola2d, domain: &BndBox2d, box_parab: &mut BndBox2d) {
        self.nb_seg = 0;
        box_parab.set_void();
        if domain.is_whole() {
            box_parab.set_whole();
            self.nb_seg = 1;
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            return;
        } else if domain.is_void() {
            return;
        }

        let nb_pi = self.inters_2d_parab(parab, domain);

        if nb_pi > 0 {
            let (mut xmin_b, mut ymin_b, mut xmax_b, mut ymax_b) =
                domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0));

            for npi in 0..nb_pi {
                xmin_b = xmin_b.min(self.xint[npi]);
                xmax_b = xmax_b.max(self.xint[npi]);
                ymin_b = ymin_b.min(self.yint[npi]);
                ymax_b = ymax_b.max(self.yint[npi]);
            }
            box_parab.update(xmin_b, ymin_b, xmax_b, ymax_b);

            for npi in 0..nb_pi {
                let mut npk = npi;
                for npj in npi + 1..nb_pi {
                    if self.parint[npj] < self.parint[npk] {
                        npk = npj;
                    }
                }
                if npk != npi {
                    let parmin = self.parint[npk];
                    self.parint[npk] = self.parint[npi];
                    self.parint[npi] = parmin;
                    let npj = self.bord[npk];
                    self.bord[npk] = self.bord[npi];
                    self.bord[npi] = npj;
                }
            }

            let mut sinan = 0.0;
            let mut out = true;

            for npi in 0..nb_pi {
                let (_, tan) = parab2d_d1(parab, self.parint[npi]);
                sinan = match self.bord[npi] {
                    1 => cross2d(DVec2::new(-1.0, 0.0), tan),
                    2 => cross2d(DVec2::new(0.0, -1.0), tan),
                    3 => cross2d(DVec2::new(1.0, 0.0), tan),
                    4 => cross2d(DVec2::new(0.0, 1.0), tan),
                    _ => 0.0,
                };
                if sinan.abs() > ANGULAR {
                    if sinan > 0.0 {
                        if self.nb_seg < 6 {
                            out = false;
                            self.begin_on_curve[self.nb_seg] = self.parint[npi];
                            self.nb_seg += 1;
                        }
                    } else {
                        if out && self.nb_seg < 6 {
                            self.begin_on_curve[self.nb_seg] = -INFINITE_VALUE;
                            self.nb_seg += 1;
                        }
                        if self.nb_seg > 0 {
                            self.end_on_curve[self.nb_seg - 1] = self.parint[npi];
                        }
                        out = true;

                        let mut ipmin = if self.begin_on_curve[self.nb_seg - 1] < -10.0 {
                            -10
                        } else {
                            self.begin_on_curve[self.nb_seg - 1] as i32
                        };
                        let mut ipmax = if self.end_on_curve[self.nb_seg - 1] > 10.0 {
                            10
                        } else {
                            self.end_on_curve[self.nb_seg - 1] as i32
                        };
                        ipmin = ipmin * 10 + 1;
                        ipmax = ipmax * 10 - 1;
                        let mut ip = ipmin;
                        let mut pas = 1;
                        while ip <= ipmax {
                            box_parab.add_point(parab2d_value(parab, ip as f64 / 10.0));
                            if ip.abs() <= 10 {
                                pas = 1;
                            } else {
                                pas = 10;
                            }
                            ip += pas;
                        }
                    }
                }
            }
            if !out && self.nb_seg > 0 {
                self.end_on_curve[self.nb_seg - 1] = INFINITE_VALUE;
            }
        } else if !domain.is_out_point(parab2d_value(parab, 0.0)) {
            *box_parab = domain.clone();
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            self.nb_seg = 1;
        }
    }

    /// OCCT Intf_Tool::LinBox (Intf_Tool.cxx L746-963).
    pub fn lin_box(&mut self, l: &Line3, domain: &BndBox, box_lin: &mut BndBox) {
        self.nb_seg = 0;
        box_lin.set_void();
        if domain.is_whole() {
            let loc = l.origin;
            let dir = l.direction.normalize_or_zero();
            box_lin.add_point(loc);
            box_lin.add_point(loc + dir);
            box_lin.add_point(loc - dir);
            self.nb_seg = 1;
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            return;
        } else if domain.is_void() {
            return;
        }

        let (xmin, ymin, zmin, xmax, ymax, zmax) =
            domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        let mut xmin_b = 0.0;
        let mut xmax_b = 0.0;
        let mut ymin_b = 0.0;
        let mut ymax_b = 0.0;
        let mut zmin_b = 0.0;
        let mut zmax_b = 0.0;
        let mut parmin = -INFINITE_VALUE;
        let mut parmax = INFINITE_VALUE;
        let (mut x_to_set, mut y_to_set, mut z_to_set) = (false, false, false);

        if l.direction.x > 0.0 {
            parmin = if domain.is_open_xmin() {
                -INFINITE_VALUE
            } else {
                (xmin - l.origin.x) / l.direction.x
            };
            parmax = if domain.is_open_xmax() {
                INFINITE_VALUE
            } else {
                (xmax - l.origin.x) / l.direction.x
            };
            x_to_set = true;
        } else if l.direction.x < 0.0 {
            parmin = if domain.is_open_xmax() {
                -INFINITE_VALUE
            } else {
                (xmax - l.origin.x) / l.direction.x
            };
            parmax = if domain.is_open_xmin() {
                INFINITE_VALUE
            } else {
                (xmin - l.origin.x) / l.direction.x
            };
            x_to_set = true;
        } else if l.origin.x < xmin || xmax < l.origin.x {
            return;
        } else {
            xmin_b = l.origin.x;
            xmax_b = l.origin.x;
            x_to_set = false;
        }

        if l.direction.y > 0.0 {
            let parcur = if domain.is_open_ymin() {
                -INFINITE_VALUE
            } else {
                (ymin - l.origin.y) / l.direction.y
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymax() {
                INFINITE_VALUE
            } else {
                (ymax - l.origin.y) / l.direction.y
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else if l.direction.y < 0.0 {
            let parcur = if domain.is_open_ymax() {
                -INFINITE_VALUE
            } else {
                (ymax - l.origin.y) / l.direction.y
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymin() {
                INFINITE_VALUE
            } else {
                (ymin - l.origin.y) / l.direction.y
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else if l.origin.y < ymin || ymax < l.origin.y {
            return;
        } else {
            ymin_b = l.origin.y;
            ymax_b = l.origin.y;
            y_to_set = false;
        }

        if l.direction.z > 0.0 {
            let parcur = if domain.is_open_zmin() {
                -INFINITE_VALUE
            } else {
                (zmin - l.origin.z) / l.direction.z
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_zmax() {
                INFINITE_VALUE
            } else {
                (zmax - l.origin.z) / l.direction.z
            };
            parmax = parmax.min(parcur);
            z_to_set = true;
        } else if l.direction.z < 0.0 {
            let parcur = if domain.is_open_zmax() {
                -INFINITE_VALUE
            } else {
                (zmax - l.origin.z) / l.direction.z
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_zmin() {
                INFINITE_VALUE
            } else {
                (zmin - l.origin.z) / l.direction.z
            };
            parmax = parmax.min(parcur);
            z_to_set = true;
        } else if l.origin.z < zmin || zmax < l.origin.z {
            return;
        } else {
            zmin_b = l.origin.z;
            zmax_b = l.origin.z;
            z_to_set = false;
        }

        self.nb_seg += 1;
        self.begin_on_curve[0] = parmin;
        self.end_on_curve[0] = parmax;

        if x_to_set {
            let par1 = l.origin.x + parmin * l.direction.x;
            let par2 = l.origin.x + parmax * l.direction.x;
            xmin_b = par1.min(par2);
            xmax_b = par1.max(par2);
        }
        if y_to_set {
            let par1 = l.origin.y + parmin * l.direction.y;
            let par2 = l.origin.y + parmax * l.direction.y;
            ymin_b = par1.min(par2);
            ymax_b = par1.max(par2);
        }
        if z_to_set {
            let par1 = l.origin.z + parmin * l.direction.z;
            let par2 = l.origin.z + parmax * l.direction.z;
            zmin_b = par1.min(par2);
            zmax_b = par1.max(par2);
        }

        box_lin.update(xmin_b, ymin_b, zmin_b, xmax_b, ymax_b, zmax_b);
    }

    /// OCCT Intf_Tool::HyprBox (Intf_Tool.cxx L967-1119) — the variant that
    /// clamps segment parameters to ±10.
    pub fn hypr_box(&mut self, hypr: &Hyperbola3, domain: &BndBox, box_hypr: &mut BndBox) {
        self.nb_seg = 0;
        box_hypr.set_void();

        if domain.is_whole() {
            box_hypr.set_whole();
            self.nb_seg = 1;
            self.begin_on_curve[0] = -100.0;
            self.end_on_curve[0] = 100.0;
            return;
        } else if domain.is_void() {
            return;
        }

        let nb_pi = self.inters_3d_hypr(hypr, domain);
        if nb_pi > 0 {
            let (mut xmin_b, mut ymin_b, mut zmin_b, mut xmax_b, mut ymax_b, mut zmax_b) =
                domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

            for npi in 0..nb_pi {
                xmin_b = xmin_b.min(self.xint[npi]);
                xmax_b = xmax_b.max(self.xint[npi]);
                ymin_b = ymin_b.min(self.yint[npi]);
                ymax_b = ymax_b.max(self.yint[npi]);
                zmin_b = zmin_b.min(self.zint[npi]);
                zmax_b = zmax_b.max(self.zint[npi]);
            }
            box_hypr.update(xmin_b, ymin_b, zmin_b, xmax_b, ymax_b, zmax_b);

            let mut sinan = 0.0;
            let mut out = true;

            for npi in 0..nb_pi {
                let (_, tan) = (hypr.point_at(self.parint[npi]), hypr.derivative_at(self.parint[npi]));
                sinan = match self.bord[npi] {
                    1 => DVec3::X.dot(tan),
                    2 => DVec3::Y.dot(tan),
                    3 => DVec3::Z.dot(tan),
                    4 => (-DVec3::X).dot(tan),
                    5 => (-DVec3::Y).dot(tan),
                    6 => (-DVec3::Z).dot(tan),
                    _ => 0.0,
                };
                if sinan.abs() > ANGULAR {
                    if sinan > 0.0 {
                        out = false;
                        self.begin_on_curve[self.nb_seg] = self.parint[npi];
                        self.end_on_curve[self.nb_seg] = 10.0;
                        self.nb_seg += 1;
                    } else {
                        if out {
                            self.begin_on_curve[self.nb_seg] = -10.0;
                            self.nb_seg += 1;
                        }
                        self.end_on_curve[self.nb_seg - 1] = self.parint[npi];
                        out = true;

                        let mut ipmin = -10.0;
                        if self.begin_on_curve[self.nb_seg - 1] > ipmin {
                            ipmin = self.begin_on_curve[self.nb_seg - 1];
                        }
                        let mut ipmax = 10.0;
                        if self.end_on_curve[self.nb_seg - 1] < ipmax {
                            ipmax = self.end_on_curve[self.nb_seg - 1];
                        }
                        ipmin = ipmin * 10.0 + 1.0;
                        ipmax = ipmax * 10.0 - 1.0;
                        let mut ip = ipmin;
                        let mut pas = 1.0;
                        while ip <= ipmax {
                            box_hypr.add_point(hypr.point_at(ip / 10.0));
                            pas = 10.0;
                            if ip.abs() <= 10.0 {
                                pas = 1.0;
                            }
                            ip += pas;
                        }
                    }
                }
            }
        } else if !domain.is_out_point(hypr.point_at(0.0)) {
            *box_hypr = domain.clone();
            self.begin_on_curve[0] = -100.0;
            self.end_on_curve[0] = 100.0;
            self.nb_seg = 1;
        }
    }

    /// OCCT Intf_Tool::ParabBox (Intf_Tool.cxx L1489-1630).
    pub fn parab_box(&mut self, parab: &Parabola3, domain: &BndBox, box_parab: &mut BndBox) {
        self.nb_seg = 0;
        box_parab.set_void();
        if domain.is_whole() {
            box_parab.set_whole();
            self.nb_seg = 1;
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            return;
        } else if domain.is_void() {
            return;
        }

        let nb_pi = self.inters_3d_parab(parab, domain);

        if nb_pi > 0 {
            let (mut xmin_b, mut ymin_b, mut zmin_b, mut xmax_b, mut ymax_b, mut zmax_b) =
                domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

            for npi in 0..nb_pi {
                xmin_b = xmin_b.min(self.xint[npi]);
                xmax_b = xmax_b.max(self.xint[npi]);
                ymin_b = ymin_b.min(self.yint[npi]);
                ymax_b = ymax_b.max(self.yint[npi]);
                zmin_b = zmin_b.min(self.zint[npi]);
                zmax_b = zmax_b.max(self.zint[npi]);
            }
            box_parab.update(xmin_b, ymin_b, zmin_b, xmax_b, ymax_b, zmax_b);

            let mut sinan = 0.0;
            let mut out = true;

            for npi in 0..nb_pi {
                let (_, tan) = (parab.point_at(self.parint[npi]), parab.derivative_at(self.parint[npi]));
                sinan = match self.bord[npi] {
                    1 => DVec3::X.dot(tan),
                    2 => DVec3::Y.dot(tan),
                    3 => DVec3::Z.dot(tan),
                    4 => (-DVec3::X).dot(tan),
                    5 => (-DVec3::Y).dot(tan),
                    6 => (-DVec3::Z).dot(tan),
                    _ => 0.0,
                };
                if sinan.abs() > ANGULAR {
                    if sinan > 0.0 {
                        if self.nb_seg < 6 {
                            out = false;
                            self.begin_on_curve[self.nb_seg] = self.parint[npi];
                            self.nb_seg += 1;
                        }
                    } else {
                        if out && self.nb_seg < 6 {
                            self.begin_on_curve[self.nb_seg] = -INFINITE_VALUE;
                            self.nb_seg += 1;
                        }
                        if self.nb_seg > 0 {
                            self.end_on_curve[self.nb_seg - 1] = self.parint[npi];
                        }
                        out = true;

                        let mut ipmin = if self.begin_on_curve[self.nb_seg - 1] < -10.0 {
                            -10
                        } else {
                            self.begin_on_curve[self.nb_seg - 1] as i32
                        };
                        let mut ipmax = if self.end_on_curve[self.nb_seg - 1] > 10.0 {
                            10
                        } else {
                            self.end_on_curve[self.nb_seg - 1] as i32
                        };
                        ipmin = ipmin * 10 + 1;
                        ipmax = ipmax * 10 - 1;
                        let mut ip = ipmin;
                        let mut pas = 1;
                        while ip <= ipmax {
                            box_parab.add_point(parab.point_at(ip as f64 / 10.0));
                            if ip.abs() <= 10 {
                                pas = 1;
                            } else {
                                pas = 10;
                            }
                            ip += pas;
                        }
                    }
                }
            }
            if !out && self.nb_seg > 0 {
                self.end_on_curve[self.nb_seg - 1] = INFINITE_VALUE;
            }
        } else if !domain.is_out_point(parab.point_at(0.0)) {
            *box_parab = domain.clone();
            self.begin_on_curve[0] = -INFINITE_VALUE;
            self.end_on_curve[0] = INFINITE_VALUE;
            self.nb_seg = 1;
        }
    }

    /// OCCT Intf_Tool::NbSegments() (L1634-1637).
    pub fn nb_segments(&self) -> usize {
        self.nb_seg
    }

    /// OCCT Intf_Tool::BeginParam(SegmentNum) (L1641-1645) — 1-based.
    pub fn begin_param(&self, segment_num: usize) -> f64 {
        assert!(
            segment_num >= 1 && segment_num <= self.nb_seg,
            "Intf_Tool::BeginParam"
        );
        self.begin_on_curve[segment_num - 1]
    }

    /// OCCT Intf_Tool::EndParam(SegmentNum) (L1649-1653) — 1-based.
    pub fn end_param(&self, segment_num: usize) -> f64 {
        assert!(
            segment_num >= 1 && segment_num <= self.nb_seg,
            "Intf_Tool::EndParam"
        );
        self.end_on_curve[segment_num - 1]
    }

    /// OCCT Intf_Tool::Inters2d(const gp_Hypr2d&, const Bnd_Box2d&)
    /// (Intf_Tool.cxx L373-473) — the hyperbola intersections with the four
    /// domain boundary lines.
    fn inters_2d_hypr(&mut self, curv: &Hyperbola2d, domain: &BndBox2d) -> usize {
        let mut nbpi = 0;
        let (xmin, ymin, xmax, ymax) = domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0));

        if !domain.is_open_ymax() {
            let l1 = Line2d::new(DVec2::new(0.0, ymax), DVec2::new(-1.0, 0.0));
            let mut inter1 = AnaIntersection2d::new();
            inter1.perform_hyperbola_conic(curv, &Conic2d::from_line(&l1));
            if inter1.is_done() && !inter1.is_empty() {
                for npi in 1..=inter1.nb_points() {
                    self.xint[nbpi] = inter1.point(npi).value().x;
                    if xmin < self.xint[nbpi] && self.xint[nbpi] <= xmax {
                        self.yint[nbpi] = ymax;
                        self.parint[nbpi] = inter1.point(npi).param_on_first();
                        self.bord[nbpi] = 1;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_xmin() {
            let l2 = Line2d::new(DVec2::new(xmin, 0.0), DVec2::new(0.0, -1.0));
            let mut inter2 = AnaIntersection2d::new();
            inter2.perform_hyperbola_conic(curv, &Conic2d::from_line(&l2));
            if inter2.is_done() && !inter2.is_empty() {
                for npi in 1..=inter2.nb_points() {
                    self.yint[nbpi] = inter2.point(npi).value().y;
                    if ymin < self.yint[nbpi] && self.yint[nbpi] <= ymax {
                        self.xint[nbpi] = xmin;
                        self.parint[nbpi] = inter2.point(npi).param_on_first();
                        self.bord[nbpi] = 2;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_ymin() {
            let l3 = Line2d::new(DVec2::new(0.0, ymin), DVec2::new(1.0, 0.0));
            let mut inter3 = AnaIntersection2d::new();
            inter3.perform_hyperbola_conic(curv, &Conic2d::from_line(&l3));
            if inter3.is_done() && !inter3.is_empty() {
                for npi in 1..=inter3.nb_points() {
                    self.xint[nbpi] = inter3.point(npi).value().x;
                    if xmin <= self.xint[nbpi] && self.xint[nbpi] < xmax {
                        self.yint[nbpi] = ymin;
                        self.parint[nbpi] = inter3.point(npi).param_on_first();
                        self.bord[nbpi] = 3;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_xmax() {
            let l4 = Line2d::new(DVec2::new(xmax, 0.0), DVec2::new(0.0, 1.0));
            let mut inter4 = AnaIntersection2d::new();
            inter4.perform_hyperbola_conic(curv, &Conic2d::from_line(&l4));
            if inter4.is_done() && !inter4.is_empty() {
                for npi in 1..=inter4.nb_points() {
                    self.yint[nbpi] = inter4.point(npi).value().y;
                    if ymin <= self.yint[nbpi] && self.yint[nbpi] < ymax {
                        self.xint[nbpi] = xmax;
                        self.parint[nbpi] = inter4.point(npi).param_on_first();
                        self.bord[nbpi] = 4;
                        nbpi += 1;
                    }
                }
            }
        }
        nbpi
    }

    /// OCCT Intf_Tool::Inters2d(const gp_Parab2d&, const Bnd_Box2d&)
    /// (Intf_Tool.cxx L642-742).
    fn inters_2d_parab(&mut self, curv: &Parabola2d, domain: &BndBox2d) -> usize {
        let mut nbpi = 0;
        let (xmin, ymin, xmax, ymax) = domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0));

        if !domain.is_open_ymax() {
            let l1 = Line2d::new(DVec2::new(0.0, ymax), DVec2::new(-1.0, 0.0));
            let mut inter1 = AnaIntersection2d::new();
            inter1.perform_parabola_conic(curv, &Conic2d::from_line(&l1));
            if inter1.is_done() && !inter1.is_empty() {
                for npi in 1..=inter1.nb_points() {
                    self.xint[nbpi] = inter1.point(npi).value().x;
                    if xmin < self.xint[nbpi] && self.xint[nbpi] <= xmax {
                        self.yint[nbpi] = ymax;
                        self.parint[nbpi] = inter1.point(npi).param_on_first();
                        self.bord[nbpi] = 1;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_xmin() {
            let l2 = Line2d::new(DVec2::new(xmin, 0.0), DVec2::new(0.0, -1.0));
            let mut inter2 = AnaIntersection2d::new();
            inter2.perform_parabola_conic(curv, &Conic2d::from_line(&l2));
            if inter2.is_done() && !inter2.is_empty() {
                for npi in 1..=inter2.nb_points() {
                    self.yint[nbpi] = inter2.point(npi).value().y;
                    if ymin < self.yint[nbpi] && self.yint[nbpi] <= ymax {
                        self.xint[nbpi] = xmin;
                        self.parint[nbpi] = inter2.point(npi).param_on_first();
                        self.bord[nbpi] = 2;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_ymin() {
            let l3 = Line2d::new(DVec2::new(0.0, ymin), DVec2::new(1.0, 0.0));
            let mut inter3 = AnaIntersection2d::new();
            inter3.perform_parabola_conic(curv, &Conic2d::from_line(&l3));
            if inter3.is_done() && !inter3.is_empty() {
                for npi in 1..=inter3.nb_points() {
                    self.xint[nbpi] = inter3.point(npi).value().x;
                    if xmin <= self.xint[nbpi] && self.xint[nbpi] < xmax {
                        self.yint[nbpi] = ymin;
                        self.parint[nbpi] = inter3.point(npi).param_on_first();
                        self.bord[nbpi] = 3;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_xmax() {
            let l4 = Line2d::new(DVec2::new(xmax, 0.0), DVec2::new(0.0, 1.0));
            let mut inter4 = AnaIntersection2d::new();
            inter4.perform_parabola_conic(curv, &Conic2d::from_line(&l4));
            if inter4.is_done() && !inter4.is_empty() {
                for npi in 1..=inter4.nb_points() {
                    self.yint[nbpi] = inter4.point(npi).value().y;
                    if ymin <= self.yint[nbpi] && self.yint[nbpi] < ymax {
                        self.xint[nbpi] = xmax;
                        self.parint[nbpi] = inter4.point(npi).param_on_first();
                        self.bord[nbpi] = 4;
                        nbpi += 1;
                    }
                }
            }
        }
        nbpi
    }

    /// OCCT Intf_Tool::Inters3d(const gp_Hypr&, const Bnd_Box&)
    /// (Intf_Tool.cxx L1123-1302) — the hyperbola intersections with the six
    /// domain boundary planes (IntAna_IntConicQuad vs plane quadrics).
    fn inters_3d_hypr(&mut self, curv: &Hyperbola3, domain: &BndBox) -> usize {
        let mut nbpi = 0;
        let (xmin, ymin, zmin, xmax, ymax, zmax) =
            domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

        macro_rules! intersect_plane {
            ($a:expr, $b:expr, $c:expr, $d:expr, $fixed:expr, $cond:expr, $border:expr, $set:expr) => {
                let plane = plane_from_coeffs($a, $b, $c, $d);
                let quad = crate::geomalgo::int_surf::quadric::Quadric::from_plane(&plane);
                if let Some((in_quadric, pts)) =
                    crate::geomalgo::int_patch::int_conic_quad::intersect_hyperbola_quadric(curv, &quad)
                {
                    if !in_quadric {
                        for (pnt, param) in pts {
                            let (c1, c2) = if $set == 0 {
                                (pnt.y, pnt.z)
                            } else if $set == 1 {
                                (pnt.x, pnt.z)
                            } else {
                                (pnt.x, pnt.y)
                            };
                            if $cond(c1, c2) {
                                if $set == 0 {
                                    self.xint[nbpi] = $fixed;
                                    self.yint[nbpi] = c1;
                                    self.zint[nbpi] = c2;
                                } else if $set == 1 {
                                    self.xint[nbpi] = c1;
                                    self.yint[nbpi] = $fixed;
                                    self.zint[nbpi] = c2;
                                } else {
                                    self.xint[nbpi] = c1;
                                    self.yint[nbpi] = c2;
                                    self.zint[nbpi] = $fixed;
                                }
                                self.parint[nbpi] = param;
                                self.bord[nbpi] = $border;
                                nbpi += 1;
                            }
                        }
                    }
                }
            };
        }

        if !domain.is_open_xmin() {
            intersect_plane!(1.0, 0.0, 0.0, -xmin, xmin,
                |c1, c2| ymin <= c1 && c1 < ymax && zmin <= c2 && c2 < zmax, 1, 0);
        }
        if !domain.is_open_ymin() {
            intersect_plane!(0.0, 1.0, 0.0, -ymin, ymin,
                |c1, c2| xmin < c1 && c1 <= xmax && zmin <= c2 && c2 < zmax, 2, 1);
        }
        if !domain.is_open_zmin() {
            intersect_plane!(0.0, 0.0, 1.0, -zmin, zmin,
                |c1, c2| xmin < c1 && c1 <= xmax && ymin < c2 && c2 <= ymax, 3, 2);
        }
        if !domain.is_open_xmax() {
            intersect_plane!(-1.0, 0.0, 0.0, xmax, xmax,
                |c1, c2| ymin < c1 && c1 <= ymax && zmin < c2 && c2 <= zmax, 4, 0);
        }
        if !domain.is_open_ymax() {
            intersect_plane!(0.0, -1.0, 0.0, ymax, ymax,
                |c1, c2| xmin <= c1 && c1 < xmax && zmin < c2 && c2 <= zmax, 5, 1);
        }
        if !domain.is_open_zmax() {
            intersect_plane!(0.0, 0.0, -1.0, zmax, zmax,
                |c1, c2| xmin <= c1 && c1 < xmax && ymin <= c2 && c2 < ymax, 6, 2);
        }

        // OCCT L1269-1301: sort parint and drop matched (duplicate) values.
        let mut a_nb_diff_points = nbpi as i64;
        let mut i = nbpi as i64 - 1;
        while i > 0 {
            for j in 0..i {
                if self.parint[i as usize] <= self.parint[j as usize] {
                    self.parint.swap(i as usize, j as usize);
                    self.zint.swap(i as usize, j as usize);
                    self.yint.swap(i as usize, j as usize);
                    self.xint.swap(i as usize, j as usize);
                    self.bord.swap(i as usize, j as usize);
                }
                if (i < nbpi as i64 - 1) && is_equal(self.parint[i as usize], self.parint[i as usize + 1]) {
                    a_nb_diff_points -= 1;
                    for k in i..a_nb_diff_points {
                        self.parint[k as usize] = self.parint[k as usize + 1];
                        self.zint[k as usize] = self.zint[k as usize + 1];
                        self.yint[k as usize] = self.yint[k as usize + 1];
                        self.xint[k as usize] = self.xint[k as usize + 1];
                        self.bord[k as usize] = self.bord[k as usize + 1];
                    }
                }
            }
            i -= 1;
        }
        a_nb_diff_points.max(0) as usize
    }

    /// OCCT Intf_Tool::Inters3d(const gp_Parab&, const Bnd_Box&)
    /// (Intf_Tool.cxx L1306-1485).
    fn inters_3d_parab(&mut self, curv: &Parabola3, domain: &BndBox) -> usize {
        let mut nbpi = 0;
        let (xmin, ymin, zmin, xmax, ymax, zmax) =
            domain.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

        macro_rules! intersect_plane {
            ($a:expr, $b:expr, $c:expr, $d:expr, $fixed:expr, $cond:expr, $border:expr, $set:expr) => {
                let plane = plane_from_coeffs($a, $b, $c, $d);
                let quad = crate::geomalgo::int_surf::quadric::Quadric::from_plane(&plane);
                if let Some((in_quadric, pts)) =
                    crate::geomalgo::int_patch::int_conic_quad::intersect_parabola_quadric(curv, &quad)
                {
                    if !in_quadric {
                        for (pnt, param) in pts {
                            let (c1, c2) = if $set == 0 {
                                (pnt.y, pnt.z)
                            } else if $set == 1 {
                                (pnt.x, pnt.z)
                            } else {
                                (pnt.x, pnt.y)
                            };
                            if $cond(c1, c2) {
                                if $set == 0 {
                                    self.xint[nbpi] = $fixed;
                                    self.yint[nbpi] = c1;
                                    self.zint[nbpi] = c2;
                                } else if $set == 1 {
                                    self.xint[nbpi] = c1;
                                    self.yint[nbpi] = $fixed;
                                    self.zint[nbpi] = c2;
                                } else {
                                    self.xint[nbpi] = c1;
                                    self.yint[nbpi] = c2;
                                    self.zint[nbpi] = $fixed;
                                }
                                self.parint[nbpi] = param;
                                self.bord[nbpi] = $border;
                                nbpi += 1;
                            }
                        }
                    }
                }
            };
        }

        if !domain.is_open_xmin() {
            intersect_plane!(1.0, 0.0, 0.0, -xmin, xmin,
                |c1, c2| ymin <= c1 && c1 < ymax && zmin <= c2 && c2 < zmax, 1, 0);
        }
        if !domain.is_open_ymin() {
            intersect_plane!(0.0, 1.0, 0.0, -ymin, ymin,
                |c1, c2| xmin < c1 && c1 <= xmax && zmin <= c2 && c2 < zmax, 2, 1);
        }
        if !domain.is_open_zmin() {
            intersect_plane!(0.0, 0.0, 1.0, -zmin, zmin,
                |c1, c2| xmin < c1 && c1 <= xmax && ymin < c2 && c2 <= ymax, 3, 2);
        }
        if !domain.is_open_xmax() {
            intersect_plane!(-1.0, 0.0, 0.0, xmax, xmax,
                |c1, c2| ymin < c1 && c1 <= ymax && zmin < c2 && c2 <= zmax, 4, 0);
        }
        if !domain.is_open_ymax() {
            intersect_plane!(0.0, -1.0, 0.0, ymax, ymax,
                |c1, c2| xmin <= c1 && c1 < xmax && zmin < c2 && c2 <= zmax, 5, 1);
        }
        if !domain.is_open_zmax() {
            intersect_plane!(0.0, 0.0, -1.0, zmax, zmax,
                |c1, c2| xmin <= c1 && c1 < xmax && ymin <= c2 && c2 < ymax, 6, 2);
        }

        // OCCT L1452-1484: sort parint and drop matched (duplicate) values.
        let mut a_nb_diff_points = nbpi as i64;
        let mut i = nbpi as i64 - 1;
        while i > 0 {
            for j in 0..i {
                if self.parint[i as usize] <= self.parint[j as usize] {
                    self.parint.swap(i as usize, j as usize);
                    self.zint.swap(i as usize, j as usize);
                    self.yint.swap(i as usize, j as usize);
                    self.xint.swap(i as usize, j as usize);
                    self.bord.swap(i as usize, j as usize);
                }
                if (i < nbpi as i64 - 1) && is_equal(self.parint[i as usize], self.parint[i as usize + 1]) {
                    a_nb_diff_points -= 1;
                    for k in i..a_nb_diff_points {
                        self.parint[k as usize] = self.parint[k as usize + 1];
                        self.zint[k as usize] = self.zint[k as usize + 1];
                        self.yint[k as usize] = self.yint[k as usize + 1];
                        self.xint[k as usize] = self.xint[k as usize + 1];
                        self.bord[k as usize] = self.bord[k as usize + 1];
                    }
                }
            }
            i -= 1;
        }
        a_nb_diff_points.max(0) as usize
    }
}

impl Default for IntfTool {
    fn default() -> Self {
        Self::new()
    }
}
