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

use glam::DVec3;
use rcad_kernel::core::precision::{COMPUTATIONAL, SQUARE_CONFUSION};

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
