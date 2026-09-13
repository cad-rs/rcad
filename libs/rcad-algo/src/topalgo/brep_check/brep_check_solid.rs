//! OCCT BRepCheck_Solid (TKTopAlgo/BRepCheck).
//!
//! Source: `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Solid.cxx`
//! (L38-330) and `BRepCheck_Solid.hxx` (L28-44).

use rcad_kernel::geom::CurveEval;
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, Orientation, ShapeType, TShape, TSolidData};
use std::sync::Arc;

use crate::topalgo::brep_class3d::solid_classifier::SolidClassifier;

use crate::brep_algo::tool::brep_tool_tolerance;
use super::brep_check_result::{
    brep_check_add, explorer, iterator_subshapes,
    BRepCheckResultBase, BRepCheckStatus,
};
use super::brep_check_wire::ShapeSet;

/// OCCT Solid.cxx L129: `aPAR_T = 0.43213918; // 10*e^(-PI)`.
const A_PAR_T: f64 = 0.43213918;

/// OCCT BRepCheck_ToolSolid (Solid.cxx L64-152) — the per-shell tool holding
/// the solid classifier.
pub struct BRepCheckToolSolid {
    /// OCCT myIsHole.
    my_is_hole: bool,
    /// OCCT myPnt.
    my_pnt: glam::DVec3,
    /// OCCT myPntTol.
    my_pnt_tol: f64,
    /// OCCT mySolid.
    my_solid: Shape,
    /// OCCT myHSC (the classifier handle). The rcad classifier is carried by
    /// value; the Option models the OCCT null-handle period before Init().
    my_hsc: Option<SolidClassifier>,
}

impl Default for BRepCheckToolSolid {
    /// OCCT BRepCheck_ToolSolid::BRepCheck_ToolSolid (Solid.cxx L70-76).
    fn default() -> Self {
        BRepCheckToolSolid {
            my_is_hole: false,
            my_pnt_tol: rcad_kernel::CONFUSION,
            my_pnt: glam::DVec3::new(-1., -1., -1.),
            my_solid: Shape::null(),
            my_hsc: None,
        }
    }
}

impl BRepCheckToolSolid {
    /// OCCT SetSolid (Solid.cxx L80-82).
    pub fn set_solid(&mut self, a_z: &Shape) {
        self.my_solid = a_z.clone();
    }

    /// OCCT InnerPoint (Solid.cxx L89-91).
    pub fn inner_point(&self) -> glam::DVec3 {
        self.my_pnt
    }

    /// OCCT CheckTol (Solid.cxx L92-94).
    pub fn check_tol(&self) -> f64 {
        self.my_pnt_tol
    }

    /// OCCT IsOut (Solid.cxx L96-108) — whether the other tool's inner point
    /// falls outside this solid (non-const in OCCT: it mutates the member
    /// classifier through the handle).
    pub fn is_out(&mut self, _brep: &BRep, a_other: &BRepCheckToolSolid) -> bool {
        let Some(a_sc) = self.my_hsc.as_mut() else {
            return false;
        };
        // OCCT L103: aSC.Perform(aOther.InnerPoint(), aOther.CheckTol()).
        a_sc.perform(a_other.inner_point(), a_other.check_tol());
        let a_state = a_sc.state();
        // OCCT L105: bFlag = (aState == TopAbs_OUT).
        a_state == 1 // TopAbs_OUT
    }

    /// OCCT Init (Solid.cxx L112-143).
    pub fn init(&mut self, brep: &BRep) {
        // OCCT L118-126: the classifier load + the infinite point.
        let mut a_sc = SolidClassifier::new();
        a_sc.load(&self.my_solid);
        // OCCT L125: aSC.PerformInfinitePoint(::RealSmall()).
        a_sc.perform_infinite_point(f64::MIN_POSITIVE);
        self.my_is_hole = a_sc.state() == 0; // TopAbs_IN
        self.my_hsc = Some(a_sc);

        // OCCT L128-142: myPnt from the first non-degenerated edge.
        for a_e in explorer(brep, &self.my_solid, ShapeType::Edge) {
            let degenerated = brep.is_edge_degenerated(&a_e);
            if !degenerated {
                let Some(ed) = a_e.as_edge() else { continue };
                if let Some(a_c3d) = &ed.curve {
                    let (a_t1, a_t2) = (ed.range[0], ed.range[1]);
                    // OCCT L137: aT = (1. - aPAR_T) * aT1 + aPAR_T * aT2.
                    let a_t = (1. - A_PAR_T) * a_t1 + A_PAR_T * a_t2;
                    // OCCT L138: myPnt = aC3D->Value(aT) (the curve value in
                    // the edge frame; the OCCT Geom_Curve value is then the
                    // located one — rcad materialises locations, and the
                    // edge location applies on read).
                    self.my_pnt = a_c3d.point_at(a_t);
                    // OCCT L139: myPntTol = BRep_Tool::Tolerance(aE).
                    self.my_pnt_tol = brep_tool_tolerance(&a_e);
                    break;
                }
            }
        }
    }
}

/// OCCT BRepCheck_Solid (Solid.hxx L28-44).
#[derive(Debug)]
pub struct BRepCheckSolid {
    /// OCCT protected base (Result.hxx L82-90).
    pub base: BRepCheckResultBase,
}

impl BRepCheckSolid {
    /// OCCT BRepCheck_Solid::BRepCheck_Solid(const TopoDS_Solid& S)
    /// (Solid.cxx L161-164).
    pub fn new(brep: &BRep, s: &Shape) -> Self {
        let mut r = BRepCheckSolid {
            base: BRepCheckResultBase::new(),
        };
        r.base.init(s);
        r.minimum(brep);
        r
    }

    /// OCCT BRepCheck_Solid::Blind (Solid.cxx L168-175).
    pub fn blind(&mut self) {
        if !self.base.my_blind {
            // nothing more than in the minimum
            self.base.my_blind = true;
        }
    }

    /// OCCT BRepCheck_Solid::InContext (Solid.cxx L179) — empty.
    pub fn in_context(&mut self, _brep: &BRep, _s: &Shape) {}

    /// OCCT BRepCheck_Solid::Minimum (Solid.cxx L183-330).
    pub fn minimum(&mut self, brep: &BRep) {
        if self.base.my_min {
            return;
        }
        // OCCT L189.
        self.base.my_min = true;
        let my_shape = self.base.my_shape.clone();

        // OCCT L201-204.
        self.base.my_map.bound(&my_shape);
        {
            let a_lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");
            a_lst.push(BRepCheckStatus::NoError);
        }

        //-------------------------------------------------
        // OCCT L207-218: 1. InvalidImbricationOfShells.
        let mut b_found = false;
        let mut a_mss = ShapeSet::new();
        for a_f in explorer(brep, &my_shape, ShapeType::Face) {
            if b_found {
                break;
            }
            if !a_mss.add(&a_f) {
                let a_lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Minimum: myShape must be bound");
                brep_check_add(a_lst, BRepCheckStatus::InvalidImbricationOfShells);
                b_found = !b_found;
            }
        }

        //-------------------------------------------------
        // OCCT L220-276: 2. too many growths / shapes out of the solid.
        let mut i_cnt_sh = 0i32;
        let mut i_cnt_sh_int = 0i32;
        let mut a_vts: Vec<BRepCheckToolSolid> = Vec::new();
        for a_sx in iterator_subshapes(brep, &my_shape) {
            if a_sx.shape_type() != ShapeType::Shell {
                // OCCT L233-238.
                let a_or = a_sx.orientation;
                if a_or != Orientation::Internal {
                    let a_lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("Minimum: myShape must be bound");
                    brep_check_add(a_lst, BRepCheckStatus::BadOrientationOfSubshape);
                }
                continue;
            }

            let a_sh = a_sx.clone();

            // OCCT L243-258: skip internal shells.
            b_found = false;
            for a_f in iterator_subshapes(brep, &a_sh) {
                if b_found {
                    break;
                }
                let a_or = a_f.orientation;
                if a_or == Orientation::Internal {
                    b_found = !b_found;
                }
            }
            if b_found {
                i_cnt_sh_int += 1;
                continue;
            }

            // OCCT L260-261.
            i_cnt_sh += 1;

            // OCCT L263-266: skip not closed shells (BRep_Tool::IsClosed).
            if !breptool_is_closed_shell(brep, &a_sh) {
                continue;
            }

            // OCCT L268-275: aBB.MakeSolid(aZ); aBB.Add(aZ, aSh).
            let a_z = Shape {
                data: Arc::new(TShape::Solid(TSolidData {
                    my_shapes: vec![a_sh.clone()],
                    flags: rcad_kernel::topods::tshape_flags::DEFAULT,
                    shells: vec![a_sh.clone()],
                    internal_vertices: Vec::new(),
                    internal_edges: Vec::new(),
                })),
                index: usize::MAX,
                location: 0,
                orientation: Orientation::Forward,
            };

            let mut a_ts = BRepCheckToolSolid::default();
            a_ts.set_solid(&a_z);
            a_vts.push(a_ts);
        }

        // OCCT L278-282.
        if i_cnt_sh == 0 && i_cnt_sh_int > 0 {
            // all shells in the solid are internal
            let a_lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");
            brep_check_add(a_lst, BRepCheckStatus::BadOrientationOfSubshape);
        }

        // OCCT L284-288.
        let a_nb_vts = a_vts.len() as i32;
        if a_nb_vts < 2 {
            return;
        }

        // OCCT L290-307.
        let mut a_nb_vts1 = 0i32;
        for i in 0..a_nb_vts as usize {
            let a_ts = &mut a_vts[i];
            a_ts.init(brep);
            let b_is_hole = a_ts.is_hole();
            if !b_is_hole {
                a_nb_vts1 += 1;
                if a_nb_vts1 > 1 {
                    // Too many growths
                    let a_lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("Minimum: myShape must be bound");
                    brep_check_add(a_lst, BRepCheckStatus::EnclosedRegion);
                    break;
                }
            }
        }

        // OCCT L309-327.
        b_found = false;
        let a_nb_vts1 = a_nb_vts - 1;
        let mut i: i32 = 0;
        while !b_found && i < a_nb_vts1 {
            let mut j = i + 1;
            while !b_found && j < a_nb_vts {
                // OCCT L319: bFlag = aTSi.IsOut(aTSj) — the split_at_mut keeps
                // the OCCT IsOut(aOther) call shape (i < j).
                let (left_part, right_part) = a_vts.split_at_mut(j as usize);
                let b_flag = left_part[i as usize].is_out(brep, &right_part[0]);
                if b_flag {
                    // smt of solid is out of solid
                    let a_lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("Minimum: myShape must be bound");
                    brep_check_add(a_lst, BRepCheckStatus::SubshapeNotInShape);
                    b_found = !b_found;
                }
                j += 1;
            }
            i += 1;
        }
        // OCCT L329: // myMin = true;
    }
}

impl BRepCheckToolSolid {
    /// OCCT IsHole (Solid.cxx L86-88).
    pub fn is_hole(&self) -> bool {
        self.my_is_hole
    }
}

/// OCCT BRep_Tool::IsClosed(theShape) for the SHELL case
/// (BRep_Tool.cxx L1707-1730): the flag is ignored for shells — the walk
/// counts each non-degenerated, non-INTERNAL/EXTERNAL edge occurrence.
pub fn breptool_is_closed_shell(brep: &BRep, the_shape: &Shape) -> bool {
    let mut a_map = ShapeSet::new();
    let mut has_bound = false;
    for e in explorer(brep, the_shape, ShapeType::Edge) {
        if brep.is_edge_degenerated(&e)
            || e.orientation == Orientation::Internal
            || e.orientation == Orientation::External
        {
            continue;
        }
        has_bound = true;
        if !a_map.add(&e) {
            a_map.remove(&e);
        }
    }
    has_bound && a_map.extent() == 0
}
