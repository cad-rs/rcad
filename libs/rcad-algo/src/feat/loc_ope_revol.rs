// OCCT LocOpe_Revol.hxx L35-83 + LocOpe_Revol.cxx L47-312 — 1:1
// translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Revol.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Revol.cxx
//
// OCCT inheritance chain (LocOpe_Revol.hxx L35): none — LocOpe_Revol is a
// standalone class in the OCCT 8.0 sources at $OCCT_SRC (no
// LocOpe_GeneratedShape virtuals to translate; see loc_ope_prism.rs).
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//    TopTools_ShapeMapHasher> (myMap) — HashMap keyed by (TShape ptr,
//    Location); never iterated.
// 2. BRepSweep_Revol (TKPrim/BRepSweep) — the sweep engine is the
//    crate::brep_sweep translation (BRepSweepRevol); IntPerf consumes it
//    through the OCCT constructor/accessor surface (the former GAP panic
//    carrier is removed), shared with loc_ope_revolution_form.rs.
// 3. BRepTools_Modifier + BRepTools_TrsfModification — imported from the
//    topalgo translations (topalgo::brep_tools_modifier /
//    topalgo::brep_tools_modification); the OCCT `BRepTools_Modifier Modif;`
//    default constructor is the rcad `new(false)`.
// 4. gp_Trsf::SetRotation(Ax, Ang) has no rcad constructor (rcad Trsf
//    carries only translation/scale/displacement forms) — the
//    trsf_set_rotation helper below carries the GAP panic; it feeds the
//    BRepTools_TrsfModification form.
// 5. gp_Circ -> geom::Circle3; the OCCT default gp_Circ CAX carrier is the
//    zeroed Circle3 (gp_circ_default) — only read after FindCircle
//    succeeds.
// 6. Geom_Circle -> Curve3::Circle (Circle3); BarycCurve returns the OCCT
//    null-handle case as Option::None.
//
// The two value constructors declared in the OCCT header (hxx L42-47) have
// no definition in LocOpe_Revol.cxx — nothing to translate.
//
// first consumer: BRepFeat_MakeRevol (3b).

use crate::brep_sweep::BRepSweepRevol;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_build_shape::LocOpeBuildShape;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::topalgo::brep_tools_modification::BRepToolsTrsfModification;
use crate::topalgo::brep_tools_modifier::BRepToolsModifier;
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::geom::{ Circle3, Curve3 };
use rcad_kernel::math::gp::{ Ax1, Ax2, Trsf };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
use std::collections::HashMap;

/// OCCT gp_Circ default constructor — the zeroed Circle3 carrier
/// (architecture difference #5; only read after FindCircle succeeds).
fn gp_circ_default() -> Circle3 {
    Circle3 {
        center: DVec3::ZERO,
        normal: DVec3::ZERO,
        x_dir: DVec3::ZERO,
        y_dir: DVec3::ZERO,
        radius: 0.0,
    }
}

/// OCCT gp_Trsf::SetRotation(Ax, Ang) (architecture difference #4 — GAP).
pub(crate) fn trsf_set_rotation(_t: &mut Trsf, _ax: &Ax1, _ang: f64) {
    panic!("GAP: gp_Trsf::SetRotation (no rcad rotation Trsf constructor)");
}

/// OCCT static FindCircle(Ax, Pt, Ci) (cxx L294-312).
fn find_circle(ax: &Ax1, pt: DVec3, ci: &mut Circle3) -> bool {
    // OCCT cxx L297.
    let dax = ax.direction;
    // OCCT cxx L298: gp_Vec OP(Ax.Location(), Pt).
    let op = pt - ax.location;

    // OCCT cxx L300.
    let prm = op.dot(dax);

    // OCCT cxx L302.
    let prj = ax.location + prm * dax;
    // OCCT cxx L303: gp_Vec axx(prj, Pt).
    let axx = pt - prj;
    // OCCT cxx L304.
    let radius = axx.length();
    // OCCT cxx L305-308.
    if radius < rcad_kernel::precision::CONFUSION {
        return false;
    }
    // OCCT cxx L309-310: Ci.SetRadius(Radius);
    // Ci.SetPosition(gp_Ax2(prj, Dax, axx)).
    *ci = Circle3::new_with_ref_dir(prj, dax, radius, axx);
    // OCCT cxx L311.
    true
}

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT LocOpe_Revol (LocOpe_Revol.hxx L35-83).
pub struct LocOpeRevol {
    my_base: Shape,        // OCCT: myBase
    my_axis: Ax1,          // OCCT: myAxis (gp_Ax1)
    my_angle: f64,         // OCCT: myAngle
    my_ang_tra: f64,       // OCCT: myAngTra
    my_is_trans: bool,     // OCCT: myIsTrans
    my_done: bool,         // OCCT: myDone
    my_res: Shape,         // OCCT: myRes
    my_first_shape: Shape, // OCCT: myFirstShape
    my_last_shape: Shape,  // OCCT: myLastShape
    // OCCT: myMap (NCollection_DataMap) — arch. diff. #1
    my_map: HashMap<(u64, u32), Vec<Shape>>,
}

impl Default for LocOpeRevol {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeRevol {
    /// OCCT LocOpe_Revol::LocOpe_Revol() (cxx L47-53).
    pub fn new() -> Self {
        LocOpeRevol {
            my_base: Shape::null(),
            my_axis: Ax1::new(DVec3::ZERO, DVec3::Z),
            my_angle: 0.0,
            my_ang_tra: 0.0,
            my_is_trans: false,
            my_done: false,
            my_res: Shape::null(),
            my_first_shape: Shape::null(),
            my_last_shape: Shape::null(),
            my_map: HashMap::new(),
        }
    }

    /// OCCT LocOpe_Revol::Perform(Base, Axis, Angle) (cxx L57-70).
    pub fn perform(&mut self, the_base: &Shape, the_axis: &Ax1, the_angle: f64) {
        self.my_map.clear();
        self.my_first_shape = Shape::null();
        self.my_last_shape = Shape::null();
        self.my_base = Shape::null();
        self.my_res = Shape::null();
        self.my_base = the_base.clone();
        self.my_angle = the_angle;
        self.my_axis = *the_axis;
        self.my_ang_tra = 0.;
        self.my_is_trans = false;
        self.int_perf();
    }

    /// OCCT LocOpe_Revol::Perform(Base, Axis, Angle, angledec)
    /// (cxx L74-90) — Rust has no overloading: the `_angle_dec` suffix
    /// carries the OCCT parameter-name distinction.
    pub fn perform_angle_dec(
        &mut self,
        the_base: &Shape,
        the_axis: &Ax1,
        the_angle: f64,
        the_angledec: f64,
    ) {
        self.my_map.clear();
        self.my_first_shape = Shape::null();
        self.my_last_shape = Shape::null();
        self.my_base = Shape::null();
        self.my_res = Shape::null();
        self.my_base = the_base.clone();
        self.my_angle = the_angle;
        self.my_axis = *the_axis;
        self.my_ang_tra = the_angledec;
        self.my_is_trans = true;
        self.int_perf();
    }

    /// OCCT LocOpe_Revol::IntPerf() (cxx L94-213).
    fn int_perf(&mut self) {
        // OCCT cxx L96-97: TopoDS_Shape theBase = myBase;
        // BRepTools_Modifier Modif; (the OCCT default constructor is the
        // theMutableInput = false form).
        let mut the_base = self.my_base.clone();
        let mut modif = BRepToolsModifier::new(false);
        if self.my_is_trans {
            // OCCT cxx L100-101 (arch. diff. #4).
            let mut t = Trsf::identity();
            trsf_set_rotation(&mut t, &self.my_axis, self.my_ang_tra);
            // OCCT cxx L102-105 (arch. diff. #3).
            let mut modbase = BRepToolsTrsfModification::new(t);
            modif.init(&the_base);
            modif.perform(&mut modbase);
            the_base = modif.modified_shape(&the_base);
        }

        // OCCT cxx L108 (arch. diff. #2): BRepSweep_Revol theRevol(theBase,
        // myAxis, myAngle) — the OCCT default argument is C=false; the
        // gp_Ax1 axis travels as the (location, direction) tuple.
        let mut the_revol = BRepSweepRevol::with_angle(
            &the_base,
            (self.my_axis.location, self.my_axis.direction),
            self.my_angle,
            false,
        );

        // OCCT cxx L110-111.
        self.my_first_shape = the_revol.first_shape();
        self.my_last_shape = the_revol.last_shape();

        // OCCT cxx L113-131: the base-FACE branch.
        if the_base.shape_type() == ShapeType::Face {
            for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L119.
                if !self.my_map.contains_key(&shape_key(&edg)) {
                    // OCCT cxx L121-122.
                    self.my_map.insert(shape_key(&edg), Vec::new());
                    // OCCT cxx L123.
                    let desc = the_revol.shape_of(&edg);
                    // OCCT cxx L124-127: if (!desc.IsNull()) — the engine
                    // null result is the Vertex-typed dummy of Shape::null();
                    // the kernel is_null() (index == usize::MAX) also fires
                    // for pool-built real shapes, so the emptiness test is
                    // the type test (a real generated shape is non-Vertex).
                    if desc.shape_type() != ShapeType::Vertex {
                        self.my_map
                            .get_mut(&shape_key(&edg))
                            .expect("myMap(edg)")
                            .push(desc);
                    }
                }
            }
            // OCCT cxx L130.
            self.my_res = the_revol.shape();
        } else {
            // Cas base != FACE
            //
            // OCCT cxx L136-140: theEFMap.
            let mut the_ef_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            map_shapes_and_ancestors(
                &the_base,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_ef_map,
            );
            // OCCT cxx L141-142.
            let mut lfaces: Vec<Shape> = Vec::new();
            let mut toremove = false;
            // OCCT cxx L143-161.
            for i in 1..=the_ef_map.len() {
                // OCCT cxx L145: edg = theEFMap.FindKey(i).
                let (_, entry) = the_ef_map.get_index(i - 1).expect("theEFMap entry");
                let edg = entry.0.clone();
                // OCCT cxx L146-147.
                self.my_map.insert(shape_key(&edg), Vec::new());
                // OCCT cxx L148.
                let desc = the_revol.shape_of(&edg);
                // OCCT cxx L149-160: if (!desc.IsNull()) — the engine null
                // result is the Vertex-typed dummy of Shape::null(); the
                // kernel is_null() (index == usize::MAX) also fires for
                // pool-built real shapes, so the emptiness test is the type
                // test (a real generated shape is non-Vertex).
                if desc.shape_type() != ShapeType::Vertex {
                    if entry.1.len() >= 2 {
                        toremove = true;
                    } else {
                        self.my_map
                            .get_mut(&shape_key(&edg))
                            .expect("myMap(edg)")
                            .push(desc.clone());
                        lfaces.push(desc);
                    }
                }
            }
            // OCCT cxx L162-176.
            if toremove {
                // Rajouter les faces de FirstShape et LastShape
                for f in explorer(&self.my_first_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }
                for f in explorer(&self.my_last_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }

                // OCCT cxx L174-175.
                let bs = LocOpeBuildShape::with_faces(&lfaces);
                self.my_res = bs.shape().cloned().unwrap_or_else(Shape::null);
            } else {
                // OCCT cxx L179-193.
                for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT cxx L182.
                    if !self.my_map.contains_key(&shape_key(&edg)) {
                        // OCCT cxx L184-185.
                        self.my_map.insert(shape_key(&edg), Vec::new());
                        // OCCT cxx L186.
                        let desc = the_revol.shape_of(&edg);
                        // OCCT cxx L187-190: if (!desc.IsNull()) — the engine
                        // null result is the Vertex-typed dummy of
                        // Shape::null(); the kernel is_null()
                        // (index == usize::MAX) also fires for pool-built
                        // real shapes, so the emptiness test is the type
                        // test (a real generated shape is non-Vertex).
                        if desc.shape_type() != ShapeType::Vertex {
                            self.my_map
                                .get_mut(&shape_key(&edg))
                                .expect("myMap(edg)")
                                .push(desc);
                        }
                    }
                }
                // OCCT cxx L193.
                self.my_res = the_revol.shape();
            }
        }

        // OCCT cxx L197-211: m-a-j des descendants.
        if self.my_is_trans {
            for an_edg in explorer(&self.my_base.clone(), ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L203-204.
                let edg = an_edg.clone();
                let edgbis = modif.modified_shape(&edg);
                // OCCT cxx L205-209.
                if !edgbis.is_same(&edg) && self.my_map.contains_key(&shape_key(&edgbis)) {
                    let list = self
                        .my_map
                        .get(&shape_key(&edgbis))
                        .expect("myMap(edgbis)")
                        .clone();
                    self.my_map.insert(shape_key(&edg), list);
                    self.my_map.remove(&shape_key(&edgbis));
                }
            }
        }
        // OCCT cxx L212.
        self.my_done = true;
    }

    /// OCCT LocOpe_Revol::Shape() (cxx L217-224).
    pub fn shape(&self) -> Shape {
        if !self.my_done {
            // OCCT cxx L219-222: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_res.clone()
    }

    /// OCCT LocOpe_Revol::FirstShape() (cxx L228-231).
    pub fn first_shape(&self) -> Shape {
        self.my_first_shape.clone()
    }

    /// OCCT LocOpe_Revol::LastShape() (cxx L235-238).
    pub fn last_shape(&self) -> Shape {
        self.my_last_shape.clone()
    }

    /// OCCT LocOpe_Revol::Shapes(S) (cxx L242-245).
    pub fn shapes(&self, s: &Shape) -> &Vec<Shape> {
        // OCCT cxx L244: return myMap(S) — operator() asserts IsBound.
        self.my_map.get(&shape_key(s)).expect("Standard_NoSuchObject")
    }

    /// OCCT LocOpe_Revol::Curves(SCurves) (cxx L249-266).
    pub fn curves(&self, s_curves: &mut Vec<Curve3>) {
        // OCCT cxx L251.
        s_curves.clear();
        // OCCT cxx L252-253.
        let mut spt: Vec<glam::DVec3> = Vec::new();
        crate::feat::loc_ope::sample_edges(&self.my_first_shape, &mut spt);
        // OCCT cxx L254-265.
        for pvt in &spt {
            let mut cax = gp_circ_default();
            // OCCT cxx L258.
            if find_circle(&self.my_axis, *pvt, &mut cax) {
                // OCCT cxx L260-261: A2 = CAX.Position(); r = CAX.Radius().
                let a2 = Ax2::new(cax.center, cax.normal, cax.x_dir);
                let r = cax.radius;
                // OCCT cxx L262-263: Ci = Geom_Circle(A2, r).
                let ci =
                    Curve3::Circle(Circle3::new_with_ref_dir(a2.location, a2.direction, r, a2.x_direction));
                s_curves.push(ci);
            }
        }
    }

    /// OCCT LocOpe_Revol::BarycCurve() (cxx L270-290) — the OCCT null
    /// handle maps to Option::None (arch. diff. #6).
    pub fn baryc_curve(&self) -> Option<Curve3> {
        // OCCT cxx L272-274.
        let mut bar = glam::DVec3::ZERO;
        let mut spt: Vec<glam::DVec3> = Vec::new();
        crate::feat::loc_ope::sample_edges(&self.my_first_shape, &mut spt);
        for pvt in &spt {
            // OCCT cxx L277-278: bar.ChangeCoord() += pvt.XYZ().
            bar += *pvt;
        }
        // OCCT cxx L280.
        bar /= spt.len() as f64;
        // OCCT cxx L281-288.
        let mut cax = gp_circ_default();
        let mut the_ci: Option<Curve3> = None;
        if find_circle(&self.my_axis, bar, &mut cax) {
            let a2 = Ax2::new(cax.center, cax.normal, cax.x_dir);
            let r = cax.radius;
            the_ci = Some(Curve3::Circle(Circle3::new_with_ref_dir(
                a2.location,
                a2.direction,
                r,
                a2.x_direction,
            )));
        }
        // OCCT cxx L289: return theCi.
        the_ci
    }
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
