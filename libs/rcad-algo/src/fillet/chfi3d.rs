//! ChFi3d_Builder / ChFi3d_FilBuilder / ChFi3d_ChBuilder — OCCT TKFillet
//! 1:1 translation.
//!
//! Sources:
//!   - ChFi3d_Builder.cxx (Compute L178-675, PerformFilletOnVertex,
//!     PerformSingularCorner, Reset L924-946, Generated L950-977)
//!   - ChFi3d_Builder_1.cxx (constructor L341-364, SetParams L367-380,
//!     SetContinuity L382-389, IsDone/Shape, Remove L1181-1199,
//!     Value L1201-1211, NbElements L1213-1229, Contains L1231-1280,
//!     Length/FirstVertex/LastVertex L1282-1310)
//!   - ChFi3d_FilBuilder.cxx (ctor L147-153, SetFilletShape L157-171,
//!     GetFilletShape L175-193, Add L197-230, radius queries L234-355,
//!     Simulate L402-435, NbSurf L439-451, Sect L455-471)
//!   - ChFi3d_ChBuilder.cxx (ctor L189-194, Add L201-262, SetDist L268-310,
//!     GetDist L312-324, Dists)
//!
//! OCCT C++ inheritance (ChFi3d_FilBuilder / ChFi3d_ChBuilder derive from
//! ChFi3d_Builder) is modeled by composition: the derived builders embed
//! the `ChFi3dBuilder` base struct as `base`.

use glam::DVec2;
use rcad_kernel::topo::topods::Orientation;
use rcad_kernel::topods::{self, Shape};

use super::chfi_ds::{
    ChFi3dFilletShape, ChFiDSSpineHandle, ChFiDS_State, ChFiDS_ChamfMode, ChFiDS_ErrorStatus,
    ChFiDSStripeMap, ChFiDSChamfSpine, ChFiDSFilSpine, ChFiDSStripe, ChFiDSMap, LawFunction,
    SharedStripe,
};
use crate::geomalgo::gtests_stubs::GeomAbsShape;

// =========================================================================
// OCCT TopOpeBRepDS_HDataStructure / TopOpeBRepBuild_HBuilder — pending
// TKTopAlgo/TKBool reconstruction subsystems.  Referenced by the builder as
// opaque handles until translated.
// =========================================================================

#[derive(Debug, Clone, Default)]
pub struct TopOpeBRepDSHDataStructure;

#[derive(Debug, Clone, Default)]
pub struct TopOpeBRepBuildHBuilder;

// =========================================================================
// OCCT ChFi3d_Builder member fields (ChFi3d_Builder.hxx L741-763, L841-846)
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFi3dBuilder {
    // OCCT: double tolappangle / tolesp / tol2d / tolapp3d / tolapp2d /
    // fleche; GeomAbs_Shape myConti
    pub tolappangle: f64,
    pub tolesp: f64,
    pub tol2d: f64,
    pub tolapp3d: f64,
    pub tolapp2d: f64,
    pub fleche: f64,
    pub my_conti: GeomAbsShape,
    // OCCT: ChFiDS_Map myEFMap / myESoMap / myEShMap / myVFMap / myVEMap
    pub my_ef_map: ChFiDSMap,
    pub my_eso_map: ChFiDSMap,
    pub my_esh_map: ChFiDSMap,
    pub my_vf_map: ChFiDSMap,
    pub my_ve_map: ChFiDSMap,
    // OCCT: occ::handle<TopOpeBRepDS_HDataStructure> myDS
    pub my_ds: Option<TopOpeBRepDSHDataStructure>,
    // OCCT: occ::handle<TopOpeBRepBuild_HBuilder> myCoup
    pub my_coup: Option<TopOpeBRepBuildHBuilder>,
    // OCCT: NCollection_List<occ::handle<ChFiDS_Stripe>> myListStripe
    pub my_list_stripe: Vec<SharedStripe>,
    // OCCT: ChFiDS_StripeMap myVDataMap
    pub my_vdata_map: ChFiDSStripeMap,
    // OCCT: NCollection_List<ChFiDS_Regul> myRegul
    pub my_regul: Vec<()>, // ChFiDS_Regul pending
    // OCCT: NCollection_List<occ::handle<ChFiDS_Stripe>> badstripes
    pub badstripes: Vec<SharedStripe>,
    // OCCT: NCollection_List<TopoDS_Shape> badvertices
    pub badvertices: Vec<Shape>,
    // OCCT: NCollection_DataMap<TopoDS_Shape, List<int>> myEVIMap
    pub my_evi_map: std::collections::HashMap<u64, Vec<i32>>,
    // OCCT: NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> myEdgeFirstFace
    pub my_edge_first_face: std::collections::HashMap<u64, Shape>,
    // OCCT: bool done / hasresult
    pub done: bool,
    pub hasresult: bool,
    // OCCT: TopoDS_Shape myShape (private section L841)
    pub my_shape: Shape,
    // OCCT: double angular
    pub angular: f64,
    // OCCT: NCollection_List<TopoDS_Shape> myGenerated (L843)
    pub my_generated: Vec<Shape>,
    // OCCT: TopoDS_Shape myShapeResult (L844)
    pub my_shape_result: Option<Shape>,
    // OCCT: TopoDS_Shape badShape (L845)
    pub bad_shape: Option<Shape>,
    /// The BRep the root shape belongs to (rcad architecture: TopoDS_Shape
    /// lives inside a BRep TShape table; OCCT has global handle graphs).
    pub my_brep: rcad_kernel::topods::BRep,
}

impl ChFi3dBuilder {
    /// OCCT ChFi3d_Builder_1.cxx L341-364.
    pub fn new(brep: &rcad_kernel::topods::BRep, s: Shape, ta: f64) -> Self {
        let mut b = ChFi3dBuilder {
            done: false,
            my_shape: s,
            my_brep: brep.clone(),
            tolappangle: 0.0,
            tolesp: 0.0,
            tol2d: 0.0,
            tolapp3d: 0.0,
            tolapp2d: 0.0,
            fleche: 0.0,
            my_conti: GeomAbsShape::C0,
            my_ef_map: ChFiDSMap::new(),
            my_eso_map: ChFiDSMap::new(),
            my_esh_map: ChFiDSMap::new(),
            my_vf_map: ChFiDSMap::new(),
            my_ve_map: ChFiDSMap::new(),
            my_ds: None,
            my_coup: None,
            my_list_stripe: Vec::new(),
            my_vdata_map: ChFiDSStripeMap::new(),
            my_regul: Vec::new(),
            badstripes: Vec::new(),
            badvertices: Vec::new(),
            my_evi_map: std::collections::HashMap::new(),
            my_edge_first_face: std::collections::HashMap::new(),
            hasresult: false,
            angular: 0.0,
            my_generated: Vec::new(),
            my_shape_result: None,
            bad_shape: None,
        };
        b.my_ds = Some(TopOpeBRepDSHDataStructure);
        b.my_coup = Some(TopOpeBRepBuildHBuilder);
        // myEFMap.Fill(S, TopAbs_EDGE, TopAbs_FACE);  (L354)
        b.my_ef_map.fill(&brep, topods::ShapeType::Edge, topods::ShapeType::Face);
        // myESoMap.Fill(S, TopAbs_EDGE, TopAbs_SOLID);  (L355)
        b.my_eso_map.fill(&brep, topods::ShapeType::Edge, topods::ShapeType::Solid);
        // myEShMap.Fill(S, TopAbs_EDGE, TopAbs_SHELL);  (L356)
        b.my_esh_map.fill(&brep, topods::ShapeType::Edge, topods::ShapeType::Shell);
        // myVFMap.Fill(S, TopAbs_VERTEX, TopAbs_FACE);  (L357)
        b.my_vf_map.fill(&brep, topods::ShapeType::Vertex, topods::ShapeType::Face);
        // myVEMap.Fill(S, TopAbs_VERTEX, TopAbs_EDGE);  (L358)
        b.my_ve_map.fill(&brep, topods::ShapeType::Vertex, topods::ShapeType::Edge);
        // SetParams(Ta, 1.0e-4, 1.e-5, 1.e-4, 1.e-5, 1.e-3);  (L359)
        b.set_params(ta, 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-5, 1.0e-3);
        // SetContinuity(GeomAbs_C1, Ta);  (L360)
        b.set_continuity(GeomAbsShape::C1, ta);
        b
    }

    /// OCCT ChFi3d_Builder_1.cxx L367-380.
    pub fn set_params(
        &mut self,
        tang: f64,
        tesp: f64,
        t2d: f64,
        tapp3d: f64,
        tolapp2d: f64,
        fleche: f64,
    ) {
        self.angular = tang;
        self.tolesp = tesp;
        self.tol2d = t2d;
        self.tolapp3d = tapp3d;
        self.tolapp2d = tolapp2d;
        self.fleche = fleche;
    }

    /// OCCT ChFi3d_Builder_1.cxx L382-389.
    pub fn set_continuity(&mut self, internal_continuity: GeomAbsShape, angular_tolerance: f64) {
        self.my_conti = internal_continuity;
        self.tolappangle = angular_tolerance;
    }

    /// OCCT ChFi3d_Builder_1.cxx L391-393.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT ChFi3d_Builder_1.cxx L393-398.
    pub fn shape(&self) -> Shape {
        assert!(self.done, "ChFi3d_Builder::Shape() - no result");
        self.my_shape_result.clone().expect("no result shape")
    }

    /// OCCT ChFi3d_Builder_1.cxx L1181-1199.
    pub fn remove(&mut self, e: &Shape) {
        let mut ic = None;
        for (i, stripe) in self.my_list_stripe.iter().enumerate() {
            let st = stripe.read().expect("stripe lock");
            let Some(sp) = st.spine() else {
                continue;
            };
            for j in 1..=sp.base().nb_edges() {
                if e.is_same(sp.base().edges(j)) {
                    ic = Some(i);
                    break;
                }
            }
            if ic.is_some() {
                break;
            }
        }
        if let Some(i) = ic {
            self.my_list_stripe.remove(i);
            return;
        }
    }

    /// OCCT ChFi3d_Builder_1.cxx L1201-1211 (Value — the stripe handle).
    pub fn value_stripe(&self, i: usize) -> SharedStripe {
        self.my_list_stripe[i - 1].clone()
    }

    /// OCCT ChFi3d_Builder_1.cxx L1201-1211 (Value — the spine handle).
    pub fn value(&self, i: usize) -> ChFiDSSpineHandle {
        self.my_list_stripe[i - 1]
            .read()
            .expect("stripe lock")
            .spine()
            .expect("null spine")
            .clone()
    }

    /// OCCT ChFi3d_Builder_1.cxx L1213-1229.
    pub fn nb_elements(&self) -> usize {
        let mut i = 0usize;
        for stripe in &self.my_list_stripe {
            let st = stripe.read().expect("stripe lock");
            match st.spine() {
                None => break,
                Some(_) => i += 1,
            }
        }
        i
    }

    /// OCCT ChFi3d_Builder_1.cxx L1231-1253.
    pub fn contains(&self, e: &Shape) -> usize {
        for (i, stripe) in self.my_list_stripe.iter().enumerate() {
            let st = stripe.read().expect("stripe lock");
            let Some(sp) = st.spine() else {
                break;
            };
            for j in 1..=sp.base().nb_edges() {
                if e.is_same(sp.base().edges(j)) {
                    return i + 1;
                }
            }
        }
        0
    }

    /// OCCT ChFi3d_Builder_1.cxx L1255-1280.
    pub fn contains_in_spine(&self, e: &Shape, index_in_spine: &mut usize) -> usize {
        *index_in_spine = 0;
        for (i, stripe) in self.my_list_stripe.iter().enumerate() {
            let st = stripe.read().expect("stripe lock");
            let Some(sp) = st.spine() else {
                break;
            };
            for j in 1..=sp.base().nb_edges() {
                if e.is_same(sp.base().edges(j)) {
                    *index_in_spine = j;
                    return i + 1;
                }
            }
        }
        0
    }

    /// OCCT ChFi3d_Builder_1.cxx L1282-1290.
    pub fn length(&self, ic: usize) -> f64 {
        if ic <= self.nb_elements() {
            let sp = self.value(ic);
            let n = sp.base().nb_edges();
            return sp.base().last_parameter_of(n);
        }
        -1.0
    }

    /// OCCT ChFi3d_Builder_1.cxx L1292-1300.
    pub fn first_vertex(&self, ic: usize) -> Shape {
        if ic <= self.nb_elements() {
            return self.value(ic).base().first_vertex();
        }
        Shape::null()
    }

    /// OCCT ChFi3d_Builder_1.cxx L1302-1310.
    pub fn last_vertex(&self, ic: usize) -> Shape {
        if ic <= self.nb_elements() {
            return self.value(ic).base().last_vertex();
        }
        Shape::null()
    }

    /// OCCT ChFi3d_Builder.cxx L178-675 (Compute).
    ///
    /// The stripe/corner numerical core (PerformSetOfSurf,
    /// PerformFilletOnVertex) and the DS reconstruction (ChFi3d_FilDS,
    /// TopOpeBRepBuild reconstruction) are pending translation; those calls
    /// follow the OCCT exception paths (catch -> badstripes/badvertices ->
    /// done = false) so the builder surfaces the same failure state.
    pub fn compute(&mut self) {
        // L223: UpdateTolesp();
        self.update_tolesp();

        // L225-228
        if self.my_list_stripe.is_empty() {
            panic!("Standard_Failure: There are no suitable edges for chamfer or fillet");
        }

        // L230-234
        self.reset();
        self.my_ds = Some(TopOpeBRepDSHDataStructure);
        self.done = true;
        self.hasresult = false;

        // L236-257: filling of myVDataMap
        for itel in self.my_list_stripe.clone() {
            let st = itel.read().expect("stripe lock");
            let sp = st.spine().expect("null spine").base().clone();
            drop(st);
            if sp.first_status() <= ChFiDS_State::BreakPoint {
                let stripe = itel.clone();
                self.my_vdata_map.add(&sp.first_vertex(), stripe);
            } else if sp.first_status() == ChFiDS_State::FreeBoundary {
                // OCCT L247: ExtentOneCorner(FirstVertex, stripe) — pending
                // ChFi3d_Builder_CnCrn translation.
            }
            if sp.last_status() <= ChFiDS_State::BreakPoint {
                let stripe = itel.clone();
                self.my_vdata_map.add(&sp.last_vertex(), stripe);
            } else if sp.last_status() == ChFiDS_State::FreeBoundary {
                // OCCT L255: ExtentOneCorner(LastVertex, stripe) — pending.
            }
        }
        // L259: preanalysis to evaluate the extensions (ExtentAnalyse).
        self.extent_analyse();

        // L266-293: Construction of the stripe of fillet on each stripe.
        for itel in self.my_list_stripe.clone() {
            {
                let mut st = itel.write().expect("stripe lock");
                if let Some(sp) = st.my_spine.as_mut() {
                    sp.base_mut().set_error_status(ChFiDS_ErrorStatus::Ok);
                }
            }
            // L273: PerformSetOfSurf(itel.ChangeValue()) — pending numerical
            // core (ChFi3d_Builder_2/6, BRepBlend_Walking).  OCCT wraps the
            // call in try/catch: the pending core raises, so the catch path
            // below is the 1:1 behavior.
            self.perform_set_of_surf_pending();
            // L281-282: badstripes.Append(itel.Value()); done = true;
            self.badstripes.push(itel.clone());
            self.done = true;
            // L283-286: if spine error is Ok, set it to ChFiDS_Error.
            {
                let mut st = itel.write().expect("stripe lock");
                if let Some(sp) = st.my_spine.as_mut() {
                    if sp.base().error_status() == ChFiDS_ErrorStatus::Ok {
                        sp.base_mut().set_error_status(ChFiDS_ErrorStatus::Error);
                    }
                }
            }
            // L288-292
            if !self.done {
                self.badstripes.push(itel.clone());
            }
            self.done = true;
        }
        // L294: done = (badstripes.IsEmpty());
        self.done = self.badstripes.is_empty();

        // L301-332: construct fillets on each vertex + feed the DS
        if self.done {
            for j in 1..=self.my_vdata_map.extent() {
                // L310: PerformFilletOnVertex(j) — pending corner machinery;
                // the OCCT catch path appends the vertex:
                self.perform_fillet_on_vertex_pending();
                self.badvertices.push(self.my_vdata_map.find_key(j).clone());
                self.hasresult = false;
                self.done = true;
                if !self.done {
                    self.badvertices.push(self.my_vdata_map.find_key(j).clone());
                }
                self.done = true;
            }
            // L328-331
            if !self.hasresult {
                self.done = self.badvertices.is_empty();
            }
        }

        // L339-354: solids/shells are registered in the DS (DStr.AddShape);
        // L354-396: stripe intersections (ChFi3d_StripeEdgeInter) +
        // ChFi3d_FilDS; L403-579: the TopOpeBRepBuild reconstruction
        // (myCoup->Perform/MergeSolid) and myShapeResult assembly; L574:
        // SetRegul.  All depend on the pending TopOpeBRepDS /
        // TopOpeBRepBuild subsystems and are unreachable with done=false.

        // L655-674: SameParameter pass over the new faces (only when done).
        if self.is_done() {
            // BRepLib::SameParameter / ShapeFix::SameParameter — pending.
        }
    }

    /// OCCT ChFi3d_Builder.cxx L924-946.
    pub fn reset(&mut self) {
        self.done = false;
        self.my_vdata_map.clear();
        self.my_regul.clear();
        self.my_evi_map.clear();
        self.badstripes.clear();
        self.badvertices.clear();

        let mut i = 0usize;
        while i < self.my_list_stripe.len() {
            let has_spine = {
                let st = self.my_list_stripe[i].read().expect("stripe lock");
                st.spine().is_some()
            };
            if has_spine {
                self.my_list_stripe[i].write().expect("stripe lock").reset();
                i += 1;
            } else {
                self.my_list_stripe.remove(i);
            }
        }
    }

    /// OCCT ChFi3d_Builder.cxx L950-977.
    pub fn generated(&mut self, eouv: &Shape) -> &Vec<Shape> {
        self.my_generated.clear();
        if eouv.is_null() {
            return &self.my_generated;
        }
        let st = eouv.shape_type();
        if st != topods::ShapeType::Edge && st != topods::ShapeType::Vertex {
            return &self.my_generated;
        }
        if let Some(l) = self.my_evi_map.get(&eouv.ptr_id()) {
            for i in l.clone() {
                // OCCT L968: myCoup->NewFaces(I) — pending reconstruction.
                let _ = i;
            }
        }
        &self.my_generated
    }

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    // Pending-subsystem boundary markers (see file header).
    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    fn update_tolesp(&mut self) {
        // OCCT ChFi3d_Builder_C2.cxx L814 — assigns to tolesp the minimal
        // spine tolesp; pending spine-parameter machinery.
    }

    fn extent_analyse(&mut self) {
        // OCCT ChFi3d_Builder.cxx L144-174 — depends on
        // ChFi3d_NumberOfSharpEdges and ExtentOne/Two/ThreeCorner (pending).
    }

    fn perform_set_of_surf_pending(&mut self) {
        // OCCT ChFi3d_Builder_2.cxx PerformSetOfSurf — pending numerical
        // core; translated callers surface the OCCT exception path.
    }

    fn perform_fillet_on_vertex_pending(&mut self) {
        // OCCT ChFi3d_Builder.cxx L759-920 — pending corner machinery.
    }
}

// =========================================================================
// OCCT ChFi3d_FilBuilder (ChFi3d_FilBuilder.hxx) — BlendFunc_Shape myShape
// field is renamed my_blend_shape to avoid the base myShape collision.
// =========================================================================

/// OCCT BlendFunc_Shape (TKGeomAlgo/BlendFunc): the fillet section shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendFuncShape {
    Rational,
    QuasiAngular,
    Polynomial,
}

#[derive(Debug, Clone)]
pub struct ChFi3dFilBuilder {
    pub base: ChFi3dBuilder,
    /// OCCT ChFi3d_FilBuilder.hxx: BlendFunc_Shape myShape.
    pub my_blend_shape: BlendFuncShape,
}

impl ChFi3dFilBuilder {
    /// OCCT ChFi3d_FilBuilder.cxx L147-153.
    pub fn new(
        brep: &rcad_kernel::topods::BRep,
        s: Shape,
        fshape: ChFi3dFilletShape,
        ta: f64,
    ) -> Self {
        let mut b = ChFi3dFilBuilder {
            base: ChFi3dBuilder::new(brep, s, ta),
            my_blend_shape: BlendFuncShape::Rational,
        };
        b.set_fillet_shape(fshape);
        b
    }

    /// OCCT ChFi3d_FilBuilder.cxx L157-171.
    pub fn set_fillet_shape(&mut self, fshape: ChFi3dFilletShape) {
        match fshape {
            ChFi3dFilletShape::Rational => self.my_blend_shape = BlendFuncShape::Rational,
            ChFi3dFilletShape::QuasiAngular => {
                self.my_blend_shape = BlendFuncShape::QuasiAngular
            }
            ChFi3dFilletShape::Polynomial => self.my_blend_shape = BlendFuncShape::Polynomial,
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L175-193.
    pub fn get_fillet_shape(&self) -> ChFi3dFilletShape {
        match self.my_blend_shape {
            BlendFuncShape::Rational => ChFi3dFilletShape::Rational,
            BlendFuncShape::QuasiAngular => ChFi3dFilletShape::QuasiAngular,
            BlendFuncShape::Polynomial => ChFi3dFilletShape::Polynomial,
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L197-218.
    pub fn add(&mut self, e: &Shape) {
        let dummy = Shape::null();

        if self.base.contains(e) == 0 && self.base.my_ef_map.contains(e) {
            let mut stripe = ChFiDSStripe::default();
            let mut sp = ChFiDSSpineHandle::Fil(ChFiDSFilSpine::with_tol(self.base.tolesp));

            let mut e_wnt = e.clone();
            e_wnt.orientation = Orientation::Forward;
            let added = {
                {
                    let Some(fsp) = sp.down_cast_fil_mut() else {
                        return;
                    };
                    fsp.base.set_edges(e_wnt);
                }
                if self.perform_element(&sp, -1.0, &dummy) {
                    self.perform_extremity_pending();
                    let Some(fsp) = sp.down_cast_fil_mut() else {
                        return;
                    };
                    fsp.base.load();
                    true
                } else {
                    false
                }
            };
            if added {
                stripe.change_spine(sp);
                self.base
                    .my_list_stripe
                    .push(std::sync::Arc::new(std::sync::RwLock::new(stripe)));
            }
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L222-230.
    pub fn add_radius(&mut self, radius: f64, e: &Shape) {
        self.add(e);
        let ic = self.base.contains(e);
        if ic > 0 {
            self.set_radius_on_edge(radius, ic, e);
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L234-241.
    pub fn set_radius_law(&mut self, c: LawFunction, ic: usize, iinc: usize) {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value_stripe(ic);
            let mut st = sp.write().expect("stripe lock");
            if let Some(ChFiDSSpineHandle::Fil(_fsp)) = st.my_spine.as_mut() {
                // OCCT: fsp->SetRadius(C, IinC) — Law_Function storage
                // pending TKMath Law package translation.
                let _ = (&c, iinc);
            }
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L245-253.
    pub fn is_constant(&self, ic: usize) -> bool {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value(ic);
            if let Some(fsp) = sp.down_cast_fil() {
                return fsp.is_constant();
            }
        }
        false
    }

    /// OCCT ChFi3d_FilBuilder.cxx L257-265.
    pub fn radius(&self, ic: usize) -> f64 {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value(ic);
            if let Some(fsp) = sp.down_cast_fil() {
                return fsp.radius();
            }
        }
        -1.0
    }

    /// OCCT ChFi3d_FilBuilder.cxx L269-276.
    pub fn reset_contour(&mut self, ic: usize) {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value_stripe(ic);
            let mut st = sp.write().expect("stripe lock");
            if let Some(fsp) = st.my_spine.as_mut().and_then(|s| s.down_cast_fil_mut()) {
                fsp.reset(true);
            }
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L280-287.
    pub fn set_radius_on_edge(&mut self, radius: f64, ic: usize, e: &Shape) {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value_stripe(ic);
            let mut st = sp.write().expect("stripe lock");
            if let Some(fsp) = st.my_spine.as_mut().and_then(|s| s.down_cast_fil_mut()) {
                fsp.set_radius_on_edge(radius, e);
            }
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L291-298.
    pub fn unset_on_edge(&mut self, ic: usize, e: &Shape) {
        if ic <= self.base.nb_elements() {
            let _ = (ic, e);
            // OCCT: fsp->UnSetRadius(E) — pending parandrad edge-parameter
            // mapping (FirstParameter(IE) over unfilled abscissa).
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L324-331.
    pub fn set_radius_uandr(&mut self, uandr: DVec2, ic: usize, iinc: usize) {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value_stripe(ic);
            let mut st = sp.write().expect("stripe lock");
            if let Some(fsp) = st.my_spine.as_mut().and_then(|s| s.down_cast_fil_mut()) {
                fsp.set_radius_uandr(uandr, iinc);
            }
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L402-435 (Simulate) — the stripe walk is
    /// real; PerformSetOfSurf(simul=true) is pending.
    pub fn simulate(&mut self, ic: usize) {
        for (i, stripe) in self.base.my_list_stripe.iter().enumerate() {
            if i + 1 == ic {
                // OCCT: PerformSetOfSurf(itel.ChangeValue(), true) — pending.
                let _ = stripe;
                break;
            }
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L439-451.
    pub fn nb_surf(&self, ic: usize) -> usize {
        for (i, stripe) in self.base.my_list_stripe.iter().enumerate() {
            if i + 1 == ic {
                let st = stripe.read().expect("stripe lock");
                return st.my_hdata.len();
            }
        }
        0
    }

    // Pending numerical-core boundary (OCCT line markers in the base).
    fn perform_element(&self, _spine: &ChFiDSSpineHandle, _offset: f64, _g: &Shape) -> bool {
        // OCCT ChFi3d_Builder_1.cxx L887 — PerformElement walks the
        // tangency-connected edge chain (ChFi3d_SameSide, FaceTangency,
        // TangentOnVertex).  Pending translation.
        false
    }

    fn perform_extremity_pending(&mut self) {
        // OCCT ChFi3d_Builder_1.cxx L714 — PerformExtremity pending.
    }
}

// =========================================================================
// OCCT ChFi3d_ChBuilder (ChFi3d_ChBuilder.cxx L189-194 ctor, L201-262 Add,
// L268-310 SetDist, L312-324 GetDist, Dists).
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFi3dChBuilder {
    pub base: ChFi3dBuilder,
    /// OCCT ChFi3d_ChBuilder.hxx: ChFiDS_ChamfMode myMode.
    pub my_mode: ChFiDS_ChamfMode,
}

impl ChFi3dChBuilder {
    /// OCCT ChFi3d_ChBuilder.cxx L189-194.
    pub fn new(brep: &rcad_kernel::topods::BRep, s: Shape, ta: f64) -> Self {
        ChFi3dChBuilder {
            base: ChFi3dBuilder::new(brep, s, ta),
            my_mode: ChFiDS_ChamfMode::ClassicChamfer,
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L201-230.
    pub fn add(&mut self, e: &Shape) {
        let dummy = Shape::null();

        if self.base.contains(e) == 0 && self.base.my_ef_map.contains(e) {
            let mut stripe = ChFiDSStripe::default();
            let mut sp = ChFiDSSpineHandle::Chamf(ChFiDSChamfSpine::with_tol(self.base.tolesp));

            let mut e_wnt = e.clone();
            e_wnt.orientation = Orientation::Forward;
            let added = {
                {
                    let Some(csp) = sp.down_cast_chamf_mut() else {
                        return;
                    };
                    csp.base.set_edges(e_wnt);
                }
                if self.perform_element(&sp, -1.0, &dummy) {
                    self.perform_extremity_pending();
                    let Some(csp) = sp.down_cast_chamf_mut() else {
                        return;
                    };
                    csp.base.load();
                    true
                } else {
                    false
                }
            };
            if added {
                stripe.change_spine(sp);
                self.base
                    .my_list_stripe
                    .push(std::sync::Arc::new(std::sync::RwLock::new(stripe)));
            }
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L232-262.
    pub fn add_dist(&mut self, dis: f64, e: &Shape) {
        if self.base.contains(e) == 0 && self.base.my_ef_map.contains(e) {
            let dummy = Shape::null();

            let mut stripe = ChFiDSStripe::default();
            let mut sp = ChFiDSSpineHandle::Chamf(ChFiDSChamfSpine::with_tol(self.base.tolesp));

            let mut e_wnt = e.clone();
            e_wnt.orientation = Orientation::Forward;

            let added = {
                {
                    let Some(csp) = sp.down_cast_chamf_mut() else {
                        return;
                    };

                    csp.set_mode(self.my_mode);

                    csp.base.set_edges(e_wnt);
                }
                if self.perform_element(&sp, -1.0, &dummy) {
                    {
                        let Some(csp) = sp.down_cast_chamf_mut() else {
                            return;
                        };
                        csp.base.load();
                        csp.set_dist(dis);
                    }

                    self.perform_extremity_pending();
                    true
                } else {
                    false
                }
            };
            if added {
                stripe.change_spine(sp);
                self.base
                    .my_list_stripe
                    .push(std::sync::Arc::new(std::sync::RwLock::new(stripe)));
            }
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L268-310.
    pub fn set_dist(&mut self, dis: f64, ic: usize, f: &Shape) {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value_stripe(ic);
            let mut st = sp.write().expect("stripe lock");
            let Some(csp) = st.my_spine.as_mut().and_then(|s| s.down_cast_chamf_mut()) else {
                return;
            };

            // Search the first edge which has a common face equal to F
            let mut i = 1usize;
            let mut found = false;
            while i <= csp.base.nb_edges() && !found {
                let (f1, f2) = search_common_faces(&self.base.my_ef_map, csp.base.edges(i));
                found = f1.is_same(f) || f2.is_same(f);
                i += 1;
            }

            if found {
                csp.set_dist(dis);
            } else {
                panic!(
                    "Standard_DomainError: the face is not common to any of edges of the contour"
                );
            }
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L312-324.
    pub fn get_dist(&self, ic: usize) -> f64 {
        let sp = self.base.value(ic);
        match sp.down_cast_chamf() {
            Some(csp) => csp.get_dist(),
            None => panic!("Standard_DomainError: not a chamfer spine"),
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx (Dists).
    pub fn dists(&self, ic: usize) -> (f64, f64) {
        let sp = self.base.value(ic);
        match sp.down_cast_chamf() {
            Some(csp) => csp.dists(),
            None => panic!("Standard_DomainError: not a chamfer spine"),
        }
    }

    fn perform_element(&self, _spine: &ChFiDSSpineHandle, _offset: f64, _g: &Shape) -> bool {
        // OCCT ChFi3d_Builder_1.cxx L887 — pending translation.
        false
    }

    fn perform_extremity_pending(&mut self) {
        // OCCT ChFi3d_Builder_1.cxx L714 — pending translation.
    }
}

/// OCCT ChFi3d_Builder_0.cxx — SearchCommonFaces(EFMap, E, F1, F2).
fn search_common_faces(efmap: &ChFiDSMap, e: &Shape) -> (Shape, Shape) {
    let list = efmap.find(e);
    let f1 = list.first().cloned().unwrap_or_else(Shape::null);
    let f2 = list.get(1).cloned().unwrap_or_else(Shape::null);
    (f1, f2)
}
