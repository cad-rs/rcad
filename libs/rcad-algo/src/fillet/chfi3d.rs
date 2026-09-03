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
use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::{Orientation, TShape};
use rcad_kernel::topo::topods::BRepTool as _;

use glam::DVec3;
use rcad_kernel::topods::{self, Shape};

use super::chfi_ds::{
    ChFi3dFilletShape, ChFiDSSpineHandle, ChFiDS_State, ChFiDS_ChamfMethod, ChFiDS_ChamfMode,
    ChFiDS_ErrorStatus,
    ChFiDSStripeMap, ChFiDSChamfSpine, ChFiDSFilSpine, ChFiDSStripe, ChFiDSMap, LawFunction, ChFiDSSurfData,
    SharedStripe,
};

// =========================================================================
// OCCT TopOpeBRepDS_HDataStructure / TopOpeBRepBuild_HBuilder — the DS
// types live in the topopebrepds module (1:1 translation of the
// TKBool/TopOpeBRepDS subset the builder uses).  The TopOpeBRepBuild
// reconstruction subsystem is pending; referenced by the builder as an
// opaque handle until translated.
// =========================================================================

pub use super::topopebrepds::{
    InterferenceRef, TopOpeBRepDSCurvePointInterference, TopOpeBRepDSCurve, TopOpeBRepDSHDataStructure,
    TopOpeBRepDSInterference, TopOpeBRepDSKind, TopOpeBRepDSPoint, TopOpeBRepDSSolidSurfaceInterference,
    TopOpeBRepDSSurface, TopOpeBRepDSSurfaceCurveInterference, TopOpeBRepDSTransition,
};

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
        b.my_ds = Some(TopOpeBRepDSHDataStructure::default());
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

    /// OCCT ChFi3d_Builder_1.cxx L1316-1328.
    pub fn abscissa(&self, ic: usize, v: &Shape) -> f64 {
        if ic <= self.nb_elements() {
            return self.value(ic).base().absc_of_vertex(v);
        }
        -1.0
    }

    /// OCCT ChFi3d_Builder_1.cxx L1330-1341.
    pub fn relative_abscissa(&self, ic: usize, v: &Shape) -> f64 {
        if ic <= self.nb_elements() {
            return self.abscissa(ic, v) / self.length(ic);
        }
        -1.0
    }

    /// OCCT ChFi3d_Builder_1.cxx L1343-1353.
    pub fn closed(&self, ic: usize) -> bool {
        if ic <= self.nb_elements() {
            return self.value(ic).base().is_closed();
        }
        false
    }

    /// OCCT ChFi3d_Builder_1.cxx L1355-1365.
    pub fn closed_and_tangent(&self, ic: usize) -> bool {
        if ic <= self.nb_elements() {
            return self.value(ic).base().is_periodic();
        }
        false
    }

    /// OCCT ChFi3d_Builder_1.cxx L398-401.
    pub fn nb_faulty_contours(&self) -> usize {
        self.badstripes.len()
    }

    /// OCCT ChFi3d_Builder_1.cxx L403-428 — the index (in myListStripe) of
    /// the I-th faulty stripe; handle identity compares Arc pointers.
    pub fn faulty_contour(&self, i: usize) -> usize {
        let st = match self.badstripes.get(i.wrapping_sub(1)) {
            Some(s) => s,
            None => return 0,
        };
        for (k, stripe) in self.my_list_stripe.iter().enumerate() {
            if std::sync::Arc::ptr_eq(st, stripe) {
                return k + 1;
            }
        }
        0
    }

    /// OCCT ChFi3d_Builder_1.cxx L430-455.
    pub fn nb_computed_surfaces(&self, ic: usize) -> usize {
        let st = match self.my_list_stripe.get(ic.wrapping_sub(1)) {
            Some(s) => s,
            None => return 0,
        };
        let guard = st.read().expect("stripe lock");
        match guard.spine() {
            None => return 0,
            Some(_) => {}
        }
        guard.my_hdata.len()
    }

    /// OCCT ChFi3d_Builder_1.cxx L457-474 — depends on the pending
    /// TopOpeBRepDS surface table (myDS->Surface(isurf)).
    pub fn computed_surface(&self, _ic: usize, _is: usize) {
        // Pending TopOpeBRepDS_HDataStructure::Surface translation.
    }

    /// OCCT ChFi3d_Builder_1.cxx L476-479.
    pub fn nb_faulty_vertices(&self) -> usize {
        self.badvertices.len()
    }

    /// OCCT ChFi3d_Builder_1.cxx L481-497.
    pub fn faulty_vertex(&self, iv: usize) -> Shape {
        match self.badvertices.get(iv.wrapping_sub(1)) {
            Some(v) => v.clone(),
            None => Shape::null(),
        }
    }

    /// OCCT ChFi3d_Builder_1.cxx L499-511.
    pub fn has_result(&self) -> bool {
        self.hasresult
    }

    /// OCCT ChFi3d_Builder_1.cxx L512-518.
    pub fn bad_shape(&self) -> Shape {
        assert!(self.hasresult, "ChFi3d_Builder::BadShape() - no result");
        self.bad_shape.clone().expect("no bad shape")
    }

    /// OCCT ChFi3d_Builder_1.cxx L520-545 (StripeStatus of the IC-th
    /// contour's spine).
    pub fn stripe_status(&self, ic: usize) -> ChFiDS_ErrorStatus {
        let st = self.value_stripe(ic);
        let guard = st.read().expect("stripe lock");
        match guard.spine() {
            Some(sp) => sp.base().error_status(),
            None => ChFiDS_ErrorStatus::Error,
        }
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
        self.my_ds = Some(TopOpeBRepDSHDataStructure::default());
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
            // L273: PerformSetOfSurf(itel.ChangeValue()) — the OCCT try/catch
            // maps to the success flag; the catch path records a bad stripe.
            let surf_ok = self.perform_set_of_surf(&itel, false);
            if !surf_ok {
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
                // L310: PerformFilletOnVertex(j) — the OCCT try/catch maps
                // to the success flag; the catch path appends the vertex:
                let ok = self.perform_fillet_on_vertex(j);
                if !ok {
                    self.badvertices.push(self.my_vdata_map.find_key(j).clone());
                    self.hasresult = false;
                    self.done = true;
                }
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
        // SetRegul.  The reconstruction subsystem is pending translation;
        // when the OCCT flow enters `if (done)` the pending raise surfaces
        // as done = false with no result shape.
        if self.done {
            self.done = false; // pending: ChFi3d_FilDS + TopOpeBRepBuild reconstruction
        }

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
                if self.base.perform_element(&mut sp, -1.0, &dummy) {
                    self.base.perform_extremity(&mut sp);
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

    /// OCCT ChFi3d_FilBuilder.cxx L302-309.
    pub fn set_radius_at_vertex(&mut self, radius: f64, ic: usize, v: &Shape) {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value_stripe(ic);
            let mut st = sp.write().expect("stripe lock");
            if let Some(fsp) = st.my_spine.as_mut().and_then(|s| s.down_cast_fil_mut()) {
                fsp.set_radius_at_vertex(radius, v);
            }
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L313-320 — fsp->UnSetRadius(V); the
    /// UnSetRadius body needs Absc(V) (pending BRepAdaptor_Curve).
    pub fn unset_at_vertex(&mut self, ic: usize, v: &Shape) {
        if ic <= self.base.nb_elements() {
            let _ = (ic, v);
            // Pending: fsp->UnSetRadius(V).
        }
    }

    /// OCCT ChFi3d_FilBuilder.cxx L335-343.
    pub fn is_constant_on_edge(&self, ic: usize, e: &Shape) -> bool {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value(ic);
            if let Some(fsp) = sp.down_cast_fil() {
                return fsp.is_constant_on(fsp.base.index_of_edge(e));
            }
        }
        false
    }

    /// OCCT ChFi3d_FilBuilder.cxx L347-355.
    pub fn radius_on_edge(&self, ic: usize, e: &Shape) -> f64 {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value(ic);
            if let Some(fsp) = sp.down_cast_fil() {
                return fsp.radius_on_edge(e);
            }
        }
        -1.0
    }

    /// OCCT ChFi3d_FilBuilder.cxx L359-372 — GetBounds via ChangeLaw(E);
    /// the Law_Function Bounds query is pending the Law package.
    pub fn get_bounds(&self, _ic: usize, _e: &Shape, _f: &mut f64, _l: &mut f64) -> bool {
        false
    }

    /// OCCT ChFi3d_FilBuilder.cxx L376-384 — pending Law package.
    pub fn get_law(&self, _ic: usize, _e: &Shape) -> Option<LawFunction> {
        None
    }

    /// OCCT ChFi3d_FilBuilder.cxx L388-398 — pending Law package.
    pub fn set_law(&mut self, _ic: usize, _e: &Shape, _l: LawFunction) {}

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
                if self.base.perform_element(&mut sp, -1.0, &dummy) {
                    self.base.perform_extremity(&mut sp);
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
                if self.base.perform_element(&mut sp, -1.0, &dummy) {
                    {
                        let Some(csp) = sp.down_cast_chamf_mut() else {
                            return;
                        };
                        csp.base.load();
                        csp.set_dist(dis);
                    }

                    self.base.perform_extremity(&mut sp);
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

    /// OCCT ChFi3d_ChBuilder.cxx L368-441 (SetDists).
    ///
    /// The face-search loop is translated 1:1; the ConcaveSide choice
    /// (ChFi3d::ConcaveSide over BRepAdaptor_Surface) is a pending
    /// boundary, so the symmetric/asymmetric branch decision keeps the
    /// OCCT structure and falls to the SetDists(Dis1, Dis2) arm.
    pub fn set_dists(&mut self, dis1: f64, dis2: f64, ic: usize, f: &Shape) {
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
                // OCCT L404-412: Sb1/Sb2 Initialize + ConcaveSide choice —
                // pending ChFi3d::ConcaveSide translation; the branch it
                // selects is csp->SetDists(Dis2, Dis1) vs SetDists(Dis1, Dis2).
                csp.set_dists(dis1, dis2);
            } else {
                panic!(
                    "Standard_DomainError: the face is not common to any of edges of the contour"
                );
            }
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L445-477 (AddDA).
    pub fn add_da(&mut self, dis1: f64, angle: f64, e: &Shape, f: &Shape) {
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
                if self.base.perform_element(&mut sp, -1.0, f) {
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
                let shared = std::sync::Arc::new(std::sync::RwLock::new(stripe));
                // OCCT: myListStripe.Append(Stripe);
                self.base.my_list_stripe.push(shared.clone());
                {
                    let mut st = shared.write().expect("stripe lock");
                    let Some(spine_ref) = st.my_spine.as_mut() else {
                        return;
                    };
                    // OCCT: Spine->SetDistAngle(Dis1, Angle);
                    if let Some(csp) = spine_ref.down_cast_chamf_mut() {
                        csp.set_dist_angle(dis1, angle);
                    }
                    // OCCT: PerformExtremity(Spine);
                    self.base.perform_extremity(spine_ref);
                }
            }
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L480-528 (SetDistAngle).
    pub fn set_dist_angle(&mut self, dis: f64, angle: f64, ic: usize, f: &Shape) {
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
                csp.set_dist_angle(dis, angle);
            } else {
                panic!("Standard_DomainError: the face is not common to any edges of the contour");
            }
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L530-537 (GetDistAngle).
    pub fn get_dist_angle(&self, ic: usize) -> (f64, f64) {
        let sp = self.base.value(ic);
        match sp.down_cast_chamf() {
            Some(csp) => csp.get_dist_angle(),
            None => panic!("Standard_DomainError: not a chamfer spine"),
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L539-544 (SetMode).
    pub fn set_mode(&mut self, the_mode: ChFiDS_ChamfMode) {
        self.my_mode = the_mode;
    }

    /// OCCT ChFi3d_ChBuilder.cxx L546-555 (IsChamfer).
    pub fn is_chamfer(&self, ic: usize) -> ChFiDS_ChamfMethod {
        let sp = self.base.value(ic);
        match sp.down_cast_chamf() {
            Some(csp) => csp.is_chamfer(),
            None => panic!("Standard_DomainError: not a chamfer spine"),
        }
    }

    /// OCCT ChFi3d_ChBuilder.cxx L557-562 (Mode).
    pub fn mode(&self) -> ChFiDS_ChamfMode {
        self.my_mode
    }

    /// OCCT ChFi3d_ChBuilder.cxx L564-570 (ResetContour).
    pub fn reset_contour(&mut self, ic: usize) {
        if ic <= self.base.nb_elements() {
            let sp = self.base.value_stripe(ic);
            let mut st = sp.write().expect("stripe lock");
            if let Some(csp) = st.my_spine.as_mut().and_then(|s| s.down_cast_chamf_mut()) {
                csp.base.reset(true);
            }
        }
    }
}

/// OCCT ChFi3d_Builder_0.cxx — SearchCommonFaces(EFMap, E, F1, F2).
fn search_common_faces(efmap: &ChFiDSMap, e: &Shape) -> (Shape, Shape) {
    let list = efmap.find(e);
    let f1 = list.first().cloned().unwrap_or_else(Shape::null);
    let f2 = list.get(1).cloned().unwrap_or_else(Shape::null);
    (f1, f2)
}

// =========================================================================
// OCCT ChFi3d.cxx — namespace ChFi3d functions (package-level free
// functions).  Architecture note: rcad passes the owning BRep where OCCT
// reads geometry from TopoDS handles.
// =========================================================================

use super::chfi_ds::ChFiDS_TypeOfConcavity;
use crate::geomalgo::gtests_stubs::GeomAbsShape;
use super::chfi_kpart;
use super::chfi3d_builder_0::{
    brep_tools_ori_edge_in_face, intermediate_point, topopebreptool_nt,
};

/// OCCT ChFi3d.cxx L650-695 — Correct2dPoint: nudge a UV point inside the
/// surface natural bounds (analytic surfaces return unchanged).
pub fn correct_2d_point(f: &Shape, p2d: &mut glam::DVec2) {
    let fd = match f.as_face() {
        Some(fd) => fd,
        None => return,
    };
    let Some(surf) = fd.surface.clone() else {
        return;
    };
    // OCCT: analytic surfaces (type < Bezier) return unchanged; rcad
    // mirrors this by skipping the offset / BSpline families (pending:
    // those families are not translated yet).
    let coeff = 0.01;
    let [u1, u2, v1, v2] = surf.default_domain();
    if u1.is_finite() && u2.is_finite() {
        let eps = (coeff * (u2 - u1)).max(super::chfi3d_builder_0::P_CONFUSION);
        if (p2d.x - u1).abs() < eps {
            p2d.x = u1 + eps;
        }
        if (p2d.x - u2).abs() < eps {
            p2d.x = u2 - eps;
        }
    }
    if v1.is_finite() && v2.is_finite() {
        let eps = (coeff * (v2 - v1)).max(super::chfi3d_builder_0::P_CONFUSION);
        if (p2d.y - v1).abs() < eps {
            p2d.y = v1 + eps;
        }
        if (p2d.y - v2).abs() < eps {
            p2d.y = v2 - eps;
        }
    }
}

/// OCCT ChFi3d.cxx L40-156 — ChFi3d::DefineConnectType.
pub fn define_connect_type(
    brep: &topods::BRep,
    e: &Shape,
    f1: &Shape,
    f2: &Shape,
    sin_tol: f64,
    correct_point: bool,
) -> ChFiDS_TypeOfConcavity {
    use rcad_kernel::geom::CurveEval as _;

    let c1 = brep.curve_on_surface(e, f1);
    // For the case of seam edge
    let mut ee = e.clone();
    if f1.is_same(f2) {
        ee.orientation = if ee.orientation == Orientation::Forward {
            Orientation::Reversed
        } else {
            Orientation::Forward
        };
    }
    let c2 = brep.curve_on_surface(&ee, f2);
    let (Some((c1c, _f1r, _l1r)), Some((c2c, _f2r, _l2r))) = (c1, c2) else {
        return ChFiDS_TypeOfConcavity::Other;
    };

    let ed = e.as_edge().expect("not an edge");
    let Some(curve) = ed.curve.clone() else {
        return ChFiDS_TypeOfConcavity::Other;
    };
    let f = ed.range[0];
    let l = ed.range[1];
    let mut par_on_c = 0.5 * (f + l);
    let mut t1 = curve.derivative_at(par_on_c);
    if t1.length_squared() <= f64::MIN_POSITIVE {
        par_on_c = intermediate_point(f, l);
        t1 = curve.derivative_at(par_on_c);
    }
    if t1.length_squared() > f64::MIN_POSITIVE {
        t1 = t1.normalize();
    }

    if brep_tools_ori_edge_in_face(brep, e, f1) == Orientation::Reversed {
        t1 = -t1;
    }
    if f1.orientation == Orientation::Reversed {
        t1 = -t1;
    }

    let mut p = c1c.point_at(par_on_c);

    if correct_point {
        correct_2d_point(f1, &mut p);
    }
    let Some(surf1) = f1.as_face().and_then(|fd| fd.surface.clone()) else {
        return ChFiDS_TypeOfConcavity::Other;
    };
    let (_p3, d1u, d1v) = surf1.derivatives(p.x, p.y);
    let mut dn1 = d1u.cross(d1v);
    if f1.orientation == Orientation::Reversed {
        dn1 = -dn1;
    }

    let mut p = c2c.point_at(par_on_c);
    if correct_point {
        correct_2d_point(f2, &mut p);
    }
    let Some(surf2) = f2.as_face().and_then(|fd| fd.surface.clone()) else {
        return ChFiDS_TypeOfConcavity::Other;
    };
    let (_p3, d1u, d1v) = surf2.derivatives(p.x, p.y);
    let mut dn2 = d1u.cross(d1v);
    if f2.orientation == Orientation::Reversed {
        dn2 = -dn2;
    }

    dn1 = dn1.normalize();
    dn2 = dn2.normalize();

    let pro_vec = dn1.cross(dn2);
    let norm_pro_vec = pro_vec.length();
    if norm_pro_vec < sin_tol {
        // plane
        if dn1.dot(dn2) > 0.0 {
            // Tangent
            ChFiDS_TypeOfConcavity::Tangential
        } else {
            // Mixed not finished!
            ChFiDS_TypeOfConcavity::Convex
        }
    } else {
        let pro_vec = if norm_pro_vec > f64::MIN_POSITIVE {
            pro_vec / norm_pro_vec
        } else {
            pro_vec
        };
        let prod = t1.dot(pro_vec);
        if prod > 0.0 {
            ChFiDS_TypeOfConcavity::Convex
        } else {
            // reenters
            ChFiDS_TypeOfConcavity::Concave
        }
    }
}

/// OCCT LocalAnalysis_SurfaceContinuity — pending TKGeomAlgo translation.
/// The OCCT IsTangentFaces sample loop instantiates it per sample; the
/// pending state reports IsDone() == false for every sample, which per the
/// OCCT flow (nbNotDone == aNbSamples) makes the check fail.
struct LocalAnalysisSurfaceContinuity;

impl LocalAnalysisSurfaceContinuity {
    fn new() -> Self {
        LocalAnalysisSurfaceContinuity
    }
    fn is_done(&self) -> bool {
        false
    }
}

/// OCCT ChFi3d.cxx L160-295 — ChFi3d::IsTangentFaces.
///
/// Pending boundaries inside a 1:1 control flow:
///   - BRep_Tool::Continuity records (L165): rcad stores none, so the
///     early accept is never taken (equivalent to C0);
///   - LocalAnalysis_SurfaceContinuity (L255): pending, see above.
/// The final mid-point normal comparison (L287-294) is translated and is
/// reached once the continuity machinery exists.
pub fn is_tangent_faces(
    brep: &topods::BRep,
    the_edge: &Shape,
    the_face1: &Shape,
    the_face2: &Shape,
    the_order: GeomAbsShape,
) -> bool {
    use rcad_kernel::geom::CurveEval as _;

    let ed = the_edge.as_edge().expect("not an edge");
    let tol_c0 = 0.001f64.max(1.5 * ed.tolerance);

    let a_c2d1;
    let mut a_c2d2 = None;
    let (a_first, a_last);

    let face_fwd = |f: &Shape| {
        let mut ff = f.clone();
        ff.orientation = Orientation::Forward;
        ff
    };
    let reverse = |s: &Shape| {
        let mut x = s.clone();
        x.orientation = if x.orientation == Orientation::Forward {
            Orientation::Reversed
        } else {
            Orientation::Forward
        };
        x
    };

    if !the_face1.is_same(the_face2) {
        // OCCT L177-204: seam edge closed on both faces.
        let a_face1 = face_fwd(the_face1);
        // Find the edge in the face 1 with its in-face orientation.
        let mut an_edge_in_face1 = Shape::null();
        if let Some(fd) = a_face1.as_face() {
            if let Some(ts) = brep.tshapes.get(fd.outer_wire.index) {
                if let TShape::Wire(wd) = ts.as_ref() {
                    for we in &wd.edges {
                        if we.is_same(the_edge) {
                            // Wire entries may carry no BRep index (ptr-only); rebuild from the
                            // seed edge keeping the in-face orientation.
                            an_edge_in_face1 = Shape::from_parts(
                                the_edge.data.clone(),
                                the_edge.index,
                                we.location,
                                we.orientation,
                            );
                            break;
                        }
                    }
                }
            }
        }
        if an_edge_in_face1.is_null() {
            return false;
        }
        match brep.curve_on_surface(&an_edge_in_face1, &a_face1) {
            Some(v) => {
                a_c2d1 = Some(v.0);
                a_first = v.1;
                a_last = v.2;
            }
            None => return false,
        }
        let rev = reverse(&an_edge_in_face1);
        let a_face2 = face_fwd(the_face2);
        match brep.curve_on_surface(&rev, &a_face2) {
            Some(v) => a_c2d2 = Some(v.0),
            None => return false,
        }
    } else {
        // Obtaining of pcurves of edge on two faces.
        match brep.curve_on_surface(the_edge, the_face1) {
            Some(v) => {
                a_c2d1 = Some(v.0);
                a_first = v.1;
                a_last = v.2;
            }
            None => return false,
        }
        // For the case of seam edge
        let ee = reverse(the_edge);
        match brep.curve_on_surface(&ee, the_face2) {
            Some(v) => a_c2d2 = Some(v.0),
            None => return false,
        }
    }
    let (Some(a_c2d1), Some(a_c2d2)) = (a_c2d1, a_c2d2) else {
        return false;
    };

    // Obtaining of two surfaces from adjacent faces.
    let Some(_a_surf1) = the_face1.as_face().and_then(|fd| fd.surface.clone()) else {
        return false;
    };
    let Some(_a_surf2) = the_face2.as_face().and_then(|fd| fd.surface.clone()) else {
        return false;
    };

    // Computation of the number of samples on the edge (OCCT uses
    // Adaptor3d_TopolTool::NbSamples; a fixed default stands in — pending
    // BRepTopAdaptor_TopolTool translation).
    let a_nb_samples = 23usize;

    // Computation of the continuity.
    let mut nb_not_done = 0usize;
    let a_delta = (a_last - a_first) / (a_nb_samples as f64 - 1.0);
    for i in 1..=a_nb_samples {
        let a_par = if i == a_nb_samples {
            a_last
        } else {
            a_first + ((i - 1) as f64) * a_delta
        };

        let _ = a_par;
        let a_cont = LocalAnalysisSurfaceContinuity::new();
        if !a_cont.is_done() {
            if the_order == GeomAbsShape::C2 {
                continue; // NullSecondDerivative case — pending
            }
            nb_not_done += 1;
            continue;
        }
        // OCCT IsG1/IsG2 checks — unreachable while the continuity is
        // pending.
    }

    if nb_not_done == a_nb_samples {
        return false;
    }

    // Compare normals of tangent faces in the middle point.
    let mid_par = 0.5 * (a_first + a_last);
    let uv1 = a_c2d1.point_at(mid_par);
    let uv2 = a_c2d2.point_at(mid_par);
    let Some(normal1) = topopebreptool_nt(brep, uv1, the_face1) else {
        return false;
    };
    let Some(normal2) = topopebreptool_nt(brep, uv2, the_face2) else {
        return false;
    };
    let _ = tol_c0;
    normal1.dot(normal2) >= 0.0
}

/// OCCT ChFi3d.cxx L302-530 — ChFi3d::ConcaveSide.  Returns the Choix code
/// (0 = edge not found in a face, 10 = locally fake tangent).
pub fn concave_side(
    brep: &topods::BRep,
    f1: &Shape,
    f2: &Shape,
    e: &Shape,
    or1: &mut Orientation,
    or2: &mut Orientation,
) -> i32 {
    use rcad_kernel::geom::CurveEval as _;
    *or1 = Orientation::Forward;
    *or2 = Orientation::Forward;
    let ed = e.as_edge().expect("not an edge");
    let Some(curve) = ed.curve.clone() else {
        return 0;
    };
    let first = ed.range[0];
    let last = ed.range[1];
    let par = 0.691254 * first + 0.308746 * last;

    // OCCT: CE.D1(par, pt, tgE); tgE.Normalize(); tgE2 = tgE1 = tgE;
    let mut tge = curve.derivative_at(par).normalize();
    let mut tge1 = tge;
    let mut tge2 = tge;
    if e.orientation == Orientation::Reversed {
        tge = -tge;
    }

    let e1_fwd = {
        let mut x = e.clone();
        x.orientation = Orientation::Forward;
        x
    };
    let mut e2_fwd = e1_fwd.clone();

    let closed_on_f1 = brep.is_edge_closed_on_face(e, f1);
    if f1.is_same(f2) && closed_on_f1 {
        e2_fwd.orientation = Orientation::Reversed;
        tge2 = -tge2;
    } else {
        // OCCT explores F1/F2 for the edge and reverses the tangent when
        // the in-face orientation is REVERSED.
        let in_face_reversed =
            |f: &Shape| -> Option<bool> {
                let fd = f.as_face()?;
                let ts = brep.tshapes.get(fd.outer_wire.index)?;
                let TShape::Wire(wd) = ts.as_ref() else {
                    return None;
                };
                for we in &wd.edges {
                    if we.is_same(e) {
                        return Some(we.orientation == Orientation::Reversed);
                    }
                }
                None
            };
        match in_face_reversed(f1) {
            Some(true) => tge1 = -tge1,
            Some(false) => {}
            None => return 0,
        }
        match in_face_reversed(f2) {
            Some(true) => tge2 = -tge2,
            Some(false) => {}
            None => return 0,
        }
    }

    let Some(pcurve1) = brep.curve_on_surface(&e1_fwd, f1) else {
        return 0;
    };
    let Some(pcurve2) = brep.curve_on_surface(&e2_fwd, f2) else {
        return 0;
    };
    let mut p2d1 = pcurve1.0.point_at(par);
    let mut p2d2 = pcurve2.0.point_at(par);

    let Some(surf1) = f1.as_face().and_then(|fd| fd.surface.clone()) else {
        return 0;
    };
    let Some(surf2) = f2.as_face().and_then(|fd| fd.surface.clone()) else {
        return 0;
    };
    let (_pt1, du1, dv1) = surf1.derivatives(p2d1.x, p2d1.y);
    let mut ns1 = du1.cross(dv1);
    ns1 = ns1.normalize();
    if f1.orientation == Orientation::Reversed {
        ns1 = -ns1;
    }
    let (_pt2, du2, dv2) = surf2.derivatives(p2d2.x, p2d2.y);
    let mut ns2 = du2.cross(dv2);
    ns2 = ns2.normalize();
    if f2.orientation == Orientation::Reversed {
        ns2 = -ns2;
    }

    let dint1 = ns1.cross(tge1);
    let dint2 = ns2.cross(tge2);
    let ang = ns1.cross(ns2).length();
    if ang > 0.0001 * std::f64::consts::PI {
        let scal = ns2.dot(dint1);
        if scal <= 0.0 {
            ns2 = -ns2;
            *or2 = Orientation::Reversed;
        }
        let scal = ns1.dot(dint2);
        if scal <= 0.0 {
            ns1 = -ns1;
            *or1 = Orientation::Reversed;
        }
    } else {
        // the faces are locally tangent - this is fake!
        if dint1.dot(dint2) < 0.0 {
            // This is a forgotten regularity — OCCT L447-481 re-evaluates
            // the normals with second derivatives (S1.D2/S2.D2) after
            // stepping the UV points along dint; the Surface3 second
            // derivative query is pending, so this sub-branch reports the
            // OCCT "no concave face" code 10.
            return 10;
        }
        // here it turns back, the points are taken in faces
        // neither too close nor too far as much as possible.
        // OCCT ChFi3d_Coefficient(dint, DU, DV, u, v): the (u,v) step of
        // dint expressed in the (DU,DV) frame.
        let coefficient = |dint: DVec3, du: DVec3, dv: DVec3| -> (f64, f64) {
            let a = dint.dot(du);
            let b = dint.dot(dv);
            let duu = du.dot(du);
            let dvv = dv.dot(dv);
            let dudv = du.dot(dv);
            let den = duu * dvv - dudv * dudv;
            if den.abs() <= 0.0 {
                (0.0, 0.0)
            } else {
                ((a * dvv - b * dudv) / den, (b * duu - a * dudv) / den)
            }
        };
        let (u, v) = coefficient(dint1, du1, dv1);
        p2d1.x += u;
        p2d1.y += v;
        let (u, v) = coefficient(dint2, du2, dv2);
        p2d2.x += u;
        p2d2.y += v;
        let (pt1, du1b, dv1b) = surf1.derivatives(p2d1.x, p2d1.y);
        let mut ns1 = du1b.cross(dv1b);
        if f1.orientation == Orientation::Reversed {
            ns1 = -ns1;
        }
        let (pt2, du2b, dv2b) = surf2.derivatives(p2d2.x, p2d2.y);
        let mut ns2 = du2b.cross(dv2b);
        if f2.orientation == Orientation::Reversed {
            ns2 = -ns2;
        }
        let vref = pt2 - pt1;
        if ns1.dot(vref) < 0.0 {
            *or1 = Orientation::Reversed;
        }
        if ns2.dot(vref) > 0.0 {
            *or2 = Orientation::Reversed;
        }
    }

    let mut choix_conge = match (*or1, *or2) {
        (Orientation::Forward, Orientation::Forward) => 1,
        (Orientation::Forward, Orientation::Reversed) => 7,
        (Orientation::Reversed, Orientation::Forward) => 3,
        (Orientation::Reversed, Orientation::Reversed) => 5,
        _ => 1,
    };
    if ns1.cross(ns2).dot(tge) >= 0.0 {
        choix_conge += 1;
    }
    choix_conge
}

// =========================================================================
// OCCT ChFi3d_Builder_1.cxx statics (L60-107, L148-232, L617-683) and
// class methods FaceTangency (L552-615), PerformExtremity (L714-878),
// PerformElement (L887-1176).
// =========================================================================

use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_conexfaces, chfi3d_edge_state, topexp_common_vertex,
    topexp_vertices, vec_angle, vec_is_parallel, shape_key,
};

/// OCCT ChFi3d_Builder_1.cxx L60-107 — static ReorderFaces.
fn reorder_faces(
    brep: &topods::BRep,
    efmap: &ChFiDSMap,
    f1: &mut Shape,
    f2: &mut Shape,
    first_face: &Shape,
    prev_edge_in: &Shape,
    common_vertex: &Shape,
) {
    if f1.is_same(first_face) {
        return;
    } else if f2.is_same(first_face) {
        std::mem::swap(f1, f2);
        return;
    }

    // Loop until find <theF1> or <theF2>
    let mut prev_edge = prev_edge_in.clone();
    let mut cur_edge = Shape::null();
    let mut prev_face = first_face.clone();
    let mut cur_face = Shape::null();
    loop {
        // OCCT: TopExp::MapShapesAndAncestors(PrevFace, VERTEX, EDGE, map)
        // — local vertex -> edges map of the face.
        let mut vertex_edge_map: std::collections::HashMap<u64, Vec<Shape>> =
            std::collections::HashMap::new();
        if let Some(fd) = prev_face.as_face() {
            if let Some(ts) = brep.tshapes.get(fd.outer_wire.index) {
                if let TShape::Wire(wd) = ts.as_ref() {
                    for we in &wd.edges {
                        if let Some(ets) = brep.tshapes.get(we.index) {
                            if let TShape::Edge(edd) = ets.as_ref() {
                                for v in [&edd.first, &edd.last] {
                                    vertex_edge_map
                                        .entry(v.ptr_id())
                                        .or_default()
                                        .push(we.clone());
                                }
                            }
                        }
                    }
                }
            }
        }
        let empty = Vec::new();
        let edges_at = vertex_edge_map.get(&common_vertex.ptr_id()).unwrap_or(&empty);
        for an_edge in edges_at {
            if an_edge.is_same(&prev_edge) {
                continue;
            }

            cur_edge = an_edge.clone();
            let flist = if efmap.contains(&cur_edge) {
                efmap.find(&cur_edge).clone()
            } else {
                Vec::new()
            };
            if flist.is_empty() {
                continue;
            }
            let first_f = &flist[0];
            let last_f = &flist[flist.len() - 1];
            cur_face = if prev_face.is_same(first_f) {
                last_f.clone()
            } else {
                first_f.clone()
            };
            if cur_face.is_same(f1) {
                return;
            } else if cur_face.is_same(f2) {
                std::mem::swap(f1, f2);
                return;
            }
        }

        prev_edge = cur_edge.clone();
        prev_face = cur_face.clone();
    }
}

/// OCCT ChFi3d_Builder_1.cxx L148-232 — static MakeOffsetEdge.
///
/// Pending: depends on Geom_OffsetSurface, GeomInt_IntSS (surface-surface
/// intersection) and Extrema_ExtPC, none translated yet.  Returns a null
/// edge, matching the OCCT failure path.
fn make_offset_edge(
    _brep: &topods::BRep,
    _the_edge: &Shape,
    _distance: f64,
    _s1: &Shape,
    _s2: &Shape,
) -> Shape {
    Shape::null()
}

/// OCCT ChFi3d_Builder_1.cxx L617-662 — static TangentExtremity.
fn tangent_extremity(
    brep: &topods::BRep,
    v: &Shape,
    e: &Shape,
    f1: &Shape,
    f2: &Shape,
    tang: f64,
) -> bool {
    use rcad_kernel::geom::Curve2dEval as _;

    let o1 = f1.orientation;
    let mut f1f = f1.clone();
    f1f.orientation = Orientation::Forward;
    let o2 = f2.orientation;
    let mut f2f = f2.clone();
    f2f.orientation = Orientation::Forward;
    let e1 = {
        let mut x = e.clone();
        x.orientation = Orientation::Forward;
        x
    };
    let mut e2 = e1.clone();
    // OCCT: if (f1.IsSame(f2) && BRep_Tool::IsClosed(e1, f1)) e2 REVERSED.
    if f1f.is_same(&f2f) && brep.is_edge_closed_on_face(&e1, &f1f) {
        e2.orientation = Orientation::Reversed;
    }
    let p1 = brep_tool_parameter(brep, v, &e1);
    let p2 = brep_tool_parameter(brep, v, &e2);
    let eps = 1.0e-9;
    let _ = eps;

    // OCCT: pc1 = BRep_Tool::CurveOnSurface(e1, f1); pc1->Value(p1) -> (u,v)
    let Some(pc1) = brep.curve_on_surface(&e1, &f1f) else {
        return false;
    };
    let uv1 = pc1.0.point_at(p1);
    let Some(n1_raw) = topopebreptool_nt(brep, uv1, &f1f) else {
        return false; // It is not known...
    };
    let mut n1 = n1_raw;
    if o1 == Orientation::Reversed {
        n1 = -n1;
    }

    let Some(pc2) = brep.curve_on_surface(&e2, &f2f) else {
        return false;
    };
    let uv2 = pc2.0.point_at(p2);
    let Some(n2_raw) = topopebreptool_nt(brep, uv2, &f2f) else {
        return false; // It is not known...
    };
    let mut n2 = n2_raw;
    if o2 == Orientation::Reversed {
        n2 = -n2;
    }

    vec_angle(n1, n2) < tang
}

/// OCCT ChFi3d_Builder_1.cxx L665-683 — static TangentOnVertex.
fn tangent_on_vertex(
    brep: &topods::BRep,
    v: &Shape,
    e: &Shape,
    efmap: &ChFiDSMap,
    tang: f64,
) -> bool {
    let (ff1, ff2) = chfi3d_conexfaces(e, efmap);
    if ff1.is_null() || ff2.is_null() {
        return false;
    }
    tangent_extremity(brep, v, e, &ff1, &ff2, tang)
}

impl ChFi3dBuilder {
    /// OCCT ChFi3d_Builder_1.cxx L552-615 — FaceTangency.
    fn face_tangency(&self, e0: &Shape, e1: &Shape, v: &Shape) -> bool {
        let mut f = [Shape::null(), Shape::null()];
        let mut nbf = 0usize;

        // It is checked if the connection is not on a regular edge.
        let e1_faces = if self.my_ef_map.contains(e1) {
            self.my_ef_map.find(e1).clone()
        } else {
            Vec::new()
        };
        for it in &e1_faces {
            if nbf > 1 {
                panic!("Standard_ConstructionError: ChFi3d_Builder:only 2 faces");
            }
            if nbf < 2 {
                f[nbf] = it.clone();
            }
            nbf += 1;
        }
        if nbf < 2 {
            return false;
        }
        if is_tangent_faces(&self.my_brep, e1, &f[0], &f[1], GeomAbsShape::G1) {
            return false;
        }

        let v_edges = if self.my_ve_map.contains(v) {
            self.my_ve_map.find(v).clone()
        } else {
            Vec::new()
        };
        for ec in &v_edges {
            if !ec.is_same(e0)
                && !ec.is_same(e1)
                && ec.orientation != Orientation::Internal
                && ec.orientation != Orientation::External
            {
                let ed = ec.as_edge().expect("not an edge");
                if ed.degenerated {
                    continue;
                }
                let ec_faces = if self.my_ef_map.contains(ec) {
                    self.my_ef_map.find(ec).clone()
                } else {
                    Vec::new()
                };
                let mut nbf2 = 0usize;
                for it in &ec_faces {
                    if nbf2 > 1 {
                        panic!(
                            "Standard_ConstructionError: ChFi3d_Builder:only 2 faces"
                        );
                    }
                    if nbf2 < 2 {
                        f[nbf2] = it.clone();
                    }
                    nbf2 += 1;
                }
                if nbf2 < 2 {
                    return false;
                }
                if !is_tangent_faces(&self.my_brep, ec, &f[0], &f[1], GeomAbsShape::G1) {
                    return false;
                }
            }
        }
        true
    }

    /// OCCT ChFi3d_Builder_1.cxx L714-878 — PerformExtremity.
    pub fn perform_extremity(&mut self, spine: &mut ChFiDSSpineHandle) {
        let mut nb_g1_connections = 0i32;

        for ii in 1..=2 {
            let mut e = [Shape::null(), Shape::null(), Shape::null()];
            let mut sst;
            let iedge;
            let v;
            if ii == 1 {
                sst = spine.base().first_status();
                iedge = 1;
                v = spine.base().first_vertex();
            } else {
                sst = spine.base().last_status();
                iedge = spine.base().nb_edges();
                v = spine.base().last_vertex();
            }
            // Before all it is checked if the tangency is not dead.
            e[0] = spine.base().edges(iedge).clone();
            // OCCT: ConexFaces(Spine, iedge, hs1, hs2) — the two support
            // faces of the end edge.
            let (hs1, hs2) = chfi3d_conexfaces(&e[0], &self.my_ef_map);
            if tangent_extremity(&self.my_brep, &v, &e[0], &hs1, &hs2, self.angular) {
                spine.base_mut().set_tangency_extremity(true, ii == 1);
            }

            if sst == ChFiDS_State::BreakPoint {
                let mut a_loc_nb_g1_connections = 0i32;
                let mut sommetpourri = false;
                // OCCT: Edges (map) + EdgesOfV (indexed map), keyed by
                // TShape+Location+Orientation.
                let mut edges: std::collections::HashSet<(u64, u32, u8)> =
                    std::collections::HashSet::new();
                let mut edges_of_v: Vec<Shape> = Vec::new();
                edges.insert(shape_key(&e[0]));
                edges_of_v.push(e[0].clone());
                let mut ind_of_e = 0usize;

                let v_edges = if self.my_ve_map.contains(&v) {
                    self.my_ve_map.find(&v).clone()
                } else {
                    Vec::new()
                };
                for an_edge in &v_edges {
                    let ed = an_edge.as_edge().expect("not an edge");
                    if ed.degenerated {
                        continue;
                    }
                    let (f1, f2) = chfi3d_conexfaces(an_edge, &self.my_ef_map);
                    if !f2.is_null()
                        && is_tangent_faces(&self.my_brep, an_edge, &f1, &f2, GeomAbsShape::G2)
                    {
                        // smooth edge
                        if !f1.is_same(&f2) {
                            nb_g1_connections += 1;
                            a_loc_nb_g1_connections += 1;
                        }
                        continue;
                    }

                    if edges.insert(shape_key(an_edge)) {
                        edges_of_v.push(an_edge.clone());
                        if ind_of_e < 2 {
                            ind_of_e += 1;
                            e[ind_of_e] = an_edge.clone();
                        }
                    } else {
                        // OCCT L785-805: an edge already seen — closed edge
                        // case, two ends of the edge in the vertex.
                        let (v1, v2) = topexp_vertices(an_edge);
                        if v1.is_same(&v2) {
                            let mut an_ind = edges_of_v
                                .iter()
                                .position(|s| shape_key(s) == shape_key(an_edge));
                            if an_ind.is_none() {
                                let rev = {
                                    let mut x = an_edge.clone();
                                    x.orientation =
                                        if x.orientation == Orientation::Forward {
                                            Orientation::Reversed
                                        } else {
                                            Orientation::Forward
                                        };
                                    x
                                };
                                an_ind = edges_of_v
                                    .iter()
                                    .position(|s| shape_key(s) == shape_key(&rev));
                                if let Some(idx) = an_ind {
                                    let mut kept = edges_of_v[idx].clone();
                                    kept.orientation =
                                        if kept.orientation == Orientation::Forward {
                                            Orientation::Reversed
                                        } else {
                                            Orientation::Forward
                                        };
                                    if edges.insert(shape_key(&kept)) {
                                        if ind_of_e < 2 {
                                            ind_of_e += 1;
                                            e[ind_of_e] = kept.clone();
                                        }
                                        edges_of_v.push(kept);
                                    }
                                }
                            } else {
                                // Same oriented edge present twice: the
                                // OCCT code finds it in EdgesOfV, reverses
                                // and re-adds.
                                if let Some(idx) = an_ind {
                                    let mut kept = edges_of_v[idx].clone();
                                    kept.orientation =
                                        if kept.orientation == Orientation::Forward {
                                            Orientation::Reversed
                                        } else {
                                            Orientation::Forward
                                        };
                                    if edges.insert(shape_key(&kept)) {
                                        if ind_of_e < 2 {
                                            ind_of_e += 1;
                                            e[ind_of_e] = kept.clone();
                                        }
                                        edges_of_v.push(kept);
                                    }
                                }
                            }
                        }
                    }
                }

                if edges_of_v.len() != 3 {
                    sommetpourri = true;
                }

                if !sommetpourri && a_loc_nb_g1_connections < 4 {
                    sst = chfi3d_edge_state(&e, &self.my_ef_map, &self.my_brep);
                }
                if ii == 1 {
                    spine.base_mut().set_first_status(sst);
                } else {
                    spine.base_mut().set_last_status(sst);
                }
            }
        }

        if !spine.base().is_periodic() {
            // OCCT L830-877: count the distinct faces at each end vertex
            // (IsSame-deduplicated over myVFMap) and mark BreakPoint when
            // more than 3.
            let count_distinct_faces = |vertex: &Shape| -> i32 {
                let faces = if self.my_vf_map.contains(vertex) {
                    self.my_vf_map.find(vertex).clone()
                } else {
                    Vec::new()
                };
                let mut nbf = 0i32;
                for (jf, cur) in faces.iter().enumerate() {
                    let mut seen = false;
                    for prev in faces.iter().take(jf) {
                        if cur.is_same(prev) {
                            seen = true;
                            break;
                        }
                    }
                    if !seen {
                        nbf += 1;
                    }
                }
                nbf
            };
            let mut nbf = count_distinct_faces(&spine.base().first_vertex());
            nbf -= nb_g1_connections;
            if nbf > 3 {
                spine.base_mut().set_first_status(ChFiDS_State::BreakPoint);
            }
            let mut nbf = count_distinct_faces(&spine.base().last_vertex());
            nbf -= nb_g1_connections;
            if nbf > 3 {
                spine.base_mut().set_last_status(ChFiDS_State::BreakPoint);
            }
        }
    }

    /// OCCT ChFi3d_Builder_1.cxx L887-1176 — PerformElement: find all
    /// mutually tangent edges.  Each edge has 2 opposing faces; for 2
    /// adjacent tangent edges it is required that the opposing faces were
    /// tangent.
    pub fn perform_element(
        &mut self,
        spine: &mut ChFiDSSpineHandle,
        offset: f64,
        the_first_face: &Shape,
    ) -> bool {
        use rcad_kernel::geom::CurveEval as _;
        let ta = self.angular;

        let mut ec = spine.base().edges(1).clone();
        {
            let ed = ec.as_edge().expect("not an edge");
            if ed.degenerated {
                return false;
            }
        }
        // it is checked if the edge is a cut edge
        let (mut ff1, mut ff2) = chfi3d_conexfaces(&ec, &self.my_ef_map);
        if ff1.is_null() || ff2.is_null() {
            return false;
        }
        if is_tangent_faces(&self.my_brep, &ec, &ff1, &ff2, GeomAbsShape::G1) {
            return false;
        }

        let mut first_face = ff1.clone();
        if !the_first_face.is_null() && ff2.is_same(the_first_face) {
            first_face = ff2.clone();
            std::mem::swap(&mut ff1, &mut ff2);
            ff1 = first_face.clone();
        }
        self.my_edge_first_face.insert(ec.ptr_id(), first_face.clone());

        // Define concavity
        let type_of_concavity =
            define_connect_type(&self.my_brep, &ec, &ff1, &ff2, 1.0e-5, true);
        spine.base_mut().set_type_of_concavity(type_of_concavity);

        let to_restrict = offset > 0.0;
        let _ = to_restrict; // BRepAdaptor_Surface restriction — pending
        if offset > 0.0 {
            let offset_edge = make_offset_edge(&self.my_brep, &ec, offset, &ff1, &ff2);
            // OCCT L939: OffsetEdge.Orientation(Ec.Orientation());
            let mut oe = offset_edge;
            oe.orientation = ec.orientation;
            spine.base_mut().set_offset_edges(oe);
        }

        let curor = ec.orientation;
        let (v_start, mut lvec) = topexp_vertices(&ec);

        let mut fini = false;
        let mut cur_st = ChFiDS_State::Closed;
        let edge_curve =
            |s: &Shape| -> Option<rcad_kernel::geom::Curve3> {
                s.as_edge().and_then(|ed| ed.curve.clone())
            };

        if v_start.is_same(&lvec) {
            // case if only one edge is closed
            let Some(c_ec) = edge_curve(&ec) else {
                return false;
            };
            let mut wl = brep_tool_parameter(&self.my_brep, &v_start, &ec);
            let mut v1 = c_ec.derivative_at(wl);
            wl = brep_tool_parameter(&self.my_brep, &lvec, &ec);
            let v2 = c_ec.derivative_at(wl);
            let is_face_tangency = self.face_tangency(&ec, &ec, &v_start);
            if vec_is_parallel(v1, v2, ta) || is_face_tangency {
                if is_face_tangency {
                    cur_st = ChFiDS_State::Closed;
                } else {
                    cur_st = ChFiDS_State::BreakPoint;
                }
            } else {
                cur_st = ChFiDS_State::BreakPoint;
            }
            spine.base_mut().set_last_status(cur_st);
            spine.base_mut().set_first_status(cur_st);
        } else {
            // Downstream progression
            let mut v1;
            let mut cur_or = curor;
            while !fini {
                cur_st = ChFiDS_State::FreeBoundary;
                let wl = brep_tool_parameter(&self.my_brep, &lvec, &ec);
                let degene_on_ec = tangent_on_vertex(&self.my_brep, &lvec, &ec, &self.my_ef_map, ta);
                let Some(c_ec) = edge_curve(&ec) else {
                    return false;
                };
                v1 = c_ec.derivative_at(wl);
                let nb = spine.base().nb_edges();

                let lv_edges = if self.my_ve_map.contains(&lvec) {
                    self.my_ve_map.find(&lvec).clone()
                } else {
                    Vec::new()
                };
                for ev in &lv_edges {
                    if ev.is_same(&ec) {
                        continue;
                    }
                    let edv = ev.as_edge().expect("not an edge");
                    if edv.degenerated {
                        continue;
                    }
                    let (mut fvev, mut lvev) = topexp_vertices(ev);
                    if lvec.is_same(&lvev) {
                        let ve1 = fvev;
                        fvev = lvev;
                        lvev = ve1;
                        let or1 = Orientation::Reversed;

                        let wf = brep_tool_parameter(&self.my_brep, &fvev, ev);
                        let Some(c_ev) = edge_curve(ev) else {
                            continue;
                        };
                        let v2 = c_ev.derivative_at(wf);
                        let av1v2 = vec_angle(v1, v2);
                        let rev = or1 != cur_or;
                        let mut on_ajoute = false;
                        if self.face_tangency(&ec, ev, &fvev) {
                            on_ajoute =
                                (!rev && av1v2 < std::f64::consts::PI / 2.0)
                                    || (rev && av1v2 > std::f64::consts::PI / 2.0);
                            if on_ajoute
                                && (degene_on_ec
                                    || tangent_on_vertex(
                                        &self.my_brep, &lvec, ev, &self.my_ef_map, ta,
                                    ))
                            {
                                on_ajoute = (!rev && av1v2 < ta)
                                    || (rev && (std::f64::consts::PI - av1v2) < ta);
                            }
                        }
                        if on_ajoute {
                            fini = false; // If this can be useful (Cf PRO14713)
                            let common_vertex = topexp_common_vertex(&ec, ev);
                            let prev_edge = ec.clone();
                            ec = ev.clone();
                            ec.orientation = or1;
                            lvec = lvev.clone();
                            spine.base_mut().set_edges(ec.clone());
                            let (mut cur_f1, mut cur_f2) =
                                chfi3d_conexfaces(&ec, &self.my_ef_map);
                            if let Some(cv) = common_vertex {
                                reorder_faces(
                                    &self.my_brep,
                                    &self.my_ef_map,
                                    &mut cur_f1,
                                    &mut cur_f2,
                                    &first_face,
                                    &prev_edge,
                                    &cv,
                                );
                            }
                            self.my_edge_first_face.insert(ec.ptr_id(), cur_f1.clone());
                            if offset > 0.0 {
                                let an_offset_edge =
                                    make_offset_edge(&self.my_brep, &ec, offset, &cur_f1, &cur_f2);
                                let mut oe = an_offset_edge;
                                oe.orientation = or1;
                                spine.base_mut().set_offset_edges(oe);
                            }
                            first_face = cur_f1;
                            cur_or = or1;
                            if v_start.is_same(&lvev) {
                                if self.face_tangency(ev, &spine.base().edges(1).clone(), &lvev) {
                                    cur_st = ChFiDS_State::Closed;
                                    fini = true;
                                } else {
                                    cur_st = ChFiDS_State::BreakPoint;
                                    fini = true;
                                }
                            }
                            break;
                        } else {
                            let nbface = if self.my_ef_map.contains(ev) {
                                self.my_ef_map.find(ev).len()
                            } else {
                                0
                            };
                            if nbface > 1 {
                                cur_st = ChFiDS_State::BreakPoint;
                            }
                            fini = (!rev && av1v2 < ta)
                                || (rev && (std::f64::consts::PI - av1v2) < ta);
                        }
                    } else {
                        let or1 = Orientation::Forward;

                        let wf = brep_tool_parameter(&self.my_brep, &fvev, ev);
                        let Some(c_ev) = edge_curve(ev) else {
                            continue;
                        };
                        let v2 = c_ev.derivative_at(wf);
                        let av1v2 = vec_angle(v1, v2);
                        let rev = or1 != cur_or;
                        let mut on_ajoute = false;
                        if self.face_tangency(&ec, ev, &fvev) {
                            on_ajoute =
                                (!rev && av1v2 < std::f64::consts::PI / 2.0)
                                    || (rev && av1v2 > std::f64::consts::PI / 2.0);
                            if on_ajoute
                                && (degene_on_ec
                                    || tangent_on_vertex(
                                        &self.my_brep, &lvec, ev, &self.my_ef_map, ta,
                                    ))
                            {
                                on_ajoute = (!rev && av1v2 < ta)
                                    || (rev && (std::f64::consts::PI - av1v2) < ta);
                            }
                        }
                        if on_ajoute {
                            fini = false;
                            let common_vertex = topexp_common_vertex(&ec, ev);
                            let prev_edge = ec.clone();
                            ec = ev.clone();
                            ec.orientation = or1;
                            lvec = lvev.clone();
                            spine.base_mut().set_edges(ec.clone());
                            let (mut cur_f1, mut cur_f2) =
                                chfi3d_conexfaces(&ec, &self.my_ef_map);
                            if let Some(cv) = common_vertex {
                                reorder_faces(
                                    &self.my_brep,
                                    &self.my_ef_map,
                                    &mut cur_f1,
                                    &mut cur_f2,
                                    &first_face,
                                    &prev_edge,
                                    &cv,
                                );
                            }
                            self.my_edge_first_face.insert(ec.ptr_id(), cur_f1.clone());
                            if offset > 0.0 {
                                let an_offset_edge =
                                    make_offset_edge(&self.my_brep, &ec, offset, &cur_f1, &cur_f2);
                                let mut oe = an_offset_edge;
                                oe.orientation = or1;
                                spine.base_mut().set_offset_edges(oe);
                            }
                            first_face = cur_f1;
                            cur_or = or1;
                            if v_start.is_same(&lvev) {
                                if self.face_tangency(ev, &spine.base().edges(1).clone(), &lvev) {
                                    cur_st = ChFiDS_State::Closed;
                                    fini = true;
                                } else {
                                    cur_st = ChFiDS_State::BreakPoint;
                                    fini = true;
                                }
                            }
                            break;
                        } else {
                            let nbface = if self.my_ef_map.contains(ev) {
                                self.my_ef_map.find(ev).len()
                            } else {
                                0
                            };
                            if nbface > 1 {
                                cur_st = ChFiDS_State::BreakPoint;
                            }
                            fini = (!rev && av1v2 < ta)
                                || (rev && (std::f64::consts::PI - av1v2) < ta);
                        }
                    }
                }
                fini = fini || (nb == spine.base().nb_edges());
            }
            spine.base_mut().set_last_status(cur_st);
            if cur_st == ChFiDS_State::Closed {
                spine.base_mut().set_first_status(cur_st);
            } else {
                // Upstream progression
                fini = false;
                ec = spine.base().edges(1).clone();
                first_face = self
                    .my_edge_first_face
                    .get(&ec.ptr_id())
                    .cloned()
                    .unwrap_or_else(Shape::null);
                let mut cur_or = ec.orientation;
                let mut fvec = v_start;
                while !fini {
                    cur_st = ChFiDS_State::FreeBoundary;
                    let wl = brep_tool_parameter(&self.my_brep, &fvec, &ec);
                    let degene_on_ec =
                        tangent_on_vertex(&self.my_brep, &fvec, &ec, &self.my_ef_map, ta);
                    let Some(c_ec) = edge_curve(&ec) else {
                        return false;
                    };
                    let v1 = c_ec.derivative_at(wl);
                    let nb = spine.base().nb_edges();

                    let fv_edges = if self.my_ve_map.contains(&fvec) {
                        self.my_ve_map.find(&fvec).clone()
                    } else {
                        Vec::new()
                    };
                    for ev in &fv_edges {
                        if ev.is_same(&ec) {
                            continue;
                        }
                        let edv = ev.as_edge().expect("not an edge");
                        if edv.degenerated {
                            continue;
                        }
                        let (mut fvev, mut lvev) = topexp_vertices(ev);
                        let or1;
                        if fvec.is_same(&fvev) {
                            let ve1 = fvev;
                            fvev = lvev;
                            lvev = ve1;
                            or1 = Orientation::Reversed;
                        } else {
                            or1 = Orientation::Forward;
                        }
                        let wf = brep_tool_parameter(&self.my_brep, &lvev, ev);
                        let Some(c_ev) = edge_curve(ev) else {
                            continue;
                        };
                        let v2 = c_ev.derivative_at(wf);
                        let av1v2 = vec_angle(v1, v2);
                        let rev = or1 != cur_or;
                        let mut on_ajoute = false;
                        if self.face_tangency(&ec, ev, &lvev) {
                            on_ajoute =
                                (!rev && av1v2 < std::f64::consts::PI / 2.0)
                                    || (rev && av1v2 > std::f64::consts::PI / 2.0);
                            if on_ajoute
                                && (degene_on_ec
                                    || tangent_on_vertex(
                                        &self.my_brep, &fvec, ev, &self.my_ef_map, ta,
                                    ))
                            {
                                on_ajoute = (!rev && av1v2 < ta)
                                    || (rev && (std::f64::consts::PI - av1v2) < ta);
                            }
                        }
                        if on_ajoute {
                            let common_vertex = topexp_common_vertex(&ec, ev);
                            let prev_edge = ec.clone();
                            ec = ev.clone();
                            ec.orientation = or1;
                            fvec = fvev.clone();
                            spine.base_mut().put_in_first(ec.clone());
                            let (mut cur_f1, mut cur_f2) =
                                chfi3d_conexfaces(&ec, &self.my_ef_map);
                            if let Some(cv) = common_vertex {
                                reorder_faces(
                                    &self.my_brep,
                                    &self.my_ef_map,
                                    &mut cur_f1,
                                    &mut cur_f2,
                                    &first_face,
                                    &prev_edge,
                                    &cv,
                                );
                            }
                            self.my_edge_first_face.insert(ec.ptr_id(), cur_f1.clone());
                            if offset > 0.0 {
                                let an_offset_edge =
                                    make_offset_edge(&self.my_brep, &ec, offset, &cur_f1, &cur_f2);
                                let mut oe = an_offset_edge;
                                oe.orientation = or1;
                                spine.base_mut().put_in_first_offset(oe);
                            }
                            first_face = cur_f1;
                            cur_or = or1;
                            break;
                        } else {
                            let nbface = if self.my_ef_map.contains(ev) {
                                self.my_ef_map.find(ev).len()
                            } else {
                                0
                            };
                            if nbface > 1 {
                                cur_st = ChFiDS_State::BreakPoint;
                            }
                            fini = (!rev && av1v2 < ta)
                                || (rev && (std::f64::consts::PI - av1v2) < ta);
                        }
                    }
                    fini = fini || (nb == spine.base().nb_edges());
                }
                spine.base_mut().set_first_status(cur_st);
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::brep_fillet_api::{explore_edges, explore_solids};
    use super::super::chfi_ds::ChFi3dFilletShape;

    /// OCCT BRepFilletAPI_MakeFillet::Add on a box: the seed edge builds one
    /// contour; PerformElement finds no tangential continuation on the
    /// sharp box edges and PerformExtremity refines the end statuses.
    #[test]
    fn add_box_edge_builds_contour() {
        let brep = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            20.0,
            20.0,
            20.0,
        )
        .unwrap();
        let solid = explore_solids(&brep).into_iter().next().unwrap();
        let mut fillet =
            ChFi3dFilBuilder::new(&brep, solid, ChFi3dFilletShape::Rational, 1.0e-2);

        let edges = explore_edges(&brep);
        assert!(!edges.is_empty());
        let e0 = edges[0].clone();

        fillet.add(&e0);
        assert_eq!(fillet.base.nb_elements(), 1, "one contour expected");
        assert_eq!(fillet.base.nb_faulty_contours(), 0, "no faulty contours");

        let spine = fillet.base.value(1);
        assert_eq!(
            spine.base().nb_edges(),
            1,
            "sharp box edge has no tangent neighbours to chain"
        );
        // TangentExtremity on perpendicular box faces is not tangent.
        assert!(
            !spine.base().is_tangency_extremity(true)
                && !spine.base().is_tangency_extremity(false),
            "box edge ends are not tangent-dead"
        );
        // Add twice on the same edge must not duplicate the contour.
        fillet.add(&e0);
        assert_eq!(fillet.base.nb_elements(), 1, "Contains dedup");
    }

    /// OCCT ChFi3d_EdgeState over the three edges of a box vertex.
    #[test]
    fn edge_state_on_box_vertex() {
        let brep = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            20.0,
            20.0,
            20.0,
        )
        .unwrap();
        let efmap = ChFiDSMap::new();
        let mut efmap = efmap;
        efmap.fill(&brep, topods::ShapeType::Edge, topods::ShapeType::Face);
        // Two sharp box edges sharing a vertex: FreeBoundary only when a
        // face is missing; here the faces exist so the state is one of the
        // concavity outcomes.
        let edges = explore_edges(&brep);
        let st = chfi3d_edge_state(&[edges[0].clone(), edges[1].clone(), edges[2].clone()], &efmap, &brep);
        assert!(
            matches!(
                st,
                ChFiDS_State::AllSame | ChFiDS_State::OnSame | ChFiDS_State::OnDiff | ChFiDS_State::BreakPoint
            ),
            "box vertex edge state = {st:?}"
        );
    }
}

// =========================================================================
// OCCT ChFi3d.cxx L525-612 — ChFi3d::NextSide (both overloads) and
// TopAbs::Compose helper.
// =========================================================================

/// OCCT TopAbs::Compose(a, b): FORWARD keeps b, REVERSED reverses b.
pub fn topabs_compose(a: Orientation, b: Orientation) -> Orientation {
    if a == Orientation::Reversed {
        if b == Orientation::Reversed {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    } else {
        b
    }
}

/// OCCT TopAbs::Reverse(o).
pub fn topabs_reverse(o: Orientation) -> Orientation {
    if o == Orientation::Reversed {
        Orientation::Forward
    } else {
        Orientation::Reversed
    }
}

/// OCCT ChFi3d.cxx L525-612 — ChFi3d::NextSide(Or1, Or2, OrSave1, OrSave2,
/// ChoixSave).
pub fn next_side(
    or1: &mut Orientation,
    or2: &mut Orientation,
    or_save1: Orientation,
    or_save2: Orientation,
    choix_save: i32,
) -> i32 {
    *or1 = if *or1 == Orientation::Forward {
        or_save1
    } else {
        topabs_reverse(or_save1)
    };
    *or2 = if *or2 == Orientation::Forward {
        or_save2
    } else {
        topabs_reverse(or_save2)
    };

    let mut choix_conge;
    if *or1 == Orientation::Forward {
        if *or2 == Orientation::Forward {
            choix_conge = 1;
        } else if choix_save < 0 {
            choix_conge = 3;
        } else {
            choix_conge = 7;
        }
    } else if *or2 == Orientation::Forward {
        if choix_save < 0 {
            choix_conge = 7;
        } else {
            choix_conge = 3;
        }
    } else {
        choix_conge = 5;
    }
    if choix_save.abs() % 2 == 0 {
        choix_conge += 1;
    }
    choix_conge
}

impl super::chfi_ds::ChFiDSElSpine {
    /// OCCT ChFiDS_ElSpine default construction.
    pub fn new() -> Self {
        super::chfi_ds::ChFiDSElSpine {
            firstparam: 0.0,
            lastparam: 0.0,
            firstpnt: DVec3::ZERO,
            firsttgt: DVec3::ZERO,
            lastpnt: DVec3::ZERO,
            lasttgt: DVec3::ZERO,
            periodic: false,
            next: None,
            previous: None,
        }
    }

    /// OCCT ChFiDS_ElSpine.hxx — SetFirstPointAndTgt(P, T).
    pub fn set_first_point_and_tgt(&mut self, p: DVec3, t: DVec3) {
        self.firstpnt = p;
        self.firsttgt = t;
    }

    /// OCCT ChFiDS_ElSpine.hxx — SetLastPointAndTgt(P, T).
    pub fn set_last_point_and_tgt(&mut self, p: DVec3, t: DVec3) {
        self.lastpnt = p;
        self.lasttgt = t;
    }
}

impl ChFi3dBuilder {
    /// OCCT ChFi3d_Builder_2.cxx L859-889 — ConexFaces: the two support
    /// faces of the iedge-th elementary spine, reordered against the first
    /// face.  rcad returns the face Shapes (OCCT wraps them into
    /// BRepAdaptor_Surface handles).
    pub fn conex_faces(&self, spine: &ChFiDSSpineHandle, iedge: usize) -> (Shape, Shape) {
        let mut ff1;
        let mut ff2;
        let an_edge = spine.base().edges(iedge).clone();
        let (c1, c2) = chfi3d_conexfaces(&an_edge, &self.my_ef_map);
        ff1 = c1;
        ff2 = c2;

        let first_face = self
            .my_edge_first_face
            .get(&an_edge.ptr_id())
            .cloned()
            .unwrap_or_else(Shape::null);
        if ff2.is_same(&first_face) {
            let tmp = ff1;
            ff1 = ff2;
            ff2 = tmp;
        }
        (ff1, ff2)
    }

    /// OCCT ChFi3d_Builder_2.cxx L826-855 — StripeOrientations.
    pub fn stripe_orientations(
        &self,
        spine: &ChFiDSSpineHandle,
    ) -> (Orientation, Orientation, i32) {
        let an_edge = spine.base().edges(1).clone();
        let first_face = self
            .my_edge_first_face
            .get(&an_edge.ptr_id())
            .cloned()
            .unwrap_or_else(Shape::null);
        let (mut ff1, mut ff2) = chfi3d_conexfaces(&an_edge, &self.my_ef_map);
        if ff2.is_same(&first_face) {
            let tmp = ff1;
            ff1 = ff2;
            ff2 = tmp;
        }
        let of1 = ff1.orientation;
        let mut f1f = ff1.clone();
        f1f.orientation = Orientation::Forward;
        let of2 = ff2.orientation;
        let mut f2f = ff2.clone();
        f2f.orientation = Orientation::Forward;

        let (mut or1, mut or2) = (Orientation::Forward, Orientation::Forward);
        let choix_conge =
            concave_side(&self.my_brep, &f1f, &f2f, spine.base().edges(1), &mut or1, &mut or2);
        or1 = topabs_compose(or1, of1);
        or2 = topabs_compose(or2, of2);
        (or1, or2, choix_conge)
    }

    /// OCCT ChFi3d_Builder_0.cxx L604-678 — ChFi3d_KParticular.
    pub fn chfi3d_k_particular(
        &self,
        spine: &ChFiDSSpineHandle,
        ie: usize,
        s1: &Shape,
        s2: &Shape,
    ) -> bool {
        use rcad_kernel::geom::Surface3;

        // OCCT: down_cast<FilSpine> — a variable radius disqualifies.
        if let Some(fs) = spine.down_cast_fil() {
            if !fs.is_constant_on(ie) {
                return false;
            }
        }

        let surf_kind = |s: &Shape| -> Option<&'static str> {
            let fd = s.as_face()?;
            let surf = fd.surface.as_ref()?;
            Some(match surf {
                Surface3::Plane(_) => "Plane",
                Surface3::Cylinder(_) => "Cylinder",
                Surface3::Cone(_) => "Cone",
                _ => "Other",
            })
        };
        let st1 = surf_kind(s1).unwrap_or("Other");
        let st2 = surf_kind(s2).unwrap_or("Other");
        let is_plane1 = st1 == "Plane";
        let is_plane2 = st2 == "Plane";
        if !(is_plane1 || is_plane2) {
            return false;
        }
        let (a_s1, a_s2, a_st1, a_st2) = if is_plane1 {
            (s1, s2, st1, st2)
        } else {
            (s2, s1, st2, st1)
        };

        if a_st2 != "Plane" && a_st2 != "Cylinder" && a_st2 != "Cone" {
            return false;
        }

        // OCCT: Spine->CurrentElementarySpine(IE).GetType() ∈ {Line, Circle}.
        let e = spine.base().edges(ie).clone();
        let ctyp = e
            .as_edge()
            .and_then(|ed| ed.curve.as_ref())
            .map(|c| match c {
                rcad_kernel::geom::Curve3::Line(_) => "Line",
                rcad_kernel::geom::Curve3::Circle(_) => "Circle",
                _ => "Other",
            })
            .unwrap_or("Other");
        if ctyp != "Line" && ctyp != "Circle" {
            return false;
        }

        let pa = rcad_kernel::core::precision::ANGULAR;

        let plane_axis = |s: &Shape| -> Option<DVec3> {
            match s.as_face()?.surface.as_ref()? {
                Surface3::Plane(p) => Some(p.normal.normalize()),
                _ => None,
            }
        };
        let cyl_axis = |s: &Shape| -> Option<DVec3> {
            match s.as_face()?.surface.as_ref()? {
                Surface3::Cylinder(c) => Some(c.axis.normalize()),
                _ => None,
            }
        };
        let cone_axis = |s: &Shape| -> Option<DVec3> {
            match s.as_face()?.surface.as_ref()? {
                Surface3::Cone(c) => Some(c.axis.normalize()),
                _ => None,
            }
        };

        if a_st2 == "Plane" {
            if ctyp == "Line" {
                return true;
            }
        } else if a_st2 == "Cylinder" {
            let d1 = plane_axis(a_s1).unwrap_or(DVec3::ZERO);
            let d2 = cyl_axis(a_s2).unwrap_or(DVec3::ZERO);
            let dot = d1.dot(d2).abs();
            if ctyp == "Line" && (1.0 - dot) <= pa {
                // IsNormal
                return true;
            } else if ctyp == "Circle" && dot >= 1.0 - pa {
                // IsParallel
                return true;
            }
        } else if a_st2 == "Cone" {
            let d1 = plane_axis(a_s1).unwrap_or(DVec3::ZERO);
            let d2 = cone_axis(a_s2).unwrap_or(DVec3::ZERO);
            let dot = d1.dot(d2).abs();
            if ctyp == "Circle" && dot >= 1.0 - pa {
                return true;
            }
        }
        false
    }

    /// OCCT ChFi3d_Builder_2.cxx L373-389 — static TgtKP.
    fn tgt_kp(
        cd: &ChFiDSSurfData,
        spine: &mut ChFiDSSpineHandle,
        iedge: usize,
        isfirst: bool,
    ) -> (DVec3, DVec3) {
        use rcad_kernel::geom::CurveEval as _;
        let wtg = cd.interference_on_s1().parameter(isfirst);
        let e = spine.base().edges(iedge).clone();
        let ed = e.as_edge().expect("not an edge");
        let curve = ed.curve.as_ref().expect("edge curve").clone();
        let (cf, cl) = (ed.range[0], ed.range[1]);
        let (ped, mut ded) = if e.orientation == Orientation::Forward {
            let u = wtg + cf;
            (curve.point_at(u), curve.derivative_at(u))
        } else {
            let u = -wtg + cl;
            (curve.point_at(u), -curve.derivative_at(u))
        };
        (ped, ded.normalize())
    }

    /// OCCT ChFi3d_Builder_SpKP.cxx L749-... - SplitKPart.
    ///
    /// The Geom2dHatch hatcher (2D trimming of the tangency lines against
    /// the face restrictions) is a pending TKGeomAlgo translation; the
    /// pending state follows the OCCT single-domain outcome (the SurfData
    /// is appended unsplit), exact when the tangency line crosses no face
    /// boundary.
    #[allow(clippy::too_many_arguments)]
    fn split_k_part(
        &mut self,
        data: &mut ChFiDSSurfData,
        set_data: &mut Vec<ChFiDSSurfData>,
        s1: &Shape,
        s2: &Shape,
        _spine: &mut ChFiDSSpineHandle,
        _iedge: usize,
        _intf: bool,
        _intl: bool,
    ) -> bool {
        // OCCT SpKP.cxx L856-871: F1/F2 from the adaptor handles are
        // registered in the DS and stored on the SurfData.
        let f1 = s1.clone();
        let f2 = s2.clone();
        let dstr = self.my_ds.as_mut().expect("DS");
        data.change_index_of_s1(dstr.add_shape(&f1));
        data.change_index_of_s2(dstr.add_shape(&f2));
        // OCCT L762-855: Geom2dHatch_Hatcher trim of InterferenceOnS1/S2
        // pcurves — pending; the multi-domain redistribution below it is
        // unreachable in the pending state.  The pending state follows the
        // OCCT single-domain outcome (NbDomains == 1 on both faces): the
        // SurfData is appended unsplit.
        set_data.push(data.clone());
        true
    }

    /// OCCT ChFi3d_Builder_2.cxx L3004-3280 — PerformSetOfKPart.
    pub fn perform_set_of_k_part(&mut self, stripe: &SharedStripe, simul: bool) {
        // L3013-3017: initialization of the stripe.
        let mut st = stripe.write().expect("stripe lock");
        st.reset();
        st.my_hdata = Vec::new();
        let mut spine_opt = st.my_spine.take();
        let Some(mut spine) = spine_opt.take() else {
            return;
        };

        // L3019-3022
        let (ref_or1, ref_or2, ref_choix) = self.stripe_orientations(&spine);
        st.my_or1 = ref_or1;
        st.my_or2 = ref_or2;
        st.my_choix = ref_choix;

        // L3027-3046: ElSpine bookkeeping initialization.
        let mut intf = false;
        let mut intl = false;
        let mut current_he = super::chfi_ds::ChFiDSElSpine::new();
        let mut current_offset_he = super::chfi_ds::ChFiDSElSpine::new();
        let first_parameter = spine.base().first_parameter();
        let (pfirst, tfirst) = spine.base_mut().d1(first_parameter);
        current_he.firstparam = first_parameter;
        current_he.set_first_point_and_tgt(pfirst, tfirst);
        current_offset_he.firstparam = first_parameter;
        current_offset_he.set_first_point_and_tgt(pfirst, tfirst);

        let mut ya_k_part = false;
        let mut iedgelastkpart = 0usize;

        let mut w_start_periodic = 0.0f64;
        let nb_edges = spine.base().nb_edges();
        let mut w_end_periodic = spine.base().last_parameter_of(nb_edges);
        let (p_end_periodic, t_end_periodic) = spine.base_mut().d1(w_end_periodic);
        let mut wlast_book = 0.0f64;

        // L3050-3220: Construction of particular cases.
        for iedge in 1..=spine.base().nb_edges() {
            let (hs1, hs2) = self.conex_faces(&spine, iedge);

            if self.chfi3d_k_particular(&spine, iedge, &hs1, &hs2) {
                intf = iedge == 1 && !spine.base().is_periodic();
                intl = iedge == spine.base().nb_edges() && !spine.base().is_periodic();
                let or1_0 = hs1.orientation;
                let or2_0 = hs2.orientation;
                let (mut or1, mut or2) = (or1_0, or2_0);
                next_side(&mut or1, &mut or2, ref_or1, ref_or2, ref_choix);

                let mut sd = ChFiDSSurfData::default();
                let mut lsd: Vec<ChFiDSSurfData> = Vec::new();

                let compute_ok = chfi_kpart::compute_data_compute(
                    &self.my_brep,
                    self.my_ds.as_mut().expect("DS"),
                    &mut sd,
                    &hs1,
                    &hs2,
                    or1,
                    or2,
                    &spine,
                    iedge,
                );
                if !compute_ok {
                    // OCCT L3071-3072: empty else — the SD is dropped.
                } else if !self.split_k_part(&mut sd, &mut lsd, &hs1, &hs2, &mut spine, iedge, intf, intl) {
                    lsd.clear();
                } else {
                    iedgelastkpart = iedge;
                }

                // OCCT L3081-3120: periodic SD resorting — unreachable for
                // non-periodic spines; translated with the OCCT order.
                if spine.base().is_periodic() {
                    let nbsd = lsd.len();
                    let period = {
                        let n = spine.base().nb_edges();
                        spine.base().last_parameter_of(n) - spine.base().first_parameter()
                    };
                    let mut wfp = w_start_periodic;
                    let mut wlp = w_end_periodic;
                    if !ya_k_part && nbsd > 0 {
                        let first_sd = &lsd[0];
                        let mut wwf = first_sd.first_spine_param();
                        let mut wwl = first_sd.last_spine_param();
                        wwf = chfi_kpart::chfi_kpart_in_period(wwf, wfp, wlp, self.tolesp);
                        wwl = chfi_kpart::chfi_kpart_in_period(wwl, wfp, wlp, self.tolesp);
                        if wwl <= wwf + self.tolesp {
                            wwl += period;
                        }
                        wfp = wwf;
                        wlp = wfp + period;
                    }
                    let mut j = 0usize;
                    while j + 1 < nbsd {
                        let jwf = {
                            let jwf = lsd[j].first_spine_param();
                            chfi_kpart::chfi_kpart_in_period(jwf, wfp, wlp, self.tolesp)
                        };
                        for k in j + 1..nbsd {
                            let kwf = {
                                let kwf = lsd[k].first_spine_param();
                                chfi_kpart::chfi_kpart_in_period(kwf, wfp, wlp, self.tolesp)
                            };
                            if kwf < jwf {
                                lsd.swap(j, k);
                            }
                        }
                        j += 1;
                    }
                }

                // L3121-3214
                let mut li: Vec<i32> = Vec::new();
                for lsd_j in 0..lsd.len() {
                    let wfirst;
                    let wlast;
                    {
                        let cur_sd = &mut lsd[lsd_j];
                        if simul {
                            // OCCT: SimulKPart(curSD) — pending ChFiDS_CircSection.
                        }
                        wfirst = cur_sd.first_spine_param();
                        wlast = cur_sd.last_spine_param();
                    }
                    let cur_sd = lsd[lsd_j].clone();
                    // OCCT: SeqSurf.Append(curSD)
                    st.my_hdata.push(std::sync::Arc::new(cur_sd.clone()));
                    if !simul {
                        li.push(cur_sd.surf());
                    }
                    let (mut wfirst_j, mut wlast_j) = (wfirst, wlast);
                    if spine.base().is_periodic() {
                        wfirst_j = chfi_kpart::chfi_kpart_in_period(
                            wfirst_j,
                            w_start_periodic,
                            w_end_periodic,
                            self.tolesp,
                        );
                        wlast_j = chfi_kpart::chfi_kpart_in_period(
                            wlast_j,
                            w_start_periodic,
                            w_end_periodic,
                            self.tolesp,
                        );
                        if wlast_j <= wfirst_j + self.tolesp {
                            wlast_j += spine.base().first_parameter(); // Period() below
                        }
                    }
                    let (pfirst_j, tfirst_j) =
                        Self::tgt_kp(&lsd[lsd_j], &mut spine, iedge, true);
                    let (plast_j, tlast_j) =
                        Self::tgt_kp(&lsd[lsd_j], &mut spine, iedge, false);

                    // L3149-3213: Determine the sections to approximate.
                    if !ya_k_part {
                        if spine.base().is_periodic() {
                            w_start_periodic = wfirst_j;
                            w_end_periodic = w_start_periodic + spine.base().period();
                            wlast_j = elclib_in_period_static(
                                wlast_j,
                                w_start_periodic,
                                w_end_periodic,
                            );
                            if wlast_j <= wfirst_j + self.tolesp {
                                wlast_j += spine.base().period();
                            }
                            spine.base_mut().set_first_parameter(w_start_periodic);
                            spine.base_mut().set_last_parameter(w_end_periodic);
                        } else if !intf || iedge > 1 {
                            spine.base_mut().set_first_tgt(0.0f64.min(wfirst_j));
                            current_he.lastparam = wfirst_j;
                            current_he.set_last_point_and_tgt(pfirst_j, tfirst_j);
                            spine.base_mut().elspines.push(current_he.clone());
                            current_he.next = Some(std::sync::Arc::new(cur_sd.clone()));
                            current_he = super::chfi_ds::ChFiDSElSpine::new();

                            current_offset_he.lastparam = wfirst_j;
                            current_offset_he.set_last_point_and_tgt(pfirst_j, tfirst_j);
                            spine.base_mut().offset_elspines.push(current_offset_he.clone());
                            current_offset_he.next = Some(std::sync::Arc::new(cur_sd.clone()));
                            current_offset_he = super::chfi_ds::ChFiDSElSpine::new();
                        }
                        current_he.firstparam = wlast_j;
                        current_he.set_first_point_and_tgt(plast_j, tlast_j);
                        current_he.previous = Some(std::sync::Arc::new(cur_sd.clone()));
                        current_offset_he.firstparam = wlast_j;
                        current_offset_he.set_first_point_and_tgt(plast_j, tlast_j);
                        current_offset_he.previous = Some(std::sync::Arc::new(cur_sd.clone()));
                        ya_k_part = true;
                    } else if wfirst_j - current_he.firstparam > self.tolesp {
                        // section between two KPart
                        current_he.lastparam = wfirst_j;
                        current_he.set_last_point_and_tgt(pfirst_j, tfirst_j);
                        spine.base_mut().elspines.push(current_he.clone());
                        current_he.next = Some(std::sync::Arc::new(cur_sd.clone()));
                        current_he = super::chfi_ds::ChFiDSElSpine::new();

                        current_offset_he.lastparam = wfirst_j;
                        current_offset_he.set_last_point_and_tgt(pfirst_j, tfirst_j);
                        spine.base_mut().offset_elspines.push(current_offset_he.clone());
                        current_offset_he.next = Some(std::sync::Arc::new(cur_sd.clone()));
                        current_offset_he = super::chfi_ds::ChFiDSElSpine::new();

                        current_he.firstparam = wlast_j;
                        current_he.set_first_point_and_tgt(plast_j, tlast_j);
                        current_he.previous = Some(std::sync::Arc::new(cur_sd.clone()));
                        current_offset_he.firstparam = wlast_j;
                        current_offset_he.set_first_point_and_tgt(plast_j, tlast_j);
                        current_offset_he.previous = Some(std::sync::Arc::new(cur_sd.clone()));
                    } else {
                        current_he.firstparam = wlast_j;
                        current_he.set_first_point_and_tgt(plast_j, tlast_j);
                        current_he.previous = Some(std::sync::Arc::new(cur_sd.clone()));
                        current_offset_he.firstparam = wlast_j;
                        current_offset_he.set_first_point_and_tgt(plast_j, tlast_j);
                        current_offset_he.previous = Some(std::sync::Arc::new(cur_sd.clone()));
                    }
                    wlast_book = wlast_j;
                }
                if !li.is_empty() {
                    // OCCT L3217: myEVIMap.Bind(Spine->Edges(iedge), li)
                    self.my_evi_map
                        .insert(spine.base().edges(iedge).ptr_id(), li);
                }
            }
        }

        // L3222-3263: last section -> end of the spine.
        if !intl || iedgelastkpart < spine.base().nb_edges() {
            if spine.base().is_periodic() {
                if w_end_periodic - wlast_book > self.tolesp {
                    current_he.lastparam = w_end_periodic;
                    current_he.set_last_point_and_tgt(p_end_periodic, t_end_periodic);
                    if !ya_k_part {
                        current_he.periodic = true;
                    }
                    spine.base_mut().elspines.push(current_he.clone());

                    current_offset_he.lastparam = w_end_periodic;
                    current_offset_he.set_last_point_and_tgt(p_end_periodic, t_end_periodic);
                    if !ya_k_part {
                        current_offset_he.periodic = true;
                    }
                    spine.base_mut().offset_elspines.push(current_offset_he);
                }
            } else {
                let last_param = spine.base().last_parameter();
                let (plast, tlast) = spine.base_mut().d1(last_param);
                let n = spine.base().nb_edges();
                let w = spine.base().last_parameter_of(n).max(wlast_book);
                spine.base_mut().set_last_tgt(w);
                if spine.base().last_parameter() - wlast_book > self.tolesp {
                    current_he.lastparam = spine.base().last_parameter();
                    current_he.set_last_point_and_tgt(plast, tlast);
                    spine.base_mut().elspines.push(current_he.clone());

                    current_offset_he.lastparam = spine.base().last_parameter();
                    current_offset_he.set_last_point_and_tgt(plast, tlast);
                    spine.base_mut().offset_elspines.push(current_offset_he);
                }
            }
        }

        // L3265-3278: ChFi3d_PerformElSpine over the sections — pending
        // (Builder_0.cxx ChFi3d_PerformElSpine, needs the composite
        // curve approximation machinery).
        // L3279
        spine.base_mut().set_split_done(true);
        st.my_spine = Some(spine);
    }

    /// OCCT ChFi3d_Builder_2.cxx L3882-3900 — PerformSetOfSurf.
    pub fn perform_set_of_surf(&mut self, stripe: &SharedStripe, simul: bool) -> bool {
        // L3887-3888
        let si = {
            let st = stripe.read().expect("stripe lock");
            chfi3d_solid_index(&st, &mut self.my_ds.as_mut().expect("DS"), &self.my_eso_map, &self.my_esh_map)
        };
        {
            let mut st = stripe.write().expect("stripe lock");
            st.index_of_solid = si;
        }
        let split_done = {
            let st = stripe.read().expect("stripe lock");
            st.spine().map(|s| s.base().split_done()).unwrap_or(false)
        };
        if !split_done {
            self.perform_set_of_k_part(stripe, simul);
        }

        // L3894: PerformSetOfKGen — the numerical walking core
        // (Builder_2.cxx L3298-3882, BRepBlend/Extrema based) is a pending
        // translation; with a fully-KPart spine it processes no section.
        let _ = simul;

        // L3896-3899: ChFi3d_MakeExtremities — pending (Builder_6.cxx).
        true
    }
}

/// OCCT ChFi3d_Builder_0.cxx L2300-2328 — ChFi3d_SolidIndex.
fn chfi3d_solid_index(
    stripe: &ChFiDSStripe,
    dstr: &mut TopOpeBRepDSHDataStructure,
    map_eso: &ChFiDSMap,
    map_esh: &ChFiDSMap,
) -> i32 {
    let Some(sp) = stripe.spine() else {
        panic!("Standard_Failure: SolidIndex : Spine incomplete");
    };
    if sp.base().nb_edges() == 0 {
        panic!("Standard_Failure: SolidIndex : Spine incomplete");
    }
    let edref = sp.base().edges(1).clone();
    let eso_has = map_eso.contains(&edref) && !map_eso.find(&edref).is_empty();
    let shell_ou_solid = if eso_has {
        map_eso.find(&edref)[0].clone()
    } else if map_esh.contains(&edref) && !map_esh.find(&edref).is_empty() {
        map_esh.find(&edref)[0].clone()
    } else {
        panic!("Standard_Failure: SolidIndex : Spine incomplete");
    };
    dstr.add_shape(&shell_ou_solid)
}

/// OCCT ElCLib::InPeriod static stand-in (see chfi_ds::elclib_in_period).
fn elclib_in_period_static(u: f64, ufirst: f64, ulast: f64) -> f64 {
    super::chfi_ds::elclib_in_period(u, ufirst, ulast)
}



#[cfg(test)]
mod kpart_tests {
    use super::*;
    use super::super::brep_fillet_api::{explore_edges, explore_solids};
    use super::super::chfi_ds::ChFi3dFilletShape;

    /// OCCT BRepFilletAPI_MakeFillet::Build on a box: PerformSetOfSurf
    /// detects the plane-plane KPart, ChFiKPart_MakeFillet computes the
    /// cylinder SurfData analytically, and the corner machinery is the
    /// only pending stage.
    #[test]
    fn compute_box_edge_produces_kpart_surfdata() {
        let brep = rcad_modeling::make_box_brep(
            glam::DVec3::ZERO,
            glam::DVec3::X,
            glam::DVec3::Y,
            20.0,
            20.0,
            20.0,
        )
        .unwrap();
        let solid = explore_solids(&brep).into_iter().next().unwrap();
        let mut fillet =
            ChFi3dFilBuilder::new(&brep, solid, ChFi3dFilletShape::Rational, 1.0e-2);
        let edges = explore_edges(&brep);
        let e0 = edges[0].clone();

        fillet.add_radius(2.0, &e0);
        assert_eq!(fillet.base.nb_elements(), 1);

        fillet.base.compute();

        // The stripe must carry one KPart SurfData with a cylindrical
        // fillet surface of radius 2.
        let stripe = fillet.base.value_stripe(1);
        let st = stripe.read().expect("stripe lock");
        assert_eq!(st.my_hdata.len(), 1, "one KPart SurfData expected");
        let sd = &st.my_hdata[0];
        assert!(sd.surf_index >= 1, "surface registered in the DS");
        assert!(
            sd.index_of_s1 >= 1 && sd.index_of_s2 >= 1,
            "support faces registered in the DS (SpKP L870-871)"
        );
        let surf = &fillet.base.my_ds.as_ref().unwrap().surfaces[sd.surf_index as usize - 1].surface;
        assert!(
            matches!(surf, rcad_kernel::geom::Surface3::Cylinder(c) if (c.radius - 2.0).abs() < 1e-9),
            "KPart surface is the radius-2 cylinder, got {surf:?}"
        );
        assert!(
            sd.index_of_s1 >= 1 && sd.index_of_s2 >= 1,
            "support faces registered in the DS (SpKP L870-871)"
        );

        // PerformFilletOnVertex -> PerformSingularCorner: the stripe end
        // curves (the quarter-circle arcs) are computed and stored in the
        // DS per OCCT Builder.cxx L682-755.
        assert!(
            st.index_ofcurve1 > 0 && st.index_ofcurve2 > 0,
            "both singular end arcs created, got {}/{}",
            st.index_ofcurve1,
            st.index_ofcurve2
        );
        for ic in [st.index_ofcurve1, st.index_ofcurve2] {
            let curve = &fillet.base.my_ds.as_ref().unwrap().curves[ic as usize - 1].curve;
            let Some(curve) = curve else {
                panic!("singular end arc curve is null");
            };
            assert!(
                matches!(curve, rcad_kernel::geom::Curve3::Circle(c) if (c.radius - 2.0).abs() < 1e-9),
                "singular end arc is a radius-2 circle"
            );
        }
    }
}

impl ChFi3dBuilder {
    /// OCCT ChFi3d_Builder.cxx L682-755 — PerformSingularCorner: load
    /// vertex and degenerated edges.
    pub fn perform_singular_corner(&mut self, index: usize) {
        let vtx = self.my_vdata_map.find_key(index).clone();
        let stripes = self.my_vdata_map.find_from_index(index).clone();

        let mut ivtx = 0i32;
        for (i, stripe) in stripes.iter().enumerate() {
            // SurfData concerned and its CommonPoints.
            let mut sens = 0i32;
            let num = {
                let st = stripe.read().expect("stripe lock");
                chfi3d_index_of_surf_data(&vtx, &st, &mut sens)
            };
            let isfirst = sens == 1;

            // OCCT: Fd = stripe->SetOfSurfData()->Sequence().Value(num);
            let (cv1, cv2) = {
                let guard = stripe.read().expect("stripe lock");
                let Some(fd) = guard.my_hdata.get((num - 1) as usize) else {
                    continue;
                };
                (
                    fd.vertex(isfirst, 1).clone(),
                    fd.vertex(isfirst, 2).clone(),
                )
            };
            // Is it always degenerated?
            if !(cv1.point().distance(cv2.point()) <= 0.0) {
                continue;
            }
            // if yes the vertex is stored in the stripe and the edge at end
            // is created.
            if i == 0 {
                ivtx = chfi3d_index_point_in_ds(&cv1, self.my_ds.as_mut().expect("DS"));
            }

            // OCCT: VOnS1/VOnS2 = PCurveOnSurf()->Value(First/LastParameter).
            let (von_s1, von_s2) = {
                let (fi1, fi2) = (
                    stripe.read().expect("stripe lock").my_hdata[(num - 1) as usize].intf1.clone(),
                    stripe.read().expect("stripe lock").my_hdata[(num - 1) as usize].intf2.clone(),
                );
                if isfirst {
                    (
                        fi1.pcurve_on_surf().map(|pc| pc.point_at(fi1.parameter(true))),
                        fi2.pcurve_on_surf().map(|pc| pc.point_at(fi2.parameter(true))),
                    )
                } else {
                    (
                        fi1.pcurve_on_surf().map(|pc| pc.point_at(fi1.parameter(false))),
                        fi2.pcurve_on_surf().map(|pc| pc.point_at(fi2.parameter(false))),
                    )
                }
            };
            let (Some(von_s1), Some(von_s2)) = (von_s1, von_s2) else {
                continue;
            };

            // OCCT: ChFi3d_ComputeArete(CV1, VOnS1, CV2, VOnS2,
            // DStr.Surface(Fd->Surf()).Surface(), C3d, PCurv, Pardeb,
            // Parfin, tolapp3d, tolapp2d, tolreached, 0);
            let surf_index = stripe.read().expect("stripe lock").my_hdata[(num - 1) as usize].surf_index;
            let surf = self.my_ds.as_ref().expect("DS").surfaces[surf_index as usize - 1]
                .surface
                .clone();
            let (c3d, pcurv, pardeb, parfin, tolreached) = super::chfi3d_builder_0::chfi3d_compute_arete(
                &self.my_brep,
                &cv1,
                von_s1,
                &cv2,
                von_s2,
                &surf,
                self.tolapp3d,
                self.tolapp2d,
                0,
            );

            // OCCT: Crv = TopOpeBRepDS_Curve(C3d, tolreached); Icurv =
            // DStr.AddCurve(Crv);
            let Some(c3d) = c3d else {
                continue;
            };
            let icurv = self
                .my_ds
                .as_mut()
                .expect("DS")
                .add_curve(super::topopebrepds::TopOpeBRepDSCurve::new(Some(c3d), tolreached));

            let mut stw = stripe.write().expect("stripe lock");
            stw.set_curve(icurv, isfirst);
            stw.set_parameters(isfirst, pardeb, parfin);
            stw.change_pcurve(isfirst, pcurv);
            stw.set_index_point(ivtx, isfirst, 1);
            stw.set_index_point(ivtx, isfirst, 2);
        }
    }

    /// OCCT ChFi3d_Builder.cxx L759-920 - PerformFilletOnVertex.
    /// Returns false where the OCCT code would raise on a pending boundary.
    pub fn perform_fillet_on_vertex(&mut self, index: usize) -> bool {
        let vtx = self.my_vdata_map.find_key(index).clone();
        let stripes = self.my_vdata_map.find_from_index(index).clone();

        let mut i = 0usize;
        let mut nondegenere = true;
        let mut toujoursdegenere = true;
        let mut sp_free_boundary_at_end = false;
        for stripe in &stripes {
            let st = stripe.read().expect("stripe lock");
            let Some(sp) = st.spine() else {
                continue;
            };
            // SurfData and its CommonPoints.
            let mut sens = 0i32;
            let num = chfi3d_index_of_surf_data(&vtx, &st, &mut sens);
            let isfirst = sens == 1;
            let Some(fd) = st.my_hdata.get((num as usize).saturating_sub(1)) else {
                i += 1;
                continue;
            };
            let cv1 = fd.vertex(isfirst, 1);
            let cv2 = fd.vertex(isfirst, 2);
            // Is it always degenerated?
            if cv1.point().distance(cv2.point()) <= 0.0 {
                nondegenere = false;
            } else {
                toujoursdegenere = false;
            }
            sp_free_boundary_at_end = sp.base().status(isfirst) == ChFiDS_State::FreeBoundary;
            i += 1;
        }

        // calcul du nombre de faces = nombre d'aretes (sharp edges).
        let nba = super::chfi3d_builder_0::chfi3d_number_of_sharp_edges(
            &vtx,
            &self.my_ve_map,
            &self.my_ef_map,
            &self.my_brep,
        );

        if nondegenere {
            // Normal processing.  A false return encodes the OCCT raise on
            // a pending boundary.
            match i {
                1 => {
                    if sp_free_boundary_at_end {
                        return true;
                    }
                    if nba > 3 {
                        // OCCT: PerformIntersectionAtEnd(Index) — pending.
                        self.perform_intersection_at_end_pending();
                        return false;
                    } else if self.more_surfdata(index) {
                        // OCCT: PerformMoreSurfdata(Index) — pending.
                        self.perform_more_surfdata_pending();
                        return false;
                    }
                    // OCCT: PerformOneCorner(Index) — pending
                    // (ChFi3d_Builder_C1.cxx L611).
                    self.perform_one_corner_pending();
                    false
                }
                2 => {
                    if nba > 3 {
                        self.perform_more_three_corner_pending();
                    } else {
                        self.perform_two_corner_pending();
                    }
                    false
                }
                3 => {
                    if nba > 3 {
                        self.perform_more_three_corner_pending();
                    } else {
                        self.perform_three_corner_pending();
                    }
                    false
                }
                _ => {
                    self.perform_more_three_corner_pending();
                    false
                }
            }
        } else if toujoursdegenere {
            // Single case processing
            self.perform_singular_corner(index);
            true
        } else {
            // Last chance...
            self.perform_more_three_corner_pending();
            false
        }
    }

    /// OCCT ChFi3d_Builder_C1.cxx L4601-4640 — MoreSurfdata.
    pub fn more_surfdata(&self, index: usize) -> bool {
        // intersection at end is created on several surfdata if:
        // - the number of surfdata concerning the vertex is more than 1.
        // - and if the last but one surfdata has one of commonpoints on one
        //   of the two arcs, which constitute the intersections of the face
        //   at end and of the fillet.  The FindFace-dependent tail is
        //   pending; the stripe-count check is translated.
        let stripes = self.my_vdata_map.find_from_index(index);
        let Some(stripe) = stripes.first() else {
            return false;
        };
        let st = stripe.read().expect("stripe lock");
        st.my_hdata.len() > 1
    }

    fn perform_intersection_at_end_pending(&mut self) {
        // OCCT ChFi3d_Builder_C2.cxx PerformIntersectionAtEnd — pending.
    }

    fn perform_more_surfdata_pending(&mut self) {
        // OCCT ChFi3d_Builder_C1.cxx L3771 PerformMoreSurfdata — pending.
    }

    fn perform_one_corner_pending(&mut self) {
        // OCCT ChFi3d_Builder_C1.cxx L611 PerformOneCorner — pending.
    }

    fn perform_two_corner_pending(&mut self) {
        // OCCT ChFi3d_Builder_CnCrn.cxx PerformTwoCorner — pending.
    }

    fn perform_three_corner_pending(&mut self) {
        // OCCT ChFi3d_Builder_CnCrn.cxx PerformThreeCorner — pending.
    }

    fn perform_more_three_corner_pending(&mut self) {
        // OCCT ChFi3d_Builder_CnCrn.cxx PerformMoreThreeCorner — pending.
    }
}

/// OCCT ChFi3d_Builder_0.cxx L3446-3495 — ChFi3d_IndexOfSurfData.
pub fn chfi3d_index_of_surf_data(v1: &Shape, cd: &ChFiDSStripe, sens: &mut i32) -> i32 {
    let spine = cd.spine().expect("null spine");
    let mut index = 0i32;
    *sens = 1;
    let e = spine.base().edges(1).clone();
    let vref = if e.orientation == Orientation::Reversed {
        e.as_edge().expect("edge").last.clone()
    } else {
        e.as_edge().expect("edge").first.clone()
    };
    if vref.is_same(v1) {
        index = 1;
    } else {
        let e1 = spine.base().edges(spine.base().nb_edges()).clone();
        let vref = if e1.orientation == Orientation::Reversed {
            e1.as_edge().expect("edge").first.clone()
        } else {
            e1.as_edge().expect("edge").last.clone()
        };
        *sens = -1;
        if vref.is_same(v1) {
            index = cd.my_hdata.len() as i32;
        } else {
            panic!("Standard_ConstructionError: ChFi3d_IndexOfSurfData() - wrong construction parameters");
        }
    }
    index
}

/// OCCT ChFi3d_Builder_0.cxx L2345-2366 — ChFi3d_IndexPointInDS.
pub fn chfi3d_index_point_in_ds(
    p1: &super::chfi_ds::ChFiDS_CommonPoint,
    dstr: &mut TopOpeBRepDSHDataStructure,
) -> i32 {
    if p1.is_vertex() {
        let v = p1.vertex().clone();
        dstr.add_shape(&v)
    } else {
        dstr.add_point(super::topopebrepds::TopOpeBRepDSPoint::new(p1.point(), p1.tolerance()))
    }
}
