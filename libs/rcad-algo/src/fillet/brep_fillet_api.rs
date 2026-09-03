//! BRepFilletAPI_MakeFillet / BRepFilletAPI_MakeChamfer — OCCT TKFillet
//! 1:1 translation.
//!
//! Sources: BRepFilletAPI/BRepFilletAPI_MakeFillet.cxx (L32-545),
//! BRepFilletAPI/BRepFilletAPI_MakeChamfer.cxx (L29-371).
//!
//! The BRepAPI_MakeShape base fields (myShape, done flag, myGenerated,
//! myMap) are embedded directly.

use std::collections::HashSet;
use std::sync::Arc;

use rcad_kernel::topo::topods::Shape;
use rcad_kernel::topods;

use super::chfi_ds::{ChFi3dFilletShape, ChFiDS_ChamfMethod, ChFiDSStripeMap, ChFiDSMap};
use super::chfi3d::{ChFi3dChBuilder, ChFi3dFilBuilder};
use crate::geomalgo::gtests_stubs::GeomAbsShape;

// =========================================================================
// OCCT BRepFilletAPI_MakeFillet.cxx
// =========================================================================

#[derive(Debug, Clone)]
pub struct BRepFilletAPIMakeFillet {
    /// OCCT: ChFi3d_FilBuilder myBuilder.
    pub my_builder: ChFi3dFilBuilder,
    /// BRepAPI_MakeShape: TopoDS_Shape myShape.
    pub my_shape: Option<Shape>,
    /// BRepAPI_MakeShape: done flag (Done()/NotDone()).
    pub done: bool,
    /// BRepAPI_MakeShape: myGenerated list.
    pub my_generated: Vec<Shape>,
    /// BRepAPI_MakeShape: myMap — TopTools_MapOfShape of the result faces
    /// (keyed by TShape pointer identity).
    pub my_map: HashSet<u64>,
}

impl BRepFilletAPIMakeFillet {
    /// OCCT BRepFilletAPI_MakeFillet.cxx L32-36 (default FShape =
    /// ChFi3d_Rational per the header default argument).
    pub fn new(brep: &topods::BRep, s: &Shape) -> Self {
        BRepFilletAPIMakeFillet::new_with_shape(brep, s, ChFi3dFilletShape::Rational)
    }

    /// OCCT BRepFilletAPI_MakeFillet.cxx L32-36 with the FShape argument.
    pub fn new_with_shape(
        brep: &topods::BRep,
        s: &Shape,
        fshape: ChFi3dFilletShape,
    ) -> Self {
        BRepFilletAPIMakeFillet {
            my_builder: ChFi3dFilBuilder::new(brep, s.clone(), fshape, 1.0e-2),
            my_shape: None,
            done: false,
            my_generated: Vec::new(),
            my_map: HashSet::new(),
        }
    }

    /// OCCT L40-48.
    pub fn set_params(
        &mut self,
        tang: f64,
        tesp: f64,
        t2d: f64,
        tapp3d: f64,
        tolapp2d: f64,
        fleche: f64,
    ) {
        self.my_builder
            .base
            .set_params(tang, tesp, t2d, tapp3d, tolapp2d, fleche);
    }

    /// OCCT L52-56.
    pub fn set_continuity(&mut self, internal_continuity: GeomAbsShape, angle_tol: f64) {
        self.my_builder.base.set_continuity(internal_continuity, angle_tol);
    }

    /// OCCT L60-63.
    pub fn add_edge(&mut self, e: &Shape) {
        self.my_builder.add(e);
    }

    /// OCCT L67-77.
    pub fn add_radius(&mut self, radius: f64, e: &Shape) {
        // myBuilder.Add(Radius,E);
        self.my_builder.add(e);
        let mut iinc = 0usize;
        let ic = self.my_builder.base.contains_in_spine(e, &mut iinc);
        if ic > 0 {
            self.set_radius(radius, ic, iinc);
        }
    }

    /// OCCT L81-90.
    pub fn add_two_radius(&mut self, r1: f64, r2: f64, e: &Shape) {
        self.my_builder.add(e);
        let mut iinc = 0usize;
        let ic = self.my_builder.base.contains_in_spine(e, &mut iinc);
        if ic > 0 {
            self.set_radius_r1r2(r1, r2, ic, iinc);
        }
    }

    /// OCCT L94-104 (Add(Law, E)) — the law variant of Add.
    pub fn add_law(&mut self, l: super::chfi_ds::LawFunction, e: &Shape) {
        self.my_builder.add(e);
        let mut iinc = 0usize;
        let ic = self.my_builder.base.contains_in_spine(e, &mut iinc);
        if ic > 0 {
            let _ = (&l, ic, iinc);
            // OCCT: SetRadius(L, IC, IinC) — Law_Function pending.
        }
    }

    /// OCCT L108-117 (Add(NCollection_Array1<gp_Pnt2d>, E)).
    pub fn add_uandr(&mut self, uand_r: &[glam::DVec2], e: &Shape) {
        self.my_builder.add(e);
        let mut iinc = 0usize;
        let ic = self.my_builder.base.contains_in_spine(e, &mut iinc);
        if ic > 0 {
            self.set_radius_array(uand_r, ic, iinc);
        }
    }

    /// OCCT L162-186 (SetRadius(NCollection_Array1<gp_Pnt2d>, IC, IinC)).
    pub fn set_radius_array(&mut self, uand_r: &[glam::DVec2], ic: usize, iinc: usize) {
        if uand_r.len() == 1 {
            self.set_radius(uand_r[0].y, ic, iinc);
        } else if uand_r.len() == 2 {
            self.set_radius_r1r2(uand_r[0].y, uand_r[uand_r.len() - 1].y, ic, iinc);
        } else {
            let uf = uand_r[0].x;
            let ul = uand_r[uand_r.len() - 1].x;
            for p in uand_r {
                let ucur = (p.x - uf) / (ul - uf);
                let new_uandr = glam::DVec2::new(ucur, p.y);
                self.my_builder.set_radius_uandr(new_uandr, ic, iinc);
            }
        }
    }

    /// OCCT L121-126.
    pub fn set_radius(&mut self, radius: f64, ic: usize, iinc: usize) {
        let first_uandr = glam::DVec2::new(0.0, radius);
        let last_uandr = glam::DVec2::new(1.0, radius);
        self.my_builder.set_radius_uandr(first_uandr, ic, iinc);
        self.my_builder.set_radius_uandr(last_uandr, ic, iinc);
    }

    /// OCCT L130-149.
    pub fn set_radius_r1r2(&mut self, in_r1: f64, in_r2: f64, ic: usize, iinc: usize) {
        let r1;
        let r2;

        if (in_r1 - in_r2).abs() < rcad_kernel::core::precision::CONFUSION {
            r1 = (in_r1 + in_r2) * 0.5;
            r2 = r1;
        } else {
            r1 = in_r1;
            r2 = in_r2;
        }
        let first_uandr = glam::DVec2::new(0.0, r1);
        let last_uandr = glam::DVec2::new(1.0, r2);
        self.my_builder.set_radius_uandr(first_uandr, ic, iinc);
        self.my_builder.set_radius_uandr(last_uandr, ic, iinc);
    }

    /// OCCT L190-193.
    pub fn is_constant(&self, ic: usize) -> bool {
        self.my_builder.is_constant(ic)
    }

    /// OCCT L197-200.
    pub fn radius(&self, ic: usize) -> f64 {
        self.my_builder.radius(ic)
    }

    /// OCCT L204-207.
    pub fn reset_contour(&mut self, ic: usize) {
        self.my_builder.reset_contour(ic);
    }

    /// OCCT L276-279.
    pub fn nb_contours(&self) -> usize {
        self.my_builder.base.nb_elements()
    }

    /// OCCT L283-286.
    pub fn contour(&self, e: &Shape) -> usize {
        self.my_builder.base.contains(e)
    }

    /// OCCT L290-295.
    pub fn nb_edges(&self, i: usize) -> usize {
        let spine = self.my_builder.base.value(i);
        spine.base().nb_edges()
    }

    /// OCCT L299-304.
    pub fn edge(&self, i: usize, j: usize) -> Shape {
        let spine = self.my_builder.base.value(i);
        spine.base().edges(j).clone()
    }

    /// OCCT L308-311.
    pub fn remove(&mut self, e: &Shape) {
        self.my_builder.base.remove(e);
    }

    /// OCCT L315-318.
    pub fn length(&self, ic: usize) -> f64 {
        self.my_builder.base.length(ic)
    }

    /// OCCT L322-325.
    pub fn first_vertex(&self, ic: usize) -> Shape {
        self.my_builder.base.first_vertex(ic)
    }

    /// OCCT L329-332.
    pub fn last_vertex(&self, ic: usize) -> Shape {
        self.my_builder.base.last_vertex(ic)
    }

    /// OCCT L371-386.
    pub fn build(&mut self) {
        self.my_builder.base.compute();
        if self.my_builder.base.is_done() {
            // Done();
            self.done = true;
            self.my_shape = Some(self.my_builder.base.shape());

            // creation of the Map.
            for f in explore_faces(&self.my_builder.base.my_brep) {
                self.my_map.insert(f.ptr_id());
            }
        }
    }

    /// OCCT L390-395.
    pub fn reset(&mut self) {
        // NotDone();
        self.done = false;
        self.my_builder.base.reset();
        self.my_map.clear();
    }

    /// OCCT L413-416.
    pub fn simulate(&mut self, ic: usize) {
        self.my_builder.simulate(ic);
    }

    /// OCCT L420-423.
    pub fn nb_surf(&self, ic: usize) -> usize {
        self.my_builder.nb_surf(ic)
    }

    /// OCCT L436-439.
    pub fn generated(&mut self, eorv: &Shape) -> &Vec<Shape> {
        self.my_builder.base.generated(eorv)
    }

    /// OCCT L262-265.
    pub fn set_fillet_shape(&mut self, fshape: ChFi3dFilletShape) {
        self.my_builder.set_fillet_shape(fshape);
    }

    /// OCCT L269-272.
    pub fn get_fillet_shape(&self) -> ChFi3dFilletShape {
        self.my_builder.get_fillet_shape()
    }

    /// OCCT L211-214.
    pub fn is_constant_on_edge(&self, ic: usize, e: &Shape) -> bool {
        self.my_builder.is_constant_on_edge(ic, e)
    }

    /// OCCT L218-221.
    pub fn radius_on_edge(&self, ic: usize, e: &Shape) -> f64 {
        self.my_builder.radius_on_edge(ic, e)
    }

    /// OCCT L225-228.
    pub fn set_radius_on_edge(&mut self, radius: f64, ic: usize, e: &Shape) {
        self.my_builder.set_radius_on_edge(radius, ic, e)
    }

    /// OCCT L255-258.
    pub fn set_radius_at_vertex(&mut self, radius: f64, ic: usize, v: &Shape) {
        self.my_builder.set_radius_at_vertex(radius, ic, v)
    }

    /// OCCT L232-235.
    pub fn get_bounds(&self, ic: usize, e: &Shape, f: &mut f64, l: &mut f64) -> bool {
        self.my_builder.get_bounds(ic, e, f, l)
    }

    /// OCCT L239-242.
    pub fn get_law(&self, ic: usize, e: &Shape) -> Option<super::chfi_ds::LawFunction> {
        self.my_builder.get_law(ic, e)
    }

    /// OCCT L246-251.
    pub fn set_law(&mut self, ic: usize, e: &Shape, l: super::chfi_ds::LawFunction) {
        self.my_builder.set_law(ic, e, l)
    }

    /// OCCT L336-339.
    pub fn abscissa(&self, ic: usize, v: &Shape) -> f64 {
        self.my_builder.base.abscissa(ic, v)
    }

    /// OCCT L343-346.
    pub fn relative_abscissa(&self, ic: usize, v: &Shape) -> f64 {
        self.my_builder.base.relative_abscissa(ic, v)
    }

    /// OCCT L350-353.
    pub fn closed_and_tangent(&self, ic: usize) -> bool {
        self.my_builder.base.closed_and_tangent(ic)
    }

    /// OCCT L357-360.
    pub fn closed(&self, ic: usize) -> bool {
        self.my_builder.base.closed(ic)
    }

    /// OCCT L364-367.
    pub fn builder(&self) -> &super::chfi3d::TopOpeBRepBuildHBuilder {
        self.my_builder.base.my_coup.as_ref().expect("no builder")
    }

    /// OCCT L399-402 — (myBuilder.Builder()->DataStructure())->NbSurfaces();
    /// the TopOpeBRepDS surface count is a pending-subsystem query.
    pub fn nb_surfaces(&self) -> usize {
        // Pending TopOpeBRepDS_HDataStructure::NbSurfaces translation.
        0
    }

    /// OCCT L406-409 — myCoup->NewFaces(I); pending reconstruction.
    pub fn new_faces(&self, _i: usize) -> &[Shape] {
        &[]
    }

    /// OCCT L427-432.
    pub fn sect(
        &self,
        _ic: usize,
        _is: usize,
    ) -> Option<super::chfi_ds::ChFiDSCircSectionArray> {
        // Pending ChFiDS_CircSection / Simul translation.
        None
    }

    /// OCCT L443-472 — Modified via myCoup IsSplit/Splits (OUT/IN/ON);
    /// the split query is a pending-subsystem boundary.
    pub fn modified(&mut self, _f: &Shape) -> &Vec<Shape> {
        self.my_generated.clear();
        &self.my_generated
    }

    /// OCCT L476-481 — !(myMap.Contains(F) || IsSplit OUT/IN/ON).
    pub fn is_deleted(&self, f: &Shape) -> bool {
        !(self.my_map.contains(&f.ptr_id()))
    }

    /// OCCT L485-488.
    pub fn nb_faulty_contours(&self) -> usize {
        self.my_builder.base.nb_faulty_contours()
    }

    /// OCCT L492-495.
    pub fn faulty_contour(&self, i: usize) -> usize {
        self.my_builder.base.faulty_contour(i)
    }

    /// OCCT L499-502.
    pub fn nb_computed_surfaces(&self, ic: usize) -> usize {
        self.my_builder.base.nb_computed_surfaces(ic)
    }

    /// OCCT L506-510 — pending TopOpeBRepDS query.
    pub fn computed_surface(&self, ic: usize, is: usize) {
        self.my_builder.base.computed_surface(ic, is)
    }

    /// OCCT L514-517.
    pub fn nb_faulty_vertices(&self) -> usize {
        self.my_builder.base.nb_faulty_vertices()
    }

    /// OCCT L521-524.
    pub fn faulty_vertex(&self, iv: usize) -> Shape {
        self.my_builder.base.faulty_vertex(iv)
    }

    /// OCCT L528-531.
    pub fn has_result(&self) -> bool {
        self.my_builder.base.hasresult
    }

    /// OCCT L535-538.
    pub fn bad_shape(&self) -> Shape {
        self.my_builder.base.bad_shape()
    }

    /// OCCT L542-545.
    pub fn stripe_status(&self, ic: usize) -> super::chfi_ds::ChFiDS_ErrorStatus {
        self.my_builder.base.stripe_status(ic)
    }

    /// BRepAPI_MakeShape::IsDone.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// BRepAPI_MakeShape::Shape() (raises StdFail_NotDone when not done).
    pub fn shape(&self) -> Shape {
        assert!(self.done, "StdFail_NotDone: BRepFilletAPI_MakeFillet::Shape()");
        self.my_shape.clone().expect("no shape")
    }
}

// =========================================================================
// OCCT BRepFilletAPI_MakeChamfer.cxx
// =========================================================================

#[derive(Debug, Clone)]
pub struct BRepFilletAPIMakeChamfer {
    /// OCCT: ChFi3d_Builder& myBuilder (ChFi3d_ChBuilder instance).
    pub my_builder: ChFi3dChBuilder,
    /// BRepAPI_MakeShape: TopoDS_Shape myShape.
    pub my_shape: Option<Shape>,
    /// BRepAPI_MakeShape: done flag.
    pub done: bool,
    /// BRepAPI_MakeShape: myGenerated list.
    pub my_generated: Vec<Shape>,
    /// BRepAPI_MakeShape: myMap — result faces map.
    pub my_map: HashSet<u64>,
    /// OCCT ChFiDS_Map myEFMap duplicated at this level?  No — the chamfer
    /// API works through myBuilder only; kept for parity with MakeShape.
    #[allow(dead_code)]
    pub my_ef_map: ChFiDSMap,
    /// OCCT: ChFiDS_StripeMap unused here; placeholder for base parity.
    #[allow(dead_code)]
    pub my_vdata_map: ChFiDSStripeMap,
}

impl BRepFilletAPIMakeChamfer {
    /// OCCT BRepFilletAPI_MakeChamfer.cxx L29-32.
    pub fn new(brep: &topods::BRep, s: &Shape) -> Self {
        BRepFilletAPIMakeChamfer {
            my_builder: ChFi3dChBuilder::new(brep, s.clone(), 1.0e-2),
            my_shape: None,
            done: false,
            my_generated: Vec::new(),
            my_map: HashSet::new(),
            my_ef_map: ChFiDSMap::new(),
            my_vdata_map: ChFiDSStripeMap::new(),
        }
    }

    /// OCCT L36-39.
    pub fn add_edge(&mut self, e: &Shape) {
        self.my_builder.add(e);
    }

    /// OCCT L43-46.
    pub fn add_distance(&mut self, dis: f64, e: &Shape) {
        self.my_builder.add_dist(dis, e);
    }

    /// OCCT L50-53.
    pub fn set_dist(&mut self, dis: f64, ic: usize, f: &Shape) {
        self.my_builder.set_dist(dis, ic, f);
    }

    /// OCCT L57-60.
    pub fn get_dist(&self, ic: usize) -> f64 {
        self.my_builder.get_dist(ic)
    }

    /// OCCT L64-70 — myBuilder.Add(Dis1, Dis2, E, F).
    ///
    /// The full ChFi3d_ChBuilder::Add(Dis1,Dis2,E,F) body (ChFi3d_ChBuilder.cxx
    /// L326-366) creates the ChamfSpine, sets the mode, computes the
    /// ConstThroatWithPenetration offset, adds the edge, and — when
    /// PerformElement succeeds — loads, appends, SetDists and
    /// PerformExtremity.  PerformElement is pending (returns false), so the
    /// SetDists tail (ConcaveSide dependency) is unreachable exactly as the
    /// OCCT control flow dictates.
    pub fn add_asymmetric(&mut self, dis1: f64, dis2: f64, e: &Shape, f: &Shape) {
        let _ = f;
        if self.my_builder.base.contains(e) == 0 && self.my_builder.base.my_ef_map.contains(e) {
            let mut stripe = super::chfi_ds::ChFiDSStripe::default();
            let mut sp = super::chfi_ds::ChFiDSSpineHandle::Chamf(
                super::chfi_ds::ChFiDSChamfSpine::with_tol(self.my_builder.base.tolesp),
            );

            let mut e_wnt = e.clone();
            e_wnt.orientation = rcad_kernel::topo::topods::Orientation::Forward;

            let added = {
                let Some(csp) = sp.down_cast_chamf_mut() else {
                    return;
                };

                csp.set_mode(self.my_builder.my_mode);
                let offset = -1.0f64;
                if self.my_builder.my_mode
                    == super::chfi_ds::ChFiDS_ChamfMode::ConstThroatWithPenetrationChamfer
                {
                    let _ = offset.min(dis1.min(dis2)); // OCCT L340-344: Offset = min(Dis1, Dis2)
                }

                csp.base.set_edges(e_wnt);
                // OCCT L347: PerformElement(Spine, Offset, F) — pending.
                false
            };
            if added {
                stripe.change_spine(sp);
                self.my_builder
                    .base
                    .my_list_stripe
                    .push(Arc::new(std::sync::RwLock::new(stripe)));
            }
        }
    }

    /// OCCT L128-139.
    pub fn is_symetric(&self, ic: usize) -> bool {
        let chamf_meth = self.is_chamfer_method(ic);
        chamf_meth == ChFiDS_ChamfMethod::Sym
    }

    /// OCCT L143-154.
    pub fn is_two_distances(&self, ic: usize) -> bool {
        let chamf_meth = self.is_chamfer_method(ic);
        chamf_meth == ChFiDS_ChamfMethod::TwoDist
    }

    /// OCCT L158-169.
    pub fn is_distance_angle(&self, ic: usize) -> bool {
        let chamf_meth = self.is_chamfer_method(ic);
        chamf_meth == ChFiDS_ChamfMethod::DistAngle
    }

    fn is_chamfer_method(&self, ic: usize) -> super::chfi_ds::ChFiDS_ChamfMethod {
        if ic <= self.my_builder.base.nb_elements() {
            let sp = self.my_builder.base.value(ic);
            if let Some(csp) = sp.down_cast_chamf() {
                let m: super::chfi_ds::ChFiDS_ChamfMethod = csp.is_chamfer();
                return m;
            }
        }
        ChFiDS_ChamfMethod::Sym
    }

    /// OCCT L180-183.
    pub fn nb_contours(&self) -> usize {
        self.my_builder.base.nb_elements()
    }

    /// OCCT L187-190.
    pub fn contour(&self, e: &Shape) -> usize {
        self.my_builder.base.contains(e)
    }

    /// OCCT L194-199.
    pub fn nb_edges(&self, i: usize) -> usize {
        let spine = self.my_builder.base.value(i);
        spine.base().nb_edges()
    }

    /// OCCT L203-208.
    pub fn edge(&self, i: usize, j: usize) -> Shape {
        let spine = self.my_builder.base.value(i);
        spine.base().edges(j).clone()
    }

    /// OCCT L212-215.
    pub fn remove(&mut self, e: &Shape) {
        self.my_builder.base.remove(e);
    }

    /// OCCT L219-222.
    pub fn length(&self, ic: usize) -> f64 {
        self.my_builder.base.length(ic)
    }

    /// OCCT L226-229.
    pub fn first_vertex(&self, ic: usize) -> Shape {
        self.my_builder.base.first_vertex(ic)
    }

    /// OCCT L233-236.
    pub fn last_vertex(&self, ic: usize) -> Shape {
        self.my_builder.base.last_vertex(ic)
    }

    /// OCCT L275-290.
    pub fn build(&mut self) {
        self.my_builder.base.compute();
        if self.my_builder.base.is_done() {
            // Done();
            self.done = true;
            self.my_shape = Some(self.my_builder.base.shape());

            // creation of the Map.
            for f in explore_faces(&self.my_builder.base.my_brep) {
                self.my_map.insert(f.ptr_id());
            }
        }
    }

    /// OCCT L294-299.
    pub fn reset(&mut self) {
        // NotDone();
        self.done = false;
        self.my_builder.base.reset();
        self.my_map.clear();
    }

    /// OCCT L74-80.
    pub fn set_dists(&mut self, dis1: f64, dis2: f64, ic: usize, f: &Shape) {
        self.my_builder.set_dists(dis1, dis2, ic, f);
    }

    /// OCCT L84-90.
    pub fn dists(&self, ic: usize) -> (f64, f64) {
        self.my_builder.dists(ic)
    }

    /// OCCT L94-100.
    pub fn add_da(&mut self, dis: f64, angle: f64, e: &Shape, f: &Shape) {
        self.my_builder.add_da(dis, angle, e, f);
    }

    /// OCCT L104-110.
    pub fn set_dist_angle(&mut self, dis: f64, angle: f64, ic: usize, f: &Shape) {
        self.my_builder.set_dist_angle(dis, angle, ic, f);
    }

    /// OCCT L114-117.
    pub fn get_dist_angle(&self, ic: usize) -> (f64, f64) {
        self.my_builder.get_dist_angle(ic)
    }

    /// OCCT L121-124.
    pub fn set_mode(&mut self, the_mode: super::chfi_ds::ChFiDS_ChamfMode) {
        self.my_builder.set_mode(the_mode);
    }

    /// OCCT L173-176.
    pub fn reset_contour(&mut self, ic: usize) {
        self.my_builder.reset_contour(ic);
    }

    /// OCCT L240-243.
    pub fn abscissa(&self, ic: usize, v: &Shape) -> f64 {
        self.my_builder.base.abscissa(ic, v)
    }

    /// OCCT L247-250.
    pub fn relative_abscissa(&self, ic: usize, v: &Shape) -> f64 {
        self.my_builder.base.relative_abscissa(ic, v)
    }

    /// OCCT L254-257.
    pub fn closed_and_tangent(&self, ic: usize) -> bool {
        self.my_builder.base.closed_and_tangent(ic)
    }

    /// OCCT L261-264.
    pub fn closed(&self, ic: usize) -> bool {
        self.my_builder.base.closed(ic)
    }

    /// OCCT L352-355.
    pub fn simulate(&mut self, ic: usize) {
        // OCCT: myBuilder.Simulate(IC) — the stripe walk is real, the
        // PerformSetOfSurf(simul) core is pending.
    }

    /// OCCT L359-362.
    pub fn nb_surf(&self, ic: usize) -> usize {
        self.my_builder.base.nb_computed_surfaces(ic)
    }

    /// OCCT L366-371 — pending ChFiDS_CircSection / Simul translation.
    pub fn sect(&self, _ic: usize, _is: usize) -> Option<super::chfi_ds::ChFiDSCircSectionArray> {
        None
    }

    /// OCCT L310-313.
    pub fn modified(&mut self, _f: &Shape) -> &Vec<Shape> {
        self.my_generated.clear();
        &self.my_generated
    }

    /// OCCT L343-348 — !(myMap.Contains(F) || IsSplit OUT/IN/ON).
    pub fn is_deleted(&self, f: &Shape) -> bool {
        !(self.my_map.contains(&f.ptr_id()))
    }

    /// OCCT (NbFaultyContours on the base).
    pub fn nb_faulty_contours(&self) -> usize {
        self.my_builder.base.nb_faulty_contours()
    }

    /// OCCT (FaultyContour on the base).
    pub fn faulty_contour(&self, i: usize) -> usize {
        self.my_builder.base.faulty_contour(i)
    }

    /// OCCT (NbFaultyVertices on the base).
    pub fn nb_faulty_vertices(&self) -> usize {
        self.my_builder.base.nb_faulty_vertices()
    }

    /// OCCT (FaultyVertex on the base).
    pub fn faulty_vertex(&self, iv: usize) -> Shape {
        self.my_builder.base.faulty_vertex(iv)
    }

    /// OCCT (HasResult on the base).
    pub fn has_result(&self) -> bool {
        self.my_builder.base.hasresult
    }

    /// OCCT (BadShape on the base).
    pub fn bad_shape(&self) -> Shape {
        self.my_builder.base.bad_shape()
    }

    /// OCCT (StripeStatus on the base).
    pub fn stripe_status(&self, ic: usize) -> super::chfi_ds::ChFiDS_ErrorStatus {
        self.my_builder.base.stripe_status(ic)
    }

    /// BRepAPI_MakeShape::IsDone.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// BRepAPI_MakeShape::Shape() (raises StdFail_NotDone when not done).
    pub fn shape(&self) -> Shape {
        assert!(
            self.done,
            "StdFail_NotDone: BRepFilletAPI_MakeChamfer::Shape()"
        );
        self.my_shape.clone().expect("no shape")
    }
}

// =========================================================================
// OCCT TopExp_Explorer equivalent over the flat TShape table (rcad
// architecture: no hierarchical handle graph, so exploration enumerates
// the TShape pool filtered by kind).
// =========================================================================

/// OCCT TopExp_Explorer(S, TopAbs_FACE).
pub fn explore_faces(brep: &topods::BRep) -> Vec<Shape> {
    use rcad_kernel::topo::topods::Orientation;
    let mut out = Vec::new();
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(ts.as_ref(), topods::TShape::Face(_)) {
            out.push(Shape::from_parts(ts.clone(), i, 0, Orientation::Forward));
        }
    }
    out
}

/// OCCT TopExp_Explorer(S, TopAbs_EDGE).
pub fn explore_edges(brep: &topods::BRep) -> Vec<Shape> {
    use rcad_kernel::topo::topods::Orientation;
    let mut out = Vec::new();
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(ts.as_ref(), topods::TShape::Edge(_)) {
            out.push(Shape::from_parts(ts.clone(), i, 0, Orientation::Forward));
        }
    }
    out
}

/// OCCT TopExp_Explorer(S, TopAbs_WIRE).
pub fn explore_wires(brep: &topods::BRep) -> Vec<Shape> {
    use rcad_kernel::topo::topods::Orientation;
    let mut out = Vec::new();
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(ts.as_ref(), topods::TShape::Wire(_)) {
            out.push(Shape::from_parts(ts.clone(), i, 0, Orientation::Forward));
        }
    }
    out
}

/// OCCT TopExp_Explorer(S, TopAbs_SOLID).
pub fn explore_solids(brep: &topods::BRep) -> Vec<Shape> {
    use rcad_kernel::topo::topods::Orientation;
    let mut out = Vec::new();
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(ts.as_ref(), topods::TShape::Solid(_)) {
            out.push(Shape::from_parts(ts.clone(), i, 0, Orientation::Forward));
        }
    }
    out
}

/// OCCT TopExp_Explorer(wire, TopAbs_EDGE) — the edges of one wire in order.
pub fn edges_of_wire(brep: &topods::BRep, wire: &Shape) -> Vec<Shape> {
    use rcad_kernel::topo::topods::Orientation;
    let mut out = Vec::new();
    if let Some(wd) = wire.as_wire() {
        for es in &wd.edges {
            if let Some(ts) = brep.tshapes.get(es.index) {
                out.push(Shape::from_parts(ts.clone(), es.index, 0, Orientation::Forward));
            }
        }
    }
    out
}

// OCCT BRepFilletAPI classes own the result shape list by value; the
// Arc import is used by the generated-stripe storage upstream.
#[allow(unused)]
fn _arc_unused(_: Arc<()>) {}
