//! OCCT ChFi3d_Builder Perform-flow completion (Stage 1f round 2) — the
//! untranslated main-flow segments of ChFi3d_Builder.cxx /
//! ChFi3d_FilBuilder.cxx, translated 1:1.
//!
//! Sources:
//!   - ChFi3d_Builder.cxx L471-578 (Compute tail: myCoup reconstruction,
//!     split-edge tolerance pass, myShapeResult / badShape assembly)
//!   - ChFi3d_Builder.cxx L647-675 (SameParameter pass over the new faces)
//!   - ChFi3d_FilBuilder.cxx L455-471 (Sect), L473-560 (SimulKPart),
//!     L2595-2637 (SetRegul)
//!   - TopOpeBRepBuild_HBuilder.hxx L52-111 (the myCoup call surface)
//!
//! chfi3d.rs is read-only this round (existing 3,746-line alignment work),
//! so this file carries the new `impl ChFi3dBuilder` /
//! `impl ChFi3dFilBuilder` blocks plus the inherent
//! `impl TopOpeBRepBuildHBuilder` block.  Rust resolves methods across
//! multiple inherent impl blocks in one crate, so the blocks below are
//! equivalent to living in chfi3d.rs.
//!
//! =========================================================================
//! Stage 1g wiring requirement list (the mechanical splices that connect
//! the translations below to their OCCT call points; the target files are
//! outside this round's ownership, so the exact edits are recorded here):
//!
//!   W1. chfi3d.rs L3351-3357 (ChFi3d_Builder::PerformSetOfSurf; OCCT
//!       ChFi3d_Builder_2.cxx L3894 and L3896-3899): replace the two
//!       pending comments and `let _ = simul;` with
//!             self.perform_set_of_k_gen(stripe, simul);            // L3894
//!             if !simul {                                          // L3896
//!                 let mut stg = stripe.write().expect("stripe lock");
//!                 super::chfi3d_builder_2c::chfi3d_make_extremities(
//!                     &self.my_brep,
//!                     &mut stg,
//!                     self.my_ds.as_mut().expect("DS"),            // L3898
//!                     &self.my_ef_map,
//!                     self.tolapp3d,
//!                     self.tol2d,
//!                 );
//!             }
//!
//!   W2. chfi3d.rs L575-580 (Compute tail): replace
//!       `self.done = false; // pending: TopOpeBRepBuild
//!       Perform/MergeSolid/NewFaces` with
//!             self.perform_hbuilder_reconstruction(&map_ind_so);   // L471-578
//!       (`map_ind_so` is the Vec<i32> local declared at chfi3d.rs L519).
//!
//!   W3. chfi3d.rs compute tail, directly after W2: OCCT L574 `SetRegul();`
//!       dispatches virtually to the derived override.  Over the rcad
//!       composition model the splice must issue the derived call
//!       explicitly — `fillet_builder.set_regul()` (this file, OCCT
//!       ChFi3d_FilBuilder.cxx L2595-2637) or the ChBuilder
//!       `set_regul()` (chfi3d_builder_chbuilder.rs L1650, OCCT
//!       ChFi3d_ChBuilder.cxx L2218-2276).
//!
//!   W4. chfi3d.rs L584-586 (Compute tail, inside the IsDone gate; OCCT
//!       L647-675): replace the pending comment with
//!             self.same_parameter_pass();
//!
//!   W5. chfi3d.rs L726 (Generated; OCCT ChFi3d_Builder.cxx L968): replace
//!       the `let _ = i;` stub with
//!             let coup = self.my_coup.as_ref().expect("HBuilder");
//!             self.my_generated = coup.new_faces(i);
//!
//!   W6. chfi3d.rs L3179-3181 (PerformSetOfKPart simul store; OCCT
//!       ChFi3d_Builder_2.cxx L3123-3125 `if (Simul) { SimulKPart(curSD); }`):
//!       replace the pending comment with the derived dispatch
//!       `simul_kpart(cur_sd)` (this file for the fillet builder,
//!       chfi3d_builder_chbuilder.rs L390 for the chamfer builder).  As in
//!       W3 the base KPart loop cannot dispatch over composition; the
//!       hook is a 1g architecture decision (tag on the base or move of
//!       the loop).
//!
//!   W7. chfi3d.rs L989-997 (ChFi3d_FilBuilder::Simulate; OCCT
//!       ChFi3d_FilBuilder.cxx L418-421): replace the pending comment with
//!             self.base.perform_set_of_surf(stripe, true);
//!
//!   W8. chfi3d_builder_2b.rs L1469 / L1488 (Builder_2.cxx L3612-3648):
//!       replace `complete_data_pending(...)` with
//!       `self.complete_data_surfcoin(...)` (chfi3d_builder_6.rs L674; the
//!       OCCT Builder.hxx L689 overload taking Newsurf + S1/PC1 + S2/PC2),
//!       passing `pc1.as_ref()` / `pc2.as_ref()` for the pcurve slots.
//!
//!   W9. chfi3d_builder_chbuilder.rs: L471-472 replace the SetSimul pending
//!       note with `sd.set_simul(sec);` (the simul slot now exists on
//!       ChFiDSSurfData); L376-378 (Sect down-cast) return the slot read
//!       `guard.simul()`; L273 `store_simul` body becomes
//!       `_data.set_simul(Some(_sec));`.
//!
//!   W10. brep_fillet_api.rs L349-355 / L694 (the Sect facade): route to
//!       `ChFi3dFilBuilder::sect` (this file) instead of returning None.
//!
//!   W11. hbuilder.rs: move the `impl TopOpeBRepBuildHBuilder` block below
//!       into fillet/hbuilder.rs and give the 7 methods their TKBO-backed
//!       bodies (the module doc there already reserves this surface).
//!       Until then the bodies are pending-boundary neutrals.
//!
//!   N1. Still-untranslated Perform-flow stages (future translation, not
//!       wiring): the corner tails behind chfi3d.rs L3698-3715 —
//!       PerformIntersectionAtEnd (ChFi3d_Builder_C2.cxx),
//!       PerformMoreSurfdata (ChFi3d_Builder_C1.cxx L3771),
//!       PerformTwoCorner / PerformThreeCorner / PerformMoreThreeCorner
//!       (ChFi3d_Builder_CnCrn.cxx) — plus PerformFilletOnVertex's
//!       non-degenerate corner paths and Compute's L144-174 ExtentAnalyse.
//! =========================================================================

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve2dEval as _, Surface3};
use rcad_kernel::topo::topods::{BRepBuilder, BRepTool as _};
use rcad_kernel::topods::{self, Shape};

use super::chfi3d::{ChFi3dBuilder, ChFi3dFilBuilder, TopOpeBRepDSHDataStructure};
use super::chfi3d_builder_0::{surface_type_of, GeomAbsSurfaceType};
use super::chfi3d_builder_2::TopAbsState;
use super::chfi_ds::{ChFiDSCircSection, ChFiDSCircSectionArray, ChFiDSSurfData};
use super::chfi_kpart_gp::{surface3_ax3, GpAx3, GpCirc};
use super::hbuilder::TopOpeBRepBuildHBuilder;

// =========================================================================
// OCCT TopOpeBRepBuild_HBuilder — the myCoup call surface consumed by the
// ChFi3d Perform flow.  The inherent impl block lives in this file per the
// Stage 1f ownership split (hbuilder.rs itself is read-only this round);
// Stage 1g moves it into hbuilder.rs with TKBO-backed bodies (item W11).
// Until then every body is the pending-boundary neutral: the reconstruction
// below keeps the exact OCCT flow shape, and the merges report "nothing
// merged" so the result assembly degrades to the identity compound.
// =========================================================================

impl TopOpeBRepBuildHBuilder {
    /// OCCT TopOpeBRepBuild_HBuilder.hxx L52 — Perform(HDS).
    pub fn perform(&mut self, _ds: &mut TopOpeBRepDSHDataStructure) {}

    /// OCCT TopOpeBRepBuild_HBuilder.hxx L85 — MergeSolid(S, TB).
    pub fn merge_solid(&mut self, _s: &Shape, _tb: TopAbsState) {}

    /// OCCT TopOpeBRepBuild_HBuilder.hxx L88 — IsSplit(S, ToBuild).
    pub fn is_split(&self, _s: &Shape, _to_build: TopAbsState) -> bool {
        false
    }

    /// OCCT TopOpeBRepBuild_HBuilder.hxx L91 — Splits(S, ToBuild).
    pub fn splits(&self, _s: &Shape, _to_build: TopAbsState) -> Vec<Shape> {
        Vec::new()
    }

    /// OCCT TopOpeBRepBuild_HBuilder.hxx L98 — Merged(S, ToBuild).
    pub fn merged(&self, _s: &Shape, _to_build: TopAbsState) -> Vec<Shape> {
        Vec::new()
    }

    /// OCCT TopOpeBRepBuild_HBuilder.hxx L105 — NewEdges(I).
    pub fn new_edges(&self, _i: i32) -> Vec<Shape> {
        Vec::new()
    }

    /// OCCT TopOpeBRepBuild_HBuilder.hxx L111 — NewFaces(I).
    pub fn new_faces(&self, _i: i32) -> Vec<Shape> {
        Vec::new()
    }
}

impl ChFi3dBuilder {
    /// OCCT ChFi3d_Builder.cxx L471-578 — the reconstruction tail of
    /// Compute: myCoup->Perform, the MergeSolid pass over MapIndSo, the
    /// split-edge tolerance pass, the myShapeResult / badShape assembly
    /// and the SetRegul call point.
    ///
    /// `map_ind_so` is OCCT's MapIndSo (the DS indices of myShape's
    /// solids, collected at Compute L339-346 — chfi3d.rs L519-527).
    pub fn perform_hbuilder_reconstruction(&mut self, map_ind_so: &[i32]) {
        // L471: myCoup->Perform(myDS);
        {
            let coup = self.my_coup.as_mut().expect("HBuilder");
            let dstr = self.my_ds.as_mut().expect("DS");
            coup.perform(dstr);
        }

        // L472-478: MergeSolid(curshape, TopAbs_IN) over MapIndSo.
        for indsol in map_ind_so {
            let curshape = {
                let dstr = self.my_ds.as_ref().expect("DS");
                dstr.shape(*indsol).clone()
            };
            let coup = self.my_coup.as_mut().expect("HBuilder");
            coup.merge_solid(&curshape, TopAbsState::In);
        }

        // L405 (scoping note): OCCT declares `BRep_Builder B1` once before
        // CompleteDS and reuses it in the L449 tolerance write (translated
        // in chfi3d.rs compute_tolerance_pass) and here; rcad keeps the
        // builder per method.
        let mut b1 = BRepBuilder::new();

        // L480-509: tolerance pass over the split edges of the DS.
        let n = {
            let dstr = self.my_ds.as_ref().expect("DS");
            dstr.nb_shapes()
        };
        for i in 1..=n {
            let s = {
                let dstr = self.my_ds.as_ref().expect("DS");
                dstr.shape(i).clone()
            };
            if s.shape_type() != topods::ShapeType::Edge {
                continue;
            }
            // L488: bool issplitIN = myCoup->IsSplit(S, TopAbs_IN);
            let issplit_in = {
                let coup = self.my_coup.as_ref().expect("HBuilder");
                coup.is_split(&s, TopAbsState::In)
            };
            if !issplit_in {
                continue;
            }
            // L493: for (it = myCoup->Splits(S, TopAbs_IN); it.More(); it.Next())
            let splits_list = {
                let coup = self.my_coup.as_ref().expect("HBuilder");
                coup.splits(&s, TopAbsState::In)
            };
            for new_e in splits_list {
                // L496: const TopoDS_Edge& newE = TopoDS::Edge(it.Value());
                // L497: double tole = BRep_Tool::Tolerance(newE);
                let tole = self.my_brep.tolerance(&new_e);
                // L498-499: TopExp_Explorer exv(newE, TopAbs_VERTEX).
                // Architecture: the rcad edge carries its two terminal
                // vertices (TopExp::FirstVertex / LastVertex), so the
                // explorer reduces to the pair.
                for v in [
                    self.my_brep.first_vertex(&new_e),
                    self.my_brep.last_vertex(&new_e),
                ] {
                    let tolv = self.my_brep.tolerance(&v);
                    if tole > tolv {
                        // L505: B1.UpdateVertex(v, tole);
                        b1.update_vertex_tolerance(&mut self.my_brep, v, tole);
                    }
                }
            }
        }

        // L510-567: myShapeResult / badShape assembly.
        if !self.hasresult {
            // L512: B1.MakeCompound(TopoDS::Compound(myShapeResult));
            // Architecture (A1): BRepBuilder::make_compound builds from the
            // accumulated list; the OCCT incremental B1.Add sequence maps
            // to the Vec pushes below.
            let mut added: Vec<Shape> = Vec::new();
            {
                let coup = self.my_coup.as_ref().expect("HBuilder");
                let dstr = self.my_ds.as_ref().expect("DS");
                for indsol in map_ind_so {
                    // L515-516
                    let curshape = dstr.shape(*indsol).clone();
                    // L517: its = myCoup->Merged(curshape, TopAbs_IN);
                    let its = coup.merged(&curshape, TopAbsState::In);
                    if its.is_empty() {
                        // L520: B1.Add(myShapeResult, curshape);
                        added.push(curshape);
                    } else {
                        // L524-541: if the old type of Shape is Shell, the
                        // Shell is placed instead of Solid.
                        for value in its {
                            let letype = curshape.shape_type();
                            if letype == topods::ShapeType::Shell {
                                // L531-533: TopExp_Explorer
                                // expsh2(its.Value(), TopAbs_SHELL);
                                // const TopoDS_Shape& cursh =
                                //     expsh2.Current(); B1.Add(cursh).
                                // Architecture: rcad has no per-shape
                                // sub-explorer yet (the flat TShape pool
                                // carries no parent links); in the
                                // HBuilder-backed flow the merged value of
                                // a shell is already the shell, so the
                                // explorer's first hit reduces to the value
                                // itself.  Revisit with W11.
                                added.push(value);
                            } else {
                                // L538: B1.Add(myShapeResult, its.Value());
                                added.push(value);
                            }
                        }
                    }
                }
            }
            self.my_shape_result = Some(b1.make_compound(&mut self.my_brep, added));
        } else {
            // L547: done = false;
            self.done = false;
            // L548: B1.MakeCompound(TopoDS::Compound(badShape));
            let mut added: Vec<Shape> = Vec::new();
            {
                let coup = self.my_coup.as_ref().expect("HBuilder");
                let dstr = self.my_ds.as_ref().expect("DS");
                for indsol in map_ind_so {
                    let curshape = dstr.shape(*indsol).clone();
                    let its = coup.merged(&curshape, TopAbsState::In);
                    if its.is_empty() {
                        // L556: B1.Add(badShape, curshape);
                        added.push(curshape);
                    } else {
                        for value in its {
                            // L562: B1.Add(badShape, its.Value());
                            added.push(value);
                        }
                    }
                }
            }
            self.bad_shape = Some(b1.make_compound(&mut self.my_brep, added));
        }

        // L574: SetRegul(); — the OCCT call dispatches virtually to the
        // derived override (ChFi3d_FilBuilder.cxx L2595 = `set_regul` on
        // ChFi3dFilBuilder in this file; ChFi3d_ChBuilder.cxx L2218 =
        // chfi3d_builder_chbuilder.rs L1650).  Composition has no virtual
        // dispatch from the base; the compute-tail splice (W2/W3) issues
        // the derived call after this method returns.
    }

    /// OCCT ChFi3d_Builder.cxx L647-675 — the SameParameter pass over the
    /// new faces (the OCCT body is gated by `if (IsDone())`; the compute
    /// tail keeps that gate at the W4 splice point).
    pub fn same_parameter_pass(&mut self) {
        // L649: double SameParTol = Precision::Confusion();
        let same_par_tol = CONFUSION;
        // L656: aNbSurfaces = myDS->NbSurfaces();
        // Architecture: the rcad DS keeps its surface table as
        // `side.surfaces`; TopOpeBRepDS_HDataStructure::NbSurfaces maps to
        // its length.
        let a_nb_surfaces = {
            let dstr = self.my_ds.as_ref().expect("DS");
            dstr.side.surfaces.len() as i32
        };
        // L657-664: for each DS surface, the new faces of myCoup.
        for i_f in 1..=a_nb_surfaces {
            // L661: const NCollection_List<TopoDS_Shape>& aLF = myCoup->NewFaces(iF);
            let a_lf = {
                let coup = self.my_coup.as_ref().expect("HBuilder");
                coup.new_faces(i_f)
            };
            for a_f in a_lf {
                // L668: BRepLib::SameParameter(aF, SameParTol, true);
                // L669: ShapeFix::SameParameter(aF, false, SameParTol);
                // Pending boundary: BRepLib::SameParameter and
                // ShapeFix::SameParameter have no rcad carrier yet; the
                // loop keeps the OCCT iteration shape (same neutralization
                // style as chfi3d_builder_chbuilder.rs L1694-1697).
                let _ = (&a_f, same_par_tol);
            }
        }
    }
}

impl ChFi3dFilBuilder {
    /// OCCT ChFi3d_FilBuilder.cxx L455-471 — Sect(IC, IS): the simulated
    /// sections of the IS-th SurfData of contour IC.
    pub fn sect(&self, ic: usize, is: usize) -> Option<ChFiDSCircSectionArray> {
        let mut i = 1usize;
        for stripe in &self.base.my_list_stripe {
            if i == ic {
                let st = stripe.read().expect("stripe lock");
                // OCCT: bid = itel.Value()->SetOfSurfData()->Value(IS)->Simul();
                //       res = down_cast<HArray1<ChFiDS_CircSection>>(bid).
                // Architecture: the rcad simul slot carries the concrete
                // section array, so the transient-handle down-cast is the
                // slot read itself.
                let sd = st.my_hdata.get(is.wrapping_sub(1))?;
                let guard = sd.read().expect("surfdata lock");
                return guard.simul();
            }
            i += 1;
        }
        None
    }

    /// OCCT ChFi3d_FilBuilder.cxx L473-560 — SimulKPart(SD): build the
    /// circular sections simulating the KPart fillet surface and store
    /// them through SD->SetSimul.
    pub fn simul_kpart(&self, sd: &mut ChFiDSSurfData) {
        // L474-475: TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
        //           occ::handle<Geom_Surface> S = DStr.Surface(SD->Surf()).Surface();
        let dstr = self.base.my_ds.as_ref().expect("DS");
        let s: &Surface3 = &dstr.surface(sd.surf()).surface;
        // L476-481: the four pcurve evaluations at the interference bounds.
        let fi1 = sd.interference_on_s1();
        let fi2 = sd.interference_on_s2();
        let p1f = fi1
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(fi1.parameter_first());
        let p1l = fi1
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(fi1.parameter_last());
        let p2f = fi2
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(fi2.parameter_first());
        let p2l = fi2
            .pcurve_on_surf()
            .expect("PCurveOnSurf")
            .point_at(fi2.parameter_last());
        // L482-484: GeomAdaptor_Surface AS(S); typ = AS.GetType();
        let typ = surface_type_of(s);
        let mut sec: Option<ChFiDSCircSectionArray> = None;
        match typ {
            GeomAbsSurfaceType::Cylinder => {
                // L489-492
                let u1 = p1f.x;
                let u2 = p2f.x;
                let v1 = p1f.y.max(p2f.y);
                let v2 = p1l.y.min(p2l.y);
                // L493: sec = new NCollection_HArray1<ChFiDS_CircSection>(1, 2);
                let mut arr = vec![ChFiDSCircSection::new(), ChFiDSCircSection::new()];
                // L494: gp_Cylinder Cy = AS.Cylinder();
                let pos = surface3_ax3(s);
                let radius = match s {
                    Surface3::Cylinder(c) => c.radius,
                    _ => unreachable!("GeomAbs_Cylinder branch"),
                };
                // L497-498: sec1/sec2.Set(ElSLib::CylinderVIso(Cy.Position(),
                // Cy.Radius(), v1/v2), u1, u2).  CylinderVIso (ElSLib.cxx
                // L1781-1789): the Ax2 of the cylinder frame translated by
                // v * Direction.
                let cy_iso = |v: f64| {
                    let axes = GpAx3 {
                        location: pos.location + pos.vzdir * v,
                        vxdir: pos.vxdir,
                        vydir: pos.vydir,
                        vzdir: pos.vzdir,
                    };
                    GpCirc::new(axes, radius).to_circle3()
                };
                arr[0].set_circ(cy_iso(v1), u1, u2);
                arr[1].set_circ(cy_iso(v2), u1, u2);
                sec = Some(arr);
            }
            GeomAbsSurfaceType::Torus => {
                // L505-509
                let v1 = p1f.y;
                let v2 = p2f.y;
                let u1 = p1f.x.max(p2f.x);
                let u2 = p1l.x.min(p2l.x);
                // L510-513
                let ang = u2 - u1;
                let (majr, minr) = match s {
                    Surface3::Torus(t) => (t.major_radius, t.minor_radius),
                    _ => unreachable!("GeomAbs_Torus branch"),
                };
                let mut n = (36.0 * ang / std::f64::consts::PI + 1.0) as i32;
                if n < 2 {
                    n = 2;
                }
                // L514: sec = new NCollection_HArray1<ChFiDS_CircSection>(1, n);
                let pos = surface3_ax3(s);
                let mut arr: ChFiDSCircSectionArray = Vec::new();
                for i in 1..=n {
                    // L516-518: isec.Set(ElSLib::TorusUIso(To.Position(),
                    // majr, minr, u), v1, v2).  TorusUIso (ElSLib.cxx
                    // L1751-1770): cx = cos(U)*dx + sin(U)*dy;
                    // axes = Ax2(Location, cx ^ dz, cx) translated by
                    // cx * MajorRadius; Circ(axes, MinorRadius).
                    let u = u1 + (i as f64 - 1.0) * (u2 - u1) / (n as f64 - 1.0);
                    let cx = pos.vxdir * u.cos() + pos.vydir * u.sin();
                    let n_dir = cx.cross(pos.vzdir);
                    let mut axes = GpAx3 {
                        location: pos.location,
                        vxdir: cx,
                        vydir: n_dir.cross(cx),
                        vzdir: n_dir,
                    };
                    axes.location += cx * majr;
                    let isec = arr.last_mut().expect("isec");
                    isec.set_circ(GpCirc::new(axes, minr).to_circle3(), v1, v2);
                }
                sec = Some(arr);
            }
            GeomAbsSurfaceType::Sphere => {
                // L527-531
                let v1 = p1f.y;
                let v2 = p2f.y;
                let u1 = p1f.x.max(p2f.x);
                let u2 = p1l.x.min(p2l.x);
                // L532-535
                let ang = u2 - u1;
                let rad = match s {
                    Surface3::Sphere(sp) => sp.radius,
                    _ => unreachable!("GeomAbs_Sphere branch"),
                };
                let mut n = (36.0 * ang / std::f64::consts::PI + 1.0) as i32;
                if n < 2 {
                    n = 2;
                }
                // L536: sec = new NCollection_HArray1<ChFiDS_CircSection>(1, n);
                let pos = surface3_ax3(s);
                let mut arr: ChFiDSCircSectionArray = Vec::new();
                for i in 1..=n {
                    // L538-540: isec.Set(ElSLib::SphereUIso(Sp.Position(),
                    // rad, u), v1, v2).  SphereUIso (ElSLib.cxx L1738-1748):
                    // cx = cos(U)*dx + sin(U)*dy;
                    // Circ(Ax2(Location, cx ^ dz, cx), Radius).
                    let u = u1 + (i as f64 - 1.0) * (u2 - u1) / (n as f64 - 1.0);
                    let cx = pos.vxdir * u.cos() + pos.vydir * u.sin();
                    let n_dir = cx.cross(pos.vzdir);
                    let axes = GpAx3 {
                        location: pos.location,
                        vxdir: cx,
                        vydir: n_dir.cross(cx),
                        vzdir: n_dir,
                    };
                    let isec = arr.last_mut().expect("isec");
                    isec.set_circ(GpCirc::new(axes, rad).to_circle3(), v1, v2);
                }
                sec = Some(arr);
            }
            _ => {
                // L544-546: default — sec stays null.
            }
        }
        // L559: SD->SetSimul(sec);
        sd.set_simul(sec);
    }

    /// OCCT ChFi3d_FilBuilder.cxx L2595-2637 — SetRegul: code the
    /// regularities after cutting (continuity of each regular edge against
    /// its two faces).
    pub fn set_regul(&mut self) {
        let coup = self.base.my_coup.as_ref().expect("HBuilder");
        for reg in &self.base.my_regul {
            // L2605: itc.Initialize(myCoup->NewEdges(reg.Curve()));
            let itc = coup.new_edges(reg.curve());
            if let Some(e) = itc.first() {
                // L2608: TopoDS_Edge E = TopoDS::Edge(itc.Value());
                let e = e.clone();
                let dstr = self.base.my_ds.as_ref().expect("DS");
                // L2609-2615
                let its1: Vec<Shape> = if reg.is_surface1() {
                    // L2611: its1.Initialize(myCoup->NewFaces(reg.S1()));
                    coup.new_faces(reg.s1())
                } else {
                    // L2613: its1.Initialize(myCoup->Merged(
                    //           myDS->Shape(reg.S1()), TopAbs_IN));
                    let s1shape = dstr.shape(reg.s1()).clone();
                    coup.merged(&s1shape, TopAbsState::In)
                };
                // L2616-2621
                let its2: Vec<Shape> = if reg.is_surface2() {
                    // L2618: its2.Initialize(myCoup->NewFaces(reg.S2()));
                    coup.new_faces(reg.s2())
                } else {
                    // L2620: its2.Initialize(myCoup->Merged(
                    //           myDS->Shape(reg.S2()), TopAbs_IN));
                    let s2shape = dstr.shape(reg.s2()).clone();
                    coup.merged(&s2shape, TopAbsState::In)
                };
                // L2622: if (its1.More() && its2.More())
                if let (Some(f1s), Some(f2s)) = (its1.first(), its2.first()) {
                    // L2623-2624: F1 = Face(its1.Value()); F2 = Face(its2.Value());
                    let f1 = f1s.clone();
                    let f2 = f2s.clone();
                    // L2625-2626:
                    //   GeomAbs_Shape cont = ChFi3d_evalconti(E, F1, F2);
                    //   B.Continuity(E, F1, F2, cont);
                    // Pending boundary: ChFi3d_evalconti is private to
                    // chfi3d_builder_chbuilder.rs and BRep_Builder::Continuity
                    // has no rcad carrier yet (the same neutralization as
                    // chfi3d_builder_chbuilder.rs L1694-1697).  The facade
                    // feeds stay empty until W11, so this branch is
                    // unreachable either way.
                    let _ = (&e, &f1, &f2);
                }
            }
        }
    }
}
