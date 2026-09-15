//! OCCT BRepLib::UpdateTolerances (ModelingAlgorithms/TKTopAlgo/BRepLib/
//! BRepLib.cxx) — the 1:1 translation of the tolerance harmonization:
//! - static `UpdShTol` (cxx L837-909) — the tolerance-map write-back
//!   (keep-max semantics, the mutable-input / reshaper write paths),
//! - static `InternalUpdateTolerances` (cxx L1744-1962) — the
//!   harmonization with rule Tolerance(VERTEX)>=Tolerance(EDGE)>=
//!   Tolerance(FACE), the IsVerifyTolerance face floor by surface kind,
//!   the EDGE walk (MapShapesAndAncestors EDGE/FACE) and the VERTEX walk
//!   (MapShapesAndUniqueAncestors VERTEX/EDGE),
//! - `BRepLib::UpdateTolerances(S, verifyFaceTolerance)` (cxx L1966-1970),
//! - `BRepLib::UpdateTolerances(S, theReshaper, verifyFaceTolerance)`
//!   (cxx L1974-1979).
//!
//! Architecture difference (the rcad pool model, topods::BRep): the OCCT
//! bodies mutate the shared TShapes in place through the handles; rcad
//! mirrors this through the pool (`BRep::edge_mut_inplace` / the in-place
//! raw slot write) — the functions take `&mut BRep` where OCCT mutates, and
//! the free-function form follows the `same_parameter.rs` convention.
//!
//! Encodings recorded along the way:
//! - `NCollection_DataMap<TopoDS_Shape, double>` -> [`ShToTolMap`]
//!   (insertion-ordered; the OCCT bucket order is irrelevant here because
//!   every entry addresses a distinct shape),
//! - `BRepBndLib::Add(curf, aB)` -> the reduced sampling re-host
//!   [`brep_bnd_lib_add_face`] (the same encoding as
//!   `brep_offset_inter3d.rs`),
//! - `TopExp::MapShapesAndAncestors` ->
//!   `crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors`,
//!   `TopExp::MapShapesAndUniqueAncestors` ->
//!   `crate::offset::brep_offset_make_offset::top_exp_map_shapes_and_unique_ancestors`,
//! - `BRepTools_ReShape` -> `crate::shhealing::shape_build::ShapeBuildReShape`,
//! - `Epsilon(tol)` -> `rcad_kernel::base::extrema_ext_elc::epsilon_of`
//!   (the kernel keeps the single canonical `std::nextafter` definition).

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use indexmap::IndexMap;
use rcad_kernel::base::extrema_ext_elc::epsilon_of;
use rcad_kernel::base::proj_lib::{GeomAbsSurfaceType, GeomSurfaceAdaptor};
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve2dEval, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo::topods::{
    edge_data_pool_free, shape_is_in_pool, tshape_flags, BRep, BRepTool, CurveRepresentation,
    Orientation, ShapeType, TEdgeData, TShape,
};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_parameter, brep_tool_pnt, brep_tool_surface, brep_tool_tolerance,
    empty_copied, explorer, sub_shapes,
};
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::offset::brep_offset_make_offset::top_exp_map_shapes_and_unique_ancestors;
use crate::shhealing::shape_build::ShapeBuildReShape;

/// OCCT `NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher>`
/// (cxx L1749): the shape -> new-tolerance map.  The key is the
/// TopTools_ShapeMapHasher identity (TShape + Location, orientation
/// ignored); the value carries the shape itself so the UpdShTol iterator
/// can address the OCCT map key.
type ShToTolMap = IndexMap<(u64, u32), (Shape, f64)>;

/// OCCT
/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// (cxx L1827) — the ancestors map (the same keying as ShToTolMap).
type ParentsMap = IndexMap<(u64, u32), (Shape, Vec<Shape>)>;

/// The shape-identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored) — the same encoding as the map helpers.
fn shape_key(the_s: &Shape) -> (u64, u32) {
    (the_s.ptr_id(), the_s.location)
}

// =========================================================================
// static UpdShTol (BRepLib.cxx L837-909).
// =========================================================================

/// OCCT static UpdShTol(theShToTol, IsMutableInput, theReshaper,
/// theVForceUpdate) — Update vertices/edges/faces according to ShToTol map
/// (create copies of necessary).
fn upd_sh_tol(
    the_brep: &mut BRep,
    the_sh_to_tol: &ShToTolMap,
    is_mutable_input: bool,
    the_reshaper: &mut ShapeBuildReShape,
    the_v_force_update: bool,
) {
    // OCCT L843: BRep_Builder aB — the pool write helpers below.
    // OCCT L844-846: the DataMap iterator.
    for (_a_key, (a_sh, a_tol)) in the_sh_to_tol.iter() {
        let a_tol = *a_tol;
        //
        // OCCT L851: TopoDS_Shape aNsh.
        let mut a_nsh;
        // OCCT L852: const TopoDS_Shape& aVsh = theReshaper.Value(aSh).
        let a_vsh = the_reshaper.value(the_brep, a_sh);
        // (The addendum-48 W2 engine-side normalization — "Value() answered a
        // NULL shape for a pool-free input" — is retired here: the reshape
        // IsNull gates are pool-free-aware since the addendum-49 reshape fix
        // (the null-sentinel signature match), so Value() returns the shape
        // itself for an unrecorded pool-free shape, exactly the OCCT
        // L243-253 contract.)
        // OCCT L853.
        let use_old_sh =
            is_mutable_input || the_reshaper.is_new_shape(a_sh) || !a_vsh.is_same(a_sh);
        if use_old_sh {
            // OCCT L855-857: aNsh = aVsh.
            a_nsh = a_vsh;
        } else {
            // OCCT L860: aNsh = aSh.EmptyCopied().
            a_nsh = if shape_is_in_pool(the_brep, a_sh) {
                the_brep.empty_copied(a_sh)
            } else {
                // The pool-free source (index == usize::MAX): the data-based
                // EmptyCopied re-host.  BRep::empty_copy indexes the pool
                // slot `tshapes[r.index]` and would panic on the pool-free
                // index; the copy is pool-free, the same encoding as
                // brep_algo::tool::empty_copied (OCCT has no pool — the
                // empty copy of the TShape is the TShape).
                empty_copied(a_sh)
            };
            // OCCT L861-866: add subshapes from the original shape
            // (TopoDS_Iterator sit(aSh); aB.Add(aNsh, sit.Value())).
            for a_sub in sub_shapes(a_sh) {
                builder_add_shape(the_brep, &mut a_nsh, &a_sub);
            }
            //
            // OCCT L868-873: aNsh.Free(aSh.Free()) ...
            // aNsh.Convex(aSh.Convex()).
            set_shape_flags(the_brep, &mut a_nsh, shape_flags(a_sh));
            //
        }
        //
        // OCCT L877: switch (aSh.ShapeType()).
        match a_nsh.shape_type() {
            // OCCT L879-881: aB.UpdateFace(TopoDS::Face(aNsh), aTol).
            ShapeType::Face => {
                builder_update_face_tol(the_brep, &mut a_nsh, a_tol);
            }
            // OCCT L883-885: aB.UpdateEdge(TopoDS::Edge(aNsh), aTol).
            ShapeType::Edge => {
                builder_update_edge_tol(the_brep, &mut a_nsh, a_tol);
            }
            // OCCT L887-898: the direct BRep_TVertex write path.
            ShapeType::Vertex => {
                tvertex_set_tolerance(the_brep, &mut a_nsh, a_tol, the_v_force_update);
            }
            // OCCT L900-901: default: break.
            _ => {}
        }
        //
        // OCCT L904-907.
        if !use_old_sh {
            the_reshaper.replace(the_brep, a_sh, &a_nsh);
        }
    }
}

// =========================================================================
// static InternalUpdateTolerances (BRepLib.cxx L1744-1962).
// =========================================================================

/// OCCT static InternalUpdateTolerances(theOldShape, IsVerifyTolerance,
/// IsMutableInput, theReshaper) — the tolerance harmonization (also the
/// tail of InternalSameParameter, cxx L1007).
pub fn internal_update_tolerances(
    the_brep: &mut BRep,
    the_old_shape: &Shape,
    is_verify_tolerance: bool,
    is_mutable_input: bool,
    the_reshaper: &mut ShapeBuildReShape,
) {
    // OCCT L1749: NCollection_DataMap<TopoDS_Shape, double> aShToTol.
    let mut a_sh_to_tol = ShToTolMap::new();
    // Harmonize tolerances
    // with rule Tolerance(VERTEX)>=Tolerance(EDGE)>=Tolerance(FACE)
    // OCCT L1752.
    let mut tol = 0.0f64;
    if is_verify_tolerance {
        // Set tolerance to its minimum value
        // OCCT L1756-1760: S; l; ex; aB; the extent locals.
        let mut a_b = BndBox::new();
        // OCCT L1761: for (ex.Init(theOldShape, TopAbs_FACE); ...).
        for a_curf in explorer(the_old_shape, ShapeType::Face, ShapeType::Shape) {
            // OCCT L1764: S = BRep_Tool::Surface(curf, l) — the location out
            // parameter l is not used by the OCCT body after the read.
            let a_s_val = brep_tool_surface(&a_curf);
            // OCCT L1765: if (!S.IsNull()).
            if let Some(a_s_val) = a_s_val {
                // OCCT L1767: aB.SetVoid().
                a_b.set_void();
                // OCCT L1768: BRepBndLib::Add(curf, aB).
                brep_bnd_lib_add_face(&a_curf, &mut a_b);
                // OCCT L1769-1772: the RectangularTrimmedSurface basis strip.
                let a_s = match &a_s_val {
                    Surface3::Trimmed(a_t) => a_t.basis.as_ref().clone(),
                    _ => a_s_val,
                };
                // OCCT L1773: GeomAdaptor_Surface AS(S).
                let a_as = GeomSurfaceAdaptor::new(a_s);
                // OCCT L1774-1789.
                match a_as.my_surface_type {
                    GeomAbsSurfaceType::Plane
                    | GeomAbsSurfaceType::Cylinder
                    | GeomAbsSurfaceType::Cone => {
                        tol = CONFUSION;
                    }
                    GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Torus => {
                        tol = CONFUSION * 2.0;
                    }
                    _ => {
                        tol = CONFUSION * 4.0;
                    }
                }
                // OCCT L1790-1820.
                if !a_b.is_whole() {
                    // OCCT L1792: aB.Get(...) — the rcad void box reads None
                    // (the Standard_ConstructionError state; unreachable here
                    // because a non-whole box carries finite corners).
                    if let Some((a_xmin, mut a_ymin, mut a_zmin, a_xmax, a_ymax, a_zmax)) =
                        a_b.get()
                    {
                        // OCCT L1793: dMax = 1.
                        let mut a_dmax = 1.0f64;
                        // OCCT L1794-1812 (the aYmin/aZmin reuse is the OCCT
                        // quirk — the interval extent lands back in the min
                        // local).
                        if !a_b.is_open_xmin() && !a_b.is_open_xmax() {
                            a_dmax = a_xmax - a_xmin;
                        }
                        if !a_b.is_open_ymin() && !a_b.is_open_ymax() {
                            a_ymin = a_ymax - a_ymin;
                        }
                        if !a_b.is_open_zmin() && !a_b.is_open_zmax() {
                            a_zmin = a_zmax - a_zmin;
                        }
                        if a_ymin > a_dmax {
                            a_dmax = a_ymin;
                        }
                        if a_zmin > a_dmax {
                            a_dmax = a_zmin;
                        }
                        // OCCT L1814.
                        tol = tol * a_dmax;
                        // Do not process tolerances > 1.
                        // OCCT L1816-1819.
                        if tol > 1.0 {
                            tol = 0.99;
                        }
                    }
                }
                // OCCT L1821: aShToTol.Bind(curf, tol).
                a_sh_to_tol.insert(shape_key(&a_curf), (a_curf.clone(), tol));
            }
        }
    }

    // Process edges
    // OCCT L1827-1829.
    let mut parents: ParentsMap = ParentsMap::new();
    map_shapes_and_ancestors(the_old_shape, ShapeType::Edge, ShapeType::Face, &mut parents);
    // OCCT L1832: for (iCur = 1; iCur <= parents.Extent(); iCur++).
    for (_a_key, (a_ek, a_list)) in parents.iter() {
        // OCCT L1834.
        tol = 0.0;
        // OCCT L1835: for (lConx.Initialize(parents(iCur)); ...).
        for a_ff in a_list {
            // OCCT L1839-1846.
            let a_ftol = if is_verify_tolerance && a_sh_to_tol.contains_key(&shape_key(a_ff)) {
                // first condition for speed-up
                a_sh_to_tol.get(&shape_key(a_ff)).expect("bound above").1
            } else {
                // tolerance have not been updated
                brep_tool_tolerance(a_ff)
            };
            // OCCT L1847.
            tol = tol.max(a_ftol);
        }
        // Update can only increase tolerance, so if the edge has a greater
        //  tolerance than its faces it is not concerned
        // OCCT L1851-1855.
        if tol > brep_tool_tolerance(a_ek) {
            a_sh_to_tol.insert(shape_key(a_ek), (a_ek.clone(), tol));
        }
    }

    // Vertices are processed
    // OCCT L1858-1860.
    let a_big_tol = 1.0e10f64;
    let mut parents: ParentsMap = ParentsMap::new();

    // The pool-free owning-face index for gcurve_surface (architecture
    // difference: OCCT's BRep_CurveOnSurface carries the surface handle in
    // the representation; the rcad representation carries the owning face
    // KEY, resolved against the pool.  A pool-free owning face — index ==
    // usize::MAX, its TShape Arc absent from the pool — is resolved from
    // this walk over theOldShape, the same bridge as edge_data_pool_free).
    let mut a_graph_faces: HashMap<u64, Surface3> = HashMap::new();
    for a_f in explorer(the_old_shape, ShapeType::Face, ShapeType::Shape) {
        if let Some(a_s) = brep_tool_surface(&a_f) {
            a_graph_faces.entry(a_f.ptr_id()).or_insert(a_s);
        }
    }

    // OCCT L1862.
    top_exp_map_shapes_and_unique_ancestors(
        the_old_shape,
        ShapeType::Vertex,
        ShapeType::Edge,
        &mut parents,
    );
    // OCCT L1863: NCollection_Map<handle(Standard_Transient)> Initialized —
    // keyed by the vertex TShape identity (the rcad pool pointer).
    let mut a_initialized: HashSet<u64> = HashSet::new();
    // OCCT L1864-1865: for (iCur = 1; iCur <= nbV; iCur++).
    for (_a_key, (a_v, a_list)) in parents.iter() {
        // OCCT L1867.
        tol = 0.0;
        // OCCT L1868: gp_Pnt aPV = BRep_Tool::Pnt(V) — the pool-free
        // BRep_Tool::Pnt re-host (the walk may reach pool-free vertices in
        // a mixed adopted graph; the pool-only vertex_position would
        // index out of bounds).
        let a_pv = brep_tool_pnt(a_v).unwrap_or(glam::DVec3::ZERO);
        // OCCT L1870-1871.
        let mut a_max_dist = 0.0f64;
        let mut a_p3d;
        // OCCT L1872: for (lConx.Initialize(parents(iCur)); ...).
        for a_e in a_list {
            // OCCT L1875-1876.
            let a_ntol = a_sh_to_tol.get(&shape_key(a_e)).map(|a_entry| a_entry.1);
            tol = tol.max(a_ntol.unwrap_or_else(|| brep_tool_tolerance(a_e)));
            // OCCT L1877-1880.
            if tol > a_big_tol {
                continue;
            }
            // OCCT L1881-1884.
            if !edge_data(the_brep, a_e).same_range {
                continue;
            }
            // OCCT L1885: double par = BRep_Tool::Parameter(V, E).
            let a_par = brep_tool_parameter(a_v, a_e);
            // OCCT L1886-1888: TE = E.TShape(); the Curves() list — the rcad
            // representation list.
            let a_te = edge_data(the_brep, a_e);
            for a_cr in a_te.representations.iter() {
                // For each CurveRepresentation, check the provided parameter
                // OCCT L1893-1894: L = Eloc * loc — the rcad world-space
                // encoding carries the locations in the stored geometry (the
                // same recorded no-op as same_parameter.rs L1333).
                // OCCT L1895-1908.
                if matches!(a_cr, CurveRepresentation::Curve3D { .. }) {
                    // OCCT L1897-1899: C = cr->Curve3D(); !C.IsNull() — the
                    // rcad Curve3D representation reads the edge 3d curve
                    // slot (the same arch. difference as build_curves3d.rs).
                    let a_ed = edge_data(the_brep, a_e);
                    if let Some(a_c) = &a_ed.curve {
                        // edge non degenerated
                        a_p3d = CurveEval::point_at(a_c, a_par);
                        let a_dist = a_p3d.distance_squared(a_pv);
                        if a_dist > a_max_dist {
                            a_max_dist = a_dist;
                        }
                    }
                }
                // OCCT L1909-1937.
                else if gcurve_is_curve_on_surface(a_cr) {
                    // OCCT L1911-1917: Su = cr->Surface(); PC = cr->PCurve();
                    // PC2 = IsCurveOnClosedSurface() ? PCurve2() : null — the
                    // surface is resolved from the owning face key (the
                    // gcurve_surface encoding of same_parameter.rs).
                    let a_su = gcurve_surface(the_brep, a_cr, &a_graph_faces);
                    let a_pc = gcurve_pcurve(a_cr);
                    let a_pc2 = if gcurve_is_curve_on_closed_surface(a_cr) {
                        gcurve_pcurve2(a_cr)
                    } else {
                        None
                    };
                    if let (Some(a_su), Some(a_pc)) = (a_su, a_pc) {
                        // OCCT L1918-1925.
                        let a_p2d = Curve2dEval::point_at(&a_pc, a_par);
                        a_p3d = SurfaceEval::point_at(&a_su, a_p2d.x, a_p2d.y);
                        let a_dist = a_p3d.distance_squared(a_pv);
                        if a_dist > a_max_dist {
                            a_max_dist = a_dist;
                        }
                        // OCCT L1926-1936.
                        if let Some(a_pc2) = a_pc2 {
                            let a_p2d = Curve2dEval::point_at(&a_pc2, a_par);
                            a_p3d = SurfaceEval::point_at(&a_su, a_p2d.x, a_p2d.y);
                            let a_dist = a_p3d.distance_squared(a_pv);
                            if a_dist > a_max_dist {
                                a_max_dist = a_dist;
                            }
                        }
                    }
                }
            }
        }
        // OCCT L1941-1942.
        tol = tol.max(a_max_dist.sqrt());
        tol += 2.0 * epsilon_of(tol);
        //
        // OCCT L1944-1947.
        let a_vtol = brep_tool_tolerance(a_v);
        let an_upd_tol = tol > a_vtol;
        let mut a_to_add = false;
        if is_verify_tolerance {
            // ASet minimum value of the tolerance
            // Attention to sharing of the vertex by other shapes
            // OCCT L1952.
            a_to_add = a_initialized.insert(a_v.ptr_id()) && a_vtol != tol;
        }
        //'Initialized' map is not used anywhere outside this block
        // OCCT L1955-1958.
        if an_upd_tol || a_to_add {
            a_sh_to_tol.insert(shape_key(a_v), (a_v.clone(), tol));
        }
    }

    // OCCT L1961: UpdShTol(aShToTol, IsMutableInput, theReshaper, true).
    upd_sh_tol(the_brep, &a_sh_to_tol, is_mutable_input, the_reshaper, true);
}

// =========================================================================
// BRepLib::UpdateTolerances — the public overloads (BRepLib.cxx
// L1966-1979).
// =========================================================================

/// OCCT BRepLib::UpdateTolerances(S, verifyFaceTolerance) (cxx L1966-1970).
pub fn update_tolerances(the_brep: &mut BRep, the_s: &Shape, verify_face_tolerance: bool) {
    // The pcurves-writer-migration sweep: materialize every map-only
    // pcurve entry as the OCCT-faithful CurveOnSurface representation
    // before the walk, so the engine's representation-based reads engage
    // on map-era edges too (the topo_builder::migrate_map_to_representations
    // contract; a no-op on pool-free graphs and on already-dual edges).
    rcad_kernel::topo::topo_builder::migrate_map_to_representations(the_brep, the_s);
    // OCCT L1968: BRepTools_ReShape aReshaper.
    let mut a_reshaper = ShapeBuildReShape::new();
    // OCCT L1969.
    internal_update_tolerances(the_brep, the_s, verify_face_tolerance, true, &mut a_reshaper);
}

/// OCCT BRepLib::UpdateTolerances(S, theReshaper, verifyFaceTolerance)
/// (cxx L1974-1979).
pub fn update_tolerances_with_reshaper(
    the_brep: &mut BRep,
    the_s: &Shape,
    the_reshaper: &mut ShapeBuildReShape,
    verify_face_tolerance: bool,
) {
    // OCCT L1978.
    internal_update_tolerances(the_brep, the_s, verify_face_tolerance, false, the_reshaper);
}

// =========================================================================
// File-local re-hosts (the OCCT builder / representation writes; private in
// the sibling modules, so this file carries its own copies — the same OCCT
// bodies).
// =========================================================================

/// OCCT BRep_Builder::UpdateFace(F, Tol) (BRep_Builder.cxx L597-607):
/// TF->Tolerance(Tol) — the plain set — then F.TShape()->Modified(true).
fn builder_update_face_tol(the_brep: &mut BRep, the_f: &mut Shape, the_tol: f64) {
    if shape_is_in_pool(the_brep, the_f) {
        let a_tf = the_brep.face_mut(the_f.clone());
        a_tf.tolerance = the_tol;
        a_tf.flags |= tshape_flags::MODIFIED;
    } else {
        // The pool-free arm: the SAME in-place TShape write as the pool
        // arm, through the shape's own Arc.  SAFETY: single-threaded; the
        // caller holds &mut BRep and no &TShape borrow of this Arc is
        // alive; every referencing shape observes the change, matching the
        // OCCT TShape handle mutation.  (An Arc::make_mut write would
        // clone-on-write and detach the update from the adopted mixed
        // graph — the harmonization would not engage.)
        let a_ptr = Arc::as_ptr(&the_f.data) as *mut TShape;
        let a_tf = unsafe {
            match &mut *a_ptr {
                TShape::Face(a_tf) => a_tf,
                _ => return,
            }
        };
        a_tf.tolerance = the_tol;
        a_tf.flags |= tshape_flags::MODIFIED;
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, Tol) (BRep_Builder.cxx L999-1009):
/// TE->UpdateTolerance(Tol) — the keep-max — then E.TShape()->Modified(true).
fn builder_update_edge_tol(the_brep: &mut BRep, the_e: &mut Shape, the_tol: f64) {
    if shape_is_in_pool(the_brep, the_e) {
        let a_te = the_brep.edge_mut_inplace(the_e.clone());
        a_te.tolerance = a_te.tolerance.max(the_tol);
        a_te.flags |= tshape_flags::MODIFIED;
    } else {
        // The pool-free arm: the SAME in-place TShape write as the pool
        // arm (see builder_update_face_tol for the SAFETY / engagement
        // rationale — the Arc::make_mut encoding would detach the update).
        let a_ptr = Arc::as_ptr(&the_e.data) as *mut TShape;
        let a_te = unsafe {
            match &mut *a_ptr {
                TShape::Edge(a_te) => a_te,
                _ => return,
            }
        };
        a_te.tolerance = a_te.tolerance.max(the_tol);
        a_te.flags |= tshape_flags::MODIFIED;
    }
}

/// OCCT L888-897: the `handle(BRep_TVertex)& aTV` write path — force
/// `aTV->Tolerance(aTol)` (L891) or the keep-max
/// `aTV->UpdateTolerance(aTol)` (L895), then `aTV->Modified(true)` (L897).
fn tvertex_set_tolerance(the_brep: &mut BRep, the_v: &mut Shape, the_tol: f64, the_force_set: bool) {
    if shape_is_in_pool(the_brep, the_v) {
        // SAFETY: single-threaded; the caller holds &mut BRep (the exclusive
        // borrow of the pool slot) and no other &TShape for this Arc is
        // alive; every referencing shape observes the change, matching the
        // OCCT TShape handle mutation (the same in-place encoding as
        // BRep::edge_mut_inplace — Arc::make_mut would clone-on-write and
        // split the identity).
        let a_ptr = Arc::as_ptr(&the_brep.tshapes[the_v.index]) as *mut TShape;
        let a_tv = unsafe {
            match &mut *a_ptr {
                TShape::Vertex(a_tv) => a_tv,
                _ => return,
            }
        };
        if the_force_set {
            // OCCT L891: aTV->Tolerance(aTol).
            a_tv.tolerance = the_tol;
        } else {
            // OCCT L895: aTV->UpdateTolerance(aTol) — keep-max.
            a_tv.tolerance = a_tv.tolerance.max(the_tol);
        }
        // OCCT L897: aTV->Modified(true).
        a_tv.flags |= tshape_flags::MODIFIED;
    } else {
        // The pool-free arm: the SAME in-place TShape write as the pool
        // arm, through the shape's own Arc (see builder_update_face_tol
        // for the SAFETY / engagement rationale — the Arc::make_mut
        // encoding would detach the update from the shared graph).
        let a_ptr = Arc::as_ptr(&the_v.data) as *mut TShape;
        let a_tv = unsafe {
            match &mut *a_ptr {
                TShape::Vertex(a_tv) => a_tv,
                _ => return,
            }
        };
        if the_force_set {
            // OCCT L891: aTV->Tolerance(aTol).
            a_tv.tolerance = the_tol;
        } else {
            // OCCT L895: aTV->UpdateTolerance(aTol) — keep-max.
            a_tv.tolerance = a_tv.tolerance.max(the_tol);
        }
        // OCCT L897: aTV->Modified(true).
        a_tv.flags |= tshape_flags::MODIFIED;
    }
}

/// The TopoDS_TShape flag word of a shape (the OCCT Free/Checked/Orientable/
/// Closed/Infinite/Convex reads of cxx L868-873).
fn shape_flags(the_s: &Shape) -> u16 {
    match the_s.data.as_ref() {
        TShape::Vertex(a_v) => a_v.flags,
        TShape::Edge(a_e) => a_e.flags,
        TShape::Wire(a_w) => a_w.flags,
        TShape::Face(a_f) => a_f.flags,
        TShape::Shell(a_s) => a_s.flags,
        TShape::Solid(a_s) => a_s.flags,
        _ => 0,
    }
}

/// The OCCT Free/Checked/Orientable/Closed/Infinite/Convex writes onto the
/// (empty-copied) shape (cxx L868-873).
fn set_shape_flags(the_brep: &mut BRep, the_s: &mut Shape, the_flags: u16) {
    let a_apply = |a_ts: &mut TShape| match a_ts {
        TShape::Vertex(a_v) => a_v.flags = the_flags,
        TShape::Edge(a_e) => a_e.flags = the_flags,
        TShape::Wire(a_w) => a_w.flags = the_flags,
        TShape::Face(a_f) => a_f.flags = the_flags,
        TShape::Shell(a_s) => a_s.flags = the_flags,
        TShape::Solid(a_s) => a_s.flags = the_flags,
        _ => {}
    };
    if shape_is_in_pool(the_brep, the_s) {
        // SAFETY: see tvertex_set_tolerance.
        let a_ptr = Arc::as_ptr(&the_brep.tshapes[the_s.index]) as *mut TShape;
        let a_ts = unsafe { &mut *a_ptr };
        a_apply(a_ts);
    } else {
        a_apply(Arc::make_mut(&mut the_s.data));
    }
}

/// OCCT L862-866: aB.Add(aNsh, sit.Value()) — the TopoDS_Builder::Add ->
/// virtual TShape::Add dispatch over the tolerance-map shape kinds
/// (Edge<-Vertex, Wire<-Edge, Face<-Wire).
fn builder_add_shape(the_brep: &mut BRep, the_nsh: &mut Shape, the_sub: &Shape) {
    if !shape_is_in_pool(the_brep, the_nsh) {
        // The pool-free arm: the same OCCT TShape::Add writes through the
        // shape's own Arc (the established pool-free setter pattern, the
        // encoding of brep_algo::tool::builder_add_edge_vertex et al.).
        // The raw-slot accessors below (edge_mut_inplace / wire_mut /
        // face_mut) index `tshapes[usize::MAX]` and would panic.
        match Arc::make_mut(&mut the_nsh.data) {
            TShape::Edge(a_ed) => match the_sub.orientation {
                Orientation::Reversed => a_ed.last = the_sub.clone(),
                _ => a_ed.first = the_sub.clone(),
            },
            TShape::Wire(a_wd) => {
                a_wd.edges.push(the_sub.clone());
                a_wd.my_shapes.push(the_sub.clone());
            }
            TShape::Face(a_fd) => {
                a_fd.inner_wires.push(the_sub.clone());
                a_fd.my_shapes.push(the_sub.clone());
            }
            _ => {}
        }
        return;
    }
    match the_nsh.data.as_ref() {
        TShape::Edge(_) => {
            // BRep_TEdge::Add(V) — the extremity assignment by orientation
            // (the tool::builder_add_edge_vertex encoding).
            let a_ed = the_brep.edge_mut_inplace(the_nsh.clone());
            match the_sub.orientation {
                Orientation::Reversed => a_ed.last = the_sub.clone(),
                _ => a_ed.first = the_sub.clone(),
            }
        }
        TShape::Wire(_) => {
            // BRep_TWire::Add(E) — the edge append (the pool
            // BRepBuilder::add_to_wire encoding).
            let a_wd = the_brep.wire_mut(the_nsh.clone());
            a_wd.edges.push(the_sub.clone());
            a_wd.my_shapes.push(the_sub.clone());
        }
        TShape::Face(_) => {
            // BRep_TFace::Add(W) — the inner-wire append (the pool
            // BRepBuilder::add_to_face encoding).
            let a_fd = the_brep.face_mut(the_nsh.clone());
            a_fd.inner_wires.push(the_sub.clone());
            a_fd.my_shapes.push(the_sub.clone());
        }
        _ => {}
    }
}

/// OCCT BRepBndLib::Add(curf, aB) (BRepLib.cxx L1768) — the reduced sampling
/// re-host (architecture difference #36): the surface part goes through the
/// kernel `surface_bounding_box` (the established BndLib_AddSurface::Add
/// re-host, the same encoding as bop/ds prepare_faces), then the edge curve
/// samples and the vertex points, world-space.  A face with no finite
/// sample (the OCCT infinite/open box) reads IsWhole.
fn brep_bnd_lib_add_face(the_face: &Shape, a_b: &mut BndBox) {
    let mut an_any = false;
    // The face vertex points (the vertex walk part, and the projection
    // anchors of the surface re-host).
    let mut a_verts: Vec<rcad_kernel::topo::topology::Vertex> = Vec::new();
    for a_vx in explorer(the_face, ShapeType::Vertex, ShapeType::Shape) {
        if let Some(a_p) = brep_tool_pnt(&a_vx) {
            if !a_p.x.is_finite() || !a_p.y.is_finite() || !a_p.z.is_finite() {
                continue;
            }
            a_verts.push(rcad_kernel::topo::topology::Vertex { point: a_p });
        }
    }
    // The surface part (the OCCT BndLib_AddSurface::Add over the face UV
    // bounds — the None of the unbounded plane re-host is the open box).
    if let Some(a_s) = brep_tool_surface(the_face) {
        if let Some([a_mn, a_mx]) = rcad_kernel::surface_bounding_box(&a_s, &a_verts) {
            if a_mn.x.is_finite() && a_mx.x.is_finite() {
                a_b.add_point(a_mn);
                a_b.add_point(a_mx);
                an_any = true;
            }
        }
    }
    // The edge curve samples (the OCCT wire/edge walk part).
    for a_e in explorer(the_face, ShapeType::Edge, ShapeType::Shape) {
        if let Some((a_c, a_f, a_l)) = brep_tool_curve(&a_e) {
            for a_t in [a_f, 0.5 * (a_f + a_l), a_l] {
                let a_p = CurveEval::point_at(&a_c, a_t);
                if !a_p.x.is_finite() || !a_p.y.is_finite() || !a_p.z.is_finite() {
                    continue;
                }
                a_b.add_point(a_p);
                an_any = true;
            }
        }
    }
    // The vertex points.
    for a_v in &a_verts {
        a_b.add_point(a_v.point);
        an_any = true;
    }
    if !an_any {
        // The OCCT box stayed open in every direction (IsWhole).
        a_b.set_whole();
    }
}

/// The guarded edge-data read: the pool walk for an in-pool edge,
/// `edge_data_pool_free` for the offset-engine edges living outside the
/// pool (no OCCT counterpart — OCCT has a single representation).
fn edge_data<'a>(the_brep: &'a BRep, the_e: &'a Shape) -> &'a TEdgeData {
    if shape_is_in_pool(the_brep, the_e) {
        the_brep.edge(the_e.clone())
    } else {
        edge_data_pool_free(the_e).expect("BRep_Tool: edge without a TEdgeData")
    }
}

/// OCCT BRep_CurveRepresentation::IsCurveOnSurface() — true for the
/// curve-on-surface and closed-surface kinds.
fn gcurve_is_curve_on_surface(a_cr: &CurveRepresentation) -> bool {
    matches!(
        a_cr,
        CurveRepresentation::CurveOnSurface { .. } | CurveRepresentation::CurveOnClosedSurface { .. }
    )
}

/// OCCT BRep_CurveOnSurface::IsCurveOnClosedSurface().
fn gcurve_is_curve_on_closed_surface(a_cr: &CurveRepresentation) -> bool {
    matches!(a_cr, CurveRepresentation::CurveOnClosedSurface { .. })
}

/// OCCT BRep_CurveOnSurface::PCurve().
fn gcurve_pcurve(a_cr: &CurveRepresentation) -> Option<rcad_kernel::geom::Curve2d> {
    match a_cr {
        CurveRepresentation::CurveOnSurface { pcurve, .. } => Some(pcurve.clone()),
        CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => Some(pcurve1.clone()),
        _ => None,
    }
}

/// OCCT BRep_CurveOnClosedSurface::PCurve2().
fn gcurve_pcurve2(a_cr: &CurveRepresentation) -> Option<rcad_kernel::geom::Curve2d> {
    match a_cr {
        CurveRepresentation::CurveOnClosedSurface { pcurve2, .. } => Some(pcurve2.clone()),
        _ => None,
    }
}

/// OCCT BRep_CurveOnSurface::Surface() — the support surface resolved from
/// the owning face key (the same encoding as same_parameter.rs L894-908).
/// Primary read: the pool scan (a pool-resident owning face, including one
/// outside theOldShape's tree).  Fallback: the graph index `the_graph_faces`
/// (built from the InternalUpdateTolerances walk over theOldShape) — a
/// pool-free owning face (index == usize::MAX) has its TShape Arc absent
/// from the pool, and the OCCT handle is carried by the graph; without the
/// fallback the OCCT L1909-1937 distance term would be silently skipped.
fn gcurve_surface(
    the_brep: &BRep,
    a_cr: &CurveRepresentation,
    the_graph_faces: &HashMap<u64, Surface3>,
) -> Option<Surface3> {
    let a_face_key = match a_cr {
        CurveRepresentation::CurveOnSurface { face, .. } => *face,
        CurveRepresentation::CurveOnClosedSurface { face, .. } => *face,
        _ => return None,
    };
    if let Some(a_ts) = the_brep
        .tshapes
        .iter()
        .find(|a_ts| std::sync::Arc::as_ptr(a_ts) as u64 == a_face_key.0)
    {
        if let TShape::Face(a_fd) = a_ts.as_ref() {
            return a_fd.surface.clone();
        }
    }
    the_graph_faces.get(&a_face_key.0).cloned()
}

// =========================================================================
// Tests (hand-derived; see the module docs).
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use glam::{DVec2, DVec3};
    use rcad_kernel::geom::{BSplineCurve3, Curve2d, Curve3, Line2d, Plane};
    use rcad_kernel::topo::topods::TFaceData;

    use crate::brep_algo::tool::{
        builder_add_compound_shape, builder_add_edge_vertex, builder_add_face_wire,
        builder_add_wire_edge, builder_make_compound, builder_make_edge, builder_make_vertex,
        builder_make_wire, oriented,
    };

    /// A straight unit X segment edge (degree-1 BSpline on [0, len]) built
    /// between two vertices, with the raw edge tolerance set.
    fn straight_edge(brep: &mut BRep, p0: DVec3, p1: DVec3, the_tol: f64) -> Shape {
        let a_v1 = brep.add_tvertex(p0);
        let a_v2 = brep.add_tvertex(p1);
        let a_len = (p1 - p0).length();
        let a_bs = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, a_len, a_len],
            control_points: vec![p0, p1],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        };
        let a_e = brep.add_tedge(Some(Curve3::BSpline(a_bs)), a_v1, a_v2, [0.0, a_len]);
        let a_ed = brep.edge_mut_inplace(a_e.clone());
        a_ed.tolerance = the_tol;
        a_e
    }

    /// Two plane faces sharing one edge (the EDGE walk fixture).
    fn two_faces_one_edge(brep: &mut BRep, f1_tol: f64, f2_tol: f64) -> (Shape, Shape, Shape) {
        let a_edge = straight_edge(brep, DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0), CONFUSION);
        // The Curve3D representation (the edge 3d curve slot read of the
        // VERTEX walk).
        let a_ed = brep.edge_mut_inplace(a_edge.clone());
        a_ed.representations.push(CurveRepresentation::Curve3D {
            curve: 0,
            location: 0,
        });        let a_wire1 = brep.add_twire(vec![a_edge.clone()]);
        let a_wire2 = brep.add_twire(vec![a_edge.clone()]);
        let a_f1 = brep.add_tface_tol(
            Some(Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))),
            a_wire1,
            vec![],
            None,
            None,
            vec![],
            false,
            f1_tol,
        );
        let a_f2 = brep.add_tface_tol(
            Some(Surface3::Plane(Plane::new(DVec3::new(1.0, 0.0, 0.0), DVec3::X))),
            a_wire2,
            vec![],
            None,
            None,
            vec![],
            false,
            f2_tol,
        );
        (a_edge, a_f1, a_f2)
    }

    /// The EDGE walk: an edge whose tolerance is BELOW its faces' tolerances
    /// is raised to max(face) (cxx L1832-1856 + UpdShTol).  The OCCT face
    /// table is not engaged (verifyFaceTolerance = false), so the max over
    /// the face tolerances 1e-3 / 2e-3 is 2e-3.
    #[test]
    fn edge_below_faces_is_raised_to_max_face() {
        let mut a_brep = BRep::new();
        let (a_edge, a_f1, a_f2) = two_faces_one_edge(&mut a_brep, 1.0e-3, 2.0e-3);
        let a_root = a_brep.add_tcompound(vec![a_f1, a_f2]);
        update_tolerances(&mut a_brep, &a_root, false);
        // The raw edge tolerance landed exactly on max(face) (the
        // UpdateEdge keep-max write of UpdShTol).
        assert_eq!(a_brep.edge(a_edge.clone()).tolerance, 2.0e-3);
    }

    /// The keep-max discriminant: an edge already ABOVE its faces must stay
    /// untouched ("Update can only increase tolerance", cxx L1849-1852).
    #[test]
    fn edge_above_faces_stays_untouched() {
        let mut a_brep = BRep::new();
        let (a_edge, a_f1, a_f2) = two_faces_one_edge(&mut a_brep, 1.0e-3, 2.0e-3);
        {
            let a_ed = a_brep.edge_mut_inplace(a_edge.clone());
            a_ed.tolerance = 5.0e-3;
        }
        let a_root = a_brep.add_tcompound(vec![a_f1, a_f2]);
        update_tolerances(&mut a_brep, &a_root, false);
        assert_eq!(a_brep.edge(a_edge.clone()).tolerance, 5.0e-3);
    }

    /// The VERTEX walk: a vertex shared by two edges (1e-3 / 3e-3, both
    /// SameRange with the straight 3d curve through the vertex so the
    /// Curve3D distance is 0) is raised to max(edge) + 2*Epsilon (cxx
    /// L1862-1958).
    #[test]
    fn vertex_between_edges_gets_max_edge_tol() {
        let mut a_brep = BRep::new();
        let a_e1 = straight_edge(&mut a_brep, DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0), 1.0e-3);
        let a_e2 = straight_edge(
            &mut a_brep,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            3.0e-3,
        );
        // SameRange + the Curve3D representation on both edges.
        for a_e in [&a_e1, &a_e2] {
            let a_ed = a_brep.edge_mut_inplace(a_e.clone());
            a_ed.same_range = true;
            a_ed.representations.push(CurveRepresentation::Curve3D {
                curve: 0,
                location: 0,
            });
        }
        let a_v_shared = a_brep.edge(a_e1.clone()).last.clone();
        let a_v_first = a_brep.edge(a_e1.clone()).first.clone();
        let a_root = a_brep.add_tcompound(vec![a_e1, a_e2]);
        update_tolerances(&mut a_brep, &a_root, false);

        // The shared vertex: tol = 3e-3 + 2*Epsilon(3e-3) (the force-set
        // arm of UpdShTol, theVForceUpdate = true).
        let a_expected = 3.0e-3 + 2.0 * epsilon_of(3.0e-3);
        let a_tol_shared = a_brep.vertex(a_v_shared.clone()).tolerance;
        assert!(
            (a_tol_shared - a_expected).abs() < 1.0e-18,
            "shared vertex tolerance {} != {}",
            a_tol_shared,
            a_expected
        );
        // The free end vertex: raised to its single edge 1e-3 + 2*Epsilon.
        let a_expected_first = 1.0e-3 + 2.0 * epsilon_of(1.0e-3);
        let a_tol_first = a_brep.vertex(a_v_first.clone()).tolerance;
        assert!(
            (a_tol_first - a_expected_first).abs() < 1.0e-18,
            "end vertex tolerance {} != {}",
            a_tol_first,
            a_expected_first
        );
    }

    /// The vertex keep-max discriminant: a vertex already above max(edge)
    /// stays untouched (anUpdTol false, cxx L1944-1955).
    #[test]
    fn vertex_above_edges_stays_untouched() {
        let mut a_brep = BRep::new();
        let a_e1 = straight_edge(&mut a_brep, DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0), 1.0e-3);
        let a_e2 = straight_edge(
            &mut a_brep,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            3.0e-3,
        );
        for a_e in [&a_e1, &a_e2] {
            let a_ed = a_brep.edge_mut_inplace(a_e.clone());
            a_ed.same_range = true;
        }
        let a_v_shared = a_brep.edge(a_e1.clone()).last.clone();
        // The raw in-place write (keeps the pool Arc identity shared — the
        // Arc::make_mut of BRep::vertex_mut would split it).
        let mut a_v_shared = a_v_shared;
        tvertex_set_tolerance(&mut a_brep, &mut a_v_shared, 9.0e-3, true);
        let a_root = a_brep.add_tcompound(vec![a_e1, a_e2]);
        update_tolerances(&mut a_brep, &a_root, false);
        assert_eq!(a_brep.vertex(a_v_shared.clone()).tolerance, 9.0e-3);
    }

    /// A rectangle wire face on the plane z = 0.
    fn rectangle_face(brep: &mut BRep, a_w: f64, a_h: f64, the_tol: f64) -> Shape {
        let a_e1 = straight_edge(brep, DVec3::ZERO, DVec3::new(a_w, 0.0, 0.0), CONFUSION);
        let a_e2 = straight_edge(
            brep,
            DVec3::new(a_w, 0.0, 0.0),
            DVec3::new(a_w, a_h, 0.0),
            CONFUSION,
        );
        let a_e3 = straight_edge(
            brep,
            DVec3::new(a_w, a_h, 0.0),
            DVec3::new(0.0, a_h, 0.0),
            CONFUSION,
        );
        let a_e4 = straight_edge(brep, DVec3::new(0.0, a_h, 0.0), DVec3::ZERO, CONFUSION);
        let a_wire = brep.add_twire(vec![a_e1, a_e2, a_e3, a_e4]);
        brep.add_tface_tol(
            Some(Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))),
            a_wire,
            vec![],
            None,
            None,
            vec![],
            false,
            the_tol,
        )
    }

    /// The verifyFaceTolerance floor on a PLANE face: the base tol is
    /// Confusion (the cxx L1776-1779 table) scaled by dMax — for the
    /// 2 x 1 rectangle the finite box gives dMax = 2 (the cxx L1790-1814
    /// extents), so the face tolerance is Confusion * 2 exactly (the
    /// UpdateFace plain set).
    #[test]
    fn verify_tolerance_plane_face_floor_is_confusion_times_dmax() {
        let mut a_brep = BRep::new();
        let a_face = rectangle_face(&mut a_brep, 2.0, 1.0, 1.0e-9);
        update_tolerances(&mut a_brep, &a_face, true);
        let a_expected = 2.0 * CONFUSION;
        let a_tol = a_brep.face(a_face.clone()).tolerance;
        assert!(
            (a_tol - a_expected).abs() < 1.0e-18,
            "plane face tolerance {} != Confusion*2 ({})",
            a_tol,
            a_expected
        );
    }

    /// The cxx L1816-1819 cap: dMax large enough that Confusion*dMax > 1
    /// pins the tolerance at exactly 0.99.
    #[test]
    fn verify_tolerance_clamps_at_0_99() {
        let mut a_brep = BRep::new();
        let a_face = rectangle_face(&mut a_brep, 2.0e7, 1.0, 1.0e-9);
        update_tolerances(&mut a_brep, &a_face, true);
        assert_eq!(a_brep.face(a_face.clone()).tolerance, 0.99);
    }

    /// The reshaper overload (cxx L1974-1979, IsMutableInput = false): a
    /// shape untouched by a fresh reshaper takes the empty-copied arm
    /// (cxx L860-906) — the ORIGINAL keeps its tolerance and the reshaper
    /// records the replacement carrying the update.
    #[test]
    fn reshaper_overload_copies_and_replaces() {
        let mut a_brep = BRep::new();
        let (a_edge, a_f1, a_f2) = two_faces_one_edge(&mut a_brep, 1.0e-3, 2.0e-3);
        let a_root = a_brep.add_tcompound(vec![a_f1, a_f2]);
        let mut a_reshaper = ShapeBuildReShape::new();
        update_tolerances_with_reshaper(&mut a_brep, &a_root, &mut a_reshaper, false);
        // The original edge TShape is untouched...
        assert_eq!(a_brep.edge(a_edge.clone()).tolerance, CONFUSION);
        // ...the reshaper recorded the replacement (IsRecorded on the key)
        assert!(a_reshaper.is_recorded(&a_edge));
        // ...and the new shape is flagged by IsNewShape (the recorded VALUE).
        let a_replacement = a_reshaper.value(&mut a_brep, &a_edge);
        assert!(!a_replacement.is_null());
        assert!(a_reshaper.is_new_shape(&a_replacement));
        // The replacement carries the updated tolerance (the max over the
        // faces).
        assert_eq!(a_brep.edge(a_replacement.clone()).tolerance, 2.0e-3);
    }

    /// An all-pool-free two-faces-one-edge fixture (every shape carries
    /// index == usize::MAX and the engine BRep stays EMPTY): face 1 is the
    /// z = 0 plane (the representation's support surface), face 2 carries no
    /// surface; the shared pool-free edge has SameRange and a
    /// CurveOnSurface representation on face 1 whose pcurve at the vertex
    /// parameter 0 maps to the plane point (5, 0, 0), so the VERTEX walk
    /// distance term (cxx L1909-1937) is exactly 5 from the ZERO vertex.
    fn pool_free_two_faces_one_edge(f1_tol: f64, f2_tol: f64) -> (Shape, Shape, Shape) {
        // Face 1: the z = 0 plane carrying the representation's surface.
        let a_f1 = Shape {
            data: Arc::new(TShape::Face(TFaceData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                surface: Some(Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))),
                surface_location: 0,
                outer_wire: Shape::null(),
                inner_wires: Vec::new(),
                sample_point: None,
                uv_domain: None,
                internal_vertices: Vec::new(),
                tolerance: f1_tol,
                natural_restriction: false,
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        };
        // Face 2: no surface (the EDGE walk reads only its tolerance).
        let a_f2 = Shape {
            data: Arc::new(TShape::Face(TFaceData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                surface: None,
                surface_location: 0,
                outer_wire: Shape::null(),
                inner_wires: Vec::new(),
                sample_point: None,
                uv_domain: None,
                internal_vertices: Vec::new(),
                tolerance: f2_tol,
                natural_restriction: false,
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        };
        // The vertices (ZERO point) and the shared pool-free edge.
        let a_v1 = builder_make_vertex();
        let a_v2 = builder_make_vertex();
        let mut a_edge = builder_make_edge();
        builder_add_edge_vertex(&mut a_edge, &a_v1);
        builder_add_edge_vertex(&mut a_edge, &oriented(&a_v2, Orientation::Reversed));
        {
            let a_te = Arc::make_mut(&mut a_edge.data);
            if let TShape::Edge(a_ed) = a_te {
                a_ed.tolerance = CONFUSION;
                a_ed.same_range = true;
                a_ed.representations.push(CurveRepresentation::CurveOnSurface {
                    face: shape_key(&a_f1),
                    pcurve: Curve2d::Line(Line2d::new(
                        DVec2::new(5.0, 0.0),
                        DVec2::new(1.0, 0.0),
                    )),
                    range: [0.0, 1.0],
                });
            }
        }
        // The wires and the faces (each face holds the shared edge once).
        let mut a_w1 = builder_make_wire();
        builder_add_wire_edge(&mut a_w1, &a_edge);
        let mut a_w2 = builder_make_wire();
        builder_add_wire_edge(&mut a_w2, &a_edge);
        let mut a_f1 = a_f1;
        builder_add_face_wire(&mut a_f1, &a_w1);
        let mut a_f2 = a_f2;
        builder_add_face_wire(&mut a_f2, &a_w2);
        // The root compound.
        let mut a_root = builder_make_compound();
        builder_add_compound_shape(&mut a_root, &a_f1);
        builder_add_compound_shape(&mut a_root, &a_f2);
        (a_root, a_edge, a_v1)
    }

    /// The pool-free mirror of `edge_below_faces_is_raised_to_max_face`:
    /// an all-pool-free graph (every shape built pool-free, never pooled;
    /// the engine BRep stays EMPTY) — the harmonization ENGAGES exactly as
    /// on the pool-resident fixture.  The EDGE below its faces' tolerances
    /// is raised to max(face) = 2e-3 (the EDGE walk keep-max write,
    /// observed through the same TShape — the in-place write reaches every
    /// handle), and the VERTEX gets 5 + 2*Epsilon(5) — the CurveOnSurface
    /// distance term (cxx L1909-1937) engaged through the pool-free
    /// owning-face resolution plus the force-set vertex write
    /// (theVForceUpdate = true).
    #[test]
    fn pool_free_edge_below_faces_is_raised_to_max_face() {
        let mut a_brep = BRep::new();
        let (a_root, a_edge, a_v1) = pool_free_two_faces_one_edge(1.0e-3, 2.0e-3);
        update_tolerances(&mut a_brep, &a_root, false);
        // The EDGE walk: the raw pool-free edge tolerance landed exactly on
        // max(face) = 2e-3 (the UpdateEdge keep-max write of UpdShTol).
        assert_eq!(
            edge_data_pool_free(&a_edge)
                .expect("pool-free edge")
                .tolerance,
            2.0e-3
        );
        // The VERTEX walk: the pcurve-on-plane distance term is 5 (the
        // pcurve at the vertex parameter 0 maps to the plane point
        // (5, 0, 0), the ZERO vertex), so tol = 5 + 2*Epsilon(5) — the
        // force-set arm of UpdShTol.
        let a_expected = 5.0 + 2.0 * epsilon_of(5.0);
        let a_v_tol = match a_v1.data.as_ref() {
            TShape::Vertex(a_vd) => a_vd.tolerance,
            _ => panic!("the fixture vertex is not a vertex"),
        };
        assert!(
            (a_v_tol - a_expected).abs() < 1.0e-18,
            "pool-free vertex tolerance {} != {}",
            a_v_tol,
            a_expected
        );
    }
}
