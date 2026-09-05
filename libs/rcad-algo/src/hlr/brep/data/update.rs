// OCCT HLRBRep_Data.cxx L522-1229 (the constructor / Destroy / Write /
// Update / InitBoundSort / InitEdge / MoreEdge / NextEdge / Edge /
// InitInterference group) + HLRBRep_Data.lxx L19-120 (the inline
// accessors).  Impl block of super::Data (proxy L); the interference /
// classification bodies from L1233 on live in [classify].
//
// Conventions of the Data object graph (established with the sibling
// modules):
// - the OCCT element pointers (`iFaceData`, `myLEData`, `iFaceMinMax`,
//   ...) keep the raw-pointer form and are dereferenced in small unsafe
//   faces annotated with the OCCT statement;
// - the (0,N) NCollection_Array1 members keep one unused slot in OCCT;
//   the rcad Vecs drop it (the repo 1-based -> 0-based shift), so the
//   1-based index i accesses `[(i - 1)]`;
// - the OCCT `HLRBRep_CurvePtr` / `HLRBRep_SurfacePtr` tool parameters
//   keep the raw pointers;
// - the OCCT uninitialized members / out parameters keep neutral defaults
//   documented at the declaration site.

use std::sync::Arc;

use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::precision::{INFINITE_VALUE, PCONFUSION};
use rcad_kernel::topods::{Orientation, Shape};

use crate::bop::int_tools::bean_face_intersector::GeomAbsCurveType;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::hlr::algo::edges_block::MinMaxIndices;
use crate::hlr::algo::hlr_algo::HLRAlgo;
use crate::hlr::algo::interference::Interference;
use crate::hlr::algo::projector::Projector;
use crate::hlr::algo::wires_block::WiresBlock;
use crate::hlr::brep::cl_props::CLProps;
use crate::hlr::brep::curve::Curve;
use crate::hlr::brep::edge_data::EdgeData;
use crate::hlr::brep::face_data::FaceData;
use crate::hlr::brep::face_iterator::FaceIterator;
use crate::hlr::brep::intersector::Intersector;
use crate::hlr::brep::sl_props::SLProps;
use crate::hlr::brep::surface::Surface;
use crate::topalgo::brep_top_adaptor::tool::BRepTopAdaptorTool;
use crate::topalgo::brep_top_adaptor::topol_tool_brep::BRepTopolTool;

use super::tableau_rejection::TableauRejection;
use super::Data;

/// OCCT `& 0x80008000` (cxx L997) — the int bit test of the encoded box
/// words; the literal needs the u32->i32 form (overflowing_literals).
const MASK_80008000: i32 = 0x80008000u32 as i32;

/// OCCT `myLEType = myLEGeom->GetType()` — the GeomAbs_CurveType of the
/// projected curve; rcad Curve reports the proj_lib CurveType twin (the
/// same OCCT enum carried by two module-local aliases).
fn curve_type_to_geom_abs(t: CurveType) -> GeomAbsCurveType {
    match t {
        CurveType::Line => GeomAbsCurveType::Line,
        CurveType::Circle => GeomAbsCurveType::Circle,
        CurveType::Ellipse => GeomAbsCurveType::Ellipse,
        CurveType::Hyperbola => GeomAbsCurveType::Hyperbola,
        CurveType::Parabola => GeomAbsCurveType::Parabola,
        CurveType::Bezier => GeomAbsCurveType::BezierCurve,
        CurveType::BSpline => GeomAbsCurveType::BSplineCurve,
        CurveType::Other => GeomAbsCurveType::OtherCurve,
    }
}

/// OCCT NCollection_IndexedMap::Add — the insertion-ordered Vec keyed by
/// Shape::ptr_id (the topo_brep/data.rs map precedent).
fn indexed_map_add(map: &mut Vec<Shape>, s: Shape) {
    let id = s.ptr_id();
    if !map.iter().any(|k| k.ptr_id() == id) {
        map.push(s);
    }
}

/// OCCT NCollection_DataMap::IsBound over the MST map (the ds_filler
/// precedent: the insertion-ordered Vec keyed by Shape::ptr_id).
fn mst_is_bound(mst: &[(Shape, BRepTopAdaptorTool)], s: &Shape) -> bool {
    let id = s.ptr_id();
    mst.iter().any(|(k, _)| k.ptr_id() == id)
}

/// The MST map position of the shape — OCCT MST.ChangeFind / Find.
fn mst_pos(mst: &[(Shape, BRepTopAdaptorTool)], s: &Shape) -> usize {
    let id = s.ptr_id();
    mst.iter()
        .position(|(k, _)| k.ptr_id() == id)
        .expect("Standard_NoSuchObject: MST key not bound")
}

/// OCCT MST.Bind(S, BRT) — insert the tool under the shape key.
fn mst_bind(mst: &mut Vec<(Shape, BRepTopAdaptorTool)>, s: &Shape, tool: BRepTopAdaptorTool) {
    mst.push((s.clone(), tool));
}

impl<'a> Data<'a> {
    /// OCCT HLRBRep_Data::HLRBRep_Data(NV, NE, NF) (cxx L522-537) — create
    /// an empty data structure of NV vertices, NE edges and NF faces.
    pub fn new(nb_vertices: usize, nb_edges: usize, nb_faces: usize) -> Data<'a> {
        // The member init list (cxx L523-533).  The members without an OCCT
        // initializer stay uninitialized in C++; the neutral defaults are
        // documented per field.
        let mut data = Data {
            my_nb_vertices: nb_vertices, // myNbVertices(NV)
            my_nb_edges: nb_edges,       // myNbEdges(NE)
            my_nb_faces: nb_faces,       // myNbFaces(NF)
            my_e_map: Vec::new(),        // myEMap — the empty IndexedMap
            my_f_map: Vec::new(),        // myFMap
            // myEData(0, NE) — the OCCT (0,NE) Array1 keeps an unused slot
            // 0; the rcad Vec drops it (the 1-based -> 0-based shift).
            my_e_data: (0..nb_edges).map(|_| EdgeData::new()).collect(),
            my_f_data: (0..nb_faces).map(|_| FaceData::new()).collect(), // myFData(0, NF)
            my_edge_indices: vec![0; nb_edges],                          // myEdgeIndices(0, NE)
            my_toler: 1.0e-5,            // myToler((float)1e-5)
            my_proj: Projector::new(),   // myProj — the default projector
            my_l_l_props: CLProps::new(2, f64::EPSILON), // myLLProps(2, Epsilon(1.))
            my_f_l_props: CLProps::new(2, f64::EPSILON), // myFLProps(2, Epsilon(1.))
            my_s_l_props: SLProps::new(2, f64::EPSILON), // mySLProps(2, Epsilon(1.))
            my_big_size: 0.0,            // uninitialized in OCCT
            // myFaceItr1 / myFaceItr2 — the default FaceIterator (the null
            // Wires handle; InitEdge binds fd.Wires()).
            my_face_itr1: FaceIterator::new(std::ptr::null_mut()),
            my_face_itr2: FaceIterator::new(std::ptr::null_mut()),
            i_face: 0,                         // uninitialized in OCCT
            i_face_data: std::ptr::null_mut(), // the null element pointers
            i_face_geom: std::ptr::null_mut(),
            i_face_min_max: std::ptr::null_mut(),
            i_face_type: GeomAbsSurfaceType::OtherSurface, // uninitialized enum
            i_face_back: false,
            i_face_simp: false,
            i_face_smpl: false,
            i_face_test: false,
            my_hide_count: 0,   // myHideCount(0)
            my_deca: [0.0; 16], // uninitialized in OCCT
            my_sur_d: [0.0; 16],
            my_cur_sort_ed: 0, // uninitialized in OCCT
            my_nbr_sort_ed: 0,
            my_le: 0,
            my_le_out_line: false,
            my_le_internal: false,
            my_le_double: false,
            my_le_iso_line: false,
            my_le_data: std::ptr::null_mut(),
            my_le_geom: std::ptr::null(),
            my_le_min_max: std::ptr::null_mut(),
            my_le_type: GeomAbsCurveType::OtherCurve, // uninitialized enum
            my_le_tol: 0.0,
            my_fe: 0,
            my_fe_ori: Orientation::Forward,
            my_fe_out_line: false,
            my_fe_internal: false,
            my_fe_double: false,
            my_fe_data: std::ptr::null_mut(),
            my_fe_geom: std::ptr::null_mut(),
            my_fe_type: GeomAbsCurveType::OtherCurve,
            my_fe_tol: 0.0,
            my_intersector: Intersector::new(), // myIntersector — the default
            // myClassifier — the null handle; the empty BRep stands in (the
            // BRepTopAdaptor_Tool() precedent).
            my_classifier: BRepTopolTool::new(Arc::new(rcad_kernel::BRep::new())),
            my_same_vertex: false, // uninitialized in OCCT
            my_intersected: false,
            my_nb_points: 0,
            my_nb_segments: 0,
            i_interf: 0,
            my_intf: Interference::new(), // the default interference
            my_above_intf: false,
            my_reject: std::ptr::null_mut(), // set below (new TableauRejection())
            // my_brep — no OCCT counterpart (the kernel-context deviation of
            // [Data]); the neutral default is the null context, injected by
            // Update / InitEdge.
            my_brep: None,
        };
        // myReject = new TableauRejection();
        let mut my_reject = Box::new(TableauRejection::new());
        // ((TableauRejection*)myReject)->SetDim(myNbEdges);
        my_reject.set_dim(nb_edges as i32);
        data.my_reject = Box::into_raw(my_reject);
        data
    }

    /// OCCT HLRBRep_Data::Destroy (cxx L539-544).
    pub fn destroy(&mut self) {
        // ((TableauRejection*)myReject)->Destroy();
        unsafe { (*self.my_reject).destroy() };
        // delete ((TableauRejection*)myReject);
        drop(unsafe { Box::from_raw(self.my_reject) });
        // the rcad guard of the OCCT one-shot destructor semantics.
        self.my_reject = std::ptr::null_mut();
    }

    /// OCCT HLRBRep_Data::Write (cxx L548-608) — write DS in me with a
    /// translation of dv, de, df.
    pub fn write(&mut self, ds: &mut Data<'a>, dv: i32, de: i32, df: i32) {
        let n1edge = ds.nb_edges(); // int n1edge = DS->NbEdges();
        let n1face = ds.nb_faces(); // int n1face = DS->NbFaces();

        // HLRBRep_EdgeData* ed = &(myEData.ChangeValue(de)); ed++;
        // (ChangeValue(de) is the 1-based slot de; the ++ steps to the
        // first written slot, 1-based de+1 = vec index de.)
        let mut ed: *mut EdgeData<'a> = &mut self.my_e_data[de as usize];
        // HLRBRep_EdgeData* e1 = &(DS->EDataArray().ChangeValue(0)); e1++;
        // (the OCCT (0,NE) array's unused slot 0; the ++ lands on 1-based
        // edge 1 = vec index 0 — the two steps fuse.)
        let mut e1: *mut EdgeData<'a> = &mut ds.my_e_data[0];

        // HLRBRep_FaceData* fd = &(myFData.ChangeValue(df)); fd++;
        let mut fd: *mut FaceData<'a> = &mut self.my_f_data[df as usize];
        // HLRBRep_FaceData* f1 = &(DS->FDataArray().ChangeValue(0)); f1++;
        let mut f1: *mut FaceData<'a> = &mut ds.my_f_data[0];

        // for (int iedge = 1; iedge <= n1edge; iedge++)
        let mut iedge: i32 = 1;
        while iedge <= n1edge as i32 {
            // *ed = *e1;
            unsafe { *ed = (*e1).clone() };

            if dv != 0 {
                // ed->VSta(ed->VSta() + dv); ed->VEnd(ed->VEnd() + dv);
                unsafe {
                    (*ed).set_v_sta((*ed).v_sta() + dv);
                    (*ed).set_v_end((*ed).v_end() + dv);
                }
            }

            // myEMap.Add(DS->EdgeMap().FindKey(iedge));
            let key = ds.my_e_map[(iedge - 1) as usize].clone();
            indexed_map_add(&mut self.my_e_map, key);

            ed = unsafe { ed.add(1) };
            e1 = unsafe { e1.add(1) };
            iedge += 1;
        }

        // for (int iface = 1; iface <= n1face; iface++)
        let mut iface: i32 = 1;
        while iface <= n1face as i32 {
            // *fd = *f1;
            unsafe { *fd = (*f1).clone() };

            if de != 0 {
                // const occ::handle<HLRAlgo_WiresBlock>& wb = fd->Wires();
                let wb: *mut WiresBlock = unsafe { (*fd).wires() };
                // int nw = wb->NbWires();
                let nw = unsafe { (*wb).nb_wires() };

                // for (int iw = 1; iw <= nw; iw++)
                let mut iw: i32 = 1;
                while iw <= nw as i32 {
                    // const occ::handle<HLRAlgo_EdgesBlock>& eb = wb->Wire(iw);
                    let eb = unsafe { (*wb).wire(iw as usize) };
                    // int ne = eb->NbEdges();
                    let ne = eb.nb_edges();

                    // for (int ie = 1; ie <= ne; ie++)
                    let mut ie: i32 = 1;
                    while ie <= ne as i32 {
                        // eb->Edge(ie, eb->Edge(ie) + de);
                        eb.set_edge(ie as usize, eb.edge(ie as usize) + de);
                        ie += 1;
                    }
                    iw += 1;
                }
            }

            // myFMap.Add(DS->FaceMap().FindKey(iface));
            let key = ds.my_f_map[(iface - 1) as usize].clone();
            indexed_map_add(&mut self.my_f_map, key);

            fd = unsafe { fd.add(1) };
            f1 = unsafe { f1.add(1) };
            iface += 1;
        }
    }

    /// OCCT HLRBRep_Data::Update (cxx L612-980) — end of building of the
    /// Data and updating of all the information linked to the projection.
    ///
    /// rcad deviation: the owning BRep travels alongside (the
    /// ds_filler::insert precedent) — the EdgeFaceTool::UVPoint read of the
    /// global TShape graph needs it; it is stored into [Data::my_brep] (the
    /// kernel-context slot) for the whole hiding phase.
    pub fn update(&mut self, brep: &'a rcad_kernel::BRep, p: &Projector) {
        // rcad: the kernel-context store (no OCCT statement).
        self.my_brep = Some(brep);
        // myProj = P;
        self.my_proj = p.clone();
        // const gp_Trsf& T = myProj.Transformation();
        let t = *self.my_proj.transformation();
        let mut tol_min_max: f64 = 0.0;

        // The MinMax work set (cxx L619-624) — uninitialized structs in
        // OCCT; the neutral all-zero default stands in.
        let mut face_min = MinMaxIndices::default();
        let mut face_max = MinMaxIndices::default();
        let mut min_max_face = MinMaxIndices::default();
        let mut wire_min = MinMaxIndices::default();
        let mut wire_max = MinMaxIndices::default();
        let mut min_max_wire = MinMaxIndices::default();
        let mut edge_min = MinMaxIndices::default();
        let mut edge_max = MinMaxIndices::default();
        let mut min_max_edge = MinMaxIndices::default();
        let mut tot_min = [0.0f64; 16];
        let mut tot_max = [0.0f64; 16];

        // HLRAlgo::InitMinMax(Precision::Infinite(), TotMin, TotMax);
        HLRAlgo::init_min_max(INFINITE_VALUE, &mut tot_min, &mut tot_max);

        // compute the global MinMax
        // *************************
        // int edge; for (edge = 1; edge <= myNbEdges; edge++)
        let mut edge: i32 = 1;
        while edge <= self.my_nb_edges as i32 {
            // HLRBRep_EdgeData& ed = myEData.ChangeValue(edge);
            let ed = &mut self.my_e_data[(edge - 1) as usize];
            // HLRBRep_Curve& EC = ed.ChangeGeometry();
            let ec = ed.change_geometry();
            // EC.Projector(&myProj);
            ec.projector(&self.my_proj);
            // double enl = EC.Update(TotMin, TotMax);
            let enl = ec.update(&mut tot_min, &mut tot_max);
            if enl > tol_min_max {
                tol_min_max = enl;
            }
            edge += 1;
        }
        HLRAlgo::enlarge_min_max(tol_min_max, &mut tot_min, &mut tot_max);
        let mut d = [0.0f64; 16];
        let mut precad = -INFINITE_VALUE; // -Precision::Infinite()

        let mut i: usize = 0;
        while i <= 15 {
            d[i] = tot_max[i] - tot_min[i];
            if precad < d[i] {
                precad = d[i];
            }
            i += 1;
        }
        self.my_big_size = precad;
        precad *= 0.0005;

        let mut i: usize = 0;
        while i <= 15 {
            // mySurD[i] = 0x00007fff / (d[i] + precad);
            self.my_sur_d[i] = 0x00007fff as f64 / (d[i] + precad);
            i += 1;
        }
        precad *= 0.5;

        let mut i: usize = 0;
        while i <= 15 {
            // myDeca[i] = -TotMin[i] + precad;
            self.my_deca[i] = -tot_min[i] + precad;
            i += 1;
        }

        let mut ver1: bool;
        let mut ver2: bool;

        // update the edges
        // ****************
        let mut edge: i32 = 1;
        while edge <= self.my_nb_edges as i32 {
            // HLRBRep_EdgeData& ed = myEData.ChangeValue(edge);
            let ed: *mut EdgeData<'a> = &mut self.my_e_data[(edge - 1) as usize];
            // HLRBRep_Curve& EC = ed.ChangeGeometry();
            let ec = unsafe { (*ed).change_geometry() };
            HLRAlgo::init_min_max(INFINITE_VALUE, &mut tot_min, &mut tot_max);
            // tolMinMax = EC.UpdateMinMax(TotMin, TotMax);
            let tol_min_max = ec.update_min_max(&mut tot_min, &mut tot_max);
            // tol = (double)(ed.Tolerance());
            let tol = unsafe { (*ed).tolerance() } as f64;
            // ed.Vertical(TotMax[0] - TotMin[0] < tol && ... && TotMax[6] - TotMin[6] < tol);
            let vertical = tot_max[0] - tot_min[0] < tol
                && tot_max[1] - tot_min[1] < tol
                && tot_max[2] - tot_min[2] < tol
                && tot_max[3] - tot_min[3] < tol
                && tot_max[4] - tot_min[4] < tol
                && tot_max[5] - tot_min[5] < tol
                && tot_max[6] - tot_min[6] < tol;
            unsafe { (*ed).set_vertical(vertical) };
            HLRAlgo::enlarge_min_max(tol_min_max, &mut tot_min, &mut tot_max);
            // Linux warning : assignment to `int' from `double'. Cast has been added.
            unsafe {
                edge_min.min[0] = ((self.my_deca[0] + tot_min[0]) * self.my_sur_d[0]) as i32;
                edge_max.min[0] = ((self.my_deca[0] + tot_max[0]) * self.my_sur_d[0]) as i32;
                edge_min.min[1] = ((self.my_deca[1] + tot_min[1]) * self.my_sur_d[1]) as i32;
                edge_max.min[1] = ((self.my_deca[1] + tot_max[1]) * self.my_sur_d[1]) as i32;
                edge_min.min[2] = ((self.my_deca[2] + tot_min[2]) * self.my_sur_d[2]) as i32;
                edge_max.min[2] = ((self.my_deca[2] + tot_max[2]) * self.my_sur_d[2]) as i32;
                edge_min.min[3] = ((self.my_deca[3] + tot_min[3]) * self.my_sur_d[3]) as i32;
                edge_max.min[3] = ((self.my_deca[3] + tot_max[3]) * self.my_sur_d[3]) as i32;
                edge_min.min[4] = ((self.my_deca[4] + tot_min[4]) * self.my_sur_d[4]) as i32;
                edge_max.min[4] = ((self.my_deca[4] + tot_max[4]) * self.my_sur_d[4]) as i32;
                edge_min.min[5] = ((self.my_deca[5] + tot_min[5]) * self.my_sur_d[5]) as i32;
                edge_max.min[5] = ((self.my_deca[5] + tot_max[5]) * self.my_sur_d[5]) as i32;
                edge_min.min[6] = ((self.my_deca[6] + tot_min[6]) * self.my_sur_d[6]) as i32;
                edge_max.min[6] = ((self.my_deca[6] + tot_max[6]) * self.my_sur_d[6]) as i32;
                edge_min.min[7] = ((self.my_deca[7] + tot_min[7]) * self.my_sur_d[7]) as i32;
                edge_max.min[7] = ((self.my_deca[7] + tot_max[7]) * self.my_sur_d[7]) as i32;
                edge_min.max[0] = ((self.my_deca[8] + tot_min[8]) * self.my_sur_d[8]) as i32;
                edge_max.max[0] = ((self.my_deca[8] + tot_max[8]) * self.my_sur_d[8]) as i32;
                edge_min.max[1] = ((self.my_deca[9] + tot_min[9]) * self.my_sur_d[9]) as i32;
                edge_max.max[1] = ((self.my_deca[9] + tot_max[9]) * self.my_sur_d[9]) as i32;
                edge_min.max[2] = ((self.my_deca[10] + tot_min[10]) * self.my_sur_d[10]) as i32;
                edge_max.max[2] = ((self.my_deca[10] + tot_max[10]) * self.my_sur_d[10]) as i32;
                edge_min.max[3] = ((self.my_deca[11] + tot_min[11]) * self.my_sur_d[11]) as i32;
                edge_max.max[3] = ((self.my_deca[11] + tot_max[11]) * self.my_sur_d[11]) as i32;
                edge_min.max[4] = ((self.my_deca[12] + tot_min[12]) * self.my_sur_d[12]) as i32;
                edge_max.max[4] = ((self.my_deca[12] + tot_max[12]) * self.my_sur_d[12]) as i32;
                edge_min.max[5] = ((self.my_deca[13] + tot_min[13]) * self.my_sur_d[13]) as i32;
                edge_max.max[5] = ((self.my_deca[13] + tot_max[13]) * self.my_sur_d[13]) as i32;
                edge_min.max[6] = ((self.my_deca[14] + tot_min[14]) * self.my_sur_d[14]) as i32;
                edge_max.max[6] = ((self.my_deca[14] + tot_max[14]) * self.my_sur_d[14]) as i32;
                edge_min.max[7] = ((self.my_deca[15] + tot_min[15]) * self.my_sur_d[15]) as i32;
                edge_max.max[7] = ((self.my_deca[15] + tot_max[15]) * self.my_sur_d[15]) as i32;
            }

            HLRAlgo::encode_min_max(&edge_min, &edge_max, &mut min_max_edge);
            // ed.UpdateMinMax(MinMaxEdge);
            unsafe { (*ed).update_min_max(&min_max_edge) };
            if unsafe { (*ed).vertical() } {
                ver1 = true;
                ver2 = true;
                // int vsta = ed.VSta(); int vend = ed.VEnd();
                let vsta = unsafe { (*ed).v_sta() };
                let vend = unsafe { (*ed).v_end() };
                // bool vout = ed.OutLVSta() || ed.OutLVEnd();
                let vout = unsafe { (*ed).out_lv_sta() } || unsafe { (*ed).out_lv_end() };
                // bool vcut = ed.CutAtSta() || ed.CutAtEnd();
                let vcut = unsafe { (*ed).cut_at_sta() } || unsafe { (*ed).cut_at_end() };

                // for (int ebis = 1; ebis <= myNbEdges; ebis++)
                let mut ebis: i32 = 1;
                while ebis <= self.my_nb_edges as i32 {
                    // HLRBRep_EdgeData& eb = myEData.ChangeValue(ebis);
                    let eb = &mut self.my_e_data[(ebis - 1) as usize];
                    if vsta == eb.v_sta() {
                        eb.set_v_sta(vend);
                        eb.set_out_lv_sta(vout);
                        eb.set_cut_at_sta(vcut);
                    } else if vsta == eb.v_end() {
                        eb.set_v_end(vend);
                        eb.set_out_lv_end(vout);
                        eb.set_cut_at_end(vcut);
                    }
                    ebis += 1;
                }
            } else {
                // gp_Pnt Pt; gp_Vec Tg1, Tg2;
                let mut pt = rcad_kernel::geom::Point3::ZERO;
                // EC.D1(EC.Parameter3d(EC.FirstParameter()), Pt, Tg1);
                let r1 = ec.d1_3d(ec.parameter_3d(ec.first_parameter()));
                pt = r1.0;
                let mut tg1 = r1.1;
                // EC.D1(EC.Parameter3d(EC.LastParameter()), Pt, Tg2);
                let r2 = ec.d1_3d(ec.parameter_3d(ec.last_parameter()));
                pt = r2.0;
                let mut tg2 = r2.1;
                // the OCCT Pt is written and never read (kept verbatim).
                let _ = pt;
                // Tg1.Transform(T); Tg2.Transform(T);
                tg1 = t.transform_vec(tg1);
                tg2 = t.transform_vec(tg2);
                if (tg1.x.abs() + tg1.y.abs()) < self.my_toler as f64 * 10.0 {
                    ver1 = true;
                } else {
                    // gp_Dir Dir1(Tg1);
                    let dir1 = tg1.normalize();
                    ver1 = (dir1.x.abs() + dir1.y.abs()) < self.my_toler as f64 * 10.0;
                }
                if (tg2.x.abs() + tg2.y.abs()) < self.my_toler as f64 * 10.0 {
                    ver2 = true;
                } else {
                    // gp_Dir Dir2(Tg2);
                    let dir2 = tg2.normalize();
                    ver2 = (dir2.x.abs() + dir2.y.abs()) < self.my_toler as f64 * 10.0;
                }
            }
            unsafe {
                // ed.VerAtSta(ed.Vertical() || ver1);
                (*ed).set_ver_at_sta((*ed).vertical() || ver1);
                // ed.VerAtEnd(ed.Vertical() || ver2);
                (*ed).set_ver_at_end((*ed).vertical() || ver2);
                // ed.AutoIntersectionDone(true);
                (*ed).set_auto_intersection_done(true);
                // ed.Simple(true);
                (*ed).set_simple(true);
            }
            edge += 1;
        }

        // update the faces
        // ****************

        // for (int face = 1; face <= myNbFaces; face++)
        let mut face: i32 = 1;
        while face <= self.my_nb_faces as i32 {
            // HLRBRep_FaceData& fd = myFData.ChangeValue(face);
            let fd: *mut FaceData<'a> = &mut self.my_f_data[(face - 1) as usize];
            // HLRBRep_Surface& FS = fd.Geometry(); iFaceGeom = &(fd.Geometry());
            let fs = unsafe { (*fd).geometry() };
            self.i_face_geom = fs as *mut Surface<'a>;
            // mySLProps.SetSurface(iFaceGeom);
            self.my_s_l_props.set_surface(unsafe { &*self.i_face_geom });
            // FS.Projector(&myProj);
            fs.projector(&self.my_proj);
            // iFaceType = FS.GetType();
            self.i_face_type = fs.get_type();

            // Is the face cut by an outline

            let mut cut = false;
            let mut with_out_l = false;

            // for (myFaceItr1.InitEdge(fd); myFaceItr1.MoreEdge(); myFaceItr1.NextEdge())
            self.my_face_itr1.init_edge(fd);
            while self.my_face_itr1.more_edge() {
                if self.my_face_itr1.internal() {
                    with_out_l = true;
                    cut = true;
                } else if self.my_face_itr1.out_line() {
                    with_out_l = true;
                    if self.my_face_itr1.double() {
                        cut = true;
                    }
                }
                self.my_face_itr1.next_edge();
            }
            // fd.Cut(cut); fd.WithOutL(withOutL);
            unsafe { (*fd).set_cut(cut) };
            unsafe { (*fd).set_with_out_l(with_out_l) };

            // Is the face simple = no auto-hiding
            // not cut and simple surface

            let simple = !with_out_l
                && (self.i_face_type == GeomAbsSurfaceType::Plane
                    || self.i_face_type == GeomAbsSurfaceType::Cylinder
                    || self.i_face_type == GeomAbsSurfaceType::Cone
                    || self.i_face_type == GeomAbsSurfaceType::Sphere
                    || self.i_face_type == GeomAbsSurfaceType::Torus);
            // fd.Simple(simple);
            unsafe { (*fd).set_simple(simple) };

            // fd.Plane/Cylinder/Cone/Sphere/Torus(iFaceType == ...);
            unsafe {
                (*fd).set_plane(self.i_face_type == GeomAbsSurfaceType::Plane);
                (*fd).set_cylinder(self.i_face_type == GeomAbsSurfaceType::Cylinder);
                (*fd).set_cone(self.i_face_type == GeomAbsSurfaceType::Cone);
                (*fd).set_sphere(self.i_face_type == GeomAbsSurfaceType::Sphere);
                (*fd).set_torus(self.i_face_type == GeomAbsSurfaceType::Torus);
            }
            // tol = (double)(fd.Tolerance());
            let tol = unsafe { (*fd).tolerance() } as f64;
            // fd.Side(FS.IsSide(tol, myToler * 10));
            unsafe { (*fd).set_side(fs.is_side(tol, self.my_toler as f64 * 10.0)) };
            let mut inverted = false;
            if unsafe { (*fd).with_out_l() } && !unsafe { (*fd).side() } {
                // inverted = OrientOutLine(face, fd);
                inverted = self.orient_out_line(face, unsafe { &mut *fd });
                // OrientOthEdge(face, fd);
                self.orient_oth_edge(face, unsafe { &mut *fd });
            }
            if unsafe { (*fd).side() } {
                unsafe { (*fd).set_hiding(false) };
                unsafe { (*fd).set_back(false) };
            } else if !unsafe { (*fd).with_out_l() } {
                // double p, pu, pv, r;
                let mut p: f64;
                let mut pu: f64 = 0.0; // OCCT uninitialized; neutral default 0.
                let mut pv: f64 = 0.0; // OCCT uninitialized; neutral default 0.
                unsafe { (*fd).set_back(false) };
                let mut found = false;

                // for (myFaceItr1.InitEdge(fd); myFaceItr1.MoreEdge() && !found; myFaceItr1.NextEdge())
                self.my_face_itr1.init_edge(fd);
                while self.my_face_itr1.more_edge() && !found {
                    // myFE = myFaceItr1.Edge();
                    self.my_fe = self.my_face_itr1.edge() as usize;
                    self.my_fe_ori = self.my_face_itr1.orientation();
                    self.my_fe_out_line = self.my_face_itr1.out_line();
                    self.my_fe_internal = self.my_face_itr1.internal();
                    self.my_fe_double = self.my_face_itr1.double();
                    // HLRBRep_EdgeData& EDataFE1 = myEData(myFE);
                    let edatafe1 = &mut self.my_e_data[(self.my_fe - 1) as usize];
                    if !self.my_fe_double
                        && (self.my_fe_ori == Orientation::Forward
                            || self.my_fe_ori == Orientation::Reversed)
                    {
                        // myFEGeom = &(EDataFE1.ChangeGeometry());
                        self.my_fe_geom = edatafe1.change_geometry() as *mut Curve<'a>;
                        // const HLRBRep_Curve& EC = EDataFE1.Geometry();
                        let ec = edatafe1.geometry();
                        // p = EC.Parameter3d((EC.LastParameter() + EC.FirstParameter()) / 2);
                        p = ec.parameter_3d((ec.last_parameter() + ec.first_parameter()) / 2.0);
                        // if (HLRBRep_EdgeFaceTool::UVPoint(p, myFEGeom, iFaceGeom, pu, pv))
                        if crate::hlr::brep::edge_face_tool::uv_point(
                            self.my_brep.expect("HLRBRep_Data::Update: no kernel context"),
                            p,
                            self.my_fe_geom,
                            &self.my_e_map[(self.my_fe - 1) as usize],
                            self.i_face_geom,
                            &mut pu,
                            &mut pv,
                        ) {
                            // mySLProps.SetParameters(pu, pv);
                            self.my_s_l_props.set_parameters(pu, pv);
                            // gp_Pnt Pt; Pt = EC.Value3D(p);
                            let mut pt = ec.value_3d(p);
                            if self.my_s_l_props.is_normal_defined() {
                                // gp_Vec Nm = mySLProps.Normal();
                                let mut nm = self
                                    .my_s_l_props
                                    .normal()
                                    .expect("StdFail_UndefinedNormal");
                                pt = t.apply(pt);
                                nm = t.transform_vec(nm);
                                // r = perspective ? Nm.Z()*Focus() - (Nm.X()*Pt.X() + ...) : Nm.Z();
                                let r = if self.my_proj.perspective() {
                                    nm.z * self.my_proj.focus()
                                        - (nm.x * pt.x + nm.y * pt.y + nm.z * pt.z)
                                } else {
                                    nm.z
                                };
                                if r.abs() > self.my_toler as f64 * 10.0 {
                                    // fd.Back(r < 0);
                                    unsafe { (*fd).set_back(r < 0.0) };
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                    self.my_face_itr1.next_edge();
                }

                if !found {
                    unsafe {
                        (*fd).set_side(true);
                        (*fd).set_hiding(false);
                        (*fd).set_back(false);
                    }
                } else if unsafe { (*fd).closed() } {
                    // switch (fd.Orientation())
                    match unsafe { (*fd).orientation() } {
                        Orientation::Reversed => unsafe { (*fd).set_hiding((*fd).back()) },
                        Orientation::Forward => unsafe { (*fd).set_hiding(!(*fd).back()) },
                        Orientation::External => unsafe { (*fd).set_hiding(true) },
                        Orientation::Internal => unsafe { (*fd).set_hiding(false) },
                    }
                } else {
                    unsafe { (*fd).set_hiding(true) };
                }
            } else {
                if inverted {
                    unsafe { (*fd).set_hiding(false) };
                    unsafe { (*fd).set_back(true) };
                } else {
                    unsafe { (*fd).set_hiding(true) };
                    unsafe { (*fd).set_back(false) };
                }
            }

            let mut first_time = true;

            // for (myFaceItr1.InitEdge(fd); myFaceItr1.MoreEdge(); myFaceItr1.NextEdge())
            self.my_face_itr1.init_edge(fd);
            while self.my_face_itr1.more_edge() {
                // myFE = myFaceItr1.Edge();
                self.my_fe = self.my_face_itr1.edge() as usize;
                // HLRBRep_EdgeData& EDataFE2 = myEData(myFE);
                let edatafe2 = &mut self.my_e_data[(self.my_fe - 1) as usize];
                if !unsafe { (*fd).simple() } {
                    // EDataFE2.AutoIntersectionDone(false);
                    edatafe2.set_auto_intersection_done(false);
                }
                // HLRAlgo::DecodeMinMax(EDataFE2.MinMax(), EdgeMin, EdgeMax);
                HLRAlgo::decode_min_max(edatafe2.min_max(), &mut edge_min, &mut edge_max);
                if self.my_face_itr1.beginning_of_wire() {
                    HLRAlgo::copy_min_max(&edge_min, &edge_max, &mut wire_min, &mut wire_max);
                } else {
                    HLRAlgo::add_min_max(&edge_min, &edge_max, &mut wire_min, &mut wire_max);
                }
                if self.my_face_itr1.end_of_wire() {
                    HLRAlgo::encode_min_max(&wire_min, &wire_max, &mut min_max_wire);
                    // myFaceItr1.Wire()->UpdateMinMax(MinMaxWire);
                    self.my_face_itr1.wire().update_min_max(min_max_wire);
                    if first_time {
                        first_time = false;
                        HLRAlgo::copy_min_max(&wire_min, &wire_max, &mut face_min, &mut face_max);
                    } else {
                        HLRAlgo::add_min_max(&wire_min, &wire_max, &mut face_min, &mut face_max);
                    }
                }
                self.my_face_itr1.next_edge();
            }
            HLRAlgo::encode_min_max(&face_min, &face_max, &mut min_max_face);
            // fd.Wires()->UpdateMinMax(MinMaxFace);
            let fd_wires: *mut WiresBlock = unsafe { (*fd).wires() };
            unsafe { (*fd_wires).update_min_max(min_max_face) };
            // fd.Size(HLRAlgo::SizeBox(FaceMin, FaceMax));
            let size = HLRAlgo::size_box(&face_min, &face_max);
            unsafe { (*fd).set_size(size) };
            face += 1;
        }
    }

    /// OCCT HLRBRep_Data::InitBoundSort (cxx L984-1018) — to compare with
    /// only non rejected edges.
    pub fn init_bound_sort(&mut self, min_max_tot: &MinMaxIndices, e1: usize, e2: usize) {
        self.my_nbr_sort_ed = 0;
        // const HLRAlgo_EdgesBlock::MinMaxIndices& MinMaxShap = MinMaxTot;
        let min_max_shap = min_max_tot;

        // for (int e = e1; e <= e2; e++)
        let mut e = e1;
        while e <= e2 {
            // HLRBRep_EdgeData& ed = myEData(e);
            let ed: *mut EdgeData<'a> = &mut self.my_e_data[e - 1];
            if !unsafe { &*ed }.status_ref().all_hidden() {
                // myLEMinMax = &ed.MinMax();
                self.my_le_min_max = unsafe { (*ed).min_max() } as *mut MinMaxIndices;
                let le_min_max = unsafe { &*self.my_le_min_max };
                if ((min_max_shap.max[0].wrapping_sub(le_min_max.min[0])) & MASK_80008000) == 0
                    && ((le_min_max.max[0].wrapping_sub(min_max_shap.min[0])) & MASK_80008000) == 0
                    && ((min_max_shap.max[1].wrapping_sub(le_min_max.min[1])) & MASK_80008000) == 0
                    && ((le_min_max.max[1].wrapping_sub(min_max_shap.min[1])) & MASK_80008000) == 0
                    && ((min_max_shap.max[2].wrapping_sub(le_min_max.min[2])) & MASK_80008000) == 0
                    && ((le_min_max.max[2].wrapping_sub(min_max_shap.min[2])) & MASK_80008000) == 0
                    && ((min_max_shap.max[3].wrapping_sub(le_min_max.min[3])) & MASK_80008000) == 0
                    && ((le_min_max.max[3].wrapping_sub(min_max_shap.min[3])) & MASK_80008000) == 0
                    && ((min_max_shap.max[4].wrapping_sub(le_min_max.min[4])) & MASK_80008000) == 0
                    && ((le_min_max.max[4].wrapping_sub(min_max_shap.min[4])) & MASK_80008000) == 0
                    && ((min_max_shap.max[5].wrapping_sub(le_min_max.min[5])) & MASK_80008000) == 0
                    && ((le_min_max.max[5].wrapping_sub(min_max_shap.min[5])) & MASK_80008000) == 0
                    && ((min_max_shap.max[6].wrapping_sub(le_min_max.min[6])) & MASK_80008000) == 0
                    && ((le_min_max.max[6].wrapping_sub(min_max_shap.min[6])) & MASK_80008000) == 0
                    && ((min_max_shap.max[7].wrapping_sub(le_min_max.min[7])) & MASK_80008000) == 0
                    && ((le_min_max.max[7].wrapping_sub(min_max_shap.min[7])) & MASK_80008000) == 0
                {
                    //- rejection en z
                    self.my_nbr_sort_ed += 1;
                    // myEdgeIndices(myNbrSortEd) = e;
                    self.my_edge_indices[self.my_nbr_sort_ed - 1] = e as i32;
                }
            }
            e += 1;
        }
    }

    /// OCCT HLRBRep_Data::InitEdge (cxx L1022-1073) — begin an iteration
    /// only on visible Edges crossing the face number FI.
    ///
    /// rcad deviations (the K/M ruling + the ds_filler precedent): the MST
    /// map is the ptr_id-keyed Vec; the kernel context rides in
    /// [Data::my_brep] (stored by Update — the hiding phase runs after it)
    /// and is refreshed from the bound tool at the ChangeFind / Bind sites.
    pub fn init_edge(
        &mut self,
        fi: usize,
        mst: &mut Vec<(Shape, BRepTopAdaptorTool)>,
    ) {
        self.my_hide_count += 1;
        self.my_hide_count += 1;

        self.i_face = fi;
        // iFaceData = &myFData(iFace);
        self.i_face_data = &mut self.my_f_data[fi - 1];
        let fd = unsafe { &mut *self.i_face_data };
        // iFaceGeom = &iFaceData->Geometry();
        self.i_face_geom = fd.geometry() as *mut Surface<'a>;
        // iFaceBack = iFaceData->Back();
        self.i_face_back = fd.back();
        // iFaceSimp = iFaceData->Simple();
        self.i_face_simp = fd.simple();
        // iFaceMinMax = &iFaceData->Wires()->MinMax();
        let fd_wires: *mut WiresBlock = fd.wires();
        self.i_face_min_max = unsafe { (*fd_wires).min_max() } as *mut MinMaxIndices;
        // iFaceType = iFaceGeom->GetType();
        self.i_face_type = unsafe { (*self.i_face_geom).get_type() };
        // iFaceTest = !iFaceSimp;
        self.i_face_test = !self.i_face_simp;
        // mySLProps.SetSurface(iFaceGeom);
        self.my_s_l_props.set_surface(unsafe { &*self.i_face_geom });
        // myIntersector.Load(iFaceGeom);
        self.my_intersector.load(self.i_face_geom);

        // HLRBRep_Surface* p1 = iFaceGeom;
        let p1 = self.i_face_geom;
        // const BRepAdaptor_Surface& bras = p1->Surface();
        let bras = unsafe { &*p1 }.my_surface();
        // const TopoDS_Face& topodsface = bras.Face();
        let topodsface = bras.face().clone();

        if mst_is_bound(mst, &topodsface) {
            // BRepTopAdaptor_Tool& BRT = MST.ChangeFind(topodsface);
            let pos = mst_pos(mst, &topodsface);
            let brt = &mut mst[pos].1;
            // myClassifier = BRT.GetTopolTool();
            self.my_classifier = brt.get_topol_tool().clone();
            // rcad: the kernel context rides with the bound tool (the same
            // refcounted BRep allocation; the HLR object graph keeps it
            // alive for 'a — the raw-handle precedent).
            self.my_brep = Some(unsafe { &*Arc::as_ptr(brt.brep()) });
        } else {
            // BRepTopAdaptor_Tool BRT(topodsface, Precision::PConfusion());
            let tool_brep = Arc::new(
                self.my_brep
                    .expect("HLRBRep_Data::InitEdge: no kernel context")
                    .clone(),
            );
            let mut brt = BRepTopAdaptorTool::new_face(tool_brep, &topodsface, PCONFUSION);
            // MST.Bind(topodsface, BRT);
            mst_bind(mst, &topodsface, brt.clone());
            // myClassifier = BRT.GetTopolTool();
            self.my_classifier = brt.get_topol_tool().clone();
            // rcad: the kernel context of the freshly bound tool.
            self.my_brep = Some(unsafe { &*Arc::as_ptr(brt.brep()) });
        }

        if self.i_face_test {
            // iFaceSmpl = !iFaceData->Cut();
            self.i_face_smpl = !fd.cut();
            // myFaceItr2.InitEdge(*iFaceData);
            self.my_face_itr2.init_edge(self.i_face_data);
        } else {
            // for (myFaceItr1.InitEdge(*iFaceData); myFaceItr1.MoreEdge(); myFaceItr1.NextEdge())
            self.my_face_itr1.init_edge(self.i_face_data);
            while self.my_face_itr1.more_edge() {
                // myFE = myFaceItr1.Edge();  // edges of a simple hiding
                self.my_fe = self.my_face_itr1.edge() as usize;
                // myEData(myFE).HideCount(myHideCount - 1);  // face must be jumped.
                self.my_e_data[self.my_fe - 1].set_hide_count(self.my_hide_count - 1);
                self.my_face_itr1.next_edge();
            }
            self.my_cur_sort_ed = 1;
        }
        self.next_edge(false);
    }

    /// OCCT HLRBRep_Data::MoreEdge (cxx L1077-1110).
    pub fn more_edge(&mut self) -> bool {
        if self.i_face_test {
            if self.my_face_itr2.more_edge() {
                // all edges must be tested if the face is not a simple one.
                self.my_le = self.my_face_itr2.edge() as usize;
                self.my_le_out_line = self.my_face_itr2.out_line();
                self.my_le_internal = self.my_face_itr2.internal();
                self.my_le_double = self.my_face_itr2.double();
                self.my_le_iso_line = self.my_face_itr2.iso_line();
                // myLEData = &myEData(myLE);
                self.my_le_data = &mut self.my_e_data[self.my_le - 1];
                let le = unsafe { &mut *self.my_le_data };
                // myLEGeom = &myLEData->ChangeGeometry();
                self.my_le_geom = le.change_geometry() as *const Curve<'a>;
                // myLEMinMax = &myLEData->MinMax();
                self.my_le_min_max = le.min_max() as *mut MinMaxIndices;
                // myLETol = myLEData->Tolerance();
                self.my_le_tol = le.tolerance();
                // myLEType = myLEGeom->GetType();
                self.my_le_type = curve_type_to_geom_abs(le.geometry().get_type());
                if !self.my_le_double {
                    // myLEData->HideCount(myHideCount - 1);
                    le.set_hide_count(self.my_hide_count - 1);
                }
                return true;
            } else {
                // at the end of the test we know if it is a simple face
                self.i_face_test = false;
                self.i_face_simp = self.i_face_smpl;
                // iFaceData->Simple(iFaceSimp);
                unsafe { (*self.i_face_data).set_simple(self.i_face_simp) };
                self.my_cur_sort_ed = 1;
                self.next_edge(false);
            }
        }
        self.my_cur_sort_ed <= self.my_nbr_sort_ed
    }

    /// OCCT HLRBRep_Data::NextEdge (cxx L1114-1205).  The OCCT default
    /// argument (skip = true) is explicit in rcad.
    pub fn next_edge(&mut self, skip: bool) {
        if skip {
            if self.i_face_test {
                // myFaceItr2.NextEdge();
                self.my_face_itr2.next_edge();
            } else {
                self.my_cur_sort_ed += 1;
            }
        }
        if !self.more_edge() {
            return;
        }
        if self.i_face_test {
            self.my_le = self.my_face_itr2.edge() as usize;
            self.my_le_out_line = self.my_face_itr2.out_line();
            self.my_le_internal = self.my_face_itr2.internal();
            self.my_le_double = self.my_face_itr2.double();
            self.my_le_iso_line = self.my_face_itr2.iso_line();
            // myLEData = &myEData(myLE);
            self.my_le_data = &mut self.my_e_data[self.my_le - 1];
            let le = unsafe { &mut *self.my_le_data };
            self.my_le_geom = le.change_geometry() as *const Curve<'a>;
            self.my_le_min_max = le.min_max() as *mut MinMaxIndices;
            self.my_le_tol = le.tolerance();
            self.my_le_type = curve_type_to_geom_abs(le.geometry().get_type());
            if le.vertical()
                || (self.my_le_double && le.hide_count() == self.my_hide_count - 1)
            {
                self.next_edge(true);
            }
            le.set_hide_count(self.my_hide_count - 1);
            return;
        } else {
            // myLE = Edge();
            self.my_le = self.edge();
            self.my_le_out_line = false;
            self.my_le_internal = false;
            self.my_le_double = false;
            self.my_le_iso_line = false;
            // myLEData = &myEData(myLE);
            self.my_le_data = &mut self.my_e_data[self.my_le - 1];
            let le = unsafe { &mut *self.my_le_data };
            self.my_le_geom = le.change_geometry() as *const Curve<'a>;
            self.my_le_min_max = le.min_max() as *mut MinMaxIndices;
            self.my_le_tol = le.tolerance();
            self.my_le_type = curve_type_to_geom_abs(le.geometry().get_type());
        }
        let le = unsafe { &mut *self.my_le_data };
        if le.vertical() {
            self.next_edge(true);
            return;
        }
        if le.hide_count() > self.my_hide_count - 2 {
            self.next_edge(true);
            return;
        }
        if le.status_ref().all_hidden() {
            self.next_edge(true);
            return;
        }
        let le_min_max = unsafe { &*self.my_le_min_max };
        let i_face_min_max = unsafe { &*self.i_face_min_max };
        if ((i_face_min_max.max[0].wrapping_sub(le_min_max.min[0])) & MASK_80008000) != 0
            || ((le_min_max.max[0].wrapping_sub(i_face_min_max.min[0])) & MASK_80008000) != 0
            || ((i_face_min_max.max[1].wrapping_sub(le_min_max.min[1])) & MASK_80008000) != 0
            || ((le_min_max.max[1].wrapping_sub(i_face_min_max.min[1])) & MASK_80008000) != 0
            || ((i_face_min_max.max[2].wrapping_sub(le_min_max.min[2])) & MASK_80008000) != 0
            || ((le_min_max.max[2].wrapping_sub(i_face_min_max.min[2])) & MASK_80008000) != 0
            || ((i_face_min_max.max[3].wrapping_sub(le_min_max.min[3])) & MASK_80008000) != 0
            || ((le_min_max.max[3].wrapping_sub(i_face_min_max.min[3])) & MASK_80008000) != 0
            || ((i_face_min_max.max[4].wrapping_sub(le_min_max.min[4])) & MASK_80008000) != 0
            || ((le_min_max.max[4].wrapping_sub(i_face_min_max.min[4])) & MASK_80008000) != 0
            || ((i_face_min_max.max[5].wrapping_sub(le_min_max.min[5])) & MASK_80008000) != 0
            || ((le_min_max.max[5].wrapping_sub(i_face_min_max.min[5])) & MASK_80008000) != 0
            || ((i_face_min_max.max[6].wrapping_sub(le_min_max.min[6])) & MASK_80008000) != 0
            || ((le_min_max.max[6].wrapping_sub(i_face_min_max.min[6])) & MASK_80008000) != 0
            || ((i_face_min_max.max[7].wrapping_sub(le_min_max.min[7])) & MASK_80008000) != 0
            || ((le_min_max.max[7].wrapping_sub(i_face_min_max.min[7])) & MASK_80008000) != 0
        {
            //-- rejection en z
            self.next_edge(true);
            return;
        }
        // iFaceGeom->IsAbove(iFaceBack, myLEGeom, (double)myLETol)
        let le_geom = unsafe { &*self.my_le_geom };
        if unsafe { &*self.i_face_geom }.is_above(
            self.i_face_back,
            le_geom.view(),
            self.my_le_tol as f64,
        ) {
            self.next_edge(true);
            return;
        }
        // edge is OK
    }

    /// OCCT HLRBRep_Data::Edge (cxx L1209-1219) — returns the current Edge.
    pub fn edge(&self) -> usize {
        if self.i_face_test {
            // myFaceItr2.Edge()
            self.my_face_itr2.edge() as usize
        } else {
            // myEdgeIndices(myCurSortEd)
            self.my_edge_indices[self.my_cur_sort_ed - 1] as usize
        }
    }

    /// OCCT HLRBRep_Data::InitInterference (cxx L1223-1229).
    pub fn init_interference(&mut self) {
        // myLLProps.SetCurve(myLEGeom);
        self.my_l_l_props.set_curve(unsafe { &*self.my_le_geom });
        // myFaceItr1.InitEdge(*((HLRBRep_FaceData*)iFaceData));
        self.my_face_itr1.init_edge(self.i_face_data);
        // myNbPoints = myNbSegments = iInterf = 0;
        self.my_nb_segments = 0;
        self.my_nb_points = 0;
        self.i_interf = 0;
        self.next_interference();
    }

    // ------------------------------------------------------------------
    // OCCT HLRBRep_Data.lxx L19-120 — the inline accessors.
    // ------------------------------------------------------------------

    /// OCCT EDataArray (lxx L19-22).
    pub fn e_data_array(&self) -> &[EdgeData<'a>] {
        &self.my_e_data
    }

    /// OCCT FDataArray (lxx L26-29).
    pub fn f_data_array(&self) -> &[FaceData<'a>] {
        &self.my_f_data
    }

    /// OCCT EDataArray (lxx L19-22) — the mutable `NCollection_Array1&`
    /// overload (the ChangeValue write path of the InternalAlgo /
    /// HLRToShape layers: Selected / Status / Used / HideCount updates).
    pub fn e_data_array_mut(&mut self) -> &mut [EdgeData<'a>] {
        &mut self.my_e_data
    }

    /// OCCT FDataArray (lxx L26-29) — the mutable overload (the
    /// `FDataArray().ChangeValue(iface)` explorer of InitEdgeStatus and
    /// HLRToShape::DrawFace).
    pub fn f_data_array_mut(&mut self) -> &mut [FaceData<'a>] {
        &mut self.my_f_data
    }

    /// OCCT Projector (lxx L47-50) — the mutable `HLRAlgo_Projector&`
    /// overload (the `DS->Projector().Scaled(true)` of
    /// HLRBRep_HLRToShape::InternalCompound).
    pub fn projector_mut(&mut self) -> &mut Projector {
        &mut self.my_proj
    }

    /// rcad kernel-context reader — no OCCT counterpart (the [Data::my_brep]
    /// deviation accessor): the InternalAlgo::Update bridge needs the owning
    /// BRep to feed the BRep-carrying rcad form of `myDS->Update(myProj)`.
    pub fn brep(&self) -> Option<&'a rcad_kernel::topods::BRep> {
        self.my_brep
    }

    /// rcad kernel-context installer — no OCCT counterpart (the [Data::my_brep]
    /// deviation accessor): the ShapeToHLR::Load installs the per-Load leaked
    /// view of the session kernel context so the post-Update hiding phases
    /// read the loaded graph even before `myDS->Update(myProj)` stores it.
    pub fn set_brep(&mut self, brep: &'a rcad_kernel::topods::BRep) {
        self.my_brep = Some(brep);
    }

    /// OCCT Tolerance(const float tol) (lxx L33-36) — set the tolerance for
    /// the rejections during the exploration.
    pub fn set_tolerance(&mut self, tol: f32) {
        self.my_toler = tol;
    }

    /// OCCT Tolerance (lxx L40-43) — the tolerance for the rejections
    /// during the exploration.
    pub fn tolerance(&self) -> f32 {
        self.my_toler
    }

    /// OCCT Projector (lxx L47-50).
    pub fn projector(&self) -> &Projector {
        &self.my_proj
    }

    /// OCCT NbVertices (lxx L54-57).
    pub fn nb_vertices(&self) -> usize {
        self.my_nb_vertices
    }

    /// OCCT NbEdges (lxx L61-64).
    pub fn nb_edges(&self) -> usize {
        self.my_nb_edges
    }

    /// OCCT NbFaces (lxx L68-71).
    pub fn nb_faces(&self) -> usize {
        self.my_nb_faces
    }

    /// OCCT EdgeMap (lxx L75-78).
    pub fn edge_map(&self) -> &[Shape] {
        &self.my_e_map
    }

    /// OCCT EdgeMap (lxx L75-78) — the mutable `NCollection_IndexedMap&`
    /// form (the ShapeToHLR::Load `DS->EdgeMap().Add(Edg)`, cxx L168).
    pub fn edge_map_mut(&mut self) -> &mut Vec<Shape> {
        &mut self.my_e_map
    }

    /// OCCT FaceMap (lxx L82-85).
    pub fn face_map(&self) -> &[Shape] {
        &self.my_f_map
    }

    /// OCCT FaceMap (lxx L82-85) — the mutable `NCollection_IndexedMap&`
    /// form (the ShapeToHLR::ExploreFace `DS->FaceMap().Add(theFace)`,
    /// cxx L239).
    pub fn face_map_mut(&mut self) -> &mut Vec<Shape> {
        &mut self.my_f_map
    }

    /// The BigSize read-back (no OCCT accessor — myBigSize is private in
    /// OCCT and read through the impl blocks; the ShapeToHLR anchor tests
    /// assert it from another module).
    pub fn big_size(&self) -> f64 {
        self.my_big_size
    }

    /// The SurD read-back (same private-field read-back as big_size).
    pub fn sur_d(&self) -> &[f64; 16] {
        &self.my_sur_d
    }

    /// OCCT SimpleHidingFace (lxx L89-92) — true if the current hiding face
    /// is not an auto-intersected one.
    pub fn simple_hiding_face(&self) -> bool {
        self.i_face_simp
    }

    /// OCCT HidingTheFace (lxx L96-99) — true if the current edge to be
    /// hidden belongs to the hiding face.
    pub fn hiding_the_face(&self) -> bool {
        self.i_face_test
    }

    // MoreInterference (lxx L103-107) — landed in [classify] next to the
    // interference exploration group.

    /// OCCT Interference (lxx L110-113).
    pub fn interference(&self) -> &Interference {
        &self.my_intf
    }

    /// OCCT EdgeOfTheHidingFace (lxx L117-120) — true if the Edge E belongs
    /// to the Hiding Face.
    pub fn edge_of_the_hiding_face(&self, _e: usize, ed: &EdgeData<'_>) -> bool {
        ed.hide_count() == self.my_hide_count - 1
    }
}

impl Drop for Data<'_> {
    /// OCCT ~HLRBRep_Data() { Destroy(); } (hxx L202).
    fn drop(&mut self) {
        if !self.my_reject.is_null() {
            self.destroy();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::brep::b_curve_tool::CurveView;
    use glam::DVec3;
    use rcad_kernel::geom::{Line3, Plane, Point3, Surface3, Vec3};
    use rcad_kernel::math::gp::Ax2;
    use rcad_kernel::topo::topods::BRepBuilder;

    /// A straight 3D segment adaptor (the BRepAdaptor_Curve test double).
    struct SegView {
        origin: Point3,
        dir: Vec3,
        len: f64,
    }

    impl CurveView for SegView {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.len
        }
        fn d0(&self, u: f64) -> Point3 {
            self.origin + self.dir * u
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (self.origin + self.dir * u, self.dir)
        }
        fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
            (self.origin + self.dir * u, self.dir, Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> Line3 {
            Line3 {
                origin: self.origin,
                direction: self.dir,
            }
        }
        fn circle(&self) -> rcad_kernel::geom::Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
            panic!("Standard_NoSuchObject");
        }
        fn degree(&self) -> i32 {
            0
        }
        fn nb_poles(&self) -> i32 {
            0
        }
        fn nb_knots(&self) -> i32 {
            0
        }
        fn is_closed(&self) -> bool {
            false
        }
        fn is_periodic(&self) -> bool {
            false
        }
        fn period(&self) -> f64 {
            0.0
        }
        fn resolution(&self, r3d: f64) -> f64 {
            r3d
        }
        fn parameter_3d(&self, p2d: f64) -> f64 {
            p2d
        }
        fn poles(&self) -> Vec<Point3> {
            Vec::new()
        }
    }

    /// The top-view projector (identity transform, type 1).
    fn top_view_projector() -> Projector {
        Projector::from_ax2(&Ax2::new(
            DVec3::ZERO,
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 0.0),
        ))
    }

    /// The OCCT Data Set path: the loaded projected curve inside a fresh
    /// EdgeData with the full-parameter domain.
    fn line_edge_data(
        proj: *const Projector,
        seg: &'static SegView,
        v1: i32,
        v2: i32,
    ) -> EdgeData<'static> {
        let mut c = Curve::new();
        c.projector(proj);
        c.load(seg);
        let mut e = EdgeData::new();
        e.set(
            false, false, c, 1.0e-7, v1, v2, false, false, false, false, 0.0, 1.0e-7, 10.0,
            1.0e-7,
        );
        e
    }

    /// A square face in the x = 0 plane ("seen edge-on" under the top
    /// view): four edges 1..4 around the corners.
    fn side_plane_brep() -> (&'static rcad_kernel::BRep, Shape) {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let p = [
            DVec3::new(0.0, 0.0, -5.0),
            DVec3::new(0.0, 10.0, -5.0),
            DVec3::new(0.0, 10.0, 5.0),
            DVec3::new(0.0, 0.0, 5.0),
        ];
        let dirs = [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, -1.0),
        ];
        let verts = [
            b.add_vertex(&mut brep, p[0], 1e-7),
            b.add_vertex(&mut brep, p[1], 1e-7),
            b.add_vertex(&mut brep, p[2], 1e-7),
            b.add_vertex(&mut brep, p[3], 1e-7),
        ];
        let mut edges = Vec::new();
        for i in 0..4 {
            let e = b.add_edge(
                &mut brep,
                Some(rcad_kernel::geom::Curve3::Line(Line3 {
                    origin: p[i],
                    direction: dirs[i],
                })),
                verts[i].clone(),
                verts[(i + 1) % 4].clone(),
                [0.0, 10.0],
            );
            edges.push(e);
        }
        let wire = brep.add_twire(edges);
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane::new(
                DVec3::ZERO,
                DVec3::new(1.0, 0.0, 0.0),
            ))),
            wire,
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );
        let brep: &'static rcad_kernel::BRep = Box::leak(Box::new(brep));
        (brep, face)
    }

    /// An EdgeData with an explicit encoded box (the bound-sort anchors).
    fn mk_edge_mm(min: [i32; 8], max: [i32; 8]) -> EdgeData<'static> {
        let mut e = EdgeData::new();
        let mut mm = MinMaxIndices::default();
        mm.min = min;
        mm.max = max;
        e.update_min_max(&mm);
        e.status().show_all();
        e
    }

    /// The ctor (cxx L522-537): the counts, the array lengths, the default
    /// tolerance, the EdgeData Selected(true) defaults and the
    /// TableauRejection SetDim.
    #[test]
    fn data_ctor_layout() {
        let mut d = Data::new(3, 5, 7);
        assert_eq!(d.nb_vertices(), 3);
        assert_eq!(d.nb_edges(), 5);
        assert_eq!(d.nb_faces(), 7);
        // myEData(0, NE) / myFData(0, NF) — the rcad Vecs keep the live slots.
        assert_eq!(d.e_data_array().len(), 5);
        assert_eq!(d.f_data_array().len(), 7);
        assert_eq!(d.edge_map().len(), 0);
        assert_eq!(d.face_map().len(), 0);
        // myToler((float)1e-5).
        assert_eq!(d.tolerance(), 1.0e-5);
        // the EdgeData default ctor: Selected(true).
        assert!(d.e_data_array()[0].selected());
        // Tolerance(tol) / Tolerance() round trip (lxx L33-43).
        d.set_tolerance(0.5);
        assert_eq!(d.tolerance(), 0.5);
    }

    /// InitBoundSort (cxx L984-1018): the accepted edge lands in
    /// myEdgeIndices in increasing order; the z-rejected and the
    /// all-hidden edges are excluded.
    #[test]
    fn init_bound_sort_order_and_z_rejection() {
        let mut d = Data::new(0, 3, 0);
        d.my_e_data[0] = mk_edge_mm([0; 8], [100; 8]); // inside the shape box
        d.my_e_data[1] = mk_edge_mm([5000; 8], [6000; 8]); // outside -> rejected en z
        d.my_e_data[2] = {
            let mut e = mk_edge_mm([0; 8], [100; 8]);
            e.status().hide_all(); // excluded by the AllHidden test
            e
        };
        let tot = MinMaxIndices {
            min: [0; 8],
            max: [200; 8],
        };
        d.init_bound_sort(&tot, 1, 3);
        assert_eq!(d.my_nbr_sort_ed, 1);
        assert_eq!(d.my_edge_indices[0], 1);
    }

    /// MoreEdge / NextEdge (cxx L1077-1205) on the sorted path: the edges
    /// come in the bound-sorted order and the walk exhausts.
    #[test]
    fn next_edge_walk_order_on_sorted_path() {
        let mut d = Data::new(0, 2, 0);
        // the loaded projected curves (the IsAbove argument reads the
        // CurveView of myLEGeom).
        let proj: &'static Projector = Box::leak(Box::new(top_view_projector()));
        let seg: &'static SegView = Box::leak(Box::new(SegView {
            origin: DVec3::ZERO,
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 10.0,
        }));
        d.my_e_data[0] = line_edge_data(proj, seg, 1, 2);
        d.my_e_data[1] = line_edge_data(proj, seg, 2, 3);
        let tot = MinMaxIndices {
            min: [0; 8],
            max: [100; 8],
        };
        d.init_bound_sort(&tot, 1, 2);
        assert_eq!(d.my_nbr_sort_ed, 2);
        // the InitEdge state on the non-test path.
        d.i_face_test = false;
        d.my_hide_count = 2;
        d.i_face_min_max = Box::leak(Box::new(tot));
        d.i_face_geom = Box::leak(Box::new(Surface::new()));
        d.my_cur_sort_ed = 1;
        // NextEdge(false): the first sorted edge is OK (a non-planar hiding
        // surface is never "above", cxx L1199 -> the edge survives).
        d.next_edge(false);
        assert!(d.more_edge());
        assert_eq!(d.my_le, 1);
        assert!(!d.my_le_out_line);
        assert!(!d.my_le_internal);
        assert!(!d.my_le_double);
        assert!(!d.my_le_iso_line);
        assert_eq!(d.edge(), 1);
        d.next_edge(true);
        assert!(d.more_edge());
        assert_eq!(d.my_le, 2);
        assert_eq!(d.edge(), 2);
        d.next_edge(true);
        assert!(!d.more_edge());
    }

    /// NextEdge (cxx L1165-1168): the vertical candidates reject
    /// recursively and the walk exhausts without the HideCount stamp (the
    /// cxx L1149 write lives on the iFaceTest path only).
    #[test]
    fn next_edge_skips_vertical_edges() {
        let mut d = Data::new(0, 1, 0);
        let mut e = mk_edge_mm([0; 8], [100; 8]);
        e.set_vertical(true);
        d.my_e_data[0] = e;
        let tot = MinMaxIndices {
            min: [0; 8],
            max: [100; 8],
        };
        d.init_bound_sort(&tot, 1, 1);
        assert_eq!(d.my_nbr_sort_ed, 1);
        d.i_face_test = false;
        d.my_hide_count = 2;
        d.i_face_min_max = Box::leak(Box::new(tot));
        d.i_face_geom = Box::leak(Box::new(Surface::new()));
        d.my_cur_sort_ed = 1;
        d.next_edge(false);
        assert!(!d.more_edge());
        assert_eq!(d.my_e_data[0].hide_count(), 0);
    }

    /// Write (cxx L548-608): the edges land at 1-based de+1.. with the dv
    /// vertex shift and the map copies; the face wires shift by de.
    #[test]
    fn write_translates_sub_data() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let va = b.add_vertex(&mut brep, DVec3::ZERO, 1e-7);
        let vb = b.add_vertex(&mut brep, DVec3::new(1.0, 0.0, 0.0), 1e-7);

        // DS: two edges (VSta 1, 2), one face with one wire of edges 1, 2.
        let mut ds = Data::new(0, 2, 1);
        ds.my_e_data[0].set_v_sta(1);
        ds.my_e_data[1].set_v_sta(2);
        ds.my_e_map.push(va.clone());
        ds.my_e_map.push(vb.clone());
        let mut wb = WiresBlock::new(1);
        let mut eb = crate::hlr::algo::edges_block::EdgesBlock::new(2);
        eb.set_edge(1, 1);
        eb.set_edge(2, 2);
        wb.set(1, eb);
        *ds.my_f_data[0].change_wires() = Some(Arc::new(wb));
        ds.my_f_map.push(va.clone());

        // the target has room for the merge at de = 2, df = 1.
        let mut d = Data::new(0, 4, 2);
        d.write(&mut ds, 10, 2, 1);
        // the edges land in 1-based slots 3, 4 with VSta shifted by dv = 10.
        assert_eq!(d.my_e_data[2].v_sta(), 11);
        assert_eq!(d.my_e_data[3].v_sta(), 12);
        assert_eq!(d.my_e_data[0].v_sta(), 0);
        // myEMap.Add(DS->EdgeMap().FindKey(iedge)).
        assert_eq!(d.edge_map().len(), 2);
        assert_eq!(d.edge_map()[0].ptr_id(), va.ptr_id());
        assert_eq!(d.edge_map()[1].ptr_id(), vb.ptr_id());
        // the face lands in 1-based slot 2; the wire edges shifted by de = 2.
        assert_eq!(d.face_map().len(), 1);
        let wb = d.my_f_data[1].wires();
        assert_eq!(unsafe { (*wb).wire(1).edge(1) }, 3);
        assert_eq!(unsafe { (*wb).wire(1).edge(2) }, 4);
    }

    /// Update (cxx L612-980) over a side plane under the top view: the
    /// analytic myBigSize / mySurD / myDeca (the OCCT formulas over the 14
    /// trig directions), the per-edge flags and the face Side/Hiding/Back
    /// states of the "seen edge-on" plane.
    #[test]
    fn update_side_plane_face_anchor() {
        let (brep, face) = side_plane_brep();
        let proj: &'static Projector = Box::leak(Box::new(top_view_projector()));
        let mut data = Data::new(4, 4, 1);
        // the projected-curve adaptors of the four square edges.
        let corners = [
            DVec3::new(0.0, 0.0, -5.0),
            DVec3::new(0.0, 10.0, -5.0),
            DVec3::new(0.0, 10.0, 5.0),
            DVec3::new(0.0, 0.0, 5.0),
        ];
        let dirs = [
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, -1.0),
        ];
        for i in 0..4 {
            let seg: &'static SegView = Box::leak(Box::new(SegView {
                origin: corners[i],
                dir: dirs[i],
                len: 10.0,
            }));
            data.my_e_data[i] = line_edge_data(proj, seg, i as i32 + 1, ((i + 1) % 4) as i32 + 1);
        }
        // the hiding face (Set + SetWire + SetWEdge of HLRBRep_FaceData).
        let mut fd = FaceData::new();
        fd.set(brep, &face, Orientation::Forward, true, 1);
        fd.set_wire(1, 4);
        for i in 0..4 {
            fd.set_w_edge(1, i + 1, i as i32 + 1, Orientation::Forward, false, false, false, false);
        }
        data.my_f_data[0] = fd;

        data.update(brep, proj);

        // ---- the 16-direction spans of the top view: x = 0, y in [0, 10],
        // z in [-5, 5] -> d[2k] in [0, 10 sin], d[2k+1] in [-10 cos, 0],
        // d[14] = d[15] = the z span. ----
        let mut span = [0.0f64; 16];
        let mut big = 10.0f64; // the z span
        for k in 0..7usize {
            let a = (k as f64) * std::f64::consts::PI / 14.0;
            let (c, s) = (a.cos(), a.sin());
            span[2 * k] = 10.0 * s;
            span[2 * k + 1] = 10.0 * c;
            big = big.max(span[2 * k]).max(span[2 * k + 1]);
        }
        span[14] = 10.0;
        span[15] = 10.0;
        assert!(
            (data.my_big_size - big).abs() <= 1e-9 * big,
            "bigSize {} vs {}",
            data.my_big_size,
            big
        );
        // mySurD[i] = 0x00007fff / (d[i] + precad), precad = bigSize*0.0005.
        let p1 = big * 0.0005;
        for i in 0..16 {
            let expected = 32767.0 / (span[i] + p1);
            assert!(
                (data.my_sur_d[i] - expected).abs() <= 1e-6 * expected.abs().max(1.0),
                "surD[{}]",
                i
            );
        }
        // myDeca[i] = -TotMin[i] + precad*0.5.
        let p2 = big * 0.00025;
        assert!((data.my_deca[0] - p2).abs() <= 1e-12);
        assert!((data.my_deca[1] - (10.0 + p2)).abs() <= 1e-9);
        assert!(
            ((data.my_deca[3] - (10.0 * (std::f64::consts::PI / 14.0).cos() + p2)).abs()) <= 1e-9
        );
        assert!((data.my_deca[14] - (5.0 + p2)).abs() <= 1e-9);
        assert!((data.my_deca[15] - (5.0 + p2)).abs() <= 1e-9);

        // the per-edge updates: the y-direction edges have in-plane
        // tangents; the z-direction edges are Vertical (their projected
        // image degenerates to a point) and stamp the shared vertices.
        for (i, e) in data.e_data_array().iter().enumerate() {
            assert!(e.simple());
            assert!(e.auto_intersection_done());
            if i % 2 == 0 {
                assert!(!e.vertical());
                assert!(!e.ver_at_sta());
                assert!(!e.ver_at_end());
            } else {
                assert!(e.vertical());
                assert!(e.ver_at_sta());
                assert!(e.ver_at_end());
            }
        }
        // the face updates: a side plane hides nothing.
        let fd2 = &data.my_f_data[0];
        assert!(fd2.side());
        assert!(!fd2.hiding());
        assert!(!fd2.back());
        assert!(fd2.simple());
        assert!(fd2.plane());
        assert!(!fd2.cut());
        assert!(!fd2.with_out_l());
        // the x-degenerate box: the OCCT SizeBox product carries the zero
        // x term, so the face size is exactly 0.
        assert_eq!(fd2.size(), 0.0);
    }

    /// InitEdge (cxx L1022-1073): the double hide-count increment, the
    /// iFace* state, the MST bind and the cache hit on the second call.
    #[test]
    fn init_edge_mst_cache_and_classifier() {
        let (brep, face) = side_plane_brep();
        let mut data = Data::new(0, 1, 1);
        data.my_brep = Some(brep);
        let mut fd = FaceData::new();
        fd.set(brep, &face, Orientation::Forward, false, 1);
        fd.set_wire(1, 1);
        fd.set_w_edge(1, 1, 1, Orientation::Forward, false, false, false, false);
        data.my_f_data[0] = fd;
        data.my_e_data[0] = EdgeData::new();
        let mut mst: Vec<(Shape, BRepTopAdaptorTool)> = Vec::new();

        data.init_edge(1, &mut mst);
        // the OCCT double HideCount increment.
        assert_eq!(data.my_hide_count, 2);
        assert_eq!(data.i_face, 1);
        // the face was bound in the MST map (the miss branch).
        assert_eq!(mst.len(), 1);
        // the kernel context rides with the bound tool.
        assert!(data.my_brep.is_some());
        // the non-simple face keeps the iFaceTest exploration over the wire.
        assert!(data.hiding_the_face());
        assert!(!data.simple_hiding_face());
        // InitEdge ends with NextEdge(false): the wire edge is the current
        // candidate and got the HideCount stamp of MoreEdge.
        assert!(data.more_edge());
        assert_eq!(data.my_le, 1);
        assert_eq!(data.my_e_data[0].hide_count(), 1);

        // the second InitEdge hits the MST cache (no new binding).
        data.init_edge(1, &mut mst);
        assert_eq!(data.my_hide_count, 4);
        assert_eq!(mst.len(), 1);
        assert!(data.more_edge());
        assert_eq!(data.my_le, 1);
    }

    /// The lxx bit semantics (lxx L89-120): HidingTheFace /
    /// SimpleHidingFace / EdgeOfTheHidingFace / MoreInterference /
    /// Interference.
    #[test]
    fn lxx_hiding_face_bit_semantics() {
        let mut d = Data::new(0, 1, 0);
        // EdgeOfTheHidingFace: the HideCount stamp of InitEdge.
        d.my_hide_count = 3;
        let ed = EdgeData::new();
        assert!(!d.edge_of_the_hiding_face(1, &ed));
        let mut ed2 = EdgeData::new();
        ed2.set_hide_count(2); // myHideCount - 1
        assert!(d.edge_of_the_hiding_face(1, &ed2));
        // HidingTheFace / SimpleHidingFace.
        d.i_face_test = true;
        assert!(d.hiding_the_face());
        d.i_face_simp = true;
        assert!(d.simple_hiding_face());
        // MoreInterference (lxx L103-107): iInterf <= NbPoints + 2*NbSegments.
        d.my_nb_points = 1;
        d.my_nb_segments = 2;
        d.i_interf = 5;
        assert!(d.more_interference());
        d.i_interf = 6;
        assert!(!d.more_interference());
        // Interference (lxx L110-113).
        assert_eq!(d.interference().orientation(), Orientation::Forward);
    }
}       
