// OCCT HLRBRep_InternalAlgo (TKHLR/HLRBRep/HLRBRep_InternalAlgo.hxx L1-139
// + .cxx L1-1020) — the internal HLR framework: the projector, the sequence
// of ShapeBounds and the hiding driver over the HLRBRep_Data structure.
//
// Contract消费 notes (the parallel-agent landing points this file relies on):
// - [`super::shape_bounds::ShapeBounds`] (proxy O1): the two OCCT ctors,
//   shape / shape_data / nb_of_iso / translate / sizes / bounds /
//   update_min_max / min_max accessors;
// - [`super::shape_to_hlr::load`] (proxy O1):
//   `load(&Arc<OutLiner>, &Projector, &mut MST, nb_iso) -> Box<Data<'static>>`;
// - [`super::hider::Hider`] (proxy N): `new(&mut Data)`, `own_hiding(fi)`,
//   `hide(fi, &mut MST)`;
// - [`super::algo::Algo`] (proxy N): `Algo { pub internal: InternalAlgo }`.
//
// Deviations (all reported to the coordinator):
// - `occ::handle<HLRBRep_Data>` keeps unique ownership
//   (`Option<Box<Data<'static>>>`) — the OCCT handle aliasing of the copy
//   constructor is unrepresentable (see [`InternalAlgo::new_from`]);
// - `occ::handle<Standard_Transient>` maps to [`SDataHandle`];
// - the OCCT overloads Select(I) / ShowAll(I) / HideAll(I) / Hide() /
//   Hide(I) / Hide(I, J) / Load(...) map to distinct names (Rust has no
//   overload resolution);
// - `myDS->Update(myProj)` feeds the rcad kernel-context BRep documented on
//   [Data::update] — scanned from the loaded Data (see [`InternalAlgo::update`]).

use std::panic::AssertUnwindSafe;
use std::sync::Arc;

use rcad_kernel::topods::Shape;

use crate::hlr::algo::edges_block::MinMaxIndices;
use crate::hlr::algo::projector::Projector;
use crate::hlr::brep::data::Data;
use crate::hlr::brep::edge_data::EdgeData;
use crate::hlr::brep::face_data::FaceData;
use crate::hlr::brep::face_iterator::FaceIterator;
use crate::hlr::brep::hider::Hider;
use crate::hlr::brep::shape_bounds::ShapeBounds;
use crate::hlr::brep::shape_to_hlr;
use crate::topalgo::brep_top_adaptor::tool::BRepTopAdaptorTool;

// the SDataHandle re-export (the sibling Algo / ShapeBounds consumers import
// it from this module — the OCCT Standard_Transient handle of the package;
// the definition lives in shape_bounds.rs).
pub use crate::hlr::brep::shape_bounds::SDataHandle;

/// OCCT `static int HLRBRep_InternalAlgo_TRACE = true;` (cxx L45).
const INTERNAL_ALGO_TRACE: bool = true;
/// OCCT `static int HLRBRep_InternalAlgo_TRACE10 = true;` (cxx L46).
const INTERNAL_ALGO_TRACE10: bool = true;

/// OCCT `& 0x80008000` (cxx L628) — the int bit test of the encoded box
/// words; the literal needs the u32->i32 form (overflowing_literals).
const MASK_80008000: i32 = 0x80008000u32 as i32;

impl Default for InternalAlgo {
    fn default() -> Self {
        InternalAlgo::new()
    }
}

/// OCCT HLRBRep_InternalAlgo (hxx L36-137).
pub struct InternalAlgo {
    /// OCCT myDS (hxx L132) — the unique-ownership handle form.
    my_ds: Option<Box<Data<'static>>>,
    /// OCCT myProj (hxx L133).
    my_proj: Projector,
    /// OCCT myShapes (hxx L134) — NCollection_Sequence<HLRBRep_ShapeBounds>.
    my_shapes: Vec<ShapeBounds>,
    /// OCCT myMapOfShapeTool (hxx L135) — NCollection_DataMap as the
    /// ptr_id-keyed Vec (the ds_filler MST precedent).
    my_map_of_shape_tool: Vec<(Shape, BRepTopAdaptorTool)>,
    /// OCCT myDebug (hxx L136).
    my_debug: bool,
}

impl InternalAlgo {
    /// OCCT HLRBRep_InternalAlgo() (cxx L50-53) — myDebug(false).
    pub fn new() -> Self {
        InternalAlgo {
            my_ds: None,
            my_proj: Projector::new(),
            my_shapes: Vec::new(),
            my_map_of_shape_tool: Vec::new(),
            my_debug: false,
        }
    }

    /// OCCT HLRBRep_InternalAlgo(const handle<HLRBRep_InternalAlgo>& A)
    /// (cxx L57-63).
    ///
    /// rcad deviation: the OCCT myDS handle is shared with A; the rcad
    /// unique-ownership Box cannot alias, so the copy starts with a null DS
    /// (the next Update rebuilds it).  myProj / myShapes / myDebug copy.
    pub fn new_from(a: &InternalAlgo) -> Self {
        InternalAlgo {
            my_ds: None, // myDS = A->DataStructure() — see the deviation note
            my_proj: a.projector().clone(),
            my_shapes: a.my_shapes.clone(),
            my_map_of_shape_tool: a.my_map_of_shape_tool.clone(),
            my_debug: a.my_debug,
        }
    }

    /// OCCT Projector(const HLRAlgo_Projector& P) (cxx L67-70) — set the
    /// projector (the setter overload; Rust has no overloads).
    pub fn set_projector(&mut self, p: &Projector) {
        // myProj = P;
        self.my_proj = p.clone();
    }

    /// OCCT HLRAlgo_Projector& Projector() (cxx L74-77) — the reference
    /// getter split into the const / mut pair (Rust has no overloads).
    pub fn projector(&self) -> &Projector {
        &self.my_proj
    }

    /// OCCT HLRAlgo_Projector& Projector() — the mutable form.
    pub fn projector_mut(&mut self) -> &mut Projector {
        &mut self.my_proj
    }

    /// OCCT DataStructure() (cxx L1017-1020) — the handle getter split into
    /// the const / mut pair.
    pub fn data_structure(&self) -> Option<&Data<'static>> {
        self.my_ds.as_deref()
    }

    /// OCCT DataStructure() — the mutable form (the non-const handle deref).
    pub fn data_structure_mut(&mut self) -> Option<&mut Data<'static>> {
        self.my_ds.as_deref_mut()
    }

    /// OCCT Update() (cxx L81-185) — update the DataStructure.
    pub fn update(&mut self) {
        // if (!myShapes.IsEmpty())
        if !self.my_shapes.is_empty() {
            // int n = myShapes.Length();
            let n = self.my_shapes.len();
            // occ::handle<HLRBRep_Data>* DS = new occ::handle<HLRBRep_Data>[n];
            let mut ds: Vec<Box<Data<'static>>> = Vec::with_capacity(n);

            // int i, dv, de, df, nv = 0, ne = 0, nf = 0;
            let mut dv: i32;
            let mut de: i32;
            let mut df: i32;
            let mut nv: i32 = 0;
            let mut ne: i32 = 0;
            let mut nf: i32 = 0;

            // for (i = 1; i <= n; i++)
            for i in 1..=n {
                // HLRBRep_ShapeBounds& SB = myShapes(i);
                // (the OCCT reference maps to the cloned ctor inputs; the
                // reassignment writes back below.)
                let (shape, s_data, nb_iso) = {
                    let sb = &self.my_shapes[i - 1];
                    (sb.shape().clone(), sb.shape_data().clone(), sb.nb_of_iso())
                };
                // the &mut OutLiner of the landed rcad Load form — the raw
                // pointer cast keeps the OCCT handle aliasing (the
                // Algo::Index raw-pointer precedent).
                let out_liner: *mut crate::hlr::topo_brep::out_liner::OutLiner =
                    Arc::as_ptr(&shape) as *const crate::hlr::topo_brep::out_liner::OutLiner
                        as *mut crate::hlr::topo_brep::out_liner::OutLiner;
                // try { OCC_CATCH_SIGNALS
                //   DS[i - 1] = HLRBRep_ShapeToHLR::Load(SB.Shape(), myProj,
                //                                        myMapOfShapeTool, SB.NbOfIso());
                //   dv = DS[i - 1]->NbVertices(); de = ...; df = ...;
                // }
                let loaded = std::panic::catch_unwind(AssertUnwindSafe(|| {
                    // the per-shape builder arena of the landed rcad Load form
                    let mut arena = rcad_kernel::BRep::new();
                    shape_to_hlr::load(
                        &mut arena,
                        unsafe { &mut *out_liner },
                        &self.my_proj,
                        &mut self.my_map_of_shape_tool,
                        nb_iso.max(0) as usize,
                    )
                }));
                let ds_i: Box<Data<'static>>;
                match loaded {
                    Ok(mut d) => {
                        dv = d.nb_vertices() as i32;
                        de = d.nb_edges() as i32;
                        df = d.nb_faces() as i32;
                        ds_i = d;
                    }
                    // catch (Standard_Failure const& anException)
                    Err(an_exception) => {
                        if self.my_debug {
                            println!("An exception was caught when preparing the Shape {}", i);
                            println!(" and computing its OutLines ");
                            let _ = an_exception; // the OCCT `std::cout << anException`
                        }
                        // DS[i - 1] = new HLRBRep_Data(0, 0, 0);
                        ds_i = Box::new(Data::new(0, 0, 0));
                        dv = 0;
                        de = 0;
                        df = 0;
                    }
                }
                ds.push(ds_i);

                // SB = HLRBRep_ShapeBounds(SB.Shape(), SB.ShapeData(), SB.NbOfIso(),
                //                          1, dv, 1, de, 1, df);
                self.my_shapes[i - 1] =
                    ShapeBounds::new_with_data(shape, s_data, nb_iso, 1, dv, 1, de, 1, df);
                // nv += dv; ne += de; nf += df;
                nv += dv;
                ne += de;
                nf += df;
            }

            // myDS->Update(myProj);
            // (the rcad [Data::update] carries the owning BRep — the
            // kernel-context deviation documented there; it rides in my_brep
            // of the loaded Data.  The OCCT global TShape-graph stand-in
            // scans the source array before the merge / drop.)
            let bridge = ds.iter().find_map(|d| d.brep());

            // if (n == 1) { myDS = DS[0]; }
            if n == 1 {
                self.my_ds = Some(ds.into_iter().next().expect("DS[0]"));
            } else {
                // myDS = new HLRBRep_Data(nv, ne, nf);
                let mut my_ds = Box::new(Data::new(nv as usize, ne as usize, nf as usize));
                // nv = 0; ne = 0; nf = 0;
                nv = 0;
                ne = 0;
                nf = 0;

                // for (i = 1; i <= n; i++)
                for i in 1..=n {
                    // HLRBRep_ShapeBounds& SB = myShapes(i);
                    // SB.Sizes(dv, de, df);
                    (dv, de, df) = self.my_shapes[i - 1].sizes();
                    // SB.Translate(nv, ne, nf);
                    self.my_shapes[i - 1].translate(nv, ne, nf);
                    // myDS->Write(DS[i - 1], nv, ne, nf);
                    my_ds.write(ds[i - 1].as_mut(), nv, ne, nf);
                    // nv += dv; ne += de; nf += df;
                    nv += dv;
                    ne += de;
                    nf += df;
                }
                self.my_ds = Some(my_ds);
            }

            // delete[] DS; — the Vec drop (the sources die here, as in OCCT).

            let my_ds = self.my_ds.as_mut().expect("HLRBRep_InternalAlgo::Update : null DS");
            let brep = match my_ds.brep() {
                Some(b) => b,
                None => bridge.expect("HLRBRep_InternalAlgo::Update : no kernel context"),
            };
            my_ds.update(brep, &self.my_proj);

            // HLRAlgo_EdgesBlock::MinMaxIndices ShapMin, ShapMax, MinMaxShap;
            // HLRAlgo_EdgesBlock::MinMaxIndices TheMin, TheMax;
            let mut shap_min = MinMaxIndices::default();
            let mut shap_max = MinMaxIndices::default();
            let mut min_max_shap = MinMaxIndices::default();
            let mut the_min = MinMaxIndices::default();
            let mut the_max = MinMaxIndices::default();
            // (the two OCCT array references borrow at each use — the rcad
            // single-borrow form.)
            let ds = self.my_ds.as_mut().expect("HLRBRep_InternalAlgo::Update : null DS");

            // for (i = 1; i <= n; i++)
            for i in 1..=n {
                // bool FirstTime = true;
                let mut first_time = true;
                // HLRBRep_ShapeBounds& SB = myShapes(i);
                // int v1, v2, e1, e2, f1, f2;
                // SB.Bounds(v1, v2, e1, e2, f1, f2);
                let (_v1, _v2, e1, e2, f1, f2) = self.my_shapes[i - 1].bounds();

                // for (int e = e1; e <= e2; e++)
                let mut e = e1;
                while e <= e2 {
                    // HLRBRep_EdgeData& ed = aEDataArray.ChangeValue(e);
                    let ed = &mut ds.e_data_array_mut()[(e - 1) as usize];
                    // HLRAlgo::DecodeMinMax(ed.MinMax(), TheMin, TheMax);
                    crate::hlr::algo::hlr_algo::HLRAlgo::decode_min_max(
                        ed.min_max(),
                        &mut the_min,
                        &mut the_max,
                    );
                    if first_time {
                        first_time = false;
                        // HLRAlgo::CopyMinMax(TheMin, TheMax, ShapMin, ShapMax);
                        crate::hlr::algo::hlr_algo::HLRAlgo::copy_min_max(
                            &the_min,
                            &the_max,
                            &mut shap_min,
                            &mut shap_max,
                        );
                    } else {
                        // HLRAlgo::AddMinMax(TheMin, TheMax, ShapMin, ShapMax);
                        crate::hlr::algo::hlr_algo::HLRAlgo::add_min_max(
                            &the_min,
                            &the_max,
                            &mut shap_min,
                            &mut shap_max,
                        );
                    }
                    e += 1;
                }

                // for (int f = f1; f <= f2; f++)
                let mut f = f1;
                while f <= f2 {
                    // HLRBRep_FaceData& fd = aFDataArray.ChangeValue(f);
                    let fd = &mut ds.f_data_array_mut()[(f - 1) as usize];
                    // HLRAlgo::DecodeMinMax(fd.Wires()->MinMax(), TheMin, TheMax);
                    let wb = fd.wires();
                    let mm = unsafe { (*wb).min_max() } as *mut MinMaxIndices;
                    crate::hlr::algo::hlr_algo::HLRAlgo::decode_min_max(
                        unsafe { &*mm },
                        &mut the_min,
                        &mut the_max,
                    );
                    // HLRAlgo::AddMinMax(TheMin, TheMax, ShapMin, ShapMax);
                    crate::hlr::algo::hlr_algo::HLRAlgo::add_min_max(
                        &the_min,
                        &the_max,
                        &mut shap_min,
                        &mut shap_max,
                    );
                    f += 1;
                }
                // HLRAlgo::EncodeMinMax(ShapMin, ShapMax, MinMaxShap);
                crate::hlr::algo::hlr_algo::HLRAlgo::encode_min_max(
                    &shap_min,
                    &shap_max,
                    &mut min_max_shap,
                );
                // SB.UpdateMinMax(MinMaxShap);
                self.my_shapes[i - 1].update_min_max(&min_max_shap);
            }
        }
    }

    /// OCCT Load(S, SData, nbIso) (cxx L189-195) — add the shape S with its
    /// ShapeData (the setter overload keeps the OCCT argument order).
    pub fn load_with_data(&mut self, s: &Arc<crate::hlr::topo_brep::out_liner::OutLiner>, s_data: SDataHandle, nb_iso: i32) {
        // myShapes.Append(HLRBRep_ShapeBounds(S, SData, nbIso, 0, 0, 0, 0, 0, 0));
        self.my_shapes
            .push(ShapeBounds::new_with_data(s.clone(), s_data, nb_iso, 0, 0, 0, 0, 0, 0));
        // myDS.Nullify();
        self.my_ds = None;
    }

    /// OCCT Load(S, nbIso) (cxx L199-203) — add the shape S.
    pub fn load(&mut self, s: &Arc<crate::hlr::topo_brep::out_liner::OutLiner>, nb_iso: i32) {
        // myShapes.Append(HLRBRep_ShapeBounds(S, nbIso, 0, 0, 0, 0, 0, 0));
        self.my_shapes
            .push(ShapeBounds::new(s.clone(), nb_iso, 0, 0, 0, 0, 0, 0));
        // myDS.Nullify();
        self.my_ds = None;
    }

    /// OCCT Index(S) (cxx L207-220) — return the index of the Shape S and
    /// return 0 if the Shape S is not found.
    pub fn index(&self, s: &Arc<crate::hlr::topo_brep::out_liner::OutLiner>) -> usize {
        // int n = myShapes.Length();
        let n = self.my_shapes.len();

        // for (int i = 1; i <= n; i++) { if (myShapes(i).Shape() == S) return i; }
        for i in 1..=n {
            // the OCCT handle equality is pointer identity (Arc::ptr_eq).
            if Arc::ptr_eq(self.my_shapes[i - 1].shape(), s) {
                return i;
            }
        }

        // return 0;
        0
    }

    /// OCCT Remove(I) (cxx L224-232) — remove the Shape of index I.
    pub fn remove(&mut self, i: usize) {
        // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
        //                              "HLRBRep_InternalAlgo::Remove : unknown Shape");
        if i == 0 || i > self.my_shapes.len() {
            panic!("HLRBRep_InternalAlgo::Remove : unknown Shape");
        }
        // myShapes.Remove(I);
        self.my_shapes.remove(i - 1);

        // myMapOfShapeTool.Clear();
        self.my_map_of_shape_tool.clear();
        // myDS.Nullify();
        self.my_ds = None;
    }

    /// OCCT ShapeData(I, SData) (cxx L236-242) — change the Shape Data of
    /// the Shape of index I.
    pub fn shape_data(&mut self, i: usize, s_data: SDataHandle) {
        // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
        //   "HLRBRep_InternalAlgo::ShapeData : unknown Shape");
        if i == 0 || i > self.my_shapes.len() {
            panic!("HLRBRep_InternalAlgo::ShapeData : unknown Shape");
        }

        // myShapes(I).ShapeData(SData);
        self.my_shapes[i - 1].set_shape_data(s_data);
    }

    /// OCCT SeqOfShapeBounds() (cxx L246-249).
    pub fn seq_of_shape_bounds(&mut self) -> &mut Vec<ShapeBounds> {
        &mut self.my_shapes
    }

    /// OCCT NbShapes() (cxx L253-256).
    pub fn nb_shapes(&self) -> usize {
        self.my_shapes.len()
    }

    /// OCCT ShapeBounds(I) (cxx L260-266) — the non-const reference getter.
    pub fn shape_bounds(&mut self, i: usize) -> &mut ShapeBounds {
        // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
        //   "HLRBRep_InternalAlgo::ShapeBounds : unknown Shape");
        if i == 0 || i > self.my_shapes.len() {
            panic!("HLRBRep_InternalAlgo::ShapeBounds : unknown Shape");
        }

        &mut self.my_shapes[i - 1]
    }

    /// OCCT InitEdgeStatus() (cxx L270-348) — init the status of the
    /// selected edges depending of the back faces of a closed shell.
    pub fn init_edge_status(&mut self) {
        // bool visible;
        let mut visible: bool;
        // HLRBRep_FaceIterator faceIt;
        let mut face_it = FaceIterator::new(std::ptr::null_mut());

        let ds = self.my_ds.as_mut().expect("HLRBRep_InternalAlgo::InitEdgeStatus : null DS");
        // (the raw-slice mirrors of the OCCT aEDataArray / aFDataArray member
        // references — the Data element-pointer precedent; ne / nf cached as
        // in OCCT.)
        let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
        let a_f_data: *mut FaceData<'static> = ds.f_data_array_mut().as_mut_ptr();
        // int ne = myDS->NbEdges(); int nf = myDS->NbFaces();
        let ne = ds.nb_edges() as i32;
        let nf = ds.nb_faces() as i32;

        // for (int e = 1; e <= ne; e++)
        let mut e: i32 = 1;
        while e <= ne {
            // HLRBRep_EdgeData& ed = aEDataArray.ChangeValue(e);
            let ed = unsafe { &mut *a_e_data.add((e - 1) as usize) };
            if ed.selected() {
                // ed.Status().ShowAll();
                ed.status().show_all();
            }
            e += 1;
        }
        //  for (int f = 1; f <= nf; f++) {
        let mut f: i32;
        // for (f = 1; f <= nf; f++)
        f = 1;
        while f <= nf {
            // HLRBRep_FaceData& fd = aFDataArray.ChangeValue(f);
            let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
            if fd.selected() {
                // for (faceIt.InitEdge(fd); faceIt.MoreEdge(); faceIt.NextEdge())
                let fd_ptr: *mut FaceData<'static> = fd;
                face_it.init_edge(fd_ptr);
                while face_it.more_edge() {
                    // HLRBRep_EdgeData* edf = &(myDS->EDataArray().ChangeValue(faceIt.Edge()));
                    let edf = unsafe {
                        &mut *a_e_data.add((face_it.edge() - 1) as usize)
                    };
                    if edf.selected() {
                        // edf->Status().HideAll();
                        edf.status().hide_all();
                    }
                    face_it.next_edge();
                }
            }
            f += 1;
        }

        // for (f = 1; f <= nf; f++)
        f = 1;
        while f <= nf {
            // HLRBRep_FaceData& fd = aFDataArray.ChangeValue(f);
            let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
            // visible = true;
            visible = true;
            // if (fd.Selected() && fd.Closed())
            if fd.selected() && fd.closed() {
                if fd.side() {
                    // visible = false;
                    visible = false;
                } else if !fd.with_out_l() {
                    // switch (fd.Orientation())
                    match fd.orientation() {
                        rcad_kernel::topods::Orientation::Reversed => {
                            // visible = fd.Back();
                            visible = fd.back();
                        }
                        rcad_kernel::topods::Orientation::Forward => {
                            // visible = !fd.Back();
                            visible = !fd.back();
                        }
                        // case TopAbs_EXTERNAL: case TopAbs_INTERNAL: visible = true;
                        rcad_kernel::topods::Orientation::External
                        | rcad_kernel::topods::Orientation::Internal => {
                            visible = true;
                        }
                    }
                }
            }
            // if (visible)
            if visible {
                // for (faceIt.InitEdge(fd); faceIt.MoreEdge(); faceIt.NextEdge())
                let fd_ptr: *mut FaceData<'static> = fd;
                face_it.init_edge(fd_ptr);
                while face_it.more_edge() {
                    // int E = faceIt.Edge();
                    let e_i = face_it.edge();
                    // HLRBRep_EdgeData* edf = &(myDS->EDataArray().ChangeValue(E));
                    let edf = unsafe { &mut *a_e_data.add((e_i - 1) as usize) };
                    // if (edf->Selected() && !edf->Vertical())
                    if edf.selected() && !edf.vertical() {
                        // edf->Status().ShowAll();
                        edf.status().show_all();
                    }
                    face_it.next_edge();
                }
            }
            f += 1;
        }
    }

    /// OCCT Select() (cxx L352-373) — select all the DataStructure.
    pub fn select(&mut self) {
        // if (!myDS.IsNull())
        if let Some(ds) = self.my_ds.as_mut() {
            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            let a_f_data: *mut FaceData<'static> = ds.f_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges(); int nf = myDS->NbFaces();
            let ne = ds.nb_edges() as i32;
            let nf = ds.nb_faces() as i32;

            // for (int e = 1; e <= ne; e++) { ... ed.Selected(true); }
            let mut e: i32 = 1;
            while e <= ne {
                let ed = unsafe { &mut *a_e_data.add((e - 1) as usize) };
                ed.set_selected(true);
                e += 1;
            }

            // for (int f = 1; f <= nf; f++) { ... fd.Selected(true); }
            let mut f: i32 = 1;
            while f <= nf {
                let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
                fd.set_selected(true);
                f += 1;
            }
        }
    }

    /// OCCT Select(const int I) (cxx L377-404) — select only the Shape of
    /// index I.
    pub fn select_shape(&mut self, i: usize) {
        // if (!myDS.IsNull())
        if let Some(ds) = self.my_ds.as_mut() {
            // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
            //   "HLRBRep_InternalAlgo::Select : unknown Shape");
            if i == 0 || i > self.my_shapes.len() {
                panic!("HLRBRep_InternalAlgo::Select : unknown Shape");
            }

            // int v1, v2, e1, e2, f1, f2;
            // myShapes(I).Bounds(v1, v2, e1, e2, f1, f2);
            let (_v1, _v2, e1, e2, f1, f2) = self.my_shapes[i - 1].bounds();

            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            let a_f_data: *mut FaceData<'static> = ds.f_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges(); int nf = myDS->NbFaces();
            let ne = ds.nb_edges() as i32;
            let nf = ds.nb_faces() as i32;

            // for (int e = 1; e <= ne; e++) { ed.Selected(e >= e1 && e <= e2); }
            let mut e: i32 = 1;
            while e <= ne {
                let ed = unsafe { &mut *a_e_data.add((e - 1) as usize) };
                ed.set_selected(e >= e1 && e <= e2);
                e += 1;
            }

            // for (int f = 1; f <= nf; f++) { fd.Selected(f >= f1 && f <= f2); }
            let mut f: i32 = 1;
            while f <= nf {
                let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
                fd.set_selected(f >= f1 && f <= f2);
                f += 1;
            }
        }
    }

    /// OCCT SelectEdge(const int I) (cxx L408-427) — select only the edges
    /// of the Shape S.
    pub fn select_edge(&mut self, i: usize) {
        // if (!myDS.IsNull())
        if let Some(ds) = self.my_ds.as_mut() {
            // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
            //   "HLRBRep_InternalAlgo::SelectEdge : unknown Shape");
            if i == 0 || i > self.my_shapes.len() {
                panic!("HLRBRep_InternalAlgo::SelectEdge : unknown Shape");
            }

            // int v1, v2, e1, e2, f1, f2;
            // myShapes(I).Bounds(v1, v2, e1, e2, f1, f2);
            let (_v1, _v2, e1, e2, _f1, _f2) = self.my_shapes[i - 1].bounds();

            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges();
            let ne = ds.nb_edges() as i32;

            // for (int e = 1; e <= ne; e++) { ed.Selected(e >= e1 && e <= e2); }
            let mut e: i32 = 1;
            while e <= ne {
                let ed = unsafe { &mut *a_e_data.add((e - 1) as usize) };
                ed.set_selected(e >= e1 && e <= e2);
                e += 1;
            }
        }
    }

    /// OCCT SelectFace(const int I) (cxx L431-450) — select only the faces
    /// of the Shape S.
    pub fn select_face(&mut self, i: usize) {
        // if (!myDS.IsNull())
        if let Some(ds) = self.my_ds.as_mut() {
            // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
            //   "HLRBRep_InternalAlgo::SelectFace : unknown Shape");
            if i == 0 || i > self.my_shapes.len() {
                panic!("HLRBRep_InternalAlgo::SelectFace : unknown Shape");
            }

            // int v1, v2, e1, e2, f1, f2;
            // myShapes(I).Bounds(v1, v2, e1, e2, f1, f2);
            let (_v1, _v2, _e1, _e2, f1, f2) = self.my_shapes[i - 1].bounds();

            let a_f_data: *mut FaceData<'static> = ds.f_data_array_mut().as_mut_ptr();
            // int nf = myDS->NbFaces();
            let nf = ds.nb_faces() as i32;

            // for (int f = 1; f <= nf; f++) { fd.Selected(f >= f1 && f <= f2); }
            let mut f: i32 = 1;
            while f <= nf {
                let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
                fd.set_selected(f >= f1 && f <= f2);
                f += 1;
            }
        }
    }

    /// OCCT ShowAll() (cxx L454-467) — set to visible all the edges.
    pub fn show_all(&mut self) {
        // if (!myDS.IsNull())
        if let Some(ds) = self.my_ds.as_mut() {
            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges();
            let ne = ds.nb_edges() as i32;

            // for (int ie = 1; ie <= ne; ie++) { ... ed.Status().ShowAll(); }
            let mut ie: i32 = 1;
            while ie <= ne {
                let ed = unsafe { &mut *a_e_data.add((ie - 1) as usize) };
                ed.status().show_all();
                ie += 1;
            }
        }
    }

    /// OCCT ShowAll(const int I) (cxx L471-492) — set to visible all the
    /// edges of the Shape S.
    pub fn show_all_shape(&mut self, i: usize) {
        // if (!myDS.IsNull())
        if self.my_ds.is_some() {
            // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
            //   "HLRBRep_InternalAlgo::ShowAll : unknown Shape");
            if i == 0 || i > self.my_shapes.len() {
                panic!("HLRBRep_InternalAlgo::ShowAll : unknown Shape");
            }

            // Select(I);
            self.select_shape(i);

            let ds = self.my_ds.as_mut().expect("HLRBRep_InternalAlgo::ShowAll : null DS");
            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges();
            let ne = ds.nb_edges() as i32;

            // for (int e = 1; e <= ne; e++) { if (ed.Selected()) ed.Status().ShowAll(); }
            let mut e: i32 = 1;
            while e <= ne {
                let ed = unsafe { &mut *a_e_data.add((e - 1) as usize) };
                if ed.selected() {
                    ed.status().show_all();
                }
                e += 1;
            }
        }
    }

    /// OCCT HideAll() (cxx L496-509) — set to hide all the edges.
    pub fn hide_all(&mut self) {
        // if (!myDS.IsNull())
        if let Some(ds) = self.my_ds.as_mut() {
            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges();
            let ne = ds.nb_edges() as i32;

            // for (int ie = 1; ie <= ne; ie++) { ... ed.Status().HideAll(); }
            let mut ie: i32 = 1;
            while ie <= ne {
                let ed = unsafe { &mut *a_e_data.add((ie - 1) as usize) };
                ed.status().hide_all();
                ie += 1;
            }
        }
    }

    /// OCCT HideAll(const int I) (cxx L513-534) — set to hide all the edges
    /// of the Shape S.
    pub fn hide_all_shape(&mut self, i: usize) {
        // if (!myDS.IsNull())
        if self.my_ds.is_some() {
            // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
            //   "HLRBRep_InternalAlgo::HideAll : unknown Shape");
            if i == 0 || i > self.my_shapes.len() {
                panic!("HLRBRep_InternalAlgo::HideAll : unknown Shape");
            }

            // Select(I);
            self.select_shape(i);

            let ds = self.my_ds.as_mut().expect("HLRBRep_InternalAlgo::HideAll : null DS");
            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges();
            let ne = ds.nb_edges() as i32;

            // for (int e = 1; e <= ne; e++) { if (ed.Selected()) ed.Status().HideAll(); }
            let mut e: i32 = 1;
            while e <= ne {
                let ed = unsafe { &mut *a_e_data.add((e - 1) as usize) };
                if ed.selected() {
                    ed.status().hide_all();
                }
                e += 1;
            }
        }
    }

    /// OCCT PartialHide() (cxx L538-556) — own hiding of all the shapes of
    /// the DataStructure without hiding by each other.
    pub fn partial_hide(&mut self) {
        // if (!myDS.IsNull())
        if self.my_ds.is_some() {
            // int i, n = myShapes.Length();
            let n = self.my_shapes.len();

            if self.my_debug {
                println!(" Partial hiding \n");
            }

            // for (i = 1; i <= n; i++) { Hide(i); }
            for i in 1..=n {
                self.hide_shape(i);
            }

            // Select();
            self.select();
        }
    }

    /// OCCT Hide() (cxx L560-589) — hide all the DataStructure.
    pub fn hide(&mut self) {
        // if (!myDS.IsNull())
        if self.my_ds.is_some() {
            // int i, j, n = myShapes.Length();
            let n = self.my_shapes.len();

            if self.my_debug {
                println!(" Total hiding");
            }

            // for (i = 1; i <= n; i++) { Hide(i); }
            for i in 1..=n {
                self.hide_shape(i);
            }

            // for (i = 1; i <= n; i++)
            for i in 1..=n {
                // for (j = 1; j <= n; j++)
                for j in 1..=n {
                    // if (i != j) { Hide(i, j); }
                    if i != j {
                        self.hide_shape_by_shape(i, j);
                    }
                }
            }

            // Select();
            self.select();
        }
    }

    /// OCCT Hide(const int I) (cxx L593-609) — hide the Shape S by itself.
    pub fn hide_shape(&mut self, i: usize) {
        // if (!myDS.IsNull())
        if self.my_ds.is_some() {
            // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length(),
            //   "HLRBRep_InternalAlgo::Hide : unknown Shape");
            if i == 0 || i > self.my_shapes.len() {
                panic!("HLRBRep_InternalAlgo::Hide : unknown Shape");
            }

            if self.my_debug {
                println!(" hiding the shape {} by itself", i);
            }

            // Select(I); InitEdgeStatus(); HideSelected(I, true);
            self.select_shape(i);
            self.init_edge_status();
            self.hide_selected(i, true);
        }
    }

    /// OCCT Hide(const int I, const int J) (cxx L613-654) — hide the Shape
    /// S1 by the shape S2.
    pub fn hide_shape_by_shape(&mut self, i: usize, j: usize) {
        // if (!myDS.IsNull())
        if self.my_ds.is_some() {
            // Standard_OutOfRange_Raise_if(I == 0 || I > myShapes.Length() ||
            //                              J == 0 || J > myShapes.Length(),
            //   "HLRBRep_InternalAlgo::Hide : unknown Shapes");
            if i == 0 || i > self.my_shapes.len() || j == 0 || j > self.my_shapes.len() {
                panic!("HLRBRep_InternalAlgo::Hide : unknown Shapes");
            }

            // if (I == J) { Hide(I); }
            if i == j {
                self.hide_shape(i);
            } else {
                // HLRAlgo_EdgesBlock::MinMaxIndices* MinMaxShBI = &myShapes(I).MinMax();
                // HLRAlgo_EdgesBlock::MinMaxIndices* MinMaxShBJ = &myShapes(J).MinMax();
                // (the MinMax reads copy the encoded boxes — the gate below
                // only reads them.)
                let min_max_bi = *self.my_shapes[i - 1].min_max();
                let min_max_bj = *self.my_shapes[j - 1].min_max();
                // the encoded-box overlap gate (the wrapping subtraction keeps
                // the OCCT int overflow form).
                if ((min_max_bj.max[0].wrapping_sub(min_max_bi.min[0])) & MASK_80008000) == 0
                    && ((min_max_bi.max[0].wrapping_sub(min_max_bj.min[0])) & MASK_80008000) == 0
                    && ((min_max_bj.max[1].wrapping_sub(min_max_bi.min[1])) & MASK_80008000) == 0
                    && ((min_max_bi.max[1].wrapping_sub(min_max_bj.min[1])) & MASK_80008000) == 0
                    && ((min_max_bj.max[2].wrapping_sub(min_max_bi.min[2])) & MASK_80008000) == 0
                    && ((min_max_bi.max[2].wrapping_sub(min_max_bj.min[2])) & MASK_80008000) == 0
                    && ((min_max_bj.max[3].wrapping_sub(min_max_bi.min[3])) & MASK_80008000) == 0
                    && ((min_max_bi.max[3].wrapping_sub(min_max_bj.min[3])) & MASK_80008000) == 0
                    && ((min_max_bj.max[4].wrapping_sub(min_max_bi.min[4])) & MASK_80008000) == 0
                    && ((min_max_bi.max[4].wrapping_sub(min_max_bj.min[4])) & MASK_80008000) == 0
                    && ((min_max_bj.max[5].wrapping_sub(min_max_bi.min[5])) & MASK_80008000) == 0
                    && ((min_max_bi.max[5].wrapping_sub(min_max_bj.min[5])) & MASK_80008000) == 0
                    && ((min_max_bj.max[6].wrapping_sub(min_max_bi.min[6])) & MASK_80008000) == 0
                    && ((min_max_bi.max[6].wrapping_sub(min_max_bj.min[6])) & MASK_80008000) == 0
                    && ((min_max_bj.max[7].wrapping_sub(min_max_bi.min[7])) & MASK_80008000) == 0
                    && ((min_max_bi.max[7].wrapping_sub(min_max_bj.min[7])) & MASK_80008000) == 0
                {
                    if self.my_debug {
                        print!(" hiding the shape {}", i);
                        println!(" by the shape : {}", j);
                    }
                    // SelectEdge(I); SelectFace(J); HideSelected(I, false);
                    self.select_edge(i);
                    self.select_face(j);
                    self.hide_selected(i, false);
                }
            }
        }
    }

    /// OCCT HideSelected(const int I, const bool SideFace) (cxx L658-999) —
    /// first if SideFace own hiding of the side faces.  After hiding of the
    /// selected parts of the DataStructure.
    fn hide_selected(&mut self, i: usize, side_face: bool) {
        // int e, f, j, nbVisEdges, nbSelEdges, nbSelFaces, nbCache;
        // int nbFSide, nbFSimp;
        // (the OCCT_DEBUG counter resets of cxx L663-674 are not compiled in
        // the reference build — omitted with the counters note.)

        // HLRBRep_ShapeBounds& SB = myShapes(I);
        // int v1, v2, e1, e2, f1, f2;
        // SB.Bounds(v1, v2, e1, e2, f1, f2);
        let (_v1, _v2, e1, e2, _f1, _f2) = self.my_shapes[i - 1].bounds();

        // if (e2 >= e1)
        if e2 >= e1 {
            // myDS->InitBoundSort(SB.MinMax(), e1, e2);
            let ds = self.my_ds.as_mut().expect("HideSelected : null DS");
            let sb_min_max = self.my_shapes[i - 1].min_max();
            ds.init_bound_sort(sb_min_max, e1 as usize, e2 as usize);
            // HLRBRep_Hider Cache(myDS); — the Cache keeps the DS handle
            // alias (the raw pointer of the landed Hider); the reborrow ends
            // with the ctor call and the aEDataArray reference below aliases
            // the Cache DS as in OCCT.
            let mut cache = Hider::new(&mut *ds);
            let a_e_data: *mut EdgeData<'static> = ds.e_data_array_mut().as_mut_ptr();
            let a_f_data: *mut FaceData<'static> = ds.f_data_array_mut().as_mut_ptr();
            // int ne = myDS->NbEdges(); int nf = myDS->NbFaces();
            let ne = ds.nb_edges() as i32;
            let nf = ds.nb_faces() as i32;

            // if (myDebug) { ... } (cxx L689-748 — the plain (non-ifdef)
            // statistics block.)
            if self.my_debug {
                let mut nb_vis_edges: i32 = 0;
                let mut nb_sel_edges: i32 = 0;
                let mut nb_sel_faces: i32 = 0;
                let mut nb_cache: i32 = 0;
                let mut nb_f_side: i32 = 0;
                let mut nb_f_simp: i32 = 0;

                // for (e = 1; e <= ne; e++)
                let mut e: i32 = 1;
                while e <= ne {
                    let ed = unsafe { &mut *a_e_data.add((e - 1) as usize) };
                    if ed.selected() {
                        nb_sel_edges += 1;
                        if !ed.status_ref().all_hidden() {
                            nb_vis_edges += 1;
                        }
                    }
                    e += 1;
                }

                // for (f = 1; f <= nf; f++)
                let mut f: i32 = 1;
                while f <= nf {
                    let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
                    if fd.selected() {
                        nb_sel_faces += 1;
                        if fd.hiding() {
                            nb_cache += 1;
                        }
                        if fd.side() {
                            nb_f_side += 1;
                        }
                        if fd.simple() {
                            nb_f_simp += 1;
                        }
                    }
                    f += 1;
                }

                println!();
                println!("Vertices  : {:5}", ds.nb_vertices());
                println!("Edges     : {:5} , ", ds.nb_edges());
                print!("Selected  : {:5} , ", nb_sel_edges);
                println!("Visible   : {:5}", nb_vis_edges);
                print!("Faces     : {:5} , ", ds.nb_faces());
                print!("Selected  : {:5} , ", nb_sel_faces);
                println!("Simple    : {:5}", nb_f_simp);
                if side_face {
                    print!("Side      : {:5} , ", nb_f_side);
                }
                println!("Cachantes : {:5} \n", nb_cache);
            }

            // if (nf == 0) { return; }
            if nf == 0 {
                return;
            }

            // int QWE = 0, QWEQWE; QWEQWE = nf / 10;
            let mut qwe: i32 = 0;
            let qwe_qwe: i32 = nf / 10;

            // if (SideFace)
            if side_face {
                // j = 0;
                let mut j: i32 = 0;

                // for (f = 1; f <= nf; f++)
                let mut f: i32 = 1;
                while f <= nf {
                    // HLRBRep_FaceData& fd = aFDataArray.ChangeValue(f);
                    let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
                    if fd.selected() {
                        if fd.side() {
                            if INTERNAL_ALGO_TRACE10 {
                                // if (++QWE > QWEQWE) { QWE = 0; if (myDebug) cout << "*"; }
                                qwe += 1;
                                if qwe > qwe_qwe {
                                    qwe = 0;
                                    if self.my_debug {
                                        print!("*");
                                    }
                                }
                            } else {
                                // if (myDebug && HLRBRep_InternalAlgo_TRACE)
                                if self.my_debug && INTERNAL_ALGO_TRACE {
                                    j += 1;
                                    println!(" OwnHiding {} of face : {}", j, f);
                                }
                            }
                            // Cache.OwnHiding(f);
                            cache.own_hiding(f as usize);
                        }
                    }
                    f += 1;
                }
            }

            //--
            // NCollection_Array1<int>    Val(1, nf);
            // NCollection_Array1<double> Size(1, nf);
            // NCollection_Array1<int>    Index(1, nf);
            let mut val: Vec<i32> = vec![0; nf as usize];
            let mut size: Vec<f64> = vec![0.0; nf as usize];
            let mut index: Vec<i32> = vec![0; nf as usize];

            // for (f = 1; f <= nf; f++)
            let mut f: i32 = 1;
            while f <= nf {
                // HLRBRep_FaceData& fd = aFDataArray.ChangeValue(f);
                let fd = unsafe { &mut *a_f_data.add((f - 1) as usize) };
                // the surface-type ranking
                if fd.plane() {
                    val[(f - 1) as usize] = 10;
                } else if fd.cylinder() {
                    val[(f - 1) as usize] = 9;
                } else if fd.cone() {
                    val[(f - 1) as usize] = 8;
                } else if fd.sphere() {
                    val[(f - 1) as usize] = 7;
                } else if fd.torus() {
                    val[(f - 1) as usize] = 6;
                } else {
                    val[(f - 1) as usize] = 0;
                }
                if fd.cut() {
                    val[(f - 1) as usize] -= 10;
                }
                if fd.side() {
                    val[(f - 1) as usize] -= 100;
                }
                if fd.with_out_l() {
                    val[(f - 1) as usize] -= 20;
                }

                // Size(f) = fd.Size();
                size[(f - 1) as usize] = fd.size();
                f += 1;
            }

            // for (int tt = 1; tt <= nf; tt++) { Index(tt) = tt; }
            let mut tt: i32 = 1;
            while tt <= nf {
                index[(tt - 1) as usize] = tt;
                tt += 1;
            }

            //-- ======================================================================
            // (the commented-out TriOk bubble sort of cxx L848-866 is kept
            // un-translated, as in OCCT.)
            //-- ======================================================================

            // if (nf > 2)
            if nf > 2 {
                // int i, ir, k, l; int rra;
                let mut i: i32;
                let mut ir: i32;
                let mut k: i32;
                let mut l: i32;
                let mut rra: i32;
                // l = (nf >> 1) + 1; ir = nf;
                l = (nf >> 1) + 1;
                ir = nf;
                // for (;;)
                loop {
                    if l > 1 {
                        // rra = Index(--l);
                        l -= 1;
                        rra = index[(l - 1) as usize];
                    } else {
                        // rra = Index(ir); Index(ir) = Index(1);
                        rra = index[(ir - 1) as usize];
                        index[(ir - 1) as usize] = index[0];
                        // if (--ir == 1) { Index(1) = rra; break; }
                        ir -= 1;
                        if ir == 1 {
                            index[0] = rra;
                            break;
                        }
                    }
                    // i = l; k = l + l;
                    i = l;
                    k = l + l;
                    // while (k <= ir)
                    while k <= ir {
                        if k < ir {
                            // if (Val(Index(k)) > Val(Index(k + 1)))
                            if val[(index[(k - 1) as usize] - 1) as usize]
                                > val[(index[(k + 1 - 1) as usize] - 1) as usize]
                            {
                                k += 1;
                            }
                            // else if (Val(Index(k)) == Val(Index(k + 1)))
                            else if val[(index[(k - 1) as usize] - 1) as usize]
                                == val[(index[(k + 1 - 1) as usize] - 1) as usize]
                            {
                                // if (Size(Index(k)) > Size(Index(k + 1)))
                                if size[(index[(k - 1) as usize] - 1) as usize]
                                    > size[(index[(k + 1 - 1) as usize] - 1) as usize]
                                {
                                    k += 1;
                                }
                            }
                        }
                        // if (Val(rra) > Val(Index(k)))
                        if val[(rra - 1) as usize] > val[(index[(k - 1) as usize] - 1) as usize] {
                            // Index(i) = Index(k); i = k; k <<= 1;
                            index[(i - 1) as usize] = index[(k - 1) as usize];
                            i = k;
                            k <<= 1;
                        }
                        // else if ((Val(rra) == Val(Index(k))) && (Size(rra) > Size(Index(k))))
                        else if val[(rra - 1) as usize] == val[(index[(k - 1) as usize] - 1) as usize]
                            && size[(rra - 1) as usize]
                                > size[(index[(k - 1) as usize] - 1) as usize]
                        {
                            // Index(i) = Index(k); i = k; k <<= 1;
                            index[(i - 1) as usize] = index[(k - 1) as usize];
                            i = k;
                            k <<= 1;
                        } else {
                            // k = ir + 1;
                            k = ir + 1;
                        }
                    }
                    // Index(i) = rra;
                    index[(i - 1) as usize] = rra;
                }
            }

            // j = 0;
            let mut j: i32 = 0;

            // QWE = 0;
            qwe = 0;
            // for (f = 1; f <= nf; f++)
            let mut f: i32 = 1;
            while f <= nf {
                // int fi = Index(f);
                let fi = index[(f - 1) as usize];
                // HLRBRep_FaceData& fd = aFDataArray.ChangeValue(fi);
                let fd = unsafe { &mut *a_f_data.add((fi - 1) as usize) };
                if fd.selected() {
                    if fd.hiding() {
                        // if (TRACE10 && !static_cast<bool>(TRACE))
                        if INTERNAL_ALGO_TRACE10 && !INTERNAL_ALGO_TRACE {
                            // if (++QWE > QWEQWE) { if (myDebug) cout << "."; QWE = 0; }
                            qwe += 1;
                            if qwe > qwe_qwe {
                                if self.my_debug {
                                    print!(".");
                                }
                                qwe = 0;
                            }
                        }
                        // else if (myDebug && HLRBRep_InternalAlgo_TRACE)
                        else if self.my_debug && INTERNAL_ALGO_TRACE {
                            // static int rty = 0; — the OCCT function static
                            // keeps the thread-local form (the Data counters
                            // precedent).
                            use std::cell::Cell;
                            thread_local! {
                                static RTY: Cell<i32> = const { Cell::new(0) };
                            }
                            j += 1;
                            print!("{:6}", fi);
                            // fflush(stdout); — the println flushes on the
                            // line break of the rty reset.
                            RTY.with(|rty| {
                                rty.set(rty.get() + 1);
                                if rty.get() > 25 {
                                    rty.set(0);
                                    println!();
                                }
                            });
                        }
                        // Cache.Hide(fi, myMapOfShapeTool);
                        cache.hide(fi as usize, &mut self.my_map_of_shape_tool);
                    }
                }
                f += 1;
            }

            // (the OCCT_DEBUG statistics of cxx L968-997 are not compiled in
            // the reference build — omitted with the counters note.)
        }
    }

    /// OCCT Debug(const bool deb) (cxx L1003-1006).
    pub fn set_debug(&mut self, deb: bool) {
        // myDebug = deb;
        self.my_debug = deb;
    }

    /// OCCT Debug() (cxx L1010-1013).
    pub fn debug(&self) -> bool {
        self.my_debug
    }
}

#[cfg(test)]
impl InternalAlgo {
    /// The test-only DS installation (the myDS slot is only reachable via
    /// Update in OCCT; the anchors install a hand-built Data directly).
    pub(crate) fn install_ds_for_test(&mut self, ds: Data<'static>) {
        self.my_ds = Some(Box::new(ds));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::algo::edges_block::EdgesBlock;
    use crate::hlr::algo::wires_block::WiresBlock;
    use crate::hlr::brep::face_data::FaceData;
    use crate::hlr::topo_brep::out_liner::OutLiner;
    use std::sync::Arc as StdArc;

    /// A selected EdgeData with the requested AllHidden status (the OCCT
    /// EdgeData ctor default Selected(true)).
    fn edge(all_hidden: bool) -> EdgeData<'static> {
        let mut e = EdgeData::new();
        if all_hidden {
            e.status().hide_all();
        } else {
            e.status().show_all();
        }
        e
    }

    /// A selected face with one wire of the given edge numbers and the
    /// Closed / Side flags (everything else defaults).
    fn face_with(edge_numbers: &[i32], closed: bool, side: bool) -> FaceData<'static> {
        let mut fd = FaceData::new();
        fd.set_closed(closed);
        fd.set_side(side);
        let mut wb = WiresBlock::new(1);
        let mut eb = EdgesBlock::new(edge_numbers.len());
        for (i, &en) in edge_numbers.iter().enumerate() {
            eb.set_edge(i + 1, en);
        }
        wb.set(1, eb);
        *fd.change_wires() = Some(StdArc::new(wb));
        fd
    }

    /// A Data of ne edges / nf faces from the builders.
    fn fixture_data(
        edges: Vec<EdgeData<'static>>,
        faces: Vec<FaceData<'static>>,
    ) -> Data<'static> {
        let mut d = Data::new(0, edges.len(), faces.len());
        let arr = d.e_data_array_mut();
        for (i, e) in edges.into_iter().enumerate() {
            arr[i] = e;
        }
        let arr = d.f_data_array_mut();
        for (i, f) in faces.into_iter().enumerate() {
            arr[i] = f;
        }
        d
    }

    /// The InternalAlgo with nb_shapes ShapeBounds entries of the given
    /// bounds over the fixture DS (the ShapeBounds ctor carries the bounds —
    /// the Update write-back form).
    fn algo_with(
        edges: Vec<EdgeData<'static>>,
        faces: Vec<FaceData<'static>>,
        bounds: (i32, i32, i32, i32, i32, i32),
        nb_shapes: usize,
    ) -> Box<InternalAlgo> {
        let mut ia = InternalAlgo::new();
        ia.install_ds_for_test(fixture_data(edges, faces));
        let (v1, v2, e1, e2, f1, f2) = bounds;
        for _ in 0..nb_shapes {
            ia.seq_of_shape_bounds().push(ShapeBounds::new(
                StdArc::new(OutLiner::new()),
                0,
                v1,
                v2,
                e1,
                e2,
                f1,
                f2,
            ));
        }
        Box::new(ia)
    }

    /// The ctor (cxx L50-53) / Load (L189-203) / Index (L207-220) / Remove
    /// (L224-232) / ShapeData (L236-242) / ShapeBounds (L260-266) round
    /// trip, the empty-list Index = 0, and the out-of-range raises.
    #[test]
    fn internal_algo_ctor_load_index_remove_round_trip() {
        // the ctor defaults (myDebug(false), the null DS, the empty sequence).
        let mut ia = InternalAlgo::new();
        assert!(!ia.debug());
        assert_eq!(ia.nb_shapes(), 0);
        assert!(ia.data_structure().is_none());
        // Projector(P) / Projector() round trip.
        let proj = crate::hlr::algo::projector::Projector::from_ax2(
            &rcad_kernel::math::gp::Ax2::new(
                glam::DVec3::ZERO,
                glam::DVec3::new(0.0, 0.0, 1.0),
                glam::DVec3::new(1.0, 0.0, 0.0),
            ),
        );
        ia.set_projector(&proj);
        assert!(!ia.projector().perspective());

        // Index on the empty list = 0 (cxx L207-220 tail).
        let s1 = StdArc::new(OutLiner::new());
        assert_eq!(ia.index(&s1), 0);

        // Load(S, nbIso): the ShapeBounds append with zero bounds and the DS
        // nullify.
        ia.load(&s1, 2);
        assert_eq!(ia.nb_shapes(), 1);
        assert!(ia.data_structure().is_none());
        assert_eq!(ia.index(&s1), 1);
        assert_eq!(ia.index(&StdArc::new(OutLiner::new())), 0);
        // NbOfIso round trip through ShapeBounds(I).
        assert_eq!(ia.shape_bounds(1).nb_of_iso(), 2);

        // Load(S, SData, nbIso) + ShapeData(I, SData).
        let s2 = StdArc::new(OutLiner::new());
        ia.load_with_data(&s2, Some(StdArc::new(7u32)), 0);
        assert_eq!(ia.nb_shapes(), 2);
        assert_eq!(ia.index(&s2), 2);
        ia.shape_data(2, Some(StdArc::new(11u32)));
        let sd = ia.shape_bounds(2).shape_data().clone().expect("SData");
        assert_eq!(sd.downcast_ref::<u32>(), Some(&11u32));

        // the out-of-range raises (Standard_OutOfRange_Raise_if).
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            ia.remove(0);
        }));
        assert!(result.is_err());
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            ia.remove(3);
        }));
        assert!(result.is_err());
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = ia.shape_bounds(0);
        }));
        assert!(result.is_err());

        // Remove(I): the shape disappears and the map / DS are cleared.
        ia.remove(2);
        assert_eq!(ia.nb_shapes(), 1);
        assert_eq!(ia.index(&s2), 0);
        assert_eq!(ia.index(&s1), 1);

        // Update on the empty sequence keeps the DS null (cxx L83 guard).
        let mut empty = InternalAlgo::new();
        empty.update();
        assert!(empty.data_structure().is_none());
    }

    /// InitEdgeStatus (cxx L270-348): the three phases — ShowAll on all the
    /// selected edges, HideAll of the selected face edges, and the closed
    /// back-face visibility flip.
    #[test]
    fn internal_algo_init_edge_status_bit_semantics() {
        // face 1: closed + side -> not visible -> its edges end HideAll.
        let mut ia = algo_with(
            vec![edge(false), edge(false)],
            vec![face_with(&[1, 2], true, true)],
            (0, 0, 1, 2, 1, 1),
            1,
        );
        ia.init_edge_status();
        let arr = ia.data_structure_mut().expect("DS").e_data_array_mut();
        assert!(!arr[0].status_ref().all_visible()); // the phase-B HideAll stuck
        assert!(arr[0].status_ref().all_hidden());
        assert!(arr[1].status_ref().all_hidden());

        // face 1: open (not closed) -> visible -> its edges end ShowAll.
        let mut ia = algo_with(
            vec![edge(false), edge(false)],
            vec![face_with(&[1, 2], false, false)],
            (0, 0, 1, 2, 1, 1),
            1,
        );
        ia.init_edge_status();
        let arr = ia.data_structure_mut().expect("DS").e_data_array_mut();
        assert!(arr[0].status_ref().all_visible());
        assert!(arr[1].status_ref().all_visible());

        // face 1: closed forward with Back=false -> visible -> ShowAll; the
        // reversed twin: closed + Back=true -> not visible -> HideAll.
        let mut fd = face_with(&[1, 2], true, false);
        fd.set_orientation(rcad_kernel::topods::Orientation::Forward);
        fd.set_back(false);
        let mut ia = algo_with(vec![edge(false), edge(false)], vec![fd], (0, 0, 1, 2, 1, 1), 1);
        ia.init_edge_status();
        let arr = ia.data_structure_mut().expect("DS").e_data_array_mut();
        assert!(arr[0].status_ref().all_visible());

        // the closed REVERSED faces follow visible = Back (cxx L321-323):
        // the seen-from-the-back reversed face stays visible, the
        // front-seen one is not.
        let mut fd = face_with(&[1, 2], true, false);
        fd.set_orientation(rcad_kernel::topods::Orientation::Reversed);
        fd.set_back(true);
        let mut ia = algo_with(vec![edge(false), edge(false)], vec![fd], (0, 0, 1, 2, 1, 1), 1);
        ia.init_edge_status();
        let arr = ia.data_structure_mut().expect("DS").e_data_array_mut();
        assert!(arr[0].status_ref().all_visible());

        let mut fd = face_with(&[1, 2], true, false);
        fd.set_orientation(rcad_kernel::topods::Orientation::Reversed);
        fd.set_back(false);
        let mut ia = algo_with(vec![edge(false), edge(false)], vec![fd], (0, 0, 1, 2, 1, 1), 1);
        ia.init_edge_status();
        let arr = ia.data_structure_mut().expect("DS").e_data_array_mut();
        assert!(arr[0].status_ref().all_hidden());
    }

    /// Select (cxx L352-450): the whole-DS and per-shape edge/face
    /// selection flags; ShowAll / HideAll (cxx L454-534) on the selected
    /// edges only.
    #[test]
    fn internal_algo_select_show_all_hide_all_flags() {
        let mut ia = algo_with(
            vec![edge(false), edge(false), edge(false)],
            vec![face_with(&[1], false, false), face_with(&[3], false, false)],
            (0, 0, 1, 2, 1, 2),
            1,
        );
        // Select(): everything selected.
        ia.select();
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array().iter().all(|e| e.selected()));
            assert!(ds.f_data_array().iter().all(|f| f.selected()));
        }

        // Select(1) with the ShapeBounds edges 1..2 / faces 1..2: edge 3 and
        // the bounds are the gate (face 2 lies inside 1..2).
        ia.select_shape(1);
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array()[0].selected());
            assert!(ds.e_data_array()[1].selected());
            assert!(!ds.e_data_array()[2].selected());
            assert!(ds.f_data_array()[0].selected());
            assert!(ds.f_data_array()[1].selected());
        }

        // SelectEdge(1) touches the edges only.
        ia.select_edge(1);
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array()[0].selected());
            assert!(!ds.e_data_array()[2].selected());
        }

        // SelectFace(1): the same bounds, faces only.
        ia.select_face(1);
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.f_data_array()[1].selected());
        }

        // HideAll(): every edge status flips to hidden.
        ia.hide_all();
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array().iter().all(|e| e.status_ref().all_hidden()));
        }

        // ShowAll(): every edge status flips back.
        ia.show_all();
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array().iter().all(|e| e.status_ref().all_visible()));
        }

        // ShowAll(1) / HideAll(1) apply to the whole DS here (edge 3 got
        // deselected by SelectEdge/SelectFace above, but status writes hit
        // only the selected ones: edge 3 keeps its state).
        ia.select_shape(1);
        ia.hide_all_shape(1);
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array()[0].status_ref().all_hidden());
            assert!(ds.e_data_array()[1].status_ref().all_hidden());
        }
        ia.show_all_shape(1);
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array()[0].status_ref().all_visible());
        }

        // the out-of-range raises.
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            ia.hide_all_shape(2);
        }));
        assert!(result.is_err());
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            ia.show_all_shape(0);
        }));
        assert!(result.is_err());
    }

    /// Hide / PartialHide / HideSelected control flow (cxx L538-654 +
    /// L658-999): the fixture has no Hiding face, so the Cache stays at the
    /// OwnHiding calls (no-ops); the trailing Select() of PartialHide / Hide
    /// and the encoded-box gate of Hide(I, J) are the asserted outcomes.
    #[test]
    fn internal_algo_hide_control_flow_reachability() {
        // -- Hide(I) by itself (SideFace = true): the selected side face
        //    takes the Cache.OwnHiding call; the non-hiding face skips
        //    Cache.Hide; the trailing state is the InitEdgeStatus result.
        let mut ia = algo_with(
            vec![edge(false)],
            vec![face_with(&[1], false, true)],
            (0, 0, 1, 1, 1, 1),
            1,
        );
        ia.hide_shape(1);
        {
            // the open face stays visible through the phase C ladder, so the
            // selected non-vertical edge ends ShowAll.
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array()[0].status_ref().all_visible());
        }
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            ia.hide_shape(0);
        }));
        assert!(result.is_err());

        // -- PartialHide: Hide(i) per shape then Select() -> everything
        //    selected.
        let mut ia = algo_with(
            vec![edge(false), edge(false)],
            vec![face_with(&[1], false, false), face_with(&[2], false, false)],
            (0, 0, 1, 2, 1, 2),
            2,
        );
        ia.partial_hide();
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array().iter().all(|e| e.selected()));
            assert!(ds.f_data_array().iter().all(|f| f.selected()));
        }

        // -- Hide() (total): the same trailing Select() through the pair
        //    loop.
        let mut ia = algo_with(
            vec![edge(false), edge(false)],
            vec![face_with(&[1], false, false), face_with(&[2], false, false)],
            (0, 0, 1, 2, 1, 2),
            2,
        );
        ia.hide();
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array().iter().all(|e| e.selected()));
        }

        // -- Hide(I, J): the encoded-box gate.  Two shapes with disjoint
        //    MinMax boxes skip the whole branch (SelectEdge / SelectFace /
        //    HideSelected are inside the gate, cxx L642-651) — the fixture
        //    keeps its ctor state (the EdgeData / FaceData Selected(true)).
        let mut ia = algo_with(
            vec![edge(false), edge(false)],
            vec![face_with(&[1], false, false), face_with(&[2], false, false)],
            (0, 0, 1, 1, 1, 1),
            2,
        );
        let mut mm = MinMaxIndices::default();
        mm.min = [0; 8];
        mm.max = [100; 8];
        ia.shape_bounds(1).update_min_max(&mm);
        let mut mm2 = MinMaxIndices::default();
        mm2.min = [10000; 8];
        mm2.max = [20000; 8];
        ia.shape_bounds(2).update_min_max(&mm2);
        ia.hide_shape_by_shape(1, 2);
        {
            let ds = ia.data_structure().expect("DS");
            // the gate skipped everything — the ctor Selected(true) stands.
            assert!(ds.e_data_array()[0].selected());
            assert!(ds.e_data_array()[1].selected());
            assert!(ds.f_data_array()[0].selected());
            assert!(ds.f_data_array()[1].selected());
        }

        // -- Hide(I, J) with I == J delegates to Hide(I).
        let mut ia =
            algo_with(vec![edge(false)], vec![face_with(&[1], false, true)], (0, 0, 1, 1, 1, 1), 1);
        ia.hide_shape_by_shape(1, 1);
        {
            let ds = ia.data_structure().expect("DS");
            assert!(ds.e_data_array()[0].selected());
        }

        // -- the overlapping-box pair passes the gate: SelectEdge(1) keeps
        //    the edge 1..1 of the bounds, SelectFace(2) the faces 2..2, and
        //    HideSelected(I, false) — the non-side hiding path with no
        //    Hiding face — is a no-op walk.
        let mut ia = algo_with(
            vec![edge(false), edge(false)],
            vec![face_with(&[1], false, false), face_with(&[2], false, false)],
            (0, 0, 1, 1, 2, 2),
            2,
        );
        let mm1 = MinMaxIndices {
            min: [0; 8],
            max: [100; 8],
        };
        ia.shape_bounds(1).update_min_max(&mm1);
        ia.shape_bounds(2).update_min_max(&mm1);
        ia.hide_shape_by_shape(1, 2);
        {
            let ds = ia.data_structure().expect("DS");
            // SelectEdge(1): edge 1 selected, edge 2 not.
            assert!(ds.e_data_array()[0].selected());
            assert!(!ds.e_data_array()[1].selected());
            // SelectFace(2): face 2 selected, face 1 not.
            assert!(!ds.f_data_array()[0].selected());
            assert!(ds.f_data_array()[1].selected());
        }
    }
}
