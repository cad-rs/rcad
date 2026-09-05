// OCCT HLRBRep_ShapeBounds (TKHLR/HLRBRep/HLRBRep_ShapeBounds.hxx L1-97
// + .cxx L1-96 + .lxx L1-70) — a Shape and the bounds of its vertices,
// edges and faces in the DataStructure.
//
// Translation notes:
// - OCCT handle<HLRTopoBRep_OutLiner> maps to Arc<OutLiner> (the O1/O2
//   contract encoding); the null handle of the default constructor stands
//   as an empty OutLiner (OriginalShape() null — the OCCT null-shape
//   behavior of the consumers).
// - OCCT handle<Standard_Transient> myShapeData (hxx L84) is OPAQUE at
//   this level — the rcad form is the SDataHandle
//   Option<Arc<dyn Any + Send + Sync>> (the InternalAlgo::Load(S, SData,
//   nbIso) rider; the null handle is None).
// - The OCCT out-parameter methods Sizes / Bounds (cxx L79-96, the
//   `int& NV` forms) map onto the returned tuples (the InternalAlgo
//   destructuring call shape).

use std::any::Any;
use std::sync::Arc;

use crate::hlr::algo::edges_block::MinMaxIndices;
use crate::hlr::topo_brep::out_liner::OutLiner;

/// OCCT `occ::handle<Standard_Transient>` — the opaque per-shape user data
/// riding with a ShapeBounds (the ShapeData slot).
pub type SDataHandle = Option<Arc<dyn Any + Send + Sync>>;

/// OCCT HLRBRep_ShapeBounds (hxx L31-93).
#[derive(Clone)]
pub struct ShapeBounds {
    /// OCCT myShape (hxx L83).
    my_shape: Arc<OutLiner>,
    /// OCCT myShapeData (hxx L84) — the null handle until a consumer binds
    /// it (ShapeData(SD)).
    my_shape_data: SDataHandle,
    /// OCCT myNbIso (hxx L85).
    my_nb_iso: i32,
    /// OCCT myVertStart / myVertEnd (hxx L86-87).
    my_vert_start: i32,
    my_vert_end: i32,
    /// OCCT myEdgeStart / myEdgeEnd (hxx L88-89).
    my_edge_start: i32,
    my_edge_end: i32,
    /// OCCT myFaceStart / myFaceEnd (hxx L90-91).
    my_face_start: i32,
    my_face_end: i32,
    /// OCCT myMinMax (hxx L92).
    my_min_max: MinMaxIndices,
}

impl ShapeBounds {
    /// OCCT HLRBRep_ShapeBounds(S, nbIso, V1, V2, E1, E2, F1, F2)
    /// (cxx L46-63) — myShapeData stays the null handle.
    #[allow(clippy::too_many_arguments)]
    pub fn new(s: Arc<OutLiner>, nb_iso: i32, v1: i32, v2: i32, e1: i32, e2: i32, f1: i32, f2: i32) -> Self {
        ShapeBounds {
            my_shape: s, // myShape(S)
            my_shape_data: None,
            my_nb_iso: nb_iso, // myNbIso(nbIso)
            my_vert_start: v1,
            my_vert_end: v2,
            my_edge_start: e1,
            my_edge_end: e2,
            my_face_start: f1,
            my_face_end: f2,
            my_min_max: MinMaxIndices::default(),
        }
    }

    /// OCCT HLRBRep_ShapeBounds(S, SData, nbIso, V1, V2, E1, E2, F1, F2)
    /// (cxx L23-42).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_data(
        s: Arc<OutLiner>,
        s_data: SDataHandle,
        nb_iso: i32,
        v1: i32,
        v2: i32,
        e1: i32,
        e2: i32,
        f1: i32,
        f2: i32,
    ) -> Self {
        ShapeBounds {
            my_shape: s,            // myShape(S)
            my_shape_data: s_data,  // myShapeData(SData)
            my_nb_iso: nb_iso,      // myNbIso(nbIso)
            my_vert_start: v1,      // myVertStart(V1)
            my_vert_end: v2,        // myVertEnd(V2)
            my_edge_start: e1,      // myEdgeStart(E1)
            my_edge_end: e2,        // myEdgeEnd(E2)
            my_face_start: f1,      // myFaceStart(F1)
            my_face_end: f2,        // myFaceEnd(F2)
            my_min_max: MinMaxIndices::default(),
        }
    }

    /// OCCT HLRBRep_ShapeBounds() (lxx L19-28) — every bound and NbOfIso
    /// zero, the shape handle null (the empty OutLiner stands in).
    pub fn default_ctor() -> Self {
        ShapeBounds {
            my_shape: Arc::new(OutLiner::new()),
            my_shape_data: None,
            my_nb_iso: 0, // myNbIso(0)
            my_vert_start: 0,
            my_vert_end: 0,
            my_edge_start: 0,
            my_edge_end: 0,
            my_face_start: 0,
            my_face_end: 0,
            my_min_max: MinMaxIndices::default(),
        }
    }

    /// OCCT HLRBRep_ShapeBounds::Translate (cxx L67-75).
    pub fn translate(&mut self, nv: i32, ne: i32, nf: i32) {
        self.my_vert_start += nv;
        self.my_vert_end += nv;
        self.my_edge_start += ne;
        self.my_edge_end += ne;
        self.my_face_start += nf;
        self.my_face_end += nf;
    }

    /// OCCT Shape(const handle&) (lxx L32-35).
    pub fn set_shape(&mut self, s: &Arc<OutLiner>) {
        self.my_shape = s.clone();
    }

    /// OCCT const handle& Shape() (lxx L39-42).
    pub fn shape(&self) -> &Arc<OutLiner> {
        &self.my_shape
    }

    /// OCCT ShapeData(const handle&) (lxx L46-49).
    pub fn set_shape_data(&mut self, sd: SDataHandle) {
        self.my_shape_data = sd;
    }

    /// OCCT const handle& ShapeData() (lxx L53-56) — the null handle is
    /// None.
    pub fn shape_data(&self) -> &SDataHandle {
        &self.my_shape_data
    }

    /// OCCT NbOfIso(const int) (lxx L60-63).
    pub fn set_nb_of_iso(&mut self, nb_iso: i32) {
        self.my_nb_iso = nb_iso;
    }

    /// OCCT int NbOfIso() (lxx L67-70).
    pub fn nb_of_iso(&self) -> i32 {
        self.my_nb_iso
    }

    /// OCCT HLRBRep_ShapeBounds::Sizes (cxx L79-84) — the OCCT `int& NV,
    /// int& NE, int& NF` out-parameters as the returned tuple.
    pub fn sizes(&self) -> (i32, i32, i32) {
        (
            self.my_vert_end + 1 - self.my_vert_start,
            self.my_edge_end + 1 - self.my_edge_start,
            self.my_face_end + 1 - self.my_face_start,
        )
    }

    /// OCCT HLRBRep_ShapeBounds::Bounds (cxx L88-96) — the OCCT `int&`
    /// out-parameters as the returned tuple.
    pub fn bounds(&self) -> (i32, i32, i32, i32, i32, i32) {
        (
            self.my_vert_start,
            self.my_vert_end,
            self.my_edge_start,
            self.my_edge_end,
            self.my_face_start,
            self.my_face_end,
        )
    }

    /// OCCT UpdateMinMax (hxx L75-78).
    pub fn update_min_max(&mut self, the_tot_min_max: &MinMaxIndices) {
        self.my_min_max = *the_tot_min_max;
    }

    /// OCCT MinMax (hxx L80).
    pub fn min_max(&mut self) -> &mut MinMaxIndices {
        &mut self.my_min_max
    }
}

impl Default for ShapeBounds {
    /// The lxx L19-28 default constructor.
    fn default() -> Self {
        Self::default_ctor()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An empty OutLiner (the null-shape stand-in of the default handle).
    fn empty_outliner() -> Arc<OutLiner> {
        Arc::new(OutLiner::new())
    }

    /// OCCT anchor: the default constructor (lxx L19-28) — all zeros.
    #[test]
    fn shape_bounds_default_ctor() {
        let sb = ShapeBounds::default_ctor();
        assert_eq!(sb.nb_of_iso(), 0);
        assert!(sb.shape_data().is_none());
        assert_eq!(sb.sizes(), (1, 1, 1)); // end + 1 - start
        assert_eq!(sb.bounds(), (0, 0, 0, 0, 0, 0));
    }

    /// OCCT anchor: the two value constructors (cxx L23-42 / L46-63), the
    /// lxx accessor round-trips (L32-70) and Translate (cxx L67-75).
    #[test]
    fn shape_bounds_ctors_translate_accessors() {
        let s = empty_outliner();

        // cxx L46-63: the constructor without SData.
        let mut sb = ShapeBounds::new(s.clone(), 3, 1, 4, 1, 8, 1, 2);
        assert!(sb.shape_data().is_none());
        assert!(Arc::ptr_eq(sb.shape(), &s));
        assert_eq!(sb.nb_of_iso(), 3);

        // cxx L23-42: the constructor with SData — the opaque
        // Standard_Transient payload.
        let s_data: SDataHandle = Some(Arc::new(7u32));
        let sb2 = ShapeBounds::new_with_data(s.clone(), s_data.clone(), 5, 2, 5, 3, 10, 2, 3);
        assert!(sb2.shape_data().is_some());
        // the same Transient handle (Arc pointer identity).
        assert!(Arc::ptr_eq(
            sb2.shape_data().as_ref().unwrap(),
            s_data.as_ref().unwrap()
        ));
        assert_eq!(sb2.nb_of_iso(), 5);

        // Sizes (cxx L79-84): end + 1 - start.
        assert_eq!(sb2.sizes(), (4, 8, 2));

        // Bounds (cxx L88-96): the constructor bounds round-trip.
        assert_eq!(sb2.bounds(), (2, 5, 3, 10, 2, 3));

        // Translate (cxx L67-75): the bounds shift, the sizes do not.
        sb.translate(10, 20, 30);
        assert_eq!(sb.bounds(), (11, 14, 21, 28, 31, 32));
        assert_eq!(sb.sizes(), (4, 8, 2));

        // the lxx setters: Shape / ShapeData / NbOfIso.
        let s2 = empty_outliner();
        sb.set_shape(&s2);
        assert!(Arc::ptr_eq(sb.shape(), &s2));
        sb.set_shape_data(s_data.clone());
        assert!(sb.shape_data().is_some());
        sb.set_nb_of_iso(7);
        assert_eq!(sb.nb_of_iso(), 7);

        // the OCCT handle copy — Clone.
        let sb_clone = sb.clone();
        assert!(Arc::ptr_eq(sb_clone.shape(), sb.shape()));
        assert_eq!(sb_clone.bounds(), sb.bounds());

        // UpdateMinMax / MinMax (hxx L75-80).
        let mut mm = MinMaxIndices::default();
        mm.min = [1; 8];
        mm.max = [2; 8];
        sb.update_min_max(&mm);
        assert_eq!(*sb.min_max(), mm);
        sb.min_max().min[0] = 5;
        assert_eq!(sb.min_max().min[0], 5);
    }

    /// OCCT anchor: the InternalAlgo::Update L115 re-Bind pattern — a
    /// ShapeBounds rebound over a real OutLiner with the loaded bounds
    /// (1..dv / 1..de / 1..df).
    #[test]
    fn shape_bounds_rebind_pattern() {
        let s = empty_outliner();

        // the OCCT L115 sequence: SB = ShapeBounds(SB.Shape(), SB.ShapeData(),
        // SB.NbOfIso(), 1, dv, 1, de, 1, df).
        let dv = 8;
        let de = 12;
        let df = 2;
        let sb0 = ShapeBounds::new(s.clone(), 0, 0, 0, 0, 0, 0, 0);
        let s_data = sb0.shape_data().clone();
        let nb_iso = sb0.nb_of_iso();
        let mut sb = ShapeBounds::new_with_data(sb0.shape().clone(), s_data, nb_iso, 1, dv, 1, de, 1, df);

        assert_eq!(sb.sizes(), (dv, de, df));
        assert_eq!(sb.bounds(), (1, dv, 1, de, 1, df));

        // the L125-138 merge loop: Sizes + Translate.
        let (dv2, de2, df2) = sb.sizes();
        sb.translate(0, 0, 0);
        let _ = (dv2, de2, df2);
        assert_eq!(sb.bounds(), (1, dv, 1, de, 1, df));
    }
}
