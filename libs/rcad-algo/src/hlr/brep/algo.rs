// OCCT HLRBRep_Algo (TKHLR/HLRBRep/HLRBRep_Algo.hxx L1-94 + .cxx L1-87) —
// the user-facing entry (the InternalAlgo subclass providing the methods
// with Shape from TopoDS).
//
// Rust encoding of the OCCT inheritance: composition + forwarding —
// [`Algo`] owns the [`InternalAlgo`] and `Deref`/`DerefMut` forward the
// base-class members (the O2 landed [`super::internal_algo::InternalAlgo`]);
// the four HLRBRep_Algo methods are explicit (Add x2 / Index /
// OutLinedShapeNullify).
//
// Deviation (reported): the OCCT `occ::handle<Standard_Transient> SData`
// slot maps to the O2 [`super::internal_algo::SDataHandle`].

use std::ops::{Deref, DerefMut};
use std::sync::Arc;

use rcad_kernel::topods::Shape;

use crate::hlr::topo_brep::out_liner::OutLiner;

use super::internal_algo::{InternalAlgo, SDataHandle};

/// OCCT HLRBRep_Algo (hxx L59-92) — a framework to compute a shape as seen
/// in a projection plane.
pub struct Algo {
    /// the HLRBRep_InternalAlgo base (the composition form of the OCCT
    /// inheritance).
    pub internal: InternalAlgo,
}

impl Algo {
    /// OCCT HLRBRep_Algo() (cxx L30) — constructs an empty framework for
    /// the calculation of visible and hidden lines of a shape in a
    /// projection.
    pub fn new() -> Algo {
        Algo {
            internal: InternalAlgo::new(),
        }
    }

    /// OCCT HLRBRep_Algo(const occ::handle<HLRBRep_Algo>& A) (cxx L34-37) —
    /// the base copy (see the InternalAlgo::new_from DS deviation note).
    pub fn new_from(a: &Algo) -> Algo {
        Algo {
            internal: InternalAlgo::new_from(&a.internal),
        }
    }

    /// OCCT Add(const TopoDS_Shape& S, const occ::handle<Standard_Transient>&
    /// SData, const int nbIso) (cxx L41-46) — add the Shape <S>.
    pub fn add_with_data(&mut self, s: &Shape, s_data: SDataHandle, nb_iso: i32) {
        // Load(new HLRTopoBRep_OutLiner(S), SData, nbIso);
        self.internal
            .load_with_data(&Arc::new(OutLiner::with_original_shape(s)), s_data, nb_iso);
    }

    /// OCCT Add(const TopoDS_Shape& S, const int nbIso) (cxx L50-53) — adds
    /// the shape S to this framework, and specifies the number of
    /// isoparameters nbiso desired in visualizing S.
    pub fn add(&mut self, s: &Shape, nb_iso: i32) {
        // Load(new HLRTopoBRep_OutLiner(S), nbIso);
        self.internal
            .load(&Arc::new(OutLiner::with_original_shape(s)), nb_iso);
    }

    /// OCCT Index(const TopoDS_Shape& S) (cxx L57-74) — return the index of
    /// the Shape <S> and return 0 if the Shape <S> is not found.
    pub fn index(&mut self, s: &Shape) -> usize {
        // int n = NbShapes();
        let n = self.internal.nb_shapes();

        // for (int i = 1; i <= n; i++)
        for i in 1..=n {
            // the OutLiner sits behind the ShapeBounds Arc handle; the
            // non-const OriginalShape()/OutLinedShape() getters go through
            // the raw-pointer deref (the FaceData Wires handle precedent).
            let sb = self.internal.shape_bounds(i);
            let ol: *mut OutLiner = Arc::as_ptr(sb.shape()) as *const OutLiner as *mut OutLiner;
            // if (ShapeBounds(i).Shape()->OriginalShape() == S)
            if *unsafe { (*ol).original_shape() } == *s {
                return i;
            }
            // if (ShapeBounds(i).Shape()->OutLinedShape() == S)
            if *unsafe { (*ol).out_lined_shape() } == *s {
                return i;
            }
        }

        // return 0;
        0
    }

    /// OCCT OutLinedShapeNullify() (cxx L78-87) — nullify all the results of
    /// OutLiner from HLRTopoBRep.
    pub fn out_lined_shape_nullify(&mut self) {
        // int n = NbShapes();
        let n = self.internal.nb_shapes();

        // for (int i = 1; i <= n; i++)
        for i in 1..=n {
            let sb = self.internal.shape_bounds(i);
            let ol: *mut OutLiner = Arc::as_ptr(sb.shape()) as *const OutLiner as *mut OutLiner;
            unsafe {
                // ShapeBounds(i).Shape()->OutLinedShape(TopoDS_Shape());
                (*ol).set_out_lined_shape(&Shape::null());
                // ShapeBounds(i).Shape()->DataStructure().Clear();
                (*ol).data_structure().clear();
            }
        }
    }
}

impl Default for Algo {
    /// OCCT HLRBRep_Algo() — the default-framework form.
    fn default() -> Self {
        Self::new()
    }
}

// the OCCT inheritance: the HLRBRep_InternalAlgo members are reachable on
// the HLRBRep_Algo object.
impl Deref for Algo {
    type Target = InternalAlgo;
    fn deref(&self) -> &InternalAlgo {
        &self.internal
    }
}

impl DerefMut for Algo {
    fn deref_mut(&mut self) -> &mut InternalAlgo {
        &mut self.internal
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::brep::data::Data;
    use glam::DVec3;
    use rcad_kernel::topo::topods::BRepBuilder;

    /// A face shape fixture (the smallest TopoDS_Shape stand-in for the
    /// OutLiner original shape).
    fn face_shape() -> (std::sync::Arc<rcad_kernel::BRep>, Shape) {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let va = b.add_vertex(&mut brep, DVec3::new(0.0, 0.0, 0.0), 1e-7);
        let vb = b.add_vertex(&mut brep, DVec3::new(1.0, 0.0, 0.0), 1e-7);
        let vc = b.add_vertex(&mut brep, DVec3::new(1.0, 1.0, 0.0), 1e-7);
        let e1 = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(
                rcad_kernel::geom::Line3 {
                    origin: DVec3::new(0.0, 0.0, 0.0),
                    direction: DVec3::new(1.0, 0.0, 0.0),
                },
            )),
            va.clone(),
            vb.clone(),
            [0.0, 1.0],
        );
        let e2 = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(
                rcad_kernel::geom::Line3 {
                    origin: DVec3::new(1.0, 0.0, 0.0),
                    direction: DVec3::new(0.0, 1.0, 0.0),
                },
            )),
            vb.clone(),
            vc.clone(),
            [0.0, 1.0],
        );
        let e3 = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(
                rcad_kernel::geom::Line3 {
                    origin: DVec3::new(1.0, 1.0, 0.0),
                    direction: DVec3::new(-1.0, 0.0, 0.0),
                },
            )),
            vc.clone(),
            va.clone(),
            [0.0, 1.0],
        );
        let wire = brep.add_twire(vec![e1, e2, e3]);
        let face = brep.add_tface(
            None,
            wire,
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );
        (std::sync::Arc::new(brep), face)
    }

    /// OCCT Add(S, nbIso) / Index(S) (cxx L50-53, L57-74) — the added shape
    /// is found by its original shape and an unknown shape returns 0.
    #[test]
    fn add_and_index_original_and_zero_for_unknown() {
        let (_brep, face) = face_shape();
        let mut algo = Algo::new();
        assert_eq!(algo.internal.nb_shapes(), 0);
        algo.add(&face, 0);
        assert_eq!(algo.internal.nb_shapes(), 1);
        // the original shape is found at index 1.
        assert_eq!(algo.index(&face), 1);
        // an unknown shape returns 0 (the OCCT not-found convention).
        let (_brep2, other) = face_shape();
        assert_eq!(algo.index(&other), 0);
    }

    /// OCCT Add(S, SData, nbIso) (cxx L41-46) — the SData slot rides into
    /// the ShapeBounds.
    #[test]
    fn add_with_data_forwards_sdata() {
        let (_brep, face) = face_shape();
        let mut algo = Algo::new();
        algo.add_with_data(&face, Some(Arc::new(7_u32)), 3);
        assert_eq!(algo.internal.nb_shapes(), 1);
        assert_eq!(algo.index(&face), 1);
        let sb = algo.internal.shape_bounds(1);
        assert_eq!(sb.nb_of_iso(), 3);
        assert!(sb.shape_data().is_some());
    }

    /// OCCT OutLinedShapeNullify (cxx L78-87) — the outlined shape is
    /// nulled and the HLRTopoBRep data structure cleared for every bound
    /// shape.
    #[test]
    fn out_lined_shape_nullify_clears_all() {
        let (_brep, face) = face_shape();
        let mut algo = Algo::new();
        algo.add(&face, 0);
        algo.add(&face, 0);
        // set a non-null outlined shape + a dirty data structure on both.
        for i in 1..=2 {
            let sb = algo.internal.shape_bounds(i);
            let ol: *mut OutLiner =
                Arc::as_ptr(sb.shape()) as *const OutLiner as *mut OutLiner;
            unsafe {
                (*ol).set_out_lined_shape(&face);
                (*ol).data_structure().add_spl_e(&face);
            }
        }
        algo.out_lined_shape_nullify();
        for i in 1..=2 {
            let sb = algo.internal.shape_bounds(i);
            let ol: *mut OutLiner =
                Arc::as_ptr(sb.shape()) as *const OutLiner as *mut OutLiner;
            unsafe {
                assert!( (*ol).out_lined_shape().is_null());
                assert!(!(*ol).data_structure().edge_has_spl_e(&face));
            }
        }
    }

    /// The composition forwardings: new_from copies the shapes (cxx L34-37)
    /// and the Deref reaches the InternalAlgo members (the inheritance
    /// encoding).
    #[test]
    fn new_from_and_deref_forwarding() {
        let (_brep, face) = face_shape();
        let mut algo = Algo::new();
        algo.add(&face, 2);
        let copy = Algo::new_from(&algo);
        assert_eq!(copy.internal.nb_shapes(), 1);
        // the Deref forwarding reaches NbShapes on the base.
        assert_eq!(algo.nb_shapes(), 1);
        // the DS of the copy is null until the next Update (the new_from
        // deviation note of InternalAlgo).
        assert!(copy.data_structure().is_none());
        // a scratch DS proves the composition storage is the test subject
        // of the Data type (the InstallDsForTest path of O2).
        let _ = Data::new(0, 0, 0);
    }
}
