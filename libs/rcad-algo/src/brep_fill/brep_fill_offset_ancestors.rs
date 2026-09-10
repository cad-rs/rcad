//! OCCT BRepFill_OffsetAncestors — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_OffsetAncestors.hxx (L26-56) +
//!         BRepFill_OffsetAncestors.cxx (L27-94).
//!
//! First consumer: BRepFill_Evolved::PlanarPerform / VerticalPerform (the
//! evolved batch); earlier consumer BRepFill_OffsetWire paths.
//!
//! Architecture difference: `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape,
//! TopTools_ShapeMapHasher>` maps to a Vec of (key, item) pairs keyed by the
//! shape identity (the offset_wire.rs reduction; IsBound = a linear scan).

use rcad_kernel::topo_shape::Shape;

use super::offset_wire::BRepFillOffsetWire;
use crate::brep_fill::generator::{shape_key, ShapeKey};

/// OCCT BRepFill_OffsetAncestors (hxx L26-56) — finds the generating
/// shapes of an OffsetWire.
#[derive(Debug)]
pub struct BRepFillOffsetAncestors {
    my_is_perform: bool,                    // OCCT: myIsPerform
    my_map: Vec<(ShapeKey, Shape)>,         // OCCT: myMap
}

impl Default for BRepFillOffsetAncestors {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFillOffsetAncestors {
    /// OCCT BRepFill_OffsetAncestors::BRepFill_OffsetAncestors() (cxx L29-33).
    pub fn new() -> Self {
        BRepFillOffsetAncestors {
            my_is_perform: false,
            my_map: Vec::new(),
        }
    }

    /// OCCT BRepFill_OffsetAncestors::BRepFill_OffsetAncestors(Paral)
    /// (cxx L35-39).
    pub fn new_with_paral(paral: &mut BRepFillOffsetWire, brep: &rcad_kernel::topods::BRep) -> Self {
        let mut r = BRepFillOffsetAncestors::new();
        r.perform(paral, brep);
        r
    }

    /// OCCT BRepFill_OffsetAncestors::Perform(Paral) (cxx L41-70).
    pub fn perform(&mut self, paral: &mut BRepFillOffsetWire, brep: &rcad_kernel::topods::BRep) {
        let spine = paral.spine().clone();

        // on itere sur les edges.
        for e in crate::brep_algo::tool::explorer(&spine, rcad_kernel::topods::ShapeType::Edge, rcad_kernel::topods::ShapeType::Shape)
        {
            for it in paral.generated_shapes(brep, &e) {
                self.my_map_bind(&it, &e);
            }
        }

        // on itere sur les vertex.
        for v in crate::brep_algo::tool::explorer(&spine, rcad_kernel::topods::ShapeType::Vertex, rcad_kernel::topods::ShapeType::Shape)
        {
            for it in paral.generated_shapes(brep, &v) {
                self.my_map_bind(&it, &v);
            }
        }

        self.my_is_perform = true;
    }

    /// OCCT NCollection_DataMap::Bind(it.Value(), Exp.Current()).
    fn my_map_bind(&mut self, key: &Shape, value: &Shape) {
        let k = shape_key(key);
        for (kk, vv) in self.my_map.iter_mut() {
            if *kk == k {
                *vv = value.clone();
                return;
            }
        }
        self.my_map.push((k, value.clone()));
    }

    /// OCCT NCollection_DataMap::IsBound(S1).
    fn my_map_is_bound(&self, key: &Shape) -> bool {
        let k = shape_key(key);
        self.my_map.iter().any(|(kk, _)| *kk == k)
    }

    /// OCCT NCollection_DataMap::operator()(S1).
    fn my_map_find(&self, key: &Shape) -> Shape {
        let k = shape_key(key);
        self.my_map
            .iter()
            .find(|(kk, _)| *kk == k)
            .map(|(_, v)| v.clone())
            .expect("myMap(S1)")
    }

    /// OCCT BRepFill_OffsetAncestors::IsDone() (cxx L72-76).
    pub fn is_done(&self) -> bool {
        self.my_is_perform
    }

    /// OCCT BRepFill_OffsetAncestors::HasAncestor(S1) (cxx L78-82).
    pub fn has_ancestor(&self, s1: &Shape) -> bool {
        self.my_map_is_bound(s1)
    }

    /// OCCT BRepFill_OffsetAncestors::Ancestor(S1) (cxx L84-94) — may
    /// return a Null Shape if S1 is not a subShape of Paral.
    pub fn ancestor(&self, s1: &Shape) -> Shape {
        if !self.my_is_perform {
            panic!(
                "StdFail_NotDone: BRepFill_OffsetAncestors::Ancestor() - Perform() \
                 should be called before accessing results"
            );
        }
        self.my_map_find(s1)
    }
}
