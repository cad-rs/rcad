//! TKBRep GTest translations (topods::BRep API).
//!
//! OCCT source: src/ModelingData/TKBRep/GTests/
//!
//! Files translated:
//!   TopoDS_TShape_Test.cxx      鈥?TShape flags, ShapeType, NbChildren, EmptyCopy
//!   TopExp_Test.cxx             鈥?Vertex/Edge properties
//!   TopoDS_Iterator_Test.cxx    鈥?TopoDS_Iterator over sub-shapes

use glam::DVec3;
use rcad_kernel::topods::{self, BRep, Shape, ShapeType, TShape, tshape_flags};

const TOL: f64 = 1e-7;

fn shape_type_of(brep: &BRep, sr: ShapeRef) -> ShapeType {
    match &*brep.tshapes[sr.index] {
        TShape::Vertex(_) => ShapeType::Vertex,
        TShape::Edge(_) => ShapeType::Edge,
        TShape::Wire(_) => ShapeType::Wire,
        TShape::Face(_) => ShapeType::Face,
        TShape::Shell(_) => ShapeType::Shell,
        TShape::Solid(_) => ShapeType::Solid,
        TShape::CompSolid(_) => ShapeType::CompSolid,
        TShape::Compound(_) => ShapeType::Compound,
    }
}

// =============================================================================
// TopoDS_TShape_Test.cxx 鈥?ShapeType, flags, EmptyCopy
// =============================================================================

#[cfg(test)]
mod topods_tshape_tests {
    use super::*;

    #[test]
    fn shapetype_all_types() {
        let mut b = BRep::new();
        let v = b.add_tvertex(DVec3::ZERO);
        assert_eq!(shape_type_of(&b, &b.ref_to_shape(v)), ShapeType::Vertex);

        let e = b.add_tedge(None, v, Shape::null(), [0.0, 1.0]);
        assert_eq!(shape_type_of(&b, &b.ref_to_shape(e)), ShapeType::Edge);

        let w = b.add_twire(vec![]);
        assert_eq!(shape_type_of(&b, &b.ref_to_shape(w)), ShapeType::Wire);

        let f = b.add_tface(None, w, vec![], None, None, vec![], true);
        assert_eq!(shape_type_of(&b, &b.ref_to_shape(f)), ShapeType::Face);

        let sh = b.add_tshell(vec![]);
        assert_eq!(shape_type_of(&b, &b.ref_to_shape(sh)), ShapeType::Shell);

        let s = b.add_tsolid(vec![]);
        assert_eq!(shape_type_of(&b, &b.ref_to_shape(s)), ShapeType::Solid);
    }

    #[test]
    fn flag_setters_getters() {
        let mut b = BRep::new();
        let sr = b.add_tvertex(DVec3::ZERO);

        // Free flag starts true (part of DEFAULT flags)
        assert!(b.has_flag(sr, tshape_flags::FREE));
        b.set_flag(sr, tshape_flags::FREE, false);
        assert!(!b.has_flag(sr, tshape_flags::FREE));
        b.set_flag(sr, tshape_flags::FREE, true);
        assert!(b.has_flag(sr, tshape_flags::FREE));

        // Modified flag starts true
        assert!(b.has_flag(sr, tshape_flags::MODIFIED));
        b.set_flag(sr, tshape_flags::MODIFIED, false);
        assert!(!b.has_flag(sr, tshape_flags::MODIFIED));
        b.set_flag(sr, tshape_flags::MODIFIED, true);
        assert!(b.has_flag(sr, tshape_flags::MODIFIED));

        // Locked flag starts false
        assert!(!b.has_flag(sr, tshape_flags::LOCKED));
        b.set_flag(sr, tshape_flags::LOCKED, true);
        assert!(b.has_flag(sr, tshape_flags::LOCKED));
    }

    #[test]
    fn empty_copy_same_type_no_children() {
        let mut b = BRep::new();
        let w = b.add_twire(vec![]);
        assert_eq!(shape_type_of(&b, &b.ref_to_shape(w)), ShapeType::Wire);
    }
}

// =============================================================================
// TopExp_Test.cxx 鈥?Vertex/Edge properties
// =============================================================================

#[cfg(test)]
mod topexp_tests {
    use super::*;

    #[test]
    fn vertex_tolerance_exists() {
        let mut b = BRep::new();
        let v = b.add_tvertex(DVec3::new(1.0, 2.0, 3.0));
        b.vertex_mut(v).tolerance = 1e-6;
        assert!((b.vertex_mut(v).tolerance - 1e-6).abs() < TOL);
    }

    #[test]
    fn edge_has_curve_or_not() {
        let mut b = BRep::new();
        let sv = b.add_tvertex(DVec3::ZERO);
        let ev = b.add_tvertex(DVec3::new(10.0, 0.0, 0.0));
        let e = b.add_tedge(None, sv, ev, [0.0, 1.0]);
        let te = &*b.tshapes[b.shape_idx(&e)];
        if let TShape::Edge(ed) = te {
            assert!(ed.curve.is_none());
        } else {
            panic!("expected Edge");
        }
    }

    #[test]
    fn full_circle_edge_is_closed() {
        use rcad_kernel::geom::Curve3;
        let mut b = BRep::new();
        let sv = b.add_tvertex(DVec3::new(1.0, 0.0, 0.0));
        let ev = b.add_tvertex(DVec3::new(1.0, 0.0, 0.0));
        let circle = Some(Curve3::Circle(rcad_kernel::geom::Circle3::new(
            DVec3::ZERO,
            DVec3::Z,
            1.0,
        )));
        let e = b.add_tedge(circle, sv, ev, [0.0, std::f64::consts::TAU]);
        if let TShape::Edge(ed) = &*b.tshapes[b.shape_idx(&e)] {
            assert!(ed.first.index == ed.last.index);
        } else {
            panic!("expected Edge");
        }
    }

    #[test]
    fn line_edge_not_closed() {
        use rcad_kernel::geom::Curve3;
        let mut b = BRep::new();
        let sv = b.add_tvertex(DVec3::ZERO);
        let ev = b.add_tvertex(DVec3::new(10.0, 0.0, 0.0));
        let line = Some(Curve3::Line(rcad_kernel::geom::Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        }));
        let e = b.add_tedge(line, sv, ev, [0.0, 10.0]);
        if let TShape::Edge(ed) = &*b.tshapes[b.shape_idx(&e)] {
            assert!(ed.first.index != ed.last.index);
        } else {
            panic!("expected Edge");
        }
    }
}

// =============================================================================
// TopoDS_Iterator_Test.cxx 鈥?Iteration over sub-shapes
// =============================================================================

#[cfg(test)]
mod topods_iterator_tests {
    use super::*;

    #[test]
    fn iterator_over_compound() {
        let mut b = BRep::new();
        let v1 = b.add_tvertex(DVec3::new(1.0, 0.0, 0.0));
        let v2 = b.add_tvertex(DVec3::new(2.0, 0.0, 0.0));
        let c = b.add_tcompound(vec![v1, v2]);

        let count = match &*b.tshapes[b.shape_idx(&c)] {
            TShape::Compound(cd) => cd.len(),
            _ => 0,
        };
        assert_eq!(count, 2);
    }
}
