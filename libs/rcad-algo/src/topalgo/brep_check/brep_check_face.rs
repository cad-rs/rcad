//! OCCT BRepCheck_Face (TKTopAlgo/BRepCheck).
//!
//! Source: `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Face.cxx`
//! (L59-954) and `BRepCheck_Face.hxx` (L31-95).
//!
//! `Intersect` (Face.cxx L626-798) runs the general 2D curve/curve
//! intersector `Geom2dInt_GInter` — the rcad real body is
//! `crate::geomalgo::geom2d_int::GInter` (the `IntCurve_IntCurveCurveGen`
//! instantiation).

use std::collections::HashMap;

use rcad_kernel::geom::{Curve2d, Curve2dEval, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{tshape_flags, BRep, BRepTool, Orientation, ShapeType, State, TShape, TFaceData};
use std::sync::Arc;

use crate::geomalgo::geom2d_int::GInter;
use crate::geomalgo::int_res2d::Domain as Res2dDomain;
use crate::topalgo::brep_class::bnd_lib_add2d_curve::add_2d_curve;
use crate::topalgo::brep_top_adaptor::fclass2d::FClass2d;
use crate::topalgo::shape_source::FaceShapeSource;

use super::brep_check_result::{
    brep_check_add, explorer, oriented, ShapeKey, BRepCheckResultBase, BRepCheckStatus,
};
use super::brep_check_wire::uv_points;

/// OCCT Face.cxx L59: `DataMapOfShapeBox2d`.
type DataMapOfShapeBox2d = HashMap<ShapeKey, (BndBox2d, Shape)>;

/// Insertion-ordered `NCollection_DataMap<TopoDS_Shape, List<TopoDS_Shape>>`
/// (myMapImb).
#[derive(Debug, Default)]
pub struct ImbricationMap {
    entries: Vec<(Shape, Vec<Shape>)>,
    index: HashMap<ShapeKey, usize>,
}

impl ImbricationMap {
    pub fn new() -> Self {
        ImbricationMap::default()
    }

    /// OCCT IsBound.
    pub fn is_bound(&self, s: &Shape) -> bool {
        self.index.contains_key(&ShapeKey::of(s))
    }

    /// OCCT Bind.
    pub fn bind(&mut self, s: &Shape, list: Vec<Shape>) {
        let key = ShapeKey::of(s);
        if let Some(&i) = self.index.get(&key) {
            self.entries[i].1 = list;
        } else {
            self.index.insert(key, self.entries.len());
            self.entries.push((s.clone(), list));
        }
    }

    /// OCCT UnBind.
    pub fn un_bind(&mut self, s: &Shape) {
        let key = ShapeKey::of(s);
        if let Some(i) = self.index.remove(&key) {
            self.entries.remove(i);
            // Reindex the shifted entries.
            for (j, (k, _)) in self.entries.iter().enumerate() {
                self.index.insert(ShapeKey::of(k), j);
            }
        }
    }

    /// OCCT operator() — the list of `s`.
    pub fn value_mut(&mut self, s: &Shape) -> Option<&mut Vec<Shape>> {
        match self.index.get(&ShapeKey::of(s)) {
            Some(&i) => Some(&mut self.entries[i].1),
            None => None,
        }
    }

    pub fn value(&self, s: &Shape) -> Option<&Vec<Shape>> {
        self.index
            .get(&ShapeKey::of(s))
            .map(|&i| &self.entries[i].1)
    }

    pub fn extent(&self) -> usize {
        self.entries.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&Shape, &Vec<Shape>)> {
        self.entries.iter().map(|(s, l)| (s, l))
    }

    pub fn clear(&mut self) {
        self.entries.clear();
        self.index.clear();
    }
}

/// OCCT BRepCheck_Face (Face.hxx L31-95).
#[derive(Debug)]
pub struct BRepCheckFace {
    /// OCCT protected base (Result.hxx L82-90).
    pub base: BRepCheckResultBase,
    /// OCCT myIntdone.
    pub my_intdone: bool,
    /// OCCT myIntres.
    pub my_intres: BRepCheckStatus,
    /// OCCT myImbdone.
    pub my_imbdone: bool,
    /// OCCT myImbres.
    pub my_imbres: BRepCheckStatus,
    /// OCCT myOridone.
    pub my_oridone: bool,
    /// OCCT myOrires.
    pub my_orires: BRepCheckStatus,
    /// OCCT myGctrl.
    pub my_gctrl: bool,
    /// OCCT myMapImb (Face.hxx private member).
    pub my_map_imb: ImbricationMap,
}

impl BRepCheckFace {
    /// OCCT BRepCheck_Face::BRepCheck_Face(const TopoDS_Face& F)
    /// (Face.cxx L75-85).
    pub fn new(brep: &BRep, f: &Shape) -> Self {
        let mut r = BRepCheckFace {
            base: BRepCheckResultBase::new(),
            my_intdone: false,
            my_intres: BRepCheckStatus::NoError,
            my_imbdone: false,
            my_imbres: BRepCheckStatus::NoError,
            my_oridone: false,
            my_orires: BRepCheckStatus::NoError,
            my_gctrl: true,
            my_map_imb: ImbricationMap::new(),
        };
        r.base.init(f);
        // OCCT L81-84.
        r.my_intdone = false;
        r.my_imbdone = false;
        r.my_oridone = false;
        r.my_gctrl = true;
        r.minimum(brep);
        r
    }

    /// OCCT BRepCheck_Face::Minimum (Face.cxx L89-112).
    pub fn minimum(&mut self, _brep: &BRep) {
        if !self.base.my_min {
            // OCCT L93-95.
            self.base.my_map.bound(&self.base.my_shape);
            let my_shape = self.base.my_shape.clone();
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");

            // OCCT L97-105: the face surface must exist.
            let has_surface = my_shape
                .as_face()
                .map(|fd| fd.surface.is_some())
                .unwrap_or(false);
            if !has_surface {
                brep_check_add(lst, BRepCheckStatus::NoSurface);
            } else {
                // Flag natural restriction???
            }
            // OCCT L106-110.
            if lst.is_empty() {
                lst.push(BRepCheckStatus::NoError);
            }
            self.base.my_min = true;
        }
    }

    /// OCCT BRepCheck_Face::InContext (Face.cxx L116-154).
    pub fn in_context(&mut self, brep: &BRep, s: &Shape) {
        // OCCT L118-133: bound check under the (parallel) lock.
        if self.base.my_map.is_bound(s) {
            return;
        }
        self.base.my_map.bind(s.clone(), Vec::new());
        let lst = self
            .base
            .my_map
            .find_mut(s)
            .expect("InContext: the context must be bound");

        // OCCT L136-148.
        let mut found = false;
        let mut more = false;
        for cur in explorer(brep, s, ShapeType::Face) {
            more = true;
            if cur.is_same(&self.base.my_shape) {
                found = true;
                break;
            }
        }
        if !more || !found {
            brep_check_add(lst, BRepCheckStatus::SubshapeNotInShape);
            return;
        }

        // OCCT L150-153.
        if lst.is_empty() {
            lst.push(BRepCheckStatus::NoError);
        }
    }

    /// OCCT BRepCheck_Face::Blind (Face.cxx L158-165).
    pub fn blind(&mut self) {
        if !self.base.my_blind {
            // nothing more than in the minimum
            self.base.my_blind = true;
        }
    }

    /// OCCT BRepCheck_Face::IntersectWires (Face.cxx L169-298).
    pub fn intersect_wires(&mut self, brep: &BRep, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        if self.my_intdone {
            // OCCT L182-189.
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("IntersectWires: myShape must be bound");
                brep_check_add(lst, self.my_intres);
            }
            return self.my_intres;
        }

        // OCCT L191-192.
        self.my_intdone = true;
        self.my_intres = BRepCheckStatus::NoError;
        // This method has to be called by an analyzer. It is assumed that
        // each edge has a correct 2d representation on the face.

        // OCCT L199-217: the wires are mapped.
        let fwd_face = oriented(&my_shape, Orientation::Forward);
        for w in explorer(brep, &fwd_face, ShapeType::Wire) {
            if !self.my_map_imb.is_bound(&w) {
                self.my_map_imb.bind(&w, Vec::new());
            } else {
                // the same wire is met twice...
                self.my_intres = BRepCheckStatus::RedundantWire;
                if update {
                    let lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("IntersectWires: myShape must be bound");
                    brep_check_add(lst, self.my_intres);
                }
                return self.my_intres;
            }
        }

        // OCCT L219-246: create 2d boxes for all edges of all wires.
        let mut a_map_shape_box2d: DataMapOfShapeBox2d = HashMap::new();
        for w in explorer(brep, &my_shape, ShapeType::Wire) {
            let mut a_box_w = BndBox2d::new();
            for e in explorer(brep, &w, ShapeType::Edge) {
                // OCCT L230: aC.Load(BRep_Tool::CurveOnSurface(anEdge, Face,
                // aFirst, aLast)).
                let Some((a_c, mut a_first, mut a_last)) = brep.curve_on_surface(&e, &my_shape)
                else {
                    continue;
                };
                // OCCT L232-239: avoid the exception in Segment for a BSpline.
                let dom = Curve2dEval::default_domain(&a_c);
                if dom[0] > a_first {
                    a_first = dom[0];
                }
                if dom[1] < a_last {
                    a_last = dom[1];
                }
                let mut a_box_e = BndBox2d::new();
                add_2d_curve(&a_c, a_first, a_last, 0., &mut a_box_e);
                a_box_w.add_box(&a_box_e);
                a_map_shape_box2d.insert(ShapeKey::of(&e), (a_box_e, e.clone()));
            }
            a_map_shape_box2d.insert(ShapeKey::of(&w), (a_box_w, w.clone()));
        }

        // OCCT L248-292: the pairwise intersection test.
        let nbwire = self.my_map_imb.extent() as i32;
        let wires_fwd = explorer(brep, &fwd_face, ShapeType::Wire);
        let mut index: i32 = 1;
        while index < nbwire {
            // OCCT L254-262: the Index-th wire.
            let mut indexbis = 0i32;
            let mut wir1: Option<Shape> = None;
            for w in &wires_fwd {
                indexbis += 1;
                if indexbis == index {
                    wir1 = Some(w.clone());
                    break;
                }
            }
            let wir1 = match wir1 {
                Some(w) => w,
                None => break,
            };
            // OCCT L264-268: to reduce the number of Intersect calls.
            let a_box1 = a_map_shape_box2d
                .get(&ShapeKey::of(&wir1))
                .map(|(b, _)| b.clone());
            let mut wir1_iter_pos = wires_fwd
                .iter()
                .position(|w| w.is_equal(&wir1))
                .map(|p| p + 1)
                .unwrap_or(wires_fwd.len());
            while wir1_iter_pos < wires_fwd.len() {
                let wir2 = &wires_fwd[wir1_iter_pos];
                wir1_iter_pos += 1;
                let a_box2 = a_map_shape_box2d
                    .get(&ShapeKey::of(wir2))
                    .map(|(b, _)| b.clone());
                // OCCT L277-280.
                if let (Some(b1), Some(b2)) = (&a_box1, &a_box2) {
                    if !b1.is_void() && !b2.is_void() && b1.is_out_box(b2) {
                        continue;
                    }
                }
                if intersect(brep, &wir1, wir2, &my_shape, &a_map_shape_box2d) {
                    self.my_intres = BRepCheckStatus::IntersectingWires;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("IntersectWires: myShape must be bound");
                        brep_check_add(lst, self.my_intres);
                    }
                    return self.my_intres;
                }
            }
            index += 1;
        }
        // OCCT L293-297.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("IntersectWires: myShape must be bound");
            brep_check_add(lst, self.my_intres);
        }
        self.my_intres
    }

    /// OCCT BRepCheck_Face::ClassifyWires (Face.cxx L302-433).
    pub fn classify_wires(&mut self, brep: &BRep, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        if self.my_imbdone {
            // OCCT L317-324.
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("ClassifyWires: myShape must be bound");
                brep_check_add(lst, self.my_imbres);
            }
            return self.my_imbres;
        }

        // OCCT L326-335.
        self.my_imbdone = true;
        self.my_imbres = self.intersect_wires(brep, false);
        if self.my_imbres != BRepCheckStatus::NoError {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("ClassifyWires: myShape must be bound");
                brep_check_add(lst, self.my_imbres);
            }
            return self.my_imbres;
        }

        // OCCT L337-345.
        let nbwire = self.my_map_imb.extent() as i32;
        if nbwire < 1 {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("ClassifyWires: myShape must be bound");
                brep_check_add(lst, self.my_imbres);
            }
            return self.my_imbres;
        }

        // OCCT L347-383: the classification of each wire in a fresh face.
        let fwd_face = oriented(&my_shape, Orientation::Forward);
        let wires_fwd = explorer(brep, &fwd_face, ShapeType::Wire);
        for wir1 in &wires_fwd {
            // OCCT L354-356: newFace = myShape.EmptyCopied() + Add(wir1).
            let Some(new_face) = empty_copied_face(&my_shape, wir1) else {
                continue;
            };
            let new_face = oriented(&new_face, Orientation::Forward);

            // OCCT L361-369: BRepTopAdaptor_FClass2d over the fresh face.
            // (The fresh face shares the TFace surface of myShape; the world
            // surface is taken from the original — the fresh copy carries no
            // pool index.)
            let Some(surf) = brep.face_surface_world(&my_shape) else {
                continue;
            };
            let Some(fclass2d_ds) = face_shape_source(brep, &new_face, surf) else {
                continue;
            };
            let f_class2d = FClass2d::new(&fclass2d_ds, 0, rcad_kernel::precision::PCONFUSION);
            let mut wire_bien_oriente = false;
            if f_class2d.perform_infinite_point(&fclass2d_ds) != State::Out {
                wire_bien_oriente = true;
                // the given wire defines a hole
                self.my_map_imb.un_bind(wir1);
                let wir1_rev = oriented(wir1, Orientation::Reversed);
                self.my_map_imb.bind(&wir1_rev, Vec::new());
            }

            for wir2 in &wires_fwd {
                if !wir2.is_same(wir1) {
                    // OCCT L377-380.
                    if is_inside(brep, wir2, wire_bien_oriente, &f_class2d, &fclass2d_ds, &new_face)
                    {
                        // OCCT myMapImb(wir1) — the wir1 key was re-bound
                        // under its REVERSED orientation when it defines a
                        // hole; IsSame matches either.
                        let wir1_key = if self.my_map_imb.is_bound(wir1) {
                            wir1.clone()
                        } else {
                            oriented(wir1, Orientation::Reversed)
                        };
                        if let Some(lst) = self.my_map_imb.value_mut(&wir1_key) {
                            lst.push(wir2.clone());
                        }
                    }
                }
            }
        }
        // OCCT L385-410: exactly one wire contains all the others.
        let mut wext: Option<Shape> = None;
        for (key, value) in self.my_map_imb.iter() {
            if !value.is_empty() {
                if wext.is_none() {
                    wext = Some(key.clone());
                } else {
                    self.my_imbres = BRepCheckStatus::InvalidImbricationOfWires;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("ClassifyWires: myShape must be bound");
                        brep_check_add(lst, self.my_imbres);
                    }
                    return self.my_imbres;
                }
            }
        }

        // OCCT L412-424.
        if let Some(wext) = &wext {
            if let Some(list) = self.my_map_imb.value(wext) {
                if list.len() as i32 != nbwire - 1 {
                    self.my_imbres = BRepCheckStatus::InvalidImbricationOfWires;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("ClassifyWires: myShape must be bound");
                        brep_check_add(lst, self.my_imbres);
                    }
                    return self.my_imbres;
                }
            }
        }

        // OCCT L427-432: quit without errors.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("ClassifyWires: myShape must be bound");
            brep_check_add(lst, self.my_imbres);
        }
        self.my_imbres
    }

    /// OCCT BRepCheck_Face::OrientationOfWires (Face.cxx L437-557).
    pub fn orientation_of_wires(&mut self, brep: &BRep, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        // WARNING : it is assumed that the edges of a wire are correctly oriented
        // OCCT L451: bool Infinite = myShape.Infinite().
        let infinite = my_shape
            .as_face()
            .map(|fd| fd.flags & tshape_flags::INFINITE != 0)
            .unwrap_or(false);
        if self.my_oridone {
            // OCCT L452-459.
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("OrientationOfWires: myShape must be bound");
                brep_check_add(lst, self.my_orires);
            }
            return self.my_orires;
        }

        // OCCT L461-470.
        self.my_oridone = true;
        self.my_orires = self.classify_wires(brep, false);
        if self.my_orires != BRepCheckStatus::NoError {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("OrientationOfWires: myShape must be bound");
                brep_check_add(lst, self.my_orires);
            }
            return self.my_orires;
        }

        // OCCT L472-492.
        let nbwire = self.my_map_imb.extent() as i32;
        let mut wext: Option<Shape> = None;
        if nbwire == 1 {
            if !infinite {
                if let Some((key, _)) = self.my_map_imb.iter().next() {
                    wext = Some(key.clone());
                }
            }
        } else {
            for (key, value) in self.my_map_imb.iter() {
                if !value.is_empty() {
                    wext = Some(key.clone());
                }
            }
        }

        // OCCT L494-505.
        if wext.is_none() && !infinite {
            if nbwire > 0 {
                self.my_orires = BRepCheckStatus::InvalidImbricationOfWires;
            }
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("OrientationOfWires: myShape must be bound");
                brep_check_add(lst, self.my_orires);
            }
            return self.my_orires;
        }

        // OCCT L508-550.
        let fwd_face = oriented(&my_shape, Orientation::Forward);
        for wir in explorer(brep, &fwd_face, ShapeType::Wire) {
            if let Some(wext_s) = &wext {
                if wir.is_same(wext_s) {
                    if wir.orientation != wext_s.orientation {
                        // the exterior wire defines a hole
                        // OCCT L517-520.
                        if check_thin(brep, &wir, &fwd_face) {
                            return self.my_orires;
                        }
                        self.my_orires = BRepCheckStatus::BadOrientationOfSubshape;
                        if update {
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("OrientationOfWires: myShape must be bound");
                            brep_check_add(lst, self.my_orires);
                        }
                        return self.my_orires;
                    }
                } else {
                    // OCCT L531-548: find the wire key in myMapImb.
                    let mut found: Option<Shape> = None;
                    for (key, _) in self.my_map_imb.iter() {
                        if key.is_same(&wir) {
                            found = Some(key.clone());
                            break;
                        }
                    }
                    // No control on More()
                    if let Some(key) = found {
                        if key.orientation == wir.orientation {
                            // the given wire does not define a hole
                            self.my_orires = BRepCheckStatus::BadOrientationOfSubshape;
                            if update {
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(&my_shape)
                                    .expect("OrientationOfWires: myShape must be bound");
                                brep_check_add(lst, self.my_orires);
                            }
                            return self.my_orires;
                        }
                    }
                }
            } else {
                // OCCT L531-548 with a null Wext — unreachable when Wext is
                // null and the face is infinite: the OCCT loop runs the same
                // key lookup.
                let mut found: Option<Shape> = None;
                for (key, _) in self.my_map_imb.iter() {
                    if key.is_same(&wir) {
                        found = Some(key.clone());
                        break;
                    }
                }
                if let Some(key) = found {
                    if key.orientation == wir.orientation {
                        self.my_orires = BRepCheckStatus::BadOrientationOfSubshape;
                        if update {
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("OrientationOfWires: myShape must be bound");
                            brep_check_add(lst, self.my_orires);
                        }
                        return self.my_orires;
                    }
                }
            }
        }
        // OCCT L552-556: quit without error.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("OrientationOfWires: myShape must be bound");
            brep_check_add(lst, self.my_orires);
        }
        self.my_orires
    }

    /// OCCT BRepCheck_Face::SetUnorientable (Face.cxx L561-569).
    pub fn set_unorientable(&mut self) {
        let lst = self
            .base
            .my_map
            .find_mut(&self.base.my_shape)
            .expect("SetUnorientable: myShape must be bound");
        brep_check_add(lst, BRepCheckStatus::UnorientableShape);
    }

    /// OCCT BRepCheck_Face::SetStatus (Face.cxx L573-581).
    pub fn set_status(&mut self, the_status: BRepCheckStatus) {
        let lst = self
            .base
            .my_map
            .find_mut(&self.base.my_shape)
            .expect("SetStatus: myShape must be bound");
        brep_check_add(lst, the_status);
    }

    /// OCCT BRepCheck_Face::IsUnorientable (Face.cxx L585-599).
    pub fn is_unorientable(&self) -> bool {
        if self.my_oridone {
            return self.my_orires != BRepCheckStatus::NoError;
        }
        if let Some(list) = self.base.my_map.find(&self.base.my_shape) {
            for itl in list {
                if *itl == BRepCheckStatus::UnorientableShape {
                    return true;
                }
            }
        }
        false
    }

    /// OCCT BRepCheck_Face::GeometricControls(B) (Face.cxx L603-615).
    pub fn set_geometric_controls(&mut self, b: bool) {
        if self.my_gctrl != b {
            if b {
                self.my_intdone = false;
                self.my_imbdone = false;
                self.my_oridone = false;
            }
            self.my_gctrl = b;
        }
    }

    /// OCCT BRepCheck_Face::GeometricControls() (Face.cxx L619-622).
    pub fn geometric_controls(&self) -> bool {
        self.my_gctrl
    }
}

/// OCCT Face.cxx L626-798: the static `Intersect(wir1, wir2, F, theMapEdgeBox)`.
///
/// Runs the real 2D curve/curve intersector — `crate::geomalgo::geom2d_int::GInter`,
/// the `Geom2dInt_GInter` real body (`Geom2dInt_GInter_0.cxx` +
/// `IntCurve_IntCurveCurveGen.gxx`).
///
/// `c1` / `c2` / `myDomain1` / `myDomain2` mirror the OCCT locals declared at
/// L665-669 and assigned inside the loop heads; Rust needs an initializer the
/// first assignment always overwrites.
#[allow(unused_assignments)]
pub fn intersect(
    brep: &BRep,
    wir1: &Shape,
    wir2: &Shape,
    f: &Shape,
    the_map_edge_box: &DataMapOfShapeBox2d,
) -> bool {
    // OCCT L631: Inter2dTol = 1.e-10.
    let inter2d_tol = 1e-10;

    // OCCT L636-649: find the common vertices of the two wires
    // (non-manifold case).
    let mut map_w1 = super::brep_check_wire::ShapeSet::new();
    let mut common_vertices: Vec<Shape> = Vec::new();
    for exp1 in explorer(brep, wir1, ShapeType::Vertex) {
        map_w1.add(&exp1);
    }
    for exp2 in explorer(brep, wir2, ShapeType::Vertex) {
        let v = exp2;
        if map_w1.contains(&v) {
            common_vertices.push(v);
        }
    }

    // OCCT L653: BRepAdaptor_Surface Surf(F, false) — the pure surface
    // adaptor (no UVBounds computation).
    let Some(surf) = brep.face_surface_world(f) else {
        return false;
    };

    // OCCT L655-663: PntSeq — the common vertices projected on the surface.
    // `BRep_Tool::Parameters(V, F)` raises Standard_NoSuchObject when the
    // vertex carries no parameter on the face; the rcad form then leaves the
    // vertex out of PntSeq.
    let mut pnt_seq: Vec<glam::DVec3> = Vec::new();
    for i in 1..=common_vertices.len() {
        // OCCT L658-661.
        let v = &common_vertices[i - 1];
        if let Some(p2d) = brep_vertex_parameters(brep, v, f) {
            let p = surf.point_at(p2d.x, p2d.y);
            pnt_seq.push(p);
        }
    }

    // OCCT L665-670: C1 / C2, the UV points, the ranges, Inter, the domains
    // and Box1 / Box2.
    let mut c1: Option<Curve2d> = None;
    let mut c2: Option<Curve2d> = None;
    let mut my_domain1 = Res2dDomain::infinite();
    let mut my_domain2 = Res2dDomain::infinite();
    let mut box1 = BndBox2d::new();
    let mut box2 = BndBox2d::new();
    let mut inter = GInter::new();

    // OCCT L672: for (exp1.Init(wir1, TopAbs_EDGE); exp1.More(); exp1.Next())
    for edg1 in explorer(brep, wir1, ShapeType::Edge) {
        // OCCT L676: C1.Load(BRep_Tool::CurveOnSurface(edg1, F, first1, last1)).
        let Some((pc1, first1_raw, last1_raw)) = brep.curve_on_surface(&edg1, f) else {
            // OCCT: Geom2dAdaptor_Curve::Load dereferences the null handle
            // (`C->DynamicType()`), i.e. undefined behaviour — the rcad form
            // skips the edge instead.
            continue;
        };
        c1 = Some(pc1.clone());
        // OCCT L678-687: clamp onto the adaptor range (to avoid exception in
        // Segment if C1 is BSpline - IFV).
        let mut first1 = first1_raw;
        let mut last1 = last1_raw;
        if Curve2dEval::default_domain(&pc1)[0] > first1 {
            first1 = Curve2dEval::default_domain(&pc1)[0];
        }
        if Curve2dEval::default_domain(&pc1)[1] < last1 {
            last1 = Curve2dEval::default_domain(&pc1)[1];
        }

        // OCCT L689-698: Box1.
        box1.set_void();
        if let Some((b, _s)) = the_map_edge_box.get(&ShapeKey::of(&edg1)) {
            box1 = b.clone();
        }
        if box1.is_void() {
            // OCCT L698: BndLib_Add2dCurve::Add(C1, first1, last1, 0., Box1).
            add_2d_curve(&pc1, first1, last1, 0., &mut box1);
        }

        // OCCT L699: for (exp2.Init(wir2, TopAbs_EDGE); exp2.More(); exp2.Next())
        for edg2 in explorer(brep, wir2, ShapeType::Edge) {
            // OCCT L702.
            if edg1.is_same(&edg2) {
                continue;
            }
            // OCCT L705: C2.Load(BRep_Tool::CurveOnSurface(edg2, F, first2, last2)).
            let Some((pc2, first2_raw, last2_raw)) = brep.curve_on_surface(&edg2, f) else {
                // Same Geom2dAdaptor_Curve::Load(null) architecture note as
                // for C1 above.
                continue;
            };
            c2 = Some(pc2.clone());
            // OCCT L707-716: clamp onto the adaptor range.
            let mut first2 = first2_raw;
            let mut last2 = last2_raw;
            if Curve2dEval::default_domain(&pc2)[0] > first2 {
                first2 = Curve2dEval::default_domain(&pc2)[0];
            }
            if Curve2dEval::default_domain(&pc2)[1] < last2 {
                last2 = Curve2dEval::default_domain(&pc2)[1];
            }

            // OCCT L718-728: Box2.
            box2.set_void();
            if let Some((b, _s)) = the_map_edge_box.get(&ShapeKey::of(&edg2)) {
                box2 = b.clone();
            }
            if box2.is_void() {
                add_2d_curve(&pc2, first2, last2, 0., &mut box2);
            }

            // OCCT L729.
            if box1.is_out_box(&box2) {
                continue;
            }

            // OCCT L731-737: UVPoints + the domains.
            let (pfirst1, plast1) = uv_points(brep, &edg1, f, &pc1, first1_raw, last1_raw);
            my_domain1 = Res2dDomain::bounded(pfirst1, first1, inter2d_tol, plast1, last1, inter2d_tol);
            let (pfirst2, plast2) = uv_points(brep, &edg2, f, &pc2, first2_raw, last2_raw);
            my_domain2 = Res2dDomain::bounded(pfirst2, first2, inter2d_tol, plast2, last2, inter2d_tol);

            // OCCT L738: Inter.Perform(C1, myDomain1, C2, myDomain2,
            // Inter2dTol, Inter2dTol).
            let c1_ref = c1.as_ref().expect("Intersect: C1 must be loaded");
            let c2_ref = c2.as_ref().expect("Intersect: C2 must be loaded");
            inter.perform_cd_cd(c1_ref, &my_domain1, c2_ref, &my_domain2, inter2d_tol, inter2d_tol);

            // OCCT L739-742.
            if !inter.is_done() {
                return true;
            }

            // OCCT L743-778: the intersection segments.
            if inter.nb_segments() > 0 {
                // OCCT L745-748.
                if pnt_seq.is_empty() {
                    return true;
                }
                // OCCT L750-777.
                let mut nb_coinc = 0;
                for i in 1..=inter.nb_segments() {
                    // OCCT L753-756.
                    let seg = inter.segment(i);
                    if !seg.has_first_point() || !seg.has_last_point() {
                        return true;
                    }
                    // OCCT L757-762.
                    let first_p2d = seg.first_point().value();
                    let last_p2d = seg.last_point().value();
                    let first_p = surf.point_at(first_p2d.x, first_p2d.y);
                    let last_p = surf.point_at(last_p2d.x, last_p2d.y);
                    // OCCT L763-771.
                    for j in 1..=pnt_seq.len() {
                        let tolv = super::brep_check_result::brep_tool_tolerance_vertex(
                            brep,
                            &common_vertices[j - 1],
                        );
                        if first_p.distance(pnt_seq[j - 1]) <= tolv
                            || last_p.distance(pnt_seq[j - 1]) <= tolv
                        {
                            nb_coinc += 1;
                            break;
                        }
                    }
                }
                // OCCT L776.
                return nb_coinc != inter.nb_segments();
            }

            // OCCT L779-795: the intersection points.
            if inter.nb_points() > 0 {
                // OCCT L781-784.
                if pnt_seq.is_empty() {
                    return true;
                }
                // OCCT L786-793.
                let mut nb_coinc = 0;
                for i in 1..=inter.nb_points() {
                    let p2d = inter.point(i).value();
                    let p = surf.point_at(p2d.x, p2d.y);
                    for j in 1..=pnt_seq.len() {
                        let mut tolv = super::brep_check_result::brep_tool_tolerance_vertex(
                            brep,
                            &common_vertices[j - 1],
                        );
                        // OCCT L789: possible tolerance of intersection point.
                        tolv += 1.0e-8;
                        let dd = p.distance_squared(pnt_seq[j - 1]);
                        if dd <= tolv * tolv {
                            nb_coinc += 1;
                            break;
                        }
                    }
                }
                // OCCT L794.
                return nb_coinc != inter.nb_points();
            }
        }
    }
    // OCCT L797.
    false
}

/// `BRep_Tool::Parameters(V, F)` — the vertex UV parameters on the face
/// (resolved through the face-point representation; falls back to the
/// projection of the 3D point on the surface).
fn brep_vertex_parameters(brep: &BRep, v: &Shape, f: &Shape) -> Option<glam::DVec2> {
    let vd = v.as_vertex()?;
    for pr in &vd.points {
        if let rcad_kernel::topods::PointRepresentation::PointOnSurface { face: fidx, u, v: vv, .. } = pr
        {
            let _ = brep;
            if *fidx == f.index {
                return Some(glam::DVec2::new(*u, *vv));
            }
        }
    }
    None
}

/// OCCT Face.cxx L802-875: the static `IsInside` — whether the wire lies
/// inside the current face region per the fresh-face classifier.
pub fn is_inside(
    brep: &BRep,
    the_wire: &Shape,
    wire_bien_oriente: bool,
    f_class2d: &FClass2d,
    ds: &dyn crate::topalgo::shape_source::ShapeSource,
    the_face: &Shape,
) -> bool {
    for an_explorer in explorer(brep, the_wire, ShapeType::Edge) {
        let an_edge = an_explorer;
        // OCCT L813: aCurve2D = BRep_Tool::CurveOnSurface(anEdge, theFace, aFirst, aLast).
        let Some((a_curve2d, a_first, a_last)) = brep.curve_on_surface(&an_edge, the_face) else {
            continue;
        };
        // OCCT L816-858: the parameter selection.
        let a_parameter: f64;
        if !rcad_kernel::precision::is_negative_infinite_value(a_first)
            && !rcad_kernel::precision::is_positive_infinite_value(a_last)
        {
            a_parameter = (a_first + a_last) * 0.5;

            // OCCT L821-824: skip when the parametric range is too small.
            if (a_parameter - a_first).abs() < rcad_kernel::precision::PCONFUSION {
                continue;
            }

            // OCCT L827-842: skip when the edge length is too small.
            let Some(ed) = an_edge.as_edge() else { continue };
            let Some(a_curve) = &ed.curve else { continue };
            let dom = rcad_kernel::geom::CurveEval::default_domain(a_curve);
            let a_first3d = dom[0];
            let a_last3d = dom[1];
            let a_points0 = a_curve.point_at(a_first);
            let a_points1 = a_curve.point_at((a_first3d + a_last3d) / 2.);
            if a_points0.distance(a_points1) < rcad_kernel::CONFUSION {
                continue;
            }
        } else if rcad_kernel::precision::is_negative_infinite_value(a_first)
            && rcad_kernel::precision::is_positive_infinite_value(a_last)
        {
            a_parameter = 0.;
        } else if rcad_kernel::precision::is_negative_infinite_value(a_first) {
            a_parameter = a_last - 1.;
        } else {
            a_parameter = a_first + 1.;
        }

        // OCCT L861-863: the classification of the point.
        let a_point2d = a_curve2d.point_at(a_parameter);
        let a_state = f_class2d.perform(ds, a_point2d, false);

        // OCCT L865-872.
        if wire_bien_oriente {
            return a_state == State::Out;
        } else {
            return a_state == State::In;
        }
    }
    // OCCT L874.
    false
}

/// OCCT Face.cxx L877-954: `CheckThin` — whether the two-edge wire is a thin
/// slit (the two edges run in opposite directions between the same
/// vertices).
pub fn check_thin(brep: &BRep, w: &Shape, f: &Shape) -> bool {
    let a_f = f;
    let a_w = w;

    // OCCT L883-891: the wire must have exactly 2 edges.
    let l_e = explorer(brep, a_w, ShapeType::Edge);
    let nb_e = l_e.len() as i32;
    if nb_e != 2 {
        return false;
    }
    let e1 = &l_e[0];
    let e2 = &l_e[l_e.len() - 1];

    // OCCT L899-910: TopExp::Vertices.
    let (Some(ed1), Some(ed2)) = (e1.as_edge(), e2.as_edge()) else {
        return false;
    };
    let v1 = &ed1.first;
    let v2 = &ed1.last;
    let v3 = &ed2.first;
    let v4 = &ed2.last;

    if v1.is_null() || v2.is_null() || v3.is_null() || v4.is_null() {
        return false;
    }

    if v1.is_same(v2) || v3.is_same(v4) {
        return false;
    }

    // OCCT L913-926.
    let mut s_f = false;
    let mut s_l = false;
    if v1.is_same(v3) || v1.is_same(v4) {
        s_f = true;
    }
    if v2.is_same(v3) || v2.is_same(v4) {
        s_l = true;
    }
    if !s_f || !s_l {
        return false;
    }

    // OCCT L928-933.
    let e1or = e1.orientation;
    let e2or = e2.orientation;

    let pc1: Option<(rcad_kernel::geom::Curve2d, f64, f64)> = brep.curve_on_surface(e1, a_f);
    let pc2: Option<(rcad_kernel::geom::Curve2d, f64, f64)> = brep.curve_on_surface(e2, a_f);
    let (Some((pc1, f1, l1)), Some((pc2, f2, l2))) = (pc1, pc2) else {
        return false;
    };

    // OCCT L940-948.
    let d1 = (l1 - f1).abs() / 100.;
    let d2 = (l2 - f2).abs() / 100.;
    let m1 = (l1 + f1) * 0.5;
    let m2 = (l2 + f2) * 0.5;

    let p1f = pc1.point_at(m1 - d1);
    let p1l = pc1.point_at(m1 + d1);
    let p2f = pc2.point_at(m2 - d2);
    let p2l = pc2.point_at(m2 + d2);

    // OCCT L950-953: (vc1 * vc2) < 0 || e1or != e2or.
    let vc1 = p1l - p1f;
    let vc2 = p2l - p2f;
    vc1.dot(vc2) < 0. || e1or != e2or
}

/// OCCT TopoDS_Shape::EmptyCopied (BRep_TFace::EmptyCopy) + B.Add(newFace,
/// wir1) + Orientation(FORWARD) — the fresh single-wire face of
/// ClassifyWires (Face.cxx L354-358).
fn empty_copied_face(my_shape: &Shape, wir1: &Shape) -> Option<Shape> {
    let fd = my_shape.as_face()?;
    let new_data = TFaceData {
        my_shapes: vec![wir1.clone()],
        flags: fd.flags,
        surface: fd.surface.clone(),
        surface_location: fd.surface_location,
        outer_wire: wir1.clone(),
        inner_wires: Vec::new(),
        sample_point: fd.sample_point,
        uv_domain: fd.uv_domain,
        internal_vertices: Vec::new(),
        tolerance: fd.tolerance,
        natural_restriction: fd.natural_restriction,
    };
    Some(Shape {
        data: Arc::new(TShape::Face(new_data)),
        index: usize::MAX,
        location: my_shape.location,
        orientation: Orientation::Forward,
    })
}

/// A FaceShapeSource over `face` with the DS-convention location table.
/// `surf` is the world surface of the face (computed by the caller — a fresh
/// EmptyCopied face carries no pool index to resolve it from).
fn face_shape_source<'a>(
    brep: &'a BRep,
    face: &'a Shape,
    surf: Surface3,
) -> Option<FaceShapeSource<'a>> {
    let mut locations: Vec<glam::DAffine3> = vec![glam::DAffine3::IDENTITY];
    locations.extend(brep.locations.iter().copied());
    Some(FaceShapeSource::new(face, surf, &locations))
}

/// Unused-import guard (Surface3 is referenced by doc links above).
const _: Option<Surface3> = None;
