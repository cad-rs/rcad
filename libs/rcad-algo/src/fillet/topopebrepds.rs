//! OCCT TopOpeBRepDS package (TKBool/TopOpeBRepDS) — 1:1 translation of the
//! subset used by the ChFi3d builder: the DataStructure tables (surfaces /
//! curves / points / shapes), the TopOpeBRepDS_Curve / _Surface / _Point
//! descriptors and the interference bookkeeping
//! (TopOpeBRepDS_Interference + SurfaceCurveInterference /
//! CurvePointInterference / SolidSurfaceInterference).
//!
//! Sources:
//!   - TopOpeBRepDS_DataStructure.cxx (AddSurface L64, AddCurve L116,
//!     AddPoint L187, AddShape L246, interference accessors L384-496,
//!     NbShapes L1036, HasGeometry L1208)
//!   - TopOpeBRepDS_Interference.cxx / _SurfaceCurveInterference /
//!     _CurvePointInterference / _SolidSurfaceInterference
//!   - TopOpeBRepDS_Curve.hxx (SetRange / Tolerance / Nullify / SetSCI)
//!   - TopOpeBRepDS_Surface.hxx, TopOpeBRepDS_Point.hxx

use std::collections::HashMap;

use glam::DVec3;
use rcad_kernel::geom::{Curve2d, Curve3, Surface3};
use rcad_kernel::topo::topods::{Orientation, Shape};

// =========================================================================
// OCCT TopOpeBRepDS_Kind (TopOpeBRepDS_Kind.hxx).
// =========================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopOpeBRepDSKind {
    Unknown,
    Vertex,
    Edge,
    Wire,
    Face,
    Shell,
    Solid,
    CompSolid,
    Compound,
    Surface,
    Curve,
    Point,
}

// =========================================================================
// OCCT TopOpeBRepDS_Transition — the fillet flow reads only the IN-state
// orientation (Transition().Orientation(TopAbs_IN)), carried by the plain
// TopAbs_Orientation the callers pass in.
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct TopOpeBRepDSTransition {
    pub orientation: Orientation,
}

impl TopOpeBRepDSTransition {
    pub fn new(orientation: Orientation) -> Self {
        TopOpeBRepDSTransition { orientation }
    }

    /// OCCT TopOpeBRepDS_Transition::Orientation(TopAbs_IN).
    pub fn orientation_in(&self) -> Orientation {
        self.orientation
    }
}

// =========================================================================
// OCCT TopOpeBRepDS_Surface (TopOpeBRepDS_Surface.hxx).
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepDSSurface {
    pub surface: Surface3,
    pub tolerance: f64,
}

impl TopOpeBRepDSSurface {
    pub fn new(surface: Surface3, tolerance: f64) -> Self {
        TopOpeBRepDSSurface { surface, tolerance }
    }

    /// OCCT TopOpeBRepDS_Surface::Surface().
    pub fn surface(&self) -> &Surface3 {
        &self.surface
    }

    /// OCCT TopOpeBRepDS_Surface::Tolerance() / ChangeTolerance.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }
}

// =========================================================================
// OCCT TopOpeBRepDS_Curve (TopOpeBRepDS_Curve.hxx L36-140).  The curve
// handle is nullable (Nullify marks the curve as "used via SCI only").
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepDSCurve {
    /// OCCT: occ::handle<Geom_Curve> myCurve.
    pub curve: Option<Curve3>,
    /// OCCT: double myTolerance.
    pub tolerance: f64,
    /// OCCT: double myFirst / myLast (SetRange).
    pub first: f64,
    pub last: f64,
    /// OCCT: occ::handle<TopOpeBRepDS_Interference> mySCI1 / mySCI2 (SetSCI).
    pub sci1: Option<InterferenceRef>,
    pub sci2: Option<InterferenceRef>,
    /// OCCT TopOpeBRepDS_Curve.hxx: int myMother (Mother()); bool myKeep.
    pub mother: i32,
    pub keep: bool,
}

/// OCCT stores the two "curve/surface-curve" interference handles on the
/// curve; rcad copies the interference payload (the pcurve + indices).
#[derive(Debug, Clone)]
pub struct InterferenceRef {
    pub pcurve: Option<Curve2d>,
    pub index_s: i32,
    pub index_g: i32,
}

impl TopOpeBRepDSCurve {
    pub fn new(curve: Option<Curve3>, tolerance: f64) -> Self {
        TopOpeBRepDSCurve {
            curve,
            tolerance,
            first: 0.0,
            last: 0.0,
            sci1: None,
            sci2: None,
            mother: 0,
            keep: false,
        }
    }

    /// OCCT TopOpeBRepDS_Curve::Curve().
    pub fn curve(&self) -> Option<&Curve3> {
        self.curve.as_ref()
    }

    /// OCCT TopOpeBRepDS_Curve::ChangeCurve() + Nullify().
    pub fn change_curve(&mut self) -> &mut Option<Curve3> {
        &mut self.curve
    }

    pub fn nullify(&mut self) {
        self.curve = None;
    }

    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    /// OCCT TopOpeBRepDS_Curve::SetRange(First, Last).
    pub fn set_range(&mut self, first: f64, last: f64) {
        self.first = first;
        self.last = last;
    }

    /// OCCT TopOpeBRepDS_Curve::SetSCI(I, S) — stores the interference pair
    /// used by the TopOpeBRepBuild reconstruction.
    pub fn set_sci(
        &mut self,
        sci1: TopOpeBRepDSSurfaceCurveInterference,
        _sci2: Option<()>,
    ) {
        self.sci1 = Some(InterferenceRef {
            pcurve: sci1.pcurve.clone(),
            index_s: sci1.index_s,
            index_g: sci1.index_g,
        });
    }
}

// =========================================================================
// OCCT TopOpeBRepDS_Point (TopOpeBRepDS_Point.hxx).
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepDSPoint {
    pub point: DVec3,
    pub tolerance: f64,
}

impl TopOpeBRepDSPoint {
    pub fn new(point: DVec3, tolerance: f64) -> Self {
        TopOpeBRepDSPoint { point, tolerance }
    }

    pub fn point(&self) -> DVec3 {
        self.point
    }

    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// OCCT TopOpeBRepDS_Point::Tolerance(Tol) (assignable slot).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }
}

// =========================================================================
// OCCT TopOpeBRepDS_Interference hierarchy.  The three concrete
// interferences used by the fillet builder share the
// (transition, support-kind/index, geometry-kind/index) payload; each adds
// its own data (pcurve for SurfaceCurve, parameter for CurvePoint).
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepDSSurfaceCurveInterference {
    pub transition: TopOpeBRepDSTransition,
    pub kind_s: TopOpeBRepDSKind,
    pub index_s: i32,
    pub kind_g: TopOpeBRepDSKind,
    pub index_g: i32,
    pub pcurve: Option<Curve2d>,
}

impl TopOpeBRepDSSurfaceCurveInterference {
    /// OCCT constructor (Transition, KindS, IndexS, KindG, IndexG, PC).
    pub fn new(
        transition: Orientation,
        kind_s: TopOpeBRepDSKind,
        index_s: i32,
        kind_g: TopOpeBRepDSKind,
        index_g: i32,
        pcurve: Option<Curve2d>,
    ) -> Self {
        TopOpeBRepDSSurfaceCurveInterference {
            transition: TopOpeBRepDSTransition::new(transition),
            kind_s,
            index_s,
            kind_g,
            index_g,
            pcurve,
        }
    }

    /// OCCT TopOpeBRepDS_Interference::GKGSKS(GK, G, SK, S).
    pub fn gkgsks(&self) -> (TopOpeBRepDSKind, i32, TopOpeBRepDSKind, i32) {
        (self.kind_g, self.index_g, self.kind_s, self.index_s)
    }
}

#[derive(Debug, Clone)]
pub struct TopOpeBRepDSCurvePointInterference {
    pub transition: TopOpeBRepDSTransition,
    pub kind_s: TopOpeBRepDSKind,
    pub index_s: i32,
    pub kind_g: TopOpeBRepDSKind,
    pub index_g: i32,
    pub parameter: f64,
}

impl TopOpeBRepDSCurvePointInterference {
    /// OCCT constructor (Transition, KindS, IndexS, KindG, IndexG, Par).
    pub fn new(
        transition: Orientation,
        kind_s: TopOpeBRepDSKind,
        index_s: i32,
        kind_g: TopOpeBRepDSKind,
        index_g: i32,
        parameter: f64,
    ) -> Self {
        TopOpeBRepDSCurvePointInterference {
            transition: TopOpeBRepDSTransition::new(transition),
            kind_s,
            index_s,
            kind_g,
            index_g,
            parameter,
        }
    }

    pub fn gkgsks(&self) -> (TopOpeBRepDSKind, i32, TopOpeBRepDSKind, i32) {
        (self.kind_g, self.index_g, self.kind_s, self.index_s)
    }
}

#[derive(Debug, Clone)]
pub struct TopOpeBRepDSSolidSurfaceInterference {
    pub transition: TopOpeBRepDSTransition,
    pub kind_s: TopOpeBRepDSKind,
    pub index_s: i32,
    pub kind_g: TopOpeBRepDSKind,
    pub index_g: i32,
}

impl TopOpeBRepDSSolidSurfaceInterference {
    pub fn new(transition: Orientation, kind_s: TopOpeBRepDSKind, index_s: i32, kind_g: TopOpeBRepDSKind, index_g: i32) -> Self {
        TopOpeBRepDSSolidSurfaceInterference {
            transition: TopOpeBRepDSTransition::new(transition),
            kind_s,
            index_s,
            kind_g,
            index_g,
        }
    }

    pub fn gkgsks(&self) -> (TopOpeBRepDSKind, i32, TopOpeBRepDSKind, i32) {
        (self.kind_g, self.index_g, self.kind_s, self.index_s)
    }
}

/// OCCT occ::handle<TopOpeBRepDS_Interference> — the heterogeneous handle
/// carried by the interference lists.
#[derive(Debug, Clone)]
pub enum TopOpeBRepDSInterference {
    SurfaceCurve(TopOpeBRepDSSurfaceCurveInterference),
    CurvePoint(TopOpeBRepDSCurvePointInterference),
    SolidSurface(TopOpeBRepDSSolidSurfaceInterference),
}

impl TopOpeBRepDSInterference {
    pub fn gkgsks(&self) -> (TopOpeBRepDSKind, i32, TopOpeBRepDSKind, i32) {
        match self {
            TopOpeBRepDSInterference::SurfaceCurve(i) => i.gkgsks(),
            TopOpeBRepDSInterference::CurvePoint(i) => i.gkgsks(),
            TopOpeBRepDSInterference::SolidSurface(i) => i.gkgsks(),
        }
    }

    /// OCCT TopOpeBRepDS_Interference::Transition().
    pub fn transition(&self) -> &TopOpeBRepDSTransition {
        match self {
            TopOpeBRepDSInterference::SurfaceCurve(i) => &i.transition,
            TopOpeBRepDSInterference::CurvePoint(i) => &i.transition,
            TopOpeBRepDSInterference::SolidSurface(i) => &i.transition,
        }
    }

    pub fn parameter(&self) -> f64 {
        match self {
            TopOpeBRepDSInterference::CurvePoint(i) => i.parameter,
            _ => 0.0,
        }
    }
}

// =========================================================================
// OCCT TopOpeBRepDS_DataStructure — the tables themselves.  Indices are
// 1-based exactly as in OCCT.
// =========================================================================
#[derive(Debug, Clone, Default)]
pub struct TopOpeBRepDSHDataStructure {
    /// OCCT: TopOpeBRepDS_Surface table.
    pub surfaces: Vec<TopOpeBRepDSSurface>,
    /// OCCT: TopOpeBRepDS_Curve table.
    pub curves: Vec<TopOpeBRepDSCurve>,
    /// OCCT: TopOpeBRepDS_Point table.
    pub points: Vec<TopOpeBRepDSPoint>,
    /// OCCT: TopOpeBRepDS_ShapeInfo table (the shapes themselves).
    pub shapes: Vec<Shape>,
    /// OCCT: shape -> index lookup by TShape pointer identity.
    pub shape_index: HashMap<u64, i32>,
    /// OCCT: per-shape interference lists (ChangeShapeInterferences(I)).
    pub shape_interferences: HashMap<i32, Vec<TopOpeBRepDSInterference>>,
    /// OCCT: per-curve interference lists (ChangeCurveInterferences(I)).
    pub curve_interferences: HashMap<i32, Vec<TopOpeBRepDSInterference>>,
    /// OCCT: per-surface interference lists (ChangeSurfaceInterferences(I)).
    pub surface_interferences: HashMap<i32, Vec<TopOpeBRepDSInterference>>,
}

impl TopOpeBRepDSHDataStructure {
    /// OCCT TopOpeBRepDS_DataStructure::AddSurface(S) — 1-based index.
    pub fn add_surface(&mut self, s: TopOpeBRepDSSurface) -> i32 {
        self.surfaces.push(s);
        self.surfaces.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::AddCurve(C) — 1-based index.
    pub fn add_curve(&mut self, c: TopOpeBRepDSCurve) -> i32 {
        self.curves.push(c);
        self.curves.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::AddPoint(PDS) — 1-based index.
    pub fn add_point(&mut self, p: TopOpeBRepDSPoint) -> i32 {
        self.points.push(p);
        self.points.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::AddShape(S) — an already-present
    /// shape keeps its index.
    pub fn add_shape(&mut self, s: &Shape) -> i32 {
        if let Some(i) = self.shape_index.get(&s.ptr_id()) {
            return *i;
        }
        self.shapes.push(s.clone());
        let i = self.shapes.len() as i32;
        self.shape_index.insert(s.ptr_id(), i);
        i
    }

    /// OCCT TopOpeBRepDS_DataStructure::NbShapes().
    pub fn nb_shapes(&self) -> i32 {
        self.shapes.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::NbCurves().
    pub fn nb_curves(&self) -> i32 {
        self.curves.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::Shape(I) — 1-based.
    pub fn shape(&self, i: i32) -> &Shape {
        &self.shapes[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeShape(I).
    pub fn change_shape(&mut self, i: i32) -> &mut Shape {
        &mut self.shapes[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::Curve(I).
    pub fn curve(&self, i: i32) -> &TopOpeBRepDSCurve {
        &self.curves[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeCurve(I).
    pub fn change_curve(&mut self, i: i32) -> &mut TopOpeBRepDSCurve {
        &mut self.curves[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::Surface(I).
    pub fn surface(&self, i: i32) -> &TopOpeBRepDSSurface {
        &self.surfaces[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeSurface(I).
    pub fn change_surface(&mut self, i: i32) -> &mut TopOpeBRepDSSurface {
        &mut self.surfaces[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::Point(I).
    pub fn point(&self, i: i32) -> &TopOpeBRepDSPoint {
        &self.points[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangePoint(I).
    pub fn change_point(&mut self, i: i32) -> &mut TopOpeBRepDSPoint {
        &mut self.points[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::HasGeometry(S) — the shape carries
    /// interferences (DataStructure.cxx L1208-1213: has =
    /// !ShapeInterferences(S).IsEmpty()).
    pub fn has_geometry(&self, s: &Shape) -> bool {
        if let Some(i) = self.shape_index.get(&s.ptr_id()) {
            self.shape_interferences
                .get(i)
                .is_some_and(|l| !l.is_empty())
        } else {
            false
        }
    }

    /// OCCT TopOpeBRepDS_DataStructure::ShapeInterferences(I) (const).
    pub fn shape_interferences(&self, i: i32) -> &[TopOpeBRepDSInterference] {
        self.shape_interferences
            .get(&i)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeShapeInterferences(I).
    pub fn change_shape_interferences(&mut self, i: i32) -> &mut Vec<TopOpeBRepDSInterference> {
        self.shape_interferences.entry(i).or_default()
    }

    /// OCCT TopOpeBRepDS_DataStructure::ShapeInterferences(S, FindKeep).
    pub fn shape_interferences_of(&self, s: &Shape) -> &[TopOpeBRepDSInterference] {
        match self.shape_index.get(&s.ptr_id()) {
            Some(i) => self.shape_interferences(*i),
            None => &[],
        }
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeShapeInterferences(S).
    pub fn change_shape_interferences_of(&mut self, s: &Shape) -> &mut Vec<TopOpeBRepDSInterference> {
        let i = self.add_shape(s);
        self.change_shape_interferences(i)
    }

    /// OCCT TopOpeBRepDS_DataStructure::CurveInterferences(I) (const).
    pub fn curve_interferences(&self, i: i32) -> &[TopOpeBRepDSInterference] {
        self.curve_interferences
            .get(&i)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeCurveInterferences(I).
    pub fn change_curve_interferences(&mut self, i: i32) -> &mut Vec<TopOpeBRepDSInterference> {
        self.curve_interferences.entry(i).or_default()
    }

    /// OCCT TopOpeBRepDS_DataStructure::SurfaceInterferences(I) (const).
    pub fn surface_interferences(&self, i: i32) -> &[TopOpeBRepDSInterference] {
        self.surface_interferences
            .get(&i)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeSurfaceInterferences(I).
    pub fn change_surface_interferences(&mut self, i: i32) -> &mut Vec<TopOpeBRepDSInterference> {
        self.surface_interferences.entry(i).or_default()
    }
}
