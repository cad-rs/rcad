//! ChFi3d data-structure facade (Stage 0.4, user ruling D6).
//!
//! OCCT 8.0 ChFi3d (TKFillet) records its intersection results (fillet
//! surfaces / spine curves / intersection points / interference
//! quadruplets) into `TopOpeBRepDS_HDataStructure`, consumed by the legacy
//! boolean reconstruction engine `TopOpeBRepBuild_HBuilder`.  Per the D6
//! ruling the old TKBool `TopOpeBRepDS` carrier is retired in rcad:
//!
//!   - shape registration routes to the rcad TKBO pipeline's BOPDS data
//!     structure (`crate::bop::ds::DS`) — the same carrier the boolean
//!     PaveFiller uses (HBuilder itself is replaced by the TKBO pipeline
//!     equivalent, Stage 1g);
//!   - the geometry payloads BOPDS has no equivalent for (Kind-indexed
//!     Surface / Curve / Point tables) live in the ChFi3d side tables
//!     [`ChFi3dDSSideTables`];
//!   - the per-shape / per-curve / per-surface interference quadruplet
//!     lists stay on the facade: the BOPDS interference model is a set of
//!     typed global tables (InterfVV / EE / EF / FF / ...) with no
//!     per-shape quadruplet-list equivalent — D6 architecture difference;
//!     the interference quadruplets are ChFi3d-private records.
//!
//! Routing map (OCCT TopOpeBRepDS_DataStructure -> rcad):
//!   AddShape / Shape / ChangeShape / NbShapes -> `bopds` (BOPDS shape table)
//!   AddSurface / AddCurve / AddPoint + accessors -> `side` ([`ChFi3dDSSideTables`])
//!   ShapeInterferences / CurveInterferences / SurfaceInterferences -> facade fields
//!
//! The type name [`TopOpeBRepDSHDataStructure`] is kept as the ChFi3d call
//! surface name (precedent: `hbuilder::TopOpeBRepBuildHBuilder`).  The
//! value types below are the form carriers for the line-by-line ChFi3d
//! alignment; per D6 they are no longer considered a TopOpeBRepDS
//! translation.
//!
//! OCCT anchors carried over from the retired
//! `fillet/topopebrepds.rs`:
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

use crate::bop::ds::DS;

// =========================================================================
// OCCT TopOpeBRepDS_Kind (TopOpeBRepDS_Kind.hxx).
// Value type = form carrier for the ChFi3d line-by-line alignment; per D6
// no longer treated as a TopOpeBRepDS translation.
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
// ChFi3d side tables — the geometry payloads the BOPDS carrier has no
// equivalent for.  OCCT TopOpeBRepDS kept Kind-indexed Surface / Curve /
// Point tables alongside the shape table; BOPDS only carries shapes.
// Indices are 1-based exactly as in OCCT (AddSurface / AddCurve / AddPoint
// return 1-based indices; the accessors subtract 1 internally).
// =========================================================================
#[derive(Debug, Clone, Default)]
pub struct ChFi3dDSSideTables {
    /// OCCT: TopOpeBRepDS_Surface table.
    pub surfaces: Vec<TopOpeBRepDSSurface>,
    /// OCCT: TopOpeBRepDS_Curve table.
    pub curves: Vec<TopOpeBRepDSCurve>,
    /// OCCT: TopOpeBRepDS_Point table.
    pub points: Vec<TopOpeBRepDSPoint>,
}

// =========================================================================
// ChFi3d DS facade — the ChFi3d call surface of the retired
// TopOpeBRepDS_HDataStructure, carried by BOPDS + side tables (D6).
// =========================================================================
#[derive(Debug, Default)]
pub struct TopOpeBRepDSHDataStructure {
    /// rcad TKBO BOPDS data structure — the shape registration carrier
    /// (D6).  Populated only through `DS::append_shape`; the facade never
    /// runs the boolean `DS::Init` pass, so no ShapeInfo enrichment is
    /// needed beyond what `append_shape` fills (shape_type + the
    /// (ptr_id, location) identity map).
    pub bopds: DS,
    /// Geometry payloads (surface / curve / point tables), ChFi3d-private.
    pub side: ChFi3dDSSideTables,
    /// OCCT: per-shape interference lists (ChangeShapeInterferences(I)).
    /// Stays on the facade: the BOPDS interference model is typed global
    /// tables with no per-shape quadruplet-list equivalent — D6
    /// architecture difference; the interference quadruplets are ChFi3d
    /// private records.
    pub shape_interferences: HashMap<i32, Vec<TopOpeBRepDSInterference>>,
    /// OCCT: per-curve interference lists (ChangeCurveInterferences(I)).
    pub curve_interferences: HashMap<i32, Vec<TopOpeBRepDSInterference>>,
    /// OCCT: per-surface interference lists (ChangeSurfaceInterferences(I)).
    pub surface_interferences: HashMap<i32, Vec<TopOpeBRepDSInterference>>,
}

impl Clone for TopOpeBRepDSHDataStructure {
    fn clone(&self) -> Self {
        // BOPDS DS is not Clone (it carries the full boolean engine state:
        // pave block pools, typed interference tables, ...).  The facade
        // only ever fills the DS shape table through `DS::append_shape`, so
        // cloning re-appends the registered shapes in order, reproducing
        // the facade-visible DS state exactly (same indices, same identity
        // map keys).
        let mut bopds = DS::new();
        for si in &self.bopds.shapes {
            bopds.append_shape(si.shape.clone());
        }
        TopOpeBRepDSHDataStructure {
            bopds,
            side: self.side.clone(),
            shape_interferences: self.shape_interferences.clone(),
            curve_interferences: self.curve_interferences.clone(),
            surface_interferences: self.surface_interferences.clone(),
        }
    }
}

impl TopOpeBRepDSHDataStructure {
    /// OCCT TopOpeBRepDS_DataStructure::AddSurface(S) L64 — 1-based index.
    /// Routes to the ChFi3d side table (BOPDS has no surface payload).
    pub fn add_surface(&mut self, s: TopOpeBRepDSSurface) -> i32 {
        self.side.surfaces.push(s);
        self.side.surfaces.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::AddCurve(C) L116 — 1-based index.
    /// Routes to the ChFi3d side table (BOPDS has no curve descriptor
    /// payload; its BOPDS_Curve is an intersection-result record).
    pub fn add_curve(&mut self, c: TopOpeBRepDSCurve) -> i32 {
        self.side.curves.push(c);
        self.side.curves.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::AddPoint(PDS) L187 — 1-based index.
    /// Routes to the ChFi3d side table (BOPDS has no point payload).
    pub fn add_point(&mut self, p: TopOpeBRepDSPoint) -> i32 {
        self.side.points.push(p);
        self.side.points.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::AddShape(S) L246 — an
    /// already-present shape keeps its index.
    /// Routes to BOPDS `DS::index` / `DS::append_shape` (D6: shape
    /// registration lives in the BOPDS shape table).  Minimal registration:
    /// `append_shape` fills shape_type and the (ptr_id, location) identity
    /// map; no `DS::Init` / InitShape enrichment — ChFi3d only needs
    /// index <-> shape lookup.  The facade keeps OCCT's 1-based index
    /// semantics on top of the 0-based BOPDS index.  Identity-key note:
    /// BOPDS indexes by (TShape ptr_id, location) where the legacy carrier
    /// keyed by ptr_id alone.
    pub fn add_shape(&mut self, s: &Shape) -> i32 {
        let found = self.bopds.index(s);
        if found >= 0 {
            return found as i32 + 1;
        }
        let i = self.bopds.append_shape(s.clone());
        i as i32 + 1
    }

    /// OCCT TopOpeBRepDS_DataStructure::NbShapes() L1036.
    /// Routes to BOPDS `DS::nb_shapes` (shape table lives in BOPDS, D6).
    pub fn nb_shapes(&self) -> i32 {
        self.bopds.nb_shapes() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::NbCurves().
    /// Routes to the ChFi3d side table.
    pub fn nb_curves(&self) -> i32 {
        self.side.curves.len() as i32
    }

    /// OCCT TopOpeBRepDS_DataStructure::Shape(I) — 1-based.
    /// Routes to BOPDS `DS::shape` (shape table lives in BOPDS, D6).
    pub fn shape(&self, i: i32) -> &Shape {
        self.bopds.shape((i - 1) as usize)
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeShape(I).
    /// Routes to BOPDS `DS::change_shape_info(..).shape` (D6).
    pub fn change_shape(&mut self, i: i32) -> &mut Shape {
        &mut self.bopds.change_shape_info((i - 1) as usize).shape
    }

    /// OCCT TopOpeBRepDS_DataStructure::Curve(I).
    /// Routes to the ChFi3d side table.
    pub fn curve(&self, i: i32) -> &TopOpeBRepDSCurve {
        &self.side.curves[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeCurve(I).
    /// Routes to the ChFi3d side table.
    pub fn change_curve(&mut self, i: i32) -> &mut TopOpeBRepDSCurve {
        &mut self.side.curves[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::Surface(I).
    /// Routes to the ChFi3d side table.
    pub fn surface(&self, i: i32) -> &TopOpeBRepDSSurface {
        &self.side.surfaces[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeSurface(I).
    /// Routes to the ChFi3d side table.
    pub fn change_surface(&mut self, i: i32) -> &mut TopOpeBRepDSSurface {
        &mut self.side.surfaces[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::Point(I).
    /// Routes to the ChFi3d side table.
    pub fn point(&self, i: i32) -> &TopOpeBRepDSPoint {
        &self.side.points[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangePoint(I).
    /// Routes to the ChFi3d side table.
    pub fn change_point(&mut self, i: i32) -> &mut TopOpeBRepDSPoint {
        &mut self.side.points[(i - 1) as usize]
    }

    /// OCCT TopOpeBRepDS_DataStructure::HasGeometry(S) L1208-1213 (has =
    /// !ShapeInterferences(S).IsEmpty()).
    /// No direct BOPDS equivalent (BOPDS carries no per-shape interference
    /// list) — the legacy logic is kept on the facade; only the index
    /// lookup routes through BOPDS `DS::index`.
    pub fn has_geometry(&self, s: &Shape) -> bool {
        let i = self.bopds.index(s);
        if i >= 0 {
            self.shape_interferences
                .get(&(i as i32 + 1))
                .is_some_and(|l| !l.is_empty())
        } else {
            false
        }
    }

    /// OCCT TopOpeBRepDS_DataStructure::ShapeInterferences(I) L384-496
    /// (const).  Stays on the facade (D6: ChFi3d-private quadruplet list).
    pub fn shape_interferences(&self, i: i32) -> &[TopOpeBRepDSInterference] {
        self.shape_interferences
            .get(&i)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeShapeInterferences(I).
    /// Stays on the facade (D6: ChFi3d-private quadruplet list).
    pub fn change_shape_interferences(&mut self, i: i32) -> &mut Vec<TopOpeBRepDSInterference> {
        self.shape_interferences.entry(i).or_default()
    }

    /// OCCT TopOpeBRepDS_DataStructure::ShapeInterferences(S, FindKeep).
    /// Interference list stays on the facade; the shape lookup routes
    /// through BOPDS `DS::index`.
    pub fn shape_interferences_of(&self, s: &Shape) -> &[TopOpeBRepDSInterference] {
        let i = self.bopds.index(s);
        if i >= 0 {
            self.shape_interferences(i as i32 + 1)
        } else {
            &[]
        }
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeShapeInterferences(S).
    /// Interference list stays on the facade; shape registration routes
    /// through `add_shape` (BOPDS).
    pub fn change_shape_interferences_of(&mut self, s: &Shape) -> &mut Vec<TopOpeBRepDSInterference> {
        let i = self.add_shape(s);
        self.change_shape_interferences(i)
    }

    /// OCCT TopOpeBRepDS_DataStructure::CurveInterferences(I) L384-496
    /// (const).  Stays on the facade (D6: ChFi3d-private quadruplet list).
    pub fn curve_interferences(&self, i: i32) -> &[TopOpeBRepDSInterference] {
        self.curve_interferences
            .get(&i)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeCurveInterferences(I).
    /// Stays on the facade (D6: ChFi3d-private quadruplet list).
    pub fn change_curve_interferences(&mut self, i: i32) -> &mut Vec<TopOpeBRepDSInterference> {
        self.curve_interferences.entry(i).or_default()
    }

    /// OCCT TopOpeBRepDS_DataStructure::SurfaceInterferences(I) L384-496
    /// (const).  Stays on the facade (D6: ChFi3d-private quadruplet list).
    pub fn surface_interferences(&self, i: i32) -> &[TopOpeBRepDSInterference] {
        self.surface_interferences
            .get(&i)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// OCCT TopOpeBRepDS_DataStructure::ChangeSurfaceInterferences(I).
    /// Stays on the facade (D6: ChFi3d-private quadruplet list).
    pub fn change_surface_interferences(&mut self, i: i32) -> &mut Vec<TopOpeBRepDSInterference> {
        self.surface_interferences.entry(i).or_default()
    }
}
