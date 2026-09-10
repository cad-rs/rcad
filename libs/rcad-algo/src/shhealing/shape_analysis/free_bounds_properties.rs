//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_FreeBoundsProperties`
//! (`ShapeAnalysis_FreeBoundsProperties.hxx` L17-160 + `.lxx` L19-120 +
//! `.cxx` L1-423).
//!
//! Builds the free bounds of the shape (through ShapeAnalysis_FreeBounds)
//! and on each free bound computes its properties: area and perimeter of
//! the contour, ratio of average length to average width, average width,
//! and the notches with their maximum widths.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — `BRep_Builder` reads and mutates the TShape
//!    graph through `rcad_kernel::BRep` (the edge.rs bridge #1); the OCCT
//!    static ContourProperties receives the pool alongside (bridge #1).
//! 2. `occ::handle<NCollection_HSequence<...>>` -> `Vec<...>`; the OCCT
//!    shared `handle<ShapeAnalysis_FreeBoundData>&` argument (mutated
//!    through the handle) -> an owned value cloned back into the sequence
//!    at the call sites (the value-model bridge).
//! 3. `ShapeExtend_WireData` -> `shape_extend::wire_data::WireData` (the
//!    W1-3 translation); `ShapeExtend_Explorer` ->
//!    `shape_extend::explorer::ShapeExtendExplorer` (the W1-3
//!    translation).  OCCT `saw->Load(wdt)` shares the WireData handle; the
//!    rcad value model rebuilds the same deterministic WireData for the
//!    analyzer (bridge #2).
//! 4. `ShapeAnalysis_FreeBounds` / `ShapeAnalysis_Wire` / `ShapeAnalysis_Edge`
//!    -> the 1:1 translations in this module tree.
//! 5. `GeomAPI_ProjectPointOnCurve ppc(pntCurr, c3d2, p1, p2)` ->
//!    `base::geom_api::project::closest_point_on_curve_range` (the landed
//!    GeomAPI projection port; `ppc.NbPoints()` is never 0 for a valid
//!    curve, so the `(ppc.NbPoints() ? ppc.LowerDistance() : 0)` pick maps
//!    to the returned distance).
//!
//! GAP carriers (untranslated other-package dependencies):
//! - `BRepTools_WireExplorer` (TKTopAlgo BRepTools, untranslated): the
//!   local [`BRepToolsWireExplorer`] re-host keeps the Init/More/Next/
//!   Current walk over the wire's stored edge order — the connected-order
//!   reordering of the OCCT explorer is the documented reduction (the same
//!   reduction as the offset-module carrier).  GAP: closes with the
//!   BRepTools batch.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval as _};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType};

use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::free_bound_data::ShapeAnalysisFreeBoundData;
use crate::shhealing::shape_analysis::free_bounds::ShapeAnalysisFreeBounds;
use crate::shhealing::shape_analysis::wire::ShapeAnalysisWire;
use crate::shhealing::shape_build::brep_tool::{builder_add, raw_subshapes};
use crate::shhealing::shape_extend::explorer::ShapeExtendExplorer;
use crate::shhealing::shape_extend::wire_data::WireData;

/// OCCT `#define NbControl 23` (cxx L41).
const NB_CONTROL: i32 = 23;

/// GAP re-host of OCCT `BRepTools_WireExplorer` (see the module GAP note):
/// Init/More/Next/Current over the wire's stored edge order.
pub struct BRepToolsWireExplorer {
    my_edges: Vec<Shape>,
    my_index: usize,
}

impl BRepToolsWireExplorer {
    /// OCCT BRepTools_WireExplorer(wire) (the direct-wire constructor; the
    /// face argument is unused by the stored-order walk).
    pub fn new(brep: &BRep, the_w: &Shape) -> Self {
        BRepToolsWireExplorer {
            my_edges: raw_subshapes(brep, the_w),
            my_index: 0,
        }
    }

    /// OCCT More().
    pub fn more(&self) -> bool {
        self.my_index < self.my_edges.len()
    }

    /// OCCT Current().
    pub fn current(&self) -> Shape {
        self.my_edges[self.my_index].clone()
    }

    /// OCCT Next().
    pub fn next(&mut self) {
        self.my_index += 1;
    }
}

/// OCCT static ContourProperties (cxx L43-87): the area and the length of
/// the contour from the NbControl-point polygon on each edge.
fn contour_properties(brep: &mut BRep, wire: &Shape, countour_area: &mut f64, countour_length: &mut f64) {
    let mut nbe = 0;
    let mut length = 0.0;
    let mut area = DVec3::ZERO;
    let mut prev = DVec3::ZERO;
    let mut cont = DVec3::ZERO;

    let mut exp = BRepToolsWireExplorer::new(brep, wire);
    while exp.more() {
        let edge = exp.current();
        nbe += 1;

        let sae = ShapeAnalysisEdge::new();
        let mut first = 0.0;
        let mut last = 0.0;
        let mut c3d: Option<Curve3> = None;
        if !sae.curve3d(brep, &edge, &mut c3d, &mut first, &mut last, true) {
            exp.next();
            continue;
        }
        let c3d = match c3d {
            Some(c) => c,
            None => {
                exp.next();
                continue;
            }
        };

        let mut ibeg = 0;
        if nbe == 1 {
            let pnt_ini = c3d.point_at(first);
            prev = pnt_ini;
            cont = prev;
            ibeg = 1;
        }

        for i in ibeg..NB_CONTROL {
            let prm = ((NB_CONTROL - 1 - i) as f64 * first + i as f64 * last)
                / (NB_CONTROL - 1) as f64;
            let pnt_curr = c3d.point_at(prm);
            let curr = pnt_curr;
            let delta = curr - prev;
            length += delta.length();
            area += curr.cross(prev);
            prev = curr;
        }
        exp.next();
    }

    area += cont.cross(prev);
    *countour_area = area.length() / 2.0;
    *countour_length = length;
}

/// OCCT ShapeAnalysis_FreeBoundsProperties (hxx L32-160).
pub struct ShapeAnalysisFreeBoundsProperties {
    /// OCCT `myShape`.
    my_shape: Shape,
    /// OCCT `myTolerance`.
    my_tolerance: f64,
    /// OCCT `mySplitClosed`.
    my_split_closed: bool,
    /// OCCT `mySplitOpen`.
    my_split_open: bool,
    /// OCCT `myClosedFreeBounds`.
    my_closed_free_bounds: Vec<ShapeAnalysisFreeBoundData>,
    /// OCCT `myOpenFreeBounds`.
    my_open_free_bounds: Vec<ShapeAnalysisFreeBoundData>,
}

impl Default for ShapeAnalysisFreeBoundsProperties {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisFreeBoundsProperties {
    /// OCCT ShapeAnalysis_FreeBoundsProperties() (cxx L91-96).
    pub fn new() -> Self {
        ShapeAnalysisFreeBoundsProperties {
            my_shape: Shape::null(),
            my_tolerance: 0.,
            my_split_closed: false,
            my_split_open: false,
            my_closed_free_bounds: Vec::new(),
            my_open_free_bounds: Vec::new(),
        }
    }

    /// OCCT ShapeAnalysis_FreeBoundsProperties(shape, tolerance,
    /// splitclosed = false, splitopen = false) (cxx L104-112).
    pub fn new_with_tolerance(
        brep: &mut BRep,
        shape: &Shape,
        tolerance: f64,
        splitclosed: bool,
        splitopen: bool,
    ) -> Self {
        let mut this = Self::new();
        this.init(brep, shape, tolerance, splitclosed, splitopen);
        this
    }

    /// OCCT ShapeAnalysis_FreeBoundsProperties(shape, splitclosed = false,
    /// splitopen = false) (cxx L120-128).
    pub fn new_with_split(
        brep: &mut BRep,
        shape: &Shape,
        splitclosed: bool,
        splitopen: bool,
    ) -> Self {
        let mut this = Self::new();
        this.init_actual(brep, shape, splitclosed, splitopen);
        this
    }

    /// OCCT Init(shape, tolerance, splitclosed, splitopen) (cxx L136-143).
    pub fn init(
        &mut self,
        brep: &mut BRep,
        shape: &Shape,
        tolerance: f64,
        splitclosed: bool,
        splitopen: bool,
    ) {
        self.init_actual(brep, shape, splitclosed, splitopen);
        self.my_tolerance = tolerance;
    }

    /// OCCT Init(shape, splitclosed, splitopen) (cxx L151-158).
    pub fn init_actual(
        &mut self,
        _brep: &mut BRep,
        shape: &Shape,
        splitclosed: bool,
        splitopen: bool,
    ) {
        self.my_shape = shape.clone();
        self.my_split_closed = splitclosed;
        self.my_split_open = splitopen;
    }

    /// OCCT Perform() (cxx L176-183): builds and analyzes free bounds of
    /// the shape.
    pub fn perform(&mut self, brep: &mut BRep) -> bool {
        let mut result = false;
        result |= self.dispatch_bounds(brep);
        result |= self.check_notches(brep, 0.0);
        result |= self.check_contours(brep, 0.0);
        result
    }

    /// OCCT IsLoaded() (lxx L33-35).
    pub fn is_loaded(&self) -> bool {
        !self.my_shape.is_null()
    }

    /// OCCT Shape() (lxx L24-26).
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT Tolerance() (lxx L39-41).
    pub fn tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT NbFreeBounds() (lxx L46-48).
    pub fn nb_free_bounds(&self) -> i32 {
        self.my_closed_free_bounds.len() as i32 + self.my_open_free_bounds.len() as i32
    }

    /// OCCT NbClosedFreeBounds() (lxx L53-55).
    pub fn nb_closed_free_bounds(&self) -> i32 {
        self.my_closed_free_bounds.len() as i32
    }

    /// OCCT NbOpenFreeBounds() (lxx L60-62).
    pub fn nb_open_free_bounds(&self) -> i32 {
        self.my_open_free_bounds.len() as i32
    }

    /// OCCT ClosedFreeBounds() (lxx L66-70).
    pub fn closed_free_bounds(&self) -> &[ShapeAnalysisFreeBoundData] {
        &self.my_closed_free_bounds
    }

    /// OCCT OpenFreeBounds() (lxx L72-76).
    pub fn open_free_bounds(&self) -> &[ShapeAnalysisFreeBoundData] {
        &self.my_open_free_bounds
    }

    /// OCCT ClosedFreeBound(index) (lxx L78-82).
    pub fn closed_free_bound(&self, index: i32) -> &ShapeAnalysisFreeBoundData {
        &self.my_closed_free_bounds[(index - 1) as usize]
    }

    /// OCCT OpenFreeBound(index) (lxx L84-88).
    pub fn open_free_bound(&self, index: i32) -> &ShapeAnalysisFreeBoundData {
        &self.my_open_free_bounds[(index - 1) as usize]
    }

    /// OCCT DispatchBounds() (cxx L187-231).
    pub fn dispatch_bounds(&mut self, brep: &mut BRep) -> bool {
        if !self.is_loaded() {
            return false;
        }

        let (tmp_closed_bounds, tmp_open_bounds) = if self.my_tolerance > 0. {
            let safb = ShapeAnalysisFreeBounds::new_forecast(
                brep,
                &self.my_shape,
                self.my_tolerance,
                self.my_split_closed,
                self.my_split_open,
            );
            (safb.get_closed_wires().clone(), safb.get_open_wires().clone())
        } else {
            // OCCT L203: ShapeAnalysis_FreeBounds safb(myShape,
            // mySplitClosed, mySplitOpen) — the (shape, splitclosed,
            // splitopen) constructor is the ACTUAL bounds form (the
            // W1-5 new_actual with checkinternaledges = false).
            let safb = ShapeAnalysisFreeBounds::new_actual(
                brep,
                &self.my_shape,
                self.my_split_closed,
                self.my_split_open,
                false,
            );
            (safb.get_closed_wires().clone(), safb.get_open_wires().clone())
        };

        let shexpl = ShapeExtendExplorer;
        let tmp_seq = shexpl.seq_from_compound(brep, &tmp_closed_bounds, false);
        for i in 1..=tmp_seq.len() as i32 {
            let wire = tmp_seq[(i - 1) as usize].clone();
            // OCCT L214: TopoDS::Wire(tmpSeq->Value(i)) — the element is a
            // wire; the rcad walk keeps the shape (the FreeBoundData
            // stores it untyped).
            let mut fb_data = ShapeAnalysisFreeBoundData::new();
            fb_data.set_free_bound(&wire);
            self.my_closed_free_bounds.push(fb_data);
        }

        let tmp_seq2 = shexpl.seq_from_compound(brep, &tmp_open_bounds, false);
        for i in 1..=tmp_seq2.len() as i32 {
            let wire = tmp_seq2[(i - 1) as usize].clone();
            let mut fb_data = ShapeAnalysisFreeBoundData::new();
            fb_data.set_free_bound(&wire);
            self.my_open_free_bounds.push(fb_data);
        }

        true
    }

    /// OCCT CheckNotches(prec = 0.0) (cxx L235-250).
    pub fn check_notches(&mut self, brep: &mut BRep, prec: f64) -> bool {
        // OCCT walks the shared handles in place; the rcad value model
        // clones the entry out, mutates it and writes it back (bridge #2).
        for i in 0..self.my_closed_free_bounds.len() {
            let mut fb_data = self.my_closed_free_bounds[i].clone();
            self.check_notches_of(brep, &mut fb_data, prec);
            self.my_closed_free_bounds[i] = fb_data;
        }
        for i in 0..self.my_open_free_bounds.len() {
            let mut fb_data = self.my_open_free_bounds[i].clone();
            self.check_notches_of(brep, &mut fb_data, prec);
            self.my_open_free_bounds[i] = fb_data;
        }

        true
    }

    /// OCCT CheckNotches(fbData, prec = 0.0) (cxx L254-273).
    pub fn check_notches_of(
        &self,
        brep: &mut BRep,
        fb_data: &mut ShapeAnalysisFreeBoundData,
        prec: f64,
    ) -> bool {
        let swd = WireData::new_from_wire(brep, &fb_data.free_bound(), false, true);
        if swd.nb_edges() > 1 {
            for j in 1..=swd.nb_edges() {
                let mut notch = Shape::null();
                let mut d_max = 0.0;
                if self.check_notches_wire(
                    brep,
                    &fb_data.free_bound(),
                    j,
                    &mut notch,
                    &mut d_max,
                    prec,
                ) {
                    fb_data.add_notch(&notch, d_max);
                }
            }
        }

        true
    }

    /// OCCT CheckContours(prec = 0.0) (cxx L277-293).
    pub fn check_contours(&mut self, brep: &mut BRep, prec: f64) -> bool {
        let mut status = false;
        for i in 0..self.my_closed_free_bounds.len() {
            let mut fb_data = self.my_closed_free_bounds[i].clone();
            status |= self.fill_properties(brep, &mut fb_data, prec);
            self.my_closed_free_bounds[i] = fb_data;
        }
        for i in 0..self.my_open_free_bounds.len() {
            let mut fb_data = self.my_open_free_bounds[i].clone();
            status |= self.fill_properties(brep, &mut fb_data, prec);
            self.my_open_free_bounds[i] = fb_data;
        }

        status
    }

    /// OCCT CheckNotches(wire, num, notch, distMax, prec = 0.0)
    /// (cxx L297-387).
    #[allow(clippy::too_many_arguments)]
    pub fn check_notches_wire(
        &self,
        brep: &mut BRep,
        wire: &Shape,
        num: i32,
        notch: &mut Shape,
        dist_max: &mut f64,
        _prec: f64,
    ) -> bool {
        let tol = self.my_tolerance.max(CONFUSION);
        let wdt = WireData::new_from_wire(brep, wire, false, true);
        // OCCT L305-306: BRep_Builder B; B.MakeWire(notch).
        *notch = brep.add_twire(Vec::new());

        if (num <= 0) || (num > wdt.nb_edges()) {
            return false;
        }

        let n1 = if num > 0 { num } else { wdt.nb_edges() };
        let mut n2 = if n1 < wdt.nb_edges() { n1 + 1 } else { 1 };

        let e1 = wdt.edge(n1);
        builder_add(brep, notch, &e1);

        let mut saw = ShapeAnalysisWire::new();
        // OCCT L319-321: saw->Load(wdt); saw->SetPrecision(myTolerance) —
        // the shared handle (bridge #3: the analyzer WireData is rebuilt
        // deterministically).
        saw.load_wire_data(WireData::new_from_wire(brep, wire, false, true));
        saw.set_precision(self.my_tolerance);
        if saw.check_small(brep, n2, tol) {
            let e = wdt.edge(n2);
            builder_add(brep, notch, &e);
            n2 = if n2 < wdt.nb_edges() { n2 + 1 } else { 1 };
        }

        let e2 = wdt.edge(n2);
        builder_add(brep, notch, &e2);

        let sae = ShapeAnalysisEdge::new();
        let mut first1 = 0.0;
        let mut last1 = 0.0;
        let mut first2 = 0.0;
        let mut last2 = 0.0;
        let mut c3d1: Option<Curve3> = None;
        let mut c3d2: Option<Curve3> = None;
        // szv#4:S4163:12Mar99 optimized
        if !sae.curve3d(brep, &e1, &mut c3d1, &mut first1, &mut last1, true)
            || !sae.curve3d(brep, &e2, &mut c3d2, &mut first2, &mut last2, true)
        {
            return false;
        }
        let c3d1 = match c3d1 {
            Some(c) => c,
            None => return false,
        };
        let c3d2 = match c3d2 {
            Some(c) => c,
            None => return false,
        };

        // OCCT L342-343: c3d1->D1(Last1, pnt, vec1); c3d2->D1(First2, pnt,
        // vec2) — the points are overwritten and never read.
        let mut vec1 = c3d1.derivative_at(last1);
        let mut vec2 = c3d2.derivative_at(first2);
        if e1.orientation == Orientation::Reversed {
            vec1 = -vec1;
        }
        if e2.orientation == Orientation::Reversed {
            vec2 = -vec2;
        }

        let angl = vec1.angle_between(vec2).abs();
        if angl > 0.95 * std::f64::consts::PI {
            *dist_max = 0.0;
            for i in 0..NB_CONTROL {
                let prm = ((NB_CONTROL - 1 - i) as f64 * first1 + i as f64 * last1)
                    / (NB_CONTROL - 1) as f64;
                let pnt_curr = c3d1.point_at(prm);

                let (p1, p2) = if first2 < last2 {
                    (first2, last2)
                } else {
                    (last2, first2)
                };

                // szv#4:S4163:12Mar99 warning
                // OCCT L375-376: GeomAPI_ProjectPointOnCurve ppc(pntCurr,
                // c3d2, p1, p2); newDist = (ppc.NbPoints() ?
                // ppc.LowerDistance() : 0) (bridge #5).
                let ppc = closest_point_on_curve_range(&c3d2, pnt_curr, p1, p2, 24);
                let new_dist = ppc.distance;
                if new_dist > *dist_max {
                    *dist_max = new_dist;
                }
            }

            return true;
        }

        false
    }

    /// OCCT FillProperties(fbData, prec = 0.0) (cxx L391-423).
    pub fn fill_properties(
        &self,
        brep: &mut BRep,
        fb_data: &mut ShapeAnalysisFreeBoundData,
        _prec: f64,
    ) -> bool {
        let mut area = 0.0;
        let mut length = 0.0;
        contour_properties(brep, &fb_data.free_bound(), &mut area, &mut length);

        let mut r = 0.0;
        let mut aver = 0.0;

        if length != 0. {
            // szv#4:S4163:12Mar99 anti-exception
            let k = area / (length * length); // szv#4:S4163:12Mar99
            // szv#4:S4163:12Mar99 optimized
            if k != 0. {
                // szv#4:S4163:12Mar99 anti-exception
                let aux = 1. - 16. * k;
                if aux >= 0. {
                    r = (1. + aux.sqrt()) / (8. * k);
                    aver = length / (2. * r);
                    r -= 1.;
                }
            }
        }

        fb_data.set_area(area);
        fb_data.set_perimeter(length);
        fb_data.set_ratio(r);
        fb_data.set_width(aver);

        true
    }
}

// The ShapeType import is part of the module surface (the topexp explorers
// of the sibling classes); kept out of this file's walk.
#[allow(unused)]
fn _shape_type_marker(_t: ShapeType) {}
