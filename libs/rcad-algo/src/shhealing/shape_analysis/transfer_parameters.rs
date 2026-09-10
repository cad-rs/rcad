//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_TransferParameters`
//! (`ShapeAnalysis_TransferParameters.hxx` L17-96 + `.cxx` L1-153).
//!
//! Transfers the parameters of an edge from its 3D curve to its pcurve and
//! back by the linear mapping (shift + scale) computed in `Init`.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — the shape-mutating members (`TransferRange`
//!    through `ShapeBuild_Edge::CopyRanges`) read and mutate the TShape
//!    graph through `rcad_kernel::BRep` (the edge.rs bridge #1).
//! 2. `occ::handle<NCollection_HSequence<double>>` -> `Vec<f64>`.
//! 3. `TopoDS_Edge` / `TopoDS_Face` -> `rcad_kernel::topo_shape::Shape`
//!    (the null face of the OCCT free-edge call -> `Shape::null()`).
//!
//! The fields are `pub(crate)` because `ShapeAnalysis_TransferParametersProj`
//! (the OCCT subclass) reads and writes them.

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::BRep;

use super::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_build::edge::ShapeBuildEdge;

/// OCCT gp::Resolution() (gp.hxx: the smallest positive coordinate, 1e-12).
const GP_RESOLUTION: f64 = 1e-12;

/// OCCT ShapeAnalysis_TransferParameters (hxx L24-96).
pub struct ShapeAnalysisTransferParameters {
    /// OCCT `protected: double myFirst`.
    pub(crate) my_first: f64,
    /// OCCT `protected: double myLast`.
    pub(crate) my_last: f64,
    /// OCCT `protected: TopoDS_Edge myEdge`.
    pub(crate) my_edge: Shape,
    /// OCCT `protected: double myMaxTolerance` (uninitialised in the OCCT
    /// default constructor; the rcad zero initialiser stands in).
    pub(crate) my_max_tolerance: f64,
    /// OCCT `protected: double myShift`.
    pub(crate) my_shift: f64,
    /// OCCT `protected: double myScale`.
    pub(crate) my_scale: f64,
    /// OCCT `protected: double myFirst2d`.
    pub(crate) my_first2d: f64,
    /// OCCT `protected: double myLast2d`.
    pub(crate) my_last2d: f64,
    /// OCCT `protected: TopoDS_Face myFace`.
    pub(crate) my_face: Shape,
}

impl Default for ShapeAnalysisTransferParameters {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisTransferParameters {
    /// OCCT ShapeAnalysis_TransferParameters() (cxx L29-33): myScale = 1.,
    /// myShift = 0. (the remaining members keep the rcad zero initialiser).
    pub fn new() -> Self {
        ShapeAnalysisTransferParameters {
            my_first: 0.0,
            my_last: 0.0,
            my_edge: Shape::null(),
            my_max_tolerance: 0.0,
            my_shift: 0.0,
            my_scale: 1.0,
            my_first2d: 0.0,
            my_last2d: 0.0,
            my_face: Shape::null(),
        }
    }

    /// OCCT ShapeAnalysis_TransferParameters(E, F) (cxx L37-41).
    pub fn new_edge_face(brep: &mut BRep, e: &Shape, face: &Shape) -> Self {
        let mut this = Self::new();
        this.init(brep, e, face);
        this
    }

    /// OCCT Init(E, F) (cxx L45-75); the face parameter is `face` (the OCCT
    /// `F`), because the body declares the local `f` (the 3d curve first).
    pub fn init(&mut self, brep: &mut BRep, e: &Shape, face: &Shape) {
        self.my_scale = 1.;
        self.my_shift = 0.;
        let mut l = 0.0;
        let mut f = 0.0;
        let mut l2d = 0.0;
        let mut f2d = 0.0;
        self.my_edge = e.clone();
        let sae = ShapeAnalysisEdge::new();
        let mut curve3d: Option<rcad_kernel::geom::Curve3> = None; // = BRep_Tool::Curve (E,f,l);
        sae.curve3d(brep, e, &mut curve3d, &mut f, &mut l, false);
        self.my_first = f;
        self.my_last = l;
        let mut curve2d: Option<rcad_kernel::geom::Curve2d> = None; // = BRep_Tool::CurveOnSurface (E, F, f2d,l2d);
        // ShapeAnalysis_Edge sae;
        if !face.is_null() {
            // process free edges
            sae.pcurve_face(brep, e, face, &mut curve2d, &mut f2d, &mut l2d, false);
        }
        self.my_first2d = f2d;
        self.my_last2d = l2d;
        self.my_face = face.clone();
        if curve3d.is_none() || curve2d.is_none() {
            return;
        }

        let ln2d = l2d - f2d;
        let ln3d = l - f;
        self.my_scale = if ln3d <= GP_RESOLUTION {
            1.
        } else {
            ln2d / ln3d
        };
        self.my_shift = f2d - f * self.my_scale;
    }

    /// OCCT SetMaxTolerance(maxtol) (cxx L79-82).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.my_max_tolerance = maxtol;
    }

    /// OCCT Perform(Params, To2d) (cxx L86-96).
    pub fn perform_params(&self, params: &[f64], to2d: bool) -> Vec<f64> {
        let mut res = Vec::new();
        for &param in params {
            res.push(self.perform(param, to2d));
        }
        res
    }

    /// OCCT Perform(Param, To2d) (cxx L100-112).
    pub fn perform(&self, param: f64, to2d: bool) -> f64 {
        let new_param;
        if to2d {
            new_param = self.my_shift + param * self.my_scale;
        } else {
            new_param = -self.my_shift / self.my_scale + param * 1. / self.my_scale;
        }
        new_param
    }

    /// OCCT TransferRange(newEdge, prevPar, currPar, Is2d) (cxx L116-146).
    pub fn transfer_range(
        &self,
        brep: &mut BRep,
        new_edge: &Shape,
        prev_par: f64,
        curr_par: f64,
        is2d: bool,
    ) {
        let sbe = ShapeBuildEdge;
        if is2d {
            let span2d = self.my_last2d - self.my_first2d;
            let tmp1;
            let tmp2;
            if prev_par > curr_par {
                tmp1 = curr_par;
                tmp2 = prev_par;
            } else {
                tmp1 = prev_par;
                tmp2 = curr_par;
            }
            let alpha = (tmp1 - self.my_first2d) / span2d;
            let beta = (tmp2 - self.my_first2d) / span2d;
            sbe.copy_ranges(brep, new_edge, &self.my_edge, alpha, beta);
        } else {
            let alpha = (prev_par - self.my_first) / (self.my_last - self.my_first);
            let beta = (curr_par - self.my_first) / (self.my_last - self.my_first);
            sbe.copy_ranges(brep, new_edge, &self.my_edge, alpha, beta);
        }
    }

    /// OCCT IsSameRange() (cxx L150-153).
    pub fn is_same_range(&self) -> bool {
        self.my_shift == 0. && self.my_scale == 1.
    }
}
