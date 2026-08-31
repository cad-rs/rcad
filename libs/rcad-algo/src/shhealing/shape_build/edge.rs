//! OCCT ShapeBuild_Edge (TKShHealing ShapeBuild package).
//!
//! 1:1 translation of the methods of `ShapeBuild_Edge.cxx` needed by the
//! ShapeFix stack. The struct is the OCCT class; methods take the owning
//! `BRep` pool because rcad TShape data lives there (OCCT reaches it through
//! the `TopoDS_Shape` handle).

use rcad_kernel::topo::topods::{BRep, Shape};

/// OCCT ShapeBuild_Edge.
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeBuildEdge;

impl ShapeBuildEdge {
    /// OCCT ShapeBuild_Edge::CopyRanges (ShapeBuild_Edge.cxx L206-296): copies
    /// ranges for curve3d and all common pcurves from `fromedge` into
    /// `toedge`, scaled to `[first + alpha * len, first + beta * len]`.
    /// ReShape::applyImpl calls it with the defaults alpha = 0, beta = 1.
    ///
    /// Architecture note: OCCT walks two parallel representation lists and
    /// matches entries by (surface, location); rcad stores the 3D range in
    /// `TEdgeData::range` and pcurve ranges in `TEdgeData::pcurves` keyed by
    /// `(face ptr, location)`, so the matching key is the pcurve map key.
    pub fn copy_ranges(
        &self,
        brep: &mut BRep,
        toedge: &Shape,
        fromedge: &Shape,
        alpha: f64,
        beta: f64,
    ) {
        let from = brep.edge(fromedge.clone()).clone();

        // fromGC IsCurve3D branch: skip when the 3D curve is null.
        if from.curve.is_some() {
            let first = from.range[0];
            let last = from.range[1];
            let len = last - first;
            // toGC IsCurve3D match: set the 3D representation range.
            let to = brep.edge_mut(toedge.clone());
            to.range = [first + alpha * len, first + beta * len];
        }

        // pcurve representations: skip when the pcurve is null (rcad map
        // entries always carry a curve); match by (surface, location) key.
        for (key, (_, pfirst, plast)) in from.pcurves.iter() {
            let first = *pfirst;
            let last = *plast;
            let len = last - first;
            let new_first = first + alpha * len;
            let new_last = first + beta * len;
            let to = brep.edge_mut(toedge.clone());
            if let Some(entry) = to.pcurves.get_mut(key) {
                entry.1 = new_first;
                entry.2 = new_last;
            }
        }
    }
}
