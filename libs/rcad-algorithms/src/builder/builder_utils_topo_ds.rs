/// Convert Vec<WireSegment> to Vec<WireSegmentTopoDS>.
///
/// Maps DS indices to ShapeRef handles so the output feeds into
/// walk_path_extract_wires_topoDS. The DSAsBRep adaptor maps back
/// by ShapeRef.index when queried.
use rcad_kernel::topods::{Orientation, ShapeRef};
use crate::bopds::ds::DS;
use super::types::{WireSegment, WireEdgeSource, WireOrientation, WireSegmentTopoDS, WireEdgeSourceTopoDS};

/// Convert existing WireSegments to WireSegmentTopoDS for the BRepTool pipeline.
pub(crate) fn segments_to_topo_ds(
    segments: &[WireSegment],
    _ds: &DS,
    face_idx: usize,
    face_refs: &[ShapeRef],
    ic_edge_map: &[Option<ShapeRef>],
) -> Vec<WireSegmentTopoDS> {
    let face_ref = face_refs[face_idx];
    let e_base = _ds.vertices.len();
    segments.iter().map(|seg| {
        let (edge_ref, source) = match &seg.source {
            WireEdgeSource::DsEdge(ei) => {
                (ShapeRef::new(e_base + *ei), WireEdgeSourceTopoDS::DsEdge(ShapeRef::new(e_base + *ei)))
            }
            WireEdgeSource::IntersectionCurve(ci) => {
                if let Some(Some(edge_ref)) = ic_edge_map.get(*ci) {
                    // IC mapped to DSEdge or retained IC TEdge (A2 dedup)
                    (*edge_ref, WireEdgeSourceTopoDS::IntersectionCurve(*edge_ref))
                } else {
                    // Fallback: should not happen
                    (ShapeRef::new(0), WireEdgeSourceTopoDS::IntersectionCurve(ShapeRef::new(0)))
                }
            }
            WireEdgeSource::SeamEdge => {
                (ShapeRef::new(0), WireEdgeSourceTopoDS::SeamEdge)
            }
        };
        let orientation = match seg.orientation {
            WireOrientation::Forward => Orientation::Forward,
            WireOrientation::Reversed => Orientation::Reversed,
            WireOrientation::Internal => Orientation::Internal,
            WireOrientation::External => Orientation::External,
        };
        WireSegmentTopoDS {
            edge: edge_ref,
            face: face_ref,
            start_vertex: ShapeRef::new(seg.start_vertex),
            end_vertex: ShapeRef::new(seg.end_vertex),
            source,
            orientation,
            is_closed_on_face: seg.is_closed_on_face,
            first_pcurve: seg.first_pcurve.clone(),
            second_pcurve: seg.second_pcurve.clone(),
            t_range: seg.t_range,
        }
    }).collect()
}

/// Convert WireSegmentTopoDS back to WireSegment (for emit_wire_face compatibility).
pub(crate) fn topo_ds_to_segments(
    topo: &[WireSegmentTopoDS],
) -> Vec<WireSegment> {
    topo.iter().map(|s| {
        let src = match &s.source {
            WireEdgeSourceTopoDS::DsEdge(e) => WireEdgeSource::DsEdge(e.index),
            WireEdgeSourceTopoDS::IntersectionCurve(c) => WireEdgeSource::IntersectionCurve(c.index),
            WireEdgeSourceTopoDS::SeamEdge => WireEdgeSource::SeamEdge,
        };
        let ori = match s.orientation {
            rcad_kernel::topods::Orientation::Forward => WireOrientation::Forward,
            rcad_kernel::topods::Orientation::Reversed => WireOrientation::Reversed,
            rcad_kernel::topods::Orientation::Internal => WireOrientation::Internal,
            rcad_kernel::topods::Orientation::External => WireOrientation::External,
        };
        WireSegment {
            start_vertex: s.start_vertex.index,
            end_vertex: s.end_vertex.index,
            source: src,
            orientation: ori,
            is_closed_on_face: s.is_closed_on_face,
            first_pcurve: s.first_pcurve.clone(),
            second_pcurve: s.second_pcurve.clone(),
            t_range: s.t_range,
        }
    }).collect()
}
