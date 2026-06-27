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
) -> Vec<WireSegmentTopoDS> {
    let face_ref = ShapeRef::new(face_idx);
    segments.iter().map(|seg| {
        let (edge_ref, source) = match &seg.source {
            WireEdgeSource::DsEdge(ei) => {
                (ShapeRef::new(*ei), WireEdgeSourceTopoDS::DsEdge(ShapeRef::new(*ei)))
            }
            WireEdgeSource::IntersectionCurve(ci) => {
                (ShapeRef::new(*ci), WireEdgeSourceTopoDS::IntersectionCurve(ShapeRef::new(*ci)))
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
            is_seam: seg.is_seam,
            tangent_start: seg.tangent_start,
            tangent_end: seg.tangent_end,
            first_pcurve: seg.first_pcurve.clone(),
            second_pcurve: seg.second_pcurve.clone(),
            t_range: seg.t_range,
        }
    }).collect()
}
