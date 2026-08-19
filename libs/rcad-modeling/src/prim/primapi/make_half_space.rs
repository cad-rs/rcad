// OCCT BRepPrimAPI_MakeHalfSpace 1:1 translation.
// OCCT ref: BRepPrimAPI_MakeHalfSpace.cxx L40-210.
//
// Builds an open solid whose single shell is the input face.  When the
// reference point lies on the outside of the face (isOutside), the shell is
// reversed so the solid interior faces the reference point.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::base::geom_api::project::closest_point_on_surface;
use rcad_kernel::geom::SurfaceEval;
use rcad_kernel::topods::{Orientation, Shape};

/// OCCT getNormalOnFace + FindExtrema (BRepPrimAPI_MakeHalfSpace.cxx L41-133):
/// project the reference point onto the face surface and return the closest
/// point plus the surface normal (reversed for a REVERSED face).
fn find_extrema(brep: &BRep, face: &Shape, ref_pnt: DVec3) -> Option<(DVec3, DVec3)> {
    let fd = match brep.tshapes.get(face.index).map(|ts| ts.as_ref())? {
        rcad_kernel::topods::TShape::Face(fd) => fd,
        _ => return None,
    };
    let surface = fd.surface.clone()?;
    let proj = closest_point_on_surface(&surface, ref_pnt, 23);
    if !proj.point.is_finite() || !proj.params.0.is_finite() || !proj.params.1.is_finite() {
        return None;
    }
    let mut normal = surface.normal_at(proj.params.0, proj.params.1);
    if face.orientation == Orientation::Reversed {
        normal = -normal;
    }
    Some((proj.point, normal))
}

/// OCCT BRepPrimAPI_MakeHalfSpace(face, refPnt) (L152-183): open solid with a
/// single-face shell; the shell is reversed when isOutside is true (L138-150).
pub fn make_half_space_brep(
    brep: &BRep,
    face: &Shape,
    ref_pnt: DVec3,
) -> Result<BRep, crate::BuildError> {
    let (min_pnt, normal) = find_extrema(brep, face, ref_pnt)
        .ok_or(crate::BuildError::DegenerateGeometry(
            "half-space: no face projection",
        ))?;
    let to_reverse = (ref_pnt - min_pnt).dot(normal) > 0.0;
    let mut result = BRep::new();
    let shell_face = if to_reverse {
        Shape {
            orientation: Orientation::Reversed,
            ..face.clone()
        }
    } else {
        face.clone()
    };
    let shell = result.add_tshell(vec![shell_face]);
    result.add_tsolid(vec![shell]);
    Ok(result)
}
