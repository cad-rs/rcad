use rcad_kernel::geom::*;
use rcad_kernel::BRep;

/// Populates `brep.geom` with analytic geometry for a box BRep.
///
/// After this call, every edge has a `Curve3::Line` and every face has a `Surface3::Plane`.
/// Precondition: brep was created by `BRep::from_primitive(Box{..})`.
pub fn populate_box_geom(brep: &mut BRep) {
    let geom = &mut brep.geom;
    geom.curves.clear();
    geom.edge_curve.clear();
    geom.surfaces.clear();
    geom.face_surface.clear();

    // Edges → Line3
    for edge in &brep.edges {
        let p0 = brep.vertices[edge.start].point;
        let p1 = brep.vertices[edge.end].point;
        let dir = (p1 - p0).normalize();
        let curve_idx = geom.curves.len();
        geom.curves.push(Curve3::Line(Line3 {
            origin: p0,
            direction: dir,
        }));
        geom.edge_curve.push(Some(curve_idx));
    }

    // Faces → Plane
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let origin = brep.vertices[face.triangles[0][0]].point;
                let surf_idx = geom.surfaces.len();
                geom.surfaces.push(Surface3::Plane(Plane {
                    origin,
                    normal: face.normal,
                }));
                geom.face_surface.push(Some(surf_idx));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn box_geom_populated() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        assert_eq!(brep.geom.edge_curve.len(), 12);
        assert!(brep.geom.edge_curve.iter().all(|c| c.is_some()));
        assert_eq!(brep.geom.face_surface.len(), 6);
        assert!(brep.geom.face_surface.iter().all(|s| s.is_some()));

        // All curves should be lines
        for c in &brep.geom.curves {
            assert!(matches!(c, Curve3::Line(_)));
        }
        // All surfaces should be planes
        for s in &brep.geom.surfaces {
            assert!(matches!(s, Surface3::Plane(_)));
        }
    }
}
