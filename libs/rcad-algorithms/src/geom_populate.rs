use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::*;

/// Populates `brep.geom` with analytic geometry for a box BRep.
///
/// After this call, every edge has a `Curve3::Line` and every face has a `Surface3::Plane`.
/// Precondition: brep was created by `BRep::from_primitive(Box{..})`.
pub fn populate_box_geom(brep: &mut BRep) {
    let geom = &mut brep.geom;
    geom.curves.clear();
    geom.edge_curve.clear();
    geom.edge_curve_range.clear();
    geom.edge_degenerated.clear();
    geom.surfaces.clear();
    geom.face_surface.clear();

    // Edges → Line3
    for edge in &brep.edges {
        let p0 = brep.vertices[edge.start].point;
        let p1 = brep.vertices[edge.end].point;
        let delta = p1 - p0;
        let len = delta.length();
        let dir = if len > 1e-12 { delta / len } else { DVec3::X };
        let curve_idx = geom.curves.len();
        geom.curves.push(Curve3::Line(Line3 {
            origin: p0,
            direction: dir,
        }));
        geom.edge_curve.push(Some(curve_idx));
        // t_range: project endpoints onto the line
        let t0 = 0.0_f64;
        let t1 = (p1 - p0).dot(dir); // = len
        geom.edge_curve_range.push(Some([t0, t1]));
        geom.edge_degenerated.push(len <= 1e-12);
    }

    // Faces → Plane
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Use the first wire vertex rather than face.triangles (triangles
                // are rendering metadata and must not be used in modeling code).
                let origin = face
                    .outer_wire
                    .edges
                    .first()
                    .and_then(|we| brep.edges.get(we.idx))
                    .map(|e| brep.vertices[e.start].point)
                    .unwrap_or(DVec3::ZERO);
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
