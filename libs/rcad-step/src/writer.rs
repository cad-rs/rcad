use rcad_kernel::{BRep, Curve3, Face, Surface3};
use std::collections::{BTreeSet, HashMap};

pub struct ExportSelection<'a> {
    pub selected_faces: &'a [usize],
    pub selected_edges: &'a [usize],
}

pub struct StepWriter;

struct FaceExportResult {
    face_ids: Vec<u64>,
    used_triangle_fallback: bool,
}

impl StepWriter {
    pub fn write_string(brep: &BRep, selection: ExportSelection<'_>) -> String {
        let mut writer = Part21Writer::new();
        writer.write_brep(brep, selection);
        writer.finish()
    }
}

struct Part21Writer {
    next_id: u64,
    records: Vec<String>,
    vertex_point_ids: HashMap<usize, u64>,
    edge_curve_ids: HashMap<usize, u64>,
}

impl Part21Writer {
    fn new() -> Self {
        Self {
            next_id: 1,
            records: Vec::new(),
            vertex_point_ids: HashMap::new(),
            edge_curve_ids: HashMap::new(),
        }
    }

    fn finish(self) -> String {
        let mut out = String::new();
        out.push_str("ISO-10303-21;\n");
        out.push_str("HEADER;\n");
        out.push_str("FILE_DESCRIPTION(('RCAD exported geometry'),'2;1');\n");
        out.push_str("FILE_NAME('rcad_export.step','2026-04-02T00:00:00',(''),(''),'RCAD','RCAD','');\n");
        out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
        out.push_str("ENDSEC;\n");
        out.push_str("DATA;\n");
        for record in self.records {
            out.push_str(&record);
            out.push('\n');
        }
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        out
    }

    fn write_brep(&mut self, brep: &BRep, selection: ExportSelection<'_>) {
        let selected_face_set: BTreeSet<usize> = selection.selected_faces.iter().copied().collect();
        let selected_edge_set: BTreeSet<usize> = selection.selected_edges.iter().copied().collect();
        let export_all = selected_face_set.is_empty() && selected_edge_set.is_empty();

        let mut face_items = Vec::new();
        let mut solid_items = Vec::new();
        let mut shell_face_groups: Vec<Vec<u64>> = Vec::new();
        let mut has_triangle_fallback = false;

        let mut face_index = 0usize;
        for (solid_index, solid) in brep.solids.iter().enumerate() {
            for (shell_index, shell) in solid.shells.iter().enumerate() {
                let mut shell_faces = Vec::new();
                for face in &shell.faces {
                    if export_all || selected_face_set.contains(&face_index) {
                        let face_surface = brep
                            .geom
                            .face_surface
                            .get(face_index)
                            .and_then(|v| *v)
                            .and_then(|sid| brep.geom.surfaces.get(sid))
                            .copied();
                        let export = self.write_face(brep, face, face_surface);
                        if export.used_triangle_fallback {
                            has_triangle_fallback = true;
                        }
                        face_items.extend(export.face_ids.iter().copied());
                        shell_faces.extend(export.face_ids);
                    }
                    face_index += 1;
                }
                if export_all && !shell_faces.is_empty() {
                    let shell_id = self.closed_shell(
                        &format!("closed_shell_{}_{}", solid_index, shell_index),
                        &shell_faces,
                    );
                    let solid_id = self.manifold_solid_brep(
                        &format!("solid_{}_{}", solid_index, shell_index),
                        shell_id,
                    );
                    solid_items.push(solid_id);
                }
                if !shell_faces.is_empty() {
                    shell_face_groups.push(shell_faces);
                }
            }
        }

        if has_triangle_fallback {
            // Triangulated fallback faces may not form a topologically valid manifold solid.
            // Export as open shell representation to maximize interoperability.
            solid_items.clear();
        }

        // Collect edge indices that belong to face boundaries — these are
        // already part of the solid/shell representation and must NOT be
        // duplicated into the wireframe.
        let mut face_edge_set: BTreeSet<usize> = BTreeSet::new();
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    face_edge_set.extend(face.outer_wire.edges.iter().copied());
                    for inner in &face.inner_wires {
                        face_edge_set.extend(inner.edges.iter().copied());
                    }
                }
            }
        }

        // Only export standalone edges (1D geometry not belonging to any face)
        // into the wireframe. When the user explicitly selected edges, include
        // those regardless.
        let mut edge_items = Vec::new();
        for (edge_index, _edge) in brep.edges.iter().enumerate() {
            if export_all && face_edge_set.contains(&edge_index) {
                continue;
            }
            if export_all || selected_edge_set.contains(&edge_index) {
                edge_items.push(self.write_edge_curve_by_index(brep, edge_index));
            }
        }

        let app_context = self.application_context("automotive_design");
        let _protocol = self.application_protocol_definition(
            "international standard",
            "automotive_design",
            2000,
            app_context,
        );
        let product_context = self.product_context("part definition", app_context, "mechanical");
        let product = self.product("rcad_export", "rcad_export", "", &[product_context]);
        let formation = self.product_definition_formation("", "", product);
        let definition_context =
            self.product_definition_context("part definition", app_context, "design");
        let definition = self.product_definition("", "", formation, definition_context);
        let product_shape = self.product_definition_shape("", "", definition);

        let length_unit = self.length_unit_meter();
        let angle_unit = self.plane_angle_unit_degree();
        let solid_angle_unit = self.solid_angle_unit_steradian();
        let uncertainty = self.uncertainty_measure_with_unit(length_unit);
        let context = self.geometric_representation_context(
            3,
            uncertainty,
            &[length_unit, angle_unit, solid_angle_unit],
            "Context #1",
            "3D Context with UNIT and UNCERTAINTY",
        );

        let base_rep = self.shape_representation("rcad_export", &[], context);
        let _shape_def = self.shape_definition_representation(product_shape, base_rep);

        let mut primary_rep = None;
        if export_all && !solid_items.is_empty() {
            let brep_rep = self.advanced_brep_shape_representation("rcad_export", &solid_items, context);
            self.shape_representation_relationship("", "", base_rep, brep_rep);
            primary_rep = Some(brep_rep);
        } else if !face_items.is_empty() {
            let mut shell_model_items = Vec::new();
            for (i, shell_faces) in shell_face_groups.iter().enumerate() {
                if shell_faces.is_empty() {
                    continue;
                }
                let shell_id = self.open_shell(&format!("export_shell_{}", i), shell_faces);
                let model_id = self.shell_based_surface_model(
                    &format!("export_shell_model_{}", i),
                    &[shell_id],
                );
                shell_model_items.push(model_id);
            }
            if !shell_model_items.is_empty() {
                let surface_rep =
                    self.manifold_surface_shape_representation("rcad_export", &shell_model_items, context);
                self.shape_representation_relationship("", "", base_rep, surface_rep);
                primary_rep = Some(surface_rep);
            }
        }

        if !edge_items.is_empty() {
            let curve_set = self.geometric_curve_set("wireframe", &edge_items);
            let wire_rep = self.geometrically_bounded_wireframe_shape_representation(
                "rcad_export",
                &[curve_set],
                context,
            );
            self.shape_representation_relationship("", "", base_rep, wire_rep);
            if let Some(surface_rep) = primary_rep {
                self.shape_representation_relationship("", "", surface_rep, wire_rep);
            }
        }
    }

    fn write_face(&mut self, brep: &BRep, face: &Face, face_surface: Option<Surface3>) -> FaceExportResult {
        // Only fall back to triangles when there is no analytic surface AND the
        // wire is degenerate (cannot form a valid edge loop).
        if face_surface.is_none() && is_degenerate_face_wire(brep, face) {
            return FaceExportResult {
                face_ids: self.write_triangle_faces(brep, face),
                used_triangle_fallback: true,
            };
        }

        let oriented_edges = oriented_face_edges(brep, face);
        if oriented_edges.is_empty() && face_surface.is_some() {
            // Seam face with no usable edge loop — fall back to triangles.
            return FaceExportResult {
                face_ids: self.write_triangle_faces(brep, face),
                used_triangle_fallback: true,
            };
        }

        let loop_points: Vec<glam::DVec3> = oriented_edges
            .iter()
            .filter_map(|edge| brep.vertices.get(edge.start).map(|v| v.point))
            .collect();
        let origin_point = loop_points.first().copied().unwrap_or(glam::DVec3::ZERO);

        // For seam faces on closed surfaces, the loop points may be collinear,
        // so compute_face_normal can fail. Use the surface's own axis instead.
        let normal = compute_face_normal(&loop_points)
            .or_else(|| surface_normal(face_surface))
            .map(dvec3_to_array)
            .unwrap_or([0.0, 0.0, 1.0]);

        let origin = self.cartesian_point("face_origin", dvec3_to_array(origin_point));
        let axis = self.direction("face_normal", normal);
        let ref_dir = self.direction("face_ref", orthogonal_dir(normal));
        let fallback_placement = self.axis2_placement_3d("face_axis", origin, axis, ref_dir);
        let surface = self.write_surface(face_surface, fallback_placement);

        // Detect seam edges: same edge_idx appearing multiple times
        let seam_edge_indices = detect_seam_edge_indices(face);

        let mut oriented_ids = Vec::new();
        for edge in &oriented_edges {
            let edge_curve = if seam_edge_indices.contains(&edge.edge_idx) {
                // Seam edge: write with a reconstructed curve lying on the surface.
                // Don't cache — the same topological edge gets two distinct STEP
                // EDGE_CURVE entities (one per orientation) so OCCT can build the
                // seam correctly.
                self.write_seam_edge_curve(brep, edge.edge_idx, face_surface)
            } else {
                self.write_edge_curve_by_index(brep, edge.edge_idx)
            };
            oriented_ids.push(self.oriented_edge("face_edge", edge_curve, edge.forward));
        }

        let edge_loop = self.edge_loop("outer_loop", &oriented_ids);
        let face_bound = self.face_outer_bound("outer_bound", edge_loop, true);
        FaceExportResult {
            face_ids: vec![self.advanced_face("face", &[face_bound], surface, true)],
            used_triangle_fallback: false,
        }
    }

    fn write_triangle_faces(&mut self, brep: &BRep, face: &Face) -> Vec<u64> {
        let mut faces = Vec::new();
        for tri in &face.triangles {
            let Some(a) = brep.vertices.get(tri[0]).map(|v| v.point) else {
                continue;
            };
            let Some(b) = brep.vertices.get(tri[1]).map(|v| v.point) else {
                continue;
            };
            let Some(c) = brep.vertices.get(tri[2]).map(|v| v.point) else {
                continue;
            };

            let n = (b - a).cross(c - a).normalize_or_zero();
            if n.length_squared() < 1e-12 {
                continue;
            }

            let origin = self.cartesian_point("tri_origin", dvec3_to_array(a));
            let axis = self.direction("tri_normal", dvec3_to_array(n));
            let ref_dir = self.direction("tri_ref", orthogonal_dir(dvec3_to_array(n)));
            let placement = self.axis2_placement_3d("tri_axis", origin, axis, ref_dir);
            let plane = self.plane("tri_plane", placement);

            let e0 = self.write_edge_curve_from_points(a, b);
            let e1 = self.write_edge_curve_from_points(b, c);
            let e2 = self.write_edge_curve_from_points(c, a);
            let o0 = self.oriented_edge("tri_edge", e0, true);
            let o1 = self.oriented_edge("tri_edge", e1, true);
            let o2 = self.oriented_edge("tri_edge", e2, true);
            let loop_id = self.edge_loop("tri_loop", &[o0, o1, o2]);
            let bound = self.face_outer_bound("tri_bound", loop_id, true);
            faces.push(self.advanced_face("tri_face", &[bound], plane, true));
        }
        faces
    }

    fn write_edge_curve_from_points(&mut self, a: glam::DVec3, b: glam::DVec3) -> u64 {
        let p0 = self.cartesian_point("tri_p0", dvec3_to_array(a));
        let p1 = self.cartesian_point("tri_p1", dvec3_to_array(b));
        let v0 = self.vertex_point("tri_v0", p0);
        let v1 = self.vertex_point("tri_v1", p1);
        let delta = dvec3_to_array(b - a);
        let dir = self.direction("tri_dir", normalize(delta));
        let vec = self.vector("tri_vec", dir, vector_length(delta).max(1e-9));
        let line = self.line("tri_line", p0, vec);
        self.edge_curve("tri_edge", v0, v1, line, true)
    }

    /// Write an EDGE_CURVE for a seam edge, synthesizing a proper 3D curve
    /// that lies on the analytic surface.  This is needed because our BRep
    /// may have lost the original curve (e.g. B-spline) during import.
    ///
    /// OCCT / FreeCAD refuse to import a face whose edge curve does not lie
    /// on the face surface, so we reconstruct a geometrically valid curve:
    ///   - Sphere: the seam is a great circle (meridian)
    ///   - Cylinder / Cone: the seam is a line along the slant/axis
    fn write_seam_edge_curve(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        face_surface: Option<Surface3>,
    ) -> u64 {
        let Some(edge) = brep.edges.get(edge_idx) else {
            return self.write_edge_curve_by_index(brep, edge_idx);
        };
        let start_pt = brep.vertices.get(edge.start).map(|v| v.point).unwrap_or(glam::DVec3::ZERO);
        let end_pt = brep.vertices.get(edge.end).map(|v| v.point).unwrap_or(glam::DVec3::ZERO);

        let v0 = self.vertex_point_by_index(brep, edge.start);
        let v1 = self.vertex_point_by_index(brep, edge.end);

        let basis_curve = match face_surface {
            Some(Surface3::Sphere(sphere)) => {
                // The seam of a sphere is a great circle (meridian).
                // Its centre is the sphere centre, radius = sphere.radius,
                // and its normal is perpendicular to the plane containing
                // the two endpoints and the sphere centre.
                let a = (start_pt - sphere.center).normalize_or_zero();
                let b = (end_pt - sphere.center).normalize_or_zero();
                let mut circle_normal = a.cross(b);
                if circle_normal.length_squared() < 1e-12 {
                    // start and end are antipodal — pick a perpendicular to the axis
                    circle_normal = any_perpendicular_dvec3(sphere.axis);
                }
                let circle_normal = circle_normal.normalize_or_zero();
                let placement = self.axis2_from_origin_axis("seam_axis", sphere.center, circle_normal);
                self.circle("seam_circle", placement, sphere.radius.max(1e-9))
            }
            Some(Surface3::Cone(_cone)) => {
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            Some(Surface3::Cylinder(_)) => {
                // Cylinder seam is a line along the axis.
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            _ => {
                // Fallback: straight line.
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
        };

        self.edge_curve("seam_edge", v0, v1, basis_curve, true)
    }

    fn write_surface(&mut self, face_surface: Option<Surface3>, fallback_placement: u64) -> u64 {
        match face_surface {
            Some(Surface3::Plane(plane)) => {
                let placement = self.axis2_from_origin_axis("plane_axis", plane.origin, plane.normal);
                self.plane("face_plane", placement)
            }
            Some(Surface3::Cylinder(cyl)) => {
                let placement = self.axis2_from_origin_axis("cyl_axis", cyl.origin, cyl.axis);
                self.cylindrical_surface("face_cylinder", placement, cyl.radius.max(1e-9))
            }
            Some(Surface3::Sphere(sphere)) => {
                let placement = self.axis2_from_origin_axis(
                    "sphere_axis",
                    sphere.center,
                    sphere.axis,
                );
                self.spherical_surface("face_sphere", placement, sphere.radius.max(1e-9))
            }
            Some(Surface3::Cone(cone)) => {
                let placement = self.axis2_from_origin_axis("cone_axis", cone.apex, cone.axis);
                self.conical_surface(
                    "face_cone",
                    placement,
                    cone.radius,
                    cone.half_angle_rad.to_degrees(),
                )
            }
            Some(Surface3::Torus(torus)) => {
                let placement = self.axis2_from_origin_axis("torus_axis", torus.center, torus.axis);
                self.toroidal_surface(
                    "face_torus",
                    placement,
                    torus.major_radius.max(1e-9),
                    torus.minor_radius.max(1e-9),
                )
            }
            None => self.plane("face_plane", fallback_placement),
        }
    }

    fn axis2_from_origin_axis(
        &mut self,
        name: &str,
        origin: glam::DVec3,
        axis: glam::DVec3,
    ) -> u64 {
        let origin_id = self.cartesian_point("surface_origin", dvec3_to_array(origin));
        let axis_arr = normalize(dvec3_to_array(axis));
        let axis_id = self.direction("surface_axis", axis_arr);
        let ref_id = self.direction("surface_ref", orthogonal_dir(axis_arr));
        self.axis2_placement_3d(name, origin_id, axis_id, ref_id)
    }

    fn axis2_from_origin_axis_ref(
        &mut self,
        name: &str,
        origin: glam::DVec3,
        axis: glam::DVec3,
        ref_dir: glam::DVec3,
    ) -> u64 {
        let origin_id = self.cartesian_point("curve_origin", dvec3_to_array(origin));
        let axis_arr = normalize(dvec3_to_array(axis));
        let axis_id = self.direction("curve_axis", axis_arr);
        let ref_arr = normalize(project_to_plane(dvec3_to_array(ref_dir), axis_arr));
        let ref_id = self.direction("curve_ref", ref_arr);
        self.axis2_placement_3d(name, origin_id, axis_id, ref_id)
    }

    fn write_edge_curve_by_index(&mut self, brep: &BRep, edge_idx: usize) -> u64 {
        if let Some(existing) = self.edge_curve_ids.get(&edge_idx) {
            return *existing;
        }

        let Some(edge) = brep.edges.get(edge_idx) else {
            let p0 = self.cartesian_point("edge_p0", [0.0, 0.0, 0.0]);
            let p1 = self.cartesian_point("edge_p1", [0.0, 0.0, 0.0]);
            let v0 = self.vertex_point("v0", p0);
            let v1 = self.vertex_point("v1", p1);
            let origin = self.cartesian_point("edge_origin", [0.0, 0.0, 0.0]);
            let dir = self.direction("edge_dir", [1.0, 0.0, 0.0]);
            let vec = self.vector("edge_vec", dir, 1.0);
            let basis = self.line("edge_line", origin, vec);
            return self.edge_curve("edge", v0, v1, basis, true);
        };

        let start_point = brep
            .vertices
            .get(edge.start)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let end_point = brep
            .vertices
            .get(edge.end)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let v0 = self.vertex_point_by_index(brep, edge.start);
        let v1 = self.vertex_point_by_index(brep, edge.end);
        let basis_curve = self.write_basis_curve_for_edge(brep, edge_idx, start_point, end_point);
        let edge_curve = self.edge_curve("edge", v0, v1, basis_curve, true);
        self.edge_curve_ids.insert(edge_idx, edge_curve);
        edge_curve
    }

    fn write_basis_curve_for_edge(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        start_point: [f64; 3],
        end_point: [f64; 3],
    ) -> u64 {
        let curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .and_then(|v| *v)
            .and_then(|curve_idx| brep.geom.curves.get(curve_idx))
            .copied();

        match curve {
            Some(Curve3::Line(line)) => {
                let origin = self.cartesian_point("line_origin", dvec3_to_array(line.origin));
                let dir = normalize(dvec3_to_array(line.direction));
                let dir_id = self.direction("line_dir", dir);
                let len = vector_length([
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ])
                .max(1e-9);
                let vec_id = self.vector("line_vec", dir_id, len);
                self.line("edge_line", origin, vec_id)
            }
            Some(Curve3::Circle(circle)) => {
                let placement = self.axis2_from_origin_axis("circle_axis", circle.center, circle.normal);
                self.circle("edge_circle", placement, circle.radius.max(1e-9))
            }
            Some(Curve3::Ellipse(ellipse)) => {
                let placement = self.axis2_from_origin_axis_ref(
                    "ellipse_axis",
                    ellipse.center,
                    ellipse.normal,
                    ellipse.major_dir,
                );
                self.ellipse(
                    "edge_ellipse",
                    placement,
                    ellipse.major_radius.max(1e-9),
                    ellipse.minor_radius.max(1e-9),
                )
            }
            None => {
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
        }
    }

    fn vertex_point_by_index(&mut self, brep: &BRep, vertex_idx: usize) -> u64 {
        if let Some(existing) = self.vertex_point_ids.get(&vertex_idx) {
            return *existing;
        }

        let point = brep
            .vertices
            .get(vertex_idx)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let cartesian = self.cartesian_point("vertex_point", point);
        let vertex = self.vertex_point("vertex", cartesian);
        self.vertex_point_ids.insert(vertex_idx, vertex);
        vertex
    }

    fn application_context(&mut self, name: &str) -> u64 {
        self.push(format!("APPLICATION_CONTEXT('{}')", name))
    }

    fn application_protocol_definition(
        &mut self,
        status: &str,
        schema: &str,
        year: i32,
        context: u64,
    ) -> u64 {
        self.push(format!(
            "APPLICATION_PROTOCOL_DEFINITION('{}','{}',{},#{})",
            status, schema, year, context
        ))
    }

    fn product_context(&mut self, name: &str, frame: u64, discipline: &str) -> u64 {
        self.push(format!(
            "PRODUCT_CONTEXT('{}',#{},'{}')",
            name, frame, discipline
        ))
    }

    fn product(&mut self, id: &str, name: &str, description: &str, contexts: &[u64]) -> u64 {
        self.push(format!(
            "PRODUCT('{}','{}','{}',({}))",
            id,
            name,
            description,
            refs(contexts)
        ))
    }

    fn product_definition_formation(&mut self, id: &str, description: &str, product: u64) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_FORMATION('{}','{}',#{})",
            id, description, product
        ))
    }

    fn product_definition_context(&mut self, name: &str, frame: u64, life_cycle: &str) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_CONTEXT('{}',#{},'{}')",
            name, frame, life_cycle
        ))
    }

    fn product_definition(
        &mut self,
        id: &str,
        description: &str,
        formation: u64,
        frame: u64,
    ) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION('{}','{}',#{},#{})",
            id, description, formation, frame
        ))
    }

    fn product_definition_shape(&mut self, name: &str, description: &str, definition: u64) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_SHAPE('{}','{}',#{})",
            name, description, definition
        ))
    }

    fn shape_definition_representation(&mut self, shape: u64, representation: u64) -> u64 {
        self.push(format!(
            "SHAPE_DEFINITION_REPRESENTATION(#{},#{})",
            shape, representation
        ))
    }

    fn length_unit_meter(&mut self) -> u64 {
        self.push("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string())
    }

    fn plane_angle_unit_degree(&mut self) -> u64 {
        let radian_unit =
            self.push("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string());
        let measure = self.push(format!(
            "PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#{})",
            radian_unit
        ));
        let dim_exp =
            self.push("DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)".to_string());
        self.push(format!(
            "( CONVERSION_BASED_UNIT('DEGREE',#{}) NAMED_UNIT(#{}) PLANE_ANGLE_UNIT() )",
            measure, dim_exp
        ))
    }

    fn solid_angle_unit_steradian(&mut self) -> u64 {
        self.push("( NAMED_UNIT(*) SOLID_ANGLE_UNIT() SI_UNIT($,.STERADIAN.) )".to_string())
    }

    fn uncertainty_measure_with_unit(&mut self, length_unit: u64) -> u64 {
        self.push(format!(
            "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-6),#{},'distance_accuracy_value','confusion accuracy')",
            length_unit
        ))
    }

    fn geometric_representation_context(
        &mut self,
        dimension: i32,
        uncertainty: u64,
        units: &[u64],
        name: &str,
        description: &str,
    ) -> u64 {
        self.push(format!(
            "( GEOMETRIC_REPRESENTATION_CONTEXT({}) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT(({})) REPRESENTATION_CONTEXT('{}','{}') )",
            dimension,
            uncertainty,
            refs(units),
            name,
            description
        ))
    }

    fn cartesian_point(&mut self, name: &str, coords: [f64; 3]) -> u64 {
        self.push(format!(
            "CARTESIAN_POINT('{}',({:.9},{:.9},{:.9}))",
            name, coords[0], coords[1], coords[2]
        ))
    }

    fn direction(&mut self, name: &str, coords: [f64; 3]) -> u64 {
        self.push(format!(
            "DIRECTION('{}',({:.9},{:.9},{:.9}))",
            name, coords[0], coords[1], coords[2]
        ))
    }

    fn vector(&mut self, name: &str, direction: u64, magnitude: f64) -> u64 {
        self.push(format!("VECTOR('{}',#{},{:.9})", name, direction, magnitude))
    }

    fn axis2_placement_3d(&mut self, name: &str, origin: u64, axis: u64, ref_dir: u64) -> u64 {
        self.push(format!(
            "AXIS2_PLACEMENT_3D('{}',#{},#{},#{})",
            name, origin, axis, ref_dir
        ))
    }

    fn line(&mut self, name: &str, origin: u64, vector: u64) -> u64 {
        self.push(format!("LINE('{}',#{},#{})", name, origin, vector))
    }

    fn circle(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!("CIRCLE('{}',#{},{:.9})", name, placement, radius))
    }

    fn ellipse(&mut self, name: &str, placement: u64, major: f64, minor: f64) -> u64 {
        self.push(format!(
            "ELLIPSE('{}',#{},{:.9},{:.9})",
            name, placement, major, minor
        ))
    }

    fn plane(&mut self, name: &str, placement: u64) -> u64 {
        self.push(format!("PLANE('{}',#{})", name, placement))
    }

    fn cylindrical_surface(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!(
            "CYLINDRICAL_SURFACE('{}',#{},{:.9})",
            name, placement, radius
        ))
    }

    fn spherical_surface(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!(
            "SPHERICAL_SURFACE('{}',#{},{:.9})",
            name, placement, radius
        ))
    }

    fn conical_surface(&mut self, name: &str, placement: u64, radius: f64, semi_angle_deg: f64) -> u64 {
        self.push(format!(
            "CONICAL_SURFACE('{}',#{},{:.9},{:.9})",
            name, placement, radius, semi_angle_deg
        ))
    }

    fn toroidal_surface(
        &mut self,
        name: &str,
        placement: u64,
        major_radius: f64,
        minor_radius: f64,
    ) -> u64 {
        self.push(format!(
            "TOROIDAL_SURFACE('{}',#{},{:.9},{:.9})",
            name, placement, major_radius, minor_radius
        ))
    }

    fn vertex_point(&mut self, name: &str, point: u64) -> u64 {
        self.push(format!("VERTEX_POINT('{}',#{})", name, point))
    }

    fn edge_curve(&mut self, name: &str, start: u64, end: u64, curve: u64, same_sense: bool) -> u64 {
        self.push(format!(
            "EDGE_CURVE('{}',#{},#{},#{},{})",
            name,
            start,
            end,
            curve,
            bool_token(same_sense)
        ))
    }

    fn oriented_edge(&mut self, name: &str, edge_curve: u64, orientation: bool) -> u64 {
        self.push(format!(
            "ORIENTED_EDGE('{}',*,*,#{},{})",
            name,
            edge_curve,
            bool_token(orientation)
        ))
    }

    fn edge_loop(&mut self, name: &str, oriented_edges: &[u64]) -> u64 {
        self.push(format!("EDGE_LOOP('{}',({}))", name, refs(oriented_edges)))
    }

    fn face_outer_bound(&mut self, name: &str, edge_loop: u64, orientation: bool) -> u64 {
        self.push(format!(
            "FACE_OUTER_BOUND('{}',#{},{})",
            name,
            edge_loop,
            bool_token(orientation)
        ))
    }

    fn advanced_face(&mut self, name: &str, bounds: &[u64], surface: u64, orientation: bool) -> u64 {
        self.push(format!(
            "ADVANCED_FACE('{}',({}),#{},{})",
            name,
            refs(bounds),
            surface,
            bool_token(orientation)
        ))
    }

    fn open_shell(&mut self, name: &str, faces: &[u64]) -> u64 {
        self.push(format!("OPEN_SHELL('{}',({}))", name, refs(faces)))
    }

    fn closed_shell(&mut self, name: &str, faces: &[u64]) -> u64 {
        self.push(format!("CLOSED_SHELL('{}',({}))", name, refs(faces)))
    }

    fn shell_based_surface_model(&mut self, name: &str, shells: &[u64]) -> u64 {
        self.push(format!("SHELL_BASED_SURFACE_MODEL('{}',({}))", name, refs(shells)))
    }

    fn manifold_solid_brep(&mut self, name: &str, outer: u64) -> u64 {
        self.push(format!("MANIFOLD_SOLID_BREP('{}',#{})", name, outer))
    }

    fn geometric_curve_set(&mut self, name: &str, curves: &[u64]) -> u64 {
        self.push(format!("GEOMETRIC_CURVE_SET('{}',({}))", name, refs(curves)))
    }

    fn shape_representation(&mut self, name: &str, items: &[u64], context: u64) -> u64 {
        self.push(format!(
            "SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn manifold_surface_shape_representation(&mut self, name: &str, items: &[u64], context: u64) -> u64 {
        self.push(format!(
            "MANIFOLD_SURFACE_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn advanced_brep_shape_representation(&mut self, name: &str, items: &[u64], context: u64) -> u64 {
        self.push(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn geometrically_bounded_wireframe_shape_representation(
        &mut self,
        name: &str,
        items: &[u64],
        context: u64,
    ) -> u64 {
        self.push(format!(
            "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn shape_representation_relationship(
        &mut self,
        name: &str,
        description: &str,
        rep_1: u64,
        rep_2: u64,
    ) -> u64 {
        self.push(format!(
            "SHAPE_REPRESENTATION_RELATIONSHIP('{}','{}',#{},#{})",
            name,
            description,
            rep_1,
            rep_2
        ))
    }

    fn push(&mut self, body: String) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        self.records.push(format!("#{}={};", id, body));
        id
    }
}

#[derive(Clone, Copy)]
struct OrientedEdgeExport {
    edge_idx: usize,
    start: usize,
    end: usize,
    forward: bool,
}

fn oriented_face_edges(brep: &BRep, face: &Face) -> Vec<OrientedEdgeExport> {
    let mut pending: Vec<(usize, usize, usize)> = face
        .outer_wire
        .edges
        .iter()
        .filter_map(|&edge_idx| brep.edges.get(edge_idx).map(|edge| (edge_idx, edge.start, edge.end)))
        .collect();
    if pending.is_empty() {
        return Vec::new();
    }

    let (first_idx, first_start, first_end) = pending.remove(0);
    let mut result = vec![OrientedEdgeExport {
        edge_idx: first_idx,
        start: first_start,
        end: first_end,
        forward: true,
    }];
    let mut current = first_end;

    while !pending.is_empty() {
        if let Some((idx, oriented)) = pending.iter().enumerate().find_map(|(idx, &(edge_idx, a, b))| {
            if a == current {
                Some((idx, OrientedEdgeExport { edge_idx, start: a, end: b, forward: true }))
            } else if b == current {
                Some((idx, OrientedEdgeExport { edge_idx, start: b, end: a, forward: false }))
            } else {
                None
            }
        }) {
            pending.remove(idx);
            current = oriented.end;
            result.push(oriented);
        } else {
            let (edge_idx, a, b) = pending.remove(0);
            current = b;
            result.push(OrientedEdgeExport { edge_idx, start: a, end: b, forward: true });
        }
    }

    result
}

fn compute_face_normal(points: &[glam::DVec3]) -> Option<glam::DVec3> {
    if points.len() < 3 {
        return None;
    }
    let origin = points[0];
    for i in 1..points.len().saturating_sub(1) {
        let a = points[i] - origin;
        let b = points[i + 1] - origin;
        let n = a.cross(b);
        if n.length_squared() > 1e-12 {
            return Some(n.normalize());
        }
    }
    None
}

fn refs(items: &[u64]) -> String {
    items
        .iter()
        .map(|id| format!("#{}", id))
        .collect::<Vec<_>>()
        .join(",")
}

fn bool_token(value: bool) -> &'static str {
    if value {
        ".T."
    } else {
        ".F."
    }
}

fn dvec3_to_array(v: glam::DVec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn vector_length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = vector_length(v);
    if len <= 1e-12 {
        [1.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn project_to_plane(v: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let dot = v[0] * normal[0] + v[1] * normal[1] + v[2] * normal[2];
    [
        v[0] - normal[0] * dot,
        v[1] - normal[1] * dot,
        v[2] - normal[2] * dot,
    ]
}

fn orthogonal_dir(normal: [f64; 3]) -> [f64; 3] {
    let helper = if normal[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    normalize(cross(normal, helper))
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn any_perpendicular_dvec3(v: glam::DVec3) -> glam::DVec3 {
    let helper = if v.dot(glam::DVec3::Y).abs() < 0.9 {
        glam::DVec3::Y
    } else {
        glam::DVec3::X
    };
    v.cross(helper).normalize_or_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StepReader;
    use glam::DVec3;
    use rcad_modeling::make_box_brep;
    const HFSS_STEP: &str = include_str!("../../../assets/hfss.step");

    #[test]
    fn exports_full_box_and_reimports() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("exported STEP should parse");
        assert!(!reparsed.edges.is_empty());
        assert!(!reparsed.solids.is_empty());
        assert!(step.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"));
        assert!(step.contains("MANIFOLD_SOLID_BREP"));
        assert!(step.contains("CLOSED_SHELL"));
    }

    #[test]
    fn exports_selected_edges_without_faces() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[0, 1],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("edge-only export should parse");
        assert!(reparsed.solids.is_empty());
        assert_eq!(reparsed.edges.len(), 2);
    }

    #[test]
    fn exports_selected_faces_via_shell_based_surface_model() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[0],
                selected_edges: &[],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("selected-face export should parse");
        assert!(!reparsed.solids.is_empty());
        assert!(step.contains("OPEN_SHELL"));
        assert!(step.contains("SHELL_BASED_SURFACE_MODEL"));
        assert!(step.contains("MANIFOLD_SURFACE_SHAPE_REPRESENTATION"));
    }

    #[test]
    fn exports_analytic_surfaces_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        // All analytic surfaces should now be exported properly as ADVANCED_FACE
        // with their respective surface types, including seam faces on
        // spheres and cones.
        assert!(step.contains("SPHERICAL_SURFACE"));
        assert!(step.contains("CYLINDRICAL_SURFACE"));
        assert!(step.contains("TOROIDAL_SURFACE"));
        assert!(step.contains("CONICAL_SURFACE"));

        // Standalone 1D curves (GEOMETRIC_CURVE_SET) must also be exported
        // alongside the solid geometry.
        assert!(
            step.contains("GEOMETRIC_CURVE_SET"),
            "standalone wireframe edges should be exported"
        );
        assert!(
            step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"),
            "wireframe shape representation should be present"
        );
    }

    #[test]
    fn round_trips_sphere_and_cone_surfaces() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        // Find the original cone half-angle and radius for comparison
        let mut orig_cone_angle = 0.0f64;
        let mut orig_cone_radius = 0.0f64;
        for surface in &brep.geom.surfaces {
            if let Surface3::Cone(c) = surface {
                orig_cone_angle = c.half_angle_rad;
                orig_cone_radius = c.radius;
            }
        }
        assert!(orig_cone_angle > 0.0, "should find a cone in hfss.step");

        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("re-exported STEP should parse");

        // Count faces with each surface type and verify cone parameters survive round-trip
        let mut sphere_count = 0usize;
        let mut cone_count = 0usize;
        for surface_binding in &reparsed.geom.face_surface {
            if let Some(sid) = surface_binding {
                match reparsed.geom.surfaces.get(*sid) {
                    Some(Surface3::Sphere(_)) => sphere_count += 1,
                    Some(Surface3::Cone(c)) => {
                        cone_count += 1;
                        assert!(
                            (c.half_angle_rad - orig_cone_angle).abs() < 1e-6,
                            "cone half-angle drifted: original={} reparsed={}",
                            orig_cone_angle,
                            c.half_angle_rad,
                        );
                        assert!(
                            (c.radius - orig_cone_radius).abs() < 1e-6,
                            "cone radius drifted: original={} reparsed={}",
                            orig_cone_radius,
                            c.radius,
                        );
                    }
                    _ => {}
                }
            }
        }
        assert!(sphere_count >= 1, "expected at least 1 sphere face after round-trip, got {}", sphere_count);
        assert!(cone_count >= 1, "expected at least 1 cone face after round-trip, got {}", cone_count);
    }
}

/// Detect which edge indices appear more than once in the face's outer wire.
/// These are seam edges on periodic surfaces.
fn detect_seam_edge_indices(face: &Face) -> BTreeSet<usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for &idx in &face.outer_wire.edges {
        *counts.entry(idx).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .filter(|&(_, c)| c >= 2)
        .map(|(idx, _)| idx)
        .collect()
}

/// Extract a representative normal/axis from an analytic surface, used as a
/// fallback when boundary loop points are collinear (e.g. seam faces).
fn surface_normal(face_surface: Option<Surface3>) -> Option<glam::DVec3> {
    match face_surface? {
        Surface3::Plane(p) => Some(p.normal),
        Surface3::Cylinder(c) => Some(c.axis),
        Surface3::Sphere(s) => Some(s.axis),
        Surface3::Cone(c) => Some(c.axis),
        Surface3::Torus(t) => Some(t.axis),
    }
}

fn is_degenerate_face_wire(brep: &BRep, face: &Face) -> bool {
    if face.outer_wire.edges.len() < 3 {
        return true;
    }

    let unique_edges: BTreeSet<usize> = face.outer_wire.edges.iter().copied().collect();
    if unique_edges.len() < 3 {
        return true;
    }

    let mut verts = BTreeSet::new();
    for &edge_idx in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(edge_idx) {
            verts.insert(edge.start);
            verts.insert(edge.end);
        }
    }
    verts.len() < 3
}