//! Python bindings for the open-source RCAD kernel (B-rep, primitives, booleans, STEP/IGES).

use glam::{DAffine3, DQuat, DVec2, DVec3};
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;
use rcad_algorithms::{
    boolean_op, boolean_op_with_options, brep_proj_cylindrical, linear_sweep_wire, offset_solid,
    pipe_sweep_wire, repair, resolved_boolean_fuzzy_tol_for_ds, tolerance, BrepProjOptions,
    BooleanOpType, BooleanOptions,
};
use rcad_kernel::properties::{centroid, inertia_tensor, signed_volume, surface_area, volume};
use rcad_kernel::BRep;
use rcad_modeling::{
    chamfer_edge, cone_brep, cylinder_brep, extrude, fillet_edge, fillet_edges, loft, make_box_brep,
    make_conical_frustum_brep, project_wire_onto_surface, revolve, sphere_brep, sweep_pipe,
    torus_brep,
};
use rcad_step::{
    ExportSelection, IgesReader, IgesWriter, StepReader,
    StepWriter,
};

fn vec3_tuple(v: (f64, f64, f64)) -> DVec3 {
    DVec3::new(v.0, v.1, v.2)
}

fn tuple3(v: DVec3) -> (f64, f64, f64) {
    (v.x, v.y, v.z)
}

fn face_count(brep: &BRep) -> usize {
    brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Flattened face index for ``face_idx`` in the first solid's first shell (matches ``extrude``).
fn face_flat_index_first_shell(brep: &BRep, face_idx: usize) -> PyResult<usize> {
    let mut flat = 0usize;
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            if si == 0 && shi == 0 {
                if face_idx >= shell.faces.len() {
                    return Err(PyValueError::new_err(format!(
                        "face_idx {face_idx} out of range (first shell has {} faces)",
                        shell.faces.len()
                    )));
                }
                return Ok(flat + face_idx);
            }
            flat += shell.faces.len();
        }
    }
    Err(PyValueError::new_err(
        "BRep has no geometry in the first solid/shell",
    ))
}

fn step_metadata_to_py(
    py: Python<'_>,
    meta: &rcad_step::StepDocumentMetadata,
) -> PyResult<Py<PyAny>> {
    let s = serde_json::to_string(meta).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let json_mod = py.import("json")?;
    let obj = json_mod.call_method1("loads", (s,))?;
    Ok(obj.unbind())
}

fn py_bool_err(e: rcad_algorithms::BooleanError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

fn py_build_err(e: rcad_modeling::BuildError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn py_step_err(e: rcad_step::StepError) -> PyErr {
    PyIOError::new_err(e.to_string())
}

fn py_iges_err(e: rcad_step::IgesError) -> PyErr {
    PyIOError::new_err(e.to_string())
}

fn py_offset_err(e: rcad_algorithms::OffsetError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

fn py_sweep_err(e: rcad_algorithms::SweepError) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

/// Boolean pipeline options matching ``rcad_algorithms::BooleanOptions`` (commonly tuned subset).
///
/// Keyword-only constructor, e.g. ``BooleanOptions(fuzzy_tol=1e-5, run_propagate_geom_tolerances=True)``.
#[pyclass(name = "BooleanOptions")]
#[derive(Clone, Copy)]
pub struct PyBooleanOptions {
    inner: BooleanOptions,
}

#[pymethods]
impl PyBooleanOptions {
    #[new]
    #[pyo3(signature = (
        *,
        use_bvh = true,
        run_healing = false,
        run_simplify = false,
        include_history = false,
        run_make_connected = false,
        make_connected_tolerance = None,
        fuzzy_tol = 0.0,
        use_glue = false,
        glue_tolerance = None,
        run_propagate_geom_tolerances = false,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        use_bvh: bool,
        run_healing: bool,
        run_simplify: bool,
        include_history: bool,
        run_make_connected: bool,
        make_connected_tolerance: Option<f64>,
        fuzzy_tol: f64,
        use_glue: bool,
        glue_tolerance: Option<f64>,
        run_propagate_geom_tolerances: bool,
    ) -> Self {
        let mut inner = BooleanOptions::default();
        inner.use_bvh = use_bvh;
        inner.run_healing = run_healing;
        inner.run_simplify = run_simplify;
        inner.include_history = include_history;
        inner.run_make_connected = run_make_connected;
        if let Some(t) = make_connected_tolerance {
            inner.make_connected_tolerance = t;
        }
        inner.fuzzy_tol = fuzzy_tol;
        inner.use_glue = use_glue;
        if let Some(t) = glue_tolerance {
            inner.glue_tolerance = t;
        }
        inner.run_propagate_geom_tolerances = run_propagate_geom_tolerances;
        Self { inner }
    }
}

/// Boundary representation of a solid model (vertices, edges, faces, solids).
#[pyclass(name = "BRep")]
#[derive(Clone)]
pub struct PyBRep {
    inner: BRep,
}

#[pymethods]
impl PyBRep {
    /// Unit box from origin, axis directions, and sizes along x/y/z.
    #[staticmethod]
    #[pyo3(signature = (origin, x_dir, y_dir, width, height, depth))]
    fn box_(
        origin: (f64, f64, f64),
        x_dir: (f64, f64, f64),
        y_dir: (f64, f64, f64),
        width: f64,
        height: f64,
        depth: f64,
    ) -> PyResult<Self> {
        let brep = make_box_brep(
            vec3_tuple(origin),
            vec3_tuple(x_dir),
            vec3_tuple(y_dir),
            width,
            height,
            depth,
        )
        .map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Sphere centered at `center` with `radius`.
    #[staticmethod]
    fn sphere(center: (f64, f64, f64), radius: f64) -> PyResult<Self> {
        let brep = sphere_brep(vec3_tuple(center), radius).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Cylinder: `axis` and `ref_dir` define the local frame (see RCAD modeling docs).
    #[staticmethod]
    #[pyo3(signature = (center, axis, ref_dir, radius, height))]
    fn cylinder(
        center: (f64, f64, f64),
        axis: (f64, f64, f64),
        ref_dir: (f64, f64, f64),
        radius: f64,
        height: f64,
    ) -> PyResult<Self> {
        let brep = cylinder_brep(
            vec3_tuple(center),
            vec3_tuple(axis),
            vec3_tuple(ref_dir),
            radius,
            height,
        )
        .map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Right circular cone.
    #[staticmethod]
    #[pyo3(signature = (center, axis, ref_dir, base_radius, height))]
    fn cone(
        center: (f64, f64, f64),
        axis: (f64, f64, f64),
        ref_dir: (f64, f64, f64),
        base_radius: f64,
        height: f64,
    ) -> PyResult<Self> {
        let brep = cone_brep(
            vec3_tuple(center),
            vec3_tuple(axis),
            vec3_tuple(ref_dir),
            base_radius,
            height,
        )
        .map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Torus with `major_radius` and `minor_radius`.
    #[staticmethod]
    #[pyo3(signature = (center, axis, ref_dir, major_radius, minor_radius))]
    fn torus(
        center: (f64, f64, f64),
        axis: (f64, f64, f64),
        ref_dir: (f64, f64, f64),
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        let brep = torus_brep(
            vec3_tuple(center),
            vec3_tuple(axis),
            vec3_tuple(ref_dir),
            major_radius,
            minor_radius,
        )
        .map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Truncated cone (frustum): `center` is midpoint between end faces; `axis` from bottom toward top.
    #[staticmethod]
    #[pyo3(signature = (center, axis, ref_dir, r_bottom, r_top, height))]
    fn conical_frustum(
        center: (f64, f64, f64),
        axis: (f64, f64, f64),
        ref_dir: (f64, f64, f64),
        r_bottom: f64,
        r_top: f64,
        height: f64,
    ) -> PyResult<Self> {
        let brep = make_conical_frustum_brep(
            vec3_tuple(center),
            vec3_tuple(axis),
            vec3_tuple(ref_dir),
            r_bottom,
            r_top,
            height,
        )
        .map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Loft between parallel cross-sections. Each profile is a list of ``(x, y, z)`` vertices; all profiles must have the same length (≥ 3).
    #[staticmethod]
    fn loft(profiles: Vec<Vec<(f64, f64, f64)>>) -> PyResult<Self> {
        let v: Vec<Vec<DVec3>> = profiles
            .into_iter()
            .map(|p| p.into_iter().map(vec3_tuple).collect())
            .collect();
        let brep = loft(&v).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Pipe sweep: 2D profile in the local XY plane (≥3 points) along a 3D spine polyline (≥2 points). Uses lofted cross-sections (``rcad_modeling::sweep_pipe``).
    #[staticmethod]
    fn sweep_pipe(profile_2d: Vec<(f64, f64)>, spine: Vec<(f64, f64, f64)>) -> PyResult<Self> {
        let p2: Vec<DVec2> = profile_2d
            .iter()
            .map(|(x, y)| DVec2::new(*x, *y))
            .collect();
        let sp: Vec<DVec3> = spine.into_iter().map(vec3_tuple).collect();
        let brep = sweep_pipe(&p2, &sp).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Pipe sweep wire (``rcad_algorithms::pipe_sweep_wire``): alternate sweep kernel with pipe mode; same arguments as ``sweep_pipe``.
    #[staticmethod]
    fn pipe_sweep_wire(profile_2d: Vec<(f64, f64)>, spine: Vec<(f64, f64, f64)>) -> PyResult<Self> {
        let p2: Vec<DVec2> = profile_2d
            .iter()
            .map(|(x, y)| DVec2::new(*x, *y))
            .collect();
        let sp: Vec<DVec3> = spine.into_iter().map(vec3_tuple).collect();
        let brep = pipe_sweep_wire(&p2, &sp).map_err(py_sweep_err)?;
        Ok(Self { inner: brep })
    }

    /// Load a STEP file (``.step`` / ``.stp``).
    #[staticmethod]
    fn read_step(path: &str) -> PyResult<Self> {
        let brep = StepReader::read_file(path).map_err(py_step_err)?;
        Ok(Self { inner: brep })
    }

    /// Load STEP geometry plus document metadata as a ``dict`` (JSON round-trip of ``StepDocumentMetadata``: products, file schema, tolerances, PMI-related lists, etc.).
    #[staticmethod]
    fn read_step_with_metadata(py: Python<'_>, path: &str) -> PyResult<(Self, Py<PyAny>)> {
        let (brep, meta) = StepReader::read_file_with_metadata(path).map_err(py_step_err)?;
        let d = step_metadata_to_py(py, &meta)?;
        Ok((Self { inner: brep }, d))
    }

    /// Write AP214 STEP to a file using OCCT-style interchange (``SURFACE_CURVE`` / ``PCURVE``, si_metre, radians).
    #[pyo3(signature = (path))]
    fn write_step(&self, path: &str) -> PyResult<()> {
        let s = StepWriter::write_string(
            &self.inner,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );
        std::fs::write(path, s).map_err(|e| PyIOError::new_err(e.to_string()))
    }

    /// Load an IGES B-rep file.
    #[staticmethod]
    fn read_iges(path: &str) -> PyResult<Self> {
        let brep = IgesReader::read_file(path).map_err(py_iges_err)?;
        Ok(Self { inner: brep })
    }

    /// Write IGES B-rep to a file.
    fn write_iges(&self, path: &str) -> PyResult<()> {
        IgesWriter::write_file(&self.inner, path).map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(())
    }

    /// Boolean union ``self | other``.
    fn union(&self, other: &PyBRep) -> PyResult<Self> {
        let out = boolean_op(BooleanOpType::Union, &self.inner, &other.inner).map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Boolean union ``self | other`` with execution options (see [`BooleanOptions`]).
    fn union_with_options(&self, other: &PyBRep, options: &PyBooleanOptions) -> PyResult<Self> {
        let (out, _) = boolean_op_with_options(
            BooleanOpType::Union,
            &self.inner,
            &other.inner,
            options.inner,
        )
        .map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Boolean intersection ``self & other``.
    fn intersection(&self, other: &PyBRep) -> PyResult<Self> {
        let out =
            boolean_op(BooleanOpType::Intersection, &self.inner, &other.inner).map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Boolean intersection with execution options (see [`BooleanOptions`]).
    fn intersection_with_options(
        &self,
        other: &PyBRep,
        options: &PyBooleanOptions,
    ) -> PyResult<Self> {
        let (out, _) = boolean_op_with_options(
            BooleanOpType::Intersection,
            &self.inner,
            &other.inner,
            options.inner,
        )
        .map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Boolean difference ``self - other``.
    fn difference(&self, other: &PyBRep) -> PyResult<Self> {
        let out =
            boolean_op(BooleanOpType::Difference, &self.inner, &other.inner).map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Boolean difference ``self - other`` with execution options (see [`BooleanOptions`]).
    fn difference_with_options(
        &self,
        other: &PyBRep,
        options: &PyBooleanOptions,
    ) -> PyResult<Self> {
        let (out, _) = boolean_op_with_options(
            BooleanOpType::Difference,
            &self.inner,
            &other.inner,
            options.inner,
        )
        .map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Extrude face ``face_idx`` (index in the first solid's outer shell) along ``direction`` by ``distance``.
    #[pyo3(signature = (face_idx, direction, distance))]
    fn extrude(&self, face_idx: usize, direction: (f64, f64, f64), distance: f64) -> PyResult<Self> {
        let brep =
            extrude(&self.inner, face_idx, vec3_tuple(direction), distance).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Revolve face ``face_idx`` around an axis; ``angle`` is in **radians**.
    #[pyo3(signature = (face_idx, axis_origin, axis_dir, angle))]
    fn revolve(
        &self,
        face_idx: usize,
        axis_origin: (f64, f64, f64),
        axis_dir: (f64, f64, f64),
        angle: f64,
    ) -> PyResult<Self> {
        let brep = revolve(
            &self.inner,
            face_idx,
            vec3_tuple(axis_origin),
            vec3_tuple(axis_dir),
            angle,
        )
        .map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Fillet a single edge (by global edge index).
    fn fillet_edge(&self, edge_idx: usize, radius: f64) -> PyResult<Self> {
        let brep = fillet_edge(&self.inner, edge_idx, radius).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Fillet multiple edges; ``edges`` is ``[(edge_index, radius), ...]``.
    fn fillet_edges(&self, edges: Vec<(usize, f64)>) -> PyResult<Self> {
        let brep = fillet_edges(&self.inner, &edges).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Chamfer a single edge by distance along adjacent faces.
    fn chamfer_edge(&self, edge_idx: usize, distance: f64) -> PyResult<Self> {
        let brep = chamfer_edge(&self.inner, edge_idx, distance).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Repair mesh / B-rep connectivity (merge vertices, etc.); ``tolerance`` is the merge tolerance.
    fn repair(&self, tolerance: f64) -> PyResult<Self> {
        let (brep, _report) = repair(&self.inner, tolerance);
        Ok(Self { inner: brep })
    }

    /// Offset the first solid's shell by ``distance`` (along face normals; see ``rcad_algorithms::offset_solid``).
    fn offset_solid(&self, distance: f64) -> PyResult<Self> {
        let solid = self
            .inner
            .solids
            .first()
            .ok_or_else(|| PyValueError::new_err("BRep has no solids"))?;
        let brep = offset_solid(solid, &self.inner, distance).map_err(py_offset_err)?;
        Ok(Self { inner: brep })
    }

    /// Linear sweep of the outer wire of face ``face_idx`` (first solid, first shell — same indexing as ``extrude``).
    #[pyo3(signature = (face_idx, direction, distance))]
    fn sweep_wire_linear(
        &self,
        face_idx: usize,
        direction: (f64, f64, f64),
        distance: f64,
    ) -> PyResult<Self> {
        let brep = linear_sweep_wire(
            &self.inner,
            face_idx,
            vec3_tuple(direction),
            distance,
        )
        .map_err(py_sweep_err)?;
        Ok(Self { inner: brep })
    }

    /// Project the outer wire of face ``face_idx`` onto that face's underlying surface; returns projected vertex positions ``[(x,y,z), ...]``.
    #[pyo3(signature = (face_idx, n_samples=32))]
    fn project_wire_to_face(
        &self,
        face_idx: usize,
        n_samples: usize,
    ) -> PyResult<Vec<(f64, f64, f64)>> {
        let flat = face_flat_index_first_shell(&self.inner, face_idx)?;
        let surf_id = self
            .inner
            .geom
            .face_surface
            .get(flat)
            .copied()
            .flatten()
            .ok_or_else(|| {
                PyValueError::new_err("face has no associated surface in GeomStore")
            })?;
        let surface = self
            .inner
            .geom
            .surfaces
            .get(surf_id)
            .ok_or_else(|| PyValueError::new_err("surface index out of range"))?;
        let face = self
            .inner
            .solids
            .first()
            .and_then(|s| s.shells.first())
            .and_then(|sh| sh.faces.get(face_idx))
            .ok_or_else(|| PyValueError::new_err("face_idx out of range"))?;
        let (_wire, pts) = project_wire_onto_surface(
            &self.inner,
            &face.outer_wire,
            surface,
            n_samples,
        )
        .map_err(py_build_err)?;
        Ok(pts.into_iter().map(tuple3).collect())
    }

    /// Translate in-place by ``(dx, dy, dz)``.
    fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        self.inner
            .apply_transform(DAffine3::from_translation(DVec3::new(dx, dy, dz)));
    }

    /// Rotate in-place around ``origin`` and ``axis`` (axis direction; need not be unit); ``angle`` in radians.
    fn rotate_axis_angle(
        &mut self,
        origin: (f64, f64, f64),
        axis: (f64, f64, f64),
        angle: f64,
    ) -> PyResult<()> {
        let o = vec3_tuple(origin);
        let a = vec3_tuple(axis);
        if a.length_squared() < 1e-30 {
            return Err(PyValueError::new_err("axis must be non-zero"));
        }
        let a = a.normalize();
        let q = DQuat::from_axis_angle(a, angle);
        let t0 = DAffine3::from_translation(-o);
        let r = DAffine3::from_quat(q);
        let t1 = DAffine3::from_translation(o);
        self.inner.apply_transform(t1 * r * t0);
        Ok(())
    }

    /// Uniform scale in-place about ``center`` (default world origin).
    #[pyo3(signature = (factor, center=None))]
    fn scale_uniform(
        &mut self,
        factor: f64,
        center: Option<(f64, f64, f64)>,
    ) -> PyResult<()> {
        if !factor.is_finite() || factor <= 0.0 {
            return Err(PyValueError::new_err("factor must be finite and positive"));
        }
        let c = center.map(vec3_tuple).unwrap_or(DVec3::ZERO);
        let s = DAffine3::from_scale(DVec3::splat(factor));
        let t0 = DAffine3::from_translation(-c);
        let t1 = DAffine3::from_translation(c);
        self.inner.apply_transform(t1 * s * t0);
        Ok(())
    }

    /// Axis-aligned bounding box ``((min_x,min_y,min_z), (max_x,max_y,max_z))`` or ``None`` if empty.
    fn bounding_box(&self) -> Option<((f64, f64, f64), (f64, f64, f64))> {
        self.inner
            .bounding_box()
            .map(|[mn, mx]| (tuple3(mn), tuple3(mx)))
    }

    /// Volume centroid (polyhedral approximation); falls back to vertex average if volume ~ 0.
    fn centroid(&self) -> (f64, f64, f64) {
        tuple3(centroid(&self.inner))
    }

    /// Signed volume (can be negative depending on face orientation).
    fn signed_volume(&self) -> f64 {
        signed_volume(&self.inner)
    }

    /// Average of all vertex positions.
    fn center_of_vertices(&self) -> (f64, f64, f64) {
        tuple3(self.inner.center())
    }

    /// Inertia tensor (unit density) about the origin as a 3×3 row-major nested list.
    fn inertia_tensor(&self) -> [[f64; 3]; 3] {
        inertia_tensor(&self.inner).to_matrix()
    }

    fn vertex_count(&self) -> usize {
        self.inner.vertices.len()
    }

    fn edge_count(&self) -> usize {
        self.inner.edges.len()
    }

    fn solid_count(&self) -> usize {
        self.inner.solids.len()
    }

    fn face_count(&self) -> usize {
        face_count(&self.inner)
    }

    /// Total surface area (model units²).
    fn surface_area(&self) -> f64 {
        surface_area(&self.inner)
    }

    /// Sum of solid volumes (model units³); compound shells may contribute zero.
    fn volume(&self) -> f64 {
        volume(&self.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "BRep(solids={}, edges={}, vertices={})",
            self.inner.solids.len(),
            self.inner.edges.len(),
            self.inner.vertices.len()
        )
    }
}

/// Effective fuzzy tolerance used inside the boolean DS (`max(configured, TOLERANCE_ABS)`).
#[pyfunction]
#[pyo3(name = "resolved_boolean_fuzzy_tol")]
fn py_resolved_boolean_fuzzy_tol(configured: f64) -> f64 {
    resolved_boolean_fuzzy_tol_for_ds(configured)
}

/// Maximum stored face tolerance in flat face order, or `TOLERANCE_ABS` if unset.
#[pyfunction]
#[pyo3(name = "max_face_tolerance")]
fn py_max_face_tolerance(brep: &PyBRep) -> f64 {
    tolerance::max_face_tolerance_or_abs(&brep.inner)
}

/// Cylindrical projection of wire-like `shape` onto `target` along `direction`.
///
/// Returns a list of `BRep` results (one per intersection chain). Options match
/// `BrepProjOptions`: `tolerance` clamps to at least `TOLERANCE_ABS` for triangle intersection.
#[pyfunction]
#[pyo3(name = "brep_proj_cylindrical")]
#[pyo3(signature = (shape, target, direction, *, tolerance=None, samples_per_edge=None))]
fn py_brep_proj_cylindrical(
    py: Python<'_>,
    shape: &PyBRep,
    target: &PyBRep,
    direction: (f64, f64, f64),
    tolerance: Option<f64>,
    samples_per_edge: Option<usize>,
) -> PyResult<Vec<Py<PyBRep>>> {
    let mut opts = BrepProjOptions::default();
    if let Some(t) = tolerance {
        opts.tolerance = t;
    }
    if let Some(n) = samples_per_edge {
        opts.samples_per_edge = n.max(2);
    }
    let dir = vec3_tuple(direction);
    if !dir.x.is_finite() || !dir.y.is_finite() || !dir.z.is_finite() {
        return Err(PyValueError::new_err("direction must be finite"));
    }
    if dir.length_squared() < 1e-30 {
        return Err(PyValueError::new_err("direction must be non-zero"));
    }
    let parts = brep_proj_cylindrical(&shape.inner, &target.inner, dir, &opts);
    let mut out = Vec::with_capacity(parts.len());
    for b in parts {
        out.push(Py::new(py, PyBRep { inner: b })?);
    }
    Ok(out)
}

#[pymodule]
fn _rcad(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBRep>()?;
    m.add_class::<PyBooleanOptions>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("TOLERANCE_ABS", tolerance::TOLERANCE_ABS)?;
    m.add_function(wrap_pyfunction!(py_resolved_boolean_fuzzy_tol, m)?)?;
    m.add_function(wrap_pyfunction!(py_max_face_tolerance, m)?)?;
    m.add_function(wrap_pyfunction!(py_brep_proj_cylindrical, m)?)?;
    Ok(())
}
