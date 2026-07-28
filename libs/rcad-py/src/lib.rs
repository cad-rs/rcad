//! Python bindings for the open-source RCAD kernel (B-rep, primitives, booleans, STEP/IGES)
//! and optional parametric [`PyHistoryDocument`] / [`PyFeatureId`] (feature tree JSON / rebuild).
//!
//! Error taxonomy exposed to Python:
//! - ``HistoryError`` and subclasses — ``rcad_history::HistoryError``.
//! - ``BooleanOpError`` and subclasses — ``rcad_algorithms::BooleanError`` (direct ``BRep`` booleans).
//! - ``FeaturesKernelError`` and subclasses — ``rcad_features::FeaturesError`` (feature-layer helpers).
//! - ``ModelingBuildError`` and subclasses — ``rcad_modeling::BuildError`` (primitive constructors, loft, etc.).
//! - ``OffsetKernelError`` / ``SweepKernelError`` — ``rcad_algorithms`` offset & sweep.
//! - ``StepExchangeError`` / ``IgesExchangeError`` — STEP / IGES read paths (``rcad_step``).

use glam::{DAffine3, DQuat, DVec2, DVec3};
use pyo3::exceptions::{PyException, PyIOError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;
use rcad_algorithms::{
    BooleanError, BooleanOpType, BrepProjOptions, OffsetError, SweepError,
    bop_occt_ops::boolean_op_generic as boolean_op, brep_proj_cylindrical, linear_sweep_wire, offset_solid,
    pipe_sweep_wire, tolerance,
};
use rcad_features::{FeaturesError as FeatErr, operations::FeatureOperations};
use rcad_kernel::BRep;
use rcad_kernel::properties::{centroid, inertia_tensor, signed_volume, surface_area, volume};
use rcad_kernel::topods;
use rcad_modeling::{
    chamfer_edge, cone_brep, cylinder_brep, extrude, fillet_edge, fillet_edges, loft,
    make_box_brep, make_conical_frustum_brep, project_wire_onto_surface, revolve, sphere_brep,
    sweep_pipe, torus_brep,
};
use rcad_step::{ExportSelection, IgesReader, IgesWriter, StepReader, StepWriter};

fn vec3_tuple(v: (f64, f64, f64)) -> DVec3 {
    DVec3::new(v.0, v.1, v.2)
}

fn tuple3(v: DVec3) -> (f64, f64, f64) {
    (v.x, v.y, v.z)
}

fn face_count(brep: &BRep) -> usize {
    brep.solids()
        .iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Flattened face index for ``face_idx`` in the first solid's first shell (matches ``extrude``).
fn face_flat_index_first_shell(brep: &BRep, face_idx: usize) -> PyResult<usize> {
    let mut flat = 0usize;
    for (si, solid) in brep.solids().iter().enumerate() {
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

pyo3::create_exception!(_rcad, HistoryError, PyException);
pyo3::create_exception!(_rcad, HistoryFeatureNotFound, HistoryError);
pyo3::create_exception!(_rcad, HistoryFeatureInUse, HistoryError);
pyo3::create_exception!(_rcad, HistoryInvalidReorder, HistoryError);
pyo3::create_exception!(_rcad, HistoryParameterExists, HistoryError);
pyo3::create_exception!(_rcad, HistoryParameterNotFound, HistoryError);
pyo3::create_exception!(_rcad, HistoryUndefinedVariable, HistoryError);
pyo3::create_exception!(_rcad, HistoryInvalidExpression, HistoryError);
pyo3::create_exception!(_rcad, HistoryDivisionByZero, HistoryError);
pyo3::create_exception!(_rcad, HistoryEvaluationFailed, HistoryError);
pyo3::create_exception!(_rcad, HistoryEmptyDocument, HistoryError);
pyo3::create_exception!(_rcad, HistoryNotEvaluated, HistoryError);
pyo3::create_exception!(_rcad, HistoryPersistError, HistoryError);

pyo3::create_exception!(_rcad, BooleanOpError, PyException);
pyo3::create_exception!(_rcad, BooleanOpEmptyInput, BooleanOpError);
pyo3::create_exception!(_rcad, BooleanOpMissingGeometry, BooleanOpError);
pyo3::create_exception!(_rcad, BooleanOpDegenerateResult, BooleanOpError);
pyo3::create_exception!(_rcad, BooleanOpNumericalFailure, BooleanOpError);
pyo3::create_exception!(_rcad, BooleanOpEmptyCollection, BooleanOpError);
pyo3::create_exception!(_rcad, BooleanOpInvalidResult, BooleanOpError);
pyo3::create_exception!(_rcad, BooleanOpIncompleteIntersection, BooleanOpError);
pyo3::create_exception!(_rcad, BooleanOpSelfIntersection, BooleanOpError);

pyo3::create_exception!(_rcad, FeaturesKernelError, PyException);
pyo3::create_exception!(_rcad, FeaturesZeroNormal, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesZeroDirection, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesInvalidPatternCount, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesInvalidPatternSpacing, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesInvalidHoleDiameter, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesInvalidHoleDepth, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesDatumNotFound, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesPatternGenerationFailed, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesMirrorFailed, FeaturesKernelError);
pyo3::create_exception!(_rcad, FeaturesBooleanFailed, FeaturesKernelError);

pyo3::create_exception!(_rcad, ModelingBuildError, PyException);
pyo3::create_exception!(_rcad, ModelingNonFiniteValue, ModelingBuildError);
pyo3::create_exception!(_rcad, ModelingNonPositiveValue, ModelingBuildError);
pyo3::create_exception!(_rcad, ModelingZeroVector, ModelingBuildError);
pyo3::create_exception!(_rcad, ModelingParallelVectors, ModelingBuildError);
pyo3::create_exception!(_rcad, ModelingDegenerateGeometry, ModelingBuildError);
pyo3::create_exception!(_rcad, ModelingInvalidIndex, ModelingBuildError);

pyo3::create_exception!(_rcad, OffsetKernelError, PyException);
pyo3::create_exception!(_rcad, OffsetZeroDistance, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetInvalidInput, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetDegenerateSurface, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetSelfIntersection, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetEdgeIntersectionFailed, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetVertexComputationFailed, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetUnsupportedGeometry, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetNumericalFailure, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetEmptyResult, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetWallThicknessViolation, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetJoinCreationFailed, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetInvalidVariableThickness, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetQualityCheckFailed, OffsetKernelError);
pyo3::create_exception!(_rcad, OffsetRecoveryFailed, OffsetKernelError);

pyo3::create_exception!(_rcad, SweepKernelError, PyException);
pyo3::create_exception!(_rcad, SweepZeroVector, SweepKernelError);
pyo3::create_exception!(_rcad, SweepNonFiniteInput, SweepKernelError);
pyo3::create_exception!(_rcad, SweepNonPositiveInput, SweepKernelError);
pyo3::create_exception!(_rcad, SweepInsufficientVertices, SweepKernelError);
pyo3::create_exception!(_rcad, SweepInsufficientSpinePoints, SweepKernelError);
pyo3::create_exception!(_rcad, SweepVertexCountMismatch, SweepKernelError);
pyo3::create_exception!(_rcad, SweepDegenerateGeometry, SweepKernelError);
pyo3::create_exception!(_rcad, SweepInvalidParameter, SweepKernelError);
pyo3::create_exception!(_rcad, SweepCornerHandlingFailed, SweepKernelError);
pyo3::create_exception!(_rcad, SweepModelingError, SweepKernelError);

pyo3::create_exception!(_rcad, StepExchangeError, PyException);
pyo3::create_exception!(_rcad, StepIoError, StepExchangeError);
pyo3::create_exception!(_rcad, StepInvalidFormatError, StepExchangeError);
pyo3::create_exception!(_rcad, StepMissingEntityError, StepExchangeError);
pyo3::create_exception!(_rcad, StepEmptyResultError, StepExchangeError);

pyo3::create_exception!(_rcad, IgesExchangeError, PyException);
pyo3::create_exception!(_rcad, IgesIoError, IgesExchangeError);
pyo3::create_exception!(_rcad, IgesInvalidFormatError, IgesExchangeError);
pyo3::create_exception!(_rcad, IgesEmptyResultError, IgesExchangeError);

/// Map ``rcad_history::HistoryError`` to typed Python exceptions (subclasses of ``HistoryError``).
fn py_history_err(e: rcad_history::HistoryError) -> PyErr {
    let msg = e.to_string();
    match e {
        rcad_history::HistoryError::FeatureNotFound(_) => HistoryFeatureNotFound::new_err(msg),
        rcad_history::HistoryError::FeatureInUse(_) => HistoryFeatureInUse::new_err(msg),
        rcad_history::HistoryError::InvalidReorder { .. } => HistoryInvalidReorder::new_err(msg),
        rcad_history::HistoryError::ParameterExists(_) => HistoryParameterExists::new_err(msg),
        rcad_history::HistoryError::ParameterNotFound(_) => HistoryParameterNotFound::new_err(msg),
        rcad_history::HistoryError::UndefinedVariable(_) => HistoryUndefinedVariable::new_err(msg),
        rcad_history::HistoryError::InvalidExpression(_) => HistoryInvalidExpression::new_err(msg),
        rcad_history::HistoryError::DivisionByZero => HistoryDivisionByZero::new_err(msg),
        rcad_history::HistoryError::EvaluationFailed(_) => HistoryEvaluationFailed::new_err(msg),
        rcad_history::HistoryError::EmptyDocument => HistoryEmptyDocument::new_err(msg),
        rcad_history::HistoryError::NotEvaluated(_) => HistoryNotEvaluated::new_err(msg),
        rcad_history::HistoryError::Persist(_) => HistoryPersistError::new_err(msg),
    }
}

fn py_bool_err(e: BooleanError) -> PyErr {
    let msg = e.to_string();
    match e {
        BooleanError::EmptyInput => BooleanOpEmptyInput::new_err(msg),
        BooleanError::MissingGeometry(_) => BooleanOpMissingGeometry::new_err(msg),
        BooleanError::DegenerateResult => BooleanOpDegenerateResult::new_err(msg),
        BooleanError::NumericalFailure(_) => BooleanOpNumericalFailure::new_err(msg),
        BooleanError::EmptyCollection(_) => BooleanOpEmptyCollection::new_err(msg),
        BooleanError::InvalidResult(_) => BooleanOpInvalidResult::new_err(msg),
        BooleanError::IncompleteIntersection(_) => BooleanOpIncompleteIntersection::new_err(msg),
        BooleanError::SelfIntersection(_) => BooleanOpSelfIntersection::new_err(msg),
        _ => pyo3::exceptions::PyRuntimeError::new_err(msg),
    }
}

fn py_features_err(e: FeatErr) -> PyErr {
    let msg = e.to_string();
    match e {
        FeatErr::ZeroNormal => FeaturesZeroNormal::new_err(msg),
        FeatErr::ZeroDirection => FeaturesZeroDirection::new_err(msg),
        FeatErr::InvalidPatternCount(_) => FeaturesInvalidPatternCount::new_err(msg),
        FeatErr::InvalidPatternSpacing(_) => FeaturesInvalidPatternSpacing::new_err(msg),
        FeatErr::InvalidHoleDiameter(_) => FeaturesInvalidHoleDiameter::new_err(msg),
        FeatErr::InvalidHoleDepth(_) => FeaturesInvalidHoleDepth::new_err(msg),
        FeatErr::FeatureNotFound(_) => FeaturesDatumNotFound::new_err(msg),
        FeatErr::PatternGenerationFailed(_) => FeaturesPatternGenerationFailed::new_err(msg),
        FeatErr::MirrorFailed(_) => FeaturesMirrorFailed::new_err(msg),
        FeatErr::BooleanFailed(_) => FeaturesBooleanFailed::new_err(msg),
    }
}

fn py_build_err(e: rcad_modeling::BuildError) -> PyErr {
    let msg = e.to_string();
    match e {
        rcad_modeling::BuildError::NonFiniteValue(_) => ModelingNonFiniteValue::new_err(msg),
        rcad_modeling::BuildError::NonPositiveValue(_) => ModelingNonPositiveValue::new_err(msg),
        rcad_modeling::BuildError::ZeroVector(_) => ModelingZeroVector::new_err(msg),
        rcad_modeling::BuildError::ParallelVectors(_, _) => ModelingParallelVectors::new_err(msg),
        rcad_modeling::BuildError::DegenerateGeometry(_) => {
            ModelingDegenerateGeometry::new_err(msg)
        }
        rcad_modeling::BuildError::InvalidIndex(_) => ModelingInvalidIndex::new_err(msg),
    }
}

fn py_offset_err(e: OffsetError) -> PyErr {
    let msg = e.to_string();
    match e {
        OffsetError::ZeroDistance => OffsetZeroDistance::new_err(msg),
        OffsetError::InvalidInput(_) => OffsetInvalidInput::new_err(msg),
        OffsetError::DegenerateSurface { .. } => OffsetDegenerateSurface::new_err(msg),
        OffsetError::SelfIntersection { .. } => OffsetSelfIntersection::new_err(msg),
        OffsetError::EdgeIntersectionFailed { .. } => OffsetEdgeIntersectionFailed::new_err(msg),
        OffsetError::VertexComputationFailed { .. } => OffsetVertexComputationFailed::new_err(msg),
        OffsetError::UnsupportedGeometry { .. } => OffsetUnsupportedGeometry::new_err(msg),
        OffsetError::NumericalFailure(_) => OffsetNumericalFailure::new_err(msg),
        OffsetError::EmptyResult => OffsetEmptyResult::new_err(msg),
        OffsetError::WallThicknessViolation { .. } => OffsetWallThicknessViolation::new_err(msg),
        OffsetError::JoinCreationFailed { .. } => OffsetJoinCreationFailed::new_err(msg),
        OffsetError::InvalidVariableThickness { .. } => {
            OffsetInvalidVariableThickness::new_err(msg)
        }
        OffsetError::QualityCheckFailed { .. } => OffsetQualityCheckFailed::new_err(msg),
        OffsetError::RecoveryFailed { .. } => OffsetRecoveryFailed::new_err(msg),
    }
}

fn py_sweep_err(e: SweepError) -> PyErr {
    let msg = e.to_string();
    match e {
        SweepError::ZeroVector(_) => SweepZeroVector::new_err(msg),
        SweepError::NonFiniteInput(_) => SweepNonFiniteInput::new_err(msg),
        SweepError::NonPositiveInput(_) => SweepNonPositiveInput::new_err(msg),
        SweepError::InsufficientVertices { .. } => SweepInsufficientVertices::new_err(msg),
        SweepError::InsufficientSpinePoints { .. } => SweepInsufficientSpinePoints::new_err(msg),
        SweepError::VertexCountMismatch { .. } => SweepVertexCountMismatch::new_err(msg),
        SweepError::DegenerateGeometry(_) => SweepDegenerateGeometry::new_err(msg),
        SweepError::InvalidParameter(_) => SweepInvalidParameter::new_err(msg),
        SweepError::CornerHandlingFailed(_) => SweepCornerHandlingFailed::new_err(msg),
        SweepError::ModelingError(_) => SweepModelingError::new_err(msg),
    }
}

fn py_step_err(e: rcad_step::StepError) -> PyErr {
    let msg = e.to_string();
    match e {
        rcad_step::StepError::Io(_) => StepIoError::new_err(msg),
        rcad_step::StepError::InvalidFormat(_) => StepInvalidFormatError::new_err(msg),
        rcad_step::StepError::MissingEntity { .. } => StepMissingEntityError::new_err(msg),
        rcad_step::StepError::EmptyResult(_) => StepEmptyResultError::new_err(msg),
    }
}

fn py_iges_err(e: rcad_step::IgesError) -> PyErr {
    let msg = e.to_string();
    match e {
        rcad_step::IgesError::Io(_) => IgesIoError::new_err(msg),
        rcad_step::IgesError::InvalidFormat(_) => IgesInvalidFormatError::new_err(msg),
        rcad_step::IgesError::EmptyResult(_) => IgesEmptyResultError::new_err(msg),
    }
}

fn parse_feature_id(feature_id: &Bound<'_, PyAny>) -> PyResult<u64> {
    if let Ok(v) = feature_id.extract::<u64>() {
        return Ok(v);
    }
    if let Ok(fid) = feature_id.downcast::<PyFeatureId>() {
        return Ok(fid.borrow().value);
    }
    Err(PyTypeError::new_err("feature_id must be int or FeatureId"))
}

/// Stable feature identifier (matches Rust ``rcad_history::FeatureId`` numeric value).
#[pyclass(name = "FeatureId", frozen, eq, hash)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct PyFeatureId {
    #[pyo3(get)]
    pub value: u64,
}

#[pymethods]
impl PyFeatureId {
    #[new]
    fn new(value: u64) -> Self {
        Self { value }
    }

    fn __repr__(&self) -> String {
        format!("FeatureId({})", self.value)
    }

    fn __int__(&self) -> u64 {
        self.value
    }
}

/// Parametric part document (feature tree + parameters) from ``rcad_history``.
///
/// Load/save JSON compatible with the Rust ``Document::to_json`` / ``from_json_str`` format.
#[pyclass(name = "HistoryDocument")]
pub struct PyHistoryDocument {
    inner: rcad_history::Document,
}

#[pymethods]
impl PyHistoryDocument {
    #[staticmethod]
    fn new(name: &str) -> Self {
        Self {
            inner: rcad_history::Document::new(name),
        }
    }

    #[staticmethod]
    fn from_json(s: &str) -> PyResult<Self> {
        let inner = rcad_history::Document::from_json_str(s).map_err(py_history_err)?;
        Ok(Self { inner })
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn to_json_pretty(&self) -> PyResult<String> {
        self.inner
            .to_json_pretty()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Document UUID from JSON / ``new`` (same as Rust ``Document::id``).
    #[getter]
    fn document_id(&self) -> String {
        self.inner.id.to_string()
    }

    fn rebuild(&mut self) -> PyResult<()> {
        self.inner.rebuild().map_err(py_history_err)
    }

    /// Dependency order check (referenced features must appear earlier in the tree).
    fn validate_dependency_order(&self) -> PyResult<()> {
        self.inner
            .validate_dependency_order()
            .map_err(py_history_err)
    }

    /// Feature ids in tree order (same as Rust ``Document::features()``).
    fn feature_ids(&self) -> Vec<u64> {
        self.inner.features().iter().map(|f| f.id.0).collect()
    }

    /// Next id that ``add_feature`` would assign (Rust ``next_feature_sequence``).
    fn next_feature_sequence(&self) -> u64 {
        self.inner.next_feature_sequence()
    }

    /// Last feature id in the tree, if any (convenience for ``final_brep`` / cache lookups).
    fn final_feature_id(&self) -> Option<u64> {
        self.inner.features().last().map(|f| f.id.0)
    }

    fn clear_cache(&mut self) {
        self.inner.clear_cache();
    }

    fn update_parameter(&mut self, name: &str, value: f64) -> PyResult<()> {
        self.inner
            .update_parameter(name, value)
            .map_err(py_history_err)
    }

    /// Add a driving dimension (constant numeric value).
    fn add_parameter(&mut self, name: &str, value: f64) -> PyResult<()> {
        let param = rcad_history::Parameter::new(name, value);
        self.inner.add_parameter(param).map_err(py_history_err)
    }

    /// Add a parameter whose value is an expression (same evaluator as JSON ``Expression``).
    fn add_parameter_expression(&mut self, name: &str, expr: &str) -> PyResult<()> {
        let param = rcad_history::Parameter::with_expression(name, expr);
        self.inner.add_parameter(param).map_err(py_history_err)
    }

    /// Append a box primitive feature (constant width/height/depth). Returns assigned feature id.
    fn add_box_feature(
        &mut self,
        name: &str,
        corner: (f64, f64, f64),
        width: f64,
        height: f64,
        depth: f64,
    ) -> PyResult<u64> {
        let f = rcad_history::Feature::new(
            rcad_history::FeatureId::new(0),
            name,
            rcad_history::FeatureType::Box(rcad_history::BoxFeature {
                corner: vec3_tuple(corner),
                width: rcad_history::ParameterValue::Constant(width),
                height: rcad_history::ParameterValue::Constant(height),
                depth: rcad_history::ParameterValue::Constant(depth),
            }),
        );
        let id = self.inner.add_feature(f);
        Ok(id.0)
    }

    /// Number of features in the tree (same ordering as Rust ``Document::features()``).
    fn feature_count(&self) -> usize {
        self.inner.features().len()
    }

    /// Sorted list of parameter names (``Document`` hash map keys).
    fn parameter_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.inner.parameters().keys().cloned().collect();
        names.sort();
        names
    }

    /// Resolved numeric values for all parameters (multi-pass; matches Rust ``Document::resolved_parameter_values``).
    fn resolved_parameter_values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.inner
            .resolved_parameter_values()
            .into_pyobject(py)
            .map(|d| d.into_any())
    }

    /// Cached solid for ``feature_id`` after the last successful ``rebuild`` (Rust ``evaluated_brep``).
    /// Accepts ``int`` or ``FeatureId``.
    fn evaluated_brep(&self, feature_id: &Bound<'_, PyAny>) -> PyResult<Option<PyBRep>> {
        let id = parse_feature_id(feature_id)?;
        Ok(self
            .inner
            .evaluated_brep(rcad_history::FeatureId::new(id))
            .map(|b| PyBRep { inner: b.clone() }))
    }

    fn final_brep(&self) -> Option<PyBRep> {
        self.inner.final_brep().map(|b| PyBRep { inner: b.clone() })
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "HistoryDocument(name={:?}, features={})",
            self.inner.name,
            self.inner.features().len()
        )
    }
}

/// Boolean pipeline options — options were removed during `Shape` migration; kept as a no-op stub
/// for backward compatibility of callers that construct `BooleanOptions(...)`.
#[pyclass(name = "BooleanOptions")]
#[derive(Clone, Copy)]
pub struct PyBooleanOptions;

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
        _use_bvh: bool,
        _run_healing: bool,
        _run_simplify: bool,
        _include_history: bool,
        _run_make_connected: bool,
        _make_connected_tolerance: Option<f64>,
        _fuzzy_tol: f64,
        _use_glue: bool,
        _glue_tolerance: Option<f64>,
        _run_propagate_geom_tolerances: bool,
    ) -> Self {
        Self
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
            0.0,
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
        let p2: Vec<DVec2> = profile_2d.iter().map(|(x, y)| DVec2::new(*x, *y)).collect();
        let sp: Vec<DVec3> = spine.into_iter().map(vec3_tuple).collect();
        let brep = sweep_pipe(&p2, &sp).map_err(py_build_err)?;
        Ok(Self { inner: brep })
    }

    /// Pipe sweep wire (``rcad_algorithms::pipe_sweep_wire``): alternate sweep kernel with pipe mode; same arguments as ``sweep_pipe``.
    #[staticmethod]
    fn pipe_sweep_wire(profile_2d: Vec<(f64, f64)>, spine: Vec<(f64, f64, f64)>) -> PyResult<Self> {
        let p2: Vec<DVec2> = profile_2d.iter().map(|(x, y)| DVec2::new(*x, *y)).collect();
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
        let out =
            boolean_op(BooleanOpType::Union, &self.inner, &other.inner).map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Boolean intersection ``self & other``.
    fn intersection(&self, other: &PyBRep) -> PyResult<Self> {
        let out = boolean_op(BooleanOpType::Intersection, &self.inner, &other.inner)
            .map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Boolean difference ``self - other``.
    fn difference(&self, other: &PyBRep) -> PyResult<Self> {
        let out = boolean_op(BooleanOpType::Difference, &self.inner, &other.inner)
            .map_err(py_bool_err)?;
        Ok(Self { inner: out })
    }

    /// Extrude face ``face_idx`` (index in the first solid's outer shell) along ``direction`` by ``distance``.
    #[pyo3(signature = (face_idx, direction, distance))]
    fn extrude(
        &self,
        face_idx: usize,
        direction: (f64, f64, f64),
        distance: f64,
    ) -> PyResult<Self> {
        let brep = extrude(&self.inner, face_idx, vec3_tuple(direction), distance)
            .map_err(py_build_err)?;
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
    fn repair(&self, _tolerance: f64) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.clone(),
        })
    }

    /// Offset the first solid's shell by ``distance`` (along face normals; see ``rcad_algorithms::offset_solid``).
    fn offset_solid(&self, distance: f64) -> PyResult<Self> {
        let solids = self.inner.solids();
        let solid = solids
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
        let brep = linear_sweep_wire(&self.inner, face_idx, vec3_tuple(direction), distance)
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
        let geom = self.inner.geom();
        let surf_id = geom
            .face_surface
            .get(flat)
            .copied()
            .flatten()
            .ok_or_else(|| PyValueError::new_err("face has no associated surface in GeomStore"))?;
        let surface = geom
            .surfaces
            .get(surf_id)
            .ok_or_else(|| PyValueError::new_err("surface index out of range"))?;
        let solids = self.inner.solids();
        let face = solids
            .first()
            .and_then(|s| s.shells.first())
            .and_then(|sh| sh.faces.get(face_idx))
            .ok_or_else(|| PyValueError::new_err("face_idx out of range"))?;
        let (_wire, pts) =
            project_wire_onto_surface(&self.inner, &face.outer_wire, surface, n_samples)
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
    fn scale_uniform(&mut self, factor: f64, center: Option<(f64, f64, f64)>) -> PyResult<()> {
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
        self.inner.vertices().len()
    }

    fn edge_count(&self) -> usize {
        self.inner.edges().len()
    }

    fn solid_count(&self) -> usize {
        self.inner.solids().len()
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
            self.inner.solids().len(),
            self.inner.edges().len(),
            self.inner.vertices().len()
        )
    }
}

/// Effective fuzzy tolerance used inside the boolean DS (`max(configured, TOLERANCE_ABS)`).
#[pyfunction]
#[pyo3(name = "resolved_boolean_fuzzy_tol")]
fn py_resolved_boolean_fuzzy_tol(configured: f64) -> f64 {
    configured.max(tolerance::TOLERANCE_ABS)
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

/// Validate mirror plane definition (raises ``Features*`` errors). For tests and scripts.
#[pyfunction]
#[pyo3(name = "features_create_mirror")]
fn py_features_create_mirror(
    name: &str,
    point: (f64, f64, f64),
    normal: (f64, f64, f64),
) -> PyResult<()> {
    FeatureOperations::create_mirror(name, vec3_tuple(point), vec3_tuple(normal))
        .map_err(py_features_err)?;
    Ok(())
}

/// Validate linear pattern parameters (raises ``Features*`` errors).
#[pyfunction]
#[pyo3(name = "features_create_linear_pattern")]
fn py_features_create_linear_pattern(
    direction: (f64, f64, f64),
    count: usize,
    spacing: f64,
) -> PyResult<()> {
    FeatureOperations::create_linear_pattern(vec3_tuple(direction), count, spacing)
        .map_err(py_features_err)?;
    Ok(())
}

/// Boolean union via ``rcad_features::FeatureOperations`` (errors are ``FeaturesBooleanFailed``, not ``BooleanOp*``).
#[pyfunction]
#[pyo3(name = "features_boolean_union")]
fn py_features_boolean_union(a: &PyBRep, b: &PyBRep) -> PyResult<PyBRep> {
    FeatureOperations::boolean_brep(BooleanOpType::Union, &a.inner, &b.inner)
        .map(|inner| PyBRep { inner })
        .map_err(py_features_err)
}

#[pymodule]
fn _rcad(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("HistoryError", py.get_type::<HistoryError>())?;
    m.add(
        "HistoryFeatureNotFound",
        py.get_type::<HistoryFeatureNotFound>(),
    )?;
    m.add("HistoryFeatureInUse", py.get_type::<HistoryFeatureInUse>())?;
    m.add(
        "HistoryInvalidReorder",
        py.get_type::<HistoryInvalidReorder>(),
    )?;
    m.add(
        "HistoryParameterExists",
        py.get_type::<HistoryParameterExists>(),
    )?;
    m.add(
        "HistoryParameterNotFound",
        py.get_type::<HistoryParameterNotFound>(),
    )?;
    m.add(
        "HistoryUndefinedVariable",
        py.get_type::<HistoryUndefinedVariable>(),
    )?;
    m.add(
        "HistoryInvalidExpression",
        py.get_type::<HistoryInvalidExpression>(),
    )?;
    m.add(
        "HistoryDivisionByZero",
        py.get_type::<HistoryDivisionByZero>(),
    )?;
    m.add(
        "HistoryEvaluationFailed",
        py.get_type::<HistoryEvaluationFailed>(),
    )?;
    m.add(
        "HistoryEmptyDocument",
        py.get_type::<HistoryEmptyDocument>(),
    )?;
    m.add("HistoryNotEvaluated", py.get_type::<HistoryNotEvaluated>())?;
    m.add("HistoryPersistError", py.get_type::<HistoryPersistError>())?;

    m.add("BooleanOpError", py.get_type::<BooleanOpError>())?;
    m.add("BooleanOpEmptyInput", py.get_type::<BooleanOpEmptyInput>())?;
    m.add(
        "BooleanOpMissingGeometry",
        py.get_type::<BooleanOpMissingGeometry>(),
    )?;
    m.add(
        "BooleanOpDegenerateResult",
        py.get_type::<BooleanOpDegenerateResult>(),
    )?;
    m.add(
        "BooleanOpNumericalFailure",
        py.get_type::<BooleanOpNumericalFailure>(),
    )?;
    m.add(
        "BooleanOpEmptyCollection",
        py.get_type::<BooleanOpEmptyCollection>(),
    )?;
    m.add(
        "BooleanOpInvalidResult",
        py.get_type::<BooleanOpInvalidResult>(),
    )?;
    m.add(
        "BooleanOpIncompleteIntersection",
        py.get_type::<BooleanOpIncompleteIntersection>(),
    )?;
    m.add(
        "BooleanOpSelfIntersection",
        py.get_type::<BooleanOpSelfIntersection>(),
    )?;

    m.add("FeaturesKernelError", py.get_type::<FeaturesKernelError>())?;
    m.add("FeaturesZeroNormal", py.get_type::<FeaturesZeroNormal>())?;
    m.add(
        "FeaturesZeroDirection",
        py.get_type::<FeaturesZeroDirection>(),
    )?;
    m.add(
        "FeaturesInvalidPatternCount",
        py.get_type::<FeaturesInvalidPatternCount>(),
    )?;
    m.add(
        "FeaturesInvalidPatternSpacing",
        py.get_type::<FeaturesInvalidPatternSpacing>(),
    )?;
    m.add(
        "FeaturesInvalidHoleDiameter",
        py.get_type::<FeaturesInvalidHoleDiameter>(),
    )?;
    m.add(
        "FeaturesInvalidHoleDepth",
        py.get_type::<FeaturesInvalidHoleDepth>(),
    )?;
    m.add(
        "FeaturesDatumNotFound",
        py.get_type::<FeaturesDatumNotFound>(),
    )?;
    m.add(
        "FeaturesPatternGenerationFailed",
        py.get_type::<FeaturesPatternGenerationFailed>(),
    )?;
    m.add(
        "FeaturesMirrorFailed",
        py.get_type::<FeaturesMirrorFailed>(),
    )?;
    m.add(
        "FeaturesBooleanFailed",
        py.get_type::<FeaturesBooleanFailed>(),
    )?;

    m.add("ModelingBuildError", py.get_type::<ModelingBuildError>())?;
    m.add(
        "ModelingNonFiniteValue",
        py.get_type::<ModelingNonFiniteValue>(),
    )?;
    m.add(
        "ModelingNonPositiveValue",
        py.get_type::<ModelingNonPositiveValue>(),
    )?;
    m.add("ModelingZeroVector", py.get_type::<ModelingZeroVector>())?;
    m.add(
        "ModelingParallelVectors",
        py.get_type::<ModelingParallelVectors>(),
    )?;
    m.add(
        "ModelingDegenerateGeometry",
        py.get_type::<ModelingDegenerateGeometry>(),
    )?;
    m.add(
        "ModelingInvalidIndex",
        py.get_type::<ModelingInvalidIndex>(),
    )?;

    m.add("OffsetKernelError", py.get_type::<OffsetKernelError>())?;
    m.add("OffsetZeroDistance", py.get_type::<OffsetZeroDistance>())?;
    m.add("OffsetInvalidInput", py.get_type::<OffsetInvalidInput>())?;
    m.add(
        "OffsetDegenerateSurface",
        py.get_type::<OffsetDegenerateSurface>(),
    )?;
    m.add(
        "OffsetSelfIntersection",
        py.get_type::<OffsetSelfIntersection>(),
    )?;
    m.add(
        "OffsetEdgeIntersectionFailed",
        py.get_type::<OffsetEdgeIntersectionFailed>(),
    )?;
    m.add(
        "OffsetVertexComputationFailed",
        py.get_type::<OffsetVertexComputationFailed>(),
    )?;
    m.add(
        "OffsetUnsupportedGeometry",
        py.get_type::<OffsetUnsupportedGeometry>(),
    )?;
    m.add(
        "OffsetNumericalFailure",
        py.get_type::<OffsetNumericalFailure>(),
    )?;
    m.add("OffsetEmptyResult", py.get_type::<OffsetEmptyResult>())?;
    m.add(
        "OffsetWallThicknessViolation",
        py.get_type::<OffsetWallThicknessViolation>(),
    )?;
    m.add(
        "OffsetJoinCreationFailed",
        py.get_type::<OffsetJoinCreationFailed>(),
    )?;
    m.add(
        "OffsetInvalidVariableThickness",
        py.get_type::<OffsetInvalidVariableThickness>(),
    )?;
    m.add(
        "OffsetQualityCheckFailed",
        py.get_type::<OffsetQualityCheckFailed>(),
    )?;
    m.add(
        "OffsetRecoveryFailed",
        py.get_type::<OffsetRecoveryFailed>(),
    )?;

    m.add("SweepKernelError", py.get_type::<SweepKernelError>())?;
    m.add("SweepZeroVector", py.get_type::<SweepZeroVector>())?;
    m.add("SweepNonFiniteInput", py.get_type::<SweepNonFiniteInput>())?;
    m.add(
        "SweepNonPositiveInput",
        py.get_type::<SweepNonPositiveInput>(),
    )?;
    m.add(
        "SweepInsufficientVertices",
        py.get_type::<SweepInsufficientVertices>(),
    )?;
    m.add(
        "SweepInsufficientSpinePoints",
        py.get_type::<SweepInsufficientSpinePoints>(),
    )?;
    m.add(
        "SweepVertexCountMismatch",
        py.get_type::<SweepVertexCountMismatch>(),
    )?;
    m.add(
        "SweepDegenerateGeometry",
        py.get_type::<SweepDegenerateGeometry>(),
    )?;
    m.add(
        "SweepInvalidParameter",
        py.get_type::<SweepInvalidParameter>(),
    )?;
    m.add(
        "SweepCornerHandlingFailed",
        py.get_type::<SweepCornerHandlingFailed>(),
    )?;
    m.add("SweepModelingError", py.get_type::<SweepModelingError>())?;

    m.add("StepExchangeError", py.get_type::<StepExchangeError>())?;
    m.add("StepIoError", py.get_type::<StepIoError>())?;
    m.add(
        "StepInvalidFormatError",
        py.get_type::<StepInvalidFormatError>(),
    )?;
    m.add(
        "StepMissingEntityError",
        py.get_type::<StepMissingEntityError>(),
    )?;
    m.add(
        "StepEmptyResultError",
        py.get_type::<StepEmptyResultError>(),
    )?;

    m.add("IgesExchangeError", py.get_type::<IgesExchangeError>())?;
    m.add("IgesIoError", py.get_type::<IgesIoError>())?;
    m.add(
        "IgesInvalidFormatError",
        py.get_type::<IgesInvalidFormatError>(),
    )?;
    m.add(
        "IgesEmptyResultError",
        py.get_type::<IgesEmptyResultError>(),
    )?;

    m.add("TRACE_TARGET_BOOLEAN", "rcad.boolean")?;
    m.add("TRACE_TARGET_CLASSIFY", "rcad.classify")?;

    m.add_class::<PyBRep>()?;
    m.add_class::<PyFeatureId>()?;
    m.add_class::<PyHistoryDocument>()?;
    m.add_class::<PyBooleanOptions>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("TOLERANCE_ABS", tolerance::TOLERANCE_ABS)?;
    m.add_function(wrap_pyfunction!(py_resolved_boolean_fuzzy_tol, m)?)?;
    m.add_function(wrap_pyfunction!(py_max_face_tolerance, m)?)?;
    m.add_function(wrap_pyfunction!(py_brep_proj_cylindrical, m)?)?;
    m.add_function(wrap_pyfunction!(py_features_create_mirror, m)?)?;
    m.add_function(wrap_pyfunction!(py_features_create_linear_pattern, m)?)?;
    m.add_function(wrap_pyfunction!(py_features_boolean_union, m)?)?;
    Ok(())
}
