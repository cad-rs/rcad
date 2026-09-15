//! OCCT TKMath/Poly package — the Poly triangulation data layer (A0 core subset).
//!
//! OCCT layering note: the Poly package lives in TKMath (FoundationClasses),
//! BELOW TKBRep in the OCCT dependency graph — `BRep_TFace` stores
//! `Poly_Triangulation` handles (BRep_TFace.hxx L22, L139-140) and
//! `BRep_PolygonOnTriangulation` stores `Poly_PolygonOnTriangulation`
//! (BRep_PolygonOnTriangulation.hxx L24-25). rcad mirrors the same direction:
//! these value types live in rcad-kernel (which plays both the
//! FoundationClasses role and the TKBRep TShape-pool role) so that `TFaceData`
//! and `CurveRepresentation` can carry them by value/pool index;
//! `rcad_brep::poly` is the TKBRep-facing module (re-exports, Poly_Connect,
//! BRep mount accessors, persistence).
//!
//! Value-type mapping (AGENTS.md 1:1 rules):
//! - `gp_Pnt`                -> `glam::DVec3`
//! - `gp_Pnt2d`              -> `glam::DVec2`
//! - `NCollection_Vec3<float>` -> `[f32; 3]` (OCCT normals are single precision)
//! - `occ::handle<T>`        -> value / pool index (handle identity = index identity)
//!
//! Not carried over from the OCCT package (documented omissions, see the
//! per-item notes below): `Poly_CoherentTriangulation` family,
//! `Poly_MakeLoops`, `Poly_MergeNodesTool` (plan BREPMESH_RMSH_PLAN §1.4
//! non-goals), `Poly_Polygon2D` (2D polygon-on-surface mount, outside A0
//! scope), `DumpJson` (debug dump, no OCCT semantic effect), the
//! single-precision node arrays (`Poly_ArrayOfNodes::SetDoublePrecision` is a
//! memory-layout knob; rcad node arrays are always f64) and the cached
//! min-max machinery (`myCachedMinMax`, a pure optimization layer; rcad has no
//! `Bnd_Box` type — `PolyTriangulation::min_max` always computes accurately).

use glam::{DVec2, DVec3};
use serde::{Deserialize, Serialize};

// =============================================================================
// Poly_MeshPurpose (OCCT Poly_MeshPurpose.hxx L18-34)
// =============================================================================

/// Purpose of triangulation using. OCCT is `typedef unsigned int Poly_MeshPurpose`
/// plus an anonymous enum; the Rust equivalent is a plain alias with the same
/// constants (Poly_MeshPurpose.hxx L18-34).
pub type PolyMeshPurpose = u32;

/// No special use (default). Poly_MeshPurpose.hxx L23.
pub const POLY_MESH_PURPOSE_NONE: PolyMeshPurpose = 0;
/// Mesh for algorithms. Poly_MeshPurpose.hxx L24.
pub const POLY_MESH_PURPOSE_CALCULATION: PolyMeshPurpose = 0x0001;
/// Mesh for presentation (LODs usage). Poly_MeshPurpose.hxx L25.
pub const POLY_MESH_PURPOSE_PRESENTATION: PolyMeshPurpose = 0x0002;
/// Mesh marked as currently active in a list. Poly_MeshPurpose.hxx L27.
pub const POLY_MESH_PURPOSE_ACTIVE: PolyMeshPurpose = 0x0004;
/// Mesh has currently loaded data. Poly_MeshPurpose.hxx L28.
pub const POLY_MESH_PURPOSE_LOADED: PolyMeshPurpose = 0x0008;
/// Special flag for BRep_Tools::Triangulation() to return any other defined
/// mesh if none matching other criteria was found. Poly_MeshPurpose.hxx L29-32.
pub const POLY_MESH_PURPOSE_ANY_FALLBACK: PolyMeshPurpose = 0x0010;
/// Application-defined flags. Poly_MeshPurpose.hxx L33.
pub const POLY_MESH_PURPOSE_USER: PolyMeshPurpose = 0x0020;

// =============================================================================
// Poly_Triangle (OCCT Poly_Triangle.hxx L30-97)
// =============================================================================

/// Describes a component triangle of a triangulation (Poly_Triangulation
/// object). A triangle is defined by a triplet of node indices within
/// [1, Poly_Triangulation::NbNodes()] range. OCCT Poly_Triangle.hxx L26-29.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolyTriangle {
    /// OCCT `int myNodes[3]` (Poly_Triangle.hxx L96).
    pub my_nodes: [i32; 3],
}

impl PolyTriangle {
    /// Constructs a triangle and sets all indices to zero.
    /// OCCT Poly_Triangle.hxx L36.
    pub fn zeroed() -> Self {
        PolyTriangle {
            my_nodes: [0, 0, 0],
        }
    }

    /// Constructs a triangle and sets its three indices, where these node
    /// values are indices in the table of nodes specific to an existing
    /// triangulation of a shape. OCCT Poly_Triangle.hxx L41-46.
    pub fn new(the_n1: i32, the_n2: i32, the_n3: i32) -> Self {
        PolyTriangle {
            my_nodes: [the_n1, the_n2, the_n3],
        }
    }

    /// Sets the value of the three nodes of this triangle.
    /// OCCT Poly_Triangle.hxx L49-54.
    pub fn set(&mut self, the_n1: i32, the_n2: i32, the_n3: i32) {
        self.my_nodes[0] = the_n1;
        self.my_nodes[1] = the_n2;
        self.my_nodes[2] = the_n3;
    }

    /// Sets the value of node with specified index of this triangle.
    /// Raises Standard_OutOfRange if index is not in 1,2,3.
    /// OCCT Poly_Triangle.hxx L58-63.
    pub fn set_value(&mut self, the_index: usize, the_node: i32) {
        assert!(
            (1..=3).contains(&the_index),
            "Poly_Triangle::Set(), invalid index"
        );
        self.my_nodes[the_index - 1] = the_node;
    }

    /// Returns the node indices of this triangle. OCCT Poly_Triangle.hxx L66-71.
    pub fn get(&self) -> (i32, i32, i32) {
        (self.my_nodes[0], self.my_nodes[1], self.my_nodes[2])
    }

    /// Get the node of given index. Raises OutOfRange if Index is not in 1,2,3.
    /// OCCT Poly_Triangle.hxx L75-80.
    pub fn value(&self, the_index: usize) -> i32 {
        assert!(
            (1..=3).contains(&the_index),
            "Poly_Triangle::Value(), invalid index"
        );
        self.my_nodes[the_index - 1]
    }
}

// =============================================================================
// Poly_TriangulationParameters (OCCT Poly_TriangulationParameters.hxx L23-70)
// =============================================================================

/// Represents initial set of parameters triangulation is built for.
/// OCCT Poly_TriangulationParameters.hxx L22-23.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PolyTriangulationParameters {
    /// OCCT `double myDeflection` (hxx L67).
    pub my_deflection: f64,
    /// OCCT `double myAngle` (hxx L68).
    pub my_angle: f64,
    /// OCCT `double myMinSize` (hxx L69).
    pub my_min_size: f64,
}

impl PolyTriangulationParameters {
    /// Constructor. Initializes object with the given parameters.
    /// OCCT Poly_TriangulationParameters.hxx L31-38 (defaults -1.).
    pub fn new(
        the_deflection: f64, // linear deflection
        the_angle: f64,      // angular deflection
        the_min_size: f64,   // minimum size
    ) -> Self {
        PolyTriangulationParameters {
            my_deflection: the_deflection,
            my_angle: the_angle,
            my_min_size: the_min_size,
        }
    }

    /// Creates a copy of current triangulation parameters.
    /// OCCT Poly_TriangulationParameters.cxx Copy() — `new` with the same
    /// three values; the derived `Clone` impl is the direct equivalent.
    pub fn copy(&self) -> Self {
        self.clone()
    }

    /// Returns true if linear deflection is defined.
    /// OCCT Poly_TriangulationParameters.hxx L47.
    pub fn has_deflection(&self) -> bool {
        !(self.my_deflection < 0.)
    }

    /// Returns true if angular deflection is defined.
    /// OCCT Poly_TriangulationParameters.hxx L49.
    pub fn has_angle(&self) -> bool {
        !(self.my_angle < 0.)
    }

    /// Returns true if minimum size is defined.
    /// OCCT Poly_TriangulationParameters.hxx L52.
    pub fn has_min_size(&self) -> bool {
        !(self.my_min_size < 0.)
    }

    /// Returns linear deflection or -1 if undefined. OCCT hxx L55.
    pub fn deflection(&self) -> f64 {
        self.my_deflection
    }

    /// Returns angular deflection or -1 if undefined. OCCT hxx L58.
    pub fn angle(&self) -> f64 {
        self.my_angle
    }

    /// Returns minimum size or -1 if undefined. OCCT hxx L61.
    pub fn min_size(&self) -> f64 {
        self.my_min_size
    }
}

// =============================================================================
// Poly_Triangulation (OCCT Poly_Triangulation.hxx L62-402 / .cxx L32-523)
// =============================================================================

/// Provides a triangulation for a surface, a set of surfaces, or more
/// generally a shape (OCCT Poly_Triangulation.hxx L38-61).
///
/// A triangulation comprises:
/// - A table of 3D nodes (3D points on the surface).
/// - A table of triangles (each a triplet of node indices).
/// - An optional table of 2D (UV) nodes, parallel to the table of 3D nodes.
/// - An optional table of 3D normals, parallel to the table of 3D nodes.
/// - An optional deflection maximizing the distance from the surface to the
///   triangulation.
///
/// Index convention: all public accessors keep the OCCT 1-based index range
/// ([1, NbNodes()] / [1, NbTriangles()]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolyTriangulation {
    /// OCCT `double myDeflection` (Poly_Triangulation.hxx L394).
    pub my_deflection: f64,
    /// OCCT `Poly_ArrayOfNodes myNodes` (hxx L395). rcad keeps plain f64
    /// storage (no single-precision variant — documented module omission).
    pub my_nodes: Vec<DVec3>,
    /// OCCT `NCollection_Array1<Poly_Triangle> myTriangles` (hxx L396).
    pub my_triangles: Vec<PolyTriangle>,
    /// OCCT `Poly_ArrayOfUVNodes myUVNodes` (hxx L397). Empty = no UV nodes
    /// (OCCT HasUVNodes() tests the array emptiness, hxx L137).
    pub my_uv_nodes: Vec<DVec2>,
    /// OCCT `NCollection_Array1<NCollection_Vec3<float>> myNormals` (hxx L398).
    /// Empty = no normals (hxx L140). Normals stay single precision like OCCT.
    pub my_normals: Vec<[f32; 3]>,
    /// OCCT `Poly_MeshPurpose myPurpose` (hxx L399).
    pub my_purpose: PolyMeshPurpose,
    /// OCCT `occ::handle<Poly_TriangulationParameters> myParams` (hxx L401).
    pub my_params: Option<PolyTriangulationParameters>,
    // OCCT hxx L392-393 (myCachedMinMax / myCachedMinMaxMutex) — omitted: the
    // cached min-max is a mutex-guarded optimization cache without semantic
    // effect; rcad has no Bnd_Box type (see module docs). min_max() below
    // always takes the accurate computeBoundingBox path (cxx L406-409).
}

impl Default for PolyTriangulation {
    fn default() -> Self {
        Self::new()
    }
}

impl PolyTriangulation {
    /// Constructs an empty triangulation. OCCT Poly_Triangulation.cxx L32-36.
    pub fn new() -> Self {
        PolyTriangulation {
            my_deflection: 0.,
            my_nodes: Vec::new(),
            my_triangles: Vec::new(),
            my_uv_nodes: Vec::new(),
            my_normals: Vec::new(),
            my_purpose: POLY_MESH_PURPOSE_NONE,
            my_params: None,
        }
    }

    /// Constructs a triangulation capable of containing the specified number
    /// of nodes and triangles. OCCT Poly_Triangulation.cxx L40-57.
    pub fn with_sizes(
        the_nb_nodes: usize,
        the_nb_triangles: usize,
        the_has_uv_nodes: bool,
        the_has_normals: bool,
    ) -> Self {
        let mut t = PolyTriangulation {
            my_deflection: 0.,
            my_nodes: vec![DVec3::ZERO; the_nb_nodes],
            my_triangles: vec![PolyTriangle::zeroed(); the_nb_triangles],
            my_uv_nodes: Vec::new(),
            my_normals: Vec::new(),
            my_purpose: POLY_MESH_PURPOSE_NONE,
            my_params: None,
        };
        if the_has_uv_nodes {
            // cxx L49-52: myUVNodes.Resize(theNbNodes, false)
            t.my_uv_nodes = vec![DVec2::ZERO; the_nb_nodes];
        }
        if the_has_normals {
            // cxx L53-56: myNormals.Resize(0, theNbNodes - 1, false)
            t.my_normals = vec![[0.0; 3]; the_nb_nodes];
        }
        t
    }

    /// Constructs a triangulation from a set of triangles: 3D points from
    /// Nodes and triangles from Triangles. OCCT Poly_Triangulation.cxx L61-71.
    pub fn from_nodes_triangles(the_nodes: &[DVec3], the_triangles: &[PolyTriangle]) -> Self {
        PolyTriangulation {
            my_deflection: 0.,
            my_nodes: the_nodes.to_vec(),
            my_triangles: the_triangles.to_vec(),
            my_uv_nodes: Vec::new(),
            my_normals: Vec::new(),
            my_purpose: POLY_MESH_PURPOSE_NONE,
            my_params: None,
        }
    }

    /// Constructs a triangulation from nodes, UV nodes and triangles.
    /// OCCT Poly_Triangulation.cxx L75-89.
    pub fn from_nodes_uv_triangles(
        the_nodes: &[DVec3],
        the_uv_nodes: &[DVec2],
        the_triangles: &[PolyTriangle],
    ) -> Self {
        PolyTriangulation {
            my_deflection: 0.,
            my_nodes: the_nodes.to_vec(),
            my_triangles: the_triangles.to_vec(),
            my_uv_nodes: the_uv_nodes.to_vec(),
            my_normals: Vec::new(),
            my_purpose: POLY_MESH_PURPOSE_NONE,
            my_params: None,
        }
    }

    /// Creates full copy of current triangulation. OCCT Poly_Triangulation.cxx
    /// L100-103 / copy constructor L107-121 (deep copy of all arrays and
    /// parameters) — the derived `Clone` impl is the direct equivalent.
    pub fn copy(&self) -> Self {
        self.clone()
    }

    /// Returns the deflection of this triangulation. OCCT hxx L109.
    pub fn deflection(&self) -> f64 {
        self.my_deflection
    }

    /// Sets the deflection of this triangulation. OCCT hxx L113.
    pub fn set_deflection(&mut self, the_deflection: f64) {
        self.my_deflection = the_deflection;
    }

    /// Returns initial set of parameters used to generate this triangulation.
    /// OCCT hxx L116.
    pub fn parameters(&self) -> Option<&PolyTriangulationParameters> {
        self.my_params.as_ref()
    }

    /// Updates initial set of parameters used to generate this triangulation.
    /// OCCT hxx L119-122.
    pub fn set_parameters(&mut self, the_params: Option<PolyTriangulationParameters>) {
        self.my_params = the_params;
    }

    /// Clears internal arrays of nodes and all attributes.
    /// OCCT Poly_Triangulation.cxx L125-140.
    pub fn clear(&mut self) {
        if !self.my_nodes.is_empty() {
            self.my_nodes = Vec::new();
        }
        if !self.my_triangles.is_empty() {
            self.my_triangles = Vec::new();
        }
        self.remove_uv_nodes();
        self.remove_normals();
    }

    /// Returns TRUE if triangulation has some geometry. OCCT hxx L128.
    pub fn has_geometry(&self) -> bool {
        !self.my_nodes.is_empty() && !self.my_triangles.is_empty()
    }

    /// Returns the number of nodes for this triangulation. OCCT hxx L131.
    pub fn nb_nodes(&self) -> usize {
        self.my_nodes.len()
    }

    /// Returns the number of triangles for this triangulation. OCCT hxx L134.
    pub fn nb_triangles(&self) -> usize {
        self.my_triangles.len()
    }

    /// Returns true if 2D nodes are associated with 3D nodes. OCCT hxx L137.
    pub fn has_uv_nodes(&self) -> bool {
        !self.my_uv_nodes.is_empty()
    }

    /// Returns true if nodal normals are defined. OCCT hxx L140.
    pub fn has_normals(&self) -> bool {
        !self.my_normals.is_empty()
    }

    /// Returns a node at the given index (1-based). OCCT hxx L145.
    pub fn node(&self, the_index: usize) -> DVec3 {
        self.my_nodes[the_index - 1]
    }

    /// Sets a node coordinates (1-based index). OCCT hxx L150.
    pub fn set_node(&mut self, the_index: usize, the_pnt: DVec3) {
        self.my_nodes[the_index - 1] = the_pnt;
    }

    /// Returns UV-node at the given index (1-based). OCCT hxx L155.
    pub fn uv_node(&self, the_index: usize) -> DVec2 {
        self.my_uv_nodes[the_index - 1]
    }

    /// Sets an UV-node coordinates (1-based index). OCCT hxx L160.
    pub fn set_uv_node(&mut self, the_index: usize, the_pnt: DVec2) {
        self.my_uv_nodes[the_index - 1] = the_pnt;
    }

    /// Returns triangle at the given index (1-based). OCCT hxx L165.
    pub fn triangle(&self, the_index: usize) -> &PolyTriangle {
        &self.my_triangles[the_index - 1]
    }

    /// Sets a triangle (1-based index). OCCT hxx L171-174.
    pub fn set_triangle(&mut self, the_index: usize, the_triangle: PolyTriangle) {
        self.my_triangles[the_index - 1] = the_triangle;
    }

    /// Returns normal at the given index (1-based), as the OCCT single
    /// precision triple. OCCT hxx L179-191.
    pub fn normal(&self, the_index: usize) -> [f32; 3] {
        self.my_normals[the_index - 1]
    }

    /// Changes normal at the given index (1-based). OCCT hxx L196-199.
    pub fn set_normal(&mut self, the_index: usize, the_normal: [f32; 3]) {
        self.my_normals[the_index - 1] = the_normal;
    }

    /// Changes normal at the given index from a double-precision direction —
    /// the gp_Dir overload converting through float (OCCT hxx L204-209:
    /// `NCollection_Vec3<float>(float(theNormal.X()), ...)`).
    pub fn set_normal_dir(&mut self, the_index: usize, the_normal: DVec3) {
        self.set_normal(the_index, [
            the_normal.x as f32,
            the_normal.y as f32,
            the_normal.z as f32,
        ]);
    }

    /// Returns mesh purpose bits. OCCT hxx L212.
    pub fn mesh_purpose(&self) -> PolyMeshPurpose {
        self.my_purpose
    }

    /// Sets mesh purpose bits. OCCT hxx L215.
    pub fn set_mesh_purpose(&mut self, the_purpose: PolyMeshPurpose) {
        self.my_purpose = the_purpose;
    }

    /// Returns TRUE if node positions are defined with double precision.
    /// OCCT hxx L259. rcad node arrays are always f64 (module note on
    /// Poly_ArrayOfNodes precision switching), so this is constant TRUE —
    /// the same value OCCT returns by default.
    pub fn is_double_precision(&self) -> bool {
        true
    }

    /// Method resizing internal arrays of nodes (synchronously for all
    /// attributes). OCCT Poly_Triangulation.cxx L286-297.
    pub fn resize_nodes(&mut self, the_nb_nodes: usize, the_to_copy_old: bool) {
        let shrink = |v: &mut Vec<DVec2>| {
            if the_to_copy_old {
                v.resize(the_nb_nodes, DVec2::ZERO);
            } else {
                v.clear();
                v.resize(the_nb_nodes, DVec2::ZERO);
            }
        };
        if the_to_copy_old {
            self.my_nodes.resize(the_nb_nodes, DVec3::ZERO);
        } else {
            self.my_nodes.clear();
            self.my_nodes.resize(the_nb_nodes, DVec3::ZERO);
        }
        if !self.my_uv_nodes.is_empty() {
            shrink(&mut self.my_uv_nodes);
        }
        if !self.my_normals.is_empty() {
            if the_to_copy_old {
                self.my_normals.resize(the_nb_nodes, [0.0; 3]);
            } else {
                self.my_normals.clear();
                self.my_normals.resize(the_nb_nodes, [0.0; 3]);
            }
        }
    }

    /// Method resizing an internal array of triangles.
    /// OCCT Poly_Triangulation.cxx L301-304.
    pub fn resize_triangles(&mut self, the_nb_triangles: usize, the_to_copy_old: bool) {
        if the_to_copy_old {
            self.my_triangles.resize(the_nb_triangles, PolyTriangle::zeroed());
        } else {
            self.my_triangles.clear();
            self.my_triangles
                .resize(the_nb_triangles, PolyTriangle::zeroed());
        }
    }

    /// If an array for UV coordinates is not allocated yet, do it now.
    /// OCCT Poly_Triangulation.cxx L308-314.
    pub fn add_uv_nodes(&mut self) {
        if self.my_uv_nodes.is_empty() || self.my_uv_nodes.len() != self.my_nodes.len() {
            self.my_uv_nodes = vec![DVec2::ZERO; self.my_nodes.len()];
        }
    }

    /// Deallocates the UV nodes array. OCCT Poly_Triangulation.cxx L144-152.
    pub fn remove_uv_nodes(&mut self) {
        if !self.my_uv_nodes.is_empty() {
            self.my_uv_nodes = Vec::new();
        }
    }

    /// If an array for normals is not allocated yet, do it now.
    /// OCCT Poly_Triangulation.cxx L318-324.
    pub fn add_normals(&mut self) {
        if self.my_normals.is_empty() || self.my_normals.len() != self.my_nodes.len() {
            self.my_normals = vec![[0.0; 3]; self.my_nodes.len()];
        }
    }

    /// Deallocates the normals array. OCCT Poly_Triangulation.cxx L156-163.
    pub fn remove_normals(&mut self) {
        if !self.my_normals.is_empty() {
            self.my_normals = Vec::new();
        }
    }

    /// Compute smooth normals by averaging triangle normals.
    /// OCCT Poly_Triangulation.cxx L442-476.
    pub fn compute_normals(&mut self) {
        // zero values (cxx L444-446)
        self.add_normals();
        self.my_normals = vec![[0.0; 3]; self.my_normals.len()];

        for a_tri_iter in 1..=self.nb_triangles() {
            let (n0i, n1i, n2i) = self.triangle(a_tri_iter).get();
            let a_node0 = self.node(n0i as usize);
            let a_node1 = self.node(n1i as usize);
            let a_node2 = self.node(n2i as usize);

            // cxx L457-459: aVec01 = n1 - n0; aVec02 = n2 - n0; cross product
            let a_vec01 = a_node1 - a_node0;
            let a_vec02 = a_node2 - a_node0;
            let a_tri_norm = a_vec01.cross(a_vec02);
            let a_norm3f: [f32; 3] = [
                a_tri_norm.x as f32,
                a_tri_norm.y as f32,
                a_tri_norm.z as f32,
            ];
            for &a_node in &[n0i, n1i, n2i] {
                // cxx L464: myNormals.ChangeValue(...) += aNorm3f
                let dst = &mut self.my_normals[(a_node - 1) as usize];
                dst[0] += a_norm3f[0];
                dst[1] += a_norm3f[1];
                dst[2] += a_norm3f[2];
            }
        }

        // Normalize all vectors (cxx L468-475); zero-length normals fall back
        // to +Z like OCCT.
        for a_norm3f in self.my_normals.iter_mut() {
            let a_mod =
                (a_norm3f[0] * a_norm3f[0] + a_norm3f[1] * a_norm3f[1] + a_norm3f[2] * a_norm3f[2])
                    .sqrt();
            if a_mod == 0.0 {
                *a_norm3f = [0.0, 0.0, 1.0];
            } else {
                *a_norm3f = [a_norm3f[0] / a_mod, a_norm3f[1] / a_mod, a_norm3f[2] / a_mod];
            }
        }
    }

    /// Extends the passed box with the bounding box of this triangulation and
    /// returns TRUE when data was added.
    /// OCCT Poly_Triangulation.cxx L388-416 (MinMax) — rcad has no Bnd_Box
    /// cache, so the cached branch (cxx L394-405) is always skipped and the
    /// accurate computeBoundingBox path (cxx L406-409) is taken, exactly the
    /// behavior of an uncached OCCT triangulation. Returns the box as
    /// (min, max); None when the box is void (no nodes).
    pub fn min_max(
        &self,
        the_trsf: glam::DAffine3,
        _the_is_accurate: bool,
    ) -> Option<(DVec3, DVec3)> {
        self.compute_bounding_box(the_trsf)
    }

    /// Calculates the bounding box of nodal data.
    /// OCCT Poly_Triangulation.cxx L420-438 (computeBoundingBox).
    pub fn compute_bounding_box(&self, the_trsf: glam::DAffine3) -> Option<(DVec3, DVec3)> {
        if self.nb_nodes() == 0 {
            return None; // Bnd_Box stays void (cxx L410 IsVoid -> false)
        }
        let identity = the_trsf == glam::DAffine3::IDENTITY;
        let mut a_min = DVec3::INFINITY;
        let mut a_max = DVec3::NEG_INFINITY;
        for a_node_idx in 1..=self.nb_nodes() {
            // cxx L423-437: identity transform adds the raw node, otherwise
            // the transformed point.
            let a_pnt = if identity {
                self.node(a_node_idx)
            } else {
                the_trsf.transform_point3(self.node(a_node_idx))
            };
            a_min = a_min.min(a_pnt);
            a_max = a_max.max(a_pnt);
        }
        Some((a_min, a_max))
    }

    /// Returns number of deferred nodes that can be loaded using
    /// LoadDeferredData(). OCCT hxx L344 (base-class default 0 — rcad has no
    /// deferred/late-load triangulation subclasses).
    pub fn nb_deferred_nodes(&self) -> usize {
        0
    }

    /// Returns number of deferred triangles. OCCT hxx L349 (base default 0).
    pub fn nb_deferred_triangles(&self) -> usize {
        0
    }

    /// Returns TRUE if there is some triangulation data that can be loaded.
    /// OCCT hxx L352 (base: NbDeferredTriangles() > 0).
    pub fn has_deferred_data(&self) -> bool {
        self.nb_deferred_triangles() > 0
    }

    /// Loads triangulation data into itself from deferred storage.
    /// OCCT Poly_Triangulation.cxx L480-492 — with the constant base-class
    /// HasDeferredData()==false this always returns false, byte-for-byte the
    /// OCCT base behavior (the purpose-bit update on success is unreachable).
    pub fn load_deferred_data(&mut self) -> bool {
        if !self.has_deferred_data() {
            return false;
        }
        self.set_mesh_purpose(self.my_purpose | POLY_MESH_PURPOSE_LOADED);
        true
    }

    /// Releases triangulation data if it has connected deferred storage.
    /// OCCT Poly_Triangulation.cxx L514-523 — always false for the same
    /// reason as load_deferred_data.
    pub fn unload_deferred_data(&mut self) -> bool {
        if self.has_deferred_data() {
            self.clear();
            self.set_mesh_purpose(self.my_purpose & !POLY_MESH_PURPOSE_LOADED);
            return true;
        }
        false
    }
}

// =============================================================================
// Poly_Polygon3D (OCCT Poly_Polygon3D.hxx L30-94 / .cxx L23-88)
// =============================================================================

/// This class provides a polygon in 3D space, generally an approximate
/// representation of a curve. If the polygon is closed, the point of closure
/// is repeated at the end of the table of nodes. Each node of an approximate
/// representation of a curve can carry the parameter of the corresponding
/// point on the curve. OCCT Poly_Polygon3D.hxx L25-29.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolyPolygon3D {
    /// OCCT `double myDeflection` (hxx L91).
    pub my_deflection: f64,
    /// OCCT `NCollection_Array1<gp_Pnt> myNodes` (hxx L92).
    pub my_nodes: Vec<DVec3>,
    /// OCCT `occ::handle<NCollection_HArray1<double>> myParameters` (hxx L93).
    /// None = null handle (HasParameters() == false).
    pub my_parameters: Option<Vec<f64>>,
}

impl PolyPolygon3D {
    /// Constructs a 3D polygon with specific number of nodes.
    /// OCCT Poly_Polygon3D.cxx L23-31.
    pub fn new(the_nb_nodes: usize, the_has_params: bool) -> Self {
        PolyPolygon3D {
            my_deflection: 0.0,
            my_nodes: vec![DVec3::ZERO; the_nb_nodes],
            my_parameters: if the_has_params {
                Some(vec![0.0; the_nb_nodes])
            } else {
                None
            },
        }
    }

    /// Constructs a 3D polygon defined by the table of points, Nodes.
    /// OCCT Poly_Polygon3D.cxx L35-44.
    pub fn from_nodes(the_nodes: &[DVec3]) -> Self {
        PolyPolygon3D {
            my_deflection: 0.,
            my_nodes: the_nodes.to_vec(),
            my_parameters: None,
        }
    }

    /// Constructs a 3D polygon defined by the table of points, Nodes, and the
    /// parallel table of parameters. OCCT Poly_Polygon3D.cxx L48-62.
    pub fn from_nodes_params(the_nodes: &[DVec3], the_parameters: &[f64]) -> Self {
        PolyPolygon3D {
            my_deflection: 0.,
            my_nodes: the_nodes.to_vec(),
            my_parameters: Some(the_parameters.to_vec()),
        }
    }

    /// Returns the deflection of this polygon. OCCT hxx L54.
    pub fn deflection(&self) -> f64 {
        self.my_deflection
    }

    /// Sets the deflection of this polygon. OCCT hxx L57.
    pub fn set_deflection(&mut self, the_defl: f64) {
        self.my_deflection = the_defl;
    }

    /// Returns the number of nodes in this polygon (on a closed triangle the
    /// point of closure is repeated, so NbNodes returns 4). OCCT hxx L63.
    pub fn nb_nodes(&self) -> usize {
        self.my_nodes.len()
    }

    /// Returns the table of nodes for this polygon. OCCT hxx L66.
    pub fn nodes(&self) -> &[DVec3] {
        &self.my_nodes
    }

    /// Returns the mutable table of nodes (OCCT ChangeNodes, hxx L69).
    pub fn change_nodes(&mut self) -> &mut Vec<DVec3> {
        &mut self.my_nodes
    }

    /// Returns true if parameters are associated with the nodes.
    /// OCCT hxx L73.
    pub fn has_parameters(&self) -> bool {
        self.my_parameters.is_some()
    }

    /// Returns the table of the parameters associated with each node.
    /// OCCT hxx L77 (null-handle dereference -> panic, same as OCCT UB).
    pub fn parameters(&self) -> &[f64] {
        self.my_parameters
            .as_deref()
            .expect("Poly_Polygon3D::Parameters : parameters is NULL")
    }

    /// Returns the mutable table of the parameters (OCCT ChangeParameters,
    /// hxx L83).
    pub fn change_parameters(&mut self) -> &mut Vec<f64> {
        self.my_parameters
            .as_mut()
            .expect("Poly_Polygon3D::ChangeParameters : parameters is NULL")
    }
}

// =============================================================================
// Poly_PolygonOnTriangulation
// (OCCT Poly_PolygonOnTriangulation.hxx L39-144 / .cxx L24-95)
// =============================================================================

/// This class provides a polygon in 3D space, based on the triangulation of a
/// surface. Each node is an index in the table of nodes specific to a
/// triangulation. If the polygon is closed, the index of the point of closure
/// is repeated at the end of the table of nodes. Each node of an approximate
/// representation of a curve can carry the parameter of the corresponding
/// point on the curve. OCCT Poly_PolygonOnTriangulation.hxx L27-38.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolyPolygonOnTriangulation {
    /// OCCT `double myDeflection` (hxx L141).
    pub my_deflection: f64,
    /// OCCT `NCollection_Array1<int> myNodes` (hxx L142).
    pub my_nodes: Vec<i32>,
    /// OCCT `occ::handle<NCollection_HArray1<double>> myParameters` (hxx L143).
    pub my_parameters: Option<Vec<f64>>,
}

impl PolyPolygonOnTriangulation {
    /// Constructs a 3D polygon on the triangulation of a shape with specified
    /// size of nodes. OCCT Poly_PolygonOnTriangulation.cxx L24-33.
    pub fn new(the_nb_nodes: usize, the_has_params: bool) -> Self {
        PolyPolygonOnTriangulation {
            my_deflection: 0.0,
            my_nodes: vec![0; the_nb_nodes],
            my_parameters: if the_has_params {
                Some(vec![0.0; the_nb_nodes])
            } else {
                None
            },
        }
    }

    /// Constructs a 3D polygon on the triangulation of a shape, defined by
    /// the table of nodes. OCCT Poly_PolygonOnTriangulation.cxx L37-42.
    pub fn from_nodes(the_nodes: &[i32]) -> Self {
        PolyPolygonOnTriangulation {
            my_deflection: 0.0,
            my_nodes: the_nodes.to_vec(),
            my_parameters: None,
        }
    }

    /// Constructs a 3D polygon on the triangulation of a shape, defined by
    /// the table of nodes and the table of parameters.
    /// OCCT Poly_PolygonOnTriangulation.cxx L46-55.
    pub fn from_nodes_params(the_nodes: &[i32], the_parameters: &[f64]) -> Self {
        PolyPolygonOnTriangulation {
            my_deflection: 0.0,
            my_nodes: the_nodes.to_vec(),
            my_parameters: Some(the_parameters.to_vec()),
        }
    }

    /// Returns the deflection of this polygon. OCCT hxx L68.
    pub fn deflection(&self) -> f64 {
        self.my_deflection
    }

    /// Sets the deflection of this polygon. OCCT hxx L72.
    pub fn set_deflection(&mut self, the_defl: f64) {
        self.my_deflection = the_defl;
    }

    /// Returns the number of nodes for this polygon. OCCT hxx L78.
    pub fn nb_nodes(&self) -> usize {
        self.my_nodes.len()
    }

    /// Returns node at the given index (1-based). OCCT hxx L81.
    pub fn node(&self, the_index: usize) -> i32 {
        self.my_nodes[the_index - 1]
    }

    /// Returns mutable node-index array (OCCT ChangeNodeArray, hxx L84).
    pub fn change_node_array(&mut self) -> &mut Vec<i32> {
        &mut self.my_nodes
    }

    /// Sets node at the given index (1-based). OCCT hxx L87.
    pub fn set_node(&mut self, the_index: usize, the_node: i32) {
        self.my_nodes[the_index - 1] = the_node;
    }

    /// Returns true if parameters are associated with the nodes.
    /// OCCT hxx L90.
    pub fn has_parameters(&self) -> bool {
        self.my_parameters.is_some()
    }

    /// Returns parameter at the given index (1-based).
    /// OCCT hxx L93-98 (Standard_NullObject raise on null handle).
    pub fn parameter(&self, the_index: usize) -> f64 {
        self.my_parameters.as_deref().expect(
            "Poly_PolygonOnTriangulation::Parameter : parameters is NULL",
        )[the_index - 1]
    }

    /// Sets parameter at the given index (1-based).
    /// OCCT hxx L101-106 (Standard_NullObject raise on null handle).
    pub fn set_parameter(&mut self, the_index: usize, the_value: f64) {
        self.my_parameters
            .as_mut()
            .expect("Poly_PolygonOnTriangulation::Parameter : parameters is NULL")[the_index - 1] =
            the_value;
    }

    /// Sets the table of the parameters associated with each node. Raises
    /// exception if array size doesn't match number of polygon nodes.
    /// OCCT Poly_PolygonOnTriangulation.cxx L74-83.
    pub fn set_parameters(&mut self, the_parameters: Option<Vec<f64>>) {
        if let Some(ref params) = the_parameters {
            assert!(
                params.len() == self.my_nodes.len(),
                "Poly_PolygonOnTriangulation::SetParameters() - invalid array size"
            );
        }
        self.my_parameters = the_parameters;
    }
}
