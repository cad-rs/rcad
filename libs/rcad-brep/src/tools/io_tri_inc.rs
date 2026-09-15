// =============================================================================
// BRepTools triangulation persistence (BRepTools.cxx Write/Clean +
// BRepTools_ShapeSet.cxx v2/v3 triangulation record set)
// =============================================================================
// OCCT writes the Poly data through BRepTools_ShapeSet in three geometry
// sections (WriteGeometry order, BRepTools_ShapeSet.cxx L257-294):
//   1. "Polygon3D"                    (WritePolygon3D,           L1419-1514)
//   2. "PolygonOnTriangulations"      (WritePolygonOnTriangulation, L1264-1349)
//   3. "Triangulations"               (WriteTriangulation,       L1575-1767)
// and binds them from the topology records:
//   - face record tail  "2 <tri-index>"   (WriteGeometry(face),  L750-785)
//   - edge record kinds "6 pt t l" (PolygonOnTriangulation) and
//     "7 pt pt2 t l" (PolygonOnClosedTriangulation)      (L709-743)
//
// rcad persistence is JSON, but the RECORD SET below is 1:1 with the OCCT
// v2/v3 triangulation section: every field of every record carries the OCCT
// anchor that produced it. The `shape_set` blob is the rcad superset (it
// additionally keeps Poly_MeshPurpose / Poly_TriangulationParameters, which
// the OCCT brep format never persists — WriteTriangulation has no fields for
// them). The read path rebuilds from `shape_set` and cross-checks the
// sections, so a file whose sections disagree with the topology is rejected.

use rcad_kernel::topo::poly::{
    PolyPolygon3D, PolyPolygonOnTriangulation, PolyTriangulation,
};
use rcad_kernel::topods::{CurveRepresentation, TShape};

/// Defined TopTools format version. OCCT TopTools_FormatVersion.hxx L17-31:
/// V1 = base, V2 = stores CurveOnSurface UV points, V3 = stores per-vertex
/// normals for triangulation-only faces; CURRENT = V3.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepFormatVersion {
    /// TopTools_FormatVersion_VERSION_1.
    Version1,
    /// TopTools_FormatVersion_VERSION_2.
    Version2,
    /// TopTools_FormatVersion_VERSION_3 (TopTools_FormatVersion_CURRENT).
    Version3,
}

impl BrepFormatVersion {
    /// The numeric value written into the file header (OCCT FormatNb()).
    pub fn nb(self) -> u32 {
        match self {
            BrepFormatVersion::Version1 => 1,
            BrepFormatVersion::Version2 => 2,
            BrepFormatVersion::Version3 => 3,
        }
    }

    /// TopTools_FormatVersion_CURRENT.
    pub fn current() -> Self {
        BrepFormatVersion::Version3
    }
}

/// One "Polygon3D" section record. OCCT WritePolygon3D compact record
/// (BRepTools_ShapeSet.cxx L1442-1512): `<nbNodes> <hasParameters>`,
/// `<deflection>`, the node triples, and (when HasParameters) the parameter
/// list. Read back by ReadPolygon3D (L1525-1571).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Polygon3DRecord {
    /// OCCT `P->NbNodes()` (L1444).
    pub nb_nodes: usize,
    /// OCCT `((P->HasParameters()) ? "1" : "0")` (L1445).
    pub has_parameters: u8,
    /// OCCT `P->Deflection()` (L1458).
    pub deflection: f64,
    /// OCCT `Nodes(j).X()/Y()/Z()` triples (L1468-1497).
    pub nodes: Vec<[f64; 3]>,
    /// OCCT `Param(i1)` list (L1506-1511); present iff has_parameters == 1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parameters: Option<Vec<f64>>,
}

/// One "PolygonOnTriangulations" section record. OCCT
/// WritePolygonOnTriangulation compact record (BRepTools_ShapeSet.cxx
/// L1284-1348): `<Nodes.Length()>`, the node index list, `"p"`, the
/// deflection, then `"1"` + parameter list or `"0"`. Read back by
/// ReadPolygonOnTriangulation (L1360-1415).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct PolygonOnTriangulationRecord {
    /// OCCT `Nodes.Length()` (L1294).
    pub nb_nodes: usize,
    /// OCCT `Nodes.Value(j)` list (L1300-1304) — indices into the owning
    /// triangulation's node table.
    pub nodes: Vec<i32>,
    /// OCCT `Poly->Deflection()` (L1318).
    pub deflection: f64,
    /// OCCT `"1 "` / `"0 \n"` after `"p "` (L1332/1346).
    pub has_parameters: u8,
    /// OCCT `Param->Value(j)` list (L1338-1342); present iff has_parameters.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parameters: Option<Vec<f64>>,
}

/// One "Triangulations" section record. OCCT WriteTriangulation compact
/// record (BRepTools_ShapeSet.cxx L1596-1767):
/// `<nbNodes> <nbTriangles> <hasUVNodes> [<hasNormals> when v3]`,
/// `<deflection>`, the 3D nodes, the UV nodes (when hasUVNodes), the
/// triangles, and the normals (when v3 && HasNormals && toWriteNormals).
/// Read back by ReadTriangulation (L1779-1853).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TriangulationRecord {
    /// OCCT `T->NbNodes()` (L1603).
    pub nb_nodes: usize,
    /// OCCT `T->NbTriangles()` (L1603).
    pub nb_triangles: usize,
    /// OCCT `((T->HasUVNodes()) ? "1" : "0")` (L1604).
    pub has_uv_nodes: u8,
    /// OCCT `((T->HasNormals() && toWriteNormals) ? "1" : "0")` — v3 only
    /// (L1607); v1/v2 files carry no normals. Note the OCCT asymmetry kept
    /// verbatim: written when HasNormals && toWriteNormals, but on read
    /// (L1838-1848) the flag directly drives SetNormal — the writer's flag is
    /// the single source of truth in the file.
    #[serde(default)]
    pub has_normals: u8,
    /// OCCT `T->Deflection()` (L1628).
    pub deflection: f64,
    /// OCCT `T->Node(j)` triples (L1638-1668).
    pub nodes: Vec<[f64; 3]>,
    /// OCCT `T->UVNode(j)` pairs (L1670-1702); present iff has_uv_nodes == 1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub uv_nodes: Option<Vec<[f64; 2]>>,
    /// OCCT `T->Triangle(j).Get(n1, n2, n3)` triples (L1708-1739).
    pub triangles: Vec<[i32; 3]>,
    /// OCCT `T->Normal(j)` triples — SINGLE precision like OCCT
    /// NCollection_Vec3<float> (L1741-1764); present iff has_normals == 1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normals: Option<Vec<[f32; 3]>>,
}

/// The topology-to-geometry binding records (the OCCT shape-record tails).
/// OCCT face record (BRepTools_ShapeSet.cxx L775-784): after the surface
/// fields the writer appends `2 <triangulation-index>` for the ACTIVE face
/// triangulation when `myWithTriangles || Surface is NULL`. OCCT edge record
/// (L709-743): representation kinds `6` (PolygonOnTriangulation) and `7`
/// (PolygonOnClosedTriangulation, the seam double-string) followed by
/// `<polygon> [<polygon2>] <triangulation> <location>`.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ShapeBindingRecords {
    /// Face records carrying `2 <tri-index>` (L779-783); `triangulation` is
    /// the 1-based section index (0 = the face record has no `2` tail).
    pub faces: Vec<FaceTriangulationBinding>,
    /// Edge records of kind 6/7 (L722-742).
    pub edges: Vec<EdgePolygonBinding>,
}

/// OCCT face-record `2 <tri-index>` tail (BRepTools_ShapeSet.cxx L775-784).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FaceTriangulationBinding {
    /// The face TShape pool index (the rcad stand-in for the TopoDS_Face
    /// occurrence being written).
    pub face: usize,
    /// 1-based index into the `triangulations` section; 0 = no tail written
    /// (the `myWithTriangles || Surface is NULL` + null-triangulation guard,
    /// L775-777).
    pub triangulation: usize,
}

/// OCCT edge-record kind 6/7 (BRepTools_ShapeSet.cxx L722-742).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EdgePolygonBinding {
    /// The edge TShape pool index.
    pub edge: usize,
    /// OCCT kind: 6 = Polygon on triangulation, 7 = Polygon on CLOSED
    /// triangulation (L726-733).
    pub kind: u8,
    /// OCCT `myNodes.FindIndex(PT->PolygonOnTriangulation())` (L734).
    pub polygon: usize,
    /// OCCT `myNodes.FindIndex(PT->PolygonOnTriangulation2())` (L737) —
    /// present only for kind 7 (the seam double string).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub polygon2: Option<usize>,
    /// OCCT `myTriangulations.FindIndex(PT->Triangulation())` (L739),
    /// 1-based section index.
    pub triangulation: usize,
    /// OCCT `Locations().Index(CR->Location())` (L740).
    pub location: u32,
}

/// The full OCCT v2/v3 mesh record set of one shape set, in the
/// WriteGeometry emission order (BRepTools_ShapeSet.cxx L257-294:
/// Polygon3D, PolygonOnTriangulations, ..., Triangulations).
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct BrepMeshSections {
    /// The OCCT "Polygon3D <count>" section records.
    pub polygon3d: Vec<Polygon3DRecord>,
    /// The OCCT "PolygonOnTriangulations <count>" section records.
    pub polygon_on_triangulations: Vec<PolygonOnTriangulationRecord>,
    /// The OCCT "Triangulations <count>" section records.
    pub triangulations: Vec<TriangulationRecord>,
    /// The face/edge topology record tails binding the sections.
    pub references: ShapeBindingRecords,
}

/// OCCT BRepTools_ShapeSet geometry collection state (BRepTools_ShapeSet.hxx):
/// `myPolygons3D` / `myNodes` (handle-indexed maps) and
/// `myTriangulations` (handle -> toWriteNormals flag map). The rcad maps key
/// by pool index / value equality — the stand-in for handle identity.
#[derive(Debug, Default)]
struct ShapeSetMeshState {
    /// OCCT myPolygons3D (TColStd_IndexedMapOfTransient).
    my_polygons3d: Vec<PolyPolygon3D>,
    /// OCCT myNodes (TColStd_IndexedMapOfTransient).
    my_nodes: Vec<PolyPolygonOnTriangulation>,
    /// OCCT myTriangulations — the flag is `toWriteNormals`
    /// (Add(Tr, needNormals), BRepTools_ShapeSet.cxx L235/L192).
    my_triangulations: Vec<(usize, bool)>,
}

impl ShapeSetMeshState {
    /// NCollection IndexedMap::Add — appends when absent, returns the
    /// 1-based index either way (first insertion wins).
    fn add_polygon3d(&mut self, p: &PolyPolygon3D) -> usize {
        if let Some(i) = self.my_polygons3d.iter().position(|x| x == p) {
            return i + 1;
        }
        self.my_polygons3d.push(p.clone());
        self.my_polygons3d.len()
    }

    fn add_node(&mut self, p: &PolyPolygonOnTriangulation) -> usize {
        if let Some(i) = self.my_nodes.iter().position(|x| x == p) {
            return i + 1;
        }
        self.my_nodes.push(p.clone());
        self.my_nodes.len()
    }

    fn add_triangulation(&mut self, tri: usize, need_normals: bool) -> usize {
        if let Some(i) = self.my_triangulations.iter().position(|x| x.0 == tri) {
            return i + 1;
        }
        self.my_triangulations.push((tri, need_normals));
        self.my_triangulations.len()
    }
}

/// Walk a brep and collect the mesh sections exactly as
/// BRepTools_ShapeSet::AddGeometry does for the Poly data
/// (BRepTools_ShapeSet.cxx L106-241). OCCT processes shapes from complex to
/// elementary, so the TopAbs_FACE branch runs before the TopAbs_EDGE branch
/// (L187-190) — the rcad walk visits all face TShapes first (pool order),
/// then all edge TShapes (pool order), preserving the section index order.
fn collect_mesh_sections(
    brep: &rcad_kernel::topods::BRep,
    the_with_triangles: bool,
    the_with_normals: bool,
    the_version: BrepFormatVersion,
) -> BrepMeshSections {
    let mut st = ShapeSetMeshState::default();
    let mut out = BrepMeshSections::default();

    // ---- TopAbs_FACE branch (L216-240) ----
    for (fi, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Face(tf) = &**ts else { continue };
        // L220: bool needNormals(myWithNormals);
        let mut need_normals = the_with_normals;
        // L222-229: if the surface is null, normals are always required
        // ("triangulation-only Faces ... no analytical geometry to restore
        // normals", TopTools_FormatVersion.hxx L27-29).
        if tf.surface.is_none() {
            need_normals = true;
        }
        // L230-237: if (myWithTriangles || TF->Surface().IsNull())
        if the_with_triangles || tf.surface.is_none() {
            if let Some(tri) = tf.active_triangulation {
                let handle = tf.triangulations[tri];
                let section_index = st.add_triangulation(handle, need_normals);
                // L775-784: the face record `2 <tri-index>` tail
                out.references.faces.push(FaceTriangulationBinding {
                    face: fi,
                    triangulation: section_index,
                });
            }
        }
    }

    // ---- TopAbs_EDGE branch (L140-214 + record writer L709-743) ----
    for (ei, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Edge(te) = &**ts else { continue };
        for cr in &te.representations {
            if cr.is_polygon3d() {
                // L177-184 / L711-721: only when myWithTriangles
                if !the_with_triangles {
                    continue;
                }
                let CurveRepresentation::Polygon3D { polygon, location } = cr else {
                    continue;
                };
                let pt = st.add_polygon3d(polygon);
                out.references.edges.push(EdgePolygonBinding {
                    edge: ei,
                    kind: 5, // OCCT -5- Polygon3D edge record kind (L716)
                    polygon: pt,
                    polygon2: None,
                    triangulation: 0,
                    location: *location,
                });
            } else if cr.is_polygon_on_triangulation() {
                // L185-200 / L722-742: only when myWithTriangles
                if !the_with_triangles {
                    continue;
                }
                let (polygon, polygon2, tri_handle, location) = match cr {
                    CurveRepresentation::PolygonOnTriangulation {
                        polygon,
                        triangulation,
                        location,
                    } => (polygon, None, *triangulation, *location),
                    CurveRepresentation::PolygonOnClosedTriangulation {
                        polygon,
                        polygon2,
                        triangulation,
                        location,
                    } => (polygon, Some(polygon2), *triangulation, *location),
                    _ => continue,
                };
                // L192: edge triangulation does not need normals
                let t_idx = st.add_triangulation(tri_handle, false);
                let pt = st.add_node(polygon);
                let pt2 = polygon2.map(|p| st.add_node(p));
                out.references.edges.push(EdgePolygonBinding {
                    edge: ei,
                    kind: if pt2.is_some() { 7 } else { 6 },
                    polygon: pt,
                    polygon2: pt2,
                    triangulation: t_idx,
                    location,
                });
            }
        }
    }

    // ---- Emit the section records in OCCT compact-record field order ----
    for p in &st.my_polygons3d {
        out.polygon3d.push(Polygon3DRecord {
            nb_nodes: p.nb_nodes(),
            has_parameters: if p.has_parameters() { 1 } else { 0 },
            deflection: p.deflection(),
            nodes: p.nodes().iter().map(|n| [n.x, n.y, n.z]).collect(),
            parameters: p
                .has_parameters()
                .then(|| p.parameters().to_vec()),
        });
    }
    for p in &st.my_nodes {
        out.polygon_on_triangulations.push(PolygonOnTriangulationRecord {
            nb_nodes: p.nb_nodes(),
            nodes: p.my_nodes.clone(),
            deflection: p.deflection(),
            has_parameters: if p.has_parameters() { 1 } else { 0 },
            parameters: p
                .has_parameters()
                .then(|| p.my_parameters.clone().unwrap_or_default()),
        });
    }
    for (handle, to_write_normals) in &st.my_triangulations {
        let t: &PolyTriangulation = &brep.triangulations[*handle];
        // L1601: const bool toWriteNormals = myTriangulations(i);
        // L1605-1608: the hasNormals flag exists from v3 only — a v2 record
        // never carries normals.
        let write_normals = the_version.nb() >= BrepFormatVersion::Version3.nb()
            && t.has_normals()
            && *to_write_normals;
        out.triangulations.push(TriangulationRecord {
            nb_nodes: t.nb_nodes(),
            nb_triangles: t.nb_triangles(),
            has_uv_nodes: if t.has_uv_nodes() { 1 } else { 0 },
            has_normals: if write_normals { 1 } else { 0 },
            deflection: t.deflection(),
            nodes: (1..=t.nb_nodes())
                .map(|j| {
                    let n = t.node(j);
                    [n.x, n.y, n.z]
                })
                .collect(),
            uv_nodes: if t.has_uv_nodes() {
                Some(
                    (1..=t.nb_nodes())
                        .map(|j| {
                            let uv = t.uv_node(j);
                            [uv.x, uv.y]
                        })
                        .collect(),
                )
            } else {
                None
            },
            triangles: (1..=t.nb_triangles())
                .map(|j| {
                    let (n1, n2, n3) = t.triangle(j).get();
                    [n1, n2, n3]
                })
                .collect(),
            normals: if write_normals {
                Some((1..=t.nb_nodes()).map(|j| t.normal(j)).collect())
            } else {
                None
            },
        });
    }
    out
}

/// The persisted document — OCCT `BRepTools::Write(theShape, stream,
/// withTriangles, withNormals, version)` (BRepTools.cxx L619-631) as a JSON
/// record set. Field anchors: `brep_format_version` = SetFormatNb
/// (TopTools_FormatVersion), `with_triangles`/`with_normals` = the
/// BRepTools_ShapeSet constructor flags (BRepTools_ShapeSet.cxx L69-73),
/// `shape_set` = the rcad superset shape dump (carries
/// MeshPurpose/TriangulationParameters, which the OCCT text format omits),
/// `sections` = the 1:1 OCCT v2/v3 mesh record set.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BrepWithTrianglesDocument {
    /// OCCT TopTools_FormatVersion (TopTools_FormatVersion.hxx L17-31).
    pub brep_format_version: u32,
    /// OCCT BRepTools_ShapeSet::myWithTriangles.
    pub with_triangles: bool,
    /// OCCT BRepTools_ShapeSet::myWithNormals.
    pub with_normals: bool,
    /// The shape set (rcad topods::BRep, lossless superset of the OCCT text
    /// format).
    pub shape_set: rcad_kernel::topods::BRep,
    /// The OCCT v2/v3 triangulation record set.
    pub sections: BrepMeshSections,
}

/// OCCT BRepTools::Write(theShape, theStream, theWithTriangles,
/// theWithNormals, theVersion) (BRepTools.cxx L619-631): builds the shape
/// set, collects the geometry and serializes. The returned string is the
/// JSON document; the OCCT "DBRep_DrawableShape\n" header line (L672) is a
/// DRAW-reader convenience and has no record content.
pub fn write_brep_with_triangles_to_string(
    brep: &rcad_kernel::topods::BRep,
    the_with_triangles: bool,
    the_with_normals: bool,
    the_version: BrepFormatVersion,
) -> Result<String, BRepToolsError> {
    let sections = collect_mesh_sections(brep, the_with_triangles, the_with_normals, the_version);
    let doc = BrepWithTrianglesDocument {
        brep_format_version: the_version.nb(),
        with_triangles: the_with_triangles,
        with_normals: the_with_normals,
        shape_set: brep.clone(),
        sections,
    };
    serde_json::to_string_pretty(&doc)
        .map_err(|e| BRepToolsError::SerializationError(e.to_string()))
}

/// File overload — OCCT BRepTools::Write(theShape, theFile, ...) open-stream
/// wrapper (BRepTools.cxx L647-687).
pub fn write_brep_with_triangles_to_file<P: AsRef<Path>>(
    brep: &rcad_kernel::topods::BRep,
    path: P,
    the_with_triangles: bool,
    the_with_normals: bool,
    the_version: BrepFormatVersion,
) -> Result<(), BRepToolsError> {
    let s = write_brep_with_triangles_to_string(brep, the_with_triangles, the_with_normals, the_version)?;
    let mut file = File::create(path.as_ref())
        .map_err(|e| BRepToolsError::IoError(e.to_string()))?;
    file.write_all(s.as_bytes())
        .map_err(|e| BRepToolsError::IoError(e.to_string()))
}

/// OCCT BRepTools::Read(Sh, S, B) (BRepTools.cxx L635-643): parse the
/// document, rebuild the shape set and cross-check the mesh sections against
/// the topology mounts (the section guards of ReadPolygon3D /
/// ReadPolygonOnTriangulation / ReadTriangulation L1111-1129 validate index
/// ranges; here the sections must also agree with the shape set they
/// accompany, otherwise the file is inconsistent).
pub fn read_brep_with_triangles_from_string(
    s: &str,
) -> Result<BrepWithTrianglesDocument, BRepToolsError> {
    let doc: BrepWithTrianglesDocument = serde_json::from_str(s)
        .map_err(|e| BRepToolsError::DeserializationError(e.to_string()))?;
    verify_mesh_sections(&doc)?;
    Ok(doc)
}

/// File overload — OCCT BRepTools::Read(Sh, File, B) (BRepTools.cxx L691-710).
pub fn read_brep_with_triangles_from_file<P: AsRef<Path>>(
    path: P,
) -> Result<BrepWithTrianglesDocument, BRepToolsError> {
    let mut s = String::new();
    File::open(path.as_ref())
        .and_then(|mut f| f.read_to_string(&mut s))
        .map_err(|e| BRepToolsError::IoError(e.to_string()))?;
    read_brep_with_triangles_from_string(&s)
}

/// Cross-check the mesh sections against the shape set: every section record
/// must reproduce the mounted Poly values (count guards first, like the
/// OCCT `t > 0 && t <= myTriangulations.Extent()` index checks at
/// BRepTools_ShapeSet.cxx L1111-1129, then the values).
fn verify_mesh_sections(doc: &BrepWithTrianglesDocument) -> Result<(), BRepToolsError> {
    let version = match doc.brep_format_version {
        1 => BrepFormatVersion::Version1,
        2 => BrepFormatVersion::Version2,
        3 => BrepFormatVersion::Version3,
        v => {
            return Err(BRepToolsError::DeserializationError(format!(
                "unknown brep format version {}",
                v
            )))
        }
    };
    let recomputed = collect_mesh_sections(
        &doc.shape_set,
        doc.with_triangles,
        doc.with_normals,
        version,
    );
    let mism = |what: &str| BRepToolsError::DeserializationError(format!(
        "brep mesh sections inconsistent with shape set: {}",
        what
    ));
    if doc.sections.polygon3d.len() != recomputed.polygon3d.len() {
        return Err(mism("Polygon3D count"));
    }
    if doc.sections.polygon_on_triangulations.len() != recomputed.polygon_on_triangulations.len()
    {
        return Err(mism("PolygonOnTriangulations count"));
    }
    if doc.sections.triangulations.len() != recomputed.triangulations.len() {
        return Err(mism("Triangulations count"));
    }
    if doc.sections.references.faces.len() != recomputed.references.faces.len()
        || doc.sections.references.edges.len() != recomputed.references.edges.len()
    {
        return Err(mism("reference record count"));
    }
    for (a, b) in doc
        .sections
        .polygon3d
        .iter()
        .zip(recomputed.polygon3d.iter())
    {
        if a != b {
            return Err(mism("Polygon3D record"));
        }
    }
    for (a, b) in doc
        .sections
        .polygon_on_triangulations
        .iter()
        .zip(recomputed.polygon_on_triangulations.iter())
    {
        if a != b {
            return Err(mism("PolygonOnTriangulations record"));
        }
    }
    for (a, b) in doc
        .sections
        .triangulations
        .iter()
        .zip(recomputed.triangulations.iter())
    {
        if a != b {
            return Err(mism("Triangulations record"));
        }
    }
    Ok(())
}

/// OCCT BRepTools::Clean(theShape, theForce) (BRepTools.cxx L714-817):
/// removes the triangulations from all geometric faces and the polygon
/// representations from their edges. With `theForce` the
/// PolygonOnTriangulation representations of ALL edges are removed even
/// without a mounted triangulation (L796-812).
///
/// The OCCT exploration (TopExp_Explorer over theShape + the location-less
/// shape map `aShapeMap`, L725-737) maps to a wrapper walk over the pool
/// (solid -> shell -> face -> wire -> edge) deduplicated by TShape index —
/// the location-less map keys by TShape identity. First occurrence wins, so
/// the wrappers keep their locations exactly like the OCCT explorer.
pub fn brep_tools_clean(brep: &mut rcad_kernel::topods::BRep, the_force: bool) {
    use std::collections::HashSet;
    use rcad_kernel::topo::topo_shape::Shape;

    // OCCT L721-723: the null handles used to nullify.
    let a_null_triangulation: Option<usize> = None;
    let a_null_poly: Option<PolyPolygonOnTriangulation> = None;

    // Collect the (face wrapper, edge wrappers) occurrences: solids first,
    // then free shells, then free faces (the same coverage a full
    // TopExp_Explorer walk of the document has). Dedup faces by TShape index
    // (L733: aShapeMap.Add(aFaceNoLoc) — second occurrences skipped).
    let mut seen_faces: HashSet<usize> = HashSet::new();
    let mut face_occurrences: Vec<(Shape, Vec<Shape>)> = Vec::new();
    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };
        for shell in &sd.shells {
            if let TShape::Shell(shd) = &*brep.tshapes[shell.index] {
                for face in &shd.faces {
                    if seen_faces.insert(face.index) {
                        face_occurrences
                            .push((face.clone(), face_wire_edges(brep, face)));
                    }
                }
            }
        }
    }
    for ts in &brep.tshapes {
        let TShape::Shell(shd) = &**ts else { continue };
        for face in &shd.faces {
            if seen_faces.insert(face.index) {
                face_occurrences.push((face.clone(), face_wire_edges(brep, face)));
            }
        }
    }
    for (fi, ts) in brep.tshapes.iter().enumerate() {
        if matches!(&**ts, TShape::Face(_)) && seen_faces.insert(fi) {
            let face = Shape::from_parts(ts.clone(), fi, 0, rcad_kernel::topods::Orientation::Forward);
            face_occurrences.push((face.clone(), face_wire_edges(brep, &face)));
        }
    }

    // L728-767: face loop
    for (a_face, edges) in &face_occurrences {
        // L740-744: BRep_Tool::IsGeometric(aFace) — no surface means no way
        // to recompute the triangulation; keep it.
        let is_geometric = match &*brep.tshapes[a_face.index] {
            TShape::Face(tf) => tf.surface.is_some(),
            _ => false,
        };
        if !is_geometric {
            continue;
        }

        // L746-747: BRep_Tool::Triangulation(aFace, aLoc) — the active
        // triangulation pool handle; aLoc receives the FACE location
        // (BRep_Tool.cxx L114) and flows into the UpdateEdge Predivided
        // composition.
        let a_triangulation = match &*brep.tshapes[a_face.index] {
            TShape::Face(tf) => tf.active_triangulation.map(|i| tf.triangulations[i]),
            _ => continue,
        };
        let Some(a_triangulation) = a_triangulation else {
            // L749-752: continue when the face carries no triangulation.
            continue;
        };

        // L754-764: nullify the PolygonOnTriangulation of the face's edges.
        for an_edge in edges {
            brep.update_edge_polygon_on_triangulation(
                an_edge,
                a_null_poly.clone(),
                a_triangulation,
                a_face.location,
            );
        }

        // L766: aBuilder.UpdateFace(aFace, aNullTriangulation) — reset the
        // face triangulation list.
        brep.update_face_triangulation(a_face, a_null_triangulation, true);
    }

    // L769-816: edge loop — 3d polygons and (forced) triangulation polygons,
    // deduplicated by TShape index (L777: the same shape map).
    let mut seen_edges: HashSet<usize> = HashSet::new();
    let mut edge_occurrences: Vec<Shape> = Vec::new();
    for (_, edges) in &face_occurrences {
        for e in edges {
            if seen_edges.insert(e.index) {
                edge_occurrences.push(e.clone());
            }
        }
    }
    for (ei, ts) in brep.tshapes.iter().enumerate() {
        if matches!(&**ts, TShape::Edge(_)) && seen_edges.insert(ei) {
            edge_occurrences.push(Shape::from_parts(
                ts.clone(),
                ei,
                0,
                rcad_kernel::topods::Orientation::Forward,
            ));
        }
    }
    for an_edge_ref in &edge_occurrences {
        let mut an_edge = an_edge_ref.clone();
        // L775: anEdge.Location(anEmptyLoc) — the location is dropped before
        // the shape-map add; the representation reads below use the wrapper.
        an_edge.location = 0;
        let is_geometric = match &*brep.tshapes[an_edge.index] {
            // L783: BRep_Tool::IsGeometric(anEdge) — a 3D curve or a curve on
            // surface (BRep_Tool.cxx L238-258).
            TShape::Edge(te) => te.curve.is_some() || !te.pcurves.is_empty(),
            _ => false,
        };
        if is_geometric {
            // L786-791: remove the 3D polygon when present.
            let (a_poly3d, _loc) = brep.polygon3d(&an_edge);
            if a_poly3d.is_some() {
                brep.update_edge_polygon3d(&an_edge, None, 0);
            }
        }
        if the_force {
            // L793-812: find and remove ALL PolygonOnTriangulation
            // representations.
            let arc = std::sync::Arc::make_mut(&mut brep.tshapes[an_edge.index]);
            if let TShape::Edge(te) = &mut *arc {
                te.representations.retain(|cr| !cr.is_polygon_on_triangulation());
                // L811: aTE->Modified(true) — flag-only effect.
            }
        }
    }
}

/// The wire-ordered edge wrappers of a face wrapper (the
/// TopExp_Explorer(aFace, TopAbs_EDGE) stand-in: outer wire first, then the
/// inner wires, wrapper locations preserved).
fn face_wire_edges(
    brep: &rcad_kernel::topods::BRep,
    the_face: &rcad_kernel::topo::topo_shape::Shape,
) -> Vec<rcad_kernel::topo::topo_shape::Shape> {
    let mut out = Vec::new();
    let TShape::Face(tf) = &*brep.tshapes[the_face.index] else {
        return out;
    };
    let mut push_wire =
        |w: &rcad_kernel::topo::topo_shape::Shape, out: &mut Vec<rcad_kernel::topo::topo_shape::Shape>| {
            if let TShape::Wire(wd) = &*brep.tshapes[w.index] {
                for we in &wd.edges {
                    out.push(we.clone());
                }
            }
        };
    push_wire(&tf.outer_wire, &mut out);
    for w in &tf.inner_wires {
        push_wire(w, &mut out);
    }
    out
}
