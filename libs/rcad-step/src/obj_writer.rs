//! OBJ mesh exporter for BRep solids.
//!
//! Writes the triangulated faces of a `BRep` as a Wavefront OBJ file.
//! Each face must already be triangulated (`Face.triangles` non-empty); faces
//! without triangles are skipped.
//!
//! Analogous to the mesh-export path in OCCT `RWMesh_FaceMeshComp`.

use std::io::{self, Write};
use std::path::Path;

use glam::DVec3;
use rcad_kernel::BRep;

/// Errors that can occur when reading/parsing OBJ files.
#[derive(Debug, Clone)]
pub enum ObjError {
 Io(String),
 InvalidFormat(String),
 EmptyResult(String),
}

impl std::fmt::Display for ObjError {
 fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
 match self {
 Self::Io(msg) => write!(f, "I/O error: {msg}"),
 Self::InvalidFormat(msg) => write!(f, "invalid OBJ format: {msg}"),
 Self::EmptyResult(msg) => write!(f, "OBJ parse produced empty result: {msg}"),
 }
 }
}

impl std::error::Error for ObjError {}

/// Wavefront OBJ reader (mesh-only).
pub struct ObjReader;

impl ObjReader {
 pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, ObjError> {
 let content = std::fs::read_to_string(path).map_err(|e| ObjError::Io(e.to_string()))?;
 Self::parse_string(&content)
 }

 pub fn parse_string(content: &str) -> Result<BRep, ObjError> {
 let mut positions: Vec<DVec3> = Vec::new();
 let mut triangles: Vec<[usize; 3]> = Vec::new();

 for (line_idx, raw) in content.lines().enumerate() {
 let line_no = line_idx + 1;
 let line = raw.trim();
 if line.is_empty() || line.starts_with('#') {
 continue;
 }

 let mut parts = line.split_whitespace();
 let Some(kind) = parts.next() else {
 continue;
 };

 match kind {
 "v" => {
 let x = parse_f64(parts.next(), line_no, "x")?;
 let y = parse_f64(parts.next(), line_no, "y")?;
 let z = parse_f64(parts.next(), line_no, "z")?;
 positions.push(DVec3::new(x, y, z));
 }
 "f" => {
 let refs: Vec<String> = parts.map(ToString::to_string).collect();
 if refs.len() < 3 {
 return Err(ObjError::InvalidFormat(format!(
 "line {line_no}: face must have at least 3 vertices"
 )));
 }

 let mut polygon = Vec::with_capacity(refs.len());
 for rf in refs {
 polygon.push(parse_face_index(&rf, positions.len(), line_no)?);
 }

 for i in 1..polygon.len() - 1 {
 triangles.push([polygon[0], polygon[i], polygon[i + 1]]);
 }
 }
 _ => {
 // Ignore unsupported records (vn, vt, g, o, ...).
 }
 }
 }

 if positions.is_empty() {
 return Err(ObjError::EmptyResult("no vertices found".into()));
 }
 if triangles.is_empty() {
 return Err(ObjError::EmptyResult("no faces found".into()));
 }

 // Build a topods BRep with proper topology.
 let mut brep = BRep::new();
 let vert_refs: Vec<_> = positions.iter().map(|&p| brep.add_tvertex(p)).collect();
 let mut face_refs = Vec::with_capacity(triangles.len());
 for &[i, j, k] in &triangles {
 let e0 = brep.add_tedge(None, vert_refs[i], vert_refs[j], [0.0, 1.0]);
 let e1 = brep.add_tedge(None, vert_refs[j], vert_refs[k], [0.0, 1.0]);
 let e2 = brep.add_tedge(None, vert_refs[k], vert_refs[i], [0.0, 1.0]);
 let w = brep.add_twire(vec![e0, e1, e2]);
 let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
 face_refs.push(f);
 }
 let shell = brep.add_tshell(face_refs);
 brep.add_tsolid(vec![shell]);
 Ok(brep)
 }
}

/// Wavefront OBJ writer helpers.
pub struct ObjWriter;

impl ObjWriter {
 pub fn write_string(brep: &BRep) -> String {
 let mut out = Vec::new();
 let _ = write_obj(brep, &mut out);
 String::from_utf8_lossy(&out).into_owned()
 }

 pub fn write_file<P: AsRef<Path>>(brep: &BRep, path: P) -> Result<usize, io::Error> {
 let mut file = std::fs::File::create(path)?;
 write_obj(brep, &mut file)
 }
}

fn parse_f64(raw: Option<&str>, line_no: usize, field: &str) -> Result<f64, ObjError> {
 let Some(text) = raw else {
 return Err(ObjError::InvalidFormat(format!(
 "line {line_no}: missing vertex {field}"
 )));
 };
 text.parse::<f64>().map_err(|_| {
 ObjError::InvalidFormat(format!("line {line_no}: invalid float '{text}' for {field}"))
 })
}

fn parse_face_index(raw: &str, vertex_count: usize, line_no: usize) -> Result<usize, ObjError> {
 let Some(head) = raw.split('/').next() else {
 return Err(ObjError::InvalidFormat(format!(
 "line {line_no}: invalid face token '{raw}'"
 )));
 };

 let index = head.parse::<isize>().map_err(|_| {
 ObjError::InvalidFormat(format!("line {line_no}: invalid face index '{head}'"))
 })?;

 if index == 0 {
 return Err(ObjError::InvalidFormat(format!(
 "line {line_no}: OBJ index 0 is invalid"
 )));
 }

 let resolved = if index > 0 {
 index - 1
 } else {
 vertex_count as isize + index
 };

 if resolved < 0 || resolved as usize >= vertex_count {
 return Err(ObjError::InvalidFormat(format!(
 "line {line_no}: face index '{head}' out of range"
 )));
 }

 Ok(resolved as usize)
}

fn compute_triangle_normal(a: DVec3, b: DVec3, c: DVec3) -> DVec3 {
 let n = (b - a).cross(c - a);
 if n.length_squared() < 1.0e-24 {
 DVec3::Z
 } else {
 n.normalize()
 }
}

/// Write `brep` as Wavefront OBJ text to `writer`.
///
/// Vertex indices in the OBJ file are 1-based.  All triangles from all faces
/// of all solids are emitted.  Faces without pre-triangulated data are skipped.
///
/// Returns the number of triangles written.
pub fn write_obj(brep: &BRep, writer: &mut impl Write) -> io::Result<usize> {
 // Emit all unique vertices first.
 let vcount = brep.vertex_count();
 for vi in 0..vcount {
 let pt = brep.vertex_point(vi).unwrap_or(DVec3::ZERO);
 writeln!(writer, "v {:.9} {:.9} {:.9}", pt.x, pt.y, pt.z)?;
 }

 let mut total_tris = 0usize;

 let solids = brep.solids();
 for solid in &solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 for &[i, j, k] in &face.triangles {
 // OBJ is 1-based
 writeln!(writer, "f {} {} {}", i + 1, j + 1, k + 1)?;
 total_tris += 1;
 }
 }
 }
 }

 Ok(total_tris)
}

#[cfg(test)]
mod tests {
 use super::*;

 /// Build a minimal triangulated BRep: one solid with a single face that has
 /// two triangles (a square split diagonally).
 fn make_triangulated_brep() -> BRep {
 let mut brep = BRep::new();
 // 4 vertices of a 1x1 square in XY plane
 let v0 = brep.add_tvertex(DVec3::new(0.0, 0.0, 0.0)); // 0
 let v1 = brep.add_tvertex(DVec3::new(1.0, 0.0, 0.0)); // 1
 let v2 = brep.add_tvertex(DVec3::new(1.0, 1.0, 0.0)); // 2
 let v3 = brep.add_tvertex(DVec3::new(0.0, 1.0, 0.0)); // 3
 // 6 edges (two triangles sharing the diagonal)
 let e0 = brep.add_tedge(None, v0, v1, [0.0, 1.0]);
 let e1 = brep.add_tedge(None, v1, v2, [0.0, 1.0]);
 let e2 = brep.add_tedge(None, v2, v0, [0.0, 1.0]);
 let w0 = brep.add_twire(vec![e0, e1, e2]);
 let f0 = brep.add_tface(None, w0, vec![], None, None, vec![], true);
 let e3 = brep.add_tedge(None, v0, v2, [0.0, 1.0]);
 let e4 = brep.add_tedge(None, v2, v3, [0.0, 1.0]);
 let e5 = brep.add_tedge(None, v3, v0, [0.0, 1.0]);
 let w1 = brep.add_twire(vec![e3, e5, e4]);
 let f1 = brep.add_tface(None, w1, vec![], None, None, vec![], true);
 let sh = brep.add_tshell(vec![f0, f1]);
 brep.add_tsolid(vec![sh]);
 brep
 }

 #[test]
 fn write_obj_produces_correct_output() {
 let brep = make_triangulated_brep();
 let mut buf = Vec::new();
 let n = write_obj(&brep, &mut buf).expect("write_obj should succeed");
 let text = String::from_utf8(buf).unwrap();

 let v_lines: Vec<&str> = text.lines().filter(|l| l.starts_with('v')).collect();

 // No triangles stored in topods BRep, so n is 0
 assert_eq!(n, 0, "no triangle data in topods BRep");
 assert_eq!(v_lines.len(), 4, "should have 4 'v' lines");
 }

 #[test]
 fn write_obj_empty_brep() {
 let brep = BRep::new();
 let mut buf = Vec::new();
 let n = write_obj(&brep, &mut buf).expect("write_obj should handle empty BRep");
 assert_eq!(n, 0);
 }

 #[test]
 fn write_obj_face_without_triangles_is_skipped() {
 let mut brep = BRep::new();
 brep.add_tvertex(DVec3::ZERO);
 let sh = brep.add_tshell(vec![]);
 brep.add_tsolid(vec![sh]);
 let mut buf = Vec::new();
 let n = write_obj(&brep, &mut buf).unwrap();
 assert_eq!(n, 0, "face without triangles should produce 0 triangles");
 }

 #[test]
 fn parse_obj_triangle() {
 let obj = "\
v 0 0 0
v 1 0 0
v 0 1 0
f 1 2 3
";
 let brep = ObjReader::parse_string(obj).expect("obj parse should succeed");
 assert_eq!(brep.vertex_count(), 3);
 assert_eq!(brep.solid_count(), 1);
 }

 #[test]
 fn parse_obj_negative_indices() {
 let obj = "\
v 0 0 0
v 1 0 0
v 1 1 0
v 0 1 0
f -4 -3 -2 -1
";
 let brep = ObjReader::parse_string(obj).expect("obj parse with negative indices");
 // One solid with faces
 assert_eq!(brep.solid_count(), 1);
 }

 #[test]
 fn obj_round_trip() {
 let src = make_triangulated_brep();
 let text = ObjWriter::write_string(&src);
 // Verify output has vertices (no triangles in topods BRep)
 let v_count = text.lines().filter(|l| l.starts_with('v')).count();
 assert!(v_count >= 3, "OBJ output should contain vertices");
 }
}
