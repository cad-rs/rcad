//! Example: OBJ / IGES mesh exchange.
//!
//! Demonstrates:
//!   1. Export a primitive BRep to OBJ and IGES strings/files
//!   2. Read those files back into BRep
//!   3. Parse bundled sample assets from the repository
//!
//! Run:
//!   cargo run -p rcad-examples --example obj_iges_io

use std::path::{Path, PathBuf};

use rcad_kernel::{BRep, PrimitiveSolid};
use rcad_step::{IgesReader, IgesWriter, ObjReader, ObjWriter};

fn triangle_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|solid| &solid.shells)
        .flat_map(|shell| &shell.faces)
        .map(|face| face.triangles.len())
        .sum()
}

fn summary(label: &str, brep: &BRep) {
    println!(
        "  {label}: vertices={} solids={} triangles={}",
        brep.vertices.len(),
        brep.solids.len(),
        triangle_count(brep)
    );
}

fn asset_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../assets")
        .join(name)
}

fn write_and_read_obj(brep: &BRep, out_dir: &Path) {
    println!("\n=== 1. OBJ export/import ===");

    let obj_text = ObjWriter::write_string(brep);
    println!("  OBJ text bytes={}", obj_text.len());

    let obj_path = out_dir.join("rcad_obj_demo.obj");
    let tri_count = ObjWriter::write_file(brep, &obj_path).expect("write OBJ file");
    println!("  wrote {} triangles -> {}", tri_count, obj_path.display());

    let reparsed = ObjReader::read_file(&obj_path).expect("read OBJ file");
    summary("OBJ round-trip", &reparsed);
}

fn write_and_read_iges(brep: &BRep, out_dir: &Path) {
    println!("\n=== 2. IGES export/import ===");

    let iges_text = IgesWriter::write_string(brep);
    println!("  IGES text bytes={}", iges_text.len());

    let iges_path = out_dir.join("rcad_obj_demo.igs");
    let tri_count = IgesWriter::write_file(brep, &iges_path).expect("write IGES file");
    println!("  wrote {} triangles -> {}", tri_count, iges_path.display());

    let reparsed = IgesReader::read_file(&iges_path).expect("read IGES file");
    summary("IGES round-trip", &reparsed);
}

fn parse_bundled_assets() {
    println!("\n=== 3. Bundled sample assets ===");

    let obj_asset = asset_path("sample_mesh.obj");
    let obj_brep = ObjReader::read_file(&obj_asset).expect("parse bundled OBJ asset");
    summary("asset OBJ", &obj_brep);

    let iges_asset = asset_path("sample_mesh.igs");
    let iges_brep = IgesReader::read_file(&iges_asset).expect("parse bundled IGES asset");
    summary("asset IGES", &iges_brep);
}

fn main() {
    println!("OBJ / IGES mesh exchange demo");

    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 1.5,
        depth: 1.0,
    });
    summary("source", &brep);

    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../target/tmp");
    std::fs::create_dir_all(&out_dir).expect("create output directory");

    write_and_read_obj(&brep, &out_dir);
    write_and_read_iges(&brep, &out_dir);
    parse_bundled_assets();

    println!("\nArtifacts written under {}", out_dir.display());
}
