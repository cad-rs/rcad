// GENERATED FILE — do not edit by hand.
//
// occt-test-gen: --occt-root "//?/C:/Users/lilu/works/OCCT" --group boolean --grid bcommon_2d --case A1
// OCCT case: tests/boolean/bcommon_2d/A1
//
// Set OCCT_SRC_ROOT to your OCCT repository root when running `cargo test`.
// (Same tree as --occt-root at generation time.)

use rcad_step::OcctBrepReader;
use std::path::{Path, PathBuf};

fn occt_src_root() -> PathBuf {
    std::env::var_os("OCCT_SRC_ROOT")
        .map(PathBuf::from)
        .expect("set OCCT_SRC_ROOT to OCCT repository root (contains tests/ and data/)")
}

fn locate_occt_data_file(basename: &str) -> PathBuf {
    let data_root = occt_src_root().join("data");
    find_file_named(&data_root, basename).unwrap_or_else(|| {
        panic!(
            "locate_data_file {} under {}",
            basename,
            data_root.display()
        );
    })
}

fn find_file_named(dir: &Path, basename: &str) -> Option<PathBuf> {
    let entries = std::fs::read_dir(dir).unwrap_or_else(|e| {
        panic!("read data directory {}: {}", dir.display(), e);
    });
    for entry in entries {
        let entry =
            entry.unwrap_or_else(|e| panic!("read directory entry in {}: {}", dir.display(), e));
        let path = entry.path();
        if path.file_name().and_then(|n| n.to_str()) == Some(basename) {
            return Some(path);
        }
        if path.is_dir() {
            if let Some(found) = find_file_named(&path, basename) {
                return Some(found);
            }
        }
    }
    None
}

fn read_occt_path(p: &Path) -> String {
    std::fs::read_to_string(&p).unwrap_or_else(|e| {
        panic!("read {}: {}", p.display(), e);
    })
}

#[test]
#[ignore = "requires OCCT_SRC_ROOT and data/case_1_edge1.brep, case_1_solid.brep from OCCT"]
fn occt_boolean_bcommon_2d_a1_geometry_loads() {
    {
        let rel = locate_occt_data_file("case_1_edge1.brep");
        let text = read_occt_path(&rel);
        OcctBrepReader::parse_string(&text).unwrap_or_else(|e| {
            panic!(
                "OcctBrepReader {} ({}): {:?}",
                rel.display(),
                "case_1_edge1.brep",
                e
            );
        });
    }
    {
        let rel = locate_occt_data_file("case_1_solid.brep");
        let text = read_occt_path(&rel);
        OcctBrepReader::parse_string(&text).unwrap_or_else(|e| {
            panic!(
                "OcctBrepReader {} ({}): {:?}",
                rel.display(),
                "case_1_solid.brep",
                e
            );
        });
    }
}

#[test]
#[ignore = "port Draw commands (bcommon, checknbshapes, …) to rcad"]
fn occt_boolean_bcommon_2d_a1_draw_script_port_pending() {
    // --- Original OCCT Draw script (verbatim lines) ---
    // restore [locate_data_file case_1_solid.brep] a
    // restore [locate_data_file case_1_edge1.brep] b
    //
    // bcommon result b a
    //
    // checkprops result -l 100.002
    // checksection result
    // checknbshapes result -vertex 4 -edge 2 -t
    //
    todo!("implement rcad equivalents of the Draw script above");
}
