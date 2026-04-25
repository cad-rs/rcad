//! Parse every fixture listed in `tools/occt-test-importer/manifest.json` (paths relative to workspace root).
//! Run from repo root: `cargo test -p rcad-step --test occt_imported_fixtures`.

use rcad_step::{OcctBrepReader, StepReader};
use std::fs;
use std::path::PathBuf;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

#[test]
fn occt_imported_manifest_fixtures_parse() {
    let root = workspace_root();
    let manifest_path = root.join("tools/occt-test-importer/manifest.json");
    let raw = match fs::read_to_string(&manifest_path) {
        Ok(s) => s,
        Err(e) => panic!(
            "read manifest {}: {e}. Importer manifest is expected at workspace root.",
            manifest_path.display()
        ),
    };
    let v: serde_json::Value =
        serde_json::from_str(&raw).expect("tools/occt-test-importer/manifest.json must be valid JSON");
    let cases = v["cases"]
        .as_array()
        .expect("manifest must contain a \"cases\" array");

    for case in cases {
        let kind = case["kind"].as_str().expect("case.kind");
        let rel = case["fixture"].as_str().expect("case.fixture");
        let path = root.join(rel);
        let content = fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read fixture {}: {e}", path.display()));
        match kind {
            "brep" => {
                OcctBrepReader::parse_string(&content).unwrap_or_else(|e| {
                    panic!("OcctBrepReader parse {}: {e:?}", path.display());
                });
            }
            "step" => {
                StepReader::parse_string(&content).unwrap_or_else(|e| {
                    panic!("StepReader parse {}: {e:?}", path.display());
                });
            }
            other => panic!("unknown manifest case kind: {other}"),
        }
    }
}
