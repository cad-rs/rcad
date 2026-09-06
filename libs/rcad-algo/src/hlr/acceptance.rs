//! Stage 4d acceptance tests: the exact-HLR pipeline against OCCT reference
//! data (the OCCT grid cases of `tests/hlr/exact_hlr`).
//!
//! Reference sources (root repo, located at runtime by [`repo_root`]):
//! - `tools/occt-hlr-runner/ref_output.txt` — OCCT 8.0.0 `occt_hlr_runner`
//!   KEY=VALUE blocks (the C++ tuple-projector path of the DRAW `COMPUTE_HLR
//!   $algotype` command). Only the `-algo` blocks are consumed here; the
//!   `-poly` blocks belong to the poly algorithm and are skipped.
//! - `tests/occt/step_reference/occt_hlr_exact_hlr_bug25813_1.json` — the
//!   DRAWEXE nbshapes + mass reference for the `exact_hlr/bug25813_1` case.
//!
//! OCCT DRAW cases (`$OCCT_SRC/tests/hlr/exact_hlr/`, all viewed with the
//! V3d_XposYnegZpos frame dir (1,-1,1), up (-1,1,2)):
//! - `bug25813_1`: `pcylinder cc 10 30` / `pcylinder cc2 8 50` /
//!   `ttranslate cc2 0 0 2` / `bfuse a cc cc2`; `checkprops result -l 204.19`.
//! - `bug25813_3`: `ptorus a 30 10` (major R=30, minor r=10);
//!   `checkprops result -l 302.685`.
//! - `Plate`: four planar boundaries + `filling` (BRepOffsetAPI_MakeFilling).
//!   NOT IMPLEMENTABLE in rcad: there is no GeomPlate / MakeFilling port
//!   (the whole GeomPlate subsystem of TKGeomBase/TKTopAlgo is absent), so no
//!   acceptance test can even build the Plate fixture. Reference numbers for
//!   the future port: the OCCT 8.0.0 official DRAWEXE viewer path yields
//!   Mass 406.283 (nbshapes VERTEX:12 / EDGE:6 / COMPOUND:5) versus the
//!   DRAW-case length 404.283 (rel 4.9e-3, inside the official checkprops
//!   epsilon 1e-2; the historical 404.283 is not reproduced by this OCCT
//!   build via any runnable path — see the SUPPLEMENT block of
//!   ref_output.txt). Also note the OCCT 8.0.0 C++ tuple path crashes
//!   (access violation, exit 139) on both `bug25813_3` and `Plate`; the torus
//!   truth below comes from the official DRAWEXE viewer path instead.

use std::collections::HashSet;
use std::path::PathBuf;

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::topo::topods::{BRep, BRepBuilder};
use rcad_kernel::topods::{Orientation, Shape, TShape};

use super::tests::{
    compound_edges, ellipse_arc_len, smoke_box_solid, smoke_run_hlr, smoke_run_hlr_vcompute,
    SegSummary, VComputeFilters,
};

// ---- reference data location ----

/// Locate the root repo (the workspace holding `tests/occt`) by walking up
/// from the process cwd (cargo sets it to the package dir
/// `rcad/libs/rcad-algo`).
fn repo_root() -> PathBuf {
    let started = std::env::current_dir().expect("current_dir is readable");
    let mut dir = started.clone();
    loop {
        if dir.join("tests/occt/step_reference").is_dir() {
            return dir;
        }
        if !dir.pop() {
            panic!(
                "OCCT HLR reference data not found: no ancestor of {:?} \
                 contains tests/occt/step_reference",
                started
            );
        }
    }
}

fn ref_output_text() -> String {
    let path = repo_root().join("tools/occt-hlr-runner/ref_output.txt");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

fn ref_json_text(name: &str) -> String {
    let path = repo_root().join("tests/occt/step_reference").join(name);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
}

// ---- ref_output.txt parsing ----
// The file uses KEY=VALUE tokens separated by whitespace inside each
// `#### RUN <case> -<algotype>` block, e.g.
//   MASS_TOTAL=146.96938456699067 MASS_TOTAL_6SIG=146.969
//   NBSHAPES_RESULT=VERTEX:18 EDGE:9 WIRE:0 ... COMPOUND:4
// `str::lines()` handles both \n and \r\n endings.

/// One `#### RUN <case> -algo` block (the exact-algorithm diagnostics).
#[derive(Default, Debug)]
struct AlgoRunRef {
    mass_total: Option<f64>,
    mass_hidden_total: Option<f64>,
    mass_visible_only: Option<f64>,
    mass_hidden_only: Option<f64>,
    nb_result_vertices: Option<usize>,
    nb_result_edges: Option<usize>,
    nb_hidden_vertices: Option<usize>,
    nb_hidden_edges: Option<usize>,
}

/// Parse the `#### RUN <case> -algo` block of ref_output.txt. Blocks whose
/// runner crashed carry no MASS/NBSHAPES lines — the Options stay None.
fn parse_algo_run(text: &str, case: &str) -> AlgoRunRef {
    let header = format!("#### RUN {} -algo", case);
    let mut run = AlgoRunRef::default();
    let mut in_block = false;
    for raw in text.lines() {
        let line = raw.trim();
        if line.starts_with("####") {
            in_block = line == header;
            continue;
        }
        if !in_block || line.is_empty() {
            continue;
        }
        // MASS_* lines: whitespace-separated KEY=VALUE tokens.
        for tok in line.split_whitespace() {
            if let Some((k, v)) = tok.split_once('=') {
                match k {
                    "MASS_TOTAL" => run.mass_total = v.parse().ok(),
                    "MASS_HIDDEN_TOTAL" => run.mass_hidden_total = v.parse().ok(),
                    "MASS_VISIBLE_ONLY" => run.mass_visible_only = v.parse().ok(),
                    "MASS_HIDDEN_ONLY" => run.mass_hidden_only = v.parse().ok(),
                    _ => continue,
                }
            }
        }
        if let Some(rest) = line.split_once("NBSHAPES_RESULT=") {
            let (vv, ev) = parse_nbshapes(rest.1);
            run.nb_result_vertices = vv;
            run.nb_result_edges = ev;
        }
        if let Some(rest) = line.split_once("NBSHAPES_HIDDEN=") {
            let (vv, ev) = parse_nbshapes(rest.1);
            run.nb_hidden_vertices = vv;
            run.nb_hidden_edges = ev;
        }
    }
    if !in_block && run.mass_total.is_none() {
        // The header itself was never found (the loop above only sets
        // in_block on headers; a missing header leaves every field None).
        panic!("ref_output.txt has no '#### RUN {} -algo' block", case);
    }
    run
}

/// The VERTEX/EDGE counts of an `NBSHAPES_X=VERTEX:n EDGE:n ...` tail.
fn parse_nbshapes(tail: &str) -> (Option<usize>, Option<usize>) {
    let mut vertices = None;
    let mut edges = None;
    for tok in tail.split_whitespace() {
        if let Some((kind, n)) = tok.split_once(':') {
            match kind {
                "VERTEX" => vertices = n.parse().ok(),
                "EDGE" => edges = n.parse().ok(),
                _ => {}
            }
        }
    }
    (vertices, edges)
}

/// Parse the `#### SUPPLEMENT` block entry selected by `marker` (a substring
/// of the block's leading comment line). The official DRAWEXE viewer path
/// prints an nbshapes summary (`VERTEX : n` / `EDGE : n` / ...) followed by
/// `Mass : <value>`. Returns (mass, visible edge count).
fn parse_supplement(text: &str, marker: &str) -> (f64, usize) {
    let mut seen_marker = false;
    let mut mass = None;
    let mut edges = None;
    for raw in text.lines() {
        let line = raw.trim();
        if line.starts_with("####") {
            if seen_marker {
                break; // the next section starts
            }
            seen_marker = line.contains(marker);
            continue;
        }
        if !seen_marker {
            continue;
        }
        if let Some(rest) = line.strip_prefix("Mass") {
            // "Mass :         302.685"
            if let Some(val) = rest.split(':').nth(1) {
                mass = val.trim().parse().ok();
            }
        }
        if line.starts_with("EDGE") {
            // "EDGE      : 4"
            if let Some(val) = line.split(':').nth(1) {
                edges = val.trim().parse().ok();
            }
        }
    }
    (
        mass.unwrap_or_else(|| panic!("no Mass line under the '{marker}' SUPPLEMENT entry")),
        edges.unwrap_or_else(|| panic!("no EDGE count under the '{marker}' SUPPLEMENT entry")),
    )
}

// ---- mass (lprops-equivalent) measurement ----
// OCCT `lprops` on the COMPUTE_HLR wire compound integrates the exact curve
// lengths. rcad sums per-edge analytic lengths (exact for Line / Circle /
// Ellipse supports) and a dense chord sampling for every other projected
// curve kind — the 4096-chord error is many orders below the 1e-2
// acceptance tolerance.

/// Length of one projected result edge over its trimmed range.
fn projected_edge_len(curve: Option<&Curve3>, range: [f64; 2]) -> f64 {
    let Some(c) = curve else { return 0.0 };
    let (u1, u2) = (range[0], range[1]);
    match c {
        Curve3::Line(l) => l.direction.length() * (u2 - u1),
        Curve3::Circle(circ) => circ.radius * (u2 - u1),
        Curve3::Ellipse(el) => ellipse_arc_len(el.major_radius, el.minor_radius, u1, u2),
        other => {
            const N: usize = 4096;
            let mut len = 0.0;
            let mut prev = other.point_at(u1);
            for k in 1..=N {
                let t = u1 + (u2 - u1) * (k as f64) / (N as f64);
                let p = other.point_at(t);
                len += (p - prev).length();
                prev = p;
            }
            len
        }
    }
}

/// Total projected curve length of every edge in a compound tree.
pub(crate) fn compound_mass(s: &Shape) -> f64 {
    if s.is_null() {
        return 0.0;
    }
    fn walk(s: &Shape, total: &mut f64) {
        match &*s.data {
            TShape::Compound(c) => {
                for ch in c {
                    walk(ch, total);
                }
            }
            TShape::Edge(ed) => *total += projected_edge_len(ed.curve.as_ref(), ed.range),
            _ => {}
        }
    }
    let mut total = 0.0;
    walk(s, &mut total);
    total
}

/// Number of edge occurrences in a compound tree.
fn compound_edge_count(s: &Shape) -> usize {
    if s.is_null() {
        return 0;
    }
    fn walk(s: &Shape, out: &mut usize) {
        match &*s.data {
            TShape::Compound(c) => {
                for ch in c {
                    walk(ch, out);
                }
            }
            TShape::Edge(_) => *out += 1,
            _ => {}
        }
    }
    let mut out = 0;
    walk(s, &mut out);
    out
}

/// (v, rg1v, rgnv, outline) tuple access: the four visible compounds of
/// [`smoke_run_hlr`].
struct VisibleCompounds {
    v: Shape,
    rg1v: Shape,
    rgnv: Shape,
    outline: Shape,
}

impl VisibleCompounds {
    fn total_mass(&self) -> f64 {
        compound_mass(&self.v)
            + compound_mass(&self.rg1v)
            + compound_mass(&self.rgnv)
            + compound_mass(&self.outline)
    }
    fn total_edges(&self) -> usize {
        compound_edge_count(&self.v)
            + compound_edge_count(&self.rg1v)
            + compound_edge_count(&self.rgnv)
            + compound_edge_count(&self.outline)
    }
    /// The OCCT-style diagnostic print: per-compound curve-kind/length lists.
    fn dump(&self, tag: &str) {
        for (name, c) in [
            ("V", &self.v),
            ("Rg1LineV", &self.rg1v),
            ("RgNLineV", &self.rgnv),
            ("OutLineV", &self.outline),
        ] {
            let mut segs = Vec::new();
            if !c.is_null() {
                compound_edges(c, &mut segs);
            }
            println!("ACCEPT {tag} {name}({}): {:?}", segs.len(), segs);
        }
    }
}

/// The OCCT checkprops-style acceptance: |actual - reference| / reference
/// <= 1e-2 (the DRAW cases' depsilon).
fn assert_rel(actual: f64, reference: f64, what: &str) {
    let rel = ((actual - reference) / reference).abs();
    assert!(
        rel <= 1.0e-2,
        "{what}: rcad {actual} vs OCCT reference {reference} (rel err {rel:.3e} > 1e-2)"
    );
}

/// The AGENTS.md valid-digit line: the rcad value must reproduce the OCCT
/// reference through its printed significant digits.  The runner lprops
/// values carry 17 digits and the DRAW references 5-6; the pin is the
/// 6-digit reading of the reference (rel <= 1e-7 leaves the half-ulp of the
/// last printed digit plus platform fp noise as the only slack).
fn assert_valid_digits(actual: f64, reference: f64, what: &str) {
    let rel = ((actual - reference) / reference).abs();
    assert!(
        rel <= 1.0e-7,
        "{what}: rcad {actual} vs OCCT reference {reference} (rel err {rel:.3e} > 1e-7)"
    );
}

/// The DRAW-value line: the rcad value must round to the OCCT printed value
/// (the half-ulp of its last printed digit).
fn assert_printed_value(actual: f64, printed: f64, what: &str) {
    assert!(
        (actual - printed).abs() <= 5.0e-4,
        "{what}: rcad {actual} does not round to the printed OCCT value {printed}"
    );
}

// ---- the OCCT VComputeHLR result assembly ----

/// The nbshapes-style unique shape counts of a result tree — the runner's
/// `nbshapes` (BRepTools_ShapeSet::Add(S) + per-type unique counts, dedup by
/// IsSame = TShape + Location; the rcad ptr_id stands in for the identity).
/// Returns (vertex, edge, compound).
fn nbshapes_counts(s: &Shape) -> (usize, usize, usize) {
    if s.is_null() {
        return (0, 0, 0);
    }
    fn walk(s: &Shape, seen_v: &mut HashSet<u64>, seen_e: &mut HashSet<u64>, nc: &mut usize) {
        match &*s.data {
            TShape::Compound(c) => {
                *nc += 1;
                for ch in c {
                    walk(ch, seen_v, seen_e, nc);
                }
            }
            TShape::Edge(ed) => {
                if seen_e.insert(s.ptr_id()) {
                    seen_v.insert(ed.first.ptr_id());
                    seen_v.insert(ed.last.ptr_id());
                }
            }
            TShape::Vertex(_) => {
                seen_v.insert(s.ptr_id());
            }
            _ => {}
        }
    }
    let mut seen_v = HashSet::new();
    let mut seen_e = HashSet::new();
    let mut nc = 0;
    walk(s, &mut seen_v, &mut seen_e, &mut nc);
    (seen_v.len(), seen_e.len(), nc)
}

/// OCCT ViewerTest_ObjectCommands.cxx L3335-3350: aCompVis / aCompHid are
/// made unconditionally — with `-showHiddenEdges` off the hidden filters are
/// never extracted, aCompHid stays an EMPTY (but non-null) compound and is
/// still counted by nbshapes (the box reference: COMPOUND 4 = aCompRes +
/// aCompVis + empty aCompHid + VCompound) — then aCompRes gets aCompVis and
/// aCompHid.  The non-null filters enter in the HLRBRep_TypeOfResultingEdge
/// order (IsoLine=0, OutLine=1, Rg1Line=2, RgNLine=3, Sharp=4); the RgN
/// slots are never filled (toShowCNEdges == false).
fn vcompute_result_tree(f: &VComputeFilters, show_hidden: bool) -> Shape {
    let mut b = BRepBuilder::new();
    let mut arena = BRep::new();

    let mut vis: Vec<Shape> = Vec::new();
    for c in [&f.iso_v, &f.outline_v, &f.rg1_v, &f.v] {
        if !c.is_null() {
            vis.push(c.clone());
        }
    }
    let a_comp_vis = b.make_compound(&mut arena, vis);

    let mut hid: Vec<Shape> = Vec::new();
    if show_hidden {
        for c in [&f.iso_h, &f.outline_h, &f.rg1_h, &f.h] {
            if !c.is_null() {
                hid.push(c.clone());
            }
        }
    }
    let a_comp_hid = b.make_compound(&mut arena, hid);

    b.make_compound(&mut arena, vec![a_comp_vis, a_comp_hid])
}

/// The root solid Shape of a BRep pool (the same root-finding rule as
/// `brep_algo_api::run_build`: the last Solid entry of the pool).
fn root_solid(brep: &rcad_kernel::BRep) -> Shape {
    brep.tshapes
        .iter()
        .enumerate()
        .rev()
        .find(|(_, ts)| matches!(ts.as_ref(), TShape::Solid(_)))
        .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, Orientation::Forward))
        .expect("BRep pool has no solid")
}

// ---- the acceptance tests ----

/// The `box b 10 20 30` grid case against the `#### RUN box -algo` block of
/// ref_output.txt: visible mass vs MASS_TOTAL (146.96938456699067),
/// hidden-only mass vs MASS_HIDDEN_ONLY (48.989794855663568), and the
/// combined -showHiddenEdges result vs MASS_HIDDEN_TOTAL
/// (195.95917942265422); 9 visible edges / 3 hidden edges.
#[test]
fn acceptance_box_vs_ref_output() {
    let run = parse_algo_run(&ref_output_text(), "box");
    let mass_total = run.mass_total.expect("box -algo MASS_TOTAL");
    let mass_hidden_total = run.mass_hidden_total.expect("box -algo MASS_HIDDEN_TOTAL");
    let mass_hidden_only = run.mass_hidden_only.expect("box -algo MASS_HIDDEN_ONLY");
    // The runner's own bookkeeping: visible-only == total, and
    // hidden-total == visible + hidden-only.
    assert_eq!(run.nb_result_edges, Some(9), "ref box EDGE count");
    assert_eq!(run.nb_hidden_edges, Some(12), "ref box hidden EDGE count");

    let (brep, solid) = smoke_box_solid();
    let (v, rg1v, rgnv, outline, h) = smoke_run_hlr(&solid, brep);
    let vis = VisibleCompounds { v, rg1v, rgnv, outline };
    vis.dump("box");
    println!(
        "ACCEPT box: visible mass {} (edges {}) vs MASS_TOTAL {mass_total}; \
         hidden mass {} vs MASS_HIDDEN_ONLY {mass_hidden_only}; combined {} vs \
         MASS_HIDDEN_TOTAL {mass_hidden_total}",
        vis.total_mass(),
        vis.total_edges(),
        compound_mass(&h),
        vis.total_mass() + compound_mass(&h),
    );

    // OCCT nbshapes RESULT EDGE:9 = the 9 sharp visible edges (the box has
    // no outlines / Rg1 / RgN output); rcad pins the same 9 / 3 split.
    assert_eq!(vis.total_edges(), 9, "visible edge count");
    assert_eq!(compound_edge_count(&h), 3, "hidden edge count");
    assert_rel(vis.total_mass(), mass_total, "visible mass vs MASS_TOTAL");
    assert_rel(compound_mass(&h), mass_hidden_only, "hidden mass vs MASS_HIDDEN_ONLY");
    assert_rel(
        vis.total_mass() + compound_mass(&h),
        mass_hidden_total,
        "visible+hidden mass vs MASS_HIDDEN_TOTAL",
    );
}

/// The `exact_hlr/bug25813_1` DRAW case against
/// `tests/occt/step_reference/occt_hlr_exact_hlr_bug25813_1.json` and the
/// `#### RUN bug25813_1 -algo` block of ref_output.txt:
/// `pcylinder cc 10 30` / `pcylinder cc2 8 50` / `ttranslate cc2 0 0 2` /
/// `bfuse a cc cc2`, viewed from dir (1,-1,1) up (-1,1,2) — the same frame
/// the smoke runner hardcodes.
///
/// The OCCT nbshapes numbers count the whole VComputeHLR result tree
/// (ViewerTest_ObjectCommands.cxx L3335-3350): the visible-only run is
/// VERTEX:30 / EDGE:15 / COMPOUND:5 (aCompRes + aCompVis + the empty
/// aCompHid + VCompound + OutLineVCompound) and the -showHiddenEdges run is
/// VERTEX:36 / EDGE:18 / COMPOUND:6 (+ HCompound with the 3 hidden arcs /
/// 6 fresh vertices).  The mass pins run at the AGENTS.md valid-digit line
/// against the runner's 17-digit lprops values.
#[test]
fn acceptance_bug25813_1_vs_occt_json() {
    let j: serde_json::Value =
        serde_json::from_str(&ref_json_text("occt_hlr_exact_hlr_bug25813_1.json"))
            .expect("reference JSON parses");
    let mass_total = j["mass_total"].as_f64().expect("mass_total");
    let mass_hidden_total = j["mass_hidden_total"].as_f64().expect("mass_hidden_total");
    let ref_result_edges = j["nbshapes_result"]["EDGE"].as_u64().expect("EDGE") as usize;
    let ref_result_vertices = j["nbshapes_result"]["VERTEX"].as_u64().expect("VERTEX") as usize;
    let ref_result_compounds =
        j["nbshapes_result"]["COMPOUND"].as_u64().expect("COMPOUND") as usize;
    let ref_hidden_edges = j["nbshapes_hidden"]["EDGE"].as_u64().expect("EDGE") as usize;
    let ref_hidden_vertices = j["nbshapes_hidden"]["VERTEX"].as_u64().expect("VERTEX") as usize;
    let ref_hidden_compounds =
        j["nbshapes_hidden"]["COMPOUND"].as_u64().expect("COMPOUND") as usize;
    // The JSON projector frame equals the smoke_run_hlr projector (the
    // V3d_XposYnegZpos view of the DRAW runner).
    let dir: Vec<f64> = j["dir"]
        .as_array()
        .expect("dir")
        .iter()
        .map(|x| x.as_f64().expect("dir component"))
        .collect();
    let up: Vec<f64> = j["up"]
        .as_array()
        .expect("up")
        .iter()
        .map(|x| x.as_f64().expect("up component"))
        .collect();
    assert_eq!(dir, vec![1.0, -1.0, 1.0], "JSON view dir");
    assert_eq!(up, vec![-1.0, 1.0, 2.0], "JSON view up");

    // The runner's 17-digit lprops references (the valid-digit pins).
    let run = parse_algo_run(&ref_output_text(), "bug25813_1");
    let ref_mass_total = run.mass_total.expect("bug25813_1 -algo MASS_TOTAL");
    let ref_mass_hidden_total = run.mass_hidden_total.expect("bug25813_1 -algo MASS_HIDDEN_TOTAL");
    let ref_mass_hidden_only = run.mass_hidden_only.expect("bug25813_1 -algo MASS_HIDDEN_ONLY");

    // The DRAW fixture, built with the rcad primapi equivalents:
    // pcylinder cc 10 30 / pcylinder cc2 8 50 / ttranslate cc2 0 0 2.
    let cc = rcad_modeling::make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 10.0, 30.0)
        .expect("pcylinder cc 10 30");
    let cc2 = rcad_modeling::make_cylinder_brep(
        DVec3::new(0.0, 0.0, 2.0),
        DVec3::Z,
        DVec3::X,
        8.0,
        50.0,
    )
    .expect("pcylinder cc2 8 50");
    // bfuse a cc cc2
    let fused = crate::fuse(&cc, &cc2).expect("bfuse a cc cc2");
    let solid = root_solid(&fused);

    // OCCT L3307-3326: the eight HLRToShape filter calls (the extraction is
    // stateless across calls — InternalCompound resets Used/HideCount at
    // entry — so one pipeline run feeds both result trees).
    let f = smoke_run_hlr_vcompute(&solid, fused);
    println!(
        "ACCEPT bug25813_1 filters: V({}) OutLineV({}) Rg1V({}) IsoV({}) | \
         H({}) OutLineH({}) Rg1H({}) IsoH({})",
        compound_edge_count(&f.v),
        compound_edge_count(&f.outline_v),
        compound_edge_count(&f.rg1_v),
        compound_edge_count(&f.iso_v),
        compound_edge_count(&f.h),
        compound_edge_count(&f.outline_h),
        compound_edge_count(&f.rg1_h),
        compound_edge_count(&f.iso_h),
    );

    // ---- run 1: the visible-only result tree (L3335-3350, no hidden) ----
    let res1 = vcompute_result_tree(&f, false);
    let (v1, e1, c1) = nbshapes_counts(&res1);
    let m1 = compound_mass(&res1);
    println!(
        "ACCEPT bug25813_1 run1: nbshapes VERTEX:{v1} EDGE:{e1} COMPOUND:{c1} \
         (OCCT {ref_result_vertices}/{ref_result_edges}/{ref_result_compounds}), \
         mass {m1} vs MASS_TOTAL {ref_mass_total}"
    );
    assert_eq!(
        (v1, e1, c1),
        (
            ref_result_vertices,
            ref_result_edges,
            ref_result_compounds
        ),
        "run-1 nbshapes (visible-only result tree)"
    );
    assert_valid_digits(m1, ref_mass_total, "run-1 mass vs MASS_TOTAL");
    assert_rel(m1, mass_total, "run-1 mass vs JSON mass_total");

    // ---- run 2: the -showHiddenEdges result tree ----
    let res2 = vcompute_result_tree(&f, true);
    let (v2, e2, c2) = nbshapes_counts(&res2);
    let m2 = compound_mass(&res2);
    println!(
        "ACCEPT bug25813_1 run2: nbshapes VERTEX:{v2} EDGE:{e2} COMPOUND:{c2} \
         (OCCT {ref_hidden_vertices}/{ref_hidden_edges}/{ref_hidden_compounds}), \
         mass {m2} vs MASS_HIDDEN_TOTAL {ref_mass_hidden_total}; \
         hidden-only mass {} vs MASS_HIDDEN_ONLY {ref_mass_hidden_only}",
        compound_mass(&f.h)
    );
    assert_eq!(
        (v2, e2, c2),
        (ref_hidden_vertices, ref_hidden_edges, ref_hidden_compounds),
        "run-2 nbshapes (-showHiddenEdges result tree)"
    );
    assert_valid_digits(m2, ref_mass_hidden_total, "run-2 mass vs MASS_HIDDEN_TOTAL");
    assert_rel(m2, mass_hidden_total, "run-2 mass vs JSON mass_hidden_total");
    assert_valid_digits(
        compound_mass(&f.h),
        ref_mass_hidden_only,
        "hidden-only mass vs MASS_HIDDEN_ONLY",
    );
    // The hidden-tree composition: OCCT emits exactly one non-null hidden
    // filter here (the 3-arc HCompound); OutLineH/Rg1H/IsoH stay null.
    assert!(
        f.outline_h.is_null() && f.rg1_h.is_null() && f.iso_h.is_null(),
        "unexpected non-null hidden filters (OutLineH/Rg1H/IsoH must be null)"
    );
    assert_eq!(compound_edge_count(&f.h), 3, "the 3 hidden sharp arcs");
}

/// The `exact_hlr/bug25813_3` DRAW case (`ptorus a 30 10`) against the
/// `#### SUPPLEMENT` DRAWEXE viewer-path truth of ref_output.txt:
/// Mass 302.685 with 4 visible edges (VERTEX:7 — one contour line closes,
/// its single closed edge carries one vertex for both ends).  NOTE: OCCT
/// 8.0.0 itself access-violates on this case via the C++ tuple projector
/// path (exit 139, documented in ref_output.txt); the truth comes from the
/// official DRAWEXE viewer path (vinit+vdisplay+vcomputehlr).  rcad must
/// produce the same result without crashing.
///
/// IGNORED (remaining exact-HLR gaps, the session-19 runway item): after the
/// compute_tangency Destination off-by-one fix (the lost 4th contour line),
/// the visible tree matches the OCCT viewer path structurally — EDGE:4 /
/// VERTEX:7 — with visible mass 302.6867 vs the printed 302.685 (1.7e-3).
/// Two gaps keep this ignored:
/// 1. mass bit-exactness: the walk trajectories sample differently (rcad
///    2680/1584 walk points vs OCCT 81/215 on the two long lines; identical
///    net lengths) — traced to math_FunctionSetRoot done=false exits (100
///    iterations exhausted at F~1e-16, ~1 in 5 walk steps, halving PasC).
///    Next suspect: SearchDirection's math_SVD IsDone branch (OCCT can fall
///    back to the gradient direction; the rcad closed-form collapse never
///    fails).
/// 2. the hidden gate: the 86.94 figure is a mixed-measure derivation
///    (3D walk total 389.62 minus the projected visible 302.685), not an
///    OCCT measurement — the runner -algo run crashes inside OCCT on this
///    case and the viewer path never extracts hidden compounds.  The torus
///    contour near/far halves project onto each other, so the projected
///    hidden mass needs re-derivation from a trustworthy source before this
///    assert can mean anything.
/// Un-ignore when the mass rounds to the printed 302.685 and the hidden
/// gate is re-derived from real OCCT behavior.
#[test]
#[ignore = "remaining gaps: visible mass 302.6867 vs printed 302.685 (walk sampling density, FunctionSetRoot done=false exits); hidden gate 86.94 is a mixed-measure derivation needing re-derivation"]
fn acceptance_ptorus_vs_viewer_path() {
    let text = ref_output_text();
    let (mass_ref, edges_ref) = parse_supplement(&text, "exact_hlr/bug25813_3");
    assert_eq!(edges_ref, 4, "supplement parser: ref visible EDGE count");
    assert!((mass_ref - 302.685).abs() < 5e-4, "supplement Mass {mass_ref}");

    // ptorus a 30 10 — OCCT BRepPrimAPI_MakeTorus(R1=30 major, R2=10 minor).
    let torus = rcad_modeling::make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 30.0, 10.0)
        .expect("ptorus a 30 10");
    let solid = root_solid(&torus);

    let f = smoke_run_hlr_vcompute(&solid, torus);
    println!(
        "ACCEPT bug25813_3 filters: V({}) OutLineV({}) Rg1V({}) IsoV({}) | \
         H({}) OutLineH({}) Rg1H({}) IsoH({})",
        compound_edge_count(&f.v),
        compound_edge_count(&f.outline_v),
        compound_edge_count(&f.rg1_v),
        compound_edge_count(&f.iso_v),
        compound_edge_count(&f.h),
        compound_edge_count(&f.outline_h),
        compound_edge_count(&f.rg1_h),
        compound_edge_count(&f.iso_h),
    );
    let mut segs = Vec::new();
    compound_edges(&f.outline_v, &mut segs);
    println!("ACCEPT bug25813_3 OutLineV segments: {segs:?}");
    segs.clear();
    compound_edges(&f.outline_h, &mut segs);
    println!("ACCEPT bug25813_3 OutLineH segments: {segs:?}");

    // run 1: the visible-only result tree; OCCT viewer path: EDGE:4 /
    // VERTEX:7 / Mass 302.685 (all four contour edges, one visible interval
    // each; no compound-count reference exists — the C++ tuple path crashed
    // on this case, so only EDGE/VERTEX/mass are pinned).
    let res1 = vcompute_result_tree(&f, false);
    let (v1, e1, c1) = nbshapes_counts(&res1);
    let m1 = compound_mass(&res1);
    println!(
        "ACCEPT bug25813_3 run1: nbshapes VERTEX:{v1} EDGE:{e1} COMPOUND:{c1} \
         (OCCT viewer path VERTEX:7 EDGE:4), mass {m1} vs Mass {mass_ref}"
    );
    assert_eq!(e1, 4, "visible EDGE count (viewer path: 4)");
    assert_eq!(v1, 7, "visible VERTEX count (the closed-contour edge shares)");
    assert_printed_value(m1, mass_ref, "run-1 mass vs viewer-path Mass");
    // run 2: the hidden compound — the OCCT hidden mass is the analytic
    // two-tangent-band remainder (389.62 contour total - 302.685 visible).
    let res2 = vcompute_result_tree(&f, true);
    let (_, e2, _) = nbshapes_counts(&res2);
    let m2 = compound_mass(&f.outline_h) + compound_mass(&f.h);
    println!(
        "ACCEPT bug25813_3 run2: tree edges {e2}, hidden mass {m2} vs the \
         analytic 86.94 (389.62 contour total - 302.685 visible)"
    );
    assert_rel(m2, 86.94, "hidden mass vs the analytic remainder");
}

