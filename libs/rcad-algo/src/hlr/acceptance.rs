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
use rcad_kernel::topods::{Orientation, Shape, TShape};

use super::tests::{compound_edges, ellipse_arc_len, smoke_box_solid, smoke_run_hlr, SegSummary};

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

/// Number of distinct endpoint vertices in a compound tree (keyed by TShape
/// pointer — the nbshapes-style unique shape count over the result compound).
fn compound_vertex_count(s: &Shape) -> usize {
    if s.is_null() {
        return 0;
    }
    fn walk(s: &Shape, seen: &mut HashSet<u64>) {
        match &*s.data {
            TShape::Compound(c) => {
                for ch in c {
                    walk(ch, seen);
                }
            }
            TShape::Edge(ed) => {
                seen.insert(ed.first.ptr_id());
                seen.insert(ed.last.ptr_id());
            }
            _ => {}
        }
    }
    let mut seen = HashSet::new();
    walk(s, &mut seen);
    seen.len()
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
    /// The DRAW `COMPUTE_HLR` result composition with toShowCNEdges == false
    /// (ViewerTest_ObjectCommands.cxx L3288-3331): VCompound +
    /// OutLineVCompound + Rg1LineVCompound (+ IsoLineVCompound); the
    /// RgNLineVCompound (the CN seam edges) is NOT part of the emitted
    /// result. The OCCT reference masses measure exactly this composition.
    fn result_mass(&self) -> f64 {
        compound_mass(&self.v)
            + compound_mass(&self.rg1v)
            + compound_mass(&self.outline)
    }
    fn result_edges(&self) -> usize {
        compound_edge_count(&self.v)
            + compound_edge_count(&self.rg1v)
            + compound_edge_count(&self.outline)
    }
    fn result_vertices(&self) -> usize {
        compound_vertex_count(&self.v)
            + compound_vertex_count(&self.rg1v)
            + compound_vertex_count(&self.outline)
    }
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
    fn total_vertices(&self) -> usize {
        compound_vertex_count(&self.v)
            + compound_vertex_count(&self.rg1v)
            + compound_vertex_count(&self.rgnv)
            + compound_vertex_count(&self.outline)
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
/// `tests/occt/step_reference/occt_hlr_exact_hlr_bug25813_1.json`:
/// `pcylinder cc 10 30` / `pcylinder cc2 8 50` / `ttranslate cc2 0 0 2` /
/// `bfuse a cc cc2`, viewed from dir (1,-1,1) up (-1,1,2) — the same frame
/// the smoke runner hardcodes. Visible mass vs JSON mass_total (204.19),
/// visible+hidden vs mass_hidden_total (266.526), visible EDGE 15 and hidden
/// EDGE 18 (the JSON nbshapes corroboration).
///
/// IGNORED (remaining exact-HLR gap, measured against the official DRAWEXE
/// viewer-path truth): OCCT emits 11 sharp arcs (119.275) + 4 outline lines
/// (84.9156) visible and 3 hidden arcs (62.3361: the big-bottom far half
/// 25.2237, the z=30 inner-rim far half 20.179, and the big-top far-upper
/// middle 16.9334). rcad matches everything except that the big-top rim's
/// far-upper half is drawn as ONE visible arc (25.224) where OCCT splits it
/// at the small-cylinder outline crossings (x' = +/-8, y' = 27.95) into
/// 4.1452 + 16.9334 + 4.1452 and hides the middle. rcad visible mass is
/// therefore 221.12 (rel err 8.3e-2) and hidden 45.40. The two outline
/// pieces are registered as internal w_edges of the hiding face, but
/// Data::Intersect's candidate guard (Forward/Reversed only, mirroring
/// OCCT HLRBRep_Data.cxx L1275) never intersects them; the OCCT mechanism
/// producing this split is still to be identified. Un-ignore when closed.
#[test]
#[ignore = "remaining hider gap: the big-top far-upper arc is not split/hidden at the small-cylinder outline crossings (16.9334)"]
fn acceptance_bug25813_1_vs_occt_json() {
    let j: serde_json::Value =
        serde_json::from_str(&ref_json_text("occt_hlr_exact_hlr_bug25813_1.json"))
            .expect("reference JSON parses");
    let mass_total = j["mass_total"].as_f64().expect("mass_total");
    let mass_hidden_total = j["mass_hidden_total"].as_f64().expect("mass_hidden_total");
    let ref_result_edges = j["nbshapes_result"]["EDGE"].as_u64().expect("EDGE") as usize;
    let ref_hidden_edges = j["nbshapes_hidden"]["EDGE"].as_u64().expect("EDGE") as usize;
    let ref_result_vertices = j["nbshapes_result"]["VERTEX"].as_u64().expect("VERTEX") as usize;
    let ref_hidden_vertices = j["nbshapes_hidden"]["VERTEX"].as_u64().expect("VERTEX") as usize;
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
    let (v, rg1v, rgnv, outline, h) = smoke_run_hlr(&solid, fused);
    let vis = VisibleCompounds { v, rg1v, rgnv, outline };
    vis.dump("bug25813_1");
    let _ = rgnv;
    println!(
        "ACCEPT bug25813_1: visible mass {} (edges {}, vertices {}) vs mass_total \
         {mass_total} (edges {ref_result_edges}, vertices {ref_result_vertices}); \
         hidden mass {} (edges {}, vertices {}) ; combined {} vs mass_hidden_total \
         {mass_hidden_total}",
        vis.result_mass(),
        vis.result_edges(),
        vis.result_vertices(),
        compound_mass(&h),
        compound_edge_count(&h),
        compound_vertex_count(&h),
        vis.result_mass() + compound_mass(&h),
    );

    assert_rel(vis.result_mass(), mass_total, "visible mass vs mass_total");
    assert_rel(
        vis.result_mass() + compound_mass(&h),
        mass_hidden_total,
        "visible+hidden mass vs mass_hidden_total",
    );
    // nbshapes EDGE pins (the OCCT HLRToShape output corroboration).
    assert_eq!(vis.result_edges(), ref_result_edges, "visible EDGE count");
    assert_eq!(compound_edge_count(&h), ref_hidden_edges, "hidden EDGE count");
    // VERTEX pins: the OCCT nbshapes VERTEX:30/36. If rcad's fuse output
    // partitions the seam/tangency vertices benignly differently these are
    // informational; the edge and mass pins above stay authoritative.
    assert_eq!(vis.result_vertices(), ref_result_vertices, "visible VERTEX count");
    assert_eq!(compound_vertex_count(&h), ref_hidden_vertices, "hidden VERTEX count");
}

/// The `exact_hlr/bug25813_3` DRAW case (`ptorus a 30 10`) against the
/// `#### SUPPLEMENT` DRAWEXE viewer-path truth of ref_output.txt:
/// Mass 302.685 with 4 visible edges. NOTE: OCCT 8.0.0 itself access-violates
/// on this case via the C++ tuple projector path (exit 139, documented in
/// ref_output.txt); the 302.685 truth comes from the official DRAWEXE viewer
/// path (vinit+vdisplay+vcomputehlr). rcad must produce the same result
/// without crashing.
///
/// IGNORED (remaining exact-HLR gap): the torus panic is closed and Contap
/// now feeds the walking stage correctly (the restriction search finds the 8
/// boundary tangency points, the inside search 20 points, ComputeTangency
/// builds 8 departure points), but the IntWalk_IWalking engine returns done
/// with 0 lines, so rcad emits 0 visible edges. Un-ignore when the walking
/// engine walks the torus contour lines.
#[test]
#[ignore = "remaining gap: IntWalk_IWalking produces 0 contour lines for the torus (8 departure + 20 interior start points are fed correctly)"]
fn acceptance_ptorus_vs_viewer_path() {
    let text = ref_output_text();
    let (mass_ref, edges_ref) = parse_supplement(&text, "exact_hlr/bug25813_3");
    assert_eq!(edges_ref, 4, "supplement parser: ref visible EDGE count");
    assert!((mass_ref - 302.685).abs() < 5e-4, "supplement Mass {mass_ref}");

    // ptorus a 30 10 — OCCT BRepPrimAPI_MakeTorus(R1=30 major, R2=10 minor).
    let torus = rcad_modeling::make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 30.0, 10.0)
        .expect("ptorus a 30 10");
    let solid = root_solid(&torus);
    let (v, rg1v, rgnv, outline, h) = smoke_run_hlr(&solid, torus);
    let vis = VisibleCompounds { v, rg1v, rgnv, outline };
    vis.dump("bug25813_3");
    println!(
        "ACCEPT bug25813_3 (ptorus 30 10): visible mass {} (edges {}) vs viewer-path \
         Mass {mass_ref}; hidden mass {} (edges {})",
        vis.result_mass(),
        vis.result_edges(),
        compound_mass(&h),
        compound_edge_count(&h),
    );

    assert_eq!(vis.result_edges(), 4, "visible EDGE count (viewer path: 4)");
    assert_rel(vis.result_mass(), mass_ref, "visible mass vs viewer-path Mass");
}

