//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_WireVertex`
//! (`ShapeAnalysis_WireVertex.hxx` L17-130 + `.cxx` L1-338).
//!
//! Analyzes the status of the connection of consecutive vertices in a wire:
//! same vertex, same coordinates, close coordinates, one edge ending on the
//! other (start/end to re-limit), intersection or disjunction.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — `BRep_Tool::Pnt/Tolerance` and
//!    `ShapeAnalysis_Edge` read the TShape graph through
//!    `rcad_kernel::BRep` (the edge.rs bridge #1).
//! 2. `ShapeExtend_WireData` -> `shape_extend::wire_data::WireData` (the
//!    W1-3 translation).
//! 3. `NCollection_HArray1<int>` / `<double>` / `<gp_XYZ>` -> `Vec<i32>` /
//!    `Vec<f64>` / `Vec<DVec3>` (1-based OCCT index -> index - 1).  The
//!    OCCT null-array test (`myStat.IsNull()`) maps to `my_wire.is_none()`
//!    (the arrays are created exactly when the wire data is stored).
//! 4. `gp_XYZ` -> `DVec3`; `gp_Pnt` -> `DVec3`.
//! 5. `V1 == V2` (TopoDS_Shape operator== / IsEqual) -> same TShape and
//!    location AND same orientation (the local `shape_is_equal`).

use glam::DVec3;
use rcad_kernel::geom::CurveEval as _;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, TShape};

use super::curve::ShapeAnalysisCurve;
use super::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_extend::wire_data::WireData;

/// OCCT BRep_Tool::Pnt(V).
fn brep_tool_pnt(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(V).
fn brep_tool_tolerance(v: &Shape) -> f64 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        _ => 0.0,
    }
}

/// OCCT TopoDS_Shape::operator== (IsEqual): same TShape + location, and
/// same orientation.
fn shape_is_equal(a: &Shape, b: &Shape) -> bool {
    a.is_same(b) && a.orientation == b.orientation
}

/// OCCT ShapeAnalysis_WireVertex (hxx L23-130).
pub struct ShapeAnalysisWireVertex {
    /// OCCT `myDone`.
    my_done: bool,
    /// OCCT `myPreci`.
    my_preci: f64,
    /// OCCT `myWire` (bridge #2: the WireData value).
    my_wire: Option<WireData>,
    /// OCCT `myStat`.
    my_stat: Vec<i32>,
    /// OCCT `myPos`.
    my_pos: Vec<DVec3>,
    /// OCCT `myUPre`.
    my_u_pre: Vec<f64>,
    /// OCCT `myUFol`.
    my_u_fol: Vec<f64>,
}

impl Default for ShapeAnalysisWireVertex {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisWireVertex {
    /// OCCT ShapeAnalysis_WireVertex() (cxx L32-36).
    pub fn new() -> Self {
        ShapeAnalysisWireVertex {
            my_done: false,
            my_preci: CONFUSION,
            my_wire: None,
            my_stat: Vec::new(),
            my_pos: Vec::new(),
            my_u_pre: Vec::new(),
            my_u_fol: Vec::new(),
        }
    }

    /// OCCT Init(wire, preci) (cxx L40-43).
    pub fn init(&mut self, brep: &mut BRep, wire: &Shape, preci: f64) {
        let sbwd = WireData::new_from_wire(brep, wire, false, true);
        self.init_wire_data(sbwd, preci);
    }

    /// OCCT Init(sbwd, preci) (cxx L47-64) — the precision is ignored
    /// (the OCCT `/*preci*/` parameter).
    pub fn init_wire_data(&mut self, sbwd: WireData, _preci: f64) {
        let nb = sbwd.nb_edges();
        if nb == 0 {
            return;
        }
        self.my_done = false;
        self.my_stat = vec![0; nb as usize]; // myStat->Init(0)
        self.my_pos = vec![DVec3::ZERO; nb as usize];
        self.my_u_pre = vec![0.0; nb as usize]; // myUPre->Init(0.0)
        self.my_u_fol = vec![0.0; nb as usize]; // myUFol->Init(0.0)
        self.my_wire = Some(sbwd);
    }

    /// OCCT Load(wire) (cxx L68-71).
    pub fn load(&mut self, brep: &mut BRep, wire: &Shape) {
        let preci = self.my_preci;
        self.init(brep, wire, preci);
    }

    /// OCCT Load(sbwd) (cxx L75-78).
    pub fn load_wire_data(&mut self, sbwd: WireData) {
        let preci = self.my_preci;
        self.init_wire_data(sbwd, preci);
    }

    /// OCCT SetPrecision(preci) (cxx L82-86).
    pub fn set_precision(&mut self, preci: f64) {
        self.my_preci = preci;
        self.my_done = false;
    }

    /// OCCT Analyze() (cxx L90-170).
    pub fn analyze(&mut self, brep: &mut BRep) {
        if self.my_wire.is_none() {
            return;
        }
        self.my_done = true;
        //  Analyse des vertex qui se suivent
        let mut c1: Option<rcad_kernel::geom::Curve3> = None;
        let mut c2: Option<rcad_kernel::geom::Curve3> = None;
        // OCCT declares `double cf, cl, upre, ufol;` outside the loop; the
        // rcad zero initialisers stand in (the uses are guarded by the
        // curve-null check exactly as in OCCT).
        let mut cf = 0.0;
        let mut cl = 0.0;
        let mut upre = 0.0;
        let mut ufol = 0.0;
        let nb = self.my_stat.len() as i32;
        let mut stat;
        let ea = ShapeAnalysisEdge::new();
        for i in 1..=nb {
            stat = -1; // au depart

            let j = if i == nb { 1 } else { i + 1 };
            // OCCT L107-108: E1/E2 declared and fetched (unused by the OCCT
            // walk besides the fetch itself).
            let _e1 = self.my_wire.as_ref().unwrap().edge(i);
            let _e2 = self.my_wire.as_ref().unwrap().edge(j);
            let v1 = ea.last_vertex(brep, &self.my_wire.as_ref().unwrap().edge(i));
            let v2 = ea.first_vertex(brep, &self.my_wire.as_ref().unwrap().edge(j));
            let pv1 = brep_tool_pnt(&v1);
            let pv2 = brep_tool_pnt(&v2);
            let tol1 = brep_tool_tolerance(&v1);
            let tol2 = brep_tool_tolerance(&v2);
            ea.curve3d(
                brep,
                &self.my_wire.as_ref().unwrap().edge(i),
                &mut c1,
                &mut cf,
                &mut upre,
                true,
            );
            ea.curve3d(
                brep,
                &self.my_wire.as_ref().unwrap().edge(j),
                &mut c2,
                &mut ufol,
                &mut cl,
                true,
            );
            let (Some(c1), Some(c2)) = (c1.as_ref(), c2.as_ref()) else {
                // c1.IsNull() || c2.IsNull() — on ne peut rien faire ...
                continue;
            };
            let p1 = c1.point_at(upre);
            let p2 = c2.point_at(ufol);

            //   Est-ce que le jeu de vertex convient ? (meme si V1 == V2, on verifie)
            let d1 = pv1.distance(p1);
            let d2 = pv2.distance(p2);
            let dd = pv1.distance(pv2);
            if d1 <= tol1 && d2 <= tol2 && dd <= (tol1 + tol2) {
                stat = 1;
            } else if d1 <= self.my_preci && d2 <= self.my_preci && dd <= self.my_preci {
                stat = 2;
            }
            self.my_stat[(i - 1) as usize] = -1; // par defaut
            if stat > 0 {
                if shape_is_equal(&v1, &v2) {
                    stat = 0;
                }
            }
            if stat >= 0 {
                self.my_stat[(i - 1) as usize] = stat;
                continue;
            }
            //    Restent les autres cas !

            //    Une edge se termine sur l autre : il faudra simplement relimiter
            //    Projection calculee sur une demi-edge (pour eviter les pbs de couture)
            let mut pj1 = DVec3::ZERO;
            let mut pj2 = DVec3::ZERO;
            let mut u1 = 0.0;
            let mut u2 = 0.0;
            let sac = ShapeAnalysisCurve;
            let dj1 = sac.project_cf_cl(
                c1,
                p2,
                self.my_preci,
                &mut pj1,
                &mut u1,
                (cf + upre) / 2.,
                upre,
                true,
            );
            let dj2 = sac.project_cf_cl(
                c2,
                p1,
                self.my_preci,
                &mut pj2,
                &mut u2,
                ufol,
                (ufol + cl) / 2.,
                true,
            );
            if dj1 <= self.my_preci {
                self.set_start(i, pj1, u1);
                continue;
            } else if dj2 <= self.my_preci {
                self.set_end(i, pj2, u2);
                continue;
            }

            //    Restent a verifier les intersections et prolongations !
        }
    }

    /// OCCT SetSameVertex(num) (cxx L174-177).
    pub fn set_same_vertex(&mut self, num: i32) {
        self.my_stat[(num - 1) as usize] = 0;
    }

    /// OCCT SetSameCoords(num) (cxx L181-184).
    pub fn set_same_coords(&mut self, num: i32) {
        self.my_stat[(num - 1) as usize] = 1;
    }

    /// OCCT SetClose(num) (cxx L188-191).
    pub fn set_close(&mut self, num: i32) {
        self.my_stat[(num - 1) as usize] = 2;
    }

    /// OCCT SetEnd(num, pos, ufol) (cxx L195-200).
    pub fn set_end(&mut self, num: i32, pos: DVec3, ufol: f64) {
        self.my_stat[(num - 1) as usize] = 3;
        self.my_pos[(num - 1) as usize] = pos;
        self.my_u_fol[(num - 1) as usize] = ufol;
    }

    /// OCCT SetStart(num, pos, upre) (cxx L204-209).
    pub fn set_start(&mut self, num: i32, pos: DVec3, upre: f64) {
        self.my_stat[(num - 1) as usize] = 4;
        self.my_pos[(num - 1) as usize] = pos;
        self.my_u_fol[(num - 1) as usize] = upre;
    }

    /// OCCT SetInters(num, pos, upre, ufol) (cxx L213-222).
    pub fn set_inters(&mut self, num: i32, pos: DVec3, upre: f64, ufol: f64) {
        self.my_stat[(num - 1) as usize] = 5;
        self.my_pos[(num - 1) as usize] = pos;
        self.my_u_pre[(num - 1) as usize] = upre;
        self.my_u_fol[(num - 1) as usize] = ufol;
    }

    /// OCCT SetDisjoined(num) (cxx L226-229).
    pub fn set_disjoined(&mut self, num: i32) {
        self.my_stat[(num - 1) as usize] = -1;
    }

    /// OCCT IsDone() (cxx L233-236).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT Precision() (cxx L240-243).
    pub fn precision(&self) -> f64 {
        self.my_preci
    }

    /// OCCT NbEdges() (cxx L247-250).
    pub fn nb_edges(&self) -> i32 {
        self.my_wire.as_ref().unwrap().nb_edges()
    }

    /// OCCT WireData() (cxx L254-257) — the stored handle; the rcad value
    /// model returns a clone (bridge #2).
    pub fn wire_data(&self) -> Option<&WireData> {
        self.my_wire.as_ref()
    }

    /// OCCT WireData() — the mutable form (the OCCT handle is shared and
    /// mutable through it; the rcad value model needs the explicit
    /// accessor, bridge #2).
    pub fn wire_data_mut(&mut self) -> Option<&mut WireData> {
        self.my_wire.as_mut()
    }

    /// OCCT Status(num) (cxx L261-264).
    pub fn status(&self, num: i32) -> i32 {
        self.my_stat[(num - 1) as usize]
    }

    /// OCCT Position(num) (cxx L268-271).
    pub fn position(&self, num: i32) -> DVec3 {
        self.my_pos[(num - 1) as usize]
    }

    /// OCCT UPrevious(num) (cxx L276-279) — szv#4:S4163:12Mar99 was bug:
    /// returned Integer.
    pub fn uprevious(&self, num: i32) -> f64 {
        self.my_u_pre[(num - 1) as usize]
    }

    /// OCCT UFollowing(num) (cxx L284-287) — szv#4:S4163:12Mar99 was bug:
    /// returned Integer.
    pub fn ufollowing(&self, num: i32) -> f64 {
        self.my_u_fol[(num - 1) as usize]
    }

    /// OCCT Data(num, pos, upre, ufol) (cxx L291-297).
    pub fn data(&self, num: i32, pos: &mut DVec3, upre: &mut f64, ufol: &mut f64) -> i32 {
        *pos = self.my_pos[(num - 1) as usize];
        *upre = self.my_u_pre[(num - 1) as usize];
        *ufol = self.my_u_fol[(num - 1) as usize];
        self.my_stat[(num - 1) as usize]
    }

    /// OCCT NextStatus(stat, num) (cxx L301-316).
    pub fn next_status(&self, stat: i32, num: i32) -> i32 {
        // szv#4:S4163:12Mar99 optimized
        if !self.my_wire.is_none() {
            let nb = self.my_stat.len() as i32;
            for i in (num + 1)..=nb {
                if self.my_stat[(i - 1) as usize] == stat {
                    return i;
                }
            }
        }
        0
    }

    /// OCCT NextCriter(crit, num) (cxx L320-338).
    pub fn next_criter(&self, crit: i32, num: i32) -> i32 {
        // szv#4:S4163:12Mar99 optimized
        if !self.my_wire.is_none() {
            let nb = self.my_stat.len() as i32;
            for i in (num + 1)..=nb {
                let stat = self.my_stat[(i - 1) as usize];
                if (crit == -1 && stat < 0)
                    || (crit == 0 && stat == 0)
                    || (crit == 1 && stat > 0)
                    || (crit == 2 && (stat >= 0 && stat <= 2))
                    || (crit == 3 && (stat == 1 || stat == 2))
                    || (crit == 4 && stat > 2)
                {
                    return i;
                }
            }
        }
        0
    }
}
