//! ChFiDS gap-fill — methods of `ChFiDS_Spine`, `ChFiDS_FilSpine` and
//! `ChFiDS_ElSpine` that were missing from `chfi_ds.rs`, translated 1:1
//! from OCCT TKFillet/ChFiDS.  Kept in a dedicated file because
//! `chfi_ds.rs` exceeds the 2000-line guideline.

use glam::DVec3;
use rcad_kernel::topo::topods::{Orientation, Shape};
use rcad_kernel::geom::CurveEval as _;

use super::chfi_ds::{ChFiDSElSpine, ChFiDSFilSpine, ChFiDSSpine, elclib_in_period};

// =========================================================================
// OCCT ChFiDS_Spine — missing method translations
// =========================================================================

impl ChFiDSSpine {
    /// OCCT ChFiDS_Spine.lxx L117-120 — OffsetEdges(I) (1-based).
    pub fn offset_edges(&self, i: usize) -> &Shape {
        &self.offsetspine[i - 1]
    }

    /// OCCT ChFiDS_Spine.cxx L89-93 — AppendElSpine(Els).
    pub fn append_el_spine(&mut self, els: ChFiDSElSpine) {
        self.elspines.push(els);
    }

    /// OCCT ChFiDS_Spine.cxx L96-99 — AppendOffsetElSpine(Els).
    pub fn append_offset_el_spine(&mut self, els: ChFiDSElSpine) {
        self.offset_elspines.push(els);
    }

    /// OCCT ChFiDS_Spine.cxx L103-106 — ElSpine(const TopoDS_Edge& E).
    pub fn el_spine_of_edge(&self, e: &Shape) -> Option<&ChFiDSElSpine> {
        self.el_spine_of_index(self.index_of_edge(e))
    }

    /// OCCT ChFiDS_Spine.cxx L108-116 — ElSpine(const int IE) (1-based).
    pub fn el_spine_of_index(&self, ie: usize) -> Option<&ChFiDSElSpine> {
        let mut wmil = 0.5 * (self.first_parameter_of(ie) + self.last_parameter_of(ie));
        if self.is_periodic() {
            wmil = elclib_in_period(wmil, self.first_parameter(), self.last_parameter());
        }
        self.el_spine_of_param(wmil)
    }

    /// OCCT ChFiDS_Spine.cxx L118-139 — ElSpine(const double W).
    pub fn el_spine_of_param(&self, w: f64) -> Option<&ChFiDSElSpine> {
        if self.elspines.len() == 1 {
            return self.elspines.first();
        }
        for cur in &self.elspines {
            let uf = cur.first_parameter();
            let ul = cur.last_parameter();
            if uf <= w && w <= ul {
                return Some(cur);
            }
        }
        None
    }

    /// OCCT ChFiDS_Spine.cxx L143-146 — ChangeElSpines().
    pub fn change_el_spines(&mut self) -> &mut Vec<ChFiDSElSpine> {
        &mut self.elspines
    }

    /// OCCT ChFiDS_Spine.cxx L150-153 — ChangeOffsetElSpines().
    pub fn change_offset_el_spines(&mut self) -> &mut Vec<ChFiDSElSpine> {
        &mut self.offset_elspines
    }

    /// OCCT ChFiDS_Spine.cxx L345-348 — Resolution(R3d).
    pub fn resolution(&self, r3d: f64) -> f64 {
        r3d
    }

    /// OCCT ChFiDS_Spine.cxx L396-399 — HasFirstTgt().
    pub fn has_first_tgt(&self) -> bool {
        self.hasfirsttgt
    }

    /// OCCT ChFiDS_Spine.cxx L403-406 — HasLastTgt().
    pub fn has_last_tgt(&self) -> bool {
        self.haslasttgt
    }

    /// OCCT ChFiDS_Spine.cxx L496-499 — UnsetReference().
    pub fn unset_reference(&mut self) {
        self.hasref = false;
    }

    /// OCCT ChFiDS_Spine.cxx L743-746 — D0(AbsC, P).
    pub fn d0(&mut self, absc: f64) -> DVec3 {
        self.value_at(absc)
    }

    /// OCCT ChFiDS_Spine.cxx L798-851 — D2(AbsC, P, V1, V2): point, tangent
    /// and second derivative (normalized, orientation-adjusted) on the
    /// composite spine.
    pub fn d2(&mut self, absc: f64) -> (DVec3, DVec3, DVec3) {
        let mut l = absc;
        let index = self.prepare(&mut l);

        if index == -1 {
            let p = self.firstori + self.firsttgt * l;
            return (p, self.firsttgt, DVec3::ZERO);
        } else if index as i32 == self.abscissa.as_ref().map_or(0, |a| a.len()) as i32 + 1 {
            let p = self.lastori + self.lasttgt * l;
            return (p, self.lasttgt, DVec3::ZERO);
        }
        let index = index as usize;
        self.indexofcurve = index as i32;
        let e = &self.spine[index - 1];
        let ed = e.as_edge().expect("not an edge");
        let curve = ed.curve.as_ref().expect("edge curve").clone();
        let t = l / self.length_of(index);
        let (cf, cl) = (ed.range[0], ed.range[1]);
        let uapp = (1.0 - t) * cf + t * cl;
        let u = rcad_kernel::base::gcpnts::abscissa_point::abscissa_point_parameter(
            &curve, cf, cl, l, uapp,
        );
        let p = curve.point_at(u);
        let mut v1 = curve.derivative_at(u);
        let mut v2 = curve.derivative2_at(u);
        // OCCT: N1 = V1.SquareMagnitude()
        let mut n1 = v1.dot(v1);
        // OCCT: D2 = -(V1.Dot(V2)) * (1./N1) * (1./N1)
        let d2c = -(v1.dot(v2)) * (1.0 / n1) * (1.0 / n1);
        v2 = v2 * (1.0 / n1);
        n1 = n1.sqrt();
        let va = v1 * d2c;
        v2 = v2 + va;
        let mut d1 = 1.0 / n1;
        if e.orientation == Orientation::Reversed {
            d1 = -d1;
        }
        v1 = v1 * d1;
        (p, v1, v2)
    }

    /// OCCT ChFiDS_Spine.cxx L855-862 — SetCurrent(Index).
    /// The OCCT myCurve.Initialize(TopoDS::Edge(...)) re-initialization is a
    /// pending boundary: rcad tracks the current elementary spine by
    /// `indexofcurve` alone (BRepAdaptor_Curve pending translation).
    pub fn set_current(&mut self, index: i32) {
        if index != self.indexofcurve {
            self.indexofcurve = index;
        }
    }

    /// OCCT ChFiDS_Spine.lxx L152-155 — CurrentIndexOfElementarySpine().
    pub fn current_index_of_elementary_spine(&self) -> i32 {
        self.indexofcurve
    }

    /// OCCT ChFiDS_Spine.lxx L159-163 — Mode().
    pub fn mode(&self) -> super::chfi_ds::ChFiDS_ChamfMode {
        self.my_mode
    }

    /// OCCT ChFiDS_Spine.lxx L167-170 — GetTolesp().
    pub fn get_tolesp(&self) -> f64 {
        self.tolesp
    }
}

// =========================================================================
// OCCT ChFiDS_FilSpine — missing method translations
// =========================================================================

impl ChFiDSFilSpine {
    /// OCCT ChFiDS_FilSpine.cxx L96-119 — UnSetRadius(const TopoDS_Edge& E)
    /// (1-based IE).
    pub fn unset_radius_on_edge(&mut self, e: &Shape) {
        self.base.splitdone = false;
        let ie = self.base.index_of_edge(e);

        let uf = self.base.first_parameter_of(ie);
        let ul = self.base.last_parameter_of(ie);
        let mut ifirst = 0usize;
        let mut ilast = 0usize;
        for i in 1..=self.parandrad.len() {
            if (self.parandrad[i - 1].x - uf).abs() <= f64::MIN_POSITIVE {
                ifirst = i;
            }
            if (self.parandrad[i - 1].x - ul).abs() <= f64::MIN_POSITIVE {
                ilast = i;
            }
        }
        if ifirst != 0 && ilast != 0 {
            // OCCT: parandrad.Remove(ifirst, ilast) — the inclusive range.
            self.parandrad.drain((ifirst - 1)..ilast);
        }
    }

    /// OCCT ChFiDS_FilSpine.cxx L220-231 — UnSetRadius(const TopoDS_Vertex& V).
    pub fn unset_radius_at_vertex(&mut self, v: &Shape) {
        let npar = self.base.absc_of_vertex(v);
        for i in 1..=self.parandrad.len() {
            if self.parandrad[i - 1].x == npar {
                self.parandrad.remove(i - 1);
                break;
            }
        }
    }
}

// =========================================================================
// OCCT ChFiDS_ElSpine — missing accessor translations.
// (SetFirstPointAndTgt / SetLastPointAndTgt / the default constructor live
// in chfi3d.rs — pre-existing.)
// =========================================================================

impl ChFiDSElSpine {
    /// OCCT ChFiDS_ElSpine.cxx L76-79 — FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.firstparam
    }

    /// OCCT ChFiDS_ElSpine.cxx L83-86 — LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.lastparam
    }

    /// OCCT ChFiDS_ElSpine.cxx L269-273 — FirstPointAndTgt(P, T).
    pub fn first_point_and_tgt(&self) -> (DVec3, DVec3) {
        (self.firstpnt, self.firsttgt)
    }

    /// OCCT ChFiDS_ElSpine.cxx L277-281 — LastPointAndTgt(P, T).
    pub fn last_point_and_tgt(&self) -> (DVec3, DVec3) {
        (self.lastpnt, self.lasttgt)
    }
}
