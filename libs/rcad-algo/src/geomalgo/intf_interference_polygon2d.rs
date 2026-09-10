//! OCCT Intf_InterferencePolygon2d (TKGeomAlgo Intf package) — interference
//! between two 2D polygons of segments (or the auto-interference of one).
//!
//! 1:1 translation of `Intf_InterferencePolygon2d.cxx` (L17-819) over the
//! [`IntfPolygon2d`] trait (the OCCT Intf_Polygon2d abstract base).  The
//! consumer instantiation is IntCurve_ThePolygon2d ==
//! IntCurve_Polygon2dGen ([`super::int_curve_generics::Polygon2dGen`]).
//! Direct consumers: `IntCurve_IntPolyPolyGen.gxx` (2a-3 item 3).

use glam::{DVec2, DVec3};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::direct_polynomial_roots::epsilon;

use super::intf::{IntfPIType, IntfSectionPoint};
use super::intf_interference::{Interference, IntfPolygon2d};
use super::intf_tangent_zone::TangentZone;

/// OCCT cxx L31: angular precision (sinus) below which two right segments
/// are considered as having a potential zone of tangency
/// (Precision::Angular()).
const PRCANG: f64 = 1.0e-12;

/// OCCT gxx use: `Epsilon(value)` (Standard_Real.hxx L242-247).
fn epsilon_of(value: f64) -> f64 {
    epsilon(value)
}

/// OCCT Intf_InterferencePolygon2d — the 2D polygon interference engine.
pub struct InterferencePolygon2d {
    /// OCCT base subobject Intf_Interference.
    pub interf: Interference,
    /// OCCT oClos.
    o_clos: bool,
    /// OCCT tClos.
    t_clos: bool,
    /// OCCT nbso.
    nbso: i32,
}

impl InterferencePolygon2d {
    /// OCCT Intf_InterferencePolygon2d() (cxx L36-42).
    pub fn new() -> Self {
        InterferencePolygon2d {
            interf: Interference::with_self(false),
            o_clos: false,
            t_clos: false,
            nbso: 0,
        }
    }

    /// OCCT Intf_InterferencePolygon2d(Obje1, Obje2) (cxx L49-69) —
    /// interference between two polygons.
    pub fn new_polygon_polygon<O1: IntfPolygon2d, O2: IntfPolygon2d>(
        obje1: &O1,
        obje2: &O2,
    ) -> Self {
        let mut r = InterferencePolygon2d::new();
        if !obje1.bounding().is_out_box(obje2.bounding()) {
            let mut tolerance =
                obje1.deflection_over_estimation() + obje2.deflection_over_estimation();
            if tolerance == 0.0 {
                tolerance = epsilon_of(1000.0);
            }
            r.interf.set_tolerance(tolerance);
            r.nbso = obje1.nb_segments();
            r.o_clos = obje1.closed();
            r.t_clos = obje2.closed();
            r.interference_two(obje1, obje2);
            r.clean();
        }
        r
    }

    /// OCCT Intf_InterferencePolygon2d(Obje) (cxx L76-91) — auto
    /// interference of one polygon.
    pub fn new_polygon_auto<O: IntfPolygon2d>(obje: &O) -> Self {
        let mut r = InterferencePolygon2d::new();
        r.interf = Interference::with_self(true);
        let mut tolerance = obje.deflection_over_estimation() * 2.0;
        if tolerance == 0.0 {
            tolerance = epsilon_of(1000.0);
        }
        r.interf.set_tolerance(tolerance);
        r.o_clos = obje.closed();
        r.t_clos = r.o_clos;
        r.interference_auto(obje);
        r.clean();
        r
    }

    /// OCCT Perform(Obje1, Obje2) (cxx L95-111).
    pub fn perform_polygon_polygon<O1: IntfPolygon2d, O2: IntfPolygon2d>(
        &mut self,
        obje1: &O1,
        obje2: &O2,
    ) {
        self.interf.self_interference(false);
        if !obje1.bounding().is_out_box(obje2.bounding()) {
            let mut tolerance =
                obje1.deflection_over_estimation() + obje2.deflection_over_estimation();
            if tolerance == 0.0 {
                tolerance = epsilon_of(1000.0);
            }
            self.interf.set_tolerance(tolerance);
            self.nbso = obje1.nb_segments();
            self.o_clos = obje1.closed();
            self.t_clos = obje2.closed();
            self.interference_two(obje1, obje2);
            self.clean();
        }
    }

    /// OCCT Perform(Obje) (cxx L115-127).
    pub fn perform_polygon_auto<O: IntfPolygon2d>(&mut self, obje: &O) {
        self.interf.self_interference(true);
        let mut tolerance = obje.deflection_over_estimation() * 2.0;
        if tolerance == 0.0 {
            tolerance = epsilon_of(1000.0);
        }
        self.interf.set_tolerance(tolerance);
        self.o_clos = obje.closed();
        self.t_clos = self.o_clos;
        self.interference_auto(obje);
        self.clean();
    }

    /// OCCT Pnt2dValue(Index) (cxx L134-137) — 1-based.
    pub fn pnt_2d_value(&self, index: usize) -> DVec2 {
        let p = self.interf.pnt_value(index).pnt();
        DVec2::new(p.x, p.y)
    }

    /// OCCT Interference(Obje1, Obje2) (cxx L141-174).
    fn interference_two<O1: IntfPolygon2d, O2: IntfPolygon2d>(&mut self, obje1: &O1, obje2: &O2) {
        let mut b_so = BndBox2d::new();
        let mut b_st = BndBox2d::new();

        let n1 = self.nbso;
        let n2 = obje2.nb_segments();
        let d1 = obje1.deflection_over_estimation();
        let d2 = obje2.deflection_over_estimation();

        for i_obje1 in 1..=n1 {
            b_so.set_void();
            let (p1b, p1e) = obje1.segment(i_obje1);
            b_so.add_point(p1b);
            b_so.add_point(p1e);
            b_so.enlarge(d1);
            if !obje2.bounding().is_out_box(&b_so) {
                for i_obje2 in 1..=n2 {
                    b_st.set_void();
                    let (p2b, p2e) = obje2.segment(i_obje2);
                    b_st.add_point(p2b);
                    b_st.add_point(p2e);
                    b_st.enlarge(d2);
                    if !b_so.is_out_box(&b_st) {
                        self.intersect(i_obje1, i_obje2, p1b, p1e, p2b, p2e);
                    }
                }
            }
        }
    }

    /// OCCT Interference(Obje) (cxx L178-210) — the auto-interference pass
    /// over the upper triangle of the segment-pair matrix.
    fn interference_auto<O: IntfPolygon2d>(&mut self, obje: &O) {
        let mut b_so = BndBox2d::new();
        let mut b_st = BndBox2d::new();

        let n = obje.nb_segments();
        let d = obje.deflection_over_estimation();

        for i_obje1 in 1..=n {
            b_so.set_void();
            let (p1b, p1e) = obje.segment(i_obje1);
            b_so.add_point(p1b);
            b_so.add_point(p1e);
            b_so.enlarge(d);
            if !obje.bounding().is_out_box(&b_so) {
                for i_obje2 in i_obje1 + 1..=n {
                    b_st.set_void();
                    let (p2b, p2e) = obje.segment(i_obje2);
                    b_st.add_point(p2b);
                    b_st.add_point(p2e);
                    b_st.enlarge(d);
                    if !b_so.is_out_box(&b_st) {
                        self.intersect(i_obje1, i_obje2, p1b, p1e, p2b, p2e);
                    }
                }
            }
        }
    }

    /// OCCT Clean() (cxx L214-305): the tangent zones that concern only one
    /// couple of segments are conserved as section points when the angle
    /// between the segments is less than <PRCANG> and when there is no real
    /// EDGE/EDGE point; then the section points located inside a tangency
    /// zone are removed from the list.
    fn clean(&mut self) {
        let nb_it = self.interf.nb_tangent_zones();
        let mut decal = 0usize;
        // Only1Seg is declared outside the zone loop in the OCCT source and
        // is not reinitialised per zone (kept verbatim).
        let mut only1seg = false;

        for ltz in 1..=nb_it {
            let mut tsp = 0usize;
            let mut tsps = 0usize;
            let (pr1mi, pr1ma) = self.interf.my_t_zones_mut()[ltz - 1 - decal].param_on_first();
            let delta1 = pr1ma - pr1mi;
            let (pr2mi, pr2ma) = self.interf.my_t_zones_mut()[ltz - 1 - decal].param_on_second();
            let delta2 = pr2ma - pr2mi;
            if delta1 < 1.0 && delta2 < 1.0 {
                only1seg = true;
            }
            if delta1 == 0.0 || delta2 == 0.0 {
                only1seg = true;
            }

            let nb_points = self.interf.my_t_zones_mut()[ltz - 1 - decal].number_of_points();
            for lpi in 1..=nb_points {
                let pi1 = self.interf.my_t_zones_mut()[ltz - 1 - decal].get_point(lpi);
                if pi1.incidence() <= PRCANG {
                    tsp = 0;
                    tsps = 0;
                    break;
                }
                let (dim1, addr1, _par) = pi1.info_first();
                let (dim2, addr2, _par) = pi1.info_second();
                let _ = (addr1, addr2);
                if dim1 == IntfPIType::Edge && dim2 == IntfPIType::Edge {
                    tsps = 0;
                    if tsp > 0 {
                        tsp = 0;
                        only1seg = false;
                        break;
                    }
                    tsp = lpi;
                } else if dim1 != IntfPIType::External && dim2 != IntfPIType::External {
                    tsps = lpi;
                }
            }
            if tsp > 0 {
                let sp = self.interf.my_t_zones_mut()[ltz - 1 - decal].get_point(tsp);
                self.interf.my_s_poins_mut().push(sp);
                self.interf.my_t_zones_mut().remove(ltz - 1 - decal);
                decal += 1;
            } else if only1seg && tsps != 0 {
                let sp = self.interf.my_t_zones_mut()[ltz - 1 - decal].get_point(tsps);
                self.interf.my_s_poins_mut().push(sp);
                self.interf.my_t_zones_mut().remove(ltz - 1 - decal);
                decal += 1;
            }
        }

        // The points of intersection located in the tangency zone are
        // removed from the list.
        let nb_it = self.interf.nb_section_points();
        let mut decal = 0usize;

        for lpi in 1..=nb_it {
            let mut ltz = 1;
            while ltz <= self.interf.nb_tangent_zones() {
                let sp = self.interf.pnt_value(lpi - decal).clone();
                let contains = self.interf.my_t_zones_mut()[ltz - 1].range_contains(&sp);
                if contains {
                    self.interf.my_s_poins_mut().remove(lpi - 1 - decal);
                    decal += 1;
                    break;
                }
                ltz += 1;
            }
        }
    }
}

impl InterferencePolygon2d {
    /// OCCT Intersect(iObje1, iObje2, BegO, EndO, BegT, EndT)
    /// (cxx L309-819) — the segment x segment interference with vertex,
    /// edge and tangency-zone classification.
    #[allow(clippy::too_many_arguments)]
    fn intersect(
        &mut self,
        i_obje1: i32,
        i_obje2: i32,
        beg_o: DVec2,
        end_o: DVec2,
        beg_t: DVec2,
        end_t: DVec2,
    ) {
        if self.interf.self_intf() && (i_obje1 - i_obje2).abs() <= 1 {
            return; //-- Ajout du 15 jan 98
        }

        let mut nbpi = 0usize;
        // OCCT uses 1-based slots of double parO[8] / parT[8].
        let mut par_o = [0.0f64; 9];
        let mut par_t = [0.0f64; 9];
        let mut the_pi: Vec<IntfSectionPoint> = Vec::new();
        let seg_t = end_t - beg_t;
        let seg_o = end_o - beg_o;

        // If the length of segment is zero, nothing is done.
        let lg_t = seg_t.dot(seg_t).sqrt();
        if lg_t <= 0.0 {
            return;
        }
        let lg_o = seg_o.dot(seg_o).sqrt();
        if lg_o <= 0.0 {
            return;
        }

        // Direction of parsing of segments.
        let sig_ps = if seg_o.dot(seg_t) > 0.0 { 1.0 } else { -1.0 };

        // Precision of calculation.
        let floatgap = epsilon_of(lg_o + lg_t);

        // Angle between two straight lines and radius of interference.
        let sin_teta = (cross_magnitude(seg_o, seg_t) / lg_o) / lg_t;
        let mut ray_intf = 0.0;
        if sin_teta > 0.0 {
            ray_intf = self.interf.get_tolerance() / sin_teta;
        }

        let tolerance = self.interf.get_tolerance();
        let d_b_o_t = cross_scalar(beg_o - beg_t, seg_t) / lg_t;
        let d_b_o_b_t = beg_o.distance(beg_t);
        let d_b_o_e_t = beg_o.distance(end_t);
        if d_b_o_t.abs() <= tolerance {
            if d_b_o_b_t <= tolerance {
                nbpi += 1;
                par_o[nbpi] = 0.0;
                par_t[nbpi] = 0.0;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(beg_o),
                    IntfPIType::Vertex,
                    0,
                    i_obje1,
                    0.0,
                    IntfPIType::Vertex,
                    0,
                    i_obje2,
                    0.0,
                    sin_teta,
                ));
            }
            if d_b_o_e_t <= tolerance {
                nbpi += 1;
                par_o[nbpi] = 0.0;
                par_t[nbpi] = 1.0;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(beg_o),
                    IntfPIType::Vertex,
                    0,
                    i_obje1,
                    0.0,
                    IntfPIType::Vertex,
                    0,
                    i_obje2 + 1,
                    0.0,
                    sin_teta,
                ));
            }
            if d_b_o_b_t > tolerance && d_b_o_e_t > tolerance && d_b_o_b_t + d_b_o_e_t <= (lg_t + tolerance)
            {
                nbpi += 1;
                par_o[nbpi] = 0.0;
                par_t[nbpi] = d_b_o_b_t / lg_t;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(beg_o),
                    IntfPIType::Vertex,
                    0,
                    i_obje1,
                    0.0,
                    IntfPIType::Edge,
                    0,
                    i_obje2,
                    par_t[nbpi],
                    sin_teta,
                ));
            }
        }

        // Interference <endO> <segT>.
        let d_e_o_t = cross_scalar(end_o - beg_t, seg_t) / lg_t;
        let d_e_o_b_t = end_o.distance(beg_t);
        let d_e_o_e_t = end_o.distance(end_t);
        if d_e_o_t.abs() <= tolerance {
            if d_e_o_b_t <= tolerance {
                nbpi += 1;
                par_o[nbpi] = 1.0;
                par_t[nbpi] = 0.0;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(end_o),
                    IntfPIType::Vertex,
                    0,
                    i_obje1 + 1,
                    0.0,
                    IntfPIType::Vertex,
                    0,
                    i_obje2,
                    0.0,
                    sin_teta,
                ));
            }
            if d_e_o_e_t <= tolerance {
                nbpi += 1;
                par_o[nbpi] = 1.0;
                par_t[nbpi] = 1.0;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(end_o),
                    IntfPIType::Vertex,
                    0,
                    i_obje1 + 1,
                    0.0,
                    IntfPIType::Vertex,
                    0,
                    i_obje2 + 1,
                    0.0,
                    sin_teta,
                ));
            }
            if d_e_o_b_t > tolerance && d_e_o_e_t > tolerance && d_e_o_b_t + d_e_o_e_t <= (lg_t + tolerance)
            {
                nbpi += 1;
                par_o[nbpi] = 1.0;
                par_t[nbpi] = d_e_o_b_t / lg_t;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(end_o),
                    IntfPIType::Vertex,
                    0,
                    i_obje1 + 1,
                    0.0,
                    IntfPIType::Edge,
                    0,
                    i_obje2,
                    par_t[nbpi],
                    sin_teta,
                ));
            }
        }

        // Interference <begT> <segO>.
        let d_b_t_o = cross_scalar(beg_t - beg_o, seg_o) / lg_o;
        if d_b_t_o.abs() <= tolerance {
            if d_b_o_b_t > tolerance && d_e_o_b_t > tolerance && d_b_o_b_t + d_e_o_b_t <= (lg_o + tolerance)
            {
                nbpi += 1;
                par_o[nbpi] = d_b_o_b_t / lg_o;
                par_t[nbpi] = 0.0;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(beg_t),
                    IntfPIType::Edge,
                    0,
                    i_obje1,
                    par_o[nbpi],
                    IntfPIType::Vertex,
                    0,
                    i_obje2,
                    0.0,
                    sin_teta,
                ));
            }
        }

        // Interference <endT> <segO>.
        let d_e_t_o = cross_scalar(end_t - beg_o, seg_o) / lg_o;
        if d_e_t_o.abs() <= tolerance {
            if d_b_o_e_t > tolerance && d_e_o_e_t > tolerance && d_b_o_e_t + d_e_o_e_t <= (lg_o + tolerance)
            {
                nbpi += 1;
                par_o[nbpi] = d_b_o_e_t / lg_o;
                par_t[nbpi] = 1.0;
                the_pi.push(IntfSectionPoint::new(
                    pnt2d_to_pnt(end_t),
                    IntfPIType::Edge,
                    0,
                    i_obje1,
                    par_o[nbpi],
                    IntfPIType::Vertex,
                    0,
                    i_obje2 + 1,
                    0.0,
                    sin_teta,
                ));
            }
        }

        let mut edge_sp = false;
        let mut par_osp = 0.0;
        let mut par_tsp = 0.0;

        if (d_b_o_t - d_e_o_t).abs() > floatgap && (d_b_t_o - d_e_t_o).abs() > floatgap {
            par_osp = d_b_o_t / (d_b_o_t - d_e_o_t);
            par_tsp = d_b_t_o / (d_b_t_o - d_e_t_o);
            if d_b_o_t * d_e_o_t <= 0.0 && d_b_t_o * d_e_t_o <= 0.0 {
                edge_sp = true;
            } else if nbpi == 0 {
                return;
            }

            // If there is no interference it is necessary to take the points
            // segment by segment.
            if nbpi == 0 && sin_teta > PRCANG {
                nbpi += 1;
                par_o[nbpi] = par_osp;
                par_t[nbpi] = par_tsp;
                the_pi.push(IntfSectionPoint::new(
                    DVec3::new(
                        beg_o.x + (seg_o.x * par_osp),
                        beg_o.y + (seg_o.y * par_osp),
                        0.0,
                    ),
                    IntfPIType::Edge,
                    0,
                    i_obje1,
                    par_osp,
                    IntfPIType::Edge,
                    0,
                    i_obje2,
                    par_tsp,
                    sin_teta,
                ));
            }
            // Otherwise it is required to check if there is no other.
            else if ray_intf >= tolerance {
                let delta_o = ray_intf / lg_o;
                let delta_t = ray_intf / lg_t;
                let mut x;
                let mut y;
                let mut par_odeb = par_osp - delta_o;
                let mut par_ofin = par_osp + delta_o;
                let mut par_tdeb = par_tsp - sig_ps * delta_t;
                let mut par_tfin = par_tsp + sig_ps * delta_t;
                if nbpi == 0 {
                    par_o[1] = par_odeb;
                    par_o[2] = par_ofin;
                    par_t[1] = par_tdeb;
                    par_t[2] = par_tfin;
                    while nbpi < 2 {
                        nbpi += 1;
                        x = beg_o.x + (seg_o.x * par_o[nbpi]);
                        y = beg_o.y + (seg_o.y * par_o[nbpi]);
                        the_pi.push(IntfSectionPoint::new(
                            DVec3::new(x, y, 0.0),
                            IntfPIType::External,
                            0,
                            i_obje1,
                            par_o[nbpi],
                            IntfPIType::External,
                            0,
                            i_obje2,
                            par_t[nbpi],
                            sin_teta,
                        ));
                    }
                } else {
                    // nbpi>0
                    if nbpi == 1 {
                        let mut ok = true;
                        if 0.0 < par_odeb && par_odeb < 1.0 && 0.0 < par_tdeb && par_tdeb < 1.0 {
                            par_o[nbpi + 1] = par_odeb;
                            par_t[nbpi + 1] = par_tdeb;
                        } else if 0.0 < par_ofin && par_ofin < 1.0 && 0.0 < par_tfin && par_tfin < 1.0
                        {
                            par_o[nbpi + 1] = par_ofin;
                            par_t[nbpi + 1] = par_tfin;
                        } else {
                            ok = false;
                        }

                        if ok {
                            x = beg_o.x + (seg_o.x * par_o[nbpi + 1]);
                            y = beg_o.y + (seg_o.y * par_o[nbpi + 1]);
                            if the_pi[0].pnt().distance(DVec3::new(x, y, 0.0)) >= (tolerance / 4.0)
                            {
                                nbpi += 1;
                                the_pi.push(IntfSectionPoint::new(
                                    DVec3::new(x, y, 0.0),
                                    IntfPIType::External,
                                    0,
                                    i_obje1,
                                    par_o[nbpi],
                                    IntfPIType::External,
                                    0,
                                    i_obje2,
                                    par_t[nbpi],
                                    sin_teta,
                                ));
                            }
                        }
                    } else {
                        // plus d une singularite
                        let mut par_omin = par_o[1];
                        let mut par_omax = par_o[1];
                        let mut par_tmin = par_t[1];
                        let mut par_tmax = par_t[1];
                        for item in 2..=nbpi {
                            par_omin = par_omin.min(par_o[item]);
                            par_omax = par_omax.max(par_o[item]);
                            par_tmin = par_tmin.min(par_t[item]);
                            par_tmax = par_tmax.max(par_t[item]);
                        }

                        let mut delta;
                        if par_odeb < 0.0 {
                            delta = -par_odeb;
                            par_odeb = 0.0;
                            par_tdeb = par_tdeb + sig_ps * (delta * (delta_t / delta_o));
                        }
                        if par_ofin > 1.0 {
                            delta = par_ofin - 1.0;
                            par_ofin = 1.0;
                            par_tfin = par_tfin - sig_ps * (delta * (delta_t / delta_o));
                        }
                        if sig_ps > 0.0 {
                            if par_tdeb < 0.0 {
                                delta = -par_tdeb;
                                par_tdeb = 0.0;
                                par_odeb = par_odeb + delta * (delta_o / delta_t);
                            }
                            if par_tfin > 1.0 {
                                delta = par_tfin - 1.0;
                                par_tfin = 1.0;
                                par_ofin = par_ofin - delta * (delta_o / delta_t);
                            }
                        } else {
                            if par_tdeb > 1.0 {
                                delta = par_tdeb - 1.0;
                                par_tdeb = 1.0;
                                par_odeb = par_odeb + delta * (delta_o / delta_t);
                            }
                            if par_tfin < 0.0 {
                                delta = -par_tfin;
                                par_tfin = 0.0;
                                par_ofin = par_ofin - delta * (delta_o / delta_t);
                            }
                        }

                        if (par_odeb < par_omin && par_omin > 0.0)
                            || (sig_ps > 0.0 && par_tdeb < par_tmin && par_tmin > 0.0)
                            || (sig_ps < 0.0 && par_tdeb > par_tmax && par_tmax < 1.0)
                        {
                            nbpi += 1;
                            par_o[nbpi] = 0.0f64.max(1.0f64.min(par_odeb));
                            par_t[nbpi] = 0.0f64.max(1.0f64.min(par_tdeb));
                            x = beg_o.x + (seg_o.x * par_o[nbpi]);
                            y = beg_o.y + (seg_o.y * par_o[nbpi]);
                            the_pi.push(IntfSectionPoint::new(
                                DVec3::new(x, y, 0.0),
                                IntfPIType::External,
                                0,
                                i_obje1,
                                par_o[nbpi],
                                IntfPIType::External,
                                0,
                                i_obje2,
                                par_t[nbpi],
                                sin_teta,
                            ));
                        }

                        if (par_ofin > par_omax && par_omax < 1.0)
                            || (sig_ps < 0.0 && par_tfin < par_tmin && par_tmin > 0.0)
                            || (sig_ps > 0.0 && par_tfin > par_tmax && par_tmax < 1.0)
                        {
                            nbpi += 1;
                            par_o[nbpi] = 1.0f64.min(0.0f64.max(par_ofin));
                            par_t[nbpi] = 1.0f64.min(0.0f64.max(par_tfin));
                            x = beg_o.x + (seg_o.x * par_o[nbpi]);
                            y = beg_o.y + (seg_o.y * par_o[nbpi]);
                            the_pi.push(IntfSectionPoint::new(
                                DVec3::new(x, y, 0.0),
                                IntfPIType::External,
                                0,
                                i_obje1,
                                par_o[nbpi],
                                IntfPIType::External,
                                0,
                                i_obje2,
                                par_t[nbpi],
                                sin_teta,
                            ));
                        }
                    }
                }
            }
        }

        //-- lbr : The points too close to each other are suspended.
        loop {
            let mut suppr = false;
            let mut i = 2usize;
            while !suppr && i <= nbpi {
                let p_im1 = the_pi[i - 2].pnt();
                let p_i = the_pi[i - 1].pnt();
                let mut d = p_i.distance(p_im1);
                d *= 50.0;
                if d < lg_t && d < lg_o {
                    for j in i..nbpi {
                        the_pi[j - 1] = the_pi[j].clone();
                    }
                    nbpi -= 1;
                    suppr = true;
                }
                i += 1;
            }
            if !suppr {
                break;
            }
        }

        if nbpi == 1 {
            if edge_sp {
                the_pi[0] = IntfSectionPoint::new(
                    DVec3::new(
                        beg_o.x + (seg_o.x * par_osp),
                        beg_o.y + (seg_o.y * par_osp),
                        0.0,
                    ),
                    IntfPIType::Edge,
                    0,
                    i_obje1,
                    par_osp,
                    IntfPIType::Edge,
                    0,
                    i_obje2,
                    par_tsp,
                    sin_teta,
                );
                par_o[1] = par_osp;
                par_t[1] = par_tsp;
            }
            if !self.interf.self_intf() {
                let mut contains = false;
                for i in 1..=self.interf.nb_section_points() {
                    if the_pi[0].is_equal(self.interf.pnt_value(i)) {
                        contains = true;
                        break;
                    }
                }
                if !contains {
                    let sp = the_pi[0].clone();
                    self.interf.my_s_poins_mut().push(sp);
                }
            } else if i_obje2 - i_obje1 != 1
                && (!self.o_clos || (i_obje1 != 1 && i_obje2 != self.nbso))
            {
                let sp = the_pi[0].clone();
                self.interf.my_s_poins_mut().push(sp);
            }
        } else if nbpi >= 2 {
            let mut the_tz = TangentZone::new();
            if nbpi == 2 {
                the_tz.polygon_insert(&the_pi[0]);
                the_tz.polygon_insert(&the_pi[1]);
            } else {
                let mut lmin = 1usize;
                let mut lmax = 1usize;
                for lpj in 2..=nbpi {
                    if par_o[lpj] < par_o[lmin] {
                        lmin = lpj;
                    } else if par_o[lpj] > par_o[lmax] {
                        lmax = lpj;
                    }
                }
                the_tz.polygon_insert(&the_pi[lmin - 1]);
                the_tz.polygon_insert(&the_pi[lmax - 1]);

                let mut ltmin = 1usize;
                let mut ltmax = 1usize;
                for lpj in 2..=nbpi {
                    if par_t[lpj] < par_t[ltmin] {
                        ltmin = lpj;
                    } else if par_t[lpj] > par_t[ltmax] {
                        ltmax = lpj;
                    }
                }
                if ltmin != lmin && ltmin != lmax {
                    the_tz.polygon_insert(&the_pi[ltmin - 1]);
                }
                if ltmax != lmin && ltmax != lmax {
                    the_tz.polygon_insert(&the_pi[ltmax - 1]);
                }
            }

            if edge_sp {
                the_tz.polygon_insert(&IntfSectionPoint::new(
                    DVec3::new(
                        beg_o.x + (seg_o.x * par_osp),
                        beg_o.y + (seg_o.y * par_osp),
                        0.0,
                    ),
                    IntfPIType::Edge,
                    0,
                    i_obje1,
                    par_osp,
                    IntfPIType::Edge,
                    0,
                    i_obje2,
                    par_tsp,
                    sin_teta,
                ));
            }

            let nb_tz = self.interf.nb_tangent_zones();
            let mut l_index: Vec<usize> = Vec::new();
            for ltz in 1..=nb_tz {
                if the_tz.has_common_range(&self.interf.my_t_zones_mut()[ltz - 1]) {
                    l_index.push(ltz);
                }
            }
            //----------------------------------------------------------------------
            //--   The list is parsed in ascending order by index, zone and tg
            //--
            if l_index.is_empty() {
                self.interf.my_t_zones_mut().push(the_tz);
            } else {
                let indexfirst = l_index.remove(0);
                let mut decal = 0usize;
                {
                    let zone = &mut self.interf.my_t_zones_mut()[indexfirst - 1];
                    zone.append_zone(&the_tz);
                }
                while !l_index.is_empty() {
                    let index = l_index.remove(0);
                    let other = self.interf.my_t_zones_mut()[index - 1 - decal].clone();
                    self.interf.my_t_zones_mut()[indexfirst - 1].append_zone(&other);
                    self.interf.my_t_zones_mut().remove(index - 1 - decal);
                    decal += 1;
                }
            }
        }
    }
}

impl Default for InterferencePolygon2d {
    /// OCCT Intf_InterferencePolygon2d() (cxx L36-42).
    fn default() -> Self {
        InterferencePolygon2d::new()
    }
}

/// OCCT gp_XY::CrossMagnitude — |X1.Y2 - X2.Y1|.
fn cross_magnitude(a: DVec2, b: DVec2) -> f64 {
    (a.x * b.y - a.y * b.x).abs()
}

/// OCCT `(gp_XY ^ gp_XY)` — the scalar 2D cross product.
fn cross_scalar(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// OCCT 2D SectionPoint ctor projection: myPnt(Where.X(), Where.Y(), 0.).
fn pnt2d_to_pnt(p: DVec2) -> DVec3 {
    DVec3::new(p.x, p.y, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::math::bnd::BndBox2d;

    /// A test polygon over explicit points (OCCT Intf_Polygon2d subclass
    /// shape).  1-based Segment(i) wraps on a closed polygon.
    struct TestPoly {
        pts: Vec<DVec2>,
        closed: bool,
        bnd: BndBox2d,
    }

    impl TestPoly {
        fn new(pts: &[DVec2], closed: bool) -> Self {
            let mut bnd = BndBox2d::new();
            for p in pts {
                bnd.add_point(*p);
            }
            TestPoly {
                pts: pts.to_vec(),
                closed,
                bnd,
            }
        }
    }

    impl IntfPolygon2d for TestPoly {
        fn bounding(&self) -> &BndBox2d {
            &self.bnd
        }
        fn bounding_mut(&mut self) -> &mut BndBox2d {
            &mut self.bnd
        }
        fn closed(&self) -> bool {
            self.closed
        }
        fn deflection_over_estimation(&self) -> f64 {
            0.0
        }
        fn nb_segments(&self) -> i32 {
            if self.closed {
                self.pts.len() as i32
            } else {
                self.pts.len() as i32 - 1
            }
        }
        fn segment(&self, the_index: i32) -> (DVec2, DVec2) {
            let n = self.pts.len() as i32;
            let i = the_index - 1;
            let j = (the_index) % n;
            (self.pts[i as usize], self.pts[j as usize])
        }
    }

    /// Open figure-eight polygon: segments 1 and 3 cross at (1, 1).
    /// (An open polygon keeps every non-adjacent crossing.)
    #[test]
    fn auto_interference_open_figure_eight() {
        let poly = TestPoly::new(
            &[
                DVec2::new(0.0, 0.0),
                DVec2::new(2.0, 2.0),
                DVec2::new(2.0, 0.0),
                DVec2::new(0.0, 2.0),
            ],
            false,
        );
        let itf = InterferencePolygon2d::new_polygon_auto(&poly);
        assert!(itf.interf.self_intf());
        assert_eq!(itf.interf.nb_section_points(), 1, "the (1,3) crossing expected");
        let sp = itf.interf.pnt_value(1);
        let p = sp.pnt();
        assert!((p.x - 1.0).abs() < 1.0e-9 && (p.y - 1.0).abs() < 1.0e-9, "got {p:?}");
        let (d1, a1, par1) = sp.info_first();
        let (d2, a2, par2) = sp.info_second();
        assert_eq!(d1, IntfPIType::Edge);
        assert_eq!(d2, IntfPIType::Edge);
        assert_eq!(a1, 1);
        assert_eq!(a2, 3);
        assert!((par1 - 0.5).abs() < 1.0e-9);
        assert!((par2 - 0.5).abs() < 1.0e-9);
        assert_eq!(itf.interf.nb_tangent_zones(), 0);
    }

    /// OCCT cxx L718 guard: on a CLOSED auto-interference the crossings
    /// involving the seam-adjacent segments (iObje1 == 1 or iObje2 == nbso)
    /// are suppressed, so the closed figure-eight (crossing on segment 1)
    /// yields no point.
    #[test]
    fn auto_interference_closed_seam_suppression() {
        let poly = TestPoly::new(
            &[
                DVec2::new(0.0, 0.0),
                DVec2::new(2.0, 2.0),
                DVec2::new(2.0, 0.0),
                DVec2::new(0.0, 2.0),
            ],
            true,
        );
        let itf = InterferencePolygon2d::new_polygon_auto(&poly);
        assert_eq!(itf.interf.nb_section_points(), 0, "closed seam guard drops it");
    }

    /// Two overlapping squares: the edges cross at (4, 2) and (2, 4).
    #[test]
    fn two_polygon_interference_crossing() {
        let s1 = TestPoly::new(
            &[
                DVec2::new(0.0, 0.0),
                DVec2::new(4.0, 0.0),
                DVec2::new(4.0, 4.0),
                DVec2::new(0.0, 4.0),
            ],
            true,
        );
        let s2 = TestPoly::new(
            &[
                DVec2::new(2.0, 2.0),
                DVec2::new(6.0, 2.0),
                DVec2::new(6.0, 6.0),
                DVec2::new(2.0, 6.0),
            ],
            true,
        );
        let itf = InterferencePolygon2d::new_polygon_polygon(&s1, &s2);
        assert!(!itf.interf.self_intf());
        assert_eq!(itf.interf.nb_section_points(), 2, "got {}", itf.interf.nb_section_points());
        let mut n_at_42 = 0;
        let mut n_at_24 = 0;
        for i in 1..=itf.interf.nb_section_points() {
            let p = itf.interf.pnt_value(i).pnt();
            if (p.x - 4.0).abs() < 1.0e-9 && (p.y - 2.0).abs() < 1.0e-9 {
                n_at_42 += 1;
            } else if (p.x - 2.0).abs() < 1.0e-9 && (p.y - 4.0).abs() < 1.0e-9 {
                n_at_24 += 1;
            }
        }
        assert_eq!((n_at_42, n_at_24), (1, 1));
    }

    /// Disjoint polygons produce no result.
    #[test]
    fn two_polygon_interference_disjoint() {
        let s1 = TestPoly::new(
            &[
                DVec2::new(0.0, 0.0),
                DVec2::new(1.0, 0.0),
                DVec2::new(1.0, 1.0),
                DVec2::new(0.0, 1.0),
            ],
            true,
        );
        let s2 = TestPoly::new(
            &[
                DVec2::new(5.0, 5.0),
                DVec2::new(6.0, 5.0),
                DVec2::new(6.0, 6.0),
                DVec2::new(5.0, 6.0),
            ],
            true,
        );
        let itf = InterferencePolygon2d::new_polygon_polygon(&s1, &s2);
        assert_eq!(itf.interf.nb_section_points(), 0);
        assert_eq!(itf.interf.nb_tangent_zones(), 0);
        // The bounding-box rejection leaves the sequences empty.
        assert!(itf.interf.nb_section_points() == 0);
    }

    /// Perform reuse: a second perform resets the sequences
    /// (SelfInterference).
    #[test]
    fn perform_resets_state() {
        let s1 = TestPoly::new(
            &[
                DVec2::new(0.0, 0.0),
                DVec2::new(4.0, 0.0),
                DVec2::new(4.0, 4.0),
                DVec2::new(0.0, 4.0),
            ],
            true,
        );
        let s2 = TestPoly::new(
            &[
                DVec2::new(2.0, 2.0),
                DVec2::new(6.0, 2.0),
                DVec2::new(6.0, 6.0),
                DVec2::new(2.0, 6.0),
            ],
            true,
        );
        let mut itf = InterferencePolygon2d::new_polygon_polygon(&s1, &s2);
        assert_eq!(itf.interf.nb_section_points(), 2);
        itf.perform_polygon_auto(&s1);
        // A simple closed square has no self-crossing; the seam-adjacent
        // suppression (cxx L718) drops the closing vertex too.
        assert_eq!(itf.interf.nb_section_points(), 0);
        let open8 = TestPoly::new(
            &[
                DVec2::new(0.0, 0.0),
                DVec2::new(2.0, 2.0),
                DVec2::new(2.0, 0.0),
                DVec2::new(0.0, 2.0),
            ],
            false,
        );
        itf.perform_polygon_auto(&open8);
        assert_eq!(itf.interf.nb_section_points(), 1, "reset + fresh auto result");
    }
}
