//! OCCT Intf_Interference + Intf_Polygon2d (TKGeomAlgo Intf package).
//!
//! 1:1 translation of `Intf_Interference.hxx` (L37-98) +
//! `Intf_Interference.cxx` (L26-294) + `Intf_Interference.lxx` (L22-77),
//! and `Intf_Polygon2d.hxx` (L29-55) + `Intf_Polygon2d.lxx` — the abstract
//! base class maps to a Rust trait over the shared `BndBox2d`.

use glam::DVec2;
use rcad_kernel::math::bnd::BndBox2d;

use super::intf::IntfSectionPoint;
use super::intf_section_line::SectionLine;
use super::intf_tangent_zone::TangentZone;

/// OCCT Intf_Polygon2d — the polygon information required to compute
/// interferences (abstract base).  The protected `myBox` maps to the
/// trait's implementor-owned box returned through [`Self::bounding_mut`].
pub trait IntfPolygon2d {
    /// OCCT Bounding() — lxx L17-24: returns the bounding box of the
    /// polygon.
    fn bounding(&self) -> &BndBox2d;

    /// Mutable access to the implementor's `myBox` (OCCT subclasses fill
    /// the protected member directly).
    fn bounding_mut(&mut self) -> &mut BndBox2d;

    /// OCCT Closed() — hxx L38: returns True if the polyline is closed.
    fn closed(&self) -> bool;

    /// OCCT DeflectionOverEstimation() — hxx L43.
    fn deflection_over_estimation(&self) -> f64;

    /// OCCT NbSegments() — hxx L46.
    fn nb_segments(&self) -> i32;

    /// OCCT Segment(theIndex, theBegin, theEnd) — hxx L49-51.
    fn segment(&self, the_index: i32) -> (DVec2, DVec2);
}

/// OCCT Intf_Interference — the interference computation result between
/// polygon2d / polygon3d / polyhedron: sequences of intersection points,
/// intersection polylines and tangent zones.
#[derive(Debug, Clone)]
pub struct Interference {
    /// NCollection_Sequence<Intf_SectionPoint> (1-based in OCCT).
    my_s_poins: Vec<IntfSectionPoint>,
    my_s_lines: Vec<SectionLine>,
    my_t_zones: Vec<TangentZone>,
    self_intf: bool,
    tolerance: f64,
}

impl Interference {
    /// OCCT Intf_Interference(Self) — cxx L26-30 (protected constructor:
    /// initialize for a deferred interference).
    pub fn with_self(self_: bool) -> Self {
        Interference {
            my_s_poins: Vec::new(),
            my_s_lines: Vec::new(),
            my_t_zones: Vec::new(),
            self_intf: self_,
            tolerance: 0.0,
        }
    }

    /// OCCT SelfInterference(Self) — cxx L37-43: reset before a new
    /// perform.
    pub fn self_interference(&mut self, self_: bool) {
        self.self_intf = self_;
        self.my_s_poins.clear();
        self.my_s_lines.clear();
        self.my_t_zones.clear();
    }

    /// OCCT SelfIntf flag read (protected member in OCCT; the subclasses
    /// test it inside Perform).
    pub fn self_intf(&self) -> bool {
        self.self_intf
    }

    /// Mutable access to Tolerance (protected member; the concrete
    /// Interference* subclasses assign it during Perform).
    pub fn set_tolerance(&mut self, t: f64) {
        self.tolerance = t;
    }

    /// OCCT NbSectionPoints() — lxx L22-25.
    pub fn nb_section_points(&self) -> usize {
        self.my_s_poins.len()
    }

    /// OCCT PntValue(Index) — lxx L31-34 (1-based).
    pub fn pnt_value(&self, index: usize) -> &IntfSectionPoint {
        &self.my_s_poins[index - 1]
    }

    /// OCCT NbSectionLines() — lxx L40-43.
    pub fn nb_section_lines(&self) -> usize {
        self.my_s_lines.len()
    }

    /// OCCT LineValue(Index) — lxx L49-52 (1-based).
    pub fn line_value(&self, index: usize) -> &SectionLine {
        &self.my_s_lines[index - 1]
    }

    /// OCCT NbTangentZones() — lxx L58-61.
    pub fn nb_tangent_zones(&self) -> usize {
        self.my_t_zones.len()
    }

    /// OCCT ZoneValue(Index) — lxx L67-70 (1-based).
    pub fn zone_value(&self, index: usize) -> &TangentZone {
        &self.my_t_zones[index - 1]
    }

    /// OCCT GetTolerance() — lxx L74-77.
    pub fn get_tolerance(&self) -> f64 {
        self.tolerance
    }

    /// OCCT Insert(LaZone) — cxx L50-173: inserts a tangent zone into the
    /// list, merging it with a connected existing zone when possible.
    pub fn insert_zone(&mut self, la_zone: &TangentZone) -> bool {
        if self.my_t_zones.len() == 0 {
            return false;
        }
        let mut lzin = 0usize; // Index in the list of the zone of interest.
        let mut lunp: i32 = 0; // Index of the 1st stop point in this zone.
        let mut lotp: i32 = 0; // Index of the 2nd stop point in this zone.
        let mut lunl: i32 = 0; // Index of the 1st point of the new zone.
        let mut lotl: i32 = 0; // Index of the 2nd point of the new zone.
        let mut same = false; // Search direction of the stop of the new zone.
        let mut inserted = true;
        let nplz = la_zone.number_of_points() as i32;

        // Loop on TangentZone:
        'zones: for iz in 1..=self.my_t_zones.len() {
            // Loop on edges of the TangentZone:
            let npcz = self.my_t_zones[iz - 1].number_of_points() as i32;
            for ipz1 in 1..=npcz {
                let ipz0 = if ipz1 - 1 <= 0 { npcz } else { ipz1 - 1 };
                let ipz2 = (ipz1 % npcz) + 1;

                // Loop on edges of the new TangentZone and search of the
                // corresponding point or edge:
                for ilz1 in 1..=nplz {
                    let ilz2 = (ilz1 % nplz) + 1;

                    if self.my_t_zones[iz - 1]
                        .get_point(ipz1 as usize)
                        .is_equal(&la_zone.get_point(ilz1 as usize))
                    {
                        if self.my_t_zones[iz - 1]
                            .get_point(ipz0 as usize)
                            .is_equal(&la_zone.get_point(ilz2 as usize))
                        {
                            lzin = iz;
                            lunp = ipz0;
                            lotp = ipz1;
                            lunl = ilz1;
                            lotl = ilz2;
                            same = false;
                            break;
                        } else if self.my_t_zones[iz - 1]
                            .get_point(ipz2 as usize)
                            .is_equal(&la_zone.get_point(ilz2 as usize))
                        {
                            lzin = iz;
                            lunp = ipz1;
                            lotp = ipz2;
                            lunl = ilz1;
                            lotl = ilz2;
                            same = true;
                            break;
                        } else {
                            lzin = iz;
                            lunp = ipz1;
                            lunl = ilz1;
                        }
                    }
                }
                if lotp != 0 {
                    break;
                }
            }
            if lotp != 0 {
                break 'zones;
            }
        }

        if lotp != 0 {
            let mut lotp = lotp;
            let mut ilc = lotl + 1;
            while ((ilc - 1) % nplz) + 1 != lunl {
                let p = la_zone.get_point((((ilc - 1) % nplz) + 1) as usize);
                self.my_t_zones[lzin - 1].insert_before(lotp as usize, &p);
                if !same {
                    lotp += 1;
                }
                ilc += 1;
            }
        } else if lunp > 0 {
            let mut lunp = lunp;
            let mut ilc = lunl;
            let mut looped = false;
            loop {
                let p = la_zone.get_point((((ilc - 1) % nplz) + 1) as usize);
                self.my_t_zones[lzin - 1].insert_before(lunp as usize, &p);
                lunp += 1;
                ilc += 1;
                if looped && ((ilc - 1 - 1) % nplz) + 1 == lunl {
                    break;
                }
                looped = true;
            }
        } else {
            inserted = false;
        }

        if inserted {
            let the_new = self.my_t_zones[lzin - 1].clone();
            self.my_t_zones.remove(lzin - 1);
            if !self.insert_zone(&the_new) {
                self.my_t_zones.push(the_new);
            }
        }
        inserted
    }

    /// OCCT Insert(pdeb, pfin) — cxx L177-281: inserts a new segment of
    /// intersection into the list of polylines, joining connected lines.
    pub fn insert_segment(&mut self, pdeb: &IntfSectionPoint, pfin: &IntfSectionPoint) {
        let mut inserted = false;
        let mut the_ls: usize = 0;
        let mut begin = false;
        let mut the_bout = *pfin;

        for ils in 1..=self.my_s_lines.len() {
            // OCCT iterates the sequence by mutable reference; mirror with a
            // temporary borrow scope.
            let (nd, nf) = {
                let sl = &mut self.my_s_lines[ils - 1];
                (sl.is_end(pdeb), sl.is_end(pfin))
            };
            if nd == 1 {
                if nf > 1 {
                    self.my_s_lines[ils - 1].close();
                }
                inserted = true;
                the_ls = ils;
                begin = true;
                break;
            } else if nd > 1 {
                if nf == 1 {
                    self.my_s_lines[ils - 1].close();
                }
                inserted = true;
                the_ls = ils;
                begin = false;
                break;
            } else if nf == 1 {
                inserted = true;
                the_ls = ils;
                begin = true;
                the_bout = *pdeb;
                break;
            } else if nf > 1 {
                inserted = true;
                the_ls = ils;
                begin = false;
                the_bout = *pdeb;
                break;
            }
        }

        if !inserted {
            let mut la_ls = SectionLine::new();
            la_ls.append(pdeb);
            la_ls.append(pfin);
            self.my_s_lines.push(la_ls);
        } else {
            let mut nd: i32 = 0;
            for ils in 1..=self.my_s_lines.len() {
                if ils != the_ls {
                    nd = self.my_s_lines[ils - 1].is_end(&the_bout);
                    if nd == 1 {
                        if begin {
                            self.my_s_lines[the_ls - 1].reverse();
                        }
                        let taken = self.my_s_lines[the_ls - 1].clone();
                        self.my_s_lines[ils - 1].prepend_line(&mut taken.clone());
                        break;
                    } else if nd > 1 {
                        if !begin {
                            self.my_s_lines[the_ls - 1].reverse();
                        }
                        let taken = self.my_s_lines[the_ls - 1].clone();
                        self.my_s_lines[ils - 1].append_line(&mut taken.clone());
                        break;
                    }
                }
            }
            if nd > 0 {
                self.my_s_lines.remove(the_ls - 1);
            } else if begin {
                self.my_s_lines[the_ls - 1].prepend(&the_bout);
            } else {
                self.my_s_lines[the_ls - 1].append(&the_bout);
            }
        }
    }

    /// OCCT Contains(LePnt) — cxx L287-306.
    pub fn contains(&self, le_pnt: &IntfSectionPoint) -> bool {
        for l in 1..=self.my_s_lines.len() {
            if self.my_s_lines[l - 1].contains(le_pnt) {
                return true;
            }
        }
        for t in 1..=self.my_t_zones.len() {
            if self.my_t_zones[t - 1].contains(le_pnt) {
                return true;
            }
        }
        false
    }
}
