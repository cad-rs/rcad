// OCCT IntSurf_LineOn2S — a sequence of PntOn2S points forming an
// intersection line, with 1D/2D bounding boxes for fast rejection.
//
// OCCT IntSurf_LineOn2S.hxx / .cxx. rcad data-model notes:
// - mySeq: NCollection_Sequence<IntSurf_PntOn2S> -> Vec<PntOn2S>.
// - myBuv1/myBuv2 (Bnd_Box2d) and myBxyz (Bnd_Box) are rebuilt lazily from
//   the points (the OCCT box Add/IsWhole/Enlarge/IsOut semantics are
//   translated inline).

use glam::{DVec2, DVec3};

use crate::geomalgo::int_surf::PntOn2S;

/// A 2D box with void/whole states, translated from the Bnd_Box2d
/// operations used by IntSurf_LineOn2S.
#[derive(Debug, Clone)]
struct BndBox2d {
    // true = whole (infinite), false = normal box; void is a flag on top.
    whole: bool,
    void: bool,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    gap: f64,
}

impl BndBox2d {
    fn new() -> Self {
        BndBox2d { whole: false, void: true, x_min: 0.0, x_max: 0.0, y_min: 0.0, y_max: 0.0, gap: 0.0 }
    }
    fn set_whole(&mut self) {
        self.whole = true;
        self.void = false;
    }
    fn set_void(&mut self) {
        self.whole = false;
        self.void = true;
    }
    fn is_whole(&self) -> bool {
        self.whole
    }
    fn add(&mut self, p: DVec2) {
        if self.void {
            self.x_min = p.x;
            self.x_max = p.x;
            self.y_min = p.y;
            self.y_max = p.y;
            self.void = false;
            self.whole = false;
        } else {
            if p.x < self.x_min { self.x_min = p.x; }
            if p.x > self.x_max { self.x_max = p.x; }
            if p.y < self.y_min { self.y_min = p.y; }
            if p.y > self.y_max { self.y_max = p.y; }
        }
    }
    fn get(&self) -> (f64, f64, f64, f64) {
        (self.x_min, self.y_min, self.x_max, self.y_max)
    }
    fn enlarge(&mut self, t: f64) {
        self.x_min -= t;
        self.x_max += t;
        self.y_min -= t;
        self.y_max += t;
    }
    fn is_out(&self, p: DVec2) -> bool {
        if self.void {
            return true;
        }
        let g = self.gap;
        p.x < self.x_min - g || p.x > self.x_max + g || p.y < self.y_min - g || p.y > self.y_max + g
    }
}

/// A 3D box with void/whole states — translated from Bnd_Box.
#[derive(Debug, Clone)]
struct BndBox3 {
    whole: bool,
    void: bool,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    z_min: f64,
    z_max: f64,
    gap: f64,
}

impl BndBox3 {
    fn new() -> Self {
        BndBox3 { whole: false, void: true, x_min: 0.0, x_max: 0.0, y_min: 0.0, y_max: 0.0, z_min: 0.0, z_max: 0.0, gap: 0.0 }
    }
    fn set_whole(&mut self) {
        self.whole = true;
        self.void = false;
    }
    fn set_void(&mut self) {
        self.whole = false;
        self.void = true;
    }
    fn is_whole(&self) -> bool {
        self.whole
    }
    fn add(&mut self, p: DVec3) {
        if self.void {
            self.x_min = p.x; self.x_max = p.x;
            self.y_min = p.y; self.y_max = p.y;
            self.z_min = p.z; self.z_max = p.z;
            self.void = false;
            self.whole = false;
        } else {
            if p.x < self.x_min { self.x_min = p.x; }
            if p.x > self.x_max { self.x_max = p.x; }
            if p.y < self.y_min { self.y_min = p.y; }
            if p.y > self.y_max { self.y_max = p.y; }
            if p.z < self.z_min { self.z_min = p.z; }
            if p.z > self.z_max { self.z_max = p.z; }
        }
    }
    fn get(&self) -> (f64, f64, f64, f64, f64, f64) {
        (self.x_min, self.y_min, self.z_min, self.x_max, self.y_max, self.z_max)
    }
    fn enlarge(&mut self, t: f64) {
        self.x_min -= t;
        self.x_max += t;
        self.y_min -= t;
        self.y_max += t;
        self.z_min -= t;
        self.z_max += t;
    }
    fn is_out(&self, p: DVec3) -> bool {
        if self.void {
            return true;
        }
        let g = self.gap;
        p.x < self.x_min - g || p.x > self.x_max + g
            || p.y < self.y_min - g || p.y > self.y_max + g
            || p.z < self.z_min - g || p.z > self.z_max + g
    }
}

/// OCCT IntSurf_LineOn2S — a sequence of points on two surfaces.
#[derive(Debug, Clone)]
pub struct LineOn2S {
    seq: Vec<PntOn2S>,
    buv1: BndBox2d,
    buv2: BndBox2d,
    bxyz: BndBox3,
}

impl LineOn2S {
    /// OCCT IntSurf_LineOn2S() — myBuv1/2/3 SetWhole().
    pub fn new() -> Self {
        let mut buv1 = BndBox2d::new();
        let mut buv2 = BndBox2d::new();
        let mut bxyz = BndBox3::new();
        buv1.set_whole();
        buv2.set_whole();
        bxyz.set_whole();
        LineOn2S { seq: Vec::new(), buv1, buv2, bxyz }
    }

    /// OCCT Add(P).
    pub fn add(&mut self, p: &PntOn2S) {
        self.seq.push(p.clone());
        if !self.bxyz.is_whole() {
            self.bxyz.add(p.value());
        }
        if !self.buv1.is_whole() {
            self.buv1.add(p.value_on_surface(true));
        }
        if !self.buv2.is_whole() {
            self.buv2.add(p.value_on_surface(false));
        }
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.seq.len()
    }

    /// OCCT Value(Index) — 1-based in OCCT; rcad uses 0-based index.
    pub fn value(&self, index: usize) -> &PntOn2S {
        &self.seq[index]
    }

    /// OCCT Reverse().
    pub fn reverse(&mut self) {
        self.seq.reverse();
    }

    /// OCCT Split(Index) — keep points 0..index-1, return index..end.
    pub fn split(&mut self, index: usize) -> LineOn2S {
        let tail: Vec<PntOn2S> = self.seq.split_off(index);
        let mut ns = LineOn2S::new();
        for p in &tail {
            ns.add(p);
        }
        ns
    }

    /// OCCT Value(Index, P) — replace point at index.
    pub fn set_value(&mut self, index: usize, p: &PntOn2S) {
        self.seq[index] = p.clone();
    }

    /// OCCT SetPoint(Index, thePnt).
    pub fn set_point(&mut self, index: usize, pt: DVec3) {
        self.seq[index].set_value_pt(pt);
    }

    /// OCCT SetUV(Index, OnFirst, U, V).
    pub fn set_uv(&mut self, index: usize, on_first: bool, u: f64, v: f64) {
        self.seq[index].set_value_uv(on_first, u, v);
        if on_first && !self.buv1.is_whole() {
            self.buv1.add(DVec2::new(u, v));
        } else if !on_first && !self.buv2.is_whole() {
            self.buv2.add(DVec2::new(u, v));
        }
    }

    /// OCCT Clear().
    pub fn clear(&mut self) {
        self.seq.clear();
        self.buv1.set_whole();
        self.buv2.set_whole();
        self.bxyz.set_whole();
    }

    /// OCCT InsertBefore(I, P) — 1-based; rcad uses 0-based.
    pub fn insert_before(&mut self, index: usize, p: &PntOn2S) {
        if index > self.seq.len() {
            self.seq.push(p.clone());
        } else {
            self.seq.insert(index, p.clone());
        }
        if !self.bxyz.is_whole() {
            self.bxyz.add(p.value());
        }
        if !self.buv1.is_whole() {
            self.buv1.add(p.value_on_surface(true));
        }
        if !self.buv2.is_whole() {
            self.buv2.add(p.value_on_surface(false));
        }
    }

    /// OCCT RemovePoint(I) — 1-based; rcad uses 0-based.
    pub fn remove_point(&mut self, index: usize) {
        self.seq.remove(index);
        self.buv1.set_whole();
        self.buv2.set_whole();
        self.bxyz.set_whole();
    }

    /// OCCT IsOutBox(Pxyz).
    pub fn is_out_box(&self, p: DVec3) -> bool {
        if self.seq.is_empty() {
            return false;
        }
        if self.bxyz.is_whole() {
            let n = self.seq.len();
            let mut bx = BndBox3::new();
            bx.set_void();
            for i in 0..n {
                bx.add(self.seq[i].value());
            }
            let (x0, y0, z0, x1, y1, z1) = bx.get();
            let (x1, y1, z1) = (x1 - x0, y1 - y0, z1 - z0);
            if x1 > y1 {
                if x1 > z1 {
                    bx.enlarge(x1 * 0.01);
                } else {
                    bx.enlarge(z1 * 0.01);
                }
            } else if y1 > z1 {
                bx.enlarge(y1 * 0.01);
            } else {
                bx.enlarge(z1 * 0.01);
            }
            return bx.is_out(p);
        }
        self.bxyz.is_out(p)
    }

    /// OCCT IsOutSurf1Box(P1uv).
    pub fn is_out_surf1_box(&self, p1uv: DVec2) -> bool {
        if self.seq.is_empty() {
            return false;
        }
        if self.buv1.is_whole() {
            let n = self.seq.len();
            let mut b = BndBox2d::new();
            b.set_void();
            for i in 0..n {
                let (pu1, pv1, _, _) = self.seq[i].parameters();
                b.add(DVec2::new(pu1, pv1));
            }
            let (pu1, pv1, pu2, pv2) = b.get();
            let (pu2, pv2) = (pu2 - pu1, pv2 - pv1);
            if pu2 > pv2 {
                b.enlarge(pu2 * 0.01);
            } else {
                b.enlarge(pv2 * 0.01);
            }
            return b.is_out(p1uv);
        }
        self.buv1.is_out(p1uv)
    }

    /// OCCT IsOutSurf2Box(P2uv).
    pub fn is_out_surf2_box(&self, p2uv: DVec2) -> bool {
        if self.seq.is_empty() {
            return false;
        }
        if self.buv2.is_whole() {
            let n = self.seq.len();
            let mut b = BndBox2d::new();
            b.set_void();
            for i in 0..n {
                let (_, _, pu2, pv2) = self.seq[i].parameters();
                b.add(DVec2::new(pu2, pv2));
            }
            let (pu1, pv1, pu2, pv2) = b.get();
            let (pu2, pv2) = (pu2 - pu1, pv2 - pv1);
            if pu2 > pv2 {
                b.enlarge(pu2 * 0.01);
            } else {
                b.enlarge(pv2 * 0.01);
            }
            return b.is_out(p2uv);
        }
        self.buv2.is_out(p2uv)
    }
}

impl Default for LineOn2S {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pnt(x: f64, y: f64, z: f64, u1: f64, v1: f64, u2: f64, v2: f64) -> PntOn2S {
        let mut p = PntOn2S::new();
        p.set_value(DVec3::new(x, y, z), true, u1, v1);
        p.set_value_uv(false, u2, v2);
        p
    }

    #[test]
    fn add_and_value() {
        let mut l = LineOn2S::new();
        l.add(&pnt(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        l.add(&pnt(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0));
        assert_eq!(l.nb_points(), 2);
        assert_eq!(l.value(1).value(), DVec3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn split_reverse() {
        let mut l = LineOn2S::new();
        l.add(&pnt(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        l.add(&pnt(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0));
        l.add(&pnt(2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0));
        let tail = l.split(1);
        assert_eq!(l.nb_points(), 1);
        assert_eq!(tail.nb_points(), 2);
        l.reverse();
        assert_eq!(l.value(0).value(), DVec3::ZERO);
    }

    #[test]
    fn is_out_box() {
        let mut l = LineOn2S::new();
        l.add(&pnt(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        l.add(&pnt(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0));
        // points near the segment are inside the enlarged box
        assert!(!l.is_out_box(DVec3::new(0.5, 0.0, 0.0)));
        assert!(l.is_out_box(DVec3::new(10.0, 10.0, 10.0)));
    }
}
