// OCCT IntSurf_PathPoint (IntSurf_PathPoint.hxx/.cxx/.lxx) and
// IntSurf_InteriorPoint (IntSurf_InteriorPoint.hxx/.lxx) 1:1 Rust
// translations.  These are the point types carried by the
// IntStart_SearchInside / IntWalk_IWalking walking infrastructure.

use glam::{DVec2, DVec3};

/// OCCT IntSurf_PathPoint — a point of the parametric surface that lies on
/// the intersection curve (a starting point for an open line, or a passing
/// point along it).
#[derive(Clone, Debug)]
pub struct PathPoint {
    // OCCT: pt.
    pt: DVec3,
    // OCCT: ispass.
    ispass: bool,
    // OCCT: istgt.
    istgt: bool,
    // OCCT: vectg.
    vectg: DVec3,
    // OCCT: dirtg.
    dirtg: DVec2,
    // OCCT: sequv (NCollection_HSequence<gp_XY>) — first entry is the current
    // UV, the remaining entries are the multiplicity parameters.
    sequv: Vec<DVec2>,
}

impl PathPoint {
    /// OCCT IntSurf_PathPoint() (IntSurf_PathPoint.cxx L20-24).
    pub fn new() -> Self {
        PathPoint {
            pt: DVec3::ZERO,
            ispass: true,
            istgt: true,
            vectg: DVec3::ZERO,
            dirtg: DVec2::ZERO,
            sequv: Vec::new(),
        }
    }

    /// OCCT IntSurf_PathPoint(P, U, V) (L26-33).
    pub fn new_uv(pt: DVec3, u: f64, v: f64) -> Self {
        PathPoint {
            pt,
            ispass: true,
            istgt: true,
            vectg: DVec3::ZERO,
            dirtg: DVec2::ZERO,
            sequv: vec![DVec2::new(u, v)],
        }
    }

    /// OCCT SetValue(P, U, V) (L35-40).
    pub fn set_value(&mut self, pt: DVec3, u: f64, v: f64) {
        self.pt = pt;
        self.sequv.clear();
        self.sequv.push(DVec2::new(u, v));
    }

    /// OCCT AddUV(U, V) (lxx L20-23).
    pub fn add_uv(&mut self, u: f64, v: f64) {
        self.sequv.push(DVec2::new(u, v));
    }

    /// OCCT SetDirections(V, D) (lxx L25-31).
    pub fn set_directions(&mut self, v: DVec3, d: DVec2) {
        self.istgt = false;
        self.vectg = v;
        self.dirtg = d;
    }

    /// OCCT SetTangency(T) (lxx L33-36).
    pub fn set_tangency(&mut self, tang: bool) {
        self.istgt = tang;
    }

    /// OCCT SetPassing(P) (lxx L38-41).
    pub fn set_passing(&mut self, pass_: bool) {
        self.ispass = pass_;
    }

    /// OCCT Value() (lxx L43-45).
    pub fn value(&self) -> DVec3 {
        self.pt
    }

    /// OCCT Value2d(U, V) (lxx L47-52).
    pub fn value_2d(&self) -> DVec2 {
        self.sequv.first().copied().unwrap_or(DVec2::ZERO)
    }

    /// OCCT IsPassingPnt() (lxx L54-56).
    pub fn is_passing_pnt(&self) -> bool {
        self.ispass
    }

    /// OCCT IsTangent() (lxx L58-60).
    pub fn is_tangent(&self) -> bool {
        self.istgt
    }

    /// OCCT Direction3d() (lxx L62-68) — throws when tangent.
    pub fn direction_3d(&self) -> DVec3 {
        if self.istgt {
            return DVec3::ZERO;
        }
        self.vectg
    }

    /// OCCT Direction2d() (lxx L70-76).
    pub fn direction_2d(&self) -> DVec2 {
        if self.istgt {
            return DVec2::ZERO;
        }
        self.dirtg
    }

    /// OCCT Multiplicity() (lxx L78-80) — number of extra parameter pairs.
    pub fn multiplicity(&self) -> i32 {
        (self.sequv.len() - 1) as i32
    }

    /// OCCT Parameters(Index, U, V) (lxx L82-85) — 1-based index.
    pub fn parameters(&self, index: i32, u: &mut f64, v: &mut f64) {
        let uv = self.sequv[index as usize + 1];
        *u = uv.x;
        *v = uv.y;
    }
}

impl Default for PathPoint {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT IntSurf_InteriorPoint — a point solution of the intersection between
/// an implicit and a parametrised surface, inside the domain.  These are the
/// starting points for closed intersection lines.
#[derive(Clone, Debug)]
pub struct InteriorPoint {
    point: DVec3,
    paramu: f64,
    paramv: f64,
    direc: DVec3,
    direc2d: DVec2,
}

impl InteriorPoint {
    /// OCCT IntSurf_InteriorPoint() — default.
    pub fn new() -> Self {
        InteriorPoint {
            point: DVec3::ZERO,
            paramu: 0.0,
            paramv: 0.0,
            direc: DVec3::ZERO,
            direc2d: DVec2::ZERO,
        }
    }

    /// OCCT IntSurf_InteriorPoint(P, U, V, Direc, Direc2d).
    pub fn new_full(point: DVec3, u: f64, v: f64, direc: DVec3, direc2d: DVec2) -> Self {
        InteriorPoint {
            point,
            paramu: u,
            paramv: v,
            direc,
            direc2d,
        }
    }

    /// OCCT SetValue(P, U, V, Direc, Direc2d).
    pub fn set_value(&mut self, point: DVec3, u: f64, v: f64, direc: DVec3, direc2d: DVec2) {
        self.point = point;
        self.paramu = u;
        self.paramv = v;
        self.direc = direc;
        self.direc2d = direc2d;
    }

    /// OCCT Value().
    pub fn value(&self) -> DVec3 {
        self.point
    }

    /// OCCT Parameters(U, V).
    pub fn parameters(&self, u: &mut f64, v: &mut f64) {
        *u = self.paramu;
        *v = self.paramv;
    }

    /// OCCT UParameter().
    pub fn u_parameter(&self) -> f64 {
        self.paramu
    }

    /// OCCT VParameter().
    pub fn v_parameter(&self) -> f64 {
        self.paramv
    }

    /// OCCT Direction() — 3D tangent of the intersection.
    pub fn direction(&self) -> DVec3 {
        self.direc
    }

    /// OCCT Direction2d() — 2D tangent in the parametric space.
    pub fn direction_2d(&self) -> DVec2 {
        self.direc2d
    }
}

impl Default for InteriorPoint {
    fn default() -> Self {
        Self::new()
    }
}
