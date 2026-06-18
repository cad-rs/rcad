use glam::DVec2;
use crate::bopds::ds::DS;

/// OCCT-aligned: IntTools_FClass2d -- 2D point-in-face classifier.
/// OCCT IntTools_FClass2d.hxx / IntTools_FClass2d.cxx
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State { In, On, Out }

pub struct FClass2d {
    uv_boundary: Vec<Vec<DVec2>>,
}

impl FClass2d {
    /// Build from DS face index. Uses ds.faces[fi].uv_boundary.
    pub fn from_ds_face(ds: &DS, fi: usize) -> Self {
        let uv_boundary = match ds.faces.get(fi) {
            Some(f) => f.uv_boundary.clone().map(|b| vec![b]).unwrap_or_default(),
            None => vec![],
        };
        FClass2d { uv_boundary }
    }

    /// Classify a UV point: IN/ON/OUT.
    /// Uses winding number for point-in-polygon on the outer boundary.
    /// If inside outer and not on/inner boundary -> In.
    /// If on any boundary edge -> On.
    /// Otherwise -> Out.
    pub fn perform(&self, uv: DVec2) -> State {
        for boundary in &self.uv_boundary {
            if boundary.len() < 3 { continue; }
            // Check ON boundary
            for i in 0..boundary.len() {
                let j = (i + 1) % boundary.len();
                let a = boundary[i];
                let b = boundary[j];
                let ab = b - a;
                let ap = uv - a;
                let t = ap.dot(ab) / (ab.length_squared() + 1e-30);
                if t >= -1e-8 && t <= 1.0 + 1e-8 {
                    let proj = a + ab * t.clamp(0.0, 1.0);
                    if (uv - proj).length_squared() < 1e-14 {
                        return State::On;
                    }
                }
            }
            // Winding number for point-in-polygon
            let mut winding = 0i32;
            let n = boundary.len();
            for i in 0..n {
                let j = (i + 1) % n;
                let a = boundary[i];
                let b = boundary[j];
                if a.y <= uv.y {
                    if b.y > uv.y && ((b.x - a.x) * (uv.y - a.y) - (b.y - a.y) * (uv.x - a.x)) > 0.0 {
                        winding += 1;
                    }
                } else if b.y <= uv.y && ((b.x - a.x) * (uv.y - a.y) - (b.y - a.y) * (uv.x - a.x)) < 0.0 {
                    winding -= 1;
                }
            }
            if winding == 0 { return State::Out; }
        }
        State::In
    }

    /// OCCT: IsHole() -- check if the face's outer boundary encloses no interior.
    /// Returns true if the outer wire is a hole (interior to the infinite face).
    /// Simplified: returns true if the winding number of the outer boundary
    /// is negative (CW in UV space = hole).
    pub fn is_hole(&self) -> bool {
        if self.uv_boundary.is_empty() { return true; }
        // Compute signed area: negative = CW = hole
        let outer = &self.uv_boundary[0];
        let mut area = 0.0;
        for i in 0..outer.len() {
            let j = (i + 1) % outer.len();
            area += outer[i].x * outer[j].y - outer[j].x * outer[i].y;
        }
        area < 0.0
    }
}
