// OCCT BRepClass3d_SolidClassifier (BRepClass3d_SolidClassifier.hxx / .cxx)
// Provides an algorithm to classify a point in a solid.
// Inherits from BRepClass3d_SClassifier.

use crate::topalgo::brep_class3d::solid_explorer::SolidExplorer;
use crate::topalgo::brep_class3d::s_classifier::SClassifier;
use rcad_kernel::geom::{Surface3, SurfaceEval};
use rcad_kernel::topods::Orientation;
use rcad_kernel::topo_shape::Shape;
use glam::DVec3;

/// OCCT BRepClass3d_SolidClassifier — classifies a point relative to a solid.
pub struct SolidClassifier {
    // BRepClass3d_SClassifier base
    pub my_state: u8,          // 0=unknown, 1=faulty, 2=ON, 3=IN, 4=OUT
    // BRepClass3d_SolidClassifier own
    pub a_solid_loaded: bool,
    pub explorer: SolidExplorer,
    pub is_a_hole_in_space: bool,
}

impl SolidClassifier {
    /// Empty constructor.
    pub fn new() -> Self {
        SolidClassifier {
            my_state: 0,
            a_solid_loaded: false,
            explorer: SolidExplorer::new(),
            is_a_hole_in_space: false,
        }
    }

    /// OCCT: Load(const TopoDS_Shape& S) — initialize explorer for the given solid.
    pub fn load(&mut self, s: &Shape) {
        if self.a_solid_loaded {
            self.explorer.clear();
        }
        self.explorer.init_shape(s);
        self.a_solid_loaded = true;
    }

    /// Constructor from a Shape.
    pub fn from_shape(s: &Shape) -> Self {
        let clsf = SolidClassifier {
            my_state: 0,
            a_solid_loaded: true,
            explorer: SolidExplorer::from_shape(s),
            is_a_hole_in_space: false,
        };
        clsf
    }

    /// Constructor to classify the point P with tolerance Tol on the solid S.
    pub fn from_shape_point(s: &Shape, p: DVec3, tol: f64) -> Self {
        let mut clsf = SolidClassifier {
            my_state: 0,
            a_solid_loaded: true,
            explorer: SolidExplorer::from_shape(s),
            is_a_hole_in_space: false,
        };
        clsf.perform(p, tol);
        clsf
    }

    /// OCCT: Perform(P, Tol) — classify point P relative to loaded solid.
    pub fn perform(&mut self, p: DVec3, tol: f64) {
        if !self.a_solid_loaded { return; }
        // OCCT L171-191: check bounding box first (fast rejection)
        if !self.is_a_hole_in_space {
            if self.explorer.reject(p) {
                self.my_state = 4; // OUT
                return;
            }
            // OCCT BRepClass3d_SClassifier::Perform — a point within the
            // tolerance of a face of the solid is ON the boundary.
            if self.explorer.point_on_face(p, tol) {
                self.my_state = 2; // ON
                return;
            }
            // OCCT L190: BRepClass3d_SClassifier::Perform(explorer, P, Tol)
            self.perform_classify(p, tol);
        } else {
            if self.explorer.reject(p) {
                self.my_state = 3; // IN
                return;
            }
            self.perform_classify(p, tol);
        }
    }

    /// OCCT: PerformInfinitePoint(Tol) — classify point at infinity.
    /// BRepClass3d_SolidClassifier::PerformInfinitePoint (SolidClassifier.cxx)
    /// delegates to BRepClass3d_SClassifier::PerformInfinitePoint and stores
    /// isaholeinspace = State() != OUT. The SClassifier algorithm (SClassifier.cxx
    /// L125-198): take a normal to a random inner point of each face and
    /// intersect this reversed normal with all the faces of the solid. If the
    /// min.par. intersection point is an inner point of a face and the
    /// transition is not TANGENT, set myState to IN or OUT according to the
    /// transition (Out -> IN, In -> OUT); otherwise try the next probing point.
    /// Used by IsHole (BOPAlgo_BuilderSolid.cxx L823-831): State() == IN.
    pub fn perform_infinite_point(&mut self, _tol: f64) {
        // OCCT L146-148: aSE.Reject(gp_Pnt(0, 0, 0)) — a solid without faces
        // is treated as a void solid (whole space) -> IN.
        if self.explorer.face_surfaces.is_empty() {
            self.my_state = 3; // IN
            return;
        }
        // OCCT: myState = 2 (ON) — the default when no definitive answer.
        self.my_state = 2;
        let faces = self.explorer.face_surfaces.clone();
        // OCCT: up to 10 random probing parameters per face; aParam in
        // [0.1, 0.9] (math_BullardGenerator). rcad uses a deterministic
        // pseudo-random sequence — the probing point must simply lie inside
        // the face, the exact position does not matter.
        for itry in 0..10 {
            // The first probe is the face center (param 0.5) — the OCCT random
            // inner point must not degenerate on edges/vertices; then spread
            // over [0.1, 0.9].
            let param = if itry == 0 {
                0.5
            } else {
                0.1 + 0.8 * (((itry as f64 + 1.0) * 0.6180339887498949).fract())
            };
            for f in &faces {
                // OCCT L158-159: FindAPointInTheFace(aF, aPoint, aU, aV, aParam).
                let Some((a_point, a_u, a_v)) = SolidExplorer::face_point(f, param) else {
                    continue;
                };
                // OCCT L160-162: FaceNormal(aF, aU, aV, aDN).
                let Some(a_dn) = SolidExplorer::face_outward_normal_at(f, a_u, a_v) else {
                    continue;
                };
                if a_dn.length_squared() < 1e-12 {
                    continue;
                }
                // OCCT L163: gp_Lin aLin(aPoint, -aDN).
                let ray_dir = -a_dn;
                // OCCT L165-181: intersect the ray with all faces; the minimal
                // WParameter wins (strict <, so among equal parameters the first
                // one wins). parmin limits the search range of the later faces
                // (Intersector3d.Perform(aLin, -RealLast(), parmin)).
                let mut parmin = f64::MAX;
                let mut best = 0u8; // 0 = no valid transition, 1 = In, 2 = Out
                for g in &faces {
                    // Outward normal of g for the transition sign. OCCT reads the
                    // intersection state/transition from IntCurvesFace; rcad
                    // derives the transition from the surface normal at the
                    // probing point of g (constant on planar faces).
                    let g_n = match SolidExplorer::face_point(g, param) {
                        Some((_, gu, gv)) => SolidExplorer::face_outward_normal_at(g, gu, gv),
                        None => None,
                    };
                    let Some(g_n) = g_n else { continue };
                    let w = match &g.surf {
                        Surface3::Plane(pl) => {
                            let normal = if g.ori == Orientation::Reversed {
                                -pl.normal
                            } else {
                                pl.normal
                            };
                            let denom = ray_dir.dot(normal);
                            if denom.abs() < 1e-12 {
                                continue;
                            }
                            let w = (pl.origin - a_point).dot(normal) / denom;
                            if w < 0.0 {
                                continue;
                            }
                            // OCCT IntCurvesFace_Intersector.cxx L291: accept
                            // only intersections whose UV lies inside the face
                            // domain (Classify(Puv) == IN/ON).
                            if let Some([umin, umax, vmin, vmax]) = g.uv_bounds {
                                let pint = a_point + ray_dir * w;
                                let rel = pint - pl.origin;
                                let u = rel.dot(pl.u_dir);
                                let v = rel.dot(pl.v_dir);
                                if u < umin || u > umax || v < vmin || v > vmax {
                                    continue;
                                }
                            }
                            Some(w)
                        }
                        _ => SolidExplorer::ray_face_param(a_point, ray_dir, g),
                    };
                    let Some(w) = w else { continue };
                    if w < parmin {
                        parmin = w;
                        // OCCT int_cs transition: cos_dir = nSurf . dirCurve;
                        // < 0 -> In, > 0 -> Out (IntCurveSurface_InterUtils.pxx
                        // L856-895).
                        best = if ray_dir.dot(g_n) > 0.0 { 2 } else { 1 };
                    }
                }
                // OCCT L184-195: a definitive transition decides the state of
                // the point at infinity: Out -> IN (hole in space), In -> OUT.
                if best == 2 {
                    self.my_state = 3; // IN
                    return;
                } else if best == 1 {
                    self.my_state = 4; // OUT
                    return;
                }
            }
        }
    }

    /// OCCT: State() — returns my_state (inherited from SClassifier).
    pub fn state(&self) -> u8 { self.my_state }

    /// Internal: perform classification using ray casting.
    /// OCCT BRepClass3d_SClassifier::Perform (L203+). Delegates the actual
    /// ray/face intersection to the explorer (BRepClass3d_SolidExplorer),
    /// which resolves face geometry from the shape tree directly.
    fn perform_classify(&mut self, p: DVec3, _tol: f64) {
        if !self.explorer.has_faces() {
            self.my_state = 3; // IN (void solid — whole space)
            return;
        }
        self.my_state = self.explorer.classify_point(p);
    }
}
