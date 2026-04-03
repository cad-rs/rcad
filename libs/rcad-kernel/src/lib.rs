use glam::DVec3;
use serde::{Deserialize, Serialize};

/// Geometric (analytic) model types: position, curve, surface, primitive descriptors.
///
/// This module describes *what shape is*.
pub mod geom;

/// Topology model types: vertex/edge/face/shell/solid incidence relationships.
///
/// This module describes *how things are connected*.
pub mod topology;

/// Shape properties: surface area, volume, centroid.
///
/// Analogous to OCCT `GProp_GProps` + `BRepGProp`.
pub mod properties;

/// Topology query helpers: edge adjacency, vertex adjacency, shape counts.
///
/// Analogous to OCCT `TopExp_Explorer` and `TopExp::MapShapesAndAncestors`.
pub mod topo_query;

/// Differential geometry: principal curvatures, Gaussian curvature, mean curvature.
///
/// Analogous to OCCT `GeomLProp_SLProps`.
pub mod curvature;

/// Curve arc-length computation.
///
/// Analogous to OCCT `GCPnts_AbscissaPoint` / `CPnts_AbscissaPoint::Length`.
pub mod arc_length;

/// Visual appearance: per-face/solid RGB color and basic material.
///
/// Analogous to OCCT `XCAFDoc_ColorTool`.
pub mod appearance;

/// Precision constants and per-entity tolerance query helpers.
///
/// Analogous to OCCT `Precision` class and `BRep_Tool::Tolerance`.
pub mod tolerance;

pub use geom::PrimitiveSolid;
pub use geom::{Curve2d, Curve3, Surface3};
pub use geom::{any_perpendicular, Curve2dEval, CurveEval, SurfaceEval};
pub use geom::BSplineCurve2;
pub use topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
pub use properties::{centroid, inertia_tensor, surface_area, volume, InertiaTensor};
pub use topo_query::{edge_adjacent_faces, edge_count, face_count, face_edges,
                     vertex_adjacent_edges, vertex_count};
pub use curvature::{gaussian_curvature, mean_curvature, principal_curvatures};
pub use arc_length::arc_length;
pub use appearance::{Color, FaceColor, StepColor};
pub use tolerance::{
    ANGULAR, APPROXIMATION, CONFUSION,
    edge_tolerance, face_tolerance, model_tolerance, vertex_tolerance,
};

/// A parameter-space curve binding that ties a 3D edge to an adjacent face's
/// surface parameter domain (u, v).  Analogous to OCCT `BRep_CurveOnSurface`.
///
/// `surface_idx` indexes into `GeomStore.surfaces`.
/// `curve2d_idx` indexes into `GeomStore.curve2ds`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PCurve {
    pub surface_idx: usize,
    pub curve2d_idx: usize,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GeomStore {
    /// Pool of 3D analytic curves.
    pub curves: Vec<Curve3>,
    /// Pool of analytic surfaces.
    pub surfaces: Vec<Surface3>,
    /// Pool of 2D parameter-space curves used by PCurves.
    pub curve2ds: Vec<Curve2d>,
    /// Indexed by `BRep.edges` index; value is index into `curves`.
    pub edge_curve: Vec<Option<usize>>,
    /// Flattened face order across solids/shells; value is index into `surfaces`.
    pub face_surface: Vec<Option<usize>>,
    /// Indexed by `BRep.edges` index; each entry is the list of PCurves for
    /// that edge on its adjacent faces (usually 1, seam edges have 2).
    pub edge_pcurves: Vec<Vec<PCurve>>,
    /// Parallel to `edge_curve`: the parameter range [t1, t2] of the edge on its
    /// 3D curve. `None` = unknown (algorithms fall back to `CurveEval::default_domain`).
    /// Analogous to `BRep_Edge::Range()` in OCCT.
    #[serde(default)]
    pub edge_curve_range: Vec<Option<[f64; 2]>>,
    /// Parallel to `BRep.edges`: `true` if this is a degenerate edge (zero-length,
    /// e.g. a polar singularity). Analogous to `BRep_Edge::Degenerated()` in OCCT.
    #[serde(default)]
    pub edge_degenerated: Vec<bool>,
    /// Per-vertex tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to `BRep.vertices`. Analogous to `BRep_Tool::Tolerance(vertex)` in OCCT.
    #[serde(default)]
    pub vertex_tolerance: Vec<f64>,
    /// Per-edge tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to `BRep.edges`. Analogous to `BRep_Tool::Tolerance(edge)` in OCCT.
    #[serde(default)]
    pub edge_tolerance: Vec<f64>,
    /// Per-face tolerance (falls back to `tolerance::CONFUSION` when absent or zero).
    /// Parallel to the flattened face order (same indexing as `face_surface`).
    /// Analogous to `BRep_Tool::Tolerance(face)` in OCCT.
    #[serde(default)]
    pub face_tolerance: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub solids: Vec<Solid>,
    #[serde(default)]
    pub geom: GeomStore,
}

impl Default for BRep {
    fn default() -> Self {
        Self::new()
    }
}

impl BRep {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            solids: Vec::new(),
            geom: GeomStore::default(),
        }
    }

    /// Creates a unit box B-Rep.
    ///
    /// Vertex layout:
    ///   0:(0,0,0)  1:(w,0,0)  2:(w,h,0)  3:(0,h,0)   <- front face (z=0)
    ///   4:(0,0,d)  5:(w,0,d)  6:(w,h,d)  7:(0,h,d)   <- back face  (z=d)
    fn create_box(width: f64, height: f64, depth: f64) -> Self {
        let (w, h, d) = (width, height, depth);

        let vertices = vec![
            Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0
            Vertex { point: DVec3::new(w,   0.0, 0.0) }, // 1
            Vertex { point: DVec3::new(w,   h,   0.0) }, // 2
            Vertex { point: DVec3::new(0.0, h,   0.0) }, // 3
            Vertex { point: DVec3::new(0.0, 0.0, d  ) }, // 4
            Vertex { point: DVec3::new(w,   0.0, d  ) }, // 5
            Vertex { point: DVec3::new(w,   h,   d  ) }, // 6
            Vertex { point: DVec3::new(0.0, h,   d  ) }, // 7
        ];

        // 12 edges: 4 front + 4 back + 4 lateral
        let edges = vec![
            Edge { start: 0, end: 1 }, // 0  front-bottom
            Edge { start: 1, end: 2 }, // 1  front-right
            Edge { start: 2, end: 3 }, // 2  front-top
            Edge { start: 3, end: 0 }, // 3  front-left
            Edge { start: 4, end: 5 }, // 4  back-bottom
            Edge { start: 5, end: 6 }, // 5  back-right
            Edge { start: 6, end: 7 }, // 6  back-top
            Edge { start: 7, end: 4 }, // 7  back-left
            Edge { start: 0, end: 4 }, // 8  lateral-bl
            Edge { start: 1, end: 5 }, // 9  lateral-br
            Edge { start: 2, end: 6 }, // 10 lateral-tr
            Edge { start: 3, end: 7 }, // 11 lateral-tl
        ];

        let faces = vec![
            // Front  (z=0, normal -Z)
            Face { outer_wire: Wire { edges: vec![WireEdge::fwd(0),WireEdge::fwd(1),WireEdge::fwd(2),WireEdge::fwd(3)] }, inner_wires: vec![], normal: DVec3::new(0.0, 0.0, -1.0), triangles: vec![[0,1,2],[0,2,3]] },
            // Back   (z=d, normal +Z)
            Face { outer_wire: Wire { edges: vec![WireEdge::fwd(4),WireEdge::fwd(5),WireEdge::fwd(6),WireEdge::fwd(7)] }, inner_wires: vec![], normal: DVec3::new(0.0, 0.0,  1.0), triangles: vec![[5,4,7],[5,7,6]] },
            // Bottom (y=0, normal -Y)
            Face { outer_wire: Wire { edges: vec![WireEdge::fwd(0),WireEdge::fwd(9),WireEdge::rev(4),WireEdge::rev(8)] }, inner_wires: vec![], normal: DVec3::new(0.0,-1.0, 0.0), triangles: vec![[0,1,5],[0,5,4]] },
            // Top    (y=h, normal +Y)
            Face { outer_wire: Wire { edges: vec![WireEdge::rev(2),WireEdge::fwd(10),WireEdge::fwd(6),WireEdge::rev(11)] }, inner_wires: vec![], normal: DVec3::new(0.0, 1.0, 0.0), triangles: vec![[3,2,6],[3,6,7]] },
            // Left   (x=0, normal -X)
            Face { outer_wire: Wire { edges: vec![WireEdge::rev(3),WireEdge::fwd(11),WireEdge::fwd(7),WireEdge::rev(8)] }, inner_wires: vec![], normal: DVec3::new(-1.0,0.0, 0.0), triangles: vec![[0,3,7],[0,7,4]] },
            // Right  (x=w, normal +X)
            Face { outer_wire: Wire { edges: vec![WireEdge::fwd(1),WireEdge::fwd(10),WireEdge::rev(5),WireEdge::rev(9)] }, inner_wires: vec![], normal: DVec3::new( 1.0,0.0, 0.0), triangles: vec![[1,2,6],[1,6,5]] },
        ];

        BRep {
            vertices,
            edges,
            solids: vec![Solid { shells: vec![Shell { faces }] }],
            geom: GeomStore::default(),
        }
    }
    ///
    /// Topology (OCCT-compatible single-seam representation):
    ///   Vertices: north (0, r, 0), south (0, -r, 0)
    ///   Edge E0:  seam meridian (north → south), Circle3 in XZ plane (normal = +Z)
    ///   Face F0:  SphericalSurface, outer_wire = [E0_fwd, E0_rev] (seam edge repeated)
    ///   PCurves:  E0 forward  → Line2d u=0,  v: 0 → π
    ///             E0 reversed → Line2d u=2π, v: π → 0
    fn create_sphere(radius: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;

        let r = radius;
        // Vertices
        let north = DVec3::new(0.0, r, 0.0);
        let south = DVec3::new(0.0, -r, 0.0);
        let vertices = vec![Vertex { point: north }, Vertex { point: south }];

        // Edge E0: seam meridian (north→south) — Circle3 in XZ plane
        // The seam lies at theta=0, i.e. x>0, z=0 plane.
        let edges = vec![Edge { start: 0, end: 1 }]; // E0

        // Face F0: outer_wire uses E0 twice (forward then reversed = seam)
        let face = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::rev(0)] },
            inner_wires: vec![],
            normal: DVec3::X, // outward, approximate
            triangles: vec![],
        };
        let shell = Shell { faces: vec![face] };
        let solid = Solid { shells: vec![shell] };

        // GeomStore
        let seam_curve = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: r,
        });
        let sphere_surf = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: r,
        });
        // PCurves for the seam edge on the sphere:
        //   Sphere param: u = longitude [0, 2π], v = colatitude [0, π] (phi from north pole)
        //   Forward half (north→south at u=0): Line2d origin=(0,0) dir=(0,1) extent π
        //   Reversed half (south→north at u=2π): Line2d origin=(2π,π) dir=(0,-1) extent π
        let pc_fwd = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        });
        let pc_rev = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(2.0 * PI, PI),
            direction: glam::DVec2::new(0.0, -1.0),
        });
        let geom = GeomStore {
            curves: vec![seam_curve],
            surfaces: vec![sphere_surf],
            curve2ds: vec![pc_fwd, pc_rev],
            edge_curve: vec![Some(0)],       // E0 → seam_curve
            face_surface: vec![Some(0)],     // F0 → sphere_surf
            edge_pcurves: vec![vec![         // E0: two pcurves (fwd and rev sides of seam)
                PCurve { surface_idx: 0, curve2d_idx: 0 },
                PCurve { surface_idx: 0, curve2d_idx: 1 },
            ]],
            // E0 is the half-meridian: t ∈ [0, π] on Circle3 (north→south)
            edge_curve_range: vec![Some([0.0, PI])],
            edge_degenerated: vec![false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
        };

        Self { vertices, edges, solids: vec![solid], geom }
    }

    /// Creates an analytic cylinder BRep along +Y axis, centered at origin.
    ///
    /// Topology:
    ///   Vertices: top_p (r, h/2, 0), bot_p (r, -h/2, 0)
    ///   Edges:
    ///     E0: top circle (Circle3, top_p → top_p seam, center=(0,h/2,0), normal=+Y)
    ///     E1: bot circle (Circle3, bot_p → bot_p seam, center=(0,-h/2,0), normal=-Y)
    ///     E2: seam line  (Line3,   top_p → bot_p)
    ///   Faces:
    ///     F0: CylindricalSurface, wire=[E2, E1_rev, E2_rev, E0]
    ///     F1: Plane +Y cap,       wire=[E0]
    ///     F2: Plane -Y cap,       wire=[E1_rev]  (stored as E1 with wire handling orientation)
    fn create_cylinder(radius: f64, height: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;

        let r = radius;
        let h = height;
        let half_h = h * 0.5;

        // Vertices: seam points where the circle meets the seam line (theta=0 → x=r, z=0)
        let top_p = DVec3::new(r, half_h, 0.0);
        let bot_p = DVec3::new(r, -half_h, 0.0);
        let vertices = vec![Vertex { point: top_p }, Vertex { point: bot_p }];

        // Edges
        // E0: top circle (seam: start == end == top_p = vertex 0)
        // E1: bot circle (seam: start == end == bot_p = vertex 1)
        // E2: seam line top→bot
        let edges = vec![
            Edge { start: 0, end: 0 }, // E0 top circle seam
            Edge { start: 1, end: 1 }, // E1 bot circle seam
            Edge { start: 0, end: 1 }, // E2 seam line
        ];

        // Faces
        // F0 lateral: outer_wire = [E2_fwd, E1_rev, E2_rev, E0_fwd]
        //   Cylinder lateral face: E2 fwd (top→bot seam down), E1 rev (bot circle CCW),
        //   E2 rev (bot→top seam up), E0 fwd (top circle CW seam).
        let f0 = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(2), WireEdge::rev(1), WireEdge::rev(2), WireEdge::fwd(0)] },
            inner_wires: vec![],
            normal: DVec3::X,
            triangles: vec![],
        };
        // F1 top cap: wire = [E0] forward
        let f1 = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(0)] },
            inner_wires: vec![],
            normal: DVec3::Y,
            triangles: vec![],
        };
        // F2 bottom cap: wire = [E1] reversed (bottom circle traversed CW from -Y view)
        let f2 = Face {
            outer_wire: Wire { edges: vec![WireEdge::rev(1)] },
            inner_wires: vec![],
            normal: -DVec3::Y,
            triangles: vec![],
        };
        let shell = Shell { faces: vec![f0, f1, f2] };
        let solid = Solid { shells: vec![shell] };

        // GeomStore — curves
        let top_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, half_h, 0.0),
            normal: DVec3::Y,
            radius: r,
        });
        let bot_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
            radius: r,
        });
        let seam_line = Curve3::Line(Line3 {
            origin: top_p,
            direction: (bot_p - top_p).normalize(),
        });

        // GeomStore — surfaces
        // F0: CylindricalSurface (origin at bottom center, axis=+Y)
        let cyl_surf = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(0.0, -half_h, 0.0),
            axis: DVec3::Y,
            radius: r,
        });
        let top_plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, half_h, 0.0),
            normal: DVec3::Y,
        });
        let bot_plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
        });

        // PCurves
        // Cylinder param: u=azimuth [0,2π], v=height [0,h] (v=0 at bottom, v=h at top)
        // E0 (top circle) on F0 (cyl): iso-line v=h, u from 0 to 2π (but seam: u=2π→0)
        //   We use the convention u goes from 2π to 0 for the "reversed" seam orientation.
        //   Forward: Line2d (2π, h) dir=(-1, 0)  [u decreases at v=h]
        // E0 (top circle) on F1 (top plane): circle in plane param space
        let e0_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(2.0 * PI, h),
            direction: glam::DVec2::new(-1.0, 0.0),
        });
        let e0_on_f1 = Curve2d::Circle(Circle2d {
            center: glam::DVec2::ZERO,
            radius: r,
        });
        // E1 (bot circle) on F0 (cyl): iso-line v=0, u from 0 to 2π
        let e1_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(1.0, 0.0),
        });
        let e1_on_f2 = Curve2d::Circle(Circle2d {
            center: glam::DVec2::ZERO,
            radius: r,
        });
        // E2 (seam line) on F0 (cyl): iso-line u=0, v from h to 0
        let e2_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, h),
            direction: glam::DVec2::new(0.0, -1.0),
        });

        let geom = GeomStore {
            curves: vec![top_circle, bot_circle, seam_line], // curve 0,1,2
            surfaces: vec![cyl_surf, top_plane, bot_plane],   // surface 0,1,2
            curve2ds: vec![e0_on_f0, e0_on_f1, e1_on_f0, e1_on_f2, e2_on_f0], // 0..4
            edge_curve: vec![Some(0), Some(1), Some(2)],      // E0→0, E1→1, E2→2
            face_surface: vec![Some(0), Some(1), Some(2)],    // F0→cyl, F1→top, F2→bot
            edge_pcurves: vec![
                // E0: top circle → on cyl (idx 0) and top plane (idx 1)
                vec![
                    PCurve { surface_idx: 0, curve2d_idx: 0 }, // E0 on cyl
                    PCurve { surface_idx: 1, curve2d_idx: 1 }, // E0 on top plane
                ],
                // E1: bot circle → on cyl (idx 0) and bot plane (idx 2)
                vec![
                    PCurve { surface_idx: 0, curve2d_idx: 2 }, // E1 on cyl
                    PCurve { surface_idx: 2, curve2d_idx: 3 }, // E1 on bot plane
                ],
                // E2: seam line → on cyl only
                vec![
                    PCurve { surface_idx: 0, curve2d_idx: 4 }, // E2 on cyl
                ],
            ],
            // E0/E1: full circles [0, 2π]; E2: seam line [0, h]
            edge_curve_range: vec![
                Some([0.0, 2.0 * PI]),  // E0 top circle
                Some([0.0, 2.0 * PI]),  // E1 bot circle
                Some([0.0, h]),         // E2 seam line
            ],
            edge_degenerated: vec![false, false, false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
        };

        Self { vertices, edges, solids: vec![solid], geom }
    }

    /// Creates an analytic cone BRep along +Y axis, apex at +Y, centered at origin.
    ///
    /// Topology:
    ///   Vertices: apex (0, h/2, 0), base_p (R, -h/2, 0)
    ///   Edges:
    ///     E0: base circle (Circle3, base_p → base_p seam, normal=-Y)
    ///     E1: slant line  (Line3,   apex → base_p)
    ///   Faces:
    ///     F0: ConicalSurface, wire=[E1, E0, E1_rev]  (seam)
    ///     F1: Plane -Y cap,  wire=[E0]
    fn create_cone(base_radius: f64, height: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;
        let r = base_radius;
        let h = height;
        let half_h = h * 0.5;

        let apex_pt  = DVec3::new(0.0,  half_h, 0.0);
        let base_pt  = DVec3::new(r,   -half_h, 0.0);
        let vertices = vec![Vertex { point: apex_pt }, Vertex { point: base_pt }];

        let edges = vec![
            Edge { start: 1, end: 1 }, // E0 base circle seam
            Edge { start: 0, end: 1 }, // E1 slant line apex→base_p
        ];

        // F0 conical lateral: E1 fwd (apex→base seam), E0 fwd (base circle), E1 rev (base→apex seam)
        let f0 = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(1), WireEdge::fwd(0), WireEdge::rev(1)] },
            inner_wires: vec![],
            normal: DVec3::X,
            triangles: vec![],
        };
        // F1 base cap: E0 reversed (base circle CW from -Y view)
        let f1 = Face {
            outer_wire: Wire { edges: vec![WireEdge::rev(0)] },
            inner_wires: vec![],
            normal: -DVec3::Y,
            triangles: vec![],
        };
        let shell = Shell { faces: vec![f0, f1] };
        let solid = Solid { shells: vec![shell] };

        // Curves
        let base_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
            radius: r,
        });
        let slant = Curve3::Line(Line3 {
            origin: apex_pt,
            direction: (base_pt - apex_pt).normalize(),
        });

        // half-angle = atan(R / h)
        let half_angle = (r / h).atan();
        // ConicalSurface: apex at top, axis pointing down (-Y), radius at apex = 0
        let cone_surf = Surface3::Cone(geom::ConicalSurface {
            apex: apex_pt,
            axis: -DVec3::Y,
            radius: 0.0, // radius at apex
            half_angle_rad: half_angle,
        });
        let base_plane = Surface3::Plane(Plane {
            origin: DVec3::new(0.0, -half_h, 0.0),
            normal: -DVec3::Y,
        });

        // PCurves
        // Cone param: u=azimuth [0,2π], v=slant distance from apex
        let slant_len = (r * r + h * h).sqrt();
        // E0 (base circle) on F0 (cone): iso-line at v=slant_len, u from 0 to 2π
        let e0_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, slant_len),
            direction: glam::DVec2::new(1.0, 0.0),
        });
        let e0_on_f1 = Curve2d::Circle(Circle2d {
            center: glam::DVec2::ZERO,
            radius: r,
        });
        // E1 (slant seam) on F0: iso-line u=0, v from 0 to slant_len
        let e1_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        });

        let geom = GeomStore {
            curves: vec![base_circle, slant],
            surfaces: vec![cone_surf, base_plane],
            curve2ds: vec![e0_on_f0, e0_on_f1, e1_on_f0],
            edge_curve: vec![Some(0), Some(1)],
            face_surface: vec![Some(0), Some(1)],
            edge_pcurves: vec![
                // E0: base circle → on cone and base plane
                vec![
                    PCurve { surface_idx: 0, curve2d_idx: 0 },
                    PCurve { surface_idx: 1, curve2d_idx: 1 },
                ],
                // E1: slant line → on cone only
                vec![
                    PCurve { surface_idx: 0, curve2d_idx: 2 },
                ],
            ],
            // E0: full base circle [0, 2π]; E1: slant from apex to base [0, slant_len]
            edge_curve_range: vec![
                Some([0.0, 2.0 * PI]),    // E0 base circle
                Some([0.0, slant_len]),   // E1 slant line
            ],
            edge_degenerated: vec![false, false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
        };

        Self { vertices, edges, solids: vec![solid], geom }
    }

    /// Creates an analytic torus BRep around +Y axis, centered at origin.
    ///
    /// Topology (double-seam representation):
    ///   Vertex: seam_pt = (R+r, 0, 0) — intersection of major and minor seams
    ///   Edges:
    ///     E0: major seam circle (major radius R, center=origin, normal=+Y)
    ///     E1: minor seam circle (minor radius r, in XY plane at (R,0,0))
    ///   Face:
    ///     F0: ToroidalSurface, outer_wire=[E0, E1, E0_rev, E1_rev]
    ///   PCurves:
    ///     E0 on F0: Line2d (0,0)→(2π,0)   [major seam: u from 0..2π at v=0]
    ///     E1 on F0: Line2d (0,0)→(0,2π)   [minor seam: v from 0..2π at u=0]
    fn create_torus(major_radius: f64, minor_radius: f64) -> Self {
        use geom::*;
        use std::f64::consts::PI;

        let big_r = major_radius;
        let small_r = minor_radius;

        // Single vertex at the seam intersection
        let seam_pt = DVec3::new(big_r + small_r, 0.0, 0.0);
        let vertices = vec![Vertex { point: seam_pt }];

        // E0: major seam (full circle of radius R in XZ plane)
        // E1: minor seam (full circle of radius r in a plane through the tube)
        let edges = vec![
            Edge { start: 0, end: 0 }, // E0 major circle seam
            Edge { start: 0, end: 0 }, // E1 minor circle seam
        ];

        // F0: outer_wire — E0 fwd (major seam), E1 fwd (minor seam),
        //     E0 rev (major seam reversed), E1 rev (minor seam reversed)
        let face = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::rev(0), WireEdge::rev(1)] },
            inner_wires: vec![],
            normal: DVec3::X,
            triangles: vec![],
        };
        let shell = Shell { faces: vec![face] };
        let solid = Solid { shells: vec![shell] };

        // Curves
        let major_circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Y,
            radius: big_r,
        });
        // Minor circle: centered at (R,0,0), in the YZ plane (normal = +X)
        let minor_circle = Curve3::Circle(Circle3 {
            center: DVec3::new(big_r, 0.0, 0.0),
            normal: DVec3::X,
            radius: small_r,
        });

        let torus_surf = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: big_r,
            minor_radius: small_r,
        });

        // PCurves — torus param: u=major angle [0,2π], v=minor angle [0,2π]
        let e0_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(1.0, 0.0),
        });
        let e0_on_f0_rev = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(2.0 * PI, 0.0),
            direction: glam::DVec2::new(-1.0, 0.0),
        });
        let e1_on_f0 = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        });
        let e1_on_f0_rev = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 2.0 * PI),
            direction: glam::DVec2::new(0.0, -1.0),
        });

        let geom = GeomStore {
            curves: vec![major_circle, minor_circle],
            surfaces: vec![torus_surf],
            curve2ds: vec![e0_on_f0, e0_on_f0_rev, e1_on_f0, e1_on_f0_rev],
            edge_curve: vec![Some(0), Some(1)],
            face_surface: vec![Some(0)],
            edge_pcurves: vec![
                // E0: major seam — two pcurves (forward and reverse passes)
                vec![
                    PCurve { surface_idx: 0, curve2d_idx: 0 },
                    PCurve { surface_idx: 0, curve2d_idx: 1 },
                ],
                // E1: minor seam — two pcurves
                vec![
                    PCurve { surface_idx: 0, curve2d_idx: 2 },
                    PCurve { surface_idx: 0, curve2d_idx: 3 },
                ],
            ],
            // Both seams are full circles [0, 2π]
            edge_curve_range: vec![
                Some([0.0, 2.0 * PI]),  // E0 major seam circle
                Some([0.0, 2.0 * PI]),  // E1 minor seam circle
            ],
            edge_degenerated: vec![false, false],
            vertex_tolerance: Vec::new(),
            edge_tolerance: Vec::new(),
            face_tolerance: Vec::new(),
        };

        Self { vertices, edges, solids: vec![solid], geom }
    }

    /// Materializes a primitive solid descriptor into an analytic B-Rep.
    ///
    /// The resulting BRep has fully populated GeomStore entries (edge_curve,
    /// face_surface, edge_pcurves). Triangles are NOT pre-populated; the render
    /// layer tessellates on demand.
    ///
    /// User-facing code should prefer `rcad-modeling` construction helpers.
    pub fn from_primitive(primitive: PrimitiveSolid) -> Self {
        match primitive {
            PrimitiveSolid::Box { width, height, depth } => Self::create_box(width, height, depth),
            PrimitiveSolid::Sphere { radius } => Self::create_sphere(radius),
            PrimitiveSolid::Cylinder { radius, height } => Self::create_cylinder(radius, height),
            PrimitiveSolid::Cone { base_radius, height } => Self::create_cone(base_radius, height),
            PrimitiveSolid::Torus { major_radius, minor_radius } => {
                Self::create_torus(major_radius, minor_radius)
            }
        }
    }

    pub fn center(&self) -> DVec3 {
        if self.vertices.is_empty() {
            return DVec3::ZERO;
        }
        let mut sum = DVec3::ZERO;
        for v in &self.vertices {
            sum += v.point;
        }
        sum / self.vertices.len() as f64
    }

    /// Returns the axis-aligned bounding box of all vertices as `[min, max]`,
    /// or `None` if the BRep has no vertices.
    pub fn bounding_box(&self) -> Option<[DVec3; 2]> {
        if self.vertices.is_empty() {
            return None;
        }
        let mut mn = DVec3::splat(f64::INFINITY);
        let mut mx = DVec3::splat(f64::NEG_INFINITY);
        for v in &self.vertices {
            mn = mn.min(v.point);
            mx = mx.max(v.point);
        }
        Some([mn, mx])
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    fn geom_populated(brep: &BRep) -> bool {
        !brep.geom.face_surface.is_empty()
            && brep.geom.face_surface.iter().all(|s| s.is_some())
            && !brep.geom.edge_pcurves.is_empty()
    }

    #[test]
    fn creates_sphere_with_analytic_geom() {
        let brep = BRep::create_sphere(1.0);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(matches!(
            brep.geom.surfaces.first(),
            Some(Surface3::Sphere(_))
        ));
        // triangles are empty — render layer tessellates on demand
        assert!(brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).all(|f| f.triangles.is_empty()));
    }

    #[test]
    fn creates_cylinder_with_analytic_geom() {
        let brep = BRep::create_cylinder(1.0, 2.0);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(brep.geom.surfaces.iter().any(|s| matches!(s, Surface3::Cylinder(_))));
    }

    #[test]
    fn creates_cone_with_analytic_geom() {
        let brep = BRep::create_cone(1.0, 2.0);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(brep.geom.surfaces.iter().any(|s| matches!(s, Surface3::Cone(_))));
    }

    #[test]
    fn creates_torus_with_analytic_geom() {
        let brep = BRep::create_torus(1.0, 0.3);
        assert!(!brep.vertices.is_empty());
        assert!(geom_populated(&brep));
        assert!(brep.geom.surfaces.iter().any(|s| matches!(s, Surface3::Torus(_))));
    }
}
