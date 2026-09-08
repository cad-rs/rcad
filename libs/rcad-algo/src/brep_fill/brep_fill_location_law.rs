//! OCCT BRepFill_LocationLaw (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_LocationLaw.hxx (L37-142) + BRepFill_LocationLaw.cxx (whole file
//! L45-745, including the file statics Norm / ToG0).
//!
//! Architecture differences:
//! - `BRepTools_WireExplorer` maps to the wire-ordered `TWireData::edges`
//!   list; `BRep_Tool::Degenerated / Pnt / Tolerance` map to the stored
//!   `TEdgeData / TVertexData` fields — the owning `BRep` pool is passed as
//!   the leading `brep` argument of the methods that touch topology.
//! - `NCollection_HArray1` (1-based) maps to `Vec` with index-1 access
//!   (`Value(i)` -> `vec[i-1]`).
//! - `handle(GeomFill_LocationLaw)` maps to
//!   `Rc<RefCell<dyn geomalgo::geomfill LocationLaw>>` (the laws are shared
//!   with the callers of `Law(i)` and mutated in place through the handle).
//! - `myPath.Closed()` maps to the `tshape_flags::CLOSED` flag of the wire.
//! - `CurvilinearBounds` is declared const in OCCT but mutates `myLength`
//!   through the array handle; the rcad form takes `&mut self`.
//! - `gp_Trsf` maps to `glam::DAffine3` (the D0 carrier construction goes
//!   through the re-hosted `gp_trsf_set_values` with the OCCT
//!   orthogonalization).
//! - GAP carrier: `BRepBuilderAPI_Transform(W, fila, true)` (TKTopAlgo) at
//!   the D0 call site (cxx L713) keeps the OCCT failure path.

use std::cell::RefCell;
use std::rc::Rc;

use glam::{DVec2, DVec3, DAffine3};

use rcad_kernel::base::gcpnts::abscissa_point::{abscissa_point_parameter, arc_length};
use rcad_kernel::topo::topods::{
    tshape_flags, BRep, BRepBuilder, Orientation, Shape,
};

use crate::brep_fill::compatible_wires::{
    angle_with_ref, shape_tolerance, top_exp_first_vertex_stored, top_exp_last_vertex_stored,
    wire_edges,
};
use crate::brep_fill::generator::top_exp_wire_vertices;
use crate::geomalgo::geomfill::gp_mat::GpMat;
use crate::geomalgo::geomfill::location_law::LocationLaw;
use crate::geomalgo::geomfill::trihedron_law::{
    curve_first_parameter, curve_last_parameter, PipeError,
};

// ---------------------------------------------------------------------------
// Pure-math re-hosts (gp statics consumed by this file)
// ---------------------------------------------------------------------------

/// OCCT gp_Mat::Invert (gp_Mat.cxx, optimized adjugate / determinant form).
fn gp_mat_invert(m: &mut GpMat) {
    let a00 = m.mat[0][0];
    let a01 = m.mat[0][1];
    let a02 = m.mat[0][2];
    let a10 = m.mat[1][0];
    let a11 = m.mat[1][1];
    let a12 = m.mat[1][2];
    let a20 = m.mat[2][0];
    let a21 = m.mat[2][1];
    let a22 = m.mat[2][2];

    // Compute adjugate matrix (transpose of cofactor matrix).
    let adj00 = a11 * a22 - a12 * a21;
    let adj10 = a12 * a20 - a10 * a22;
    let adj20 = a10 * a21 - a11 * a20;
    let adj01 = a02 * a21 - a01 * a22;
    let adj11 = a00 * a22 - a02 * a20;
    let adj21 = a01 * a20 - a00 * a21;
    let adj02 = a01 * a12 - a02 * a11;
    let adj12 = a02 * a10 - a00 * a12;
    let adj22 = a00 * a11 - a01 * a10;

    let a_det = a00 * adj00 + a01 * adj10 + a02 * adj20;
    assert!(
        a_det.abs() > f64::MIN_POSITIVE,
        "gp_Mat::Invert() - matrix has zero determinant"
    );
    let inv_det = 1.0 / a_det;

    m.mat[0][0] = adj00 * inv_det;
    m.mat[1][0] = adj10 * inv_det;
    m.mat[2][0] = adj20 * inv_det;
    m.mat[0][1] = adj01 * inv_det;
    m.mat[1][1] = adj11 * inv_det;
    m.mat[2][1] = adj21 * inv_det;
    m.mat[0][2] = adj02 * inv_det;
    m.mat[1][2] = adj12 * inv_det;
    m.mat[2][2] = adj22 * inv_det;
}

/// OCCT gp_Mat::Inverted().
pub(super) fn gp_mat_inverted(m: &GpMat) -> GpMat {
    let mut new_mat = *m;
    gp_mat_invert(&mut new_mat);
    new_mat
}

/// OCCT gp_Mat operator- (coefficient-wise difference).
fn gp_mat_minus(a: &GpMat, b: &GpMat) -> GpMat {
    let mut r = *a;
    for i in 0..3 {
        for j in 0..3 {
            r.mat[i][j] -= b.mat[i][j];
        }
    }
    r
}

/// OCCT static Norm (BRepFill_LocationLaw.cxx L47-67) — the maximal
/// sum-of-abs row norm.
fn norm_mat(m: &GpMat) -> f64 {
    let mut norme;
    let mut coord = DVec3::new(m.mat[0][0], m.mat[0][1], m.mat[0][2]);
    norme = coord.x.abs() + coord.y.abs() + coord.z.abs();
    coord = DVec3::new(m.mat[1][0], m.mat[1][1], m.mat[1][2]);
    let r = coord.x.abs() + coord.y.abs() + coord.z.abs();
    if r > norme {
        norme = r;
    }
    coord = DVec3::new(m.mat[2][0], m.mat[2][1], m.mat[2][2]);
    let r = coord.x.abs() + coord.y.abs() + coord.z.abs();
    if r > norme {
        norme = r;
    }
    norme
}

/// OCCT static ToG0 (L69-78) — Calculate transformation T such as T.M2 = M1.
pub(super) fn to_g0(m1: &GpMat, m2: &GpMat, t: &mut GpMat) {
    *t = gp_mat_inverted(m2);
    *t = t.multiplied_mat(m1);
}

/// OCCT gp_Mat::Column(theCol) — the column read back as an XYZ.
pub(super) fn gp_mat_column(m: &GpMat, the_col: usize) -> DVec3 {
    DVec3::new(
        m.mat[0][the_col - 1],
        m.mat[1][the_col - 1],
        m.mat[2][the_col - 1],
    )
}

/// OCCT gp_Dir::IsParallel(Other, AngularTolerance) — the angle is 0 or PI
/// within the angular tolerance (gp_Dir.cxx).
fn dir_is_parallel(a: DVec3, b: DVec3, ang_tol: f64) -> bool {
    let la = a.length();
    let lb = b.length();
    if la < 1e-300 || lb < 1e-300 {
        return false;
    }
    let cross = a.cross(b).length() / (la * lb);
    cross <= ang_tol.abs().sin() + 1e-18
}

/// OCCT gp_Dir::IsOpposite(Other, AngularTolerance) — parallel with a
/// negative dot product.
fn dir_is_opposite(a: DVec3, b: DVec3, ang_tol: f64) -> bool {
    dir_is_parallel(a, b, ang_tol) && a.dot(b) < 0.0
}

/// OCCT gp_Vec::Rotate (the gp_Trsf rotation form) — Rodrigues around the
/// unit axis direction.
fn vec_rotated(v: DVec3, axis_dir: DVec3, ang: f64) -> DVec3 {
    let u = axis_dir.normalize_or_zero();
    let c = ang.cos();
    let s = ang.sin();
    let d = v.dot(u);
    v * c + u.cross(v) * s + u * (d * (1.0 - c))
}

/// OCCT gp_Trsf::Orthogonalize (gp_Trsf.cxx L863-945) — Gram-Schmidt over
/// the columns then over the rows (the geomfill/sweep.rs re-host).
fn trsf_orthogonalize(m: &GpMat) -> GpMat {
    let mut a_tm = *m;

    let mut v1 = gp_mat_column(&a_tm, 1);
    let mut v2 = gp_mat_column(&a_tm, 2);
    let mut v3 = gp_mat_column(&a_tm, 3);

    v1 = v1.normalize_or_zero();
    v2 -= v1 * (v2.dot(v1));
    v2 = v2.normalize_or_zero();
    v3 -= v1 * (v3.dot(v1)) + v2 * (v3.dot(v2));
    v3 = v3.normalize_or_zero();
    a_tm.set_cols(v1, v2, v3);

    let row = |t: &GpMat, i: usize| DVec3::new(t.mat[i - 1][0], t.mat[i - 1][1], t.mat[i - 1][2]);
    let mut r1 = row(&a_tm, 1);
    let mut r2 = row(&a_tm, 2);
    let mut r3 = row(&a_tm, 3);

    r1 = r1.normalize_or_zero();
    r2 -= r1 * (r2.dot(r1));
    r2 = r2.normalize_or_zero();
    r3 -= r1 * (r3.dot(r1)) + r2 * (r3.dot(r2));
    r3 = r3.normalize_or_zero();
    a_tm.mat[0] = [r1.x, r1.y, r1.z];
    a_tm.mat[1] = [r2.x, r2.y, r2.z];
    a_tm.mat[2] = [r3.x, r3.y, r3.z];

    a_tm
}

/// OCCT gp_Trsf::SetValues(a11..a34) (gp_Trsf.cxx L346-385) — the null
/// determinant raise is the Err branch; the vectorial part is scaled by
/// det^(1/3) and orthogonalized, loc = col4 (the geomfill/sweep.rs re-host).
fn gp_trsf_set_values(a: &[f64; 12]) -> Result<DAffine3, ()> {
    let m = GpMat::from_rows(a[0], a[1], a[2], a[4], a[5], a[6], a[8], a[9], a[10]);
    let mut s = m.mat[0][0] * (m.mat[1][1] * m.mat[2][2] - m.mat[1][2] * m.mat[2][1])
        - m.mat[0][1] * (m.mat[1][0] * m.mat[2][2] - m.mat[1][2] * m.mat[2][0])
        + m.mat[0][2] * (m.mat[1][0] * m.mat[2][1] - m.mat[1][1] * m.mat[2][0]);
    if s.abs() < f64::MIN_POSITIVE {
        return Err(()); // "gp_Trsf::SetValues, null determinant"
    }
    if s > 0.0 {
        s = s.powf(1.0 / 3.0);
    } else {
        s = -((-s).powf(1.0 / 3.0));
    }
    let mut scaled = m.multiplied_scalar(1.0 / s);
    scaled = trsf_orthogonalize(&scaled);

    let mut trsf = DAffine3::from_mat3(glam::DMat3::from_cols(
        gp_mat_column(&scaled, 1),
        gp_mat_column(&scaled, 2),
        gp_mat_column(&scaled, 3),
    ));
    trsf.translation = DVec3::new(a[3], a[7], a[11]);
    Ok(trsf)
}

/// OCCT gp_Vec::Angle (gp_XYZ::Angle) — the angle between two vectors.
fn gv_angle(v1: DVec3, v2: DVec3) -> f64 {
    let an_norm = v1.length();
    let a_no_norm = v2.length();
    let mut value = v1.dot(v2) / (an_norm * a_no_norm);
    if value > 1.0 {
        value = 1.0;
    } else if value < -1.0 {
        value = -1.0;
    }
    value.acos()
}

/// GAP carrier: OCCT BRepBuilderAPI_Transform(S, T, Copy) (TKTopAlgo,
/// BRepBuilderAPI_Transform.cxx) — the shape-transform engine is not
/// translated (plan D3); the D0 call site (BRepFill_LocationLaw.cxx L713)
/// keeps the OCCT failure path.
fn brep_builder_api_transform_copy(_w: &Shape, _fila: &DAffine3, _copy: bool) -> Shape {
    panic!(
        "GAP: BRepBuilderAPI_Transform (TKTopAlgo) is not translated — \
         BRepFill_LocationLaw::D0 (BRepFill_LocationLaw.cxx L713)"
    )
}

// ---------------------------------------------------------------------------
// BRepFill_LocationLaw
// ---------------------------------------------------------------------------

/// OCCT `handle<BRepFill_LocationLaw>` slot trait.  The OCCT derived classes
/// (BRepFill_Edge3DLaw / BRepFill_ACRLaw / BRepFill_EdgeOnSurfLaw /
/// BRepFill_DraftLaw) override no base method; the rcad derived structs embed
/// the [`BRepFillLocationLaw`] base and forward this accessor, so the
/// consumers (BRepFill_Sweep / BRepFill_SectionPlacement) stay 1:1 with the
/// OCCT handle usage (`myLoc->Law(i)` -> `h.borrow().base().law(i)`).
pub trait BRepFillLocationLawOps {
    /// The embedded base sub-object.
    fn base(&self) -> &BRepFillLocationLaw;
    /// The embedded base sub-object (mutable).
    fn base_mut(&mut self) -> &mut BRepFillLocationLaw;
}

/// OCCT BRepFill_LocationLaw (hxx L37-142): Location Law on a Wire.
pub struct BRepFillLocationLaw {
    /// OCCT myPath (TopoDS_Wire).
    pub my_path: Shape,
    /// OCCT myTol.
    pub my_tol: f64,
    /// OCCT myLaws (HArray1(1, NbEdge) of handle(GeomFill_LocationLaw)) —
    /// the 1-based array maps to a Vec (`Value(i)` -> `my_laws[i-1]`).
    pub my_laws: Vec<Rc<RefCell<dyn LocationLaw>>>,
    /// OCCT myLength (HArray1(1, NbEdge + 1)).
    pub my_length: Vec<f64>,
    /// OCCT myEdges (HArray1(1, NbEdge) of TopoDS_Shape).
    pub my_edges: Vec<Shape>,
    /// OCCT myDisc (HArray1 of int), null when not computed.
    pub my_disc: Option<Vec<i32>>,
    /// OCCT myType (private).
    pub my_type: i32,
}

impl BRepFillLocationLaw {
    /// OCCT Init (cxx L82-110) — Initialize all the fields, this method has
    /// to be called by the constructors of the inherited classes.
    pub fn init(&mut self, brep: &BRep, path: &Shape) {
        // L94-101: count the non-degenerated edges.
        let mut nb_edge = 0usize;
        for e in wire_edges(brep, path) {
            let degenerated = brep.edge(e.clone()).degenerated;
            if !degenerated {
                nb_edge += 1;
            }
        }

        self.my_path = path.clone();
        self.my_tol = 1.0e-4;

        self.my_laws = Vec::new();
        // myLength = new HArray1(1, NbEdge + 1); Init(-1.); SetValue(1, 0.)
        self.my_length = vec![-1.0; nb_edge + 1];
        self.my_length[0] = 0.0;
        self.my_edges = vec![Shape::null(); nb_edge];
        self.my_disc = None;
        self.tangent_is_main();
    }

    /// OCCT myPath.Closed() — the wire CLOSED flag.
    fn path_is_closed(&self, brep: &BRep) -> bool {
        brep.has_flag(self.my_path.clone(), tshape_flags::CLOSED)
    }

    /// OCCT GetStatus (cxx L114-123).
    pub fn get_status(&self) -> PipeError {
        let mut status = PipeError::PipeOk;
        for law in &self.my_laws {
            if status != PipeError::PipeOk {
                break;
            }
            status = law.borrow().error_status();
        }
        status
    }

    /// OCCT TangentIsMain (cxx L127-130).
    pub fn tangent_is_main(&mut self) {
        self.my_type = 1;
    }

    /// OCCT NormalIsMain (cxx L134-137).
    pub fn normal_is_main(&mut self) {
        self.my_type = 2;
    }

    /// OCCT BiNormalIsMain (cxx L141-144).
    pub fn binormal_is_main(&mut self) {
        self.my_type = 3;
    }

    /// OCCT TransformInCompatibleLaw (cxx L148-198) — Apply a linear
    /// transformation on each law, to reduce the discontinuities of law at
    /// one rotation.
    pub fn transform_in_compatible_law(&mut self, tol_angular: f64) {
        let mut first = 0.0;
        let mut last = 0.0;
        let mut angle;
        let mut trsf = GpMat::identity();
        let mut m1 = GpMat::identity();
        let mut m2 = GpMat::identity();
        let mut v = DVec3::ZERO;
        let mut t1;
        let mut t2;
        let mut n1;
        let mut n2;
        let oz = DVec3::new(0.0, 0.0, 1.0); // OCCT gp_XYZ OZ(0, 0, 1)

        self.my_laws[0].borrow().get_domain(&mut first, &mut last);

        let mut ipath = 2usize;
        while ipath <= self.my_laws.len() {
            self.my_laws[ipath - 2].borrow().d0(last, &mut m1, &mut v);
            self.my_laws[ipath - 1]
                .borrow()
                .get_domain(&mut first, &mut last);
            self.my_laws[ipath - 1].borrow().d0(first, &mut m2, &mut v);
            t1 = gp_mat_column(&m1, 3); // T1.SetXYZ(M1.Column(3))
            t2 = gp_mat_column(&m2, 3);
            n1 = gp_mat_column(&m1, 1);
            n2 = gp_mat_column(&m2, 1);
            if dir_is_parallel(t1, t2, tol_angular) && !dir_is_opposite(t1, t2, tol_angular) {
                // Correction G0
                to_g0(&m1, &m2, &mut trsf);
            } else {
                let cross = t1.cross(t2); // gp_Vec cross(T1); cross.Cross(T2)
                let alpha = angle_with_ref(t2, t1, cross);
                // N2.Rotate(axe, alpha) with axe(gp::Origin(), cross.XYZ())
                n2 = vec_rotated(n2, cross, alpha);

                angle = angle_with_ref(n2, n1, t1);
                trsf.set_rotation(oz, angle);
            }
            self.my_laws[ipath - 1].borrow_mut().set_trsf(trsf);
            ipath += 1;
        }
    }

    /// OCCT TransformInG0Law (cxx L202-226).
    pub fn transform_in_g0_law(&mut self, brep: &BRep) {
        let mut first = 0.0;
        let mut last = 0.0;
        let mut aux = GpMat::identity();
        let mut m1 = GpMat::identity();
        let mut m2 = GpMat::identity();
        let mut v = DVec3::ZERO;

        self.my_laws[0].borrow().get_domain(&mut first, &mut last);
        let mut ipath = 2usize;
        while ipath <= self.my_laws.len() {
            self.my_laws[ipath - 2].borrow().d0(last, &mut m1, &mut v);
            self.my_laws[ipath - 1]
                .borrow()
                .get_domain(&mut first, &mut last);
            self.my_laws[ipath - 1].borrow().d0(first, &mut m2, &mut v);
            to_g0(&m1, &m2, &mut aux);
            self.my_laws[ipath - 1].borrow_mut().set_trsf(aux);
            ipath += 1;
        }

        // Is the law periodical ?
        // (the closed-path junction evaluation of cxx L219-225)
        if self.path_is_closed(brep) {
            let n = self.my_laws.len() - 1;
            self.my_laws[n].borrow().d0(last, &mut m1, &mut v);
            self.my_laws[0].borrow().get_domain(&mut first, &mut last);
            self.my_laws[0].borrow().d0(first, &mut m2, &mut v);
        }
    }

    /// OCCT DeleteTransform (cxx L232-241) — Remove the setting in
    /// continuity of law.
    pub fn delete_transform(&mut self) {
        let id = GpMat::identity();
        for law in &self.my_laws {
            law.borrow_mut().set_trsf(id);
        }
        self.my_disc = None;
    }

    /// OCCT NbHoles (cxx L247-275) — Find "Holes".
    /// (the `brep` pool argument is the rcad BRep_Tool mapping — see the
    /// file-header architecture note)
    pub fn nb_holes(&mut self, brep: &BRep, tol: f64) -> i32 {
        if self.my_disc.is_none() {
            let mut seq: Vec<i32> = Vec::new();
            let mut ii = 2i32;
            while ii <= self.my_laws.len() as i32 + 1 {
                if self.is_g1(brep, ii - 1, tol, 1.0e-12) == -1 {
                    seq.push(ii);
                }
                ii += 1;
            }
            let nb_disc = seq.len() as i32;
            if nb_disc > 0 {
                self.my_disc = Some(seq);
            }
        }
        match &self.my_disc {
            None => 0,
            Some(disc) => disc.len() as i32,
        }
    }

    /// OCCT Holes (cxx L279-288).
    pub fn holes(&self, disc: &mut [i32]) {
        if let Some(my_disc) = &self.my_disc {
            for (ii, value) in my_disc.iter().enumerate() {
                disc[ii] = *value;
            }
        }
    }

    /// OCCT NbLaw (cxx L292-295) — Return the number of elementary Law.
    pub fn nb_law(&self) -> i32 {
        self.my_laws.len() as i32
    }

    /// OCCT Law (cxx L299-302) — Return the elementary Law of rank <Index>
    /// (<Index> have to be in [1, NbLaw()]).
    pub fn law(&self, index: i32) -> Rc<RefCell<dyn LocationLaw>> {
        self.my_laws[(index - 1) as usize].clone()
    }

    /// OCCT Wire (cxx L306-309) — return the path.
    pub fn wire(&self) -> Shape {
        self.my_path.clone()
    }

    /// OCCT Edge (cxx L313-316) — Return the Edge of rank <Index> in the
    /// path (TopoDS::Edge cast).
    pub fn edge(&self, index: i32) -> Shape {
        self.my_edges[(index - 1) as usize].clone()
    }

    /// OCCT Vertex (cxx L320-349) — Return the vertex of rank <Index> in the
    /// path (<Index> have to be in [0, NbLaw()]).
    pub fn vertex(&self, brep: &BRep, index: i32) -> Shape {
        let mut v = Shape::null();
        if index <= self.my_edges.len() as i32 {
            let e = self.my_edges[(index - 1) as usize].clone();
            if e.orientation == Orientation::Reversed {
                v = top_exp_last_vertex_stored(brep, &e);
            } else {
                v = top_exp_first_vertex_stored(brep, &e);
            }
        } else if index == self.my_edges.len() as i32 + 1 {
            let e = self.my_edges[(index - 2) as usize].clone();
            if e.orientation == Orientation::Reversed {
                v = top_exp_first_vertex_stored(brep, &e);
            } else {
                v = top_exp_last_vertex_stored(brep, &e);
            }
        }
        v
    }

    /// OCCT PerformVertex (cxx L356-452) — Calculate a vertex of sweeping
    /// from a vertex of section and the index of the edge in the trajectory.
    pub fn perform_vertex(
        &self,
        brep: &mut BRep,
        index: i32,
        input: &Shape,
        tol_min: f64,
        output: &mut Shape,
        iloc: i32,
    ) {
        let mut b = BRepBuilder::new();
        let mut is_bary = iloc == 0;
        let mut first = 0.0;
        let mut last = 0.0;
        let mut m1 = GpMat::identity();
        let mut m2 = GpMat::identity();
        let mut v1 = DVec3::ZERO;
        let mut v2 = DVec3::ZERO;

        if index > 0 && index < self.my_laws.len() as i32 {
            if iloc <= 0 {
                self.my_laws[(index - 1) as usize]
                    .borrow()
                    .get_domain(&mut first, &mut last);
                self.my_laws[(index - 1) as usize].borrow().d0(last, &mut m1, &mut v1);
            }

            if iloc >= 0 {
                self.my_laws[index as usize]
                    .borrow()
                    .get_domain(&mut first, &mut last);
                if iloc == 0 {
                    self.my_laws[index as usize].borrow().d0(first, &mut m2, &mut v2);
                } else {
                    self.my_laws[index as usize].borrow().d0(first, &mut m1, &mut v1);
                }
            }
        }

        if index == 0 || index == self.my_laws.len() as i32 {
            if !self.path_is_closed(brep) || self.is_g1(brep, index, tol_min, 1.0e-4) != 1 {
                is_bary = false;
                if index == 0 {
                    self.my_laws[0].borrow().get_domain(&mut first, &mut last);
                    self.my_laws[0].borrow().d0(first, &mut m1, &mut v1);
                } else {
                    let n = self.my_laws.len() - 1;
                    self.my_laws[n].borrow().get_domain(&mut first, &mut last);
                    self.my_laws[n].borrow().d0(last, &mut m1, &mut v1);
                }
            } else {
                if iloc <= 0 {
                    let n = self.my_laws.len() - 1;
                    self.my_laws[n].borrow().get_domain(&mut first, &mut last);
                    self.my_laws[n].borrow().d0(last, &mut m1, &mut v1);
                }

                if iloc >= 0 {
                    self.my_laws[0].borrow().get_domain(&mut first, &mut last);
                    if iloc == 0 {
                        self.my_laws[0].borrow().d0(first, &mut m2, &mut v2);
                    } else {
                        self.my_laws[0].borrow().d0(first, &mut m1, &mut v1);
                    }
                }
            }
        }

        let mut p = brep.vertex(input.clone()).point; // BRep_Tool::Pnt(Input)

        if is_bary {
            // gp_XYZ P1(P.XYZ()), P2(P.XYZ()); P1 *= M1; P1 += V1.XYZ();
            // P2 *= M2; P2 += V2.XYZ();
            let mut p1 = GpMat::multiply_xyz_row(p, &m1) + v1;
            let p2 = GpMat::multiply_xyz_row(p, &m2) + v2;

            // P.ChangeCoord().SetLinearForm(0.5, P1, 0.5, P2);
            p = p1 * 0.5 + p2 * 0.5;
            p1 -= p2;
            let mut tol = p1.length() / 2.0;
            tol += tol_min;
            *output = b.add_vertex(brep, p, tol);
        } else {
            // P.ChangeCoord() *= M1; P.ChangeCoord() += V1.XYZ();
            p = GpMat::multiply_xyz_row(p, &m1) + v1;
            *output = b.add_vertex(brep, p, tol_min);
        }
    }

    /// OCCT CurvilinearBounds (cxx L456-476) — Return the Curvilinear Bounds
    /// of the <Index> Law.
    pub fn curvilinear_bounds(&mut self, index: i32, first: &mut f64, last: &mut f64) {
        *first = self.my_length[(index - 1) as usize];
        *last = self.my_length[index as usize];
        if *last < 0.0 {
            // It is required to carry out the calculation
            let nb_e = self.my_edges.len() as i32;
            let mut length = 0.0f64;
            let mut f = 0.0;
            let mut l = 0.0;
            for ii in 1i32..=nb_e {
                self.my_laws[(ii - 1) as usize]
                    .borrow()
                    .get_domain(&mut f, &mut l);
                // GCPnts_AbscissaPoint::Length(*myLaws->Value(ii)->GetCurve(),
                //                              myTol)
                let curve = self.my_laws[(ii - 1) as usize].borrow().get_curve();
                if let Some(curve) = curve {
                    length += arc_length(
                        &curve,
                        curve_first_parameter(&curve),
                        curve_last_parameter(&curve),
                    )
                    .abs();
                }
                self.my_length[ii as usize] = length;
            }

            *first = self.my_length[(index - 1) as usize];
            *last = self.my_length[index as usize];
        }
    }

    /// OCCT IsClosed (cxx L478-488).
    pub fn is_closed(&self, brep: &BRep) -> bool {
        if self.path_is_closed(brep) {
            return true;
        }

        // TopExp::Vertices(myPath, V1, V2); return (V1.IsSame(V2));
        let (v1, v2) = top_exp_wire_vertices(brep, &self.my_path);
        v1.ptr_id() == v2.ptr_id()
    }

    /// OCCT IsG1 (cxx L494-623) — Evaluate the continuity of the law by a
    /// vertex: -1 not connex, 0 G0, 1 tangent (G1).
    pub fn is_g1(
        &self,
        brep: &BRep,
        index: i32,
        spatial_tolerance: f64,
        angular_tolerance: f64,
    ) -> i32 {
        let mut v1 = DVec3::ZERO;
        let mut dv1 = DVec3::ZERO;
        let mut v2 = DVec3::ZERO;
        let mut dv2 = DVec3::ZERO;
        let mut m1 = GpMat::identity();
        let mut m2 = GpMat::identity();
        let mut dm1 = GpMat::identity();
        let mut dm2 = GpMat::identity();
        let mut first = 0.0;
        let mut last = 0.0;
        let eps_nul = 1.0e-12;
        let mut tol_eps = spatial_tolerance;
        let mut ok_d1 = false;
        let mut v = Shape::null();
        let mut e = Shape::null();
        // NCollection_Array1<gp_Pnt2d> Bid1(1, 1); / gp_Vec2d Bid2(1, 1);
        let mut bid1 = [DVec2::ZERO; 1];
        let mut bid2 = [DVec2::ZERO; 1];

        if index > 0 && index < self.my_laws.len() as i32 {
            self.my_laws[(index - 1) as usize]
                .borrow()
                .get_domain(&mut first, &mut last);
            ok_d1 = self.my_laws[(index - 1) as usize].borrow().d1(
                last,
                &mut m1,
                &mut v1,
                &mut dm1,
                &mut dv1,
                &mut bid1,
                &mut bid2,
            );
            if !ok_d1 {
                self.my_laws[(index - 1) as usize].borrow().d0(last, &mut m1, &mut v1);
            }

            self.my_laws[index as usize]
                .borrow()
                .get_domain(&mut first, &mut last);
            if ok_d1 {
                ok_d1 = self.my_laws[index as usize].borrow().d1(
                    first,
                    &mut m2,
                    &mut v2,
                    &mut dm2,
                    &mut dv2,
                    &mut bid1,
                    &mut bid2,
                );
            }
            if !ok_d1 {
                self.my_laws[index as usize].borrow().d0(first, &mut m2, &mut v2);
            }

            e = self.my_edges[index as usize].clone();
        }
        if index == 0 || index == self.my_laws.len() as i32 {
            if !self.path_is_closed(brep) {
                return -1;
            }
            let n = self.my_laws.len() - 1;
            self.my_laws[n].borrow().get_domain(&mut first, &mut last);
            ok_d1 = self.my_laws[n].borrow().d1(
                last,
                &mut m1,
                &mut v1,
                &mut dm1,
                &mut dv1,
                &mut bid1,
                &mut bid2,
            );
            if !ok_d1 {
                self.my_laws[n].borrow().d0(last, &mut m1, &mut v1);
            }

            self.my_laws[0].borrow().get_domain(&mut first, &mut last);
            if ok_d1 {
                self.my_laws[0].borrow().d1(
                    first,
                    &mut m2,
                    &mut v2,
                    &mut dm2,
                    &mut dv2,
                    &mut bid1,
                    &mut bid2,
                );
            }
            if !ok_d1 {
                self.my_laws[0].borrow().d0(first, &mut m2, &mut v2);
            }

            e = self.my_edges[0].clone();
        }

        if e.orientation == Orientation::Reversed {
            v = top_exp_last_vertex_stored(brep, &e);
        } else {
            v = top_exp_first_vertex_stored(brep, &e);
        }

        tol_eps += 2.0 * shape_tolerance(brep, &v);

        let mut is_g0 = true;
        let mut is_g1 = true;

        if (v1 - v2).length() > tol_eps {
            is_g0 = false;
        }
        // if (Norm(M1 - M2) > SpatialTolerance)
        if norm_mat(&gp_mat_minus(&m1, &m2)) > spatial_tolerance {
            is_g0 = false;
        }

        if !is_g0 {
            return -1;
        }
        if !ok_d1 {
            return 0; // No control of the derivative
        }

        // if ((DV1.Magnitude() > EpsNul) && (DV2.Magnitude() > EpsNul)
        //     && (DV1.Angle(DV2) > AngularTolerance))
        if dv1.length() > eps_nul && dv2.length() > eps_nul && gv_angle(dv1, dv2) > angular_tolerance
        {
            is_g1 = false;
        }

        // For the next, the tests are mostly empirical
        let norm1 = norm_mat(&dm1);
        let norm2 = norm_mat(&dm2);
        // It two 2 norms are null, it is good
        if norm1 > eps_nul || norm2 > eps_nul {
            // otherwise the normalized matrices are compared
            if norm1 > eps_nul && norm2 > eps_nul {
                dm1 = dm1.multiplied_scalar(1.0 / norm1); // DM1 /= Norm1
                dm2 = dm2.multiplied_scalar(1.0 / norm2);
                if norm_mat(&gp_mat_minus(&dm1, &dm2)) > angular_tolerance {
                    is_g1 = false;
                }
            } else {
                is_g1 = false; // 1 Null the other is not
            }
        }

        if is_g1 {
            1
        } else {
            0
        }
    }

    /// OCCT Parameter (cxx L627-680) — Find the index Law and the parameter,
    /// for a given Curvilinear abscissa.
    pub fn parameter(&mut self, abcissa: f64, index: &mut i32, u: &mut f64) {
        let nb_e = self.my_edges.len() as i32;
        let mut trouve = false;

        // Control that the lengths are calculated
        if self.my_length[nb_e as usize] < 0.0 {
            let mut f = 0.0;
            let mut l = 0.0;
            self.curvilinear_bounds(nb_e, &mut f, &mut l);
        }

        // Find the interval
        let mut iedge = 1i32;
        while iedge <= nb_e && !trouve {
            if self.my_length[iedge as usize] >= abcissa {
                trouve = true;
            } else {
                iedge += 1;
            }
        }

        if trouve {
            let mut f = 0.0;
            let mut l = 0.0;
            let law = self.my_laws[(iedge - 1) as usize].clone();
            law.borrow().get_domain(&mut f, &mut l);

            if abcissa == self.my_length[iedge as usize] {
                *u = l;
            } else if abcissa == self.my_length[(iedge - 1) as usize] {
                *u = f;
            } else {
                // GCPnts_AbscissaPoint(myTol, *GetCurve(),
                //                      Abcissa - myLength(iedge), f)
                let curve = law.borrow().get_curve();
                if let Some(curve) = curve {
                    *u = abscissa_point_parameter(
                        &curve,
                        f,
                        curve_last_parameter(&curve),
                        abcissa - self.my_length[(iedge - 1) as usize],
                        f,
                    );
                }
            }
            *index = iedge;
        } else {
            *index = 0;
        }
    }

    /// OCCT D0 (cxx L686-723) — Position of a section, with a given
    /// curvilinear abscissa.
    pub fn d0(&mut self, abcissa: f64, w: &mut Shape) {
        let mut u = 0.0;
        let mut ind = 0i32;
        let mut m = GpMat::identity();
        let mut v = DVec3::ZERO;

        self.parameter(abcissa, &mut ind, &mut u);
        if ind != 0 {
            // Positionement
            self.my_laws[(ind - 1) as usize].borrow().d0(u, &mut m, &mut v);
            // gp_Trsf fila;
            // fila.SetValues(M(1,1), M(1,2), M(1,3), V.X(), M(2,1), ... V.Z())
            match gp_trsf_set_values(&[
                m.mat[0][0], m.mat[0][1], m.mat[0][2], v.x, m.mat[1][0], m.mat[1][1], m.mat[1][2],
                v.y, m.mat[2][0], m.mat[2][1], m.mat[2][2], v.z,
            ]) {
                Ok(fila) => {
                    // W = BRepBuilderAPI_Transform(W, fila, true); // copy
                    *w = brep_builder_api_transform_copy(w, &fila, true);
                }
                Err(()) => panic!("gp_Trsf::SetValues, null determinant"),
            }
        } else {
            // W.Nullify();
            *w = Shape::null();
        }
    }

    /// OCCT Abscissa (cxx L729-744) — Calculate the abscissa of a point.
    pub fn abscissa(&mut self, index: i32, param: f64) -> f64 {
        let mut length = self.my_length[(index - 1) as usize];
        if length < 0.0 {
            let mut bid = 0.0;
            self.curvilinear_bounds(index, &mut bid, &mut length);
        }

        // GCPnts_AbscissaPoint::Length(*GetCurve(),
        //                              GetCurve()->FirstParameter(), Param,
        //                              myTol)
        let curve = self.my_laws[(index - 1) as usize].borrow().get_curve();
        if let Some(curve) = curve {
            length += arc_length(&curve, curve_first_parameter(&curve), param);
        }
        length
    }
}

/// OCCT BRepFill_LocationLaw is its own handle target (the base-class
/// upcast form).
impl BRepFillLocationLawOps for BRepFillLocationLaw {
    fn base(&self) -> &BRepFillLocationLaw {
        self
    }
    fn base_mut(&mut self) -> &mut BRepFillLocationLaw {
        self
    }
}
