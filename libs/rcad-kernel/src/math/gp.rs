//! OCCT gp package (TKMath) — axes, lines and rigid transformations used by
//! the helix algorithms.
//!
//! 1:1 subset of:
//! - `gp_Ax1`   (location + direction)
//! - `gp_Ax2`   (right-handed 2-axis frame; Y = N ^ Vx)
//! - `gp_Ax3`   (coordinate system with sense flag)
//! - `gp_Lin`   (infinite line)
//! - `gp_Trsf`  (rigid transformation, `SetDisplacement`)
//!
//! Only the members consumed by HelixGeom / HelixBRep are provided; the
//! operation semantics (normalization, cross products, matrix layout and
//! floating-point operation order) follow the OCCT sources exactly.
//! Matrix convention: `matrix[i][j]` == OCCT `gp_Trsf::Value(i+1, j+1)`.

use crate::core::precision::ANGULAR;
use glam::DVec3;

/// OCCT gp_Ax1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ax1 {
    pub location: DVec3,
    pub direction: DVec3,
}

impl Ax1 {
    /// OCCT gp_Ax1(P, V).
    pub fn new(location: DVec3, direction: DVec3) -> Self {
        Ax1 {
            location,
            direction: direction.normalize_or_zero(),
        }
    }

    /// OCCT gp_Ax1::SetLocation.
    pub fn set_location(&mut self, p: DVec3) {
        self.location = p;
    }

    /// OCCT gp_Ax1::SetDirection.
    pub fn set_direction(&mut self, v: DVec3) {
        self.direction = v.normalize_or_zero();
    }
}

/// OCCT gp_Ax2 (always right-handed: YDirection = N ^ Vx).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ax2 {
    pub location: DVec3,
    pub direction: DVec3,
    pub x_direction: DVec3,
    pub y_direction: DVec3,
}

impl Ax2 {
    /// OCCT gp_Ax2(P, N, Vx): main direction N; the X direction is
    /// re-orthogonalized against N (`Vx -> N ^ (Vx ^ N)`); Y = N ^ X.
    pub fn new(location: DVec3, direction: DVec3, x_direction: DVec3) -> Self {
        let n = direction.normalize_or_zero();
        let vx = x_direction.normalize_or_zero();
        // OCCT gp_Ax2.cxx: XDirection = N ^ (Vx ^ N), normalized.
        let xd = n.cross(vx).cross(n);
        let xd = if xd.length_squared() < f64::EPSILON {
            // Degenerate (Vx parallel to N) — OCCT raises ConstructionError;
            // never hit with valid input.
            DVec3::ZERO
        } else {
            xd.normalize_or_zero()
        };
        let yd = n.cross(xd).normalize_or_zero();
        Ax2 {
            location,
            direction: n,
            x_direction: xd,
            y_direction: yd,
        }
    }

    /// OCCT gp_Ax2::Axis.
    pub fn axis(&self) -> Ax1 {
        Ax1::new(self.location, self.direction)
    }

    /// OCCT gp_Ax2(P, V) two-argument constructor (gp_Ax2.cxx L27-85): the X
    /// direction is the unit vector perpendicular to V having a zero in the
    /// coordinate of the smallest |component| of V, applied through
    /// SetXDirection (gp_Ax2.hxx L143-147).  With D already perpendicular to
    /// V, SetXDirection(D) equals Ax2::new(P, V, D).
    pub fn from_direction(location: DVec3, direction: DVec3) -> Self {
        let a = direction.x;
        let b = direction.y;
        let c = direction.z;
        let aabs = a.abs();
        let babs = b.abs();
        let cabs = c.abs();
        let d = if babs <= aabs && babs <= cabs {
            if aabs > cabs {
                DVec3::new(-c, 0.0, a)
            } else {
                DVec3::new(c, 0.0, -a)
            }
        } else if aabs <= babs && aabs <= cabs {
            if babs > cabs {
                DVec3::new(0.0, -c, b)
            } else {
                DVec3::new(0.0, c, -b)
            }
        } else if aabs > babs {
            DVec3::new(-b, a, 0.0)
        } else {
            DVec3::new(b, -a, 0.0)
        };
        Self::new(location, direction, d)
    }
}

/// OCCT gp_Ax3 (right- or left-handed coordinate system).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ax3 {
    /// OCCT gp_Ax3::Axis (main direction + location).
    pub axis: Ax1,
    pub y_direction: DVec3,
    pub x_direction: DVec3,
}

impl Ax3 {
    /// OCCT gp_Ax3() — default right-handed OXYZ.
    pub fn new() -> Self {
        Ax3 {
            axis: Ax1::new(DVec3::ZERO, DVec3::Z),
            y_direction: DVec3::Y,
            x_direction: DVec3::X,
        }
    }

    /// OCCT gp_Ax3(P, N, Vx) — right-handed system (via gp_Ax2).
    pub fn from_pnt_n_vx(location: DVec3, direction: DVec3, x_direction: DVec3) -> Self {
        let a2 = Ax2::new(location, direction, x_direction);
        Ax3 {
            axis: a2.axis(),
            y_direction: a2.y_direction,
            x_direction: a2.x_direction,
        }
    }

    /// OCCT gp_Ax3(const gp_Ax2& theA) — right-handed.
    pub fn from_ax2(a: &Ax2) -> Self {
        Ax3 {
            axis: a.axis(),
            y_direction: a.y_direction,
            x_direction: a.x_direction,
        }
    }

    /// OCCT gp_Ax3::Location.
    pub fn location(&self) -> DVec3 {
        self.axis.location
    }

    /// OCCT gp_Ax3::Direction.
    pub fn direction(&self) -> DVec3 {
        self.axis.direction
    }

    /// OCCT gp_Ax3::Direct (gp_Ax3.hxx L208): DERIVED from the frame, not a
    /// stored flag — `(vxdir ^ vydir) . axis.Direction() > 0`.  rcad carried a
    /// `sense` field that every constructor set to `true`, which cannot agree
    /// with OCCT once a setter leaves the frame indirect.
    pub fn direct(&self) -> bool {
        self.x_direction
            .cross(self.y_direction)
            .dot(self.axis.direction)
            > 0.0
    }

    /// OCCT gp_Ax3::SetLocation.
    pub fn set_location(&mut self, p: DVec3) {
        self.axis.location = p;
    }

    /// OCCT gp_Ax3::SetDirection (gp_Ax3.hxx L506-534) — 1:1, including the
    /// parallel branch that the previous version was missing (it produced a
    /// null X direction instead of swapping the axes).
    pub fn set_direction(&mut self, v: DVec3) {
        let v = v.normalize_or_zero();
        // OCCT L508: double aDot = theV.Dot(vxdir).
        let a_dot = v.dot(self.x_direction);
        if 1.0 - a_dot.abs() <= ANGULAR {
            // OCCT L509-521: theV is (anti)parallel to the current X direction,
            // so `theV ^ (vxdir ^ theV)` would degenerate.  OCCT relabels the
            // axes instead of producing a null vector.
            if a_dot > 0.0 {
                // OCCT L511-515.
                self.x_direction = self.y_direction;
                self.y_direction = self.axis.direction;
            } else {
                // OCCT L516-519.
                self.x_direction = self.axis.direction;
            }
            // OCCT L520: axis.SetDirection(theV).
            self.axis.direction = v;
        } else {
            // OCCT L523-533.  `direct` is read BEFORE the main direction moves.
            let direct = self.direct();
            self.axis.direction = v;
            // vxdir = theV.CrossCrossed(vxdir, theV): `gp_XYZ::CrossCrossed`
            // computes `<me> ^ (theCoord1 ^ theCoord2)` (gp_XYZ.hxx L238-246).
            self.x_direction = v.cross(self.x_direction.cross(v)).normalize_or_zero();
            self.y_direction = if direct {
                v.cross(self.x_direction)
            } else {
                self.x_direction.cross(v)
            }
            .normalize_or_zero();
        }
    }

    /// OCCT gp_Ax3::SetXDirection (gp_Ax3.hxx L540-570) — 1:1, including the
    /// parallel branch.
    pub fn set_x_direction(&mut self, vx: DVec3) {
        let vx = vx.normalize_or_zero();
        let n = self.axis.direction;
        // OCCT L542: double aDot = theVx.Dot(axis.Direction()).
        let a_dot = vx.dot(n);
        if 1.0 - a_dot.abs() <= ANGULAR {
            // OCCT L543-555.
            if a_dot > 0.0 {
                // OCCT L545-549: axis.SetDirection(vxdir); vydir = -vydir.
                self.axis.direction = self.x_direction;
                self.y_direction = -self.y_direction;
            } else {
                // OCCT L550-553: axis.SetDirection(vxdir).
                self.axis.direction = self.x_direction;
            }
            // OCCT L554: vxdir = theVx.
            self.x_direction = vx;
        } else {
            // OCCT L556-568.
            let direct = self.direct();
            // vxdir = axis.Direction().CrossCrossed(theVx, axis.Direction()).
            self.x_direction = n.cross(vx.cross(n)).normalize_or_zero();
            self.y_direction = if direct {
                n.cross(self.x_direction)
            } else {
                self.x_direction.cross(n)
            }
            .normalize_or_zero();
        }
    }
}

impl Default for Ax3 {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT gp_Lin.
#[derive(Debug, Clone, Copy)]
pub struct Lin {
    pub pos: DVec3,
    pub dir: DVec3,
}

impl Lin {
    /// OCCT gp_Lin(theA1) — line through the axis location along its direction.
    pub fn from_ax1(a: &Ax1) -> Self {
        Lin {
            pos: a.location,
            dir: a.direction,
        }
    }

    /// OCCT gp_Lin::Distance(theP) — perpendicular distance to the line.
    pub fn distance(&self, p: DVec3) -> f64 {
        let d = p - self.pos;
        d.cross(self.dir).length()
    }

    /// OCCT gp_Lin(gp_Pnt, gp_Dir) — the direction is normalized.
    pub fn from_pnt_dir(p: DVec3, dir: DVec3) -> Self {
        Lin {
            pos: p,
            dir: dir.normalize_or_zero(),
        }
    }

    /// OCCT gp_Lin::Transform(theT) — gp_Lin.hxx L178: `pos.Transform(theT)`
    /// where `pos` is a gp_Ax1 (gp_Lin.hxx L213); gp_Ax1::Transform
    /// (gp_Ax1.hxx L201-205) moves the location AND transforms the
    /// direction (a gp_Dir: normalized).
    pub fn transform(&mut self, t: &Trsf) {
        self.pos = t.apply(self.pos);
        self.dir = t.transform_dir(self.dir);
    }
}

/// OCCT gp_TrsfForm (gp_TrsfForm.hxx) — the shape tag driving the
/// shape-special-cased gp_Trsf paths (SetScaleFactor / SetTranslationPart /
/// Invert / Multiply).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrsfForm {
    Identity,
    Rotation,
    Translation,
    PntMirror,
    Ax1Mirror,
    Ax2Mirror,
    Scale,
    CompoundTrsf,
    Other,
}

/// OCCT gp_Trsf — transformation with a separate scale factor and form tag
/// (gp_Trsf.hxx L379-385 member order: scale, shape, matrix, loc).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trsf {
    /// Row-major 3x3: `matrix[i][j]` == OCCT `matrix.Value(i+1, j+1)`
    /// (WITHOUT the scale folded in — OCCT keeps `scale` separate).
    pub matrix: [[f64; 3]; 3],
    /// Translation part.
    pub loc: DVec3,
    /// OCCT gp_Trsf scale factor (gp_Trsf.hxx L227).
    pub scale: f64,
    /// OCCT gp_TrsfForm tag (gp_Trsf.hxx member `shape`).
    pub form: TrsfForm,
}

/// OCCT gp::Resolution() — gp.hxx L60 = RealSmall() = DBL_MIN
/// (Standard_Real.hxx L132-134).
pub const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT gp_Mat::SetRotation(theAxis, theAng) — gp_Mat.cxx L122-159.
/// Rodrigues' rotation formula: R = I + sin(theta)K + (1-cos(theta))K^2,
/// where K is the skew-symmetric matrix of the normalized axis.
fn mat_set_rotation(matrix: &mut [[f64; 3]; 3], axis: DVec3, ang: f64) {
    let a_v = axis.normalize_or_zero(); // OCCT L126: theAxis.Normalized().

    let a = a_v.x; // OCCT L128.
    let b = a_v.y; // OCCT L129.
    let c = a_v.z; // OCCT L130.

    let a_cos = ang.cos(); // OCCT L133.
    let a_sin = ang.sin(); // OCCT L134.
    let a_om_cos = 1.0 - a_cos; // OCCT L135: one minus cosine.

    let a2 = a * a; // OCCT L138-143: precomputed terms.
    let b2 = b * b;
    let c2 = c * c;
    let ab = a * b;
    let ac = a * c;
    let bc = b * c;

    // OCCT L148-158: R = I + sin(theta)K + (1-cos(theta))K^2.
    matrix[0][0] = 1.0 + a_om_cos * (-(b2 + c2));
    matrix[0][1] = a_om_cos * ab - a_sin * c;
    matrix[0][2] = a_om_cos * ac + a_sin * b;

    matrix[1][0] = a_om_cos * ab + a_sin * c;
    matrix[1][1] = 1.0 + a_om_cos * (-(a2 + c2));
    matrix[1][2] = a_om_cos * bc - a_sin * a;

    matrix[2][0] = a_om_cos * ac - a_sin * b;
    matrix[2][1] = a_om_cos * bc + a_sin * a;
    matrix[2][2] = 1.0 + a_om_cos * (-(a2 + b2));
}

impl Trsf {
    /// OCCT gp_Trsf default constructor (gp_Trsf.hxx L379-385):
    /// scale = 1, shape = gp_Identity, identity matrix, zero location.
    pub fn identity() -> Self {
        Trsf {
            matrix: [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            loc: DVec3::ZERO,
            scale: 1.0,
            form: TrsfForm::Identity,
        }
    }

    /// OCCT gp_Trsf::SetRotation(A1, Ang) — gp_Trsf.cxx L90-101.
    pub fn set_rotation(&mut self, ax1: &Ax1, ang: f64) {
        self.form = TrsfForm::Rotation; // OCCT L92: shape = gp_Rotation.
        self.scale = 1.0; // OCCT L93: scale = 1.
        self.loc = ax1.location; // OCCT L94: loc = A1.Location().XYZ().
        // OCCT L95: matrix.SetRotation(A1.Direction().XYZ(), Ang).
        mat_set_rotation(&mut self.matrix, ax1.direction, ang);
        self.loc = -self.loc; // OCCT L96: loc.Reverse().
        // OCCT L97: loc.Multiply(matrix).  gp_XYZ::Multiply(const gp_Mat&) is
        // documented as "<me> = theMatrix * <me>" (gp_XYZ.hxx L308-309), i.e.
        // the matrix acts on the left.  The previous row-vector form computed
        // the transpose (M^T * loc) and therefore produced a wrong offset for
        // every axis that does not pass through the origin.
        let l = self.loc;
        self.loc = DVec3::new(
            self.matrix[0][0] * l.x + self.matrix[0][1] * l.y + self.matrix[0][2] * l.z,
            self.matrix[1][0] * l.x + self.matrix[1][1] * l.y + self.matrix[1][2] * l.z,
            self.matrix[2][0] * l.x + self.matrix[2][1] * l.y + self.matrix[2][2] * l.z,
        );
        self.loc += ax1.location; // OCCT L98: loc.Add(A1.Location().XYZ()).
    }

    /// OCCT gp_Trsf::SetDisplacement(FromA1, ToA2) — gp_Trsf.cxx L218-240,
    /// with the same floating-point operation order.
    pub fn set_displacement(from_a1: &Ax3, to_a2: &Ax3) -> Self {
        // OCCT L219-220: shape = gp_CompoundTrsf; scale = 1.0.
        // matrix from ToA2 to XOY: SetCol(1..3, ToA2 X / Y / Z directions).
        let tx = to_a2.x_direction;
        let ty = to_a2.y_direction;
        let tz = to_a2.direction();
        let mut matrix = [[0.0f64; 3]; 3];
        // SetCol(j, xyz): column j gets the vector.
        matrix[0][0] = tx.x;
        matrix[1][0] = tx.y;
        matrix[2][0] = tx.z;
        matrix[0][1] = ty.x;
        matrix[1][1] = ty.y;
        matrix[2][1] = ty.z;
        matrix[0][2] = tz.x;
        matrix[1][2] = tz.y;
        matrix[2][2] = tz.z;
        let mut loc = to_a2.location();

        // matrix XOY to FromA1: gp_Mat MA1(xDir, yDir, zDir) sets the COLUMNS,
        // then MA1.Transpose() — so MA1 rows are the FromA1 directions.
        let x_dir = from_a1.x_direction;
        let y_dir = from_a1.y_direction;
        let z_dir = from_a1.direction();
        let ma1 = [
            [x_dir.x, x_dir.y, x_dir.z],
            [y_dir.x, y_dir.y, y_dir.z],
            [z_dir.x, z_dir.y, z_dir.z],
        ];
        // MA1loc = FromA1.Location(); MA1loc.Multiply(MA1).  gp_XYZ::Multiply
        // (const gp_Mat&) is "<me> = theMatrix * <me>" (gp_XYZ.hxx L308-309):
        // the matrix acts on the LEFT.  The transposed form used here before
        // made every FromA1 with a non-zero location wrong.
        let loc_from = from_a1.location();
        let mut ma1_loc = [0.0f64; 3];
        for i in 0..3 {
            ma1_loc[i] =
                ma1[i][0] * loc_from.x + ma1[i][1] * loc_from.y + ma1[i][2] * loc_from.z;
        }
        // MA1loc.Reverse().
        for v in ma1_loc.iter_mut() {
            *v = -*v;
        }
        // MA1loc.Multiply(matrix) — again the matrix on the LEFT.
        let mut mult = [0.0f64; 3];
        for i in 0..3 {
            mult[i] =
                matrix[i][0] * ma1_loc[0] + matrix[i][1] * ma1_loc[1] + matrix[i][2] * ma1_loc[2];
        }
        // loc.Add(MA1loc).
        loc += DVec3::from_array(mult);
        // matrix.Multiply(MA1) — full 3x3 product.
        let mut prod = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0f64;
                for k in 0..3 {
                    s += matrix[i][k] * ma1[k][j];
                }
                prod[i][j] = s;
            }
        }
        Trsf {
            matrix: prod,
            loc,
            scale: 1.0,
            form: TrsfForm::CompoundTrsf,
        }
    }

    /// OCCT gp_Trsf::SetTransformation(A3) — gp_Trsf.cxx L196-206.
    /// The transformation maps points of the reference frame to the A3 frame
    /// (rows of the matrix are the A3 directions).
    pub fn set_transformation_ax3(a3: &Ax3) -> Self {
        // matrix.SetRows(A3.XDirection(), A3.YDirection(), A3.Direction()).
        let x_dir = a3.x_direction;
        let y_dir = a3.y_direction;
        let z_dir = a3.direction();
        let matrix = [
            [x_dir.x, x_dir.y, x_dir.z],
            [y_dir.x, y_dir.y, y_dir.z],
            [z_dir.x, z_dir.y, z_dir.z],
        ];
        // loc = A3.Location(); loc.Multiply(matrix) — matrix times the
        // location (column-vector convention); loc.Reverse().
        let p = a3.location();
        let mut loc = DVec3::new(
            matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2] * p.z,
            matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2] * p.z,
            matrix[2][0] * p.x + matrix[2][1] * p.y + matrix[2][2] * p.z,
        );
        loc = -loc;
        Trsf {
            matrix,
            loc,
            scale: 1.0,
            form: TrsfForm::CompoundTrsf,
        }
    }

    /// OCCT gp_Trsf::SetTransformation(const gp_Ax3& FromA1, const gp_Ax3&
    /// ToA2) — gp_Trsf.cxx L172-194.  This is the two-argument overload, a
    /// DIFFERENT function from `set_displacement`: it maps the COORDINATES of
    /// a point expressed in FromA1 into the equivalent coordinates expressed
    /// in ToA2, whereas SetDisplacement carries the FromA1 frame ONTO the ToA2
    /// frame.  OCCT states the relation directly (gp_Trsf.hxx L136-137): "the
    /// vectorial part of one is the inverse of the vectorial part of the
    /// other".
    pub fn set_transformation_from_to(from_a1: &Ax3, to_a2: &Ax3) -> Self {
        // OCCT L175-179: matrix.SetRows(ToA2 XDirection / YDirection /
        // Direction) — SetRows puts the vectors in the ROWS.
        let tx = to_a2.x_direction;
        let ty = to_a2.y_direction;
        let tz = to_a2.direction();
        let matrix = [
            [tx.x, tx.y, tx.z],
            [ty.x, ty.y, ty.z],
            [tz.x, tz.y, tz.z],
        ];
        // OCCT L177-179: loc = ToA2.Location(); loc.Multiply(matrix);
        // loc.Reverse().  gp_XYZ::Multiply(const gp_Mat&) is
        // "<me> = theMatrix * <me>" (gp_XYZ.hxx L308-309) — matrix on the LEFT.
        let l2 = to_a2.location();
        let mut loc = -DVec3::new(
            matrix[0][0] * l2.x + matrix[0][1] * l2.y + matrix[0][2] * l2.z,
            matrix[1][0] * l2.x + matrix[1][1] * l2.y + matrix[1][2] * l2.z,
            matrix[2][0] * l2.x + matrix[2][1] * l2.y + matrix[2][2] * l2.z,
        );

        // OCCT L181-185: gp_Mat MA1(xDir, yDir, zDir) sets the COLUMNS —
        // note there is NO Transpose() here, unlike SetDisplacement.
        let x_dir = from_a1.x_direction;
        let y_dir = from_a1.y_direction;
        let z_dir = from_a1.direction();
        let ma1 = [
            [x_dir.x, y_dir.x, z_dir.x],
            [x_dir.y, y_dir.y, z_dir.y],
            [x_dir.z, y_dir.z, z_dir.z],
        ];
        // OCCT L186-189: MA1loc = FromA1.Location(); MA1loc.Multiply(matrix) —
        // this multiplies by `matrix`, NOT by MA1; then loc.Add(MA1loc).
        let l1 = from_a1.location();
        loc += DVec3::new(
            matrix[0][0] * l1.x + matrix[0][1] * l1.y + matrix[0][2] * l1.z,
            matrix[1][0] * l1.x + matrix[1][1] * l1.y + matrix[1][2] * l1.z,
            matrix[2][0] * l1.x + matrix[2][1] * l1.y + matrix[2][2] * l1.z,
        );
        // OCCT L190: matrix.Multiply(MA1) — gp_Mat::Multiply(B) is `this = this
        // * B`, i.e. FromA1 to XOY composed after XOY to ToA2.
        let mut prod = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0f64;
                for k in 0..3 {
                    s += matrix[i][k] * ma1[k][j];
                }
                prod[i][j] = s;
            }
        }
        Trsf {
            matrix: prod,
            loc,
            scale: 1.0,
            form: TrsfForm::CompoundTrsf,
        }
    }

    /// OCCT gp_Trsf::SetTranslation(V) — gp_Trsf.hxx L400-406.
    pub fn set_translation(&mut self, v: DVec3) {
        self.form = TrsfForm::Translation; // OCCT L402: shape = gp_Translation.
        self.scale = 1.0; // OCCT L403: scale = 1.
        self.matrix = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]; // OCCT L404: matrix.SetIdentity().
        self.loc = v; // OCCT L405: loc = theV.XYZ().
    }

    /// OCCT gp_Trsf::SetTranslationPart(V) — gp_Trsf.cxx L244-281.
    pub fn set_translation_part(&mut self, v: DVec3) {
        self.loc = v;
        let loc_null = self.loc.length_squared() < GP_RESOLUTION;
        match self.form {
            TrsfForm::Identity => {
                if !loc_null {
                    self.form = TrsfForm::Translation;
                }
            }
            TrsfForm::Translation => {
                if loc_null {
                    self.form = TrsfForm::Identity;
                }
            }
            // gp_Rotation, gp_PntMirror, gp_Ax1Mirror, gp_Ax2Mirror,
            // gp_Scale, gp_CompoundTrsf, gp_Other (OCCT L264-280).
            _ => {
                if !loc_null {
                    self.form = TrsfForm::CompoundTrsf;
                }
            }
        }
    }

    /// OCCT gp_Trsf::IsNegative() — gp_Trsf.hxx L218: scale < 0.0.
    pub fn is_negative(&self) -> bool {
        self.scale < 0.0
    }

    /// OCCT gp_Trsf::SetScaleFactor(S) — gp_Trsf.cxx L284-330.
    pub fn set_scale_factor(&mut self, s: f64) {
        if s.abs() <= GP_RESOLUTION {
            panic!("gp_Trsf::SetScaleFactor");
        }
        self.scale = s;
        let unit = (self.scale - 1.0).abs() <= GP_RESOLUTION; // = (scale == 1)
        let munit = (self.scale + 1.0).abs() <= GP_RESOLUTION; // = (scale == -1)
        match self.form {
            TrsfForm::Identity | TrsfForm::Translation => {
                if !unit {
                    self.form = TrsfForm::Scale;
                }
                if munit {
                    self.form = TrsfForm::PntMirror;
                }
            }
            TrsfForm::Rotation => {
                if !unit {
                    self.form = TrsfForm::CompoundTrsf;
                }
            }
            TrsfForm::PntMirror | TrsfForm::Ax1Mirror | TrsfForm::Ax2Mirror => {
                if !munit {
                    self.form = TrsfForm::Scale;
                }
                if unit {
                    self.form = TrsfForm::Identity;
                }
            }
            TrsfForm::Scale => {
                if unit {
                    self.form = TrsfForm::Identity;
                }
                if munit {
                    self.form = TrsfForm::PntMirror;
                }
            }
            TrsfForm::CompoundTrsf | TrsfForm::Other => {}
        }
    }

    /// OCCT gp_Trsf::Invert() — gp_Trsf.cxx L396-428.
    pub fn invert(&mut self) {
        //  X' = scale * R * X + T  =>  X = (R^t / scale) * (X' - T)
        // The scale is kept outside the matrix, so R^-1 = R transposed.
        match self.form {
            TrsfForm::Identity => {}
            TrsfForm::Translation | TrsfForm::PntMirror => {
                self.loc = -self.loc;
            }
            TrsfForm::Scale => {
                if self.scale.abs() <= GP_RESOLUTION {
                    panic!("gp_Trsf::Invert() - transformation has zero scale");
                }
                self.scale = 1.0 / self.scale;
                self.loc *= -self.scale;
            }
            _ => {
                if self.scale.abs() <= GP_RESOLUTION {
                    panic!("gp_Trsf::Invert() - transformation has zero scale");
                }
                self.scale = 1.0 / self.scale;
                // matrix.Transpose().
                let m = self.matrix;
                self.matrix = [
                    [m[0][0], m[1][0], m[2][0]],
                    [m[0][1], m[1][1], m[2][1]],
                    [m[0][2], m[1][2], m[2][2]],
                ];
                // loc.Multiply(matrix) — matrix times the location.
                let p = self.loc;
                self.loc = DVec3::new(
                    self.matrix[0][0] * p.x + self.matrix[0][1] * p.y + self.matrix[0][2] * p.z,
                    self.matrix[1][0] * p.x + self.matrix[1][1] * p.y + self.matrix[1][2] * p.z,
                    self.matrix[2][0] * p.x + self.matrix[2][1] * p.y + self.matrix[2][2] * p.z,
                );
                self.loc *= -self.scale;
            }
        }
    }

    /// OCCT gp_Trsf::VectorialPart() — gp_Trsf.hxx L464-471: the matrix with
    /// the scale folded in.
    pub fn vectorial_part(&self) -> [[f64; 3]; 3] {
        if self.scale == 1.0 {
            return self.matrix;
        }
        let mut m = self.matrix;
        for row in m.iter_mut() {
            for v in row.iter_mut() {
                *v *= self.scale;
            }
        }
        m
    }

    /// OCCT gp_Trsf::Value(Row, Col) — gp_Trsf.hxx L420-434: Col < 4 returns
    /// scale * matrix[Row-1][Col-1], Col 4 returns the location.
    pub fn value(&self, r: usize, c: usize) -> f64 {
        if c < 4 {
            self.scale * self.matrix[r - 1][c - 1]
        } else {
            [self.loc.x, self.loc.y, self.loc.z][r - 1]
        }
    }

    /// OCCT gp_Trsf::Multiplied(theT) — `this * theT`.
    ///
    /// Note: OCCT special-cases the form combinations in Multiply; this
    /// implementation is the compound general branch. All current consumers
    /// combine CompoundTrsf forms with scale 1 (SetDisplacement /
    /// SetTransformation results), where the branches agree.
    pub fn multiplied(&self, t: &Trsf) -> Trsf {
        let mut matrix = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0f64;
                for k in 0..3 {
                    s += self.matrix[i][k] * t.matrix[k][j];
                }
                matrix[i][j] = s;
            }
        }
        // Combined transform — OCCT gp_Trsf::Multiply's general branch
        // (gp_Trsf.cxx L544-559): Tloc = matrix * T.loc; if (scale != 1.0)
        // Tloc *= scale; loc += Tloc.  The scale factor on the translated
        // location was missing here before.
        let mut loc = self.loc;
        for i in 0..3 {
            let mut tloc = self.matrix[i][0] * t.loc.x
                + self.matrix[i][1] * t.loc.y
                + self.matrix[i][2] * t.loc.z;
            if self.scale != 1.0 {
                tloc *= self.scale;
            }
            loc[i] += tloc;
        }
        Trsf {
            matrix,
            loc,
            scale: self.scale * t.scale,
            form: TrsfForm::CompoundTrsf,
        }
    }

    /// OCCT gp_Trsf::Transforms(theCoord) — gp_Trsf.hxx L452-461:
    /// matrix multiply, then scale (only if != 1), then translate.
    pub fn apply(&self, p: DVec3) -> DVec3 {
        let mut x = self.matrix[0][0] * p.x + self.matrix[0][1] * p.y + self.matrix[0][2] * p.z;
        let mut y = self.matrix[1][0] * p.x + self.matrix[1][1] * p.y + self.matrix[1][2] * p.z;
        let mut z = self.matrix[2][0] * p.x + self.matrix[2][1] * p.y + self.matrix[2][2] * p.z;
        if self.scale != 1.0 {
            x *= self.scale;
            y *= self.scale;
            z *= self.scale;
        }
        DVec3::new(x + self.loc.x, y + self.loc.y, z + self.loc.z)
    }

    /// OCCT gp_Vec::Transform(theT) — gp_Vec.cxx L120-138: vectors feel no
    /// translation; the form decides the operation.
    pub fn transform_vec(&self, v: DVec3) -> DVec3 {
        match self.form {
            TrsfForm::Identity | TrsfForm::Translation => v,
            TrsfForm::PntMirror => -v,
            TrsfForm::Scale => v * self.scale,
            _ => {
                // coord.Multiply(theTransformation.VectorialPart()).
                let m = self.vectorial_part();
                DVec3::new(
                    m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                    m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                    m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
                )
            }
        }
    }

    /// OCCT gp_Dir::Transform(T) — gp_Dir.cxx L129-154: like the vector
    /// transform, followed by normalization (and a reversal for a negative
    /// scale).
    pub fn transform_dir(&self, v: DVec3) -> DVec3 {
        match self.form {
            TrsfForm::Identity | TrsfForm::Translation => v,
            TrsfForm::PntMirror => -v,
            TrsfForm::Scale => {
                if self.scale < 0.0 {
                    -v
                } else {
                    v
                }
            }
            _ => {
                // coord.Multiply(T.HVectorialPart()); coord.Divide(Modulus).
                let m = self.matrix;
                let mut d = DVec3::new(
                    m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                    m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                    m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
                );
                let len = d.length();
                d /= len;
                if self.scale < 0.0 {
                    d = -d;
                }
                d
            }
        }
    }

    /// The transformation as a glam affine (for `transform_curve`).
    pub fn to_daffine3(&self) -> glam::DAffine3 {
        let m = glam::DMat3::from_cols(
            DVec3::new(
                self.matrix[0][0],
                self.matrix[1][0],
                self.matrix[2][0],
            ),
            DVec3::new(
                self.matrix[0][1],
                self.matrix[1][1],
                self.matrix[2][1],
            ),
            DVec3::new(
                self.matrix[0][2],
                self.matrix[1][2],
                self.matrix[2][2],
            ),
        ) * self.scale;
        glam::DAffine3::from_mat3_translation(m, self.loc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT gp_Trsf::SetRotation(Z axis, PI/2) applied to (10,0,0) is
    /// (0,10,0): the axis location stays fixed and the rotation matrix is
    /// the Rodrigues form.
    #[test]
    fn set_rotation_about_z_axis() {
        let mut t = Trsf::identity();
        let axis = Ax1::new(DVec3::ZERO, DVec3::Z);
        t.set_rotation(&axis, std::f64::consts::FRAC_PI_2);
        assert_eq!(t.form, TrsfForm::Rotation);
        assert_eq!(t.scale, 1.0);
        let p = t.apply(DVec3::new(10.0, 0.0, 0.0));
        let expected = DVec3::new(0.0, 10.0, 0.0);
        assert!(
            (p - expected).length() < 1e-12,
            "rotation moved the point to {p:?}"
        );
    }

    /// OCCT gp_Trsf::SetRotation(A1, Ang) with a non-origin axis and an angle
    /// whose matrix is *not* symmetric: the result must be the Rodrigues
    /// rotation about the axis line.  The `PI` case above cannot detect a
    /// transposed `loc.Multiply(matrix)` (gp_XYZ::Multiply(const gp_Mat&) is
    /// `<me> = theMatrix * <me>`, gp_XYZ.hxx L308-309), because `M(PI)` is
    /// symmetric and `M^T L == M L` there.
    #[test]
    fn set_rotation_about_off_origin_axis_matches_rodrigues() {
        let loc = DVec3::new(0.0, 0.0, 50.0);
        let dir = DVec3::new(0.0, 1.0, 0.0);
        let axis = Ax1::new(loc, dir);
        let ang = 0.4_f64;
        let mut t = Trsf::identity();
        t.set_rotation(&axis, ang);
        for &p in &[
            DVec3::new(50.0, 20.3, 40.0),
            DVec3::new(-13.0, 7.5, 3.0),
            DVec3::new(30.0, -2.0, 99.0),
        ] {
            let d = p - loc;
            let d_par = dir * d.dot(dir);
            let d_perp = d - d_par;
            let expected = loc + d_par + d_perp * ang.cos() + dir.cross(d_perp) * ang.sin();
            let got = t.apply(p);
            assert!(
                (got - expected).length() < 1e-12,
                "p={p:?} got={got:?} expected={expected:?}"
            );
        }
    }

    /// OCCT gp_Trsf::SetRotation(A1, Ang) with a non-origin axis: the axis
    /// location is invariant (gp_Trsf.cxx L94-98 location handling).
    #[test]
    fn set_rotation_keeps_axis_location_invariant() {
        let mut t = Trsf::identity();
        let axis = Ax1::new(DVec3::new(1.0, 1.0, 0.0), DVec3::Z);
        t.set_rotation(&axis, std::f64::consts::PI);
        let p = t.apply(DVec3::new(1.0, 1.0, 5.0));
        let expected = DVec3::new(1.0, 1.0, 5.0);
        assert!(
            (p - expected).length() < 1e-12,
            "a point on the axis must stay fixed, got {p:?}"
        );
    }

    /// OCCT gp_Trsf::SetTranslation(V) — identity matrix, translation form.
    #[test]
    fn set_translation_moves_points() {
        let mut t = Trsf::identity();
        t.set_translation(DVec3::new(3.0, -4.0, 5.0));
        assert_eq!(t.form, TrsfForm::Translation);
        assert_eq!(t.scale, 1.0);
        assert_eq!(t.matrix, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let p = t.apply(DVec3::new(1.0, 1.0, 1.0));
        assert_eq!(p, DVec3::new(4.0, -3.0, 6.0));
    }

    /// OCCT gp_Trsf::SetDisplacement(FromA1, ToA2) — gp_Trsf.cxx L218-240:
    /// the result carries the FromA1 frame onto the ToA2 frame, i.e.
    /// `F * FromA1.Location() == ToA2.Location()` and `F * FromA1.XDir()` is
    /// `ToA2.XDir()` (and likewise for Y and Z).  Both `MA1loc.Multiply(...)`
    /// calls put the matrix on the LEFT (`<me> = theMatrix * <me>`,
    /// gp_XYZ.hxx L308-309); the transposed form is invisible when FromA1 is
    /// the default frame (identity axes, zero location) — which is the only
    /// call site in the repo — so this test moves AND tilts FromA1.
    #[test]
    fn set_displacement_carries_from_frame_onto_to_frame() {
        let ang = 0.7_f64;
        let from = Ax3::from_pnt_n_vx(
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(0.0, ang.cos(), ang.sin()),
            DVec3::X,
        );
        let to = Ax3::from_pnt_n_vx(DVec3::new(5.0, 7.0, 11.0), DVec3::Z, DVec3::Y);
        let t = Trsf::set_displacement(&from, &to);
        assert_eq!(t.form, TrsfForm::CompoundTrsf);
        assert_eq!(t.scale, 1.0);
        let got = t.apply(from.location());
        let expected = to.location();
        assert!(
            (got - expected).length() < 1e-12,
            "FromA1 origin must map to ToA2 origin: got={got:?} expected={expected:?}"
        );
        for (f, x) in [
            (from.x_direction, to.x_direction),
            (from.y_direction, to.y_direction),
            (from.direction(), to.direction()),
        ] {
            let g = t.transform_dir(f);
            assert!(
                (g - x).length() < 1e-12,
                "direction {f:?} maps to {g:?}, expected {x:?}"
            );
        }
    }

    /// OCCT gp_Trsf::Multiplied — the general branch of gp_Trsf::Multiply
    /// (gp_Trsf.cxx L544-559): `Tloc = matrix * T.loc; if (scale != 1.0)
    /// Tloc *= scale; loc += Tloc`.  The scale factor on the translated
    /// location was missing here before, and a non-unit scale alone exposes
    /// it: `A(B(p))` with `A` a pure scale of 3 and `B` a translation by
    /// (1,2,3) is `3p + (3,6,9)`, not `3p + (1,2,3)`.
    #[test]
    fn multiplied_scales_the_left_translation() {
        let mut a = Trsf::identity();
        a.set_scale_factor(3.0);
        let mut b = Trsf::identity();
        b.set_translation(DVec3::new(1.0, 2.0, 3.0));
        let c = a.multiplied(&b);
        let got = c.apply(DVec3::new(1.0, 1.0, 1.0));
        let expected = DVec3::new(6.0, 9.0, 12.0);
        assert!(
            (got - expected).length() < 1e-12,
            "got={got:?} expected={expected:?}"
        );
    }

    /// OCCT gp_Trsf::Value(Row, Col) (gp_Trsf.hxx L420-434): Col < 4 is the
    /// SCALED matrix entry, Col 4 the location coordinate.  Both halves are
    /// checked against a state built independently of `value`.
    #[test]
    fn value_folds_in_the_scale_and_reads_the_location() {
        let ang = 0.4_f64;
        let loc = DVec3::new(3.0, -4.0, 5.0);
        let mut t = Trsf::identity();
        t.set_rotation(&Ax1::new(DVec3::ZERO, DVec3::Z), ang);
        t.set_scale_factor(2.5);
        t.set_translation_part(loc);
        let (s, c) = ang.sin_cos();
        assert!((t.value(1, 1) - 2.5 * c).abs() < 1e-12);
        assert!((t.value(1, 2) - 2.5 * -s).abs() < 1e-12);
        assert!((t.value(2, 1) - 2.5 * s).abs() < 1e-12);
        assert!((t.value(3, 3) - 2.5).abs() < 1e-12);
        assert_eq!(t.value(1, 4), loc.x);
        assert_eq!(t.value(2, 4), loc.y);
        assert_eq!(t.value(3, 4), loc.z);
    }

    /// OCCT gp_Trsf::VectorialPart (gp_Trsf.hxx L464-471): the matrix with the
    /// scale FOLDED IN.  A non-unit scale is the discriminating input — the
    /// `scale == 1.0` shortcut must not swallow it.
    #[test]
    fn vectorial_part_folds_in_the_scale() {
        let ang = 0.4_f64;
        let (s, c) = ang.sin_cos();
        let mut t = Trsf::identity();
        t.set_rotation(&Ax1::new(DVec3::ZERO, DVec3::Z), ang);
        t.set_scale_factor(3.0);
        let m = t.vectorial_part();
        assert!((m[0][0] - 3.0 * c).abs() < 1e-12);
        assert!((m[0][1] - 3.0 * -s).abs() < 1e-12);
        assert!((m[1][0] - 3.0 * s).abs() < 1e-12);
        assert!((m[1][1] - 3.0 * c).abs() < 1e-12);
        // At scale 1 the matrix is returned unchanged.
        let mut r = Trsf::identity();
        r.set_rotation(&Ax1::new(DVec3::ZERO, DVec3::Z), ang);
        assert_eq!(r.vectorial_part(), r.matrix);
    }

    /// OCCT gp_Trsf::IsNegative (gp_Trsf.hxx L218): `scale < 0.0`.  The
    /// discriminating case is a plain negative SCALE (form gp_Scale), which an
    /// implementation testing `form == PntMirror` would report as positive.
    #[test]
    fn is_negative_reads_the_scale_not_the_form() {
        let mut positive = Trsf::identity();
        positive.set_scale_factor(2.0);
        assert!(!positive.is_negative());

        let mut negative = Trsf::identity();
        negative.set_scale_factor(-2.0);
        assert_eq!(negative.form, TrsfForm::Scale);
        assert!(negative.is_negative());

        let mut mirror = Trsf::identity();
        mirror.set_scale_factor(-1.0);
        assert_eq!(mirror.form, TrsfForm::PntMirror);
        assert!(mirror.is_negative());
    }

    /// OCCT gp_Trsf::SetTranslationPart(V) (gp_Trsf.cxx L244-281): loc = V and
    /// then the form transition.  The discriminating branch is a Rotation with
    /// a ZERO location — OCCT keeps gp_Rotation, while an "always demote to
    /// CompoundTrsf" implementation would fail.
    #[test]
    fn set_translation_part_follows_the_occt_form_transitions() {
        // Identity + non-zero V -> Translation.
        let mut t = Trsf::identity();
        t.set_translation_part(DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(t.form, TrsfForm::Translation);

        // Translation + zero V -> Identity.
        let mut t = Trsf::identity();
        t.set_translation(DVec3::new(1.0, 2.0, 3.0));
        t.set_translation_part(DVec3::ZERO);
        assert_eq!(t.form, TrsfForm::Identity);

        // Rotation + non-zero V -> CompoundTrsf.
        let mut t = Trsf::identity();
        t.set_rotation(&Ax1::new(DVec3::ZERO, DVec3::Z), 0.4);
        t.set_translation_part(DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(t.form, TrsfForm::CompoundTrsf);

        // Rotation + ZERO V keeps gp_Rotation.
        let mut t = Trsf::identity();
        t.set_rotation(&Ax1::new(DVec3::ZERO, DVec3::Z), 0.4);
        t.set_translation_part(DVec3::ZERO);
        assert_eq!(t.form, TrsfForm::Rotation);
    }

    /// OCCT gp_Trsf::SetScaleFactor (gp_Trsf.cxx L284-330): the form transition
    /// as a function of the starting form and of S in {1, -1, other}.  Built
    /// from the documented transition table rather than from the code.
    #[test]
    fn set_scale_factor_follows_the_occt_form_transitions() {
        fn rotation() -> Trsf {
            let mut t = Trsf::identity();
            t.set_rotation(&Ax1::new(DVec3::ZERO, DVec3::Z), 0.4);
            t
        }
        fn translation() -> Trsf {
            let mut t = Trsf::identity();
            t.set_translation(DVec3::new(1.0, 0.0, 0.0));
            t
        }
        fn scale() -> Trsf {
            let mut t = Trsf::identity();
            t.set_scale_factor(2.0);
            t
        }
        fn pnt_mirror() -> Trsf {
            let mut t = Trsf::identity();
            t.set_scale_factor(-1.0);
            t
        }
        fn compound() -> Trsf {
            let mut t = rotation();
            t.set_scale_factor(3.0);
            t
        }

        let cases: [(&str, fn() -> Trsf, f64, TrsfForm); 12] = [
            ("identity", Trsf::identity, 3.0, TrsfForm::Scale),
            ("identity", Trsf::identity, 1.0, TrsfForm::Identity),
            ("identity", Trsf::identity, -1.0, TrsfForm::PntMirror),
            ("translation", translation, 3.0, TrsfForm::Scale),
            ("translation", translation, -1.0, TrsfForm::PntMirror),
            ("rotation", rotation, 3.0, TrsfForm::CompoundTrsf),
            ("rotation", rotation, 1.0, TrsfForm::Rotation),
            ("scale", scale, 1.0, TrsfForm::Identity),
            ("scale", scale, -1.0, TrsfForm::PntMirror),
            ("pntMirror", pnt_mirror, 1.0, TrsfForm::Identity),
            ("pntMirror", pnt_mirror, 3.0, TrsfForm::Scale),
            ("compound", compound, 1.0, TrsfForm::CompoundTrsf),
        ];
        for (name, build, s, expected) in cases {
            let mut t = build();
            t.set_scale_factor(s);
            assert_eq!(t.form, expected, "start={name} S={s}");
            assert_eq!(t.scale, s, "start={name} S={s}");
        }
    }

    /// A free vector rotated by `ang` about the Z axis, straight from the
    /// Rodrigues formula (independent of every helper under test here).
    fn rot_free_about_z(v: DVec3, ang: f64) -> DVec3 {
        let par = DVec3::Z * v.dot(DVec3::Z);
        let perp = v - par;
        par + perp * ang.cos() + DVec3::Z.cross(perp) * ang.sin()
    }

    /// OCCT gp_Vec::Transform (gp_Vec.cxx L120-138) multiplies by
    /// `VectorialPart()` — the matrix with the scale FOLDED IN — while
    /// gp_Dir::Transform (gp_Dir.cxx L129-154) multiplies by
    /// `HVectorialPart()` — the RAW matrix — before normalizing.
    ///
    /// Swapping the two is only observable through the SIGN: normalization
    /// cancels a positive scale, so a non-unit-but-positive scale cannot tell
    /// them apart at all.  The discriminating case is a NEGATIVE scale, where
    /// the scaled matrix would flip the direction and then be flipped again by
    /// the `scale < 0` reversal at the end of gp_Dir::Transform.  This test
    /// therefore runs both signs.
    #[test]
    fn transform_vec_folds_in_scale_while_transform_dir_does_not() {
        let ang = 0.4_f64;
        for scale in [3.0_f64, -3.0_f64] {
            let mut t = Trsf::identity();
            t.set_rotation(&Ax1::new(DVec3::new(1.0, 0.0, 0.0), DVec3::Z), ang);
            t.set_scale_factor(scale);
            assert_eq!(t.form, TrsfForm::CompoundTrsf);
            assert_eq!(t.scale, scale);

            for &v in &[
                DVec3::new(1.0, 2.0, -0.5),
                DVec3::new(-3.0, 0.25, 7.0),
                DVec3::X,
            ] {
                let r = rot_free_about_z(v, ang);
                let got_vec = t.transform_vec(v);
                assert!(
                    (got_vec - r * scale).length() < 1e-12,
                    "transform_vec must fold in the scale: scale={scale} v={v:?} \
                     got={got_vec:?} expected={:?}",
                    r * scale
                );
                let expected_dir = if scale < 0.0 { -r } else { r }.normalize();
                let got_dir = t.transform_dir(v);
                assert!(
                    (got_dir - expected_dir).length() < 1e-12,
                    "transform_dir must stay a unit direction: scale={scale} v={v:?} \
                     got={got_dir:?} expected={expected_dir:?}"
                );
                assert!(
                    (got_dir.length() - 1.0).abs() < 1e-12,
                    "transform_dir returned a non-unit vector {got_dir:?}"
                );
            }
        }
    }

    /// OCCT gp_Trsf::Transforms (gp_Trsf.hxx L452-461) is `matrix*p`, then
    /// `*= scale` (only if scale != 1.0), then `+= loc`, in that order.  A
    /// non-unit scale distinguishes the order: the scale must touch the
    /// matrix term only, never the location.
    #[test]
    fn apply_scales_the_matrix_term_only() {
        let ang = 0.4_f64;
        let loc = DVec3::new(3.0, -4.0, 5.0);
        let mut t = Trsf::identity();
        t.set_rotation(&Ax1::new(DVec3::new(1.0, 0.0, 0.0), DVec3::Z), ang);
        t.set_scale_factor(2.5);
        t.set_translation_part(loc);

        let p = DVec3::new(2.0, -1.0, 0.5);
        let expected = rot_free_about_z(p, ang) * 2.5 + loc;
        let got = t.apply(p);
        assert!(
            (got - expected).length() < 1e-12,
            "got={got:?} expected={expected:?}"
        );
    }

    /// OCCT gp_Trsf::Invert (gp_Trsf.cxx L396-428), general branch:
    /// `scale = 1/scale; matrix.Transpose(); loc = matrix * loc; loc *= -scale`.
    /// The in-place transpose is legitimate here (OCCT does exactly that),
    /// unlike the `set_rotation` defect.  A non-unit scale plus an off-origin
    /// axis fails the round trip if any of those steps is wrong.
    #[test]
    fn invert_round_trips_scale_rotation_and_translation() {
        let mut t = Trsf::identity();
        t.set_rotation(&Ax1::new(DVec3::new(1.0, 0.0, 0.0), DVec3::Z), 0.4_f64);
        t.set_scale_factor(2.5);
        t.set_translation_part(DVec3::new(3.0, -4.0, 5.0));
        let mut back = t;
        back.invert();
        assert_eq!(back.scale, 1.0 / 2.5);
        for &p in &[
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(-7.0, 0.5, 11.0),
            DVec3::ZERO,
        ] {
            let q = back.apply(t.apply(p));
            assert!((q - p).length() < 1e-9, "p={p:?} round trip gave {q:?}");
        }
    }

    /// OCCT gp_Trsf::SetTransformation(const gp_Ax3&) (gp_Trsf.cxx L196-206)
    /// takes the Ax3 frame onto the standard frame: the Ax3 origin maps to
    /// (0,0,0) and origin + XDirection maps to +X (likewise Y and Z).  The
    /// Ax3 is moved and tilted so that a wrong row/column convention in the
    /// matrix would show up.
    #[test]
    fn set_transformation_ax3_maps_its_frame_to_the_standard_frame() {
        let ang = 0.7_f64;
        let a3 = Ax3::from_pnt_n_vx(
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(0.0, ang.cos(), ang.sin()),
            DVec3::X,
        );
        let t = Trsf::set_transformation_ax3(&a3);
        assert_eq!(t.form, TrsfForm::CompoundTrsf);
        let o = t.apply(a3.location());
        assert!(
            o.length() < 1e-12,
            "the Ax3 origin must map to the standard origin, got {o:?}"
        );
        let origin = a3.location();
        for (dir, axis) in [
            (a3.x_direction, DVec3::X),
            (a3.y_direction, DVec3::Y),
            (a3.direction(), DVec3::Z),
        ] {
            let got = t.apply(origin + dir);
            assert!(
                (got - axis).length() < 1e-12,
                "origin + {dir:?} must map to {axis:?}, got {got:?}"
            );
        }
    }

    /// OCCT gp_Trsf::SetTransformation(FromA1, ToA2) (gp_Trsf.cxx L172-194) is
    /// the COORDINATE-change transform: a point whose coordinates in FromA1 are
    /// `p` has coordinates `T(p)` in ToA2.  The expected values here are
    /// obtained by projecting onto each frame's own directions, which is
    /// independent of the implementation.  Note this is the INVERSE relation to
    /// `set_displacement` (gp_Trsf.hxx L136-137), so a test shape that suits
    /// one does NOT suit the other.
    #[test]
    fn set_transformation_from_to_changes_between_frames() {
        let ang = 0.7_f64;
        let from = Ax3::from_pnt_n_vx(
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(0.0, ang.cos(), ang.sin()),
            DVec3::X,
        );
        let to = Ax3::from_pnt_n_vx(DVec3::new(5.0, 7.0, 11.0), DVec3::Z, DVec3::Y);
        let t = Trsf::set_transformation_from_to(&from, &to);
        assert_eq!(t.form, TrsfForm::CompoundTrsf);
        assert_eq!(t.scale, 1.0);

        let coords = |a3: &Ax3, p: DVec3| {
            let d = p - a3.location();
            DVec3::new(
                d.dot(a3.x_direction),
                d.dot(a3.y_direction),
                d.dot(a3.direction()),
            )
        };

        for &p in &[
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(-4.0, 0.5, 8.0),
            DVec3::ZERO,
        ] {
            let got = t.apply(coords(&from, p));
            let expected = coords(&to, p);
            assert!(
                (got - expected).length() < 1e-12,
                "p={p:?} got={got:?} expected={expected:?}"
            );
        }

        // FromA1 == ToA2 must reduce to the identity.
        let same = Trsf::set_transformation_from_to(&from, &from);
        for &p in &[DVec3::new(3.0, -1.0, 0.25), DVec3::ZERO] {
            let got = same.apply(p);
            assert!(
                (got - p).length() < 1e-12,
                "the same-frame case gave {got:?} for {p:?}"
            );
        }
    }

    /// The glam affine carries the same transform as `apply`: the `DMat3` is
    /// built from the matrix COLUMNS, scaled, then translated.
    #[test]
    fn to_daffine3_agrees_with_apply() {
        let mut t = Trsf::identity();
        t.set_rotation(&Ax1::new(DVec3::new(1.0, 0.0, 0.0), DVec3::Z), 0.4_f64);
        t.set_scale_factor(2.5);
        t.set_translation_part(DVec3::new(3.0, -4.0, 5.0));
        let a = t.to_daffine3();
        for &p in &[DVec3::new(1.0, 2.0, 3.0), DVec3::new(-7.0, 0.5, 11.0)] {
            let got = a.transform_point3(p);
            let expected = t.apply(p);
            assert!(
                (got - expected).length() < 1e-9,
                "p={p:?} affine={got:?} apply={expected:?}"
            );
        }
    }

    /// Unit length and mutual orthogonality of an `Ax3` frame, plus the
    /// agreement of `direct()` with the frame it is derived from.
    fn assert_ax3_frame(a3: &Ax3) {
        for (name, v) in [
            ("x", a3.x_direction),
            ("y", a3.y_direction),
            ("axis", a3.axis.direction),
        ] {
            assert!((v.length() - 1.0).abs() < 1e-12, "{name} is not unit: {v:?}");
        }
        assert!(a3.x_direction.dot(a3.y_direction).abs() < 1e-12);
        assert!(a3.x_direction.dot(a3.axis.direction).abs() < 1e-12);
        assert!(a3.y_direction.dot(a3.axis.direction).abs() < 1e-12);
        let derived = a3.x_direction.cross(a3.y_direction).dot(a3.axis.direction) > 0.0;
        assert_eq!(a3.direct(), derived);
    }

    /// OCCT gp_Ax3::SetDirection, the PARALLEL branch (gp_Ax3.hxx L509-521):
    /// theV is (anti)parallel to the current X direction, so the
    /// `theV ^ (vxdir ^ theV)` re-orthogonalization would degenerate and OCCT
    /// relabels the axes instead.  The previous rcad version produced a NULL X
    /// direction here, so the frame assertions are the discriminating part.
    #[test]
    fn set_direction_handles_the_parallel_case() {
        // aDot > 0 (OCCT L511-515): X <- old Y, Y <- old axis.
        let mut a3 = Ax3::new();
        a3.set_direction(DVec3::X);
        assert_eq!(a3.axis.direction, DVec3::X);
        assert_eq!(a3.x_direction, DVec3::Y);
        assert_eq!(a3.y_direction, DVec3::Z);
        assert_ax3_frame(&a3);

        // aDot < 0 (OCCT L516-519): only X is relabelled.
        let mut a3 = Ax3::new();
        a3.set_direction(-DVec3::X);
        assert_eq!(a3.axis.direction, -DVec3::X);
        assert_eq!(a3.x_direction, DVec3::Z);
        assert_eq!(a3.y_direction, DVec3::Y);
        assert_ax3_frame(&a3);
    }

    /// OCCT gp_Ax3::SetDirection, the general branch (gp_Ax3.hxx L523-533):
    /// `vxdir = theV ^ (vxdir ^ theV)`, then Y = theV ^ vxdir for a direct frame.
    #[test]
    fn set_direction_reorthogonalizes_x_in_the_general_case() {
        let mut a3 = Ax3::new();
        let v = DVec3::new(1.0, 1.0, 1.0).normalize();
        a3.set_direction(v);
        let expected_x = v.cross(DVec3::X.cross(v)).normalize();
        assert!((a3.x_direction - expected_x).length() < 1e-12);
        assert!((a3.y_direction - v.cross(expected_x)).length() < 1e-12);
        assert!(a3.direct());
        assert_ax3_frame(&a3);
    }

    /// OCCT gp_Ax3::Direct (gp_Ax3.hxx L208) is DERIVED from the frame, not a
    /// stored flag: an indirect frame must report `false`, and SetDirection must
    /// then take the `vxdir ^ theV` ordering (gp_Ax3.hxx L531-533).  rcad's
    /// removed `sense` field was always `true`, so it could not reach this.
    #[test]
    fn direct_is_derived_and_drives_the_indirect_branch() {
        let mut a3 = Ax3::new();
        // The frame fields are public, so an indirect frame is built directly.
        a3.y_direction = -DVec3::Y;
        assert!(!a3.direct());
        let v = DVec3::new(1.0, 1.0, 0.0).normalize();
        a3.set_direction(v);
        let expected_x = v.cross(DVec3::X.cross(v)).normalize();
        assert!((a3.x_direction - expected_x).length() < 1e-12);
        // Indirect: Y = vxdir ^ theV instead of theV ^ vxdir.
        assert!((a3.y_direction - expected_x.cross(v)).length() < 1e-12);
        assert_ax3_frame(&a3);
    }

    /// OCCT gp_Ax3::SetXDirection, the PARALLEL branch (gp_Ax3.hxx L543-555).
    #[test]
    fn set_x_direction_handles_the_parallel_case() {
        // aDot > 0 (OCCT L545-549): axis <- old X, Y flips, X <- theVx.
        let mut a3 = Ax3::new();
        a3.set_x_direction(DVec3::Z);
        assert_eq!(a3.axis.direction, DVec3::X);
        assert_eq!(a3.x_direction, DVec3::Z);
        assert_eq!(a3.y_direction, -DVec3::Y);
        assert_ax3_frame(&a3);

        // aDot < 0 (OCCT L550-553): Y keeps its direction.
        let mut a3 = Ax3::new();
        a3.set_x_direction(-DVec3::Z);
        assert_eq!(a3.axis.direction, DVec3::X);
        assert_eq!(a3.x_direction, -DVec3::Z);
        assert_eq!(a3.y_direction, DVec3::Y);
        assert_ax3_frame(&a3);
    }
}
