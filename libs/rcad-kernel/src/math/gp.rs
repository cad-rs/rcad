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
    /// OCCT gp_Ax3::sense — true for right-handed ("direct") systems.
    sense: bool,
}

impl Ax3 {
    /// OCCT gp_Ax3() — default right-handed OXYZ.
    pub fn new() -> Self {
        Ax3 {
            axis: Ax1::new(DVec3::ZERO, DVec3::Z),
            y_direction: DVec3::Y,
            x_direction: DVec3::X,
            sense: true,
        }
    }

    /// OCCT gp_Ax3(P, N, Vx) — right-handed system (via gp_Ax2).
    pub fn from_pnt_n_vx(location: DVec3, direction: DVec3, x_direction: DVec3) -> Self {
        let a2 = Ax2::new(location, direction, x_direction);
        Ax3 {
            axis: a2.axis(),
            y_direction: a2.y_direction,
            x_direction: a2.x_direction,
            sense: true,
        }
    }

    /// OCCT gp_Ax3(const gp_Ax2& theA) — right-handed.
    pub fn from_ax2(a: &Ax2) -> Self {
        Ax3 {
            axis: a.axis(),
            y_direction: a.y_direction,
            x_direction: a.x_direction,
            sense: true,
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

    /// OCCT gp_Ax3::Direct.
    pub fn direct(&self) -> bool {
        self.sense
    }

    /// OCCT gp_Ax3::SetLocation.
    pub fn set_location(&mut self, p: DVec3) {
        self.axis.location = p;
    }

    /// OCCT gp_Ax3::SetDirection — keeps sense, recomputes X (`V ^ (X ^ V)`)
    /// and Y (`N ^ X`).
    pub fn set_direction(&mut self, v: DVec3) {
        let v = v.normalize_or_zero();
        let old_x = self.x_direction;
        let xd = v.cross(old_x.cross(v));
        let xd = if xd.length_squared() < f64::EPSILON {
            DVec3::ZERO
        } else {
            xd.normalize_or_zero()
        };
        self.axis.direction = v;
        self.x_direction = xd;
        self.y_direction = v.cross(xd).normalize_or_zero();
    }

    /// OCCT gp_Ax3::SetXDirection — Y recomputed as N ^ new X.
    pub fn set_x_direction(&mut self, vx: DVec3) {
        let n = self.axis.direction;
        let xd = n.cross(vx.cross(n));
        let xd = if xd.length_squared() < f64::EPSILON {
            DVec3::ZERO
        } else {
            xd.normalize_or_zero()
        };
        self.x_direction = xd;
        self.y_direction = n.cross(xd).normalize_or_zero();
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

    /// OCCT gp_Lin::Transform(theT) — gp_Lin.hxx L178: only the location is
    /// transformed (this OCCT version does not touch the direction here).
    pub fn transform(&mut self, t: &Trsf) {
        self.pos = t.apply(self.pos);
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
        // MA1loc = FromA1.Location(); MA1loc.Multiply(MA1) — row-vector times
        // matrix: r[j] = v . (column j of MA1).
        let mut ma1_loc = [0.0f64; 3];
        for j in 0..3 {
            ma1_loc[j] = from_a1.location().dot(DVec3::new(
                ma1[0][j],
                ma1[1][j],
                ma1[2][j],
            ));
        }
        // MA1loc.Reverse().
        for v in ma1_loc.iter_mut() {
            *v = -*v;
        }
        // MA1loc.Multiply(matrix) — again row-vector times matrix, in place.
        let mut mult = [0.0f64; 3];
        for i in 0..3 {
            mult[i] =
                ma1_loc[0] * matrix[0][i] + ma1_loc[1] * matrix[1][i] + ma1_loc[2] * matrix[2][i];
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
        // Combined transform: A(B(p)) = (A*B)*p + A.matrix*B.loc + A.loc.
        let mut loc = self.loc;
        for i in 0..3 {
            loc[i] += self.matrix[i][0] * t.loc.x
                + self.matrix[i][1] * t.loc.y
                + self.matrix[i][2] * t.loc.z;
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
