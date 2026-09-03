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
}

/// OCCT gp_Trsf — rigid transformation as built by `gp_Trsf::SetDisplacement`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Trsf {
    /// Row-major 3x3: `matrix[i][j]` == OCCT `Value(i+1, j+1)`.
    pub matrix: [[f64; 3]; 3],
    /// Translation part.
    pub loc: DVec3,
}

impl Trsf {
    /// OCCT gp_Trsf::SetDisplacement(FromA1, ToA2) — gp_Trsf.cxx L218-240,
    /// with the same floating-point operation order.
    pub fn set_displacement(from_a1: &Ax3, to_a2: &Ax3) -> Self {
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
        }
    }

    /// OCCT gp_Trsf::Multiplied(theT) — `this * theT`.
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
        Trsf { matrix, loc }
    }

    /// OCCT gp_XYZ::Transform(theT): `r[i] = Sum_j Value(i, j) * p[j] + loc[i]`.
    pub fn apply(&self, p: DVec3) -> DVec3 {
        DVec3::new(
            self.matrix[0][0] * p.x + self.matrix[0][1] * p.y + self.matrix[0][2] * p.z
                + self.loc.x,
            self.matrix[1][0] * p.x + self.matrix[1][1] * p.y + self.matrix[1][2] * p.z
                + self.loc.y,
            self.matrix[2][0] * p.x + self.matrix[2][1] * p.y + self.matrix[2][2] * p.z
                + self.loc.z,
        )
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
        );
        glam::DAffine3::from_mat3_translation(m, self.loc)
    }
}
