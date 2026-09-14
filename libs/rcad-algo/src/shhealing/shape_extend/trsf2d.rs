//! Local carrier for OCCT `gp_Trsf2d` (TKMath, `gp_Trsf2d.hxx` /
//! `gp_Trsf2d.cxx`), translated 1:1.
//!
//! This block used to be the tail of `composite_surface.rs`; it moved to its
//! own module to keep both files under the 2000-line guideline.  The carrier
//! is exercised by `ShapeExtendCompositeSurface::global_to_local_transformation`
//! (`super::composite_surface`), which reads the form tag and composes the
//! scale and shift transformations.
//!
//! Translation notes:
//! - `gp_Mat2d` maps to the row-major `[[f64; 2]; 2]` in which `m[i][j]` is
//!   `gp_Mat2d::Value(i + 1, j + 1)`.
//! - `gp_XY` maps to `glam::DVec2`.
//! - `gp_Ax2d` maps to the (origin, direction) pair, as in `set_mirror_ax2d`.
//! - The OCCT `Standard_ConstructionError` raises panic here.

use glam::DVec2;
use rcad_kernel::math::gp::GP_RESOLUTION;

/// OCCT gp_TrsfForm (gp_TrsfForm.hxx L21-30), restricted to the values a
/// `gp_Trsf2d` can actually hold.  The 2D class never produces gp_Ax2Mirror
/// (bilateral symmetry exists only in the 3D class), so that value has no
/// counterpart here; `Other` carries the OCCT `gp_Other` fallback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrsfForm {
    Identity,
    Rotation,
    Translation,
    PntMirror,
    Ax1Mirror,
    Scale,
    CompoundTrsf,
    Other,
}

/// 1:1 carrier for OCCT `gp_Trsf2d` (TKMath, `gp_Trsf2d.hxx` / `gp_Trsf2d.cxx`):
/// a 2D transformation with a separate scale factor and form tag, following the
/// gp_Trsf member layout (scale, shape, matrix, loc).  The member names follow
/// OCCT's private members: `matrix` is the 2x2 part WITHOUT the scale folded in
/// (`myMat`), `loc` the translation part, `scale` the gp_Trsf scale factor and
/// `form` the gp_TrsfForm tag.
///
/// Point transformation follows gp_Trsf2d::Transforms (gp_Trsf2d.hxx
/// L316-338): `P' = matrix * P * scale + loc`.
/// OCCT gp_XY::Multiply(const gp_Mat2d&) (gp_XY.hxx L401-406):
/// `<me> = theMatrix * <me>` — the matrix acts on the LEFT, not the transpose.
fn mat2_mul_xy(m: [[f64; 2]; 2], v: DVec2) -> DVec2 {
    DVec2::new(m[0][0] * v.x + m[0][1] * v.y, m[1][0] * v.x + m[1][1] * v.y)
}

/// OCCT gp_Mat2d::Multiply(const gp_Mat2d&) (gp_Mat2d.hxx L326-334):
/// `this = this * theOther`.
fn mat2_mul(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

/// OCCT gp_Mat2d::Transpose() (gp_Mat2d.hxx L394-399): `A(j, i) -> A(i, j)`.
fn mat2_transpose(m: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [[m[0][0], m[1][0]], [m[0][1], m[1][1]]]
}

/// OCCT gp_XY::Normalize() (gp_XY.hxx L410-417): divides both coordinates by
/// `Modulus()`; the `Standard_ConstructionError` of a modulus below
/// gp::Resolution() becomes a panic.
fn gp_xy_normalize(v: DVec2) -> DVec2 {
    let a_d = v.length(); // gp_XY::Modulus().
    if a_d <= GP_RESOLUTION {
        panic!("gp_XY::Normalize() - vector has zero norm");
    }
    DVec2::new(v.x / a_d, v.y / a_d)
}

#[derive(Debug, Clone, Copy)]
pub struct Trsf2d {
    /// Row-major 2x2 matrix (WITHOUT the scale folded in, like OCCT).
    matrix: [[f64; 2]; 2],
    /// Translation part.
    loc: DVec2,
    /// gp_Trsf scale factor.
    scale: f64,
    /// gp_TrsfForm tag.
    form: TrsfForm,
}

impl Default for Trsf2d {
    /// OCCT gp_Trsf2d default constructor: scale = 1, gp_Identity, identity
    /// matrix, zero location.
    fn default() -> Self {
        Trsf2d {
            matrix: [[1.0, 0.0], [0.0, 1.0]],
            loc: DVec2::ZERO,
            scale: 1.0,
            form: TrsfForm::Identity,
        }
    }
}

impl Trsf2d {
    /// OCCT gp_Trsf::SetTranslation(gp_Vec) (gp_Trsf.cxx): scale = 1,
    /// identity matrix, translation part = V, form = gp_Translation.
    pub fn set_translation(&mut self, v: DVec2) {
        self.matrix = [[1.0, 0.0], [0.0, 1.0]];
        self.loc = v;
        self.scale = 1.0;
        self.form = TrsfForm::Translation;
    }

    /// OCCT gp_Trsf2d::SetScale(const gp_Pnt2d&, double) (gp_Trsf2d.hxx
    /// L270-277): form = gp_Scale, scale = S, identity matrix,
    /// loc = P * (1 - S). The 2D class raises for no scale value: the
    /// Standard_ConstructionError of gp_Trsf::SetScale (gp_Trsf.cxx L164-165)
    /// belongs to the 3D class only.
    pub fn set_scale(&mut self, p: DVec2, s: f64) {
        self.form = TrsfForm::Scale;
        self.scale = s;
        self.matrix = [[1.0, 0.0], [0.0, 1.0]];
        self.loc = p;
        self.loc *= 1.0 - s;
    }

    /// OCCT gp_Trsf2d::SetMirror(const gp_Pnt2d&) (gp_Trsf2d.hxx L259-266):
    /// the central symmetry about `p` — scale = -1, identity matrix,
    /// loc = P, loc.Multiply(2.0).
    pub fn set_mirror_pnt(&mut self, p: DVec2) {
        self.form = TrsfForm::PntMirror;
        self.scale = -1.0;
        self.matrix = [[1.0, 0.0], [0.0, 1.0]];
        self.loc = p;
        self.loc *= 2.0;
    }

    /// OCCT gp_Trsf2d::SetMirror(const gp_Ax2d&) (gp_Trsf2d.cxx L31-46): the
    /// rotational symmetry about the line `A` — scale = -1 and the reflection
    /// matrix `I - 2*V*V^T` written column by column.
    pub fn set_mirror_ax2d(&mut self, origin: DVec2, dir: DVec2) {
        self.form = TrsfForm::Ax1Mirror;
        self.scale = -1.0;
        let vx = dir.x;
        let vy = dir.y;
        let x0 = origin.x;
        let y0 = origin.y;
        // matrix.SetCol(1, (1 - 2 vx^2, -2 vx vy));
        // matrix.SetCol(2, (-2 vx vy, 1 - 2 vy^2)).
        self.matrix = [
            [1.0 - 2.0 * vx * vx, -2.0 * vx * vy],
            [-2.0 * vx * vy, 1.0 - 2.0 * vy * vy],
        ];
        // loc.SetCoord(-2*((vx^2 - 1) x0 + (vx vy y0)),
        //              -2*((vx vy x0) + (vy^2 - 1) y0)).
        self.loc = DVec2::new(
            -2.0 * ((vx * vx - 1.0) * x0 + (vx * vy * y0)),
            -2.0 * ((vx * vy * x0) + (vy * vy - 1.0) * y0),
        );
    }

    /// OCCT gp_Trsf2d::SetRotation(const gp_Pnt2d&, double) (gp_Trsf2d.hxx
    /// L246-257): the rotation of `ang` radians about the centre `p` —
    /// `loc = P; loc.Reverse(); matrix.SetRotation(ang); loc.Multiply(matrix);
    /// loc.Add(P)`.
    pub fn set_rotation(&mut self, p: DVec2, ang: f64) {
        self.form = TrsfForm::Rotation;
        self.scale = 1.0;
        self.loc = -p;
        // gp_Mat2d::SetRotation (gp_Mat2d.hxx L271-278).
        let (sin, cos) = ang.sin_cos();
        self.matrix = [[cos, -sin], [sin, cos]];
        self.loc = mat2_mul_xy(self.matrix, self.loc);
        self.loc += p;
    }

    /// OCCT gp_Trsf2d::Value(Row, Col) (gp_Trsf2d.hxx L301-314): Col < 3
    /// returns scale * myMat[Row-1][Col-1], Col 3 the location coordinate.
    pub fn value(&self, r: usize, c: usize) -> f64 {
        if c < 3 {
            self.scale * self.matrix[r - 1][c - 1]
        } else {
            [self.loc.x, self.loc.y][r - 1]
        }
    }

    /// OCCT gp_Trsf::Form().
    pub fn form(&self) -> TrsfForm {
        self.form
    }

    /// OCCT gp_Trsf2d::Multiply(const gp_Trsf2d& T) (gp_Trsf2d.cxx L253-382) —
    /// `this = this * T`, so T applies first.  All fifteen branches are
    /// reproduced in OCCT's order and the form tags follow the C++ code
    /// exactly.
    pub fn multiplied(&self, t: &Trsf2d) -> Trsf2d {
        let mut res = *self;
        if t.form == TrsfForm::Identity {
            // OCCT L255-257: T identity -> this unchanged.
        } else if res.form == TrsfForm::Identity {
            // OCCT L258-264: this identity -> copy T.
            res.form = t.form;
            res.scale = t.scale;
            res.loc = t.loc;
            res.matrix = t.matrix;
        } else if res.form == TrsfForm::Rotation && t.form == TrsfForm::Rotation {
            // OCCT L265-272.
            if res.loc.x != 0.0 || res.loc.y != 0.0 {
                // loc.Add(T.loc.Multiplied(matrix)).
                res.loc += mat2_mul_xy(res.matrix, t.loc);
            }
            // matrix.Multiply(T.matrix).
            res.matrix = mat2_mul(res.matrix, t.matrix);
        } else if res.form == TrsfForm::Translation && t.form == TrsfForm::Translation {
            // OCCT L273-276: loc.Add(T.loc).
            res.loc += t.loc;
        } else if res.form == TrsfForm::Scale && t.form == TrsfForm::Scale {
            // OCCT L277-281: loc.Add(T.loc.Multiplied(scale));
            // scale = scale * T.scale.
            res.loc += t.loc * res.scale;
            res.scale = res.scale * t.scale;
        } else if res.form == TrsfForm::PntMirror && t.form == TrsfForm::PntMirror {
            // OCCT L282-287: scale = 1.0; shape = gp_Translation;
            // loc.Add(T.loc.Reversed()).
            res.scale = 1.0;
            res.form = TrsfForm::Translation;
            res.loc += -t.loc;
        } else if res.form == TrsfForm::Ax1Mirror && t.form == TrsfForm::Ax1Mirror {
            // OCCT L288-297.
            res.form = TrsfForm::Rotation;
            let mut tloc = mat2_mul_xy(res.matrix, t.loc);
            tloc *= res.scale;
            res.scale = res.scale * t.scale;
            res.loc += tloc;
            res.matrix = mat2_mul(res.matrix, t.matrix);
        } else if (res.form == TrsfForm::CompoundTrsf
            || res.form == TrsfForm::Rotation
            || res.form == TrsfForm::Ax1Mirror)
            && t.form == TrsfForm::Translation
        {
            // OCCT L298-308.
            let mut tloc = mat2_mul_xy(res.matrix, t.loc);
            if res.scale != 1.0 {
                tloc *= res.scale;
            }
            res.loc += tloc;
        } else if (res.form == TrsfForm::Scale || res.form == TrsfForm::PntMirror)
            && t.form == TrsfForm::Translation
        {
            // OCCT L309-314: Tloc = T.loc * scale.
            res.loc += t.loc * res.scale;
        } else if res.form == TrsfForm::Translation
            && (t.form == TrsfForm::CompoundTrsf
                || t.form == TrsfForm::Rotation
                || t.form == TrsfForm::Ax1Mirror)
        {
            // OCCT L315-322.
            res.form = TrsfForm::CompoundTrsf;
            res.scale = t.scale;
            res.loc += t.loc;
            res.matrix = t.matrix;
        } else if res.form == TrsfForm::Translation
            && (t.form == TrsfForm::Scale || t.form == TrsfForm::PntMirror)
        {
            // OCCT L323-329: shape = T.shape.
            res.form = t.form;
            res.loc += t.loc;
            res.scale = t.scale;
        } else if (res.form == TrsfForm::PntMirror || res.form == TrsfForm::Scale)
            && (t.form == TrsfForm::PntMirror || t.form == TrsfForm::Scale)
        {
            // OCCT L330-339.
            res.form = TrsfForm::CompoundTrsf;
            res.loc += t.loc * res.scale;
            res.scale = res.scale * t.scale;
        } else if (res.form == TrsfForm::CompoundTrsf
            || res.form == TrsfForm::Rotation
            || res.form == TrsfForm::Ax1Mirror)
            && (t.form == TrsfForm::Scale || t.form == TrsfForm::PntMirror)
        {
            // OCCT L340-356.
            res.form = TrsfForm::CompoundTrsf;
            let mut tloc = mat2_mul_xy(res.matrix, t.loc);
            if res.scale == 1.0 {
                res.scale = t.scale;
            } else {
                tloc *= res.scale;
                res.scale = res.scale * t.scale;
            }
            res.loc += tloc;
        } else if (t.form == TrsfForm::CompoundTrsf
            || t.form == TrsfForm::Rotation
            || t.form == TrsfForm::Ax1Mirror)
            && (res.form == TrsfForm::Scale || res.form == TrsfForm::PntMirror)
        {
            // OCCT L357-368.
            res.form = TrsfForm::CompoundTrsf;
            res.loc += t.loc * res.scale;
            res.scale = res.scale * t.scale;
            res.matrix = t.matrix;
        } else {
            // OCCT L369-382: the general branch.
            res.form = TrsfForm::CompoundTrsf;
            let mut tloc = mat2_mul_xy(res.matrix, t.loc);
            if res.scale != 1.0 {
                tloc *= res.scale;
                res.scale = res.scale * t.scale;
            } else {
                res.scale = t.scale;
            }
            res.loc += tloc;
            res.matrix = mat2_mul(res.matrix, t.matrix);
        }
        res
    }

    /// OCCT gp_Trsf2d::Transformed(P): `P' = matrix * P * scale + loc`.
    pub fn transformed(&self, p: DVec2) -> DVec2 {
        let m = self.matrix;
        DVec2::new(
            (m[0][0] * p.x + m[0][1] * p.y) * self.scale + self.loc.x,
            (m[1][0] * p.x + m[1][1] * p.y) * self.scale + self.loc.y,
        )
    }

    /// OCCT gp_Trsf2d::SetTransformation(const gp_Ax2d& FromA1,
    /// const gp_Ax2d& ToA2) (gp_Trsf2d.cxx L48-70): the transformation
    /// changing from the coordinate system FromA1 to ToA2.  As everywhere in
    /// this carrier, a gp_Ax2d is the (origin, direction) pair.
    pub fn set_transformation_from_to(
        &mut self,
        from_origin: DVec2,
        from_dir: DVec2,
        to_origin: DVec2,
        to_dir: DVec2,
    ) {
        // OCCT L50-51: shape = gp_CompoundTrsf; scale = 1.0.
        self.form = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        // matrix from XOY to A2 :
        let v1 = to_dir;
        let v2 = DVec2::new(-v1.y, v1.x);
        // matrix.SetCol(1, V1); matrix.SetCol(2, V2).
        self.matrix = [[v1.x, v2.x], [v1.y, v2.y]];
        self.loc = to_origin;
        self.matrix = mat2_transpose(self.matrix);
        self.loc = mat2_mul_xy(self.matrix, self.loc);
        self.loc = -self.loc;
        // matrix FromA1 to XOY
        let v3 = from_dir;
        let v4 = DVec2::new(-v3.y, v3.x);
        // gp_Mat2d MA1(V3, V4) — the constructor takes the two columns.
        let ma1 = [[v3.x, v4.x], [v3.y, v4.y]];
        let mut ma1loc = from_origin;
        // matrix * MA1 => FromA1 ToA2
        ma1loc = mat2_mul_xy(self.matrix, ma1loc);
        self.loc += ma1loc;
        self.matrix = mat2_mul(self.matrix, ma1);
    }

    /// OCCT gp_Trsf2d::SetTransformation(const gp_Ax2d& A) (gp_Trsf2d.cxx
    /// L72-84): the transformation from the basic coordinate system to the
    /// local coordinate system A.
    pub fn set_transformation_ax2d(&mut self, origin: DVec2, dir: DVec2) {
        // OCCT L74-75: shape = gp_CompoundTrsf; scale = 1.0.
        self.form = TrsfForm::CompoundTrsf;
        self.scale = 1.0;
        let v1 = dir;
        let v2 = DVec2::new(-v1.y, v1.x);
        // matrix.SetCol(1, V1); matrix.SetCol(2, V2).
        self.matrix = [[v1.x, v2.x], [v1.y, v2.y]];
        self.loc = origin;
        self.matrix = mat2_transpose(self.matrix);
        self.loc = mat2_mul_xy(self.matrix, self.loc);
        self.loc = -self.loc;
    }

    /// OCCT gp_Trsf2d::SetTranslationPart(const gp_Vec2d& V) (gp_Trsf2d.cxx
    /// L86-117): replaces the translation vector and re-tags the form.
    pub fn set_translation_part(&mut self, v: DVec2) {
        self.loc = v;
        if self.loc.x.abs() <= GP_RESOLUTION && self.loc.y.abs() <= GP_RESOLUTION {
            if self.form == TrsfForm::Identity
                || self.form == TrsfForm::PntMirror
                || self.form == TrsfForm::Scale
                || self.form == TrsfForm::Rotation
                || self.form == TrsfForm::Ax1Mirror
            {
                // OCCT L91-94.
            } else if self.form == TrsfForm::Translation {
                // OCCT L95-98.
                self.form = TrsfForm::Identity;
            } else {
                // OCCT L99-102.
                self.form = TrsfForm::CompoundTrsf;
            }
        } else if self.form == TrsfForm::Translation
            || self.form == TrsfForm::Scale
            || self.form == TrsfForm::PntMirror
        {
            // OCCT L106-108.
        } else if self.form == TrsfForm::Identity {
            // OCCT L109-112.
            self.form = TrsfForm::Translation;
        } else {
            // OCCT L113-116.
            self.form = TrsfForm::CompoundTrsf;
        }
    }

    /// OCCT gp_Trsf2d::SetScaleFactor(const double S) (gp_Trsf2d.cxx
    /// L120-196): changes the scale factor and re-tags the form.  Every
    /// branch and its order follow the C++ source; the tag transitions are
    /// deliberately not "cleaned up" — OCCT leaves some states tagged with a
    /// form that no longer describes the map.
    pub fn set_scale_factor(&mut self, s: f64) {
        if s == 1.0 {
            let mut x = self.loc.x;
            if x < 0.0 {
                x = -x;
            }
            let mut y = self.loc.y;
            if y < 0.0 {
                y = -y;
            }
            if x <= GP_RESOLUTION && y <= GP_RESOLUTION {
                if self.form == TrsfForm::Identity || self.form == TrsfForm::Rotation {
                    // OCCT L136-138.
                } else if self.form == TrsfForm::Scale {
                    // OCCT L139-142.
                    self.form = TrsfForm::Identity;
                } else if self.form == TrsfForm::PntMirror {
                    // OCCT L143-146.
                    self.form = TrsfForm::Translation;
                } else {
                    // OCCT L147-150.
                    self.form = TrsfForm::CompoundTrsf;
                }
            } else if self.form == TrsfForm::Identity
                || self.form == TrsfForm::Rotation
                || self.form == TrsfForm::Scale
            {
                // OCCT L154-156.
            } else if self.form == TrsfForm::PntMirror {
                // OCCT L157-160.
                self.form = TrsfForm::Translation;
            } else {
                // OCCT L161-164.
                self.form = TrsfForm::CompoundTrsf;
            }
        } else if s == -1.0 {
            if self.form == TrsfForm::PntMirror || self.form == TrsfForm::Ax1Mirror {
                // OCCT L169-171.
            } else if self.form == TrsfForm::Identity || self.form == TrsfForm::Scale {
                // OCCT L172-175.
                self.form = TrsfForm::PntMirror;
            } else {
                // OCCT L176-179.
                self.form = TrsfForm::CompoundTrsf;
            }
        } else if self.form == TrsfForm::Scale {
            // OCCT L183-185.
        } else if self.form == TrsfForm::Identity
            || self.form == TrsfForm::Translation
            || self.form == TrsfForm::PntMirror
        {
            // OCCT L186-189.
            self.form = TrsfForm::Scale;
        } else {
            // OCCT L190-193.
            self.form = TrsfForm::CompoundTrsf;
        }
        self.scale = s;
    }

    /// OCCT gp_Trsf2d::VectorialPart() (gp_Trsf2d.cxx L198-214): the 2x2
    /// matrix with the scale factor folded in.
    pub fn vectorial_part(&self) -> [[f64; 2]; 2] {
        if self.scale == 1.0 {
            // OCCT L200-203.
            return self.matrix;
        }
        let mut m = self.matrix;
        if self.form == TrsfForm::Scale || self.form == TrsfForm::PntMirror {
            // OCCT L205-208: M.SetDiagonal(matrix.Value(1, 1) * scale,
            // matrix.Value(2, 2) * scale) — only the main diagonal.
            m[0][0] = self.matrix[0][0] * self.scale;
            m[1][1] = self.matrix[1][1] * self.scale;
        } else {
            // OCCT L209-212: M.Multiply(scale), in gp_Mat2d order.
            m[0][0] *= self.scale;
            m[0][1] *= self.scale;
            m[1][0] *= self.scale;
            m[1][1] *= self.scale;
        }
        m
    }

    /// OCCT gp_Trsf2d::RotationPart() (gp_Trsf2d.cxx L216-219): the angle of
    /// the rotation, read from matrix.Value(2, 1) / matrix.Value(1, 1).
    pub fn rotation_part(&self) -> f64 {
        self.matrix[1][0].atan2(self.matrix[0][0])
    }

    /// OCCT gp_Trsf2d::Invert() (gp_Trsf2d.cxx L221-251).
    pub fn invert(&mut self) {
        //                                    -1
        //  X' = scale * R * X + T  =>  X = (R  / scale)  * ( X' - T)
        //
        // For a gp_Trsf2d, since the scale is extracted from the matrix R,
        // determinant (R) = 1 always and R-1 = R transposed.
        if self.form == TrsfForm::Identity {
            // OCCT L228-230.
        } else if self.form == TrsfForm::Translation || self.form == TrsfForm::PntMirror {
            // OCCT L231-234: loc.Reverse().
            self.loc = -self.loc;
        } else if self.form == TrsfForm::Scale {
            // OCCT L235-241.
            if self.scale.abs() <= GP_RESOLUTION {
                panic!("gp_Trsf2d::Invert() - transformation has zero scale");
            }
            self.scale = 1.0 / self.scale;
            self.loc *= -self.scale;
        } else {
            // OCCT L242-250.
            if self.scale.abs() <= GP_RESOLUTION {
                panic!("gp_Trsf2d::Invert() - transformation has zero scale");
            }
            self.scale = 1.0 / self.scale;
            self.matrix = mat2_transpose(self.matrix);
            self.loc = mat2_mul_xy(self.matrix, self.loc);
            self.loc *= -self.scale;
        }
    }

    /// OCCT gp_Trsf2d::Power(const int N) (gp_Trsf2d.cxx L384-548): the
    /// N-fold composition of <me> with itself.  The form dispatch happens
    /// after the optional inversion because Invert() itself can re-tag the
    /// form (gp_Ax1Mirror becomes gp_Rotation).
    pub fn power(&mut self, n: i32) {
        if self.form == TrsfForm::Identity {
            // OCCT L386-389.
        } else if n == 0 {
            // OCCT L391-397.
            self.scale = 1.0;
            self.form = TrsfForm::Identity;
            self.matrix = [[1.0, 0.0], [0.0, 1.0]]; // matrix.SetIdentity().
            self.loc = DVec2::new(0.0, 0.0);
        } else if n == 1 {
            // OCCT L398-400.
        } else if n == -1 {
            // OCCT L401-404.
            self.invert();
        } else {
            // OCCT L405-411.
            if n < 0 {
                self.invert();
            }
            if self.form == TrsfForm::Translation {
                // OCCT L411-433.
                let mut npower = n;
                if npower < 0 {
                    npower = -npower;
                }
                npower -= 1;
                let mut temploc = self.loc;
                loop {
                    if npower % 2 == 1 {
                        self.loc += temploc;
                    }
                    if npower == 1 {
                        break;
                    }
                    temploc += temploc;
                    npower = npower / 2;
                }
            } else if self.form == TrsfForm::Scale {
                // OCCT L434-459.
                let mut npower = n;
                if npower < 0 {
                    npower = -npower;
                }
                npower -= 1;
                let mut temploc = self.loc;
                let mut tempscale = self.scale;
                loop {
                    if npower % 2 == 1 {
                        self.loc += temploc * self.scale;
                        self.scale = self.scale * tempscale;
                    }
                    if npower == 1 {
                        break;
                    }
                    temploc += temploc * tempscale;
                    tempscale = tempscale * tempscale;
                    npower = npower / 2;
                }
            } else if self.form == TrsfForm::Rotation {
                // OCCT L460-504.
                let mut npower = n;
                if npower < 0 {
                    npower = -npower;
                }
                npower -= 1;
                let mut tempmatrix = self.matrix;
                if self.loc.x == 0.0 && self.loc.y == 0.0 {
                    // OCCT L469-484.
                    loop {
                        if npower % 2 == 1 {
                            self.matrix = mat2_mul(self.matrix, tempmatrix);
                        }
                        if npower == 1 {
                            break;
                        }
                        tempmatrix = mat2_mul(tempmatrix, tempmatrix);
                        npower = npower / 2;
                    }
                } else {
                    // OCCT L485-503.
                    let mut temploc = self.loc;
                    loop {
                        if npower % 2 == 1 {
                            self.loc += mat2_mul_xy(self.matrix, temploc);
                            self.matrix = mat2_mul(self.matrix, tempmatrix);
                        }
                        if npower == 1 {
                            break;
                        }
                        temploc += mat2_mul_xy(tempmatrix, temploc);
                        tempmatrix = mat2_mul(tempmatrix, tempmatrix);
                        npower = npower / 2;
                    }
                }
            } else if self.form == TrsfForm::PntMirror || self.form == TrsfForm::Ax1Mirror {
                // OCCT L505-514.
                if n % 2 == 0 {
                    self.form = TrsfForm::Identity;
                    self.scale = 1.0;
                    self.matrix = [[1.0, 0.0], [0.0, 1.0]];
                    self.loc = DVec2::new(0.0, 0.0);
                }
            } else {
                // OCCT L515-545.
                self.form = TrsfForm::CompoundTrsf;
                let mut npower = n;
                if npower < 0 {
                    npower = -npower;
                }
                npower -= 1;
                // matrix.SetDiagonal(scale * matrix.Value(1, 1),
                //                    scale * matrix.Value(2, 2)).
                let d1 = self.scale * self.matrix[0][0];
                let d2 = self.scale * self.matrix[1][1];
                self.matrix[0][0] = d1;
                self.matrix[1][1] = d2;
                let mut temploc = self.loc;
                let mut tempscale = self.scale;
                let mut tempmatrix = self.matrix;
                loop {
                    if npower % 2 == 1 {
                        self.loc += mat2_mul_xy(self.matrix, temploc) * self.scale;
                        self.scale = self.scale * tempscale;
                        self.matrix = mat2_mul(self.matrix, tempmatrix);
                    }
                    if npower == 1 {
                        break;
                    }
                    tempscale = tempscale * tempscale;
                    temploc += mat2_mul_xy(tempmatrix, temploc) * tempscale;
                    tempmatrix = mat2_mul(tempmatrix, tempmatrix);
                    npower = npower / 2;
                }
            }
        }
    }

    /// OCCT gp_Trsf2d::PreMultiply(const gp_Trsf2d& T) (gp_Trsf2d.cxx
    /// L550-671) — `<me> = T * <me>`, so T applies last.  Kept in OCCT's
    /// in-place form (gp_Trsf2d.hxx L163); gp_Trsf2d has no `PreMultiplied`
    /// helper, unlike `Multiplied` which `Trsf2d::multiplied` mirrors.
    pub fn pre_multiply(&mut self, t: &Trsf2d) {
        if t.form == TrsfForm::Identity {
            // OCCT L552-554.
        } else if self.form == TrsfForm::Identity {
            // OCCT L555-561.
            self.form = t.form;
            self.scale = t.scale;
            self.loc = t.loc;
            self.matrix = t.matrix;
        } else if self.form == TrsfForm::Rotation && t.form == TrsfForm::Rotation {
            // OCCT L562-567.
            self.loc = mat2_mul_xy(t.matrix, self.loc); // loc.Multiply(T.matrix).
            self.loc += t.loc;
            self.matrix = mat2_mul(t.matrix, self.matrix); // matrix.PreMultiply(T.matrix).
        } else if self.form == TrsfForm::Translation && t.form == TrsfForm::Translation {
            // OCCT L568-571.
            self.loc += t.loc;
        } else if self.form == TrsfForm::Scale && t.form == TrsfForm::Scale {
            // OCCT L572-577.
            self.loc *= t.scale;
            self.loc += t.loc;
            self.scale = self.scale * t.scale;
        } else if self.form == TrsfForm::PntMirror && t.form == TrsfForm::PntMirror {
            // OCCT L578-584.
            self.scale = 1.0;
            self.form = TrsfForm::Translation;
            self.loc = -self.loc; // loc.Reverse().
            self.loc += t.loc;
        } else if self.form == TrsfForm::Ax1Mirror && t.form == TrsfForm::Ax1Mirror {
            // OCCT L585-593.
            self.form = TrsfForm::Rotation;
            self.loc = mat2_mul_xy(t.matrix, self.loc); // loc.Multiply(T.matrix).
            self.loc *= t.scale; // loc.Multiply(T.scale).
            self.scale = self.scale * t.scale;
            self.loc += t.loc;
            self.matrix = mat2_mul(t.matrix, self.matrix); // matrix.PreMultiply(T.matrix).
        } else if (self.form == TrsfForm::CompoundTrsf
            || self.form == TrsfForm::Rotation
            || self.form == TrsfForm::Ax1Mirror)
            && t.form == TrsfForm::Translation
        {
            // OCCT L594-598.
            self.loc += t.loc;
        } else if (self.form == TrsfForm::Scale || self.form == TrsfForm::PntMirror)
            && t.form == TrsfForm::Translation
        {
            // OCCT L599-602.
            self.loc += t.loc;
        } else if self.form == TrsfForm::Translation
            && (t.form == TrsfForm::CompoundTrsf
                || t.form == TrsfForm::Rotation
                || t.form == TrsfForm::Ax1Mirror)
        {
            // OCCT L603-619.
            self.form = TrsfForm::CompoundTrsf;
            self.matrix = t.matrix;
            if t.scale == 1.0 {
                self.loc = mat2_mul_xy(t.matrix, self.loc); // loc.Multiply(T.matrix).
            } else {
                self.scale = t.scale;
                self.loc = mat2_mul_xy(self.matrix, self.loc); // loc.Multiply(matrix).
                self.loc *= self.scale; // loc.Multiply(scale).
            }
            self.loc += t.loc;
        } else if (t.form == TrsfForm::Scale || t.form == TrsfForm::PntMirror)
            && self.form == TrsfForm::Translation
        {
            // OCCT L620-626.
            self.loc *= t.scale;
            self.loc += t.loc;
            self.scale = t.scale;
            self.form = t.form;
        } else if (self.form == TrsfForm::PntMirror || self.form == TrsfForm::Scale)
            && (t.form == TrsfForm::PntMirror || t.form == TrsfForm::Scale)
        {
            // OCCT L627-634.
            self.form = TrsfForm::CompoundTrsf;
            self.loc *= t.scale;
            self.loc += t.loc;
            self.scale = self.scale * t.scale;
        } else if (self.form == TrsfForm::CompoundTrsf
            || self.form == TrsfForm::Rotation
            || self.form == TrsfForm::Ax1Mirror)
            && (t.form == TrsfForm::Scale || t.form == TrsfForm::PntMirror)
        {
            // OCCT L635-642.
            self.form = TrsfForm::CompoundTrsf;
            self.loc *= t.scale;
            self.loc += t.loc;
            self.scale = self.scale * t.scale;
        } else if (t.form == TrsfForm::CompoundTrsf
            || t.form == TrsfForm::Rotation
            || t.form == TrsfForm::Ax1Mirror)
            && (self.form == TrsfForm::Scale || self.form == TrsfForm::PntMirror)
        {
            // OCCT L643-659.
            self.form = TrsfForm::CompoundTrsf;
            self.matrix = t.matrix;
            if t.scale == 1.0 {
                self.loc = mat2_mul_xy(t.matrix, self.loc); // loc.Multiply(T.matrix).
            } else {
                self.loc = mat2_mul_xy(self.matrix, self.loc); // loc.Multiply(matrix).
                self.loc *= t.scale; // loc.Multiply(T.scale).
                self.scale = t.scale * self.scale;
            }
            self.loc += t.loc;
        } else {
            // OCCT L660-671.
            self.form = TrsfForm::CompoundTrsf;
            self.loc = mat2_mul_xy(t.matrix, self.loc); // loc.Multiply(T.matrix).
            if t.scale != 1.0 {
                self.loc *= t.scale; // loc.Multiply(T.scale).
                self.scale = self.scale * t.scale;
            }
            self.loc += t.loc;
            self.matrix = mat2_mul(t.matrix, self.matrix); // matrix.PreMultiply(T.matrix).
        }
    }

    /// OCCT gp_Trsf2d::SetValues(a11, a12, a13, a21, a22, a23) (gp_Trsf2d.cxx
    /// L676-710): sets the coefficients of the transformation, splitting the
    /// 2x2 part into a uniform scale and an orthogonal matrix.
    pub fn set_values(&mut self, a11: f64, a12: f64, a13: f64, a21: f64, a22: f64, a23: f64) {
        let col1 = DVec2::new(a11, a21);
        let col2 = DVec2::new(a12, a22);
        let col3 = DVec2::new(a13, a23);
        // compute the determinant
        // gp_Mat2d M(col1, col2) — the constructor takes the two columns.
        let m = [[col1.x, col2.x], [col1.y, col2.y]];
        let mut s = m[0][0] * m[1][1] - m[1][0] * m[0][1]; // M.Determinant().
        if s.abs() < GP_RESOLUTION {
            panic!("gp_Trsf2d::SetValues, null determinant");
        }

        if s > 0.0 {
            s = s.sqrt();
        } else {
            s = (-s).sqrt();
        }

        // M.Divide(s).
        let m = [[m[0][0] / s, m[0][1] / s], [m[1][0] / s, m[1][1] / s]];

        self.scale = s;
        self.form = TrsfForm::CompoundTrsf;

        self.matrix = m;
        self.orthogonalize();

        self.loc = col3;
    }

    /// OCCT gp_Trsf2d::Orthogonalize() (gp_Trsf2d.cxx L721-748).
    ///
    /// ATTENTION: orthogonalization is not an equivalent transformation,
    /// therefore a transformation built from the source matrix and one built
    /// from the orthogonalized matrix can give different results for one
    /// shape; the source matrix must already be close to orthogonal.
    fn orthogonalize(&mut self) {
        let mut atm = self.matrix;

        let mut a_v1 = DVec2::new(atm[0][0], atm[1][0]); // aTM.Column(1).
        let mut a_v2 = DVec2::new(atm[0][1], atm[1][1]); // aTM.Column(2).

        a_v1 = gp_xy_normalize(a_v1);

        // aV2 -= aV1 * (aV2.Dot(aV1)).
        a_v2 -= a_v1 * a_v2.dot(a_v1);
        a_v2 = gp_xy_normalize(a_v2);

        // aTM.SetCols(aV1, aV2).
        atm = [[a_v1.x, a_v2.x], [a_v1.y, a_v2.y]];

        let mut a_v1 = DVec2::new(atm[0][0], atm[0][1]); // aTM.Row(1).
        let mut a_v2 = DVec2::new(atm[1][0], atm[1][1]); // aTM.Row(2).

        a_v1 = gp_xy_normalize(a_v1);

        // aV2 -= aV1 * (aV2.Dot(aV1)).
        a_v2 -= a_v1 * a_v2.dot(a_v1);
        a_v2 = gp_xy_normalize(a_v2);

        // aTM.SetRows(aV1, aV2).
        atm = [[a_v1.x, a_v1.y], [a_v2.x, a_v2.y]];

        self.matrix = atm;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An arbitrary gp_CompoundTrsf state: a non-identity HVectorialPart, an
    /// explicit loc and scale.  The fields are written directly instead of
    /// through `set_values` because `SetValues` orthogonalizes the matrix it
    /// is given (gp_Trsf2d.cxx L721-749) and the counterexamples below need a
    /// deliberately non-orthogonal matrix.
    fn compound_trsf(m: [[f64; 2]; 2], loc: DVec2, scale: f64) -> Trsf2d {
        Trsf2d {
            matrix: m,
            loc,
            scale,
            form: TrsfForm::CompoundTrsf,
        }
    }

    /// Counterexample for the general branch of `Trsf2d::multiplied`
    /// (OCCT gp_Trsf2d::Multiply final else, gp_Trsf2d.cxx L365-381):
    /// `Tloc` is `matrix * T.loc`, the left scale folds into it afterwards,
    /// and the vectorial parts compose as `matrix * T.matrix`.
    #[test]
    fn multiplied_general_branch_puts_the_matrix_on_the_translation() {
        let m = [[1.0, 2.0], [3.0, 5.0]];
        let t_m = [[2.0, 0.0], [0.0, 3.0]];
        let self_t = compound_trsf(m, DVec2::ZERO, 2.0);
        let t = compound_trsf(t_m, DVec2::new(1.0, 1.0), 1.0);

        let r = self_t.multiplied(&t);

        // OCCT L368-380: Tloc = matrix * T.loc = (1 + 2, 3 + 5) = (3, 8);
        // scale = 2 != 1 -> Tloc *= 2 -> (6, 16) and scale = 2 * 1 = 2;
        // loc = (0, 0) + (6, 16); matrix = matrix * T.matrix.
        assert_eq!(r.form(), TrsfForm::CompoundTrsf);
        assert_eq!(r.scale, 2.0);
        assert_eq!(r.loc, DVec2::new(6.0, 16.0));
        // Both transformations applied in order (t first, then self):
        // t(1, 0) = T.matrix * (1, 0) * 1 + (1, 1) = (2, 0) + (1, 1) = (3, 1);
        // self: matrix * (3, 1) * 2 = (3 + 2, 9 + 5) * 2 = (10, 28).
        let got = r.transformed(DVec2::new(1.0, 0.0));
        assert!(
            (got.x - 10.0).abs() < 1e-12 && (got.y - 28.0).abs() < 1e-12,
            "got={got:?} expected=(10, 28)"
        );
    }

    /// The composition performed by
    /// ShapeExtend_CompositeSurface::GlobalToLocalTransformation
    /// (ShapeExtend_CompositeSurface.cxx L382-391): `Scale * Shift` with a
    /// non-unit scale, a non-zero scale center and a non-zero translation.
    /// OCCT gp_Trsf2d::Multiply (Scale || PntMirror) && Translation branch
    /// (gp_Trsf2d.cxx L309-314) folds the scale into the copied translation.
    #[test]
    fn scale_translation_composition_folds_the_scale_into_the_translation() {
        let mut scale_t = Trsf2d::default();
        // loc = P * (1 - S) = (1, 1) * (1 - 3) = (-2, -2).
        scale_t.set_scale(DVec2::new(1.0, 1.0), 3.0);
        let mut shift_t = Trsf2d::default();
        shift_t.set_translation(DVec2::new(4.0, 4.0));

        let trsf = scale_t.multiplied(&shift_t);

        assert_eq!(trsf.form(), TrsfForm::Scale);
        assert_eq!(trsf.scale, 3.0);
        // OCCT L311-313: Tloc = T.loc * scale = (12, 12); loc = (-2, -2) +
        // (12, 12) = (10, 10).
        assert_eq!(trsf.loc, DVec2::new(10.0, 10.0));
        // Shift first, then the scale about (1, 1): (p + (4, 4) - (1, 1)) * 3
        // + (1, 1) = 3 * p + (10, 10).
        let got = trsf.transformed(DVec2::new(1.0, 1.0));
        assert!(
            (got.x - 13.0).abs() < 1e-12 && (got.y - 13.0).abs() < 1e-12,
            "got={got:?} expected=(13, 13)"
        );
        assert_eq!(trsf.transformed(DVec2::ZERO), DVec2::new(10.0, 10.0));
    }

    /// gp_Trsf2d::SetScale (gp_Trsf2d.hxx L270-277) raises for no scale
    /// value, so S = 0 must yield the degenerate gp_Scale transformation
    /// instead of aborting.
    #[test]
    fn set_scale_accepts_a_zero_scale_like_the_2d_class() {
        let mut t = Trsf2d::default();
        t.set_scale(DVec2::new(2.0, 3.0), 0.0);
        assert_eq!(t.form(), TrsfForm::Scale);
        assert_eq!(t.scale, 0.0);
        // loc = P * (1 - S) = (2, 3); `matrix * P * scale` vanishes.
        assert_eq!(t.transformed(DVec2::new(1.0, 1.0)), DVec2::new(2.0, 3.0));
    }

    /// gp_Trsf2d::SetRotation(P, Ang) (gp_Trsf2d.hxx L246-257) rotates about
    /// the centre P: P' = P_centre + R(Ang) * (P - P_centre).  The centre and
    /// the angle are both off the degeneracies (0, 0) / 0 / pi so that a
    /// transposed `loc.Multiply(matrix)` or a sign slip cannot pass.
    #[test]
    fn set_rotation_turns_about_the_centre() {
        let centre = DVec2::new(1.5, -0.5);
        let ang = 0.7_f64;
        let mut t = Trsf2d::default();
        t.set_rotation(centre, ang);
        assert_eq!(t.form(), TrsfForm::Rotation);
        assert_eq!(t.scale, 1.0);
        let (s, c) = ang.sin_cos();
        for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
            let d = p - centre;
            let expected = centre + DVec2::new(c * d.x - s * d.y, s * d.x + c * d.y);
            let got = t.transformed(p);
            assert!(
                (got - expected).length() < 1e-12,
                "p={p:?} got={got:?} expected={expected:?}"
            );
        }
        // The centre itself is invariant.
        assert!((t.transformed(centre) - centre).length() < 1e-12);
    }

    /// gp_Trsf2d::SetMirror(const gp_Pnt2d&) (gp_Trsf2d.hxx L259-266) is the
    /// central symmetry: P' = 2 * centre - P.
    #[test]
    fn set_mirror_pnt_is_central_symmetry() {
        let centre = DVec2::new(-1.0, 4.0);
        let mut t = Trsf2d::default();
        t.set_mirror_pnt(centre);
        assert_eq!(t.form(), TrsfForm::PntMirror);
        assert_eq!(t.scale, -1.0);
        for &p in &[DVec2::new(1.0, 2.0), DVec2::new(0.0, 0.0)] {
            let got = t.transformed(p);
            let expected = centre * 2.0 - p;
            assert!(
                (got - expected).length() < 1e-12,
                "p={p:?} got={got:?} expected={expected:?}"
            );
        }
    }

    /// gp_Trsf2d::SetMirror(const gp_Ax2d&) (gp_Trsf2d.cxx L31-46) reflects
    /// across the LINE through `origin` with direction `dir`.  The group law
    /// applied here is the Householder reflection `2*v*v^T - I` about the
    /// line's own origin, derived independently of the OCCT matrix/offset
    /// expressions.
    #[test]
    fn set_mirror_ax2d_reflects_across_the_line() {
        let origin = DVec2::new(1.0, 1.0);
        let dir = DVec2::new(0.6, 0.8); // unit
        let mut t = Trsf2d::default();
        t.set_mirror_ax2d(origin, dir);
        assert_eq!(t.form(), TrsfForm::Ax1Mirror);
        assert_eq!(t.scale, -1.0);
        for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25), origin] {
            let d = p - origin;
            let refl = dir * (2.0 * d.dot(dir)) - d;
            let expected = origin + refl;
            let got = t.transformed(p);
            assert!(
                (got - expected).length() < 1e-12,
                "p={p:?} got={got:?} expected={expected:?}"
            );
        }
    }

    /// The composition law of gp_Trsf2d::Multiply: `A.multiplied(B)` applied to
    /// a point equals `A(B(P))`.  Every branch of the OCCT switch computes its
    /// own loc/scale/matrix combination, so this single assertion covers the
    /// algebra of all of them without restating the formulas.
    #[test]
    fn multiplied_composes_like_apply_after_apply() {
        let mut rot = Trsf2d::default();
        rot.set_rotation(DVec2::new(1.5, -0.5), 0.7);
        let mut trans = Trsf2d::default();
        trans.set_translation(DVec2::new(2.0, 3.0));
        let mut scale = Trsf2d::default();
        scale.set_scale(DVec2::new(0.5, 0.5), 2.5);
        let mut pnt_mirror = Trsf2d::default();
        pnt_mirror.set_mirror_pnt(DVec2::new(-1.0, 4.0));
        let mut ax_mirror = Trsf2d::default();
        ax_mirror.set_mirror_ax2d(DVec2::new(1.0, 1.0), DVec2::new(0.6, 0.8));

        let all = [rot, trans, scale, pnt_mirror, ax_mirror];
        let names = ["rotation", "translation", "scale", "pntMirror", "ax1Mirror"];
        for (i, a) in all.iter().enumerate() {
            for (j, b) in all.iter().enumerate() {
                let c = a.multiplied(b);
                for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
                    let got = c.transformed(p);
                    let expected = a.transformed(b.transformed(p));
                    assert!(
                        (got - expected).length() < 1e-9,
                        "{}({}) p={p:?}: got={got:?} expected={expected:?}",
                        names[i],
                        names[j]
                    );
                }
            }
        }
    }

    /// The form tags must follow gp_Trsf2d::Multiply's branches exactly — the
    /// earlier partial body collapsed the ones it did not implement to
    /// CompoundTrsf.
    #[test]
    fn multiplied_follows_the_occt_form_tags() {
        let mut rot = Trsf2d::default();
        rot.set_rotation(DVec2::new(1.5, -0.5), 0.7);
        let mut rot2 = Trsf2d::default();
        rot2.set_rotation(DVec2::ZERO, 0.3);
        // OCCT L265-272 keeps gp_Rotation.
        assert_eq!(rot.multiplied(&rot2).form(), TrsfForm::Rotation);

        let mut t1 = Trsf2d::default();
        t1.set_translation(DVec2::new(1.0, 0.0));
        let mut t2 = Trsf2d::default();
        t2.set_translation(DVec2::new(0.0, 1.0));
        // OCCT L273-276 keeps gp_Translation.
        assert_eq!(t1.multiplied(&t2).form(), TrsfForm::Translation);

        let mut s = Trsf2d::default();
        s.set_scale(DVec2::ZERO, 2.0);
        // OCCT L277-281 keeps gp_Scale; OCCT L323-329 takes shape = T.shape.
        assert_eq!(s.multiplied(&s).form(), TrsfForm::Scale);
        assert_eq!(t1.multiplied(&s).form(), TrsfForm::Scale);

        let mut m1 = Trsf2d::default();
        m1.set_mirror_pnt(DVec2::new(1.0, 0.0));
        let mut m2 = Trsf2d::default();
        m2.set_mirror_pnt(DVec2::new(0.0, 1.0));
        // OCCT L282-287: two central symmetries are a translation.
        let mm = m1.multiplied(&m2);
        assert_eq!(mm.form(), TrsfForm::Translation);
        assert_eq!(mm.scale, 1.0);

        let mut a1 = Trsf2d::default();
        a1.set_mirror_ax2d(DVec2::ZERO, DVec2::new(1.0, 0.0));
        let mut a2 = Trsf2d::default();
        a2.set_mirror_ax2d(DVec2::ZERO, DVec2::new(0.0, 1.0));
        // OCCT L288-297: two line mirrors are a rotation.
        assert_eq!(a1.multiplied(&a2).form(), TrsfForm::Rotation);
    }

    /// gp_Trsf2d::Value(Row, Col) (gp_Trsf2d.hxx L301-314): Col < 3 is the
    /// scaled matrix entry, Col 3 the location coordinate.  Note that
    /// gp_Trsf2d::SetScale resets the matrix to the identity (hxx L270-277),
    /// so the rotation and the scale are read from two separate states.
    #[test]
    fn value_folds_in_the_scale_and_reads_the_location() {
        let centre = DVec2::new(1.5, -0.5);
        let ang = 0.7_f64;
        let (s, c) = ang.sin_cos();
        let mut rot = Trsf2d::default();
        rot.set_rotation(centre, ang);
        assert!((rot.value(1, 1) - c).abs() < 1e-12);
        assert!((rot.value(1, 2) - -s).abs() < 1e-12);
        assert!((rot.value(2, 1) - s).abs() < 1e-12);
        assert!((rot.value(2, 2) - c).abs() < 1e-12);
        // Col 3 is the location coordinate, independent of the scale.
        assert_eq!(rot.value(1, 3), rot.loc.x);
        assert_eq!(rot.value(2, 3), rot.loc.y);

        let mut sc = Trsf2d::default();
        sc.set_scale(centre, 2.0);
        assert_eq!(sc.value(1, 1), 2.0);
        assert_eq!(sc.value(2, 2), 2.0);
        assert_eq!(sc.value(1, 2), 0.0);
        assert_eq!(sc.value(1, 3), sc.loc.x);
    }

    /// Assert that two points coincide within 1e-12 (tight enough to catch a
    /// transposed matrix, a missing scale factor or a wrong branch).
    fn assert_near(got: DVec2, expected: DVec2, ctx: &str) {
        assert!(
            (got - expected).length() < 1e-12,
            "{ctx}: got={got:?} expected={expected:?}"
        );
    }

    /// gp_Trsf2d::SetTransformation(const gp_Ax2d& FromA1, const gp_Ax2d&
    /// ToA2) (gp_Trsf2d.cxx L48-70) is a frame change: it rewrites the
    /// coordinates of one and the same world point from the FromA1 frame into
    /// the ToA2 frame.  The property asserted here — the world position
    /// rebuilt from the FromA1 coordinates equals the world position rebuilt
    /// from the ToA2 coordinates — is the definition of a frame change and
    /// does not restate the matrix/offset expressions of the source.
    #[test]
    fn set_transformation_from_to_changes_the_frame() {
        let o1 = DVec2::new(2.0, 1.0);
        let d1 = DVec2::new(0.6, 0.8); // unit
        let o2 = DVec2::new(-1.0, 3.0);
        let d2 = DVec2::new(0.0, 1.0); // unit, not parallel to d1
        let p1 = DVec2::new(-d1.y, d1.x);
        let p2 = DVec2::new(-d2.y, d2.x);

        let mut t = Trsf2d::default();
        t.set_transformation_from_to(o1, d1, o2, d2);
        assert_eq!(t.form(), TrsfForm::CompoundTrsf);
        assert_eq!(t.scale, 1.0);

        // world(A1, x) == world(A2, t(x)).
        let world_a1 = |x: DVec2| d1 * x.x + p1 * x.y + o1;
        let world_a2 = |x: DVec2| d2 * x.x + p2 * x.y + o2;
        for local in [DVec2::ZERO, DVec2::new(1.0, 0.0), DVec2::new(-2.0, 3.0)] {
            assert_near(world_a2(t.transformed(local)), world_a1(local), "frame change");
        }

        // FromA1 == ToA2 is the identity frame change.
        let mut id = Trsf2d::default();
        id.set_transformation_from_to(o1, d1, o1, d1);
        for x in [DVec2::ZERO, DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
            assert_near(id.transformed(x), x, "identity frame change");
        }
    }

    /// gp_Trsf2d::SetTransformation(const gp_Ax2d& A) (gp_Trsf2d.cxx L72-84)
    /// maps the basic frame onto A, so the image of a point is the projection
    /// of its offset from the origin of A onto the two axes of A: for a unit
    /// direction d, `T(x) == (dot(x - o, d), dot(x - o, perp(d)))`.
    #[test]
    fn set_transformation_ax2d_rewrites_in_the_axis_frame() {
        let o = DVec2::new(1.5, -0.5);
        let d = DVec2::new(0.6, 0.8); // unit
        let q = DVec2::new(-d.y, d.x);
        let mut t = Trsf2d::default();
        t.set_transformation_ax2d(o, d);
        assert_eq!(t.form(), TrsfForm::CompoundTrsf);
        assert_eq!(t.scale, 1.0);
        for &x in &[DVec2::ZERO, DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25), o] {
            let v = x - o;
            let expected = DVec2::new(v.dot(d), v.dot(q));
            assert_near(t.transformed(x), expected, "axis frame change");
        }
        // The origin of A is the origin of the new frame.
        assert_near(t.transformed(o), DVec2::ZERO, "axis origin");
    }

    /// gp_Trsf2d::SetTranslationPart(V) (gp_Trsf2d.cxx L86-117) replaces the
    /// translation part.  The asserted content is the geometry plus the tag
    /// that geometry forces: a gp_Identity with a non-zero V is a genuine
    /// translation by V, and a gp_Translation whose V is set back to zero is
    /// the identity.
    #[test]
    fn set_translation_part_retags_identity_and_translation() {
        let v = DVec2::new(1.5, -2.0);
        let mut t = Trsf2d::default();
        t.set_translation_part(v);
        assert_eq!(t.form(), TrsfForm::Translation);
        for &p in &[DVec2::ZERO, DVec2::new(1.0, 2.0)] {
            assert_near(t.transformed(p), p + v, "identity -> translation");
        }

        t.set_translation_part(DVec2::ZERO);
        assert_eq!(t.form(), TrsfForm::Identity);
        for &p in &[DVec2::ZERO, DVec2::new(1.0, 2.0)] {
            assert_near(t.transformed(p), p, "translation -> identity");
        }

        // From a rotation the vectorial part is untouched and the offset
        // becomes V.  The tag does not stay gp_Rotation: OCCT falls through to
        // the catch-all `shape = gp_CompoundTrsf` (cxx L113-116), which is
        // still correct for the caller, whose only test is `!= gp_Identity`.
        // `T(p) - T(0)` is the linear part of the map.
        let mut before = Trsf2d::default();
        before.set_rotation(DVec2::new(1.5, -0.5), 0.7);
        let mut r = before;
        r.set_translation_part(v);
        assert_eq!(r.form(), TrsfForm::CompoundTrsf);
        assert_eq!(r.loc, v);
        for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
            let linear_part = before.transformed(p) - before.transformed(DVec2::ZERO);
            assert_near(r.transformed(p), linear_part + v, "rotation keeps its linear part");
        }
    }

    /// gp_Trsf2d::SetScaleFactor(S) (gp_Trsf2d.cxx L120-196).  The asserted
    /// transitions are the ones whose resulting tag is forced by the map
    /// alone, so the expectations are not a restatement of the branch table:
    /// - a central symmetry whose scale becomes 1 is the translation by
    ///   2 * centre, and must be tagged gp_Translation (cxx L143-146);
    /// - a scale by 1 about the origin is the identity (cxx L139-142);
    /// - the identity with S = -1 is the central symmetry about the origin
    ///   (cxx L172-175).
    #[test]
    fn set_scale_factor_retags_to_the_form_the_map_has() {
        let centre = DVec2::new(-1.0, 4.0);

        let mut pm = Trsf2d::default();
        pm.set_mirror_pnt(centre); // p -> 2c - p
        pm.set_scale_factor(1.0); // the map becomes p -> p + 2c
        assert_eq!(pm.form(), TrsfForm::Translation);
        assert_eq!(pm.scale, 1.0);
        assert_near(pm.transformed(DVec2::ZERO), centre * 2.0, "mirror point scaled by 1");

        let mut sc = Trsf2d::default();
        sc.set_scale(DVec2::ZERO, 3.0); // p -> 3p, loc == 0
        sc.set_scale_factor(1.0); // the map becomes p -> p
        assert_eq!(sc.form(), TrsfForm::Identity);
        assert_eq!(sc.scale, 1.0);
        assert_near(sc.transformed(DVec2::new(1.0, 2.0)), DVec2::new(1.0, 2.0), "scale by 1");

        let mut id = Trsf2d::default();
        id.set_scale_factor(-1.0); // the map becomes p -> -p
        assert_eq!(id.form(), TrsfForm::PntMirror);
        assert_eq!(id.scale, -1.0);
        assert_near(id.transformed(DVec2::new(1.0, 2.0)), DVec2::new(-1.0, -2.0), "scale by -1");

        // Any other S leaves the matrix and the offset alone and multiplies
        // the linear part: T'(p) = S * (T(p) - T(0)) + T(0).
        let mut before = Trsf2d::default();
        before.set_rotation(DVec2::new(1.5, -0.5), 0.7);
        let mut t = before;
        t.set_scale_factor(2.5);
        assert_eq!(t.form(), TrsfForm::CompoundTrsf);
        assert_eq!(t.scale, 2.5);
        for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
            let expected = before.transformed(p) * 2.5 + before.loc * (1.0 - 2.5);
            assert_near(t.transformed(p), expected, "scale factor only rescales the linear part");
        }
    }

    /// gp_Trsf2d::VectorialPart() (gp_Trsf2d.cxx L198-214) is, by its own
    /// documentation (gp_Trsf2d.hxx L114-116), the 2x2 matrix which applied to
    /// a point gives the linear part of the map:
    /// `VectorialPart * p == T(p) - T(0)`.  That identity is independent of
    /// the diagonal/multiply split of the source.
    #[test]
    fn vectorial_part_is_the_linear_part_of_the_map() {
        let mut rot = Trsf2d::default();
        rot.set_rotation(DVec2::new(1.5, -0.5), 0.7);
        let mut scale = Trsf2d::default();
        scale.set_scale(DVec2::new(1.0, 1.0), 3.0);
        // a gp_CompoundTrsf with scale != 1: the M.Multiply(scale) branch.
        let comp = compound_trsf([[0.6, -0.8], [0.8, 0.6]], DVec2::new(2.0, -1.0), 2.5);

        for (name, t) in [("rotation", rot), ("scale", scale), ("compound", comp)] {
            let v = t.vectorial_part();
            for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
                let got = DVec2::new(v[0][0] * p.x + v[0][1] * p.y, v[1][0] * p.x + v[1][1] * p.y);
                let expected = t.transformed(p) - t.transformed(DVec2::ZERO);
                assert_near(got, expected, name);
            }
        }

        // The scale == 1 short-circuit returns the raw HVectorialPart.
        let raw = compound_trsf([[1.0, 2.0], [3.0, 5.0]], DVec2::ZERO, 1.0);
        assert_eq!(raw.vectorial_part(), [[1.0, 2.0], [3.0, 5.0]]);
    }

    /// gp_Trsf2d::RotationPart() (gp_Trsf2d.cxx L216-219) is documented as the
    /// operation opposite to SetRotation, so it must return the angle that
    /// built the state.  A swapped atan2 argument would return pi/2 - ang.
    #[test]
    fn rotation_part_inverts_set_rotation() {
        for ang in [-2.5_f64, -0.4, 0.0, 0.7, 2.5] {
            let mut t = Trsf2d::default();
            t.set_rotation(DVec2::new(1.5, -0.5), ang);
            assert!(
                (t.rotation_part() - ang).abs() < 1e-12,
                "ang={ang} got={}",
                t.rotation_part()
            );
        }
    }

    /// gp_Trsf2d::Invert() (gp_Trsf2d.cxx L221-251).  The defining property of
    /// an inverse — `t_inv(t(p)) == p` — is independent of the branch
    /// structure; the gp_Scale and gp_Rotation states additionally pin their
    /// branch geometrically (a rotation's inverse is the rotation by -ang
    /// about the same centre, a scale's inverse is the scale by 1/S about the
    /// same centre).
    #[test]
    fn invert_undoes_the_transformation() {
        let centre = DVec2::new(1.5, -0.5);
        let ang = 0.7_f64;
        let probes = [DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25), centre];

        // gp_Rotation, with a centre off the origin (the general else branch).
        let mut rot = Trsf2d::default();
        rot.set_rotation(centre, ang);
        let mut rot_inv = rot;
        rot_inv.invert();
        assert_eq!(rot_inv.form(), TrsfForm::Rotation);
        let (s, c) = (-ang).sin_cos();
        for &p in &probes {
            let d = p - centre;
            let expected = centre + DVec2::new(c * d.x - s * d.y, s * d.x + c * d.y);
            assert_near(rot_inv.transformed(p), expected, "inverse rotation");
            assert_near(rot_inv.transformed(rot.transformed(p)), p, "rotation round trip");
            assert_near(rot.transformed(rot_inv.transformed(p)), p, "rotation round trip");
        }

        // gp_Scale (its own branch).
        let mut sc = Trsf2d::default();
        sc.set_scale(centre, 3.0);
        let mut sc_inv = sc;
        sc_inv.invert();
        assert_eq!(sc_inv.form(), TrsfForm::Scale);
        assert_eq!(sc_inv.scale, 1.0 / 3.0);
        for &p in &probes {
            assert_near(sc_inv.transformed(p), centre + (p - centre) / 3.0, "inverse scale");
            assert_near(sc_inv.transformed(sc.transformed(p)), p, "scale round trip");
        }

        // gp_Translation / gp_PntMirror branch: the offset is just reversed.
        let mut tr = Trsf2d::default();
        tr.set_translation(DVec2::new(2.0, 3.0));
        let mut tr_inv = tr;
        tr_inv.invert();
        assert_eq!(tr_inv.form(), TrsfForm::Translation);
        for &p in &probes {
            assert_near(tr_inv.transformed(p), p - DVec2::new(2.0, 3.0), "inverse translation");
        }
        // gp_PntMirror: OCCT only reverses the offset here (cxx L231-234) and
        // keeps the tag, so the central symmetry about c becomes the central
        // symmetry about -c.  That is NOT the inverse of the map — the 3D
        // gp_Trsf::Invert does exactly the same (gp_Trsf.cxx L406-409) — so
        // this assertion pins OCCT's behaviour instead of a round trip.
        let mut pm = Trsf2d::default();
        pm.set_mirror_pnt(centre); // p -> 2c - p
        let mut pm_inv = pm;
        pm_inv.invert();
        assert_eq!(pm_inv.form(), TrsfForm::PntMirror);
        for &p in &probes {
            assert_near(pm_inv.transformed(p), -p - centre * 2.0, "reversed point mirror");
        }

        // An arbitrary gp_CompoundTrsf whose matrix is orthogonal (det = -1):
        // R-1 == R transposed, so the round trip closes.
        let mut comp = compound_trsf([[0.6, -0.8], [0.8, 0.6]], DVec2::new(0.5, -1.0), 2.0);
        let comp_fwd = comp;
        comp.invert();
        assert_eq!(comp.form(), TrsfForm::CompoundTrsf);
        for &p in &probes {
            assert_near(comp.transformed(comp_fwd.transformed(p)), p, "compound round trip");
        }
    }

    /// gp_Trsf2d::Invert() raises Standard_ConstructionError when the scale is
    /// not above gp::Resolution() (gp_Trsf2d.cxx L237-238, L244-245); the
    /// carrier panics on the same condition.
    #[test]
    #[should_panic(expected = "zero scale")]
    fn invert_panics_on_a_zero_scale() {
        let mut t = Trsf2d::default();
        t.set_scale(DVec2::new(1.0, 1.0), 0.0);
        t.invert();
    }

    /// gp_Trsf2d::PreMultiply(T) (gp_Trsf2d.cxx L550-671) composes
    /// `<me> = T * <me>`, so for a point the result is `T(<me>(p))` — the
    /// definition of pre-multiplication, independent of the branch structure.
    /// The states include two non-commuting rotations with non-zero offsets,
    /// so a wrong factor order (or a transposed matrix) cannot pass.
    #[test]
    fn pre_multiply_applies_its_argument_last() {
        let mut rot = Trsf2d::default();
        rot.set_rotation(DVec2::new(1.5, -0.5), 0.7);
        let mut rot2 = Trsf2d::default();
        rot2.set_rotation(DVec2::new(-0.5, 2.0), 0.3);
        let mut trans = Trsf2d::default();
        trans.set_translation(DVec2::new(2.0, 3.0));
        let mut scale = Trsf2d::default();
        scale.set_scale(DVec2::new(0.5, 0.5), 2.5);
        let mut pnt_mirror = Trsf2d::default();
        pnt_mirror.set_mirror_pnt(DVec2::new(-1.0, 4.0));
        let mut ax_mirror = Trsf2d::default();
        ax_mirror.set_mirror_ax2d(DVec2::new(1.0, 1.0), DVec2::new(0.6, 0.8));
        let identity = Trsf2d::default();

        let all = [rot, rot2, trans, scale, pnt_mirror, ax_mirror, identity];
        let names = ["rot", "rot2", "trans", "scale", "pntMirror", "ax1Mirror", "identity"];
        for (i, a) in all.iter().enumerate() {
            for (j, b) in all.iter().enumerate() {
                let mut c = *a;
                c.pre_multiply(b);
                for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
                    let expected = b.transformed(a.transformed(p));
                    assert_near(c.transformed(p), expected, &format!("{}({})", names[j], names[i]));
                }
                // The caller of the carrier (GlobalToLocalTransformation, cxx
                // L382-391) reads the tag as `form() != Identity`: gp_Identity
                // must mean the identity map.
                if c.form() == TrsfForm::Identity {
                    for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
                        assert_near(c.transformed(p), p, "identity tag");
                    }
                    assert_eq!(c.scale, 1.0, "identity tag");
                    assert_eq!(c.loc, DVec2::ZERO, "identity tag");
                }
            }
        }
    }

    /// gp_Trsf2d::SetValues(a11..a23) (gp_Trsf2d.cxx L676-710) documents
    /// `x' = a11 x + a12 y + a13`, `y' = a21 x + a22 y + a23` and
    /// `Value(i, j) == aij` (gp_Trsf2d.hxx L186-201).  For an input whose 2x2
    /// part is already a uniform scale times an orthogonal matrix — the case
    /// the orthogonalization pass leaves alone — the map must reproduce those
    /// equations exactly.
    #[test]
    fn set_values_reproduces_the_documented_coefficients() {
        let ang = 0.7_f64;
        let (s, c) = ang.sin_cos();
        let k = 2.5; // uniform scale of the 2x2 part
        let (a11, a12, a21, a22) = (k * c, -k * s, k * s, k * c);
        let (a13, a23) = (1.5, -0.25);
        let mut t = Trsf2d::default();
        t.set_values(a11, a12, a13, a21, a22, a23);
        assert_eq!(t.form(), TrsfForm::CompoundTrsf);
        assert_eq!(t.scale, k);

        for (r, col, expected) in [(1, 1, a11), (1, 2, a12), (2, 1, a21), (2, 2, a22)] {
            assert!(
                (t.value(r, col) - expected).abs() < 1e-12,
                "Value({r}, {col}) got={} expected={expected}",
                t.value(r, col)
            );
        }
        assert_eq!(t.value(1, 3), a13);
        assert_eq!(t.value(2, 3), a23);

        for &p in &[DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25)] {
            let expected = DVec2::new(a11 * p.x + a12 * p.y + a13, a21 * p.x + a22 * p.y + a23);
            assert_near(t.transformed(p), expected, "documented map");
        }
    }

    /// `Orthogonalize` (gp_Trsf2d.cxx L721-749) replaces the 2x2 part by an
    /// orthogonal matrix, dividing out the determinant sign through the scale.
    /// The characterizing properties — independent of the source statements —
    /// are that the result is orthogonal and that the first column keeps its
    /// direction.
    #[test]
    fn orthogonalize_makes_the_matrix_orthogonal_keeping_the_first_column() {
        let mut t = Trsf2d::default();
        // Not orthogonal, and the first column is not axis aligned.
        t.set_values(1.0, 1.0, 0.0, 1.0, 0.0, 0.0);
        assert_eq!(t.scale, 1.0);
        let v = t.vectorial_part();
        let (a, b, c, d) = (v[0][0], v[0][1], v[1][0], v[1][1]);
        assert!((a * a + c * c - 1.0).abs() < 1e-12, "column 1 not unit: {v:?}");
        assert!((b * b + d * d - 1.0).abs() < 1e-12, "column 2 not unit: {v:?}");
        assert!((a * b + c * d).abs() < 1e-12, "columns not orthogonal: {v:?}");
        // The first column direction (1, 1) is kept.
        assert!(
            (a - c).abs() < 1e-12 && a > 0.0,
            "first column direction changed: {v:?}"
        );
        // The determinant sign of the source survives (det < 0 here).
        assert!((a * d - b * c + 1.0).abs() < 1e-12, "determinant sign lost: {v:?}");
    }

    /// gp_Trsf2d::SetValues raises Standard_ConstructionError for a null
    /// determinant (gp_Trsf2d.cxx L689-690); the carrier panics.
    #[test]
    #[should_panic(expected = "null determinant")]
    fn set_values_panics_on_a_null_determinant() {
        let mut t = Trsf2d::default();
        t.set_values(1.0, 2.0, 0.0, 2.0, 4.0, 0.0);
    }

    /// gp_Trsf2d::Power(N) (gp_Trsf2d.cxx L384-548) is the N-fold composition,
    /// so the result applied to a point must equal N successive applications
    /// of the original — the definition of a power.  The states cover every
    /// arm of the form dispatch (translation, scale, rotation with and without
    /// an offset, point mirror, negative N).
    #[test]
    fn power_repeats_the_transformation() {
        let centre = DVec2::new(1.5, -0.5);
        let ang = 0.4_f64;
        let probes = [DVec2::new(1.0, 2.0), DVec2::new(-3.0, 0.25), centre];

        // gp_Rotation about a centre off the origin (the loc != 0 arm).
        let mut rot = Trsf2d::default();
        rot.set_rotation(centre, ang);
        let mut p3 = rot;
        p3.power(3);
        assert_eq!(p3.form(), TrsfForm::Rotation);
        let (s, c) = (3.0 * ang).sin_cos();
        for &p in &probes {
            assert_near(
                p3.transformed(p),
                rot.transformed(rot.transformed(rot.transformed(p))),
                "rotation power 3",
            );
            let d = p - centre;
            let expected = centre + DVec2::new(c * d.x - s * d.y, s * d.x + c * d.y);
            assert_near(p3.transformed(p), expected, "rotation power 3 angle");
        }

        // gp_Rotation about the origin (the loc == 0 arm).
        let mut rot0 = Trsf2d::default();
        rot0.set_rotation(DVec2::ZERO, ang);
        let mut p3_0 = rot0;
        p3_0.power(3);
        let (s, c) = (3.0 * ang).sin_cos();
        assert!((p3_0.value(1, 1) - c).abs() < 1e-12, "origin rotation power 3 cos");
        assert!((p3_0.value(2, 1) - s).abs() < 1e-12, "origin rotation power 3 sin");

        // gp_Translation.
        let mut tr = Trsf2d::default();
        tr.set_translation(DVec2::new(1.0, -2.0));
        let mut tr4 = tr;
        tr4.power(4);
        assert_eq!(tr4.form(), TrsfForm::Translation);
        assert_near(tr4.transformed(DVec2::ZERO), DVec2::new(4.0, -8.0), "translation power 4");
        assert_near(
            tr4.transformed(DVec2::new(3.0, 5.0)),
            DVec2::new(7.0, -3.0),
            "translation power 4",
        );

        // gp_Scale: S^3 is the scale about the same centre by S^3.
        let mut sc = Trsf2d::default();
        sc.set_scale(centre, 1.5);
        let mut sc3 = sc;
        sc3.power(3);
        assert_eq!(sc3.form(), TrsfForm::Scale);
        for &p in &probes {
            assert_near(sc3.transformed(p), centre + (p - centre) * 3.375, "scale power 3");
        }

        // gp_PntMirror is an involution: an even power is the identity, an odd
        // one is the mirror itself.
        let mut pm = Trsf2d::default();
        pm.set_mirror_pnt(centre);
        let mut pm2 = pm;
        pm2.power(2);
        assert_eq!(pm2.form(), TrsfForm::Identity);
        let mut pm3 = pm;
        pm3.power(3);
        assert_eq!(pm3.form(), TrsfForm::PntMirror);
        for &p in &probes {
            assert_near(pm2.transformed(p), p, "point mirror power 2");
            assert_near(pm3.transformed(p), pm.transformed(p), "point mirror power 3");
        }

        // A negative power is the inverse applied |N| times.
        let mut rotm = rot;
        rotm.power(-2);
        let mut expected = Trsf2d::default();
        expected.set_rotation(centre, -2.0 * ang);
        for &p in &probes {
            assert_near(rotm.transformed(p), expected.transformed(p), "rotation power -2");
        }

        // N == 0 and N == 1 are the documented shortcuts (gp_Trsf2d.hxx
        // L167-173), and N == -1 is Invert().
        let mut p0 = rot;
        p0.power(0);
        assert_eq!(p0.form(), TrsfForm::Identity);
        assert_near(p0.transformed(DVec2::new(1.0, 2.0)), DVec2::new(1.0, 2.0), "power 0");
        let mut p1 = rot;
        p1.power(1);
        assert_near(p1.transformed(DVec2::new(1.0, 2.0)), rot.transformed(DVec2::new(1.0, 2.0)), "power 1");
        let mut pm1 = rot;
        pm1.power(-1);
        let mut inv = rot;
        inv.invert();
        assert_near(pm1.transformed(DVec2::new(1.0, 2.0)), inv.transformed(DVec2::new(1.0, 2.0)), "power -1");
    }
}
