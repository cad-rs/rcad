//! OCCT ShapeExtend_CompositeSurface (TKShHealing): `.hxx` L67-307 and
//! `.cxx` L32-758 — a composite surface represented by a grid of surface
//! patches with a linear local-to-global parametrisation.
//!
//! Architecture mapping:
//! - `handle(NCollection_HArray2(handle(Geom_Surface)))` (always built as
//!   `(1, NbU, 1, NbV)`) -> `Option<Vec<Vec<Surface3>>>` (outer index = U
//!   rank, inner index = V rank, 1-based ranks mapped to `i-1`/`j-1`).
//! - `handle(NCollection_HArray1(double))` (`(1, N+1)`) ->
//!   `Option<Vec<f64>>`.
//! - `handle(Geom_Surface)` -> the rcad value type
//!   [`rcad_kernel::geom::Surface3`]; evaluation uses the `SurfaceEval`
//!   trait (`point_at` = Value/D0, `derivatives` = D1, `derivatives2` = D2,
//!   `default_domain` = Bounds).
//! - `gp_Trsf` -> `glam::DAffine3` (the `transform_surface` kernel bridge).
//!
//! GAP carrier: `gp_Trsf2d` (TKMath) lives at the bottom of this file and is
//! only exercised by [`ShapeExtendCompositeSurface::global_to_local_transformation`].

use super::status::ShapeExtendParametrisation;
use rcad_kernel::geom::{transform_surface, Surface3, SurfaceEval};
use rcad_kernel::{CONFUSION, PCONFUSION};
use glam::{DVec2, DVec3};

/// OCCT Geom_Surface::ResD1 (Geom_Surface.hxx L54-59): point and first
/// partial derivatives.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceResD1 {
    pub point: DVec3,
    pub d1u: DVec3,
    pub d1v: DVec3,
}

/// OCCT Geom_Surface::ResD2 (Geom_Surface.hxx L62-70): point and partial
/// derivatives up to 2nd order.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceResD2 {
    pub point: DVec3,
    pub d1u: DVec3,
    pub d1v: DVec3,
    pub d2u: DVec3,
    pub d2v: DVec3,
    pub d2uv: DVec3,
}

/// OCCT Geom_Surface::ResD3 (Geom_Surface.hxx L73-83): point and partial
/// derivatives up to 3rd order.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceResD3 {
    pub point: DVec3,
    pub d1u: DVec3,
    pub d1v: DVec3,
    pub d2u: DVec3,
    pub d2v: DVec3,
    pub d2uv: DVec3,
    pub d3u: DVec3,
    pub d3v: DVec3,
    pub d3uuv: DVec3,
}

/// OCCT ShapeExtend_CompositeSurface (ShapeExtend_CompositeSurface.hxx
/// L67-307).
#[derive(Debug, Clone, Default)]
pub struct ShapeExtendCompositeSurface {
    /// OCCT myPatches (hxx L302).
    my_patches: Option<Vec<Vec<Surface3>>>,
    /// OCCT myUJointValues (hxx L303).
    my_u_joint_values: Option<Vec<f64>>,
    /// OCCT myVJointValues (hxx L304).
    my_v_joint_values: Option<Vec<f64>>,
    /// OCCT myUClosed (hxx L305).
    my_u_closed: bool,
    /// OCCT myVClosed (hxx L306).
    my_v_closed: bool,
}

impl ShapeExtendCompositeSurface {
    /// OCCT ShapeExtend_CompositeSurface() (cxx L32): empty constructor.
    pub fn new() -> Self {
        ShapeExtendCompositeSurface {
            my_patches: None,
            my_u_joint_values: None,
            my_v_joint_values: None,
            my_u_closed: false,
            my_v_closed: false,
        }
    }

    /// OCCT ShapeExtend_CompositeSurface(GridSurf, param) (cxx L36-41):
    /// initializes by a grid of surfaces (calls Init()).
    pub fn new_with_grid_param(
        grid_surf: Vec<Vec<Surface3>>,
        param: ShapeExtendParametrisation,
    ) -> Self {
        let mut surf = ShapeExtendCompositeSurface::new();
        surf.init_param(grid_surf, param);
        surf
    }

    /// OCCT ShapeExtend_CompositeSurface(GridSurf, UJoints, VJoints)
    /// (cxx L45-51): initializes by a grid of surfaces (calls Init()).
    pub fn new_with_grid_joints(
        grid_surf: Vec<Vec<Surface3>>,
        ujoints: &[f64],
        vjoints: &[f64],
    ) -> Self {
        let mut surf = ShapeExtendCompositeSurface::new();
        surf.init_joints(grid_surf, ujoints, vjoints);
        surf
    }

    /// OCCT Init(GridSurf, param) (cxx L55-66): initializes by a grid of
    /// surfaces; returns False when geometrical connectivity is not
    /// satisfied (the class is initialized even in that case).
    pub fn init_param(
        &mut self,
        grid_surf: Vec<Vec<Surface3>>,
        param: ShapeExtendParametrisation,
    ) -> bool {
        // OCCT L59-62: GridSurf.IsNull() -> false; myPatches = GridSurf.
        self.my_patches = Some(grid_surf);
        self.compute_joint_values(param);
        self.check_connectivity(CONFUSION)
    }

    /// OCCT Init(GridSurf, UJoints, VJoints) (cxx L70-92).
    pub fn init_joints(
        &mut self,
        grid_surf: Vec<Vec<Surface3>>,
        ujoints: &[f64],
        vjoints: &[f64],
    ) -> bool {
        self.my_patches = Some(grid_surf);

        let mut ok = true;
        if !self.set_u_joint_values(ujoints) || !self.set_v_joint_values(vjoints) {
            ok = false;
            self.compute_joint_values(ShapeExtendParametrisation::Natural);
            // OCCT_DEBUG warning (cxx L86-88) is compiled out.
        }

        if self.check_connectivity(CONFUSION) {
            ok
        } else {
            false
        }
    }

    /// OCCT NbUPatches() (cxx L96-99): returns number of patches in U
    /// direction (the HArray2 ColLength).
    pub fn nb_u_patches(&self) -> i32 {
        self.my_patches.as_ref().map(|p| p.len()).unwrap_or(0) as i32
    }

    /// OCCT NbVPatches() (cxx L103-106): returns number of patches in V
    /// direction (the HArray2 RowLength).
    pub fn nb_v_patches(&self) -> i32 {
        self.my_patches
            .as_ref()
            .and_then(|p| p.first())
            .map(|row| row.len())
            .unwrap_or(0) as i32
    }

    /// OCCT Patch(i, j) (cxx L110-113): returns one surface patch.
    pub fn patch_ij(&self, i: i32, j: i32) -> &Surface3 {
        &self.my_patches.as_ref().expect("myPatches")[(i - 1) as usize][(j - 1) as usize]
    }

    /// OCCT Patches() (cxx L117-121): returns grid of surfaces.
    pub fn patches(&self) -> &Option<Vec<Vec<Surface3>>> {
        &self.my_patches
    }

    /// OCCT UJointValues() (cxx L125-128): returns the array of U values
    /// corresponding to joint points.
    pub fn u_joint_values(&self) -> Option<&Vec<f64>> {
        self.my_u_joint_values.as_ref()
    }

    /// OCCT VJointValues() (cxx L132-135).
    pub fn v_joint_values(&self) -> Option<&Vec<f64>> {
        self.my_v_joint_values.as_ref()
    }

    /// OCCT UJointValue(i) (cxx L139-143): returns i-th joint value in U
    /// direction.
    pub fn u_joint_value(&self, i: i32) -> f64 {
        self.my_u_joint_values.as_ref().expect("myUJointValues")[(i - 1) as usize]
    }

    /// OCCT VJointValue(j) (cxx L146-149).
    pub fn v_joint_value(&self, j: i32) -> f64 {
        self.my_v_joint_values.as_ref().expect("myVJointValues")[(j - 1) as usize]
    }

    /// OCCT SetUJointValues(UJoints) (cxx L153-173): sets the array of U
    /// values corresponding to joint points; does nothing and returns False
    /// when the length is not NbUPatches()+1 or the values are not sorted.
    pub fn set_u_joint_values(&mut self, ujoints: &[f64]) -> bool {
        let nb_u = self.nb_u_patches();
        if ujoints.len() as i32 != nb_u + 1 {
            return false;
        }

        // OCCT L161-170: UJointValues = new HArray1(1, NbU + 1); the cursor
        // j starts at UJoints.Lower() (rcad slices are 0-based, Lower() = 0).
        let mut ujoint_values = vec![0.0f64; (nb_u + 1) as usize];
        let mut j = 0usize;
        for i in 1..=nb_u + 1 {
            ujoint_values[(i - 1) as usize] = ujoints[j];
            if i > 1 && ujoints[j] - ujoints[j - 1] < PCONFUSION {
                return false;
            }
            j += 1;
        }
        self.my_u_joint_values = Some(ujoint_values);
        true
    }

    /// OCCT SetVJointValues(VJoints) (cxx L177-197).
    pub fn set_v_joint_values(&mut self, vjoints: &[f64]) -> bool {
        let nb_v = self.nb_v_patches();
        if vjoints.len() as i32 != nb_v + 1 {
            return false;
        }

        let mut vjoint_values = vec![0.0f64; (nb_v + 1) as usize];
        let mut j = 0usize;
        for i in 1..=nb_v + 1 {
            vjoint_values[(i - 1) as usize] = vjoints[j];
            if i > 1 && vjoints[j] - vjoints[j - 1] < PCONFUSION {
                return false;
            }
            j += 1;
        }
        self.my_v_joint_values = Some(vjoint_values);
        true
    }

    /// OCCT SetUFirstValue(UFirst) (cxx L201-214): changes starting value for
    /// global U parametrisation (all other joint values are shifted
    /// accordingly).
    pub fn set_u_first_value(&mut self, ufirst: f64) {
        let Some(values) = self.my_u_joint_values.as_mut() else {
            return;
        };
        let shift = ufirst - values[0];
        let nb_u = values.len() as i32;
        for i in 1..=nb_u {
            values[(i - 1) as usize] += shift;
        }
    }

    /// OCCT SetVFirstValue(VFirst) (cxx L218-231).
    pub fn set_v_first_value(&mut self, vfirst: f64) {
        let Some(values) = self.my_v_joint_values.as_mut() else {
            return;
        };
        let shift = vfirst - values[0];
        let nb_v = values.len() as i32;
        for i in 1..=nb_v {
            values[(i - 1) as usize] += shift;
        }
    }

    /// OCCT LocateUParameter(U) (cxx L235-246): returns number of col that
    /// contains given (global) parameter.
    pub fn locate_u_parameter(&self, u: f64) -> i32 {
        let nb_patch = self.nb_u_patches();
        let values = self.my_u_joint_values.as_ref().expect("myUJointValues");
        for i in 2..=nb_patch {
            if u < values[(i - 1) as usize] {
                return i - 1;
            }
        }
        nb_patch
    }

    /// OCCT LocateVParameter(V) (cxx L250-261).
    pub fn locate_v_parameter(&self, v: f64) -> i32 {
        let nb_patch = self.nb_v_patches();
        let values = self.my_v_joint_values.as_ref().expect("myVJointValues");
        for i in 2..=nb_patch {
            if v < values[(i - 1) as usize] {
                return i - 1;
            }
        }
        nb_patch
    }

    /// OCCT LocateUVPoint(pnt, i, j) (cxx L265-269): returns number of row
    /// and col of surface that contains given point.  The C++ output
    /// parameters map to the returned tuple.
    pub fn locate_uv_point(&self, pnt: DVec2) -> (i32, i32) {
        let i = self.locate_u_parameter(pnt.x);
        let j = self.locate_v_parameter(pnt.y);
        (i, j)
    }

    /// OCCT Patch(U, V) (cxx L273-277): returns one surface patch that
    /// contains given (global) parameters.
    pub fn patch_uv(&self, u: f64, v: f64) -> &Surface3 {
        let i = self.locate_u_parameter(u);
        let j = self.locate_v_parameter(v);
        self.patch_ij(i, j)
    }

    /// OCCT Patch(pnt) (cxx L281-284): returns one surface patch that
    /// contains given point.
    pub fn patch_point(&self, pnt: DVec2) -> &Surface3 {
        let i = self.locate_u_parameter(pnt.x);
        let j = self.locate_v_parameter(pnt.y);
        self.patch_ij(i, j)
    }

    /// OCCT ULocalToGlobal(i, j, u) (cxx L288-296): converts local parameter
    /// u on patch i,j to global parameter U.
    pub fn u_local_to_global(&self, i: i32, j: i32, u: f64) -> f64 {
        let [u1, u2, _v1, _v2] = self.patch_ij(i, j).default_domain();
        let values = self.my_u_joint_values.as_ref().expect("myUJointValues");
        let scale = (values[i as usize] - values[(i - 1) as usize]) / (u2 - u1);
        // ! this formula is stable if u1 is infinite (cxx L294).
        u * scale + (values[(i - 1) as usize] - u1 * scale)
    }

    /// OCCT VLocalToGlobal(i, j, v) (cxx L300-308).
    pub fn v_local_to_global(&self, i: i32, j: i32, v: f64) -> f64 {
        let [_u1, _u2, v1, v2] = self.patch_ij(i, j).default_domain();
        let values = self.my_v_joint_values.as_ref().expect("myVJointValues");
        let scale = (values[j as usize] - values[(j - 1) as usize]) / (v2 - v1);
        // ! this formula is stable if v1 is infinite (cxx L306).
        v * scale + (values[(j - 1) as usize] - v1 * scale)
    }

    /// OCCT LocalToGlobal(i, j, uv) (cxx L312-324): converts local parameters
    /// uv on patch i,j to global parameters UV.
    pub fn local_to_global(&self, i: i32, j: i32, uv: DVec2) -> DVec2 {
        let [u1, u2, v1, v2] = self.patch_ij(i, j).default_domain();
        let uvalues = self.my_u_joint_values.as_ref().expect("myUJointValues");
        let vvalues = self.my_v_joint_values.as_ref().expect("myVJointValues");
        let scaleu = (uvalues[i as usize] - uvalues[(i - 1) as usize]) / (u2 - u1);
        let scalev = (vvalues[j as usize] - vvalues[(j - 1) as usize]) / (v2 - v1);
        DVec2::new(
            // ! this formula is stable if u1 or v1 is infinite (cxx L322).
            uv.x * scaleu + (uvalues[(i - 1) as usize] - u1 * scaleu),
            uv.y * scalev + (vvalues[(j - 1) as usize] - v1 * scalev),
        )
    }

    /// OCCT UGlobalToLocal(i, j, U) (cxx L328-336).
    pub fn u_global_to_local(&self, i: i32, j: i32, u: f64) -> f64 {
        let [u1, u2, _v1, _v2] = self.patch_ij(i, j).default_domain();
        let values = self.my_u_joint_values.as_ref().expect("myUJointValues");
        let scale = (u2 - u1) / (values[i as usize] - values[(i - 1) as usize]);
        // ! this formula is stable if u1 is infinite (cxx L334).
        u * scale + (u1 - values[(i - 1) as usize] * scale)
    }

    /// OCCT VGlobalToLocal(i, j, V) (cxx L340-348).
    pub fn v_global_to_local(&self, i: i32, j: i32, v: f64) -> f64 {
        let [_u1, _u2, v1, v2] = self.patch_ij(i, j).default_domain();
        let values = self.my_v_joint_values.as_ref().expect("myVJointValues");
        let scale = (v2 - v1) / (values[j as usize] - values[(j - 1) as usize]);
        // ! this formula is stable if v1 is infinite (cxx L346).
        v * scale + (v1 - values[(j - 1) as usize] * scale)
    }

    /// OCCT GlobalToLocal(i, j, UV) (cxx L352-365).
    pub fn global_to_local(&self, i: i32, j: i32, uv: DVec2) -> DVec2 {
        let [u1, u2, v1, v2] = self.patch_ij(i, j).default_domain();
        let uvalues = self.my_u_joint_values.as_ref().expect("myUJointValues");
        let vvalues = self.my_v_joint_values.as_ref().expect("myVJointValues");
        let scaleu = (u2 - u1) / (uvalues[i as usize] - uvalues[(i - 1) as usize]);
        let scalev = (v2 - v1) / (vvalues[j as usize] - vvalues[(j - 1) as usize]);
        DVec2::new(
            // ! this formula is stable if u1 or v1 is infinite (cxx L363).
            uv.x * scaleu + (u1 - uvalues[(i - 1) as usize] * scaleu),
            uv.y * scalev + (v1 - vvalues[(j - 1) as usize] * scalev),
        )
    }

    /// OCCT GlobalToLocalTransformation(i, j, uFact, Trsf) (cxx L369-393):
    /// computes transformation operator and uFactor describing the affine
    /// transformation required to convert global parameters on the composite
    /// surface to local parameters on patch (i,j):
    /// uv = ( uFactor, 1. ) X Trsf * UV.  Returns True if the transformation
    /// is not an identity.
    pub fn global_to_local_transformation(&self, i: i32, j: i32) -> (bool, f64, Trsf2d) {
        let [u1, u2, v1, v2] = self.patch_ij(i, j).default_domain();
        let uvalues = self.my_u_joint_values.as_ref().expect("myUJointValues");
        let vvalues = self.my_v_joint_values.as_ref().expect("myVJointValues");

        let scaleu = (u2 - u1) / (uvalues[i as usize] - uvalues[(i - 1) as usize]);
        let scalev = (v2 - v1) / (vvalues[j as usize] - vvalues[(j - 1) as usize]);
        let shift = DVec2::new(
            u1 / scaleu - uvalues[(i - 1) as usize],
            v1 / scalev - vvalues[(j - 1) as usize],
        );

        let ufact = scaleu / scalev;
        let mut shift_trsf = Trsf2d::default();
        let mut scale_trsf = Trsf2d::default();
        if shift.x != 0. || shift.y != 0. {
            shift_trsf.set_translation(shift);
        }
        if scalev != 1. {
            scale_trsf.set_scale(DVec2::new(0., 0.), scalev);
        }
        let trsf = scale_trsf.multiplied(&shift_trsf);
        let not_identity = ufact != 1. || trsf.form() != TrsfForm::Identity;
        (not_identity, ufact, trsf)
    }

    // ------------------------------------------------------------------
    // Inherited methods (from Geom_Geometry and Geom_Surface), cxx L395-534
    // ------------------------------------------------------------------

    /// OCCT Transform(T) (cxx L401-414): applies transformation to all the
    /// patches.
    pub fn transform(&mut self, t: &glam::DAffine3) {
        let Some(patches) = self.my_patches.as_mut() else {
            return;
        };
        for row in patches.iter_mut() {
            for patch in row.iter_mut() {
                *patch = transform_surface(patch, t);
            }
        }
    }

    /// OCCT Copy() (cxx L418-437): returns a copy of the surface.
    pub fn copy(&self) -> ShapeExtendCompositeSurface {
        let mut surf = ShapeExtendCompositeSurface::new();
        let Some(patches) = self.my_patches.as_ref() else {
            return surf;
        };
        // OCCT L426-434: each patch is copied through Geom_Geometry::Copy()
        // (the rcad value clone) and the copy is initialized with the default
        // (Natural) parametrisation.
        let patches: Vec<Vec<Surface3>> =
            patches.iter().map(|row| row.iter().map(|p| p.clone()).collect()).collect();
        surf.init_param(patches, ShapeExtendParametrisation::Natural);
        surf
    }

    /// OCCT UReverse() (cxx L441): NOT IMPLEMENTED (does nothing).
    pub fn u_reverse(&mut self) {}

    /// OCCT UReversedParameter(U) (cxx L445-448): returns U.
    pub fn u_reversed_parameter(&self, u: f64) -> f64 {
        u
    }

    /// OCCT VReverse() (cxx L452): NOT IMPLEMENTED (does nothing).
    pub fn v_reverse(&mut self) {}

    /// OCCT VReversedParameter(V) (cxx L456-459): returns V.
    pub fn v_reversed_parameter(&self, v: f64) -> f64 {
        v
    }

    /// OCCT Bounds(U1, U2, V1, V2) (cxx L463-469): returns the parametric
    /// bounds of grid.  The C++ output parameters map to the returned
    /// `[U1, U2, V1, V2]`.
    pub fn bounds(&self) -> [f64; 4] {
        [
            self.u_joint_value(1),
            self.u_joint_value(self.nb_u_patches() + 1),
            self.v_joint_value(1),
            self.v_joint_value(self.nb_v_patches() + 1),
        ]
    }

    /// OCCT IsUPeriodic() (cxx L473-476): returns False.
    pub fn is_u_periodic(&self) -> bool {
        false
    }

    /// OCCT IsVPeriodic() (cxx L480-483): returns False.
    pub fn is_v_periodic(&self) -> bool {
        false
    }

    /// OCCT UIso(U) (cxx L487-491): NOT IMPLEMENTED (returns Null curve).
    pub fn u_iso(&self, _u: f64) -> Option<rcad_kernel::geom::Curve3> {
        None
    }

    /// OCCT VIso(V) (cxx L495-499): NOT IMPLEMENTED (returns Null curve).
    pub fn v_iso(&self, _v: f64) -> Option<rcad_kernel::geom::Curve3> {
        None
    }

    /// OCCT Continuity() (cxx L503-506): returns C0.
    pub fn continuity(&self) -> rcad_kernel::math::GeomAbsShape {
        rcad_kernel::math::GeomAbsShape::C0
    }

    /// OCCT IsCNu(N) (cxx L510-513): returns True if N <= 0.
    pub fn is_cn_u(&self, n: i32) -> bool {
        n <= 0
    }

    /// OCCT IsCNv(N) (cxx L517-520): returns True if N <= 0.
    pub fn is_cn_v(&self, n: i32) -> bool {
        n <= 0
    }

    /// OCCT IsUClosed() (cxx L524-527): True when the grid is closed in U
    /// direction (checked by CheckConnectivity with Precision::Confusion).
    pub fn is_u_closed(&self) -> bool {
        self.my_u_closed
    }

    /// OCCT IsVClosed() (cxx L531-534).
    pub fn is_v_closed(&self) -> bool {
        self.my_v_closed
    }

    /// OCCT EvalD0(U, V) (cxx L538-544): computes the point of parameter U,V
    /// on the grid.
    pub fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        let i = self.locate_u_parameter(u);
        let j = self.locate_v_parameter(v);
        let uv = self.global_to_local(i, j, DVec2::new(u, v));
        self.patch_ij(i, j).point_at(uv.x, uv.y)
    }

    /// OCCT EvalD1(U, V) (cxx L548-554).
    pub fn eval_d1(&self, u: f64, v: f64) -> SurfaceResD1 {
        let i = self.locate_u_parameter(u);
        let j = self.locate_v_parameter(v);
        let uv = self.global_to_local(i, j, DVec2::new(u, v));
        let (point, d1u, d1v) = self.patch_ij(i, j).derivatives(uv.x, uv.y);
        SurfaceResD1 { point, d1u, d1v }
    }

    /// OCCT EvalD2(U, V) (cxx L558-564).
    pub fn eval_d2(&self, u: f64, v: f64) -> SurfaceResD2 {
        let i = self.locate_u_parameter(u);
        let j = self.locate_v_parameter(v);
        let uv = self.global_to_local(i, j, DVec2::new(u, v));
        let (point, d1u, d1v, d2u, d2v, d2uv) = self.patch_ij(i, j).derivatives2(uv.x, uv.y);
        SurfaceResD2 {
            point,
            d1u,
            d1v,
            d2u,
            d2v,
            d2uv,
        }
    }

    /// OCCT EvalD3(U, V) (cxx L568-574).
    ///
    /// Architecture note: the rcad `Surface3` general derivative dispatcher
    /// is `Surface3::dn(u, v, Nu, Nv)` (the GeomAdaptor_DN re-host); the
    /// third-order fields are assembled from it in the Geom_Surface::ResD3
    /// layout (D3U, D3V, D3UUV).
    pub fn eval_d3(&self, u: f64, v: f64) -> SurfaceResD3 {
        let i = self.locate_u_parameter(u);
        let j = self.locate_v_parameter(v);
        let uv = self.global_to_local(i, j, DVec2::new(u, v));
        let patch = self.patch_ij(i, j);
        let (point, d1u, d1v, d2u, d2v, d2uv) = patch.derivatives2(uv.x, uv.y);
        let d3u = patch.dn(uv.x, uv.y, 3, 0);
        let d3v = patch.dn(uv.x, uv.y, 0, 3);
        let d3uuv = patch.dn(uv.x, uv.y, 2, 1);
        SurfaceResD3 {
            point,
            d1u,
            d1v,
            d2u,
            d2v,
            d2uv,
            d3u,
            d3v,
            d3uuv,
        }
    }

    /// OCCT EvalDN(U, V, Nu, Nv) (cxx L578-587).
    pub fn eval_dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
        let i = self.locate_u_parameter(u);
        let j = self.locate_v_parameter(v);
        let uv = self.global_to_local(i, j, DVec2::new(u, v));
        self.patch_ij(i, j).dn(uv.x, uv.y, nu, nv)
    }

    /// OCCT Value(pnt) (cxx L591-599): computes the point of parameter pnt on
    /// the grid (the patch D0 call of cxx L597).
    pub fn value(&self, pnt: DVec2) -> DVec3 {
        let i = self.locate_u_parameter(pnt.x);
        let j = self.locate_v_parameter(pnt.y);
        let uv = self.global_to_local(i, j, pnt);
        self.patch_ij(i, j).point_at(uv.x, uv.y)
    }

    /// OCCT ComputeJointValues(param) (cxx L603-653): computes Joint values
    /// according to parameter.
    pub fn compute_joint_values(&mut self, param: ShapeExtendParametrisation) {
        let nb_u = self.nb_u_patches();
        let nb_v = self.nb_v_patches();
        self.my_u_joint_values = Some(vec![0.0f64; (nb_u + 1) as usize]);
        self.my_v_joint_values = Some(vec![0.0f64; (nb_v + 1) as usize]);

        if param == ShapeExtendParametrisation::Natural {
            let mut u = 0.0f64;
            let mut v = 0.0f64;
            let patches = self.my_patches.as_ref().expect("myPatches");
            let uvalues = self.my_u_joint_values.as_mut().expect("myUJointValues");
            // OCCT L614-623: for i = 1..NbU over the 1st row: Bounds of patch
            // (i, 1); i == 1 stores U = U1 at rank 1; U += (U2 - U1).
            for i in 1..=nb_u {
                let [u1, u2, _v1, _v2] = patches[(i - 1) as usize][0].default_domain();
                if i == 1 {
                    u = u1;
                    uvalues[0] = u;
                }
                u += u2 - u1;
                uvalues[i as usize] = u;
            }
            let vvalues = self.my_v_joint_values.as_mut().expect("myVJointValues");
            // OCCT L624-633: for i = 1..NbV over the 1st column.
            for i in 1..=nb_v {
                let [u1, u2, v1, v2] = patches[0][(i - 1) as usize].default_domain();
                let _ = (u1, u2);
                if i == 1 {
                    v = v1;
                    vvalues[0] = v;
                }
                v += v2 - v1;
                vvalues[i as usize] = v;
            }
        } else {
            // OCCT L637: suppose param == ShapeExtend_Uniform.
            let mut stepu = 1.0f64;
            let mut stepv = 1.0f64;
            if param == ShapeExtendParametrisation::Unitary {
                stepu /= nb_u as f64;
                stepv /= nb_v as f64;
            }
            let uvalues = self.my_u_joint_values.as_mut().expect("myUJointValues");
            // OCCT L644-647: for i = 0..NbU (inclusive).
            for i in 0..=nb_u {
                uvalues[i as usize] = i as f64 * stepu;
            }
            let vvalues = self.my_v_joint_values.as_mut().expect("myVJointValues");
            for i in 0..=nb_v {
                vvalues[i as usize] = i as f64 * stepv;
            }
        }
    }

    /// OCCT CheckConnectivity(Prec) (cxx L675-758): checks geometrical
    /// connectivity of the patches, including closedness (sets fields
    /// myUClosed and myVClosed).
    pub fn check_connectivity(&mut self, prec: f64) -> bool {
        const NPOINTS: i32 = 23;
        let mut ok = true;
        let nb_u = self.nb_u_patches();
        let nb_v = self.nb_v_patches();

        // check in u direction (cxx L682-716): for (i = 1, j = NbU; i <= NbU;
        // j = i++) — j takes the previous i (NbU when i == 1).
        for i in 1..=nb_u {
            let j = if i == 1 { nb_u } else { i - 1 };
            let mut maxdist2 = 0.0f64;
            for k in 1..=nb_v {
                let sj = self.patch_ij(j, k);
                let si = self.patch_ij(i, k);
                let [_uj1, uj2, vj1, vj2] = get_limited_bounds(sj);
                let [ui1, _ui2, vi1, vi2] = get_limited_bounds(si);
                let stepj = (vj2 - vj1) / (NPOINTS - 1) as f64;
                let stepi = (vi2 - vi1) / (NPOINTS - 1) as f64;
                for isample in 0..NPOINTS {
                    let parj = vj1 + stepj * isample as f64;
                    let pari = vi1 + stepi * isample as f64;
                    let dist2 = sj.point_at(uj2, parj).distance_squared(si.point_at(ui1, pari));
                    if maxdist2 < dist2 {
                        maxdist2 = dist2;
                    }
                }
            }
            if i == 1 {
                self.my_u_closed = maxdist2 <= prec * prec;
            } else if maxdist2 > prec * prec {
                ok = false;
            }
        }

        // check in v direction (cxx L718-751).
        for i in 1..=nb_v {
            let j = if i == 1 { nb_v } else { i - 1 };
            let mut maxdist2 = 0.0f64;
            for k in 1..=nb_u {
                let sj = self.patch_ij(k, j);
                let si = self.patch_ij(k, i);
                let [uj1, uj2, _vj1, vj2] = get_limited_bounds(sj);
                let [ui1, ui2, vi1, _vi2] = get_limited_bounds(si);
                let stepj = (uj2 - uj1) / (NPOINTS - 1) as f64;
                let stepi = (ui2 - ui1) / (NPOINTS - 1) as f64;
                for isample in 0..NPOINTS {
                    let parj = uj1 + stepj * isample as f64;
                    let pari = ui1 + stepi * isample as f64;
                    let dist2 = sj.point_at(parj, vj2).distance_squared(si.point_at(pari, vi1));
                    if maxdist2 < dist2 {
                        maxdist2 = dist2;
                    }
                }
            }
            if i == 1 {
                self.my_v_closed = maxdist2 <= prec * prec;
            } else if maxdist2 > prec * prec {
                ok = false;
            }
        }

        // OCCT_DEBUG warning (cxx L753-756) is compiled out.
        ok
    }
}

/// OCCT LimitValue (cxx L657-660): an infinite parameter clamps to
/// -10000./10000.
fn limit_value(par: f64) -> f64 {
    if rcad_kernel::is_infinite_value(par) {
        if par < 0.0 {
            -10000.0
        } else {
            10000.0
        }
    } else {
        par
    }
}

/// OCCT GetLimitedBounds (cxx L662-673): the patch bounds with infinite
/// values clamped.
fn get_limited_bounds(surf: &Surface3) -> [f64; 4] {
    let [u1, u2, v1, v2] = surf.default_domain();
    [limit_value(u1), limit_value(u2), limit_value(v1), limit_value(v2)]
}

/// OCCT gp_TrsfForm (gp_Trsf.hxx) — the subset of forms the gp_Trsf2d GAP
/// carrier distinguishes (Identity, Translation, Scale, CompoundTrsf);
/// other forms are reported as Other and are pending the TKMath 2D batch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrsfForm {
    Identity,
    Translation,
    Scale,
    CompoundTrsf,
    Other,
}

/// GAP carrier for OCCT `gp_Trsf2d` (TKMath, `gp_Trsf.hxx` / `gp_Trsf.cxx`):
/// a 2D transformation with a separate scale factor and form tag, following
/// the gp_Trsf member layout (scale, shape, matrix, loc).  Only the forms
/// exercised by
/// [`ShapeExtendCompositeSurface::global_to_local_transformation`]
/// (Identity / Translation / Scale and their products) are carried; the
/// general 2D transformation set is pending the TKMath batch.
///
/// Point transformation follows gp_Trsf::Transformed: `P' = matrix * P *
/// scale + loc`.
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

    /// OCCT gp_Trsf::SetScale(gp_Pnt, S) (gp_Trsf.cxx L159-168): scale = S,
    /// identity matrix, loc = P * (1 - S), form = gp_Scale (with the
    /// Standard_ConstructionError raise for |S| <= Resolution kept as a
    /// panic).
    pub fn set_scale(&mut self, p: DVec2, s: f64) {
        self.form = TrsfForm::Scale;
        self.scale = s;
        self.loc = p;
        assert!(self.scale.abs() > f64::MIN_POSITIVE, "gp_Trsf::SetScaleFactor");
        self.matrix = [[1.0, 0.0], [0.0, 1.0]];
        self.loc *= 1.0 - s;
    }

    /// OCCT gp_Trsf::Form().
    pub fn form(&self) -> TrsfForm {
        self.form
    }

    /// OCCT gp_Trsf::Multiplied(T) restricted to the (Identity | Scale) x
    /// (Identity | Translation | Scale) combinations
    /// (gp_Trsf.cxx L430-530 branches): `A * B` applies B first, then A.
    pub fn multiplied(&self, t: &Trsf2d) -> Trsf2d {
        let mut res = *self;
        if t.form == TrsfForm::Identity {
            // OCCT L432-434: T identity -> this unchanged.
        } else if res.form == TrsfForm::Identity {
            // OCCT L435-441: this identity -> copy T.
            res.form = t.form;
            res.scale = t.scale;
            res.loc = t.loc;
            res.matrix = t.matrix;
        } else if res.form == TrsfForm::Scale && t.form == TrsfForm::Translation {
            // OCCT L501-509: (Scale || PntMirror) && T Translation:
            // Tloc = T.loc * scale; loc += Tloc.
            let tloc = t.loc * res.scale;
            res.loc += tloc;
        } else if res.form == TrsfForm::Scale && t.form == TrsfForm::Scale {
            // OCCT L475-479: Scale && Scale: loc += T.loc * scale;
            // scale = scale * T.scale.
            let tloc = t.loc * res.scale;
            res.loc += tloc;
            res.scale = res.scale * t.scale;
        } else {
            // General forms pending the TKMath 2D batch.
            res.form = TrsfForm::CompoundTrsf;
            let tloc = t.loc * res.scale;
            res.loc += tloc;
            res.scale = res.scale * t.scale;
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
}
