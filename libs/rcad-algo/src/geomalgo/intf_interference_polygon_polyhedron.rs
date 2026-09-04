//! OCCT Intf_InterferencePolygonPolyhedron (TKGeomAlgo Intf package) —
//! interference between a polygon of segments and a polyhedron of triangles
//! (or between straight lines and a polyhedron).
//!
//! 1:1 translation of `Intf_InterferencePolygonPolyhedron.gxx` (L1-1368)
//! over the template parameters `Polygon3d`/`ToolPolygon3d` and
//! `Polyhedron`/`ToolPolyh`, mapped to the [`ToolPolygon3d`]/[`ToolPolyh`]
//! traits with 1-based index semantics like the OCCT tools:
//! - IntCurveSurface_TheInterferenceOfHInter — Polygon3d =
//!   IntCurveSurface_ThePolygonOfHInter, Polyhedron =
//!   IntCurveSurface_ThePolyhedronOfHInter (consumed by
//!   IntCurveSurface_HInter, 2a-3);
//! - HLRBRep_TheInterferenceOfInterCSurf — the HLRBRep polygon/polyhedron
//!   wrappers (Stage 3b supplies the tool impls over the same traits).
//!
//! The Intf.cxx package helper `Intf::PlaneEquation` is included here.

use glam::DVec3;
use rcad_kernel::geom::Line3;
use rcad_kernel::math::bnd::{BndBox, BoundSortBox};
use rcad_kernel::math::direct_polynomial_roots::epsilon;

use super::intf::{IntfPIType, IntfSectionPoint, IntfTool};
use super::intf_interference::Interference;

/// OCCT gp::Resolution() = 1e-15.
const GP_RESOLUTION: f64 = 1e-15;
/// OCCT Precision::Computational() (Precision.hxx).
const PRECISION_COMPUTATIONAL: f64 = 1e-14;

/// OCCT gxx L34: `static const int Pourcent3[4] = {0, 1, 2, 0};`
const POURCENT3: [usize; 4] = [0, 1, 2, 0];

/// OCCT gxx L770/L1085: `Epsilon(1000.)` — the ULP gap at 1000.0.
fn epsilon_1000() -> f64 {
    epsilon(1000.0)
}

/// OCCT gxx L36-56: `static bool IsInSegment(...)`.
fn is_in_segment(
    p1_p2: DVec3,
    p1_p: DVec3,
    n_p1_p2: f64,
    param: &mut f64,
    tolerance: f64,
) -> bool {
    if n_p1_p2 < PRECISION_COMPUTATIONAL {
        return false;
    }
    *param = p1_p2.dot(p1_p);
    *param /= n_p1_p2;
    if *param > (n_p1_p2 + tolerance) {
        return false;
    }
    if *param < -tolerance {
        return false;
    }
    *param /= n_p1_p2;
    if *param < 0.0 {
        *param = 0.0;
    }
    if *param > 1.0 {
        *param = 1.0;
    }
    true
}

/// OCCT Intf::PlaneEquation (Intf.cxx L20-42).
pub fn plane_equation(p1: DVec3, p2: DVec3, p3: DVec3) -> (DVec3, f64) {
    let v1 = p2 - p1;
    let v2 = p3 - p2;
    let v3 = p1 - p3;
    let mut normal_vector = v1.cross(v2) + v2.cross(v3) + v3.cross(v1);
    let a_norm_len = normal_vector.length();
    let polar_distance = if a_norm_len < GP_RESOLUTION {
        0.0
    } else {
        normal_vector /= a_norm_len;
        normal_vector.dot(p1)
    };
    (normal_vector, polar_distance)
}

/// OCCT `ToolPolygon3d` template parameter — the static accessors over the
/// polygon of segments (IntCurveSurface_ThePolygonToolOfHInter for the HInter
/// instantiation, HLRBRep_ThePolygonToolOfInterCSurf for the CInter one).
pub trait ToolPolygon3d<P> {
    /// OCCT ToolPolygon3d::DeflectionOverEstimation(thePolyg).
    fn deflection_over_estimation(the_polyg: &P) -> f64;
    /// OCCT ToolPolygon3d::Bounding(thePolyg).
    fn bounding(the_polyg: &P) -> &BndBox;
    /// OCCT ToolPolygon3d::Closed(thePolyg).
    fn closed(the_polyg: &P) -> bool;
    /// OCCT ToolPolygon3d::NbSegments(thePolyg).
    fn nb_segments(the_polyg: &P) -> usize;
    /// OCCT ToolPolygon3d::BeginOfSeg(thePolyg, iLin) — iLin is 1-based.
    fn begin_of_seg(the_polyg: &P, i_lin: usize) -> DVec3;
    /// OCCT ToolPolygon3d::EndOfSeg(thePolyg, iLin) — iLin is 1-based.
    fn end_of_seg(the_polyg: &P, i_lin: usize) -> DVec3;
}

/// OCCT `ToolPolyh` template parameter — the static accessors over the
/// polyhedron of triangles (IntCurveSurface_ThePolyhedronToolOfHInter for
/// the HInter instantiation, HLRBRep_ThePolyhedronToolOfInterCSurf for the
/// CInter one).  All index arguments are 1-based like the OCCT tools.
pub trait ToolPolyh<H> {
    /// OCCT ToolPolyh::DeflectionOverEstimation(thePolyh).
    fn deflection_over_estimation(the_polyh: &H) -> f64;
    /// OCCT ToolPolyh::Bounding(thePolyh).
    fn bounding(the_polyh: &H) -> &BndBox;
    /// OCCT ToolPolyh::ComponentsBounding(thePolyh).
    fn components_bounding(the_polyh: &H) -> &[BndBox];
    /// OCCT ToolPolyh::Triangle(thePolyh, TTri, pTri0, pTri1, pTri2).
    fn triangle(the_polyh: &H, ttri: usize) -> [usize; 3];
    /// OCCT ToolPolyh::Point(thePolyh, I).
    fn point(the_polyh: &H, i: usize) -> DVec3;
    /// OCCT ToolPolyh::TriConnex(thePolyh, Triang, Pivot, Pedge, TriCon,
    /// PEdg) — returns (TriCon, PEdg).
    fn tri_connex(the_polyh: &H, triang: usize, pivot: usize, pedge: usize) -> (i32, i32);
    /// OCCT ToolPolyh::IsOnBound(thePolyh, indP1, indP2).
    fn is_on_bound(the_polyh: &H, ind_p1: usize, ind_p2: usize) -> bool;
    /// OCCT ToolPolyh::GetBorderDeflection(thePolyh).
    fn get_border_deflection(the_polyh: &H) -> f64;
}

/// OCCT Intf_InterferencePolygonPolyhedron — the interference engine.  The
/// class stores no template-typed state (only the base sequences,
/// BeginOfClosedPolygon and iLin), so the P/H tool pairs are generic method
/// parameters like in the other IntCurve ports.
pub struct InterferencePolygonPolyhedron {
    /// OCCT base subobject Intf_Interference(false).
    pub interf: Interference,
    /// OCCT BeginOfClosedPolygon.
    begin_of_closed_polygon: bool,
    /// OCCT iLin.
    i_lin: i32,
}

impl InterferencePolygonPolyhedron {
    /// OCCT Intf_InterferencePolygonPolyhedron() (gxx L60-65).
    pub fn new() -> Self {
        InterferencePolygonPolyhedron {
            interf: Interference::with_self(false),
            begin_of_closed_polygon: false,
            i_lin: 0,
        }
    }

    /// OCCT Intf_InterferencePolygonPolyhedron(thePolyg, thePolyh)
    /// (gxx L73-88).
    pub fn new_polygon_polyhedron<P, H, TP, TH>(the_polyg: &P, the_polyh: &H) -> Self
    where
        TP: ToolPolygon3d<P>,
        TH: ToolPolyh<H>,
    {
        let mut r = InterferencePolygonPolyhedron::new();
        let mut tolerance = TP::deflection_over_estimation(the_polyg)
            + TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        r.interf.set_tolerance(tolerance);

        if !TP::bounding(the_polyg).is_out_box(TH::bounding(the_polyh)) {
            r.interference::<P, H, TP, TH>(the_polyg, the_polyh);
        }
        r
    }

    /// OCCT Intf_InterferencePolygonPolyhedron(thePolyg, thePolyh, PolyhGrid)
    /// (gxx L90-106).
    pub fn new_polygon_polyhedron_grid<P, H, TP, TH>(
        the_polyg: &P,
        the_polyh: &H,
        polyh_grid: &mut BoundSortBox,
    ) -> Self
    where
        TP: ToolPolygon3d<P>,
        TH: ToolPolyh<H>,
    {
        let mut r = InterferencePolygonPolyhedron::new();
        let mut tolerance = TP::deflection_over_estimation(the_polyg)
            + TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        r.interf.set_tolerance(tolerance);

        if !TP::bounding(the_polyg).is_out_box(TH::bounding(the_polyh)) {
            r.interference_grid::<P, H, TP, TH>(the_polyg, the_polyh, polyh_grid);
        }
        r
    }

    /// OCCT Intf_InterferencePolygonPolyhedron(theLin, thePolyh)
    /// (gxx L114-147).
    pub fn new_line_polyhedron<H, TH>(the_lin: &Line3, the_polyh: &H) -> Self
    where
        TH: ToolPolyh<H>,
    {
        let mut r = InterferencePolygonPolyhedron::new();
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        r.interf.set_tolerance(tolerance);

        r.begin_of_closed_polygon = false;

        let mut polyh_grid = BoundSortBox::new();
        polyh_grid.initialize_with_enclosing(
            TH::bounding(the_polyh).clone(),
            TH::components_bounding(the_polyh).to_vec(),
        );
        r.i_lin = 0;

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        btoo.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

        let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
        for &ind_tri in &liste {
            r.intersect::<H, TH>(
                the_lin.origin,
                the_lin.origin + the_lin.direction,
                true,
                ind_tri,
                the_polyh,
            );
        }
        r
    }

    /// OCCT Intf_InterferencePolygonPolyhedron(theLins, thePolyh)
    /// (gxx L155-193).
    pub fn new_lines_polyhedron<H, TH>(the_lins: &[Line3], the_polyh: &H) -> Self
    where
        TH: ToolPolyh<H>,
    {
        let mut r = InterferencePolygonPolyhedron::new();
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        r.interf.set_tolerance(tolerance);

        let mut bof_lin = BndBox::new();
        let mut b_too = IntfTool::new();
        r.begin_of_closed_polygon = false;

        let mut polyh_grid = BoundSortBox::new();
        polyh_grid.initialize_with_enclosing(
            TH::bounding(the_polyh).clone(),
            TH::components_bounding(the_polyh).to_vec(),
        );

        r.i_lin = 0;
        while (r.i_lin as usize) < the_lins.len() {
            r.i_lin += 1;
            let the_lin = &the_lins[(r.i_lin - 1) as usize];

            b_too.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

            let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
            for &ind_tri in &liste {
                r.intersect::<H, TH>(
                    the_lin.origin,
                    the_lin.origin + the_lin.direction,
                    true,
                    ind_tri,
                    the_polyh,
                );
            }
        }
        r
    }

    /// OCCT Perform(thePolyg, thePolyh) (gxx L197-210).
    pub fn perform_polygon_polyhedron<P, H, TP, TH>(&mut self, the_polyg: &P, the_polyh: &H)
    where
        TP: ToolPolygon3d<P>,
        TH: ToolPolyh<H>,
    {
        self.interf.self_interference(false);
        let mut tolerance = TP::deflection_over_estimation(the_polyg)
            + TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        self.interf.set_tolerance(tolerance);

        if !TP::bounding(the_polyg).is_out_box(TH::bounding(the_polyh)) {
            self.interference::<P, H, TP, TH>(the_polyg, the_polyh);
        }
    }

    /// OCCT Perform(theLin, thePolyh) (gxx L214-245).
    pub fn perform_line<H, TH>(&mut self, the_lin: &Line3, the_polyh: &H)
    where
        TH: ToolPolyh<H>,
    {
        self.interf.self_interference(false);
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        self.interf.set_tolerance(tolerance);

        self.begin_of_closed_polygon = false;

        let mut polyh_grid = BoundSortBox::new();
        polyh_grid.initialize_with_enclosing(
            TH::bounding(the_polyh).clone(),
            TH::components_bounding(the_polyh).to_vec(),
        );

        self.i_lin = 0;

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        btoo.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

        let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
        for &ind_tri in &liste {
            self.intersect::<H, TH>(
                the_lin.origin,
                the_lin.origin + the_lin.direction,
                true,
                ind_tri,
                the_polyh,
            );
        }
    }

    /// OCCT Perform(theLins, thePolyh) (gxx L253-288).
    pub fn perform_lines<H, TH>(&mut self, the_lins: &[Line3], the_polyh: &H)
    where
        TH: ToolPolyh<H>,
    {
        self.interf.self_interference(false);
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        self.interf.set_tolerance(tolerance);

        let mut bof_lin = BndBox::new();
        let mut b_too = IntfTool::new();
        self.begin_of_closed_polygon = false;

        let mut polyh_grid = BoundSortBox::new();
        polyh_grid.initialize_with_enclosing(
            TH::bounding(the_polyh).clone(),
            TH::components_bounding(the_polyh).to_vec(),
        );

        self.i_lin = 0;
        while (self.i_lin as usize) < the_lins.len() {
            self.i_lin += 1;
            let the_lin = &the_lins[(self.i_lin - 1) as usize];

            b_too.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

            let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
            for &ind_tri in &liste {
                self.intersect::<H, TH>(
                    the_lin.origin,
                    the_lin.origin + the_lin.direction,
                    true,
                    ind_tri,
                    the_polyh,
                );
            }
        }
    }

    /// OCCT Interference(thePolyg, thePolyh) (gxx L296-350): compare the
    /// boundings between the segments of the polygon and the facets of the
    /// polyhedron.
    fn interference<P, H, TP, TH>(&mut self, the_polyg: &P, the_polyh: &H)
    where
        TP: ToolPolygon3d<P>,
        TH: ToolPolyh<H>,
    {
        let mut bof_seg = BndBox::new();

        let mut polyh_grid = BoundSortBox::new();
        polyh_grid.initialize_with_enclosing(
            TH::bounding(the_polyh).clone(),
            TH::components_bounding(the_polyh).to_vec(),
        );

        self.begin_of_closed_polygon = TP::closed(the_polyg);

        let def_ph = TH::deflection_over_estimation(the_polyh);

        self.i_lin = 0;
        while (self.i_lin as usize) < TP::nb_segments(the_polyg) {
            self.i_lin += 1;

            bof_seg.set_void();
            bof_seg.add_point(TP::begin_of_seg(the_polyg, self.i_lin as usize));
            bof_seg.add_point(TP::end_of_seg(the_polyg, self.i_lin as usize));
            bof_seg.enlarge(TP::deflection_over_estimation(the_polyg));

            let maliste: Vec<i32> = polyh_grid.compare(&bof_seg).clone();
            for &ind_tri in &maliste {
                let p1 = TP::begin_of_seg(the_polyg, self.i_lin as usize);
                let p2 = TP::end_of_seg(the_polyg, self.i_lin as usize);
                let p_tri = TH::triangle(the_polyh, ind_tri as usize);
                let pa = TH::point(the_polyh, p_tri[0]);
                let pb = TH::point(the_polyh, p_tri[1]);
                let pc = TH::point(the_polyh, p_tri[2]);
                let pa_pb = pb - pa;
                let pa_pc = pc - pa;
                let mut normale = pa_pb.cross(pa_pc);
                let norm_normale = normale.length();
                if norm_normale < 1e-14 {
                    continue;
                }
                normale *= def_ph / norm_normale;
                let p1m = p1 - normale;
                let p1p = p1 + normale;
                let p2m = p2 - normale;
                let p2p = p2 + normale;
                self.intersect::<H, TH>(p1m, p2p, false, ind_tri, the_polyh);
                self.intersect::<H, TH>(p1p, p2m, false, ind_tri, the_polyh);
            }
            self.begin_of_closed_polygon = false;
        }
    }

    /// OCCT Intf_InterferencePolygonPolyhedron(theLin, thePolyh, PolyhGrid)
    /// (gxx L358-388).
    pub fn new_line_polyhedron_grid<H, TH>(
        the_lin: &Line3,
        the_polyh: &H,
        polyh_grid: &mut BoundSortBox,
    ) -> Self
    where
        TH: ToolPolyh<H>,
    {
        let mut r = InterferencePolygonPolyhedron::new();
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        r.interf.set_tolerance(tolerance);

        r.begin_of_closed_polygon = false;

        r.i_lin = 0;

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        btoo.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

        let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
        for &ind_tri in &liste {
            r.intersect::<H, TH>(
                the_lin.origin,
                the_lin.origin + the_lin.direction,
                true,
                ind_tri,
                the_polyh,
            );
        }
        r
    }

    /// OCCT Intf_InterferencePolygonPolyhedron(theLins, thePolyh, PolyhGrid)
    /// (gxx L396-430).
    pub fn new_lines_polyhedron_grid<H, TH>(
        the_lins: &[Line3],
        the_polyh: &H,
        polyh_grid: &mut BoundSortBox,
    ) -> Self
    where
        TH: ToolPolyh<H>,
    {
        let mut r = InterferencePolygonPolyhedron::new();
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        r.interf.set_tolerance(tolerance);

        let mut bof_lin = BndBox::new();
        let mut b_too = IntfTool::new();
        r.begin_of_closed_polygon = false;

        r.i_lin = 0;
        while (r.i_lin as usize) < the_lins.len() {
            r.i_lin += 1;
            let the_lin = &the_lins[(r.i_lin - 1) as usize];

            b_too.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

            let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
            for &ind_tri in &liste {
                r.intersect::<H, TH>(
                    the_lin.origin,
                    the_lin.origin + the_lin.direction,
                    true,
                    ind_tri,
                    the_polyh,
                );
            }
        }
        r
    }

    /// OCCT Perform(thePolyg, thePolyh, PolyhGrid) (gxx L434-448).
    pub fn perform_polygon_polyhedron_grid<P, H, TP, TH>(
        &mut self,
        the_polyg: &P,
        the_polyh: &H,
        polyh_grid: &mut BoundSortBox,
    ) where
        TP: ToolPolygon3d<P>,
        TH: ToolPolyh<H>,
    {
        self.interf.self_interference(false);
        let mut tolerance = TP::deflection_over_estimation(the_polyg)
            + TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        self.interf.set_tolerance(tolerance);

        if !TP::bounding(the_polyg).is_out_box(TH::bounding(the_polyh)) {
            self.interference_grid::<P, H, TP, TH>(the_polyg, the_polyh, polyh_grid);
        }
    }

    /// OCCT Perform(theLin, thePolyh, PolyhGrid) (gxx L452-482).
    pub fn perform_line_grid<H, TH>(
        &mut self,
        the_lin: &Line3,
        the_polyh: &H,
        polyh_grid: &mut BoundSortBox,
    ) where
        TH: ToolPolyh<H>,
    {
        self.interf.self_interference(false);
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        self.interf.set_tolerance(tolerance);

        self.begin_of_closed_polygon = false;

        self.i_lin = 0;

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        btoo.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

        let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
        for &ind_tri in &liste {
            self.intersect::<H, TH>(
                the_lin.origin,
                the_lin.origin + the_lin.direction,
                true,
                ind_tri,
                the_polyh,
            );
        }
    }

    /// OCCT Perform(theLins, thePolyh, PolyhGrid) (gxx L490-523).
    pub fn perform_lines_grid<H, TH>(
        &mut self,
        the_lins: &[Line3],
        the_polyh: &H,
        polyh_grid: &mut BoundSortBox,
    ) where
        TH: ToolPolyh<H>,
    {
        self.interf.self_interference(false);
        let mut tolerance = TH::deflection_over_estimation(the_polyh);
        if tolerance == 0.0 {
            tolerance = epsilon_1000();
        }
        self.interf.set_tolerance(tolerance);

        let mut bof_lin = BndBox::new();
        let mut b_too = IntfTool::new();
        self.begin_of_closed_polygon = false;

        self.i_lin = 0;
        while (self.i_lin as usize) < the_lins.len() {
            self.i_lin += 1;
            let the_lin = &the_lins[(self.i_lin - 1) as usize];

            b_too.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

            let liste: Vec<i32> = polyh_grid.compare(&bof_lin).clone();
            for &ind_tri in &liste {
                self.intersect::<H, TH>(
                    the_lin.origin,
                    the_lin.origin + the_lin.direction,
                    true,
                    ind_tri,
                    the_polyh,
                );
            }
        }
    }

    /// OCCT Interference(thePolyg, thePolyh, PolyhGrid) (gxx L531-623) — the
    /// grid-driven variant with the MKK boundary enlargement of 2007-10-25.
    fn interference_grid<P, H, TP, TH>(
        &mut self,
        the_polyg: &P,
        the_polyh: &H,
        polyh_grid: &mut BoundSortBox,
    ) where
        TP: ToolPolygon3d<P>,
        TH: ToolPolyh<H>,
    {
        let mut bof_seg = BndBox::new();

        self.begin_of_closed_polygon = TP::closed(the_polyg);

        self.i_lin = 0;
        let nb_segments = TP::nb_segments(the_polyg);
        while (self.i_lin as usize) < nb_segments {
            self.i_lin += 1;

            bof_seg.set_void();
            bof_seg.add_point(TP::begin_of_seg(the_polyg, self.i_lin as usize));
            bof_seg.add_point(TP::end_of_seg(the_polyg, self.i_lin as usize));
            bof_seg.enlarge(TP::deflection_over_estimation(the_polyg));

            //  Modified by MKK - Thu Oct  25 12:40:11 2007
            let def_ph = TH::deflection_over_estimation(the_polyh);
            let maliste: Vec<i32> = polyh_grid.compare(&bof_seg).clone();
            let mut p1 = DVec3::ZERO;
            let mut beg0 = DVec3::ZERO;
            let mut p2 = DVec3::ZERO;
            let mut end0 = DVec3::ZERO;
            if !maliste.is_empty() {
                p1 = TP::begin_of_seg(the_polyg, self.i_lin as usize);
                p2 = TP::end_of_seg(the_polyg, self.i_lin as usize);
                beg0 = p1;
                end0 = p2;
            }
            for &ind_tri in &maliste {
                let p_tri = TH::triangle(the_polyh, ind_tri as usize);
                // Vecteur normal; distance polaire.
                let (tri_nor, tri_dp) = plane_equation(
                    TH::point(the_polyh, p_tri[0]),
                    TH::point(the_polyh, p_tri[1]),
                    TH::point(the_polyh, p_tri[2]),
                );

                // enlarge boundary segment
                if self.i_lin == 1 {
                    let mut dif = p1 - p2;
                    let dist = dif.length();
                    if dist > GP_RESOLUTION {
                        dif /= dist;
                        let mut a_cos = dif.dot(tri_nor);
                        a_cos = a_cos.abs();
                        if a_cos > GP_RESOLUTION {
                            let shift = def_ph / a_cos;
                            beg0 = p1 + dif * shift;
                        }
                    }
                } else if self.i_lin as usize == nb_segments {
                    let mut dif = p2 - p1;
                    let dist = dif.length();
                    if dist > GP_RESOLUTION {
                        dif /= dist;
                        let mut a_cos = dif.dot(tri_nor);
                        a_cos = a_cos.abs();
                        if a_cos > GP_RESOLUTION {
                            let shift = def_ph / a_cos;
                            end0 = p2 + dif * shift;
                        }
                    }
                }
                let d_beg_tri = tri_nor.dot(beg0) - tri_dp; // Distance <p1> plane
                let d_end_tri = tri_nor.dot(end0) - tri_dp; // Distance <p2> plane

                self.intersect_with_normal::<H, TH>(
                    beg0,
                    end0,
                    false,
                    ind_tri,
                    the_polyh,
                    tri_nor,
                    d_beg_tri,
                    d_end_tri,
                );
            }
            self.begin_of_closed_polygon = false;
        }
    }
}

impl InterferencePolygonPolyhedron {
    /// OCCT Extrema_ExtElC(theLin1, theLin2, 0.00000001) for the line x line
    /// case: None = IsParallel(), Some((SquareDistance, U1, U2)) otherwise
    /// (kernel line_line_extrema keeps the OCCT parametrization).
    fn extrema_line_line(
        lin_pol: &Line3,
        lin_tri: &Line3,
    ) -> Option<(f64, f64, f64)> {
        let a_d1 = lin_pol.direction.normalize_or_zero();
        let a_d2 = lin_tri.direction.normalize_or_zero();
        let a_cos_a = a_d1.dot(a_d2);
        let a_sq_sin_a = 1.0 - a_cos_a * a_cos_a;
        if a_sq_sin_a < 1e-30 || a_d1.cross(a_d2).length() < 1e-12 {
            None
        } else {
            let mut v = rcad_kernel::base::extrema::line_line_extrema(lin_pol, lin_tri);
            if v.is_empty() {
                None
            } else {
                Some(v.remove(0))
            }
        }
    }

    /// OCCT Intersect(BegO, EndO, Infinite, TTri, thePolyh) — the active
    /// #else branch (gxx L741-1053).  Computes the triangle plane equation
    /// itself, then classifies the segment/line x triangle intersection.
    #[allow(clippy::too_many_arguments)]
    fn intersect<H, TH>(
        &mut self,
        beg_o: DVec3,
        end_o: DVec3,
        infinite: bool,
        ttri: i32,
        the_polyh: &H,
    ) where
        TH: ToolPolyh<H>,
    {
        let mut typ_on_g = IntfPIType::Edge;
        let p_tri = TH::triangle(the_polyh, ttri as usize);

        let (tri_nor, tri_dp) = plane_equation(
            TH::point(the_polyh, p_tri[0]),
            TH::point(the_polyh, p_tri[1]),
            TH::point(the_polyh, p_tri[2]),
        );

        let d_beg_tri = tri_nor.dot(beg_o) - tri_dp; // Distance <BegO> plan
        let d_end_tri = tri_nor.dot(end_o) - tri_dp; // Distance <EndO> plan
        let mut no_intersection_with_triangle = false;

        let mut param;
        let t = d_beg_tri - d_end_tri;
        if t >= 1.0e-16 || t <= -1.0e-16 {
            param = d_beg_tri / t;
        } else {
            param = d_beg_tri;
        }
        let floatgap = epsilon_1000();

        if !infinite {
            if d_beg_tri <= floatgap && d_beg_tri >= -floatgap {
                param = 0.0;
                typ_on_g = IntfPIType::Vertex;
                if self.begin_of_closed_polygon {
                    no_intersection_with_triangle = false;
                }
            } else if d_end_tri <= floatgap && d_end_tri >= -floatgap {
                param = 1.0;
                typ_on_g = IntfPIType::Vertex;
                no_intersection_with_triangle = false;
            }
            if param < 0.0 || param > 1.0 {
                no_intersection_with_triangle = true;
            }
        }
        if !no_intersection_with_triangle {
            let sp_lieu = beg_o + (end_o - beg_o) * param;
            let mut d_pi_e = [0.0f64; 3];
            let mut d_pt_pi = [0.0f64; 3];
            let mut sigd;
            let mut is = 0usize;
            let mut s_edge: i32 = -1;
            let mut s_vertex: i32 = -1;
            let mut tbreak = 0;
            { //-- is = 0
                let seg_t = TH::point(the_polyh, p_tri[1]) - TH::point(the_polyh, p_tri[0]);
                let vec_p = sp_lieu - TH::point(the_polyh, p_tri[0]);
                d_pt_pi[0] = vec_p.length();
                if d_pt_pi[0] <= floatgap {
                    s_vertex = 0;
                    is = 0;
                    tbreak = 1;
                } else {
                    let seg_t_x_vec_p = seg_t.cross(vec_p);
                    let modulus_seg_t_x_vec_p = seg_t_x_vec_p.length();
                    sigd = seg_t_x_vec_p.dot(tri_nor);
                    if sigd > floatgap {
                        sigd = 1.0;
                    } else if sigd < -floatgap {
                        sigd = -1.0;
                    } else {
                        sigd = 0.0;
                    }
                    d_pi_e[0] = sigd * (modulus_seg_t_x_vec_p / seg_t.length());
                    if d_pi_e[0] <= floatgap && d_pi_e[0] >= -floatgap {
                        s_edge = 0;
                        is = 0;
                        tbreak = 1;
                    }
                }
            }

            if tbreak == 0 {
                { //-- is = 1
                    let seg_t = TH::point(the_polyh, p_tri[2]) - TH::point(the_polyh, p_tri[1]);
                    let vec_p = sp_lieu - TH::point(the_polyh, p_tri[1]);
                    d_pt_pi[1] = vec_p.length();
                    if d_pt_pi[1] <= floatgap {
                        s_vertex = 1;
                        is = 1;
                        tbreak = 1;
                    } else {
                        let seg_t_x_vec_p = seg_t.cross(vec_p);
                        let modulus_seg_t_x_vec_p = seg_t_x_vec_p.length();
                        sigd = seg_t_x_vec_p.dot(tri_nor);
                        if sigd > floatgap {
                            sigd = 1.0;
                        } else if sigd < -floatgap {
                            sigd = -1.0;
                        } else {
                            sigd = 0.0;
                        }
                        d_pi_e[1] = sigd * (modulus_seg_t_x_vec_p / seg_t.length());
                        if d_pi_e[1] <= floatgap && d_pi_e[1] >= -floatgap {
                            s_edge = 1;
                            is = 1;
                            tbreak = 1;
                        }
                    }
                }
            }
            if tbreak == 0 {
                { //-- is = 2
                    let seg_t = TH::point(the_polyh, p_tri[0]) - TH::point(the_polyh, p_tri[2]);
                    let vec_p = sp_lieu - TH::point(the_polyh, p_tri[2]);
                    d_pt_pi[2] = vec_p.length();
                    if d_pt_pi[2] <= floatgap {
                        s_vertex = 2;
                        is = 2;
                    }
                    let seg_t_x_vec_p = seg_t.cross(vec_p);
                    let modulus_seg_t_x_vec_p = seg_t_x_vec_p.length();
                    sigd = seg_t_x_vec_p.dot(tri_nor);
                    if sigd > floatgap {
                        sigd = 1.0;
                    } else if sigd < -floatgap {
                        sigd = -1.0;
                    } else {
                        sigd = 0.0;
                    }
                    d_pi_e[2] = sigd * (modulus_seg_t_x_vec_p / seg_t.length());
                    if d_pi_e[2] <= floatgap && d_pi_e[2] >= -floatgap {
                        s_edge = 2;
                        is = 2;
                    }
                }
            }
            //-- fin for i=0 to 2

            if s_vertex > -1 {
                // triCon/pedg kept for the commented-out while loop of the
                // OCCT source (gxx L1220-1224); they are not consumed.
                let _tri_con = ttri;
                let _pedg = p_tri[POURCENT3[s_vertex as usize + 1]];
                let sp = IntfSectionPoint::new(
                    sp_lieu,
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Vertex,
                    p_tri[is] as i32,
                    0,
                    0.0,
                    1.0,
                );
                self.interf.my_s_poins_mut().push(sp);
            } else if s_edge > -1 {
                let (_tri_con, _pedg) = TH::tri_connex(
                    the_polyh,
                    ttri as usize,
                    p_tri[s_edge as usize],
                    p_tri[POURCENT3[s_edge as usize + 1]],
                );
                let sp = IntfSectionPoint::new(
                    sp_lieu,
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Edge,
                    p_tri[s_edge as usize]
                        .min(p_tri[POURCENT3[s_edge as usize + 1]])
                        as i32,
                    p_tri[s_edge as usize]
                        .max(p_tri[POURCENT3[s_edge as usize + 1]])
                        as i32,
                    0.0,
                    1.0,
                );
                self.interf.my_s_poins_mut().push(sp);
            } else if d_pi_e[0] > 0.0 && d_pi_e[1] > 0.0 && d_pi_e[2] > 0.0 {
                let sp = IntfSectionPoint::new(
                    sp_lieu,
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Face,
                    ttri,
                    0,
                    0.0,
                    1.0,
                );
                self.interf.my_s_poins_mut().push(sp);
            }
            //  Modified by Sergey KHROMOV - Fri Dec  7 14:40:11 2001
            // Sometimes triangulation doesn't cover whole the face. In this
            // case it is necessary to take into account the deflection
            // between boundary isolines of the surface and boundary
            // triangles.
            else {
                for i in 1..=3usize {
                    let ind_p1 = if i == 3 { p_tri[0] } else { p_tri[i] };
                    let ind_p2 = p_tri[i - 1];

                    if TH::is_on_bound(the_polyh, ind_p1, ind_p2) {
                        // For boundary line it is necessary to check the
                        // border deflection.
                        let deflection = TH::get_border_deflection(the_polyh);
                        let beg_p = TH::point(the_polyh, ind_p1);
                        let end_p = TH::point(the_polyh, ind_p2);
                        let vec_tri = end_p - beg_p;
                        let lin_tri = Line3::new(beg_p, vec_tri);
                        let a_p_on_e = sp_lieu;
                        let a_dist = lin_tri.distance(a_p_on_e);

                        if a_dist <= deflection {
                            let a_v_loc_p_on_e = a_p_on_e - beg_p;
                            let a_vec_dir_tri = lin_tri.direction;
                            let a_par = a_v_loc_p_on_e.dot(a_vec_dir_tri);
                            let a_max_par = vec_tri.length();

                            if a_par >= 0.0 && a_par <= a_max_par {
                                let sp = IntfSectionPoint::new(
                                    sp_lieu,
                                    typ_on_g,
                                    0,
                                    self.i_lin,
                                    param,
                                    IntfPIType::Face,
                                    ttri,
                                    0,
                                    0.0,
                                    1.0,
                                );
                                self.interf.my_s_poins_mut().push(sp);
                            }
                        }
                    }
                }
            }
            //  Modified by Sergey KHROMOV - Fri Dec  7 14:40:29 2001 End
        } //---- if(NoIntersectionWithTriangle == false)

        //-------------------------------------------------------------------
        //-- On teste la distance entre les cotes du triangle et le polygone
        //--
        //-- Si cette distance est inferieure a Tolerance, on cree un SP.
        //-------------------------------------------------------------------
        {
            let vec_pol = end_o - beg_o;
            let n_vec_pol = vec_pol.length();
            if n_vec_pol < PRECISION_COMPUTATIONAL {
                return;
            }
            let lin_pol = Line3::new(beg_o, vec_pol);

            for i in 0..3usize {
                let p_tri_ip1pc3 = p_tri[POURCENT3[i + 1]];
                let p_tri_i = p_tri[i];
                let beg_t = TH::point(the_polyh, p_tri_ip1pc3);
                let end_t = TH::point(the_polyh, p_tri_i);
                let vec_tri = end_t - beg_t;
                let n_vec_tri = vec_tri.length();
                if n_vec_tri < PRECISION_COMPUTATIONAL {
                    continue;
                }
                let lin_tri = Line3::new(beg_t, vec_tri);
                // OCCT Extrema_ExtElC Extrema(LinPol, LinTri, 0.00000001).
                if let Some((dist2, u_o, u_t)) = Self::extrema_line_line(&lin_pol, &lin_tri) {
                    // Extrema.IsParallel() == false && Extrema.NbExt() != 0.
                    if dist2 <= self.interf.get_tolerance() * self.interf.get_tolerance() {
                        let po = lin_pol.origin + lin_pol.direction * u_o;
                        let pt = lin_tri.origin + lin_tri.direction * u_t;
                        let mut param_on_o = 0.0;
                        if is_in_segment(vec_pol, po - beg_o, n_vec_pol, &mut param_on_o, self.interf.get_tolerance())
                        {
                            let mut param_on_t = 0.0;
                            if is_in_segment(vec_tri, pt - beg_t, n_vec_tri, &mut param_on_t, self.interf.get_tolerance())
                            {
                                let sp_lieu = beg_t + (end_t - beg_t) * param_on_t;
                                let (tmin, tmax) = if p_tri_i > p_tri_ip1pc3 {
                                    (p_tri_ip1pc3, p_tri_i)
                                } else {
                                    (p_tri_i, p_tri_ip1pc3)
                                };
                                let sp = IntfSectionPoint::new(
                                    sp_lieu,
                                    typ_on_g,
                                    0,
                                    self.i_lin,
                                    param_on_o,
                                    IntfPIType::Edge,
                                    tmin as i32,
                                    tmax as i32,
                                    0.0,
                                    1.0,
                                );
                                self.interf.my_s_poins_mut().push(sp);
                            }
                        }
                    }
                }
            }
        }
    }

    /// OCCT Intersect(BegO, EndO, Infinite, TTri, thePolyh, TriNormal,
    /// TriDp, dBegTri, dEndTri) — the 9-argument variant (gxx L1057-1368)
    /// reusing the caller's plane equation (TriDp is unused there too).
    #[allow(clippy::too_many_arguments)]
    fn intersect_with_normal<H, TH>(
        &mut self,
        beg_o: DVec3,
        end_o: DVec3,
        infinite: bool,
        ttri: i32,
        the_polyh: &H,
        tri_normal: DVec3,
        d_beg_tri: f64,
        d_end_tri: f64,
    ) where
        TH: ToolPolyh<H>,
    {
        let mut typ_on_g = IntfPIType::Edge;
        let p_tri = TH::triangle(the_polyh, ttri as usize);
        let tri_nor = tri_normal; // Vecteur normal.

        let mut no_intersection_with_triangle = false;

        let mut param;
        let t = d_beg_tri - d_end_tri;
        if t >= 1.0e-16 || t <= -1.0e-16 {
            param = d_beg_tri / t;
        } else {
            param = d_beg_tri;
        }
        let floatgap = epsilon_1000();

        if !infinite {
            if d_beg_tri <= floatgap && d_beg_tri >= -floatgap {
                param = 0.0;
                typ_on_g = IntfPIType::Vertex;
                if self.begin_of_closed_polygon {
                    no_intersection_with_triangle = false;
                }
            } else if d_end_tri <= floatgap && d_end_tri >= -floatgap {
                param = 1.0;
                typ_on_g = IntfPIType::Vertex;
                no_intersection_with_triangle = false;
            }
            if param < 0.0 || param > 1.0 {
                no_intersection_with_triangle = true;
            }
        }
        if !no_intersection_with_triangle {
            let sp_lieu = beg_o + (end_o - beg_o) * param;
            let mut d_pi_e = [0.0f64; 3];
            let mut d_pt_pi = [0.0f64; 3];
            let mut sigd;
            let mut is = 0usize;
            let mut s_edge: i32 = -1;
            let mut s_vertex: i32 = -1;
            let mut tbreak = 0;
            { //-- is = 0
                let seg_t = TH::point(the_polyh, p_tri[1]) - TH::point(the_polyh, p_tri[0]);
                let vec_p = sp_lieu - TH::point(the_polyh, p_tri[0]);
                d_pt_pi[0] = vec_p.length();
                if d_pt_pi[0] <= floatgap {
                    s_vertex = 0;
                    is = 0;
                    tbreak = 1;
                } else {
                    let seg_t_x_vec_p = seg_t.cross(vec_p);
                    let modulus_seg_t_x_vec_p = seg_t_x_vec_p.length();
                    sigd = seg_t_x_vec_p.dot(tri_nor);
                    if sigd > floatgap {
                        sigd = 1.0;
                    } else if sigd < -floatgap {
                        sigd = -1.0;
                    } else {
                        sigd = 0.0;
                    }
                    d_pi_e[0] = sigd * (modulus_seg_t_x_vec_p / seg_t.length());
                    if d_pi_e[0] <= floatgap && d_pi_e[0] >= -floatgap {
                        s_edge = 0;
                        is = 0;
                        tbreak = 1;
                    }
                }
            }

            if tbreak == 0 {
                { //-- is = 1
                    let seg_t = TH::point(the_polyh, p_tri[2]) - TH::point(the_polyh, p_tri[1]);
                    let vec_p = sp_lieu - TH::point(the_polyh, p_tri[1]);
                    d_pt_pi[1] = vec_p.length();
                    if d_pt_pi[1] <= floatgap {
                        s_vertex = 1;
                        is = 1;
                        tbreak = 1;
                    } else {
                        let seg_t_x_vec_p = seg_t.cross(vec_p);
                        let modulus_seg_t_x_vec_p = seg_t_x_vec_p.length();
                        sigd = seg_t_x_vec_p.dot(tri_nor);
                        if sigd > floatgap {
                            sigd = 1.0;
                        } else if sigd < -floatgap {
                            sigd = -1.0;
                        } else {
                            sigd = 0.0;
                        }
                        d_pi_e[1] = sigd * (modulus_seg_t_x_vec_p / seg_t.length());
                        if d_pi_e[1] <= floatgap && d_pi_e[1] >= -floatgap {
                            s_edge = 1;
                            is = 1;
                            tbreak = 1;
                        }
                    }
                }
            }
            if tbreak == 0 {
                { //-- is = 2
                    let seg_t = TH::point(the_polyh, p_tri[0]) - TH::point(the_polyh, p_tri[2]);
                    let vec_p = sp_lieu - TH::point(the_polyh, p_tri[2]);
                    d_pt_pi[2] = vec_p.length();
                    if d_pt_pi[2] <= floatgap {
                        s_vertex = 2;
                        is = 2;
                    }
                    let seg_t_x_vec_p = seg_t.cross(vec_p);
                    let modulus_seg_t_x_vec_p = seg_t_x_vec_p.length();
                    sigd = seg_t_x_vec_p.dot(tri_nor);
                    if sigd > floatgap {
                        sigd = 1.0;
                    } else if sigd < -floatgap {
                        sigd = -1.0;
                    } else {
                        sigd = 0.0;
                    }
                    d_pi_e[2] = sigd * (modulus_seg_t_x_vec_p / seg_t.length());
                    if d_pi_e[2] <= floatgap && d_pi_e[2] >= -floatgap {
                        s_edge = 2;
                        is = 2;
                    }
                }
            }
            //-- fin for i=0 to 2

            if s_vertex > -1 {
                let _tri_con = ttri;
                let _pedg = p_tri[POURCENT3[s_vertex as usize + 1]];
                let sp = IntfSectionPoint::new(
                    sp_lieu,
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Vertex,
                    p_tri[is] as i32,
                    0,
                    0.0,
                    1.0,
                );
                self.interf.my_s_poins_mut().push(sp);
            } else if s_edge > -1 {
                let (_tri_con, _pedg) = TH::tri_connex(
                    the_polyh,
                    ttri as usize,
                    p_tri[s_edge as usize],
                    p_tri[POURCENT3[s_edge as usize + 1]],
                );
                let sp = IntfSectionPoint::new(
                    sp_lieu,
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Edge,
                    p_tri[s_edge as usize]
                        .min(p_tri[POURCENT3[s_edge as usize + 1]])
                        as i32,
                    p_tri[s_edge as usize]
                        .max(p_tri[POURCENT3[s_edge as usize + 1]])
                        as i32,
                    0.0,
                    1.0,
                );
                self.interf.my_s_poins_mut().push(sp);
            } else if d_pi_e[0] > 0.0 && d_pi_e[1] > 0.0 && d_pi_e[2] > 0.0 {
                let sp = IntfSectionPoint::new(
                    sp_lieu,
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Face,
                    ttri,
                    0,
                    0.0,
                    1.0,
                );
                self.interf.my_s_poins_mut().push(sp);
            }
            //  Modified by Sergey KHROMOV - Fri Dec  7 14:40:11 2001
            else {
                for i in 1..=3usize {
                    let ind_p1 = if i == 3 { p_tri[0] } else { p_tri[i] };
                    let ind_p2 = p_tri[i - 1];

                    if TH::is_on_bound(the_polyh, ind_p1, ind_p2) {
                        let deflection = TH::get_border_deflection(the_polyh);
                        let beg_p = TH::point(the_polyh, ind_p1);
                        let end_p = TH::point(the_polyh, ind_p2);
                        let vec_tri = end_p - beg_p;
                        let lin_tri = Line3::new(beg_p, vec_tri);
                        let a_p_on_e = sp_lieu;
                        let a_dist = lin_tri.distance(a_p_on_e);

                        if a_dist <= deflection {
                            let a_v_loc_p_on_e = a_p_on_e - beg_p;
                            let a_vec_dir_tri = lin_tri.direction;
                            let a_par = a_v_loc_p_on_e.dot(a_vec_dir_tri);
                            let a_max_par = vec_tri.length();

                            if a_par >= 0.0 && a_par <= a_max_par {
                                let sp = IntfSectionPoint::new(
                                    sp_lieu,
                                    typ_on_g,
                                    0,
                                    self.i_lin,
                                    param,
                                    IntfPIType::Face,
                                    ttri,
                                    0,
                                    0.0,
                                    1.0,
                                );
                                self.interf.my_s_poins_mut().push(sp);
                            }
                        }
                    }
                }
            }
            //  Modified by Sergey KHROMOV - Fri Dec  7 14:40:29 2001 End
        } //---- if(NoIntersectionWithTriangle == false)

        //-------------------------------------------------------------------
        //-- On teste la distance entre les cotes du triangle et le polygone
        //--
        //-- Si cette distance est inferieure a Tolerance, on cree un SP.
        //-------------------------------------------------------------------
        {
            let vec_pol = end_o - beg_o;
            let n_vec_pol = vec_pol.length();
            if n_vec_pol < PRECISION_COMPUTATIONAL {
                return;
            }
            let lin_pol = Line3::new(beg_o, vec_pol);

            for i in 0..3usize {
                let p_tri_ip1pc3 = p_tri[POURCENT3[i + 1]];
                let p_tri_i = p_tri[i];
                let beg_t = TH::point(the_polyh, p_tri_ip1pc3);
                let end_t = TH::point(the_polyh, p_tri_i);
                let vec_tri = end_t - beg_t;
                let n_vec_tri = vec_tri.length();
                if n_vec_tri < PRECISION_COMPUTATIONAL {
                    continue;
                }
                let lin_tri = Line3::new(beg_t, vec_tri);
                // OCCT Extrema_ExtElC Extrema(LinPol, LinTri, 0.00000001).
                if let Some((dist2, u_o, u_t)) = Self::extrema_line_line(&lin_pol, &lin_tri) {
                    // Extrema.IsParallel() == false && Extrema.NbExt() != 0.
                    if dist2 <= self.interf.get_tolerance() * self.interf.get_tolerance() {
                        let po = lin_pol.origin + lin_pol.direction * u_o;
                        let pt = lin_tri.origin + lin_tri.direction * u_t;
                        let mut param_on_o = 0.0;
                        if is_in_segment(vec_pol, po - beg_o, n_vec_pol, &mut param_on_o, self.interf.get_tolerance())
                        {
                            let mut param_on_t = 0.0;
                            if is_in_segment(vec_tri, pt - beg_t, n_vec_tri, &mut param_on_t, self.interf.get_tolerance())
                            {
                                let sp_lieu = beg_t + (end_t - beg_t) * param_on_t;
                                let (tmin, tmax) = if p_tri_i > p_tri_ip1pc3 {
                                    (p_tri_ip1pc3, p_tri_i)
                                } else {
                                    (p_tri_i, p_tri_ip1pc3)
                                };
                                let sp = IntfSectionPoint::new(
                                    sp_lieu,
                                    typ_on_g,
                                    0,
                                    self.i_lin,
                                    param_on_o,
                                    IntfPIType::Edge,
                                    tmin as i32,
                                    tmax as i32,
                                    0.0,
                                    1.0,
                                );
                                self.interf.my_s_poins_mut().push(sp);
                            }
                        }
                    }
                }
            }
        }
    }
}

impl Default for InterferencePolygonPolyhedron {
    /// OCCT Intf_InterferencePolygonPolyhedron() (gxx L60-65).
    fn default() -> Self {
        InterferencePolygonPolyhedron::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::int_curv_surf::{
        TheInterferenceOfHInter, ThePolyhedronOfHInter, ThePolygonOfHInter,
        ThePolyhedronToolOfHInter, ThePolygonToolOfHInter,
    };
    use rcad_kernel::geom::{CurveEval, Line3, Plane};

    /// The z=0 plane sampled as a 3x3 polyhedron over [0,10]x[0,10]
    /// (OCCT IntCurveSurface_ThePolyhedronOfHInter over a plane).
    fn plane_polyhedron() -> ThePolyhedronOfHInter {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        ThePolyhedronOfHInter::new(&plane, 3, 3, 0.0, 0.0, 10.0, 10.0)
    }

    /// OCCT Intf_InterferencePolygonPolyhedron(theLin, thePolyh)
    /// (gxx L114-147) via the HInter instantiation: a vertical line through
    /// the plane polyhedron gives one face section point at the crossing.
    #[test]
    fn line_through_plane_polyhedron() {
        let polyh = plane_polyhedron();
        let lin = Line3::new(DVec3::new(5.3, 4.7, -1.0), DVec3::Z);
        let mut itf =
            TheInterferenceOfHInter::new_line_polyhedron::<_, ThePolyhedronToolOfHInter>(&lin, &polyh);
        assert!(itf.interf.nb_section_points() >= 1, "expected a face crossing");
        // Exactly one triangle contains (5.3, 4.7); the neighbours yield no
        // point (their dPiE classification rejects them).
        assert_eq!(itf.interf.nb_section_points(), 1);
        let sp = itf.interf.pnt_value(1);
        let p = sp.pnt();
        assert!((p.x - 5.3).abs() < 1.0e-9 && (p.y - 4.7).abs() < 1.0e-9 && p.z.abs() < 1.0e-9, "got {p:?}");
        // The tool side records the crossed triangle and the FACE kind.
        let (tool_kind, tool_t1, _, _) = sp.info_second_full();
        assert_eq!(tool_kind, IntfPIType::Face);
        assert!(tool_t1 >= 1, "1-based triangle address, got {}", tool_t1);
        // OCCT quirk: with Infinite=true the object-side kind stays EDGE.
        let (obj_kind, _, obj_lin, obj_param) = sp.info_first_full();
        assert_eq!(obj_kind, IntfPIType::Edge);
        assert_eq!(obj_lin, 0, "iLin stays 0 on the line ctors");
        assert!((obj_param - 1.0).abs() < 1.0e-9);
    }

    /// OCCT Interference(thePolyg, thePolyh) (gxx L296-350): a polygon
    /// segment crossing the plane gives the +/-normal offset pair of
    /// section points at the crossing.
    #[test]
    fn polygon_segment_through_plane_polyhedron() {
        let polyh = plane_polyhedron();
        let seg = Line3::new(DVec3::new(5.3, 4.7, -1.0), DVec3::Z);
        let polyg = ThePolygonOfHInter::new_params(&seg, &[0.0, 2.0]);
        assert_eq!(polyg.nb_segments(), 1);
        let mut itf = TheInterferenceOfHInter::new_polygon_polyhedron::<
            _,
            _,
            ThePolygonToolOfHInter,
            ThePolyhedronToolOfHInter,
        >(&polyg, &polyh);
        let n = itf.interf.nb_section_points();
        assert!(n >= 2, "the +/-normal offset pair both cross: got {n}");
        // Both offset lines pass through the same plane point.
        for i in 1..=n {
            let p = itf.interf.pnt_value(i).pnt();
            assert!((p.x - 5.3).abs() < 1.0e-6 && (p.y - 4.7).abs() < 1.0e-6 && p.z.abs() < 1.0e-6, "got {p:?}");
        }
        assert_eq!(itf.interf.nb_section_lines(), 0);
        assert_eq!(itf.interf.nb_tangent_zones(), 0);
        assert!(itf.interf.get_tolerance() > 0.0);
    }

    /// OCCT Perform(theLin, thePolyh) reuse (gxx L214-245): a second
    /// perform resets the result (SelfInterference).
    #[test]
    fn perform_resets_state() {
        let polyh = plane_polyhedron();
        let lin = Line3::new(DVec3::new(5.3, 4.7, -1.0), DVec3::Z);
        let mut itf =
            TheInterferenceOfHInter::new_line_polyhedron::<_, ThePolyhedronToolOfHInter>(&lin, &polyh);
        assert_eq!(itf.interf.nb_section_points(), 1);
        let miss = Line3::new(DVec3::new(50.0, 50.0, -1.0), DVec3::Z);
        itf.perform_line::<_, ThePolyhedronToolOfHInter>(&miss, &polyh);
        assert_eq!(itf.interf.nb_section_points(), 0, "previous points cleared");
    }

    /// OCCT Intf::PlaneEquation (Intf.cxx L20-42): unit normal and polar
    /// distance d = N . P1.
    #[test]
    fn plane_equation_anchor() {
        let (n, d) = plane_equation(DVec3::ZERO, DVec3::X, DVec3::Y);
        assert!((n - DVec3::Z).length() < 1.0e-12 || (n + DVec3::Z).length() < 1.0e-12);
        assert!(d.abs() < 1.0e-12);
        let (n, d) = plane_equation(DVec3::new(0.0, 0.0, 5.0), DVec3::new(1.0, 0.0, 5.0), DVec3::new(0.0, 1.0, 5.0));
        assert!((n.length() - 1.0).abs() < 1.0e-12);
        assert!((d - 5.0).abs() < 1.0e-12 || (d + 5.0).abs() < 1.0e-12);
    }

    /// OCCT gxx L36-56 IsInSegment: projection with clamping.
    #[test]
    fn is_in_segment_anchor() {
        let mut param = 0.0;
        // p at half of a unit segment.
        assert!(is_in_segment(DVec3::X, DVec3::new(0.5, 0.0, 0.0), 1.0, &mut param, 1.0e-6));
        assert!((param - 0.5).abs() < 1.0e-12);
        // p beyond the end: rejected.
        assert!(!is_in_segment(DVec3::X, DVec3::new(1.5, 0.0, 0.0), 1.0, &mut param, 1.0e-6));
        // Degenerate segment: rejected.
        assert!(!is_in_segment(DVec3::ZERO, DVec3::new(0.5, 0.0, 0.0), 0.0, &mut param, 1.0e-6));
    }
}
