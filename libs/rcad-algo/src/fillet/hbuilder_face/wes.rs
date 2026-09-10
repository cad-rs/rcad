//! OCCT TopOpeBRepBuild_WireEdgeSet (WireEdgeSet.cxx L51-575 + .hxx) —
//! continuation file of `fillet::hbuilder_face` (the SplitFace1 D6
//! carrier): a bound is a wire, a boundelement is an edge.  The ShapeSet
//! stores a list of wires (bounds), a list of edges (boundelements) to
//! start reconstructions, and a map of vertex giving the list of edges
//! incident to a vertex (WireEdgeSet.hxx L36-40).

#![allow(dead_code)]

use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _, SurfaceEval as _};
use rcad_kernel::topods::{BRep, BRepTool as _, Orientation, Shape, ShapeType, TShape};

use super::classify::{
    brep_tool_parameter, fc2d_curve_on_surface, face_surface_value, shape_explore_children,
    shape_key, shapes_same, ShapeSet, ShapeSetBase,
};

// =========================================================================
// OCCT TopOpeBRepBuild_WireEdgeSet (WireEdgeSet.cxx L51-575 + .hxx).
// =========================================================================

/// OCCT TopOpeBRepBuild_WireEdgeSet — a bound is a wire, a boundelement
/// is an edge.  The ShapeSet stores a list of wires (bounds), a list of
/// edges (boundelements) to start reconstructions, and a map of vertex
/// giving the list of edges incident to a vertex (WireEdgeSet.hxx L36-40).
pub(crate) struct WireEdgeSet {
    pub(crate) base: ShapeSetBase,
    pub(crate) my_face: Shape,
}

impl WireEdgeSet {
    /// WireEdgeSet.cxx L51-55.
    pub(crate) fn new(f: &Shape) -> Self {
        WireEdgeSet {
            base: ShapeSetBase::new(ShapeType::Vertex),
            my_face: f.clone(),
        }
    }

    /// WireEdgeSet.cxx L124-127.
    pub(crate) fn face(&self) -> &Shape {
        &self.my_face
    }

    /// WireEdgeSet.cxx L292-318.
    pub(crate) fn vertex_connects_edges(
        &self,
        brep: &BRep,
        v: &Shape,
        e1: &Shape,
        e2: &Shape,
        o1: &mut Orientation,
        o2: &mut Orientation,
    ) -> bool {
        let ex1 = shape_explore_children(brep, e1, ShapeType::Vertex);
        for c1 in &ex1 {
            if shapes_same(v, c1) {
                let ex2 = shape_explore_children(brep, e2, ShapeType::Vertex);
                for c2 in &ex2 {
                    if shapes_same(v, c2) {
                        *o1 = c1.orientation;
                        *o2 = c2.orientation;
                        if o1 != o2 {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    /// WireEdgeSet.cxx L322-388.
    pub(crate) fn vertex_connects_edges_closing(
        &self,
        brep: &BRep,
        v: &Shape,
        e1: &Shape,
        e2: &Shape,
    ) -> bool {
        // bool VertexConnectsEdgesClosing (WireEdgeSet.cxx L327-346):
        // if E1 and E2 are not closed : edges are NOT connected
        // if E1 or E2 is/are closed :
        //   if V changes of relative orientation between E1,E2 : connected
        //   else : NOT connected
        let c1 = self.is_closed(brep, e1);
        let c2 = self.is_closed(brep, e2);

        let mut testconnect = c1 || c2;
        let resu;
        let mut o1 = Orientation::Forward;
        let mut o2 = Orientation::Forward;

        if c1 && c2 {
            let u1 = if c1 { self.is_u_closed(brep, e1) } else { false };
            let v1 = if c1 { self.is_v_closed(brep, e1) } else { false };
            let u2 = if c2 { self.is_u_closed(brep, e2) } else { false };
            let v2 = if c2 { self.is_v_closed(brep, e2) } else { false };
            let uvdiff = (u1 && v2) || (u2 && v1);
            testconnect = uvdiff;
        }

        if testconnect {
            resu = self.vertex_connects_edges(brep, v, e1, e2, &mut o1, &mut o2);
        } else {
            // cto 012 O2 arete de couture de face cylindrique
            let oe1 = e1.orientation;
            let oe2 = e2.orientation;
            // OCCT: bool iseq = E1.IsEqual(E2) — IsSame + same orientation.
            let iseq = shapes_same(e1, e2) && e1.orientation == e2.orientation;
            if (c1 && c2) && (oe1 == oe2) && (!iseq) {
                resu = self.vertex_connects_edges(brep, v, e1, e2, &mut o1, &mut o2);
            } else {
                resu = false;
            }
        }
        resu
    }

    /// WireEdgeSet.cxx L392-404.
    pub(crate) fn nb_closing_shapes(&self, brep: &BRep, l: &[Shape]) -> i32 {
        let mut n = 0;
        for s in l {
            if self.is_closed(brep, s) {
                n += 1;
            }
        }
        n
    }

    /// WireEdgeSet.cxx L408-439.
    pub(crate) fn local_d1(
        &self,
        brep: &BRep,
        sf: &Shape,
        se: &Shape,
        sv: &Shape,
        pe: &mut glam::DVec2,
        d1e: &mut glam::DVec2,
    ) {
        let f = sf;
        let e = se;
        let v = sv;
        let par_e = brep_tool_parameter(brep, v, e);

        // OCCT: CE = BRep_Tool::Curve(E, Loc, fiE, laE);
        // CE = CE->Transformed(Loc.Transformation()) — the world curve.
        let Some((ce, _range)) = brep.edge_curve_world(e) else {
            *d1e = glam::DVec2::ZERO;
            return;
        };
        let p3de;
        let d3de;
        // CE->D1(parE, p3dE, d3dE).
        p3de = ce.point_at(par_e);
        d3de = ce.derivative_at(par_e);

        // S = BRep_Tool::Surface(F); proj(p3dE, S); u, v =
        // proj.LowerDistanceParameters().
        let Some(s) = face_surface_value(brep, f) else {
            *d1e = glam::DVec2::ZERO;
            return;
        };
        let domain = s.default_domain();
        // OCCT GeomAPI_ProjectPointOnSurf proj(p3dE, S) — the natural
        // surface domain projection.
        let mut proj = crate::bop::int_tools::context::ProjectOnSurface::new_init(
            s.clone(),
            domain,
            rcad_kernel::core::precision::CONFUSION,
        );
        proj.perform(p3de);
        let (mut u, mut v) = proj.lower_distance_parameters();
        pe.x = u;
        pe.y = v;
        // S->D1(u, v, bid, d1u, d1v).
        let (_bid, d1u, d1v) = s.derivatives(u, v);
        u = d3de.dot(d1u);
        v = d3de.dot(d1v);
        d1e.x = u;
        d1e.y = v;
    }

    /// WireEdgeSet.cxx L443-455 — IsClosed(E): BRep_Tool::IsClosed(EE,
    /// myFace).
    pub(crate) fn is_closed(&self, brep: &BRep, e: &Shape) -> bool {
        let ee = e;
        let closed = brep.is_edge_closed_on_face(ee, &self.my_face);
        closed
    }

    /// WireEdgeSet.cxx L459-490 (static IsUVISO).
    pub(crate) fn is_uviso(
        brep: &BRep,
        e: &Shape,
        f: &Shape,
        uiso: &mut bool,
        viso: &mut bool,
    ) {
        *uiso = false;
        *viso = false;
        let Some((pc, _fe, _le, _tolpc)) = fc2d_curve_on_surface(brep, e, f, true) else {
            // WireEdgeSet.cxx L471: throw Standard_ProgramError.
            panic!("TopOpeBRepBuild_WireEdgeSet::IsUVISO");
        };

        // OCCT: if the pcurve is a Geom2d_Line, check the direction
        // parallelism with the Y axis (uiso) / X axis (viso) at
        // Precision::Angular().
        if let rcad_kernel::geom::Curve2d::Line(hl) = &pc {
            let d = hl.direction;
            let tol = rcad_kernel::core::precision::ANGULAR;
            // gp_Dir2d::IsParallel(D, tol) — the angle with the axis is 0
            // or PI; for unit directions |cross| = |sin(angle)| <= tol.
            let par_y = (d.x * 0.0 - d.y * 1.0).abs() <= tol;
            let par_x = (d.x * 1.0 - d.y * 0.0).abs() <= tol;
            if par_y {
                *uiso = true;
            } else if par_x {
                *viso = true;
            }
        }
    }

    /// WireEdgeSet.cxx L494-500.
    pub(crate) fn is_u_closed(&self, brep: &BRep, e: &Shape) -> bool {
        let ee = e;
        let mut bid = false;
        let mut closed = false;
        Self::is_uviso(brep, ee, &self.my_face, &mut closed, &mut bid);
        closed
    }

    /// WireEdgeSet.cxx L504-510.
    pub(crate) fn is_v_closed(&self, brep: &BRep, e: &Shape) -> bool {
        let ee = e;
        let mut bid = false;
        let mut closed = false;
        Self::is_uviso(brep, ee, &self.my_face, &mut bid, &mut closed);
        closed
    }

    // WireEdgeSet.cxx L514-575: the SNameVEE / SNameVEL / DumpSS / SName /
    // SNameori debug-name helpers — the OCCT bodies return empty strings
    // (the debug dump machinery is compiled out); carried as stubs with
    // the same shape.

    /// WireEdgeSet.cxx L514-520.
    pub(crate) fn s_name_vee(&self, _v: &Shape, _e1: &Shape, _e2: &Shape) -> String {
        String::new()
    }

    /// WireEdgeSet.cxx L524-531.
    pub(crate) fn s_name_vel(&self, _v: &Shape, _e: &Shape, _l: &[Shape]) -> String {
        String::new()
    }

    /// WireEdgeSet.cxx L539-545.
    pub(crate) fn s_name(&self, _s: &Shape, sb: &str, _sa: &str) -> String {
        String::from(sb)
    }

    /// WireEdgeSet.cxx L549-555.
    pub(crate) fn s_nameori(&self, _s: &Shape, sb: &str, _sa: &str) -> String {
        String::from(sb)
    }

    /// WireEdgeSet.cxx L559-565.
    pub(crate) fn s_name_list(&self, _s: &[Shape], _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// WireEdgeSet.cxx L569-575.
    pub(crate) fn s_nameori_list(&self, _s: &[Shape], _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// WireEdgeSet.cxx L535.
    pub(crate) fn dump_ss(&self) {}
}

impl ShapeSet for WireEdgeSet {
    /// WireEdgeSet.cxx L59-84.
    fn add_shape(&mut self, brep: &BRep, s: &Shape) {
        let mut tocheck = true;
        let iswire = s.shape_type() == ShapeType::Wire;
        if iswire {
            // OCCT: BRepAdaptor_Surface bas(myFace, false); uc =
            // bas.IsUClosed(); vc = bas.IsVClosed().
            let bas = crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface::initialize_face(
                brep, &self.my_face, false,
            );
            let uc = bas.is_u_closed();
            let vc = bas.is_v_closed();
            if uc || vc {
                tocheck = false;
            }
        }
        let mut chk = true;
        if tocheck {
            chk = self.base.check_shape(s);
        }

        if !chk {
            return;
        }
        self.base.process_add_shape(s);
    }

    /// WireEdgeSet.cxx L88-113.
    fn add_start_element(&mut self, brep: &BRep, s: &Shape) {
        let mut tocheck = true;
        let isedge = s.shape_type() == ShapeType::Edge;
        if isedge {
            // OCCT: BRepAdaptor_Curve cac(TopoDS::Edge(S));
            // t = cac.GetType(); b = (t == GeomAbs_BSplineCurve ||
            // t == GeomAbs_BezierCurve); tocheck = !b.
            let b = match &*s.data {
                TShape::Edge(ed) => {
                    let t = ed.curve.as_ref().map(|c| match c {
                        rcad_kernel::geom::Curve3::BSpline(_) => GeomAbsCurveType::BSplineCurve,
                        rcad_kernel::geom::Curve3::Bezier(_) => GeomAbsCurveType::BezierCurve,
                        _ => GeomAbsCurveType::OtherCurve,
                    });
                    matches!(
                        t,
                        Some(GeomAbsCurveType::BSplineCurve) | Some(GeomAbsCurveType::BezierCurve)
                    )
                }
                _ => false,
            };
            tocheck = !b;
        }
        let mut chk = true;
        if tocheck {
            chk = self.base.check_shape(s);
        }

        if !chk {
            return;
        }
        self.base.process_add_start_element(brep, s);
    }

    /// WireEdgeSet.cxx L117-120 — TopOpeBRepBuild_ShapeSet::AddElement
    /// (the explicit base-class call; the base AddElement body is
    /// ShapeSet.cxx L94-106).
    fn add_element(&mut self, brep: &BRep, s: &Shape) {
        let chk = self.base.check_shape(s);
        if !chk {
            return;
        }
        self.base.process_add_element(brep, s);
    }

    /// ShapeSet.cxx L155-158.
    fn start_elements(&self) -> &Vec<Shape> {
        self.base.start_elements()
    }

    /// ShapeSet.cxx L162-165.
    fn init_shapes(&mut self) {
        self.base.init_shapes();
    }

    /// ShapeSet.cxx L169-173.
    fn more_shapes(&self) -> bool {
        self.base.more_shapes()
    }

    /// ShapeSet.cxx L177-180.
    fn next_shape(&mut self) {
        self.base.next_shape();
    }

    /// ShapeSet.cxx L184-188.
    fn shape(&self) -> Shape {
        self.base.shape()
    }

    /// ShapeSet.cxx L192-195.
    fn init_start_elements(&mut self) {
        self.base.init_start_elements();
    }

    /// ShapeSet.cxx L199-203.
    fn more_start_elements(&self) -> bool {
        self.base.more_start_elements()
    }

    /// ShapeSet.cxx L207-210.
    fn next_start_element(&mut self) {
        self.base.next_start_element();
    }

    /// ShapeSet.cxx L214-218.
    fn start_element(&self) -> Shape {
        self.base.start_element()
    }

    /// WireEdgeSet.cxx L131-137.
    fn init_neighbours(&mut self, brep: &BRep, e: &Shape) {
        self.base.base_init_neighbours(brep, e);
        self.find_neighbours(brep);
    }

    /// ShapeSet.cxx L231-235.
    fn more_neighbours(&self) -> bool {
        self.base.more_neighbours()
    }

    /// ShapeSet.cxx L239-252 (the base NextNeighbour — not overridden by
    /// WireEdgeSet).
    fn next_neighbour(&mut self, brep: &BRep) {
        self.base.my_incident_shapes_iter += 1;
        let noisimore = !self.more_neighbours();
        if noisimore {
            let ssemore =
                self.base.my_sub_shape_explorer_pos < self.base.my_sub_shape_explorer.len();
            if ssemore {
                self.base.my_sub_shape_explorer_pos += 1;
                self.find_neighbours(brep);
            }
        }
    }

    /// ShapeSet.cxx L256-260.
    fn neighbour(&self) -> Shape {
        self.base.neighbour()
    }

    /// WireEdgeSet.cxx L141-166.
    fn find_neighbours(&mut self, brep: &BRep) {
        while self.base.my_sub_shape_explorer_pos < self.base.my_sub_shape_explorer.len() {
            // l = list of edges neighbour of edge myCurrentShape through
            // the vertex mySubShapeExplorer.Current().
            let v = self.base.my_sub_shape_explorer[self.base.my_sub_shape_explorer_pos].clone();
            let cur = self.base.my_current_shape.clone();
            let l = self.make_neighbours_list(brep, &cur, &v);

            self.base.my_incident_shapes = l;
            self.base.my_incident_shapes_iter = 0;
            if self.more_neighbours() {
                break;
            } else {
                self.base.my_sub_shape_explorer_pos += 1;
            }
        }
    }

    /// WireEdgeSet.cxx L172-288.
    fn make_neighbours_list(&mut self, brep: &BRep, earg: &Shape, varg: &Shape) -> Vec<Shape> {
        let e = earg.clone();
        let v = varg.clone();
        let l: Vec<Shape> = self
            .base
            .my_sub_shape_map
            .get(&shape_key(&v))
            .cloned()
            .unwrap_or_default();

        let nclosing = self.nb_closing_shapes(brep, &l);

        if nclosing != 0 {
            // build myCurrentShapeNeighbours = edge list made of connected
            // shapes to Earg through Varg
            self.base.my_current_shape_neighbours.clear();
            for curn in &l {
                // current neighbour
                let k = self.vertex_connects_edges_closing(brep, &v, &e, curn);
                if k {
                    self.base.my_current_shape_neighbours.push(curn.clone());
                }
            }

            let newn = self.nb_closing_shapes(brep, &self.base.my_current_shape_neighbours);

            if newn >= 2 {
                let f = self.my_face.clone();

                // plusieurs aretes de couture connexes a E par V et telles
                // que : orientation de V dans E # orientation de V dans ces
                // aretes.  on ne garde, parmi les aretes de couture
                // connexes, que l'arete A qui verifie tg(E) ^ tg(A) > 0
                let mut d1e = glam::DVec2::ZERO;
                let mut pe = glam::DVec2::ZERO;
                let par_e = brep_tool_parameter(brep, &v, &e);
                // trim3d = true; PCE = FC2D_CurveOnSurface(E, F, fiE, laE,
                // tolpc, trim3d).
                let pce = fc2d_curve_on_surface(brep, &e, &f, true);

                match pce {
                    Some((pce_c2d, _fie, _lae, _tolpc)) => {
                        // PCE->D1(parE, pE, d1E) — the 2D point + derivative
                        // at parE.
                        pe = pce_c2d.point_at(par_e);
                        d1e = pce_c2d.tangent_at(par_e);
                    }
                    None => {
                        self.local_d1(brep, &f, &e, &v, &mut pe, &mut d1e);
                    }
                }

                let e_ori = e.orientation;
                if e_ori == Orientation::Reversed {
                    d1e = -d1e;
                }

                let mut lclo: Vec<Shape> = std::mem::take(&mut self.base.my_current_shape_neighbours);
                let mut idx = 0usize;
                while idx < lclo.len() {
                    if !self.is_closed(brep, &lclo[idx]) {
                        idx += 1;
                        continue;
                    }

                    let ee = lclo[idx].clone();
                    let mut d1ee = glam::DVec2::ZERO;
                    let mut pee = glam::DVec2::ZERO;
                    let par_ee = brep_tool_parameter(brep, &v, &ee);
                    let pcee = fc2d_curve_on_surface(brep, &ee, &f, true);

                    match pcee {
                        Some((pcee_c2d, _fiee, _laee, _tolpc1)) => {
                            pee = pcee_c2d.point_at(par_ee);
                            d1ee = pcee_c2d.tangent_at(par_ee);
                        }
                        None => {
                            self.local_d1(brep, &f, &ee, &v, &mut pee, &mut d1ee);
                        }
                    }

                    let ee_ori = ee.orientation;
                    if ee_ori == Orientation::Reversed {
                        d1ee = -d1ee;
                    }

                    let cross = d1e.x * d1ee.y - d1e.y * d1ee.x;
                    let mut ove = Orientation::Forward;
                    let mut ovee = Orientation::Forward;
                    self.vertex_connects_edges(brep, &v, &e, &ee, &mut ove, &mut ovee);

                    let t2 = (cross > 0.0 && ove == Orientation::Reversed)
                        || (cross < 0.0 && ove == Orientation::Forward);

                    if t2 {
                        //-- t1 : c'est la bonne IsClosed, on ne garde
                        // qu'elle parmi les IsClosed
                        idx += 1;
                    } else {
                        // on vire l'arete IsClosed
                        lclo.remove(idx);
                    }
                }
                self.base.my_current_shape_neighbours = lclo;
            }
            return self.base.my_current_shape_neighbours.clone();
        } else {
            return l;
        }
    }

    /// ShapeSet.cxx L312-334.
    fn max_number_sub_shape(&mut self, brep: &BRep, shape: &Shape) -> i32 {
        self.base.max_number_sub_shape(brep, shape)
    }
}

/// The BRepAdaptor_Curve GetType result the AddStartElement check needs
/// (the OCCT GeomAbs_CurveType values consumed).  The canonical nine-variant
/// enum lives in rcad_kernel::math (GeomAbs_CurveType.hxx L23-31); the
/// former private three-variant copy (BSplineCurve/BezierCurve/Other) was
/// deleted (Rule 4).  WireEdgeSet.cxx L88-113 tests only BSpline/Bezier, so
/// every other canonical type takes the branch the local Other took.
use rcad_kernel::math::GeomAbsCurveType;
