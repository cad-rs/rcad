// OCCT LocOpe_SplitDrafts.cxx — the file-local statics (split out of
// loc_ope_split_drafts.rs for the single-file <2000-line rule):
//   NewPlane  cxx L1461-1538
//   MakeFace  cxx L1542-1712
//   Contains  cxx L1716-1727
//   NewEdge   cxx L1731-1976
// plus the API-surface carriers for the deferred dependencies
// (LocOpeSplitShape, GeomFillPipe, GeomIntIntSS, BRepToolsSubstitution —
// architecture differences #2/#3/#4/#5 of the module header).
//
// 1:1 translation; the architecture differences are numbered in the
// loc_ope_split_drafts.rs header.

use crate::brep_algo::tool::brep_tool_tolerance;
use crate::feat::loc_ope_split_drafts::{
    brep_tool_curve_on_surface, brep_tool_pnt, brep_tool_range, brep_tool_surface,
    builder_add_face_wire, geom_rectangular_trimmed_basis_surface,
    gp_pln_axis, gp_pln_direct, gp_pln_rotated, shape_key, top_abs_reverse, with_orientation,
};
use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Line3, Plane, Surface3, SurfaceEval,
};
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::precision::{ANGULAR, CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use rcad_kernel::topo::topods::{BRep, BRepBuilder};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::Orientation;
use crate::fillet::chfi3d_builder_0::topexp_vertices;
use std::collections::{HashMap, HashSet};

/// OCCT static NewPlane (cxx L1461-1538) — the new inclined plane.
#[allow(clippy::too_many_arguments)]
pub(crate) fn new_plane(
    the_f: &Shape,
    the_extr: DVec3,
    the_neutr: &Plane,
    the_ang: f64,
    the_newpl: &mut Plane,
    the_normal_f: &mut Ax1,
    the_modify: bool,
) -> bool {
    // Determination du nouveau plan incline
    // OCCT L1471-1480.
    let s0 = brep_tool_surface(the_f).expect("BRep_Tool::Surface");
    let s = geom_rectangular_trimmed_basis_surface(&s0);
    let plorig = match s {
        Surface3::Plane(p) => p,
        _ => return false, // occ::down_cast<Geom_Plane>(S).IsNull()
    };

    // OCCT L1483-1493.
    if !the_modify {
        *the_newpl = plorig;
        *the_normal_f = gp_pln_axis(the_newpl);
        if (gp_pln_direct(the_newpl) && the_f.orientation == Orientation::Reversed)
            || (!gp_pln_direct(the_newpl) && the_f.orientation == Orientation::Forward)
        {
            // OCCT L1490: NormalF.Reverse().
            the_normal_f.direction = -the_normal_f.direction;
        }
        return true;
    }

    // OCCT L1495-1496: gp_Ax1 Axe; double Theta; (Theta is assigned in the
    // intersection branch below).
    let mut axe = Ax1::new(DVec3::ZERO, DVec3::Z);

    // OCCT L1498-1535.
    let i2pl = intersect_plane_plane(the_neutr, &plorig);
    if let Some(lin_inters) = i2pl {
        // OCCT L1503-1505.
        let nx = lin_inters.direction;
        *the_normal_f = gp_pln_axis(&plorig);
        let ny = the_normal_f.direction.cross(nx);
        let a = the_extr.dot(nx);
        if a.abs() <= 1.0 - ANGULAR {
            // OCCT L1509-1512.
            let mut b = the_extr.dot(ny);
            let mut c = the_extr.dot(the_normal_f.direction);
            let direct = gp_pln_direct(&plorig);
            let oris = the_f.orientation;
            if (direct && oris == Orientation::Reversed)
                || (!direct && oris == Orientation::Forward)
            {
                b = -b;
                c = -c;
                // OCCT L1517: NormalF.Reverse().
                the_normal_f.direction = -the_normal_f.direction;
            }
            // OCCT L1519-1533.
            let denom = (1.0 - a * a).sqrt();
            let sina = the_ang.sin();
            if denom > sina.abs() {
                let phi = (b / denom).atan2(c / denom);
                let theta0 = (sina / denom).acos();
                let mut theta = theta0 - phi;
                if theta.cos() < 0.0 {
                    theta = -theta0 - phi;
                }
                axe = Ax1::new(lin_inters.origin, lin_inters.direction);
                *the_newpl = gp_pln_rotated(&plorig, &axe, theta);
                return true;
            }
        }
    }
    println!("fin newplane return standard_false");
    false
}

/// OCCT IntAna_QuadQuadGeo(Plorig, Neutr, Precision::Angular(),
/// Precision::Confusion()) with TypeInter() == IntAna_Line and Line(1) —
/// the plane-plane case (architecture difference #9).
fn intersect_plane_plane(the_p1: &Plane, the_p2: &Plane) -> Option<Line3> {
    match rcad_kernel::base::int_ana::intersect_plane_plane_intana(the_p1, the_p2) {
        rcad_kernel::base::int_ana::PlnPlnResult::Line(l) => Some(l),
        _ => None,
    }
}

/// OCCT static MakeFace (cxx L1542-1712) — builds the face from the edge
/// list (completes the missing p-curves, assembles the wires and fixes the
/// orientation through the surface properties).
pub(crate) fn make_face(
    the_f: &mut Shape,
    the_ledg: &mut Vec<Shape>,
    pool: &mut BRep,
    b: &mut BRepBuilder,
) {
    // ledg est une liste d'edge

    // Verification de l`existence des p-curves. Celles qui manquent
    // correspondent necessairement a des isos (et meme des iso u).

    // OCCT L1552-1561.
    let mut f = 0.0f64;
    let mut l = 0.0f64;
    let mut itl = 0usize;
    while itl < the_ledg.len() {
        let edg = the_ledg[itl].clone(); // TopoDS::Edge(itl.Value())
        let c2d = brep_tool_curve_on_surface(&edg, the_f);
        // OCCT L1559-1611.
        if let Some((_c2d, cf, cl)) = c2d {
            let _ = (cf, cl);
        } else {
            // OCCT L1561: BRep_Tool::Range(edg, f, l).
            let (rf, rl) = brep_tool_range(&edg);
            f = rf;
            l = rl;
            // OCCT L1562-1563.
            let (v1, v2) = topexp_vertices(&edg);
            // OCCT L1564-1606.
            let mut itl2 = 0usize;
            while itl2 < the_ledg.len() {
                let edg2 = the_ledg[itl2].clone(); // TopoDS::Edge(itl2.Value())
                if edg2.is_same(&edg) {
                    itl2 += 1;
                    continue;
                }
                let (vp1, vp2) = topexp_vertices(&edg2);
                if vp1.is_same(&v1) || vp2.is_same(&v1) || vp1.is_same(&v2) || vp2.is_same(&v2) {
                    // OCCT L1576-1577.
                    let c22d = brep_tool_curve_on_surface(&edg2, the_f);
                    if let Some((c22d, f2, l2)) = c22d {
                        // OCCT L1580-1600.
                        let mut pt2d = DVec2::ZERO;
                        if vp1.is_same(&v1) {
                            pt2d = c22d.point_at(f2);
                            pt2d.y = pt2d.y - f;
                        } else if vp2.is_same(&v1) {
                            pt2d = c22d.point_at(l2);
                            pt2d.y = pt2d.y - f;
                        } else if vp1.is_same(&v2) {
                            pt2d = c22d.point_at(f2);
                            pt2d.y = pt2d.y - l;
                        } else if vp2.is_same(&v2) {
                            pt2d = c22d.point_at(l2);
                            pt2d.y = pt2d.y - l;
                        }
                        // OCCT L1601-1603.
                        let c2d_new = Curve2d::Line(Line2d::new(pt2d, DVec2::new(0.0, 1.0)));
                        b.update_edge_pcurve(
                            pool,
                            edg.clone(),
                            c2d_new,
                            the_f.clone(),
                            brep_tool_tolerance(&edg),
                        );
                        break;
                    }
                }
                itl2 += 1;
            }
            // OCCT L1607-1610.
            let still_null = brep_tool_curve_on_surface(&edg, the_f).is_none();
            if still_null {
                println!("Ca merde violemment");
            }
        }
        itl += 1;
    }

    // OCCT L1614-1688.
    let mut lwires: Vec<Shape> = Vec::new();
    let mut alldone = the_ledg.is_empty();
    while !alldone {
        // OCCT L1618-1619.
        let wnew = b.make_wire(pool);
        // OCCT L1620-1631.
        let edg = the_ledg[0].clone(); // TopoDS::Edge(ledg.First())
        let (mut vfirst, mut vlast) = if edg.orientation == Orientation::Forward {
            topexp_vertices(&edg)
        } else {
            let (a, bl) = topexp_vertices(&edg);
            (bl, a)
        };
        b.add_to_wire(pool, wnew.clone(), edg);
        the_ledg.remove(0);
        // on suppose VFirst et VLast non nuls
        // OCCT L1635.
        let mut wdone = the_ledg.is_empty() || vfirst.is_same(&vlast);
        while !wdone {
            // OCCT L1638-1640.
            let mut vf = Shape::null();
            let mut vl = Shape::null();
            let mut oredg = Orientation::Forward;
            let mut hit: Option<usize> = None;

            // OCCT L1642-1672.
            for (i2, s2) in the_ledg.iter().enumerate() {
                let edg2 = s2; // TopoDS::Edge(itl.Value())
                let fwd = with_orientation(edg2, Orientation::Forward);
                let (a, bpt) = topexp_vertices(&fwd);
                vf = a;
                vl = bpt;
                if vf.is_same(&vlast) {
                    vlast = vl.clone();
                    oredg = Orientation::Forward;
                    hit = Some(i2);
                    break;
                } else if vl.is_same(&vfirst) {
                    vfirst = vf.clone();
                    oredg = Orientation::Forward;
                    hit = Some(i2);
                    break;
                } else if vf.is_same(&vfirst) {
                    vfirst = vl.clone();
                    oredg = Orientation::Reversed;
                    hit = Some(i2);
                    break;
                } else if vl.is_same(&vlast) {
                    vlast = vf.clone();
                    oredg = Orientation::Reversed;
                    hit = Some(i2);
                    break;
                }
            }

            // OCCT L1673-1684.
            if hit.is_none() {
                wdone = true;
            } else {
                let i2 = hit.expect("itl");
                let a_local_shape = with_orientation(&the_ledg[i2], oredg);
                b.add_to_wire(pool, wnew.clone(), a_local_shape);
                the_ledg.remove(i2);
                wdone = the_ledg.is_empty() || vfirst.is_same(&vlast);
            }
        }
        // OCCT L1686-1687.
        lwires.push(wnew);
        alldone = the_ledg.is_empty();
    }

    // OCCT L1690-1711.
    the_f.orientation = Orientation::Forward;
    for i in 0..lwires.len() {
        // OCCT L1693-1696.
        let mut new_face = pool.empty_copied(the_f);
        builder_add_face_wire(pool, &mut new_face, lwires[i].clone());
        // OCCT L1697-1699.
        let gp = crate::feat::loc_ope_split_drafts::brep_gprop_surface_properties_mass(pool, &new_face);
        if gp < 0.0 {
            // OCCT L1701: itl.ChangeValue().Reverse().
            lwires[i] = with_orientation(&lwires[i], top_abs_reverse(lwires[i].orientation));
        }
    }
    // OCCT L1704-1711.
    if lwires.len() == 1 {
        builder_add_face_wire(pool, the_f, lwires[0].clone());
    } else {
        println!("Not yet implemented : nbwire >= 2");
    }
}

/// OCCT static Contains (cxx L1716-1727).
pub(crate) fn contains(the_ll: &[Shape], the_s: &Shape) -> bool {
    for itl in the_ll {
        if itl.is_same(the_s) {
            return true;
        }
    }
    false
}

/// OCCT static NewEdge (cxx L1731-1976) — the edge of the draft surface
/// through the given vertices; a null Shape carries the OCCT null
/// TopoDS_Edge.
#[allow(clippy::too_many_arguments)]
pub(crate) fn new_edge(
    the_edg: &Shape,
    the_f: &Shape,
    new_s: &Surface3,
    the_v1: &Shape,
    the_v2: &Shape,
    pool: &mut BRep,
    b: &mut BRepBuilder,
) -> Shape {
    // OCCT L1737-1743.
    let s1 = brep_tool_surface(the_f).expect("BRep_Tool::Surface");
    let mut app_s1 = false;
    if !matches!(s1, Surface3::Plane(_)) {
        app_s1 = true;
    }

    // OCCT L1745-1749.
    let mut i2s = GeomIntIntSS::new();
    i2s.perform(&s1, new_s, CONFUSION, true, app_s1, false);
    if !i2s.is_done() || i2s.nb_lines() == 0 {
        return Shape::null();
    }

    // OCCT L1751-1756.
    let mut prmf = 0.0f64;
    let mut prml = 0.0f64;

    // OCCT L1758-1760.
    let pvf = brep_tool_pnt(the_v1);
    let pvl = brep_tool_pnt(the_v2);
    let nb_lines = i2s.nb_lines();
    let mut i = 1usize;
    while i <= nb_lines {
        let crv = i2s.line(i).expect("GeomInt_IntSS::Line");
        // OCCT L1762-1763: TheCurve.Load(i2s.Line(i)); Extrema_ExtPC
        // myExtPC(pvf, TheCurve) — the two-arg ctor over the full domain,
        // the default theTolF is 1.0e-10.
        let dom = crv.default_domain();
        let a_adaptor = GeomCurveAdaptor::new(crv.clone());
        let a_tool = CurveToolHandle::for_curve3(&crv, &a_adaptor, &a_adaptor);
        let mut my_ext_pc = ExtremaExtPC::new_point_curve(pvf, &a_tool, 1.0e-10);

        if my_ext_pc.is_done() {
            // OCCT L1767-1769.
            let mut thepmin = dom[0]; // TheCurve.FirstParameter()
            // OCCT L1768: myExtPC.TrimmedSquareDistances(Dist2Min, Dist2,
            // p1b, p2b).
            let (d2first, d2last, _p1b, _p2b) = my_ext_pc.trimmed_square_distances();
            let mut dist2_min = d2first;
            let mut dist2 = d2last;
            // OCCT L1770-1774.
            if dist2 < dist2_min && !crv.is_periodic() {
                dist2_min = dist2;
                thepmin = dom[1]; // TheCurve.LastParameter()
            }
            // OCCT L1775-1783.
            for k in 1..=my_ext_pc.nb_ext() {
                dist2 = my_ext_pc.square_distance(k);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    thepmin = my_ext_pc.point(k).param;
                }
            }

            // OCCT L1785-1814.
            if dist2_min <= SQUARE_CONFUSION {
                prmf = thepmin;
                // OCCT L1788: myExtPC.Perform(pvl).
                my_ext_pc.perform(pvl);
                if my_ext_pc.is_done() {
                    // OCCT L1791-1796.
                    let mut thepmin = dom[1];
                    // OCCT L1789: myExtPC.TrimmedSquareDistances(Dist2,
                    // Dist2Min, p1b, p2b) — the first output lands in the
                    // Dist2 variable (the FIRST endpoint square distance),
                    // the second in Dist2Min (the LAST endpoint).
                    let (d2_out, d2min_out, _p1b2, _p2b2) =
                        my_ext_pc.trimmed_square_distances();
                    let mut dist2_min2 = d2min_out;
                    let mut dist22 = d2_out;
                    // OCCT: `if (Dist2 < Dist2Min && !TheCurve.IsClosed())`
                    // compares the FIRST endpoint square distance against the
                    // LAST one.
                    if dist22 < dist2_min2 && !crv.is_closed() {
                        dist2_min2 = dist22;
                        thepmin = dom[0];
                    }
                    // OCCT L1798-1806.
                    for k in 1..=my_ext_pc.nb_ext() {
                        dist22 = my_ext_pc.square_distance(k);
                        if dist22 < dist2_min2 {
                            dist2_min2 = dist22;
                            thepmin = my_ext_pc.point(k).param;
                        }
                    }

                    // OCCT L1808-1812.
                    if dist2_min2 <= SQUARE_CONFUSION {
                        prml = thepmin;
                        break;
                    }
                }
            }
        }
        i += 1;
    }

    // OCCT L1818-1974.
    if i <= nb_lines {
        // OCCT L1820-1828.
        let mut rev = false;
        let mut vf = the_v1.clone();
        let mut vl = the_v2.clone();
        let cimg = i2s.line(i).expect("GeomInt_IntSS::Line");
        let mut cimg2d: Option<Curve2d> = None;
        if app_s1 {
            cimg2d = i2s.line_on_s1(i);
        }

        // OCCT L1830-1916.
        if cimg.is_periodic() {
            let dom = cimg.default_domain();
            // OCCT: period = Cimg->Period() — the intrinsic period equals
            // LastParameter - FirstParameter for the bounded periodic
            // curves (Geom_Curve::Period()).
            let period = dom[1] - dom[0];
            let imf = dom[0];
            let iml = dom[1];

            // OCCT L1837-1841.
            let (f, l) = brep_tool_range(the_edg);
            let delt = l - f;
            let delt1 = (prml - prmf).abs();
            let delt2 = (period - delt1).abs();

            // OCCT L1843-1908.
            if delt1 == 0.0 || delt2 == 0.0 {
                prmf = imf;
                prml = iml;
            } else if (delt1 - delt).abs() > (delt2 - delt).abs() {
                // le bon ecart est delt2...
                if prml > prmf {
                    if prml < iml {
                        prmf += period;
                    } else {
                        prml -= period;
                    }
                } else if prmf < iml {
                    prml += period;
                } else {
                    prmf -= period;
                }
            } else if (delt1 - delt).abs() < (delt2 - delt).abs() {
                if prmf >= iml && prml >= iml {
                    prmf -= period;
                    prml -= period;
                } else if prmf <= imf && prml <= imf {
                    prmf += period;
                    prml += period;
                }
            } else {
                // egalite; on priveligie l'ordre f,l
                if prmf > prml {
                    prmf -= period;
                }
                if prmf >= iml && prml >= iml {
                    prmf -= period;
                    prml -= period;
                } else if prmf <= imf && prml <= imf {
                    prmf += period;
                    prml += period;
                }
            }
        }

        // OCCT L1918-1946.
        let s1_eval = &s1;
        if s1_eval.is_u_periodic() {
            // OCCT: speriod = S1->UPeriod() — the U period equals the
            // u-domain span for the rcad periodic surfaces (GAP with the
            // GeomInt_IntSS translation for the non-standard domains).
            let sdom = s1_eval.default_domain();
            let speriod = sdom[1] - sdom[0];
            // OCCT L1923-1926.
            let c2 = cimg2d.clone().expect("Cimg2d");
            let pf = c2.point_at(prmf);
            let pl = c2.point_at(prml);

            // OCCT L1928-1941.
            let uf = pf.x;
            let ul = pl.x;
            let mut ptra = 0.0f64;

            let mut ustart = uf.min(ul);
            while ustart < -PCONFUSION {
                ustart += speriod;
                ptra += speriod;
            }
            while ustart > speriod - PCONFUSION {
                ustart -= speriod;
                ptra -= speriod;
            }
            // OCCT L1942-1945: Cimg2d->Translate(gp_Vec2d(ptra, 0.)) — the
            // rcad Curve2d carries a translation vehicle for the line
            // variant only (GAP with the GeomInt_IntSS translation).
            if ptra != 0.0 {
                if let Some(Curve2d::Line(ref mut l2)) = cimg2d {
                    l2.origin += DVec2::new(ptra, 0.0);
                }
            }
        }
        // OCCT L1947-1957.
        if prmf < prml {
            vf.orientation = Orientation::Forward;
            vl.orientation = Orientation::Reversed;
        } else {
            vf.orientation = Orientation::Reversed;
            vl.orientation = Orientation::Forward;
            rev = true;
        }

        // OCCT L1959-1968: MakeEdge(Cimg) + Add(Vf, Vl) + UpdateVertex
        // (prmf, prml) + UpdateEdge(Cimg2d).
        let new_edg = b.add_edge(pool, Some(cimg.clone()), vf.clone(), vl.clone(), [prmf, prml]);
        b.set_vertex_param(pool, new_edg.clone(), vf.clone(), prmf);
        b.set_vertex_param(pool, new_edg.clone(), vl.clone(), prml);
        if let Some(c2) = &cimg2d {
            b.update_edge_pcurve(pool, new_edg.clone(), c2.clone(), the_f.clone(), CONFUSION);
        }

        // OCCT L1970-1973.
        if rev {
            let mut out = new_edg;
            out.orientation = Orientation::Reversed;
            return out;
        }
        return new_edg;
    }
    // OCCT L1976: return NewEdg; (null when the line loop ran out).
    Shape::null()
}


// ---------------------------------------------------------------------------
// API-surface carriers for the deferred dependencies (architecture
// differences #2/#3/#4/#5).
// ---------------------------------------------------------------------------

/// OCCT LocOpe_SplitShape (LocOpe_SplitShape.hxx L38-91) — deferred body;
/// the member set and the surface consumed by LocOpe_SplitDrafts::Perform
/// are carried below (architecture difference #2).
pub(crate) struct LocOpeSplitShape {
    my_done: bool,                           // OCCT: myDone
    my_shape: Shape,                         // OCCT: myShape
    my_map: HashMap<(u64, u32), Vec<Shape>>, // OCCT: myMap
    my_dble: HashSet<(u64, u32)>,            // OCCT: myDblE
    my_left: Vec<Shape>,                     // OCCT: myLeft
}

impl LocOpeSplitShape {
    /// OCCT LocOpe_SplitShape::LocOpe_SplitShape(S) (hxx L47) — deferred
    /// body; the shape is carried.
    pub(crate) fn new_with_shape(the_s: &Shape) -> Self {
        LocOpeSplitShape {
            my_done: false,
            my_shape: the_s.clone(),
            my_map: HashMap::new(),
            my_dble: HashSet::new(),
            my_left: Vec::new(),
        }
    }

    /// OCCT LocOpe_SplitShape::Add(V, P, E) (hxx L56) — deferred body; the
    /// bookkeeping lands with the LocOpe_SplitShape translation.
    pub(crate) fn add(&mut self, _the_v: &Shape, _the_p: f64, _the_e: &Shape) {
        // deferred: the map stays unbound until the body lands.
    }

    /// OCCT LocOpe_SplitShape::DescendantShapes(S) (hxx L68) — the myMap
    /// lookup; an absent key carries the OCCT absent DataMap entry.
    pub(crate) fn descendant_shapes(&mut self, the_s: &Shape) -> Option<&Vec<Shape>> {
        self.my_map.get(&shape_key(the_s))
    }
}

/// OCCT GeomFill_Pipe (GeomFill_Pipe.hxx) — the real TKGeomAlgo translation
/// (geomalgo/geomfill/pipe.rs, sweep-closure batch 2).  The deferred-body
/// stub is gone; the LocOpe_SplitDrafts call sites carry the real class.
pub(crate) use crate::geomalgo::geomfill::pipe::Pipe as GeomFillPipe;

/// OCCT GeomInt_IntSS (GeomInt_IntSS.hxx L30-120) — deferred body; the
/// member surface consumed by LocOpe_SplitDrafts (myDone, myLines,
/// myLinesOnS1, myLinesOnS2 and the AppS1/AppS2 Perform overloads) is
/// carried below (architecture difference #4).
pub(crate) struct GeomIntIntSS {
    my_done: bool,              // OCCT: myDone
    my_s1: Option<Surface3>,    // OCCT: myS1 (carried argument)
    my_s2: Option<Surface3>,    // OCCT: myS2 (carried argument)
    my_lines: Vec<Curve3>,      // OCCT: myLines
    my_lines_on_s1: Vec<Curve2d>, // OCCT: myLinesOnS1
    my_lines_on_s2: Vec<Curve2d>, // OCCT: myLinesOnS2
}

impl GeomIntIntSS {
    /// OCCT GeomInt_IntSS::GeomInt_IntSS() — empty.
    pub(crate) fn new() -> Self {
        GeomIntIntSS {
            my_done: false,
            my_s1: None,
            my_s2: None,
            my_lines: Vec::new(),
            my_lines_on_s1: Vec::new(),
            my_lines_on_s2: Vec::new(),
        }
    }

    /// OCCT GeomInt_IntSS::GeomInt_IntSS(S1, S2, Tol) — the constructor
    /// performs (GeomInt_IntSS.cxx L59-63); deferred body (architecture
    /// difference #4).
    pub(crate) fn new_with_surfaces(the_s1: &Surface3, the_s2: &Surface3, _the_tol: f64) -> Self {
        let mut i2s = GeomIntIntSS::new();
        i2s.my_s1 = Some(the_s1.clone());
        i2s.my_s2 = Some(the_s2.clone());
        // deferred: IsDone stays false until the GeomInt_IntSS body lands.
        i2s
    }

    /// OCCT GeomInt_IntSS::Perform(S1, S2, Tol, Deflection, AppS1, AppS2)
    /// (GeomInt_IntSS.cxx) — deferred body (architecture difference #4).
    pub(crate) fn perform(
        &mut self,
        the_s1: &Surface3,
        the_s2: &Surface3,
        _the_tol: f64,
        _the_deflection: bool,
        _the_app_s1: bool,
        _the_app_s2: bool,
    ) {
        self.my_s1 = Some(the_s1.clone());
        self.my_s2 = Some(the_s2.clone());
        // deferred: IsDone stays false until the GeomInt_IntSS body lands.
    }

    /// OCCT GeomInt_IntSS::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT GeomInt_IntSS::NbLines().
    pub(crate) fn nb_lines(&self) -> usize {
        self.my_lines.len()
    }

    /// OCCT GeomInt_IntSS::Line(I) — None carries the OCCT out-of-range
    /// raise.
    pub(crate) fn line(&self, the_i: usize) -> Option<Curve3> {
        self.my_lines.get(the_i - 1).cloned()
    }

    /// OCCT GeomInt_IntSS::LineOnS1(I).
    pub(crate) fn line_on_s1(&self, the_i: usize) -> Option<Curve2d> {
        self.my_lines_on_s1.get(the_i - 1).cloned()
    }

    /// OCCT GeomInt_IntSS::LineOnS2(I).
    pub(crate) fn line_on_s2(&self, the_i: usize) -> Option<Curve2d> {
        self.my_lines_on_s2.get(the_i - 1).cloned()
    }
}

/// OCCT BRepTools_Substitution (BRepTools_Substitution.cxx) — the
/// Substitute/Copy/IsCopied map operations are the real OCCT ones
/// (BRepTools_Substitution.cxx L30-58); Build(S) is deferred (architecture
/// difference #5).
pub(crate) struct BRepToolsSubstitution {
    my_map: HashMap<(u64, u32), Vec<Shape>>, // OCCT: myMap
}

impl BRepToolsSubstitution {
    /// OCCT BRepTools_Substitution default state — an empty copy map.
    pub(crate) fn new() -> Self {
        BRepToolsSubstitution {
            my_map: HashMap::new(),
        }
    }

    /// OCCT BRepTools_Substitution::Substitute(S, L) (cxx L30-40).
    pub(crate) fn substitute(&mut self, the_s: &Shape, the_l: Vec<Shape>) {
        self.my_map.insert(shape_key(the_s), the_l);
    }

    /// OCCT BRepTools_Substitution::Build(S) (cxx L45-110) — deferred body:
    /// the ancestor rebuild walk needs the TopoDS rebuild vehicle.
    pub(crate) fn build(&mut self, _the_s: &Shape) {
        // deferred: the substitution map is carried; the rebuild is a no-op
        // until the body lands.
    }

    /// OCCT BRepTools_Substitution::IsCopied(S) (cxx).
    pub(crate) fn is_copied(&self, the_s: &Shape) -> bool {
        self.my_map.contains_key(&shape_key(the_s))
    }

    /// OCCT BRepTools_Substitution::Copy(S) (cxx) — None carries the OCCT
    /// absent-entry raise.
    pub(crate) fn copy(&self, the_s: &Shape) -> Option<&Vec<Shape>> {
        self.my_map.get(&shape_key(the_s))
    }
}
