// OCCT Contap_Contour.cxx — the file-static helpers.
//
//   Recadre                     (cxx L229-281)
//   LineConstructor             (cxx L283-511)
//   KeepInsidePoints            (cxx L515-554)
//   ComputeTangency             (cxx L556-891)
//   ComputeTransitionOnLine     (cxx L893-956)
//   ProcessSegments             (cxx L958-1118)
//   ComputeInternalPointsOnRstr (cxx L1120-1284)
//   ComputeInternalPoints       (cxx L1286-1539)
//   FindLine                    (cxx L1973-2011)
//   PutPointsOnLine             (cxx L2013-2091)
//   ComputeTransitionOngpLine   (cxx L2100-2126)
//   ComputeTransitionOngpCircle (cxx L2128-2154)
//
// `Surf` is the `occ::handle<Adaptor3d_Surface>` → `&dyn SurfaceAdapter`.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Circle3, Line3, Point3, Vec3};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{Orientation, State};

use crate::geomalgo::int_patch::imp_prm::function_set_root::FunctionSetRoot;
use crate::geomalgo::int_patch::imp_prm::path_point::{InteriorPoint, PathPoint};

use crate::geomalgo::int_patch::transitions::{make_transition, Transition, TypeTrans};
use crate::geomalgo::int_surf::{LineOn2S, PntOn2S};
use crate::geomalgo::top_trans::CurveTransition;
use crate::hlr::contap::contour::TOLE;
use crate::hlr::contap::domain::ContapDomain;
use crate::hlr::contap::h_cont_tool as hcont;
use crate::hlr::contap::h_curve2d_tool as hcurve2d;
use crate::hlr::contap::i_type::IType;
use crate::hlr::contap::line::Line;
use crate::hlr::contap::point::Point as ContapPoint;
use crate::hlr::contap::surf_function::SurfFunction;
use crate::hlr::contap::surf_props;
use crate::hlr::contap::surface_adaptor::SurfaceAdapter;
use crate::hlr::contap::the_path_point_of_the_search::ThePathPointOfTheSearch;
use crate::hlr::contap::the_search::TheSearch;
use crate::hlr::contap::the_search_inside::TheSearchInside;
use crate::hlr::contap::t_function::TFunction;

/// OCCT RealEpsilon().
const REAL_EPSILON: f64 = 2.220446049250313e-16;
/// OCCT gp::Resolution().
const GP_RESOLUTION: f64 = 1e-15;

// ---------------------------------------------------------------------------
// Recherche des portions utiles sur les lignes
// ---------------------------------------------------------------------------

/// OCCT Recadre (cxx L229-281).
pub(crate) fn recadre(my_hs1: &dyn SurfaceAdapter, u1: &mut f64, v1: &mut f64) {
    let typs1 = my_hs1.get_type();

    let (my_hs1_is_u_periodic, my_hs1_is_v_periodic);
    match typs1 {
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere => {
            my_hs1_is_u_periodic = true;
            my_hs1_is_v_periodic = false;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Torus => {
            my_hs1_is_u_periodic = true;
            my_hs1_is_v_periodic = true;
        }
        _ => {
            my_hs1_is_u_periodic = false;
            my_hs1_is_v_periodic = false;
        }
    }
    if my_hs1_is_u_periodic {
        let lmf = std::f64::consts::PI + std::f64::consts::PI; //-- myHS1->UPeriod();
        let f = my_hs1.first_u_parameter();
        let l = my_hs1.last_u_parameter();
        while *u1 < f {
            *u1 += lmf;
        }
        while *u1 > l {
            *u1 -= lmf;
        }
    }
    if my_hs1_is_v_periodic {
        let lmf = std::f64::consts::PI + std::f64::consts::PI; //-- myHS1->VPeriod();
        let f = my_hs1.first_v_parameter();
        let l = my_hs1.last_v_parameter();
        while *v1 < f {
            *v1 += lmf;
        }
        while *v1 > l {
            *v1 -= lmf;
        }
    }
}

/// OCCT LineConstructor (cxx L283-511) — decoupe la ligne en portions
/// entre 2 vertex.
pub(crate) fn line_constructor(
    slin: &mut Vec<Line>,
    domain: &mut dyn ContapDomain,
    l: &mut Line,
    surf: &dyn SurfaceAdapter,
) {
    let tol = rcad_kernel::precision::PCONFUSION;
    let typl = l.type_contour();
    if typl == IType::Walking {
        let nbvtx = l.nb_vertex();
        for i in 1..nbvtx {
            let firstp = l.vertex(i).parameter_on_line() as usize;
            let lastp = l.vertex(i + 1).parameter_on_line() as usize;
            if firstp != lastp {
                let pmid = (firstp + lastp) / 2; //-- entiers
                let pmid_v = l.point(pmid);
                let (mut u1, mut v1, mut u2, mut v2) = pmid_v.parameters();
                recadre(surf, &mut u2, &mut v2);
                let in2 = domain.classify(DVec2::new(u2, v2), tol, true);
                if in2 == State::Out {
                    // (empty branch in OCCT)
                } else {
                    let mut line_on_2s = LineOn2S::new();
                    line_on_2s.add(l.point(firstp));
                    for j in (firstp + 1)..=lastp {
                        let a_sq_dist =
                            l.point(j).value().distance_squared(l.point(j - 1).value());
                        if a_sq_dist > GP_RESOLUTION {
                            line_on_2s.add(l.point(j));
                        }
                    }
                    if line_on_2s.nb_points() < 2 {
                        continue;
                    }
                    let mut line = Line::new();
                    line.set_line_on_2s(&line_on_2s);
                    let mut pvtx = l.vertex(i).clone();
                    pvtx.set_parameter(1.0);
                    line.add(pvtx);

                    let mut pvtx = l.vertex(i + 1).clone();
                    pvtx.set_parameter(line_on_2s.nb_points() as f64);
                    line.add(pvtx);
                    line.set_transition_on_s(l.transition_on_s());
                    slin.push(line);
                }
            }
        }
    } else if typl == IType::Lin {
        let nbvtx = l.nb_vertex();
        for i in 1..nbvtx {
            let firstp = l.vertex(i).parameter_on_line();
            let lastp = l.vertex(i + 1).parameter_on_line();
            if firstp != lastp {
                let pmid = (firstp + lastp) * 0.5;
                let lin = l.line();
                let pmid_p = crate::geomalgo::int_patch::elclib::line_value(&lin, pmid);
                let mut u2 = 0.0;
                let mut v2 = 0.0;
                let cyl_cone =
                    elslib_parameters(surf, pmid_p, &mut u2, &mut v2);
                let _ = cyl_cone;
                let _ = (firstp, lastp);

                recadre(surf, &mut u2, &mut v2);
                let in2 = domain.classify(DVec2::new(u2, v2), tol, true);
                if in2 == State::Out {
                    // (empty branch in OCCT)
                } else {
                    let mut line = Line::new();
                    line.set_value_line(&lin);
                    let pvtx = l.vertex(i).clone();
                    line.add(pvtx);

                    let pvtx = l.vertex(i + 1).clone();
                    line.add(pvtx);
                    line.set_transition_on_s(l.transition_on_s());
                    slin.push(line);
                }
            }
        }
    } else if typl == IType::Circle {
        let nbvtx = l.nb_vertex();
        let mut novtx = true;
        if nbvtx > 0 {
            novtx = false;
        }
        let mut i = 1;
        while i < nbvtx || novtx {
            let (mut firstp, mut lastp) = (0.0, std::f64::consts::PI + std::f64::consts::PI);
            if !novtx {
                firstp = l.vertex(i).parameter_on_line();
                lastp = l.vertex(i + 1).parameter_on_line();
            }
            if (firstp - lastp).abs() > 0.000000001 {
                let pmid = (firstp + lastp) * 0.5;
                let cir = l.circle();
                let pmid_p = crate::geomalgo::int_patch::elclib::circle_value(&cir, pmid);
                let mut u2 = 0.0;
                let mut v2 = 0.0;
                elslib_parameters(surf, pmid_p, &mut u2, &mut v2);

                recadre(surf, &mut u2, &mut v2);
                let in2 = domain.classify(DVec2::new(u2, v2), tol, true);
                if in2 == State::Out {
                    // (empty branch in OCCT)
                } else {
                    let mut line = Line::new();
                    line.set_value_circle(&cir);
                    if !novtx {
                        let pvtx = l.vertex(i).clone();
                        line.add(pvtx);
                        let pvtx = l.vertex(i + 1).clone();
                        line.add(pvtx);
                    }
                    line.set_transition_on_s(l.transition_on_s());
                    slin.push(line);
                }
            }
            novtx = false;
            i += 1;
        }
        if nbvtx > 0 {
            let firstp = l.vertex(nbvtx).parameter_on_line();
            let lastp = l.vertex(1).parameter_on_line()
                + std::f64::consts::PI
                + std::f64::consts::PI;
            if (firstp - lastp).abs() > 0.0000000001 {
                let pmid = (firstp + lastp) * 0.5;
                let cir = l.circle();
                let pmid_p = crate::geomalgo::int_patch::elclib::circle_value(&cir, pmid);
                let mut u2 = 0.0;
                let mut v2 = 0.0;
                elslib_parameters(surf, pmid_p, &mut u2, &mut v2);

                recadre(surf, &mut u2, &mut v2);
                let in2 = domain.classify(DVec2::new(u2, v2), tol, true);
                if in2 == State::Out {
                    // (empty branch in OCCT)
                } else {
                    let mut line = Line::new();
                    line.set_value_circle(&cir);
                    let pvtx = l.vertex(nbvtx).clone();
                    line.add(pvtx);

                    let mut pvtx = l.vertex(1).clone();
                    pvtx.set_parameter(
                        pvtx.parameter_on_line() + std::f64::consts::PI + std::f64::consts::PI,
                    );
                    line.add(pvtx);
                    line.set_transition_on_s(l.transition_on_s());
                    slin.push(line);
                }
            }
        }
    } else {
        //-- std::cout<<" ni WLine ni Lin ni Circ "<<std::endl;
        slin.push(l.clone());
    }
}

/// The ElSLib::Parameters(Cylinder/Cone/Sphere, P, u, v) switch of
/// LineConstructor (cxx L359-370, L416-431, L465-480).
fn elslib_parameters(surf: &dyn SurfaceAdapter, p: Point3, u2: &mut f64, v2: &mut f64) -> bool {
    match surf.get_type() {
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder => {
            let cy = surf.cylinder();
            let (u, v) = rcad_kernel::math::el::elslib_cylinder_parameters(
                p,
                cy.origin,
                cy.ref_dir.normalize_or_zero(),
                cy.y_axis(),
                cy.axis,
                cy.radius,
            );
            *u2 = u;
            *v2 = v;
            true
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone => {
            let co = surf.cone();
            let (u, v) = rcad_kernel::math::el::elslib_cone_parameters(
                p,
                co.apex,
                co.ref_dir.normalize_or_zero(),
                co.axis.cross(co.ref_dir).normalize_or_zero(),
                co.axis,
                co.radius,
                co.half_angle_rad,
            );
            *u2 = u;
            *v2 = v;
            true
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere => {
            let sp = surf.sphere();
            let (u, v) = rcad_kernel::math::el::elslib_sphere_parameters(
                p,
                sp.center,
                sp.ref_dir.normalize(),
                sp.axis.cross(sp.ref_dir).normalize(),
                sp.axis,
            );
            *u2 = u;
            *v2 = v;
            true
        }
        _ => false,
    }
}

/// OCCT KeepInsidePoints (cxx L515-554).
pub(crate) fn keep_inside_points(
    solins: &TheSearchInside,
    solrst: &TheSearch,
    func: &SurfFunction,
    seqpins: &mut Vec<InteriorPoint>,
) {
    let nba = solrst.nb_segments();
    let surf = func.surface();

    let nbp = solins.nb_points();
    for indp in 1..=nbp {
        let mut tokeep = true;
        let pti = solins.value(indp);
        let (mut u, mut v) = (0.0, 0.0);
        pti.parameters(&mut u, &mut v);
        let toproj = DVec2::new(u, v);
        for inda in 1..=nba {
            let thearc = solrst.segment(inda).curve();
            if let Some((_paramproj, ptproj)) = hcont::project(thearc.as_ref(), toproj) {
                let pprojete = surf.value(ptproj.x, ptproj.y);
                if pti.value().distance(pprojete) <= CONFUSION {
                    tokeep = false;
                    break;
                }
            }
        }
        if tokeep {
            seqpins.push(pti.clone());
        }
    }
}

/// OCCT ComputeTangency (cxx L556-891).
pub(crate) fn compute_tangency(
    solrst: &TheSearch,
    domain: &mut dyn ContapDomain,
    func: &mut SurfFunction,
    seqpdep: &mut Vec<PathPoint>,
    destination: &mut [i32],
) {
    let nb_points = solrst.nb_points();
    let mut seqlength = 0usize;

    let surf = func.surface().clone();

    for i in 1..=nb_points {
        if destination[i] == 0 {
            let pstart = solrst.point(i);
            let thearc = pstart.arc();
            let theparam = pstart.parameter();
            let ptoproj = hcurve2d::value(thearc.as_ref(), theparam);
            //-- lbr le 15 mai 97
            //-- On elimine les points qui sont egalement present sur une
            //   restriction solution
            let mut sur_une_restriction_solution = false;
            let mut restriction = 1usize;
            while !sur_une_restriction_solution && restriction <= solrst.nb_segments() {
                let thearcsol = solrst.segment(restriction).curve();
                if let Some((_paramproj, pproj)) = hcont::project(thearcsol.as_ref(), ptoproj) {
                    // gp_Pnt pprojete = Value(Surf, pproj.X(), pproj.Y()); (IFV)
                    let pprojete = surf.value(pproj.x, pproj.y);
                    if pstart.value().distance(pprojete) <= CONFUSION {
                        sur_une_restriction_solution = true;
                    }
                }
                restriction += 1;
            }
            if !sur_une_restriction_solution {
                let arcorien = domain.orientation_arc(&thearc);
                let mut ispassing =
                    arcorien == Orientation::Internal || arcorien == Orientation::External;

                let (pt2d, tg2drst) = hcurve2d::d1(thearc.as_ref(), theparam);
                let mut x = [pt2d.x, pt2d.y];
                let mut ppoint = PathPoint::new();
                ppoint.set_value(pstart.value(), x[0], x[1]);

                let values_ok = func.values(&x).is_some();
                if values_ok && func.is_tangent() {
                    ppoint.set_tangency(true);
                    destination[i] = seqlength as i32 + 1;
                    if !pstart.is_new() {
                        let vtx = pstart.vertex().clone();
                        for k in (i + 1)..=nb_points {
                            if destination[k] == 0 {
                                let pstart2 = solrst.point(k);
                                if !pstart2.is_new() {
                                    let vtx2 = pstart2.vertex().clone();
                                    if domain.identical(&vtx, &vtx2) {
                                        let thearc2 = pstart2.arc();
                                        let theparam = pstart2.parameter();
                                        let arcorien = domain.orientation_arc(&thearc2);
                                        ispassing = ispassing
                                            && (arcorien == Orientation::Internal
                                                || arcorien == Orientation::External);

                                        let (pt2d, _tg2d) = hcurve2d::d1(thearc2.as_ref(), theparam);
                                        x[0] = pt2d.x;
                                        x[1] = pt2d.y;
                                        ppoint.add_uv(x[0], x[1]);
                                        destination[k] = seqlength as i32 + 1;
                                    }
                                }
                            }
                        }
                    }
                    ppoint.set_passing(ispassing);
                    seqpdep.push(ppoint);
                    seqlength += 1;
                } else if values_ok {
                    // on a un point de depart potentiel
                    let mut vectg = func.direction_3d();
                    let mut dirtg = func.direction_2d();

                    let mut normale = Vec3::ZERO;
                    let mut ptbid = Point3::ZERO;
                    let (mut v1, mut v2) = (Vec3::ZERO, Vec3::ZERO);
                    surf_props::deriv_and_norm(
                        &*surf,
                        x[0],
                        x[1],
                        &mut ptbid,
                        &mut v1,
                        &mut v2,
                        &mut normale,
                    );
                    let mut tg3drst = tg2drst.x * v1 + tg2drst.y * v2;
                    if normale.length_squared() < REAL_EPSILON {
                        //-- Normale Nulle en U/V
                    } else {
                        let mut test = vectg.dot(normale.cross(tg3drst));

                        if pstart.is_new() {
                            let tbis = vectg.normalize().dot(tg3drst.normalize());
                            if tbis.abs() < 1.0 - TOLE {
                                if (test < 0.0 && arcorien == Orientation::Forward)
                                    || (test > 0.0 && arcorien == Orientation::Reversed)
                                {
                                    vectg = -vectg;
                                    dirtg = -dirtg;
                                }
                                ppoint.set_directions(vectg, dirtg);
                            } else {
                                // on garde le point comme point d'arret (tangent)
                                ppoint.set_tangency(true);
                            }
                            ppoint.set_passing(ispassing);
                            destination[i] = seqlength as i32 + 1;
                            seqpdep.push(ppoint);
                            seqlength += 1;
                        } else {
                            // traiter la transition complexe
                            let bidnorm = DVec3::new(1.0, 1.0, 1.0);

                            let mut tobeverified = false;
                            let mut loc_trans;
                            let mut comptrans = CurveTransition::<DVec3>::new();
                            comptrans.reset_3d(vectg, bidnorm, 0.0);
                            let mut tg3drst = tg3drst;
                            let mut test = test;
                            if arcorien != Orientation::Internal
                                && arcorien != Orientation::External
                            {
                                // pour essai
                                let vtx = pstart.vertex().clone();
                                let vtxorien = domain.orientation_vertex(&vtx);
                                test /= vectg.length();
                                test /= normale.cross(tg3drst).length();

                                if test.abs() <= TOLE {
                                    tobeverified = true;
                                    loc_trans = Orientation::External; // et pourquoi pas INTERNAL
                                } else {
                                    if (test > 0.0 && arcorien == Orientation::Forward)
                                        || (test < 0.0 && arcorien == Orientation::Reversed)
                                    {
                                        loc_trans = Orientation::Forward;
                                    } else {
                                        loc_trans = Orientation::Reversed;
                                    }
                                    if arcorien == Orientation::Reversed {
                                        tg3drst = -tg3drst;
                                    } // pas deja fait ???
                                }

                                comptrans.compare(
                                    TOLE,
                                    tg3drst,
                                    bidnorm,
                                    0.0,
                                    loc_trans,
                                    vtxorien,
                                );
                            } else {
                                loc_trans = Orientation::Forward; // (not reached in OCCT)
                            }
                            destination[i] = seqlength as i32 + 1;
                            for k in (i + 1)..=nb_points {
                                if destination[k] == 0 {
                                    let pstart2 = solrst.point(k);
                                    if !pstart2.is_new() {
                                        let vtx2 = pstart2.vertex().clone();
                                        if domain.identical(pstart.vertex(), &vtx2) {
                                            let thearc2 = pstart2.arc();
                                            let theparam = pstart2.parameter();
                                            let arcorien = domain.orientation_arc(&thearc2);

                                            let (pt2d, tg2drst2) =
                                                hcurve2d::d1(thearc2.as_ref(), theparam);
                                            x[0] = pt2d.x;
                                            x[1] = pt2d.y;
                                            ppoint.add_uv(x[0], x[1]);

                                            if arcorien != Orientation::Internal
                                                && arcorien != Orientation::External
                                            {
                                                ispassing = false;
                                                tg3drst = tg2drst2.x * v1 + tg2drst2.y * v2;
                                                test = vectg.dot(normale.cross(tg3drst));
                                                test /= vectg.length();
                                                test /= normale.cross(tg3drst).length();

                                                let vtxorien =
                                                    domain.orientation_vertex(&vtx2);
                                                if test.abs() <= TOLE {
                                                    tobeverified = true;
                                                    loc_trans = Orientation::External; // et pourquoi pas INTERNAL
                                                } else {
                                                    if (test > 0.0
                                                        && arcorien == Orientation::Forward)
                                                        || (test < 0.0
                                                            && arcorien == Orientation::Reversed)
                                                    {
                                                        loc_trans = Orientation::Forward;
                                                    } else {
                                                        loc_trans = Orientation::Reversed;
                                                    }
                                                    if arcorien == Orientation::Reversed {
                                                        tg3drst = -tg3drst;
                                                    } // deja fait????
                                                }

                                                comptrans.compare(
                                                    TOLE,
                                                    tg3drst,
                                                    bidnorm,
                                                    0.0,
                                                    loc_trans,
                                                    vtxorien,
                                                );
                                            }
                                            destination[k] = seqlength as i32 + 1;
                                        }
                                    }
                                }
                            }
                            let mut fairpt = true;
                            if !ispassing {
                                let before = comptrans.state_before();
                                let after = comptrans.state_after();
                                if before == State::Unknown || after == State::Unknown {
                                    fairpt = false;
                                } else if before == State::In {
                                    if after == State::In {
                                        ispassing = true;
                                    } else {
                                        vectg = -vectg;
                                        dirtg = -dirtg;
                                    }
                                } else if after != State::In {
                                    fairpt = false;
                                }
                            }

                            // evite de partir le long d une restriction solution

                            if fairpt && tobeverified {
                                for k in i..=nb_points {
                                    if destination[k] == seqlength as i32 + 1 {
                                        let theparam = solrst.point(k).parameter();
                                        let thearc2 = solrst.point(k).arc();
                                        let arcorien = domain.orientation_arc(&thearc2);

                                        if arcorien == Orientation::Forward
                                            || arcorien == Orientation::Reversed
                                        {
                                            let (pt2d, tg2drst3) =
                                                hcurve2d::d1(thearc2.as_ref(), theparam);
                                            tg3drst = tg2drst3.x * v1 + tg2drst3.y * v2;
                                            let vtxorien = domain
                                                .orientation_vertex(solrst.point(k).vertex());
                                            if (arcorien == Orientation::Forward
                                                && vtxorien == Orientation::Reversed)
                                                || (arcorien == Orientation::Reversed
                                                    && vtxorien == Orientation::Forward)
                                            {
                                                tg3drst = -tg3drst;
                                            }
                                            let test = vectg.normalize().dot(tg3drst.normalize());
                                            if test >= 1.0 - TOLE {
                                                fairpt = false;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }

                            if fairpt {
                                ppoint.set_directions(vectg, dirtg);
                                ppoint.set_passing(ispassing);
                                seqpdep.push(ppoint);
                                seqlength += 1;
                            } else {
                                // il faut remettre en "ordre" si on ne garde pas le point.
                                for k in i..=nb_points {
                                    if destination[k] == seqlength as i32 + 1 {
                                        destination[k] = -destination[k];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/// OCCT ComputeTransitionOnLine (cxx L893-956).
pub(crate) fn compute_transition_on_line(
    sfunc: &mut SurfFunction,
    u: f64,
    v: f64,
    tgline: DVec3,
) -> TypeTrans {
    let surf = sfunc.surface().clone();
    let (pntbid, d1u, d1v) = surf.d1(u, v);
    let _ = pntbid;

    //------------------------------------------------------
    //--   Calcul de la tangente dans l espace uv        ---
    //------------------------------------------------------

    let d1ut = d1u.dot(tgline);
    let d1vt = d1v.dot(tgline);
    let normu2 = d1u.dot(d1u);
    let normv2 = d1v.dot(d1v);
    let d1ud1v = d1u.dot(d1v);
    let det = normu2 * normv2 - d1ud1v * d1ud1v;
    if det < REAL_EPSILON {
        //-- On ne doit pas passer ici !!
        return TypeTrans::Undecided;
    }

    let alpha = (d1ut * normv2 - d1vt * d1ud1v) / det;
    let beta = (normu2 * d1vt - d1ud1v * d1ut) / det;
    //-----------------------------------------------------
    //--  Calcul du Gradient de la fonction Utilisee     --
    //--  pour le contour apparent                       --
    //-----------------------------------------------------

    let x = [u, v];
    let df = sfunc.derivatives(&x).unwrap_or([0.0, 0.0]);
    let v1 = df[0];
    let v2 = df[1];

    //-----------------------------------------------------
    //-- On calcule si la fonction                       --
    //--        F(.) = Normale . Dir_Regard              --
    //-- Croit Losrque l on se deplace sur la Gauche     --
    //--  de la direction de deplacement sur la ligne.   --
    //-----------------------------------------------------

    let det = -v1 * beta + v2 * alpha;

    if det < REAL_EPSILON {
        // revoir le test jag 940620
        return TypeTrans::Undecided;
    }
    if det > 0.0 {
        return TypeTrans::Out;
    }
    TypeTrans::In
}

/// OCCT ProcessSegments (cxx L958-1118).
pub(crate) fn process_segments(
    solrst: &TheSearch,
    slin: &mut Vec<Line>,
    tol_arc: f64,
    sfunc: &mut SurfFunction,
    domain: &mut dyn ContapDomain,
) {
    let nbedg = solrst.nb_segments();
    let surf = sfunc.surface().clone();

    for i in 1..=nbedg {
        let thesegsol = solrst.segment(i);
        let mut theline = Line::new();
        theline.set_value_arc(thesegsol.curve());

        // Traitement des points debut/fin du segment solution.

        let mut dofirst = false;
        let mut dolast = false;
        let (mut procf, mut procl) = (false, false);
        let mut pstartf: Option<ThePathPointOfTheSearch> = None;
        let mut pstartl: Option<ThePathPointOfTheSearch> = None;
        if thesegsol.has_first_point() {
            dofirst = true;
            pstartf = Some(thesegsol.first_point().clone());
        }
        if thesegsol.has_last_point() {
            dolast = true;
            pstartl = Some(thesegsol.last_point().clone());
        }
        let paramf = pstartf.as_ref().map_or(0.0, |p| p.parameter());
        let paraml = pstartl.as_ref().map_or(0.0, |p| p.parameter());

        // determination de la transition
        let u;
        if dofirst && dolast {
            u = (paramf + paraml) / 2.0;
        } else if dofirst {
            u = paramf + 1.0;
        } else if dolast {
            u = paraml - 1.0;
        } else {
            u = 0.0;
        }

        let (p2d, d2d) = hcurve2d::d1(thesegsol.curve().as_ref(), u);
        let (valpt, d1u, d1v) = surf.d1(p2d.x, p2d.y);
        let _ = valpt;
        let tgline = d2d.x * d1u + d2d.y * d1v;
        let tral = compute_transition_on_line(sfunc, p2d.x, p2d.y, tgline);

        theline.set_transition_on_s(tral);

        if dofirst || dolast {
            let nblines = slin.len();
            for j in 1..=nblines {
                let nbpts = slin[j - 1].nb_vertex();
                for k in 1..=nbpts {
                    let ptvtx = slin[j - 1].vertex(k).clone();
                    if dofirst {
                        if let Some(psf) = &pstartf {
                            if ptvtx.value().distance(psf.value()) <= tol_arc {
                                slin[j - 1].vertex_mut(k).set_multiple();
                                let mut ptvtx = ptvtx.clone();
                                ptvtx.set_multiple();
                                ptvtx.set_parameter(paramf);
                                theline.add(ptvtx);
                                procf = true;
                            }
                        }
                    }
                    if dolast {
                        if let Some(psl) = &pstartl {
                            if ptvtx.value().distance(psl.value()) <= tol_arc {
                                slin[j - 1].vertex_mut(k).set_multiple();
                                let mut ptvtx = ptvtx.clone();
                                ptvtx.set_multiple();
                                ptvtx.set_parameter(paraml);
                                theline.add(ptvtx);
                                procl = true;
                            }
                        }
                    }
                }
                // Si on a traite le pt debut et/ou fin, on ne doit pas
                // recommencer si il (ils) correspond(ent) a un point multiple.

                if procf {
                    dofirst = false;
                }
                if procl {
                    dolast = false;
                }
            }
        }

        // Si on n a pas trouve le point debut et./ou fin sur une des lignes
        // d intersection, il faut quand-meme le placer sur la restriction
        // solution

        if dofirst {
            if let Some(psf) = &pstartf {
                let p2d = hcurve2d::value(thesegsol.curve().as_ref(), paramf);
                let mut ptvtx = ContapPoint::new();
                ptvtx.set_value(psf.value(), p2d.x, p2d.y);
                ptvtx.set_parameter(paramf);
                if !psf.is_new() {
                    ptvtx.set_vertex(psf.vertex().clone());
                }
                theline.add(ptvtx);
            }
        }
        if dolast {
            if let Some(psl) = &pstartl {
                let p2d = hcurve2d::value(thesegsol.curve().as_ref(), paraml);
                let mut ptvtx = ContapPoint::new();
                ptvtx.set_value(psl.value(), p2d.x, p2d.y);
                ptvtx.set_parameter(paraml);
                if !psl.is_new() {
                    ptvtx.set_vertex(psl.vertex().clone());
                }
                theline.add(ptvtx);
            }
        }

        // il faut chercher le points internal sur les restrictions solutions.
        if thesegsol.has_first_point() && thesegsol.has_last_point() {
            compute_internal_points_on_rstr(&mut theline, paramf, paraml, sfunc);
        }
        line_constructor(slin, domain, &mut theline, &*surf); //-- lbr
        //-- slin.Append(theline);
        theline.clear();
    }
}

/// OCCT ComputeInternalPointsOnRstr (cxx L1120-1284).
pub(crate) fn compute_internal_points_on_rstr(
    line: &mut Line,
    paramf: f64,
    paraml: f64,
    sfunc: &mut SurfFunction,
) {
    // On recherche les points ou la tangente a la ligne de contour et
    // la direction sont alignees.
    // 1ere etape : recherche de changement de signe.
    // 2eme etape : localisation de la solution par dichotomie

    let surf = sfunc.surface().clone();
    let type_func = sfunc.function_type();

    if line.type_contour() != IType::Restriction {
        return;
    }

    let thearc = line.arc();

    let nbpnts = hcont::nb_samples_on_arc(thearc.as_ref());
    let mut indexinf = 1i32;
    let mut vecregard = sfunc.direction();
    let toler = hcurve2d::resolution(thearc.as_ref(), CONFUSION);
    let mut found = false;

    let mut vecref = DVec3::ZERO;
    let mut paraminf = 0.0;

    loop {
        paraminf = ((nbpnts - indexinf) as f64 * paramf + (indexinf - 1) as f64 * paraml)
            / (nbpnts - 1) as f64;
        let (p2d, d2d) = hcurve2d::d1(thearc.as_ref(), paraminf);
        let (pcour, d1u, d1v) = surf.d1(p2d.x, p2d.y);
        let tgt = d2d.x * d1u + d2d.y * d1v;

        if tgt.length() > GP_RESOLUTION {
            if type_func == TFunction::ContourPrs || type_func == TFunction::DraftPrs {
                vecregard = pcour - sfunc.eye();
            }
            vecref = vecregard.cross(tgt);

            if vecref.length() <= GP_RESOLUTION {
                indexinf += 1;
            } else {
                found = true;
            }
        } else {
            indexinf += 1;
        }
        if !((indexinf <= nbpnts) && !found) {
            break;
        }
    }

    let mut indexsup = indexinf + 1;
    let mut toutvu = indexsup > nbpnts;
    while !toutvu {
        let paramsup = ((nbpnts - indexsup) as f64 * paramf + (indexsup - 1) as f64 * paraml)
            / (nbpnts - 1) as f64;
        let (p2d, d2d) = hcurve2d::d1(thearc.as_ref(), paramsup);
        let (pcour, d1u, d1v) = surf.d1(p2d.x, p2d.y);
        let tgt = d2d.x * d1u + d2d.y * d1v;

        let vectest;
        if tgt.length() > GP_RESOLUTION {
            if type_func == TFunction::ContourPrs || type_func == TFunction::DraftPrs {
                vecregard = pcour - sfunc.eye();
            }
            vectest = vecregard.cross(tgt);
        } else {
            vectest = DVec3::ZERO;
        }
        if vectest.length() <= GP_RESOLUTION {
            // On cherche un vrai changement de signe
            indexsup += 1;
        } else {
            if vectest.dot(vecref) < 0.0 {
                // Essayer de converger
                let mut pinf = paraminf;
                let mut psup = paramsup;
                let mut solution = false;
                let mut ok = false;
                let mut paramp = 0.0;
                while !solution {
                    paramp = (pinf + psup) / 2.0;
                    let (p2d, d2d) = hcurve2d::d1(thearc.as_ref(), paramp);
                    let (pcour, d1u, d1v) = surf.d1(p2d.x, p2d.y);
                    let tgt = d2d.x * d1u + d2d.y * d1v;

                    let vtestb;
                    if tgt.length() > GP_RESOLUTION {
                        if type_func == TFunction::ContourPrs || type_func == TFunction::DraftPrs {
                            vecregard = pcour - sfunc.eye();
                        }
                        vtestb = vecregard.cross(tgt);
                    } else {
                        vtestb = DVec3::ZERO;
                    }

                    if vtestb.length() <= GP_RESOLUTION
                        || (paramp - pinf).abs() <= toler
                        || (paramp - psup).abs() <= toler
                    {
                        // on est a la solution
                        solution = true;
                        ok = true;
                    } else if vtestb.dot(vecref) < 0.0 {
                        psup = paramp;
                    } else {
                        pinf = paramp;
                    }
                }

                if ok {
                    // On verifie que le point trouve ne correspond pas a un ou
                    // des vertex deja existant(s). On teste sur paramp.
                    for i in 1..=line.nb_vertex() {
                        let thevtx = line.vertex_mut(i);
                        if (thevtx.parameter_on_line() - paramp).abs() <= toler {
                            thevtx.set_internal();
                            ok = false; // on a correspondence
                        }
                    }
                    if ok {
                        // il faut alors rajouter le point
                        let mut internalp = ContapPoint::with_uv(pcour, p2d.x, p2d.y);
                        internalp.set_parameter(paramp);
                        internalp.set_internal();
                        line.add(internalp);
                    }
                }
                // cxx L1275 — paramsup re-evaluated for the current indexsup.
                let paramsup = ((nbpnts - indexsup) as f64 * paramf
                    + (indexsup - 1) as f64 * paraml)
                    / (nbpnts - 1) as f64;
                vecref = vectest;
                indexinf = indexsup;
                indexsup += 1;
                paraminf = paramsup;
            } else {
                vecref = vectest;
                indexinf = indexsup;
                indexsup += 1;
                paraminf = paramsup;
            }
        }
        toutvu = indexsup > nbpnts;
    }
}

/// OCCT ComputeInternalPoints (cxx L1286-1539).
pub(crate) fn compute_internal_points(
    line: &mut Line,
    sfunc: &mut SurfFunction,
    ureso: f64,
    vreso: f64,
) {
    // On recherche les points ou la tangente a la ligne de contour et
    // la direction sont alignees.
    // 1ere etape : recherche de changement de signe.
    // 2eme etape : localisation de la solution par simili dichotomie

    let surf = sfunc.surface().clone();
    let type_func = sfunc.function_type();

    if line.type_contour() != IType::Walking {
        return;
    }

    let mut nbpnts = line.nb_pnts() as i32;

    let toler = [ureso, vreso]; //-- Trop long !!! (cxx L1317-1318)
    let infb = [surf.first_u_parameter(), surf.first_v_parameter()];
    let supb = [surf.last_u_parameter(), surf.last_v_parameter()];

    let mut rsnld = FunctionSetRoot::new(sfunc, toler);

    let mut indexinf = 1i32;
    let mut vecregard = sfunc.direction();

    let mut found = false;
    let mut xinf = [0.0f64; 2];
    let mut vecref = DVec3::ZERO;
    loop {
        let (uu, vv) = line.point(indexinf as usize).parameters_on_surface(false);
        xinf = [uu, vv];
        let _ = sfunc.values(&xinf);
        if !sfunc.is_tangent() {
            let tgt = sfunc.direction_3d();
            if type_func == TFunction::ContourPrs || type_func == TFunction::DraftPrs {
                vecregard = line.point(indexinf as usize).value() - sfunc.eye();
            }
            vecref = vecregard.cross(tgt);

            if vecref.length() <= GP_RESOLUTION {
                indexinf += 1;
            } else {
                found = true;
            }
        } else {
            indexinf += 1;
        }
        if !((indexinf <= nbpnts) && !found) {
            break;
        }
    }

    let mut indexsup = indexinf + 1;
    let mut toutvu = indexsup > nbpnts;
    let mut x = [0.0f64; 2];
    let mut xsup = [0.0f64; 2];
    while !toutvu {
        let (uu, vv) = line.point(indexsup as usize).parameters_on_surface(false);
        xsup = [uu, vv];
        let _ = sfunc.values(&xsup);
        let vectest;
        if !sfunc.is_tangent() {
            let tgt = sfunc.direction_3d();

            if type_func == TFunction::ContourPrs || type_func == TFunction::DraftPrs {
                vecregard = line.point(indexsup as usize).value() - sfunc.eye();
            }
            vectest = vecregard.cross(tgt);
        } else {
            vectest = DVec3::ZERO;
        }
        if vectest.length() <= GP_RESOLUTION {
            // On cherche un vrai changement de signe
            indexsup += 1;
        } else {
            if vectest.dot(vecref) < 0.0 {
                // Essayer de converger
                let mut solution = false;
                let mut ok = false;
                while !solution {
                    // Selecting the middle point between XInf and XSup leads
                    // situation, where X values almost do not change. To
                    // prevent this situation, select shifted point instead of
                    // middle.
                    let a_coef = 2.0 / 3.0;
                    x[0] = xinf[0] + a_coef * (xsup[0] - xinf[0]);
                    x[1] = xinf[1] + a_coef * (xsup[1] - xinf[1]);
                    rsnld.perform(sfunc, x, infb, supb);

                    if !rsnld.is_done() {
                        // cxx L1402: std::cout << "Echec recherche internal points"
                        solution = true;
                        ok = false;
                    } else {
                        x = rsnld.root();
                        let _ = sfunc.values(&x);
                        if sfunc.root().abs() <= sfunc.tolerance() {
                            let vtestb;
                            if !sfunc.is_tangent() {
                                let tgt = sfunc.direction_3d();
                                if type_func == TFunction::ContourPrs
                                    || type_func == TFunction::DraftPrs
                                {
                                    vecregard = sfunc.point() - sfunc.eye();
                                }
                                vtestb = vecregard.cross(tgt);
                            } else {
                                vtestb = DVec3::ZERO;
                            }
                            if vtestb.length() <= GP_RESOLUTION
                                || ((x[0] - xinf[0]).abs() <= toler[0]
                                    && (x[1] - xinf[1]).abs() <= toler[1])
                                || ((x[0] - xsup[0]).abs() <= toler[0]
                                    && (x[1] - xsup[1]).abs() <= toler[1])
                            {
                                // on est a la solution
                                solution = true;
                                ok = true;
                            } else if vtestb.dot(vecref) < 0.0 {
                                xsup = x;
                            } else {
                                xinf = x;
                            }
                        } else {
                            // on n est pas sur une solution
                            solution = true;
                            ok = false;
                        }
                    }
                }

                if ok {
                    let mut newpoint = false;
                    let (u0, v0) = line.point(indexinf as usize).parameters_on_surface(false);
                    let mut vinf = DVec2::new(x[0] - u0, x[1] - v0);
                    let mut paramp = 0.0;
                    if vinf.x.abs() <= toler[0] && vinf.y.abs() <= toler[1] {
                        paramp = indexinf as f64;
                    } else {
                        for index in (indexinf + 1)..=indexsup {
                            let (u1, v1) =
                                line.point(index as usize).parameters_on_surface(false);
                            let vsup = DVec2::new(x[0] - u1, x[1] - v1);
                            if vsup.x.abs() <= toler[0] && vsup.y.abs() <= toler[1] {
                                paramp = index as f64;
                                break;
                            } else if vinf.dot(vsup) < 0.0 {
                                // on est entre les 2 points
                                paramp = index as f64;
                                let mut pt2s = PntOn2S::new();
                                pt2s.set_value(sfunc.point(), false, x[0], x[1]);
                                line.line_on_2s_mut().insert_before(index as usize, &pt2s);

                                //-- Il faut decaler les parametres des vertex
                                for v in 1..=line.nb_vertex() {
                                    let vertex = line.vertex_mut(v);
                                    if vertex.parameter_on_line() >= index as f64 {
                                        vertex.set_parameter(vertex.parameter_on_line() + 1.0);
                                    }
                                }

                                nbpnts += 1;
                                indexsup += 1;
                                newpoint = true;
                                break;
                            } else {
                                vinf = vsup;
                            }
                        }
                    }

                    if !newpoint {
                        // on est sur un point de cheminement. On regarde alors
                        // la correspondance avec un vertex existant.
                        newpoint = true;
                        for v in 1..=line.nb_vertex() {
                            let vertex = line.vertex_mut(v);
                            if vertex.parameter_on_line() == paramp {
                                vertex.set_internal();
                                newpoint = false;
                            }
                        }
                    }

                    if newpoint && paramp > 1.0 && paramp < nbpnts as f64 {
                        // on doit creer un nouveau vertex.
                        let mut internalp = ContapPoint::with_uv(sfunc.point(), x[0], x[1]);
                        internalp.set_parameter(paramp);
                        internalp.set_internal();
                        line.add(internalp);
                    }
                }
                let (uu, vv) = line.point(indexsup as usize).parameters_on_surface(false);
                xsup = [uu, vv];
            }
            vecref = vectest;
            indexinf = indexsup;
            indexsup += 1;
            xinf = xsup;
        }
        toutvu = indexsup > nbpnts;
    }
}

/// OCCT FindLine (cxx L1973-2011).
pub(crate) fn find_line(
    line: &Line,
    surf: &dyn SurfaceAdapter,
    pt2d: DVec2,
    ptref: &mut Point3,
    paramin: &mut f64,
    tgmin: &mut DVec3,
    norm: &mut DVec3,
) -> bool {
    let mut pt = Point3::ZERO;
    let mut tg = DVec3::ZERO;

    surf_props::normale(surf, pt2d.x, pt2d.y, ptref, norm);

    let para;
    if line.type_contour() == IType::Lin {
        let lin = line.line();
        para = crate::geomalgo::int_patch::elclib::line_parameter(&lin, *ptref);
        let (pt_v, tg_v) = crate::geomalgo::int_patch::elclib::line_d1(&lin, para);
        pt = pt_v;
        tg = tg_v;
        let _dist = pt.distance(*ptref) + norm.dot(lin.direction).abs();
    } else {
        // Contap_Circle
        let cir = line.circle();
        para = crate::geomalgo::int_patch::elclib::circle_parameter(&cir, *ptref);
        let (pt_v, tg_v) = crate::geomalgo::int_patch::elclib::circle_d1(&cir, para);
        pt = pt_v;
        tg = tg_v;
        let _dist = pt.distance(*ptref) + norm.dot(tg / cir.radius).abs();
    }
    // OCCT: the single (dist < dismin) pass assigns Paramin/ptmin/Tgmin
    // unconditionally (dismin starts at RealLast()).
    *paramin = para;
    *tgmin = tg;
    pt.distance_squared(*ptref) <= super::TOLPETIT
}

/// OCCT PutPointsOnLine (cxx L2013-2091).
pub(crate) fn put_points_on_line(
    solrst: &TheSearch,
    surf: &dyn SurfaceAdapter,
    slin: &mut Vec<Line>,
) {
    let nb_points = solrst.nb_points();

    let nb_lin = slin.len();
    for l in 1..=nb_lin {
        for i in 1..=nb_points {
            let pstart = solrst.point(i);
            let thearc = pstart.arc();
            let theparam = pstart.parameter();

            let (pt2d, d2d) = hcurve2d::d1(thearc.as_ref(), theparam);
            let mut ptonsurf = Point3::ZERO;
            let mut vectg = DVec3::ZERO;
            let mut normale = DVec3::ZERO;
            let mut paramlin = 0.0;
            let goon = find_line(
                &slin[l - 1],
                surf,
                pt2d,
                &mut ptonsurf,
                &mut paramlin,
                &mut vectg,
                &mut normale,
            );

            let mut ppoint = ContapPoint::new();

            if goon {
                let (_bidpt, d1u, d1v) = surf.d1(pt2d.x, pt2d.y);
                ppoint.set_value(ptonsurf, pt2d.x, pt2d.y);
                let mut tline = Transition::new();
                let mut tarc = Transition::new();
                if normale.length() < REAL_EPSILON {
                    // TLine.SetValue(); TArc.SetValue();
                    tline.set_value_undecided();
                    tarc.set_value_undecided();
                } else {
                    // Petit test qui devrait permettre de bien traiter les
                    // pointes des cones, et les sommets d'une sphere.
                    let tgtrst;
                    if d2d.y.abs() <= rcad_kernel::precision::PCONFUSION {
                        let mut t = d1v.cross(normale);
                        if d2d.x < 0.0 {
                            t = -t;
                        }
                        tgtrst = t;
                    } else {
                        tgtrst = d2d.x * d1u + d2d.y * d1v;
                    }
                    make_transition(vectg, tgtrst, normale, &mut tline, &mut tarc);
                }

                ppoint.set_arc(thearc.clone(), theparam, tline, tarc);
                ppoint.set_parameter(paramlin);
                if !pstart.is_new() {
                    ppoint.set_vertex(pstart.vertex().clone());
                }
                slin[l - 1].add(ppoint);
            }
        }
    }
}

/// OCCT ComputeTransitionOngpLine (cxx L2100-2126).
pub(crate) fn compute_transition_on_gp_line(sfunc: &mut SurfFunction, l: &Line3) -> TypeTrans {
    let surf = sfunc.surface().clone();
    let typ_s = surf.get_type();
    let (p, t) = crate::geomalgo::int_patch::elclib::line_d1(l, 0.0);
    let mut u = 0.0;
    let mut v = 0.0;
    match typ_s {
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder => {
            let cy = surf.cylinder();
            let (uu, vv) = rcad_kernel::math::el::elslib_cylinder_parameters(
                p,
                cy.origin,
                cy.ref_dir.normalize_or_zero(),
                cy.y_axis(),
                cy.axis,
                cy.radius,
            );
            u = uu;
            v = vv;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone => {
            let co = surf.cone();
            let (uu, vv) = rcad_kernel::math::el::elslib_cone_parameters(
                p,
                co.apex,
                co.ref_dir.normalize_or_zero(),
                co.axis.cross(co.ref_dir).normalize_or_zero(),
                co.axis,
                co.radius,
                co.half_angle_rad,
            );
            u = uu;
            v = vv;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere => {
            let sp = surf.sphere();
            let (uu, vv) = rcad_kernel::math::el::elslib_sphere_parameters(
                p,
                sp.center,
                sp.ref_dir.normalize(),
                sp.axis.cross(sp.ref_dir).normalize(),
                sp.axis,
            );
            u = uu;
            v = vv;
        }
        _ => {}
    }
    compute_transition_on_line(sfunc, u, v, t)
}

/// OCCT ComputeTransitionOngpCircle (cxx L2128-2154).
pub(crate) fn compute_transition_on_gp_circle(sfunc: &mut SurfFunction, c: &Circle3) -> TypeTrans {
    let surf = sfunc.surface().clone();
    let typ_s = surf.get_type();
    let (p, t) = crate::geomalgo::int_patch::elclib::circle_d1(c, 0.0);
    let mut u = 0.0;
    let mut v = 0.0;
    match typ_s {
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder => {
            let cy = surf.cylinder();
            let (uu, vv) = rcad_kernel::math::el::elslib_cylinder_parameters(
                p,
                cy.origin,
                cy.ref_dir.normalize_or_zero(),
                cy.y_axis(),
                cy.axis,
                cy.radius,
            );
            u = uu;
            v = vv;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone => {
            let co = surf.cone();
            let (uu, vv) = rcad_kernel::math::el::elslib_cone_parameters(
                p,
                co.apex,
                co.ref_dir.normalize_or_zero(),
                co.axis.cross(co.ref_dir).normalize_or_zero(),
                co.axis,
                co.radius,
                co.half_angle_rad,
            );
            u = uu;
            v = vv;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere => {
            let sp = surf.sphere();
            let (uu, vv) = rcad_kernel::math::el::elslib_sphere_parameters(
                p,
                sp.center,
                sp.ref_dir.normalize(),
                sp.axis.cross(sp.ref_dir).normalize(),
                sp.axis,
            );
            u = uu;
            v = vv;
        }
        _ => {}
    }
    compute_transition_on_line(sfunc, u, v, t)
}
