//! OCCT ShapeUpgrade_UnifySameDomain.cxx L4560-4687 — `SplitWire` (the body
//! of the forward declaration at cxx L142-145).

use rcad_kernel::geom::Curve2dEval;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::BRepTool;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType};

use super::statics_a::REAL_LAST;
use super::topexp::{first_vertex, last_vertex, occt_is_same_shape};
use super::{map_add, shape_key, IndexedDataMapOfShapeListOfShape, MapOfShape};
use crate::shhealing::shape_build::brep_tool::iter_subshapes;

/// OCCT static SplitWire (cxx L4560-4687): splits the wire at the splitting
/// vertices into the sequence of wires.
pub fn split_wire(
    brep: &mut BRep,
    the_wire: &Shape,
    the_face: &Shape,
    the_vmap: &MapOfShape,
    the_wire_seq: &mut Vec<Shape>,
) {
    // OCCT L4565: the vertex -> edges map.
    let mut a_vemap: std::collections::HashMap<(u64, u32), (Shape, Vec<Shape>)> =
        std::collections::HashMap::new();

    // OCCT L4567-4592.
    let mut a_emap = MapOfShape::new();
    for itw in iter_subshapes(brep, the_wire, true, true) {
        let an_edge = itw;
        if !map_add(&mut a_emap, &an_edge) {
            continue;
        }
        if an_edge.orientation != Orientation::Forward
            && an_edge.orientation != Orientation::Reversed
        {
            continue;
        }
        let a_vertex = first_vertex(brep, &an_edge, true);
        let key = shape_key(&a_vertex);
        match a_vemap.get_mut(&key) {
            Some((_, list)) => list.push(an_edge),
            None => {
                a_vemap.insert(key, (a_vertex, vec![an_edge]));
            }
        }
    }

    // OCCT L4594-4685.
    let mut a_bb = BRepBuilder::new();
    for ii in 1..=the_vmap.len() {
        let (_, an_origin) = the_vmap.get_index(ii - 1).unwrap();
        let an_origin = an_origin.clone();
        let Some((_, a_branches)) = a_vemap.get_mut(&shape_key(&an_origin)) else {
            continue;
        };
        let mut a_branches = a_branches.clone();
        a_vemap.insert(shape_key(&an_origin), (an_origin.clone(), Vec::new()));
        let mut bi = 0usize;
        while bi < a_branches.len() {
            let mut cur_edge = a_branches[bi].clone();
            a_branches.remove(bi);

            // OCCT L4605-4607.
            let a_new_wire = a_bb.make_wire(brep);
            loop {
                // OCCT L4609.
                a_bb.add_to_wire(brep, a_new_wire.clone(), cur_edge.clone());

                // OCCT L4612-4617.
                let a_vertex = last_vertex(brep, &cur_edge, true);
                if occt_is_same_shape(&a_vertex, &an_origin) {
                    break;
                }
                // OCCT L4619-4628.
                let Some((_, a_elist)) = a_vemap.get_mut(&shape_key(&a_vertex)) else {
                    break;
                };
                if a_elist.is_empty() {
                    break;
                }

                if a_elist.len() == 1 {
                    // OCCT L4630-4634.
                    cur_edge = a_elist[0].clone();
                    a_elist.clear();
                } else {
                    // OCCT L4636-4671: the rightest-direction choice.
                    let Some((cur_pc, fpar, lpar)) = brep.curve_on_surface(&cur_edge, the_face)
                    else {
                        break;
                    };
                    let a_param = if cur_edge.orientation == Orientation::Forward {
                        lpar
                    } else {
                        fpar
                    };
                    let (_a_point, cur_dir_v) =
                        (cur_pc.point_at(a_param), cur_pc.derivative_at(a_param));
                    let mut cur_dir: glam::DVec2 = cur_dir_v.normalize_or_zero();
                    if cur_edge.orientation == Orientation::Reversed {
                        cur_dir = -cur_dir;
                    }
                    let mut min_angle = REAL_LAST;
                    let mut next_edge = Shape::null();
                    for an_edge in a_elist.iter() {
                        let Some((a_pc, fpar, lpar)) = brep.curve_on_surface(an_edge, the_face)
                        else {
                            continue;
                        };
                        let a_param = if an_edge.orientation == Orientation::Forward {
                            fpar
                        } else {
                            lpar
                        };
                        let a_dir_v = a_pc.derivative_at(a_param);
                        let mut a_dir: glam::DVec2 = a_dir_v.normalize_or_zero();
                        if an_edge.orientation == Orientation::Reversed {
                            a_dir = -a_dir;
                        }
                        let an_angle = super::gap_deps::dir_angle_2d(cur_dir, a_dir);
                        if an_angle < min_angle {
                            min_angle = an_angle;
                            next_edge = an_edge.clone();
                        }
                    }
                    cur_edge = next_edge;
                    // OCCT L4673-4681: remove CurEdge from the list.
                    if let Some((_, a_elist)) = a_vemap.get_mut(&shape_key(&a_vertex)) {
                        for k in 0..a_elist.len() {
                            if occt_is_same_shape(&cur_edge, &a_elist[k]) {
                                a_elist.remove(k);
                                break;
                            }
                        }
                    }
                } // else (more than one edge)
            } // for (;;)
            the_wire_seq.push(a_new_wire);
        }
        let _ = a_branches;
    }
    let _ = (ShapeType::Edge, IndexedDataMapOfShapeListOfShape::new());
}
