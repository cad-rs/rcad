//! OCCT BRepBuilderAPI_Sewing.cxx — the output group:
//! `CreateSewedShape` (L5030-5255), `CreateOutputInformations`
//! (L5257-5365) and `SameParameterShape` (L5882-5910).

use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;

use super::{idx_list_get, idx_shape_get, set_add, set_contains};
use super::BRepBuilderAPISewing;

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::CreateSewedShape() (cxx L5030-5255) —
    /// creates the sewed shape through the BRepTools_Quilt.
    pub(crate) fn create_sewed_shape(&mut self, brep: &mut BRep) {
        // ---------------------
        // create the new shapes
        // ---------------------
        // OCCT L5033: BRepTools_Quilt aQuilt;
        let mut a_quilt = crate::topalgo::brep_tools_quilt::BRepToolsQuilt::new();
        // OCCT L5034: bool isLocal = !myShape.IsNull();
        let is_local = !self.my_shape.is_null();
        if is_local {
            // Local sewing
            // OCCT L5037-5039.
            let ns = self.my_re_shape.apply(brep, &self.my_shape, ShapeType::Shape);
            a_quilt.add(brep, &ns);
        }
        // OCCT L5040-5052.
        for i in 0..self.my_old_shapes.len() {
            let (k, mut sh) = super::idx_shape_get(&self.my_old_shapes, i);
            let _ = &k;
            if !sh.is_null() {
                sh = self.my_re_shape.apply(brep, &sh, ShapeType::Shape);
                if let Some(entry) = self.my_old_shapes.get_index_mut(i) {
                    entry.1.1 = sh.clone();
                }
                if !is_local {
                    a_quilt.add(brep, &sh);
                }
            }
        }
        // OCCT L5053-5054.
        let a_new_shape = a_quilt.shells(brep);
        let mut numsh = 0i32;

        // OCCT L5056-5062.
        let mut old_shells: super::ShapeSet = super::ShapeSet::new();

        let mut a_b = BRepBuilder::new();
        let mut a_comp = a_b.make_compound(brep, vec![]);
        for a_exp_sh in bat::sub_shapes(&bat::oriented(&a_new_shape, Orientation::Forward)) {
            let mut sh = a_exp_sh;
            let mut has_edges = false;
            if sh.shape_type() == ShapeType::Shell {
                if self.my_nonmanifold {
                    has_edges = !set_contains(&old_shells, &sh);
                } else {
                    // OCCT L5072-5086.
                    let mut face = Shape::null();
                    let mut numf = 0i32;
                    for a_exp_f in
                        bat::explorer(&sh, ShapeType::Face, ShapeType::Shape)
                    {
                        if numf >= 2 {
                            break;
                        }
                        face = a_exp_f;
                        numf += 1;
                    }
                    if numf == 1 {
                        a_b.add_to_compound(brep, a_comp.clone(), face.clone());
                    } else if numf > 1 {
                        a_b.add_to_compound(brep, a_comp.clone(), sh.clone());
                    }
                    if numf != 0 {
                        numsh += 1;
                    }
                }
            } else if sh.shape_type() == ShapeType::Face {
                if self.my_nonmanifold {
                    // OCCT L5090-5097.
                    let ss = a_b.make_shell(brep);
                    a_b.add_to_shell(brep, ss.clone(), sh.clone());
                    sh = ss;
                    has_edges = true;
                } else {
                    // OCCT L5099-5103.
                    a_b.add_to_compound(brep, a_comp.clone(), sh.clone());
                    numsh += 1;
                }
            } else {
                // OCCT L5105-5108.
                a_b.add_to_compound(brep, a_comp.clone(), sh.clone());
                numsh += 1;
            }
            // OCCT L5110-5113.
            if has_edges {
                set_add(&mut old_shells, &sh);
            }
        }
        // Process collected shells
        // OCCT L5116-5236.
        if self.my_nonmanifold {
            let nb_old_shells = old_shells.len();
            if nb_old_shells == 1 {
                // Single shell - check for single face
                // OCCT L5120-5141.
                let sh = super::set_get(&old_shells, 0);
                let mut face = Shape::null();
                let mut numf = 0i32;
                for a_exp_f in bat::explorer(&sh, ShapeType::Face, ShapeType::Shape) {
                    if numf >= 2 {
                        break;
                    }
                    face = a_exp_f;
                    numf += 1;
                }
                if numf == 1 {
                    a_b.add_to_compound(brep, a_comp.clone(), face.clone());
                } else if numf > 1 {
                    a_b.add_to_compound(brep, a_comp.clone(), sh.clone());
                }
                if numf != 0 {
                    numsh += 1;
                }
            } else if nb_old_shells != 0 {
                // Several shells should be merged
                // OCCT L5144-5222.
                let mut index_merged: std::collections::HashSet<usize> =
                    std::collections::HashSet::new();
                while index_merged.len() < nb_old_shells {
                    let mut new_shell = Shape::null();
                    let mut new_edges: super::ShapeSet = super::ShapeSet::new();
                    for i in 1..=nb_old_shells {
                        if index_merged.contains(&i) {
                            continue;
                        }
                        let shell = super::set_get(&old_shells, i - 1);
                        if new_shell.is_null() {
                            // OCCT L5158-5169.
                            let mut a_b2 = BRepBuilder::new();
                            let ns = a_b2.make_shell(brep);
                            for a_it_ss in bat::sub_shapes(&shell) {
                                a_b2.add_to_shell(brep, ns.clone(), a_it_ss.clone());
                            }
                            new_shell = ns;
                            // Fill map of edges
                            for eexp in bat::explorer(&shell, ShapeType::Edge, ShapeType::Shape)
                            {
                                let edge = eexp;
                                set_add(&mut new_edges, &edge);
                            }
                            index_merged.insert(i);
                        } else {
                            // OCCT L5172-5199.
                            let mut has_shared_edge = false;
                            for eexp in
                                bat::explorer(&shell, ShapeType::Edge, ShapeType::Shape)
                            {
                                if has_shared_edge {
                                    break;
                                }
                                has_shared_edge = set_contains(&new_edges, &eexp);
                            }
                            if has_shared_edge {
                                // Add edges to the map
                                for eexp1 in
                                    bat::explorer(&shell, ShapeType::Edge, ShapeType::Shape)
                                {
                                    let edge = eexp1;
                                    set_add(&mut new_edges, &edge);
                                }
                                // Add faces to the shell
                                for fexp in
                                    bat::explorer(&shell, ShapeType::Face, ShapeType::Shape)
                                {
                                    let face = fexp;
                                    let mut a_b2 = BRepBuilder::new();
                                    a_b2.add_to_shell(brep, new_shell.clone(), face.clone());
                                }
                                index_merged.insert(i);
                            }
                        }
                    }
                    // Process new shell
                    // OCCT L5201-5221.
                    let mut face = Shape::null();
                    let mut numf = 0i32;
                    for a_exp_f in
                        bat::explorer(&new_shell, ShapeType::Face, ShapeType::Shape)
                    {
                        if numf >= 2 {
                            break;
                        }
                        face = a_exp_f;
                        numf += 1;
                    }
                    if numf == 1 {
                        a_b.add_to_compound(brep, a_comp.clone(), face.clone());
                    } else if numf > 1 {
                        a_b.add_to_compound(brep, a_comp.clone(), new_shell.clone());
                    }
                    if numf != 0 {
                        numsh += 1;
                    }
                }
            }
        }
        // OCCT L5237-5253.
        if numsh == 1 {
            // Extract single component
            // OCCT L5240-5242: TopoDS_Iterator aIt(aComp, false).
            let comps = super::iter_no_cumori(&a_comp);
            if let Some(first) = comps.first() {
                self.my_sewed_shape = first.clone();
            } else {
                self.my_sewed_shape = a_comp.clone();
            }
        } else {
            self.my_sewed_shape = a_comp.clone();
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::CreateOutputInformations() (cxx
    /// L5257-5365) — constructs the free/contiguous/multiple edge maps.
    pub(crate) fn create_output_informations(&mut self, brep: &mut BRep) {
        // Construct edgeSections
        // OCCT L5262-5264.
        let mut edge_sections: super::IdxListMap = super::IdxListMap::new();
        // (the IndexedMap carrier regulates the free edges)
        for i in 0..self.my_bound_faces.len() {
            let (bound, _) = idx_list_get(&self.my_bound_faces, i);
            // OCCT L5265-5268.
            let lsect: Vec<Shape> = self
                .my_bound_sections
                .get(&bat::shape_key(&bound))
                .cloned()
                .unwrap_or_default();
            // OCCT L5269-5271.
            let applied = self.my_re_shape.apply(brep, &bound, ShapeType::Shape);
            for a_exp in bat::explorer(&applied, ShapeType::Edge, ShapeType::Shape) {
                let edge = a_exp;
                // OCCT L5272-5281.
                let mut sec = bound.clone();
                for a_i in &lsect {
                    let section = a_i.clone();
                    let applied_section =
                        self.my_re_shape.apply(brep, &section, ShapeType::Shape);
                    if edge.is_same(&applied_section) {
                        sec = section;
                        break;
                    }
                }
                // OCCT L5282-5292.
                if edge_sections.contains_key(&bat::shape_key(&edge)) {
                    if let Some(entry) = edge_sections.get_mut(&bat::shape_key(&edge)) {
                        entry.1.push(sec);
                    }
                } else {
                    let list_sec = vec![sec];
                    super::idx_list_add(&mut edge_sections, &edge, list_sec);
                }
            }
        }

        // Fill maps of Free, Contiguous and Multiple edges
        // OCCT L5295-5317.
        for i in 0..edge_sections.len() {
            let (edge, list_section) = idx_list_get(&edge_sections, i);
            if list_section.len() == 1 {
                if super::brep_tool_degenerated(&edge) {
                    set_add(&mut self.my_degenerated, &edge);
                } else {
                    set_add(&mut self.my_free_edges, &edge);
                }
            } else if list_section.len() == 2 {
                super::idx_list_add(&mut self.my_contigous_edges, &edge, list_section);
            } else {
                set_add(&mut self.my_multiple_edges, &edge);
            }
        }

        // constructs myContigSectBound
        // OCCT L5320-5351.
        let a_edge_map: super::DataListMap = super::DataListMap::new(); // gka
        let _ = &a_edge_map;
        for i in 0..self.my_bound_faces.len() {
            let (bound, _) = idx_list_get(&self.my_bound_faces, i);
            if let Some(sections) = self.my_bound_sections.get(&bat::shape_key(&bound)).cloned() {
                for iter in &sections {
                    let section = iter.clone();
                    if !set_contains(&self.my_merged_edges, &section) {
                        continue;
                    }
                    // OCCT L5335-5341.
                    let nedge = self.my_re_shape.apply(brep, &section, ShapeType::Shape);
                    if nedge.is_null() {
                        continue; // szv debug
                    }
                    if !bound.is_same(&section)
                        && self.my_contigous_edges.contains_key(&bat::shape_key(&nedge))
                    {
                        self.my_contig_sec_bound
                            .insert(bat::shape_key(&section), bound.clone());
                    }
                }
            }
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::SameParameterShape() (cxx L5882-5910) —
    /// makes all edges from the sewed shape same parameter if
    /// mySameParameterMode.
    pub(crate) fn same_parameter_shape(&mut self, brep: &mut BRep) {
        // OCCT L5883-5885.
        if !self.my_same_parameter_mode {
            return;
        }
        // Le flag sameparameter est a false pour chaque edge cousue
        // OCCT L5886-5893.
        for exp in bat::explorer(&self.my_sewed_shape, ShapeType::Edge, ShapeType::Shape) {
            let sec = exp;
            // OCCT L5894-5906: BRepLib::SameParameter(sec, Tol(sec)) inside
            // try/catch — the rcad carrier does not throw.  The sewed-shape
            // edges are pool-resident in `brep` (my_sewed_shape was created
            // by create_sewed_shape over the same pool), so the engine's
            // in-place TShape writes reach the sec handle the way the OCCT
            // engine mutates the shared TShape through the handle.
            crate::topalgo::brep_lib::same_parameter::same_parameter(
                brep,
                &sec,
                bat::brep_tool_tolerance(&sec),
            );
        }
    }
}

// The unused import anchors (parity with the OCCT include set).
#[allow(unused_imports)]
use idx_shape_get as _idx_shape_get_anchor;
