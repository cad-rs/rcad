//! OCCT BRepTools_Substitution (TKBRep/BRepTools; TKTopAlgo landing zone per
//! docs/module-map.md) — 1:1 translation of BRepTools_Substitution.hxx
//! (L24-77) + BRepTools_Substitution.cxx (L27-176).
//!
//! Architecture differences:
//! - `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
//!   maps to `HashMap<ShapeKey, Vec<Shape>>` (TopTools_ShapeMapHasher).
//! - `BRep_Builder` maps to the pool-leading `rcad_kernel::topo::topods::
//!   BRepBuilder` + `BRep` pair (the owning pool is the leading `brep`
//!   argument, the same convention as the BRepFill sweep part A).
//! - OCCT `TopoDS_Shape::Nullify()` maps to `Shape::null()` (the OCCT null
//!   handle; the shape disappears from its ancestors through the empty
//!   replacement list).
//! - `TopoDS_Iterator` maps to the local [`topods_iterator`] children walk
//!   (the direct sub-shapes stored by each TShape kind).
//! - `NewS.EmptyCopy()` maps to the pool `BRep::empty_copy` (a fresh TShape
//!   with the same kind and data, no children).

use std::collections::HashMap;

use rcad_kernel::topo::topods::{Orientation, Shape, TShape};

use crate::brep_fill::generator::ShapeKey;

/// OCCT TopoDS_Iterator — the direct sub-shapes of a shape (BRep_Builder
/// storage order: Edge -> [first, last] vertices; Wire -> edges; Face ->
/// [outer, inner...]; Shell -> faces; Solid -> shells; CompSolid -> solids;
/// Compound -> shapes).
fn topods_iterator(s: &Shape) -> Vec<Shape> {
    if s.is_null() {
        return Vec::new();
    }
    match s.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => {
            let mut out = Vec::new();
            if !ed.first.is_null() {
                out.push(ed.first.clone());
            }
            if !ed.last.is_null() {
                out.push(ed.last.clone());
            }
            out
        }
        TShape::Wire(wd) => wd.edges.clone(),
        TShape::Face(fd) => {
            let mut out = Vec::new();
            if !fd.outer_wire.is_null() {
                out.push(fd.outer_wire.clone());
            }
            out.extend(fd.inner_wires.iter().cloned());
            out
        }
        TShape::Shell(sd) => sd.faces.clone(),
        TShape::Solid(sd) => sd.shells.clone(),
        TShape::CompSolid(sd) => sd.clone(),
        TShape::Compound(cd) => cd.clone(),
    }
}

/// OCCT BRep_Builder::Add(S, SS) over the pool — the per-kind container
/// append (BRep_Builder.cxx Add).
fn builder_add(brep: &mut BRep, s: &Shape, ss: Shape) {
    match s.data.as_ref() {
        TShape::Edge(_) => {
            rcad_kernel::topo::topods::BRepBuilder::new().add_to_edge(brep, s.clone(), ss)
        }
        TShape::Wire(_) => {
            rcad_kernel::topo::topods::BRepBuilder::new().add_to_wire(brep, s.clone(), ss)
        }
        TShape::Face(_) => {
            rcad_kernel::topo::topods::BRepBuilder::new().add_to_face(brep, s.clone(), ss)
        }
        TShape::Shell(_) => {
            rcad_kernel::topo::topods::BRepBuilder::new().add_to_shell(brep, s.clone(), ss)
        }
        TShape::Solid(_) => {
            // OCCT B.Add(solid, shell) — the shell list is appended in place
            // (same TShape identity, the OCCT builder edits the TShape).
            brep.solid_mut(s.clone()).shells.push(ss);
        }
        TShape::CompSolid(_) | TShape::Compound(_) => {
            rcad_kernel::topo::topods::BRepBuilder::new().add_to_compound(brep, s.clone(), ss)
        }
        TShape::Vertex(_) => {}
    }
}

use rcad_kernel::topo::topods::BRep;

/// OCCT BRepTools_Substitution (hxx L40-76): a tool to substitute subshapes
/// by other shapes.
pub struct BRepToolsSubstitution {
    /// OCCT myMap (NCollection_DataMap<TopoDS_Shape,
    /// NCollection_List<TopoDS_Shape>>).
    my_map: HashMap<ShapeKey, Vec<Shape>>,
}

impl Default for BRepToolsSubstitution {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepToolsSubstitution {
    /// OCCT BRepTools_Substitution() (cxx L31).
    pub fn new() -> Self {
        BRepToolsSubstitution {
            my_map: HashMap::new(),
        }
    }

    /// OCCT Clear (cxx L35-38) — reset all the fields.
    pub fn clear(&mut self) {
        self.my_map.clear();
    }

    /// OCCT Substitute (cxx L42-47) — <OldShape> will be replaced by
    /// <NewShapes>.
    pub fn substitute(&mut self, os: &Shape, ns: Vec<Shape>) {
        // Standard_ConstructionError_Raise_if(IsCopied(OS), ...)
        if self.is_copied(os) {
            panic!("BRepTools_CutClue::Substitute");
        }
        self.my_map.insert(ShapeKey(os.ptr_id()), ns);
    }

    /// OCCT Build (cxx L53-141) — build the new shape from S if its
    /// subshapes have been modified.
    pub fn build(&mut self, brep: &mut BRep, s: &Shape) {
        if self.is_copied(s) {
            return;
        }

        // TopoDS_Iterator iteS(S.Oriented(TopAbs_FORWARD));
        let ite_s = topods_iterator(s);
        let mut is_modified = false;
        let mut has_sub_shape = false;

        // ------------------------------------------
        // look S is modified and build subshapes.
        // ------------------------------------------
        for ss in &ite_s {
            self.build(brep, ss);
            if self.is_copied(ss) {
                is_modified = true;
            }
        }

        // TopoDS_Shape NewS = S.Oriented(TopAbs_FORWARD);
        let mut new_s = {
            let mut ns = s.clone();
            ns.orientation = Orientation::Forward;
            ns
        };
        if is_modified {
            // ----------------------------------------
            // Rebuild S.
            // ----------------------------------------
            // NewS.EmptyCopy();
            new_s = brep.empty_copy(s.clone());
            new_s.orientation = Orientation::Forward;

            if matches!(new_s.data.as_ref(), TShape::Edge(_)) {
                // BRep_Tool::Range(TopoDS::Edge(S), f, l);
                // B.Range(TopoDS::Edge(NewS), f, l);
                let range = brep.edge(s.clone()).range;
                rcad_kernel::topo::topods::BRepBuilder::new()
                    .set_edge_range(brep, new_s.clone(), range[0], range[1]);
            }

            // ------------------------------------------
            // Add the copy of subshapes of S to NewS.
            // ------------------------------------------
            for ss in &ite_s {
                // TopAbs_Orientation OS = iteS.Value().Orientation();
                let os = ss.orientation;
                // NCollection_List<TopoDS_Shape> L = myMap(iteS.Value());
                let l = self
                    .my_map
                    .get(&ShapeKey(ss.ptr_id()))
                    .cloned()
                    .unwrap_or_default();

                for nss in l {
                    // ------------------------------------------
                    // Rebuild NSS and add its copy to NewS.
                    // ------------------------------------------
                    self.build(brep, &nss);

                    // const NCollection_List<TopoDS_Shape>& NL = myMap(NSS);
                    let nl = self
                        .my_map
                        .get(&ShapeKey(nss.ptr_id()))
                        .cloned()
                        .unwrap_or_default();
                    // TopAbs_Orientation NewOr = TopAbs::Compose(OS,
                    // NSS.Orientation());
                    let new_or = os.compose(nss.orientation);

                    for value in nl {
                        // B.Add(NewS, iteNL.Value().Oriented(NewOr));
                        let mut oriented = value;
                        oriented.orientation = new_or;
                        builder_add(brep, &new_s, oriented);
                        has_sub_shape = true;
                    }
                }
            }
            if !has_sub_shape {
                if matches!(
                    new_s.data.as_ref(),
                    TShape::Wire(_) | TShape::Shell(_) | TShape::Solid(_) | TShape::Compound(_)
                ) {
                    // -----------------------------------------------------------------
                    // Wire,Solid,Shell,Compound must have subshape else they
                    // disappear
                    // -----------------------------------------------------------------
                    new_s = Shape::null();
                }
            }
        }
        // -------------------------------------------------------
        // NewS has the same orientation than S in its ancestors
        // so NewS is bound with orientation FORWARD.
        // -------------------------------------------------------
        let mut l: Vec<Shape> = Vec::new();
        if !new_s.is_null() {
            let mut forward = new_s;
            forward.orientation = Orientation::Forward;
            l.push(forward);
        }
        self.substitute(s, l);
    }

    /// OCCT IsCopied (cxx L147-164) — returns true if S has been replaced.
    pub fn is_copied(&self, s: &Shape) -> bool {
        if let Some(list) = self.my_map.get(&ShapeKey(s.ptr_id())) {
            if list.is_empty() {
                true
            } else {
                !s.is_same(&list[0])
            }
        } else {
            false
        }
    }

    /// OCCT Copy (cxx L170-175) — the set of shapes substituted to S.
    pub fn copy(&self, s: &Shape) -> Vec<Shape> {
        if !self.is_copied(s) {
            panic!("BRepTools_Substitution::Copy");
        }
        self.my_map
            .get(&ShapeKey(s.ptr_id()))
            .cloned()
            .unwrap_or_default()
    }
}
