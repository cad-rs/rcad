//! OCCT AdvApp2Var_Framework (AdvApp2Var_Framework.hxx +
//! AdvApp2Var_Framework.cxx).
//!
//! Encoding notes (architecture):
//! - `NCollection_Sequence<occ::handle<AdvApp2Var_Node>>` ->
//!   `Vec<Node>`; `NCollection_Sequence<NCollection_Sequence<handle<Iso>>>`
//!   -> `Vec<Vec<Iso>>` (1-based Value(i) -> vec[i - 1], Append -> push,
//!   InsertBefore -> vec.insert(i - 1, ..), InsertAfter -> vec.insert(i, ..)).
//! - The shared-handle `ChangeIso` write-back and the `*Node(ind) = N1`
//!   assignment are reproduced by writing the mutated local value back into
//!   the Vec slot (the only readers are FirstNotApprox / accessors, so the
//!   value-copy-back is behaviorally identical to the OCCT handle store).

use rcad_kernel::base::proj_lib::proj_lib_projected_curve::IsoType;

use super::iso::Iso;
use super::node::Node;
use glam::DVec2;

/// OCCT findStripByBounds (AdvApp2Var_Framework.cxx L31-44) - anonymous
/// namespace helper.
fn find_strip_by_bounds(the_constraints: &[Vec<Iso>], the_first: f64, the_last: f64) -> i32 {
    let mut a_strip_index: i32 = 1;
    while a_strip_index < the_constraints.len() as i32
        && (the_constraints[(a_strip_index - 1) as usize].first().expect("empty strip").t0()
            != the_first
            || the_constraints[(a_strip_index - 1) as usize]
                .first()
                .expect("empty strip")
                .t1()
                != the_last)
    {
        a_strip_index += 1;
    }
    a_strip_index
}

/// OCCT findIsoByConstant (AdvApp2Var_Framework.cxx L46-56) - anonymous
/// namespace helper.
fn find_iso_by_constant(the_strip: &[Iso], the_nb_iso: i32, the_const: f64) -> i32 {
    let mut an_iso_index: i32 = 1;
    while an_iso_index <= the_nb_iso
        && the_strip[(an_iso_index - 1) as usize].constante() != the_const
    {
        an_iso_index += 1;
    }
    an_iso_index
}

/// OCCT AdvApp2Var_Framework (AdvApp2Var_Framework.hxx L36-95).
#[derive(Debug, Clone)]
pub struct Framework {
    /// hxx L92: myNodeConstraints.
    my_node_constraints: Vec<Node>,
    /// hxx L93: myUConstraints.
    my_u_constraints: Vec<Vec<Iso>>,
    /// hxx L94: myVConstraints.
    my_v_constraints: Vec<Vec<Iso>>,
}

impl Framework {
    /// OCCT AdvApp2Var_Framework() (AdvApp2Var_Framework.cxx L61).
    pub fn new() -> Self {
        Framework {
            my_node_constraints: Vec::new(),
            my_u_constraints: Vec::new(),
            my_v_constraints: Vec::new(),
        }
    }

    /// OCCT AdvApp2Var_Framework(Frame, UFrontier, VFrontier)
    /// (AdvApp2Var_Framework.cxx L65-73).
    pub fn new_from(frame: Vec<Node>, u_frontier: Vec<Vec<Iso>>, v_frontier: Vec<Vec<Iso>>) -> Self {
        Framework {
            my_node_constraints: frame,
            my_u_constraints: u_frontier,
            my_v_constraints: v_frontier,
        }
    }

    /// OCCT FirstNotApprox(IndexIso, IndexStrip) (AdvApp2Var_Framework.cxx
    /// L77-107) - search the Index of the first Iso not approximated, if all
    /// Isos are approximated None is returned.  The returned handle is a
    /// value copy of the stored Iso (written back through change_iso()).
    pub fn first_not_approx(&self, index_iso: &mut i32, index_strip: &mut i32) -> Option<Iso> {
        for an_uv_iter in 0..2 {
            let a_seq: &Vec<Vec<Iso>> = if an_uv_iter == 0 {
                &self.my_u_constraints
            } else {
                &self.my_v_constraints
            };
            let mut i: i32 = 1;
            for s in a_seq.iter() {
                let mut j: i32 = 1;
                for an_iso in s.iter() {
                    if !an_iso.is_approximated() {
                        *index_iso = j;
                        *index_strip = i;
                        return Some(an_iso.clone());
                    }
                    j += 1;
                }
                i += 1;
            }
        }
        None // OCCT: return occ::handle<AdvApp2Var_Iso>() (null handle).
    }

    /// OCCT FirstNode(Type, IndexIso, IndexStrip) (AdvApp2Var_Framework.cxx
    /// L111-121).
    pub fn first_node(&self, the_type: IsoType, index_iso: i32, index_strip: i32) -> i32 {
        let a_nb_iso_in_v = self.my_u_constraints.len() as i32 + 1;
        if the_type == IsoType::IsoU {
            return a_nb_iso_in_v * (index_strip - 1) + index_iso;
        }
        a_nb_iso_in_v * (index_iso - 1) + index_strip
    }

    /// OCCT LastNode(Type, IndexIso, IndexStrip) (AdvApp2Var_Framework.cxx
    /// L125-135).
    pub fn last_node(&self, the_type: IsoType, index_iso: i32, index_strip: i32) -> i32 {
        let a_nb_iso_in_v = self.my_u_constraints.len() as i32 + 1;
        if the_type == IsoType::IsoU {
            return a_nb_iso_in_v * index_strip + index_iso;
        }
        a_nb_iso_in_v * (index_iso - 1) + index_strip + 1
    }

    /// OCCT ChangeIso(IndexIso, IndexStrip, theIso) (AdvApp2Var_Framework.cxx
    /// L139-147) - stores theIso in the selected strip (the handle store
    /// becomes a slot write).
    pub fn change_iso(&mut self, index_iso: i32, index_strip: i32, the_iso: Iso) {
        let a_strip = if the_iso.type_() == IsoType::IsoV {
            &mut self.my_u_constraints
        } else {
            &mut self.my_v_constraints
        };
        a_strip[(index_strip - 1) as usize][(index_iso - 1) as usize] = the_iso;
    }

    /// OCCT Node(IndexNode) (hxx L64-67) - const handle returned as a
    /// reference to the stored node.
    pub fn node_index(&self, index_node: i32) -> &Node {
        &self.my_node_constraints[(index_node - 1) as usize]
    }

    /// OCCT Node(IndexNode) mutable slot access (the OCCT
    /// `*myConstraints.Node(indN1) = N1` assignment target).
    pub fn node_index_mut(&mut self, index_node: i32) -> &mut Node {
        &mut self.my_node_constraints[(index_node - 1) as usize]
    }

    /// OCCT Node(U, V) (AdvApp2Var_Framework.cxx L151-164).
    pub fn node(&self, u: f64, v: f64) -> &Node {
        for a_node in self.my_node_constraints.iter() {
            if a_node.coord().x == u && a_node.coord().y == v {
                return a_node;
            }
        }
        // OCCT: return myNodeConstraints.Last().
        self.my_node_constraints.last().expect("empty node constraints")
    }

    /// OCCT IsoU(U, V0, V1) (AdvApp2Var_Framework.cxx L168-177).
    pub fn iso_u(&self, u: f64, v0: f64, v1: f64) -> &Iso {
        let a_strip_index = find_strip_by_bounds(&self.my_v_constraints, v0, v1);
        let a_strip = &self.my_v_constraints[(a_strip_index - 1) as usize];
        let a_iso_index = find_iso_by_constant(a_strip, self.my_u_constraints.len() as i32, u);
        &a_strip[(a_iso_index - 1) as usize]
    }

    /// OCCT IsoV(U0, U1, V) (AdvApp2Var_Framework.cxx L181-190).
    pub fn iso_v(&self, u0: f64, u1: f64, v: f64) -> &Iso {
        let a_strip_index = find_strip_by_bounds(&self.my_u_constraints, u0, u1);
        let a_strip = &self.my_u_constraints[(a_strip_index - 1) as usize];
        let a_iso_index = find_iso_by_constant(a_strip, self.my_v_constraints.len() as i32, v);
        &a_strip[(a_iso_index - 1) as usize]
    }

    /// OCCT UpdateInU(CuttingValue) (AdvApp2Var_Framework.cxx L194-288).
    pub fn update_in_u(&mut self, cutting_value: f64) {
        let mut a_u_strip_index: i32 = 1;
        for a_u_const in self.my_u_constraints.iter() {
            let first = a_u_const.first().expect("empty U strip");
            if first.u0() <= cutting_value && first.u1() >= cutting_value {
                break;
            }
            a_u_strip_index += 1;
        }

        {
            // The OCCT loop mutates the shared handles of the strip in
            // place; rcad works on a value copy of the strip and writes it
            // back (same final state).
            let mut a_strip = self.my_u_constraints[(a_u_strip_index - 1) as usize].clone();
            let a_u_first = a_strip.first().expect("empty U strip").u0();
            let a_u_last = a_strip.first().expect("empty U strip").u1();

            // Modify the V isos of the U strip at aUStripIndex.
            for k in 0..a_strip.len() {
                a_strip[k].change_domain(a_u_first, cutting_value);
                a_strip[k].reset_approx();
            }

            // Insert a new U strip after aUStripIndex.
            let mut a_new_strip: Vec<Iso> = Vec::new();
            for an_iso in a_strip.iter() {
                let mut a_new_iso = Iso::new_full(
                    an_iso.type_(),
                    an_iso.constante(),
                    cutting_value,
                    a_u_last,
                    an_iso.v0(),
                    an_iso.v1(),
                    0,
                    an_iso.u_order(),
                    an_iso.v_order(),
                );
                a_new_iso.reset_approx();
                a_new_strip.push(a_new_iso);
            }
            // OCCT: myUConstraints.InsertAfter(aUStripIndex, aNewStrip);
            // and the modified strip is stored back.
            self.my_u_constraints.insert(a_u_strip_index as usize, a_new_strip);
            self.my_u_constraints[(a_u_strip_index - 1) as usize] = a_strip;
        }

        //  Insert a new Iso U=U* in each V strip after aUStripIndex
        //  and restrict the domains of the adjacent Isos
        for a_v_strip_index in 1..=(self.my_v_constraints.len() as i32) {
            let a_strip = &mut self.my_v_constraints[(a_v_strip_index - 1) as usize];
            let mut an_iso = a_strip[(a_u_strip_index - 1) as usize].clone();
            an_iso.change_domain_full(an_iso.u0(), cutting_value, an_iso.v0(), an_iso.v1());
            a_strip[(a_u_strip_index - 1) as usize] = an_iso.clone();

            let mut a_new_iso = Iso::new_full(
                an_iso.type_(),
                cutting_value,
                an_iso.u0(),
                cutting_value,
                an_iso.v0(),
                an_iso.v1(),
                0,
                an_iso.u_order(),
                an_iso.v_order(),
            );
            a_new_iso.reset_approx();
            // OCCT: aStrip.InsertAfter(aUStripIndex, aNewIso).
            a_strip.insert(a_u_strip_index as usize, a_new_iso);

            let mut an_iso = a_strip[(a_u_strip_index + 1) as usize].clone();
            an_iso.change_domain_full(cutting_value, an_iso.u1(), an_iso.v0(), an_iso.v1());
            a_strip[(a_u_strip_index + 1) as usize] = an_iso;
        }

        //  Insert the new nodes (U*,Vj)
        let mut a_next: Node;
        let mut a_prev = self.my_node_constraints.first().expect("empty nodes").clone();
        let mut a_node_index: i32 = 1;
        while a_node_index < self.my_node_constraints.len() as i32 {
            a_next = self.my_node_constraints[(a_node_index) as usize].clone();
            if a_prev.coord().x < cutting_value
                && a_next.coord().x > cutting_value
                && a_prev.coord().y == a_next.coord().y
            {
                let a_new_uv = DVec2::new(cutting_value, a_prev.coord().y);
                let a_new_node = Node::new_with_coord(a_new_uv, a_prev.u_order(), a_prev.v_order());
                // OCCT: myNodeConstraints.InsertAfter(aNodeIndex, aNewNode).
                self.my_node_constraints.insert(a_node_index as usize, a_new_node);
            }
            a_prev = a_next;
            a_node_index += 1;
        }
    }

    /// OCCT UpdateInV(CuttingValue) (AdvApp2Var_Framework.cxx L292-381).
    pub fn update_in_v(&mut self, cutting_value: f64) {
        let mut a_v_strip_index: i32 = 1;
        while self.my_v_constraints[(a_v_strip_index - 1) as usize]
            .first()
            .expect("empty V strip")
            .v0()
            > cutting_value
            || self.my_v_constraints[(a_v_strip_index - 1) as usize]
                .first()
                .expect("empty V strip")
                .v1()
                < cutting_value
        {
            a_v_strip_index += 1;
        }

        {
            // The OCCT loop mutates the shared handles of the strip in
            // place; rcad works on a value copy of the strip and writes it
            // back (same final state).
            let mut a_strip = self.my_v_constraints[(a_v_strip_index - 1) as usize].clone();
            let a_v_first = a_strip.first().expect("empty V strip").v0();
            let a_v_last = a_strip.first().expect("empty V strip").v1();

            // Modify the U isos of the V strip at aVStripIndex.
            for k in 0..a_strip.len() {
                a_strip[k].change_domain(a_v_first, cutting_value);
                a_strip[k].reset_approx();
            }

            // Insert a new V strip after aVStripIndex.
            let mut a_new_strip: Vec<Iso> = Vec::new();
            for an_iso in a_strip.iter() {
                let mut a_new_iso = Iso::new_full(
                    an_iso.type_(),
                    an_iso.constante(),
                    an_iso.u0(),
                    an_iso.u1(),
                    cutting_value,
                    a_v_last,
                    0,
                    an_iso.u_order(),
                    an_iso.v_order(),
                );
                a_new_iso.reset_approx();
                a_new_strip.push(a_new_iso);
            }
            // OCCT: myVConstraints.InsertAfter(aVStripIndex, aNewStrip);
            // and the modified strip is stored back.
            self.my_v_constraints.insert(a_v_strip_index as usize, a_new_strip);
            self.my_v_constraints[(a_v_strip_index - 1) as usize] = a_strip;
        }

        // Insert a new Iso V=V* in each U strip after aVStripIndex
        // and restrict the domains of the adjacent Isos
        for k in 0..self.my_u_constraints.len() {
            let a_strip = &mut self.my_u_constraints[k];
            let mut an_iso = a_strip[(a_v_strip_index - 1) as usize].clone();
            an_iso.change_domain_full(an_iso.u0(), an_iso.u1(), an_iso.v0(), cutting_value);

            let mut a_new_iso = Iso::new_full(
                an_iso.type_(),
                cutting_value,
                an_iso.u0(),
                an_iso.u1(),
                an_iso.v0(),
                cutting_value,
                0,
                an_iso.u_order(),
                an_iso.v_order(),
            );
            a_new_iso.reset_approx();
            // OCCT: aStrip.InsertAfter(aVStripIndex, aNewIso).
            a_strip.insert(a_v_strip_index as usize, a_new_iso);

            let mut an_iso = a_strip[(a_v_strip_index + 1) as usize].clone();
            an_iso.change_domain_full(an_iso.u0(), an_iso.u1(), cutting_value, an_iso.v1());
            a_strip[(a_v_strip_index + 1) as usize] = an_iso;
        }

        //  Insert the new nodes (Ui,V*)
        let mut a_node_start_index: i32 = 1;
        while a_node_start_index <= self.my_node_constraints.len() as i32
            && self.my_node_constraints[(a_node_start_index - 1) as usize].coord().y
                < cutting_value
        {
            a_node_start_index += self.my_u_constraints.len() as i32 + 1;
        }
        for an_iso_index in 1..=(self.my_u_constraints.len() as i32 + 1) {
            let a_j_node = self.my_node_constraints[(an_iso_index - 1) as usize].clone();
            let new_uv = DVec2::new(a_j_node.coord().x, cutting_value);
            let a_new_node = Node::new_with_coord(new_uv, a_j_node.u_order(), a_j_node.v_order());
            // OCCT: myNodeConstraints.InsertAfter(
            //   aNodeStartIndex + anIsoIndex - 2, aNewNode).
            self.my_node_constraints
                .insert((a_node_start_index + an_iso_index - 2) as usize, a_new_node);
        }
    }

    /// OCCT UEquation(IndexIso, IndexStrip) (AdvApp2Var_Framework.cxx
    /// L385-390) - the const handle& return becomes a clone of the Polynom
    /// array (None when the iso has no result).
    pub fn u_equation(&self, index_iso: i32, index_strip: i32) -> Option<Vec<f64>> {
        self.my_v_constraints[(index_strip - 1) as usize][(index_iso - 1) as usize].polynom()
    }

    /// OCCT VEquation(IndexIso, IndexStrip) (AdvApp2Var_Framework.cxx
    /// L394-399).
    pub fn v_equation(&self, index_iso: i32, index_strip: i32) -> Option<Vec<f64>> {
        self.my_u_constraints[(index_strip - 1) as usize][(index_iso - 1) as usize].polynom()
    }
}
