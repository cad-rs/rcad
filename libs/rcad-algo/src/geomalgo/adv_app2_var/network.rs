//! OCCT AdvApp2Var_Network (AdvApp2Var_Network.hxx + AdvApp2Var_Network.cxx).
//!
//! Encoding notes (architecture):
//! - `NCollection_Sequence<occ::handle<AdvApp2Var_Patch>>` -> `Vec<Patch>`
//!   (each patch appears exactly once in the sequence, so the shared-handle
//!   semantics collapse to slot access; Value(i) -> vec[i - 1], Append ->
//!   push, InsertBefore(i, x) -> vec.insert(i - 1, x), InsertAfter(i, x) ->
//!   vec.insert(i, x)).
//! - `NCollection_Sequence<double>` -> `Vec<f64>`.

use super::patch::Patch;

/// OCCT AdvApp2Var_Network (AdvApp2Var_Network.hxx L29-78).
#[derive(Debug, Clone)]
pub struct Network {
    /// hxx L75: myNet.
    my_net: Vec<Patch>,
    /// hxx L76: myUParameters.
    my_u_parameters: Vec<f64>,
    /// hxx L77: myVParameters.
    my_v_parameters: Vec<f64>,
}

impl Network {
    /// OCCT AdvApp2Var_Network() (AdvApp2Var_Network.cxx L27).
    pub fn new() -> Self {
        Network {
            my_net: Vec::new(),
            my_u_parameters: Vec::new(),
            my_v_parameters: Vec::new(),
        }
    }

    /// OCCT AdvApp2Var_Network(Net, TheU, TheV) (AdvApp2Var_Network.cxx
    /// L31-39).
    pub fn new_from(net: Vec<Patch>, the_u: Vec<f64>, the_v: Vec<f64>) -> Self {
        Network {
            my_net: net,
            my_u_parameters: the_u,
            my_v_parameters: the_v,
        }
    }

    /// OCCT FirstNotApprox(Index) (AdvApp2Var_Network.cxx L43-58) - search
    /// the Index of the first Patch not approximated, if all Patches are
    /// approximated false is returned.
    pub fn first_not_approx(&self, the_index: &mut i32) -> bool {
        let mut an_index: i32 = 1;
        for a_patch in self.my_net.iter() {
            if !a_patch.is_approximated() {
                *the_index = an_index;
                return true;
            }
            an_index += 1;
        }
        false
    }

    /// OCCT ChangePatch(Index) (hxx L44) - mutable slot access (the
    /// `myResult(FirstNA)` operator() target).
    pub fn change_patch(&mut self, index: i32) -> &mut Patch {
        &mut self.my_net[(index - 1) as usize]
    }

    /// OCCT UpdateInU(CuttingValue) (AdvApp2Var_Network.cxx L62-91).
    pub fn update_in_u(&mut self, cutting_value: f64) {
        // Insert the new cutting parameter.
        let mut i: i32 = 1;
        while self.my_u_parameters[(i - 1) as usize] < cutting_value {
            i += 1;
        }
        // OCCT: myUParameters.InsertBefore(i, CuttingValue).
        self.my_u_parameters.insert((i - 1) as usize, cutting_value);

        let a_nb_patch_in_u = self.my_u_parameters.len() as i32 - 1;
        let mut j: i32 = 1;
        while j < self.my_v_parameters.len() as i32 {
            // Modify the patch impacted by the cut.
            let a_patch_index = a_nb_patch_in_u * (j - 1) + i - 1;
            // OCCT evaluation order: the arguments aPat->U0()/V0()/V1() are
            // read before ChangeDomain mutates the patch.
            let (a_pat_u0, a_pat_v0, a_pat_v1, a_pat_uorder, a_pat_vorder) = {
                let a_pat = &self.my_net[(a_patch_index - 1) as usize];
                (
                    a_pat.u0(),
                    a_pat.v0(),
                    a_pat.v1(),
                    a_pat.u_order(),
                    a_pat.v_order(),
                )
            };
            {
                let a_pat = &mut self.my_net[(a_patch_index - 1) as usize];
                a_pat.change_domain(a_pat_u0, cutting_value, a_pat_v0, a_pat_v1);
                a_pat.reset_approx();
            }

            // Insert the right-side patch.
            let a_new_pat = Patch::new_with_domain(
                cutting_value,
                self.my_u_parameters[i as usize],
                self.my_v_parameters[(j - 1) as usize],
                self.my_v_parameters[j as usize],
                a_pat_uorder,
                a_pat_vorder,
            );
            let mut a_new_pat = a_new_pat;
            a_new_pat.reset_approx();
            // OCCT: myNet.InsertAfter(aPatchIndex, aNewPat).
            self.my_net.insert(a_patch_index as usize, a_new_pat);
            j += 1;
        }
    }

    /// OCCT UpdateInV(CuttingValue) (AdvApp2Var_Network.cxx L95-130).
    pub fn update_in_v(&mut self, cutting_value: f64) {
        // Insert the new cutting parameter.
        let mut j: i32 = 1;
        while self.my_v_parameters[(j - 1) as usize] < cutting_value {
            j += 1;
        }
        // OCCT: myVParameters.InsertBefore(j, CuttingValue).
        self.my_v_parameters.insert((j - 1) as usize, cutting_value);

        let a_nb_patch_in_u = self.my_u_parameters.len() as i32 - 1;

        // Modify the patches affected by the cut.
        let mut i: i32 = 1;
        while i <= a_nb_patch_in_u {
            let a_patch_index = a_nb_patch_in_u * (j - 2) + i;
            let (a_patch_u0, a_patch_u1, a_patch_v0) = {
                let a_patch = &self.my_net[(a_patch_index - 1) as usize];
                (a_patch.u0(), a_patch.u1(), a_patch.v0())
            };
            let a_patch = &mut self.my_net[(a_patch_index - 1) as usize];
            a_patch.change_domain(a_patch_u0, a_patch_u1, a_patch_v0, cutting_value);
            a_patch.reset_approx();
            i += 1;
        }

        // Insert the top patches.
        let mut i: i32 = 1;
        while i <= a_nb_patch_in_u {
            let a_patch_index = a_nb_patch_in_u * (j - 1) + i - 1;
            let (a_source_uorder, a_source_vorder) = {
                let a_source_patch = &self.my_net[(a_nb_patch_in_u * (j - 2) + i - 1) as usize];
                (a_source_patch.u_order(), a_source_patch.v_order())
            };
            let a_new_pat = Patch::new_with_domain(
                self.my_u_parameters[(i - 1) as usize],
                self.my_u_parameters[i as usize],
                cutting_value,
                self.my_v_parameters[j as usize],
                a_source_uorder,
                a_source_vorder,
            );
            let mut a_new_pat = a_new_pat;
            a_new_pat.reset_approx();
            // OCCT: myNet.InsertAfter(aPatchIndex, aNewPat).
            self.my_net.insert(a_patch_index as usize, a_new_pat);
            i += 1;
        }
    }

    /// OCCT SameDegree(iu, iv, ncfu, ncfv) (AdvApp2Var_Network.cxx
    /// L134-156) - a non-const member (it mutates the patches through the
    /// shared handles).
    pub fn same_degree(&mut self, iu: i32, iv: i32, ncfu: &mut i32, ncfv: &mut i32) {
        //  Compute the max coefficients, initialized according to the
        //  continuity order
        *ncfu = 2 * iu + 2;
        *ncfv = 2 * iv + 2;
        for a_pat in self.my_net.iter() {
            *ncfu = (*ncfu).max(a_pat.nb_coeff_in_u());
            *ncfv = (*ncfv).max(a_pat.nb_coeff_in_v());
        }

        //  Increase the number of coefficients
        for a_pat in self.my_net.iter_mut() {
            a_pat.change_nb_coeff(*ncfu, *ncfv);
        }
    }

    /// OCCT NbPatch() (AdvApp2Var_Network.cxx L160-163).
    pub fn nb_patch(&self) -> i32 {
        self.my_net.len() as i32
    }

    /// OCCT NbPatchInU() (AdvApp2Var_Network.cxx L167-170).
    pub fn nb_patch_in_u(&self) -> i32 {
        self.my_u_parameters.len() as i32 - 1
    }

    /// OCCT NbPatchInV() (AdvApp2Var_Network.cxx L174-177).
    pub fn nb_patch_in_v(&self) -> i32 {
        self.my_v_parameters.len() as i32 - 1
    }

    /// OCCT UParameter(Index) (AdvApp2Var_Network.cxx L181-184).
    pub fn u_parameter(&self, index: i32) -> f64 {
        self.my_u_parameters[(index - 1) as usize]
    }

    /// OCCT VParameter(Index) (AdvApp2Var_Network.cxx L188-191).
    pub fn v_parameter(&self, index: i32) -> f64 {
        self.my_v_parameters[(index - 1) as usize]
    }

    /// OCCT Patch(UIndex, VIndex) (hxx L64-67) - const access.
    pub fn patch(&self, u_index: i32, v_index: i32) -> &Patch {
        &self.my_net[((v_index - 1) * (self.my_u_parameters.len() as i32 - 1) + u_index - 1)
            as usize]
    }
}
