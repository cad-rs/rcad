//! 1:1 translation of OCCT `ShapeFix_Shell` — the class
//! (`TKShHealing/ShapeFix/ShapeFix_Shell.hxx` L17-133 +
//! `ShapeFix_Shell.lxx` L17-36 + `ShapeFix_Shell.cxx` L66-177 and
//! L1425-1727).  The file statics live in [`super::shell_statics`].
//!
//! Function-count equation (OCCT Shell.cxx class members + lxx = rcad
//! shell.rs): `ShapeFix_Shell()` (cxx L68-76) / `ShapeFix_Shell(shape)`
//! (cxx L80-88) / `Init` (cxx L92-97) / `Perform` (cxx L101-174) /
//! `FixFaceOrientation` (cxx L1425-1653) / `Status` (cxx L1657-1660) /
//! `Shell` (cxx L1664-1667) / `Shape` (cxx L1671-1674) / `ErrorFaces` (cxx
//! L1678-1681) / `SetMsgRegistrator` (cxx L1685-1689) / `SetPrecision` (cxx
//! L1693-1697) / `SetMinTolerance` (cxx L1701-1705) / `SetMaxTolerance` (cxx
//! L1709-1713) / `NbShells` (cxx L1717-1720) / `SetNonManifoldFlag` (cxx
//! L1724-1727) — 15 members = 15 rcad methods; the lxx inlines
//! `FixFaceTool` (L19-22) / `FixFaceMode` (L26-29) / `FixOrientationMode`
//! (L33-36) — 3 = 3.
//!
//! Architecture bridges:
//! 1. `ShapeFix_Root` inheritance — the rcad composition field `base` (the
//!    ShapeFix_Wire precedent).
//! 2. `BRep` pool argument — the rcad equivalents take `brep: &mut BRep`.
//! 3. OCCT overload disambiguation: `Shape()` -> `shape_result` (the member
//!    read; `Shape::shape_type` is the kernel method), `Shell()` ->
//!    `shell_result`.
//! 4. `TopoDS_Shell::Closed()` — the TShape CLOSED flag read;
//!    `TopoDS_Shell::Closed(bool)` — [`super::shell_statics::
//!    set_shell_closed`].

use std::collections::HashMap;

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType, TShape, tshape_flags};

use crate::shhealing::shape_analysis::shell::ShapeAnalysisShell;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes, topexp_explorer};
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::face_a::ShapeFixFace;
use crate::shhealing::shape_fix::root::{MsgRegistratorHandle, ShapeFixRoot};
use crate::shhealing::shape_fix::shape_fix::MessageProgressRange;

use super::shell_statics::{
    add_multi_conexity_faces, create_closed_shell, create_non_manifold_shells, get_shells,
    set_shell_closed, ShapeSet,
};

/// OCCT ShapeFix_Shell (hxx L37-129): fixing orientation of faces in shell.
pub struct ShapeFixShell {
    /// The OCCT ShapeFix_Root base subobject (bridge #1).
    pub base: ShapeFixRoot,
    /// OCCT myShell (hxx L121).
    pub(crate) my_shell: Shape,
    /// OCCT myErrFaces (hxx L122).
    pub(crate) my_err_faces: Shape,
    /// OCCT myStatus (hxx L123).
    pub(crate) my_status: i32,
    /// OCCT myFixFace (hxx L124).
    pub(crate) my_fix_face: ShapeFixFace,
    /// OCCT myFixFaceMode (hxx L125).
    pub(crate) my_fix_face_mode: i32,
    /// OCCT myFixOrientationMode (hxx L126).
    pub(crate) my_fix_orientation_mode: i32,
    /// OCCT myNbShells (hxx L127).
    pub(crate) my_nb_shells: i32,
    /// OCCT myNonManifold (hxx L128).
    pub(crate) my_non_manifold: bool,
}

impl Default for ShapeFixShell {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixShell {
    // OCCT ShapeFix_Shell.cxx L68-76 — ShapeFix_Shell().
    /// OCCT ShapeFix_Shell::ShapeFix_Shell() (cxx L68-76): empty constructor.
    pub fn new() -> Self {
        ShapeFixShell {
            base: ShapeFixRoot::new(),
            my_shell: Shape::null(),
            my_err_faces: Shape::null(),
            my_status: encode_status(ShapeExtendStatus::Ok),
            my_fix_face: ShapeFixFace::new(),
            my_fix_face_mode: -1,
            my_fix_orientation_mode: -1,
            my_nb_shells: 0,
            my_non_manifold: false,
        }
    }

    // OCCT ShapeFix_Shell.cxx L80-88 — ShapeFix_Shell(shape).
    /// OCCT ShapeFix_Shell::ShapeFix_Shell(shape) (cxx L80-88): initializes
    /// by a shell.
    pub fn with_shell(brep: &mut BRep, shape: &Shape) -> Self {
        let mut this = ShapeFixShell {
            base: ShapeFixRoot::new(),
            my_shell: Shape::null(),
            my_err_faces: Shape::null(),
            my_status: encode_status(ShapeExtendStatus::Ok),
            my_fix_face: ShapeFixFace::new(),
            my_fix_face_mode: -1,
            my_fix_orientation_mode: -1,
            my_nb_shells: 0,
            my_non_manifold: false,
        };
        this.init(brep, shape);
        this.my_non_manifold = false;
        this
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shell.lxx L19-36 — the inline accessors.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shell::FixFaceTool (lxx L19-22): returns the tool for
    /// fixing faces.
    pub fn fix_face_tool(&mut self) -> &mut ShapeFixFace {
        &mut self.my_fix_face
    }

    /// OCCT ShapeFix_Shell::FixFaceMode (lxx L26-29).
    pub fn fix_face_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_face_mode
    }

    /// OCCT ShapeFix_Shell::FixOrientationMode (lxx L33-36).
    pub fn fix_orientation_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_orientation_mode
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shell.cxx L92-97 — Init.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shell::Init (cxx L92-97): initializes by a shell.
    pub fn init(&mut self, _brep: &mut BRep, shell: &Shape) {
        self.base.my_shape = shell.clone();
        self.my_shell = shell.clone();
        self.my_nb_shells = 0;
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shell.cxx L101-174 — Perform.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shell::Perform (cxx L101-174): iterates on subshapes
    /// and performs fixes (for each face calls ShapeFix_Face::Perform and
    /// then calls FixFaceOrientation).
    pub fn perform(&mut self, brep: &mut BRep, the_progress: MessageProgressRange) -> bool {
        let mut status = false;
        // OCCT L104-107.
        if self.base.my_context.is_none() {
            self.base.set_context(ShapeBuildReShape::new());
        }
        self.my_fix_face.base.my_context = self.base.my_context.clone();

        // OCCT L110-137.
        if ShapeFixRoot::need_fix(self.my_fix_face_mode, true) {
            let my_shell = self.my_shell.clone();
            let s = self.context_apply(brep, &my_shell);

            // OCCT L114-118: the progress scope (the no-abort bridge).
            let a_ps = MessageProgressScope::new();

            for value in iter_subshapes(brep, &s, true, true) {
                if !a_ps.more() {
                    break;
                }
                let sh = value;
                let tmp_face = sh;
                self.my_fix_face.init_face(brep, &tmp_face);
                if self.my_fix_face.perform(brep, the_progress) {
                    status = true;
                    self.my_status |= encode_status(ShapeExtendStatus::Done1);
                }
                a_ps.next();
            }

            // Halt algorithm in case of user's abort (the no-abort bridge:
            // the loop always runs to completion).
        }

        // OCCT L139.
        let my_shell = self.my_shell.clone();
        let newsh = self.context_apply(brep, &my_shell);
        // OCCT L140-143.
        if ShapeFixRoot::need_fix(self.my_fix_orientation_mode, true) {
            self.fix_face_orientation(brep, &newsh, true, self.my_non_manifold);
        }

        // OCCT L145-163.
        let a_newsh = self.context_apply(brep, &newsh);
        let mut a_sas = ShapeAnalysisShell::new();
        for a_shell_exp in topexp_explorer(brep, &a_newsh, ShapeType::Shell) {
            let a_cur_shell = a_shell_exp;
            if shell_closed(brep, &a_cur_shell) {
                a_sas.load_shells(brep, &a_cur_shell);
                a_sas.check_oriented_shells(brep, &a_cur_shell, true, false);
                if a_sas.has_free_edges() {
                    set_shell_closed(brep, &a_cur_shell, false);
                    // Shell has incorrect flag isClosed
                    self.base
                        .send_warning_own(&MessageMsg::from_key("FixAdvShell.FixClosedFlag.MSG0"));
                }
                a_sas.clear();
            }
        }

        // OCCT L165-173.
        if status {
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
        }
        if self.status(ShapeExtendStatus::Done2) {
            status = true;
        }
        status
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shell.cxx L1425-1653 — FixFaceOrientation.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shell::FixFaceOrientation (cxx L1425-1653): fixes the
    /// orientation of faces in the shell; changes the orientation of a face
    /// opposite to the neighbouring faces; the faces that cannot be oriented
    /// are stored in the Error compound.
    pub fn fix_face_orientation(
        &mut self,
        brep: &mut BRep,
        shell: &Shape,
        is_account_multi_conex: bool,
        non_manifold: bool,
    ) -> bool {
        // OCCT L1429: myStatus = EncodeStatus(OK) (commented out in OCCT).
        let mut done = false;
        let mut a_seq_shells: Vec<Shape> = Vec::new();
        let mut a_err_faces: Vec<Shape> = Vec::new(); // Compound of faces like to Mebius leaf.
        let mut lface: Vec<Shape> = Vec::new();
        let mut a_map_face_shells: HashMap<(u64, u32), Shape> = HashMap::new();
        self.my_shell = shell.clone();
        self.base.my_shape = shell.clone();
        let mut a_num_mult_shell = 0i32;
        let mut nb_f = 0i32;
        let mut a_map_added = ShapeSet::new();
        for value in iter_subshapes(brep, shell, true, true) {
            nb_f += 1;
            if a_map_added.add(&value) {
                lface.push(value);
            }
        }
        if (lface.len() as i32) < nb_f {
            done = true;
        }

        // OCCT L1452-1454.
        let mut a_map_edge_faces = super::shell_statics::IndexedShapeListMap::new();
        map_shapes_and_ancestors(brep, &self.my_shell, ShapeType::Edge, ShapeType::Face, &mut a_map_edge_faces);
        // OCCT L1455-1473.
        let mut a_map_multi_connect_edges = ShapeSet::new();
        let mut is_free_boundaries = false;
        for k in 1..=(a_map_edge_faces.extent()) {
            let a_face_count = a_map_edge_faces.find_from_index(k).len() as i32;
            if !is_free_boundaries && a_face_count == 1 {
                let e = a_map_edge_faces.find_key(k).clone();
                if !crate::shhealing::shape_fix::wire::brep_tool_degenerated(&e) {
                    is_free_boundaries = true;
                }
            }
            // Finds multishared edges
            else if is_account_multi_conex && a_face_count > 2 {
                a_map_multi_connect_edges.add(a_map_edge_faces.find_key(k));
            }
        }
        // OCCT L1474-1480.
        let my_shell_is_closed = crate::shhealing::shape_build::brep_tool::brep_tool_is_closed(brep, &self.my_shell);
        if my_shell_is_closed == is_free_boundaries {
            set_shell_closed(brep, &self.my_shell, !is_free_boundaries);
            // Shell has incorrect flag isClosed
            self.base
                .send_warning_own(&MessageMsg::from_key("FixAdvShell.FixClosedFlag.MSG0"));
        }
        // OCCT L1481-1495: gets possible shells with taking in account of
        // multiconnexity.
        let mut is_get_shells = true;
        while is_get_shells && !lface.is_empty() {
            let mut a_tmp_seq_shells: Vec<Shape> = Vec::new();
            if get_shells(
                brep,
                &mut lface,
                &a_map_multi_connect_edges,
                &mut a_tmp_seq_shells,
                &mut a_map_face_shells,
                &mut a_err_faces,
            ) {
                done = true;
            }
            is_get_shells = !a_tmp_seq_shells.is_empty();
            if is_get_shells {
                a_seq_shells.append(&mut a_tmp_seq_shells);
            }
        }
        // OCCT L1496-1499.
        if !done {
            done = a_seq_shells.len() > 1;
        }
        // OCCT L1500-1520.
        let mut a_is_done = false;
        if !lface.is_empty() && !a_seq_shells.is_empty() {
            for jj in 1..=(lface.len()) {
                let key = (lface[jj - 1].ptr_id(), lface[jj - 1].location);
                if a_map_face_shells.contains_key(&key) {
                    a_map_face_shells.remove(&key);
                }
            }

            // Addition of faces having only multiconnexity boundary to shells
            // having holes containing only the multiconnexity edges
            a_is_done = add_multi_conexity_faces(
                brep,
                &mut lface,
                &a_map_multi_connect_edges,
                &mut a_seq_shells,
                &a_map_face_shells,
                &a_map_edge_faces,
                &mut a_err_faces,
                non_manifold,
            );
        }
        // OCCT L1521.
        a_num_mult_shell = a_seq_shells.len() as i32;
        // OCCT L1522-1571.
        if !a_err_faces.is_empty() {
            // if Shell contains of Mebius faces one shell will be created from
            // each those face.
            let my_err_faces = brep.add_tcompound(Vec::new());
            self.my_err_faces = my_err_faces.clone();
            let a_comp_shells = brep.add_tcompound(Vec::new());
            for n in 1..=(a_err_faces.len()) {
                builder_add(brep, &my_err_faces, &a_err_faces[n - 1]);
            }
            if a_num_mult_shell != 0 {
                if a_num_mult_shell == 1 {
                    builder_add(brep, &a_comp_shells, &a_seq_shells[0]);
                    for n1 in 1..=(a_err_faces.len()) {
                        let a_sh = brep.add_tshell(Vec::new());
                        builder_add(brep, &a_sh, &a_err_faces[n1 - 1]);
                        builder_add(brep, &a_comp_shells, &a_sh);
                    }
                    self.base.my_shape = a_comp_shells.clone();
                } else {
                    for i in 1..=(a_seq_shells.len()) {
                        builder_add(brep, &a_comp_shells, &a_seq_shells[i - 1]);
                    }
                    for n1 in 1..=(a_err_faces.len()) {
                        let a_sh = brep.add_tshell(Vec::new());
                        builder_add(brep, &a_sh, &a_err_faces[n1 - 1]);
                        builder_add(brep, &a_comp_shells, &a_sh);
                    }
                    self.base.my_shape = a_comp_shells.clone();
                }
            }

            done = true;
            self.my_status = encode_status(ShapeExtendStatus::Fail);
            // Impossible to orient faces in shell, several shells created
            self.base
                .send_warning_own(&MessageMsg::from_key("FixAdvShell.FixOrientation.MSG20"));
            return true;
        }
        // OCCT L1572-1590.
        if a_num_mult_shell > 1 {
            let mut open_shells: Vec<Shape> = Vec::new();
            let mut i1 = 0usize;
            while i1 < a_seq_shells.len() {
                let a_shell = a_seq_shells[i1].clone();
                if !crate::shhealing::shape_build::brep_tool::brep_tool_is_closed(brep, &a_shell) {
                    open_shells.push(a_shell);
                    a_seq_shells.remove(i1);
                } else {
                    i1 += 1;
                }
            }
            if open_shells.len() > 1 {
                // Attempt of creation closed shell from open shells with
                // taking into account multiconnexity.
                create_closed_shell(brep, &mut open_shells, &a_map_multi_connect_edges);
            }
            a_seq_shells.append(&mut open_shells);
        }

        // OCCT L1592-1605.
        if !lface.is_empty() {
            for i in 1..=(lface.len()) {
                let one_shell = brep.add_tshell(Vec::new());
                builder_add(brep, &one_shell, &lface[i - 1]);
                a_seq_shells.push(one_shell);
            }
        }
        // OCCT L1606-1609.
        if non_manifold && a_seq_shells.len() > 1 {
            create_non_manifold_shells(brep, &mut a_seq_shells, &a_map_multi_connect_edges);
        }
        // OCCT L1610-1631.
        if !done {
            done = a_seq_shells.len() > 1 || a_is_done;
        }
        if a_seq_shells.len() == 1 {
            self.my_shell = a_seq_shells[0].clone();
            let my_shell = self.my_shell.clone();
            self.base.my_shape = my_shell;
            self.my_nb_shells = 1;
        } else {
            let a_comp_shells = brep.add_tcompound(Vec::new());
            for i in 1..=(a_seq_shells.len()) {
                builder_add(brep, &a_comp_shells, &a_seq_shells[i - 1]);
            }
            self.base.my_shape = a_comp_shells.clone();
            self.my_nb_shells = a_seq_shells.len() as i32;
        }
        // OCCT L1632-1652.
        if done {
            self.my_status = encode_status(ShapeExtendStatus::Done2);
            if self.base.my_context.is_some() {
                let ctx = self.base.my_context.as_mut().unwrap();
                ctx.replace(brep, shell, &self.base.my_shape);
            }
            if self.my_nb_shells == 1 {
                // Faces were incorrectly oriented in the shell, corrected
                self.base
                    .send_warning_own(&MessageMsg::from_key("FixAdvShell.FixOrientation.MSG0"));
            } else {
                // Improperly connected shell split into parts
                self.base
                    .send_warning_own(&MessageMsg::from_key("FixAdvShell.FixOrientation.MSG30"));
            }
            return true;
        } else {
            return false;
        }
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shell.cxx L1657-1727 — the accessors and the setters.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shell::Status (cxx L1657-1660): the status of the last
    /// Fix.
    pub fn status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }

    /// OCCT ShapeFix_Shell::Shell (cxx L1664-1667): the fixed shell (or the
    /// subset of oriented faces).
    pub fn shell_result(&self) -> Shape {
        self.my_shell.clone()
    }

    /// OCCT ShapeFix_Shell::Shape (cxx L1671-1674): in case of multiconnexity
    /// returns the compound of the fixed shells, else the one shell.
    pub fn shape_result(&self) -> Shape {
        self.base.my_shape.clone()
    }

    /// OCCT ShapeFix_Shell::ErrorFaces (cxx L1678-1681): the not oriented
    /// subset of faces.
    pub fn error_faces(&self) -> Shape {
        self.my_err_faces.clone()
    }

    /// OCCT ShapeFix_Shell::SetMsgRegistrator (cxx L1685-1689).
    pub fn set_msg_registrator(&mut self, msgreg: MsgRegistratorHandle) {
        self.base.set_msg_registrator(msgreg.clone());
        self.my_fix_face.set_msg_registrator(msgreg);
    }

    /// OCCT ShapeFix_Shell::SetPrecision (cxx L1693-1697).
    pub fn set_precision(&mut self, preci: f64) {
        self.base.set_precision(preci);
        self.my_fix_face.set_precision(preci);
    }

    /// OCCT ShapeFix_Shell::SetMinTolerance (cxx L1701-1705).
    pub fn set_min_tolerance(&mut self, mintol: f64) {
        self.base.set_min_tolerance(mintol);
        self.my_fix_face.base.set_min_tolerance(mintol);
    }

    /// OCCT ShapeFix_Shell::SetMaxTolerance (cxx L1709-1713).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.base.set_max_tolerance(maxtol);
        self.my_fix_face.base.set_max_tolerance(maxtol);
    }

    /// OCCT ShapeFix_Shell::NbShells (cxx L1717-1720): the number of the
    /// obtained shells.
    pub fn nb_shells(&self) -> i32 {
        self.my_nb_shells
    }

    /// OCCT ShapeFix_Shell::SetNonManifoldFlag (cxx L1724-1727): sets the
    /// NonManifold flag.
    pub fn set_non_manifold_flag(&mut self, is_non_manifold: bool) {
        self.my_non_manifold = is_non_manifold;
    }

    // -----------------------------------------------------------------------
    // Internal helpers shared with the statics submodule.
    // -----------------------------------------------------------------------

    /// OCCT `Context()->Apply(shape)` — the root-context application; a None
    /// context returns the shape unchanged (the callers guard with
    /// `Context().IsNull()` exactly as OCCT does).
    pub(crate) fn context_apply(&mut self, brep: &mut BRep, shape: &Shape) -> Shape {
        match self.base.my_context.as_mut() {
            Some(ctx) => ctx.apply(brep, shape, ShapeType::Shape),
            None => shape.clone(),
        }
    }
}

/// OCCT `Message_ProgressScope` — the no-abort bridge (the shape_fix.rs
/// precedent).
pub struct MessageProgressScope;

impl MessageProgressScope {
    pub fn new() -> Self {
        MessageProgressScope
    }
    pub fn more(&self) -> bool {
        true
    }
    pub fn next(&self) {}
}

impl Default for MessageProgressScope {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT `TopoDS_Shell::Closed()` — the TShape CLOSED flag read.
pub(crate) fn shell_closed(brep: &BRep, shell: &Shape) -> bool {
    match brep.tshapes.get(shell.index).map(|ts| ts.as_ref()) {
        Some(TShape::Shell(sd)) => sd.flags & tshape_flags::CLOSED != 0,
        _ => false,
    }
}

/// OCCT `TopExp::MapShapesAndAncestors(S, TS, TA, M)` (TopExp.cxx L336-368):
/// for each ancestor of type TA, the subshapes of type TS are keyed with the
/// ancestor appended (the duplicates of the OCCT exploration preserved).
pub(crate) fn map_shapes_and_ancestors(
    brep: &mut BRep,
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut super::shell_statics::IndexedShapeListMap,
) {
    // OCCT L339-350: visit the ancestors.
    for anc in topexp_explorer(brep, s, ta) {
        for sub in topexp_explorer(brep, &anc, ts) {
            let key = (sub.ptr_id(), sub.location);
            match m.seek(&sub) {
                Some(_) => {
                    // append to the existing entry
                    let idx = m.index_of(&key);
                    m.append(idx, anc.clone());
                }
                None => {
                    m.add(sub, vec![anc.clone()]);
                }
            }
        }
    }
    let _ = Orientation::Forward;
}
