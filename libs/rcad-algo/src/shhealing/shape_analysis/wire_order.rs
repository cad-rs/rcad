//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_WireOrder` (`ShapeAnalysis_WireOrder.hxx` L17-181 +
//! `ShapeAnalysis_WireOrder.cxx` L1-868).
//!
//! This class is intended to control and, if possible, redefine the order of
//! a list of edges which define a wire.  Edges are not given directly, but
//! as their bounds (start, end).  It can work either in 2D, in 3D, or in
//! miscible mode.
//!
//! Architecture bridges:
//! 1. `gp_XYZ` / `gp_XY` -> `DVec3` / `DVec2`; `NCollection_HArray1<int>` /
//!    `NCollection_HSequence<T>` -> `Vec<T>` (null handle = `None` for the
//!    result arrays, per the IsDone/IsNull checks).
//! 2. OCCT arrays and sequences are 1-based; the rcad Vec indices are the
//!    OCCT indices minus one (kept visible at each access).
//! 3. Output parameters use `&mut` (the AGENTS 1:1 rule); OCCT overloads of
//!    `Add` are disambiguated with `_3d` / `_2d` / `_both` suffixes.

use glam::{DVec2, DVec3};
use rcad_kernel::precision::{square_p_confusion, CONFUSION, SQUARE_CONFUSION};

// OCCT Standard_Real.hxx L179-186: RealLast() - the biggest representable real.
const REAL_LAST: f64 = f64::MAX;
// OCCT Standard_Real.hxx L132-140: RealSmall() - the smallest positive real.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;

/// OCCT ModeType (hxx L161-166): the mode in which the algorithm works.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ModeType {
    Mode2D,
    Mode3D,
    ModeBoth,
}

/// OCCT ShapeAnalysis_WireOrder (hxx L51-178).
pub struct ShapeAnalysisWireOrder {
    // occ::handle<NCollection_HArray1<int>> myOrd (hxx L168); null = None.
    my_ord: Option<Vec<i32>>,
    // occ::handle<NCollection_HArray1<int>> myChains (hxx L169).
    my_chains: Option<Vec<i32>>,
    // occ::handle<NCollection_HArray1<int>> myCouples (hxx L170).
    my_couples: Option<Vec<i32>>,
    // occ::handle<NCollection_HSequence<gp_XYZ>> myXYZ (hxx L171).
    my_xyz: Vec<DVec3>,
    // occ::handle<NCollection_HSequence<gp_XY>> myXY (hxx L172).
    my_xy: Vec<DVec2>,
    my_tol: f64,
    my_gap: f64,
    my_stat: i32,
    my_keep_loops: bool,
    my_mode: ModeType,
}

/// OCCT static DISTABS (cxx L170-173).
fn distabs(v1: DVec3, v2: DVec3) -> f64 {
    (v1.x - v2.x).abs() + (v1.y - v2.y).abs() + (v1.z - v2.z).abs()
}

/// OCCT gp_XYZ::IsEqual(theOther, theTolerance) (gp_XYZ.hxx L164-169): the
/// per-component absolute difference below the tolerance.
fn gp_xyz_is_equal(a: DVec3, b: DVec3, the_tolerance: f64) -> bool {
    (a.x - b.x).abs() < the_tolerance
        && (a.y - b.y).abs() < the_tolerance
        && (a.z - b.z).abs() < the_tolerance
}

impl Default for ShapeAnalysisWireOrder {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisWireOrder {
    /// OCCT ShapeAnalysis_WireOrder() (cxx L32-40): empty constructor.
    pub fn new() -> Self {
        let mut w = ShapeAnalysisWireOrder {
            my_ord: None,
            my_chains: None,
            my_couples: None,
            my_xyz: Vec::new(),
            my_xy: Vec::new(),
            my_gap: 0.0,
            my_stat: 0,
            my_keep_loops: false,
            my_mode: ModeType::Mode3D,
            my_tol: 0.0,
        };
        w.my_tol = CONFUSION;
        w.clear();
        w
    }

    /// OCCT ShapeAnalysis_WireOrder(theMode3D, theTolerance, theModeBoth)
    /// (cxx L44-68): creates a WireOrder.  Flag theMode3D defines 3D or 2d
    /// mode; theModeBoth defines miscible mode (theMode3D ignored).
    /// Warning: theTolerance is not used in algorithm.
    pub fn new_with_mode(the_mode3d: bool, the_tolerance: f64, the_mode_both: bool) -> Self {
        let mut w = ShapeAnalysisWireOrder {
            my_ord: None,
            my_chains: None,
            my_couples: None,
            my_xyz: Vec::new(),
            my_xy: Vec::new(),
            my_tol: the_tolerance,
            my_gap: 0.0,
            my_stat: 0,
            my_keep_loops: false,
            my_mode: ModeType::Mode3D,
        };
        if the_mode_both {
            w.my_mode = ModeType::ModeBoth;
        } else {
            if the_mode3d {
                w.my_mode = ModeType::Mode3D;
            } else {
                w.my_mode = ModeType::Mode2D;
            }
        }
        w.clear();
        w
    }

    /// OCCT SetMode (cxx L72-102): sets new values; clears the edge list if
    /// the mode changes; clears the connexion list.
    pub fn set_mode(&mut self, the_mode3d: bool, the_tolerance: f64, the_mode_both: bool) {
        let a_new_mode;
        if the_mode_both {
            a_new_mode = ModeType::ModeBoth;
        } else {
            if the_mode3d {
                a_new_mode = ModeType::Mode3D;
            } else {
                a_new_mode = ModeType::Mode2D;
            }
        }
        if self.my_mode != a_new_mode {
            self.clear();
        }
        self.my_mode = a_new_mode;
        self.my_ord = None; // OCCT L98: myOrd.Nullify()
        self.my_stat = 0;
        self.my_gap = 0.0;
        self.my_tol = if the_tolerance > 0.0 {
            the_tolerance
        } else {
            1.0e-08
        };
    }

    /// OCCT Tolerance (cxx L106-109): returns the working tolerance.
    pub fn tolerance(&self) -> f64 {
        self.my_tol
    }

    /// OCCT Clear (cxx L113-119): clears the list of edges, but not mode
    /// and tol.
    pub fn clear(&mut self) {
        self.my_xyz = Vec::new(); // OCCT L115: myXYZ = new HSequence<gp_XYZ>()
        self.my_xy = Vec::new(); // OCCT L116: myXY = new HSequence<gp_XY>()
        self.my_stat = 0;
        self.my_gap = 0.0;
    }

    /// OCCT Add(theStart3d, theEnd3d) (cxx L123-131): adds a couple of
    /// points 3D (start, end).
    pub fn add_3d(&mut self, the_start3d: DVec3, the_end3d: DVec3) {
        if self.my_mode == ModeType::Mode3D {
            self.my_xyz.push(the_start3d);
            self.my_xyz.push(the_end3d);
        }
    }

    /// OCCT Add(theStart2d, theEnd2d) (cxx L134-144): adds a couple of
    /// points 2D (start, end).
    pub fn add_2d(&mut self, the_start2d: DVec2, the_end2d: DVec2) {
        if self.my_mode == ModeType::Mode2D {
            // OCCT L139: val.SetCoord(theStart2d.X(), theStart2d.Y(), 0.0).
            self.my_xyz.push(DVec3::new(the_start2d.x, the_start2d.y, 0.0));
            // OCCT L141: val.SetCoord(theEnd2d.X(), theEnd2d.Y(), 0.0).
            self.my_xyz.push(DVec3::new(the_end2d.x, the_end2d.y, 0.0));
        }
    }

    /// OCCT Add(theStart3d, theEnd3d, theStart2d, theEnd2d) (cxx L148-161):
    /// adds a couple of points 3D and 2D (start, end).
    pub fn add_both(
        &mut self,
        the_start3d: DVec3,
        the_end3d: DVec3,
        the_start2d: DVec2,
        the_end2d: DVec2,
    ) {
        if self.my_mode == ModeType::ModeBoth {
            self.my_xyz.push(the_start3d);
            self.my_xyz.push(the_end3d);

            self.my_xy.push(the_start2d);
            self.my_xy.push(the_end2d);
        }
    }

    /// OCCT NbEdges (cxx L165-168): the count of added couples of points
    /// (one per edges).
    pub fn nb_edges(&self) -> i32 {
        (self.my_xyz.len() / 2) as i32
    }

    /// OCCT KeepLoopsMode (cxx L194-197): if this mode is True method
    /// perform does not sort edges of different loops; the resulting order
    /// is first loop, second one etc...
    pub fn keep_loops_mode(&mut self) -> &mut bool {
        &mut self.my_keep_loops
    }

    /// OCCT Perform (cxx L205-690): make wire order analysis and propose the
    /// better order of the edges taking into account the gaps between edges.
    /// Warning: parameter `closed` not used.
    pub fn perform(&mut self, _closed: bool) {
        self.my_stat = 0;
        let a_nb_edges = self.nb_edges();
        // no edges loaded, nothing to do -- return with status OK
        if a_nb_edges == 0 {
            return;
        }
        // OCCT L214-215: myOrd = new HArray1<int>(1, aNbEdges); Init(0).
        let mut my_ord = vec![0i32; a_nb_edges as usize];

        // sequence of the edge nums in the right order
        let mut an_edge_seq: Vec<i32> = Vec::new();
        let mut a_loops: Vec<Vec<i32>> = Vec::new();

        // the beginnings and ends of the edges (OCCT L222-225: the 1-based
        // arrays aBegins3D/anEnds3D/aBegins2D/anEnds2D).
        let mut a_begins3d = vec![DVec3::ZERO; a_nb_edges as usize];
        let mut an_ends3d = vec![DVec3::ZERO; a_nb_edges as usize];
        let mut a_begins2d = vec![DVec2::ZERO; a_nb_edges as usize];
        let mut an_ends2d = vec![DVec2::ZERO; a_nb_edges as usize];
        for i in 1..=a_nb_edges {
            a_begins3d[(i - 1) as usize] = self.my_xyz[(2 * i - 2) as usize];
            an_ends3d[(i - 1) as usize] = self.my_xyz[(2 * i - 1) as usize];
            if self.my_mode == ModeType::ModeBoth {
                a_begins2d[(i - 1) as usize] = self.my_xy[(2 * i - 2) as usize];
                an_ends2d[(i - 1) as usize] = self.my_xy[(2 * i - 1) as usize];
            }
        }
        // the flags that the edges was considered
        let mut is_edge_used = vec![false; a_nb_edges as usize];

        let a_tol2 = SQUARE_CONFUSION;
        let a_tol_p2 = square_p_confusion();

        // take the first edge to the constructed chain
        is_edge_used[0] = true;
        let mut a_first_pnt3d = a_begins3d[0];
        let mut a_last_pnt3d = an_ends3d[0];
        let mut a_first_pnt2d = DVec2::ZERO;
        let mut a_last_pnt2d = DVec2::ZERO;
        if self.my_mode == ModeType::ModeBoth {
            a_first_pnt2d = a_begins2d[0];
            a_last_pnt2d = an_ends2d[0];
        }
        an_edge_seq.push(1);

        // cycle until all edges are considered
        loop {
            // joint type
            // 0 - the start of the best edge to the end of constructed sequence (nothing to do)
            // 1 - the end of the best edge to the start of constructed sequence (need move the edge)
            // 2 - the end of the best edge to the end of constructed sequence (need to reverse)
            // 3 - the start of the best edge to the start of constructed sequence (need to reverse
            // and move the edge)
            let mut a_best_joint_type = 3i32;
            // the best minimum distance between constructed sequence and the best edge
            let mut a_best_min3d = REAL_LAST;
            // number of the best edge
            let mut a_best_edge_num = 0i32;
            // the best edge was found
            let mut is_found = false;
            let mut is_connected = false;
            // loop to find the best edge among all the remaining
            for i in 1..=a_nb_edges {
                if is_edge_used[(i - 1) as usize] {
                    continue;
                }

                // find minimum distance and joint type for 3D and 2D (if necessary) modes
                let mut a_cur_joint_type;
                let a_cur_min;
                // distance for four possible cases
                let a_seq_tail_edge_head = a_last_pnt3d.distance_squared(a_begins3d[(i - 1) as usize]);
                let a_seq_tail_edge_tail = a_last_pnt3d.distance_squared(an_ends3d[(i - 1) as usize]);
                let a_seq_head_edge_tail = a_first_pnt3d.distance_squared(an_ends3d[(i - 1) as usize]);
                let a_seq_head_edge_head = a_first_pnt3d.distance_squared(a_begins3d[(i - 1) as usize]);
                // the best distances for joints with head and tail of sequence
                let a_min_dist_to_tail: f64;
                let a_min_dist_to_head: f64;
                let a_tail_join_type: i32;
                let a_head_joint_type: i32;
                if a_seq_tail_edge_head <= a_seq_tail_edge_tail {
                    a_tail_join_type = 0;
                    a_min_dist_to_tail = a_seq_tail_edge_head;
                } else {
                    a_tail_join_type = 2;
                    a_min_dist_to_tail = a_seq_tail_edge_tail;
                }
                if a_seq_head_edge_tail <= a_seq_head_edge_head {
                    a_head_joint_type = 1;
                    a_min_dist_to_head = a_seq_head_edge_tail;
                } else {
                    a_head_joint_type = 3;
                    a_min_dist_to_head = a_seq_head_edge_head;
                }
                // comparing the head and the tail cases
                // if distances are close enough then we use rule for joint type: 0 < 1 < 2 < 3
                if (a_min_dist_to_tail - a_min_dist_to_head).abs() < a_tol2 {
                    if a_tail_join_type < a_head_joint_type {
                        a_cur_joint_type = a_tail_join_type;
                        a_cur_min = a_min_dist_to_tail;
                    } else {
                        a_cur_joint_type = a_head_joint_type;
                        a_cur_min = a_min_dist_to_head;
                    }
                } else {
                    if a_min_dist_to_tail <= a_min_dist_to_head {
                        a_cur_joint_type = a_tail_join_type;
                        a_cur_min = a_min_dist_to_tail;
                    } else {
                        a_cur_joint_type = a_head_joint_type;
                        a_cur_min = a_min_dist_to_head;
                    }
                }
                // update for the best values
                if self.my_mode == ModeType::ModeBoth {
                    // distances in 2D
                    let mut a_joint_mask3d = 0u32;
                    let mut a_joint_mask2d = 0u32;
                    if a_seq_tail_edge_head < a_tol2 {
                        a_joint_mask3d |= 1 << 0;
                    }
                    if a_seq_tail_edge_tail < a_tol2 {
                        a_joint_mask3d |= 1 << 2;
                    }
                    if a_seq_head_edge_tail < a_tol2 {
                        a_joint_mask3d |= 1 << 1;
                    }
                    if a_seq_head_edge_head < a_tol2 {
                        a_joint_mask3d |= 1 << 3;
                    }
                    let a_seq_tail_edge_head2d = a_last_pnt2d.distance_squared(a_begins2d[(i - 1) as usize]);
                    let a_seq_tail_edge_tail2d = a_last_pnt2d.distance_squared(an_ends2d[(i - 1) as usize]);
                    let a_seq_head_edge_tail2d = a_first_pnt2d.distance_squared(an_ends2d[(i - 1) as usize]);
                    let a_seq_head_edge_head2d = a_first_pnt2d.distance_squared(a_begins2d[(i - 1) as usize]);
                    if a_seq_tail_edge_head2d < a_tol_p2 {
                        a_joint_mask2d |= 1 << 0;
                    }
                    if a_seq_tail_edge_tail2d < a_tol_p2 {
                        a_joint_mask2d |= 1 << 2;
                    }
                    if a_seq_head_edge_tail2d < a_tol_p2 {
                        a_joint_mask2d |= 1 << 1;
                    }
                    if a_seq_head_edge_head2d < a_tol_p2 {
                        a_joint_mask2d |= 1 << 3;
                    }
                    // new approch for detecting best edge connection, for all other cases used
                    // old 3D algorithm
                    let a_full_mask = a_joint_mask3d & a_joint_mask2d;
                    if a_full_mask != 0 {
                        // find the best current joint type
                        a_cur_joint_type = 3;
                        for j in 0..4 {
                            if a_full_mask & (1 << j) != 0 {
                                a_cur_joint_type = j;
                                break;
                            }
                        }
                        if !is_connected || a_cur_joint_type < a_best_joint_type {
                            is_found = true;
                            is_connected = true;
                            match a_cur_joint_type {
                                0 => a_best_min3d = a_seq_tail_edge_head,
                                1 => a_best_min3d = a_seq_head_edge_tail,
                                2 => a_best_min3d = a_seq_tail_edge_tail,
                                3 => a_best_min3d = a_seq_head_edge_head,
                                _ => {}
                            }
                            a_best_joint_type = a_cur_joint_type;
                            a_best_edge_num = i;
                        }
                    }
                    // if there is still no connection, continue to use ald 3D algorithm
                    if is_connected {
                        continue;
                    }
                }
                // if the best distance is still not reached (aBestMin3D > aTol2) or we found a
                // better joint type
                if a_best_min3d > a_tol2 || a_cur_joint_type < a_best_joint_type {
                    // make a decision that this edge is good enough:
                    // - it gets the best distance but there is fabs(aCurMin3d - aBestMin3d) < aTol2 &&
                    // (aCurJointType < aBestJointType) ?
                    // - it gets the best joint in some cases
                    if a_cur_min < a_best_min3d
                        || ((a_cur_min == a_best_min3d || a_cur_min < a_tol2)
                            && (a_cur_joint_type < a_best_joint_type))
                    {
                        is_found = true;
                        a_best_min3d = a_cur_min;
                        a_best_joint_type = a_cur_joint_type;
                        a_best_edge_num = i;
                    }
                }
            }

            // check that we found edge for connecting
            if is_found {
                // distance between first and last point in sequence
                let a_close_dist = a_first_pnt3d.distance_squared(a_last_pnt3d);
                // if it's better to insert the edge than to close the loop, just insert the edge
                // according to joint type
                if a_best_min3d <= REAL_SMALL || a_best_min3d < a_close_dist {
                    match a_best_joint_type {
                        0 => {
                            an_edge_seq.push(a_best_edge_num);
                            a_last_pnt3d = an_ends3d[(a_best_edge_num - 1) as usize];
                        }
                        1 => {
                            an_edge_seq.insert(0, a_best_edge_num);
                            a_first_pnt3d = a_begins3d[(a_best_edge_num - 1) as usize];
                        }
                        2 => {
                            an_edge_seq.push(-a_best_edge_num);
                            a_last_pnt3d = a_begins3d[(a_best_edge_num - 1) as usize];
                        }
                        3 => {
                            an_edge_seq.insert(0, -a_best_edge_num);
                            a_first_pnt3d = an_ends3d[(a_best_edge_num - 1) as usize];
                        }
                        _ => {}
                    }
                    if self.my_mode == ModeType::ModeBoth {
                        match a_best_joint_type {
                            0 => {
                                a_last_pnt2d = an_ends2d[(a_best_edge_num - 1) as usize];
                            }
                            1 => {
                                a_first_pnt2d = a_begins2d[(a_best_edge_num - 1) as usize];
                            }
                            2 => {
                                a_last_pnt2d = a_begins2d[(a_best_edge_num - 1) as usize];
                            }
                            3 => {
                                a_first_pnt2d = an_ends2d[(a_best_edge_num - 1) as usize];
                            }
                            _ => {}
                        }
                    }
                }
                // closing loop and creating new one
                else {
                    a_loops.push(std::mem::take(&mut an_edge_seq));
                    a_first_pnt3d = a_begins3d[(a_best_edge_num - 1) as usize];
                    a_last_pnt3d = an_ends3d[(a_best_edge_num - 1) as usize];
                    if self.my_mode == ModeType::ModeBoth {
                        a_first_pnt2d = a_begins2d[(a_best_edge_num - 1) as usize];
                        a_last_pnt2d = an_ends2d[(a_best_edge_num - 1) as usize];
                    }
                    an_edge_seq.push(a_best_edge_num);
                }
                // mark the edge as used
                is_edge_used[(a_best_edge_num - 1) as usize] = true;
            } else {
                // the only condition under which we can't find an edge is when all edges are done
                break;
            }
        }
        // append the last loop
        a_loops.push(an_edge_seq);

        // handling with constructed loops
        let mut a_main_loop: Vec<i32>;
        if self.my_keep_loops {
            // keeping the loops, adding one after another.
            a_main_loop = Vec::new();
            for i in 1..=a_loops.len() {
                let a_cur_loop = &a_loops[i - 1];
                a_main_loop.extend_from_slice(a_cur_loop);
            }
        } else {
            // connecting loops
            a_main_loop = a_loops.remove(0);
            while !a_loops.is_empty() {
                // iterate over all loops to find the closest one
                let mut a_min_dist1 = REAL_LAST;
                let mut a_loop_num1 = 0usize;
                let mut a_cur_loop_it1 = 0i32;
                let mut a_direct1 = false;
                let mut a_main_loop_it1 = 0i32;
                for a_loop_it in 1..=a_loops.len() {
                    let a_cur_loop = &a_loops[a_loop_it - 1];
                    // iterate over all gaps between edges in current loop
                    let mut a_cur_loop_it2 = 0i32;
                    let mut a_main_loop_it2 = 0i32;
                    let mut a_direct2 = false;
                    let mut a_min_dist2 = REAL_LAST;
                    let a_cur_loop_length = a_cur_loop.len() as i32;
                    for a_cur_edge_it in 1..=a_cur_loop_length {
                        // get the distance between the current edge and the previous edge taking
                        // into account the edge's orientation
                        let a_prev_edge_it = if a_cur_edge_it == 1 {
                            a_cur_loop_length
                        } else {
                            a_cur_edge_it - 1
                        };
                        let a_cur_edge_idx = a_cur_loop[(a_cur_edge_it - 1) as usize];
                        let a_prev_edge_idx = a_cur_loop[(a_prev_edge_it - 1) as usize];
                        let a_cur_loop_first = if a_cur_edge_idx > 0 {
                            a_begins3d[(a_cur_edge_idx - 1) as usize]
                        } else {
                            an_ends3d[(-a_cur_edge_idx - 1) as usize]
                        };
                        let a_cur_loop_last = if a_prev_edge_idx > 0 {
                            an_ends3d[(a_prev_edge_idx - 1) as usize]
                        } else {
                            a_begins3d[(-a_prev_edge_idx - 1) as usize]
                        };
                        // iterate over all gaps between edges in main loop
                        let mut a_min_dist3 = REAL_LAST;
                        let mut a_main_loop_it3 = 0i32;
                        let mut a_direct3 = false;
                        let a_main_loop_length = a_main_loop.len() as i32;
                        let mut a_cur_edge_it2 = 1i32;
                        while a_cur_edge_it2 <= a_main_loop_length && a_min_dist3 != 0.0 {
                            // get the distance between the current edge and the next edge taking
                            // into account the edge's orientation
                            let a_next_edge_it2 = if a_cur_edge_it2 == a_main_loop_length {
                                1
                            } else {
                                a_cur_edge_it2 + 1
                            };
                            let a_cur_edge_idx2 = a_main_loop[(a_cur_edge_it2 - 1) as usize];
                            let a_next_edge_idx2 = a_main_loop[(a_next_edge_it2 - 1) as usize];
                            let a_main_loop_first = if a_cur_edge_idx2 > 0 {
                                an_ends3d[(a_cur_edge_idx2 - 1) as usize]
                            } else {
                                a_begins3d[(-a_cur_edge_idx2 - 1) as usize]
                            };
                            let a_main_loop_last = if a_next_edge_idx2 > 0 {
                                a_begins3d[(a_next_edge_idx2 - 1) as usize]
                            } else {
                                an_ends3d[(-a_next_edge_idx2 - 1) as usize]
                            };
                            // getting the sum of square distances if we try to sew the current loop
                            // with the main loop in current positions
                            let a_direct_dist = a_cur_loop_first.distance_squared(a_main_loop_first)
                                + a_cur_loop_last.distance_squared(a_main_loop_last);
                            let mut a_reverse_dist = a_cur_loop_first.distance_squared(a_main_loop_last)
                                + a_cur_loop_last.distance_squared(a_main_loop_first);
                            // take the best result
                            let a_join_dist;
                            if (a_direct_dist < a_tol2) || (a_direct_dist < 2.0 * a_reverse_dist) {
                                a_join_dist = a_direct_dist;
                                a_reverse_dist = a_direct_dist;
                            } else {
                                a_join_dist = a_reverse_dist;
                            }
                            // check if we found a better distance
                            if a_join_dist < a_min_dist3
                                && (a_min_dist3 - a_join_dist).abs() > a_tol2
                            {
                                a_min_dist3 = a_join_dist;
                                a_direct3 = a_direct_dist <= a_reverse_dist;
                                a_main_loop_it3 = a_cur_edge_it2;
                            }
                            a_cur_edge_it2 += 1;
                        }
                        // check if we found a better distance
                        if a_min_dist3 < a_min_dist2 && (a_min_dist2 - a_min_dist3).abs() > a_tol2
                        {
                            a_min_dist2 = a_min_dist3;
                            a_direct2 = a_direct3;
                            a_main_loop_it2 = a_main_loop_it3;
                            a_cur_loop_it2 = a_cur_edge_it;
                        }
                    }
                    // check if we found a better distance
                    if a_min_dist2 < a_min_dist1 && (a_min_dist1 - a_min_dist2).abs() > a_tol2 {
                        a_min_dist1 = a_min_dist2;
                        a_loop_num1 = a_loop_it;
                        a_direct1 = a_direct2;
                        a_main_loop_it1 = a_main_loop_it2;
                        a_cur_loop_it1 = a_cur_loop_it2;
                    }
                }
                // insert the found loop into main loop
                let a_loop = &a_loops[a_loop_num1 - 1];
                let a_factor = if a_direct1 { 1 } else { -1 };
                for i in 0..a_loop.len() {
                    // OCCT L627-628: anIdx = (aCurLoopIt1 + i > aLoop->Length()
                    // ? aCurLoopIt1 + i - aLoop->Length() : aCurLoopIt1 + i)
                    // - the 1-based index into the loop (OCCT Value() of 0
                    // would be out of bounds; the mirrored arithmetic keeps
                    // the OCCT form).
                    let an_idx = if a_cur_loop_it1 + (i as i32) > a_loop.len() as i32 {
                        a_cur_loop_it1 + (i as i32) - a_loop.len() as i32
                    } else {
                        a_cur_loop_it1 + (i as i32)
                    };
                    let inserted = a_loop[(an_idx - 1) as usize] * a_factor;
                    // OCCT L629: InsertAfter(aMainLoopIt1 + i, value).
                    a_main_loop.insert((a_main_loop_it1 + (i as i32)) as usize, inserted);
                }
                a_loops.remove(a_loop_num1 - 1);
            }
        }

        // checking the new order of the edges
        //  0 - order is the same
        //  1 - some edges were reordered
        // -1 - some edges were reversed
        let mut a_temp_status = 0i32;
        for i in 1..=(a_main_loop.len() as i32) {
            if i != a_main_loop[(i - 1) as usize] && a_temp_status >= 0 {
                a_temp_status = if a_main_loop[(i - 1) as usize] > 0 {
                    1
                } else {
                    -1
                };
            }
            my_ord[(i - 1) as usize] = a_main_loop[(i - 1) as usize];
        }
        if a_temp_status == 0 {
            self.my_stat = a_temp_status;
            self.my_ord = Some(my_ord);
            return;
        } else {
            // check if edges were only shifted in reverse or forward, not reordered
            let mut is_shift_reverse = true;
            let mut is_shift_forward = true;
            let mut a_first_idx;
            let mut a_second_idx;
            let a_length = a_main_loop.len() as i32;
            for i in 1..=(a_length - 1) {
                a_first_idx = a_main_loop[(i - 1) as usize];
                a_second_idx = a_main_loop[i as usize];
                if a_second_idx - a_first_idx != 1
                    && (a_first_idx != a_length || a_second_idx != 1)
                {
                    is_shift_forward = false;
                }
                if a_first_idx - a_second_idx != 1
                    && (a_second_idx != a_length || a_first_idx != 1)
                {
                    is_shift_reverse = false;
                }
            }
            a_first_idx = a_main_loop[(a_length - 1) as usize];
            a_second_idx = a_main_loop[0];
            if a_second_idx - a_first_idx != 1 && (a_first_idx != a_length || a_second_idx != 1) {
                is_shift_forward = false;
            }
            if a_first_idx - a_second_idx != 1 && (a_second_idx != a_length || a_first_idx != 1) {
                is_shift_reverse = false;
            }
            if is_shift_forward || is_shift_reverse {
                a_temp_status = 3;
            }
            self.my_stat = a_temp_status;
            self.my_ord = Some(my_ord);
            return;
        }
    }

    /// OCCT IsDone (cxx L694-697): tells if Perform has been done.
    pub fn is_done(&self) -> bool {
        self.my_ord.is_some()
    }

    /// OCCT Status (cxx L701-704): the status of the order (0 if not done):
    /// 0 : all edges are direct and in sequence
    /// 1 : all edges are direct but some are not in sequence
    /// -1 : some edges are reversed, but no gap remain
    /// 3 : edges in sequence are just shifted in forward or reverse manner
    pub fn status(&self) -> i32 {
        self.my_stat
    }

    /// OCCT Ordered (cxx L708-716): the number of original edge which
    /// correspond to the newly ordered number theIdx.  Warning: the returned
    /// value is NEGATIVE if edge should be reversed.
    pub fn ordered(&self, the_idx: i32) -> i32 {
        let my_ord = match self.my_ord.as_ref() {
            Some(o) if (o.len() as i32) >= the_idx => o,
            _ => return the_idx, // OCCT: myOrd.IsNull() || myOrd->Upper() < theIdx
        };
        let an_old_idx = my_ord[(the_idx - 1) as usize];
        if an_old_idx == 0 {
            the_idx
        } else {
            an_old_idx
        }
    }

    /// OCCT XYZ (cxx L720-724): the values of the couple theIdx, as 3D
    /// values.
    pub fn xyz(&self, the_idx: i32, the_start3d: &mut DVec3, the_end3d: &mut DVec3) {
        *the_start3d = self.my_xyz[(if the_idx > 0 {
            2 * the_idx - 1
        } else {
            -2 * the_idx
        } - 1) as usize];
        *the_end3d = self.my_xyz[(if the_idx > 0 {
            2 * the_idx
        } else {
            -2 * the_idx - 1
        } - 1) as usize];
    }

    /// OCCT XY (cxx L728-742): the values of the couple theIdx, as 2D
    /// values.
    pub fn xy(&self, the_idx: i32, the_start2d: &mut DVec2, the_end2d: &mut DVec2) {
        if self.my_mode == ModeType::ModeBoth {
            *the_start2d = self.my_xy[(if the_idx > 0 {
                2 * the_idx - 1
            } else {
                -2 * the_idx
            } - 1) as usize];
            *the_end2d = self.my_xy[(if the_idx > 0 {
                2 * the_idx
            } else {
                -2 * the_idx - 1
            } - 1) as usize];
        } else {
            let a_start3d = self.my_xyz[(if the_idx > 0 {
                2 * the_idx - 1
            } else {
                -2 * the_idx
            } - 1) as usize];
            the_start2d.x = a_start3d.x;
            the_start2d.y = a_start3d.y;
            let an_end3d = self.my_xyz[(if the_idx > 0 {
                2 * the_idx
            } else {
                -2 * the_idx - 1
            } - 1) as usize];
            the_end2d.x = an_end3d.x;
            the_end2d.y = an_end3d.y;
        }
    }

    /// OCCT Gap (cxx L746-758): the gap between a couple and its preceding;
    /// num is considered ordered; num = 0 (the OCCT default) returns the
    /// greatest gap found.
    pub fn gap(&self, num: i32) -> f64 {
        if num == 0 {
            return self.my_gap;
        }
        let n1 = self.ordered(num);
        let n0 = self.ordered(if num == 1 { self.nb_edges() } else { num - 1 });
        //  Distance entre fin (n0) et debut (n1)
        distabs(
            self.my_xyz[(if n0 > 0 {
                2 * n0
            } else {
                -2 * n0 - 1
            } - 1) as usize],
            self.my_xyz[(if n1 > 0 {
                2 * n1 - 1
            } else {
                -2 * n1
            } - 1) as usize],
        )
        ////  return (myXYZ->Value(2*n0)).Distance (myXYZ->Value(2*n1-1));
    }

    /// OCCT SetChains (cxx L762-802): determines the chains inside which
    /// successive edges have a gap less than a given value.  Queried by
    /// NbChains and Chain.
    pub fn set_chains(&mut self, gap: f64) {
        let mut n0;
        let mut n1;
        let mut n2;
        let mut nb = self.nb_edges(); // szv#4:S4163:12Mar99 o0,o1,o2 not needed
        if nb == 0 {
            return;
        }
        let mut chain: Vec<i32> = Vec::new();
        n0 = 0;
        chain.push(1); // On demarre la partie
        // OCCT L772: gp_XYZ f3d, l3d, f13d, l13d; - declared uninitialized;
        // the first loop iteration (n0 == 0) fills f13d/l13d before use, so
        // the ZERO pre-initialisation is observation-equivalent.
        let mut f3d = DVec3::ZERO;
        let mut l3d = DVec3::ZERO;
        let mut f13d = DVec3::ZERO;
        let mut l13d = DVec3::ZERO; // szv#4:S4163:12Mar99 f03d,l03d unused
        n1 = 1;
        while n1 <= nb {
            if n0 == 0 {
                // nouvelle boucle
                n0 = n1;
                // szv#4:S4163:12Mar99 optimized
                self.xyz(self.ordered(n0), &mut f13d, &mut l13d);
            }
            // szv#4:S4163:12Mar99 optimized
            n2 = if n1 == nb { n0 } else { n1 + 1 };
            self.xyz(self.ordered(n2), &mut f3d, &mut l3d);
            if !gp_xyz_is_equal(f3d, l13d, gap) {
                chain.push(n2);
                n0 = 0;
            }
            f13d = f3d;
            l13d = l3d;
            n1 += 1;
        }
        nb = chain.len() as i32;
        if nb == 0 {
            return;
        }
        // OCCT L797-801: myChains = new HArray1<int>(1, nb); the copy loop.
        self.my_chains = Some(vec![0i32; nb as usize]);
        let my_chains = self.my_chains.as_mut().unwrap();
        for n1 in 1..=nb {
            my_chains[(n1 - 1) as usize] = chain[(n1 - 1) as usize];
        }
    }

    /// OCCT NbChains (cxx L806-809): the count of computed chains.
    pub fn nb_chains(&self) -> i32 {
        match self.my_chains.as_ref() {
            Some(c) => c.len() as i32,
            None => 0,
        }
    }

    /// OCCT Chain (cxx L813-834): for the chain n0 num, starting and ending
    /// numbers of edges.  In the list of ordered edges (see Ordered for
    /// originals).
    pub fn chain(&self, num: i32, n1: &mut i32, n2: &mut i32) {
        *n1 = 0;
        *n2 = 0;
        let my_chains = match self.my_chains.as_ref() {
            Some(c) => c,
            None => return,
        };
        let nb = my_chains.len() as i32;
        if num == 0 || num > nb {
            return;
        }
        *n1 = my_chains[(num - 1) as usize];
        if num == nb {
            *n2 = self.nb_edges();
        } else {
            *n2 = my_chains[num as usize] - 1;
        }
    }

    /// OCCT SetCouples (cxx L838-843): determines the couples of edges for
    /// which end and start fit inside a given gap.  Queried by NbCouples and
    /// Couple.  Warning: function isn't implemented (the OCCT_DEBUG print is
    /// omitted).
    pub fn set_couples(&mut self, _gap: f64) {}

    /// OCCT NbCouples (cxx L847-850): the count of computed couples.
    pub fn nb_couples(&self) -> i32 {
        match self.my_couples.as_ref() {
            Some(c) => c.len() as i32,
            None => 0,
        }
    }

    /// OCCT Couple (cxx L854-868): for the couple n0 num, the two implied
    /// edges.  In the list of ordered edges.
    pub fn couple(&self, num: i32, n1: &mut i32, n2: &mut i32) {
        *n1 = 0;
        *n2 = 0;
        let my_couples = match self.my_couples.as_ref() {
            Some(c) => c,
            None => return,
        };
        let nb = my_couples.len() as i32;
        if num == 0 || num * 2 > nb {
            return;
        }
        *n1 = my_couples[(2 * num - 1 - 1) as usize];
        *n2 = my_couples[(2 * num - 1) as usize];
    }
}
