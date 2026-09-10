//! OCCT HLRAlgo_EdgeIterator (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_EdgeIterator.hxx` (L26-65) +
//! `HLRAlgo_EdgeIterator.cxx` (L26-94) + `HLRAlgo_EdgeIterator.lxx`
//! (L21-70).
//!
//! The OCCT raw status pointers (`EVis` / `EHid`) map to borrows of the
//! edge status (Rust cannot hold raw member pointers); the iterator caches
//! the running hidden interval exactly as OCCT does.

use super::edge_status::EdgeStatus;

/// OCCT HLRAlgo_EdgeIterator — iterator on the visible or hidden parts of
/// an edge status.
#[derive(Debug, Clone)]
pub struct EdgeIterator<'a> {
    my_nb_vis: i32,
    my_nb_hid: i32,
    e_vis: Option<&'a EdgeStatus>,
    e_hid: Option<&'a EdgeStatus>,
    i_vis: i32,
    i_hid: i32,
    my_hid_start: f64,
    my_hid_end: f64,
    my_hid_tol_start: f32,
    my_hid_tol_end: f32,
}

impl Default for EdgeIterator<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> EdgeIterator<'a> {
    /// OCCT HLRAlgo_EdgeIterator() — cxx L26-38.
    pub fn new() -> Self {
        EdgeIterator {
            my_nb_vis: 0,
            my_nb_hid: 0,
            e_vis: None,
            e_hid: None,
            i_vis: 0,
            i_hid: 0,
            my_hid_start: 0.0,
            my_hid_end: 0.0,
            my_hid_tol_start: 0.0,
            my_hid_tol_end: 0.0,
        }
    }

    /// OCCT InitHidden(status) — cxx L42-64.
    pub fn init_hidden(&mut self, status: &'a EdgeStatus) {
        self.e_hid = Some(status);
        self.i_hid = 1;
        if status.all_hidden() {
            let (s, ts, e, te) = status.bounds();
            self.my_hid_start = s;
            self.my_hid_tol_start = ts;
            self.my_hid_end = e;
            self.my_hid_tol_end = te;
            self.my_nb_hid = 0;
        } else {
            self.my_nb_hid = status.nb_visible_part() as i32;
            let (s, ts, _b1, _b2) = status.bounds();
            self.my_hid_start = s;
            self.my_hid_tol_start = ts;
            let (e, te, _, _) = status.visible_part(self.i_hid as usize);
            self.my_hid_end = e;
            self.my_hid_tol_end = te;
        }
        if self.my_hid_start + self.my_hid_tol_start as f64 >= self.my_hid_end - self.my_hid_tol_end as f64
            && self.my_hid_end + self.my_hid_tol_end as f64 >= self.my_hid_start - self.my_hid_tol_start as f64
        {
            self.next_hidden();
        }
    }

    /// OCCT MoreHidden() — lxx L21-24.
    pub fn more_hidden(&self) -> bool {
        self.i_hid <= self.my_nb_hid + 1
    }

    /// OCCT NextHidden() — cxx L68-94.
    pub fn next_hidden(&mut self) {
        if self.i_hid >= self.my_nb_hid + 1 {
            self.i_hid += 1;
        } else {
            let status = self.e_hid.expect("InitHidden not called");
            let (b1, b2, s, ts) = status.visible_part(self.i_hid as usize);
            self.my_hid_start = s;
            self.my_hid_tol_start = ts;
            self.i_hid += 1;
            if self.i_hid == self.my_nb_hid + 1 {
                let (_b1, _b2, e, te) = status.bounds();
                self.my_hid_end = e;
                self.my_hid_tol_end = te;
                if self.my_hid_start + self.my_hid_tol_start as f64
                    >= self.my_hid_end - self.my_hid_tol_end as f64
                    && self.my_hid_end + self.my_hid_tol_end as f64
                        >= self.my_hid_start - self.my_hid_tol_start as f64
                {
                    self.i_hid += 1;
                }
            } else {
                let (e, te, _, _) = status.visible_part(self.i_hid as usize);
                self.my_hid_end = e;
                self.my_hid_tol_end = te;
            }
            let _ = (b1, b2);
        }
    }

    /// OCCT Hidden(Start, TolStart, End, TolEnd) — lxx L28-37: the bounds
    /// and tolerances of the current hidden interval.
    pub fn hidden(&self) -> (f64, f32, f64, f32) {
        (
            self.my_hid_start,
            self.my_hid_tol_start,
            self.my_hid_end,
            self.my_hid_tol_end,
        )
    }

    /// OCCT InitVisible(status) — lxx L41-46.
    pub fn init_visible(&mut self, status: &'a EdgeStatus) {
        self.e_vis = Some(status);
        self.i_vis = 1;
        self.my_nb_vis = status.nb_visible_part() as i32;
    }

    /// OCCT MoreVisible() — lxx L50-53.
    pub fn more_visible(&self) -> bool {
        self.i_vis <= self.my_nb_vis
    }

    /// OCCT NextVisible() — lxx L57-60.
    pub fn next_visible(&mut self) {
        self.i_vis += 1;
    }

    /// OCCT Visible(Start, TolStart, End, TolEnd) — lxx L64-70: the bounds
    /// and tolerances of the current visible interval.
    pub fn visible(&self) -> (f64, f32, f64, f32) {
        self.e_vis
            .expect("InitVisible not called")
            .visible_part(self.i_vis as usize)
    }
}
