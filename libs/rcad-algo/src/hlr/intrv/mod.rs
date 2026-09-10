//! OCCT Intrv package (TKHLR) — interval algebra on the real line.
//!
//! 1:1 translation of `Intrv_Position.hxx` + `Intrv_Interval.hxx/.cxx/.lxx`
//! + `Intrv_Intervals.hxx/.cxx/.lxx`.

/// OCCT Intrv_Position (Intrv_Position.hxx L17-31) — relative position of
/// one interval with respect to another.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Position {
    Before,
    JustBefore,
    OverlappingAtStart,
    JustEnclosingAtEnd,
    Enclosing,
    JustOverlappingAtStart,
    Similar,
    JustEnclosingAtStart,
    Inside,
    JustOverlappingAtEnd,
    OverlappingAtEnd,
    JustAfter,
    After,
}

/// OCCT Standard::Epsilon(theValue) — Standard_Real.hxx L239-246: the
/// absolute difference between the value and its nearest neighbour in the
/// direction of infinity with the same sign.
pub fn epsilon(the_value: f64) -> f64 {
    if the_value >= 0.0 {
        next_after(the_value, f64::MAX) - the_value
    } else {
        the_value - next_after(the_value, -f64::MAX)
    }
}

/// `std::nextafter` (not yet stable in Rust std) — the representable value
/// adjacent to `x` in the direction of `to`.
fn next_after(x: f64, to: f64) -> f64 {
    if x.is_nan() || to.is_nan() {
        return f64::NAN;
    }
    if x == to {
        return x;
    }
    if x == 0.0 {
        return if to > 0.0 {
            f64::MIN_POSITIVE
        } else {
            -f64::MIN_POSITIVE
        };
    }
    let bits = x.to_bits();
    let next = if (to > x) == (x > 0.0) {
        bits + 1
    } else {
        bits - 1
    };
    f64::from_bits(next)
}

/// OCCT Standard::RealFirst() (Standard_Real.hxx L170-173).
pub const REAL_FIRST: f64 = -f64::MAX;

/// OCCT Standard::RealLast() (Standard_Real.hxx L181-184).
pub const REAL_LAST: f64 = f64::MAX;

/// OCCT Intrv_Interval — a real interval with independent start/end
/// tolerances (Intrv_Interval.hxx L37-188 + .cxx L35-151 + .lxx L14-325).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    my_start: f64,
    my_end: f64,
    my_tol_start: f32,
    my_tol_end: f32,
}

impl Default for Interval {
    fn default() -> Self {
        Self::new()
    }
}

impl Interval {
    /// OCCT Intrv_Interval() — Intrv_Interval.cxx L35-41: the whole real
    /// line with infinite tolerances.
    pub fn new() -> Self {
        Interval {
            my_start: REAL_FIRST,
            my_end: REAL_LAST,
            // myTolStart = (float)Epsilon(RealFirst()); — the nextafter
            // distance past -DBL_MAX is +inf in both OCCT and Rust.
            my_tol_start: epsilon(REAL_FIRST) as f32,
            my_tol_end: epsilon(REAL_LAST) as f32,
        }
    }

    /// OCCT Intrv_Interval(Start, End) — Intrv_Interval.cxx L45-51.
    pub fn from_bounds(start: f64, end: f64) -> Self {
        Interval {
            my_start: start,
            my_end: end,
            my_tol_start: epsilon(start) as f32,
            my_tol_end: epsilon(end) as f32,
        }
    }

    /// OCCT Intrv_Interval(Start, TolStart, End, TolEnd) —
    /// Intrv_Interval.cxx L55-74: tolerances are raised to at least the
    /// representable epsilon of their bound.
    pub fn from_bounds_tol(start: f64, tol_start: f32, end: f64, tol_end: f32) -> Self {
        let mut interval = Interval {
            my_start: start,
            my_end: end,
            my_tol_start: tol_start,
            my_tol_end: tol_end,
        };
        let eps_start = epsilon(interval.my_start) as f32;
        let eps_end = epsilon(interval.my_end) as f32;
        if interval.my_tol_start < eps_start {
            interval.my_tol_start = eps_start;
        }
        if interval.my_tol_end < eps_end {
            interval.my_tol_end = eps_end;
        }
        interval
    }

    /// OCCT Start() — Intrv_Interval.lxx L24-31.
    pub fn start(&self) -> f64 {
        self.my_start
    }

    /// OCCT End() — Intrv_Interval.lxx L36-43.
    pub fn end(&self) -> f64 {
        self.my_end
    }

    /// OCCT TolStart() — Intrv_Interval.lxx L48-55.
    pub fn tol_start(&self) -> f32 {
        self.my_tol_start
    }

    /// OCCT TolEnd() — Intrv_Interval.lxx L60-67.
    pub fn tol_end(&self) -> f32 {
        self.my_tol_end
    }

    /// OCCT Bounds() — Intrv_Interval.lxx L72-80.
    pub fn bounds(&self) -> (f64, f32, f64, f32) {
        (
            self.my_start,
            self.my_tol_start,
            self.my_end,
            self.my_tol_end,
        )
    }

    /// OCCT SetStart() — Intrv_Interval.lxx L85-92.
    pub fn set_start(&mut self, start: f64, tol_start: f32) {
        self.my_start = start;
        self.my_tol_start = tol_start;
    }

    /// OCCT FuseAtStart() — Intrv_Interval.lxx L97-118.
    pub fn fuse_at_start(&mut self, start: f64, tol_start: f32) {
        if self.my_start != REAL_FIRST {
            let a = (self.my_start - self.my_tol_start as f64).min(start - tol_start as f64);
            let b = (self.my_start + self.my_tol_start as f64).min(start + tol_start as f64);
            self.my_start = (a + b) / 2.0;
            self.my_tol_start = ((b - a) / 2.0) as f32;
        }
    }

    /// OCCT CutAtStart() — Intrv_Interval.lxx L123-144.
    pub fn cut_at_start(&mut self, start: f64, tol_start: f32) {
        if self.my_start != REAL_FIRST {
            let a = (self.my_start - self.my_tol_start as f64).max(start - tol_start as f64);
            let b = (self.my_start + self.my_tol_start as f64).max(start + tol_start as f64);
            self.my_start = (a + b) / 2.0;
            self.my_tol_start = ((b - a) / 2.0) as f32;
        }
    }

    /// OCCT SetEnd() — Intrv_Interval.lxx L149-156.
    pub fn set_end(&mut self, end: f64, tol_end: f32) {
        self.my_end = end;
        self.my_tol_end = tol_end;
    }

    /// OCCT FuseAtEnd() — Intrv_Interval.lxx L161-182.
    pub fn fuse_at_end(&mut self, end: f64, tol_end: f32) {
        if self.my_end != REAL_LAST {
            let a = (self.my_end - self.my_tol_end as f64).max(end - tol_end as f64);
            let b = (self.my_end + self.my_tol_end as f64).max(end + tol_end as f64);
            self.my_end = (a + b) / 2.0;
            self.my_tol_end = ((b - a) / 2.0) as f32;
        }
    }

    /// OCCT CutAtEnd() — Intrv_Interval.lxx L187-208.
    pub fn cut_at_end(&mut self, end: f64, tol_end: f32) {
        if self.my_end != REAL_LAST {
            let a = (self.my_end - self.my_tol_end as f64).min(end - tol_end as f64);
            let b = (self.my_end + self.my_tol_end as f64).min(end + tol_end as f64);
            self.my_end = (a + b) / 2.0;
            self.my_tol_end = ((b - a) / 2.0) as f32;
        }
    }

    /// OCCT IsProbablyEmpty() — Intrv_Interval.lxx L213-223.
    pub fn is_probably_empty(&self) -> bool {
        are_fused(
            self.my_start,
            self.my_tol_start,
            self.my_end,
            self.my_tol_end,
        )
    }

    /// OCCT Position(Other) — Intrv_Interval.cxx L78-151.
    pub fn position(&self, other: &Interval) -> Position {
        let my_s_min = self.my_start - self.my_tol_start as f64;
        let my_s_max = self.my_start + self.my_tol_start as f64;
        let my_e_min = self.my_end - self.my_tol_end as f64;
        let my_e_max = self.my_end + self.my_tol_end as f64;
        let ot_s_min = other.my_start - other.my_tol_start as f64;
        let ot_s_max = other.my_start + other.my_tol_start as f64;
        let ot_e_min = other.my_end - other.my_tol_end as f64;
        let ot_e_max = other.my_end + other.my_tol_end as f64;
        let p;
        if my_s_max < ot_s_min {
            if my_e_max < ot_s_min {
                p = Position::Before;
            } else if ot_s_max >= my_e_min {
                p = Position::JustBefore;
            } else if my_e_max < ot_e_min {
                p = Position::OverlappingAtStart;
            } else if ot_e_max >= my_e_min {
                p = Position::JustEnclosingAtEnd;
            } else {
                p = Position::Enclosing;
            }
        } else if ot_s_max >= my_s_min {
            if my_e_max < ot_e_min {
                p = Position::JustOverlappingAtStart;
            } else if ot_e_max >= my_e_min {
                p = Position::Similar;
            } else {
                p = Position::JustEnclosingAtStart;
            }
        } else if my_s_max < ot_e_min {
            if my_e_max < ot_e_min {
                p = Position::Inside;
            } else if ot_e_max >= my_e_min {
                p = Position::JustOverlappingAtEnd;
            } else {
                p = Position::OverlappingAtEnd;
            }
        } else if ot_e_max >= my_s_min {
            p = Position::JustAfter;
        } else {
            p = Position::After;
        }
        p
    }

    /// OCCT IsBefore(Other) — Intrv_Interval.lxx L228-238.
    pub fn is_before(&self, other: &Interval) -> bool {
        (self.my_tol_end as f64 + other.my_tol_start as f64) < other.my_start - self.my_end
    }

    /// OCCT IsAfter(Other) — Intrv_Interval.lxx L243-253.
    pub fn is_after(&self, other: &Interval) -> bool {
        (self.my_tol_start as f64 + other.my_tol_end as f64) < self.my_start - other.my_end
    }

    /// OCCT IsInside(Other) — Intrv_Interval.lxx L258-269.
    pub fn is_inside(&self, other: &Interval) -> bool {
        (self.my_tol_start as f64 + other.my_tol_start as f64) < self.my_start - other.my_start
            && (self.my_tol_end as f64 + other.my_tol_end as f64) < other.my_end - self.my_end
    }

    /// OCCT IsEnclosing(Other) — Intrv_Interval.lxx L274-285.
    pub fn is_enclosing(&self, other: &Interval) -> bool {
        (self.my_tol_start as f64 + other.my_tol_start as f64) < other.my_start - self.my_start
            && (self.my_tol_end as f64 + other.my_tol_end as f64) < self.my_end - other.my_end
    }

    /// OCCT IsJustEnclosingAtStart(Other) — Intrv_Interval.lxx L290-301.
    pub fn is_just_enclosing_at_start(&self, other: &Interval) -> bool {
        are_fused(
            self.my_start,
            self.my_tol_start,
            other.my_start,
            other.my_tol_start,
        ) && (self.my_tol_end as f64 + other.my_tol_end as f64) < self.my_end - other.my_end
    }

    /// OCCT IsJustEnclosingAtEnd(Other) — Intrv_Interval.lxx L306-317.
    pub fn is_just_enclosing_at_end(&self, other: &Interval) -> bool {
        (self.my_tol_start as f64 + other.my_tol_start as f64) < other.my_start - self.my_start
            && are_fused(
                other.my_end,
                other.my_tol_end,
                self.my_end,
                self.my_tol_end,
            )
    }

    /// OCCT IsJustBefore(Other) — Intrv_Interval.lxx L322-332.
    pub fn is_just_before(&self, other: &Interval) -> bool {
        are_fused(self.my_end, self.my_tol_end, other.my_start, other.my_tol_start)
    }

    /// OCCT IsJustAfter(Other) — Intrv_Interval.lxx L337-347.
    pub fn is_just_after(&self, other: &Interval) -> bool {
        are_fused(other.my_end, other.my_tol_end, self.my_start, self.my_tol_start)
    }

    /// OCCT IsOverlappingAtStart(Other) — Intrv_Interval.lxx L352-363.
    pub fn is_overlapping_at_start(&self, other: &Interval) -> bool {
        (self.my_tol_start as f64 + other.my_tol_start as f64) < other.my_start - self.my_start
            && (self.my_tol_end as f64 + other.my_tol_start as f64) < self.my_end - other.my_start
            && (self.my_tol_end as f64 + other.my_tol_end as f64) < other.my_end - self.my_end
    }

    /// OCCT IsOverlappingAtEnd(Other) — Intrv_Interval.lxx L368-379.
    pub fn is_overlapping_at_end(&self, other: &Interval) -> bool {
        (self.my_tol_start as f64 + other.my_tol_start as f64) < self.my_start - other.my_start
            && (self.my_tol_start as f64 + other.my_tol_end as f64) < other.my_end - self.my_start
            && (self.my_tol_end as f64 + other.my_tol_end as f64) < self.my_end - other.my_end
    }

    /// OCCT IsJustOverlappingAtStart(Other) — Intrv_Interval.lxx L384-395.
    pub fn is_just_overlapping_at_start(&self, other: &Interval) -> bool {
        are_fused(
            self.my_start,
            self.my_tol_start,
            other.my_start,
            other.my_tol_start,
        ) && (self.my_tol_end as f64 + other.my_tol_end as f64) < other.my_end - self.my_end
    }

    /// OCCT IsJustOverlappingAtEnd(Other) — Intrv_Interval.lxx L400-411.
    pub fn is_just_overlapping_at_end(&self, other: &Interval) -> bool {
        (self.my_tol_start as f64 + other.my_tol_start as f64) < self.my_start - other.my_start
            && are_fused(other.my_end, other.my_tol_end, self.my_end, self.my_tol_end)
    }

    /// OCCT IsSimilar(Other) — Intrv_Interval.lxx L416-429.
    pub fn is_similar(&self, other: &Interval) -> bool {
        let b1 = are_fused(
            self.my_start,
            self.my_tol_start,
            other.my_start,
            other.my_tol_start,
        );
        let b2 = are_fused(self.my_end, self.my_tol_end, other.my_end, other.my_tol_end);
        b1 && b2
    }
}

/// OCCT AreFused(c1, t1, c2, t2) — Intrv_Interval.lxx L212-218.
fn are_fused(c1: f64, t1: f32, c2: f64, t2: f32) -> bool {
    t1 as f64 + t2 as f64 >= (c1 - c2).abs()
}

/// OCCT Intrv_Intervals — a sorted sequence of non-overlapping intervals
/// (Intrv_Intervals.hxx L30-68 + .cxx L39-248 + .lxx L20-33).
#[derive(Debug, Clone, Default)]
pub struct Intervals {
    /// NCollection_Sequence<Intrv_Interval> — 1-based in OCCT, stored as a
    /// plain Vec here.
    my_inter: Vec<Interval>,
}

impl Intervals {
    /// OCCT Intrv_Intervals() — Intrv_Intervals.cxx L39.
    pub fn new() -> Self {
        Intervals { my_inter: Vec::new() }
    }

    /// OCCT Intrv_Intervals(Int) — Intrv_Intervals.cxx L43-46.
    pub fn from_interval(int: Interval) -> Self {
        Intervals {
            my_inter: vec![int],
        }
    }

    /// OCCT Intersect(Tool) — Intrv_Intervals.cxx L50-54.
    pub fn intersect(&mut self, tool: &Interval) {
        let mut inter = Intervals::from_interval(*tool);
        self.intersect_intervals(&inter);
    }

    /// OCCT Intersect(Tool) — Intrv_Intervals.cxx L58-64.
    pub fn intersect_intervals(&mut self, tool: &Intervals) {
        let mut x_uni = self.clone();
        x_uni.x_unite_intervals(tool);
        self.unite_intervals(tool);
        self.subtract_intervals(&x_uni);
    }

    /// OCCT Subtract(Tool) — Intrv_Intervals.cxx L68-122.
    pub fn subtract(&mut self, tool: &Interval) {
        let mut index: i32 = 1;

        while index <= self.my_inter.len() as i32 {
            match tool.position(&self.my_inter[(index - 1) as usize]) {
                Position::Before => {
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::JustBefore => {
                    // modifier le debut
                    self.my_inter[(index - 1) as usize]
                        .cut_at_start(tool.end(), tool.tol_end());
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::OverlappingAtStart | Position::JustOverlappingAtStart => {
                    // garder la fin
                    self.my_inter[(index - 1) as usize]
                        .set_start(tool.end(), tool.tol_end());
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::JustEnclosingAtEnd
                | Position::Enclosing
                | Position::Similar
                | Position::JustEnclosingAtStart => {
                    self.my_inter.remove((index - 1) as usize); // detruire et
                    index -= 1; // continuer
                }
                Position::Inside => {
                    let copy = self.my_inter[(index - 1) as usize];
                    self.my_inter.insert(index as usize, copy);
                    // garder le debut
                    self.my_inter[(index - 1) as usize]
                        .set_end(tool.start(), tool.tol_start());
                    // garder la fin
                    self.my_inter[index as usize].set_start(tool.end(), tool.tol_end());
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::JustOverlappingAtEnd | Position::OverlappingAtEnd => {
                    // garder le debut
                    self.my_inter[(index - 1) as usize]
                        .set_end(tool.start(), tool.tol_start());
                    // continuer
                }
                Position::JustAfter => {
                    // modifier la fin
                    self.my_inter[(index - 1) as usize]
                        .cut_at_end(tool.start(), tool.tol_start());
                    // continuer
                }
                Position::After => {
                    // continuer
                }
            }
            index += 1;
        }
    }

    /// OCCT Subtract(Tool) — Intrv_Intervals.cxx L126-133.
    pub fn subtract_intervals(&mut self, tool: &Intervals) {
        for index in 1..=tool.my_inter.len() as i32 {
            self.subtract(&tool.my_inter[(index - 1) as usize]);
        }
    }

    /// OCCT Unite(Tool) — Intrv_Intervals.cxx L137-219.
    pub fn unite(&mut self, tool: &Interval) {
        let mut inserted = false;
        let mut tins = *tool;
        let mut index: i32 = 1;

        while index <= self.my_inter.len() as i32 {
            match tins.position(&self.my_inter[(index - 1) as usize]) {
                Position::Before => {
                    inserted = true;
                    // inserer avant et
                    self.my_inter.insert((index - 1) as usize, tins);
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::JustBefore | Position::OverlappingAtStart => {
                    inserted = true;
                    // changer le debut
                    let (s, ts) = (tins.start(), tins.tol_start());
                    self.my_inter[(index - 1) as usize].set_start(s, ts);
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::Similar => {
                    // modifier le debut
                    let (s, ts) = (
                        self.my_inter[(index - 1) as usize].start(),
                        self.my_inter[(index - 1) as usize].tol_start(),
                    );
                    tins.fuse_at_start(s, ts);
                    // modifier la fin (OCCT fallthrough into JustEnclosingAtEnd)
                    let (e, te) = (
                        self.my_inter[(index - 1) as usize].end(),
                        self.my_inter[(index - 1) as usize].tol_end(),
                    );
                    tins.fuse_at_end(e, te);
                    self.my_inter.remove((index - 1) as usize); // detruire et
                    index -= 1; // continuer
                }
                Position::JustEnclosingAtEnd => {
                    // modifier la fin
                    let (e, te) = (
                        self.my_inter[(index - 1) as usize].end(),
                        self.my_inter[(index - 1) as usize].tol_end(),
                    );
                    tins.fuse_at_end(e, te);
                    self.my_inter.remove((index - 1) as usize); // detruire et
                    index -= 1; // continuer
                }
                Position::Enclosing => {
                    self.my_inter.remove((index - 1) as usize); // detruire et
                    index -= 1; // continuer
                }
                Position::JustOverlappingAtEnd => {
                    // changer le debut
                    let (s, ts) = (
                        self.my_inter[(index - 1) as usize].start(),
                        self.my_inter[(index - 1) as usize].tol_start(),
                    );
                    tins.set_start(s, ts);
                    // modifier la fin
                    let (e, te) = (
                        self.my_inter[(index - 1) as usize].end(),
                        self.my_inter[(index - 1) as usize].tol_end(),
                    );
                    tins.fuse_at_end(e, te);
                    self.my_inter.remove((index - 1) as usize); // detruire et
                    index -= 1; // continuer
                }
                Position::JustOverlappingAtStart => {
                    inserted = true;
                    // modifier le debut
                    let (s, ts) = (tins.start(), tins.tol_start());
                    self.my_inter[(index - 1) as usize].fuse_at_start(s, ts);
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::JustEnclosingAtStart => {
                    // modifier le debut
                    let (s, ts) = (
                        self.my_inter[(index - 1) as usize].start(),
                        self.my_inter[(index - 1) as usize].tol_start(),
                    );
                    tins.fuse_at_start(s, ts);
                    self.my_inter.remove((index - 1) as usize); // detruire et
                    index -= 1; // continuer
                }
                Position::Inside => {
                    inserted = true;
                    index = self.my_inter.len() as i32; // sortir
                }
                Position::OverlappingAtEnd | Position::JustAfter => {
                    // changer le debut
                    let (s, ts) = (
                        self.my_inter[(index - 1) as usize].start(),
                        self.my_inter[(index - 1) as usize].tol_start(),
                    );
                    tins.set_start(s, ts);
                    self.my_inter.remove((index - 1) as usize); // detruire et
                    index -= 1; // continuer
                }
                Position::After => {
                    // continuer
                }
            }
            index += 1;
        }
        if !inserted {
            self.my_inter.push(tins);
        }
    }

    /// OCCT Unite(Tool) — Intrv_Intervals.cxx L223-230.
    pub fn unite_intervals(&mut self, tool: &Intervals) {
        for index in 1..=tool.my_inter.len() as i32 {
            self.unite(&tool.my_inter[(index - 1) as usize]);
        }
    }

    /// OCCT XUnite(Tool) — Intrv_Intervals.cxx L234-238.
    pub fn x_unite(&mut self, tool: &Interval) {
        let mut inter = Intervals::from_interval(*tool);
        self.x_unite_intervals(&inter);
    }

    /// OCCT XUnite(Tool) — Intrv_Intervals.cxx L242-248.
    pub fn x_unite_intervals(&mut self, tool: &Intervals) {
        let mut sub2 = tool.clone();
        sub2.subtract_intervals(self);
        self.subtract_intervals(tool);
        self.unite_intervals(&sub2);
    }

    /// OCCT NbIntervals() — Intrv_Intervals.lxx L20-28.
    pub fn nb_intervals(&self) -> usize {
        self.my_inter.len()
    }

    /// OCCT Value(Index) — Intrv_Intervals.lxx L30-38 (1-based).
    pub fn value(&self, index: usize) -> &Interval {
        &self.my_inter[index - 1]
    }
}
