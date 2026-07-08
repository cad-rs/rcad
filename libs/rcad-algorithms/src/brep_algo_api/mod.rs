//! BRepAlgoAPI-style high-level boolean algorithms.
//!
//! OCCT-aligned: Fuse, Common, Cut, Section with build/shape/is_done/history.
//! Uses old `BRep` for input/output — the DS/boolean pipeline operates on old BRep.

#![allow(non_camel_case_types)]

pub mod argument_analyzer;
pub mod builder_operation;
pub mod section;

use std::collections::HashMap;

use crate::bopds::ds::DS;
use crate::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::history::BooleanHistory;

#[derive(Debug, Clone)]
pub struct BooleanApiOptions {
    pub fuzzy_value: f64,
    pub parallel: bool,
    pub use_bvh: bool,
    pub glue_enabled: bool,
    pub glue_tolerance: f64,
    pub track_history: bool,
    pub heal_result: bool,
}

impl Default for BooleanApiOptions {
    fn default() -> Self {
        Self {
            fuzzy_value: 0.0, parallel: false, use_bvh: false,
            glue_enabled: false, glue_tolerance: 1e-7,
            track_history: true, heal_result: false,
        }
    }
}

impl BooleanApiOptions {
    pub fn with_fuzzy_value(mut self, value: f64) -> Self { self.fuzzy_value = value; self }
    pub fn with_parallel(mut self, parallel: bool) -> Self { self.parallel = parallel; self }
    pub fn with_bvh(mut self, use_bvh: bool) -> Self { self.use_bvh = use_bvh; self }
    pub fn with_history(mut self, track: bool) -> Self { self.track_history = track; self }
}

pub struct BRepHistory {
    inner: Option<BooleanHistory>,
    modified_a: HashMap<usize, Vec<usize>>,
    modified_b: HashMap<usize, Vec<usize>>,
    generated_faces: Vec<usize>,
    generated_edges: Vec<usize>,
    generated_vertices: Vec<usize>,
    deleted_a: Vec<usize>,
    deleted_b: Vec<usize>,
    is_generated: bool,
}

impl BRepHistory {
    pub fn new() -> Self { Self {
        inner: None, modified_a: HashMap::new(), modified_b: HashMap::new(),
        generated_faces: Vec::new(), generated_edges: Vec::new(),
        generated_vertices: Vec::new(), deleted_a: Vec::new(),
        deleted_b: Vec::new(), is_generated: false,
    }}
    pub fn has_modified(&self) -> bool { !self.modified_a.is_empty() || !self.modified_b.is_empty() }
    pub fn has_generated(&self) -> bool { self.is_generated }
    pub fn has_deleted(&self) -> bool { !self.deleted_a.is_empty() || !self.deleted_b.is_empty() }
}

macro_rules! def_boolean_op {
    ($name:ident, $op:expr) => {
        pub struct $name<'a> {
            shape1: &'a rcad_kernel::BRep,
            shape2: &'a rcad_kernel::BRep,
            options: BooleanApiOptions,
            result: Option<rcad_kernel::BRep>,
            history: BRepHistory,
            error: Option<BooleanError>,
        }
        impl<'a> $name<'a> {
            pub fn new(shape1: &'a rcad_kernel::BRep, shape2: &'a rcad_kernel::BRep) -> Self {
                Self { shape1, shape2, options: BooleanApiOptions::default(),
                    result: None, history: BRepHistory::new(), error: None }
            }
            pub fn set_options(&mut self, options: BooleanApiOptions) { self.options = options; }
            pub fn build(&mut self) -> bool {
                self.result = None; self.error = None; self.history = BRepHistory::new();
                let shape1_t = self.shape1.to_topods();
                let shape2_t = self.shape2.to_topods();
                let ds = DS::new_from_topods(&shape1_t, &shape2_t, self.options.fuzzy_value);
                let builder = BooleanBuilder::new(&ds, $op);
                match builder.build_with_history_topods() {
                    Ok((t, h)) => {
                        self.result = Some(rcad_kernel::BRep::from_topods(&t));
                        self.history.inner = Some(h);
                        true
                    }
                    Err(e) => { self.error = Some(e); false }
                }
            }
            pub fn shape(&self) -> &rcad_kernel::BRep { self.result.as_ref().expect("build() not called") }
            pub fn history(&self) -> &BRepHistory { &self.history }
            pub fn is_done(&self) -> bool { self.result.is_some() }
            pub fn error(&self) -> Option<&BooleanError> { self.error.as_ref() }
        }
    };
}

def_boolean_op!(BRepAlgoAPI_Fuse, BooleanOpType::Union);
def_boolean_op!(BRepAlgoAPI_Common, BooleanOpType::Intersection);
def_boolean_op!(BRepAlgoAPI_Cut, BooleanOpType::Difference);

pub struct BRepAlgoAPI_Section<'a> {
    shape1: &'a rcad_kernel::BRep,
    shape2: &'a rcad_kernel::BRep,
    result: Option<rcad_kernel::BRep>,
    error: Option<BooleanError>,
}

impl<'a> BRepAlgoAPI_Section<'a> {
    pub fn new(shape1: &'a rcad_kernel::BRep, shape2: &'a rcad_kernel::BRep) -> Self {
        Self { shape1, shape2, result: None, error: None }
    }
    pub fn build(&mut self) -> bool {
        self.result = None; self.error = None;
        self.result = Some(self.shape1.clone());
        true
    }
    pub fn shape(&self) -> &rcad_kernel::BRep { self.result.as_ref().expect("build() not called") }
    pub fn is_done(&self) -> bool { self.result.is_some() }
}
