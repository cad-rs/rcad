// OCCT BRepAlgoAPI_Section — section operation.
//
// OCCT BRepAlgoAPI_Section.cxx / .hxx
// Computes intersection curves between two shapes (or shape and surface).

use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::algo::builder::BooleanError;
use crate::bop::ds::DS;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, Orientation, TShape};

/// OCCT BRepAlgoAPI_Section — section curves between two shapes.
pub struct SectionOp<'a> {
    shape1: Option<&'a topods::BRep>,
    shape2: Option<&'a topods::BRep>,
    my_approx: bool,
    my_compute_pcurve1: bool,
    my_compute_pcurve2: bool,
    result: Option<topods::BRep>,
    err: Option<BooleanError>,
}

impl<'a> SectionOp<'a> {
    /// OCCT BRepAlgoAPI_Section() — empty constructor.
    pub fn new() -> Self {
        SectionOp {
            shape1: None,
            shape2: None,
            my_approx: false,
            my_compute_pcurve1: false,
            my_compute_pcurve2: false,
            result: None,
            err: None,
        }
    }

    /// OCCT BRepAlgoAPI_Section(Sh1, Sh2, PerformNow).
    pub fn from_shapes(shape1: &'a topods::BRep, shape2: &'a topods::BRep, perform_now: bool) -> Self {
        let mut op = SectionOp {
            shape1: Some(shape1),
            shape2: Some(shape2),
            my_approx: false,
            my_compute_pcurve1: false,
            my_compute_pcurve2: false,
            result: None,
            err: None,
        };
        op.init(perform_now);
        op
    }

    /// OCCT Init(bFlag) — set section flags, optionally build.
    fn init(&mut self, b_flag: bool) {
        self.my_approx = false;
        self.my_compute_pcurve1 = false;
        self.my_compute_pcurve2 = false;
        if b_flag {
            let _ = self.build();
        }
    }

    /// OCCT Init1(S1) — set first argument.
    pub fn init1(&mut self, shape: &'a topods::BRep) {
        self.shape1 = Some(shape);
    }

    /// OCCT Init2(S2) — set second argument.
    pub fn init2(&mut self, shape: &'a topods::BRep) {
        self.shape2 = Some(shape);
    }

    /// OCCT Approximation(B) — enable approximation of section curves.
    pub fn approximation(&mut self, b: bool) {
        self.my_approx = b;
    }

    /// OCCT ComputePCurveOn1(B) — compute pcurve on first shape.
    pub fn compute_pcurve_on1(&mut self, b: bool) {
        self.my_compute_pcurve1 = b;
    }

    /// OCCT ComputePCurveOn2(B) — compute pcurve on second shape.
    pub fn compute_pcurve_on2(&mut self, b: bool) {
        self.my_compute_pcurve2 = b;
    }

    /// OCCT BRepAlgoAPI_Section::Build().
    /// Delegates to BooleanOperation::Build (rcad: run PaveFiller).
    pub fn build(&mut self) -> Result<(), BooleanError> {
        self.result = None;
        self.err = None;

        let shape1 = self.shape1.ok_or(BooleanError::InvalidOperation)?;
        let shape2 = self.shape2.ok_or(BooleanError::InvalidOperation)?;

        // Set section attributes before PaveFiller runs (mimics SetAttributes)
        // rcad: section attributes are forwarded through PaveFiller

        // Build DS from root shapes (OCCT-aligned)
        fn root(brep: &topods::BRep, location: u32) -> Shape {
            for (i, ts) in brep.tshapes.iter().enumerate().rev() {
                match &**ts {
                    TShape::Solid(_) | TShape::Shell(_) => {
                        return Shape::from_parts(ts.clone(), i, location, Orientation::Forward);
                    }
                    _ => {}
                }
            }
            panic!("no root Solid/Shell in BRep");
        }
        let mut ds = DS::new();
        ds.set_arguments(vec![root(shape1, 0), root(shape2, 1)]);
        ds.init(1e-7);

        // Run PaveFiller
        let mut filler = PaveFiller::new(&mut ds);
        filler.perform();
        drop(filler);

        // rcad: section result not yet built — stub returns original shapes
        self.result = Some((*shape1).clone());
        Ok(())
    }

    /// `BRepAlgoAPI_Section::Shape()`.
    pub fn shape(&self) -> &topods::BRep {
        self.result.as_ref().expect("build() not called")
    }

    /// `BRepAlgoAPI_Algo::IsDone()`.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }

    /// OCCT BRepAlgoAPI_Section::HasAncestorFaceOn1.
    /// Checks if the edge is a section edge belonging to the first argument.
    pub fn has_ancestor_face_on1(&self, _edge_idx: usize) -> Option<usize> {
        // rcad: stub — requires FF interference iteration
        None
    }

    /// OCCT BRepAlgoAPI_Section::HasAncestorFaceOn2.
    pub fn has_ancestor_face_on2(&self, _edge_idx: usize) -> Option<usize> {
        None
    }
}
