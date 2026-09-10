//! OCCT MAT_Bisector — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_Bisector.hxx L17-111, MAT_Bisector.cxx L25-205

use std::sync::{Arc, RwLock};

use rcad_kernel::core::precision::INFINITE_VALUE; // OCCT Precision::Infinite()

use super::mat_edge::HandleMatEdge;
use super::mat_list_of_bisector_0::{HandleMatListOfBisector, MatListOfBisector};

/// OCCT `occ::handle<MAT_Bisector>`
pub type HandleMatBisector = Arc<RwLock<MatBisector>>;

/// OCCT MAT_Bisector.hxx L29-109
pub struct MatBisector {
    thebisectornumber: i32,
    theindexnumber: i32,
    thefirstedge: Option<HandleMatEdge>,
    thesecondedge: Option<HandleMatEdge>,
    thelistofbisectors: HandleMatListOfBisector,
    theissuepoint: i32,
    theendpoint: i32,
    thefirstvector: i32,
    thesecondvector: i32,
    thesense: f64,
    thefirstparameter: f64,
    thesecondparameter: f64,
    distissuepoint: f64,
}

impl MatBisector {
    /// OCCT MAT_Bisector.cxx L25-31
    pub fn new() -> Self {
        MatBisector {
            thebisectornumber: -1,
            thefirstparameter: INFINITE_VALUE,
            thesecondparameter: INFINITE_VALUE,
            thelistofbisectors: Arc::new(RwLock::new(MatListOfBisector::new())),
            theindexnumber: 0,
            thefirstedge: None,
            thesecondedge: None,
            theissuepoint: 0,
            theendpoint: 0,
            thefirstvector: 0,
            thesecondvector: 0,
            thesense: 0.0,
            distissuepoint: 0.0,
        }
    }

    /// OCCT MAT_Bisector.cxx L33-36
    pub fn add_bisector(&self, abisector: &HandleMatBisector) {
        self.thelistofbisectors.write().unwrap().back_add(abisector);
    }

    /// OCCT MAT_Bisector.cxx L38-41
    pub fn list(&self) -> HandleMatListOfBisector {
        self.thelistofbisectors.clone()
    }

    /// OCCT MAT_Bisector.cxx L43-46
    pub fn first_bisector(&self) -> Option<HandleMatBisector> {
        self.thelistofbisectors.read().unwrap().first_item()
    }

    /// OCCT MAT_Bisector.cxx L48-51
    pub fn last_bisector(&self) -> Option<HandleMatBisector> {
        self.thelistofbisectors.read().unwrap().last_item()
    }

    /// OCCT MAT_Bisector.cxx L53-56 (setter overload of BisectorNumber)
    pub fn set_bisector_number(&mut self, anumber: i32) {
        self.thebisectornumber = anumber;
    }

    /// OCCT MAT_Bisector.cxx L58-61 (setter overload of IndexNumber)
    pub fn set_index_number(&mut self, anumber: i32) {
        self.theindexnumber = anumber;
    }

    /// OCCT MAT_Bisector.cxx L63-66 (setter overload of FirstEdge)
    pub fn set_first_edge(&mut self, anedge: &HandleMatEdge) {
        self.thefirstedge = Some(anedge.clone());
    }

    /// OCCT MAT_Bisector.cxx L68-71 (setter overload of SecondEdge)
    pub fn set_second_edge(&mut self, anedge: &HandleMatEdge) {
        self.thesecondedge = Some(anedge.clone());
    }

    // OCCT ~MAT2d_Mat2d (cxx L1878-1910): the destructor nulls the edge
    // handles — the clear entry points the MAT2d Drop path calls.
    pub fn clear_first_edge(&mut self) {
        self.thefirstedge = None;
    }

    pub fn clear_second_edge(&mut self) {
        self.thesecondedge = None;
    }

    /// OCCT MAT_Bisector.cxx L73-76 (setter overload of IssuePoint)
    pub fn set_issue_point(&mut self, apoint: i32) {
        self.theissuepoint = apoint;
    }

    /// OCCT MAT_Bisector.cxx L78-81 (setter overload of EndPoint)
    pub fn set_end_point(&mut self, apoint: i32) {
        self.theendpoint = apoint;
    }

    /// OCCT MAT_Bisector.cxx L83-86 (setter overload of DistIssuePoint)
    pub fn set_dist_issue_point(&mut self, areal: f64) {
        self.distissuepoint = areal;
    }

    /// OCCT MAT_Bisector.cxx L88-91 (setter overload of FirstVector)
    pub fn set_first_vector(&mut self, avector: i32) {
        self.thefirstvector = avector;
    }

    /// OCCT MAT_Bisector.cxx L93-96 (setter overload of SecondVector)
    pub fn set_second_vector(&mut self, avector: i32) {
        self.thesecondvector = avector;
    }

    /// OCCT MAT_Bisector.cxx L98-101 (setter overload of Sense)
    pub fn set_sense(&mut self, asense: f64) {
        self.thesense = asense;
    }

    /// OCCT MAT_Bisector.cxx L103-106 (setter overload of FirstParameter)
    pub fn set_first_parameter(&mut self, aparameter: f64) {
        self.thefirstparameter = aparameter;
    }

    /// OCCT MAT_Bisector.cxx L108-111 (setter overload of SecondParameter)
    pub fn set_second_parameter(&mut self, aparameter: f64) {
        self.thesecondparameter = aparameter;
    }

    /// OCCT MAT_Bisector.cxx L113-116 (getter overload of BisectorNumber)
    pub fn bisector_number(&self) -> i32 {
        self.thebisectornumber
    }

    /// OCCT MAT_Bisector.cxx L118-121 (getter overload of IndexNumber)
    pub fn index_number(&self) -> i32 {
        self.theindexnumber
    }

    /// OCCT MAT_Bisector.cxx L123-126 (getter overload of FirstEdge)
    pub fn first_edge(&self) -> Option<HandleMatEdge> {
        self.thefirstedge.clone()
    }

    /// OCCT MAT_Bisector.cxx L128-131 (getter overload of SecondEdge)
    pub fn second_edge(&self) -> Option<HandleMatEdge> {
        self.thesecondedge.clone()
    }

    /// OCCT MAT_Bisector.cxx L133-136 (getter overload of IssuePoint)
    pub fn issue_point(&self) -> i32 {
        self.theissuepoint
    }

    /// OCCT MAT_Bisector.cxx L138-141 (getter overload of EndPoint)
    pub fn end_point(&self) -> i32 {
        self.theendpoint
    }

    /// OCCT MAT_Bisector.cxx L143-146 (getter overload of DistIssuePoint)
    pub fn dist_issue_point(&self) -> f64 {
        self.distissuepoint
    }

    /// OCCT MAT_Bisector.cxx L148-151 (getter overload of FirstVector)
    pub fn first_vector(&self) -> i32 {
        self.thefirstvector
    }

    /// OCCT MAT_Bisector.cxx L153-156 (getter overload of SecondVector)
    pub fn second_vector(&self) -> i32 {
        self.thesecondvector
    }

    /// OCCT MAT_Bisector.cxx L158-161 (getter overload of Sense)
    pub fn sense(&self) -> f64 {
        self.thesense
    }

    /// OCCT MAT_Bisector.cxx L163-166 (getter overload of FirstParameter)
    pub fn first_parameter(&self) -> f64 {
        self.thefirstparameter
    }

    /// OCCT MAT_Bisector.cxx L168-171 (getter overload of SecondParameter)
    pub fn second_parameter(&self) -> f64 {
        self.thesecondparameter
    }

    /// OCCT MAT_Bisector.cxx L173-205
    pub fn dump(&self, ashift: i32, alevel: i32) {
        for _i in 0..ashift {
            print!("  ");
        }
        println!(" BISECTOR : {}", self.thebisectornumber);
        for _i in 0..ashift {
            print!("  ");
        }
        println!(
            "   First edge     : {}",
            self.thefirstedge
                .as_ref()
                .expect("MAT_Bisector::Dump")
                .read()
                .unwrap()
                .edge_number()
        );
        for _i in 0..ashift {
            print!("  ");
        }
        println!(
            "   Second edge    : {}",
            self.thesecondedge
                .as_ref()
                .expect("MAT_Bisector::Dump")
                .read()
                .unwrap()
                .edge_number()
        );
        for _i in 0..ashift {
            print!("  ");
        }
        if alevel != 0 {
            if !self.thelistofbisectors.read().unwrap().more() {
                println!("   Bisectors List : ");
                self.thelistofbisectors.write().unwrap().dump(ashift + 1, 1);
            }
        }
        println!();
    }
}

impl Default for MatBisector {
    fn default() -> Self {
        Self::new()
    }
}
