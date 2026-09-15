//! OCCT BRepBuilderAPI_BndBoxTreeSelector
//! (BRepBuilderAPI_BndBoxTreeSelector.hxx L20-59),
//! BRepBuilderAPI_VertexInspector (BRepBuilderAPI_VertexInspector.hxx
//! L20-66 + BRepBuilderAPI_Sewing.cxx L5912-5940 Inspect) and the
//! `NCollection_UBTree<int, Bnd_Box>` / `NCollection_UBTreeFiller<int,
//! Bnd_Box>` / `NCollection_CellFilter<BRepBuilderAPI_VertexInspector>`
//! carriers they ride (TKFoundation).
//!
//! GAP: the TKFoundation NCollection_UBTree / NCollection_CellFilter are
//! untranslated; the re-hosts keep the OCCT Add/Fill/Select/Inspect call
//! form.  `Select` walks the filled entries in insertion order invoking
//! Reject/Accept — the same documented reduction as the ShapeAnalysis
//! box_bnd_tree.rs carrier (the traversal order is the observable
//! difference; the selector semantics are unchanged).  The CellFilter keeps
//! the OCCT cell decomposition over a hash map keyed by the cell index
//! (the extrema_math_opt.rs kernel re-host convention, local to this
//! subtree).

use std::collections::HashMap;

use glam::DVec3;
use rcad_kernel::math::bnd::BndBox;

/// OCCT `NCollection_UBTree<int, Bnd_Box>` (the int/Bnd_Box instantiation
/// used by Cutting, cxx L4455).
pub struct BndBoxTree {
    /// The (object, box) pairs in insertion order.
    entries: Vec<(i32, BndBox)>,
}

/// OCCT `NCollection_UBTreeFiller<int, Bnd_Box>` — Add queues the (object,
/// box) pair; Fill commits them into the tree.
pub struct BndBoxTreeFiller {
    tree: BndBoxTree,
}

/// OCCT BRepBuilderAPI_BndBoxTreeSelector (hxx L25-58) — used to select
/// overlapping boxes; maintains the selection condition and retrieves the
/// selected objects after the search.
pub struct BndBoxTreeSelector {
    /// OCCT `NCollection_List<int> myResInd`.
    my_res_ind: Vec<i32>,
    /// OCCT `Bnd_Box myBox`.
    my_box: BndBox,
}

impl BndBoxTree {
    /// OCCT NCollection_UBTree<int, Bnd_Box>().
    pub fn new() -> Self {
        BndBoxTree { entries: Vec::new() }
    }
}

impl Default for BndBoxTree {
    fn default() -> Self {
        BndBoxTree::new()
    }
}

impl BndBoxTreeFiller {
    /// OCCT NCollection_UBTreeFiller<int, Bnd_Box>(aTree).
    pub fn new(a_tree: BndBoxTree) -> Self {
        BndBoxTreeFiller { tree: a_tree }
    }

    /// OCCT Add(theObj, theBnd).
    pub fn add(&mut self, the_obj: i32, the_bnd: BndBox) {
        self.tree.entries.push((the_obj, the_bnd));
    }

    /// OCCT Fill() — commits the queued pairs (the carrier tree stores them
    /// in insertion order).
    pub fn fill(&mut self) {}

    /// OCCT aTree.Select(aSelector) — the tree traversal invoking
    /// Reject/Accept; the accepted count is not consumed by the callers.
    pub fn select(&mut self, a_selector: &mut BndBoxTreeSelector) {
        for (the_obj, the_bnd) in &self.tree.entries {
            // OCCT hxx L38: Reject(theBox) = myBox.IsOut(theBox).
            if a_selector.my_box.is_out_box(the_bnd) {
                continue;
            }
            // OCCT hxx L45: Accept(theObj) appends and returns true.
            a_selector.my_res_ind.push(*the_obj);
        }
    }
}

impl BndBoxTreeSelector {
    /// OCCT BRepBuilderAPI_BndBoxTreeSelector() (hxx L30-32).
    pub fn new() -> Self {
        BndBoxTreeSelector {
            my_res_ind: Vec::new(),
            my_box: BndBox::new(),
        }
    }

    /// OCCT ClearResList() (hxx L51-52).
    pub fn clear_res_list(&mut self) {
        self.my_res_ind.clear();
    }

    /// OCCT SetCurrent(theBox) (hxx L54-55).
    pub fn set_current(&mut self, the_box: BndBox) {
        self.my_box = the_box;
    }

    /// OCCT ResInd() (hxx L57-58) — the indexes of the boxes intersecting
    /// the current box.
    pub fn res_ind(&self) -> &Vec<i32> {
        &self.my_res_ind
    }
}

impl Default for BndBoxTreeSelector {
    fn default() -> Self {
        BndBoxTreeSelector::new()
    }
}

// ---------------------------------------------------------------------------
// VertexInspector + NCollection_CellFilter (GlueVertices, cxx L3310+).
// ---------------------------------------------------------------------------

/// OCCT `NCollection_CellFilter_Action` (NCollection_CellFilter.hxx L26-30).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum CellFilterAction {
    /// CellFilter_Keep.
    Keep,
    /// CellFilter_Purge.
    Purge,
}

/// OCCT `BRepBuilderAPI_VertexInspector` (hxx L30-66 + cxx L5912-5940) —
/// inspector for the CellFilter algorithm working with gp_XYZ points in 3d
/// space, used in the search of coincidence points with a certain
/// tolerance.
pub struct VertexInspector {
    /// OCCT `double myTol` (the squared tolerance, cxx ctor).
    my_tol: f64,
    /// OCCT `NCollection_List<int> myResInd`.
    my_res_ind: Vec<i32>,
    /// OCCT `VectorOfPoint myPoints` (NCollection_DynamicArray<gp_XYZ>).
    my_points: Vec<DVec3>,
    /// OCCT `gp_XYZ myCurrent`.
    my_current: DVec3,
}

impl VertexInspector {
    /// OCCT BRepBuilderAPI_VertexInspector::Coord(i, thePnt) (hxx L41-42).
    #[allow(dead_code)]
    pub fn coord(i: usize, the_pnt: DVec3) -> f64 {
        the_pnt[i]
    }

    /// OCCT BRepBuilderAPI_VertexInspector::Shift(thePnt, theTol)
    /// (hxx L44-47).
    pub fn shift(the_pnt: DVec3, the_tol: f64) -> DVec3 {
        DVec3::new(the_pnt.x + the_tol, the_pnt.y + the_tol, the_pnt.z + the_tol)
    }

    /// OCCT BRepBuilderAPI_VertexInspector(theTol) (hxx L50-52): myTol is
    /// the squared tolerance.
    pub fn new(the_tol: f64) -> Self {
        VertexInspector {
            my_tol: the_tol * the_tol,
            my_res_ind: Vec::new(),
            my_points: Vec::new(),
            my_current: DVec3::ZERO,
        }
    }

    /// OCCT Add(thePnt) (hxx L55-56) — keeps the points used for comparison.
    pub fn add(&mut self, the_pnt: DVec3) {
        self.my_points.push(the_pnt);
    }

    /// OCCT ClearResList() (hxx L59-60).
    pub fn clear_res_list(&mut self) {
        self.my_res_ind.clear();
    }

    /// OCCT SetCurrent(theCurPnt) (hxx L63-64).
    pub fn set_current(&mut self, the_cur_pnt: DVec3) {
        self.my_current = the_cur_pnt;
    }

    /// OCCT ResInd() (hxx L66) — the indexes of the points adjacent with the
    /// current one.
    pub fn res_ind(&self) -> &Vec<i32> {
        &self.my_res_ind
    }

    /// OCCT Inspect(theTarget) (BRepBuilderAPI_Sewing.cxx L5912-5940) —
    /// selection and storage of coinciding points.
    pub fn inspect(&mut self, the_target: i32) -> CellFilterAction {
        // OCCT L5921: const gp_XYZ& aPnt = myPoints.Value(theTarget - 1);
        let a_pnt = self.my_points[(the_target - 1) as usize];
        // OCCT L5922-5925.
        let a_dx = self.my_current.x - a_pnt.x;
        let a_dy = self.my_current.y - a_pnt.y;
        let a_dz = self.my_current.z - a_pnt.z;

        // OCCT L5927-5931.
        if (a_dx * a_dx <= self.my_tol) && (a_dy * a_dy <= self.my_tol) && (a_dz * a_dz <= self.my_tol)
        {
            self.my_res_ind.push(the_target);
        }
        CellFilterAction::Keep
    }
}

/// OCCT `NCollection_CellFilter<Inspector>` (NCollection_CellFilter.hxx
/// L112-506) over the (Target = int, Point = DVec3, Dimension = 3) form the
/// VertexInspector instantiates; the cell map is a hash map keyed by the
/// cell index, every cell holds its targets most recently added first (hxx
/// L355-359) and Inspect walks them in that order (hxx L452-469).
pub struct CellFilter {
    /// hxx L503: NCollection_Array1<double> myCellSize.
    my_cell_size: f64,
    /// hxx L502: CellMap myCells.
    my_cells: HashMap<(i64, i64, i64), Vec<i32>>,
}

impl CellFilter {
    /// OCCT NCollection_CellFilter(theCellSize) (hxx L140-146).
    pub fn new(the_cell_size: f64) -> Self {
        CellFilter {
            my_cell_size: the_cell_size,
            my_cells: HashMap::new(),
        }
    }

    /// OCCT Cell(thePnt, theCellSize) (hxx L252-267) — the cell index of a
    /// point (the scalar cell-size form).
    fn cell_index(&self, the_pnt: DVec3) -> (i64, i64, i64) {
        const INT_MAX: f64 = i32::MAX as f64;
        const INT_MIN: f64 = i32::MIN as f64;
        let clamp = |a_val: f64| -> i64 {
            if a_val > INT_MAX - 1.0 {
                (a_val % INT_MAX) as i64
            } else if a_val < INT_MIN + 1.0 {
                (a_val % INT_MIN) as i64
            } else {
                a_val as i64
            }
        };
        (
            clamp(the_pnt.x / self.my_cell_size),
            clamp(the_pnt.y / self.my_cell_size),
            clamp(the_pnt.z / self.my_cell_size),
        )
    }

    /// OCCT Add(theTarget, thePnt) (hxx L165-169) — into the single cell
    /// that contains the point; the new target is prepended to the cell's
    /// list (hxx L355-359).
    pub fn add(&mut self, the_target: i32, the_pnt: DVec3) {
        let index = self.cell_index(the_pnt);
        let cell = self.my_cells.entry(index).or_default();
        cell.insert(0, the_target);
    }

    /// OCCT Inspect(thePntMin, thePntMax, theInspector) (hxx L217-225 + the
    /// iterateInspect walk hxx L477-497) — the inspector's Inspect is
    /// invoked on every target of every cell in the range, in the OCCT
    /// (newest-first) order.
    pub fn inspect(&mut self, the_pnt_min: DVec3, the_pnt_max: DVec3, the_inspector: &mut VertexInspector) {
        let cell_min = self.cell_index(the_pnt_min);
        let cell_max = self.cell_index(the_pnt_max);
        // OCCT iterateInspect: the 3-dimensional nested loop over the index
        // ranges (hxx L477-497).
        let mut ix = cell_min.0;
        while ix <= cell_max.0 {
            let mut iy = cell_min.1;
            while iy <= cell_max.1 {
                let mut iz = cell_min.2;
                while iz <= cell_max.2 {
                    if let Some(targets) = self.my_cells.get(&(ix, iy, iz)) {
                        for the_target in targets.clone() {
                            the_inspector.inspect(the_target);
                        }
                    }
                    iz += 1;
                }
                iy += 1;
            }
            ix += 1;
        }
    }
}
