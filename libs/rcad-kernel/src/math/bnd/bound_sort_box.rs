//! OCCT Bnd_BoundSortBox (FoundationClasses TKMath Bnd package).
//!
//! 1:1 translation of `Bnd_BoundSortBox.hxx` (L51-148) +
//! `Bnd_BoundSortBox.cxx` (L102-624), including the file-local
//! `Bnd_VoxelGrid` helper class and `getBnd_VoxelGridResolution`.
//!
//! `NCollection_HArray1<Bnd_Box>` (1-based) maps to `Vec<BndBox>`,
//! `NCollection_List<int>` (Compare result) to `Vec<i32>`, and the OCCT
//! incremental-allocation vectors to plain `Vec<i32>`.
//!
//! OCCT `Compare(const gp_Pln&)` (cxx L508-521) is deferred: it requires
//! `Bnd_Box::IsOut(gp_Pln)`, which rcad has not ported yet; no HLR consumer
//! needs it (HLRBRep_InterCSurf only calls the box overload).

use glam::DVec3;

use super::BndBox;

/// OCCT static getBnd_VoxelGridResolution (cxx L102-121).
fn get_bnd_voxel_grid_resolution(the_boxes_count: i32) -> i32 {
    if the_boxes_count > 40000 {
        return 128;
    }
    if the_boxes_count > 10000 {
        return 64;
    }
    if the_boxes_count > 1000 {
        return 32;
    }
    if the_boxes_count > 100 {
        return 16;
    }
    8
}

/// OCCT class Bnd_VoxelGrid (cxx L129-307): a 3D voxel grid for spatial
/// partitioning — per-axis slices of box-index lists.
#[derive(Debug, Clone, Default)]
struct VoxelGrid {
    /// SliceArray mySlicesX/Y/Z over voxel indices [0, resolution-1].
    my_slices_x: Vec<Vec<i32>>,
    my_slices_y: Vec<Vec<i32>>,
    my_slices_z: Vec<Vec<i32>>,
}

impl VoxelGrid {
    /// OCCT Bnd_VoxelGrid(theResolution, theExpectedBoxCount) — cxx
    /// L216-231 (the presized vectors are plain Vec here).
    fn new(the_resolution: usize, _the_expected_box_count: i32) -> Self {
        VoxelGrid {
            my_slices_x: vec![Vec::new(); the_resolution],
            my_slices_y: vec![Vec::new(); the_resolution],
            my_slices_z: vec![Vec::new(); the_resolution],
        }
    }

    /// OCCT AddBox(theBoxIndex, theVoxelBox) — cxx L235-247.
    fn add_box(&mut self, the_box_index: i32, the_voxel_box: &[i32; 6]) {
        let (a_min_x, a_min_y, a_min_z) = (the_voxel_box[0], the_voxel_box[1], the_voxel_box[2]);
        let (a_max_x, a_max_y, a_max_z) = (the_voxel_box[3], the_voxel_box[4], the_voxel_box[5]);

        // AppendSliceX/Y/Z (cxx L275-307).
        for i in a_min_x..=a_max_x {
            self.my_slices_x[i as usize].push(the_box_index);
        }
        for i in a_min_y..=a_max_y {
            self.my_slices_y[i as usize].push(the_box_index);
        }
        for i in a_min_z..=a_max_z {
            self.my_slices_z[i as usize].push(the_box_index);
        }
    }

    /// OCCT GetSliceX — cxx L251-255 (None when the slice is empty).
    fn get_slice_x(&self, the_voxel_index: i32) -> Option<&Vec<i32>> {
        let s = &self.my_slices_x[the_voxel_index as usize];
        if s.is_empty() { None } else { Some(s) }
    }

    /// OCCT GetSliceY — cxx L259-263.
    fn get_slice_y(&self, the_voxel_index: i32) -> Option<&Vec<i32>> {
        let s = &self.my_slices_y[the_voxel_index as usize];
        if s.is_empty() { None } else { Some(s) }
    }

    /// OCCT GetSliceZ — cxx L267-271.
    fn get_slice_z(&self, the_voxel_index: i32) -> Option<&Vec<i32>> {
        let s = &self.my_slices_z[the_voxel_index as usize];
        if s.is_empty() { None } else { Some(s) }
    }
}

/// OCCT Bnd_BoundSortBox — compares a bounding box against a sorted set of
/// bounding boxes via a voxel grid over the enclosing box.
#[derive(Debug, Clone, Default)]
pub struct BoundSortBox {
    my_enclosing_box: BndBox,
    /// NCollection_HArray1<Bnd_Box> (1-based in OCCT; index 0 unused here).
    my_boxes: Option<Vec<BndBox>>,
    my_coeff_x: f64,
    my_coeff_y: f64,
    my_coeff_z: f64,
    my_resolution: i32,
    my_last_result: Vec<i32>,
    my_large_boxes: Vec<i32>,
    my_voxel_grid: Option<VoxelGrid>,
}

impl BoundSortBox {
    /// OCCT Bnd_BoundSortBox() — cxx L311-320.
    pub fn new() -> Self {
        BoundSortBox {
            my_enclosing_box: BndBox::new(),
            my_boxes: None,
            my_coeff_x: 0.0,
            my_coeff_y: 0.0,
            my_coeff_z: 0.0,
            my_resolution: 0,
            my_last_result: Vec::new(),
            my_large_boxes: Vec::new(),
            my_voxel_grid: None,
        }
    }

    /// OCCT Initialize(theSetOfBoxes) — cxx L324-348.
    pub fn initialize(&mut self, the_set_of_boxes: Vec<BndBox>) {
        // myBoxes = theSetOfBoxes (the Vec IS the 1-based array shifted).
        let mut boxes = the_set_of_boxes;
        // The OCCT array is 1-based; keep a leading unused slot so indices
        // stay identical to OCCT.
        boxes.insert(0, BndBox::new());
        self.my_boxes = Some(boxes);

        let upper = self.my_boxes.as_ref().unwrap().len() - 1;
        for a_box_index in 1..=upper {
            let a_box = &self.my_boxes.as_ref().unwrap()[a_box_index];
            if !a_box.is_void() {
                self.my_enclosing_box.add_box(a_box);
            }
        }

        self.my_resolution = get_bnd_voxel_grid_resolution(upper as i32);

        if self.my_enclosing_box.is_void() {
            return;
        }

        self.calculate_coefficients();
        self.reset_voxel_grid(upper as i32);
        self.sort_boxes();
    }

    /// OCCT Initialize(theEnclosingBox, theSetOfBoxes) — cxx L352-368.
    pub fn initialize_with_enclosing(&mut self, the_enclosing_box: BndBox, the_set_of_boxes: Vec<BndBox>) {
        let mut boxes = the_set_of_boxes;
        boxes.insert(0, BndBox::new());
        self.my_boxes = Some(boxes);
        self.my_enclosing_box = the_enclosing_box;
        let upper = self.my_boxes.as_ref().unwrap().len() - 1;
        self.my_resolution = get_bnd_voxel_grid_resolution(upper as i32);

        if self.my_enclosing_box.is_void() {
            return;
        }

        self.calculate_coefficients();
        self.reset_voxel_grid(upper as i32);
        self.sort_boxes();
    }

    /// OCCT Initialize(theEnclosingBox, theNbBoxes) — cxx L372-390: boxes
    /// are added later through [`Self::add`].
    pub fn initialize_incremental(&mut self, the_enclosing_box: BndBox, the_nb_boxes: i32) {
        if the_nb_boxes <= 0 {
            panic!("Unexpected: theNbBoxes <= 0");
        }
        // myBoxes = new HArray1(1, theNbBoxes) Init(emptyBox).
        let mut boxes = vec![BndBox::new(); (the_nb_boxes + 1) as usize];
        boxes[0] = BndBox::new();
        self.my_boxes = Some(boxes);
        self.my_enclosing_box = the_enclosing_box;
        self.my_resolution = get_bnd_voxel_grid_resolution(the_nb_boxes);

        if self.my_enclosing_box.is_void() {
            return;
        }

        self.calculate_coefficients();
        self.reset_voxel_grid(the_nb_boxes);
        // OCCT stops here: sortBoxes() runs inside Add() through addBox.
    }

    /// OCCT Add(theBox, theIndex) — cxx L394-406.
    pub fn add(&mut self, the_box: &BndBox, the_index: i32) {
        let boxes = self.my_boxes.as_mut().expect("Initialize not called");
        if !boxes[the_index as usize].is_void() {
            panic!(" This box is already defined !");
        }
        if the_box.is_void() {
            return;
        }

        boxes[the_index as usize] = the_box.clone();

        self.add_box(the_box, the_index);
    }

    /// OCCT Compare(theBox) — cxx L410-504: returns the list of box indices
    /// intersecting (or inside) theBox.  The result lives in the sorter
    /// (`myLastResult`), matching the OCCT reference semantics.
    pub fn compare(&mut self, the_box: &BndBox) -> &Vec<i32> {
        self.my_last_result.clear();

        if the_box.is_void() || self.my_enclosing_box.is_out_box(the_box) {
            return &self.my_last_result;
        }

        // Processing the large boxes (cxx L420-429).
        let large = self.my_large_boxes.clone();
        for a_box_index in &large {
            let a_box = &self.my_boxes.as_ref().unwrap()[*a_box_index as usize];
            if !a_box.is_out_box(the_box) {
                self.my_last_result.push(*a_box_index);
            }
        }

        // Obtaining the box voxel coordinates (cxx L432-433).
        let [a_min_voxel_x, a_min_voxel_y, a_min_voxel_z, a_max_voxel_x, a_max_voxel_y, a_max_voxel_z] =
            self.get_bounding_voxels(the_box);

        // aResultIndices bit masks (cxx L442-446).
        let upper = self.my_boxes.as_ref().unwrap().len() - 1;
        let mut a_result_indices = vec![0u8; upper + 1];
        const AN_OCCUPIED_X: u8 = 0b01;
        const AN_OCCUPIED_Y: u8 = 0b10;
        const AN_OCCUPIED_XY: u8 = 0b11;

        let grid = self.my_voxel_grid.as_ref().unwrap();

        // Checking the voxels along X-axis (cxx L448-461).
        for a_voxel_x in a_min_voxel_x..=a_max_voxel_x {
            let Some(a_box_indices) = grid.get_slice_x(a_voxel_x) else {
                continue;
            };
            for &a_box_index in a_box_indices {
                a_result_indices[a_box_index as usize] |= AN_OCCUPIED_X;
            }
        }

        // Checking the voxels along Y-axis (cxx L463-476).
        for a_voxel_y in a_min_voxel_y..=a_max_voxel_y {
            let Some(a_box_indices) = grid.get_slice_y(a_voxel_y) else {
                continue;
            };
            for &a_box_index in a_box_indices {
                a_result_indices[a_box_index as usize] |= AN_OCCUPIED_Y;
            }
        }

        // Checking the voxels along Z-axis (cxx L478-501).
        for a_voxel_z in a_min_voxel_z..=a_max_voxel_z {
            let Some(a_box_indices) = grid.get_slice_z(a_voxel_z) else {
                continue;
            };
            for &a_box_index in a_box_indices {
                if a_result_indices[a_box_index as usize] == AN_OCCUPIED_XY {
                    a_result_indices[a_box_index as usize] = 0;
                    let a_box = &self.my_boxes.as_ref().unwrap()[a_box_index as usize];
                    if !a_box.is_out_box(the_box) {
                        self.my_last_result.push(a_box_index);
                    }
                }
            }
        }

        &self.my_last_result
    }

    /// OCCT calculateCoefficients() — cxx L525-532.
    fn calculate_coefficients(&mut self) {
        let Some((a_xmin, a_ymin, a_zmin, a_xmax, a_ymax, a_zmax)) = self.my_enclosing_box.get()
        else {
            return;
        };
        let r = self.my_resolution as f64;
        self.my_coeff_x = if a_xmax - a_xmin == 0.0 { 0.0 } else { r / (a_xmax - a_xmin) };
        self.my_coeff_y = if a_ymax - a_ymin == 0.0 { 0.0 } else { r / (a_ymax - a_ymin) };
        self.my_coeff_z = if a_zmax - a_zmin == 0.0 { 0.0 } else { r / (a_zmax - a_zmin) };
    }

    /// OCCT resetVoxelGrid() — cxx L536-542.
    fn reset_voxel_grid(&mut self, box_count: i32) {
        self.my_voxel_grid = Some(VoxelGrid::new(self.my_resolution as usize, box_count));
        self.my_large_boxes.clear();
    }

    /// OCCT sortBoxes() — cxx L546-552.
    fn sort_boxes(&mut self) {
        let upper = self.my_boxes.as_ref().unwrap().len() - 1;
        for a_box_index in 1..=upper {
            let a_box = self.my_boxes.as_ref().unwrap()[a_box_index].clone();
            self.add_box(&a_box, a_box_index as i32);
        }
    }

    /// OCCT getBoundingVoxels(theBox) — cxx L556-596: voxel indices of the
    /// box min/max corners, clamped into the grid.
    ///
    /// The C++ expression is `std::clamp(static_cast<int>(v) - 1, 0,
    /// myResolution - 1)`. The `static_cast<int>` of an out-of-range double
    /// saturates to INT_MIN/INT_MAX, and the `- 1` / `+ 1` are plain C++ int
    /// operations (they wrap; the clamp then discards the wrapped value). A
    /// box that reaches the grid from an unbounded segment (the sampled
    /// polygon of an infinite line, e.g. the ±Infinite y of an axis-aligned
    /// segment) therefore produces exactly this saturation, so the Rust
    /// arithmetic must be the wrapping one to keep the same result instead of
    /// panicking on overflow.
    fn get_bounding_voxels(&self, the_box: &BndBox) -> [i32; 6] {
        // Start point of the voxel grid (cxx L559).  Note: `unwrap` mirrors
        // the OCCT CornerMin() on a non-void enclosing box.
        let a_grid_start = self.my_enclosing_box.corner_min().unwrap_or(DVec3::ZERO);

        let Some((a_x_min, a_y_min, a_z_min, a_x_max, a_y_max, a_z_max)) = the_box.get() else {
            return [0; 6];
        };

        let clamp_res = |v: i32| -> i32 { v.clamp(0, self.my_resolution - 1) };
        // The C++ shape is `clamp((int)v - 1, 0, res - 1)`: the `- 1`/`+ 1`
        // sits INSIDE the clamp argument and is a plain C++ int operation.
        let idx_min = |v: f64| clamp_res((v as i32).wrapping_sub(1));
        let idx_max = |v: f64| clamp_res((v as i32).wrapping_add(1));
        let a_x_min_index = idx_min((a_x_min - a_grid_start.x) * self.my_coeff_x);
        let a_y_min_index = idx_min((a_y_min - a_grid_start.y) * self.my_coeff_y);
        let a_z_min_index = idx_min((a_z_min - a_grid_start.z) * self.my_coeff_z);
        let a_x_max_index = idx_max((a_x_max - a_grid_start.x) * self.my_coeff_x);
        let a_y_max_index = idx_max((a_y_max - a_grid_start.y) * self.my_coeff_y);
        let a_z_max_index = idx_max((a_z_max - a_grid_start.z) * self.my_coeff_z);


        [a_x_min_index, a_y_min_index, a_z_min_index, a_x_max_index, a_y_max_index, a_z_max_index]
    }

    /// OCCT addBox(theBox, theIndex) — cxx L600-624: boxes spanning a large
    /// portion of the grid go to `myLargeBoxes` instead of the grid.
    fn add_box(&mut self, the_box: &BndBox, the_index: i32) {
        if the_box.is_void() {
            return;
        }

        let voxel_box = self.get_bounding_voxels(the_box);

        let a_box_min_side = (voxel_box[3] - voxel_box[0])
            .min(voxel_box[4] - voxel_box[1])
            .min(voxel_box[5] - voxel_box[2]);

        if a_box_min_side * 4 > self.my_resolution {
            self.my_large_boxes.push(the_index);
        } else {
            self.my_voxel_grid
                .as_mut()
                .unwrap()
                .add_box(the_index, &voxel_box);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn box_at(a: [f64; 3], b: [f64; 3]) -> BndBox {
        BndBox::from_corners(a[0], a[1], a[2], b[0], b[1], b[2])
    }

    /// OCCT Bnd_BoundSortBox behavior anchor: Compare returns exactly the
    /// boxes overlapping the query box.
    #[test]
    fn compare_finds_overlapping_boxes() {
        let mut bsb = BoundSortBox::new();
        let boxes = vec![
            box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_at([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]),
            box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]),
            box_at([10.0, 10.0, 10.0], [11.0, 11.0, 11.0]),
        ];
        bsb.initialize(boxes);
        let query = box_at([0.8, 0.8, 0.8], [1.2, 1.2, 1.2]);
        let res = bsb.compare(&query).clone();
        assert_eq!(res.len(), 2);
        assert!(res.contains(&1));
        assert!(res.contains(&3));
        // Far-away query intersects nothing.
        let far = box_at([50.0, 50.0, 50.0], [51.0, 51.0, 51.0]);
        assert_eq!(bsb.compare(&far).len(), 0);
    }
}
