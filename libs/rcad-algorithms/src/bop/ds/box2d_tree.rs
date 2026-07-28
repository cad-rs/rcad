/// BOPTools_Box2dTree / BOPTools_BoxSelector / Bnd_Box2d.
use glam::DVec2;

#[derive(Debug, Clone, Copy)]
pub struct BndBox2d {
    pub min: DVec2,
    pub max: DVec2,
}

impl BndBox2d {
    pub fn new(min: DVec2, max: DVec2) -> Self {
        BndBox2d { min, max }
    }
    pub fn is_out(&self, other: &BndBox2d) -> bool {
        self.max.x < other.min.x
            || self.min.x > other.max.x
            || self.max.y < other.min.y
            || self.min.y > other.max.y
    }
    pub fn is_out_pt(&self, pt: DVec2) -> bool {
        pt.x < self.min.x || pt.x > self.max.x || pt.y < self.min.y || pt.y > self.max.y
    }
}

/// OCCT BOPTools_Box2dTree (BVH over 2D bboxes, linear-scan fallback).
pub struct Box2dTree {
    objects: Vec<(i32, BndBox2d)>,
}

impl Box2dTree {
    pub fn new() -> Self {
        Box2dTree {
            objects: Vec::new(),
        }
    }
    pub fn add(&mut self, id: i32, bbox: BndBox2d) {
        self.objects.push((id, bbox));
    }
    pub fn build(&mut self) {}
    pub fn len(&self) -> usize {
        self.objects.len()
    }
    pub fn element_box(&self, idx: usize) -> &BndBox2d {
        &self.objects[idx].1
    }
    pub fn element_id(&self, idx: usize) -> i32 {
        self.objects[idx].0
    }
}

/// OCCT BOPTools_BoxSelector<2>.
pub struct Box2dTreeSelector {
    my_box: BndBox2d,
    my_indices: Vec<i32>,
}

impl Box2dTreeSelector {
    pub fn new() -> Self {
        Box2dTreeSelector {
            my_box: BndBox2d::new(DVec2::ZERO, DVec2::ZERO),
            my_indices: Vec::new(),
        }
    }
    pub fn set_box(&mut self, bbox: BndBox2d) {
        self.my_box = bbox;
    }
    pub fn clear(&mut self) {
        self.my_indices.clear();
    }
    pub fn indices(&self) -> &[i32] {
        &self.my_indices
    }
    pub fn perform(&mut self, tree: &Box2dTree) {
        self.my_indices.clear();
        for (id, bbox) in &tree.objects {
            if !self.my_box.is_out(bbox) {
                self.my_indices.push(*id);
            }
        }
    }
}
