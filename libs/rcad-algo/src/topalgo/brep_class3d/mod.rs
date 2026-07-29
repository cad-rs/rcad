// OCCT BRepClass3d — Solid classification algorithms.
//
// OCCT ref: TKTopAlgo/BRepClass3d/
//
// Provides point-in-solid classification via BRepClass3d_SolidClassifier.
// The classifier uses ray casting with BVH-accelerated face intersection.

pub mod solid_classifier;
pub mod s_classifier;
pub mod solid_explorer;
pub mod bnd_box_tree;
pub mod intersector3d;
pub mod passive_classifier;
