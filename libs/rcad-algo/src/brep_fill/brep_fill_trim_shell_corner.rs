//! OCCT BRepFill_TrimShellCorner (TKBool/BRepFill —
//! BRepFill_TrimShellCorner.hxx L25-100 + BRepFill_TrimShellCorner.cxx,
//! 2855 lines) — trims sets of faces in the corner to make proper parts of
//! pipe.
//!
//! D6 / GAP carrier: the real body runs the boolean machinery over the
//! corner faces (BOPDS_PDS — BOPAlgo_PaveFiller F-F intersections,
//! MakeFacesSec / MakeFacesNonSec / ChooseSection).  BOPDS is the boolean
//! DS owned by the bop/ module (outside this dispatch's owned file set);
//! per the dispatch protocol the class is carried at the exact OCCT API
//! surface and Perform() preserves the OCCT failure path.  Consumer:
//! BRepFill_Sweep::PerformCorner (part B, cxx L3714-3767).

use std::collections::HashMap;

use glam::DVec3;

use rcad_kernel::topo::topods::Shape;

use crate::brep_fill::brep_fill_pipe_shell_b::{BRepFillTransitionStyle, ShapeHArray2};
use crate::brep_fill::generator::ShapeKey;

/// OCCT BRepFill_TrimShellCorner (hxx L27-100).
pub struct BRepFillTrimShellCorner {
    /// OCCT myTransition.
    my_transition: BRepFillTransitionStyle,
    /// OCCT myAxeOfBisPlane (gp_Ax2: location + normal + x direction).
    my_axe_of_bis_plane: (DVec3, DVec3, DVec3),
    /// OCCT myIntPointCrossDir.
    my_int_point_cross_dir: DVec3,
    /// OCCT myShape1.
    my_shape1: Shape,
    /// OCCT myShape2.
    my_shape2: Shape,
    /// OCCT myBounds.
    my_bounds: Option<ShapeHArray2>,
    /// OCCT myUEdges.
    my_u_edges: Option<ShapeHArray2>,
    /// OCCT myVEdges (NCollection_HArray1<TopoDS_Shape>).
    my_v_edges: Option<Vec<Shape>>,
    /// OCCT myFaces.
    my_faces: Option<ShapeHArray2>,
    /// OCCT myDone.
    my_done: bool,
    /// OCCT myHasSection.
    my_has_section: bool,
    /// OCCT myHistMap
    /// (NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>).
    my_hist_map: HashMap<ShapeKey, Vec<Shape>>,
}

impl BRepFillTrimShellCorner {
    /// OCCT BRepFill_TrimShellCorner(theFaces, theTransition,
    /// theAxeOfBisPlane, theIntPointCrossDir) (cxx L44-64).
    pub fn new(
        the_faces: &ShapeHArray2,
        the_transition: BRepFillTransitionStyle,
        the_axe_of_bis_plane: (DVec3, DVec3, DVec3),
        the_int_point_cross_dir: DVec3,
    ) -> Self {
        // OCCT: myFaces(theFaces->Array2()), the rest defaulted.
        BRepFillTrimShellCorner {
            my_transition: the_transition,
            my_axe_of_bis_plane: the_axe_of_bis_plane,
            my_int_point_cross_dir: the_int_point_cross_dir,
            my_shape1: Shape::null(),
            my_shape2: Shape::null(),
            my_bounds: None,
            my_u_edges: None,
            my_v_edges: None,
            my_faces: Some(the_faces.clone()),
            my_done: false,
            my_has_section: false,
            my_hist_map: HashMap::new(),
        }
    }

    /// OCCT AddBounds (cxx L70-75).
    pub fn add_bounds(&mut self, bounds: &ShapeHArray2) {
        self.my_bounds = Some(bounds.clone());
    }

    /// OCCT AddUEdges (cxx L81-86).
    pub fn add_u_edges(&mut self, the_u_edges: &ShapeHArray2) {
        self.my_u_edges = Some(the_u_edges.clone());
    }

    /// OCCT AddVEdges (cxx L92-106) — the column `the_index` of the 2d array
    /// becomes the 1d array.
    pub fn add_v_edges(&mut self, the_v_edges: &ShapeHArray2, the_index: i32) {
        let col = (the_index - 1) as usize;
        self.my_v_edges = Some(
            the_v_edges
                .iter()
                .map(|row| row[col].clone())
                .collect::<Vec<Shape>>(),
        );
    }

    /// OCCT Perform (cxx L112-260) — the boolean corner trim
    /// (BOPDS_PDS-driven; see the file-header GAP note).
    pub fn perform(&mut self) {
        panic!(
            "GAP: BRepFill_TrimShellCorner::Perform (TKBool/BRepFill, \
             BRepFill_TrimShellCorner.cxx L112-260; BOPDS_PDS boolean \
             machinery owned by the bop/ module) is not translated — \
             see file header (D6)"
        );
    }

    /// OCCT IsDone (cxx L266-269).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT HasSection (cxx L275-278).
    pub fn has_section(&self) -> bool {
        self.my_has_section
    }

    /// OCCT Modified (cxx L284-296) — the history lookup.
    pub fn modified(&self, s: &Shape) -> Vec<Shape> {
        self.my_hist_map
            .get(&ShapeKey(s.ptr_id()))
            .cloned()
            .unwrap_or_default()
    }
}
