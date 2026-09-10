//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_ShapeTolerance` (`ShapeAnalysis_ShapeTolerance.hxx`
//! L17-72 + `.cxx` L1-332).
//!
//! Queries and accumulates the tolerances of the sub-shapes of a shape
//! (minimum / average / maximum).
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — `BRep_Tool::Tolerance` reads the TShape data
//!    through `rcad_kernel::BRep` (the edge.rs bridge #1).
//! 2. `TopExp_Explorer` -> `brep_tool::topexp_explorer` (the depth-first
//!    walk collecting every occurrence of the type).
//! 3. `occ::handle<NCollection_HSequence<TopoDS_Shape>>` -> `Vec<Shape>`.
//! 4. `TopAbs_ShapeEnum` -> `rcad_kernel::topods::ShapeType`.

use std::collections::HashSet;

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, ShapeType, TShape};

use crate::shhealing::shape_build::brep_tool::topexp_explorer;

/// OCCT BRep_Tool::Tolerance(S) — the vertex/edge/face tolerance.
fn brep_tool_tolerance(the_s: &Shape) -> f64 {
    match the_s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT static AddTol (cxx L32-65): accumulates one tolerance into the
/// running (nbt, cmin, cmoy, cmax) aggregates.
fn add_tol(tol: f64, nbt: &mut i32, cmin: &mut f64, cmoy: &mut f64, cmax: &mut f64) {
    *nbt += 1;
    if *nbt == 1 {
        *cmin = tol;
        *cmoy = tol;
        *cmax = tol;
    } else {
        if *cmin > tol {
            *cmin = tol;
        }
        if *cmax < tol {
            *cmax = tol;
        }
        //    cmoy += tol;
        //  Calcul en moyenne geometrique  entre 1 et 1e-7
        let mult = 1;
        // #76 rln 11.03.99 S4135: compute average without weights according to tolerances
        /*    if      (tol < 1.e-07) mult = 10000;
            else if (tol < 1.e-06) mult =  3000;
            else if (tol < 1.e-05) mult =  1000;
            else if (tol < 1.e-04) mult =   300;
            else if (tol < 1.e-03) mult =   100;
            else if (tol < 1.e-02) mult =    30;
            else if (tol < 1.e-01) mult =    10;
            else if (tol < 1.    ) mult =     3;
        */
        *nbt += mult - 1;
        *cmoy += tol * mult as f64;
    }
}

/// OCCT ShapeAnalysis_ShapeTolerance (hxx L24-72).
pub struct ShapeAnalysisShapeTolerance {
    /// OCCT `myNbTol`.
    my_nb_tol: i32,
    /// OCCT `double myTols[3]` — [0]: min, [1]: average sum, [2]: max.
    my_tols: [f64; 3],
}

impl Default for ShapeAnalysisShapeTolerance {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisShapeTolerance {
    /// OCCT ShapeAnalysis_ShapeTolerance() (cxx L27-30).
    pub fn new() -> Self {
        ShapeAnalysisShapeTolerance {
            my_nb_tol: 0,
            my_tols: [0.0; 3],
        }
    }

    /// OCCT Tolerance(shape, mode, type = TopAbs_SHAPE) (cxx L69-76).
    pub fn tolerance(
        &mut self,
        brep: &mut BRep,
        shape: &Shape,
        mode: i32,
        the_type: ShapeType,
    ) -> f64 {
        self.init_tolerance();
        self.add_tolerance(brep, shape, the_type);
        self.global_tolerance(mode)
    }

    /// OCCT OverTolerance(shape, value, type) (cxx L80-93).
    pub fn over_tolerance(
        &self,
        brep: &mut BRep,
        shape: &Shape,
        value: f64,
        the_type: ShapeType,
    ) -> Vec<Shape> {
        if value >= 0. {
            self.in_tolerance(brep, shape, value, 0., the_type)
        } else {
            self.in_tolerance(brep, shape, 0., value, the_type)
        }
    }

    /// OCCT InTolerance(shape, valmin, valmax, type) (cxx L97-226).
    pub fn in_tolerance(
        &self,
        brep: &mut BRep,
        shape: &Shape,
        valmin: f64,
        valmax: f64,
        the_type: ShapeType,
    ) -> Vec<Shape> {
        let mut tol;
        // pas de liminite max
        let over = valmax < valmin;
        let mut sl: Vec<Shape> = Vec::new();

        // Iteration sur les Faces

        if the_type == ShapeType::Face || the_type == ShapeType::Shape {
            for cur in topexp_explorer(brep, shape, ShapeType::Face) {
                tol = brep_tool_tolerance(&cur);
                if tol >= valmin && (over || tol <= valmax) {
                    sl.push(cur);
                }
            }
        }

        // Iteration sur les Edges

        if the_type == ShapeType::Edge || the_type == ShapeType::Shape {
            for cur in topexp_explorer(brep, shape, ShapeType::Edge) {
                tol = brep_tool_tolerance(&cur);
                if tol >= valmin && (over || tol <= valmax) {
                    sl.push(cur);
                }
            }
        }

        // Iteration sur les Vertex

        if the_type == ShapeType::Vertex || the_type == ShapeType::Shape {
            for cur in topexp_explorer(brep, shape, ShapeType::Vertex) {
                tol = brep_tool_tolerance(&cur);
                // OCCT L149: the vertex branch tests `tol >= valmax` (the
                // OCCT source literally reads `>=` where the face/edge
                // branches read `<=`); kept as-is.
                if tol >= valmin && (over || tol >= valmax) {
                    sl.push(cur);
                }
            }
        }

        // Iteration combinee (cumul) SHELL+FACE+EDGE+VERTEX, on retourne SHELL+FACE

        if the_type == ShapeType::Shell {
            //  Exploration des shells
            // OCCT NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>
            // mapface — keyed by the IsSame identity (TShape pointer +
            // location).
            let mut mapface: HashSet<(u64, u32)> = HashSet::new();
            for ash in topexp_explorer(brep, shape, ShapeType::Shell) {
                let mut iashell = false;
                for face in topexp_explorer(brep, &ash, ShapeType::Face) {
                    mapface
                        .insert((std::sync::Arc::as_ptr(&face.data) as u64, face.location));
                    let fc = self.in_tolerance(brep, &face, valmin, valmax, the_type);
                    if !fc.is_empty() {
                        sl.extend(fc);
                        iashell = true;
                    }
                }
                if iashell {
                    sl.push(ash);
                }
            }

            //  Les faces (libres ou sous shell)
            for cur in topexp_explorer(brep, shape, ShapeType::Face) {
                let mut iaface = false;
                if mapface
                    .contains(&(std::sync::Arc::as_ptr(&cur.data) as u64, cur.location))
                {
                    continue;
                }
                tol = brep_tool_tolerance(&cur);
                if tol >= valmin && (over || tol <= valmax) {
                    iaface = true;
                } else {
                    // les edges contenues ?
                    let fl = self.in_tolerance(brep, &cur, valmin, valmax, ShapeType::Edge);
                    if !fl.is_empty() {
                        iaface = true;
                    } else {
                        let fl = self.in_tolerance(brep, &cur, valmin, valmax, ShapeType::Vertex);
                        if !fl.is_empty() {
                            iaface = true;
                        }
                    }
                }
                if iaface {
                    sl.push(cur);
                }
            }
        }

        sl
    }

    /// OCCT InitTolerance() (cxx L230-234).
    pub fn init_tolerance(&mut self) {
        self.my_nb_tol = 0;
        self.my_tols[1] = 0.0;
    }

    /// OCCT AddTolerance(shape, type) (cxx L238-300).
    pub fn add_tolerance(&mut self, brep: &mut BRep, shape: &Shape, the_type: ShapeType) {
        let mut nbt = 0;
        let mut tol;
        let mut cmin = 0.;
        let mut cmoy = 0.;
        let mut cmax = 0.;

        // Iteration sur les Faces

        if the_type == ShapeType::Face || the_type == ShapeType::Shape {
            for cur in topexp_explorer(brep, shape, ShapeType::Face) {
                tol = brep_tool_tolerance(&cur);
                add_tol(tol, &mut nbt, &mut cmin, &mut cmoy, &mut cmax);
            }
        }

        // Iteration sur les Edges

        if the_type == ShapeType::Edge || the_type == ShapeType::Shape {
            for cur in topexp_explorer(brep, shape, ShapeType::Edge) {
                tol = brep_tool_tolerance(&cur);
                add_tol(tol, &mut nbt, &mut cmin, &mut cmoy, &mut cmax);
            }
        }

        // Iteration sur les Vertices

        if the_type == ShapeType::Vertex || the_type == ShapeType::Shape {
            for cur in topexp_explorer(brep, shape, ShapeType::Vertex) {
                tol = brep_tool_tolerance(&cur);
                add_tol(tol, &mut nbt, &mut cmin, &mut cmoy, &mut cmax);
            }
        }

        //  Resultat : attention en mode cumul
        if nbt == 0 {
            return;
        }
        if self.my_nb_tol == 0 || self.my_tols[0] > cmin {
            self.my_tols[0] = cmin;
        }
        if self.my_nb_tol == 0 || self.my_tols[2] < cmax {
            self.my_tols[2] = cmax;
        }
        self.my_tols[1] += cmoy;
        self.my_nb_tol += nbt;
    }

    /// OCCT GlobalTolerance(mode) (cxx L304-332).
    pub fn global_tolerance(&self, mode: i32) -> f64 {
        // szv#4:S4163:12Mar99 optimized
        let mut result = 0.;
        // OCCT L308: `if (myNbTol != 0.)` — the int member compared against
        // the double literal; kept as the zero comparison.
        if self.my_nb_tol != 0 {
            if mode < 0 {
                result = self.my_tols[0];
            } else if mode == 0 {
                if self.my_tols[0] == self.my_tols[2] {
                    result = self.my_tols[0];
                } else {
                    result = self.my_tols[1] / self.my_nb_tol as f64;
                }
            } else {
                result = self.my_tols[2];
            }
        }

        result
    }
}
