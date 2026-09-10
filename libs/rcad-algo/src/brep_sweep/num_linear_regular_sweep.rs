//! OCCT BRepSweep_NumLinearRegularSweep (TKPrim/BRepSweep) — the generic
//! sweep engine building swept primitives from a generating shape and a
//! directing line.
//!
//! Sources:
//! - BRepSweep_NumLinearRegularSweep.hxx L67-263
//! - BRepSweep_NumLinearRegularSweep.cxx L32-621
//!
//! Architecture difference (Rust has no inheritance): the C++
//! `BRepSweep_NumLinearRegularSweep` base class with its 16 pure-virtual
//! slots maps to the [`NumLinearRegularSweepCore`] state struct plus the
//! [`NumLinearRegularSweepSlots`] trait; the engine (the non-virtual Shape /
//! FirstShape / ... methods) are trait default methods — the OCCT call form
//! `sweep.Shape(g, d)` maps to `sl.shape(g, d)` on the concrete composition
//! root.

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};

use super::brep_sweep_iterator::BRepSweepIterator;
use super::brep_sweep_tool::BRepSweepTool;
use super::sweep_num_shape_iterator::SweepNumShapeIterator;
use super::sweep_num_shape_tool::SweepNumShapeTool;
use super::sweep_num_shape::SweepNumShape;
use super::tool_rehost::{brep_tool_is_closed_shape, set_closed_flag};
use super::BRepSweepBuilder;

/// OCCT `NCollection_Array2<T>` with the sweep's 1-based ranges
/// (1, myGenShapeTool.NbShapes(), 1, myDirShapeTool.NbShapes()).
pub struct Array2<T> {
    #[allow(dead_code)] // the OCCT Array2 row extent (kept for form)
    n_row: usize,
    n_col: usize,
    data: Vec<T>,
}

impl<T: Clone> Array2<T> {
    /// OCCT `NCollection_Array2<T>(1, nRow, 1, nCol)` filled with `value`.
    pub fn new_filled(n_row: usize, n_col: usize, value: T) -> Self {
        Array2 {
            n_row,
            n_col,
            data: vec![value; n_row * n_col],
        }
    }

    /// OCCT `A(i, j)` read (1-based).
    pub fn get(&self, i: i32, j: i32) -> &T {
        &self.data[((i - 1) as usize) * self.n_col + (j - 1) as usize]
    }

    /// OCCT `A(i, j)` write (1-based).
    pub fn set(&mut self, i: i32, j: i32, value: T) {
        self.data[((i - 1) as usize) * self.n_col + (j - 1) as usize] = value;
    }
    /// OCCT `A.UpperCol()` — the number of columns.
    pub fn upper_col(&self) -> i32 {
        self.n_col as i32
    }
}

/// The base-class members (BRepSweep_NumLinearRegularSweep.hxx L255-262).
pub struct NumLinearRegularSweepCore {
    /// OCCT: myBuilder.
    pub my_builder: BRepSweepBuilder,
    /// OCCT: myGenShape.
    pub my_gen_shape: Shape,
    /// OCCT: myDirWire.
    pub my_dir_wire: SweepNumShape,
    /// OCCT: myGenShapeTool.
    pub my_gen_shape_tool: BRepSweepTool,
    /// OCCT: myDirShapeTool.
    pub my_dir_shape_tool: SweepNumShapeTool,
    /// OCCT: myShapes (the shape grid).
    pub my_shapes: Array2<Shape>,
    /// OCCT: myBuiltShapes.
    pub my_built_shapes: Array2<bool>,
    /// OCCT: myUsedShapes.
    pub my_used_shapes: Array2<bool>,
}

impl NumLinearRegularSweepCore {
    /// OCCT BRepSweep_NumLinearRegularSweep::BRepSweep_NumLinearRegularSweep
    /// (cxx L32-54): the member tools and the three grids (the built/used
    /// grids initialized to false).
    pub fn new(a_builder: BRepSweepBuilder, a_gen_shape: &Shape, a_dir_wire: &SweepNumShape) -> Self {
        let my_gen_shape_tool = BRepSweepTool::new(a_gen_shape);
        let my_dir_shape_tool = SweepNumShapeTool::new(a_dir_wire);
        let n_gen = my_gen_shape_tool.nb_shapes();
        let n_dir = my_dir_shape_tool.nb_shapes();
        NumLinearRegularSweepCore {
            my_builder: a_builder,
            my_gen_shape: a_gen_shape.clone(),
            my_dir_wire: *a_dir_wire,
            my_gen_shape_tool,
            my_dir_shape_tool,
            my_shapes: Array2::new_filled(n_gen as usize, n_dir as usize, Shape::null()),
            my_built_shapes: Array2::new_filled(n_gen as usize, n_dir as usize, false),
            my_used_shapes: Array2::new_filled(n_gen as usize, n_dir as usize, false),
        }
    }
}

/// The pure-virtual slots of OCCT BRepSweep_NumLinearRegularSweep
/// (BRepSweep_NumLinearRegularSweep.hxx L74-208) plus the engine (the
/// non-virtual methods) as default methods.
pub trait NumLinearRegularSweepSlots {
    /// The base-class member access (the composition form of the C++
    /// inheritance); mutable.
    fn core(&mut self) -> &mut NumLinearRegularSweepCore;
    /// The base-class member access (const form).
    fn core_ref(&self) -> &NumLinearRegularSweepCore;

    /// Builds the vertex addressed by [aGenV,aDirV], with its geometric
    /// part, but without subcomponents (hxx L74-75).
    fn make_empty_vertex(&mut self, a_gen_v: &Shape, a_dir_v: &SweepNumShape) -> Shape;

    /// Builds the edge addressed by [aGenV,aDirE], with its geometric part,
    /// but without subcomponents (hxx L80-81).
    fn make_empty_directing_edge(&mut self, a_gen_v: &Shape, a_dir_e: &SweepNumShape) -> Shape;

    /// Builds the edge addressed by [aGenE,aDirV], with its geometric part,
    /// but without subcomponents (hxx L85-86).
    fn make_empty_generating_edge(&mut self, a_gen_e: &Shape, a_dir_v: &SweepNumShape) -> Shape;

    /// Sets the parameters of the new vertex on the new face (hxx L91-95).
    fn set_parameters(
        &mut self,
        a_new_face: &Shape,
        a_new_vertex: &mut Shape,
        a_gen_f: &Shape,
        a_gen_v: &Shape,
        a_dir_v: &SweepNumShape,
    );

    /// Sets the parameter of the new vertex on the new edge (hxx L100-104).
    fn set_directing_parameter(
        &mut self,
        a_new_edge: &Shape,
        a_new_vertex: &mut Shape,
        a_gen_v: &Shape,
        a_dir_e: &SweepNumShape,
        a_dir_v: &SweepNumShape,
    );

    /// Sets the parameter of the new vertex on the new edge (hxx L109-113).
    fn set_generating_parameter(
        &mut self,
        a_new_edge: &Shape,
        a_new_vertex: &mut Shape,
        a_gen_e: &Shape,
        a_gen_v: &Shape,
        a_dir_v: &SweepNumShape,
    );

    /// Builds the face addressed by [aGenS,aDirS], with its geometric part,
    /// but without subcomponents (hxx L120-121).
    fn make_empty_face(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> Shape;

    /// Sets the PCurve for a new edge on a new face (hxx L126-131).
    fn set_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        a_gen_f: &Shape,
        a_gen_e: &Shape,
        a_dir_v: &SweepNumShape,
        orien: Orientation,
    );

    /// Sets the PCurve for a new edge on a new face (hxx L136-141).
    fn set_generating_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        a_gen_e: &Shape,
        a_dir_e: &SweepNumShape,
        a_dir_v: &SweepNumShape,
        orien: Orientation,
    );

    /// Sets the PCurve for a new edge on a new face (hxx L146-151).
    fn set_directing_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        a_gen_e: &Shape,
        a_gen_v: &Shape,
        a_dir_e: &SweepNumShape,
        orien: Orientation,
    );

    /// Returns the Orientation of the shell in the solid generated by the
    /// face aGenS with the edge aDirS (hxx L157-158).
    fn direct_solid(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> Orientation;

    /// Returns true if aNewSubShape (addressed by aSubGenS and aDirS) must
    /// be added in aNewShape (addressed by aGenS and aDirS) (hxx L163-167).
    fn ggd_shape_is_to_add(
        &self,
        a_new_shape: &Shape,
        a_new_sub_shape: &Shape,
        a_gen_s: &Shape,
        a_sub_gen_s: &Shape,
        a_dir_s: &SweepNumShape,
    ) -> bool;

    /// Returns true if aNewSubShape (addressed by aGenS and aSubDirS) must
    /// be added in aNewShape (addressed by aGenS and aDirS) (hxx L172-176).
    fn gdd_shape_is_to_add(
        &self,
        a_new_shape: &Shape,
        a_new_sub_shape: &Shape,
        a_gen_s: &Shape,
        a_dir_s: &SweepNumShape,
        a_sub_dir_s: &SweepNumShape,
    ) -> bool;

    /// True when the generated face topology must be composed of independent
    /// closed wires (hxx L182-186).
    fn separated_wires(
        &self,
        a_new_shape: &Shape,
        a_new_sub_shape: &Shape,
        a_gen_s: &Shape,
        a_sub_gen_s: &Shape,
        a_dir_s: &SweepNumShape,
    ) -> bool;

    /// True when the generated Shell must be composed of independent closed
    /// Shells; returns a Compound of independent Shells (hxx L192; the base
    /// implementation cxx L543-549 wraps the shell in a compound).
    fn split_shell(&self, a_new_shape: &Shape) -> Shape {
        // OCCT BRepSweep_NumLinearRegularSweep::SplitShell (cxx L543-549).
        let core = self.core_ref();
        let comp = core.my_builder.my_builder.make_compound();
        core.my_builder.add(&comp, a_new_shape);
        comp
    }

    /// Propagates the continuity of every vertex between two edges of the
    /// generating wire aGenS on the generated edge and faces (hxx L197-198).
    fn set_continuity(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape);

    /// Returns true if aDirS and aGenS addresses a resulting Shape
    /// (hxx L204-205).
    fn has_shape(&self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> bool;

    /// Returns true if aGenS cannot be transformed (hxx L208).
    fn is_invariant(&self, a_gen_s: &Shape) -> bool;

    /// OCCT BRepSweep_NumLinearRegularSweep::Shape() (cxx L62-73) — returns
    /// the resulting Shape indexed by myDirWire and myGenShape.
    fn shape_full(&mut self) -> Shape {
        let core = self.core_ref();
        let a_gen = core.my_gen_shape.clone();
        let a_dir = core.my_dir_wire;
        if self.has_shape(&a_gen, &a_dir) {
            self.shape(&a_gen, &a_dir)
        } else {
            // OCCT: TopoDS_Shape bidon; return bidon.
            Shape::null()
        }
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::Shape(aGenS) (cxx L80-91) —
    /// returns the Shape generated with aGenS.
    fn shape_of_gen(&mut self, a_gen_s: &Shape) -> Shape {
        let core = self.core_ref();
        let a_dir = core.my_dir_wire;
        if core.my_gen_shape_tool.index(a_gen_s) != 0 && self.has_shape(a_gen_s, &a_dir) {
            self.shape(a_gen_s, &a_dir)
        } else {
            Shape::null()
        }
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::Shape(aGenS, aDirS)
    /// (cxx L98-460) — the engine: returns the resulting Shape indexed by
    /// the arguments.
    fn shape(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> Shape {
        let i_gen_s = self.core_ref().my_gen_shape_tool.index(a_gen_s);
        let i_dir_s = self.core_ref().my_dir_shape_tool.index(a_dir_s);
        eprintln!("[ENTRY] shape(i={},j={}) built={}", i_gen_s, i_dir_s, self.core_ref().my_built_shapes.get(i_gen_s, i_dir_s));
        if !*self.core_ref().my_built_shapes.get(i_gen_s, i_dir_s) {
            // OCCT L105-111: the newShape / bGenS / cGenS / subGenS /
            // subsubGenS / bDirS / subDirS / It / Kt / Lt / Or / Pr locals
            // (the OCCT default constructors are null shapes).
            let mut new_shape = Shape::null();
            let mut b_gen_s = Shape::null();
            if self.core_ref().my_dir_shape_tool.type_of(a_dir_s) == ShapeType::Vertex {
                // Ici on construit les "planchers" du Shape.
                let a_gen_s_type = self.core_ref().my_gen_shape_tool.type_of(a_gen_s);
                match a_gen_s_type {
                    ShapeType::Vertex => {
                        let made = self.make_empty_vertex(a_gen_s, a_dir_s);
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Edge => {
                        let made = self.make_empty_generating_edge(a_gen_s, a_dir_s);
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Wire => {
                        let made = self.core_ref().my_builder.my_builder.make_wire();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Face => {
                        let made = self.make_empty_face(a_gen_s, a_dir_s);
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Shell => {
                        let made = self.core_ref().my_builder.my_builder.make_shell();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Solid | ShapeType::CompSolid => {
                        // OCCT: throw Standard_NoSuchObject("Solids are not
                        // Processed").
                        panic!("Standard_NoSuchObject: Solids are not Processed");
                    }
                    ShapeType::Compound => {
                        let made = self.core_ref().my_builder.my_builder.make_compound();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    _ => {
                        // OCCT: throw Standard_NoSuchObject("Unknown Shape").
                        panic!("Standard_NoSuchObject: Unknown Shape");
                    }
                }
                b_gen_s = a_gen_s.clone();
                self.core_ref()
                    .my_gen_shape_tool
                    .set_orientation(&mut b_gen_s, Orientation::Forward);
                let mut it = BRepSweepIterator::default();
                it.init(&b_gen_s);
                while it.more() {
                    let sub_gen_s = it.value();
                    let or = it.orientation();
                    if self.has_shape(&sub_gen_s, a_dir_s) {
                        new_shape = self.shape(&sub_gen_s, a_dir_s);
                        let i_new_gen_s = self.core_ref().my_gen_shape_tool.index(&sub_gen_s);
                        let i_new_dir_s = i_dir_s;
                        let a_new_shape = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                        if self.ggd_shape_is_to_add(&a_new_shape, &new_shape, a_gen_s, &sub_gen_s, a_dir_s)
                        {
                            // Les "planchers" doivent etre construits par les
                            // fonctions de construcion geometrique
                            // identiquement au shape generateur.  On leur
                            // recolle juste une orientation pour etre bien
                            // sur.
                            let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                            self.core_ref()
                                .my_builder
                                .add_oriented(&target, &new_shape, or);
                            self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                            let sub_gen_s_type = self.core_ref().my_gen_shape_tool.type_of(&sub_gen_s);
                            if a_gen_s_type == ShapeType::Face {
                                if sub_gen_s_type == ShapeType::Vertex {
                                    let a_new_face = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                    self.set_parameters(&a_new_face, &mut new_shape, a_gen_s, &sub_gen_s, a_dir_s);
                                } else if sub_gen_s_type == ShapeType::Edge {
                                    let a_new_face = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                    self.set_pcurve(&a_new_face, &mut new_shape, a_gen_s, &sub_gen_s, a_dir_s, or);
                                } else if sub_gen_s_type == ShapeType::Wire {
                                    let mut c_gen_s = sub_gen_s.clone();
                                    self.core_ref()
                                        .my_gen_shape_tool
                                        .set_orientation(&mut c_gen_s, Orientation::Forward);
                                    let mut jt = BRepSweepIterator::default();
                                    jt.init(&c_gen_s);
                                    while jt.more() {
                                        let subsub_gen_s = jt.value();
                                        let pr = jt.orientation();
                                        if self.has_shape(&subsub_gen_s, a_dir_s) {
                                            let mut newsub_edge = self.shape(&subsub_gen_s, a_dir_s);
                                            let a_new_face = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                            self.set_pcurve(
                                                &a_new_face,
                                                &mut newsub_edge,
                                                a_gen_s,
                                                &subsub_gen_s,
                                                a_dir_s,
                                                pr,
                                            );
                                        }
                                        jt.next();
                                    }
                                }
                            } else if a_gen_s_type == ShapeType::Edge {
                                let a_new_edge = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.set_generating_parameter(&a_new_edge, &mut new_shape, &b_gen_s, &sub_gen_s, a_dir_s);
                            }
                        }
                    }
                    it.next();
                }
            } else if self.core_ref().my_dir_shape_tool.type_of(a_dir_s) == ShapeType::Edge {
                // Ici on construit les murs du Shape.
                let a_gen_s_type = self.core_ref().my_gen_shape_tool.type_of(a_gen_s);
                let mut new_wire = Shape::null();
                let mut new_shell = Shape::null();
                let mut wire_seq: Vec<Shape> = Vec::new();
                let mut sepwires = false;
                match a_gen_s_type {
                    ShapeType::Vertex => {
                        let made = self.make_empty_directing_edge(a_gen_s, a_dir_s);
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Edge => {
                        // On cree un wire intermediaire qui contient tous les
                        // edges du montant (face) du Shape pour le cas
                        // standard, et une sequence de wires pour les cas
                        // merdiques necessitant des wires independants.
                        new_wire = self.core_ref().my_builder.my_builder.make_wire();
                        let made = self.make_empty_face(a_gen_s, a_dir_s);
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Wire => {
                        let made = self.core_ref().my_builder.my_builder.make_shell();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Face => {
                        // On cree un shell intermediaire dans lequel on jette
                        // toutes les faces en direct, pour eviter les
                        // empilages compliques de shells et sous shells dans
                        // la structure du solide.
                        new_shell = self.core_ref().my_builder.my_builder.make_shell();
                        let made = self.core_ref().my_builder.my_builder.make_solid();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Shell => {
                        let made = self.core_ref().my_builder.my_builder.make_compsolid();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Solid | ShapeType::CompSolid => {
                        panic!("Standard_NoSuchObject: Solids are not Processed");
                    }
                    ShapeType::Compound => {
                        let made = self.core_ref().my_builder.my_builder.make_compound();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    _ => {
                        panic!("Standard_NoSuchObject: Unknown Shape");
                    }
                }
                b_gen_s = a_gen_s.clone();
                self.core_ref()
                    .my_gen_shape_tool
                    .set_orientation(&mut b_gen_s, Orientation::Forward);
                let mut it = BRepSweepIterator::default();
                it.init(&b_gen_s);
                while it.more() {
                    let sub_gen_s = it.value();
                    if self.has_shape(&sub_gen_s, a_dir_s) {
                        new_shape = self.shape(&sub_gen_s, a_dir_s);
                        let i_new_gen_s = self.core_ref().my_gen_shape_tool.index(&sub_gen_s);
                        let i_new_dir_s = i_dir_s;
                        let a_new_shape = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                        if self.ggd_shape_is_to_add(&a_new_shape, &new_shape, a_gen_s, &sub_gen_s, a_dir_s)
                        {
                            let sub_gen_s_type = self.core_ref().my_gen_shape_tool.type_of(&sub_gen_s);
                            if a_gen_s_type == ShapeType::Edge {
                                let or = it.orientation();
                                let a_new_shape2 = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                if self.separated_wires(&a_new_shape2, &new_shape, a_gen_s, &sub_gen_s, a_dir_s)
                                {
                                    sepwires = true;
                                    let wi = self.core_ref().my_builder.my_builder.make_wire();
                                    self.core_ref().my_builder.add_oriented(&wi, &new_shape, or);
                                    self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                                    // OCCT: wi.Closed(BRep_Tool::IsClosed(wi)).
                                    let is_closed = brep_tool_is_closed_shape(&wi);
                                    set_closed_flag(&wi, is_closed);
                                    wire_seq.push(wi);
                                } else {
                                    self.core_ref().my_builder.add_oriented(&new_wire, &new_shape, or);
                                    self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                                }
                                let a_new_face = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.set_directing_pcurve(
                                    &a_new_face,
                                    &mut new_shape,
                                    &b_gen_s,
                                    &sub_gen_s,
                                    a_dir_s,
                                    or,
                                );
                            } else if a_gen_s_type == ShapeType::Wire {
                                let or = it.orientation();
                                let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.core_ref().my_builder.add_oriented(&target, &new_shape, or);
                                self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                            } else if a_gen_s_type == ShapeType::Face {
                                let or = it.orientation();
                                if sub_gen_s_type == ShapeType::Wire {
                                    let mut lt = BRepSweepIterator::default();
                                    lt.init(&new_shape);
                                    while lt.more() {
                                        let lt_value = lt.value();
                                        // OCCT: TopAbs::Compose(Lt.Orientation(), Or).
                                        let comp = lt.orientation().compose(or);
                                        self.core_ref().my_builder.add_oriented(&new_shell, &lt_value, comp);
                                        lt.next();
                                    }
                                } else if sub_gen_s_type == ShapeType::Edge {
                                    self.core_ref().my_builder.add_oriented(&new_shell, &new_shape, or);
                                    self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                                }
                            } else if a_gen_s_type == ShapeType::Shell {
                                let or = Orientation::Forward;
                                let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.core_ref().my_builder.add_oriented(&target, &new_shape, or);
                                self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                            } else if a_gen_s_type == ShapeType::Compound {
                                let or = Orientation::Forward;
                                let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.core_ref().my_builder.add_oriented(&target, &new_shape, or);
                                self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                            } else {
                                let or = it.orientation();
                                let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.core_ref().my_builder.add_oriented(&target, &new_shape, or);
                                self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                            }
                        }
                    }
                    it.next();
                }
                // OCCT L325: bDirS = aDirS; the Kt loop over the directing
                // sub-shapes.
                let b_dir_s = *a_dir_s;
                let mut kt = SweepNumShapeIterator::default();
                kt.init(&b_dir_s);
                while kt.more() {
                    let sub_dir_s = *kt.value();
                    if self.has_shape(a_gen_s, &sub_dir_s) {
                        new_shape = self.shape(a_gen_s, &sub_dir_s);
                        let i_new_gen_s = i_gen_s;
                        let i_new_dir_s = self.core_ref().my_dir_shape_tool.index(&sub_dir_s);
                        let a_new_shape = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                        if self.gdd_shape_is_to_add(&a_new_shape, &new_shape, a_gen_s, a_dir_s, &sub_dir_s)
                        {
                            if a_gen_s_type == ShapeType::Edge {
                                // OCCT: Or = TopAbs::Reverse(Kt.Orientation()).
                                let or = top_abs_reverse(kt.orientation());
                                self.core_ref().my_builder.add_oriented(&new_wire, &new_shape, or);
                                self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                                let a_new_face = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.set_generating_pcurve(&a_new_face, &mut new_shape, a_gen_s, a_dir_s, &sub_dir_s, or);
                            } else if a_gen_s_type == ShapeType::Vertex {
                                let or = kt.orientation();
                                let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.core_ref().my_builder.add_oriented(&target, &new_shape, or);
                                self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                                let a_new_edge = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                                self.set_directing_parameter(&a_new_edge, &mut new_shape, a_gen_s, a_dir_s, &sub_dir_s);
                            } else if a_gen_s_type == ShapeType::Face {
                                let or = kt.orientation();
                                self.core_ref().my_builder.add_oriented(&new_shell, &new_shape, or);
                                self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                            }
                        }
                    }
                    kt.next();
                }
                if a_gen_s_type == ShapeType::Edge {
                    if sepwires {
                        for ij in 1..=wire_seq.len() {
                            let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                            let wi = wire_seq[ij - 1].clone();
                            self.core_ref().my_builder.add(&target, &wi);
                        }
                    } else {
                        // OCCT: newWire.Closed(BRep_Tool::IsClosed(newWire)).
                        let is_closed = brep_tool_is_closed_shape(&new_wire);
                        set_closed_flag(&new_wire, is_closed);
                        let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                        self.core_ref().my_builder.add(&target, &new_wire);
                    }
                    self.core().my_built_shapes.set(i_gen_s, i_dir_s, true);
                    self.set_continuity(a_gen_s, a_dir_s);
                }
                if a_gen_s_type == ShapeType::Wire {
                    self.set_continuity(a_gen_s, a_dir_s);
                }
                if a_gen_s_type == ShapeType::Face {
                    // OCCT: newShell.Closed(BRep_Tool::IsClosed(newShell)).
                    let is_closed = brep_tool_is_closed_shape(&new_shell);
                    set_closed_flag(&new_shell, is_closed);
                    let temp = self.split_shell(&new_shell);
                    let shell_ori = self.direct_solid(a_gen_s, a_dir_s);
                    let mut lt = BRepSweepIterator::default();
                    lt.init(&temp);
                    if lt.more() {
                        lt.next();
                    }
                    if lt.more() {
                        let mut lt = BRepSweepIterator::default();
                        lt.init(&temp);
                        while lt.more() {
                            let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                            let lt_value = lt.value();
                            self.core_ref().my_builder.add_oriented(&target, &lt_value, shell_ori);
                            lt.next();
                        }
                    } else {
                        let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                        self.core_ref().my_builder.add_oriented(&target, &new_shell, shell_ori);
                    }
                }
            } else if self.core_ref().my_dir_shape_tool.type_of(a_dir_s) == ShapeType::Wire {
                let a_gen_s_type = self.core_ref().my_gen_shape_tool.type_of(a_gen_s);
                match a_gen_s_type {
                    ShapeType::Vertex => {
                        let made = self.core_ref().my_builder.my_builder.make_wire();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Edge | ShapeType::Wire => {
                        let made = self.core_ref().my_builder.my_builder.make_shell();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Face | ShapeType::Shell => {
                        let made = self.core_ref().my_builder.my_builder.make_compsolid();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    ShapeType::Solid | ShapeType::CompSolid => {
                        panic!("Standard_NoSuchObject: Solids are not Processed");
                    }
                    ShapeType::Compound => {
                        let made = self.core_ref().my_builder.my_builder.make_compound();
                        self.core().my_shapes.set(i_gen_s, i_dir_s, made);
                    }
                    _ => {
                        panic!("Standard_NoSuchObject: Unknown Shape");
                    }
                }
                let b_dir_s = *a_dir_s;
                let mut kt = SweepNumShapeIterator::default();
                kt.init(&b_dir_s);
                while kt.more() {
                    let sub_dir_s = *kt.value();
                    if self.has_shape(a_gen_s, &sub_dir_s) {
                        let i_new_gen_s = i_gen_s;
                        let i_new_dir_s = self.core_ref().my_dir_shape_tool.index(&sub_dir_s);
                        let or = kt.orientation();
                        new_shape = self.shape(a_gen_s, &sub_dir_s);
                        let target = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
                        self.core_ref().my_builder.add_oriented(&target, &new_shape, or);
                        self.core().my_used_shapes.set(i_new_gen_s, i_new_dir_s, true);
                    }
                    kt.next();
                }
            }
            self.core().my_built_shapes.set(i_gen_s, i_dir_s, true);
        }
        // Change the "Closed" flag only for Wires and Shells.
        let result = self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone();
        if result.shape_type() == ShapeType::Wire || result.shape_type() == ShapeType::Shell {
            let is_closed = brep_tool_is_closed_shape(&result);
            set_closed_flag(&result, is_closed);
        }
        self.core_ref().my_shapes.get(i_gen_s, i_dir_s).clone()
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::FirstShape() (cxx L467-478).
    fn first_shape_full(&mut self) -> Shape {
        let mut result = Shape::null();
        if self.core_ref().my_dir_shape_tool.has_first_vertex() {
            let first = self.core_ref().my_dir_shape_tool.first_vertex();
            if self.has_shape(&self.core_ref().my_gen_shape, &first) {
                let the_gen = self.core_ref().my_gen_shape.clone();
                result = self.shape(&the_gen, &first);
            }
        }
        result
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::LastShape() (cxx L485-496).
    fn last_shape_full(&mut self) -> Shape {
        let mut result = Shape::null();
        if self.core_ref().my_dir_shape_tool.has_last_vertex() {
            let last = self.core_ref().my_dir_shape_tool.last_vertex();
            if self.has_shape(&self.core_ref().my_gen_shape, &last) {
                let the_gen = self.core_ref().my_gen_shape.clone();
                result = self.shape(&the_gen, &last);
            }
        }
        result
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::FirstShape(aGenS) (cxx L503-514).
    fn first_shape(&mut self, a_gen_s: &Shape) -> Shape {
        let mut result = Shape::null();
        if self.core_ref().my_dir_shape_tool.has_first_vertex() {
            let first = self.core_ref().my_dir_shape_tool.first_vertex();
            if self.has_shape(a_gen_s, &first) {
                result = self.shape(a_gen_s, &first);
            }
        }
        result
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::LastShape(aGenS) (cxx L521-532).
    fn last_shape(&mut self, a_gen_s: &Shape) -> Shape {
        let mut result = Shape::null();
        if self.core_ref().my_dir_shape_tool.has_last_vertex() {
            let last = self.core_ref().my_dir_shape_tool.last_vertex();
            if self.has_shape(a_gen_s, &last) {
                result = self.shape(a_gen_s, &last);
            }
        }
        result
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::Closed() (cxx L536-539).
    fn closed(&self) -> bool {
        self.core_ref().my_dir_wire.closed()
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::IsUsed(aGenS) (cxx L553-602).
    fn is_used(&self, a_gen_s: &Shape) -> bool {
        let core = self.core_ref();
        let i_gen_s = core.my_gen_shape_tool.index(a_gen_s);
        if i_gen_s == 0 {
            return false;
        }
        let mut is_built = false;
        let mut is_used = false;
        let mut j = 2;
        while j <= core.my_built_shapes.upper_col() {
            is_built = is_built || *core.my_built_shapes.get(i_gen_s, j);
            is_used = is_used || *core.my_used_shapes.get(i_gen_s, j);
            j += 1;
        }
        if is_used {
            if a_gen_s.shape_type() == ShapeType::Vertex && self.is_invariant(a_gen_s) {
                if *core.my_used_shapes.get(i_gen_s, 1) || !self.closed() {
                    return is_used;
                } else {
                    return false;
                }
            } else {
                return is_used;
            }
        }
        // OCCT: if (isBuilt) //&& !IsUsed
        if is_built {
            let dir_wire = core.my_dir_wire;
            if !self.has_shape(a_gen_s, &dir_wire) && !self.closed() {
                return true;
            } else if a_gen_s.shape_type() == ShapeType::Vertex && !self.closed() {
                if !*core.my_built_shapes.get(i_gen_s, 1) {
                    return true;
                }
            }
        }
        is_used
    }

    /// OCCT BRepSweep_NumLinearRegularSweep::GenIsUsed(theS) (cxx L606-621).
    fn gen_is_used(&self, the_s: &Shape) -> bool {
        let core = self.core_ref();
        let i_gen_s = core.my_gen_shape_tool.index(the_s);
        if i_gen_s == 0 {
            return false;
        }
        if i_gen_s == 1 {
            *core.my_built_shapes.get(i_gen_s, 1)
        } else {
            *core.my_built_shapes.get(i_gen_s, 1) && *core.my_used_shapes.get(i_gen_s, 1)
        }
    }
}

/// OCCT TopAbs::Reverse (TopAbs.hxx L80-91): FORWARD<->REVERSED,
/// INTERNAL/EXTERNAL unchanged.
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}
