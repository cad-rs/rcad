// OCCT HLRBRep_FaceData (TKHLR/HLRBRep/HLRBRep_FaceData.hxx L1-155
// + .cxx L1-70 + .lxx L1-286) — the per-face record of the HLRBRep_Data
// structure: the bit-packed flags, the wire block of the projected face and
// its geometry (HLRBRep_Surface).
//
// The OCCT lxx bodies are all inline; they land in this file next to the
// struct.  The EMaskFlags bit layout is preserved bit-exactly (the flags
// word is shared with the mask arithmetic of the Hider / Data).

use std::sync::Arc;

use rcad_kernel::topo::topods::{BRepTool, Orientation, Shape};

use super::surface::Surface;

use crate::hlr::algo::edges_block::EdgesBlock;
use crate::hlr::algo::wires_block::WiresBlock;

/// OCCT HLRBRep_FaceData (hxx L30-153).
//
// Clone: the OCCT class is copied by the default copy assignment in
// HLRBRep_Data::Write (`*fd = *f1`, cxx L584) — a copy of the flags, the
// Wires handle, the surface value and the scalars.
#[derive(Clone)]
pub struct FaceData<'a> {
    /// OCCT myFlags (hxx L146) — the EMaskFlags bit field (hxx L128-144),
    /// bit-exact: Orient = bits 0-3, Selected = 16, Back = 32, Side = 64,
    /// Closed = 128, Hiding = 256, Simple = 512, Cut = 1024, WithOutL =
    /// 2048, Plane = 4096, Cylinder = 8192, Cone = 16384, Sphere = 32768,
    /// Torus = 65536.
    my_flags: i32,
    /// OCCT myWires (hxx L147).
    my_wires: Option<Arc<WiresBlock>>,
    /// OCCT myGeometry (hxx L148).
    my_geometry: Surface<'a>,
    /// OCCT mySize (hxx L149).
    my_size: f64,
    /// OCCT myTolerance (hxx L150).
    my_tolerance: f32,
}

/// OCCT enum EMaskFlags (hxx L128-144) — bit-exact.
pub mod emask_flags {
    pub const EMASK_ORIENT: i32 = 15;
    pub const FMASK_SELECTED: i32 = 16;
    pub const FMASK_BACK: i32 = 32;
    pub const FMASK_SIDE: i32 = 64;
    pub const FMASK_CLOSED: i32 = 128;
    pub const FMASK_HIDING: i32 = 256;
    pub const FMASK_SIMPLE: i32 = 512;
    pub const FMASK_CUT: i32 = 1024;
    pub const FMASK_WITH_OUT_L: i32 = 2048;
    pub const FMASK_PLANE: i32 = 4096;
    pub const FMASK_CYLINDER: i32 = 8192;
    pub const FMASK_CONE: i32 = 16384;
    pub const FMASK_SPHERE: i32 = 32768;
    pub const FMASK_TORUS: i32 = 65536;

    pub use EMASK_ORIENT as EMASK_F_OCCULTATION;
}

impl<'a> FaceData<'a> {
    /// OCCT HLRBRep_FaceData::HLRBRep_FaceData() (cxx L25-30) — zero flags
    /// and zero size, then Selected(true).  myWires stays a null handle and
    /// myTolerance is left uninitialized in OCCT (the neutral 0 default
    /// stands in).
    pub fn new() -> Self {
        let mut fd = FaceData {
            my_flags: 0,
            my_wires: None,
            my_geometry: Surface::new(),
            my_size: 0.0,
            my_tolerance: 0.0,
        };
        fd.set_selected(true);
        fd
    }

    /// OCCT Set(FG, Or, Cl, NW) (cxx L34-44) — <Or> is the orientation of
    /// the face, <Cl> is true if the face belongs to a closed volume, <NW>
    /// is the number of wires (or block of edges) of the face.  The TopoDS
    /// face carries its owning BRep through the kernel-boundary `brep`
    /// reference (the HLRBRep_EdgeData::Set precedent).
    pub fn set(
        &mut self,
        brep: &'a rcad_kernel::BRep,
        fg: &Shape,
        or_: Orientation,
        cl: bool,
        nw: usize,
    ) {
        self.set_closed(cl);
        self.my_geometry.load(brep, fg); // Geometry().Surface(FG)
        // myTolerance = (float)(BRep_Tool::Tolerance(FG));
        self.my_tolerance = brep.tolerance(fg) as f32;
        self.set_orientation(or_);
        // Wires() = new HLRAlgo_WiresBlock(NW);
        self.my_wires = Some(Arc::new(WiresBlock::new(nw)));
    }

    /// OCCT SetWire(WI, NE) (cxx L48-51) — set <NE> the number of edges of
    /// the wire number <WI>.
    pub fn set_wire(&mut self, wi: usize, ne: usize) {
        // Wires()->Set(WI, new HLRAlgo_EdgesBlock(NE));
        let wb = self
            .my_wires
            .as_mut()
            .expect("HLRBRep_FaceData::SetWire: null Wires handle");
        Arc::get_mut(wb)
            .expect("HLRBRep_FaceData::SetWire: shared Wires handle")
            .set(wi, EdgesBlock::new(ne));
    }

    /// OCCT SetWEdge(WI, EWI, EI, Or, OutL, Inte, Dble, IsoL) (cxx L55-70) —
    /// set the edge number <EWI> of the wire <WI>.
    #[allow(clippy::too_many_arguments)]
    pub fn set_w_edge(
        &mut self,
        wi: usize,
        ewi: usize,
        ei: i32,
        or_: Orientation,
        out_l: bool,
        inte: bool,
        dble: bool,
        iso_l: bool,
    ) {
        // Wires()->Wire(WI) — the sole-handle mutation through Arc::get_mut
        // stands for the OCCT handle dereference.
        let wb = self
            .my_wires
            .as_mut()
            .expect("HLRBRep_FaceData::SetWEdge: null Wires handle");
        let w = Arc::get_mut(wb)
            .expect("HLRBRep_FaceData::SetWEdge: shared Wires handle")
            .wire(wi);
        w.set_edge(ewi, ei); // Edge(EWI, EI)
        w.set_orientation(ewi, or_); // Orientation(EWI, Or)
        w.set_out_line(ewi, out_l); // OutLine(EWI, OutL)
        w.set_internal(ewi, inte); // Internal(EWI, Inte)
        w.set_double(ewi, dble); // Double(EWI, Dble)
        w.set_iso_line(ewi, iso_l); // IsoLine(EWI, IsoL)
    }

    /// OCCT Selected() (lxx L19-22).
    pub fn selected(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_SELECTED) != 0
    }

    /// OCCT Selected(B) (lxx L26-32).
    pub fn set_selected(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_SELECTED;
        } else {
            self.my_flags &= !emask_flags::FMASK_SELECTED;
        }
    }

    /// OCCT Back() (lxx L36-39).
    pub fn back(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_BACK) != 0
    }

    /// OCCT Back(B) (lxx L43-49).
    pub fn set_back(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_BACK;
        } else {
            self.my_flags &= !emask_flags::FMASK_BACK;
        }
    }

    /// OCCT Side() (lxx L53-56).
    pub fn side(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_SIDE) != 0
    }

    /// OCCT Side(B) (lxx L60-66).
    pub fn set_side(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_SIDE;
        } else {
            self.my_flags &= !emask_flags::FMASK_SIDE;
        }
    }

    /// OCCT Closed() (lxx L70-73).
    pub fn closed(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_CLOSED) != 0
    }

    /// OCCT Closed(B) (lxx L77-83).
    pub fn set_closed(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_CLOSED;
        } else {
            self.my_flags &= !emask_flags::FMASK_CLOSED;
        }
    }

    /// OCCT Hiding() (lxx L87-90).
    pub fn hiding(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_HIDING) != 0
    }

    /// OCCT Hiding(B) (lxx L94-100).
    pub fn set_hiding(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_HIDING;
        } else {
            self.my_flags &= !emask_flags::FMASK_HIDING;
        }
    }

    /// OCCT Simple() (lxx L104-107).
    pub fn simple(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_SIMPLE) != 0
    }

    /// OCCT Simple(B) (lxx L111-117).
    pub fn set_simple(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_SIMPLE;
        } else {
            self.my_flags &= !emask_flags::FMASK_SIMPLE;
        }
    }

    /// OCCT Cut() (lxx L121-124).
    pub fn cut(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_CUT) != 0
    }

    /// OCCT Cut(B) (lxx L128-134).
    pub fn set_cut(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_CUT;
        } else {
            self.my_flags &= !emask_flags::FMASK_CUT;
        }
    }

    /// OCCT WithOutL() (lxx L138-141).
    pub fn with_out_l(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_WITH_OUT_L) != 0
    }

    /// OCCT WithOutL(B) (lxx L145-151).
    pub fn set_with_out_l(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_WITH_OUT_L;
        } else {
            self.my_flags &= !emask_flags::FMASK_WITH_OUT_L;
        }
    }

    /// OCCT Plane() (lxx L155-158).
    pub fn plane(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_PLANE) != 0
    }

    /// OCCT Plane(B) (lxx L162-168).
    pub fn set_plane(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_PLANE;
        } else {
            self.my_flags &= !emask_flags::FMASK_PLANE;
        }
    }

    /// OCCT Cylinder() (lxx L172-175).
    pub fn cylinder(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_CYLINDER) != 0
    }

    /// OCCT Cylinder(B) (lxx L179-185).
    pub fn set_cylinder(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_CYLINDER;
        } else {
            self.my_flags &= !emask_flags::FMASK_CYLINDER;
        }
    }

    /// OCCT Cone() (lxx L189-192).
    pub fn cone(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_CONE) != 0
    }

    /// OCCT Cone(B) (lxx L196-202).
    pub fn set_cone(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_CONE;
        } else {
            self.my_flags &= !emask_flags::FMASK_CONE;
        }
    }

    /// OCCT Sphere() (lxx L206-209).
    pub fn sphere(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_SPHERE) != 0
    }

    /// OCCT Sphere(B) (lxx L213-219).
    pub fn set_sphere(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_SPHERE;
        } else {
            self.my_flags &= !emask_flags::FMASK_SPHERE;
        }
    }

    /// OCCT Torus() (lxx L223-226).
    pub fn torus(&self) -> bool {
        (self.my_flags & emask_flags::FMASK_TORUS) != 0
    }

    /// OCCT Torus(B) (lxx L230-236).
    pub fn set_torus(&mut self, b: bool) {
        if b {
            self.my_flags |= emask_flags::FMASK_TORUS;
        } else {
            self.my_flags &= !emask_flags::FMASK_TORUS;
        }
    }

    /// OCCT Size() (lxx L240-243).
    pub fn size(&self) -> f64 {
        self.my_size
    }

    /// OCCT Size(S) (lxx L247-250).
    pub fn set_size(&mut self, s: f64) {
        self.my_size = s;
    }

    /// OCCT Orientation() (lxx L254-257) — the cast of the low EMaskOrient
    /// nibble.
    pub fn orientation(&self) -> Orientation {
        let v = self.my_flags & emask_flags::EMASK_ORIENT;
        match v {
            1 => Orientation::Reversed,
            2 => Orientation::Internal,
            3 => Orientation::External,
            _ => Orientation::Forward,
        }
    }

    /// OCCT Orientation(O) (lxx L261-265).
    pub fn set_orientation(&mut self, o: Orientation) {
        self.my_flags &= !emask_flags::EMASK_ORIENT;
        self.my_flags |= (o as i32) & emask_flags::EMASK_ORIENT;
    }

    /// OCCT Wires() (lxx L269-272) — the wire-block handle.  The OCCT
    /// `occ::handle<HLRAlgo_WiresBlock>&` dereference maps to the raw
    /// pointer (the HLR methodological exception, HLRBRep_Surface::myProj
    /// precedent); it is null before Set (the OCCT null handle).
    pub fn wires(&self) -> *mut WiresBlock {
        match &self.my_wires {
            Some(wb) => Arc::as_ptr(wb) as *const WiresBlock as *mut WiresBlock,
            None => std::ptr::null_mut(),
        }
    }

    /// OCCT Wires() (lxx L269-272) — the mutable handle reference (the
    /// `Wires() = new ...` assignment target; the shared-handle mutation
    /// goes through Arc::get_mut).
    pub fn change_wires(&mut self) -> &mut Option<Arc<WiresBlock>> {
        &mut self.my_wires
    }

    /// OCCT Geometry() (lxx L276-279) — the face surface (mutable, as in
    /// OCCT where Set loads it through Geometry().Surface(FG)).
    pub fn geometry(&mut self) -> &mut Surface<'a> {
        &mut self.my_geometry
    }

    /// The const Geometry() read (the OCCT callers bind the non-const
    /// reference; the rcad read-only contexts take this).
    pub fn geometry_ref(&self) -> &Surface<'a> {
        &self.my_geometry
    }

    /// OCCT Tolerance() (lxx L283-286).
    pub fn tolerance(&self) -> f32 {
        self.my_tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{Line3, Plane, Surface3};
    use rcad_kernel::topo::topods::BRepBuilder;

    /// A 2x2 planar face at z = `z` (the HLRBRep_Surface plane fixture).
    fn plane_face(brep: &mut rcad_kernel::BRep, z: f64) -> Shape {
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(brep, DVec3::new(0.0, 0.0, z), 1e-7);
        let v2 = b.add_vertex(brep, DVec3::new(2.0, 0.0, z), 1e-7);
        let e = b.add_edge(
            brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, z),
                direction: DVec3::new(1.0, 0.0, 0.0),
            })),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::new(0.0, 0.0, z),
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 2.0, 0.0, 2.0]),
            Vec::new(),
            true,
        )
    }

    /// OCCT anchor: the ctor sets only Selected (cxx L25-30); every FMask
    /// bit setter/getter round-trips (lxx, one assertion per bit) and the
    /// bits stay independent.
    #[test]
    fn face_data_fmask_bits_round_trip() {
        // ctor: myFlags = 0 then Selected(true) (cxx L26-29).
        let mut fd = FaceData::new();
        assert!(fd.selected());
        assert!(!fd.back());
        assert!(!fd.side());
        assert!(!fd.closed());
        assert!(!fd.hiding());
        assert!(!fd.simple());
        assert!(!fd.cut());
        assert!(!fd.with_out_l());
        assert!(!fd.plane());
        assert!(!fd.cylinder());
        assert!(!fd.cone());
        assert!(!fd.sphere());
        assert!(!fd.torus());

        // Selected(false) clears (lxx L26-32); the bit table order follows
        // the hxx L128-144 declaration order.
        fd.set_selected(false);
        assert!(!fd.selected());

        // Each FMask bit: set true, read back true.
        fd.set_back(true);
        assert!(fd.back());
        fd.set_side(true);
        assert!(fd.side());
        fd.set_closed(true);
        assert!(fd.closed());
        fd.set_hiding(true);
        assert!(fd.hiding());
        fd.set_simple(true);
        assert!(fd.simple());
        fd.set_cut(true);
        assert!(fd.cut());
        fd.set_with_out_l(true);
        assert!(fd.with_out_l());
        fd.set_plane(true);
        assert!(fd.plane());
        fd.set_cylinder(true);
        assert!(fd.cylinder());
        fd.set_cone(true);
        assert!(fd.cone());
        fd.set_sphere(true);
        assert!(fd.sphere());
        fd.set_torus(true);
        assert!(fd.torus());
        fd.set_selected(true);
        assert!(fd.selected());

        // Clearing one bit leaves the others (bit independence).
        fd.set_cut(false);
        assert!(!fd.cut());
        assert!(fd.selected());
        assert!(fd.back());
        assert!(fd.side());
        assert!(fd.closed());
        assert!(fd.hiding());
        assert!(fd.simple());
        assert!(fd.with_out_l());
        assert!(fd.plane());
        assert!(fd.cylinder());
        assert!(fd.cone());
        assert!(fd.sphere());
        assert!(fd.torus());

        // Size round-trip (lxx L240-250).
        assert_eq!(fd.size(), 0.0);
        fd.set_size(2.5);
        assert_eq!(fd.size(), 2.5);
    }

    /// OCCT anchor: the Orientation nibble holds the four TopAbs_Orientation
    /// values 0-3 (lxx L254-265) without clobbering the FMask bits above
    /// bit 3.
    #[test]
    fn face_data_orientation_nibble() {
        let mut fd = FaceData::new();
        assert_eq!(fd.orientation(), Orientation::Forward);

        fd.set_back(true);
        fd.set_orientation(Orientation::Reversed);
        assert_eq!(fd.orientation(), Orientation::Reversed);
        fd.set_orientation(Orientation::Internal);
        assert_eq!(fd.orientation(), Orientation::Internal);
        fd.set_orientation(Orientation::External);
        assert_eq!(fd.orientation(), Orientation::External);
        fd.set_orientation(Orientation::Forward);
        assert_eq!(fd.orientation(), Orientation::Forward);
        // The EMaskOrient rewrite kept the neighbouring FMask bit.
        assert!(fd.back());
    }

    /// OCCT anchor: Set (cxx L34-44) loads the geometry, the BRep_Tool
    /// tolerance and the WiresBlock; SetWire/SetWEdge (cxx L48-70) fill the
    /// blocks; Wires/Geometry/Tolerance read back (lxx L269-286).
    #[test]
    fn face_data_set_setwire_setwedge_round_trip() {
        let mut brep = rcad_kernel::BRep::new();
        let face = plane_face(&mut brep, 5.0);

        let mut fd = FaceData::new();
        fd.set(&brep, &face, Orientation::Reversed, true, 2);

        // Closed(Cl) (cxx L39).
        assert!(fd.closed());
        // Geometry().Surface(FG) — the plane classifies (cxx L40).
        assert_eq!(
            fd.geometry_ref().get_type(),
            crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane
        );
        // myTolerance = (float)(BRep_Tool::Tolerance(FG)) (cxx L41) — the
        // kernel face tolerance verbatim (the add_tface default is
        // f64::EPSILON, not Precision::Confusion()).
        assert_eq!(fd.tolerance(), brep.tolerance(&face) as f32);
        // Orientation(Or) (cxx L42).
        assert_eq!(fd.orientation(), Orientation::Reversed);
        // Wires() = new HLRAlgo_WiresBlock(NW) (cxx L43).
        let wb = unsafe { &*fd.wires() };
        assert_eq!(wb.nb_wires(), 2);

        // SetWire(WI, NE) (cxx L48-51).
        fd.set_wire(1, 3);
        fd.set_wire(2, 1);

        // SetWEdge(WI, EWI, EI, Or, OutL, Inte, Dble, IsoL) (cxx L55-70).
        fd.set_w_edge(1, 1, 10, Orientation::Reversed, true, false, true, false);
        fd.set_w_edge(1, 2, 11, Orientation::Forward, false, true, false, true);
        fd.set_w_edge(2, 1, 12, Orientation::Internal, true, true, true, true);

        // Round-trip through the wire blocks (WiresBlock::wire is the OCCT
        // mutable accessor).
        let wb = unsafe { &mut *fd.wires() };
        let w1 = &wb.wire(1);
        assert_eq!(w1.nb_edges(), 3);
        assert_eq!(w1.edge(1), 10);
        assert_eq!(w1.orientation(1), Orientation::Reversed);
        assert!(w1.out_line(1));
        assert!(!w1.internal(1));
        assert!(w1.double(1));
        assert!(!w1.iso_line(1));
        assert_eq!(w1.edge(2), 11);
        assert_eq!(w1.orientation(2), Orientation::Forward);
        assert!(w1.internal(2));
        assert!(w1.iso_line(2));
        let w2 = &wb.wire(2);
        assert_eq!(w2.nb_edges(), 1);
        assert_eq!(w2.edge(1), 12);
        assert_eq!(w2.orientation(1), Orientation::Internal);
        assert!(w2.out_line(1));
        assert!(w2.double(1));
    }
}
