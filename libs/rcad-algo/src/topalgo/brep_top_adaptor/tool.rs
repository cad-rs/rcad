// OCCT BRepTopAdaptor_Tool (TKTopAlgo/BRepTopAdaptor/BRepTopAdaptor_Tool.hxx
// L26-55 + .cxx L44-140) — the (TopolTool, surface) pair cached per face in
// the HLR MST map (NCollection_DataMap<TopoDS_Shape, BRepTopAdaptor_Tool>).
//
// Architecture notes:
// - OCCT myTopolTool is a handle created by every constructor; rcad owns the
//   BRepTopolTool value (the (brep, face)-backed adaptor per
//   topol_tool_brep.rs).
// - OCCT myHSurface is the handle<Adaptor3d_Surface> upcast of the
//   BRepAdaptor_Surface; rcad stores the restricted
//   `Arc<dyn SurfaceAdapter>` (the established `handle<Adaptor3d_Surface>`
//   mapping, hlr/contap/surface_adaptor.rs) — the null handle is `None`.
// - The OCCT TopolTool::Initialize(S) downcast to BRepAdaptor_Surface is
//   resolved at the call boundary: rcad's concrete BRepAdaptorSurface is the
//   downcast result, and the owning BRep (which BRepAdaptorSurface does not
//   carry) is passed alongside it where the surface-form API needs it.

use std::sync::Arc;

use rcad_kernel::topods::{BRep, Shape};

use crate::hlr::contap::surface_adaptor::SurfaceHandle;
use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;
use crate::topalgo::brep_top_adaptor::topol_tool_brep::BRepTopolTool;

/// OCCT BRepTopAdaptor_Tool.
// Clone: the OCCT HLR MST map (NCollection_DataMap::Bind, e.g.
// HLRTopoBRep_DSFiller.cxx L93 and HLRTopoBRep_OutLiner.cxx L301) stores a
// VALUE copy of the Tool (a copy of the two handles — the underlying
// myTopolTool / surface objects are shared); the rcad handle-copy equivalent
// is a plain Clone.
#[derive(Clone)]
pub struct BRepTopAdaptorTool {
    /// OCCT myloaded (hxx L51).
    my_loaded: bool,
    /// OCCT myTopolTool (hxx L52) — live after every constructor.
    my_topol_tool: BRepTopolTool,
    /// OCCT myHSurface (hxx L53) — the null handle until a loaded Init.
    my_h_surface: Option<SurfaceHandle>,
}

impl BRepTopAdaptorTool {
    /// OCCT BRepTopAdaptor_Tool() (cxx L44-49).
    pub fn new() -> Self {
        // OCCT myTopolTool = new BRepTopAdaptor_TopolTool(); the unloaded
        // tool wraps no surface — the empty BRep stands for the
        // uninitialized myS.
        BRepTopAdaptorTool {
            my_loaded: false,
            my_topol_tool: BRepTopolTool::new(Arc::new(BRep::new())),
            my_h_surface: None,
        }
    }

    /// OCCT BRepTopAdaptor_Tool(const TopoDS_Face& F, const double /*Tol2d*/)
    /// (cxx L51-63).
    pub fn new_face(brep: Arc<BRep>, f: &Shape, _tol2d: f64) -> Self {
        // OCCT myTopolTool = new BRepTopAdaptor_TopolTool().
        let mut my_topol_tool = BRepTopolTool::new(brep.clone());
        // OCCT surface = new BRepAdaptor_Surface(); surface->Initialize(F, true).
        let surface = BRepAdaptorSurface::initialize_face(&brep, f, true);
        // OCCT myTopolTool->Initialize(aSurf) — the downcast lands on the
        // loaded face; rcad's (brep, face)-backed TopolTool takes the face.
        my_topol_tool.initialize_surface(f);
        // OCCT myHSurface = surface (the upcast to Adaptor3d_Surface).
        let my_h_surface = Some(Arc::new(surface.adaptor_surface().clone()) as SurfaceHandle);
        BRepTopAdaptorTool {
            my_loaded: true,
            my_topol_tool,
            my_h_surface,
        }
    }

    /// OCCT BRepTopAdaptor_Tool(const occ::handle<Adaptor3d_Surface>&
    /// surface, const double /*Tol2d*/) (cxx L65-73).  rcad deviation: the
    /// owning BRep travels alongside the adaptor (BRepAdaptorSurface does
    /// not carry it) — the pair is the pre-downcast data of the OCCT
    /// signature.
    pub fn new_surface(brep: Arc<BRep>, surface: &BRepAdaptorSurface<'_>, _tol2d: f64) -> Self {
        // OCCT myTopolTool = new BRepTopAdaptor_TopolTool().
        let mut my_topol_tool = BRepTopolTool::new(brep);
        // OCCT myTopolTool->Initialize(surface).
        my_topol_tool.initialize_surface(surface.face());
        // OCCT myHSurface = surface; myloaded = true.
        let my_h_surface = Some(Arc::new(surface.adaptor_surface().clone()) as SurfaceHandle);
        BRepTopAdaptorTool {
            my_loaded: true,
            my_topol_tool,
            my_h_surface,
        }
    }

    /// OCCT Init(const TopoDS_Face& F, const double /*Tol2d*/) (cxx L77-88).
    pub fn init_face(&mut self, brep: Arc<BRep>, f: &Shape, _tol2d: f64) {
        // OCCT surface = new BRepAdaptor_Surface(); surface->Initialize(F)
        // (the Initialize Restriction flag defaults to true).
        let surface = BRepAdaptorSurface::initialize_face(&brep, f, true);
        // OCCT myHSurface = surface (the upcast to Adaptor3d_Surface).
        let my_h_surface = Some(Arc::new(surface.adaptor_surface().clone()) as SurfaceHandle);
        // OCCT myTopolTool->Initialize(aSurf) — rcad's tool carries the brep,
        // so the loaded pair is rebuilt with the same state reset as
        // Initialize(S).
        self.my_topol_tool = BRepTopolTool::new(brep);
        self.my_topol_tool.initialize_surface(f);
        self.my_h_surface = my_h_surface;
        self.my_loaded = true;
    }

    /// OCCT Init(const occ::handle<Adaptor3d_Surface>& surface, const double
    /// /*Tol2d*/) (cxx L90-97).  rcad deviation: the owning BRep travels
    /// alongside the adaptor, as in new_surface.
    pub fn init_surface(&mut self, brep: Arc<BRep>, surface: &BRepAdaptorSurface<'_>, _tol2d: f64) {
        // OCCT myHSurface = surface (the upcast to Adaptor3d_Surface).
        let my_h_surface = Some(Arc::new(surface.adaptor_surface().clone()) as SurfaceHandle);
        // OCCT myTopolTool->Initialize(surface).
        self.my_topol_tool = BRepTopolTool::new(brep);
        self.my_topol_tool.initialize_surface(surface.face());
        self.my_h_surface = my_h_surface;
        self.my_loaded = true;
    }

    /// OCCT GetTopolTool() (cxx L99-114) — myTopolTool is returned in both
    /// branches (the unloaded branch only prints the OCCT_DEBUG message).
    pub fn get_topol_tool(&mut self) -> &mut BRepTopolTool {
        &mut self.my_topol_tool
    }

    /// OCCT GetSurface() (cxx L116-131) — the surface handle (None encodes
    /// the null handle of the unloaded tool).
    pub fn get_surface(&self) -> Option<&SurfaceHandle> {
        self.my_h_surface.as_ref()
    }

    /// OCCT SetTopolTool(const occ::handle<BRepTopAdaptor_TopolTool>& TT)
    /// (cxx L133-136).
    pub fn set_topol_tool(&mut self, tt: BRepTopolTool) {
        self.my_topol_tool = tt;
    }

    /// OCCT Destroy() (cxx L138-140).
    pub fn destroy(&mut self) {}
}

impl Drop for BRepTopAdaptorTool {
    /// OCCT ~BRepTopAdaptor_Tool() { Destroy(); } (hxx L47).
    fn drop(&mut self) {
        self.destroy();
    }
}

impl Default for BRepTopAdaptorTool {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topalgo::brep_top_adaptor::topol_tool_brep::tests::square_face;

    /// OCCT anchor: the face constructor (cxx L51-63) — the TopolTool is
    /// initialized from the face (the restriction edge list is live) and
    /// GetSurface hands the restricted plane adaptor.
    #[test]
    fn tool_face_round_trip() {
        let (brep, face) = square_face();
        let mut tool = BRepTopAdaptorTool::new_face(Arc::new(brep), &face, 1e-7);
        assert!(tool.my_loaded);

        // GetTopolTool: the initialized tool iterates the face edges.
        let tt = tool.get_topol_tool();
        tt.init();
        assert!(tt.more());

        // GetSurface: the restricted plane window of the unit square.
        let surf = tool.get_surface().expect("loaded surface");
        assert_eq!(
            surf.get_type(),
            crate::geomalgo::int_curve_surface::SurfaceType::Plane
        );
        assert_eq!(surf.first_u_parameter(), 0.0);
        assert_eq!(surf.last_u_parameter(), 1.0);
        assert_eq!(surf.first_v_parameter(), 0.0);
        assert_eq!(surf.last_v_parameter(), 1.0);
    }

    /// OCCT anchor: the default constructor (cxx L44-49) — myloaded is
    /// false, the surface handle is null (GetSurface reports None) and the
    /// TopolTool is alive but empty (the OCCT_DEBUG branch of GetTopolTool).
    #[test]
    fn tool_default_unloaded() {
        let mut tool = BRepTopAdaptorTool::new();
        assert!(!tool.my_loaded);
        assert!(tool.get_surface().is_none());
        let tt = tool.get_topol_tool();
        tt.init();
        assert!(!tt.more());
    }

    /// OCCT anchor: Init(F) (cxx L77-88) reloads an existing tool and
    /// SetTopolTool (cxx L133-136) replaces the tool handle.
    #[test]
    fn tool_init_and_set_topol_tool() {
        let (brep, face) = square_face();
        let mut tool = BRepTopAdaptorTool::new();
        assert!(!tool.my_loaded);
        tool.init_face(Arc::new(brep), &face, 1e-7);
        assert!(tool.my_loaded);
        assert!(tool.get_surface().is_some());
        let tt = tool.get_topol_tool();
        tt.init();
        assert!(tt.more());

        // SetTopolTool swaps the handle (the fresh tool is empty again).
        tool.set_topol_tool(BRepTopolTool::new(Arc::new(BRep::new())));
        let tt = tool.get_topol_tool();
        tt.init();
        assert!(!tt.more());
        // myloaded is untouched by SetTopolTool (OCCT leaves it true).
        assert!(tool.my_loaded);
    }
}
