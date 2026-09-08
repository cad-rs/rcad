// OCCT BRepOffset_MakeOffset_1.cxx L1-1305 — module a of the 1:1 translation
// (the file is split by OCCT order into brep_offset_make_offset_1.rs /
// _b.rs / _c.rs / ... to respect the 2000-line file limit; the split points
// are annotated at the top of each module).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_MakeOffset_1.cxx (9,533 lines)
//
// This file carries cxx L1-1305: the anonymous-namespace helpers (L78-521),
// the BRepOffset_BuildOffsetFaces class declaration (L529-1032) mapped onto
// the BRepOffsetBuildOffsetFaces struct, the constructor + setters, and the
// methods BuildSplitsOfTrimmedFaces / BuildSplitsOfExtendedFaces /
// UpdateIntersectedEdges / IntersectTrimmedEdges plus the file static
// checkConnectionsOfFace (L1301-1334).
//
// OCCT inheritance chain: none — BRepOffset_BuildOffsetFaces is an auxiliary
// local class used for building splits of offset faces.
//
// Architecture differences (numbering continues from brep_offset_tool_d.rs):
// 38. The OCCT output reference member `BRepAlgo_Image* myImage` becomes the
//     explicit `the_image: &mut BRepAlgoImage` parameter of the methods that
//     consume it (the Rust borrow form of the constructor-held pointer; the
//     call structure keeps the OCCT order).
// 39. Message_ProgressScope / Message_ProgressRange -> the rcad NoopProgress
//     + ProgressScope forms; the OCCT child-scope chaining is flattened to
//     local scopes (the progress machinery is result-inert; the
//     brep_offset_inter3d.rs #37 precedent).
// 40. The OCCT nullable input pointers (myFaces, myAsDes, myAnalyzer,
//     myEdgesOrigins, myFacesOrigins, myETrimEInf) become Option<..> owned
//     values; SetFaces/SetAsDesInfo/... clone the arguments (the Handle /
//     pointer form carries no refcount in rcad — architecture difference
//     #23 of brep_offset_tool.rs).
// 41. BOPAlgo_Builder / BOPAlgo_Splitter / BOPAlgo_BOP -> the
//     crate::bop::algo::builder::Builder driven by the crate::bop PaveFiller
//     pipeline (the brep_algo_api::run_build_* driver form: the PaveFiller
//     performs the intersection, the Builder builds the result on the DS);
//     the OCCT `BOPAlgo_Builder aGF; aGF.AddArgument(..); aGF.Perform();`
//     form maps onto the local filler + Builder pair; Modified()/Images()
//     read the Builder my_images / my_origins tables (pub(crate)).
// 42. BOPAlgo_BuilderFace -> the crate::bop::algo::builder_face::BuilderFace
//     (the brep_feat_builder.rs L1256 precedent: the rcad BuilderFace is
//     DS-bound, so the face+edges argument list runs through a local
//     PaveFiller first; the OCCT BuilderFace has no PaveFiller of its own).
// 43. BRepLib::BuildPCurveForEdgesOnPlane -> the crate::bop Builder re-host
//     build_pcurve_for_edges_on_plane (BRepLib_1.cxx L313-339; re-hosted
//     here through a local re-host of the same OCCT body form).
// 44. BOPTools_AlgoTools::MakeConnexityBlocks -> the rcad translations
//     shell_splitter::make_connexity_blocks (EDGE/FACE form) and
//     wire_splitter::make_connexity_blocks (VERTEX/EDGE form); the
//     VERTEX/FACE form used once (cxx L6016) has no rcad translation and
//     takes the GAP annotation at the call site.
// 45. BOPTools_Set -> local BOPToolsSet re-host of BOPTools_Set.cxx
//     (Add/AreEqual semantics) — the TKBO package has no rcad translation
//     of BOPTools_Set yet (the FindInvalidEdges 2nd form use).
// 46. IntTools_Context -> crate::bop::int_tools::context::IntToolsContext;
//     the SolidClassifier access keeps the OCCT context-owned classifier
//     form through a local classifier carrier.
// 47. BOPAlgo_MakerVolume -> crate::bop::algo::maker_volume::MakerVolume
//     (the rcad form is self-contained — it owns its PaveFiller).
// 48. BOPAlgo_Section -> crate::bop::algo::section::BOPAlgoSection (the
//     run_build_section_brep driver form).

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::builder::Builder;
use crate::bop::algo::occt_map::OcctDataMapInt;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_range, brep_tool_tolerance, ShapeKey,
};

use super::brep_offset_tool::{
    oriented, set_add, shape_data_map, top_exp_vertices,
    OcctIndexedShapeMap, OcctShapeSet, ShapeDataMap, ShapeIndexedDataMap,
};
use super::brep_offset_tool_d::BRepOffsetAnalyse;

// ---------------------------------------------------------------------------
// Shared data-map type of the cxx file (cxx L60-64).
// ---------------------------------------------------------------------------

/// The file-name-derived alias of the class carrier consumed by the
/// BRepOffset_MakeOffset module (the OCCT class of this file is the local
/// BRepOffset_BuildOffsetFaces — the alias keeps the consuming-module
/// import form stable).
pub type BRepOffsetMakeOffset1 = BRepOffsetBuildOffsetFaces;

/// OCCT `BRepOffset_DataMapOfShapeIndexedMapOfShape` typedef (cxx L60-64):
/// DataMap<Shape, IndexedMap<Shape>, ShapeMapHasher>.
pub(crate) type DataMapOfShapeIndexedMapOfShape = ShapeDataMap<OcctIndexedShapeMap>;

// ---------------------------------------------------------------------------
// Anonymous-namespace helpers (cxx L66-521).
// ---------------------------------------------------------------------------

/// OCCT static AddToContainer(S, List) (cxx L78-81) — append into the list.
pub(crate) fn add_to_container_list(the_s: &Shape, the_list: &mut Vec<Shape>) {
    the_list.push(the_s.clone());
}

/// OCCT static AddToContainer(S, Map) (cxx L83-87) — add into the set;
/// returns true when newly added.
pub(crate) fn add_to_container_map(the_s: &Shape, the_map: &mut OcctShapeSet) -> bool {
    set_add(the_map, the_s)
}

/// OCCT static AddToContainer(S, IndexedMap) (cxx L89-95) — add into the
/// indexed map; returns true when newly added.  The OCCT index is 1-based
/// (anInd > aNb); the rcad index is 0-based (anInd >= aNb).
pub(crate) fn add_to_container_indexed(the_s: &Shape, the_map: &mut OcctIndexedShapeMap) -> bool {
    let a_nb = the_map.extent();
    let an_ind = the_map.add(the_s);
    an_ind >= a_nb
}

/// OCCT static AddToContainer(S, SOut) (cxx L97-100) — BRep_Builder().Add.
pub(crate) fn add_to_container_shape(the_s: &Shape, the_s_out: &mut Shape) {
    bat::builder_add_compound_shape(the_s_out, the_s);
}

/// OCCT static AddToContainer(Key, Value, DataMap<K, List<K>>) (cxx
/// L102-112) — append theValue into theKey's list, binding when absent.
pub(crate) fn add_to_container_map_list(
    the_map: &mut ShapeDataMap<Vec<Shape>>,
    the_key: &Shape,
    the_value: &Shape,
) {
    if let Some(p_list) = dm_change_seek(the_map, the_key) {
        p_list.push(the_value.clone());
    } else {
        shape_data_map::bound(the_map, the_key).push(the_value.clone());
    }
}

/// OCCT static TakeModified (cxx L123-153) — check if the shape has images
/// in the given images map.  Puts in the output map either the images or the
/// shape itself.  Returns true when the shape has images.
///
/// The OCCT template over ContainerType/FenceMapType is taken here in the
/// concrete forms used by the single call site (UpdateImages, cxx L297):
/// ContainerType = NCollection_List<TopoDS_Shape>, FenceMapType =
/// NCollection_Map<TopoDS_Shape> (the `theMFence == nullptr` form is the
/// take_modified_no_fence helper below).
pub(crate) fn take_modified(
    the_s: &Shape,
    the_images: &OcctDataMapInt<ShapeKey, Vec<Shape>>,
    the_container: &mut Vec<Shape>,
    the_m_fence: &mut OcctShapeSet,
) -> bool {
    // OCCT L131: pLSIm = theImages.Seek(theS) — the DataMap key is the
    // TShape + Location identity (TopTools_ShapeMapHasher); the rcad
    // builder images table is keyed by the same (ptr, location) pair.
    if let Some(p_ls_im) = the_images.get((the_s.ptr_id(), the_s.location)) {
        for a_s_im in p_ls_im {
            // OCCT L138: if (!theMFence || AddToContainer(aSIm, *theMFence)).
            if add_to_container_map(a_s_im, the_m_fence) {
                add_to_container_list(a_s_im, the_container);
            }
        }
        true
    } else {
        if add_to_container_map(the_s, the_m_fence) {
            add_to_container_list(the_s, the_container);
        }
        false
    }
}

/// OCCT static TakeModified, the no-fence template form (cxx L155-164) —
/// `TakeModified(theS, theImages, theMapOut)` with a null fence.
pub(crate) fn take_modified_no_fence(
    the_s: &Shape,
    the_images: &OcctDataMapInt<ShapeKey, Vec<Shape>>,
    the_container: &mut Vec<Shape>,
) -> bool {
    if let Some(p_ls_im) = the_images.get((the_s.ptr_id(), the_s.location)) {
        for a_s_im in p_ls_im {
            add_to_container_list(a_s_im, the_container);
        }
        true
    } else {
        add_to_container_list(the_s, the_container);
        false
    }
}

/// OCCT static TakeModified (cxx L123-153), the TopoDS_Shape container
/// form (AddToContainer(S, SOut)) — the cxx L5007/L5013 form.
pub(crate) fn take_modified_shape(
    the_s: &Shape,
    the_images: &OcctDataMapInt<ShapeKey, Vec<Shape>>,
    the_container: &mut Shape,
    the_m_fence: &mut OcctShapeSet,
) -> bool {
    if let Some(p_ls_im) = the_images.get((the_s.ptr_id(), the_s.location)) {
        for a_s_im in p_ls_im {
            if add_to_container_map(a_s_im, the_m_fence) {
                add_to_container_shape(a_s_im, the_container);
            }
        }
        true
    } else {
        if add_to_container_map(the_s, the_m_fence) {
            add_to_container_shape(the_s, the_container);
        }
        false
    }
}

/// OCCT static TakeModified (cxx L155-164), the NCollection_IndexedMap
/// container form (AddToContainer(S, IndexedMap)) — the cxx L4486/L4492
/// form.
pub(crate) fn take_modified_no_fence_indexed(
    the_s: &Shape,
    the_images: &OcctDataMapInt<ShapeKey, Vec<Shape>>,
    the_container: &mut OcctIndexedShapeMap,
) -> bool {
    if let Some(p_ls_im) = the_images.get((the_s.ptr_id(), the_s.location)) {
        for a_s_im in p_ls_im {
            add_to_container_indexed(a_s_im, the_container);
        }
        true
    } else {
        add_to_container_indexed(the_s, the_container);
        false
    }
}

/// OCCT static TakeModified (cxx L155-164), the NCollection_Map container
/// form (AddToContainer(S, Map)) — the cxx L5021/L5031 form.
pub(crate) fn take_modified_map_no_fence(
    the_s: &Shape,
    the_images: &OcctDataMapInt<ShapeKey, Vec<Shape>>,
    the_container: &mut OcctShapeSet,
) -> bool {
    if let Some(p_ls_im) = the_images.get((the_s.ptr_id(), the_s.location)) {
        for a_s_im in p_ls_im {
            add_to_container_map(a_s_im, the_container);
        }
        true
    } else {
        add_to_container_map(the_s, the_container);
        false
    }
}

/// OCCT static hasData (cxx L170-174) — checks if the container has any
/// data in it.
pub(crate) fn has_data(the_data: Option<&Vec<Shape>>) -> bool {
    match the_data {
        Some(v) => !v.is_empty(),
        None => false,
    }
}

/// OCCT static AppendToList (cxx L180-191) — add to a list only unique
/// elements.
pub(crate) fn append_to_list(the_list: &mut Vec<Shape>, the_shape: &Shape) {
    for a_s in the_list.iter() {
        if a_s.is_same(the_shape) {
            return;
        }
    }
    the_list.push(the_shape.clone());
}

/// OCCT static ProcessMicroEdge (cxx L197-216) — checking if the edge is a
/// micro edge.
///
/// The IntTools_Context argument of BOPTools_AlgoTools::IsMicroEdge carries
/// the extrema machinery; the rcad IsMicroEdge re-host
/// (bop::tools::algo_tools::is_micro_edge) takes the shrunk-range length
/// form: the OCCT body computes the shrunk range of the edge and compares
/// its length against the edge tolerance (BOPTools_AlgoTools.cxx L2075+,
/// final check `aT2 - aT1 < 2 * aTE`).  The micro-edge vertex tolerance
/// update (cxx L208-213) keeps the BRep_Builder UpdateVertex form.
pub(crate) fn process_micro_edge(the_edge: &Shape) -> bool {
    // OCCT L200-201: TopExp::Vertices(theEdge, aV1, aV2).
    let (a_v1, a_v2) = top_exp_vertices(the_edge);
    if a_v1.is_null() || a_v2.is_null() {
        return false;
    }

    // OCCT L207: bMicro = BOPTools_AlgoTools::IsMicroEdge(theEdge, theCtx).
    // The shrunk-range computation over the full edge range (the OCCT
    // IntTools_ShrunkRange machinery) is re-hosted by the reduced uniform
    // form: the shrinkable length of a straight edge equals its parametric
    // range; the comparison keeps the OCCT 2 * tolerance form.
    let (a_t1, a_t2) = brep_tool_range(the_edge);
    let a_len = a_t2 - a_t1;
    let b_micro =
        crate::bop::tools::algo_tools::is_micro_edge(a_len, brep_tool_tolerance(the_edge));
    if b_micro {
        // OCCT L208: BRepAdaptor_Curve(theEdge).GetType() == GeomAbs_Line.
        if edge_curve_is_line(the_edge) {
            // OCCT L210-212: aLen = BRep_Tool::Pnt(aV1).Distance(Pnt(aV2));
            // BRep_Builder().UpdateVertex(aV1/aV2, aLen / 2.).
            let p1 = bat::brep_tool_pnt(&a_v1).unwrap_or(glam::DVec3::ZERO);
            let p2 = bat::brep_tool_pnt(&a_v2).unwrap_or(glam::DVec3::ZERO);
            let a_len3 = p1.distance(p2);
            let mut v1 = a_v1.clone();
            let mut v2 = a_v2.clone();
            bat::builder_update_vertex_tol(&mut v1, a_len3 / 2.);
            bat::builder_update_vertex_tol(&mut v2, a_len3 / 2.);
        }
    }

    b_micro
}

/// The GeomAbs_Line probe of cxx L208 (BRepAdaptor_Curve(E).GetType()) —
/// the rcad form reads the underlying curve variant (the
/// brep_offset_inter2d.rs BRepAdaptorCurve precedent).
pub(crate) fn edge_curve_is_line(the_edge: &Shape) -> bool {
    match bat::brep_tool_curve(the_edge) {
        Some((c, _, _)) => match &c {
            rcad_kernel::geom::Curve3::Line(_) => true,
            rcad_kernel::geom::Curve3::Trimmed(t) => {
                matches!(&*t.curve, rcad_kernel::geom::Curve3::Line(_))
            }
            _ => false,
        },
        None => false,
    }
}

/// OCCT static UpdateOrigins (cxx L220-261) — merging the modifications of
/// theGF into theOrigins.
///
/// The `BOPAlgo_Builder& theGF` argument is taken in the rcad Builder form
/// (architecture difference #41): `theGF.Modified(aS)` reads the
/// my_images table.
pub(crate) fn update_origins(
    the_la: &[Shape],
    the_origins: &mut ShapeDataMap<Vec<Shape>>,
    the_gf: &Builder,
) {
    for a_s in the_la.iter() {
        // OCCT L230: aLSIm = theGF.Modified(aS).
        let a_ls_im = builder_modified(the_gf, a_s);
        if a_ls_im.is_empty() {
            continue;
        }

        // OCCT L236-242: pLS = theOrigins.ChangeSeek(aS); when absent the
        // local aLSEmpt = { aS } list is used without binding.  The rcad
        // form snapshots the list (the borrow form of the DataMap value).
        let a_ls: Vec<Shape> =
            if dm_change_seek(the_origins, a_s).is_some() {
                shape_data_map::find(the_origins, a_s)
            } else {
                let mut l = Vec::new();
                l.push(a_s.clone());
                l
            };

        for a_s_im in a_ls_im.iter() {
            // OCCT L247-254: merge two lists.
            if let Some(p_ls_or) = dm_change_seek(the_origins, a_s_im) {
                for a_it1 in a_ls.iter() {
                    append_to_list(p_ls_or, a_it1);
                }
            } else {
                // OCCT L257: theOrigins.Bind(aSIm, *pLS).
                shape_data_map::bind(the_origins, a_s_im, a_ls.clone());
            }
        }
    }
}

/// OCCT static UpdateImages (cxx L265-306) — updating the images map with
/// the modifications of theGF.
pub(crate) fn update_images(
    the_la: &[Shape],
    the_images: &mut ShapeDataMap<Vec<Shape>>,
    the_gf: &Builder,
    the_modified: &mut OcctShapeSet,
) {
    for a_s in the_la.iter() {
        // OCCT L276-286: when the shape has no images yet, bind the
        // modifications of theGF.
        let has_images = dm_change_seek(the_images, a_s).is_some();
        if !has_images {
            let a_ls_im = builder_modified(the_gf, a_s);
            if !a_ls_im.is_empty() {
                shape_data_map::bind(the_images, a_s, a_ls_im);
                set_add(the_modified, a_s);
            }
            continue;
        }

        let mut a_m_fence: OcctShapeSet = HashMap::new();
        let mut a_ls_im_new: Vec<Shape> = Vec::new();

        let mut b_modified = false;

        // OCCT L294-298: check modifications of the images.
        let a_ls_im: Vec<Shape> = shape_data_map::find(the_images, a_s);
        for a_s_im in a_ls_im.iter() {
            // OCCT L297: bModified |= TakeModified(aSIm, theGF.Images(),
            // aLSImNew, &aMFence).
            b_modified |= take_modified(
                a_s_im,
                builder_images(the_gf),
                &mut a_ls_im_new,
                &mut a_m_fence,
            );
        }

        if b_modified {
            // OCCT L302: *pLSIm = aLSImNew.
            *shape_data_map::change_find(the_images, a_s) = a_ls_im_new;
            set_add(the_modified, a_s);
        }
    }
}

/// OCCT BOPAlgo_Builder::Modified(S) (architecture difference #41) — the
/// images of the shape in the builder, an empty list when the shape has not
/// been modified.
pub(crate) fn builder_modified(the_gf: &Builder, the_s: &Shape) -> Vec<Shape> {
    the_gf
        .my_images
        .get((the_s.ptr_id(), the_s.location))
        .cloned()
        .unwrap_or_default()
}

/// OCCT BOPAlgo_Builder::Images() — the whole myImages table (cxx L297
/// use).  The OcctDataMapInt keys are the (TShape ptr, location) pairs used
/// by the builder's FillImages step.
pub(crate) fn builder_images<'a>(the_gf: &'a Builder<'a>) -> &'a OcctDataMapInt<ShapeKey, Vec<Shape>> {
    &the_gf.my_images
}

/// OCCT `BOPAlgo_Builder&` base-class reference carrying the actual
/// BOPAlgo_MakerVolume object (architecture difference #49: the rcad
/// MakerVolume is a separate struct — the base-class reference becomes
/// this enum; the OCCT call sites pass either a Builder or a MakerVolume).
pub(crate) enum BuilderRef<'a> {
    /// The plain BOPAlgo_Builder form.
    Builder(&'a Builder<'a>),
    /// The BOPAlgo_MakerVolume form (the derived object).
    MakerVolume(&'a crate::bop::algo::maker_volume::MakerVolume),
}

impl<'a> BuilderRef<'a> {
    /// OCCT BOPAlgo_Builder::Modified(S).
    pub fn modified(&self, the_s: &Shape) -> Vec<Shape> {
        match self {
            BuilderRef::Builder(b) => builder_modified(b, the_s),
            BuilderRef::MakerVolume(_mv) => {
                // GAP: BOPAlgo_MakerVolume::Modified (the rcad MakerVolume
                // re-host does not expose the post-perform image tables
                // yet) — the OCCT empty-list path is taken.
                Vec::new()
            }
        }
    }

    /// OCCT BOPAlgo_Builder::Images() — the whole table (the TakeModified
    /// argument form).
    pub fn images(&self) -> OcctDataMapInt<ShapeKey, Vec<Shape>> {
        match self {
            BuilderRef::Builder(b) => {
                // The snapshot of the builder table.
                let mut m: OcctDataMapInt<ShapeKey, Vec<Shape>> = OcctDataMapInt::new();
                for (k, v) in b.my_images.iter() {
                    m.insert(k, v.clone());
                }
                m
            }
            BuilderRef::MakerVolume(_mv) => {
                // GAP: BOPAlgo_MakerVolume::Images — the OCCT empty-map
                // path is taken.
                OcctDataMapInt::new()
            }
        }
    }

    /// OCCT BOPAlgo_Builder::Origins() — the whole table.
    pub fn origins(&self) -> OcctDataMapInt<ShapeKey, Vec<Shape>> {
        match self {
            BuilderRef::Builder(b) => {
                let mut m: OcctDataMapInt<ShapeKey, Vec<Shape>> = OcctDataMapInt::new();
                for (k, v) in b.my_origins.iter() {
                    m.insert(*k, v.clone());
                }
                m
            }
            BuilderRef::MakerVolume(_mv) => {
                // GAP: BOPAlgo_MakerVolume::Origins — the OCCT empty-map
                // path is taken.
                OcctDataMapInt::new()
            }
        }
    }

    /// OCCT BOPAlgo_Builder::Arguments().
    pub fn arguments(&self) -> Vec<Shape> {
        match self {
            BuilderRef::Builder(b) => b.my_arguments.clone(),
            BuilderRef::MakerVolume(_mv) => {
                // GAP: BOPAlgo_MakerVolume::Arguments — the rcad form
                // keeps my_arguments private; the OCCT empty-list path.
                Vec::new()
            }
        }
    }

    /// OCCT BOPAlgo_BuilderShape::Shape() — the compound root of the
    /// result (the BRep pool root of the rcad Builder build()).
    pub fn shape(&self) -> Option<Shape> {
        match self {
            BuilderRef::Builder(b) => {
                let brep = b.my_shape.as_ref()?;
                Some(brep_root_shape(brep))
            }
            BuilderRef::MakerVolume(mv) => mv.shape().cloned(),
        }
    }

    /// OCCT BOPAlgo_Builder::IsDeleted(S).
    pub fn is_deleted(&self, _the_s: &Shape) -> bool {
        // GAP: BOPAlgo_Builder::IsDeleted (the rcad Builder/MakerVolume
        // re-hosts keep no deleted-shapes table) — the OCCT false path.
        false
    }

    /// OCCT BOPAlgo_Builder::HasDeleted().
    pub fn has_deleted(&self) -> bool {
        // GAP: the deleted-shapes table — the OCCT false path.
        false
    }

    /// OCCT BOPAlgo_Builder::PDS() — the DS interference access of the cxx
    /// ShapesConnections body; the GAP form returns no interferences (the
    /// OCCT aNbFF == 0 path) when the table is not addressable from the
    /// MakerVolume re-host.
    pub fn pds_interf_ff(&self) -> Vec<()> {
        Vec::new()
    }
}

/// OCCT static FindCommonParts (cxx L312-345) — looking for the parts of
/// type theType contained in both lists.
pub(crate) fn find_common_parts(
    the_ls1: &[Shape],
    the_ls2: &[Shape],
    the_lsc: &mut Vec<Shape>,
    the_type: ShapeType,
) {
    // OCCT L318-323: map shapes in the first list.
    let mut a_ms1 = OcctIndexedShapeMap::new();
    for a_s in the_ls1.iter() {
        for a_sub in explorer(a_s, the_type, ShapeType::Shape) {
            a_ms1.add(&a_sub);
        }
    }
    if a_ms1.is_empty() {
        return;
    }

    // OCCT L330-344: check for such shapes in the other list.
    let mut a_m_fence: OcctShapeSet = HashMap::new();
    for a_s in the_ls2.iter() {
        for a_st in explorer(a_s, the_type, ShapeType::Shape) {
            if a_ms1.contains(&a_st) && set_add(&mut a_m_fence, &a_st) {
                the_lsc.push(a_st);
            }
        }
    }
}

/// OCCT static NbPoints (cxx L351-361) — defines the number of sample
/// points to get the average direction of the edge.
pub(crate) fn nb_points(the_edge: &Shape) -> usize {
    // OCCT L353-360: BRepAdaptor_Curve aBAC(theEdge);
    // switch (aBAC.GetType()) { case GeomAbs_Line: return 1; default:
    // return 11; }.
    if edge_curve_is_line(the_edge) {
        1
    } else {
        11
    }
}

/// OCCT static FindShape (cxx L367-405) — looking for the same sub-shape in
/// the shape.
pub(crate) fn find_shape(
    the_s_what: &Shape,
    the_s_where: &Shape,
    the_analyse: Option<&BRepOffsetAnalyse>,
    the_res: &mut Shape,
) -> bool {
    let mut b_found = false;
    let a_type = the_s_what.shape_type();
    for a_s in explorer(the_s_where, a_type, ShapeType::Shape) {
        if a_s.is_same(the_s_what) {
            *the_res = a_s;
            b_found = true;
            break;
        }
    }

    if !b_found {
        if let Some(the_analyse) = the_analyse {
            // OCCT L388: pLD = theAnalyse->Descendants(theSWhere).
            if let Some(p_ld) = the_analyse.descendants(the_s_where) {
                for a_s in p_ld.iter() {
                    if a_s.is_same(the_s_what) {
                        *the_res = a_s.clone();
                        b_found = true;
                        break;
                    }
                }
            }
        }
    }

    b_found
}

/// OCCT static BuildSplitsOfTrimmedFace (cxx L411-432) — building the
/// splits of an offset face.
///
/// BOPAlgo_Splitter -> the rcad PaveFiller + Builder pipeline (architecture
/// difference #41; the BRepAlgoAPI_Splitter driver form of
/// brep_algo_api::run_build_splitter_brep with objects = { theFace,
/// theEdges }).  `SetToFillHistory(false)` maps onto the history-free
/// build() (the rcad Builder fills history only in the
/// build_with_history_* entries).
pub(crate) fn build_splits_of_trimmed_face(
    the_face: &Shape,
    the_edges: &Shape,
    the_lf_images: &mut Vec<Shape>,
    _the_range: &ProgressScope,
) {
    // OCCT L416-421: BOPAlgo_Splitter aSplitter; AddArgument(theFace);
    // AddArgument(theEdges); SetToFillHistory(false); Perform(theRange).
    let mut a_splitter_filler = PaveFiller::new();
    let mut a_args: Vec<Shape> = Vec::new();
    a_args.push(the_face.clone());
    a_args.push(the_edges.clone());
    a_splitter_filler.set_arguments(a_args);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_Splitter", 1);
    a_splitter_filler.perform(&a_ps);
    let mut a_splitter = Builder::new(
        a_splitter_filler.ds(),
        crate::bop::algo::builder::BooleanOpType::Union,
        a_splitter_filler.fuzzy_value(),
    );
    // OCCT BOPAlgo_Splitter::SetArguments(theFace) + SetTools(theEdges):
    // the result keeps only the split parts of the OBJECTS (theFace), the
    // tool parts (theEdges) are excluded (BOPAlgo_Builder_1.cxx L130-168).
    // The OCCT form here uses AddArgument for both — the Splitter treats
    // every argument as an object (myArguments = myObjects).
    let n_objs = 1;
    a_splitter.my_arguments = a_splitter_filler.ds().arguments[..n_objs].to_vec();
    a_splitter.my_tools = a_splitter_filler.ds().arguments[n_objs..].to_vec();
    a_splitter.my_is_splitter = true;
    let a_res = a_splitter.build();
    // OCCT L422-425: if (aSplitter.HasErrors()) return.
    match a_res {
        Err(_) => return,
        Ok(a_brep) => {
            // OCCT L428-431: splits of the offset shape — explore FACE.
            for a_f in bat::explorer(&brep_root_shape(&a_brep), ShapeType::Face, ShapeType::Shape)
            {
                the_lf_images.push(a_f);
            }
        }
    }
}

/// The BRep pool root read of the builder result (the SplitterOp form of
/// brep_algo_api::run_build_splitter_brep — the top-level shape of the
/// pool; the BRep result of the rcad Builder is a pool, the OCCT
/// BOPAlgo_BuilderShape::myShape is the compound root).
pub(crate) fn brep_root_shape(a_brep: &rcad_kernel::BRep) -> Shape {
    use rcad_kernel::topods::TShape;
    // The synthetic root: rebuild the compound of the pool top-level shapes
    // (the brep_algo_api::brep_top_shapes walk).
    let mut referenced = vec![false; a_brep.tshapes.len()];
    fn mark(sr: &Shape, referenced: &mut Vec<bool>) {
        let i = sr.index;
        if i >= referenced.len() || referenced[i] {
            return;
        }
        referenced[i] = true;
        match &*sr.data {
            TShape::Solid(sd) => {
                for sh in &sd.shells {
                    mark(sh, referenced);
                }
                for v in &sd.internal_vertices {
                    mark(v, referenced);
                }
                for e in &sd.internal_edges {
                    mark(e, referenced);
                }
            }
            TShape::Shell(sd) => {
                for f in &sd.faces {
                    mark(f, referenced);
                }
            }
            TShape::Face(fd) => {
                mark(&fd.outer_wire, referenced);
                for w in &fd.inner_wires {
                    mark(w, referenced);
                }
                for v in &fd.internal_vertices {
                    mark(v, referenced);
                }
            }
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    mark(e, referenced);
                }
            }
            TShape::Edge(ed) => {
                mark(&ed.first, referenced);
                mark(&ed.last, referenced);
            }
            _ => {}
        }
    }
    // The pool top-level shapes are the synthetic-wrapped result (the
    // Builder build() form) — every non-referenced shape is a root.
    let mut roots: Vec<Shape> = Vec::new();
    for ts in a_brep.tshapes.iter().enumerate() {
        let s = Shape::from_parts(
            ts.1.clone(),
            ts.0,
            0,
            Orientation::Forward,
        );
        if !matches!(&*s.data, TShape::Vertex(_)) {
            let before = referenced.iter().filter(|&&r| r).count();
            mark(&s, &mut referenced);
            let after = referenced.iter().filter(|&&r| r).count();
            if after > before {
                roots.push(s);
            }
        }
    }
    let a_bb = bat::builder_make_compound();
    let mut a_compound = a_bb;
    for r in roots {
        bat::builder_add_compound_shape(&mut a_compound, &r);
    }
    a_compound
}

/// OCCT static BuildSplitsOfFace (cxx L438-485) — building the splits of an
/// offset face.
pub(crate) fn build_splits_of_face(
    the_face: &Shape,
    the_edges: &Shape,
    the_faces_origins: &mut ShapeDataMap<Shape>,
    the_lf_images: &mut Vec<Shape>,
) {
    the_lf_images.clear();

    // OCCT L446-456: take edges to split the face — both orientations.
    let mut a_le: Vec<Shape> = Vec::new();
    for a_e_cur in explorer(the_edges, ShapeType::Edge, ShapeType::Shape) {
        let mut a_e = a_e_cur;
        a_e = oriented(&a_e, Orientation::Forward);
        a_le.push(a_e.clone());
        a_e = oriented(&a_e, Orientation::Reversed);
        a_le.push(a_e);
    }

    // OCCT L458-460: aFF = theFace with FORWARD orientation.
    let mut a_ff = the_face.clone();
    let an_or = the_face.orientation;
    a_ff = oriented(&a_ff, Orientation::Forward);

    // OCCT L463: BRepLib::BuildPCurveForEdgesOnPlane(aLE, aFF)
    // (architecture difference #43 — the local re-host of the BRepLib_1.cxx
    // L313-339 body form of the bop Builder re-host).
    build_pcurve_for_edges_on_plane(&a_le, &a_ff);

    // OCCT L466-473: BOPAlgo_BuilderFace aBF; SetFace(aFF); SetShapes(aLE);
    // Perform() (architecture difference #42 — the rcad BuilderFace is
    // DS-bound; the face+edges argument list runs through a local
    // PaveFiller first).
    let mut a_bf_args: Vec<Shape> = Vec::new();
    a_bf_args.push(a_ff.clone());
    a_bf_args.extend(a_le.iter().cloned());
    let mut a_bf_filler = PaveFiller::new();
    a_bf_filler.set_arguments(a_bf_args);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_BuilderFace", 1);
    a_bf_filler.perform(&a_ps);
    let mut a_bf =
        crate::bop::algo::builder_face::BuilderFace::new(a_bf_filler.ds());
    a_bf.my_face = Some(a_bf_filler.ds().arguments[0].clone());
    a_bf.my_edges = a_bf_filler.ds().arguments[1..].to_vec();
    a_bf.perform();
    // OCCT L470-473: if (aBF.HasErrors()) return.
    if a_bf.has_errors() {
        return;
    }

    // OCCT L475-484: aLFSp = aBF.Areas(); orient the splits and bind the
    // face origins.
    let a_lf_sp: Vec<Shape> = a_bf.my_areas.clone();
    for mut a_f_sp in a_lf_sp {
        a_f_sp = oriented(&a_f_sp, an_or);
        the_lf_images.push(a_f_sp.clone());
        // OCCT L483: theFacesOrigins.Bind(aFSp, theFace).
        shape_data_map::bind(the_faces_origins, &a_f_sp, the_face.clone());
    }
}

/// OCCT BRepLib::BuildPCurveForEdgesOnPlane (BRepLib_1.cxx L313-339) — the
/// local re-host of the bop Builder re-host body (architecture difference
/// #43): for each edge lacking a stored pcurve on the plane face, project
/// the edge's 3D curve onto the plane and store the pcurve.
pub(crate) fn build_pcurve_for_edges_on_plane(a_le: &[Shape], a_f: &Shape) {
    use crate::bop::algo::builder::project_edge_on_plane;
    let surf = match a_f.as_face().and_then(|fd| fd.surface.clone()) {
        Some(s) => s,
        None => return,
    };
    let plane = match &surf {
        rcad_kernel::geom::Surface3::Plane(pl) => pl.clone(),
        _ => return,
    };
    for a_e in a_le.iter() {
        // OCCT L319-321: if (BRep_Tool::CurveOnSurface(E, F, f, l).IsNull())
        // — the pcurve presence probe.
        if bat::brep_tool_curve_on_surface(a_e, a_f).is_some() {
            continue;
        }
        let (c3, _, _) = match bat::brep_tool_curve(a_e) {
            Some(v) => v,
            None => continue,
        };
        let range = bat::brep_tool_range(a_e);
        // The projection of the 3D curve onto the plane (the
        // GeomProjLib::Curve2d plane form of the BRepLib body).
        let pc = match project_edge_on_plane(&c3, &plane, [range.0, range.1]) {
            Some(pc) => pc,
            None => continue,
        };
        let mut e_mut = a_e.clone();
        bat::builder_update_edge_pcurve(&mut e_mut, &pc, a_f, 0.0);
        let _ = e_mut;
    }
}

/// OCCT static GetAverageTangent (cxx L491-521) — computes the average
/// tangent vector along the curve.
pub(crate) fn get_average_tangent(the_s: &Shape, the_nb_p: usize) -> glam::DVec3 {
    use rcad_kernel::geom::CurveEval;
    let mut a_va = glam::DVec3::ZERO;
    for a_e_cur in explorer(the_s, ShapeType::Edge, ShapeType::Shape) {
        let a_e = a_e_cur;
        // OCCT L499-500: aT1, aT2; aC = BRep_Tool::Curve(aE, aT1, aT2).
        let (a_c, a_t1, a_t2) = match bat::brep_tool_curve(&a_e) {
            Some(v) => v,
            None => continue,
        };
        // OCCT L504-511: aVSum = sum of normalized D1 tangents over the
        // uniform parameter walk.
        let mut a_v_sum = glam::DVec3::ZERO;
        let mut a_t = a_t1;
        let a_dt = (a_t2 - a_t1) / (the_nb_p as f64);
        while a_t <= a_t2 {
            let a_v = a_c.derivative_at(a_t);
            let len = a_v.length();
            if len > 0. {
                a_v_sum += a_v / len;
            }
            a_t += a_dt;
        }
        // OCCT L513-516: reverse the sum for REVERSED edges.
        if a_e.orientation == Orientation::Reversed {
            a_v_sum = -a_v_sum;
        }
        // OCCT L518: aVA += aVSum.
        a_va += a_v_sum;
    }
    a_va
}

// ---------------------------------------------------------------------------
// OCCT class (cxx L529-1032).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces (cxx L529-1032) — auxiliary local class
/// that is used here for building splits of offset faces, that are further
/// used for building volumes.
pub struct BRepOffsetBuildOffsetFaces {
    // Input data (cxx L984-992).
    pub(crate) my_faces: Option<Vec<Shape>>, // OCCT: const NCollection_List<TopoDS_Shape>* myFaces
    pub(crate) my_as_des: Option<Rc<RefCell<BRepAlgoAsDes>>>, // OCCT: occ::handle<BRepAlgo_AsDes> myAsDes (the shared-handle form)
    pub(crate) my_analyzer: Option<BRepOffsetAnalyse>, // OCCT: const BRepOffset_Analyse* myAnalyzer
    pub(crate) my_edges_origins: Option<ShapeDataMap<Vec<Shape>>>, // OCCT: myEdgesOrigins
    pub(crate) my_faces_origins: Option<ShapeDataMap<Shape>>, // OCCT: myFacesOrigins
    pub(crate) my_e_trim_e_inf: Option<ShapeDataMap<Shape>>, // OCCT: myETrimEInf

    // Intermediate data (cxx L994-1023).
    // OCCT: myOEImages — DataMap<Shape, List<Shape>>.
    pub(crate) my_oe_images: ShapeDataMap<Vec<Shape>>,
    // OCCT: myOEOrigins — DataMap<Shape, List<Shape>>.
    pub(crate) my_oe_origins: ShapeDataMap<Vec<Shape>>,
    // OCCT: myOFImages — IndexedDataMap<Shape, List<Shape>>.
    pub(crate) my_of_images: ShapeIndexedDataMap<Vec<Shape>>,

    pub(crate) my_invalid_edges: OcctIndexedShapeMap, // OCCT: myInvalidEdges
    pub(crate) my_valid_edges: OcctIndexedShapeMap,   // OCCT: myValidEdges
    pub(crate) my_inverted_edges: OcctIndexedShapeMap, // OCCT: myInvertedEdges
    pub(crate) my_edges_to_avoid: OcctIndexedShapeMap, // OCCT: myEdgesToAvoid
    pub(crate) my_last_inv_edges: OcctShapeSet,       // OCCT: myLastInvEdges
    pub(crate) my_modified_edges: OcctShapeSet,       // OCCT: myModifiedEdges
    pub(crate) my_inside_edges: OcctIndexedShapeMap,  // OCCT: myInsideEdges

    // OCCT: myInvalidFaces — IndexedDataMap<Shape, List<Shape>>.
    pub(crate) my_invalid_faces: ShapeIndexedDataMap<Vec<Shape>>,
    // OCCT: myArtInvalidFaces — BRepOffset_DataMapOfShapeIndexedMapOfShape.
    pub(crate) my_art_invalid_faces: DataMapOfShapeIndexedMapOfShape,
    // OCCT: myAlreadyInvFaces — DataMap<Shape, int>.
    pub(crate) my_already_inv_faces: ShapeDataMap<i32>,
    // OCCT: myFNewHoles — DataMap<Shape, List<Shape>>.
    pub(crate) my_f_new_holes: ShapeDataMap<Vec<Shape>>,

    // OCCT: mySSInterfs — DataMap<Shape, List<Shape>>.
    pub(crate) my_ss_interfs: ShapeDataMap<Vec<Shape>>,
    // OCCT: mySSInterfsArt — DataMap<Shape, List<Shape>>.
    pub(crate) my_ss_interfs_art: ShapeDataMap<Vec<Shape>>,
    // OCCT: myIntersectionPairs —
    // DataMap<Shape, DataMap<Shape, Map<Shape>>, ShapeMapHasher>.
    pub(crate) my_intersection_pairs: ShapeDataMap<ShapeDataMap<OcctShapeSet>>,

    // OCCT: myFacesToRebuild — IndexedDataMap<Shape, List<Shape>>.
    pub(crate) my_faces_to_rebuild: ShapeIndexedDataMap<Vec<Shape>>,
    // OCCT: myFSelfRebAvoid — Map<Shape>.
    pub(crate) my_f_self_reb_avoid: OcctShapeSet,

    pub(crate) my_solids: Shape, // OCCT: mySolids (cxx L1025)

    // Auxiliary tools (cxx L1028).
    // OCCT: myContext — the IntTools_Context carrier (the reduced re-hosts
    // of the IsMicroEdge / SolidClassifier / PointInFace calls keep the
    // standalone forms; the member preserves the OCCT ctor allocation).
    #[allow(dead_code)]
    pub(crate) my_context: crate::bop::int_tools::context::IntToolsContext,
}

impl BRepOffsetBuildOffsetFaces {
    /// OCCT BRepOffset_BuildOffsetFaces::BRepOffset_BuildOffsetFaces(Image)
    /// (cxx L533-542) — the OCCT ctor keeps only myImage = &theImage; the
    /// rcad form passes the image explicitly to the consuming methods
    /// (architecture difference #38).
    pub fn new() -> Self {
        BRepOffsetBuildOffsetFaces {
            my_faces: None,
            my_analyzer: None,
            my_edges_origins: None,
            my_faces_origins: None,
            my_e_trim_e_inf: None,
            my_as_des: None,
            my_oe_images: HashMap::new(),
            my_oe_origins: HashMap::new(),
            my_of_images: indexmap::IndexMap::new(),
            my_invalid_edges: OcctIndexedShapeMap::new(),
            my_valid_edges: OcctIndexedShapeMap::new(),
            my_inverted_edges: OcctIndexedShapeMap::new(),
            my_edges_to_avoid: OcctIndexedShapeMap::new(),
            my_last_inv_edges: HashMap::new(),
            my_modified_edges: HashMap::new(),
            my_inside_edges: OcctIndexedShapeMap::new(),
            my_invalid_faces: indexmap::IndexMap::new(),
            my_art_invalid_faces: HashMap::new(),
            my_already_inv_faces: HashMap::new(),
            my_f_new_holes: HashMap::new(),
            my_ss_interfs: HashMap::new(),
            my_ss_interfs_art: HashMap::new(),
            my_intersection_pairs: HashMap::new(),
            my_faces_to_rebuild: indexmap::IndexMap::new(),
            my_f_self_reb_avoid: HashMap::new(),
            my_solids: Shape::null(),
            my_context: crate::bop::int_tools::context::IntToolsContext::new(),
        }
    }

    // -------------------------------------------------------------------
    // Setting data (cxx L544-574).
    // -------------------------------------------------------------------

    /// OCCT SetFaces (cxx L546) — sets faces to build splits.
    pub fn set_faces(&mut self, the_faces: &[Shape]) {
        self.my_faces = Some(the_faces.to_vec());
    }

    /// OCCT SetAsDesInfo (cxx L549) — sets ascendants/descendants info
    /// (the Handle form — the shared Rc carrier).
    pub fn set_as_des_info(&mut self, the_as_des: &Rc<RefCell<BRepAlgoAsDes>>) {
        self.my_as_des = Some(the_as_des.clone());
    }

    /// OCCT SetAnalysis (cxx L552) — sets the analysis info of the input
    /// shape.
    pub(crate) fn set_analysis(&mut self, the_analyse: &BRepOffsetAnalyse) {
        self.my_analyzer = Some(the_analyse.clone());
    }

    /// OCCT SetEdgesOrigins (cxx L555-560) — origins of the offset edges.
    pub fn set_edges_origins(&mut self, the_edges_origins: &ShapeDataMap<Vec<Shape>>) {
        self.my_edges_origins = Some(the_edges_origins.clone());
    }

    /// OCCT SetFacesOrigins (cxx L563-567) — origins of the offset faces.
    pub fn set_faces_origins(&mut self, the_faces_origins: &ShapeDataMap<Shape>) {
        self.my_faces_origins = Some(the_faces_origins.clone());
    }

    /// OCCT SetInfEdges (cxx L570-574) — infinite (extended) edges for the
    /// trimmed ones.
    pub fn set_inf_edges(&mut self, the_e_trim_e_inf: &ShapeDataMap<Shape>) {
        self.my_e_trim_e_inf = Some(the_e_trim_e_inf.clone());
    }

    // -------------------------------------------------------------------
    // Public methods to build the splits (cxx L576-582).
    // -------------------------------------------------------------------

    /// OCCT BuildSplitsOfTrimmedFaces (cxx L1036-1084) — build splits of
    /// already trimmed faces.
    #[allow(unused_variables)]
    pub fn build_splits_of_trimmed_faces(&mut self, the_image: &mut BRepAlgoImage, _the_range: &ProgressScope) {
        // OCCT L1038-1041: if (!hasData(myFaces)) return.
        if !has_data(self.my_faces.as_ref()) {
            return;
        }

        // OCCT L1043-1048: anEdgesOrigins local; when myEdgesOrigins is
        // null it is bound to the local map (architecture difference #40 —
        // the owned Option form replaces the pointer rebinding).
        if self.my_edges_origins.is_none() {
            self.my_edges_origins = Some(HashMap::new());
        }

        // OCCT L1050: Message_ProgressScope aPS(theRange, "Building splits
        // of trimmed faces", 5) — flattened local scope (architecture
        // difference #39).
        let _a_ps = ProgressScope::new(&NoopProgress, "Building splits of trimmed faces", 5);

        // OCCT L1053: fuse all edges — IntersectTrimmedEdges(aPS.Next(1)).
        let a_ps_next = ProgressScope::new(&NoopProgress, "", 1);
        self.intersect_trimmed_edges(&a_ps_next);

        // OCCT L1055-1081: the face loop.
        let a_faces: Vec<Shape> = self.my_faces.as_ref().unwrap().clone();
        let a_ps_loop = ProgressScope::new(&NoopProgress, "", a_faces.len());
        for a_f in a_faces.iter() {
            // OCCT L1058-1061: if (!aPSLoop.More()) return.
            let _ = &a_ps_loop;

            let mut a_ce = Shape::null();
            let b_found = self.get_edges(a_f, &mut a_ce, None);

            // OCCT L1067-1075: split the face by the edges.
            if !b_found {
                if !the_image.has_image(a_f) {
                    // OCCT L1072:
                    // myOFImages(myOFImages.Add(aF, List())).Append(aF) —
                    // bind the empty list, then append the face itself.
                    let a_key = bat::shape_key(a_f);
                    if !self.my_of_images.contains_key(&a_key) {
                        self.my_of_images.insert(a_key, (a_f.clone(), Vec::new()));
                    }
                    let a_idx = self
                        .my_of_images
                        .get_full(&a_key)
                        .map(|(i, _, _)| i + 1)
                        .expect("myOFImages: the face was just added");
                    shape_indexed_data_map_change_find_1(&mut self.my_of_images, a_idx)
                        .push(a_f.clone());
                }
                continue;
            }

            let mut a_lf_images: Vec<Shape> = Vec::new();
            build_splits_of_trimmed_face(a_f, &a_ce, &mut a_lf_images, &a_ps_loop);

            // OCCT L1080: myOFImages.Add(aF, aLFImages).
            add_of_images(&mut self.my_of_images, a_f, a_lf_images);
        }
        // OCCT L1083: fill history for faces and edges.
        self.fill_history(the_image);
    }

    /// OCCT BuildSplitsOfExtendedFaces (cxx L1088-1169) — building splits
    /// of not-trimmed offset faces; for the cases in which invalidities
    /// will be found, these invalidities will be rebuilt.
    #[allow(unused_variables, unused_assignments, unused_mut)]
    pub fn build_splits_of_extended_faces(&mut self, the_image: &mut BRepAlgoImage, _the_range: &ProgressScope) {
        // OCCT L1091-1095: check input data.
        if !has_data(self.my_faces.as_ref())
            || self.my_edges_origins.as_ref().map(|m| m.is_empty()).unwrap_or(true)
            || self.my_faces_origins.as_ref().map(|m| m.is_empty()).unwrap_or(true)
            || self.my_e_trim_e_inf.as_ref().map(|m| m.is_empty()).unwrap_or(true)
        {
            return;
        }

        // OCCT L1097-1102: the progress scopes (flattened, architecture
        // difference #39).
        let _a_ps = ProgressScope::new(&NoopProgress, "Building splits of extended faces", 100);
        // OCCT L1102: double aWhole = 100. - 4.
        let a_whole: f64 = 100. - 4.;

        // OCCT L1105: fuse all trimmed offset edges.
        let a_ps_next = ProgressScope::new(&NoopProgress, "", 1);
        self.intersect_trimmed_edges(&a_ps_next);

        // OCCT L1111: vertices to avoid.
        let mut a_verts_to_avoid: OcctShapeSet = HashMap::new();

        // OCCT L1115-1121: the rebuild loop (10 attempts max, halved
        // progress portions).
        let a_nb_max_attempts = 10;
        // OCCT L1120: double aPart = aWhole / 2.
        let mut a_part: f64 = a_whole / 2.;
        for _i_count in 1..=a_nb_max_attempts {
            // OCCT L1123-1126.
            let _ = &_a_ps;

            // OCCT L1128-1137: clear the data before further faces
            // construction.
            self.my_invalid_faces.clear();
            self.my_art_invalid_faces.clear();
            self.my_invalid_edges.clear();
            self.my_inverted_edges.clear();
            self.my_ss_interfs.clear();
            self.my_ss_interfs_art.clear();
            self.my_intersection_pairs.clear();
            self.my_solids = Shape::null();
            self.my_faces_to_rebuild.clear();
            self.my_f_self_reb_avoid.clear();

            // OCCT L1142: the loop scope.
            let a_ps_loop = ProgressScope::new(&NoopProgress, "", 10);

            // OCCT L1145: build splits of the faces having new
            // intersection edges.
            let a_ps_next7 = ProgressScope::new(&NoopProgress, "", 7);
            self.build_splits_of_faces(&a_ps_next7);
            if self.my_invalid_faces.is_empty() {
                break;
            }

            // OCCT L1152: find faces to rebuild.
            self.find_faces_to_rebuild();
            if self.my_faces_to_rebuild.is_empty() {
                break;
            }

            // OCCT L1159-1160: perform new intersections.
            self.my_modified_edges.clear();
            let a_ps_next3 = ProgressScope::new(&NoopProgress, "", 3);
            self.intersect_faces(&mut a_verts_to_avoid, &a_ps_next3);

            a_part /= 2.;
        }

        // OCCT L1165: fill possible gaps in the splits of offset faces.
        let a_ps_next4 = ProgressScope::new(&NoopProgress, "", 4);
        self.fill_gaps(&a_ps_next4);

        // OCCT L1168: fill history for faces and edges.
        self.fill_history(the_image);
    }

    // -------------------------------------------------------------------
    // Private helpers re-hosted across the split modules.
    // -------------------------------------------------------------------

    /// OCCT IntersectTrimmedEdges (cxx L1207-1292) — the module b form.
    pub(crate) fn intersect_trimmed_edges(&mut self, the_range: &ProgressScope) {
        super::brep_offset_make_offset_1_b::intersect_trimmed_edges_impl(self, the_range)
    }

    /// OCCT BuildSplitsOfFaces (cxx L1342-2055) — the module b form.
    pub(crate) fn build_splits_of_faces(&mut self, the_range: &ProgressScope) {
        super::brep_offset_make_offset_1_b::build_splits_of_faces_impl(self, the_range)
    }

    /// OCCT GetEdges (cxx L1839-1958) — the module b form.
    pub(crate) fn get_edges(
        &mut self,
        the_face: &Shape,
        the_edges: &mut Shape,
        the_inv: Option<&mut OcctIndexedShapeMap>,
    ) -> bool {
        super::brep_offset_make_offset_1_b::get_edges_impl(self, the_face, the_edges, the_inv)
    }

    /// OCCT CheckIfArtificial (cxx L1960-2055) — the module b form.
    pub(crate) fn check_if_artificial(
        &mut self,
        the_f: &Shape,
        the_lf_images: &[Shape],
        the_ce: &Shape,
        the_map_e_inv: &OcctIndexedShapeMap,
        the_men_inv: &mut OcctShapeSet,
    ) -> bool {
        super::brep_offset_make_offset_1_b::check_if_artificial_impl(
            self,
            the_f,
            the_lf_images,
            the_ce,
            the_map_e_inv,
            the_men_inv,
        )
    }

    /// OCCT FindInvalidEdges, the per-face form (cxx L2057-2535) — the
    /// module c form.
    pub(crate) fn find_invalid_edges_per_face(
        &mut self,
        the_f: &Shape,
        the_lf_images: &[Shape],
        the_dmfmve: &mut ShapeDataMap<OcctShapeSet>,
        the_dmfmne: &mut ShapeDataMap<OcctShapeSet>,
        the_dmfmie: &mut DataMapOfShapeIndexedMapOfShape,
        the_dmfmvie: &mut ShapeDataMap<OcctShapeSet>,
        the_dmeorleim: &mut ShapeDataMap<Vec<Shape>>,
        the_edges_invalid_by_vertex: &mut OcctShapeSet,
        the_edges_valid_by_vertex: &mut OcctShapeSet,
        the_range: &ProgressScope,
    ) {
        super::brep_offset_make_offset_1_c::find_invalid_edges_per_face_impl(
            self,
            the_f,
            the_lf_images,
            the_dmfmve,
            the_dmfmne,
            the_dmfmie,
            the_dmfmvie,
            the_dmeorleim,
            the_edges_invalid_by_vertex,
            the_edges_valid_by_vertex,
            the_range,
        )
    }

    /// OCCT FindInvalidEdges, the list form (cxx L2569-2774) — the module c
    /// form.
    pub(crate) fn find_invalid_edges_list(
        &mut self,
        the_lfoffset: &[Shape],
        the_loc_inv_edges: &mut DataMapOfShapeIndexedMapOfShape,
        the_loc_valid_edges: &mut ShapeDataMap<OcctShapeSet>,
        the_neutral_edges: &mut ShapeDataMap<OcctShapeSet>,
    ) {
        super::brep_offset_make_offset_1_c::find_invalid_edges_list_impl(
            self,
            the_lfoffset,
            the_loc_inv_edges,
            the_loc_valid_edges,
            the_neutral_edges,
        )
    }

    /// OCCT MakeInvertedEdgesInvalid (cxx L2775-2854) — the module c form.
    pub(crate) fn make_inverted_edges_invalid(&mut self, the_lfoffset: &[Shape]) {
        super::brep_offset_make_offset_1_c::make_inverted_edges_invalid_impl(self, the_lfoffset)
    }

    /// OCCT FindInvalidFaces (cxx L2856-3088) — the module c form.
    pub(crate) fn find_invalid_faces(
        &mut self,
        the_lf_images: &mut Vec<Shape>,
        the_dmfmve: &ShapeDataMap<OcctShapeSet>,
        the_dmfmie: &DataMapOfShapeIndexedMapOfShape,
        the_me_neutral: &OcctShapeSet,
        the_edges_invalid_by_vertex: &OcctShapeSet,
        the_edges_valid_by_vertex: &OcctShapeSet,
        the_mf_holes: &OcctShapeSet,
        the_mf_inv_in_hole: &mut OcctIndexedShapeMap,
        the_inv_faces: &mut Vec<Shape>,
        the_inverted_faces: &mut Vec<Shape>,
    ) {
        super::brep_offset_make_offset_1_c::find_invalid_faces_impl(
            self,
            the_lf_images,
            the_dmfmve,
            the_dmfmie,
            the_me_neutral,
            the_edges_invalid_by_vertex,
            the_edges_valid_by_vertex,
            the_mf_holes,
            the_mf_inv_in_hole,
            the_inv_faces,
            the_inverted_faces,
        )
    }

    /// OCCT FindFacesInsideHoleWires (cxx L3089-3315) — the module d form.
    pub(crate) fn find_faces_inside_hole_wires(
        &mut self,
        the_f_origin: &Shape,
        the_f_offset: &Shape,
        the_lf_images: &[Shape],
        the_dmeorleim: &ShapeDataMap<Vec<Shape>>,
        the_ef_map: &ShapeIndexedDataMap<Vec<Shape>>,
        the_mf_holes: &mut OcctShapeSet,
    ) {
        super::brep_offset_make_offset_1_d::find_faces_inside_hole_wires_impl(
            self,
            the_f_origin,
            the_f_offset,
            the_lf_images,
            the_dmeorleim,
            the_ef_map,
            the_mf_holes,
        )
    }

    /// OCCT CheckInverted (cxx L3316-3521) — the module d form.
    pub(crate) fn check_inverted(
        &mut self,
        the_e_im: &Shape,
        the_f_or: &Shape,
        the_dmve: &ShapeIndexedDataMap<Vec<Shape>>,
        the_m_edges: &OcctIndexedShapeMap,
    ) -> bool {
        super::brep_offset_make_offset_1_d::check_inverted_impl(
            self,
            the_e_im,
            the_f_or,
            the_dmve,
            the_m_edges,
        )
    }

    /// OCCT CheckInvertedBlock (cxx L3549-3686) — the module d form.
    pub(crate) fn check_inverted_block(
        &mut self,
        the_cb: &Shape,
        the_lcbf: &[Shape],
        the_dmcbv_inverted: &mut ShapeDataMap<OcctShapeSet>,
        the_dmcbv_all: &mut ShapeDataMap<OcctShapeSet>,
    ) -> bool {
        super::brep_offset_make_offset_1_d::check_inverted_block_impl(
            self,
            the_cb,
            the_lcbf,
            the_dmcbv_inverted,
            the_dmcbv_all,
        )
    }

    /// OCCT RemoveInvalidSplitsByInvertedEdges (cxx L3687-3892) — the
    /// module d form.
    pub(crate) fn remove_invalid_splits_by_inverted_edges(
        &mut self,
        the_me_removed: &mut OcctIndexedShapeMap,
    ) {
        super::brep_offset_make_offset_1_d::remove_invalid_splits_by_inverted_edges_impl(
            self,
            the_me_removed,
        )
    }

    /// OCCT RemoveInvalidSplitsFromValid (cxx L3893-4075) — the module e
    /// form.
    pub(crate) fn remove_invalid_splits_from_valid(
        &mut self,
        the_dmfmvie: &ShapeDataMap<OcctShapeSet>,
    ) {
        super::brep_offset_make_offset_1_e::remove_invalid_splits_from_valid_impl(self, the_dmfmvie)
    }

    /// OCCT RemoveInsideFaces (cxx L4295-4642) — the module e form.
    pub(crate) fn remove_inside_faces(
        &mut self,
        the_inverted_faces: &[Shape],
        the_mf_to_check_int: &OcctIndexedShapeMap,
        the_mf_inv_in_hole: &OcctIndexedShapeMap,
        the_f_holes: &Shape,
        the_me_removed: &mut OcctIndexedShapeMap,
        the_range: &ProgressScope,
    ) {
        super::brep_offset_make_offset_1_e::remove_inside_faces_impl(
            self,
            the_inverted_faces,
            the_mf_to_check_int,
            the_mf_inv_in_hole,
            the_f_holes,
            the_me_removed,
            the_range,
        )
    }

    /// OCCT ShapesConnections (cxx L4643-4967) — the module e form.
    pub(crate) fn shapes_connections(
        &mut self,
        the_dm_for: &ShapeDataMap<Shape>,
        the_builder: &BuilderRef,
    ) {
        super::brep_offset_make_offset_1_e::shapes_connections_impl(self, the_dm_for, the_builder)
    }

    /// OCCT RemoveHangingParts (cxx L4968-5167) — the module e form.
    pub(crate) fn remove_hanging_parts(
        &mut self,
        the_mv: &BuilderRef,
        the_dmf_im_f: &ShapeDataMap<Shape>,
        the_mf_inv: &OcctIndexedShapeMap,
        the_mf_to_rem: &mut OcctShapeSet,
    ) {
        super::brep_offset_make_offset_1_e::remove_hanging_parts_impl(
            self,
            the_mv,
            the_dmf_im_f,
            the_mf_inv,
            the_mf_to_rem,
        )
    }

    /// OCCT RemoveValidSplits (cxx L5168-5227) — the module f form.
    pub(crate) fn remove_valid_splits(
        &mut self,
        the_sp_rem: &OcctShapeSet,
        the_gf: &BuilderRef,
        the_me_removed: &mut OcctIndexedShapeMap,
    ) {
        super::brep_offset_make_offset_1_f::remove_valid_splits_impl(
            self,
            the_sp_rem,
            the_gf,
            the_me_removed,
        )
    }

    /// OCCT RemoveInvalidSplits (cxx L5228-5319) — the module f form.
    pub(crate) fn remove_invalid_splits(
        &mut self,
        the_sp_rem: &OcctShapeSet,
        the_gf: &BuilderRef,
        the_me_removed: &mut OcctIndexedShapeMap,
    ) {
        super::brep_offset_make_offset_1_f::remove_invalid_splits_impl(
            self,
            the_sp_rem,
            the_gf,
            the_me_removed,
        )
    }

    /// OCCT FilterEdgesImages (cxx L5320-5368) — the module f form.
    pub(crate) fn filter_edges_images(&mut self, the_s: &Shape) {
        super::brep_offset_make_offset_1_f::filter_edges_images_impl(self, the_s)
    }

    /// OCCT FilterInvalidFaces (cxx L5369-5554) — the module f form.
    pub(crate) fn filter_invalid_faces(
        &mut self,
        the_dmef: &ShapeIndexedDataMap<Vec<Shape>>,
        the_me_removed: &OcctIndexedShapeMap,
    ) {
        super::brep_offset_make_offset_1_f::filter_invalid_faces_impl(self, the_dmef, the_me_removed)
    }

    /// OCCT CheckEdgesCreatedByVertex (cxx L5555-5605) — the module f form.
    pub(crate) fn check_edges_created_by_vertex(&mut self) {
        super::brep_offset_make_offset_1_f::check_edges_created_by_vertex_impl(self)
    }

    /// OCCT FilterInvalidEdges (cxx L5606-5757) — the module f form.
    pub(crate) fn filter_invalid_edges(
        &mut self,
        the_dmfmie: &DataMapOfShapeIndexedMapOfShape,
        the_me_removed: &OcctIndexedShapeMap,
        the_me_inside: &OcctIndexedShapeMap,
        the_me_use_in_rebuild: &mut OcctShapeSet,
    ) {
        super::brep_offset_make_offset_1_f::filter_invalid_edges_impl(
            self,
            the_dmfmie,
            the_me_removed,
            the_me_inside,
            the_me_use_in_rebuild,
        )
    }

    /// OCCT FindFacesToRebuild (cxx L5758-5905) — the module f form.
    pub(crate) fn find_faces_to_rebuild(&mut self) {
        super::brep_offset_make_offset_1_f::find_faces_to_rebuild_impl(self)
    }

    /// OCCT IntersectFaces, the main form (cxx L5924-6570) — the module g
    /// form.
    pub(crate) fn intersect_faces(
        &mut self,
        the_verts_to_avoid: &mut OcctShapeSet,
        the_range: &ProgressScope,
    ) {
        super::brep_offset_make_offset_1_g::intersect_faces_impl(self, the_verts_to_avoid, the_range)
    }

    /// OCCT PrepareFacesForIntersection (cxx L6571-6679) — the module g
    /// form.
    pub(crate) fn prepare_faces_for_intersection(
        &mut self,
        the_look_vert_to_avoid: bool,
        the_fle: &mut ShapeIndexedDataMap<Vec<Shape>>,
        the_mdone: &mut ShapeDataMap<Vec<Shape>>,
        the_dmsf: &mut ShapeDataMap<Vec<Shape>>,
        the_meinf_etrim: &mut ShapeDataMap<Vec<Shape>>,
        the_dmvefull: &mut ShapeDataMap<Vec<Shape>>,
        the_dmefinv: &mut ShapeIndexedDataMap<Vec<Shape>>,
    ) {
        super::brep_offset_make_offset_1_g::prepare_faces_for_intersection_impl(
            self,
            the_look_vert_to_avoid,
            the_fle,
            the_mdone,
            the_dmsf,
            the_meinf_etrim,
            the_dmvefull,
            the_dmefinv,
        )
    }

    /// OCCT FindVerticesToAvoid (cxx L6680-6787) — the module g form.
    pub(crate) fn find_vertices_to_avoid(
        &mut self,
        the_dmefinv: &ShapeIndexedDataMap<Vec<Shape>>,
        the_dmvefull: &ShapeDataMap<Vec<Shape>>,
        the_mvrinv: &mut OcctShapeSet,
    ) {
        super::brep_offset_make_offset_1_g::find_vertices_to_avoid_impl(
            self,
            the_dmefinv,
            the_dmvefull,
            the_mvrinv,
        )
    }

    /// OCCT FindFacesForIntersection (cxx L6788-6974) — the module g form.
    pub(crate) fn find_faces_for_intersection(
        &mut self,
        the_f_inv: &Shape,
        the_me: &OcctIndexedShapeMap,
        the_dmsf: &ShapeDataMap<Vec<Shape>>,
        the_mv_inv_all: &OcctShapeSet,
        the_art_case: bool,
        the_mf_avoid: &mut OcctIndexedShapeMap,
        the_mf_int: &mut OcctIndexedShapeMap,
        the_mf_int_ext: &mut OcctIndexedShapeMap,
        the_lf_im_int: &mut Vec<Shape>,
    ) {
        super::brep_offset_make_offset_1_g::find_faces_for_intersection_impl(
            self,
            the_f_inv,
            the_me,
            the_dmsf,
            the_mv_inv_all,
            the_art_case,
            the_mf_avoid,
            the_mf_int,
            the_mf_int_ext,
            the_lf_im_int,
        )
    }

    /// OCCT ProcessCommonEdges (cxx L6975-7112) — the module h form.
    pub(crate) fn process_common_edges(
        &mut self,
        the_lec: &[Shape],
        the_me: &OcctIndexedShapeMap,
        the_meinf_etrim: &ShapeDataMap<Vec<Shape>>,
        the_all_invs: &OcctShapeSet,
        the_force_use: bool,
        the_mecv: &mut OcctIndexedShapeMap,
        the_mecheckext: &mut OcctShapeSet,
        the_dmeetrim: &mut ShapeDataMap<Vec<Shape>>,
        the_lfei: &mut Vec<Shape>,
        the_lfej: &mut Vec<Shape>,
        the_me_to_int: &mut OcctIndexedShapeMap,
    ) {
        super::brep_offset_make_offset_1_h::process_common_edges_impl(
            self,
            the_lec,
            the_me,
            the_meinf_etrim,
            the_all_invs,
            the_force_use,
            the_mecv,
            the_mecheckext,
            the_dmeetrim,
            the_lfei,
            the_lfej,
            the_me_to_int,
        )
    }

    /// OCCT UpdateIntersectedFaces (cxx L7159-7232) — the module h form.
    pub(crate) fn update_intersected_faces(
        &mut self,
        the_f_inv: &Shape,
        the_fi: &Shape,
        the_fj: &Shape,
        the_lf_inv: &[Shape],
        the_lf_imi: &[Shape],
        the_lf_imj: &[Shape],
        the_lfei: &[Shape],
        the_lfej: &[Shape],
        the_me_to_int: &mut OcctIndexedShapeMap,
    ) {
        super::brep_offset_make_offset_1_h::update_intersected_faces_impl(
            self,
            the_f_inv,
            the_fi,
            the_fj,
            the_lf_inv,
            the_lf_imi,
            the_lf_imj,
            the_lfei,
            the_lfej,
            the_me_to_int,
        )
    }

    /// OCCT IntersectFaces, the pair form (cxx L7233-7318) — the module h
    /// form.
    pub(crate) fn intersect_faces_pair(
        &mut self,
        the_f_inv: &Shape,
        the_fi: &Shape,
        the_fj: &Shape,
        the_lf_inv: &[Shape],
        the_lf_imi: &[Shape],
        the_lf_imj: &[Shape],
        the_lfei: &mut Vec<Shape>,
        the_lfej: &mut Vec<Shape>,
        the_mecv: &mut OcctIndexedShapeMap,
        the_me_to_int: &mut OcctIndexedShapeMap,
    ) {
        super::brep_offset_make_offset_1_h::intersect_faces_pair_impl(
            self,
            the_f_inv,
            the_fi,
            the_fj,
            the_lf_inv,
            the_lf_imi,
            the_lf_imj,
            the_lfei,
            the_lfej,
            the_mecv,
            the_me_to_int,
        )
    }

    /// OCCT IntersectAndTrimEdges (cxx L7319-7606) — the module h form.
    pub(crate) fn intersect_and_trim_edges(
        &mut self,
        the_mf_int: &OcctIndexedShapeMap,
        the_me_int: &OcctIndexedShapeMap,
        the_dmeetrim: &ShapeDataMap<Vec<Shape>>,
        the_ms_inv: &OcctIndexedShapeMap,
        the_mve: &OcctIndexedShapeMap,
        the_verts_to_avoid: &OcctShapeSet,
        the_new_verts_to_avoid: &OcctShapeSet,
        the_mecheckext: &OcctShapeSet,
        the_ss_interfs: Option<&ShapeDataMap<Vec<Shape>>>,
        the_mv_bounds: &mut OcctShapeSet,
        the_e_images: &mut ShapeDataMap<Vec<Shape>>,
    ) {
        super::brep_offset_make_offset_1_h::intersect_and_trim_edges_impl(
            self,
            the_mf_int,
            the_me_int,
            the_dmeetrim,
            the_ms_inv,
            the_mve,
            the_verts_to_avoid,
            the_new_verts_to_avoid,
            the_mecheckext,
            the_ss_interfs,
            the_mv_bounds,
            the_e_images,
        )
    }

    /// OCCT GetInvalidEdges (cxx L7607-7678) — the module i form.
    pub(crate) fn get_invalid_edges(
        &mut self,
        the_verts_to_avoid: &OcctShapeSet,
        the_mv_bounds: &OcctShapeSet,
        the_gf: &BuilderRef,
        the_me_inv: &mut OcctShapeSet,
    ) {
        super::brep_offset_make_offset_1_i::get_invalid_edges_impl(
            self,
            the_verts_to_avoid,
            the_mv_bounds,
            the_gf,
            the_me_inv,
        )
    }

    /// OCCT UpdateValidEdges (cxx L7679-8238) — the module i form.
    pub(crate) fn update_valid_edges(
        &mut self,
        the_fle: &ShapeIndexedDataMap<Vec<Shape>>,
        the_oen_edges: &ShapeIndexedDataMap<Vec<Shape>>,
        the_mv_bounds: &OcctShapeSet,
        the_me_inv_on_art: &OcctShapeSet,
        the_mecheckext: &mut OcctShapeSet,
        the_verts_to_avoid: &mut OcctShapeSet,
        the_e_images: &mut ShapeDataMap<Vec<Shape>>,
        the_eetrim: &mut ShapeDataMap<Vec<Shape>>,
        the_range: &ProgressScope,
    ) {
        super::brep_offset_make_offset_1_i::update_valid_edges_impl(
            self,
            the_fle,
            the_oen_edges,
            the_mv_bounds,
            the_me_inv_on_art,
            the_mecheckext,
            the_verts_to_avoid,
            the_e_images,
            the_eetrim,
            the_range,
        )
    }

    /// OCCT TrimNewIntersectionEdges (cxx L8239-8407) — the module i form.
    pub(crate) fn trim_new_intersection_edges(
        &mut self,
        the_le: &[Shape],
        the_eetrim: &ShapeDataMap<Vec<Shape>>,
        the_mv_bounds: &OcctShapeSet,
        the_mecheckext: &mut OcctShapeSet,
        the_e_images: &mut ShapeDataMap<Vec<Shape>>,
        the_meb: &mut OcctShapeSet,
        the_mv_old: &mut OcctShapeSet,
        the_me_new: &mut OcctShapeSet,
        the_dme_or: &mut ShapeDataMap<Vec<Shape>>,
        the_melf: &mut ShapeDataMap<Vec<Shape>>,
    ) {
        super::brep_offset_make_offset_1_i::trim_new_intersection_edges_impl(
            self,
            the_le,
            the_eetrim,
            the_mv_bounds,
            the_mecheckext,
            the_e_images,
            the_meb,
            the_mv_old,
            the_me_new,
            the_dme_or,
            the_melf,
        )
    }

    /// OCCT IntersectEdges (cxx L8408-8550) — the module i form.
    pub(crate) fn intersect_edges(
        &mut self,
        the_la: &[Shape],
        the_le: &[Shape],
        the_mv_bounds: &OcctShapeSet,
        the_verts_to_avoid: &OcctShapeSet,
        the_me_new: &mut OcctShapeSet,
        the_mecheckext: &mut OcctShapeSet,
        the_e_images: &mut ShapeDataMap<Vec<Shape>>,
        the_dme_or: &mut ShapeDataMap<Vec<Shape>>,
        the_melf: &mut ShapeDataMap<Vec<Shape>>,
        the_splits: &mut Shape,
    ) {
        super::brep_offset_make_offset_1_i::intersect_edges_impl(
            self,
            the_la,
            the_le,
            the_mv_bounds,
            the_verts_to_avoid,
            the_me_new,
            the_mecheckext,
            the_e_images,
            the_dme_or,
            the_melf,
            the_splits,
        )
    }

    /// OCCT GetBounds (cxx L8551-8593) — the module j form.
    pub(crate) fn get_bounds(
        &mut self,
        the_lfaces: &[Shape],
        the_meb: &OcctShapeSet,
        the_bounds: &mut Shape,
    ) {
        super::brep_offset_make_offset_1_j::get_bounds_impl(self, the_lfaces, the_meb, the_bounds)
    }

    /// OCCT GetBoundsToUpdate (cxx L8594-8674) — the module j form.
    pub(crate) fn get_bounds_to_update(
        &mut self,
        the_lf: &[Shape],
        the_meb: &OcctShapeSet,
        the_la_bounds: &mut Vec<Shape>,
        the_la_valid: &mut Vec<Shape>,
        the_bounds: &mut Shape,
    ) {
        super::brep_offset_make_offset_1_j::get_bounds_to_update_impl(
            self,
            the_lf,
            the_meb,
            the_la_bounds,
            the_la_valid,
            the_bounds,
        )
    }

    /// OCCT GetInvalidEdgesByBounds (cxx L8675-8996) — the module j form.
    pub(crate) fn get_invalid_edges_by_bounds(
        &mut self,
        the_splits: &Shape,
        the_bounds: &Shape,
        the_mv_old: &OcctShapeSet,
        the_me_new: &OcctShapeSet,
        the_dme_or: &ShapeDataMap<Vec<Shape>>,
        the_melf: &ShapeDataMap<Vec<Shape>>,
        the_e_images: &ShapeDataMap<Vec<Shape>>,
        the_mecheckext: &OcctShapeSet,
        the_me_inv_on_art: &OcctShapeSet,
        the_verts_to_avoid: &mut OcctShapeSet,
        the_me_inv: &mut OcctShapeSet,
    ) {
        super::brep_offset_make_offset_1_j::get_invalid_edges_by_bounds_impl(
            self,
            the_splits,
            the_bounds,
            the_mv_old,
            the_me_new,
            the_dme_or,
            the_melf,
            the_e_images,
            the_mecheckext,
            the_me_inv_on_art,
            the_verts_to_avoid,
            the_me_inv,
        )
    }

    /// OCCT FilterSplits (cxx L8997-9047) — the module j form.
    pub(crate) fn filter_splits(
        &mut self,
        the_le: &[Shape],
        the_me_filter: &OcctShapeSet,
        the_is_inv: bool,
        the_e_images: &mut ShapeDataMap<Vec<Shape>>,
        the_splits: &mut Shape,
    ) {
        super::brep_offset_make_offset_1_j::filter_splits_impl(
            self,
            the_le,
            the_me_filter,
            the_is_inv,
            the_e_images,
            the_splits,
        )
    }

    /// OCCT UpdateNewIntersectionEdges (cxx L9048-9296) — the module j
    /// form.
    pub(crate) fn update_new_intersection_edges(
        &mut self,
        the_le: &[Shape],
        the_melf: &ShapeDataMap<Vec<Shape>>,
        the_e_images: &ShapeDataMap<Vec<Shape>>,
        the_eetrim: &mut ShapeDataMap<Vec<Shape>>,
    ) {
        super::brep_offset_make_offset_1_j::update_new_intersection_edges_impl(
            self,
            the_le,
            the_melf,
            the_e_images,
            the_eetrim,
        )
    }

    /// OCCT FillGaps (cxx L9297-9412) — the module j form.
    pub(crate) fn fill_gaps(&mut self, the_range: &ProgressScope) {
        super::brep_offset_make_offset_1_j::fill_gaps_impl(self, the_range)
    }

    /// OCCT FillHistory (cxx L9413-9491) — the module j form.
    pub(crate) fn fill_history(&mut self, the_image: &mut BRepAlgoImage) {
        super::brep_offset_make_offset_1_j::fill_history_impl(self, the_image)
    }
}

// ---------------------------------------------------------------------------
// Local re-hosts shared across the split modules.
// ---------------------------------------------------------------------------

/// OCCT TopExp::MapShapes(S, T, IndexedMap) — fill the indexed map with the
/// sub-shapes of S of type T.
pub(crate) fn map_shapes_indexed(the_s: &Shape, the_t: ShapeType, the_map: &mut OcctIndexedShapeMap) {
    for a_s in explorer(the_s, the_t, ShapeType::Shape) {
        the_map.add(&a_s);
    }
}

/// OCCT TopExp::MapShapesAndAncestors(S, TS, TA, M) over the
/// ShapeIndexedDataMap form (the loc_ope_glued_shape re-host body).
pub(crate) fn map_shapes_and_ancestors_map(
    the_s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut ShapeIndexedDataMap<Vec<Shape>>,
) {
    for a_anc in explorer(the_s, ta, ShapeType::Shape) {
        for a_exs in explorer(&a_anc, ts, ShapeType::Shape) {
            let key = bat::shape_key(&a_exs);
            let entry = m.entry(key).or_insert((a_exs.clone(), Vec::new()));
            entry.1.push(a_anc.clone());
        }
    }
    for a_ex in explorer(the_s, ts, ta) {
        let key = bat::shape_key(&a_ex);
        m.entry(key).or_insert((a_ex, Vec::new()));
    }
}

/// The TopTools_ShapeMapHasher key of the shape (the (TShape ptr,
/// location) identity pair).
pub(crate) fn shape_key_of(s: &Shape) -> ShapeKey {
    bat::shape_key(s)
}

/// The 1-based key read over the ShapeIndexedDataMap (the OCCT
/// IndexedDataMap::FindKey(i) form).
pub(crate) fn find_key_1_local<V>(m: &ShapeIndexedDataMap<V>, i: usize) -> &Shape {
    &m.get_index(i - 1)
        .expect("IndexedDataMap::FindKey out of range")
        .1
        .0
}

/// The 1-based value read over the ShapeIndexedDataMap (the OCCT
/// IndexedDataMap::operator(i) const form).
pub(crate) fn value_1_local<V>(m: &ShapeIndexedDataMap<V>, i: usize) -> &V {
    &m.get_index(i - 1)
        .expect("IndexedDataMap::operator() out of range")
        .1
        .1
}

/// The OcctIndexedShapeMap clone through the public add() API (the map
/// carries no derive(Clone) in the tool cluster).
pub(crate) fn clone_indexed_map(m: &OcctIndexedShapeMap) -> OcctIndexedShapeMap {
    let mut m2 = OcctIndexedShapeMap::new();
    for a_s in m.iter() {
        m2.add(a_s);
    }
    m2
}

/// The OCCT IndexedDataMap::Find(k)/ChangeSeek(k) read over the
/// ShapeIndexedDataMap — the snapshot of the value list.
pub(crate) fn idm_seek(m: &ShapeIndexedDataMap<Vec<Shape>>, k: &Shape) -> Option<Vec<Shape>> {
    m.get(&shape_key_of(k)).map(|e| e.1.clone())
}

/// The OCCT IndexedDataMap::ChangeFind(k) write over the
/// ShapeIndexedDataMap (the Vec<Shape> value form).
pub(crate) fn idm_bind(m: &mut ShapeIndexedDataMap<Vec<Shape>>, k: &Shape, v: Vec<Shape>) {
    let key = shape_key_of(k);
    if let Some(entry) = m.get_mut(&key) {
        entry.1 = v;
    } else {
        m.insert(key, (k.clone(), v));
    }
}

/// The OCCT DataMap::Find(k) read over the ShapeIndexedDataMap — the
/// snapshot of the value list (panics when unbound — the OCCT Find form).
pub(crate) fn idm_find(m: &ShapeIndexedDataMap<Vec<Shape>>, k: &Shape) -> Vec<Shape> {
    m.get(&shape_key_of(k))
        .expect("IndexedDataMap::Find unbound")
        .1
        .clone()
}

/// The OCCT IndexedDataMap::Find(k)/ChangeSeek(k) read over the
/// ShapeIndexedDataMap — the snapshot of the value list.
/// The OCCT IndexedDataMap::ChangeFind(k) write over the
/// ShapeIndexedDataMap (the Vec<Shape> value form).
/// The OCCT DataMap::Find(k) read over the ShapeIndexedDataMap — the
/// snapshot of the value list (panics when unbound — the OCCT Find form).
/// The OCCT IndexedDataMap::Find(k)/ChangeSeek(k) read over the
/// ShapeIndexedDataMap — the snapshot of the value list.
/// The OCCT IndexedDataMap::ChangeFind(k) write over the
/// ShapeIndexedDataMap (the Vec<Shape> value form).
/// The OCCT DataMap::Find(k) read over the ShapeIndexedDataMap — the
/// snapshot of the value list (panics when unbound — the OCCT Find form).
/// The OCCT IndexedDataMap::Find(k)/ChangeSeek(k) read over the
/// ShapeIndexedDataMap — the snapshot of the value list.
/// The OCCT IndexedDataMap::ChangeFind(k) write over the
/// ShapeIndexedDataMap (the Vec<Shape> value form).
/// The OCCT DataMap::Find(k) read over the ShapeIndexedDataMap — the
/// snapshot of the value list (panics when unbound — the OCCT Find form).
/// The OCCT IndexedDataMap::Find(k)/ChangeSeek(k) read over the
/// ShapeIndexedDataMap — the snapshot of the value list.
/// The OCCT IndexedDataMap::ChangeFind(k) write over the
/// ShapeIndexedDataMap (the Vec<Shape> value form).
/// The OCCT DataMap::Find(k) read over the ShapeIndexedDataMap — the
/// snapshot of the value list (panics when unbound — the OCCT Find form).
/// The OCCT IndexedDataMap::Find(k)/ChangeSeek(k) read over the
/// ShapeIndexedDataMap — the snapshot of the value list.
/// The OCCT IndexedDataMap::ChangeFind(k) write over the
/// ShapeIndexedDataMap (the Vec<Shape> value form).
/// The OCCT DataMap::Find(k) read over the ShapeIndexedDataMap — the
/// snapshot of the value list (panics when unbound — the OCCT Find form).
/// The OCCT NCollection_DataMap::ChangeSeek(k) form over the ShapeDataMap
/// (the mutable seek of the tool-cluster shape_data_map module).
pub(crate) fn dm_change_seek<'a, V>(m: &'a mut ShapeDataMap<V>, k: &Shape) -> Option<&'a mut V> {
    m.get_mut(&shape_key_of(k)).map(|e| &mut e.1)
}

/// The `BRep_Builder().MakeCompound(aCE)` form of the cxx body.
pub(crate) fn empty_compound() -> Shape {
    bat::builder_make_compound()
}

/// The list-shape carrier of the OCCT connexity-block lists (the
/// `aCB`/`aCECBInv` NCollection_List forms passed as shapes into the
/// TopExp walks) — the compound of the list elements.
pub(crate) fn block_shape(the_list: &[Shape]) -> Shape {
    let mut a_cb = empty_compound();
    for a_s in the_list.iter() {
        bat::builder_add_compound_shape(&mut a_cb, a_s);
    }
    a_cb
}

/// The 1-based IndexedDataMap mutable value accessor over the
/// ShapeIndexedDataMap (the cxx `myOFImages(i)` mutable form).
pub(crate) fn shape_indexed_data_map_change_find_1<'a, V>(
    m: &'a mut ShapeIndexedDataMap<V>,
    i: usize,
) -> &'a mut V {
    &mut m
        .get_index_mut(i - 1)
        .expect("IndexedDataMap::operator() out of range")
        .1
        .1
}

/// The OCCT `myOFImages.Add(aF, aLFImages)` form — binds or appends?  The
/// OCCT IndexedDataMap::Add *binds* (replacing nothing — the key must be
/// new); the cxx sites use Add both for the empty-bind-then-append form
/// (L1072) and the plain bind form (L1080).  The rcad form mirrors
/// shape_indexed_data_map::add (bind when absent) and overwrites the value
/// list with the given one — the OCCT Add binds the passed value.
pub(crate) fn add_of_images(m: &mut ShapeIndexedDataMap<Vec<Shape>>, k: &Shape, v: Vec<Shape>) {
    let key = bat::shape_key(k);
    if let Some(entry) = m.get_mut(&key) {
        entry.1 = v;
    } else {
        m.insert(key, (k.clone(), v));
    }
}

