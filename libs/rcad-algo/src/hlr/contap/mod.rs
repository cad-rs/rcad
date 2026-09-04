// OCCT Contap (TKHLR) — the contour (apparent outline) engine for a
// single surface.
//
// 1:1 translations, one OCCT class per file:
//   - [`surface_adaptor`]  Adaptor3d_Surface instance interface + the
//     GeomAdaptor_Surface-equivalent concrete wrapper over the kernel
//     surface (the `occ::handle<Adaptor3d_Surface>` of OCCT).
//   - [`i_type`]           Contap_IType enum.
//   - [`t_function`]       Contap_TFunction enum.
//   - [`surf_props`]       Contap_SurfProps (normals + derivatives of the
//     normal on analytic surfaces).
//   - [`point`]            Contap_Point.
//   - [`line`]             Contap_Line.
//   - [`h_curve2d_tool`]   Contap_HCurve2dTool (statics over the 2D arc).
//   - [`h_cont_tool`]      Contap_HContTool (surface/arc sampling + the
//     EPCOfExtPC2d projection).
//   - [`surf_function`]    Contap_SurfFunction (F(u,v) = N . Dir).
//   - [`arc_function`]     Contap_ArcFunction (F(t) on a restriction arc).
//   - [`cont_ana`]         Contap_ContAna (analytic contours).
//   - [`the_path_point_of_the_search`] Contap_ThePathPointOfTheSearch.
//   - [`the_segment_of_the_search`]    Contap_TheSegmentOfTheSearch.
//   - [`domain`]           the TopolTool interface used by the engines +
//     its impl over the Adaptor3d TopolTool.
//   - [`the_search_inside`] Contap_TheSearchInside (IntStart_SearchInside.gxx).
//   - [`the_search`]       Contap_TheSearch (IntStart_SearchOnBoundaries.gxx).
//   - [`contour`]          Contap_Contour (the 2389-line main engine).
//
// The Contap_TheIWalking instantiation reuses the landed IntWalk_IWalking
// engine (geomalgo::int_patch::imp_prm::i_walking) via the SurfFunction
// interface bridge; see [`surf_function`].

pub mod surface_adaptor;
pub mod geom_tool;

pub use surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter, SurfaceHandle};
pub use geom_tool::{GeomBasisCurve, GeomBasisSurface, GeomTool};

pub mod i_type;
pub mod t_function;
pub mod surf_props;
pub mod point;
pub mod line;
pub mod h_curve2d_tool;
pub mod h_cont_tool;
pub mod cont_ana;
pub mod surf_function;
pub mod arc_function;
pub mod the_path_point_of_the_search;
pub mod the_segment_of_the_search;
pub mod domain;
pub mod the_search_inside;
pub mod the_search;
pub mod contour;

pub use i_type::IType;
pub use t_function::TFunction;
pub use surf_props::{deriv_and_norm, norm_and_dn, normale};
pub use point::Point;
pub use line::Line;
pub use cont_ana::ContAna;
pub use surf_function::SurfFunction;
pub use arc_function::ArcFunction;
pub use the_path_point_of_the_search::ThePathPointOfTheSearch;
pub use the_segment_of_the_search::TheSegmentOfTheSearch;
pub use domain::ContapDomain;
pub use the_search_inside::TheSearchInside;
pub use the_search::TheSearch;
pub use contour::Contour;
pub use h_cont_tool as hcont;
