// IntPatch_ImpPrmIntersection — analytic-parametric surface intersection
// (the walking algorithm).
//
// OCCT IntPatch_ImpPrmIntersection.cxx (3891 lines), ported into rcad-algo.
// The analytic surface (plane/cylinder/sphere/cone) is treated as an implicit
// quadric; the parametric surface is walked along the curve F(u,v) = 0 using
// IntPatch_TheSurfFunction / IntPatch_TheSearchInside / IntPatch_TheIWalking.

pub mod function_set_root;
pub mod i_walking;
pub mod imp_prm_intersection;
pub mod path_point;
pub mod search_inside;
pub mod surf_function;

pub use imp_prm_intersection::ImpPrmIntersection;
