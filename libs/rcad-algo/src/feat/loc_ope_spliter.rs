// OCCT LocOpe_Spliter.hxx L30-73 + LocOpe_Spliter.cxx L17-718 — DEFERRED
// TRANSLATION (Stage 3d batch).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Spliter.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Spliter.cxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Why deferred: Perform(PW) (cxx L63-718) depends on four pieces outside
// this batch:
// 1. LocOpe_WiresOnShape (1623 lines) — the PW handle and its
//    InitEdgeIterator/MoreEdge/NextEdge/Edge/OnVertex/OnEdge/Bind/BindAll
//    surface;
// 2. LocOpe_BuildWires (this package — deferred, see loc_ope_build_wires.rs)
//    used at cxx L25/L145;
// 3. BRepTools_Substitution (cxx L85/L117-118) — the TopTools substitution
//    map driving the vertex/edge replacement;
// 4. GeomAPI_ProjectPointOnCurve (cxx L23) for the EdgeOnEdg bookkeeping.
// Per the port discipline none of these may be invented here.
//
// Translation anchors for the next batch (already read; cxx line ranges):
// - LocOpe_Spliter::Init(S)            lxx (inline; myShape = S)
// - Perform(PW)                        cxx L63-718
//   (vertex substitution pass L81-115; theSubs.Build L117; edge/face
//    split passes L120-...; RebuildWires/Put/Select statics L52-59 and
//    their bodies; DirectLeft/Left assembly; DescendantShapes map)
// - ResultingShape() / Shape() / DirectLeft() / Left() /
//   DescendantShapes(S) — accessors over myRes/myDLeft/myLeft/myMap.
//
// first consumer: BRepFeat_Form family (3b) — LocOpe_Gluer::Perform
// (cxx L203-214) runs the splitter over the WiresOnShape bindings.
