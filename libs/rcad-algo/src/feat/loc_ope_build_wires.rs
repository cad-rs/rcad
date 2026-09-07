// OCCT LocOpe_BuildWires.hxx L28-48 + LocOpe_BuildWires.cxx L37-281 —
// DEFERRED TRANSLATION (Stage 3d batch).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_BuildWires.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_BuildWires.cxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Why deferred: the single entry point
//   Perform(L, PW)                      (cxx L63-228)
// takes the handle to LocOpe_WiresOnShape and drives its
// OnVertex(vtx, aV_border) / OnEdge(vtx, etmp, partmp) predicates
// (cxx L102-104). LocOpe_WiresOnShape (1623 lines: BRepAdaptor_Curve2d,
// Extrema_ExtCC/ExtCC2d, ShapeAnalysis_Edge, ShapeConstruct_ProjectCurveOn-
// Surface, BRepTopAdaptor_FClass2d, ...) is a later batch of the port plan,
// and per the port discipline the parameter type cannot be invented here.
//
// Translation anchors for the next batch (already read; cxx line ranges):
// - LocOpe_BuildWires()                        cxx L44-47
// - LocOpe_BuildWires(L, PW)                   cxx L52-56
// - Perform(L, PW)                             cxx L63-228
//   (compound C + theMap Add loop L70-84; MapShapesAndAncestors
//    VERTEX/EDGE L86-88; Bords via PW predicates L90-108; the wire walk
//    FindFirstEdge loop L110-162; the closed-loop reopen L164-194; the open
//    wire tail L195-203; myRes.Append + Bords/edge consumption L205-224)
// - IsDone()                                   cxx L232-235
// - Result()                                   cxx L239-246
// - static FindFirstEdge(theMapVE, theBord)    cxx L250-281
//
// first consumer: BRepFeat_Form family (3b) through LocOpe_Spliter::Perform
// (which builds the wires of the WiresOnShape edge list).
