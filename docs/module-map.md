# rcad ↔ OCCT 模块映射（Module Map）

> 目的：给出 rcad 各 crate/目录与 OCCT toolkit/package 的一一对应关系，
> 供后续 agent 按图补码。补码规则以仓库根 `AGENTS.md` 的方法论为准
> （1:1 形式翻译、逐行对照 OCCT 源码、每翻译一个函数编译一次）。
>
> OCCT 源码根：`$OCCT_SRC`（默认 `C:\tools\opencascade-8.0.0-vc14-64\opencascade-8.0.0`）。
> 各 toolkit 的源码路径：`$OCCT_SRC/src/ModelingAlgorithms/<TK>/<Package>/`、
> `$OCCT_SRC/src/ModelingData/<TK>/<Package>/`、`$OCCT_SRC/src/FoundationClasses/<TK>/<Package>/`。
>
> 状态标记：
> - ✅ 已对齐（1:1 翻译完成，有 OCCT 布尔网格测试基线覆盖）
> - ◐ 部分完成（有代码，尚未声明全链路对齐）
> - ⬜ 空占位（目录已建、mod.rs 为空、未在 lib.rs 声明——等价于"待开工"）
> - ◻ 未建目录（OCCT 有、rcad 尚无对应目录）
> - — 不映射（rcad 采用自有 Rust 实现，不逐行对照 OCCT）

## 1. crate 级映射

| rcad crate | OCCT toolkit | 说明 |
|---|---|---|
| `rcad-kernel` | TKMath + TKG2d + TKG3d + TKGeomBase + TKBRep(基础) | 几何内核：`src/geom`(Geom/Geom2d 曲线曲面)、`src/base`(TKGeomBase 各包)、`src/topo`(TopoDS/BRep/TopExp)、`src/math`(TKMath: math_*、Bnd)、`src/core`(Precision/gp 常量) |
| `rcad-brep` | TKBRep | `adaptor.rs`↔BRepAdaptor、`graph`↔BRepGraph、`lprop.rs`↔BRepLProp、`tools`↔BRepTools |
| `rcad-algo` | ModelingAlgorithms 各算法 toolkit | 详见第 2、3 节 |
| `rcad-modeling` | TKPrim + BRepBuilderAPI 级构造 | `make_cylinder/cone/box/sphere`↔BRepPrimAPI、`prism_face_solid_brep`↔BRepSweep 平移棱柱 |
| `rcad-step` / `rcad-iges` | TKDESTEP / TKDEIGES | 数据交换 |
| `rcad-scene` / `rcad-render` / `rcad-xmesh 相关显示` | — | 自有实现，不映射 |

## 2. rcad-algo 现有模块 ↔ OCCT

### 2.1 `bop/` ↔ TKBO（布尔管线，✅ 已对齐，测试基线 = 四个 boolean 网格）

| rcad | OCCT package | 状态 |
|---|---|---|
| `bop/algo/`（pave_filler*.rs, builder*.rs, wire_splitter.rs, shell_splitter.rs, checker_si.rs, argument_analyzer.rs, section_attribute.rs） | `BOPAlgo` | ✅ |
| `bop/ds/`（mod.rs, common_block.rs, face_info.rs, pave.rs, iterator.rs, tools.rs, topods_builder.rs） | `BOPDS` | ✅ |
| `bop/int_tools/`（edge_edge.rs, edge_face.rs, face_face.rs, bean_face_intersector.rs, common_prt.rs, pnt_on_face.rs, …） | `BOPInt` + `IntTools` | ✅ |
| `bop/tools/`（algo_tools.rs, box_tree.rs, bvh_tree.rs） | `BOPTools`（box_tree↔BOPTools_BoxTree；bvh_tree 为 BVH 求交加速的等价实现） | ✅ |
| `bop/brep_algo_api/` | `BRepAlgoAPI`（fuse/common/cut/cut21 即 BRepAlgoAPI_Fuse/Common/Cut/Cut21） | ✅ |
| `bop/history.rs` | `BRepTools_History` / `BOPAlgo_History`（ShapeHistory） | ◐ |

注意：`face_make_curve.rs` 在 `bop/int_tools/`，它是 `IntTools_FaceFace::MakeCurve` +
`GeomInt_LineConstructor` + `GeomInt_IntSS::BuildPCurves` 的联合翻译（布尔馈送层）。

### 2.2 `geomalgo/` ↔ TKGeomAlgo（交线/交点几何算法）

| rcad | OCCT package / 文件 | 状态 |
|---|---|---|
| `geomalgo/int_patch/`（intersection.rs, imp_imp_intersection.rs, int_cycy.rs, int_quad_quad.rs, quad_quad_geo.rs, int_xx.rs, int_cs.rs, curve_surface.rs, restriction.rs, transitions.rs, w_line_tool.rs, a_line_to_w_line.rs, …） | `IntPatch`（IntPatch_ImpImpIntersection、IntPatch_ImpPrmIntersection、IntPatch_RLine/WLine/GLine、IntPatch_PointLine 等） | ✅（布尔网格覆盖） |
| `geomalgo/int_patch/imp_prm/imp_prm_intersection.rs` | `IntPatch_ImpPrmIntersection.cxx`（Perform/ComputeTangency/Recadre/IsCoincide/DecomposeResult/GetLocalStep） | ✅ 审查过（L224-465 compute_tangency 逐行核对） |
| `geomalgo/int_patch/imp_prm/i_walking.rs` | `IntWalk_IWalking.gxx` + `IntPatch_TheIWLineOfTheIWalking_0.cxx` | ✅（86236f99 审查对齐，2c49fa1c 已切回） |
| `geomalgo/int_patch/imp_prm/surf_function.rs` | `IntImp_ZerImpFunc.gxx`（=IntPatch_TheSurfFunction，经 IntSurf_QuadricTool 委托 IntSurf_Quadric） | ✅ 逐行核对无偏差 |
| `geomalgo/int_patch/imp_prm/search_inside.rs` | `IntPatch_TheSearchInside_0.cxx` → `IntStart_SearchInside` | ✅（行走馈送） |
| `geomalgo/int_patch/imp_prm/function_set_root.rs` | `math_FunctionSetRoot`（TKMath/math）+ IntWalk 用法 | ◐ |
| `geomalgo/int_patch/imp_prm/path_point.rs` | `IntSurf_PathPoint` / `IntPatch_ThePathPointOfTheSOnBounds` | ✅ |
| `geomalgo/int_patch/so_on_bounds.rs` | `IntPatch_TheSOnBounds` → `IntStart_SOnBounds`/`IntStart_Segment` | ✅ |
| `geomalgo/int_surf/quadric.rs` | `IntSurf_Quadric`（Value/Gradient/ValAndGrad 逐面类型公式） | ✅ 逐行核对 |
| `geomalgo/int_surf/line_on_2s.rs` | `IntSurf_LineOn2S` | ✅ |
| `geomalgo/int_res2d/` | `IntRes2d` | ◐ |
| `geomalgo/top_trans/` | `TopTrans`（CurveTransition） | ✅ |
| `geomalgo/approx_int.rs` | `ApproxInt`（ApproxInt_Approx、ApproxInt_KnotTools、MultiBSpCurve）+ `AppParCurves`(TKGeomBase) | ✅（WLine 近似路径覆盖） |

### 2.3 `topalgo/` ↔ TKTopAlgo（拓扑算法层）

| rcad | OCCT package | 状态 |
|---|---|---|
| `topalgo/brep_bnd_lib/` | `BRepBndLib` | ◐ |
| `topalgo/brep_class/`（face_classifier.rs, edge.rs, intersector.rs, g_inter.rs, bnd_box2d.rs, face_explorer.rs） | `BRepClass`（FClassifier/FaceExplorer/Edge）+ `Bnd`(Bnd_Box2d) | ◐ |
| `topalgo/brep_class3d/`（mod.rs, intersector3d.rs, passive_classifier.rs, bnd_box_tree.rs） | `BRepClass3d` | ◐ |
| `topalgo/brep_extrema/dist_shape_shape.rs` | `BRepExtrema`（DistShapeShape） | ◐ |
| `topalgo/brep_int_curve_surface/inter.rs` | `BRepIntCurveSurface`（Inter） | ◐ |
| `topalgo/brep_lib/brep_lib.rs` | `BRepLib`（SameParameter/SameRange/MakeEdge/BoundingVertex…） | ✅（布尔链路使用） |
| `topalgo/brep_top_adaptor/`（fclass2d.rs, class2d.rs） | `BRepTopAdaptor`（FClass2d/Class2d） | ✅（WireSplitter/分类使用） |
| `topalgo/gcpnts/` | `GCPnts`（注意：TKGeomBase 也有 GCPnts，rcad-kernel/src/base/gcpnts 为另一份） | ◐ |
| `topalgo/shape_source.rs` | `TopExp_Explorer`/TopoDS_Iterator 的遍历辅助（自有薄封装） | — |

### 2.4 `algo_ext/`（混合层：OCCT 对应 + 自有扩展）

| rcad | OCCT 对应 | 状态 |
|---|---|---|
| `algo_ext/brep_check/` | `BRepCheck`（TKTopAlgo） | ◐ |
| `algo_ext/healing/`、`algo_ext/shape_analysis/` | `ShapeAnalysis`/`ShapeFix`（TKShHealing）早期子集 | ◐ |
| `algo_ext/shape_custom.rs` | `ShapeCustom`（TKShHealing） | ◐ |
| `algo_ext/brep_repair/` | 自有修复链（无 1:1 OCCT 对应；融合了 ShapeFix 思路） | — |
| `algo_ext/bool_ops_ext.rs`、`brep_algo.rs`、`brep_tools.rs`、`topods_ext.rs`、`tolerance.rs`、`geom_populate.rs`、`bspline_edit.rs`、`features.rs`、`revolve.rs`、`fillet.rs`、`extrude_profile.rs` | 自有扩展/遗留 API（生成测试的兼容层） | — |

## 3. 新增空占位目录 ↔ OCCT（待开工，即本次新增）

| rcad 目录 | OCCT toolkit | OCCT packages（源码路径 `$OCCT_SRC/src/ModelingAlgorithms/<TK>/<Package>/`） | 状态 |
|---|---|---|---|
| `feat/` | TKFeat | `BRepFeat`、`LocOpe` | ⬜ 空占位 |
| `fillet/` | TKFillet | `BRepFilletAPI`、`ChFi2d`、`ChFi3d`、`ChFiDS`、`ChFiKPart`、`BRepBlend`、`Blend`、`BlendFunc`、`FilletSurf` | ⬜ 空占位 |
| `helix/` | TKHelix | `HelixBRep`、`HelixGeom` | ⬜ 空占位 |
| `hlr/` | TKHLR | `HLRBRep`、`HLRAlgo`、`HLRAppli`、`HLRTopoBRep`、`Contap`、`Intrv`、`TopBas`、`TopCnx` | ⬜ 空占位 |
| `offset/` | TKOffset | `BRepOffset`、`BRepOffsetAPI`、`BiTgte`、`Draft` | ⬜ 空占位 |
| `shhealing/` | TKShHealing | `ShapeFix`、`ShapeAnalysis`、`ShapeBuild`、`ShapeExtend`、`ShapeConstruct`、`ShapeCustom`、`ShapeUpgrade`、`ShapeProcess`、`ShapeProcessAPI`、`ShapeAlgo`、`SHMessage` | ⬜ 空占位 |
| `xmesh/` | TKXMesh | `XBRepMesh` | ⬜ 空占位 |

尚未建目录的 OCCT toolkit（需要时再建，勿混入上述目录）：

| OCCT toolkit | packages | 备注 |
|---|---|---|
| TKBool | `BRepAlgo`、`BRepFill`、`BRepProj`、`TopOpeBRep`、`TopOpeBRepBuild`、`TopOpeBRepDS`、`TopOpeBRepTool` | TopOpeBRep* 为旧版布尔（rcad 用 TKBO 路线，勿移植）；`BRepFill`（sweep/section）可按需建 `boolfill/` 或并入 feat |
| TKMesh | `BRepMesh`、`IMeshData`、`IMeshTools`、`BRepMeshData` | rcad 三角剖分走自有 mesh（见 AGENTS.md「后续重构指导」），不 1:1 移植 |

## 4. 命名与翻译约定（摘要，全文见根 AGENTS.md）

1. OCCT package（PascalCase）→ rcad 目录（snake_case）：`BRepFilletAPI` → `brep_fillet_api/`、`ChFi3d` → `chfi3d/`、`HLRBRep` → `hlr_brep/`（hlr 已有外层 `hlr/` 时为 `hlr/brep/` 等，按包名逐个建子目录，勿合并多个 package 进一个文件）。
2. 类 → struct、成员函数 → impl 方法、参数名保留 OCCT 命名转 snake_case（`my`/`the`/`a`/`an` 前缀照搬）。
3. 每个文件头部标注 OCCT 源文件与行号范围：`// OCCT ChFi3d_Builder.cxx L100-200`。
4. `.gxx`/`.lxx`/`.hxx` 内联实现一并翻译；`#define` 别名（如 TheIWFunction=IntPatch_TheSurfFunction）展开为注释注明。
5. 每翻译一个函数 `cargo check -p rcad-algo`；对齐期间不跑测试看数值（方法论：阶段1 形式对齐 → 阶段2 调试）。
6. 新目录需在 `libs/rcad-algo/src/lib.rs` 声明后才能编译（当前 7 个占位目录均未声明，属正常状态）。

## 5. 验证基线（补码后必须保持）

- 生成测试（gitignored，须逐网格重生成）：
  `cargo run -p occt-test-gen -- --batch-boolean --merge-groups --batch-grid <grid>`
- 四网格当前基线（全部通过，作为补码的回归门槛）：
  - bopfuse 748/0、bopcommon 755/0、boptuc 745/0、bcut 729/0（含 g6）
  - `cargo test -p rcad-algo --lib` 70/0
- bcut 全量跑法：`cargo test -p occt-generated-tests --test generated_occt_boolean_bopcut_simple`（g6 已可直接跑，无需 --skip）。
- 参考 STEP/拓扑 JSON：`tests/occt/step_output/ref/`、`tests/occt/step_reference/`；重新生成单个用例的 JSON 可复用 `tools/gen-occt-ref/gen_ref_topology.py` 的 `get_ref_topology(grid, case)`。
