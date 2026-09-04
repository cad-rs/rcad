# TKHLR 1:1 翻译推进计划（多 session runway）

> 本文档是 TKHLR（Hidden Line Removal）移植的持久 runway 文档。每个 session
> 开场必读：先看 §7 进度勾选表，从第一个未勾选任务继续；完成一项勾一项。
> 勘察日期：2026-09-04（数字与结论基于当日 OCCT 树实测）。

## 0. 给后续 agent 的开场流程（每个 session 必读）

1. 读本文档（重点 §7 进度表）+ `rcad/docs/module-map.md` 的 hlr 行。
2. 用 `TodoWrite` 建当前 Stage 的任务清单；从第一个未勾选任务继续。
3. 环境事实：
   - `OCCT_SRC=C:\Users\lilu\works\OCCT`；TKHLR 源码 `$OCCT_SRC/src/ModelingAlgorithms/TKHLR/`。
   - 编译验证：`cd rcad && cargo check -p rcad-algo`（每改一处跑一次）。
   - rcad 是 rcad-pro 根仓库下的 git 子模块；提交只加本任务文件（子模块工作树里可能有其它 session 的遗留修改，不要混入）。
4. 方法论与纪律全文见根 AGENTS.md，要点：
   - 阶段 1 只做 1:1 语句翻译（不思考算法、不跑测试）；阶段 2 才调试。
   - 每函数顶部标注 `// OCCT <文件>.cxx L<起>-<止>`；`.gxx`/`.lxx` 内联翻译；不用 emoji。
   - 单文件 ≤2000 行（`HLRBRep_Data.cxx` 2683 行须按 OCCT 文件内的 section 注释拆 `brep/data/` 子模块）。
   - 代码注释全英文；「OCCT 没有的概念直接删」。
   - OCCT 源码普遍 `#define No_Exception` / `No_Standard_OutOfRange`（数组越界保护被关闭）——移植不得加 rcad 自创检查。
   - 每个 Stage 结束跑回归：`cargo test -p rcad-algo --lib helix` + `cargo test -p rcad-algo --test tkhelix_gtests`（不破坏既有基线）。

## 1. 规模与结构（勘察结论）

总量 ~44k 行 C++（实现层）：

| 包 | 文件数 | 行数 | 内容 |
|---|---|---|---|
| HLRBRep | 110 | 26,051 | 入口/内部DS/adaptor/相交/The*实例化 |
| Contap | 40 | ~7,300 | 轮廓线引擎（Contour 2389 行主引擎） |
| HLRAlgo | 39 | ~5,700 | Projector/EdgeStatus/编码盒/PolyData 等 |
| HLRTopoBRep | 16 | ~3,100 | DSFiller/Data/OutLiner/FaceIsoLiner |
| Intrv | 7 | ~1,060 | 区间代数（Interval+Intervals+Position） |
| TopCnx | 2 | 205 | EdgeFaceTransition |
| TopBas | 2 | ~200 | TestInterference 值类型 |
| HLRAppli | 2 | 172 | ReflectLines |

HLRBRep 内含 ~24 个 `The*` 类，全部是 TKGeomAlgo 泛型引擎的模板实例化：
IntWalk_IWalking、IntStart_SearchOnBoundaries/SearchInside、IntImp_IntCS、
IntImp_ZerCSParFunc、IntImpParGen_Intersector、IntCurve_IntCurveCurveGen、
IntCurveSurface 系 —— **翻译它们必须连带翻译泛型体**。
rcad 已有：`geomalgo/int_patch/imp_prm/i_walking.rs` = IntWalk_IWalking 1:1 ✅。

关键入口 API（精确线路）：
`HLRBRep_Algo`(Add→Load) → `HLRBRep_InternalAlgo`(Projector/Update/Hide/Select)
→ `HLRBRep_HLRToShape`（V/Rg1/RgN/OutLine/Iso × Visible/Hidden compounds，
`CompoundOfEdges(type, visible, In3d)`）。
poly 线路：`HLRBRep_PolyAlgo::Update` → `HLRBRep_PolyHLRToShape`（本计划 Stage 5，单列）。

## 2. 验收标准（诚实边界）

- **DRAW 网格**：exact_hlr 85 + poly_hlr 89 中仅 6 个自包含（各 grid 的
  `Plate` / `bug25813_1` / `bug25813_3`），其余 ~168 个依赖 `locate_data_file`
  外部数据且文件不在 `$OCCT_SRC/data`（已逐一验证，data/occ 仅 44 个文件）
  → **永久排除**（AGENTS.md 同类规则）。
- **断言机制**（已核实）：
  - `COMPUTE_HLR` 是 `tests/hlr/begin` 里的 Tcl proc：
    `vcomputehlr a result -algoType $algotype`，algo 后加 `build3d result`。
  - `vcomputehlr` 的 C++ 实现是 `ViewerTest_ObjectCommands.cxx:3141-3377` 的
    `VComputeHLR`；支持无头友好的显式 `anEye aDir anUp` 9 元组形式；
    result = 可见 compound（在前）+ 隐藏 compound。
  - 断言在 `checkprops result -l ${length}`（`CheckCommands.tcl` proc 531 行起）：
    `lprops` 取 Mass，正则提取，相对误差 `depsilon` 默认 1e-2（case 里的
    `set depsilon` 不影响该比较）。case 文件 `set length 6.34983` 是期望总长。
- **rcad 侧验收 = 三条**：
  (a) 阶段单测（数学锚点：Projector 变换、Intrv 区间 13 种 Position、
  EdgeStatus 布尔语义、编码盒编解码等）；
  (b) `tools/gen-occt-ref/gen_hlr_ref.py` 用无头 DRAWEXE（显式 eye/dir/up，
  viewname→方向按 `ViewerTest_ViewerCommands.cxx` 的 VViewProj 语义硬编码）
  生成参考 JSON（nbshapes + visible/hidden 两个 sub-compound 的 lprops Mass），
  rcad 端同投影方向跑 1:1 管线，断言 `|ref-got|/ref ≤ 1e-2`；
  (c) occt-test-gen 接入 hlr 网格翻译，生成
  `tests/occt/tests/generated_occt_boolean_hlr_exact_hlr.rs`（仅 3 个自包含用例）。
- **poly 线路单列**（HLRBRep_PolyAlgo 4451+636 行依赖 Poly_Triangulation/
  BRepMesh 1:1，与 AGENTS.md「mesh 不 1:1 移植」冲突）→ 需用户决策后另立
  runway，本计划只覆盖精确线路（Stage 0-4）。
- 无 GTests（TKHLR/GTests 目录只有 FILES.cmake）。

## 3. 目录布局（module-map L114-119 约定，外层 hlr/ 已存在）

```
rcad/libs/rcad-algo/src/hlr/
├── mod.rs          # 入口：pub mod 声明 + pub use（lib.rs 加 pub mod hlr;）
├── tests.rs        # Stage 单测
├── commands.rs     # HLRTest DRAW 命令层等价（Stage 4）
├── algo/           # HLRAlgo（Projector、EdgeStatus、编码盒、PolyData…）
├── brep/           # HLRBRep（Algo/InternalAlgo/Data/Hider/HLRToShape/adaptor…）
├── contap/         # Contap（Contour/ContAna/SurfFunction/The*…）
├── topo_brep/      # HLRTopoBRep（DSFiller/Data/OutLiner/FaceIsoLiner…）
├── intrv/          # Intrv（Interval/Intervals/Position）
├── top_bas/        # TopBas（TestInterference）
├── top_cnx/        # TopCnx（EdgeFaceTransition）
└── appli/          # HLRAppli（ReflectLines）
```

OCCT 类 → snake_case 文件，一 OCCT 类一文件（Data 例外拆子模块）。
文件头标注 OCCT 行号；BRep 层映射注释沿 `helix/helix_brep/mod.rs` 先例
（MakeVertex→add_tvertex、MakeEdge→edge TShape、UpdateEdge→TShape 替换）。

## 4. 阶段划分

### Stage 0：叶子包（~2.7k 行，1 session）—— 无外部依赖
- `top_bas/`：TopBas_TestInterference（158 hxx 全内联 + 36 _0.cxx）。
- `intrv/`：Interval（188 hxx / 151 cxx / 325 lxx）+ Intervals（70/248/38）
  + Position 枚举（37）。RealFirst/RealLast = ±DBL_MAX，Epsilon 按 OCCT
  Standard_Real 宏语义。
- `top_cnx/`：EdgeFaceTransition（71/134）——复用 rcad
  `geomalgo/top_trans::CurveTransition`，State/Orientation 用
  `rcad_kernel::topods`。
- `appli/`：HLRAppli_ReflectLines（74/98）。
- `algo/` 起步：`HLRAlgo` 静态 MinMax/编码盒 + π/14 三角表（87/205）、
  `HLRAlgo_Projector`（128/396/77：gp_Ax2/gp_Trsf 构造、Project×3、Shoot、
  Focus；复用 kernel gp::Trsf）。
- DoD：cargo check 全绿；`hlr/tests.rs`（Intrv Position 13 分支、Interval
  Fuse/Cut、Projector round-trip 数值锚点）；不破坏 tkhelix/既有单测；提交。

### Stage 1：HLRAlgo 数据结构（~1.8k 行，1 session）
EdgeStatus(119/123)、Interference(160+39)、Intersection(86/47/113)、
BiPoint(262/240)、Coincidence(82)、EdgesBlock(164/32，MinMaxIndices 4×4
编码盒——**位级保真**，culling 依赖它)、WiresBlock(67/24)、PolyMask(37)、
TriangleData(33)、PolyInternalNode(71/19)、PolyInternalSegment(33)、
PolyHidingData(73)、EdgeIterator(69/94/70)。
DoD：单测覆盖 EdgeStatus 状态机与编码盒编解码；提交。

### Stage 2：Contap 路线（2~3 session）⚠️ 最大不确定区，开工先清点
- 2a 前置审计（0.5 session）：清点并补齐泛型引擎 ——
  `IntStart_SearchOnBoundaries.gxx` / `IntStart_SearchInside.gxx`（体量未测，
  最大未知数）、`Intf_Interference/Intf_Polygon2d`（rcad intf.rs 只有 trait
  层）、`IntCurve_IntConicConic` 缺 8 个重载（int_conic_conic.rs 仅
  Ellipse-Ellipse）、`IntImp_IntCS/ZerCSParFunc/IntImpParGen_Intersector`、
  `IntCurveSurface_HInter` 组装、`Bnd_BoundSortBox`+`Bnd_Range`、
  `Extrema_EPCOfExtPC2d`、`GProp_PEquation`、`BRepAdaptor_Curve2d`。
  产出：补齐后的缺口勾选表写入本文档 §7。
- 2b Contap 包 1:1（~7.3k 行）：SurfProps → HContTool → HCurve2dTool →
  Point/Line → SurfFunction → ArcFunction → ContAna → Contour(2389) →
  The* 实例化（TheIWalking 复用 i_walking.rs；TheSearch/TheSearchInside 按
  2a 结果组装）。
- DoD：Contap 对解析面（圆柱/圆锥/球轮廓线）的单元锚点测试（参考值
  DRAWEXE 实测生成）；逐文件提交。

### Stage 3：HLRBRep 精确线路（4~6 session）
自底向上：
- 3a adaptor/tool 层：Curve(215/632/169)/Surface(197/335/271)/CurveTool/
  SurfaceTool/BCurveTool/BSurfaceTool/LineTool/CLProps/SLProps。
- 3b 相交层：Intersector(99/858)/InterCSurf(211/578)/CInter(275/85) +
  全部 The*OfInterCSurf/CInter 实例化。
- 3c 前置：BRepApprox_Approx 链路（AppDef_BSplineCompute —— 正是当前
  math runway「BSplComputeLine→BSplineCompute」的下一站，按 AppParCurves
  已 landed 的 LeastSquare/ResolConstraint/Gradient 底座继续）。
- 3d HLRTopoBRep 包：DSFiller(758)/Data(303+106)/FaceData/VData/
  FaceIsoLiner(488)/OutLiner(344)。
- 3e 干涉数据结构：EdgeIList/VertexList/EdgeBuilder(521)/
  EdgeInterferenceTool/AreaLimit/FaceIterator/BiPoint/BiPnt2D。
- 3f `Data`（2683 行，拆 `brep/data/` 子模块）。
- 3g Hider(852) → ShapeBounds → InternalAlgo(1020) → Algo → HLRToShape。
- DoD：每文件 cargo check + 提交；3g 末跑通 `HLRBRep_Algo` 对盒体/圆柱的
  烟囱测试（自建，不急对数值）。

### Stage 4：验收闭环（1~2 session）
- 4a `tools/gen-occt-ref/gen_hlr_ref.py`（模板=gen_helix_ref_topology.py；
  捕获 nbshapes + 两 sub-compound lprops Mass；先做无头 vcomputehlr spike）。
- 4b `hlr/commands.rs`（HLRTest 等价层：compute_hlr / 投影参数显式化；
  静态量语义对齐——HLRBRep_InternalAlgo TRACE 默认 true、Contap 两处
  Tolpetit 冲突值照抄）。
- 4c occt-test-gen 接入（4 个接入点：main.rs L6 mod 声明、L463 batch
  分支、L3390 translate_draw_script_inner 分发、DrawPort needs_hlr import）。
- 4d 生成并跑 exact_hlr 3 用例，长度断言 1e-2 相对误差 + nbshapes 佐证；
  增强：自建简单图元×视图方向的参考用例扩充覆盖。
- DoD：3 用例全绿，module-map.md hlr 行更新为对齐状态+基线数字。

### Stage 5（用户决策后另立）：poly 线路
HLRBRep_PolyAlgo(4451+636)、HLRAlgo_PolyAlgo/PolyData(1128)/
PolyInternalData(1039)/PolyShellData、PolyHLRToShape、
Poly_Triangulation/Polygon3D/PolygonOnTriangulation/BRepMesh 依赖 ——
与 mesh 策略冲突，默认推迟。

## 5. 已知坑（勘察实证）

- `HLRTest.cxx` 的 `static handle<HLRBRep_Algo> hider` 与
  `Contap_HContTool.cxx` 的 `static double uinf..vsup` 是文件级可变静态
  → 翻译成 thread_local（沿 helix commands.rs 的 theHelixAxis 先例）。
- 同名冲突常量照抄：Contap_ContAna Tolpetit=1.e-8 vs Contap_Contour
  Tolpetit=1.e-10；HLRBRep_Data CutLar=2.e-1 / CutBig=1.e-1 / SIZEUV=8；
  HLRAlgo.cxx π/14 三角表。
- OCCT 数组越界保护被宏关闭 → 不加 rcad 自创检查；1-based 索引语义用
  val()/下标映射沿先例。
- EdgesBlock::MinMaxIndices 的位打包必须位级保真。
- `HLRBRep_PolyAlgo` 私有矩阵副本 TMat/TTMa/TIMa（poly 线路，Stage 5）。
- IntSurf_Allocator（NCollection_BaseAllocator 共享缓冲）→ Rust 侧用
  Vec/Arc 语义等价即可，注释说明。
- 测试超时 RCAD_TEST_TIMEOUT_SECS 默认 10s；HLR 在圆柱/圆环上应毫秒级，
  超时即查对齐。

## 6. 依赖缺口总表（rcad 侧，勘察 2026-09-04）

| 组件 | 状态 | 位置/备注 |
|---|---|---|
| AppDef_MultiLine | ✅ 完整 | geomalgo/app_def.rs |
| AppParCurves/MultiBSpCurve/LeastSquare/ResolConstraint/Gradient | ✅ 完整 | geomalgo/app_par_curves.rs（近期 landed） |
| AppDef_BSplineCompute | ❌ 缺 | math runway 下一站（Stage 3c） |
| BRepApprox_Approx/ApproxLine | ❌ 缺 | Stage 3c |
| Adaptor3d_TopolTool/BRepTopAdaptor_TopolTool/HVertex | ❌ 缺 | topalgo/brep_top_adaptor 仅 class2d/fclass2d（Stage 3a 前补） |
| BRepAdaptor_Curve/Surface | ◐ | rcad-brep/adaptor.rs（EdgeAdaptor/FaceAdaptor） |
| BRepAdaptor_Curve2d | ❌ 缺 | Stage 3a 前补 |
| Geom2dHatch_Hatcher/HatchGen_Domain | ❌ 缺引擎 | hatch.rs 仅类型（视 Contap/HLR 用途定） |
| IntCurve_IntConicConic | ◐ 仅 Ellipse-Ellipse | 缺 8 个重载（Stage 2a） |
| IntCurveSurface_HInter | ◐ 组装缺 | 采样类已在 int_curv_surf.rs（Stage 2a/3b） |
| IntRes2d | ✅ 数据类型完整 | geomalgo/int_res2d |
| IntSurf Quadric/LineOn2S/PntOn2S/PathPoint/InteriorPoint | ✅ | geomalgo/int_surf + int_patch/imp_prm |
| IntWalk_IWalking | ✅ 1:1 | int_patch/imp_prm/i_walking.rs |
| IntStart_SearchOnBoundaries/SearchInside | ❓ 未清点 | Stage 2a 首任务 |
| Intf_Interference/Polygon2d | ❌ 缺 | intf.rs 仅 trait 层（Stage 2a） |
| IntAna/IntImpParGen | ✅ | kernel base/int_ana.rs + int_patch/imp_prm |
| Extrema_ExtPC/LocateExtPC/ExtPC2d/POnCurve2d | ✅ | kernel base/extrema.rs |
| Extrema_EPCOfExtPC2d | ❌ 缺 | Stage 2a |
| Bnd_Box/Box2d | ✅ | kernel math/bnd |
| Bnd_BoundSortBox/Bnd_Range | ❌ 缺 | Stage 2a |
| BndLib_AddSurface | ✅ | kernel base/bnd_lib（surface_bounding_box） |
| TopTrans_CurveTransition | ✅ | geomalgo/top_trans |
| CSLib/GeomLProp | ✅ | kernel math/cs_lib.rs、base/geom_lprop |
| GProp_PEquation | ❌ 缺 | Stage 2a |
| ProjLib/GeomProjLib | ◐ | kernel base/proj_lib、geom_proj_lib |
| Poly_*（OCCT 对齐） | ❌ 缺 | Stage 5（mesh 策略冲突） |
| TopExp/BRep_Tool/TopoDS/MakeEdge 等价 | ✅ | kernel topo/（helix_brep 先例映射） |

## 7. 进度勾选表（每 session 更新）

- [x] 计划落盘 + module-map 指向（2026-09-04）
- [x] **Stage 0** 叶子包：top_bas / intrv / top_cnx / HLRAlgo 静态+EdgesBlock+Projector（2026-09-04 完成并提交；kernel Trsf 同步补齐 OCCT scale/form 语义 + SetTransformation/SetScaleFactor/SetTranslationPart/Invert/VectorialPart/Value + gp_Vec/gp_Dir/gp_Lin::Transform；单测 12 个全绿。**HLRAppli_ReflectLines 推迟至 Stage 3g 之后**——它依赖 HLRBRep_Algo/HLRToShape/BRepLib::SameParameter）
- [x] **Stage 1** HLRAlgo 数据结构（EdgeStatus/Interference/Intersection/BiPoint/Coincidence/WiresBlock/Poly 系列小结构/EdgeIterator；2026-09-04 完成并提交。要点：EdgeStatus::Hide 的 `if (!OnFace)` guard 照抄、EdgeIterator 的隐藏区间缓存语义（当前段=上一可见段 End 至下一可见段 Start）、指针字段→借用映射；单测 15/15）
- [x] **Stage 2a（审计完成 2026-09-04）**：泛型引擎闭包已逐个实测（The* 类 `_0.cxx` 的 gxx include → 实际行数 → rcad 现状）：

| 泛型体（OCCT 实测行数） | 消费者 | rcad 状态 |
|---|---|---|
| IntWalk_IWalking.gxx（3152） | Contap_TheIWalking | ✅ int_patch/imp_prm/i_walking.rs（3037 行 1:1） |
| IntStart_SearchOnBoundaries.gxx（1232） | Contap_TheSearch | ◐ int_patch/so_on_bounds.rs 1:1，但 domain=UV 矩形近似；Contap 需真实 TopolTool（面 wires 边界弧）→ 2a-4 |
| IntStart_SearchInside.gxx（313） | Contap_TheSearchInside | ✅ imp_prm/search_inside.rs（TopolTool→网格采样适配） |
| IntImp_IntCS.gxx（198） | TheExactInterCSurf | ◐ int_patch/int_cs.rs 仅 canonic 路径；非 canonic（多面体）路径 → 2a-3 |
| IntImp_ZerCSParFunc.gxx（116） | TheCSFunctionOfInterCSurf | ❌ 新写 |
| IntImpParGen_Intersector.gxx（824） | TheIntersectorOfTheIntConicCurve | ❌ 新写（imp_prm/function_set_root.rs 是同族但不等同） |
| IntCurve_IntCurveCurveGen.gxx（1033） | HLRBRep_CInter | ✅ geomalgo/int_curve_curve_gen.rs（Geom2dInt_GInter 实例） |
| IntCurve_IntPolyPolyGen.gxx（1797） | TheIntPCurvePCurveOfCInter | ❌ 新写（int_curve_curve_gen.rs L14 明确 NOT yet） |
| IntCurve_UserIntConicCurveGen.gxx（889） | IntConicCurveOfCInter | ❌ 新写 |
| IntCurve_IntConicCurveGen.gxx（92） | TheIntConicCurveOfCInter | ❌ 新写 |
| IntCurve_Polygon2dGen.gxx（383） | ThePolygon2dOf… | ❌ 新写 |
| IntCurve_DistBetweenPCurvesGen.gxx（105） | TheDistBetweenPCurves | ❌ 新写 |
| IntCurve_ExactIntersectionPoint.gxx（271） | ExactIntersectionPoint… | ❌ 新写 |
| Intf_InterferencePolygonPolyhedron.gxx（1368） | TheInterferenceOfInterCSurf | ❌ 新写（intf.rs 仅 PIType/SectionPoint 数据类） |
| IntCurveSurface_HInter.cxx（581） | HInter 组装 | ❌ 新写（int_curv_surf.rs 仅采样类） |
| Intf_Interference 基类 + Intf_Polygon2d（488） | 上述 Intf 引擎的地基 | ❌ 新写 |
| Bnd_Range（327）→ 消费者 Contap_TheIWalking.hxx；Bnd_BoundSortBox（774）→ HLRBRep_InterCSurf | HLRBRep_Data | ❌ 本轮起补 |
| GProp_PEquation（342）→ 消费者 HLRBRep_Surface.cxx | HLRBRep_Data | ❌ 新写 |
| Extrema_EPCOfExtPC2d → 消费者 Contap_HContTool.cxx | Contap | ❌ 2a-2 |
| IntCurve_IntConicConic（缺 8/9 重载）→ 消费者 CInter/IntConicCurveOfCInter | Stage 3b | ◐ |
| Geom2dHatch_Hatcher/Intersector → 消费者 HLRTopoBRep_FaceIsoLiner.cxx（在 Stage 3d，非 Contap） | Stage 3d | ❌ |

  **Stage 2 内部补齐顺序（重排）：**
  - 2a-1 小项批次（~2300 行，低风险）：Bnd_Range、Bnd_BoundSortBox、GProp_PEquation、Intf_Interference/Polygon2d、IntCurve_IntConicCurveGen、IntCurve_Polygon2dGen、IntCurve_DistBetweenPCurvesGen、IntCurve_ExactIntersectionPoint
  - 2a-2 中项批次：IntImp_ZerCSParFunc、IntImp_IntCS 非 canonic 补全、IntImpParGen_Intersector
  - 2a-3 大项批次：Intf_InterferencePolygonPolyhedron、IntCurve_UserIntConicCurveGen、IntCurve_IntPolyPolyGen、IntCurveSurface_HInter 组装
  - 2a-4 TopolTool 真实版：Adaptor3d_TopolTool/BRepTopAdaptor_TopolTool/HVertex/BRepAdaptor_Curve2d + so_on_bounds 的 BRep 边分支保真（Contap_TheSearch 依赖）
- [x] **2a-1** 小项批次翻译（✅ Bnd_Range → `rcad-kernel/src/math/bnd/range.rs`；✅ Bnd_BoundSortBox → `bnd/bound_sort_box.rs`（voxel grid + LargeBoxes 通道 + 位掩码三轴过滤，Compare(gp_Pln) 因 Bnd_Box::IsOut(Pln) 未移植而 defer 并注释）；✅ GProp_PGProps + GProp_PEquation → `kernel base/gprop/pg_props.rs` + `pequation.rs`（**关键语义：MatrixOfInertia 经 GProp::HOperator 从原点惯性移轴到质心**，GProp_GProps.cxx L110-115；Jacobi 复用 math_Jacobi 移植）；✅ **2a-1c** IntCurve 三泛型 → `geomalgo/int_curve_generics.rs`：Polygon2dGen（ComputeWithBox 的 deflection 收敛循环逐行，`t` 与 dx/dy 在分支前声明照抄）、DistBetweenPCurvesGen（实现 FunctionSetWithDerivatives）、ExactIntersectionPoint（math_FunctionSetRoot 驱动 + 四段 bound-widening 重试循环），`TheCurveTool` 模板参数 → `ProjPCurveTool` trait；kernel BndBox2d 补 IsOut(Bnd_Box2d)（Bnd_Box2d.cxx L456-511 快慢双路径）。测试 2/2。⚠️ IntConicCurveGen(92)/UserIntConicCurveGen(889) 的 Perform 引擎体绑定 HLRBRep_Curve/CurveTool + IntImpParGen_Intersector → 并入 2a-2/3a）
- [x] **2a-1b** Intf 数据类：SectionPoint 补 ParamOnFirst/ParamOnSecond/IsEqual（SectionPoint.lxx L17-41）、Intf_TangentZone、Intf_SectionLine（IsEnd 末点返回 Length 的怪癖照抄）、Intf_Interference（Insert zone 合并 + Insert 线段拼接）、Intf_Polygon2d→trait（geomalgo/intf_tangent_zone.rs、intf_section_line.rs、intf_interference.rs）
- [ ] **2a-2** 中项批次翻译
- [ ] **2a-3** 大项批次翻译
- [ ] **2a-4** TopolTool 真实版
- [ ] **Stage 2b** Contap 包 1:1
- [ ] **Stage 3a** HLRBRep adaptor/tool 层（含 TopolTool/BRepAdaptor_Curve2d）
- [ ] **Stage 3b** HLRBRep 相交层 + The* 实例化
- [ ] **Stage 3c** BRepApprox_Approx 链路（AppDef_BSplineCompute）
- [ ] **Stage 3d** HLRTopoBRep 包
- [ ] **Stage 3e** HLRBRep 干涉数据结构
- [ ] **Stage 3f** HLRBRep Data（拆子模块）
- [ ] **Stage 3g** Hider→InternalAlgo→Algo→HLRToShape + 烟囱测试
- [ ] **Stage 4a** gen_hlr_ref.py（含无头 spike）
- [ ] **Stage 4b** hlr/commands.rs
- [ ] **Stage 4c** occt-test-gen hlr 接入
- [ ] **Stage 4d** exact_hlr 3 用例闭环 + module-map 更新
- [ ] Stage 5 poly 线路（用户决策后另立计划）
