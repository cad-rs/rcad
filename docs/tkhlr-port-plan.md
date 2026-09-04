# TKHLR 1:1 翻译推进计划（多 session runway）

> **交接快照（2026-09-05 session 结束时）**
> - 提交链（rcad 子模块分支 `sd-hash-wip`）：`7011f938` 计划落盘 → `89fdd3ec` Stage 0 → `78789bbf` Stage 1 → `6f0f12b8` 2a 审计+Bnd_Range → `34dfaa5b` 2a-1 批次（Intf 数据类/BoundSortBox/PGProps+PEquation）→ `da383601` 2a-1 完成（IntCurve 三泛型）→ `5ef9ae00` 2a-2 IntImp 部分（ZerCSParFunc + IntCS）→ `57ab4ff5` 2a-2 完成（IntImpParGen_Intersector 泛型引擎 + IntConicCurve 双泛型）→ `5cde756c` 2a-3①（Intf_InterferencePolygonPolyhedron 1368 行 1:1）→ `f528da46` 2a-3②（Intf_InterferencePolygon2d 819 行 1:1 + Polygon2dGen 实现 IntfPolygon2d trait）→ `3466868b` 2a-3③（IntCurve_IntPolyPolyGen 1797 行 1:1）→ **2a-3④ 完成（IntCurveSurface_HInter 组装，见 §8 新翻译模式）** → **2a-4 Adaptor3d 层完成**（Adaptor2d_Line2d + Adaptor3d_HVertex + Adaptor3d_TopolTool 全 1:1，HInter 的 UTrim/VTrim + BSplineSurface-TopolTool 分支接通）→ **2a-4 BRep 层完成，2a-4 全部关闭**（`62a7c5d3` 代码：BRepAdaptor_Curve2d + BRepTopAdaptor_HVertex + BRepTopAdaptor_FClass2d 859 行 + BRepTopAdaptor_TopolTool 653 行 全 1:1；`bb8a7913` 文档，见 §8）。
> - 回归基线（全部全绿）：algo lib **150 通过 + 1 ignore**（146+4 新增 BRep 层锚点；ignore = FClass2d 边界 On 锚点，根因见 §8）、kernel lib **645**/645、tkhelix_gtests 16/16、pavefiller_stage_tests 26/26、tkgeom_algo_gtests 134+1。
> - **下一步 = 2b：Contap 包 1:1**。源码 `$OCCT_SRC/src/ModelingAlgorithms/TKHLR/Contap/`（40 文件，~7.3k 行），建议顺序（按依赖）：
>   1. Contap_SurfProps（351，静态面属性：UVBounds/IsCanonic 域）；
>   2. Contap_Point（cxx 小）+ Contap_Line（120）+ Contap_IType（枚举）——数据类；
>   3. Contap_HCurve2dTool（52）+ Contap_HContTool（370，依赖 **Extrema_EPCOfExtPC2d——kernel 未移植，需先补**（现仅有 ExtPC2d 无参变体；EPCOfExtPC2d 带 Nbu/epsX/Tol 参数化构造，OCCT TKMath/Extrema/Extrema_EPCOfExtPC2d.cxx ~500 行））；
>   4. Contap_SurfFunction（287）+ Contap_ArcFunction（209）+ Contap_TFunction.hxx——math_FunctionSetWithDerivative 系；
>   5. Contap_ContAna（483，解析面轮廓：Line/Circle/椭圆抛物双曲 vs Cylinder/Cone/Sphere）；
>   6. Contap_TheIWalking_0/TheSearch_0/TheSearchInside_0 等 _0.cxx 实例化——**引擎已全部就绪**（i_walking.rs / so_on_bounds.rs / search_inside.rs），只需按 §8 marker 模式组装；
>   7. Contap_Contour（2389，主引擎）收尾。
>   注意：Contap_HContTool.cxx 的 `static double uinf..vsup` 是文件级可变静态 → thread_local（§5 已登记）；Contap 的 domain 参数用本 session 落地的 BRepTopAdaptor_TopolTool（`topalgo/brep_top_adaptor/topol_tool_brep.rs`）。
>   遗留对齐项（非阻塞 2b，已定位根因）：FClass2d 边界 On 兜底——FaceShapeSource 传 kernel BRep locations 表（无 identity 槽位 0）给 `edge_pcurve_on_face`（按 bop DS 表约定读），pcurve 键落空 → FaceExplorer::segment=None → nowires → In；对齐该表约定后 un-ignore `fclass2d_topol_perform_boundary_on`。
>   HInter 已就绪的挂载点：`geomalgo/int_curve_surface/`（mod=数据类+HCurveTool/HSurfaceTool/HInterHost trait、inter_utils=InterUtils.pxx、inter_impl=Inter.pxx、hinter=HInter 壳+4 锚点、quad_curv_exact=QuadCurvExact 引擎）；IntConicQuad 类壳（含 Line/Ellipse/Line×Plane 分支）在 `int_patch/int_conic_quad.rs`；IntLinTorus 在 `int_patch/int_lin_torus.rs`；ElSLib 反向 Parameters×5 在 kernel `math/el.rs`；采样类工具泛型构造器 `ThePolygonOfHInter::new_tool*/ThePolyhedronOfHInter::new_tool*` 在 `int_curv_surf.rs`。Stage 3b 的 HLRBRep_CInter/The*OfInterCSurf 实例化只需为 HLRBRep_Curve/HLRBRep_Surface 实现同一套 trait。
> - 已确立的翻译模式（沿例勿改）：OCCT gxx 模板参数 → Rust trait，**曲线类型本身也是 trait/struct 泛型参数 `C: ?Sized`**（`ParTool<C>`/`ProjectOnPCurveTool<C>`/`PCurveTool<C>`，见 `geomalgo/int_imp_par_gen.rs`）；具体实例化 = marker 类型上的 `impl<'a> ParTool<dyn Curve2dAdaptor + 'a>`（生命周期弹性，等价 C++ `const Adaptor2d_Curve2d&`）+ 具体薄壳结构体（TheIntersectorOfTheIntConicCurveOfGInter / TheIntConicCurveOfGInter）；gxx/lxx 内联翻译；每函数 `// OCCT <文件> L<起>-<止>` 标注；每个叶子翻译配 OCCT 解析锚点单测。
> - 两个已踩过的契约坑：**(1)** `ThePSurfaceTool::D1`/`TheCurveTool::D1` 的 P 输出就是曲面/曲线点本身（测试桩写成 ZERO 会使 F≡0、触发 f2<=eps 的 done=true 提前返回——症状是 root 原样返回起点）；**(2)** `GProp_GProps::MatrixOfInertia` 经 HOperator 移轴到质心。另有怪癖照抄清单见 §5。**新坑 (3)**：IntConicCurveGen 的 C/E ctor+Perform 在 `!D1.IsClosed()` 时只注入 SetEquivalentParameters(F, F+2π)，**域边界 first/last 不变**——弧域外的交点角参数仍会被丢弃（锚点测试 circle_ctor_quirk_non_closed_domain 记录了该行为）。
> - 子模块工作树有 3 个**非本任务**的遗留修改（`algo_ext/mod.rs`、`fillet/topopebrepbuild.rs`、`tests/tkgeom_algo_gtests.rs`）——提交时只 `git add` 本任务文件，勿混入。
> - shell 的 cwd 会被重置到 `C:\Users\lilu\works\rcad-pro`（根仓库），编译/测试前先 `cd rcad`。

>
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
| Adaptor3d_TopolTool/BRepTopAdaptor_TopolTool/HVertex | ◐ Adaptor3d 层 ✅（topalgo/adaptor3d/：TopolTool+HVertex 1:1）；BRep 层待补（2a-4 续） | topalgo/brep_top_adaptor 仅 class2d/fclass2d（Stage 3a 前补） |
| BRepAdaptor_Curve/Surface | ◐ | rcad-brep/adaptor.rs（EdgeAdaptor/FaceAdaptor） |
| BRepAdaptor_Curve2d | ❌ 缺 | Stage 3a 前补（2a-4 续；kernel BRepTool::curve_on_surface 已备） |
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
| IntImp_ZerCSParFunc.gxx（116） | TheCSFunctionOfInterCSurf | ✅ geomalgo/int_imp/zer_cs_par_func.rs（5ef9ae00） |
| IntImpParGen_Intersector.gxx（824） | TheIntersectorOfTheIntConicCurve | ✅ geomalgo/int_imp_par_gen.rs 泛型 `Intersector<C, PT, JT>`（引擎 1:1）+ geom2d_int.rs 的 GInter 具体薄壳（引擎体最初由 InterCurveCurve session 以 GInter 具体实例落地，本轮泛型化去重） |
| IntCurve_IntCurveCurveGen.gxx（1033） | HLRBRep_CInter | ✅ geomalgo/int_curve_curve_gen.rs（Geom2dInt_GInter 实例） |
| IntCurve_IntPolyPolyGen.gxx（1797） | TheIntPCurvePCurveOfCInter | ❌ 新写（int_curve_curve_gen.rs L14 明确 NOT yet） |
| IntCurve_UserIntConicCurveGen.gxx（889） | IntConicCurveOfCInter | ✅ geomalgo/user_int_conic_curve_gen.rs（5 Perform + 5 InternalPerform 泛型 1:1；conic-typed pcurve 臂调用 IntConicConic 各重载——E-E 已实现，其余 13 个为 documented unimplemented，Stage 3b） |
| IntCurve_IntConicCurveGen.gxx（92） | TheIntConicCurveOfCInter | ✅ geomalgo/int_conic_curve_gen.rs 泛型壳（gxx C/E/P/H ctor + lxx Lin ctor/Perform×6）|
| IntCurve_Polygon2dGen.gxx（383） | ThePolygon2dOf… | ❌ 新写 |
| IntCurve_DistBetweenPCurvesGen.gxx（105） | TheDistBetweenPCurves | ❌ 新写 |
| IntCurve_ExactIntersectionPoint.gxx（271） | ExactIntersectionPoint… | ❌ 新写 |
| Intf_InterferencePolygonPolyhedron.gxx（1368） | TheInterferenceOfInterCSurf | ✅ geomalgo/intf_interference_polygon_polyhedron.rs（7 ctor + 6 Perform + Interference×2 + Intersect×2 + Intf::PlaneEquation；ToolPolygon3d/ToolPolyh trait + TheInterferenceOfHInter alias；2a-3① 5cde756c） |
| Intf_InterferencePolygon2d.cxx（819） | IntPolyPolyGen 的前置 | ✅ geomalgo/intf_interference_polygon2d.rs（Polygon2dGen 实现 IntfPolygon2d trait；2a-3② f528da46） |
| IntCurveSurface_HInter.cxx（581） | HInter 组装 | ✅ geomalgo/int_curve_surface/（2a-3④：数据类+HCurveTool/HSurfaceTool/HInterHost trait、inter_utils、inter_impl、hinter、quad_curv_exact） |
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
- [x] **2a-2** 中项批次（✅ 2026-09-04 完成并提交。①IntImp 部分（5ef9ae00）：`geomalgo/int_imp/`（ZerCSParFunc 全文 + IntCS 全文含 3 次 w-restart 循环与 MarginCoef 边界扩张）；`ThePSurfaceTool`/`TheCurveTool` → `PSurfaceTool`/`CurveTool3d` trait；kernel `math_FunctionSetRoot::SetTolerance` 补齐；Precision 常量落位。解析锚点：平面 z=0 × 直线 (w,0,w) 收敛到原点 2/2。②IntImpParGen_Intersector(824)：引擎体已在 InterCurveCurve session 以 GInter 具体实例落地（geom2d_int.rs），本轮抽出为泛型 `geomalgo/int_imp_par_gen.rs`——`Intersector<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>>`（824 行 1:1：FindU/FindV/And_Domaine_Objet1_Intersections/Perform 含封闭隐式曲线周期偏移与 Calcule_Toutes_Transitions 内联）+ IntImpParGen.cxx 静态函数（NormalizeOnDomain/DeterminePosition/DetermineTransition×2）+ 泛化 MyImpParTool；geom2d_int.rs 删除硬编码副本（-1134 行），恢复为 `impl<'a> ParTool<dyn Curve2dAdaptor + 'a>` marker + 具体薄壳包装（生命周期弹性 = C++ `const Adaptor2d_Curve2d&`）。③IntConicCurve 双泛型：`int_conic_curve_gen.rs`（gxx C/E/P/H ctor + lxx Lin/Perform×6，含 !IsClosed→SetEquivalentParameters 怪癖）+ `user_int_conic_curve_gen.rs`（889 行 1:1：5 ctor + 5 Perform 含 NbIntervals>1 复合区间循环（Ok 中断怪癖照抄）+ 5 InternalPerform 按曲线类型分发，conic 臂走 IntConicConic（E-E 实装，其余 documented unimplemented 待 3b），default 臂走 IntConicCurveGen→Intersector）。锚点测试 5 个：圆×共线 Bezier 精确交点（u=(17±2√17)/34）经 GInter 壳与 User 变体双路径、弧域怪癖两态、线×Bezier 单交点、IntImpParGen 静态函数。⚠️ 工作区 3 个非本任务文件未混入。回归：algo lib 120（115+5）、kernel 645、tkhelix 16、pavefiller 26、tkgeom_algo_gtests 134+1）
- [◐] **2a-3** 大项批次（**① Intf_InterferencePolygonPolyhedron.gxx（1368）✅ 5cde756c**：7 ctor + 6 Perform（含 Bnd_BoundSortBox 网格/非网格双形态）、Interference 无网格版（±normal 偏移双穿越）、Interference 网格版（MKK 2007 边界外扩 Beg0/End0）、Intersect 5 参（活跃 #else 分支，#if 0 草稿不移植）与 9 参变体、sVertex/sEdge/dPiE 分类 + KHROMOV 2001 边界挠度分支、Extrema_ExtElC(line,line) 边贴近尾段（IsParallel 跳过）、IsInSegment、Pourcent3 表；Intf::PlaneEquation（Intf.cxx L20-42）；`ToolPolygon3d<P>`/`ToolPolyh<H>` trait + ThePolygonToolOfHInter/ThePolyhedronToolOfHInter marker + TheInterferenceOfHInter alias（引擎无模板态状态，实例化 = 调用点的工具对）；补 ThePolygonOfHInter::bounding / ThePolyhedronOfHInter::components_bounding 访问器；Interference 基类补 mySPoins pub(crate) 访问器。锚点 5/5：竖直线×平面多面体精确穿越（1-based 三角形地址、Infinite=true 的 EDGE-kind 怪癖）、多边形段穿越（±normal 偏移对都落在平面点）、perform 复位、PlaneEquation/IsInSegment 数值锚。**② Intf_InterferencePolygon2d.cxx（819）✅ f528da46**（两 ctor/两 Perform、Clean 的 Only1Seg 怪癖与 decal 回卷、Intersect 的 parO/parT 1-based 槽 + sinTeta/rayIntf 相切带 + 50 倍贴近悬挂 + L718 封闭接缝 guard；Polygon2dGen 实现 IntfPolygon2d trait。锚点 5/5：开放 8 字自交精确 E/E 点、封闭接缝抑制 0 点、双方形穿越 (4,2)+(2,4)、disjoint 拒绝、perform 复位。**③ IntCurve_IntPolyPolyGen.gxx（1797）✅ 3466868b**（Perform×3 + findIntersect + HeadOrEndPoint + GetIntersection 全 1:1，Tol/TolConf 递归换位照抄；GInterPolyTool<'a> marker + TheIntPCurvePCurveOfGInter wrapper；IntCurveCurveGen 补 intcurvcurv 成员并接通两个 unimplemented 臂。锚点 5/5：直线×直线 (1,1) 精确 √2 参数、degree-1 开放 8 字 BSpline 自交 (u=1/6, v=5/6)、X 交二次 Bezier 原点 0.5/0.5、disjoint 空、min-sample 访问器。**④ IntCurveSurface_HInter 组装 ✅**（HInter.cxx 581 + Inter.pxx 1006 + InterUtils.pxx 1637 全 1:1 → `geomalgo/int_curve_surface/`：mod=TransitionOnCurve/IntersectionPoint/IntersectionSegment/Intersection 基类（PARAMEQUAL 1e-8 去重）+HCurveTool/HSurfaceTool/HInterHost trait+Adaptor3d{Curve,Surface}Basis；inter_utils=Collect/Sort/Process/ProcessIntAna/ComputeAppendPoint/ComputeTransitions/ComputeParamsOnQuadric/DoSurface/DoNewBounds/Decompose/ClampUV + EstLimForInfSurf/Extr/Revl/Offs + ProjectIntersectAndEstLim + ProcessLinTorus + PerformCurveQuadric；inter_impl=Perform/PerformBounds/PerformPolygon/Polyhedron/PolygonPolyhedron/BSB/InternalPerform×3/PerformConicSurf×5/AppendIntAna；hinter=HInter 壳 + AppendPoint/AppendSegment；quad_curv_exact=TheQuadCurvExactHInter+TheQuadCurvFunc+QuadricCurveExactInterUtils（math_FunctionAllRoots 路径）。连带：IntConicQuad 类壳 + Line×Quadric/Line×Plane/Ellipse×Quadric 分支（int_conic_quad.rs）、IntAna_IntLinTorus（int_patch/int_lin_torus.rs）、ElSLib 反向 Parameters×5（kernel math/el.rs）、采样类工具泛型构造器 new_tool*（int_curv_surf.rs）、ZerCSParFunc/IntCS 的 S/C ?Sized 放宽。锚点 4/4：直线×平面（IntAna 直连，w=√2 单位方向参数化）、直线×圆柱（两穿越 In/Out，u=π→0）、Bezier×平面 QuadCurvExact（w=(3±√3)/6 精确根）、Bezier×平面 IntCS 精化（perform_polygon_polyhedron→干涉+ZerCSParFunc+IntCS）。UTrim/VTrim 与 InternalPerformPolygonBounds 的 BSplineSurface-TopolTool 分支 = documented unimplemented（2a-4）））
- [x] **2a-4** TopolTool 真实版——Adaptor3d 层 ✅：`topalgo/adaptor2d/line2d.rs`（Adaptor2d_Line2d 全 1:1 + Curve2dAdaptor impl + Trim/Load/Continuity/Intervals）、`topalgo/adaptor3d/hvertex.rs`（Adaptor3d_HVertex + ElCLib::Parameter(Lin2d,Pnt2d) 锚点）、`topalgo/adaptor3d/topol_tool.rs`（Adaptor3d_TopolTool 全 1:1：Initialize(S) 四 restriction 线 + 锥顶附加线、Init/More/Value/Next、Initialize(C) 端点顶点迭代、Classify/IsThePointOn 全状态机、ComputeSamplePoints 含 pole-grid Analyse + 各向异性检查、SamplePnts/BSplSamplePnts 自适应采样（NIZHNY-EMV/MinPnts 怪癖照抄）、GetConeApexParam；Analyse 为 pub(crate) 共享，OCCT 双份 file-static 注明）。连带：kernel `math/el.rs` 补 elclib_line_parameter_2d；HSurfaceTool 的 u_trim/v_trim 落地为「同一曲面收窄参数窗口」（GeomAdaptor_Surface::UTrim cxx L886-892 语义）+ bezier/bspline 访问器；Surface 关联类型去 ?Sized（handle 指向完整对象，按值返回）；HInter `internal_perform_polygon_bounds` 的 BSplineSurface 分支接通（Inter.pxx L489-510：UTrim→VTrim→TopolTool::SamplePnts→Parameters→Polyhedron(Upars,Vpars)）。锚点 8 个（restrictions 四线序/Classify IN-OUT-ON/曲线端点顶点/SamplePnts plane 6x6 min-clamp/锥顶参数/Line2d/HVertex）。
- [x] **2a-4 续** TopolTool 真实版 BRep 层（✅ 2026-09-05：`topalgo/brep_adaptor/curve2d.rs`（BRepAdaptor_Curve2d = pcurve + Edge/Face 访问器，Curve2dAdaptor 委托）；`topalgo/brep_top_adaptor/hvertex_brep.rs`（BRepTopAdaptor_HVertex：Value 的 RealFirst 怪癖照抄、Parameter=BRep_Tool::Parameter、Resolution 三段 refine 全 1:1）；`topalgo/brep_top_adaptor/fclass2d_topol.rs`（BRepTopAdaptor_FClass2d 859 行 1:1：逐 wire 采样 + 面积/周长比驱动 QuasiUniformDeflection 重离散循环 + TabOrien + Perform/TestOnRestriction 周期 recadre walk + FaceClassifier 兜底）；`topalgo/brep_top_adaptor/topol_tool_brep.rs`（BRepTopAdaptor_TopolTool 653 行 1:1：逐边 Curve2d 列表/迭代器、Classify/IsThePointOn 委托 FClass2d、buc60462 差分钳位 ComputeSamplePoints、Has3d/Tol3d/Pnt）。连带：FaceShapeSource 补 wire 注册（FaceExplorer 枚举前置）；BRepTool 语义落地 face UV 域优先（TFaceData.uv_domain，BRepAdaptor_Surface::FirstUParameter 语义）。锚点 4 组：Curve2d 参数化/edge-face、FClass2d In-Out-On(TestOnRestriction)/无限点、TopolTool 四边迭代+Classify+10x10 钳位+SamplePoint、顶点迭代器+Parameter+Tol3d/Pnt。**遗留 1 个 #[ignore]**：Perform 边界点兜底 On——根因已定位（FaceShapeSource 传 kernel locations 表 vs edge_pcurve_on_face 期望 bop DS 表），见头部快照）
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

## 8. Session 交接记录（2026-09-04，2a-2 全部 + 2a-3①②③④ 完成，2a-3 全部关闭）

**接续提示词：读本文档 §0 + §7 + §8 + 头部快照，从 §7 第一个未勾选项（2a-4 TopolTool 真实版）继续。**

### 本 session 提交链（rcad 子模块 sd-hash-wip，全部零回归）

| commit | 内容 |
|---|---|
| `57ab4ff5` | **2a-2 完成**：IntImpParGen_Intersector.gxx（824）从 geom2d_int.rs 的 GInter 硬编码副本泛型化为 `int_imp_par_gen.rs::Intersector<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>>`（+IntImpParGen.cxx 静态函数 + 泛化 MyImpParTool）；IntConicCurve 双泛型 = `int_conic_curve_gen.rs`（gxx 92 + lxx）+ `user_int_conic_curve_gen.rs`（889：5 ctor/5 Perform/5 InternalPerform）；geom2d_int.rs 净删 1134 行，GInter 实例化 = `impl<'a> ParTool<dyn Curve2dAdaptor + 'a>` marker（Geom2dCurveTool/TheProjPCurOfGInter）+ 具体薄壳包装 |
| `5cde756c` | **2a-3①**：Intf_InterferencePolygonPolyhedron.gxx（1368）→ `intf_interference_polygon_polyhedron.rs`（7 ctor/6 Perform/Interference±normal 双穿越/MKK 边界外扩/Intersect 5 参+9 参/sVertex-sEdge-dPiE 分类/KHROMOV 边界挠度分支/Extrema_ExtElC 边贴近尾段）+ Intf::PlaneEquation + ToolPolygon3d/ToolPolyh trait + TheInterferenceOfHInter alias |
| `f528da46` | **2a-3②**：Intf_InterferencePolygon2d.cxx（819）→ `intf_interference_polygon2d.rs`（自/双干涉、Clean 的 Only1Seg 怪癖与 decal 回卷、Intersect 的 parO/parT 1-based 槽 + sinTeta/rayIntf 相切带 + L718 封闭接缝 guard）；Polygon2dGen 实现 IntfPolygon2d trait |
| `3466868b` | **2a-3③**：IntCurve_IntPolyPolyGen.gxx（1797）→ `int_poly_poly_gen.rs`（Perform 双曲线公共版 + 自交版 + findIntersect + HeadOrEndPoint + GetIntersection 递归二分；**gxx L1390 的 Tol/TolConf 递归换位怪癖照抄**）；geom2d_int.rs 新增 GInterPolyTool<'a> 三合一工具 marker + TheIntPCurvePCurveOfGInter 具体壳；int_curve_curve_gen.rs 补 `intcurvcurv` 成员（GInter 第三个子引擎）并接通两个 unimplemented 臂 |
| `6b331eee`/`10ece47f`/`49bfb935`/`9a69c4a6` | 文档 runway 更新（2a-3④ 勘察清单在文档头部快照：Inter.pxx/InterUtils.pxx 函数级清单、依赖状态、④-a/④-b/④-c 切分） |

### 回归基线（全绿，以此为准）

algo lib **135**（120→125→130→135 逐批 +5 锚点）、kernel lib **645**、tkgeom_algo_gtests **134+1**、tkhelix **16**、pavefiller **26**。锚点测试全部在各自模块 `#[cfg(test)]` 内（tkgeom_algo_gtests.rs 是外部遗留修改区，勿动）。

### 本 session 新确立的翻译模式（沿例勿改）

1. **曲线类型也是模板参数**：`ParTool<C: ?Sized>` / `ProjectOnPCurveTool<C>` / `PCurveTool<C>` 把曲线类型做成 trait 泛型参数（OCCT 模板的 ParCurve/TheCurve）。**禁忌**：关联类型写 `type Curve = dyn Trait` 会固化 `+ 'static`，破坏 C++ `const Adaptor2d_Curve2d&` 的生命周期弹性（E0521）——正确做法是 marker 带 `PhantomData<&'a ()>` + `impl<'a> Trait<dyn Trait + 'a> for Marker<'a>`（见 GInterPolyTool<'a>），或引擎方法级泛型 + 具体壳内 helper-fn（`fn engine<'a>() -> Intersector<dyn Curve2dAdaptor + 'a, ...>`）。
2. **具体壳包装泛型引擎**：OCCT `_0.cxx` 具体实例化 = rcad 具体薄壳结构体（持有 base/域等无生命周期状态），方法内构造泛型引擎（helper-fn 技巧），perform 后拷回 base——递归状态在单次 perform 内保持（IntPolyPolyGen 的 RecursD1/D2 递归即此）。
3. **macro_rules 在 impl 块内的 hygiene**：宏体引用 `self`/局部变量（如 `composite`）必须在**函数体内**定义宏（int_curve_curve_gen.rs 的 finish! 先例），impl 层定义会 hygiene 失败。
4. **已踩坑**：bash 管道 `| tail` 会吞掉 python 的失败退出码 → 后续 `rm` 误删暂存文件（两次中招）——append+rm 必须分开两条命令并先验证；`python -c "..."` 内含反引号/`\U` 会转义地狱——大段脚本一律 Write 落盘再跑。

### 本 session 追加（2026-09-04，2a-3④ 完成，2a-3 全部关闭）

| commit | 内容 |
|---|---|
| （本轮代码提交） | **2a-3④**：IntCurveSurface_HInter 组装——新建 `geomalgo/int_curve_surface/`（mod：TransitionOnCurve/IntersectionPoint/IntersectionSegment/Intersection 数据类 + HCurveTool/HSurfaceTool/HInterHost trait + Adaptor3d{Curve,Surface}Basis；inter_utils.rs = InterUtils.pxx 全函数含 EstLim*/ProcessLinTorus/PerformCurveQuadric；inter_impl.rs = Inter.pxx 全函数；hinter.rs = HInter 壳 + 4 锚点；quad_curv_exact.rs = TheQuadCurvExactHInter + TheQuadCurvFunc + QuadricCurveExactInterUtils.pxx）。连带补齐：IntConicQuad 类壳 + Line×Quadric/Line×Plane/Ellipse×Quadric（int_patch/int_conic_quad.rs）、IntAna_IntLinTorus（int_patch/int_lin_torus.rs）、ElSLib 反向 Parameters（Plane/Cylinder/Cone/Sphere/Torus + normalizeAngle，kernel math/el.rs）、ThePolygonOfHInter/ThePolyhedronOfHInter 工具泛型构造器 new_tool*（int_curv_surf.rs，工具经 CurveEvalByTool/SurfaceEvalByTool 适配器走 landed 的 Init 路径）、PolyhedronLike/PolygonLike impl、ZerCSParFunc/IntCS 的 S/C 放宽 ?Sized |

**本 session 新确立的翻译模式（沿例勿改）：**

1. **回调束 = Host trait**：OCCT HInter.cxx 传给 Inter.pxx 模板的 lambda 束（ResetFunc/PerformBoundsFunc/AppendFunc/...，均捕获 this）在 Rust 编码为 `HInterHost<C, CT, S, ST>` trait（一方法一 OCCT 成员函数，`done`/`myIsParallel` 经 `done_flag()/is_parallel_flag()` 返回 &mut bool），HInter 实现之；引擎函数保持 OCCT 模板形状（泛型 + 逐条语句 1:1），经 host 回调回 HInter 成员——借用安全且无递归（每条引擎路径只回调"其它"成员）。
2. **同一工具类满足多个模板概念**：HCurveTool/HSurfaceTool 与 IntImp 的 CurveTool3d/PSurfaceTool 是同一 OCCT 工具类（TheHCurveTool/Adaptor3d_HSurfaceTool）的两个 Rust trait 面；`HInterHost` 的泛型参数直接带上双重约束，HInter 的 impl 处统一满足。engine 内所有工具调用点用 `<CT as HCurveTool>::`/`<ST as HSurfaceTool>::` 显式限定消歧。
3. **OCCT 隐式转换的落点**：Inter.pxx 的 `IntAna_IntConicQuad(theLine, gp_Cylinder)` 依赖 gp_Cylinder→IntAna_Quadric 隐式转换——Rust 侧显式 `Quadric::from_cylinder(...)`；`Perform(E, Pln, Tolang, Tol)`（Ellipse/Parab/Hypr×Plane）在 OCCT 内部转投 Quadric 路径（IntConicQuad.cxx L562-575），照抄（new_ellipse_plane 等 → perform_*_quadric）。
4. **采样类的工具泛型构造器**：landed 的 ThePolygonOfHInter/ThePolyhedronOfHInter 构造器收 `&dyn CurveEval/SurfaceEval`（等价 OCCT GeomGridEval_Surface 角色）；HInter 引擎从 adaptor 构造采样类时经 `CurveEvalByTool/SurfaceEvalByTool` 适配器（实现 CurveEval/SurfaceEval，采样走工具 trait）复用同一 Init 路径——OCCT 的 PolygonUtils/GeomGridEval 本就是模板分派。
5. **坑**：bash heredoc 在本环境再次截断（第 3 次）——大段追加一律 Write 落盘 + 单独一条 `cat >>`/`mv` 命令；`Quadric` 是 Clone 非 Copy——模板参数按 const& 传值的构造在循环内需 `.clone()`。

### 本 session 追加 2（2026-09-04，2a-4 Adaptor3d 层完成）

| commit | 内容 |
|---|---|
| （本轮代码提交） | **2a-4 Adaptor3d 层**：`topalgo/adaptor2d/line2d.rs`（Adaptor2d_Line2d 444 行 1:1，实现 Curve2dAdaptor）+ `topalgo/adaptor3d/hvertex.rs`（Adaptor3d_HVertex）+ `topalgo/adaptor3d/topol_tool.rs`（Adaptor3d_TopolTool 1697 行 1:1，含 Analyse/BSplSamplePnts/GetConeApexParam）；kernel el.rs 补 elclib_line_parameter_2d；HSurfaceTool 升级（u_trim/v_trim = 同曲面收窄窗口、bezier/bspline、Surface 去 ?Sized）；HInter BSplineSurface-TopolTool 分支接通。锚点 +7 全绿（algo lib 146） |

**2a-4 续（BRep 层）勘察结论（下一 session 直接开工）：**

- OCCT 剩余体量：BRepTopAdaptor_FClass2d.cxx（859，BRepTopAdaptor_TopolTool::Classify 的同包依赖）+ BRepTopAdaptor_TopolTool.cxx（653）+ BRepTopAdaptor_HVertex（250，Resolution 的 refine 算法）+ BRepAdaptor_Curve2d（143，= Geom2dAdaptor_Curve + pcurve + Edge/Face 访问器）。~1.9k 行。
- kernel `BRepTool` trait（topods.rs L1823 起，impl for BRep）已备齐 Shape 级访问器：parameter_on_edge（= BRep_Tool::Parameter(V,E,F)，vertex_params 按 vertex ptr_id 键）、curve_on_surface（= CurveOnSurface，pcurve 键 = face ptr + 预除位置）、is_edge_closed_on_face、is_edge_degenerated、vertex_tolerance/position、face_surface、u/v_resolution、edge_curve_world。
- BRepAdaptor_Curve2d = 持 (brep: &BRep, my_edge: Shape, my_face: Shape) + pcurve 委托到已有的 `impl Curve2dAdaptor for Curve2d`。
- BRepTopAdaptor_HVertex::Value 的 `return gp_Pnt2d(RealFirst(), RealFirst())`（do nothing 怪癖）照抄；Resolution 的三段 refine 全 1:1（面 surface 用 BRepTool::face_surface 的局部系，与 pcurve 同系自洽）。
- BRepTopAdaptor_FClass2d 的 ctor 算法与已移植的 IntTools_FClass2d（topalgo/brep_top_adaptor/fclass2d.rs）同骨架，但 (a) 以 kernel BRep+Shape 为源（非 DS 索引）、(b) 有 QuasiUniformDeflection 重离散循环（gxx L347-435，面积/周长比驱动）、(c) U1/V1/U2/V2 周期补偿。WireExplorer 需 Shape 级版本（现有 order_wire_edges 是 DS-based；可参照其 1:1 注释移植）。
- BRepTopAdaptor_TopolTool：Initialize(S) 的 down_cast<BRepAdaptor_Surface> 在 rcad = 直接持 (BRep, face_idx)；myCurves = 逐 edge 的 BRepAdaptor_Curve2d 列表；Classify/IsThePointOn 委托 FClass2d::Perform/TestOnRestriction；ComputeSamplePoints 的 Cylinder/Cone/Sphere/Torus 分支有 buc60462 的 nbsu/nbsv 差分钳位（30/15）照抄。
- 坑：bash heredoc 截断第 4 次（长 python 一律 Write 落盘再 uv run）；`tktopalgo_gtests.rs` 存在先于本 session 的编译损坏（Trsf form/scale 字段缺失，commit 4a823e11 引入面），与 TKHLR 线无关、勿混入。

### 本 session 追加（2026-09-05，2a-4 BRep 层完成，2a-4 全部关闭）

| commit | 内容 |
|---|---|
| （本轮代码提交） | **2a-4 BRep 层**：`topalgo/brep_adaptor/curve2d.rs`（BRepAdaptor_Curve2d）+ `topalgo/brep_top_adaptor/hvertex_brep.rs`（BRepTopAdaptor_HVertex，Resolution 三段 refine）+ `topalgo/brep_top_adaptor/fclass2d_topol.rs`（BRepTopAdaptor_FClass2d 859 行 1:1）+ `topalgo/brep_top_adaptor/topol_tool_brep.rs`（BRepTopAdaptor_TopolTool 653 行 1:1）。连带：FaceShapeSource 补 wire 注册（FaceExplorer 枚举前置）；BRepTool 语义落地 face UV 域优先（TFaceData.uv_domain）。锚点 4 组全过 + 1 个 #[ignore]（见下）。基线 algo lib 150+1i |

**本 session 新确立的翻译事实（沿例勿改）：**

1. **面 UV 域语义**：BRepAdaptor_Surface::FirstUParameter 等返回面的 UV 限界域（wire 限界），rcad 落地为 TFaceData.uv_domain 优先、SurfaceEval::default_domain 兜底——底层平面自然域是 ±∞，直接用会让 DomainIsInfinite/Classify 全错。
2. **fixture 顶点方向**：OCCT BRep_Builder::Add(E,V) 在边上存顶点方向（起点 FORWARD、终点 REVERSED）；BRepTools_WireExplorer 的 V1/V2 链依赖它。rcad fixture 用 add_tedge 时必须给终顶点传 Reversed 克隆，否则 wire 走断（ordered 只剩 1 条边）。
3. **遗留对齐项（唯一 #[ignore]，已定位根因）**：`fclass2d_topol_perform_boundary_on`——OCCT 边界点走 SiDans==0 → BRepClass_FaceClassifier → On。rcad 兜底返回 In 的根因：FaceShapeSource 把 kernel BRep locations 表（identity 无槽位 0）传给 `edge_pcurve_on_face`（该函数按 bop DS 表约定"槽 0 = identity"读）→ pcurve 键全部落空 → FaceExplorer::segment=None → nowires → In。对齐表约定后 un-ignore。
4. **坑**：bash heredoc 大文件追加必截断（第 4 次）；对既有子系统加临时 eprintln 探针后必须立即清除（本次 face_explorer.rs 探针已全部清除）。

### 下一 session 入口（按序）

1. **2b**：Contap 包 1:1（SurfProps→HContTool→HCurve2dTool→Point/Line→SurfFunction→ArcFunction→ContAna→Contour 2389 行）。Contap_TheSearch 的 domain 参数 = BRepTopAdaptor_TopolTool（已就绪）。
2. 顺手项：对齐 FaceShapeSource/BRepTool 的 locations 表约定，un-ignore `fclass2d_topol_perform_boundary_on`。
3. Stage 3b 的 HLRBRep_CInter/The*OfInterCSurf 实例化只需为 HLRBRep_Curve/HLRBRep_Surface 实现同一套 HCurveTool/HSurfaceTool（+ CurveTool3d/PSurfaceTool）trait。

