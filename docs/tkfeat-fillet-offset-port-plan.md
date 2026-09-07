# TKFeat / TKFillet / TKOffset 1:1 翻译推进计划（多 session runway）

> **调研快照（2026-09-06 勘察，基于 OCCT 8.0.0 源码树实测）**
> - **目标**：把 TKFeat（特征）、TKFillet（圆角/倒角）、TKOffset（抽壳/偏移）三个 toolkit 按布尔管线相同的"逐行形式对齐"方法论移植进 rcad（`feat/`、`fillet/`、`offset/`），布尔能力一律走 rcad 已对齐的 TKBO 等价层，不引入老布尔。
> - **核心事实（决定 whole 计划走向）**：OCCT 8.0 已把 `BRepAlgoAPI` 整包移入 **TKBO**（TKBO = BOPAlgo / BOPDS / BOPTools / **BRepAlgoAPI** / IntTools）。TKBool 只剩：老布尔求交器 `TopOpeBRep`、`TopOpeBRepBuild`/`TopOpeBRepDS`/`TopOpeBRepTool`、`BRepAlgo`（工具门面）、`BRepFill`（扫掠/放样，内容与布尔无关）、`BRepProj`。
> - **逐包裁决**：TKFeat 布尔调用 100% TKBO（`BRepFeat_Builder` 直接继承 `BOPAlgo_BOP`，`BRepFeat_MakeCylindricalHole` 五处 `BOPAlgo_BOP::Perform`，另用 `BRepAlgoAPI_Cut×7/Fuse×4/Section×2/Common×2` + `BOPAlgo_BuilderFace`），**零 TKBool include**；TKOffset 的 `BRepOffset` 已完全 BOPAlgo 化（`BOPAlgo_PaveFiller`（BRepOffset_Tool.cxx:1475）、`MakerVolume`（MakeOffset.cxx:5071/5159/5214）、`Splitter/Section/BOP/BuilderFace`（MakeOffset_1.cxx:42-45）、`BRepAlgoAPI_Check`）；TKFillet 的 ChFi3d 体系写在 `TopOpeBRepDS` + `TopOpeBRepBuild_HBuilder`（仅 7 方法面：Perform/MergeSolid/IsSplit/Splits/Merged/NewFaces/NewEdges，成员 `myCoup`，ChFi3d_Builder.hxx:180 `Builder()` 访问器）上，**但老布尔求交器 TopOpeBRep 零引用**（TopOpeBRepBuild/DS/Tool 内部也不引用它）。
> - **已确认决策（用户拍板，2026-09-06）**：① "不依赖 TKBool" 精确化为"不引入 TopOpeBRep 求交器与 BRepProj；TopOpeBRepDS / TopOpeBRepBuild(HBuilder 子集) / TopOpeBRepTool(子集) / BRepAlgo(工具类) / BRepFill(按需) 作为伴随翻译的非布尔支撑"；② 顺序 = **fillet → offset → feat**；③ BRepFill **按需引入、逐类对齐**；④ TKBO 缺口（BOPAlgo_Section / BOPAlgo_MakerVolume）作为 Stage 0 前置先补齐。详见 §8。
> - **规模**：三 toolkit 本体 ~148.8k OCCT 行（TKFillet ~79.7k / TKOffset ~41.4k / TKFeat ~27.7k）+ 伴随件按需。rcad `fillet/` 已有 ~17.3k 行翻译但**整个模块 TEMP-EXCLUDED**（topopebrepbuild.rs 编辑中途），`feat/`、`offset/` 空占位未声明。
> - **下一 session 入口**：§7 进度表 Stage 0.1（fillet WIP 收尾恢复编译）。
> - 回归基线（沿用 tkhlr-port-plan 口径，全绿为准）：algo lib 361/0、kernel 663、builder_stage_tests 76 + smoke 1、pavefiller_stage_tests 26、boolean 网格（bopfuse 748/0、bopcommon 755/0、boptuc 745/0、bcut 729/0 含 g6、splitter 12/12）。
> - 本文档体例沿 `tkhlr-port-plan.md`：§0 开场流程 → §7 进度勾选表从第一个未勾选项继续，完成一项勾一项。

## 0. 给后续 agent 的开场流程（每个 session 必读）

1. 通读本档 + 根 `AGENTS.md` + `rcad/docs/module-map.md`（三个模块的行在 module-map §3/§6）。
2. `cd rcad`（shell cwd 会被重置到根仓库，编译/测试前必须先 cd）；`cargo test -p rcad-algo --lib` 确认基线全绿，再动代码。
3. 看 §7 勾选表，从第一个未勾选任务继续；一次只推进一个子阶段，完成后勾选并在 §9 追加 session 记录。
4. 每个子阶段一次提交（rcad 子模块），提交信息注明 `Stage x.y：<内容摘要>`；提交前跑 §2 的回归集。
5. 纪律沿 AGENTS.md：阶段 1（翻译期）不跑测试不看结果、每函数 `cargo check -p rcad-algo`、`// OCCT <文件> L<起>-<止>` 锚点、代码注释全英文、单文件 <2000 行。

## 1. 规模与结构（勘察结论）

OCCT 根：`C:/Users/lilu/works/OCCT/src/ModelingAlgorithms/`（`$OCCT_SRC`）。

### 1.1 三 toolkit 本体（.cxx 行数）

| Toolkit | 包 | cxx | 行数 | 备注 |
|---|---|---|---|---|
| TKFeat | BRepFeat | 14 | 14,362 | 布尔面纯 TKBO；`Builder` 继承 `BOPAlgo_BOP` |
| | LocOpe | 21 | 13,342 | 用 BRepAlgo_Loop + BRepFill_Pipe/Evolved |
| TKFillet | ChFi3d | 19 | **36,532** | 老布尔 DS/HBuilder 消费方（唯一例外点） |
| | BRepBlend / BlendFunc / Blend | 57 | 27,115 | 纯曲面求解（滚动球逼近），无布尔 |
| | ChFiKPart | 14 | 5,389 | 标准例族（平面/球/柱/锥/环组合） |
| | ChFi2d | 7 | 4,645 | **2D 圆角，零布尔依赖** |
| | ChFiDS | 12 | 3,381 | Spine/Stripe 数据结构 |
| | FilletSurf | 2 | 1,266 | API 门面 |
| | Blend | 11 | 827 | |
| | BRepFilletAPI | 3 | 1,087 | MakeFillet/MakeChamfer/MakeFillet2d/LocalOperation |
| TKOffset | BRepOffset | 12 | **28,932** | 布尔面纯 TKBO（PaveFiller/MakerVolume/Splitter/Section/BOP） |
| | BRepOffsetAPI | 13 | 5,736 | 消费 BRepFill_* 七类 + BRepAlgo 工具 |
| | BiTgte | 3 | 3,218 | BRepAlgo_Loop/Image/AsDes |
| | Draft | 6 | 3,470 | Draft_Modification → BRepOffsetAPI_DraftAngle |

合计 ~148.8k 行。类清单（翻译单元枚举）：
- **BRepFeat**：Builder Form Gluer MakeCylindricalHole MakeDPrism MakeLinearForm MakePipe MakePrism MakeRevol MakeRevolutionForm PerfSelection RibSlot SplitShape Status StatusError
- **LocOpe**：BuildShape BuildWires CSIntersector CurveShapeIntersector DPrism FindEdges FindEdgesInFace GeneratedShape Generator GluedShape Gluer LinearForm Operation Pipe PntFace Prism Revol RevolutionForm SplitDrafts SplitShape Spliter WiresOnShape
- **BRepOffsetAPI**：DraftAngle FindContigousEdges MakeDraft MakeEvolved MakeFilling MakeOffset MakeOffsetShape MakePipe MakePipeShell MakeThickSolid MiddlePath NormalProjection Sewing ThruSections

### 1.2 大文件榜（拆文件规划用）

| 文件 | 行数 | 去向 |
|---|---|---|
| BRepOffset_MakeOffset_1.cxx | **9,533** | Stage 2b，按 OCCT 函数组拆多 rs 文件 |
| ChFi3d_Builder_0.cxx | 5,944 | rcad 已有 chfi3d_builder_0.rs（2,007）续 |
| BRepOffset_MakeOffset.cxx | 5,659 | Stage 2b |
| ChFi3d_Builder_C1.cxx | 5,445 | rcad 已有 chfi3d_builder_c1.rs（1,651）续 |
| BRepOffset_Tool.cxx | 4,659 | Stage 2b |
| ChFi3d_Builder_CnCrn.cxx | 3,927 | Stage 1f 新增 |
| ChFi3d_Builder_2.cxx | 3,900 | Stage 1f 新增 |
| BRepFeat_RibSlot.cxx | 2,736 | Stage 3b |
| ChFi3d_Builder_6.cxx | 2,703 | Stage 1f 新增 |
| ChFi3d_FilBuilder.cxx | 2,635 | Stage 1f |
| BRepOffset_Inter2d.cxx | 2,394 | Stage 2b |
| ChFi3d_ChBuilder.cxx | 2,318 | Stage 1f 新增 |
| BRepFeat_MakeRevolutionForm.cxx | 2,030 | Stage 3b |
| LocOpe_SplitDrafts.cxx | 1,976 | Stage 3c |
| BRepOffset_Offset.cxx | 1,814 | Stage 2b |
| LocOpe_SplitShape.cxx | 1,776 | Stage 3c |
| LocOpe_WiresOnShape.cxx | 1,623 | Stage 3c |

### 1.3 伴随件（非布尔，按需引入）

| OCCT 包（物理位置） | 规模 | 谁需要 | 内容性质 |
|---|---|---|---|
| TopOpeBRepDS（TKBool） | 58 cxx / 20,664 | ChFi3d/ChFiKPart/BRepFilletAPI/FilletSurf | 老布尔**数据结构**（ChFi3d 自己填 DS） |
| TopOpeBRepBuild（TKBool） | 66 cxx / 38,095 | ChFi3d 经 HBuilder 7 方法 | 老布尔**构建**的可达子集 |
| TopOpeBRepTool（TKBool） | 36 cxx / 19,428 | 仅 TOOL 等小子集被 ChFi3d 引用 | 工具 |
| BRepAlgo（TKBool） | 7 cxx / 3,953 | BRepOffset/BRepOffsetAPI/BiTgte/LocOpe/BRepFill/ChFi3d | 工具类：AsDes(298)/Image(353)/Loop(1151)/FaceRestrictor(454)/NormalProjection(677)，NormalProjection 内部仅依赖 BRepAlgoAPI_Section（TKBO） |
| BRepFill（TKBool） | 31 cxx / 30,567 | BRepOffsetAPI、LocOpe（Pipe/Evolved） | 扫掠/放样，内部已是现代 BOP（PaveFiller×3、Builder/MakerVolume/BuilderFace、BRepAlgoAPI_Section×2 于 Draft/TrimShellCorner） |
| **不需要**：TopOpeBRep（19,775 求交器）、BRepProj | — | 零引用 | 永不引入 |

## 2. 验收标准（诚实边界）

1. **拓扑优先**：每个用例先对齐 V/E/F/S/SOLID 数量与曲面类型分布（step-topo-diff ✅），SA 后行；数值达 OCCT `checkprops` 5-6 位有效位。
2. **网格**（OCCT `$OCCT_SRC/tests/` 顶层用例数）：blend 10、chamfer 13、fillet2d 5（Stage 1）；offset 24、thrusection 10、draft 4（Stage 2）；feat 8、evolved 5、mkface 10、pipe 5、nproject 4（Stage 2/3）。另有 bugs/modalg_* 的 featprism/featdprism/featrevol/featrf/blend/chamf/offsetshape 用例按需纳入。
3. **永久排除**（不算待办）：offset/draft 等网格中 `restore *.rle` 数据驱动用例（沿 AGENTS.md 外部数据规则）；nproject 网格如依赖外部数据同此处理。
4. **GTests**：TKFillet 2（BRepFilletAPI_MakeFillet_Test / MakeChamfer_Test）+ TKOffset 3（BRepOffset_MakeOffset_Test / BRepOffsetAPI_MakeThickSolid_Test / MakePipeShell_Test）落位 `tests/tkfillet_gtests.rs` / `tests/tkoffset_gtests.rs`；Sewing_Test 归 shhealing 范畴不翻。TKFeat 无 GTests。
5. **红线**：不引入 TopOpeBRep 求交器与 BRepProj；不把 ChFi3d 合并层重写到 BOPAlgo（违反 1:1 方法论）；OCCT 没有的步骤/概念一律不加。
6. **每子阶段回归集**（提交前必跑）：`cd rcad && cargo test -p rcad-algo --lib` 全绿；boolean 网格基线（bopfuse/bopcommon/boptuc/bcut/splitter，run_grid.ps1 必须取最新编译产物）；builder_stage_tests + pavefiller_stage_tests 不回归。

## 3. 目录布局

```
libs/rcad-algo/src/
├── fillet/                  # TKFillet 全包 + TopOpeBRepDS/Build/Tool 子集（延续现状，目录即对齐边界）
│   ├── chfi2d_*.rs          # 1a 新增（AnaFilletAlgo/FilletAlgo/ChamferAPI/Builder/FilletAPI）
│   ├── chfi_ds.rs           # 1b 补全（现 2,284 行）
│   ├── chfi_kpart.rs        # 1c 扩充（现 397 行）
│   ├── topopebrepds.rs      # 1d 按需扩充（现 514 行；超 2000 行拆 topopebrepds/ 子目录）
│   ├── brep_blend*.rs       # 1e 新增（BRepBlend/Blend/BlendFunc）
│   ├── chfi3d*.rs           # 1f 延续六件 + 新增 builder_2/_6/_cncrn/_chbuilder/_notimp
│   ├── topopebrepbuild.rs   # 1g（现 2,107 行 WIP；超限拆 topopebrepbuild/ 子目录）
│   ├── brep_fillet_api.rs   # 1h（现 830 行盘点补全）
│   └── fillet_surf.rs       # 1h 新增
├── brep_algo/               # 2a 新建：OCCT TKBool/BRepAlgo 工具包（非布尔求交，module-map 挂 TKBool/BRepAlgo 归属）
│   └── as_des.rs / image.rs / loop.rs / face_restrictor.rs / normal_projection.rs
├── brep_fill/               # 按需新建：OCCT TKBool/BRepFill（逐类引入，每类标首个消费方）
│   └── offset_wire.rs / compatible_wires.rs / generator.rs / draft.rs / trim_shell_corner.rs / pipe.rs / evolved.rs / ...
├── offset/                  # TKOffset 全包（现空占位）
│   ├── brep_offset_*.rs     # 2b（MakeOffset_1 按函数组拆多文件）
│   ├── bi_tgte.rs / draft_*.rs   # 2d
│   └── brep_offset_api_*.rs # 2e
└── feat/                    # TKFeat 全包（现空占位）
    ├── brep_feat_*.rs       # 3a/3b
    └── loc_ope_*.rs         # 3c
```

约定：新目录（brep_algo/brep_fill）与 feat/offset 在 `lib.rs` 声明后才编译；module-map.md 对应行同步更新归属与状态。

## 4. 阶段划分

### Stage 0：前置收口（~1 session）——两项硬前置

- **0.1 fillet WIP 收尾**：`fillet/topopebrepbuild.rs` 编辑收尾 → `lib.rs` L12 取消 `// TEMP-EXCLUDED`，恢复 `pub mod fillet;` 编译 → 对 fillet/ 现有 12 文件做对齐标记盘点（哪些函数已 `✅ OCCT-aligned`、哪些 WIP），缺口清单回填 §7。`algo_ext/mod.rs` L24/L33 再导出同步恢复。
- **0.2 BOPAlgo_Section 1:1**（OCCT TKBO/BOPAlgo/BOPAlgo_Section.cxx，414 行）：真 SECTION 语义（section edges + 可选 PCurve/SetApproxPCurve），替换 `brep_algo_api/mod.rs` 中 `BooleanOpType::Section → 退化为 cut` 的路径。消费方：BRepFill_Draft、BRepFill_TrimShellCorner、BRepFeat_MakeLinearForm、BRepFeat_MakeRevolutionForm + BRepAlgo_NormalProjection。配解析锚点单测。
- **0.3 BOPAlgo_MakerVolume 1:1**（OCCT TKBO/BOPAlgo/BOPAlgo_MakerVolume.cxx，414 行）：BRepOffset_MakeOffset 三处调用（MakeOffset.cxx:5071/5159/5214）。注意与 `builder_solid.rs`（BOPAlgo_BuilderSolid）是两个类，勿混淆。配锚点单测。

### Stage 1：TKFillet（~79.7k + 老布尔子集伴随，6-10 session）

按依赖序推进；每字母项 1-2 session，完成即跑 blend/chamfer/fillet2d 对应用例。

- **1a ChFi2d 热身**（4,645 行，零布尔依赖）：AnaFilletAlgo（解析 2D 圆角）/FilletAlgo/ChamferAPI/Builder/FilletAPI → `fillet/chfi2d_*.rs`。锚点 = 直线-圆/圆-圆 2D 圆角解析解单测；验收 fillet2d 网格 5 用例。
- **1b ChFiDS 盘点补全**（OCCT 3,381 行，rcad 已有 2,284）：Spine/Stripe/FilSpine/FaceInterference/FilBuilder 等逐类对齐标记，补缺。
- **1c ChFiKPart**（5,389 行，rcad 已有 397 起步）：19 个 cxx 标准例族（平面-平面/平面-柱/球/锥/环…）逐例 1:1。
- **1d TopOpeBRepDS + TopOpeBRepTool 子集**：以 ChFi3d 实际引用清单为准（HDataStructure/DataStructure/Surface/Curve/Curve/Point/Kind/Transition/Interference 族 + CurvePointInterference/SurfaceCurveInterference/SolidSurfaceInterference + PointIterator/InterferenceIterator/CurveExplorer + TopOpeBRepTool_TOOL），对 `topopebrepds.rs`（现 514 行）按清单 1:1 扩充；只翻引用到的，不整包扫。
- **1e BRepBlend / Blend / BlendFunc**（27,115 行，底层支撑数学）：滚动球曲面求解链；先盘点 `geomalgo/` 是否已有重叠实现，再按 ChFi3d 消费面引入（Blend_Function/Blend_Ruled/CSFunction/FilletPoint/BlendFunc_Chamfer/CorrdConic 等按消费清单）。
- **1f ChFi3d 主体**（36,532 行）：续现有六件（chfi3d.rs 3,745 / builder_0 2,007 / c1 1,651 / filds 1,290 / spkp 647），新增 Builder_2、Builder_6、CnCrn、ChBuilder（2,318）、NotImp/Debug 等文件；ChFi3d_Builder 主类（含 `myCoup` HBuilder 交互）补全。依赖 1b/1d/1e。
- **1g TopOpeBRepBuild HBuilder 链补全**：以 7 方法（Perform/MergeSolid/IsSplit/Splits/Merged/NewEdges/NewFaces）可达闭包为界续 topopebrepbuild.rs；只翻闭包内函数，不整包扫。
- **1h BRepFilletAPI 门面 + FilletSurf**（1,087 + 1,266 行）：MakeFillet/MakeChamfer/MakeFillet2d/LocalOperation 盘点补全（现 830 行）+ fillet_surf.rs 新增。验收：GTests 2 项 + blend/chamfer 网格。

### Stage 2：TKOffset（~41.4k + 按需 BRepFill/BRepAlgo，4-7 session）

- **2a BRepAlgo 工具包**（新建 `brep_algo/`，~4.0k 行）：AsDes → Image → Loop（1,151）→ FaceRestrictor → NormalProjection（677，依赖 0.2 的 SECTION）。
- **2b BRepOffset 核心**（28,932 行）：建议序 MakeSimpleOffset（最短闭合路径，先行烟囱）→ Offset/OffsetSurface（偏移曲面构造）→ Inter2d/Inter3d → Tool（4,659，含 PaveFiller 驱动段 Tool.cxx:1475）→ MakeOffset_1（9,533，按函数组拆多文件）→ MakeOffset（5,659）。前置盘点 `bop/int_tools/` 对 IntTools_Context/ShrunkRange/FaceFace/FClass2d/BeanFaceIntersector/Tools 的覆盖度，缺的按消费清单补（它们本身是 TKBO 层件）。
- **2c BRepFill 第一批按需**（`brep_fill/`，逐类）：OffsetWire（MakeOffset 的 wire 路径）→ CompatibleWires + Generator（ThruSections）→ Draft（BRepFill_Draft，用 0.2 SECTION）→ TrimShellCorner（MakeThickSolid）。PipeShell/Pipe/Filling/Evolved 等到对应 API 启用时再引入（2e/3c 触发）。
- **2d BiTgte + Draft 包**（3,218 + 3,470 行）：BiTgte_Blended/Contact/… → Draft_Modification 族 → 支撑 BRepOffsetAPI_DraftAngle。
- **2e BRepOffsetAPI 门面**（5,736 行，15 类按序）：MakeOffsetShape → MakeOffset → MakeThickSolid → ThruSections（rcad topalgo/thru_sections.rs ruled 路径迁移归位）→ DraftAngle → MakePipe/MakePipeShell → MakeFilling → MakeEvolved → MakeDraft → MiddlePath → FindContigousEdges → NormalProjection。Sewing 不翻（归 shhealing）。
- 验收：offset 网格非 .rle 用例 + thrusection 10 + draft 4；GTests 3 项。

### Stage 3：TKFeat（~27.7k + BRepFill_Pipe/Evolved 按需，3-5 session）

- **3a 架构映射先行 + 烟囱**：`BRepFeat_Builder` 继承 `BOPAlgo_BOP` → rcad 组合 + 委托 `Builder(my_operation)`（OCCT 继承链写进文件头注释，属"架构差异层"先消灭再对齐）；`BRepFeat_MakeCylindricalHole`（纯 BOPAlgo_BOP 驱动，5 处 Perform）作为最小闭合烟囱案例先行。
- **3b BRepFeat Form 家族**（14,362 行）：Form（基类 Perform/Propagate）→ MakeDPrism/MakePrism/MakeRevol/MakeRevolutionForm（2,030）/MakePipe/MakeLinearForm（两处用 0.2 SECTION）→ RibSlot（2,736 大头）→ Gluer/PerfSelection/SplitShape/Status/StatusError。
- **3c LocOpe**（13,342 行，21 类）：先 Operation/BuildShape/BuildWires/Gluer/GluedShape/GeneratedShape/PntFace/FindEdges/FindEdgesInFace 等小件 → Prism/DPrism/Revol/Pipe/LinearForm/RevolutionForm/Generator（用 2a 的 Loop）→ 大文件 SplitShape（1,776）/WiresOnShape（1,623）/SplitDrafts（1,976）/Spliter/CurveShapeIntersector/CSIntersector。LocOpe_Pipe/RevolutionForm 触发 `brep_fill/pipe.rs`、`evolved.rs` 引入。
- 验收：feat 8 + mkface 10 + evolved 5 网格；bugs/modalg 的 feat 命令族用例按需。

### Stage 4：验收闭环（1-2 session）

- **4.1 生成器接入**：`tools/occt-test-gen/src/main.rs` 接入 11 个网格翻译（blend/chamfer/fillet2d/feat/draft/evolved/mkface/thrusection/pipe/nproject/offset 非 .rle）。先例 = hlr 4c 的四个接入点（mod 声明/batch/分发/needs_x 判定）。⚠️ 前置：生成器当前引用已不存在的 `BooleanOp/BooleanOptions` 类型（与子模块 HEAD 漂移），接入时一并修复。
- **4.2 ref 基线**：`tools/gen-occt-ref/draw_ref_step.py` + `gen_ref_topology.py` 扩展命令族（blend/chamf/offsetshape/featprism/featdprism/featrevol/featrf/draft/thrusections/pipe/mkface…），产物入 `tests/occt/step_output/ref/` + `tests/occt/step_reference/`；DRAWEXE 注意 `-b` batch 与一次一进程（tkhlr 踩坑 17）。
- **4.3 逐用例对齐**：step-topo-diff 目录批量 → 汇总表 → 按 §2 标准逐格转绿。
- **4.4 收尾**：module-map.md 三行状态更新（feat/fillet/offset + brep_algo/brep_fill 新目录）；本档 §7 全勾，移交维护模式。

## 5. 已知坑（勘察实证）

1. **OCCT 8.0 BRepAlgoAPI 属 TKBO**——按 7.x 旧认知判依赖会把 BRepFeat/ 的调用点误记成 TKBool；EXTERNLIB.cmake 是权威。
2. **BRepFill 在 TKBool 目录但内容是扫掠/放样**，内部已走现代 BOP——按需引入不算"依赖 TKBool"，目录归属 ≠ 依赖语义。
3. **ChFi3d 求交自己算**（ChFiDS spine 数据 + ChFi3d 内部相交），老布尔只用于最终合并（HBuilder 7 方法）；TopOpeBRepBuild/DS/Tool 内部零引用 TopOpeBRep 求交器——所以老布尔"构建半边"要翻，"求交半边"永不翻。
4. **BRepFeat_Builder 继承 BOPAlgo_BOP**：Rust 无继承，用组合 + 委托映射，继承关系必须写进文件头注释；不许借机改架构。
5. **BRepOffset_MakeOffset_1.cxx 9,533 行**：按 OCCT 函数组拆多个 rs 文件（先例 chfi3d_builder_*.rs 拆法），单文件 <2000 行。
6. **rcad `BooleanOpType::Section` 当前退化为 cut**：BRepFill_Draft/TrimShellCorner/MakeLinearForm/MakeRevolutionForm/NormalProjection 五类调用点必须真 SECTION 语义——Stage 0.2 关闭前不得开 2a/3b。
7. **fillet 模块 TEMP-EXCLUDED + topopebrepbuild.rs 编辑中途**：Stage 0.1 不收口，Stage 1 全部无从谈起；子模块工作树另有非本任务遗留修改（沿 tkhlr 计划口径，提交时只 add 本任务文件）。
8. **生成器漂移**：tools/occt-test-gen 引用已不存在的 `BooleanOp/BooleanOptions`；网格测试 gitignored——Stage 4.1 一并修复，勿提前。
9. **.rle 永久排除**：offset/draft 网格大头是 restore 外部数据——见 AGENTS.md，不算待办、不要尝试生成。
10. shell cwd 重置根仓库：编译/测试前先 `cd rcad`；DRAWEXE 必须 `-b` batch、每命令族一次进程（tkhlr 踩坑 15/17 沿用）。
11. 对齐纪律沿 AGENTS.md 三不要 + 阶段 1 不跑测试；翻译模式沿 tkhlr 已确立项（gxx→trait、lxx 内联、`// OCCT <文件> L<起>-<止>`、每叶子配解析锚点、OCCT 字面 no-op 照抄）。

## 6. 依赖缺口总表（rcad 侧）

| OCCT 件 | rcad 现状 | 动作 |
|---|---|---|
| BOPAlgo_Section | 无（Section 退化为 cut） | Stage 0.2 补 1:1（414 行） |
| BOPAlgo_MakerVolume | 无（builder_solid.rs 是 BuilderSolid，另一类） | Stage 0.3 补 1:1（414 行） |
| IntTools_Context/ShrunkRange/FaceFace/FClass2d/BeanFaceIntersector/Tools | bop/int_tools/ 部分覆盖 | Stage 2b 前盘点，缺则按消费清单补 |
| TopOpeBRepDS 子集 | topopebrepds.rs 514 行 | Stage 1d 按引用清单扩 |
| TopOpeBRepBuild HBuilder 闭包 | topopebrepbuild.rs 2,107 行 WIP | Stage 0.1 收尾 + 1g 补全 |
| TopOpeBRepTool 子集（TOOL 等） | 未盘点 | Stage 1d 一并盘点 |
| ChFiDS | chfi_ds.rs 2,284 行 | Stage 1b 盘点补全 |
| ChFiKPart | chfi_kpart.rs 397 行 | Stage 1c 扩充 |
| ChFi2d | 无 | Stage 1a |
| BRepBlend/Blend/BlendFunc | 未盘点（geomalgo/ 或有重叠） | Stage 1e |
| BRepAlgo 工具 5 类 | 无 | Stage 2a 新建 brep_algo/ |
| BRepFill 按需类 | topalgo/thru_sections.rs（ruled 路径） | Stage 2c/3c 逐类引入 brep_fill/ |
| BRepFilletAPI | brep_fillet_api.rs 830 行 | Stage 1h 盘点补全 |
| FilletSurf | 无 | Stage 1h |
| BRepOffset/BRepOffsetAPI/BiTgte/Draft | offset/ 空（仅 thru_sections 沾边） | Stage 2 |
| BRepFeat/LocOpe | feat/ 空 | Stage 3 |

## 7. 进度勾选表（每 session 更新）

**Stage 0 前置收口**
- [ ] 0.1 fillet WIP 收尾：topopebrepbuild.rs 收尾 + lib.rs 恢复编译 + 12 文件对齐标记盘点
- [ ] 0.2 BOPAlgo_Section 1:1 + 锚点单测 + 替换 Section→cut 退化路径
- [ ] 0.3 BOPAlgo_MakerVolume 1:1 + 锚点单测

**Stage 1 TKFillet**
- [ ] 1a ChFi2d（4,645）热身 + fillet2d 网格 5 用例
- [ ] 1b ChFiDS 盘点补全
- [ ] 1c ChFiKPart 标准例族（5,389）
- [ ] 1d TopOpeBRepDS/Tool 子集（按 ChFi3d 引用清单）
- [ ] 1e BRepBlend/Blend/BlendFunc（27,115，按消费面）
- [ ] 1f ChFi3d 主体（36,532；续六件 + 新增 builder_2/_6/_cncrn/_chbuilder）
- [ ] 1g TopOpeBRepBuild HBuilder 7 方法闭包补全
- [ ] 1h BRepFilletAPI 门面 + FilletSurf；验收 blend 10 + chamfer 13 网格 + GTests 2

**Stage 2 TKOffset**
- [ ] 2a brep_algo/ 工具包 5 类（~4.0k）
- [ ] 2b BRepOffset 核心（28,932；MakeSimpleOffset 烟囱先行 → Offset/Inter2d/Inter3d/Tool → MakeOffset_1 → MakeOffset）
- [ ] 2c BRepFill 第一批（OffsetWire/CompatibleWires/Generator/Draft/TrimShellCorner）
- [ ] 2d BiTgte（3,218）+ Draft 包（3,470）
- [ ] 2e BRepOffsetAPI 门面 15 类按序；验收 offset 非 .rle + thrusection 10 + draft 4 + GTests 3

**Stage 3 TKFeat**
- [ ] 3a BRepFeat_Builder 架构映射 + MakeCylindricalHole 烟囱
- [ ] 3b BRepFeat Form 家族（14,362；Form → Prism 族 → RibSlot → Gluer/SplitShape）
- [ ] 3c LocOpe（13,342；小件 → Prism 族/Generator → SplitShape/WiresOnShape/SplitDrafts）；验收 feat 8 + mkface 10 + evolved 5

**Stage 4 验收闭环**
- [ ] 4.1 occt-test-gen 接入 11 网格（先修 BooleanOp/BooleanOptions 漂移）
- [ ] 4.2 ref 基线扩展（draw_ref_step/gen_ref_topology 命令族）
- [ ] 4.3 step-topo-diff 逐用例对齐转绿
- [ ] 4.4 module-map.md 更新 + 本档收尾

## 8. 决策记录

| # | 决策 | 结论 | 日期 |
|---|---|---|---|
| D1 | "不依赖 TKBool" 的精确定义 | 不引入 TopOpeBRep 求交器与 BRepProj；TopOpeBRepDS / TopOpeBRepBuild(HBuilder 闭包) / TopOpeBRepTool(子集) / BRepAlgo(工具类) / BRepFill(按需) 作为伴随翻译的非布尔支撑 | 2026-09-06 用户确认 |
| D2 | 三模块推进顺序 | fillet → offset → feat（与 fillet WIP 衔接最顺） | 2026-09-06 用户确认 |
| D3 | BRepFill 引入策略 | 按需引入、逐类对齐（每类标注首个消费方），不整包扫 | 2026-09-06 用户确认 |
| D4 | TKBO 缺口时序 | BOPAlgo_Section / BOPAlgo_MakerVolume 作为 Stage 0 前置先补齐并加单测 | 2026-09-06 用户确认 |
| D5 | 伴随件目录位置 | 新建 `brep_algo/`（OCCT TKBool/BRepAlgo 包）与 `brep_fill/`（OCCT TKBool/BRepFill 包）两个顶层目录，module-map 挂归属；ChFi3d 链的 TopOpeBRepDS/Build 子集留在 fillet/ 内（延续现状） | 建议值，开工如无异议照此执行 |

## 9. Session 交接记录

### Session 0 交接（2026-09-07：调研 + 计划落盘，零代码改动）

- **本 session 产出**：本文档全篇（§1 勘察数据 / §4 阶段划分 / §5 坑 / §6 缺口表 / §7 勾选表 / §8 决策 D1-D5）。未动任何代码、未提交；本档在 rcad 子模块为新增未跟踪文件。
- **开工 session 第一个动作**：`cd rcad && git add docs/tkfeat-fillet-offset-port-plan.md` + 提交"计划落盘"（先例 tkhlr `7011f938`）。只 add 本档——工作树还有非本任务遗留修改（见下），严禁一起提交。
- **裁决生效记录**：D1-D4 用户已拍板（2026-09-07）；D5（`brep_algo/` + `brep_fill/` 新顶层目录，TopOpeBRepDS/Build 子集留 fillet/）为建议值，Stage 2a/2c 开工时如无异议照此执行。
- **勘察硬数据存放位置**（防重复调研）：逐包规模与大文件榜 = §1.1/§1.2；伴随件与"永不引入"清单 = §1.3；TKBO 缺口（BOPAlgo_Section / BOPAlgo_MakerVolume 各 414 行）= §6 头两行；ChFi3d 的 HBuilder 7 方法面（成员 `myCoup`，`ChFi3d_Builder.hxx:180` 的 `Builder()` 访问器）= 头部快照；代表性调用点行号（MakeOffset.cxx:5071/5159/5214、BRepOffset_Tool.cxx:1475、MakeOffset_1.cxx:42-45、MakeCylindricalHole.cxx 五处 Perform）= 头部快照；OCCT 源码根 = `C:/Users/lilu/works/OCCT/src/ModelingAlgorithms/`（`$OCCT_SRC` 即 `C:\Users\lilu\works\OCCT`，注意不是 C:\tools 下的安装树）。
- **rcad 侧现状快照（2026-09-07 实测）**：`fillet/` 12 文件 ~17.3k 行（chfi3d.rs 3,745 / chfi_ds.rs 2,284 / topopebrepbuild.rs 2,107 WIP / chfi3d_builder_0.rs 2,007 / fillet.rs 1,769 遗留兼容层 / builder_c1 1,651 / filds 1,290 / brep_fillet_api.rs 830 / spkp 647 / topopebrepds.rs 514 / chfi_kpart.rs 397），无 todo!/unimplemented!；`lib.rs` L12 TEMP-EXCLUDED、`algo_ext/mod.rs` L24/L33 再导出注释中；`feat/`、`offset/`、`xmesh/` 空占位未在 lib.rs 声明；boolean 入口 `fuse/cut/common/cut21/splitter` 全在 `bop/brep_algo_api/mod.rs`（全走 PaveFiller 路径，无 Union 快速路径），`BooleanOpType::Section` 退化为 cut（Stage 0.2 关闭点）；splitter 走 `algo_ext/bool_ops_ext.rs::n_ary_partition`。
- **环境遗留（非本任务，勿动勿修）**：根仓库工作树有 `Cargo.lock` 修改、`tools/occt-test-gen/src/main.rs` 修改、`tools/occt-test-gen/src/hlr_translate.rs` 未跟踪（TKHLR session 遗留）；生成器引用已不存在的 `BooleanOp/BooleanOptions` 类型 + 网格测试 gitignored = 已知漂移（§5.8），Stage 4.1 才处理；`run_grid.ps1` 必须取最新编译产物（AGENTS.md 教训）。
- **下一 session 入口（Stage 0.1，按序）**：
  1. `cd rcad`，`cargo test -p rcad-algo --lib` 确认基线全绿（361/0）；
  2. `git add docs/tkfeat-fillet-offset-port-plan.md` 提交计划落盘；
  3. 收尾 `fillet/topopebrepbuild.rs` → `lib.rs` L12 恢复 `pub mod fillet;` + `algo_ext/mod.rs` L24/L33 再导出 → `cargo check -p rcad-algo` 恢复编译；
  4. 对 fillet/ 12 文件做逐函数对齐标记盘点（`✅ OCCT-aligned` / WIP 清单），缺口回填 §7 后转 0.2（BOPAlgo_Section）。
- **回归基线口径**：沿 tkhlr-port-plan 头部（algo lib 361/0、kernel 663、builder_stage_tests 76 + smoke 1、pavefiller_stage_tests 26）+ §2.6 的 boolean 网格基线（bopfuse 748/0、bopcommon 755/0、boptuc 745/0、bcut 729/0 含 g6、splitter 12/12）。
- 后续 session 交接格式（沿 tkhlr-port-plan §8 体例）：本 session 提交链 / 回归基线 / 新确立翻译模式 / 下一 session 入口。
