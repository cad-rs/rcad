# E3-W 交接：三域（TKFeat / TKFillet / TKOffset）翻译推进 —— 2026-09-12

> **一句话（追加 14 收尾态）**：tkfeat featlf 的 profile 有效性**五连破**——追加 13 记的"布尔分割面无效"
> **被翻案**（真因是 rcad 边顶点 tag 约定不统一，不是布尔输出质量）；`is_done` 门 **0/12 → 11/15 越过**、
> **全部 kernel panic 清零**。落地 **7 处 1:1 修复**（边顶点 tag 不变量 / `direct_children` 不改标签 /
> `MakeEdge::Init` reordonate+周期 / `MakeFace(Pln,W,true)` 的 `CheckInside` / `MakeFace(W)` 的 `FindSurface` /
> `Curve3::reversed_parameter` 10 变体分发 / `Geom_TrimmedCurve::Reverse` 与 `cc->Reverse()` 字面翻译）
> \+ **1 处架构修复**（跨池子图重编号 `renumbered_pool` / `adopt_subgraph_into`，消除 rcad 池 index 别名）。
> 八网格与六门槛**全程零回归**。剩余墙全在**下游**，且**新增一处最急的非终止**（`featrf_a1`，见 §0.2）。

## 0. 新 session 一句话提示词（直接粘贴）

> 读 `rcad/docs/e3w-handover-tkfeat-fillet-offset.md`（本交接：门槛实测值/提交链/三域队列/**§0.2 新增非终止待办**/配方/16 条坑）
> 与 `rcad/docs/tkfeat-fillet-offset-port-plan.md` §E3-W 追加 11–14（权威脉络）；
> 先 `cd rcad` 跑 6 条门槛（见 §1）确认 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**，
> 再 `cd /c/Users/lilu/works/rcad-pro && bash output/run_eight_grids.sh` 确认**八网格 8/8**，
> 然后按 §4 队列开工——**第 0 项（最急）= 用 §5 配方 3 的 `occt_bool_runner` 逐轮对拍
> `BRepFeat_MakeRevolutionForm::init` 的 `while(!FirstOK)` 循环（cxx L724-853），消掉 `featrf_a1` 的不终止**
> （rcad 结构与 OCCT 逐行一致，差异必在 `cc` 重建 / `theLastPnt` 推进 / `myTol` 判据的数值上）；
> 第 1 项 = `BOPAlgo_Section` 对共面面片对的输出（`FalseSide ×5` 的根因，见 §4.1）；
> 第 2 项 = `BRepFeat_RibSlot::ExtremeFaces`（`NoExtFace ×3`）；第 3 项 = 追加 13 批次 2 的
> `Geom2dInt_GInter` 通用批；第 4 项 = tkoffset 三件 / tkfillet a1 上游。
> **铁律（不可动摇）**：非 TKBool 模块的代码与修复**一律严格 1:1 翻译对齐**——逐行语句对照 + OCCT 行号锚点 +
> 禁载体（GAP 只许留给外部未翻依赖且保留 OCCT 失败路径）/ **禁等价替换**（本轮 `cc->Reverse()` 的"换基曲线保留 [f,l]"
> 就是反面教材）/ 禁凑结果 / 架构差异先消灭再对齐。
> 每 Edit 后 `cargo check`；**跑非门槛网格一律加 `timeout`**（见 §6 坑 16）；探针即用即清
> （提交前 `git diff | grep -c "+.*eprintln"` = 0）；完成后更新 §E3-W 追加 15 并提交两仓库（**均不推送**）。

> **⚠ 追加 15 对本节第 1–2 条的覆盖（2026-09-12）**：本节原"第 0 项 = 用 `occt_bool_runner` 逐轮对拍 `init` 的 `while(!FirstOK)`"
> **已完成，且结论翻案**——差异**不在数值**，在**一处 OCCT 没有的重算语句**（见 §0.3）。原第 1–4 项顺延；
> **新第 0 项 = `BRepFeat_RibSlot::LFPerform` 的结果装配**（见 §4.0）。其余铁律与操作协议（§5 / §6）不变。

## 0.1 追加 14 本轮落地（2026-09-12，rcad `3bf6858f` / `cc38f690` / `7dfa5884`，均未推送）

| 修复 | 位置 | OCCT 锚点 | 消除的症状 |
|------|------|-----------|-----------|
| 边顶点 **tag 不变量** | `kernel/topo/topods.rs::add_tedge` | `BRepLib_MakeEdge.cxx L771-772` + `BRepPrim_Builder.cxx L143-155` + `TopExp.cxx L214-252` | `Edge InvalidPointOnCurve`、wire `NotClosed` |
| `direct_children` 不改写标签 | `topalgo/brep_check/brep_check_result.rs` | `TopoDS_Iterator(E)` 语义 | `orv` 落错参数（analyzer 误报） |
| `MakeEdge::Init` reordonate | `feat/brep_feat_rib_slot.rs::reorder_edge_endpoints` | `BRepLib_MakeEdge.cxx L603-651` | 顶点与 range 错位 |
| `MakeFace(Pln,W,true)` 的 `CheckInside` | `feat/brep_feat_rib_slot_b.rs`（+ `..make_revolution_form.rs`） | `BRepLib_MakeFace.cxx L262-272 / L905-924` | `Face BadOrientationOfSubshape` |
| `MakeFace(W)` 曲面探测 | `feat/brep_feat_rib_slot_b.rs::make_face_wire` | `BRepLib_MakeFace.cxx L189-262` | 无曲面面 ⇒ `NoSurface`（7 例 `NoFaceProf`） |
| **跨池子图重编号** | `feat/brep_feat_form_2.rs::renumbered_pool` / `adopt_subgraph_into` | 架构差异（OCCT 按指针携带图） | `Shape N is not a Face/Edge`（12 例 panic） |

**效果**：featlf 真实断言 15 例中 `is_done` 门通过 11（此前 0）；kernel panic 0（此前 12）。
**剩余失败全部在下游**：`FalseSide ×5`、`NoExtFace ×3`、`NoFaceProf ×3`、
`BRepTools_Modifier::Perform` GAP（a3）、`draft_modification_1_b.rs:1269` GAP（b4）、
off-chain 的 Draft（depouille）e4/e5。

## 0.2 追加 14 补记：`Geom_Curve::ReversedParameter` 链（同 session，**含一处新增的非终止待办**）

**修复（1:1，OCCT 锚点齐）**：
- `impl CurveEval for Curve3` **未覆盖 `reversed_parameter`** ⇒ 落 trait 默认**恒等**。OCCT 的
  `Geom_Curve::ReversedParameter` 是逐类型虚函数：`Geom_Line` → `-U`（Geom_Line.cxx **L163**）、
  `Geom_Circle`/`Geom_Ellipse` → `2π−U`（Geom_Circle.cxx **L184**、Geom_Ellipse.cxx **L199**）、
  `Geom_TrimmedCurve` → 委托基曲线（Geom_TrimmedCurve.cxx **L88-91**）。**rcad 已有多个消费者一直拿到恒等值**
  （`brep_fill_section_law.rs:267` 的 `reversed_parameter_of`、`chfi3d_perform_elspine.rs:1069`、
  `brep_fill/generator.rs:914`、`pave_filler.rs:5415`）。已补分发（`kernel/src/geom/eval.rs`）。
- `geom_curve_reversed` 的 `Trimmed` 臂原来只做 `(basis_reversed, t.last, t.first)` 槽位互换、
  **不套 `ReversedParameter`** ⇒ 修剪区间反向（`first > last`）⇒ `extrema.rs` 数值回退臂
  `f64::clamp(t0,t1)` panic（`min = 0.0, max = -1.0`，featrf_a1）。按
  `Geom_TrimmedCurve::Reverse`（cxx **L78-84**）修正。
- `BRepFeat_MakeRevolutionForm::init` 里 cxx **L741-744** 的 `cc->Reverse()` 改为**字面翻译**
  （原实现是"换基曲线但保留 `[f,l]`"的近似——属 AGENTS.md 禁止的"等价替换"）。

**★ 新增待办（下一 session 第一优先，比 `FalseSide` 更紧急）：`featrf_a1` 从不终止**
（复现：`cd /c/Users/lilu/works/rcad-pro && timeout 120 cargo test -q -p occt-generated-tests --test generated_occt_boolean_featrf feat_featrf_a1::`）。
上条修正把 a1 从 clamp panic 推到 `while (!FirstOK)` 循环（OCCT cxx **L724-853**）**空转**。
临时探针实测（已清）：`it_idx` 序列 `0 → 1 → 0 → 0 → …`、`counter1` 单调递增、`last_ok` 每轮为真
⇒ 触发 `it.Initialize(myListOfEdges)` 回卷 ⇒ 永不收敛；`first_ok` 恒假。
**下一手 = OCCT 侧逐轮对拍**（§5 配方 3：`occt_bool_runner` + 在 `BRepFeat_MakeRevolutionForm.cxx` L724-853
插桩打 `it`/`LastOK`/`FirstOK`/`theLastPnt` 轨迹）。rcad 的循环结构与 OCCT 逐行一致，
差异只可能在 `cc` 重建 / `theLastPnt` 推进 / `myTol` 判据的**数值**上。
另记一处**无害命名偏差**：OCCT L824 `theFEdge = edg;`，rcad 写作 `the_l_edge = edg;`——OCCT 侧两者皆死存储，不影响行为。

## 0.3 追加 15 本轮落地（2026-09-12，rcad 提交链见 §3；均未推送）

**一句话**：`featrf_a1` 的**不终止已清零**——根因**不是数值**，而是 `while(!FirstOK)` 里一处 **OCCT 没有的重算语句**；
顺带把同一条链后面的**三道墙**真身化（Transform 载体 / 滑动 profile 面的 outer wire / 子形状枚举的空占位）。
**六门槛与八网格全程零回归**（415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1；八网格 375·378·379·373·12·102·83·110）。

| 修复 | 位置 | OCCT 锚点 | 消除的症状 |
|------|------|-----------|-----------|
| **`fp`/`lp` 重算删除**（本轮根因） | `feat/brep_feat_make_revolution_form.rs::init` | `BRepFeat_MakeRevolutionForm.cxx` **L734-853**（`cc->Reverse()` 在 **L761**，`fp/lp` 只在 **L746-747** 算一次） | **`featrf_a1` 挂死**（`sens==2` 的边永远置不上 `FirstOK`，链在环上回卷） |
| **Transform 载体真身化** | 同上（`brep_builder_api_transform` 载体 → `topalgo::brep_builderapi_transform`） | Perform **L1204-1205 / L1245** | `panic("GAP… rotation engine pending")` |
| **同 trsf 多次 Perform 去重**（新 `perform_shape_once`） | `topalgo/brep_builderapi_transform.rs` | 一个 `gp_Trsf T`（**L1202-1204**）喂两次 Perform（**L1205/L1245**），且滑动路径 `mySkface = myPbase = Prof`（Init **L1107-1108**） | 同一 TShape **双重旋转**（结果几何错） |
| **滑动 profile 面的 outer wire** | `feat/brep_feat_make_revolution_form.rs::init` | Init **L912-921**（`MakeFace(f, myPln, 0.)` + `BB.Add(f, w)`：空面的第一条 wire = outer） | `BRep_Tool::Parameters` 的 **`Standard_NoSuchObject`**（面被喂了 null 占位当顶点） |
| **子形状枚举跳过空占位** | `brep_algo/tool.rs::sub_shapes`（Face 分支，**类型判据**） | `TopoDS_Iterator.cxx` **L28-51**（走 `TShape::myShapes` = 只有 Add 过的形状） | 扫掠迭代器把 null 占位当子形状 |
| **face 子形状按类型分派** | `brep_sweep/brep_sweep_builder.rs::add`（Face 分支） | `BRep_Builder::Add(F, W)` / `(F, V)` 的单列表语义 | 顶点可能占/混入 wire 槽 |

**实测（`generated_occt_boolean_feat_featrf`，真实断言 5 例）**：`0/5` 不变，但**全部 2s 内失败、无挂死**。
**a1 推进最深**：init 通过 + `IsDone()` 真 + `perform` 走完 Transform/扫掠，停在**结果装配**墙
（`surface area: expected 109.511, got 0` ⇒ `BRepFeat_RibSlot::LFPerform` 产空结果）。**a4/a5/a7/a9** 仍停在 init 的 `assert!(is_done)`。

**★ 本轮新发现两处卫生问题（需单独立卡）**：
1. **featrf 的参考拓扑断言一直被静默跳过**：测试找 `step_reference/occt_boolean_featrf_a1.json`（**不存在**），
   生成器实际写的是 **`occt_boolean_feat_featrf_a1.json`**（多 `feat_` 前缀）⇒ `exists()` 恒假 ⇒ 该网格**从未校验** V/E/F/S。
   owner = `occt-test-gen` 的 grid 命名。A1 参考值：**V10/E17/F9/SHELL1/SOLID1**、PLANE 5 + CYL 3 + CONICAL 1、面积 **109.511**。
2. **已提交的遗留探针**：`brep_sweep/num_linear_regular_sweep.rs:284` 的 `eprintln!("[ENTRY] …")`（来自 `ddc8d551`）——跑扫掠的测试都会刷屏，待清。

## 1. 门槛（本轮终测，全部实测；2026-09-12）

| 基线 | 值 | 命令 |
|------|-----|------|
| rcad-algo lib | **415/0/0** | `cargo test -p rcad-algo --lib` |
| rcad-kernel lib | **689/0** | `cargo test -p rcad-kernel --lib` |
| tktopalgo_gtests | **36/36**（本轮从 35/1 转绿） | `cargo test -p rcad-algo --test tktopalgo_gtests` |
| pavefiller_stage | **26/26** | `cargo test -p rcad-algo --test pavefiller_stage_tests` |
| builder_stage | **76/76** + smoke **1/1** | `cargo test -p rcad-algo --test builder_stage_tests` / `..._smoke` |
| 八网格 | **8/8**（375 · 378 · 379 · 373 · 12 · 102 · 83 · 110） | `cd /c/Users/lilu/works/rcad-pro && bash output/run_eight_grids.sh` |

**本轮已修 / 仍挂的既有失败**：
- 已修：`tktopalgo_gtests::occ10006_loft_and_fusion`（loft VIso GAP，见 §3 第 5 条）。
- **仍挂（off-gate，需单独立卡）**：`tkg3d_gtests` 3 例（`plane_default_domain_open` /
  `cylinder_default_domain` / `plane_default_domain_infinite`，default-domain 口径）。
- 无关既有破损（stash 验证先于本轮，**不在三域范围，勿顺手改**）：`rcad-constraints`
  的 `WireEdge` 字段/陈旧 import、`rcad-gordon` 两个 `Circle3` 测试字面量、`rcad-render`
  的 `use rcad_algorithms::...`（该 crate 未列入 workspace members）。

## 2. 域网格基线（off-gate，本轮实测；口径见 §5 陷阱 6）

组文件里每个 case 有两条测试：`*_geometry_loads`（占位，恒过）与
`*_draw_script_rcad_equivalent`（真实断言）。**下列"通过数"只数真实断言**：

| 网格 | 真实断言 | 状态 |
|------|---------|------|
| fillet2d_fillet2d | **10/10** | 全绿 |
| fillet2d_chamfer2d | **2/2** | 全绿 |
| mkface_after_offset | **4/4** | 全绿 |
| mkface_after_extsurf_and_offset | **32/32** | 全绿 |
| blend_simple | **0/11** | 全败（a1/a3/q1 早期待过的三例现亦败，与 E3-U 记档一致，非本轮回归） |
| blend_complex | **0/2** | 全败 |
| feat_featlf | **0/15** | 全败；`is_done` 门**11/15 已越过**（追加 14），**kernel panic 清零**，剩余全部为下游墙（§4.1） |
| feat_featprism / featrf | **0/6** / **0/5** | 全败 |
| feat_featrevol | **0/45** | 全败 |
| offset_shape_type_a / _i / _i_c | **0/1** / **0/12** / **0/19** | 全败 |
| offset_faces_type_i | **0/8** | 全败 |
| draft_angle | **0/49** | 全败 |
| thrusection_specific | **0/26** | 全败（thrusection 另见 §4.4） |

- `mkface_after_revsurf_and_offset` **只有 per-case 文件、无 merged target**（当前生成器产物），
  需要时按 §5 配方重生成该组。
- `chamfer_dist_angle/dist_dist/equal_dist` 等 `generated_occt_chamfer_*` 是**生成器中间产物**
  （`_begin`/`_cases_list` 形态），不是可跑网格。

**测量配方**（逐组拿真实断言数）：
```bash
cd /c/Users/lilu/works/rcad-pro
for g in blend_simple blend_complex fillet2d_fillet2d fillet2d_chamfer2d \
         feat_featlf feat_featprism feat_featrevol feat_featrf \
         offset_shape_type_a offset_shape_type_i offset_shape_type_i_c offset_faces_type_i \
         draft_angle; do
  r=$(timeout 900 cargo test -q -p occt-generated-tests --test generated_occt_boolean_${g} 2>&1 | grep "test result" | head -1)
  echo "${g}: ${r:-TIMEOUT/NO-TARGET}"
done
```
（`PASS = 恒过的 geometry_loads 占位数`；**FAIL 才是真实断言失败数**。）

## 3. 提交链与落地内容（**三轮**：追加 13 / 追加 14 / **追加 15 = 本轮**；均未推送）

**追加 15（本轮，2026-09-12）**——rcad `main`：**本交接文件所在提交**（docs：追加 15 + 本交接）
← `d280a867`（code：`featrf_a1` 链四处修复）← `e161078a`（= 追加 14 链尾）。
根 `main`：**对应的 rcad pointer sync 提交**（`sync: rcad pointer (<本交接所在 rcad 提交> …)`，根 main 当时的顶尖）← `1686ccd`（= 追加 14 链尾）。**两仓库均未推送。**

**追加 14（本轮）**——rcad `main`：**本交接文件所在提交** ← `3cf8d81d` ← `26031d47` ← `3422a867` ← `7dfa5884` ← `cc38f690` ← `3bf6858f` ← `2d01a7f9`（= 追加 13 链尾）
根 `main`：**对应的 rcad pointer sync 提交** ← `3d1c389` ← `6a08885` ← `38451b3`（= 追加 13 链尾）

| 提交 | 内容 |
|------|------|
| `3bf6858f` | **featlf init-profile 有效性三缺口 1:1**：① kernel `add_tedge` 建立**边顶点 tag 不变量**（顶点↔参数由既有 `vertex_params` 派生，= `BRepLib_MakeEdge.cxx L771-772` + `BRepPrim_Builder.cxx L143-155` + `TopExp.cxx L214-252`）；② `brep_check_result::direct_children` 不再改写边子顶点标签（原"镜像交换"只对 GWedge 来源成立）；③ 新 `reorder_edge_endpoints` 补齐 `BRepLib_MakeEdge::Init` 的 **reordonate**（cxx L636-651）与**周期 AdjustPeriodic**（L628-634） |
| `cc38f690` | **`brep_algo_is_valid` scratch 池改为重编号**：原先按原 index 塞池 ⇒ 跨池等值 index **别名**到无关 TShape（12 例 panic）。新 `renumbered_pool`：**后序** DFS + 每 TShape 唯一 Arc + 子引用同时改槽位与 `data` 句柄 |
| `7dfa5884` | **`BRepLib_MakeFace` 两形态补齐**：`(Pln,W,Inside=true)` 的 `CheckInside`（cxx L262-272 / L905-924）；`(W)` 的 `BRepLib_FindSurface` 曲面探测（cxx L189-262）。并抽出 `adopt_subgraph_into(pool, s)` 把 CutVehicle 池的 wire 收养进 profile 池 |
| `3422a867` | docs：追加 14 + 本交接重写 |
| `26031d47` | docs+libs：把 `FalseSide` 根因（**section 0 边**）钉进追加 14/交接队列；探针清理 |
| （本 session 尾批，见 §0.2） | **`Geom_Curve::ReversedParameter` 链 1:1**：`Curve3::reversed_parameter` 补齐 10 变体分发（此前落 trait 默认**恒等**）+ `geom_curve_reversed` 的 `Trimmed` 臂按 `Geom_TrimmedCurve::Reverse`（cxx L78-84）套 `ReversedParameter` + `BRepFeat_MakeRevolutionForm::init` 的 `cc->Reverse()`（cxx L741-744）按字面翻译 |

**追加 13（上一轮）**——rcad 7 提交 / 根 4 提交：

rcad `main`：`a4a4b0de` ← `8b5a7b26` ← `6b6c7089` ← `04e1b5b3` ← `4dfe09ac` ← `d6730766` ← `f4fa8168` ← `2cdee66b`（上轮交接）
根 `main`：`4f5306f` ← `62356d8` ← `933deda` ← `360eef7` ← `f4dda9e`

| 提交 | 内容（追加 13 = 上一轮） |
|------|------|
| `f4fa8168` | **feat Transform 链 1:1**：kernel `Trsf::set_rotation/set_translation/is_negative` + `mat_set_rotation`（gp_Mat.cxx L122-159 Rodrigues）+ `Plane::transform`；新文件 `topalgo/brep_builderapi_transform.rs`（Perform 分支判定 cxx L48-49 + TrsfModification 附加步骤：容差 ×\|scale\|、RevFace；两分支物化几何以适配 rcad 平读架构）；featlf **居中肋 GAP 接通**（cxx L186-194）；modeling 载体删除、gtests 改指新家 |
| `4dfe09ac` | **面级 BOP 三修复**：① `bnd_lib::surface_bounding_box` 增 `Trimmed(Plane)` 分支（= `GeomBndLib_Plane::Box` hxx L52-74）——无 wire 的 `MakeFace(Pln,u,v)` 面此前 bbox VOID ⇒ 不进 FF 的 BB 树 ⇒ **FF 候选 0 对**；② `CutVehicle` root 挑选加面级回退（BuildShape L1092 compound 平铺语义）；③ `geom_proj_lib::try_project_direct` 增 (BSpline, Plane) 臂 |
| `04e1b5b3` | **`IntTools_FaceFace::perform` 解包 RectangularTrimmed**（= `GeomAdaptor_Surface::Load` cxx L423-431）——plane×plane FF 截面从采样 BSpline 变**精确 Geom_Line**，顺带跳过 Geom2dInt 墙 |
| `6b6c7089` | **并行批**：BRepCheck_Analyzer 全家桶（8 文件 6,355 行，Analyzer+Result+六子类，38 个 Status，3 GAP 中性）+ Extrema_ExtCC2d 链（5 文件 ~2,700 行，全 25 臂）+ `brep_algo_is_valid` 真身（保索引 scratch 池） |
| `8b5a7b26` | **`Geom_BSplineSurface::UIso/VIso` 真身**（+ `BSplSLib::Iso` L1617-1740 与 `BSplCLib::Eval` **in-place 角切割** L865-870 重载）→ tktopalgo 36/36 |
| `a4a4b0de` | **`BSplineSurface` 补 OCCT 四标志**：存储 `is_periodic_u/v`（ctor 入参语义）+ 访问器 `is_rational_u/v`（= static `Rational` L110-138 的相邻权重 ULP 判定）；38 构造点逐点补齐；nsections 的 `PeriodicFlags` 载体删除 |

## 4. 下一轮队列（三域，按优先级）

### 4.0 追加 15 后的队列重排（2026-09-12，**从这一节往下取**）

**已完成（追加 15）**：原第 0 项 = `featrf_a1` 不终止 ⇒ **清零**（见 §0.3）；顺带清掉同链三道墙。
**新增第 0 项（最急，a1 的直接下一墙）**：

0. **`BRepFeat_RibSlot::LFPerform` 的结果装配**（`featrf_a1`）：init/`perform` 骨架已全通，停在
   `surface area: expected 109.511, got 0` ⇒ 结果为空。入口 = `feat/brep_feat_rib_slot.rs::lf_perform()`
   → `my_gs_hape / my_map / my_glued_f` 的装配路径，与 OCCT `BRepFeat_RibSlot::LFPerform` 逐行对照。
   **顺带**：`featrf` 的 a4/a5/a7/a9 停在 init 的 `is_done`（各自的 `NoFaceProf`/`NoSlidingProfile` 面），可分头取证。
   **注意**：该网格的参考拓扑断言当前被静默跳过（§0.3 卫生问题 1）——修生成器命名前，必须先手写/对齐
   `step_reference` 里的 **`occt_boolean_feat_featrf_a1.json`** 才能拿到 V/E/F/S 判据。

**原第 1–4 项原样保留**（见 §4.1–§4.4 正文），其中第 1 项 = `BOPAlgo_Section` 对共面面片对的输出（`FalseSide ×5` 根因）。

### 4.1 tkfeat —— 第 1 优先级：featlf 下游墙（**追加 14 后重写**）

**现状（追加 14 实测）**：featlf 真实断言 15 例中 **`is_done` 门通过 11**（追加 13 时 0），
**kernel panic 清零**（追加 13 时 12 例 panic 于 `topods.rs` 的池别名）。
profile 有效性问题（`InvalidPointOnCurve` / `NotClosed` / `UnorientableShape` / `NoSurface`）
**已全部关闭**——见 §0.1 六条修复。剩余失败**全部在下游**：

| 墙 | 例数 | 入口 |
|----|------|------|
| `FalseSide` | **5**（b3/c5/d7/d8/d9） | **根因已定界到"section 0 边"**（本轮探针实测）：`Propagate` 的 `BRepAlgoAPI_Section(fac, CurrentFace, approximation=true)` 在 **b3 返回 `nedges=0`**（Compound 在、无边），而 a3 的同类 section **有 1 条边且判定完全正确**。⇒ 墙在 **`BOPAlgo_Section` 对共面面片对的输出**，不在 `Propagate`。入口 = `bop/brep_algo_api/mod.rs::SectionOp` → `run_build_section_brep` → `BOPAlgo_Section.cxx`（重点：FF 重叠面片对） |
| `NoExtFace` | **3**（b5/b6/b7） | `BRepFeat_RibSlot::ExtremeFaces`（cxx **L747-1319**）：`ChoiceOfFaces` / `LocOpe_CSIntersector` 路径 |
| `NoFaceProf` | **3**（b1/b2/c6） | `brep_feat_rib_slot_b.rs` 的 6 个 `profile_ok=false` 返回点之一（临时探针逐个区分） |
| `BRepTools_Modifier::Perform` GAP | **1**（a3） | `feat/loc_ope_prism.rs:93`（消费 = `LocOpeLinearForm::int_perf` / `perform_trans`） |
| `draft_modification_1_b.rs:1269` GAP | **1**（b4） | tkoffset Draft 前沿 |
| Draft（depouille）断言 | **2**（e4/e5） | off-chain，属 tkoffset Draft |

**上一轮（追加 13）的根因记录已被追加 14 修正**：`Edge InvalidPointOnCurve ×3` +
`Face NotClosed` + `Face UnorientableShape` 的**共同上游是 rcad 边顶点 tag 约定不统一**，
**不是**布尔输出几何质量——布尔分割面本身健康。

**该靶子的下游关联**：**featprism(0/6) / featrevol(0/45) / featrf(0/5)** 共享同一条
RibSlot/Form 链，profile 修好后它们的失败点也随之下移。

### 4.2 通用依赖：Geom2dInt_GInter 通用 2D 求交批（一石二鸟）

- **消费点 1（feat）**：`feat/brep_feat_rib_slot.rs:466` 的 `Geom2dAPIInterCurveCurve` 非 line/line 臂仍 panic
  （`IntCurve_IntConicConic` / `Geom2dInt_GInter` 未翻）。**当前 12 例 featlf 已不经过它**（Trimmed 解包后
  截面为 Line ⇒ 走解析臂），但其它几何（BSpline 截面）会回到这条墙。
- **消费点 2（analyzer）**：`brep_check_wire.rs` 的 `SelfIntersect` 与 `brep_check_face.rs` 的静态 `Intersect`
  的 GAP 臂（模块头已列）——译好即同时激活。
- 参考真身：`fillet/chfi3d_builder_c2c.rs` 的 cncrn `Geom2dIntGInter` re-host（E3-G 批次已交付，可作模板）。

### 4.3 tkfillet —— blend 前线

**现状**：blend_simple 0/11 · blend_complex 0/2（a1/a3/q1 早期待过、现亦败，与 E3-U 记档一致）。
`fillet2d_*` 两组全绿。**新立档项（E3-Q 起，全部仍开）**：

1. **a1 上游审计（最高优先）**：`ChFiKPart_ComputeData::Compute`（cxx **L94-106**）→
   `ChFiKPart_MakeFillet` / `FilPlnPln` 的 SD 构造与**支撑面 pcurve 建立**。
   E3-Q 实测：SpKP hatcher 已真身化、域已从退化 [10,10] 变为 [0,100] 两端有点，
   但 a1 可观测输出**逐位未变**（面积 30671.24、7F/16E/12V+OPEN_SHELL vs 参考 59527.9、7F/15E/10V+CLOSED_SHELL）
   ⇒ 缺陷在 `SplitKPart` 上游。
2. **★ 优先查 `SimulSurf`/`PerformSurf` 的 44 处 pending-leaf stand-in**
   （`chfi3d_builder_2.rs:840/869`、`chfi3d_builder_2b.rs:76`，**0 调用者**）—— fillet 双曲面求交路径
   不可达，**极可能是 a1 上游缺陷根源**（E3-Q 审计记录原话）。
3. **OCCT 8 次 `SetArc` 的构成**：OCCT 真值 = a1 调 `ChFiDS_CommonPoint::SetArc` 8 次
   （`FORWARD@param10.0`×4 + `REVERSED@param90.0`×4），rcad 仅 6 次且全在 param 10.0
   ⇒ rcad 从不生成第二端（param 90）的弧值。扩展尾依赖 `SearchFace` / `ChFi3d_EdgeState`
   （rcad 目前把 SpKP.cxx L1151-1203/L1287 的扩展尾降级为 flag 重置）。
4. **blend 断言余项**：q4（SD 全链健康但最终结果 = 输入未动）、q2（三边定半径多棱非 EvolRad 路径）、
   p9（`filbuilder_c2.rs:977` pivot curve）、q7（blend-on-blend 角球，面积 150 vs 133.982）、
   a2/a4/p8/x1（面积差）、g9（2199.11 vs 2104.35，4.5%）。
5. **E3-Q 方向铁律（血证，必须遵守）**：对齐**单向**——OCCT 行为即规格，
   **禁止任何反向表述**（"OCCT 应走…"），分析语言必须始终是"OCCT 在第 N 行做 A，rcad 在第 M 行做成 B，故改 rcad"。

### 4.4 tkoffset —— 三件专批 + 矩阵余项

1. **`GeomAPI_ProjectPointOnCurve` 真身**（解 `trim_edge` 的 `FindParameter` 回退；
   `make_offset_d.rs:609` 为锚点）。
2. **`ExtentEdge` 面盒延伸语义**（`brep_offset_inter2d.rs`）。
3. **`GeomInt_IntSS` 专批（~1,860 行）**：`PrepareSurfaces` / `DefineUVMaxStep` /
   `IntSS_1::MakeCurve`（L275-1096）/ `BuildPCurves` / `TreatRLine` / `TrimILine` + `WLApprox`。
   ⚠ E3-U 翻案：**所谓"真身" `loc_ope_split_drafts_b.rs:595` 本身是 deferred stub**，需从零翻译。
4. **offset IsDone 40 例矩阵**：`perform_planes` 桩**已修**（现委托 bare-face
   `bop::int_tools::face_face::perform_face_face_planes`，见 `brep_offset_tool_d.rs:209`），
   按 GAP 迭代清障模式继续推进（下一候选 = `BRepLib::SortFaces` @ `brep_offset_make_offset.rs:297`
   及其后的 `MakeLoops`/`MakerVolume` 链）。
5. `offset_shape_type_i_c` **19 例 = 生成器坏输入工件**（owner = `occt-test-gen`，**只分类不修**）。
6. `draft_angle` 0/49（TKOffset Draft）——本轮未涉，属独立前沿。

### 4.5 其它（低优先，随手）

`rst_int.rs` 的 `FindParameter` Restriction 分支（`IntPatch_HInterTool::Project`）·
`u/v_iso` 其余面型臂（LinearExtrusion / Ruled / Offset / Pipe / Coons / Trimmed 仍 GAP）·
StepWriter 周期面 seam 保真度 · prism 构造迁移至 `primapi/` ·
`tools/occt-test-gen/tests/` 遗留跟踪文件 · `rcad-kernel/src/math/gprop/` 已删（本轮）。

## 5. 操作配方（复用，勿重摸）

```bash
# --- 1) 门槛（每次改动后）---
cd /c/Users/lilu/works/rcad-pro/rcad
cargo test -p rcad-algo --lib                 # 415/0/0
cargo test -p rcad-kernel --lib               # 689/0
cargo test -p rcad-algo --test tktopalgo_gtests        # 36/36
cargo test -p rcad-algo --test pavefiller_stage_tests  # 26/26
cargo test -p rcad-algo --test builder_stage_tests     # 76/76
cargo test -p rcad-algo --test builder_stage_smoke     # 1/1

# --- 2) 八网格（exe 必须取根仓库 target/debug/deps 最新产物）---
cd /c/Users/lilu/works/rcad-pro && bash output/run_eight_grids.sh

# --- 3) OCCT 侧插桩对拍（最快的定界手段；TKBO/TKPrim/TKFillet/TKOffset 同理）---
# 1) 改 C:/Users/lilu/works/OCCT/src/...（必须用 Edit 工具，Bash 写 OCCT 源被钩子拦）
# 2) cmd.exe //c "C:\Users\lilu\works\rcad-pro\output\build_tkprim.bat"   # 换 target 即可
# 3) cp -f C:/tools/occt-debug/win64/vc14/bind/<DLL>.dll tools/occt-bool-runner/build/Debug/
#    grep -ac <探针串> tools/occt-bool-runner/build/Debug/<DLL>.dll    # 必须 >0（DLL 遮蔽陷阱）
# 4) cd tools/occt-bool-runner/build/Debug && RCAD_WS_PROBE=1 ./occt_bool_runner.exe <grid> <CASE>
# 5) 用后 git checkout -- 还原 OCCT 源（探针即用即清）

# --- 4) 单用例复跑（从根仓库，target 名用完整 merged 名）---
cd /c/Users/lilu/works/rcad-pro
cargo test -p occt-generated-tests --test generated_occt_boolean_feat_featlf "feat_featlf_a3::" -- --nocapture

# --- 5) 域网格真实断言计数（§2 配方见正文循环体）---
cargo test -q -p occt-generated-tests --test generated_occt_boolean_<grid>

# --- 6) 生成器重生成（**必须在仓库根目录跑**，按 CWD 解析输出路径）---
OCCT_SRC="C:/Users/lilu/works/OCCT" cargo run -q -p occt-test-gen -- --batch-boolean --merge-groups --batch-grid <grid>
```

## 6. 本轮固化/强化的坑（勿重复踩）

0. **（追加 15 新坑，最贵）① `Shape::is_null()` 不能当作"槽位为空"的判据**：null 形状按 `index == usize::MAX` 判定，
   而**池外 builder 形状**（`BRepSweepBRepBuilder::{make_face,make_wire,make_vertex}`）的 `index` **也是 `usize::MAX`**
   ⇒ 合法形状被判空。第一版用 `is_null()` 判 `outer_wire` ⇒ `brep_sweep::prism` 单测立刻 6 面→2 面（侧面全丢）；
   改判**子形状类型**（`TShape::Wire`）后恢复。**判"槽位有没有东西"一律看类型，不看 index。**
   ② **同一个 `gp_Trsf` 的多次 `Perform` 在 rcad 的就地物化下必须去重**：OCCT 把运动挂在**结果 location** 上、TShape 不变
   ⇒ 同一 TShape 出现在两个结果里是"各转一次"；rcad 就地改 TShape ⇒ 第二次必须跳过（`perform_shape_once` 的 `done` 集），
   否则**双重变换**（`BRepFeat_MakeRevolutionForm::Perform` 的 `mySkface == myPbase` 就是活例）。
   ③ **"对拍数值"之前先做"语句级对拍"**：`featrf_a1` 的不终止**不是**浮点/容差问题，而是一处 OCCT 没有的**重算语句**
   （追加 14 的探针把范围收窄到"cc 重建 / theLastPnt 推进 / myTol 判据的数值"，方向仍错了一格）。
   **`while(!FirstOK)` 这类带回卷（`it.Initialize`）的链式循环对任何多余语句都极度敏感**——循环体内凡见 OCCT 没有的赋值，先删再过。
   ④ **face 的子形状枚举必须与 `BRep_Builder::Add` 的语义两侧对齐**：OCCT 把 wire 与 internal vertex 放在**同一个** `myShapes`
   列表里（TopoDS_Iterator 按此枚举），rcad 用三个类型化槽表示；**生产者**（`brep_sweep_builder::add`）与**消费者**
   （`brep_algo/tool.rs::sub_shapes`）必须同时对齐，否则"顶点占 wire 槽 / 空占位被子形状枚举吃掉"这类隐患只在扫掠路径才炸。

1. **"零候选/空结果"先查 bbox**：无 wire 的 `BRepLib_MakeFace(Pln,u,v)` 面的**唯一** bbox 来源是曲面 UV 窗分支
   （`GeomBndLib_Plane::Box`）；缺分支 ⇒ VOID 盒 ⇒ **静默退出 FF 的 BB 树**（症状是"候选对为 0"而非显式报错）。
   本轮 featlf 的第一个决定性根因。
2. **scratch 池必须**重编号**（追加 14 修正本条）**：`Shape::index` **只在创建它的 pool 内有效**；
   把多个 pool 拼装出来的子图按原 index 塞进 scratch 池 ⇒ **等值 index 别名到无关 TShape**
   （症状：`Shape 14 is not a Face` / `Shape 3 is not an Edge`）。正确做法见
   `feat/brep_feat_form_2.rs::renumbered_pool`：**后序 DFS**（子图是 DAG，pre-order 不是拓扑序，会在
   `children are built before their parent` 上 panic）+ 每 TShape 恰好一个 Arc + 子引用**同时**改
   **槽位与 `data` 句柄**（只改槽位无效——`child_occurrences` 从 `s.data` 读子件，缺口会"下沉一层"）。
   `Shape::null()`（`index == usize::MAX`）必须在收集阶段就跳过。
2b. **边顶点 tag 是 rcad 的隐性数据模型约定（追加 14）**：`TEdgeData.first/last` 的**语义由 `orientation`
   标签决定，不由槽位顺序决定**——OCCT 侧 `BRepLib_MakeEdge::Init`（低→高，F/R）与
   `BRepPrim_Builder::AddEdgeVertex`（高→低，R/F）**槽位顺序相反、标签一致**，
   而 `TopExp::Vertices(E,V1,V2,CumOri)` 与 analyzer、WireExplorer **只认标签**。
   任何"first 就是低参数顶点"的位置假设都会在另一类来源的边上翻车。**收口点 = `add_tedge`**
   （用 `vertex_params` 派生标签），新写边构造器时不要绕过它。
3. **Trimmed 包装面在按 surface-type 分派处必须解包**（`GeomAdaptor_Surface::Load` cxx L423-431 语义）：
   否则 plane×plane 会落进"other"臂产出**采样 BSpline 截面**（t=[26,68] 弧长参数化特征）而非精确 Geom_Line。
   与 E3-U 的 Trimmed 曲线未解包同族，**第二次踩**。
4. **`BSplCLib::Eval` 有两个重载**：`BSplSLib::Iso` 用的是 **L865-870 的 in-place 角切割**版，
   不是 L3640 的基函数版（rcad 既有 `eval_flat` 是后者）——误用会在 2·Degree 局部窗上**索引越界**。
5. **立卡前先 grep OCCT 源码确认类名**：`Extrema_ExtCF` **在 OCCT 不存在**（`BRepExtrema_ExtCF` 的内核引擎是
   `Extrema_ExtCS`，rcad 早已翻译）——队列里的"ExtCF 真身"是**伪缺口**。同理"translator 说有 GAP"先查真身是否已在库
   （E3-W 追加 5 的 `BRepTools_Quilt` 教训）。
6. **域网格的通过数口径**：`*_geometry_loads` 是**恒过的占位**，真实断言只有
   `*_draw_script_rcad_equivalent`——报数时必须区分，否则会把 0/11 误读成 11/22。
7. **文档基线会漂移**：tktopalgo/tkg3d 的失败在交接树上就存在——**开场必须实跑**，不能照抄文档的"全绿"断言。
8. **stash 归因法**：`git stash -u → 跑测试 → git stash pop` 是区分"既有失败 vs 本轮回归"的最快手段（~2 min）。
   ⚠ 根仓库里**存在其他会话遗留的 stash**（如 `830e33e`/`9e4a818`）——stash 是 LIFO，操作前后核对 `git stash list`，
   勿动他人 stash。
9. **探针即用即清**：提交前必查 `git diff | grep -c "+.*eprintln"`（本轮全程 = 0）。探针必须 `--nocapture` 才可见。
10. **源码一律走 Edit/Write 工具**：`Bash` 直写 `libs/**/*.rs` 与 OCCT 源会被 Mimosa 钩子拦（`sed -i` 亦拦）。
11. **`cargo check --workspace` 会因 `rcad-constraints` / `rcad-render` 的既有破损报错**——按 crate 检查
    （`-p rcad-kernel/rcad-algo/rcad-step/rcad-modeling`）。
12. **`Shape::null()` 的 `index == usize::MAX`**：把它当池索引会 panic；重编号时在收集阶段按 `is_null()` 跳过。
13. **`BRepLib_MakeFace` 有**三个**形态，别串**（追加 14）：`(Pln, W, Inside)` 需要 `CheckInside` 尾巴
    （L267-270）；`(W)` 需要 `FindSurface` 曲面探测（L189-262）；`(S, U1,U2,V1,V2)` 是矩形窗。rcad 此前把
    `(W)` 翻成了"无曲面面"⇒ `BRepCheck_Face::Minimum` 的 `NoSurface`。
14. **`add_tedge` 的 `vertex_params` 是位置匹配**（`d0 <= d1`）：闭合/退化曲线上两个顶点会都落到 `range[0]`，
    这时标签派生退化为"两面都 FORWARD"——与 OCCT 的 closed-edge 分支（V1 与 V2 同一点）语义一致，不是缺陷。
15. **`reversed_parameter` 曾是**恒等****（追加 14 补记）：`impl CurveEval for Curve3` 漏了分发 ⇒ 所有
    跨类型调用静默退化。写"反向曲线/反向修剪"类代码时，**不要**用手写的位置互换（`(basis_rev, last, first)`)
    代替 OCCT 的 `Reversed()` / `ReversedParameter()`——两者的差别在 `Line`/`Circle`/`Ellipse` 上**一定是错的**
    （`f64::clamp` panic 只是其最响的症状；静默错几何更难查）。例：`brep_fill_section_law.rs:267`、
    `chfi3d_perform_elspine.rs:1069` 两处消费者因此一直在拿恒等值——**修好后要复查它们的下游**。
16. **域网格"失败"有三层**：panic < assert < **不终止**。最后一层最贵（烧掉 timeout 才暴露，且
    `RCAD_TEST_TIMEOUT_SECS` 只在 `run_grid.ps1` 侧生效，`cargo test` 直跑会一直挂着）。
    跑非门槛网格时**一律加 `timeout`**（如 `timeout 120 cargo test …`）。当前已知不终止 1 例：`featrf_a1`。

## 7. 资产位置（本轮新增/更新）

| 资产 | 路径 |
|------|------|
| 边顶点 tag 不变量（一处收口） | `libs/rcad-kernel/src/topo/topods.rs::add_tedge`（追加 14） |
| `Curve3::reversed_parameter` 10 变体分发 | `libs/rcad-kernel/src/geom/eval.rs`（追加 14 补记） |
| `Geom_TrimmedCurve::Reverse` 的 `Trimmed` 臂 | `libs/rcad-algo/src/feat/brep_feat_rib_slot.rs::geom_curve_reversed`（追加 14 补记） |
| `cc->Reverse()` 字面翻译（RevolutionForm） | `libs/rcad-algo/src/feat/brep_feat_make_revolution_form.rs`（init 的 `while !first_ok` 头部，cxx L741-744） |
| 边子顶点不再改写标签 | `libs/rcad-algo/src/topalgo/brep_check/brep_check_result.rs::direct_children`（追加 14） |
| `BRepLib_MakeEdge::Init` reordonate/周期 | `libs/rcad-algo/src/feat/brep_feat_rib_slot.rs::reorder_edge_endpoints`（追加 14） |
| `BRepLib_MakeFace` 的 `CheckInside` | `libs/rcad-algo/src/feat/brep_feat_rib_slot_b.rs::check_inside`（追加 14；`..make_revolution_form.rs` 同用） |
| `BRepLib_MakeFace(W)` 曲面探测 | `libs/rcad-algo/src/feat/brep_feat_rib_slot_b.rs::make_face_wire`（追加 14） |
| 跨池子图重编号 + 收养 | `libs/rcad-algo/src/feat/brep_feat_form_2.rs::{renumbered_pool, adopt_subgraph_into}`（追加 14） |
| BRepBuilderAPI_Transform + TrsfModification 附加步骤 | `libs/rcad-algo/src/topalgo/brep_builderapi_transform.rs`（新） |
| BRepCheck_Analyzer 全家桶（8 文件 6,355 行） | `libs/rcad-algo/src/topalgo/brep_check/{brep_check_analyzer,brep_check_result,brep_check_vertex,brep_check_edge,brep_check_wire,brep_check_face,brep_check_shell,brep_check_solid}.rs`（新） |
| Extrema_ExtCC2d 链（5 文件 ~2,700 行） | `libs/rcad-kernel/src/base/{extrema_curve2d_tool,extrema_ext_p_elc2d,extrema_ext_elc2d,extrema_gen_ext_cc2d,extrema_ext_cc2d}.rs`（新） |
| `BRepAlgo::IsValid` 真身 | `libs/rcad-algo/src/feat/brep_feat_form_2.rs::brep_algo_is_valid` |
| `BSplSLib::Iso`（flat-knots + in-place 角切割） | `libs/rcad-kernel/src/math/bspl_lib.rs::bspl_slib_iso` |
| `Geom_BSplineSurface::UIso/VIso` 臂 | `libs/rcad-algo/src/geomalgo/approx_curve_on_surface.rs::surface_u_iso/surface_v_iso` |
| `BSplineSurface` 四标志 + Rational 访问器 | `libs/rcad-kernel/src/geom/mod.rs`（`is_periodic_u/v` 字段 + `is_rational_u/v()` + `standard_epsilon`/`next_after`） |
| Trimmed 解包 | `libs/rcad-algo/src/bop/int_tools/face_face.rs::unwrap_trimmed_surface` |
| `Trimmed(Plane)` 盒 | `libs/rcad-kernel/src/base/bnd_lib/mod.rs::surface_bounding_box` |
| featlf 居中肋接线 | `libs/rcad-algo/src/feat/brep_feat_make_linear_form.rs`（L186-194 处） |
| 既有可复用探针 | `RCAD_BS_DEBUG`（builder/pave_filler 的 build_rc/结果装配）、`RCAD_FF_DEBUG`（FF 配对+曲线类型）、`RCAD_MB_DEBUG`（MakeBlocks 移位）、`RCAD_WS_DEBUG`（WireSplitter）、`RCAD_SPLIT_DEBUG`（CUT 结果） |
| 权威长档 | `rcad/docs/tkfeat-fillet-offset-port-plan.md` §E3-W 追加 11–14 |
| 模块归属表 | `rcad/docs/module-map.md` §2/§3（`feat/`↔TKFeat、`fillet/`↔TKFillet、`offset/`↔TKOffset） |
| 重复实现审计 | `tools/occt-impl-audit/audit.py`（`uv run python tools/occt-impl-audit/audit.py .`）→ `docs/occt-impl-audit.md`（域风险排序：fillet 31 · feat 25 · offset 7） |
