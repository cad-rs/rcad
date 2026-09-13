# E3-W 交接：三域（TKFeat / TKFillet / TKOffset）翻译推进 —— 2026-09-12

> **一句话（追加 19 收尾态）**：**工作口径仍是"先译后调"**——先把缺失的 body 逐一 1:1 译全、把重复/过期载体收敛到真身，等代码基本译完再回头调试与修测试（用户明确指令）。
> 本轮落地 **4 批**：① `ProjLib_HCompProjectedCurve` 在 `GeomPlate_BuildPlateSurface` 的**三处接线**（真身本就在库）② **池外 `BRep_Tool::Curve`**（`curve_pool_free` + 守卫分流）
> ③ **`GeomLib::BuildCurve3d` 家族 1:1**（`AdvApprox_PrefAndRec` / `with_cut_tool` + 子空间存储 / `GeomLib_CurveOnSurfaceEvaluator` / `isIsoLine`+`buildC3dOnIsoLine` / `GeomLib_MakeCurvefromApprox`）+ 四处 GAP 载体收敛
> ④ **池外读取续链**（`edge_data_pool_free` + `build_curves3d.rs::edge_data`，该链上十处池索引读收敛）。
> **off-gate 效果**：`blend_simple` a2/p8/p9 离开 GAP panic；**`offset_shape_type_i` 的池外 panic 清零**（a3/a4/d2/d3 全部落到测试断言）；`offset_shape_type_a` a4 推进到 `brep_algo/image.rs:159`；
> `draft_angle` 库内 panic **25 → 20**。六门槛与八网格**全程零回归**。
> **★ 本轮有三条"旧笔记/旧注释会过期"的硬教训**（坑 28/29/30）：**同一 helper 行号会掩盖层推进（必须看 backtrace）**、**旧笔记的"架构难点"要回查 OCCT 基类**、**`GetType()` 常量返回决定分支**。
> （前三轮：**追加 18** = 翻译补全轮（五处真身 + 三处重复删除 + `int` 数据模型对齐）；**追加 17** = D3 结案 + `featrf_a1` init 首次通过；**追加 16** = 0a 收尾 / a1 拓扑全等 / TKOffset 定界。）

## 0. 新 session 一句话提示词（直接粘贴 —— 追加 19 收尾态，2026-09-13）

> 读 `rcad/docs/e3w-handover-tkfeat-fillet-offset.md`（本交接：门槛实测值 / 提交链 / 三域队列 / 配方 / 坑清单）
> 与 `rcad/docs/tkfeat-fillet-offset-port-plan.md` §E3-W 追加 11–19（权威脉络，**追加 19 是当前状态**）；
> 先 `cd rcad` 跑 6 条门槛确认 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**，
> 再 `cd /c/Users/lilu/works/rcad-pro && cargo test --no-run -p occt-generated-tests`（**重编 exe，否则八网格会拿旧产物误判**）
> 后 `bash output/run_eight_grids.sh` 确认**八网格 8/8**；
> 然后**按 §4.6 的三域队列继续"先译后调"**——**先只做 1:1 翻译/接线**（逐行语句对照 + OCCT 行号锚点 +
> 禁载体/禁等价替换/禁凑结果 + 架构差异先消灭再对齐），**代码基本译完前不针对性修测试数值**；
> **涉及布尔层的代码直接复用已对齐的 `bop/**`（TKBO）实现、不要重写第二份**（`BRepAlgoAPI_*` 包装层按 OCCT 形式接线到既有真身）；
> **立卡前先按 OCCT 函数名 grep 全库**（`panic!("GAP…")`/`unimplemented!` 的文案会过期，本轮两例真身其实早已在库）；
> **失败层深度自检必须看 backtrace 而不是失败地图行号**（追加 19 实测：同一 `topods.rs:1799` 背后换了调用帧，行号完全看不出推进 —— 见坑 28）；
> 每 Edit 后 `cargo check -p rcad-algo`；**跑非门槛网格一律加 `timeout`**（§6 坑 16）；探针即用即清
> （提交前 `git diff | grep -c "+.*eprintln"` = 0，OCCT 侧探针用后 `git checkout --` 还原 + 重建 DLL）；
> 每批做完跑**六门槛 + 八网格 + 该域网格**，按坑 21 的**失败层深度**（不是通过数）自检，更新 port-plan §E3-W 追加，
> 并提交**两仓库**（rcad + 根仓库指针，rcad 推得上就推）。

### 0.0 当前状态速览（追加 19 收尾，2026-09-13；**rcad 顶尖 = 本交接文件所在提交**（用 `cd rcad && git log -1 --oneline` 即得；写就时基线为 `e841e805`）/ 根仓库指针 = 本文件所在提交，**均已推送**）

- **门槛与网格**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（375·378·379·373·12·102·83·110，**重编 exe 后**实测）。域网格**逐格通过数不变**：`fillet2d_fillet2d` 10/10 · `fillet2d_chamfer2d` 2/2 · `mkface_after_offset` 4/4 · `mkface_after_extsurf_and_offset` 32/32 · `feat_featlf` 0/15 · `feat_featprism` 0/6 · `feat_featrevol` **1/45**（a5）· `feat_featrf` 0/5 · `blend_simple` 0/11 · `blend_complex` 0/2 · `offset_shape_type_a/_i` 0/1 · 0/12 · `offset_faces_type_i` 0/8 · `draft_angle` 0/49 · `thrusection_specific` 0/26。
- **★ 工作模式（用户指令，追加 18 起生效）**：**先完成代码的等价实现（近乎 1:1 的翻译），代码基本译完再开始调试/修测试**。⇒ 队列里**先取"翻译/接线"项**，取"调试/定界"项前先确认没有未译的 body 挡在前面。
- **本轮已落地（4 批，见 §0.7）**：geomplate 的 ProjLib 三处接线（`d7beea3c`）· 池外 `BRep_Tool::Curve`（`830eb2c2`）· `GeomLib::BuildCurve3d` 家族 + 四处载体收敛（`e841e805`）· **池外读取续链收敛（`fe297616`：`offset_shape_type_i` 池外 panic 清零）**。
- **★ 域网格实测失败地图（下一轮的对照基线；口径 = 逐例 file:line）**：
  - `offset_shape_type_i`：a1/a2 → `brep_offset_inter2d.rs:1060`（`EdgeInter: E2 carries no pcurve`，OCCT 同处也 raise ⇒ **状态**）；**a3/a4/d2/d3 → 已离开库内，落到测试断言（310/422/534/646）** —— 批次 2 + 补记 1 的直接收益，**该网格池外 panic 已清零**；e1/e2/e3/e4/e6/e7 → 测试断言（758/870/948/1060/1171/1283）。
  - `offset_shape_type_a`：a4 → `brep_algo/image.rs:159`（本轮由 `brep_offset_make_offset_c.rs:52` 推到这里）。
  - `offset_faces_type_i`：a1/a2 → `brep_offset_inter2d.rs:1060`；a5/a6/g1/g2/g6/g7 → 测试断言。
  - `blend_simple`（11 例）：a1 → 测试断言 **L113**（面积，门 = 无界 pcurve）；**a2/p8/p9 → `proj_lib_h_comp_projected_curve.rs:449`**（本轮由 `geomplate/build_plate_surface.rs:1003` 的 GAP panic 推进到这里；`D0` 的 `Standard_DomainError` ⇒ **状态**类）；a3/a4 → `chfi3d_builder_2b.rs:583`；q1 → 测试断言 L963；q2 → `chfi3d_builder_c1.rs:1162`；q4 → `geom/bspline_ops.rs:329`；q7 → `brep_blend_walking.rs:143`；x1 → `base/convert/mod.rs:2008`。
  - `feat_featrf`：a1 → 测试断言 **L866**（面积 0）；a4/a5/a7/a9 → 测试断言（**未劣化**）。
  - `feat_featlf`（15 例）：全部为测试断言（本轮实测**库内 panic 已清零**，与追加 18 记的 3 例库内 panic 相比又前进一层）。
  - `feat_featrevol`：44 例全为测试断言 + a5 通过（1/45）。
  - `draft_angle`（49 例）：库内 panic **20 = 12 × `draft_modification_1_b.rs:1213` + 4 × `_1_c.rs:868` + 2 × `make_revol.rs:583` + 2 × 新层次点**（`brep_offset_api_draft_angle.rs:83`、`approx_int.rs:1071`）；**追加 18 为 25** ⇒ 本轮 5 例离开库内 panic。
- **三域下一步（详见 §4.6）**：**TKOffset** = ① 池外读取**续链**（`check_same_range` / `gcurve_range` / `brep_tool_curve_on_surface_index` / `brep_tool_range_on_surface` / `brep_tool_degenerated` / `brep_tool_tolerance` 收敛到受守卫的 edge-data 读取；**写回侧池外无池可变，需另法，勿硬凑**）→ ② `BRepTools_Quilt`/`FaceRestrictor` **产物入池**（根；入池后整条池读链一次解开，且拓扑计数随之正确；**牵涉拓扑计数 ⇒ 全套复测**）；**TKFillet** = 接力项 A 的"无界 pcurve 附着点"（a1 面积 −2e100 的门，本轮未触及）+ blend a2/p8/p9 的新墙（状态类）；**TKFeat** = `LocOpe_Generator::Perform` 的 `IsDone`（a1 粘合路径下一墙）。

### 0.7 追加 19 本轮落地（**翻译/接线轮**，2026-09-13；四批，rcad 提交链见 §3）

| 批次 | 提交 | 内容 | 实测（失败层深度） |
|------|------|------|--------------------|
| 1 | `d7beea3c` | **`ProjLib_HCompProjectedCurve` 在 `GeomPlate_BuildPlateSurface` 三处接线**（metrics 比较 cxx L1746-1802 / ProjectCurve L254-303 / ProjectedCurve L307-349）——真身与 `Adaptor3d_CurveOnSurface` 真身早已在库 | `blend_simple` **a2/p8/p9** 由 GAP panic 推进到 `proj_lib_h_comp_projected_curve.rs:449` |
| 2 | `830eb2c2` | **池外 `BRep_Tool::Curve`**：`topods.rs` 新增 `curve_pool_free`（OCCT **BRep_Tool.cxx L172-196**）+ `shape_is_in_pool`；`build_curves3d.rs::brep_tool_curve` 与 `topexp.rs::{brep_tool_curve_loc,brep_tool_range}` 按守卫分流 | `offset_shape_type_i` **a3/a4/d2/d3** 由 `brep_tool_curve` 内推进到 `BRepLib::check_same_range`（backtrace 取证） |
| 3 | `e841e805` | **`GeomLib::BuildCurve3d` 家族 1:1**（新增 `AdvApprox_PrefAndRec`、`ApproxAFunction::with_cut_tool` + 子空间存储、`kernel_curve2d` 桥、`GeomLib_CurveOnSurfaceEvaluator`、`build_curve3d`、`isIsoLine`/`buildC3dOnIsoLine`、`GeomLib_MakeCurvefromApprox`）+ **四处 GAP 载体收敛** | `offset_shape_type_a` **a4** 推进到 `brep_algo/image.rs:159`；`draft_angle` 库内 panic **25 → 20** |
| 4 | `fe297616` | **池外读取续链**：kernel `edge_data_pool_free` + `build_curves3d.rs::edge_data` 守卫访问器，该链上**十处**池索引读（`check_same_range`/`same_range`/`gcurve_range` + 四个 `brep_tool_*` re-host）全部收敛 | `offset_shape_type_i` **a3/a4/d2/d3 离开库内**（→ 测试断言 310/422/534/646）⇒ **该网格池外 panic 清零** |

**三条关键提醒（本轮实测）**：
1. **追加 18 记的"`Adaptor3d_CurveOnSurface` 以通用 `Adaptor3d_Curve` 为底"是不成立的**：`ProjLib_CompProjectedCurve` 的基类是 **`Adaptor2d_Curve2d`**（hxx **L37**），`AProj` 走的就是 kernel `CurveOnSurface` 既有的 2D 载体口径 —— **不要**因为旧笔记没落这条线。
2. **锚点勘误**：`BRep_Tool.cxx L410-452` 是 `BRep_Tool::CurveOnPlane`，**`Curve` 本体是 L172-196**。旧注释（本档与源码）引错过，批次 2 的注释已改为正确锚点。
3. **`ProjLib_CompProjectedCurve::GetType()` 恒返回 `GeomAbs_OtherCurve`**（cxx **L2108-2111**），所以 `Adaptor3d_CurveOnSurface::EvalKPart` 的 Plane 分支只可能落到 `myType = Other`、`EvalD1` 走通用链（`myCurve->D1` + `mySurface->D1` + 线性组合）——kernel `CurveOnSurface::d1` 正是这条，**无需为这个消费点加 Line/Circle 分支**。

### 0.5 追加 17 本轮落地（1:1 修复三处；rcad 提交链见 §3）

| 修复 | 位置 | OCCT 锚点 | 消除的症状 |
|------|------|-----------|-----------|
| **BndWire 走位改用真 WireExplorer 连通序**（D3 真因） | `topalgo/brep_top_adaptor/fclass2d.rs::wire_explorer_order`（新）+ `feat/brep_feat_rib_slot_b.rs::sliding_profile` 调用点 | `BRepFeat_RibSlot::SlidingProfile` **L1567-1569**（`BRepTools_WireExplorer explo(BndWire)`）+ `BRepTools_WireExplorer.cxx` L121-705（既有真身 `order_wire_edges`） | `featrf_a1` init 失败（profile wire 7 条自交 ⇒ `UnorientableShape` ⇒ `NoFaceProf`）；修后 5 条与 OCCT 逐项相同 |
| **`BoundSortBox::get_bounding_voxels` 的 C++ int 语义** | `rcad-kernel/src/math/bnd/bound_sort_box.rs` | `Bnd_BoundSortBox.cxx` **L577-590**（`std::clamp(static_cast<int>(v) - 1, 0, myResolution - 1)`：饱和 + int 回绕 + clamp 吃掉） | 无界段盒（`(-2,-1e100,5)-(-2,0,5)`，`IntCurvesFace_Intersector` 对无限直线采样）⇒ `subtract with overflow` panic |
| **解析访问器解包 `Trimmed`**（D2 的补完） | `bop/int_tools/hinter_adaptor.rs::{basis_curve_of(新), line, circle, ellipse, hyperbola, parabola, bezier, bspline}` | `GeomAdaptor_Curve::load` **L239-255**（`Load(BasisCurve, UFirst, ULast)` ⇒ `myCurveData` 全部来自基曲线） | `feat_featlf` **8 例**（a3/b3/b6/b7/c5/d7/d8/d9）死在 `hinter_adaptor.rs:218` 的 `Standard_NoSuchObject: Line`（`GetType()` 报 Line 而 `line()` 按原始变体匹配） |
| **TKFillet 端盖弧 range**（追加 17 补记 2） | `fillet/chfi3d_builder_0.rs::{elclib_parameter_circle, elclib_parameter_ellipse, reverse_curve}` | `ElCLib::CircleParameter` **L1199-1222**（`AngleWithRef` ⇒ `atan2(v·Y, v·X)`）+ `gp_Ax2::SetDirection` **hxx L548-571**（通用分支 ⇒ 圆反向保 X 翻 Y） | 端盖弧 270°（`[π,2.5π]`）⇒ 与 OCCT 一致为 90°（`[0.5π,π]`）；⚠ 落两条后 a1 面积 57328.76 → **−2e100**（暴露下一个已立档的"无界 pcurve"缺陷，见 §4.6 TKFillet） |
| **`geom_proj_lib::curve2d` 解包 `Trimmed`**（追加 17 补记 1） | `rcad-kernel/src/base/geom_proj_lib/mod.rs` 入口 | `GeomAdaptor_Curve::load` 同族（第三次：D2 → 解析访问器 → `curve2d`） | `featrf_a1` 的 `BRepFeat::IsInside` 在分类器前返回 false ⇒ `collage=false`；修后 `collage=true ope=Fuse`（与 OCCT 一致），粘合路径打通 |

**对拍通道（可复用，成本极低）**：OCCT 侧 `RCAD_WS_PROBE` 门控探针 → `output/build_tkfeat.bat` → 拷 `TKFeat.dll` 到 `tools/occt-bool-runner/build/Debug/`（`grep -ac WS-PROBE <dll>` 必须 > 0，防 DLL 遮蔽）→ `occt_bool_runner feat_featrf A1`。**featrf_a1 的 OCCT 真值已落档**（port-plan 追加 17）。用完 `git checkout --` 还原 OCCT 源并重建干净 DLL。

### 0.6 上一轮：追加 18 落地（**翻译补全轮**：先译后调，2026-09-12）

| 真身（OCCT 锚点） | rcad 落位 | 收敛掉的重复/载体 |
|------|-----------|-----------|
| `GeomLib::SameRange`（`GeomLib.cxx` **L842-970**） | kernel `geom::same_range_2d`（补 `Tolerance` 首参 + L902-907 `PConfusion` 守卫 + L924-969 periodic/非 periodic 分段 + 改调 `bspl_lib::reparametrize`） | `geomalgo/geom_lib_same_range.rs` 的 GAP → 委托；`shhealing/shape_fix/edge.rs` 的"返回输入不变"stand-in → 改名 `geom_lib_same_range` 并委托（ShapeFix_Edge.cxx L418/L447）；`bop/algo/pave_filler.rs` 的容差由硬编码 1e-7 改为 OCCT 的 `aTolPPC = PConfusion()`（1e-9，**此前是失真**） |
| `BRepCheck_Edge::Tolerance`（`BRepCheck_Edge.cxx` **L598-707**） | `topalgo/brep_check/brep_check_edge.rs::BRepCheckEdge::tolerance`（NCONTROL=23、L636-641 槽位搬迁、**L669 缝 pcurve 只用 `cr->Location()`**、逐坐标 `IsInfinite` 短路、`sqrt(max)*1.05`）；`HCurveAdaptor` 补 `value(U)` | `offset/brep_offset_make_offset.rs` 的 GAP → 委托（`update_tolerance` 串入 `&BRep`）；`brep_fill_sweep_b.rs` 的 `GeomLib_CheckCurveOnSurface` **近似载体删除**并委托 |
| `BRepCheck_Vertex::Tolerance`（`BRepCheck_Vertex.cxx` **L343-383**） | `topalgo/brep_check/brep_check_vertex.rs::BRepCheckVertex::tolerance`（收尾 **`sqrt(Tol*1.05)`**，与 Edge 版 **不同**）；已注明 rcad 只有 2 种点表示 | `offset/brep_offset_make_offset.rs` 的 GAP → 委托 |
| `ElCLib::To3d` 全集（`ElCLib.cxx` **L1339-1440**） | kernel `math/el.rs` 的 `elclib_to3d_*`（8 个：pnt/vec/ax22d/line/circle/ellipse/hyperbola/parabola；`Ax22d` 经 3 参 `gp_Ax2(P, VX×VY, VX)`） | — |
| `GeomLib::To3d`（`GeomLib.cxx` **L559-675**） | `geomalgo/geom_lib.rs::to_3d`（Trimmed 递归再裁剪 / Offset 重建 / Bezier+BSpline 极点提升 / 五类解析帧提升 / `Standard_NotImplemented`）；BSpline 周期由 `bspl::bspline_is_periodic` 导出 | `topalgo/brep_lib/build_curves3d.rs::geom_lib_to_3d` 的 GAP → 委托（⇒ **`BRepLib::build_curve3d` 的 cxx L362 平面分支可用**） |
| `GeomLib::ExtendSurfByLength`（`GeomLib.cxx` **L1485-1972**） | 真身本就在 `fillet/chfi3d_builder_c2_geomlib.rs::geom_lib_extend_surf_by_length` | `geomalgo/geom_lib_same_range.rs` 的 GAP → 委托（消费点 `brep_offset_tool_c.rs:617/683`、`brep_fill_sweep_c.rs:362/374`） |
| **数据模型对齐**：`Approx_ComputeLine` / `GeomInt_WLApprox` / `BRepApprox_Approx::SetParameters` / `IntTools_FaceFace::ApproxParameters` 的**度数全是 OCCT `int`** | `mydegremin/mydegremax`、`myDegMin/myDegMax`、`deg_min/deg_max`、两个 `nbp` 局部量 ⇒ **全部 `i32`**；以 usize 为索引的 rcad API 边界显式 `as usize` 并注明 | `panic: attempt to subtract with overflow`（`m_degmax = nbp - 5`）——OCCT 的 `nbp - 5` **允许为负**且**紧接着就被 clamp 回 `mydegremin`**（`Approx_ComputeLine.gxx` **L1336/L1344-1350**），usize 既下溢又丢符号语义 |

**关键提醒（两条，来自本轮实测）**：
1. `geomalgo/approx_curve_on_surface.rs` 里的 `is_iso_line`/`build_c3d_on_iso_line` 是 **`Approx_CurveOnSurface` 的静态副本**（OCCT 自身也有两份），`GeomLib::BuildCurve3d` 要用的是 **`GeomLib` 自己的**静态——翻译前先分清，否则会把两份混成一份。
2. **`panic!("GAP…")` / `unimplemented!` 的文案会过期**：本轮两例（`ExtendSurfByLength`、`ProjLib_HCompProjectedCurve`）的真身**早已在库**，只是没接线/没注册到消费者。**立卡前先按 OCCT 函数名 grep 全库**，再决定是"翻译"还是"接线"（与坑 7 同族）。

### 0.0.1 必读的三条硬约束（血的教训）

1. **`Shape::is_null()` 判空不可靠**（坑 0① / 坑 19）：**池外**构造的合法形状 `index == usize::MAX` 会被判成 null。**要给下游长期持有的形状一律入池**；确须池外时看**子形状类型**。
2. **DS 键数据的翻译契约**（坑 20）：主 DS **深拷贝**参数（`clone_arguments`），所以 `my_images`/历史/任何 DS 键数据都用**克隆** TShape；持**原始**形状的 API 级消费者必须经 `DS::argument_shapes` 翻译，否则**静默全空**。
3. **OCCT 的 `int` 字段不要译成 `usize`**（坑 26）：有符号中间值（如 `nbp - 5`）是**设计的一部分**；`usize` 会 panic 且丢语义。同理 C++ 的 `int` 算术在 Rust 里要 `wrapping_*`（坑 22）。搬运字段/表达式前先读 OCCT 声明。



> **★★★ 追加 16 收尾态（2026-09-12 尾，**最新权威入口**）**：三域**并行推进轮**——主代理做 0a 收尾、3 个子代理分域推进（详见 port-plan §E3-W 追加 16）。
> **TKFeat**：`featrf_a1` 的 bind 墙拆掉——Cut 走完全程（`done=true` status OK，结果为 4 面 compound）；两处架构修复 = **`BRepAlgoAPI::Shape()` 必须是 `aResult` 的 Compound**（OCCT `BOPAlgo_BOP` L1042-1050/L1106 + `BRepAlgoAPI_BuilderAlgo::Build` L157）、**DS 深拷贝参数的翻译契约**（`DS::argument_shapes` + `CutVehicle::ds_keyed`）。提交 `fa0ddc43`。
> **TKFillet**：**blend_simple_a1 拓扑首次与 OCCT 全等**（V/E/F/S/SHELL/SOLID 完全一致，`step-topo-diff` fully matches），面积 30671.24 → **57328.76**（OCCT 59527.9）。两处 1:1 修复 = `gp_Ax3::YReverse()` 的左手系（FilPlnPln）+ pcurve range 由 3D 表示覆盖。提交 `0263c85d`。**剩余缺口已精确到 fillet 端盖弧 range `[π,2.5π]` vs 参考 90°**（残差 700π 逐面核对吻合），且喂入缺陷 `chfi3d_builder_0.rs::reverse_curve` 对 Circle 是 no-op（须与弧 range 同批落）。
> **TKOffset**：本轮无代码落地，换来**逐簇证据地图 + 6 条域外任务**（port-plan 追加 16 ③）。
> 六门槛与八网格**全程零回归**（收尾时三批同时在树复跑：415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1；八网格 375·378·379·373·12·102·83·110）。**新坑 19/20/21** 见 port-plan 追加 16 尾（`is_null()` 池外误判第二次实测 / DS clone_arguments 契约 / 并行时的 A/B 污染）。

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
2. **`[ENTRY]` 遗留探针已清（追加 15 补记 3）**：`brep_sweep/num_linear_regular_sweep.rs:284` 那行无门控的 `eprintln!("[ENTRY] …")`（来自 `ddc8d551`）已删除。
   ⚠ 复核结论（追加 15 补记 3）：`bop/**` 里的 `[EF-DBG]`/`[EF-EDGE]`/`[EF-CB]`（56 处）**是 env 门控的可复用探针**（实测跑 boolean 用例命中 0 次），
   **不是垃圾、勿删**——"源码里有 eprintln"必须先跑一次确认是否真的打印再动手。

## 0.4 ★★ 追加 15 尾的两个交付批次（**已验收并提交**；新 session 只需做收尾的**拆分**活）

两个子工作流（域不相交、各自 `CARGO_TARGET_DIR`）**均已交付、已由主代理验收并提交**；六门槛与八网格在**两批同时在树**的树上实测通过。

| 批次 | 提交 | 域 | 交付内容 | 验收实测 |
|------|------|----|---------|---------|
| **A. `GeomInt_IntSS` + `IntSS_1` + `GeomInt_LineConstructor`** | `3f6f746a` | `geomalgo/**` | 新 `geom_int_int_ss.rs`(528) / `geom_int_int_ss_1.rs`(**2368 ⚠**) / `geom_int_line_constructor.rs`(1186) + `int_patch/intersection.rs`(+27，补 `IntPatch_Intersection::Perform` 起点重载 = cxx **L2002-2035**)；计数等式 **30 = 30**、**21 = 21**；7 条未译清单见提交信息 | 六门槛绿 + 八网格 8/8；**无在役消费者**（`feat/loc_ope_split_drafts*.rs` 不在其域内）⇒ 无可观测变化属预期 |
| **B. `LocOpe_Gluer::Perform` + `AddEdges` 真身** | `824d7c63` | `feat/**` | `loc_ope_gluer.rs` 的 `perform`(L471-723)/`add_edges`(L781-873) 真身（原为 DEFERRED 空体）；计数等式 **17 = 17**；唯一残留替身 = `BRepExtrema_ExtPF`（**惰性**：OCCT 自己的 `AddEdges` 丢弃全部计算结果，cxx L549-551 两个 `if (flag==1) { }` 是空体） | 六门槛绿 + 八网格 8/8；featrf_a1 **逐字不变**（见下） |

**★ 两个必须记住的结论（本批次最有价值的信息）**：
1. **`LocOpe_Gluer::Perform` 的"DEPRECATED/待译"注释是**过期**的**：`LocOpe_WiresOnShape`（745+1515 行）、`LocOpe_GluedShape`、`LocOpe_Spliter`(1286)、`LocOpe_Generator`(1259)、`LocOpe::TgtFaces`（`loc_ope.rs:322`）**早已在库**——只有 `BRepExtrema_ExtPF` 真缺且它**惰性**。⇒ **教训：文件头的 DEFERRED 注释会过期；动手前先 grep 依赖是否已在库**（与 E3-W 追加 5 的 `BRepTools_Quilt` 教训同族）。
2. **`featrf_a1` 现在只剩**一个**阻塞点**：探针实测 `glued_f=2 the_ope=1 ope=Invalid collage=true` 且**零条 gluer 日志** ⇒ `LocOpe_Gluer::Perform` **根本没被进入**（没有一次 `Bind(Face,Face)` 成功——`glued_f` 的 key 不在 `myGShape` 里）。⇒ **粘合那一半已完成，剩下的是 §4.0 第 0a 项（BOP 历史镜像同一性）**。

**新 session 唯一要在本页做的收尾活（先做再往下推进）**：**已完成**（`8c052537`，见 §0 补记 5 收尾态与 port-plan 追加 15 补记 5）——
`geomalgo/geom_int_int_ss_1.rs` 已拆为 **1736 + 665 行**两文件（`geom_int_int_ss_1_curves.rs`），验收六门槛 + 八网格全绿、探针 = 0。

**（本轮已作废的）验收/提交配方留档**（本轮实际执行过一遍，有效）：`git status` 看清在飞改动 → **禁 `git add -A`**、按域显式 add → `cargo check -p rcad-kernel`/`-p rcad-algo` 零 error → 六门槛 + **重编 exe 后**八网格 → 形式抽查（计数等式 + OCCT 行号锚点 + 英文注释 + 单文件 <2000 行 + 探针计数 = 0）→ 按域逐批提交 + 同步根仓库指针（**不推送**）→ 不合格则按域整批 `git checkout --` 回退并把原因写进 §E3-W。

**⚠ 换行坑（会毁掉 review）**：`feat/brep_feat_rib_slot*.rs` 等文件的 blob **混有 CR**，用 Edit 工具改它们可能产出**整文件换行差异**（实测 3,946 行）。发现"整文件 diff"就 `git checkout -- <file>` 回退，把改动放进**新文件**；确实必须在原文件改时，在提交信息里单列一句"含 EOL 归一化"。（纯 LF blob 的文件——如 `feat/brep_feat_make_revolution_form.rs`、`bop/**`、`geomalgo/**` 多数——不受影响。）

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
| feat_featrevol | **1/45**（a5 **已通过**，2026-09-12 尾实测；其余 44 败） | 追加 15 尾修正：旧记档"0/45"**已过期**——`feat_featrevol_a5` 的真实断言**现在通过**（作者已用 path-scoped `git stash` + 重编证明**不是**批次 B 带来的；最可能来自本 session 早期的**扫掠子形状/面构造**修复，未逐条归因）。**未逐格复核，仅此一条已确认** |
| offset_shape_type_a / _i / _i_c | **0/1** / **0/12** / **0/19** | 全败（**逐例失败地图见 §0.0**；追加 19 后 `_i` 的 a3/a4/d2/d3 池外帧前进一层，`_a` 的 a4 推进到 `brep_algo/image.rs:159`） |
| offset_faces_type_i | **0/8** | 全败 |
| draft_angle | **0/49** | 全败（**追加 19 实测库内 panic 20** = 12 × `draft_modification_1_b.rs:1213` + 4 × `_1_c.rs:868` + 2 × `make_revol.rs:583` + 2 个新层次点 `brep_offset_api_draft_angle.rs:83` / `approx_int.rs:1071`；其余测试断言。追加 18 前为 28，追加 18 后 25） |
| thrusection_specific | **0/26** | 全败（thrusection 另见 §4.4） |

> **★ 口径提醒（追加 18 强化）**：**通过数不变 ≠ 无进展**。本轮所有批次都**没有**翻转任何网格的通过数，但把 6 例 offset + 3 例 draft 从"库内 panic"推到"测试断言"（失败层深度前进）。
> 逐例 file:line 的**失败地图**见 §0.0，跑法：`timeout 900 cargo test -q -p occt-generated-tests --test generated_occt_boolean_<grid> -- --nocapture 2>&1 | grep -E "^thread" | sed -E "s/^thread '([a-z0-9_]+)::.*panicked at ([^:]+):([0-9]+).*/\1 \2:\3/" | sort`。

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

## 3. 提交链与落地内容（**七轮**：追加 13 / 14 / 15 / 16 / 17 / 18 / **追加 19 = 本轮**）

**追加 19（本轮，2026-09-13；rcad `main` 顶尖 `fe297616`，根 `main` 顶尖 = 本文件所在指针 sync，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`fe297616`（**核心 4**：池外读取续链 —— kernel `edge_data_pool_free` + `build_curves3d.rs::edge_data` 守卫，十处池索引读收敛；`offset_shape_type_i` 池外 panic 清零）
← `e841e805`（**核心 3**：`GeomLib::BuildCurve3d` 家族 1:1 + 四处 GAP 载体收敛 —— `AdvApprox_PrefAndRec`、`ApproxAFunction::with_cut_tool` + 1D/2D 子空间存储、`kernel_curve2d` 桥、`GeomLib_CurveOnSurfaceEvaluator`、`GeomLib::isIsoLine`/`buildC3dOnIsoLine`、`GeomLib_MakeCurvefromApprox`）
← `830eb2c2`（**核心 2**：池外 `BRep_Tool::Curve` —— `topods.rs::curve_pool_free` + `shape_is_in_pool`，`brep_tool_curve`/`brep_tool_curve_loc`/`brep_tool_range` 按守卫分流）
← `d7beea3c`（**核心 1**：`ProjLib_HCompProjectedCurve` 在 `GeomPlate_BuildPlateSurface` 三处接线）
← `33ad8514`（= 追加 18 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `e841e805`）← `7528964` ← `030caf6`（追加 18）。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
域网格逐格通过数不变（失败层深度见 §0.7 / 追加 19）；探针 = 0。
**本轮新坑 28/29/30** 见 §6 尾（同 helper 行号掩盖层推进 / 旧笔记的"架构难点"会过期 / `GetType()` 常量返回决定分支）。

**追加 18（上一轮，2026-09-12；rcad `main` 顶尖 `05ddd5e1`，根 `main` 顶尖 `030caf6`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`05ddd5e1`（docs：追加 18 补记 2 —— 域网格 panic 分类 + 两处"真身已在库未接线"）
← `e8beac92`（geomalgo：`ExtendSurfByLength` 载体改为委托既有真身）
← `b2f2ab2b`（approx：OCCT `int` 度数字段改 `i32`（数据模型对齐）—— 修掉 `nbp - 5` 下溢）
← `e49dd9d2`（docs：追加 18 —— 翻译补全轮 + 失败层推进 + 剩余项依赖状态）
← `bf6eb0d6`（**核心**：translation round —— 5 处 GAP 真身 + 3 处重复删除）
← `4ab35d6d`（docs：追加 17 补记 3 —— TKOffset 第 1 项证据化）
← `480785e8`（fillet：**TKFillet 端盖弧 range 结案** —— `atan2` 实参顺序 + `reverse_curve` Circle 帧）
← `536f5a24`（kernel：`geom_proj_lib::curve2d` 解包 `Trimmed` —— a1 粘合路径打通）
← `39ed72f5`（docs：追加 17）← `0c211eb1`（**核心**：**D3 结案** —— BndWire 真 WireExplorer 连通序 + `BoundSortBox` int 语义 + `hinter_adaptor` 解包 `Trimmed`）
← `193da129`（= 追加 16 链尾）。
根 `main`：`030caf6` ← `810820d` ← `0f7452f` ← `9e2cdb0` ← `6ef64ba` ← `b7247e3`（追加 17/18 各批的 rcad 指针 sync）。
**本轮验证（全部在树）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
域网格逐格通过数不变（清单与逐例失败地图见 §0.0）；探针 = 0（rcad 与 OCCT 两侧都已清，OCCT DLL 已重建）。

**追加 17（上一轮，2026-09-12；rcad `main` 顶尖 `39ed72f5`，根 `main` 顶尖 `b7247e3`，均已推送）**
— rcad `main`：`39ed72f5`（docs：追加 17）← `0c211eb1`（feat+topalgo+kernel：BndWire 真 WireExplorer 连通序（D3 真因，含 `brep_feat_rib_slot_b.rs` 的 EOL 归一化）+ `BoundSortBox` C++ int 语义 + `hinter_adaptor` 解析访问器解包 `Trimmed`）← `193da129`（= 追加 16 链尾）。
其后同轮续推三轮（均含代码）：`536f5a24`（curve2d 的 `Trimmed`）、`480785e8`（TKFillet 端盖弧）、`4ab35d6d`（TKOffset 定界）。

**追加 16（2026-09-12；rcad `main` 顶尖 `193da129`，根 `main` 顶尖 `d23c335`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`193da129`（docs：追加 16 —— IntCurvesFace 轮结论 + D1/D2 + D3 任务书）
← `fe6b2750`（topalgo+bop：**D1** `fclass2d_topol` 的 wire 朝向复合（`TopoDS_Iterator` cumOri）+ **D2** `hinter_adaptor::curve_type_of` 的 Trimmed→基曲线）
← `cae1cd88`（topalgo+feat：**`IntCurvesFace_Intersector` 真身 1:1**，874 行 / 20=20 无桩；含 `NoExtFace` 假说证伪 + D1/D2 定位）
← `19369a85`（docs：追加 16 补注 —— 假说证伪 + featprism 结果根定为真因）
← `7e786fd6`（docs：追加 16 —— 三域并行轮）
← `fa0ddc43`（feat+bop：**0a 收尾** —— `BRepAlgoAPI::Shape()` = `aResult` Compound + `DS::argument_shapes` 历史翻译）
← `0263c85d`（fillet：**`ChFiKPart` 轴系 `YReverse` 左手系 + pcurve range** ⇒ `blend_simple_a1` 拓扑与 OCCT 全等）
← `514aefe8`（docs：追加 15 补记 5）← `8cd348f2`（feat+topalgo：FaceUntil 真 `MakeFace(S,Tol)` + outer-wire 槽 + Plane iso）← `ce5f01dd`（kernel：`add_to_compound` 穿句柄就地变异）← `8c052537`（geomalgo：`geom_int_int_ss_1` 拆分）。
根 `main`：`d23c335` ← `5cef39b`（两笔 rcad 指针 sync）。
**本轮验证（全部批次在树）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（375·378·379·373·12·102·83·110）；域网格逐格不变。

**追加 15（上一轮，2026-09-12，当时未推送）**——rcad `main`（**自上而下 = 新到旧**）：`824d7c63`（feat：`LocOpe_Gluer::Perform`+`AddEdges` 真身 = 队列第 0 项 (b) 半，**且其"依赖未译"注释经查过期**）
← `3f6f746a`（geomalgo：`GeomInt_IntSS` + `IntSS_1` + `GeomInt_LineConstructor` = 队列第 4 项第 3 条；⚠ 含 1 个 2368 行文件待拆分）
← `50f24743`（docs：补记 4 = 在飞批次交接）
← `3fa09a6a`（sweep+docs：删 `[ENTRY]` 探针 + 探针审计纠正）
← `4276bacb`（docs：补记 3 并行推进轮）
← `e7af4f95`（feat：`Geom2dAPIInterCurveCurve` 的 GAP 臂接真身 `Geom2dInt_GInter`）
← `fb4c8d94`（offset+geomalgo：`GeomAPI_ProjectPointOnCurve` / `GeomAPI::To2d/To3d` / `GeomLib::ExtendCurveToPoint` 真身 + `ExtentEdge` 真尾 + 4 处重复载体删除）
← `6f4b499a`（topalgo：`BRepCheck_Wire::SelfIntersect` / `BRepCheck_Face::Intersect` 真跑 2D 求交器 + 伪成功重复本体替换）
← `c9537dd0`（docs：补记 2 —— `featrf_a1` 根因链打到 OCCT 真值）
← `24ea2f77`（docs：追加 15 + 本交接）← `d280a867`（code：`featrf_a1` 链四处修复）← `e161078a`（= 追加 14 链尾）。
**工作树在本交接写就时是干净的**（两批都已提交）。
根 `main`：`sync: rcad pointer (… addendum 15 …)`（pointer → 本轮 rcad 顶尖；含 feat 对拍资产 `cases/featrf.hpp`）← `1686ccd`（= 追加 14 链尾）。**两仓库均未推送。**

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

> **★ 取项入口 = §4.6**（唯一在役队列）。§4.0–§4.5 是历史脉络（追加 15/16 时期），其中的"已定界"结论已被追加 16/17/18 更正，**只作背景阅读，不要按它们取项**。

### 4.0 追加 15 后的队列重排（**已归档，见 §4.6**）

**第 0 项 = 先验收 §0.4 的两个在飞批次**（提交或回退），再把下面按序推进。

**已完成（追加 15 同 session，三段）**：① `featrf_a1` 不终止 **清零**（根因 = `while(!FirstOK)` 里一处 OCCT 没有的 `fp`/`lp` 重算，**不是数值**）+ 同链三道墙（Transform 载体真身化、滑动 profile 面的 outer wire、子形状枚举空占位）；② 原第 3 项（`Geom2dInt_GInter` 通用批）⇒ **引擎 + 两个 analyzer 消费点落地**（并删掉一处"伪成功"重复本体），feat 侧 GAP 臂接真身（潜在路径真身化，feat 网格逐字不变）；③ 原第 4 项第 1/2 条（`GeomAPI_ProjectPointOnCurve` 真身 + `ExtendEdge` 真尾，含新译 `GeomAPI::To2d/To3d`、`GeomLib::ExtendCurveToPoint`）⇒ **落地并删掉 4 处重复载体**（offset 网格逐格零移动）。
**两个新定界（都有实测）**：原第 0 项的下一墙 = **两批活**（BOP 历史镜像同一性 + `LocOpe_Gluer::Perform`）；原第 2 项（`NoExtFace ×3`）= **`IntCurvesFace_Intersector::Perform` 欠计数**（b5 实测 `NbPoints(1)=1`，OCCT 需 ≥2）。

**新第 0 项（最急，a1 的直接下一墙；**粘合那一半已完成，只剩 0a**）**：

0. **`BRepFeat_RibSlot::LFPerform` 的结果装配**（`featrf_a1`）——**追加 15 补记 2 已把根因链打通到"要两批活"**（探针已清，权威记录见 port-plan 追加 15 补记 2/4）：
   - OCCT 真值：`OpeType()=LocOpe_FUSE` ⇒ `theOpe=1` ⇒ **走粘合路径** ⇒ `myShape = theGlue.ResultingShape()`（F=9 / area 109.511）。
   - rcad：**一次 bind 都没发生**（`glued_f` 的 key 不在 `myGShape` 里）⇒ `ope` 停在 `Invalid` ⇒ 回落构造器路径 ⇒ 空结果。
   - 再往下一层：`CutVehicle::modified(fac)` 在 rcad 返回**空**（OCCT 同一探针读 **`fac_in_cut_shape=0 n_mod=1 mod_in_cut_shape=1`**，即 OCCT 也是副本、但历史里有 1 个镜像**且该镜像就在结果里**）。
   - 最深层（实测）：rcad 的 `bop/algo/builder.rs::prepare_history` 只在"镜像 TShape 命中 `shape_remap`"时记 `Modified`；实测该 Cut 的源面为 **`n_imgs=2 imgs(ptr,in_remap)=[(x,false),(y,false)]`** 或 **`n_imgs=0`** ⇒ **一个都没进历史** ⇒ 同一性链断。**入口 = `fill_images_faces` 填 `my_images` 的点 / `shape_remap` 的构建点（`push_shape_recursive`/`build_result`），先证"镜像 ptr ∉ remap"是构造期就错还是查表期错。**
   - **两批活**：(a) **BOP 历史的镜像同一性**（`bop/algo/builder*.rs`，**核心代码，单独评审 + 八网格逐格复测**）——**新第 0a 项 = 唯一剩下的阻塞点**（实测：`glued_f=2 the_ope=1 ope=Invalid collage=true` + **零条 gluer 日志** ⇒ 没有一次 `Bind(Face,Face)` 成功，`LocOpe_Gluer::Perform` 从未被进入）；(b) **`LocOpe_Gluer::Perform` 真身** —— **已完成并提交（`824d7c63`）**，且其"依赖未译"的注释经查**是过期的**（`LocOpe_WiresOnShape` 等早已在库）。**⇒ 只需做 (a) 就能把 a1 推进到下一层。**
   - **对拍通道已建好（可复用）**：`tools/occt-bool-runner/cases/featrf.hpp` + `output/build_tkfeat.bat` / `build_runner_feat.bat`；跑法（含 **DLL 遮蔽检查**）见 port-plan 追加 15 补记 2。

**原第 1–4 项原样保留**（见 §4.1–§4.4 正文），其中第 1 项 = `BOPAlgo_Section` 对共面面片对的输出（`FalseSide ×5` 根因）——
**本轮已把它再定界一层**（探针已清，详见 §4.1 的补记）：b3 实测 `sc_pb=0 sc_v=0` + DS 的 FF 记录 `curves=0 points=0`，
⇒ 墙**不在** `build_section`（它逐字取 `PaveBlocksSc`，而该集合的唯一生产者是 FF 曲线），而在**共面/同曲面面片对的 FF 分支**；
**下一手必须是 OCCT 侧对拍**（§5 配方 3：`BOPAlgo_PaveFiller` 的 FF 段 / `IntTools_FaceFace` 同曲面分支插桩，跑 featlf b3），
先拿"OCCT 在该对上产出了什么"（曲线 / SD / 公共边 PB）再改 rcad。

**★ 追加 15 补记 5 对第 0 项/0a 的更新（2026-09-12，本 session）**：0a 的**上游**已打通——
`featrf_a1` 的 Cut **首次真正切分**（4 处 1:1 修复，见 port-plan 补记 5 / 提交 `ce5f01dd` + `8cd348f2`）。
**(b) `LocOpe_Gluer::Perform` 真身**上轮已交付；**(a)** 的**入口换了**：不再是 `prepare_history` 的镜像同一性，
而是 **`CutVehicle::with_operation` 的根形状口径**（`feat/brep_feat_form_2.rs`）——
rcad 取「池里**最后一个** Solid/Shell TShape」当 `my_shape`，OCCT `BOPAlgo_BOP::BuildShape` L1106 是
`myShape = aResult` = **全部结果容器的 compound**。实测链（探针已清）：Cut 现在把 primitive 切成 3 个 solid，
但 `glued_f` 的 key 落在**另外两块**的 face 上 ⇒ `LocOpe_Gluer::Bind(Face,Face)` **零调用**（`[HIST-BIND]` 零行）
⇒ `ope=Invalid` ⇒ 回落构造器路径 ⇒ 空结果（`surface area ... got 0`）。
**下一手 = 让 CutVehicle 返回 `a_result` 的 compound**（`set_shape_from_shapes` 当前把容器平铺进池、丢了 compound 语义）；
CutVehicle 被 prism/revol/d_prism/pipe 共用 ⇒ **改它必须跑全套域网格复测**，属核心代码单独评审。

### 4.1 tkfeat —— 第 1 优先级：featlf 下游墙（**追加 14 后重写**）

**现状（追加 14 实测）**：featlf 真实断言 15 例中 **`is_done` 门通过 11**（追加 13 时 0），
**kernel panic 清零**（追加 13 时 12 例 panic 于 `topods.rs` 的池别名）。
profile 有效性问题（`InvalidPointOnCurve` / `NotClosed` / `UnorientableShape` / `NoSurface`）
**已全部关闭**——见 §0.1 六条修复。剩余失败**全部在下游**：

| 墙 | 例数 | 入口 |
|----|------|------|
| `FalseSide` | **5**（b3/c5/d7/d8/d9） | **根因已定界到"section 0 边"**（本轮探针实测）：`Propagate` 的 `BRepAlgoAPI_Section(fac, CurrentFace, approximation=true)` 在 **b3 返回 `nedges=0`**（Compound 在、无边），而 a3 的同类 section **有 1 条边且判定完全正确**。⇒ 墙在 **`BOPAlgo_Section` 对共面面片对的输出**，不在 `Propagate`。入口 = `bop/brep_algo_api/mod.rs::SectionOp` → `run_build_section_brep` → `BOPAlgo_Section.cxx`（重点：FF 重叠面片对） |
| `NoExtFace` | **3**（b5/b6/b7） | `BRepFeat_RibSlot::ExtremeFaces`（cxx **L747-1319**）；**追加 15 补记 3 已再定界**：b5 走**单边分支**的"少于 2 个交点"退出（OCCT **L1052-1060**），上游判据 = **L1008** 的 `ASI.IsDone() && ASI.NbPoints(1) >= 2`，实测 rcad **`NbPoints(1)=1`** ⇒ 根源在 **`IntCurvesFace_Intersector::Perform`（TKTopAlgo）欠计数**，不在 `ExtremeFaces`/`LocOpe_CSIntersector`（后两者与 OCCT 逐句一致）。下一手 = 对该类插桩 + featlf b5 对拍 |
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

### 4.6 ★★ 三域队列（2026-09-12，**新 session 从此节取项**；上一节 4.0–4.5 为历史脉络，其"已定界"结论部分已被追加 16/17/18 更正）

**★ 工作模式（用户指令，追加 18 起）**：**先完成代码的等价实现（近乎 1:1 的翻译），代码基本译完再开始调试/修测试。**
⇒ 取项顺序永远是：**先"翻译/接线"类，再"调试/定界"类**；遇到 `panic!("GAP…")` / `unimplemented!` **先按 OCCT 函数名 grep 全库**——
本轮两例（`ExtendSurfByLength`、`ProjLib_HCompProjectedCurve`）的**真身早已在库**，只是没接线（坑 7 家族：文案会过期）。

**工作方式（连续三轮验证有效，继续沿用）**：主代理做**关键路径**，同时派**文件域不相交**的子代理并行推进三域；子代理各自 `CARGO_TARGET_DIR`、**只 check/门槛、不提交、不改 docs**；主代理**统一评审批量提交 + 跑八网格**。任务书必须点明：① 该批的 OCCT 清单与落位；② **不该碰的邻域**；③ 已知边界；④ 换行坑。**验收必须比对"失败层深度"**（见坑 21）。

**当前三域"翻译/接线"类待办（按可开工度排序，全部已定界）**（**★ 追加 19 起：原第 1/2/3 项已完成**，见下表）：

| 原项 | 状态 |
|------|------|
| 原 1（池外 `BRep_Tool::Curve`） | **已完成（`830eb2c2` + 续链 `fe297616`）**：`curve_pool_free` + `edge_data_pool_free` + `build_curves3d.rs::edge_data` 守卫，该链上十处池索引读全部收敛 ⇒ **`offset_shape_type_i` 池外 panic 清零**。⚠ 仅剩：**写回侧**（`same_range` 的 `builder_*` / `edge_mut_inplace`）池外**无池可变**——本网格未触发（`check_same_range` 为真时跳过 `same_range`），但只要某 case 需要它就仍会撞上；届时按坑 17 的 `Shape::data` 就地变异另做，**不要硬凑**。 |
| 原 2（接线 `ProjLib_HCompProjectedCurve`） | **已完成（`d7beea3c`）**；三个消费点全落。blend a2/p8/p9 的新墙为 `proj_lib_h_comp_projected_curve.rs:449` 的 `D0` `Standard_DomainError`（**状态类**，排在其他翻译项之后）。 |
| 原 3（`GeomLib::BuildCurve3d`） | **已完成（`e841e805`）**。残留：`offset/brep_offset_make_offset{.rs,_e.rs}` 的自有 GAP 文案、`offset/brep_offset_offset.rs` 与 `fillet/chfi3d_perform_elspine.rs` 的 **`BRepLib::BuildCurve3d` 级**载体（**不同函数**，仍为 stand-in）。 |

**⇒ 追加 19 后的新队列**：
1. **TKOffset 第 2 项（根，最高价值）**：`BRepTools_Quilt::builder_make_shell` / `brep_algo/face_restrictor.rs` 的**产物入池** —— 入池后本轮整条池读/写链（含写回侧）一次性解开，且 `result_brep()` 的拓扑计数随之正确。**牵涉拓扑计数 ⇒ 改完全套复测**（八网格 + 全域网格）。
2. **TKOffset 第 1 项写回侧（残留）**：见上表"写回侧"一行；只在某 case 真触发 `same_range` 时才需要，**先等第 2 项（入池）**——入池后它自然消失。
3. **TKFeat 第 0 项**：`LocOpe_Generator::Perform` 的 `IsDone`（a1 的直接下一墙，见下）。
4. **TKFillet 接力项 A**：`blend_simple_a1` 的**无界 pcurve 附着点**（面积 −2e100 的门）；本轮未触及，仍是 TKFillet 唯一"翻译/定界缺口"。
5. **调试类（按用户口径排在翻译/接线项之后）**：blend a2/p8/p9 的 `D0` 域错误；`offset_shape_type_i` a1/a2 的 `EdgeInter: E2 carries no pcurve`；`draft_angle` 的 20 例库内 panic（12 + 4 + 2 同 OCCT 亦 raise 的状态类）。

**TKFeat（`feat/**` + `topalgo/int_curves_face_intersector.rs`）**
0. **★ 第 0 项（最急，a1 的直接下一墙，已定界到函数）**：**`LocOpe_Generator::Perform` 的 `IsDone`**（`feat/loc_ope_generator.rs` / `loc_ope_generator_b.rs`）。
   `geom_proj_lib::curve2d` 的 `Trimmed` 解包**已落地**（追加 17 补记 1）：`BRepFeat::IsInside` 恢复正常分类 ⇒ `LFPerform` 的粘合判定 `collage=true ope=Fuse`（与 OCCT 一致）⇒ 进入 `theOpe=1` 的粘合路径。**新墙**：`theGlue.Perform()` 后 rcad **`IsDone=false`** vs OCCT **`IsDone=1` / `resultFaces=9`**（面积 109.511）。OCCT 的 `myDone` 来自 `myDone = theGen.IsDone()`（`LocOpe_Gluer.cxx` **L230**，`theGen` = `LocOpe_Generator`）⇒ **下一手 = 查 `LocOpe_Generator::Perform` 为何不达**（对拍通道现成：`RCAD_WS_PROBE` + `output/build_tkfeat.bat`；OCCT 侧插桩点在 `LocOpe_Generator.cxx` 的 `Perform` / `myDone` 赋值处）。
1. **`Propagate` 的 `BOPAlgo_Section` 共面 FF**（`FalseSide ×5` 根因，追加 15 已定界）：`build_section` 逐字取 `PaveBlocksSc`，而该集合的唯一生产者是 **FF 曲线**；b3 实测 `sc_pb=0 sc_v=0` + DS 的 FF 记录 `curves=0 points=0` ⇒ 墙在**共面/同曲面面片对的 FF 分支**。**下一手必须是 OCCT 侧对拍**（在 `BOPAlgo_PaveFiller` 的 FF 段 / `IntTools_FaceFace` 同曲面分支插桩，跑 featlf b3），先拿"OCCT 在该对上产出了什么"再改 rcad。
2. **featprism 结果根**（`NoExtFace ×3` 的真因，④ 已证伪"欠计数"假说）：OCCT 侧 `nbshapes r` = **11 面 SOLID**（V12/E20/W12/F11/SHELL1/SOLID1），rcad 侧 `root_shape`（生成测试 helper）拿到 **3 面 Shell**。**先判定两条完全不同的修法**：是**结果池缺根**（`feat/brep_feat_make_prism.rs` / `loc_ope_*.rs` / `brep_feat_form*.rs` / `brep_feat_rib_slot*.rs`）还是**生成器 helper 取错**（`tools/occt-test-gen`）——判定依据 = 直接从结果 BRep 数出顶层 Solid 的面数是否等于 11。
3. **b6/b7**：`PtOnEdgeVertex`（`brep_feat_rib_slot_b.rs:1110-1115`）同一个 `my_sbase` 问题；D2 **已修**（追加 17 又补了 `basis_curve_of`），复测。
4. 其余既有队列（`NoFaceProf ×3`、`BRepTools_Modifier::Perform` GAP）见 §4.1；feat 侧实测失败地图见 §0.0。

**TKFillet（`fillet/**`）—— 端盖弧 range 与 `reverse_curve` **已结案**；**接力项 B（blend 侧 ProjLib 接线）已由追加 19 落地**；接力项 = ① 无界 pcurve（a1 的门）② blend 剩余失败点**

> **★ 追加 19 更新（2026-09-13）**：下面**接力项 B 的接线已完成**（`d7beea3c`：`geomplate` 三个消费点全落），
> 且其"⚠ 先核实通用底曲线 `D1` 口径"的顾虑**已证伪**（`ProjLib_CompProjectedCurve` 基类是 `Adaptor2d_Curve2d`，见 §0.7 提醒 1）。
> blend a2/p8/p9 的新墙 = `proj_lib_h_comp_projected_curve.rs:449`（`D0` 的 `Standard_DomainError`，状态类）。
> **接力项 A（无界 pcurve）仍是本轮之后 TKFillet 唯一未动的翻译/定界缺口。**
1. **已完成（追加 17 补记 2，1:1 落地 + 实测）**：`chfi3d_compute_curves` 的端盖弧 range 由 **270°（`[π, 2.5π]`）** 修正为 **90°（`[0.5π, π]`，与 OCCT 一致）**。真因**不是**追加 16 记的"`Vint.Dot(Vref) < 0` 镜像步未生效"（探针实测 `dot=−100`，镜像一直生效），而是 **`elclib_parameter_circle` 的 `atan2` 实参写反**（OCCT `ElCLib::CircleParameter` = `AngleWithRef` ⇒ `atan2(v·YDir, v·XDir)`；rcad 写成 `atan2(v·X, v·Y)` ⇒ 整体偏 π/2）。同批落地 `reverse_curve` 的 Circle `y_dir` 翻号（`gp_Ax2::SetDirection` **真身** **L548-571** 通用分支：X 保留、Y 翻号）。**两条必须同批**。
2. **★ 接力项 A（a1 的门，调试类）**：**定位无界 pcurve 的附着点**。两条修复落地后 a1 面积由 57328.76 变 **−2e100**（无界），失败断言行不变（L113）、拓扑仍全过 —— 即追加 16 §4.6 第 3 条那个**未定位**的无界来源（某边界线 pcurve 带自然 `Line2d ±2e100` 域，挂在**不属于 7 张结果面**的 face 指针上）。起点：`ChFi3d_ProjectPCurv` 的返回值 / `hbuilder_face` 之外的生产者。**只做这一项就能把 a1 推到面积断言。**
3. **★ 接力项 B（翻译/接线类，blend a2/p8/p9 的门）**：**接线 `ProjLib_HCompProjectedCurve`**。真身已在库（`geomalgo/proj_lib_h_comp_projected_curve.rs` + `_b.rs`，已注册），消费者仍是 `unimplemented!`：`geomalgo/geomplate/build_plate_surface.rs:1003`。OCCT 现场 = `GeomPlate_BuildPlateSurface.cxx` **L1746-1802**（"Comparing metrics of curves and projected curves"，~45 行）：`new ProjLib_HCompProjectedCurve(hsur, Curve, myTol3d, myTol3d)` → `Adaptor3d_CurveOnSurface AProj(ProjCurve, hsur)` → 逐参数 `D1` 比模长算 `Ratio`，越界即 `myIsLinear=false`。⚠ **`Adaptor3d_CurveOnSurface` 以通用 `Adaptor3d_Curve` 为底时的 `D1` 口径要先核实**（rcad 既有 `Adaptor3dCurveOnSurface::new` 收的是 2D 曲线载体）——**不要**用"直接用 ProjCurve 的 D1"近似蒙过去。
4. **blend 其余失败点（本轮实测，逐例）**：a3/a4 → `chfi3d_builder_2b.rs:583`；q1 → 测试断言 L963；q2 → `chfi3d_builder_c1.rs:1162`；q4 → `geom/bspline_ops.rs:329`；q7 → `brep_blend_walking.rs:143`；x1 → `base/convert/mod.rs:2008`（按用户口径，先译后调 ⇒ 这些先当"待分类"，别急着按断言修）。
5. `ChFiDS_CommonPoint::SetArc` 旧记档（"rcad 只调 6 次"）**未经运行时复核**，静态 grep 有 ≥15 个调用点 ⇒ 先做运行时计数再立卡。


**TKOffset（`offset/**`）—— ★ 追加 19 后：第 1 项前段（池外读取）与第 3 项（`GeomLib::BuildCurve3d`）**已落地**；接力项 = ① 第 2 项（入池，根）② 池外读取续链（读侧）③ 回绕/附着两簇**

> **★ 追加 19 更新（2026-09-13，勿重复立卡）**：下面第 1 项的**首个池外调用点已通**（`830eb2c2`：`curve_pool_free` + 守卫分流），
> 墙已下移到 `BRepLib::check_same_range`；第 3 项（`GeomLib::BuildCurve3d` 家族）**整批落地**（`e841e805`）。
> **下一手 = 第 2 项（`BRepTools_Quilt`/`FaceRestrictor` 产物入池）= 根**：入池后整条池读/写链一次解开。
0. **★ 已完成（追加 18，本轮；勿重复立卡）**：`GeomLib::SameRange`、`BRepCheck_Edge::Tolerance`、`BRepCheck_Vertex::Tolerance`、`ElCLib::To3d` 全集、`GeomLib::To3d` 五个真身；`ExtendSurfByLength` 接线；OCCT `int` 度数字段对齐。**实测**：`offset_shape_type_i` e1/e2（原 SameRange GAP）、e3/e4（原 Edge Tolerance GAP）、e6/e7 **全部离开库内 panic 落到测试断言**；`offset_shape_type_a` 的 a4 由 `GeomLib::To3d` GAP 推进到 `brep_offset_make_offset_c.rs:52`。**当前 `offset_shape_type_i` 仅剩 6 例库内 panic = a1/a2（`brep_offset_inter2d.rs:1057`，状态类）与 a3/a4/d2/d3（池外，第 1 项）**。
1. **★ 第 1 项（池外，行号已证据化）**：`BRepLib::build_curve3d` 的**首个池外调用点是 `brep_tool_curve`**（不是写回路径）。
   实测链：`make_offset_shape` → `encode_regularity`（`brep_offset_make_offset_e.rs:992`）→ `brep_lib_build_curve3d_edge`（`brep_offset_make_offset.rs:278`）→ `BRepLib::build_curve3d` → **`brep_tool_curve(the_brep, an_edge)`（`build_curves3d.rs:677`）** → `BRep::edge_curve_world` → `BRep::edge`（`topods.rs:1799`）→ **`index out of bounds: the len is 16 but the index is 18446744073709551615`**（池外 `usize::MAX`）；4 例 = `offset_shape_type_i` **a3/a4/d2/d3**。
   **模板已在库（勿新造）**：kernel `topods::curve_on_surface_pool_free(&Shape,&Shape)`（`topods.rs:3955`）就是为同一批 offset 池外形状写的池外适配器 ⇒ 按 `BRep_Tool::Curve(E,L,f,l)` 再加一个同族池外读取（读 `the_e.data` 的 `TEdgeData.curve` + range），**不要**在 offset 里另写一份；然后逐个确认后续池依赖（`brep_tool_curve_on_surface_index` / `same_range` / `check_same_range` / 写回）。
2. **`BRepTools_Quilt::builder_make_shell` / `brep_algo/face_restrictor.rs` 的产物入池**（第 1 项同 4 例的**根**，且会让下游 `result_brep()` 的拓扑计数正确）——牵涉拓扑计数，**改完必须跑全套域网格复测**。
3. **`GeomLib::BuildCurve3d` 的剩余分支**（`GeomLib.cxx` **L1051+**）：平面分支现已可用 `geom_lib::to_3d`；iso 分支要用 **`GeomLib::isIsoLine` / `buildC3dOnIsoLine`**（⚠ **勿**混用 `geomalgo/approx_curve_on_surface.rs` 里的 `Approx_CurveOnSurface` 同名静态副本）；尾部需 `AdvApprox_ApproxAFunction`（**已在库**）+ **`GeomLib_CurveOnSurfaceEvaluator`（缺，~40 行，需新译）**。
4. **回绕链产出 trimmed 面**（8 例 `is_done()==false`）：`SelectShells` → `Deboucle3D` 返回 null，被拒的是**从未回绕成 trimmed 面**的延伸面（边界在 r=2.4e6 / ±1e7）；`BRepOffset_MakeLoops::Build` 已被正确调用（`loops.rs:45`）⇒ 修点在 `BRepAlgo_Loop::Perform` / `BRepAlgo_FaceRestrictor`。
5. **offset/split 面的 pcurve 附着**（`EdgeInter: E2 carries no pcurve on F`，a1/a2）；OCCT 侧同一处也是失败路径（`Geom2dAdaptor_Curve` 的 null 分支）⇒ 真缺陷在上游面的 pcurve 建立。
6. `draft_angle` 的两簇（17 × `draft_modification_1_b.rs:1213` 与 6 × `_1_c.rs:868` 的 `expect("… null Geometry")`）= **EMap 缺几何**（状态类，OCCT 同处也会 raise）⇒ 属调试项，先别当翻译任务。


**跨域通用**：`BRepAlgoAPI_*` 包装层（`CutVehicle` 等）**不要重写**——按 OCCT 形式接线到 `bop/**`（TKBO）的既有真身；本轮 0a 的两处修复（`aResult` compound 根 + `DS::argument_shapes` 翻译）已经把这条链打通，后续 feat 形态（prism/revol/d_prism/pipe）历史语义随之忠实。

### 4.7 其它（低优先，随手）

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
   ⑤ **对拍时别拿两边枚举的"名字"直接比**：OCCT `LocOpe_Operation` 从 **0** 起（`FUSE=0, CUT=1, INVALID=2`），C++ 探针 `(int)ope`
   打出的 **0 是 FUSE**，而 rcad 侧打印的是 Rust enum 的 `Debug` 名 —— 本轮差点把"OCCT 走粘合、rcad 走构造器"这个**决定性结论**读反。
   对拍输出一律换成人可读名或统一映射（本轮 `[FEAT] glue ope=0 Collage=1` → 0=FUSE 是靠查 `LocOpe_Operation.hxx` 才定性的）。
   ⑥ **feat 域的对拍通道先建再用**：`occt_bool_runner` 原先只有 boolean 用例，feat 用例要自己加（`cases/featrf.hpp` 模式：
   DRAW 参数映射写在文件头 + `TKFeat/TKOffset` 进 `OCCT_LIBS`），配 `output/build_tkfeat.bat`/`build_runner_feat.bat`；
   **一次建好，后面对拍零成本**——本轮 featrf 的根因链（3 层）全靠它一次性拿到 OCCT 真值。
   ⑦ **文件头的 DEFERRED/待译注释会**过期**，动手前先 grep**：`feat/loc_ope_gluer.rs` 的头注列出 5 个"待译依赖"
   （`LocOpe_WiresOnShape` 745+1515 行、`LocOpe_Spliter` 1286 行、`LocOpe_Generator` 1259 行、`LocOpe::TgtFaces`、`BRepExtrema_ExtPF`）——
   本轮实测**前四个早已在库**，真缺的只有最后一个、而且它**惰性**（OCCT 自己的 `AddEdges` 把算出来的东西全丢掉，cxx L549-551 两个 `if (flag==1) { }` 是空体）。
   ⇒ 与 E3-W 追加 5 的 `BRepTools_Quilt` 教训**同族**：**"注释说缺" ≠ "真缺"**；先 grep 真身、再立任务书，否则会把一个 873 行的活错估成 ~2500 行。
   ⑧ **旧记档的域网格通过数会漂移**：`feat_featrevol` 记档"0/45"，本轮实测 **a5 的真实断言通过**（用 path-scoped `git stash` + 重编确认**不是**新批次的功劳）
   ⇒ **报数前先实跑**（§6 坑 7 的老话，本轮再次应验）。

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
17. **（追加 15 补记 5，最贵）`Arc::make_mut` 不能用来 re-host OCCT 的"穿句柄就地变异"语义**。OCCT 的
    `BRep_Builder::Add` 一族（`Add(Shape,Shape)` / `Add(F,W)` / `Add(Solid,Shell)` / `Add(Compound,S)`）
    都是**就地改 TShape、所有拷贝的句柄立即可见**；rcad 若用 `Arc::make_mut`，**只要调用方还留着同一
    TShape 的 Shape 句柄，它就会被克隆**，调用方的句柄此后永远看着**旧的那份**。症状极具迷惑性：
    **所有构造语句都执行过了**，但容器/面/壳**看上去是空的**——本轮 `LocOpe_BuildShape` 的 compound
    （`C` 空 ⇒ `Nbedges=0` ⇒ `BRepFeat::Tool` 返回 null solid）与 `BRepLib_MakeFace::Init` 的面
    （wire 全落 `inner_wires`、`outer_wire` 留空占位）都是这一类。**正确做法** = `Arc::as_ptr(&pool.tshapes[idx])
    as *mut TShape` 就地变异（与结果构建期的索引重指向同一套共享变异模型，单线程管线安全）。
    判"槽位有没有东西"仍按坑 0① 看**子形状类型**，不看 `is_null()`。
18. **ff（FaceFace）产出 0 曲线时，先量 `raw_lines()` 再怀疑求交**：`has_intersection()` 是 **MakeCurve 域裁剪之后**
    的结果，而 `raw_lines()` 是裁剪之前。本轮 8 对候选"0 曲线"实际是 **raw 1–2 条、被 `part_in_face_hole` 全部丢弃**；
    若只看 `has_intersection` 会误判成"求交器不支持该面型"。
19. **（追加 16）`Shape::is_null()` 的池外误判会**静默**吞掉合法形状**（坑 0① 的第二次实测）：任何**池外**构造的
    形状（`BRepSweepBRepBuilder` 产物、`BRepTools_Quilt` / `FaceRestrictor` 产物、新建的 compound 根）`index == usize::MAX`
    ⇒ 被 `is_null()` 判成 null。**要给下游长期持有的形状一律入池**（本轮 `aResult` compound 入池后 LFPerform 才 `done=true`）；
    确须池外时，判空必须看**子形状类型**。**同族症状**：`result_brep()` 返回 None、`topods.rs` 的 `BRep::edge()` 越界
    （`index = 18446744073709551615`）。
20. **（追加 16）DS 的 `clone_arguments` 契约必须被每个 DS 键消费者遵守**：主 DS **深拷贝**参数（防就地修改泄漏到调用方
    BRep），所以 `my_images` / 历史 / 任何 DS 键数据都用**克隆** TShape；持**原始**形状的 API 级消费者（`BRepAlgoAPI_*`
    包装层）必须经 `DS::argument_shapes` 翻译，否则**静默全空**——本轮表现为 `Modified()` 对**每一张**源面都返回空，
    症状与"几何根本没切"极像（差点误判为布尔输出质量）。
21. **（追加 16，最贵）并行批次验收必须比对"失败层深度"，不能只比通过数**：本轮子代理的批次让 `featrf_a1` 的失败点由
    面积断言（L866）**退到** init 断言（L77），而六门槛 / 八网格 / 三域网格的**通过数完全不变** ⇒ 按通过数自检**完全不可见**。
    另外：**主代理做 A/B 归因前先 `git status` + `grep` 确认被测路径是否正被子代理改写**（本轮两次 A/B 都被在飞改动掩盖，
    差点把自己的修复判成回退）。
22. **（追加 17）把 C++ 表达式搬进 Rust 时，先看**括号在哪**——`std::clamp(static_cast<int>(v) - 1, 0, res-1)` 是
    `clamp(x-1)`，不是 `clamp(x)-1`**：`Bnd_BoundSortBox.cxx` L577-590 的越界输入（无限直线的采样段盒）会让
    `static_cast<int>` 饱和到 `INT_MIN`、`- 1` 在 x86 上回绕成 `INT_MAX`，最后**被 clamp 吃掉**；rcad 逐字翻译会
    debug 下 panic，首版把 `- 1` 挪到 clamp 外则立刻产出 `[-1,…]` 的越界索引（症状是 `VoxelGrid::add_box` 索引 panic，
    极易误判成"网格分辨率不一致"）。正确形式 = `clamp_res((v as i32).wrapping_sub(1))`。**同族**：C++ 的 `int` 算术在
    Rust 里是 `wrapping_*`（debug 默认 panic），搬运表达式时必须逐算子核对。
23. **（追加 17）`GeomAdaptor_Curve::load` 的"载入基曲线"语义必须贯穿**全部**访问器，不能只修 `GetType()`**：
    OCCT `Load`（L239-255）对 `Geom_TrimmedCurve` **递归载入基曲线**（区间仍是被裁的那个），所以 `Line()` / `Circle()` /
    `Ellipse()` / `Bezier()` / `BSpline()` 与 `GetType()` **同源**。追加 16 的 D2 只修了 `curve_type_of`，于是
    `GetType()` 报 `Line` 而 `line()` 按**原始变体**匹配 ⇒ `Standard_NoSuchObject` panic，**8 例 featlf 死在这里**。
    ⇒ 凡"把某类型归一化到另一形态"的翻译，**同一族的每个访问器都要过一遍**；`curve2d`/`curve2d_simple` 是**同一族的第三个**
    （见 §4.6 第 0 项）。
24. **（追加 17）对拍通道是**资产**不是一次性开销**：`RCAD_WS_PROBE` 门控 + `build_<TK>.bat` + 拷 DLL + 遮蔽检查
    （`grep -ac WS-PROBE <dll>` > 0）这套流程在追加 15/16/17 之间**零成本复用**，本轮三处修复全部由它一次定位（含
    OCCT 真值 `SlidingProfile ProfileOK=1 nWireEdges=5`、`LFPerform glue IsDone=1 ope=FUSE resultFaces=9`）。
    **用完必须 `git checkout --` 还原 OCCT 源并重建干净 DLL**（本轮已做：`grep -ac WS-PROBE` = 0），否则下一个 session 会
    在"OCCT 输出里为什么有旧探针行"上浪费时间。
25. **（追加 17）`wire_explorer_edges` 这类"载体"取的是**存储序**，而 OCCT 的 `BRepTools_WireExplorer` 走**连通序**——
    两者在**布尔产出的 wire 上必然不同**（OCCT 自己的 `BW[0..3]` 存储序就不是连通序）**：任何"按 wire 顺序推进"的算法
    （`SlidingProfile` 的 BndEdge 环走、`LocOpe_*` 的 wire 环）都必须用真 WireExplorer（rcad 真身 =
    `fclass2d::order_wire_edges` + 本次的 `wire_explorer_order` 包装），**不能**直接遍历 `TWireData::edges`。
    症状极具迷惑性：wire 会多出几条边并出现**对角线**（`make_edge(LastPnt, pp)` 把绕远路的端点连起来）。

26. **（追加 18，最便宜的一类）`panic!("GAP…")` / `unimplemented!` 的**文案会过期**——立卡前先按 OCCT 函数名 grep 全库**：
    本轮两例的真身**早已在库**，缺的只是"接线"：`GeomLib::ExtendSurfByLength`（真身在 `fillet/chfi3d_builder_c2_geomlib.rs`，
    offset/BRepFill 的载体却在 panic）、`ProjLib_HCompProjectedCurve`（`geomalgo/proj_lib_h_comp_projected_curve.rs` 已注册，
    而 `geomplate/build_plate_surface.rs:1003` 仍是 `unimplemented!`）。**做法**：`grep -rn "fn <OCCT 函数名小写>" libs/`
    + 看该文件是否已注册到 `mod.rs`；命中 ⇒ 这是"接线"任务（十几行），不是"翻译"任务（几百行）。与坑 7 同族（`BRepTools_Quilt` 教训）。
27. **（追加 18）OCCT 的 `int` 字段/局部量不要译成 Rust `usize`——有符号中间值是**设计的一部分****（与坑 22 同族）：
    `Approx_ComputeLine.gxx` **L1336/L1344-1350** 的 `int nbp = lpt - fpt + 1;` 与 `Mdegmax = nbp - 5;`（短点段时**为负**，
    紧接着被 `if (Mdegmax < mydegremin) Mdegmax = mydegremin;` clamp 回来）——rcad 用 `usize` ⇒ **debug 下 panic 且丢语义**。
    同类已落地三处：`ComputeLine::{mydegremin,mydegremax}`、`GeomInt_WLApprox::{myDegMin,myDegMax}`、`BRepApprox_Approx::SetParameters` 的
    `deg_min/deg_max`（连同 `IntTools_FaceFace::ApproxParameters` 的局部量）全部改 `i32`，只在以 usize 为索引的 rcad API 边界
    显式 `as usize`。**另**：C++ 表达式搬运时别挪括号——`std::clamp(static_cast<int>(v) - 1, …)` 是 `clamp(x-1)` 不是 `clamp(x)-1`（坑 22）。

28. **（追加 19，最便宜也最容易误判）失败层推进**可以**完全不改变失败地图的行号**——必须看 backtrace**：
    批次 2（池外 `BRep_Tool::Curve`）把 `offset_shape_type_i` a3/a4/d2/d3 从 `brep_tool_curve` 推到
    `BRepLib::check_same_range`，而两者**都**panic 在同一个 `topods.rs:1799`（`BRep::edge` 是共用的越界点）
    ⇒ 只比"逐例 file:line"会得出**"毫无进展"的错误结论**。
    **做法**：对池外/共通 helper 的失败，`RUST_BACKTRACE=1 cargo test … -- --nocapture` 取调用帧
    （本轮帧序：`BRep::edge` ← `BRepLib::check_same_range` ← `BRepLib::build_curve3d` ←
    `brep_lib_build_curve3d_edge` ← `encode_regularity` ← `make_offset_shape`）。坑 21 的加强版。

29. **（追加 19）旧笔记的"架构难点"也会过期，落线前必须回查 OCCT 头文件的基类**：
    追加 18 记着 `GeomPlate` 的 `Adaptor3d_CurveOnSurface AProj(ProjCurve, hsur)`"难点在通用 `Adaptor3d_Curve` 为底"，
    实测 `ProjLib_CompProjectedCurve` 的基类是 **`Adaptor2d_Curve2d`**（`ProjLib_CompProjectedCurve.hxx` **L37**）
    ⇒ 走的是**既有的 2D 载体口径**（`Adaptor3d_CurveOnSurface.hxx` **L43-44**），零架构成本。
    **与坑 7 / 坑 26 同族，但方向相反**：那两条说"注释说缺 ≠ 真缺"，这条说"注释说难 ≠ 真难"——
    **共同做法都是：动手前先读 OCCT 声明/基类/grep 全库，别信笔记的转述。**

30. **（追加 19）OCCT 的 `GetType()` 常量返回会决定 `EvalKPart` 的分支，translation 前先确认**：
    `ProjLib_CompProjectedCurve::GetType()` **恒返回 `GeomAbs_OtherCurve`**（cxx **L2108-2111**）⇒
    `Adaptor3d_CurveOnSurface::EvalKPart` 的 Plane 分支只可能落到 `myType = Other`（既不进 Circle 也不进 Line），
    `EvalD1` 因此**恒走通用链**（`myCurve->D1` + `mySurface->D1` + `SetLinearForm`）——kernel `CurveOnSurface::d1` 正是这条。
    **不要**为了"看起来更全"给这个消费点补 Line/Circle 分支；OCCT 在这里根本不走。

## 7. 资产位置（本轮新增/更新）

| 资产 | 路径 |
|------|------|
| **追加 19 新增**：`ProjLib_HCompProjectedCurve` 三处接线（metrics 比较 / ProjectCurve / ProjectedCurve） | `libs/rcad-algo/src/geomalgo/geomplate/build_plate_surface.rs`（metrics / 非线性回退）+ `build_plate_surface_b.rs::{project_curve, projected_curve}`；真身 `geomalgo/proj_lib_h_comp_projected_curve.rs` + kernel `base/proj_lib/adaptor.rs::CurveOnSurface` |
| **追加 19 新增**：池外 `BRep_Tool::Curve`（`curve_pool_free` + `shape_is_in_pool`） | `libs/rcad-kernel/src/topo/topods.rs`（紧随 `curve_on_surface_pool_free`）；消费点 `topalgo/brep_lib/build_curves3d.rs::brep_tool_curve`、`shhealing/shape_upgrade/unify_same_domain/topexp.rs::{brep_tool_curve_loc, brep_tool_range}` |
| **追加 19 新增**：`GeomLib::BuildCurve3d` 家族 | `libs/rcad-algo/src/geomalgo/geom_lib.rs::build_curve3d` + `geom_lib_iso_line.rs`（`is_iso_line` / `build_c3d_on_iso_line`）+ `geom_lib_make_curve_from_approx.rs`；kernel `math/adv_approx/pref_and_rec.rs`、`approx_a_function.rs::with_cut_tool`、`base/proj_lib/adaptor.rs::kernel_curve2d` |
| **追加 18 新增**：`GeomLib::SameRange` 真身（含 `Tolerance` 首参/守卫/periodic 分段） | `libs/rcad-kernel/src/geom/mod.rs::same_range_2d`；OCCT 签名入口 `libs/rcad-algo/src/geomalgo/geom_lib_same_range.rs::same_range` |
| **追加 18 新增**：`BRepCheck_Edge::Tolerance` / `BRepCheck_Vertex::Tolerance` 真身 | `libs/rcad-algo/src/topalgo/brep_check/brep_check_edge.rs::BRepCheckEdge::tolerance` / `brep_check_vertex.rs::BRepCheckVertex::tolerance`（`HCurveAdaptor::value` 新） |
| **追加 18 新增**：`ElCLib::To3d` 全集（8 个） | `libs/rcad-kernel/src/math/el.rs::elclib_to3d_*` |
| **追加 18 新增**：`GeomLib::To3d` 真身 | `libs/rcad-algo/src/geomalgo/geom_lib.rs::to_3d`（消费点 `topalgo/brep_lib/build_curves3d.rs::geom_lib_to_3d`） |
| **追加 18 新增**：OCCT `int` 度数（数据模型对齐） | `geomalgo/approx_int.rs`（`ComputeLine` / `WLineApprox`）、`geomalgo/brep_approx_approx.rs`、`bop/int_tools/face_make_curve.rs::approx_parameters_for`、`hlr/topo_brep/ds_filler.rs` |
| **追加 17 新增**：真 WireExplorer 枚举载体 | `libs/rcad-algo/src/topalgo/brep_top_adaptor/fclass2d.rs::wire_explorer_order`（消费点 `feat/brep_feat_rib_slot_b.rs::sliding_profile`） |
| **追加 17 新增**：`curve2d` 的 `Trimmed` 解包 | `libs/rcad-kernel/src/base/geom_proj_lib/mod.rs`（入口归一化） |
| **追加 17 新增**：`elclib_parameter_circle/ellipse` + `reverse_curve`（Circle 帧） | `libs/rcad-algo/src/fillet/chfi3d_builder_0.rs` |
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
| 权威长档 | `rcad/docs/tkfeat-fillet-offset-port-plan.md` §E3-W 追加 11–15（含补记 1/2） |
| **feat 对拍通道（追加 15 新增）** | `tools/occt-bool-runner/cases/featrf.hpp`（featrf A1 ground-truth driver）+ `main.cpp` 注册 + `CMakeLists.txt` 的 `TKFeat TKOffset`；脚本 `output/build_tkfeat.bat`（增量编/装 TKFeat Debug）、`output/build_runner_feat.bat`（编 runner）；跑法含 **DLL 遮蔽检查**（`grep -ac <探针串> TKFeat.dll` 必须 >0） |
| **在飞批次的隔离 target（可删）** | `rcad/target_a15_r3a`（批次 A `GeomInt_IntSS`）、`rcad/target_a15_r3b`（批次 B `LocOpe_Gluer`）；上一轮的 `rcad/target_a15_g2d` / `rcad/target_a15_offset` 亦可删 |
| 模块归属表 | `rcad/docs/module-map.md` §2/§3（`feat/`↔TKFeat、`fillet/`↔TKFillet、`offset/`↔TKOffset） |
| 重复实现审计 | `tools/occt-impl-audit/audit.py`（`uv run python tools/occt-impl-audit/audit.py .`）→ `docs/occt-impl-audit.md`（域风险排序：fillet 31 · feat 25 · offset 7） |
