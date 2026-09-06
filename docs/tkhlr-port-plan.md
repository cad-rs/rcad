# TKHLR 1:1 翻译推进计划（多 session runway）

> **交接快照（2026-09-07 session 19 终（含后续追加）—— 取代下方 session 18 终快照，两份都要读）**
> - **追加（同窗口第二次提交）：ptorus 数值链审计第一轮完成**——①发现并修复：`ToroidalSurface` 未覆盖 `derivatives2`，走 trait 默认的**有限差分**（h=1e-5，无 OCCT 对应物），而 OCCT 走 ElSLib 解析 D2（Contap_SurfProps::NormAndDn 的 default 分支消费）→ 已补解析 D2（eval.rs，对照环面参数化解析式）。效果：39 碎片/258.5（原 42/246.4），方向正确但非主因。②同轮已逐行核实 1:1：`Contap_SurfProps`（Plane/Sphere/Cyl/Cone/default 全部分支，含 torus 走 default）、`Contap_SurfFunction`（Value/Derivatives/Values/IsTangent/root）、`Cadrage`（gxx L413-591 全部钳制分支）、domain 传递链（uv_window→perform 的 BornInf/BornSup）。③剩余嫌疑（按序，下轮直接开工）：`math_FunctionSetRoot::Perform`（math_FunctionSetRoot.cxx L796-1100+，1452 行 vs rcad 特化版 949 行——停止准则/迭代流程/Rsnld.IsDone 语义是走线 bad-root 塌缩的关键，碎片化最像频繁 bad-root→PasC 减半→塌缩终止）；其次 `IsTangent` 的 `d <= GP_RESOLUTION` 阈值与 Contap tol 的 Set 来源链。验证锚：修完的标准 = 4 整线/302.685 逐位/隐藏 86.94；每改一步 `cargo test -p rcad-algo --lib acceptance_ptorus -- --ignored --nocapture` 看 filters 计数。
> - **提交链（本子模块）**：`7ce7245f`（apply_transform regularity 修复 + wd1/wd2 哑槽真对齐 + acceptance VCompute 树组装 + 9 文件探针删净）+ `<追加提交>`（torus 解析 D2）；根仓库指针同步。
> - **① CLOSED**：OCCT `nbshapes_hidden`（36/18/6）实为**整棵 resulth 树**的计数（可见 15 边/30 顶点 + 隐藏 3 弧/6 顶点 + 空 aCompHid 嵌套——VComputeHLR L3335-3350 无条件 MakeCompound aCompHid，空也非 null 也被计入；box 4/5 compound 同理）。acceptance 重写为 VComputeHLR 组装：`smoke_run_hlr_vcompute`（8 次抽取 OCCT 序）+ `vcompute_result_tree` + `nbshapes_counts`（BRepTools_ShapeSet 唯一计数）+ `assert_valid_digits`（1e-7 有效位线）。实测 run1=30/15/5、run2=36/18/6 与 OCCT 全等；mass 204.19032644821684 / 266.5264217926007 / 62.336095344383885（vs runner 17 位值 rel≤3.5e-8）。`acceptance_bug25813_1_vs_occt_json` 已 un-ignore。
> - **修复先在回归（stash 实证非本 session 引入）**：`apply_transform` 不变换 `CurveOn2Surfaces.surface1/surface2` → 烘焙 fixture 融合后缝边 CN regularity 按曲面值匹配失败 → C0 → 缝线漏进 VCompound（222.153 = 204.190+17.963），`bug25813_1_hidden_compound_three_arcs` 在 HEAD+当时脏树下已红。OCCT 锚：BRepTools_Modifier.cxx L191-195（变换后用新面重建 aBB.Continuity）+ BRepTools_TrsfModification.cxx L441-449（Continuity=原始值）。注意：纯 HEAD（不含他 session 脏文件）**编译不过**（algo_ext/mod.rs 脏版 TEMP-EXCLUDED fillet）——绿基线一直依赖脏树。
> - **② 前提反转（session 16/17 的关键误判，本次亲验 gxx + git 历史定论）**：OCCT `IntWalk_IWalking::Clear()`（gxx L110-137，注释头即 "adds dummy data to maintain start index of 1"）对 wd1/wd2 追加 `etat=-10` 哑槽、nbMultiplicities 追加 `-1`——0-based LinearVector 变成 [dummy, p1..pN]，**所有 `I=1..=n` 循环对齐且在界内**（`wd1[nbPath]` 读的是最后一个真点，无越界）。session 16 的「Reserve 仅容量 + 槽 0=第 1 点 + ComputeOpenLine 越界读」是误读；session 17 的「off-by-one 复刻」（去哑槽 + etat=0 尾槽）是建立在误读上的补偿 hack。已按真 OCCT 修复：clear() 注哑槽、删 with-interior 两处尾槽、删 3 处自创 `i < len` 守卫（TestArretPassage×2/TestPassedSolutionWithNegativeState）；no-interior Perform 无 Clear（gxx L317-391 实证）→ 保留其尾槽编码新对象越界。ToFillHoles 重建哑槽段（rcad 原本就对）不动。
> - **② 现状（真对齐下的诚实输出）**：ptorus 走线 42 碎片 / 可见 246.4 / OutLineH 0（= session 16 修复前状态）；生成测试 ptorus 1e-2 门槛转红（3/4 绿）。**gxx 的 wd 逻辑已逐行审计确认文本对齐**（代理全量对照：填充内容/全部读取/阈值/PWalking 分发全 MATCH），剩余分歧在 wd 数组下游的**数值链**：`IntWalk_PWalking.cxx`（从未逐行审计）与/或 Contap_SurfFunction / math_FunctionSetRoot 的 Root/导数值——导致走线提前终止（总长 246.4 vs OCCT 389.62，缺 ~143；42 碎片 vs 4 整线）。下一 session 首任务：IntWalk_PWalking.cxx + function_set_root/contap 函数层逐行审计，配 OCCT 侧单次插桩取走线步进真相（纪律：同命令 build+cp+grep+run 后还原）。
> - **③ 完成**：9 文件在册探针删净（classify/update/hider/internal_algo/shape_to_hlr/i_walking/contap perform/ds_filler/out_liner；grep 零残留）；acceptance bug25813_1 un-ignore（ptorus 保持 ignore，理由已更新）；回归全绿：algo lib **368+1i**（364→368：+3 回归测试绿 + acceptance bug25813_1 un-ignore）/ builder **76+1** / pavefiller **26** / tktopalgo **36** / tkbo **40** / tkgeom_algo **134+1**；生成测试 exact_hlr **3/4**（bug25813_1=204.1903 逐位、box 隐含、ptorus 红=诚实态）。kernel 663 绿。
> - **验收线（沿袭）**：AGENTS.md 有效位标准（逐位）；1e-2 只是中间里程碑——ptorus 的 303.867@1e-2 曾是补偿态假阳性，快照原文「吻合属覆盖/隐藏互相补偿」今日被证实。
> - **纪律提醒（沿袭）**：4 个他 session 遗留脏文件勿 add（algo_ext/mod.rs、fillet/topopebrepbuild.rs、tests/tkgeom_algo_gtests.rs——注意纯 HEAD 缺脏版 mod.rs 不可编译）；tkfeat-fillet-offset-port-plan.md 保持未跟踪；方法论 = 静态逐行审查 + OCCT 侧插桩取真相，rcad 侧不加探针。
> - **验证命令（不变）**：生成测试 `cargo test -p occt-generated-tests --test generated_occt_boolean_hlr_exact_hlr`（根仓库跑）；acceptance `cargo test -p rcad-algo --lib acceptance -- --nocapture`；回归 lib/builder/pavefiller/tktopalgo/tkbo/tkgeom_algo（清单见上）。

> **交接快照（2026-09-06 session 18 终 —— 取代下方 session 18 中与 session 17 快照）**
> - **提交链：rcad `8f97ffd6`（AngleWithRef/IsParallel 对齐）→ `30ee75bf`（Classify 四门 15 项修复 = gap② 关闭）→ `27d842ae`（回归测试+快照）；根仓库 `1fd0bd6`/`c8b5169`/`6f1fa17`（指针 ×2 + occt-test-gen 折叠修复）。**
> - **② CLOSED（本 session 收口）：bug25813_1 visible = 204.19032644821687（OCCT 204.19032694460864，8 位一致）、hidden = 62.336095344383885（OCCT 62.336097499736837）、隐藏 3 弧（大底缘远半 25.2237 + z=30 内缘远半 20.179 + 大顶缘远上中段 16.9334）✓、V 11 弧 + 4 轮廓线 ✓、rim 在 x'=±8/y'=27.959 处分裂与 OCCT 记录一致。**
> - **真正的根因（与 session 18 中段快照所记不同，最终定论）**：OCCT `Classify` 的全部四个 Z 门（Data.cxx L2026-2040 lf=1、L2059-2073 sta、L2087-2101 end、L2159-2173 mid）**只有 15 项比较——统一缺失第 16 项 `(MinMaxVert.Max[7] − iFaceMinMax.Min[7])`**（OCCT 原码的不对称；SimplClassify L2303-2317 则是完整 16 项）。rcad 翻译时"补对称"成 16 项：大顶缘中点的量化深度 tick 0x3231 < 面 4（小圆柱侧面）盒 z-min tick 0x3478，借位恰发生在多出的那项 → rcad 判 Out → 16.93 不隐藏；OCCT 该项不存在 → 直通 CS（2 点全接受 IN lvl=2）→ 隐藏。已修：classify.rs 四个门删第 16 项（保留 SimplClassify 的 16 项），`CLASSIFY edge=1 → state=In level=2` 与 OCCT CLR#65 逐位一致。探针佐证链：X17（隐藏通道+CS 值）、X18（myOX=−0.7854、D1.Z=+0.577）、X19（Load 时 IntL w_edge ori=INTERNAL/Int=1）、X20（Hide 时已重写 F/R——OrientOutLine 机制属实但非本 bug 直接因）、X21（面盒词两边逐位相同）、X22（缝边 reg1=regn=1 rg=6=CN）。日志在 `temp/x1*.log`、`temp/x2*.log`；OCCT 源码每次探针后均已 `git checkout` 还原。
> - **新增回归测试（rcad-algo lib 367+2i，+3 新测试）**：`bug25813_1_hidden_compound_three_arcs`（全管线：隐藏 3 弧 + mass 62.3360975 + V+Rg1+OutLine mass 204.1903269 + RgN=42.4578 gap③ 见证）、`bug25813_1_intl_wedges_present_and_reoriented_after_update`（IntL wedge 注册+OrientOutLine 重写断言）、`bug25813_1_interference_loop_over_outlined_faces`（诊断：RJ 真值复现 HIT In@1.4289/3.2835 + CLASSIFY 门判打印）。注意：直接 `ds.init_edge()` 的诊断绕过 `init_bound_sort` 会使 MoreEdge phase-2（全局排序边表）为空——复现必须走 `algo.hide_shape(1)` 或补 `init_bound_sort`。
> - **缺口③ CLOSED（生成测试 4/4 全绿：bug25813_1 = 204.1903 达标）**：所谓"Iso 双算"其实是**生成 harness 的 fixture 构造差异**——`hlr_translate.rs` 对 `ttranslate` 一律发 `apply_transform`（烘焙几何），与 OCCT 的 location 语义不等价，其融合 BRep 使缝线泄漏进 Iso（222.153 = 204.190 + 17.963）。已修：`hlr_translate.rs` 记录圆柱基点、ttranslate 折进构造重发（cyl_bases map + 行替换）；`main.rs` 的通用 ttranslate 臂同步加了 Cylinder 折叠分支。实测各 compound：v=119.2747 ✓、outline=84.9156 ✓、rg1v=0 ✓、iso=0 ✓、rgnv=42.4578（**OCCT 的 DRAW result 由 toShowCNEdges=false 永不取 RgN compound（runner main.cpp L383/414），rcad 的 RgN 内容正确且不计入 DRAW 口径——"RgN 42.458 多余"的推断作废**）。harness 口径（V+OutLine+Rg1+Iso）= 204.1903 ✓。
> - **在册待修（非阻塞）**：Elips2d::Reverse 误译（需 sense 字段，43 处构造点；本例 D1.Z>0 未激活）、Parameter3d 的 `1e-15 vs DBL_MIN` 阈值差、三代理 ulp 级清单（glam recip-multiply vs 逐分量除法、Ax2 框架链、line_parameter 重复归一化等）。
> - **验证命令（不变）**：生成测试 `cargo test -p occt-generated-tests --test generated_occt_boolean_hlr_exact_hlr`（**4/4 全绿**：bug25813_1=204.1903、ptorus=302.685@1e-2）；acceptance `cargo test -p rcad-algo --lib acceptance -- --ignored --nocapture`；回归 lib 367+2i / builder 76 / pavefiller 26 / tktopalgo 36 / tkbo 40 / tkgeom_algo 134+1 全绿。OCCT 插桩纪律：同命令 build+cp+grep+run，用后还原。
> - **剩余（按序，下一 session）**：① acceptance 的隐藏/顶点 pin（OCCT 隐藏树 = 6 子 compound / 18 边实例 / 36 顶点 vs rcad 平铺 3 弧 / 6 顶点；质量已逐位一致）——HLRToShape 结果组装（6 子 compound 的分片规则、added/Used/HideCount 簿记）的静态对齐。② ptorus：RgN=5→4 边、隐藏 compound ≈86.94、302.685 逐位（现 303.867，rel 3.9e-3）。③ 删净 9 文件在册探针、un-ignore acceptance、全量回归、module-map 终版。
> - **纪律提醒（沿袭）**：4 个他 session 遗留脏文件勿 add；tkfeat-fillet-offset-port-plan.md 保持未跟踪；rcad 探针清单（9 文件）删净待收口；方法论 = 静态逐行审查 + OCCT 侧插桩取真相，rcad 侧不加新探针；验收线 = AGENTS.md 有效位标准（逐位）。

> **交接快照（2026-09-06 session 18 中 —— 取代下方 session 17 快照，两份都要读）**
> - **② 的 OCCT 机制已完整钉死（本次四轮插桩 X17/X18/X19/X20，全部 build+cp+grep+run 同命令、用后 `git checkout` 还原；日志在根仓库 `temp/x1[789]_*.log`、`temp/x20_*.log`）**：
>   - **16.9334 的隐藏通道（X17）**：E=2（大顶缘半边）× 隐藏面 FI=4（**小圆柱侧面**，非顶环域！远上弧恰从小圆柱身后经过）的两条 IntL 轮廓线段 w_edge（ie=17/18）发生 2D 求交 → `RJ#3/#4: p1=1.4288992721907325/3.2834897081939571, p2=32.242640687119284, dz=-7.3484692283495345/-...363, TolZ=5.486e-4, st=IN` → RJF: `ori=2(INTERNAL), orie=FORWARD/REVERSED, decal=2` → 路由 ILHidden（**非 ILOn！session 17 的 ON 路由假设作废**）→ HidingStartLevel=0（tolpar=全边程1%）、两个交点都被保留 → EB.Builds(IN) 配对 [1.4289,3.2835] → `Classify(lf=1, mid=2.3561944901923448)` 无 Z 门拒绝、CS 2 点全接受 `state=IN lvl=2`（TolZ=bigSize*1e-6, wLim=-9.1555422664115138, w=-31.200949951460117/-11.605032009194696 全部 < wLim）→ `ES.Hide(1.4289,3.2835,onface=0,onbnd=0)`。**整个 exact pass 唯一一次区间 ES.Hide 就是它**；其余边全走 Compare→HideAll（该分支未插桩也无需）。
>   - **①③⑤前提修正（X18）**：`myOX = -0.78539816339744828`（负 π/4；cosinus=d3·d2=0.70710678118654757 比分支阈值 0.70710678118655 小 2.4e-15 → acos 分支）；**D1.Z=+0.57735026918962584 > 0 → `El.Reverse()` 分支本例未激活**；circle→Ellipse 判定走 `|D1·DZ|=0.577` 的 Ellipse 臂。
>   - **缺失机制 = OrientOutLine 方向重写（X19+X20 对照铁证）**：ExploreFace 时小圆柱侧面 wire2=[17,18] 记录为 `ori=2(INTERNAL), Int=1`（OutLiner 把自身轮廓线段以 TopAbs_INTERNAL 加入面副本的新 wire；Int 旗标 = IsIntLFaceEdge(F 原面,E)）。`Data::Update`（Data.cxx L841-845，`fd.WithOutL() && !fd.Side()` 时）调 `OrientOutLine`（L1746-1877）：对每条 `(OutLine||Internal) && !Vertical` 的 w_edge，取边中点/端点做 `UVPoint`→`CurvatureValue`（平面恒 0 → 不翻法向）→ `r=(Nm×Tg).z` → **`eb1->Orientation(ie1, myFEOri)` 把方向重写为 FORWARD/REVERSED（L1861）**；`OrientOthEdge` 处理其余边。X20 实证 Hide 时同一 w_edge 已变 `ori=0/1(F/R), Int=1` → NextInterference 的 F/R 守卫（L1277）放行 → 求交。**前 session "Data::Intersect 候选守卫不与 internal w_edge 求交、机制待识别" 的悬案就此关闭：守卫本身 1:1，缺的是 Update 阶段的方向重写生效。**
> - **rcad 侧已修（本 session 提交）**：`angle_with_ref` 补齐 OCCT `gp_Dir::AngleWithRef`（gp_Dir.cxx L58-83）的 asin 分支 + 去掉除长度/clamp（单位向量直接点积）；`D1.IsParallel(DZ)` 判据改为 OCCT Angle（acos/asin 分支）+ `ang<=tol || π-ang<=tol` 形式（gp_Dir.hxx L191-196）。修复后 bug25813_1 数值不变（239.086652），属正确性纠偏非成因。`El.Reverse` 误译（rcad 翻长轴，OCCT 翻 YDirection；Ellipse2d 需加 sense 字段，43 处构造点）与 `Parameter3d` 的 `1e-15 vs DBL_MIN` 阈值差为在册待修（本例未激活）。
> - **下一动作（gap② 收口，静态逐行）**：rcad 链条 `ds_filler.rs`（InsertFace `IntL.Append` L538 段：自身轮廓线段是否进 IntL）→ `out_liner.rs`（FaceHasIntL → wire append）→ `shape_to_hlr.rs`（SetWEdge Int 旗标，IsIntLFaceEdge 用**原始面**键控，data.rs L116-212 已核 1:1）→ `update.rs` L559 withOutL 链（已核 1:1）→ `classify.rs` orient_out_line（结构 1:1；`curvature_value`/`uv_point` 已核 1:1）——三层代理对照全部"结构对齐"，分歧必在**这条链上的值/旗标流**：小圆柱侧面在 rcad 的 FaceData 是否真的持有 IntL w_edge（Internal 旗标）→ withOutL 是否为 true → orient_out_line 是否写回。三个 Explore 代理报告全文在 session 记录（FaceData/Intersector/Projector 三方向），全部结构对齐 + 少量在册 ulp 级待修项。
> - **缺口③ 现值（acceptance 普查实证）**：seam 17.962924780409974（x'=5.657 的缝线）同时出现在 rcad 的 **RgNLineV 与（harness 口径的）Iso** → 双算 +17.96；OCCT nbIso=0 时 IsoLineVCompound 恒空。ptorus：acceptance RgN=5（OCCT 可见=4 条轮廓边、无缝边）、OutLineV=14 条 Other——结构差距未动。
> - **验证命令（不变，根仓库跑）**：`cargo test -p occt-generated-tests --test generated_occt_boolean_hlr_exact_hlr`；rcad 内 `cargo test -p rcad-algo --lib`（364+2i）+ builder 76+1 + pavefiller 26 + tktopalgo 36 + tkbo 40 + tkgeom_algo 134+1；acceptance 普查：`cargo test -p rcad-algo --lib acceptance -- --ignored --nocapture`。OCCT 插桩纪律沿用（同命令 build+cp+grep+run；`git checkout -- src/ModelingAlgorithms/TKHLR/` 还原）。
> - **纪律提醒（沿袭）**：4 个他 session 遗留脏文件勿 add；tkfeat-fillet-offset-port-plan.md 保持未跟踪；rcad 探针在册待删（9 文件清单已核：classify/update/hider/internal_algo/shape_to_hlr/i_walking/contap perform/ds_filler/out_liner 的 RCAD_HLR_TRACE/IWALK_DEBUG/SFDBG/[PFDBG]/[CWDBG] 块）；方法论 = 静态逐行审查 + OCCT 侧插桩取真相，rcad 侧不加新探针；验收线 = AGENTS.md 有效位标准（逐位）。

> **交接快照（2026-09-06 session 17 终 —— 取代下方 session 16 快照，两份都要读）**
> - **已完成（本 session，rcad 提交链 `6f8cc300`→`07154e54`→`8a1815a7`→`f751245e`）**：① **IntWalk LinearVector 0-based off-by-one 复刻完成**（clear() 去 dummy、wd1/wd2 填充后各补一个 etat=0 尾槽编码 `I<=nbPath` 的越界读、ToFillHoles 的显式 dummy 照抄；pnts1 侧保持错位）→ **ptorus 生成测试转绿**：visible 303.863 vs 302.685（rel 3.9e-3 ≤1e-2）。全基线零回归：algo lib **364+2i** / builder **76+1** / pavefiller **26** / tktopalgo **36** / tkbo **40** / tkgeom_algo **134+1**。
> - **⚠️ 精度现状（AGENTS.md 有效位标准，防误读）**：ptorus 的 1e-2 通过**不是达标**——303.863 vs 302.685 在第 3 位有效数字分叉，且结构仍错（17 条碎片 vs 4 条边、隐藏 compound 0 vs ≈86.94，吻合属覆盖/隐藏互相补偿）。AGENTS.md 要求精确到 OCCT 参考值末位，即最终必须：ptorus 4 边 302.685 逐位、bug25813_1 204.19 逐位、box 已逐位（146.969/48.990 达标）。当前 rcad 侧所有断言（生成测试与 acceptance 的 assert_rel）都是 1e-2 官方容差层，有效位层判定依赖 acceptance.rs 的边数/隐藏 pin 与逐位对拍，收口时三者都要过。
> - **② bug25813_1 侦查进展（当前 visible 239.09 harness 口径 / 221.12 acceptance 口径 vs 204.19；差 = 16.93 大顶缘中段 + 17.96 seam 双算）**：
>   - **② 机制侦查（本 session 末；⚠️ 含一次未解的构建间证据矛盾，如实记录）**：
>     - **输出侧铁证**：OCCT 隐藏 compound 含 16.9334（[DBG-E2] len=16.9334 sta=0.6435 end=2.4981 = 3D [1.4289,3.2835] 经 Parameter2d 偏移 ox=0.7854 的像）→ **OCCT 确实隐藏了中段**，rcad 缺失的就是它。
>     - **门判定证据自相矛盾（两次构建）**：run5（DBG-C+DBG-G 构建）：`classify(E=2,iFace=4,lf=1,2.35619)` → **CS nbpts=2 → ST=0(IN) → EB-IN state=0 → ES.Hide 执行（隐藏成功）**。run6（外部 revert 后仅 DBG-G2 重建）：同调用 **t2=0x80008000（w7: vMax word − fMin word 触发借位）→ Z 门拒绝 → 返回 OUT**。两次构建仅差 cout 打印，逻辑应同——矛盾未解，疑点：① 同形 (E=2,iFace=4,lf=1) 调用有多处（ILHidden-IN 通道 + ON 通道 + …），print 与 gate 的归属可能被 grep 序列错误配对；② revert 时序造成两次构建源不同的可能性无法完全排除。
>     - **静态关键线索（已核实代码）**：跨越点 State 由 **`RejectedPoint`（Data.cxx L2330-2368）** 判定：`TolZ = myBigSize*1e-5`；`dz = LE.Z(p1) − FE.Z(p2)`；`dz ≥ TolZ → above 拒绝`；**`st = (dz ≤ −TolZ) ? IN : ON`** → IN 进 ILHidden（classify 门控），**ON 进 ILOn → Hider L753-794 `EB.Builds(TopAbs_IN)` 于 ILOn **无条件** `ES.Hide(...,true,false)`"on the Face"**（无 classify！）。大顶缘中段若走 ILOn 通道即被无条件隐藏——与 Z 门无关。rcad 的 ILOn 通道存在且 1:1（hider.rs L862-949），但 ILOn 里没有这个区间。
>   - **下一动作（收口，一次性决定性插桩）**：单次构建同时打印（a）Classify 的 top/gate 结果/CS/ST/返回值（带调用序号），（b）**全部 5 个 ES.Hide 调用点**的 (E, p1, p2, OnFace, OnBoundary)，跑 bug25813_1 一次，把 [0.6435,2.4981]（2D 口径 16.9334）的隐藏归到确切调用点与通道；同时打印 RejectedPoint 对该交点的 (dz, TolZ, st) 验证 ON/IN 路由。然后静态对照 rcad 对应段修复。预期终态：16.93 隐藏 → bug25813_1 harness 239.09−16.93−(③ seam 17.96) → 204.19。
>   - **⚠️ 环境坑（本 session 两次实测）**：OCCT 工作树的 TKHLR 源码会被外部进程还原（刚加的探针文件被抹、`git status` 变干净），导致"构建二进制与源码不一致"的假矛盾（run5 与 run6 的门判定"矛盾"极可能即此因：run5 构建自旧源+探针，run6 构建自被还原后的源）。**下 session 插桩后必须立即在同一条命令里 build+cp+run，并在分析前用 `git status` + grep 探针字符串确认 DLL 与源一致**。runner 目录下的 build/Debug/TKHLR.dll 本地副本会遮蔽 PATH，必须重新 cp。
> - **③ seam**：bug25813_1 的 harness 口径比 acceptance 口径恰好多 17.96（小圆柱缝在 Iso 里被多算；acceptance 的 VisibleCompounds 不含 Iso）。OCCT nbIso=0 时 IsoLineVCompound 恒空 → 排查 rcad 该缝边的 iso 旗标（ShapeToHLR::ExploreFace 的 SetWEdge iso 入参 / HLRTopoBRep DS 的 IsoL 列表）。ptorus 的 RgNV=5 也需按"OCCT 可见=4 条轮廓边、无缝边"重新对照（acceptance pin result_edges()==4，当前 14）。
> - **验证命令（收口用，全部从根仓库 C:\Users\lilu\works\rcad-pro 跑）**：三用例 `cargo test -p occt-generated-tests --test generated_occt_boolean_hlr_exact_hlr`（box 隐含于 bug25813_1/3 之外——box 在 rcad acceptance.rs：`cargo test -p rcad-algo --lib acceptance -- --ignored --nocapture`）；回归 `cargo test -p rcad-algo --lib`（364+2i）+ `--test builder_stage_tests`(76) + `builder_stage_smoke`(1) + `pavefiller_stage_tests`(26) + `tktopalgo_gtests`(36) + `tkbo_gtests`(40) + `tkgeom_algo_gtests`(134+1)。OCCT 侧插桩调试循环：改 TKHLR 源 → `cd %OCCT_SRC%\build_debug && cmake --build . --config Debug --target TKHLR && cmake --build . --config Debug --target INSTALL` → **必须** `cp C:/Users/lilu/works/OCCT/build_debug/win64/vc14/bind/TKHLR.dll tools/occt-hlr-runner/build/Debug/TKHLR.dll`（本地副本遮蔽 PATH）→ `PATH="/c/tools/occt-debug/win64/vc14/bind:$PATH" tools/occt-hlr-runner/build/Debug/occt_hlr_runner.exe bug25813_1`（Debug 工具配 Debug OCCT；插桩后**同一条命令**里 build+cp+run，防外部还原）。收尾后 `git checkout -- src/ModelingAlgorithms/TKHLR/`（在 OCCT 仓库）。
> - **纪律提醒（沿袭）**：4 个他 session 遗留脏文件勿 add；tkfeat-fillet-offset-port-plan.md 保持未跟踪；rcad 探针收口时删净（清单不变）；方法论 = 静态逐行审查 + OCCT 侧调试工作流取真相，rcad 侧不加新探针；**验收线 = AGENTS.md 有效位标准（逐位），1e-2 官方容差只是中间里程碑**。

> **交接快照（2026-09-06 session 16 终 —— 取代下方 session 15 快照，两份都要读）**
> - **方法论纠正（本 session 用户裁定，永久有效）**：rcad 侧禁止"加探针+跑运行时对拍"定位断点（AGENTS.md 阶段 2：不 println/不 trace）；定位靠源码逐行静态审查；OCCT 侧真相按 AGENTS.md「OCCT Boolean Runner 调试工作流」取（插桩 TKHLR→Debug 重编→runner→`git checkout` 还原）。本 session 已按此执行。
> - **已提交（本子模块）**：`6f8cc300` —— ① `HLRBRep::MakeEdge` 补齐 Bezier/BSpline 分支（HLRBRep.cxx L63-158；CurveView 增 knots/multiplicities/weights；Curve 增 Poles/PolesAndWeights/Knots/Multiplicities 投影极点；is_rational 落地）——**ptorus 轮廓边此前是 BSpline 落在空臂返回 null，任何 compound 都不出现**，这是 session 15 "全部 hidden" 结论的真正成因（EdgeStatus 实为全可见）。② `Intersector::Perform(gp_Lin,P)` polyhedron 臂的角点临时变量遮蔽了函数参数 P（wLim）——OCCT 用大小写区分 p/P，rcad 改名 `pp`。③ Walking 容差对齐：`i_walking` 的 tolerance 从 `tol/u_extent` 域宽近似改为 OCCT gxx L244-245 的 `ThePSurfaceTool::UResolution(Caro, Confusion)` 解析公式（IWFunction 增 u_resolution/v_resolution；`surface3_u_resolution/v_resolution` 按 GeomAdaptor cxx L1818-1959 分型；contap SurfFunction 委托 SurfaceAdapter）。
> - **ptorus 当前值**：visible mass 246.40（42 条碎片轮廓边，OutV=42）+ RgNV=5 + RgNH/IsoH=1；回归基线全绿：algo lib **364+2i** / builder **76+1** / pavefiller **26** / tktopalgo **36**。
> - **OCCT 真值（已按调试工作流实测并还原插桩）**：ptorus（30/10，dir=(1,-1,1)）的 Contap 产出 **4 条走线**：nbpts 30/81/215/109，折线长 31.478/88.469/149.723/119.947（总长 389.62；解析积分两切线带 = 外 239.903 + 内 149.731 = 389.63，互证）；4 条线全部 HasFirst+HasLast（8 个 departure 点每线首尾各一）；DSFiller 每线 nbVertex=2 → **IntL 恰 4 条边**。viewer-path 可见结果 = **EDGE:4 / VERTEX:7 / Mass 302.685** ⇒ 隐藏 ≈ 389.62−302.685 = 86.94。**session 15 快照"近半 visible/远半 hidden"的前提作废**：结构目标 = 4 条轮廓边（acceptance 已 pin `result_edges()==4`），无缝边入可见 compound。
> - **唯一剩余断点（已钉死到语句级，下一 session 首任务）**：**复刻 OCCT `NCollection_LinearVector` 的 0-based 越界语义（off-by-one 怪癖）**。
>   - OCCT 事实（IntWalk_IWalking.gxx + NCollection_LinearVector.hxx 实证）：wd1/wd2/nbMultiplicities 由 `Reserve(n)`（仅容量）+ `Append`（从槽 0 起）填充 ⇒ **槽 0 = 第 1 个点**；而全部算法循环从 `i=1` 开始 ⇒ **第一个 departure/interior 点永远不被检查**（不作起点、不作停靠）；ComputeOpenLine L1486 `for (I=1; I<=nbPath; I++)` 还越界读 `wd1[nbPath]`（indeterminate；实践上不满足 >11/<-11 测试即跳过）。数据错位是双层的：`wd1[i]` 携带 Pnts1(i+1) 的 etat/ustart/vstart，而 `Pnts1.Value(I)` 处（L1494/L1744 等）取的是 Pnts1(I) 的 PathPoint——OCCT 自身混用。
>   - 作者意图证据：ToFillHoles 块（gxx L286-290）在 `wd2.Clear()` 后**显式 Append 一个 dummy** 再放 hole 点——即只有那里恢复了 1-based 语义。
>   - rcad 现状：clear() 预置哑槽把语义规整成干净 1-based（检查全部 8+20 个点）→ 线装配拓扑与 OCCT 不同 → 42 条碎片 vs OCCT 4 条长线。
>   - **复刻配方（纯机械，约 110 处下标点位）**：① clear() 不再预置哑槽（wd1/wd2/nb_multiplicities 空）；② perform 填充后给 wd1 追加一个**尾部默认槽**（etat=0 —— 越界读的确定性编码，两个测试 `>11`/`<-11` 自然不命中）；wd2 不加尾槽（OCCT 循环 `i<Size()` 无越界）；③ 全部 `self.wd1[i]/wd2[i]/nb_multiplicities[i]` 读取保持不变（数组内容已移位）；④ 对 `i>=len` 的 etat 写入忽略（OCCT 写 OOB 容量区，良性）；⑤ ToFillHoles 路径照抄 `wd2.clear(); push(dummy); push(holes)`；⑥ `pnts1[i-1]`/`pnts1[n-1]` 等 PathPoint 取址保持不动（保住错位）。验证：`cargo test -p rcad-algo --lib` 364 不回归 + ptorus 4 线/302.685 + box/cylinder 逐位不回归。
> - **gap②（bug25813_1 16.93）与 gap③（seam Iso→RgN）未动**，待 ①收口后按同链审查。注意 OCCT ptorus 结果无缝边可见——ptorus 的 RgNV=5 亦属多余，须与 ③ 一并按"可见 compound 结构"对照。
> - **纪律提醒（沿袭）**：4 个他 session 遗留脏文件勿 add；tkfeat-fillet-offset-port-plan.md 保持未跟踪；探针收口时删净（rcad 侧 RCAD_HLR_TRACE/IWALK_DEBUG/SFDBG/harness 的 CDBG/IWDBG/TDEF/OPL/BADROOT/NOTDONE 块全部在册待删）；OCCT 侧插桩已全部还原（`git checkout` 已执行）。

> **交接快照（2026-09-06 session 15 终 —— 取代下方 session 11 快照，两份都要读）**
> - **当前位置：Stage 0-3 + 4a + 4c 全关；4b 被无头路径取代（见 §7 勾选）；4d 进行中——box 验收全绿，bug25813_1/ptorus 两个验收测试已写好并 #[ignore]，剩余三个窄缺口（见下）。**
> - **提交链（本子模块）**：`c436f91b`(4c+4d+6修复) → `6f0b5439`(ptorus三连解锁) → `8e489252`(8环节静态审查) → `466cdb97`(容器变异器in-place) → `23dc22d2`(SplE机制还原) → `6cbebd64`(Int旗标实证+OutLineHCompound) → `ccd2d78a`/`a88011a0`(终极探针数据)。根仓库指针已同步。注意：曾误提交他 session 的 `docs/tkfeat-fillet-offset-port-plan.md` 后已 reset 撤出，该文件保持未跟踪。
> - **回归基线（当前全绿）**：algo lib **364+2i** / builder **76+1** / pavefiller **26** / tktopalgo **36** / kernel **663**。生成测试 `cargo test -p occt-generated-tests --test generated_occt_boolean_hlr_exact_hlr`（根仓库跑）：bug25813_1 现值未复测（修复前 239.09）、ptorus 现值 0.000（目标 302.685）、box 用例在 rcad acceptance.rs 全绿（146.969/48.990）。
>
> **剩余三缺口 + 收尾（按序）：**
>
> **① ptorus 轮廓边全部误判 hidden（visible mass 0，目标 302.685）——断点已钉死到数值级：**
> - 已打通（勿重查）：IntWalking 产 52 线（d2d 归一化修复后）→ Contour 采纳 51 → insert_face 建边入 IntL（51 条全 Wlk 型 BSpline，idx 9+，range [0,1]）→ process_face 加 INTERNAL 线到 NF → explore_face 旗标 **int=true 全部正确**（SFDBG 实证）→ DrawFace/DrawEdge/InternalCompound 逐行一致。
> - 实测症状：轮廓边 EdgeStatus **全部 hidden**（visible 区间空）→ OutV 空 → mass 0。缝边 RgN 2→5（SplE 分段+提取已生效，容器变异器 in-place 修复的实效）。
> - 已排除（勿重查）：Classify 三重门（param/end/0.4·sta+0.6·end）、TolZ 三分支（iFaceTest×OutLine/Internal→bigSize·0.01）、Shoot 方向（NZ=负Z，rcad 一致）、wLim-=TolZ、周期调整、w<wLim→classifier→IN、DrawFace/DrawEdge/InternalCompound（Used 重置/HideCount 延迟/leftover）、is_int_l_face_edge/is_spl_e_edge_edge、add_edge/add_to_edge/update_edge_pcurve（均 in-place）。
> - **收口动作**：单次运行抓一个近半轮廓片段的 classify `(w, w_lim)` 实测对拍——探针已埋：classify.rs 的 `trace_cl` 块打 `[TRACE] classify hits (u,v,w,in)`（注意：该 trace 只在 CS 实际执行后打印；torus 运行中**没有任何 classify: 行**——先查 Compare 回退路径（Hider.cxx L385 `Compare→Classify(LevelFlag=false)→HideAll`）是否为隐藏来源及其 trace 条件）。候选偏差：`elclib::line_parameter`（边上点射线参数化符号）与 `perform_line` 多面体求交 w 值。若 CS 无命中且 Compare 判 IN——对照 OCCT Data.cxx L2192-2262 的 Shoot/ElCLib::Parameter 符号约定。
>
> **② bug25813_1 差 16.9334**：combined 266.5264 已精确；可见侧缺的恰是 OCCT 隐藏的大顶缘中段弧 [1.4289, 3.2835]（3D 参，隐藏区间构建正确但被 Classify 盒门在 k=7 深度维拒绝——见追加 8 段矛盾分析）。与 ① 同属 Hider 深度语义，先复测当前 mass 再定。
>
> **③ seam Iso→RgN**：nbIso=0 时 OCCT IsoLineVCompound 恒空，融合体小圆柱 seam 被 rcad 计入 Iso（harness 会多算）。cylinder smoke 的 seam 归 RgN 正确 → 对照 bfuse 输出的 regularity（CurveOn2Surfaces/CN）是否丢失（bop 层）或 Data 的 Iso 旗标判定。
>
> **收尾机械项**：删全部临时探针（RCAD_HLR_TRACE：classify.rs/update.rs/hider.rs/internal_algo.rs；RCAD_IWALK_DEBUG：i_walking.rs/function_set_root.rs/ds_filler.rs/out_liner.rs/contour perform.rs；SFDBG：shape_to_hlr.rs；生成 harness 的 hlr_debug_compounds/census）→ un-ignore 两个 acceptance 测试 → 重生成 harness → 三用例全绿（204.19/302.685/146.969+48.990）→ §7 勾选 → module-map 终版 → 提交同步。
>
> **验证命令（从各自目录）**：根仓库 `cargo test -p occt-generated-tests --test generated_occt_boolean_hlr_exact_hlr`；rcad 内 `cargo test -p rcad-algo --lib`（364+2i）+ builder 76+1 + pavefiller 26 + tktopalgo 36。重生成 harness：`cargo run -p occt-test-gen -- --batch-boolean --batch-grid exact_hlr`（生成文件 gitignored）。
> - 纪律提醒：4 个他 session 遗留脏文件勿 add（algo_ext/mod.rs、fillet/topopebrepbuild.rs、tests/tkgeom_algo_gtests.rs、topalgo/brep_class/face_explorer.rs）；tkfeat-fillet-offset-port-plan.md 保持未跟踪；代码注释全英文；探针在收口提交时必须删净。
> - 下方 session 11 快照与其余追加段（含追加 8-13 的完整诊断链）仍有效，冲突处以本快照为准。
> **交接快照（2026-09-06 session 11 终）**
> - **当前位置：Stage 0-3 全关，HLR 精确管线全链贯通且盒体/圆柱端到端 = OCCT 真值逐位一致。剩 Stage 4d 收官（exact_hlr 断言闭环 + 4c 接入 + module-map）与两条残余记录项。**
> - 提交链（rcad 子模块 sd-hash-wip）：`7011f938` 计划落盘 → Stage 0/1 → 2a 全关 → `dc9cf575` 2b Contap 全包 → 3a → 3b 全关（`fe89ac91`/`729f0366`/`2f9b2d3e`）→ `2f1ac08e` Lin×Circ → `d1dc1ad1` 3c-1 → `432bff94` approx_int 泛型化 → `a0867176` conic 全臂 + BRepApprox 层 → session 8 `5e782450`/`79807080` **3c-2 全关 + Geom2dHatch 叶子 + 3e 干涉八件** → session 9 `91e5ff1f` **3d 全关（Hatcher 主引擎 + BRepTopAdaptor_Tool + HLRTopoBRep 全包）** → session 10 `5a6adf51` **3f+3g 全关（Data/Hider/InternalAlgo/HLRToShape/ShapeToHLR/ShapeBounds + 3e EIT Data trait）** → session 11 `38936196` **烟囱闭环（盒体=OCCT 真值）+ Q/R 修复 + MakeEdge2d + gen_hlr_ref.py** → `dc7de1e0` **缺口①②③（DomainIntersection/CurveOnPlane/会话 BRep）** → `81a03db1` **缺口④⑤ + Hider 圆柱远半弧全关**。根仓库：`a04c7da`/`4b08687`（tools/occt-hlr-runner 与 occt_hlr_*.json 属根仓库 tools/tests）。
> - 回归基线（全部全绿，以此为准）：algo lib **361/0 failed/0 ignored**、kernel **663**、builder_stage_tests **76** + smoke **1**、pavefiller_stage_tests **26**、tkhelix_gtests 16、tkgeom_algo_gtests 134+1、tkbo_gtests 40、**tktopalgo_gtests 36**（Trsf 字面量已按 OCCT SetScale/SetMirror/SetRotation 成员序 matrix=identity+scale 承载修复）。
> - **session 11 完成清单（九代理批次 P/Q/R/S + T/U/V/W + X/Y/Z + 主线程裁决/生成器）**：
>   1. **3f（K/L/M）**：`hlr/brep/data/`（Data<'a> 70 字段契约 + update.rs/classify.rs 2416 行/tableau_rejection.rs）+ face_data.rs + edge_face_tool.rs（UVPoint 7 参）；uv_point 契约裁决（Data.my_brep kernel 上下文 + my_e_map 槽位）。
>   2. **3g（N/O1/O2）**：hider.rs（Hide 主循环 L41-852 全量 + EIT Data trait impl）+ algo（组合+Deref）+ internal_algo.rs（HideSelected heapsort 逐行）+ hlr_to_shape.rs（含 HLRBRep.cxx 包静态）+ shape_to_hlr.rs + shape_bounds.rs。
>   3. **Q**：hider 双锚转绿——根因是 fixture 对象寿命（栈 Data update() 后被 Arc 移动，&myProj 悬空）；**ON 锚"失败"经安装版 OCCT 二进制求证为手工对拍错误**（EdgeStatus::Hide 的 ！OnFace 短路门），0 ignore。
>   4. **R**：fclass2d vmap → IndexedVertexMap（插入序，BRepTools_WireExplorer.cxx 依据）+ kernel `topo/brep_lib/MakeEdge2d`（669 行 1:1，14 锚点）。
>   5. **S（烟囱）**：盒体 box 10x20x30 端到端 = OCCT 真值（VCompound 9 边 146.96939 / HCompound 3 边 48.98980 / 12 线段解析像 1e-6）；圆柱部分锚定 + 5 缺口暴露。
>   6. **T**：圆柱双轮廓——真缺陷在 g_inter.rs 求**交器**缺 DomainIntersection 语义（IntCurve_IntConicConic_1.cxx L641-714：域外解丢弃/域端钳位 Head-End）；C++ 工具对 OCCT 二进制求证 6 探针 IN。
>   7. **U**：CurveOnPlane 回退——FClass2d init 无回退语句（与 rcad 一致），真回退在 BRep_Tool::CurveOnSurface **内部**（BRep_Tool.cxx L367-450）→ kernel curve_on_plane 落在正确层。
>   8. **V**：InternalAlgo arena 产线化（my_brep 会话捕获，smoke 走自然 add→update→hide，断言不变）。
>   9. **W（4a 收官）**：`tools/occt-hlr-runner/`（C++ 无头复刻 VComputeHLR L3299-3362 双分支；Release 构建实测）——box/bug25813_1 与 OCCT lprops **逐位吻合**；**ptorus viewer-path 302.685 与 case length 逐位一致**；**定论**：Plate/ptorus 的 tuple 路径崩溃 = OCCT 8.0.0 自身 HLRBRep exact-algo 数值稳健性 bug（headless DRAWEXE 同崩），非 rcad 问题。
>   10. **X（缺口④）**：seam 携带 **GeomAbs_CN 而非 G1**（BRepPrim_Builder.cxx L107-118）→ 归 RgNLineVCompound（typ 4）；kernel GeomAbsShape 真实枚举序 + CurveOn2Surfaces regularity + BRepBuilder::continuity + 真实现 BRep_Tool::Continuity；四个 primitive 补 CN。
>   11. **Y（缺口⑤）**：核实已在 3f/3g 关闭（"注释 OCCT 原文"是引用注释惯例，活代码在正下方）；加成对端点唯一性 pin 断言。
>   12. **Z（Hider 远半弧）**：两处 1:1 底层偏差——kernel `Lin::transform` 丢方向旋转（gp_Ax1::Transform 同转 vdir）+ `Curve::First/LastParameter` 缺 Parameter2d 包装（HLRBRep_Curve.lxx L66-76）；修复后 **cylinder smoke 与 OCCT 真值逐位一致（V=5 弧 75.671 / H=1 弧 25.2237 / outline 2 条 / seam 在 RgN）**。
>   13. **gen_hlr_ref.py（主线程）**：`tools/gen-occt-ref/gen_hlr_ref.py`——每 vcomputehlr 一次独立进程、**必须 `-b` batch**（交互 exit 不 flush C stdio，lprops 的 printf Mass 块丢失）、lprops 走 C++ stdout 而非 Tcl 返回、algotype 剥 Tcl 引号、默认视图等效元组 dir=(1,-1,1)/up=(-1,1,2)。产出 4/6 JSON（`tests/occt/step_reference/occt_hlr_*.json`；exact bug25813_1 自检 204.19 精确 PASS + 隐藏 266.526；poly×3 已合入 runner 值）。
> - **下一 session 入口（Stage 4d 收官，按序）**：
>   1. **4d 断言闭环**：新增 exact_hlr 的 rcad 管线断言测试（消费 `tests/occt/step_reference/occt_hlr_*.json` + `tools/occt-hlr-runner/ref_output.txt`——box 146.969/48.990 已由烟囱证明；bug25813_1 204.19；ptorus 用 viewer-path 302.685；Plate 用 406.283 viewer-path 或 CI 404.283 需记录差异），相对误差 ≤1e-2。
>   2. **4c**：occt-test-gen 接入 hlr 网格翻译（main.rs L6 mod、L463 batch、L3390 分发、DrawPort needs_hlr——4 个接入点）。
>   3. **module-map.md hlr 行更新**（对齐状态 + 基线数字）。
>   4. **残余记录项（非阻塞）**：① `brep_tool_is_closed_edge_face`（ExploreFace 的 Dbl/IsReallyClosed）按 face ptr_id 匹配而非曲面值匹配，对 OutLinedShape 面副本把 seam 的 Dbl 判 false（X 报告）；② 烟囱 fixture 的 pcurve 注册注释过时（U 报告，stored 优先行为不变）；③ hlr_to_shape 的 Hyperbola/Parabola 臂 null（HLRBRep_Curve 适配器 stub + CurveView 缺 2D 投影 poles/knots 访问器）；④ 参考脚本注释里的过时描述（fixture 属主更新）。
> - **HLR 消费链现状（全链真值验证）**：Algo::Add → InternalAlgo::Update（ShapeToHLR::Load → DSFiller::Insert → BRepApprox/Contap/FaceIsoLiner/Hatcher）→ Data.Update → Hider.Hide → HLRToShape。盒体 12 棱可见/隐藏归属 + 圆柱双轮廓/椭圆弧分段/seam RgN 归类/远半弧隐藏全部与 OCCT 一致。
> - **已知遗留（非阻塞，均已记录）**：Contap quadric-exact 死分支待 Adaptor3d_CurveOnSurface；landed approx_int 引擎的非 Bezier ComputeLine 分支未落；DSFiller insert 签名带 rcad &mut BRep arena 参数（显式化备案）；FaceIsoLiner pcurve key 的 Shape-only 约定；poly 线路 = Stage 5 单列（runner 已可产 poly 参考）。
> - **并行子代理模式（session 11 实证补充）**：主代理预落 struct 字段契约 + 骨架先行（多代理并发零注册冲突）；跨代理契约冲突经 SendMessage 一轮裁决（uv_point 7 参 vs 5 参调用点）；**"修正任务前提"型发现是常态**（Q 的 fixture 寿命、U 的回退层位、X 的 CN 值、Y 的已关闭核实）——代理必须以 OCCT 源码/二进制求证后才允许改期望，禁止凭推理改测试；X 的 git checkout 误回滚靠逐字重建重放 + 协调方核对恢复。
> - 已确立的翻译模式（沿例勿改）：gxx 模板参数 → trait；实例化 = 别名 + marker；gxx/lxx 内联；`// OCCT <文件> L<起>-<止>`；每叶子配解析锚点。**HLR 特例**：裸指针对/void* → 生命周期槽（Arc<dyn T + 'a>）/map 键 = ptr_id/handle 值拷贝 → Clone/OCCT 继承 → 组合 + Deref/try-catch → catch_unwind/kernel-context 显式参数（brep/edge Shape）。
> - 踩坑累积（续）：**(12)** 固有方法遮蔽 trait 方法——显式限定；**(13)** 引擎 value 分路径走访问器；**(14)** heredoc 大段必截断——Write 落盘 + cat 拼接；**(15)** mod.rs 预注册 + 骨架先行；**(16)** 追加注册防重复（E0428）；**(17)** DRAWEXE 交互模式 exit 不 flush C stdio——**必须 `-b`**；lprops 走 cout 不走 Tcl 返回；同进程第二次 vcomputehlr 损坏首个结果——**每 vcomputehlr 一次进程**；**(18)** 对象寿命：update() 植入内部指针后不得再移动宿主（先 Arc/Box 稳定地址）；**(19)** OCCT 枚举序（GeomAbs_Shape C0<G1<C1<G2<C2<C3<CN）与 SetScale/SetMirror 的成员承载（matrix=identity+scale）要以源为准，不许按"常识"近似；**(20)** gp::Trsf 增字段后必须同步 tests/ 下的字面量。
> - 子模块工作树遗留修改（非本任务）：`algo_ext/mod.rs`、`fillet/topopebrepbuild.rs`、`tests/tkgeom_algo_gtests.rs`、`topalgo/brep_class/face_explorer.rs`——提交时只 `git add` 本任务文件。
> - shell 的 cwd 会被重置到 `C:\Users\lilu\works\rcad-pro`（根仓库），编译/测试前先 `cd rcad`；回归必含 builder 分阶段（builder_stage_tests + pavefiller_stage_tests）+ tktopalgo_gtests。

> - **session 5 完成清单（3b 剩余三项全关）**：
>   1. `fe89ac91` **ContapDomain for BRepTopAdaptor_TopolTool**：`Curve2dAdaptor::as_any` = occ::down_cast 支撑（BRepTopolTool::Initialize/Orientation 下行转换限制弧句柄）；`HVertexBehavior` trait + `HVertexHandle` = `handle<Adaptor3d_HVertex>` 的多态句柄（Value/Parameter/Resolution/Orientation/IsSame 全虚，BRepTopAdaptor_HVertex 全覆盖）；BRep 侧四类（BRepCurve2d/BRepHVertex/FClass2dTopol/BRepTopolTool）整体从借用转 **Arc\<BRep\> 所有权**（HLRBRep_Data 所有权设计）；ContapDomain::edge() → Option\<Shape\>（OCCT void* 语义）；锚点 = 经 ContapDomain trait 驱动 BRep 域全迭代。
>   2. `729f0366` **IntCurveCurveGen 泛型化**：concrete Geom2dInt 绑定 → `IntCurveCurveGen<C, PT: PCurveTool<C>, IC, IPP>` 引擎（+子引擎成员接口 `IntConicCurveMember`/`IntPCurvePCurveMember` + gxx L1023 SetMinNbSamples）；GInter = 类型别名级实例化；**Extrema 定位机制泛型化**（LocatorCurveTool：GCurveLocator/GFuncExtPC/GenLocateExtPC + 共享 FindParameter 体 proj_cur_find_parameter_{bounded,unbounded}，Geom2d 签名零回归）。
>   3. `2f9b2d3e` **CInter 家族 + Intersector**：`hlr/brep/curve_tool.rs`（HLRBRep_CurveTool 静态工具，PCurveTool/ParTool/LocatorCurveTool 三面）；`c_inter.rs`（TheProjPCurOfCInter/TheIntersectorOfCInter/TheIntConicCurveOfCInter/IntConicCurveOfCInter/CPolyTool+TheIntPCurvePCurveOfCInter/CInter 别名级实例化）；**IntConicConic::Perform(Lin,Lin) 1:1**（_1.cxx L1381-2233 + 7 个文件静态助手 + 文件局部 TOLERANCE_ANGULAIRE=1e-15）；`edge_data.rs`（HLRBRep_EdgeData 位旗标全套）；`intersector.rs`（HLRBRep_Intersector 99/858：双边 decalage 循环/SimulateOnePoint/CS 分支含 polyhedron 参数窗剔除）。锚点 6 个。
> - **session 7 完成清单**：
>   1. `d1dc1ad1` **Stage 3c-1 落地**：`app_par_curves_bsp.rs`（AppParCurves_BSpFunction.gxx L21-337 → BSpParFunction；AppParCurves_BSpGradient.gxx L22-412 → BSpGradient 双构造 + lambda 搜索 + Rogers-Fog 投影 + BFGS 收尾，实现 GradientFunction 面，BSpGradient_BFGS = landed GradientBfgs 泛型壳直用）+ `bspl_compute_line.rs`（Approx_BSplComputeLine.gxx L139-1458 → BSplineCompute，按 AppDef_BSplineCompute_0.cxx 别名表绑定：四构造器 / Perform 切割循环 / Parameters 三参数化 / Compute 度数循环 + Interpol 兜底 / Interpol C2 三次插值含两点度 1 分支 / First-LastTangencyVector 抛物线回退 / FindRealConstraints；OCCT 无定义的 ComputeCurve 死声明省略并注释）。锚点 2 个（四分之一圆 17 点逼近达容差；两点 Interpol 度 1）。
>   2. `432bff94` **approx_int 引擎接缝泛型化**（零行为变更，worktree 隔离验证 196/134 前后一致）：`ApproxIntWLine` trait（TheWLine 缝）+ `ApproxIntMultiLine` trait（TheMultiLine+LineTool 缝，含 sub_line 构造）+ `ApproxLineTrsf`；引擎体全走 trait——**BRepApprox 绑定只需实现这两个 trait，引擎零改动**。Alias 表结论：BRepApprox_TheMultiLineOfApprox 与 GeomInt_TheMultiLineOfWLApprox 是同一 ApproxInt_MultiLine.gxx，仅 TheLine 不同（BRepApprox_ApproxLine vs IntPatch_WLine）。
>   3. `a0867176` **IntConicConic 全部 13 个 dispatch 臂激活 + BRepApprox 数据/工具层（5 个并行子代理产物集成）**：Circ×Circ 闭式（int_conic_conic_circ_circ.rs 987 行：IdentCircles 修正、IndirectCircles 域反转、位级 next_after）+ Lin×Ells 闭式（int_conic_conic_lin_ells.rs 1053 行：规范系变换、Extrema_ExtElC2d、Esup=DL 字面怪癖保留）+ 10 个委托型重载（int_conic_conic.rs，各按 .cxx 的侧分配/SetAccuracy/闭域处理）+ PConic::new_parabola 焦点修正（prm1=Focal()=focal_param/2，此前死代码）+ **brep_approx.rs**（1915 行：SvSurfaces trait=ApproxInt_SvSurfaces.hxx、ApproxLine、TheMultiLineOfApprox=ApproxInt_MultiLine.gxx 全 13 成员含 MakeMLBetween 弧长重采样与 CodeErreur=1 重复行、MultiLineTool 16 静态、SurfaceTool over BRepAdaptorSurface；interval/D3/DN/Bezier-BSpline/basis-curve 为具名 unimplemented!() 延迟）。
> - **并行子代理模式（实证有效，沿用）**：单 session 吞吐 ≈3 倍。做法：主代理先做文件切分（每代理独占 1-2 个文件，mod.rs 注册由主代理完成）+ 发自包含 brief（OCCT 绝对路径 + 目标文件 + landed 机制清单 + 1:1 纪律 + 锚点要求 + 验证命令 `cargo test -p rcad-algo --lib <module>` 全绿才返回）；代理并行期间主代理不得触碰其文件，全部交卷后统一集成编译 → 全量回归 → 集成提交。**session 8 实证补充**：见头部新快照的 (a)(b)(c)。
> - **session 6 完成清单**：
>   1. `2f1ac08e` **IntConicConic::Perform(gp_Lin2d, gp_Circ2d) 1:1**（_1.cxx L2236-2652）+ **IntConicConic_Tool 机制**（Tool.hxx/.cxx：Interval/PeriodicInterval 含 FirstIntersection/SecondIntersection/Normalize/Complement、Determine_Transition_LC 的 TOUCH 曲率比较分支、NormalizeOnCircleDomain）+ LineCircleGeometricIntersection（Circle± 带、IsDirect 翻转、KHROMOV 2000 双解拆分）+ ProjectOnLAndIntersectWithLDomain + gp 助手（Coefficients/Angle/CircleD1/D2/LineD1）。CInter 的 Line×Circle 与 Circle×Line 两个 dispatch 臂激活。锚点 2 个（圆心在线上的双穿越退化为 2 点；切线 y=1 单点 TOUCH）。**IntCurveCurveGen dispatch 剩余缺口降为 11 个 conic×conic 重载（Lin×Ells/Circ×Circ/Lin×Prb/Lin×Hypr/Ells×Prb/Ells×Hypr/Prb×Prb/Circ×Prb/Circ×Hypr/Ells×Ells 已走 imp-par 委托、Prb×Hypr/Hypr×Hypr）**——注意 E-E 用的 imp-par 委托先例，其余按 _1.cxx 闭式补。
> - **Stage 3c 勘察结论（下一 session 直接开工）**：
>   - BRepApprox_Approx_0.cxx include 的是 **ApproxInt_Approx.gxx**（rcad approx_int.rs 已按 GeomInt WLApprox 实例落地过同款 gxx）→ 3c-2 = approx_int 引擎泛型化 + BRepApprox 工具绑定，先例同 2a 的 int_curve_curve_gen。
>   - BRepApprox_TheComputeLineOfApprox = **AppDef_BSplineCompute** = **Approx_BSplComputeLine.gxx**（1458 行）实例化 → 3c-1 需先译该 gxx（泛型 over MultiLine/LineTool/LSQ 门面），连带 **AppDef_MyBSplGradientOfBSplineCompute / MyGradientbisOfBSplineCompute**（hxx-only ~136 行各，包 BSpParLeastSquare+BSpParFunction+BSpGradient_BFGS / ParLeastSquare+ParFunction+Gradient_BFGS 链）。
>   - app_par_curves.rs 已有 LeastSquare（含 new_bsp BSP 构造族）/ResolConstraint/ParFunction/Gradient/GradientBfgs，但 **AppParCurves_BSpParFunction.gxx 未落地**（MyBSplGradient 需要）——3c-1 的一部分。
>   - 量级：3c-1 ≈ 1458(gxx) + 2×136(hxx) + ~600(BSpParFunction) ≈ 2400 OCCT 行，独占一个 session；3c-2（ApproxInt 泛型化 + BRepApprox_Approx/ApproxLine + SurfaceTool + PrmPrmSvSurfaces/Int2S 函数族）另占一个 session。
> - **HLR 消费链现状**：`hlr/brep/` 十三件；Contap 可经 BRepTopolTool 的 ContapDomain 直接消费真实面域；CInter/Intersector 可对投影边跑 2D 相交（直线×直线/圆已闭式激活）、Intersector 可跑线/面 CS。
> - **已知遗留（非 3b 验收项）**：IntConicConic 其余 11 个 conic×conic 重载（见上）；Contap quadric-exact 死分支仍待 Adaptor3d_CurveOnSurface（3d）；EdgeData::Set 的 TopoDS→CurveView 桥接在 3f Data 落地。
> - 已确立的翻译模式（沿例勿改）：OCCT gxx 模板参数 → Rust trait；具体实例化 = 类型别名（GInter/CInter 先例）+ marker 类型（`PhantomData<&'a ()>` 生命周期弹性）；工具静态组绑进 landed 泛型引擎（as_any = DownCast 先例）；gxx/lxx 内联翻译；每函数 `// OCCT <文件> L<起>-<止>` 标注；每个叶子翻译配 OCCT 解析锚点单测。**HLR 方法论特例**：HLRBRep_Surface 持 `my_proj: *const Projector` 裸指针、Intersector 持 `my_surface: Option<*const Surface>`（OCCT 同形，Data 拥有点对象）；`HLRBRep_Curve::d3_2d`/`dn` 是 OCCT 空函数体的照抄。
> - 踩坑累积：**(6)** `&Arc<dyn Trait>` → `&dyn Trait` 不自动 deref 强转——调用点必须 `.as_ref()`；**(7)** heredoc 大文件追加第 6 次截断——一律 Write 落盘 + `cat >>` 拼接，截断后 `sed -i 'N,Md'` 清尾再重拼；**(8)** `git checkout <file>` 恢复前先确认目标状态；**(9)** 1:1 审查点：OCCT 字面 no-op（如 Lin-Circ 的 `if (Cinf >= Csup)` 同值自赋）与未初始化 out 参数（Lin-Lin 的 Pos1a/2a）必须照抄/以中性默认初始化，不许"顺手修正"。
> - 子模块工作树遗留修改（非本任务）：`algo_ext/mod.rs`、`fillet/topopebrepbuild.rs`、`tests/tkgeom_algo_gtests.rs`、`topalgo/brep_class/face_explorer.rs`——提交时只 `git add` 本任务文件。
> - shell 的 cwd 会被重置到 `C:\Users\lilu\works\rcad-pro`（根仓库），编译/测试前先 `cd rcad`；回归必含 builder 分阶段（builder_stage_tests + pavefiller_stage_tests）。

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
| IntCurve_IntConicConic | ◐ E-E（imp-par 先例）+ Lin-Lin（_1.cxx 闭式） | 缺 12 个重载（HLR 圆/圆弧边对消费前按 E-E 先例补；Stage 3b 已落地 Lin-Lin） |
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
| IntCurve_IntCurveCurveGen（1033） | HLRBRep_CInter | ✅ `geomalgo/int_curve_curve_gen.rs` 泛型引擎 `IntCurveCurveGen<C, PT, IC, IPP>`；GInter/CInter = 别名级实例化（729f0366/2f9b2d3e） |
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
- [x] **Stage 2b** Contap 包 1:1（✅ 2026-09-05：`hlr/contap/` 全包落地。叶子层：surface_adaptor（SurfaceAdapter instance-trait = `handle<Adaptor3d_Surface>` + GeomSurfaceAdapter=GeomAdaptor_Surface，含 U/VResolution 分型公式、RealFirst/RealLast 无穷域标记、face UV window）、geom_tool（GeomTool=HSurfaceTool over GeomSurfaceAdapter，喂 TopolTool）、i_type/t_function 枚举、surf_props（Normale/DerivAndNorm/NormAndDn 解析分支 + RealEpsilon apex 处理）、point/line 数据类（Arc=Arc<dyn Curve2dAdaptor>，顶点按 ParameterOnLine 排序 Add）、h_curve2d_tool（NbSamples 分型）、h_cont_tool（NbSamplesU/V/NbSamplePoints/SamplePoint + thread_local 采样窗口 + Project=EPCOfExtPC2d）、surf_function/arc_function（4 种 TFunction 分支 + IsTangent 的 d2d/d3d 推导）、cont_ana（Sphere/Cylinder/Cone × Dir/Dir+Angle/Eye 共 9 个 Perform，Tolpetit=1.e-8 照抄）。引擎层：domain（ContapDomain trait = TopolTool 引擎接口，Adaptor3d TopolTool 已 impl；3b 时 BRepTopAdaptor 补 impl）、the_search_inside（SearchInside.gxx 1:1，双 Perform）、the_search（SearchOnBoundaries.gxx 1:1：FindVertex/BoundedArc（含 Rejection 预检、SolInfo 排序、OCC569 节点扫描、BrentMinimum 精化、PointProcess、TreatLC/IsRegularity 的 edge() 早退、ComputeBoundsfromInfinite、MinFunction）——**Contap 实例的 quadric-exact IntCS 分支为死代码**（ArcFunction.myQuad 永不赋值 → TypeQuad 恒 OtherSurface，照 OCCT 字面守卫保留结构，分支体 3b 随 HLRBRep 落地）、contour（2389 行拆 mod/functions/perform/perform_ana：Recadre、LineConstructor 三分型、ComputeTangency 复杂 transition 状态机（TopTrans_CurveTransition）、ComputeInternalPoints(OnRstr) 求符号变化+二分、KeepInsidePoints、ProcessSegments、PutPointsOnLine、ComputeTransitionOngp{Line,Circle}；Tolpetit=1.e-10 vs ContAna 1.e-8 冲突照抄）。连带泛型化：**Extrema_EPCOfExtPC2d**（`geomalgo/extrema_gen_ext_pc2d.rs`=GGenExtPC+GFuncExtPC 2D 实例化，含奇异点 DN-Taylor/三点导数分支与 SearchOfTolerance；修正 PPc 方向为 pc−p）、imp_prm `FunctionSetRoot` 抽 `FunctionSetWithDerivatives2` trait（value/values/tolerance/get_state_number）、`i_walking.rs` 抽 `IWFunction` trait（+derivatives/root/is_tangent/direction_2d/point/set_surface/surface_value）——IntPatch 调用点零破坏。锚点 24 个全绿：EPCOfExtPC2d 3（圆 min/max、直线投影、区间限制）、叶子 17、Contour 集成 1（圆柱侧影端到端：两 silhouette 线位置/方向/顶点、+Y 线 Out、−Y 线 Undecided=ComputeTransitionOnLine `det<RealEpsilon()` 的 jag 940620 怪癖字面行为）+ 叶子内含。基线 algo lib **174+1i**）
- [◐] **Stage 3a** HLRBRep adaptor/tool 层（**2026-09-05 session 3 前半 ✅**：`55549b30` 三修复 + `ccde4456` Surface/Curve 主体。① 顺手项全清：kernel `Curve2d` 枚举级 is_closed/is_periodic 按 variant 委托（原来恒 false）；BRepBuilder add_pcurve/update_edge_pcurve/update_edge_pcurve_closed/update_edge_curve3d/update_edge_tolerance/set_edge_degenerated_with_clear 从 `edge_mut`(Arc::make_mut COW) 改 `edge_mut_inplace`——COW 克隆分裂 TShape 身份，wire 引用的边看不到后补的 pcurve，这是 2a-4 边界 On ignore 的第二层根因；FaceShapeSource 持 DS 约定 locations 表（槽 0=identity），kernel-BRep 调用方前置 identity——`fclass2d_topol_perform_boundary_on` 已 un-ignore，**algo lib 179/179 零 ignore**。② `hlr/brep/` 落地：BRepAdaptorSurface（Initialize(F,R=true) 受限形，窗口=2a-4 face UV domain，GeomAdaptor 角色委托 contap::GeomSurfaceAdapter；poles_grid()）、BSurfaceTool NbSamplesU/V（负分母 range 变体照抄）、BCurveTool=CurveView trait、HLRBRep_Surface（load 分类 degree1-Bezier→Plane、SideRowsOfPoles+GProp_PEquation Z 平面臂、IsSide 四分支、IsAbove 30 点测试、Plane() Bezier 拟合、lxx 一行组、实现 SurfaceAdapter 供 Contap 直连）、HLRBRep_Curve（Parameter2d/3d 透视公式、Update 投影分类（圆→椭圆 AngleWithRef 偏移、deg1-Bezier→线、透视线系数 VX/OX/OF/VZ/OZ）、UpdateMinMax 30 点挠度、D0/D1/D2 2D 透视公式、D3/DN 空体照抄、Line/Circle/Ellipse 的 ProjLib XOY 像）。锚点 4 个。**3a 后半 ✅（`8eb52a47`）**：kernel 新增 `geom_lprop/cl_props_base.rs`（GeomLProp_CLPropsBase + LProp_CurveUtils 泛型引擎：EnsureDeriv/IsTangentDefined 显著导数搜索/ComputeTangent 弦号校正/Curvature/Normal/CentreOfCurvature，Accessor=ToolAccess trait；LPropStatus 重排为 OCCT 判别序 Undecided/Undefined/Defined/Computed，rcad-only Zero 追加尾部保旧实现）；`hlr/brep/cl_props.rs` HLRBRep_CLProps=ClPropsBase over Curve（CLPropsATool 委托，D3 空体照抄）+ `Curve::tangent()` 已接通（Epsilon(1.) 容差）；`hlr/brep/line_tool.rs` LineTool gp_Lin 静态组。锚点 +3。**3a 全部关闭**。3b 时补：SLPropsATool/HLRBRep_SLProps（消费方 InterCSurf，需 LProp_SurfaceUtils）+ BRepTopAdaptor_TopolTool 的 ContapDomain impl（需 Arc<BRep> 所有权设计））
- [◐] **Stage 3b** HLRBRep 相交层（**2026-09-05 session 4 首项 ✅ `e226a8f9`**：kernel `geom_lprop/sl_props_base.rs` = GeomLProp_SLPropsBase + GeomLProp_SurfaceUtils 泛型引擎（FindSurfTangentOrder/ComputeSurfTangent U-V 弦号校正/ComputeSurfNormal→已落地的 CSLib::Normal/ComputeSurfCurvatures 第一二基本形式 + math_DirectPolynomialRoots 二次分支 + 脐点捷径；SLPropsSurface trait = ToolAccess 策略）；`hlr/brep/sl_props.rs` HLRBRep_SLProps = SlPropsBase over HLRBRep_Surface（SLPropsATool 委托 + Bounds=First/Last 参数）。锚点：kernel 圆柱主曲率 0/-1 + 球脐点 K=1；HLR 圆柱经 Surface 适配器（kernel trait-default D2 为 h=1e-5 有限差分，容差 1e-6）。**3b 后续 ✅（`d92659f9`）**：`hlr/brep/surface_tool.rs` = HLRBRep_LineTool（HCurveTool+CurveTool3d over gp_Lin；无限域/CN/ElCLib/空 Intervals 体/零高阶导/NbSamples=3）+ HLRBRep_SurfaceTool（HSurfaceTool+PSurfaceTool over HLRBRep_Surface；全成员一行转发；UTrim/VTrim=窗口收窄；BasisCurve/BasisSurface/Direction/OffsetValue/Bezier/BSpline=NoSuchObject 臂）；`hlr/brep/inter_csurf.rs` HLRBRep_InterCSurf = IntCurveSurface Inter.pxx 组装绑定固定类型（TheCurve=gp_Lin/TheCurveTool=LineTool/TheSurface=HLRBRep_Surface/TheSurfaceTool=SurfaceTool），Perform/PerformBounds/Polygon/Polyhedron/BSB/InternalPerform x4 与 OCCT cxx L107-378 成员一一对应；BRepAdaptorSurface/Surface 补 u_trim/v_trim。锚点：线/平面 w=1.5 UV(0,0.5)、线/单位圆柱 w=1/3 两次穿越。**3b 收尾 ✅（2026-09-05 session 5：`fe89ac91`+`729f0366`+`2f9b2d3e`）**：① ContapDomain for BRepTopAdaptor_TopolTool——`Curve2dAdaptor::as_any`（occ::down_cast）+ `HVertexBehavior`/`HVertexHandle`（handle<Adaptor3d_HVertex> 多态句柄：Value/Parameter/Resolution/Orientation/IsSame 全虚）+ BRep 侧四类（BRepCurve2d/BRepHVertex/FClass2dTopol/BRepTopolTool）转 **Arc\<BRep\> 所有权** + ContapDomain::edge()→Option\<Shape\>；锚点 = trait 驱动 BRep 域全迭代（arc 下行转换/顶点迭代/Identical 拓扑 IsSame/Edge/Classify/采样格）。② IntCurveCurveGen 泛型化——`IntCurveCurveGen<C, PT: PCurveTool<C>, IC: IntConicCurveMember<C>, IPP: IntPCurvePCurveMember<C>>` 引擎 + GInter 别名级实例化（inter_cc.rs 持 GInter）+ gxx L1023 SetMinNbSamples；连带 **Extrema 定位机制泛型化**（LocatorCurveTool：locate_on_curve/GFuncExtPC/gen_locate_ext_pc/proj_cur_find_parameter_{bounded,unbounded} 共享体，Geom2d 侧签名不变）。③ CInter 家族——`hlr/brep/curve_tool.rs`（HLRBRep_CurveTool hxx/cxx/lxx 全量：NbSamples 双重载 Line=2/Bezier=3+NbPoles/BSpline 钳 50、EpsX=1e-10、单 C1 区间约定）+ `c_inter.rs`（TheProjPCurOfCInter/TheIntersectorOfCInter/TheIntConicCurveOfCInter/IntConicCurveOfCInter/CPolyTool+TheIntPCurvePCurveOfCInter/`CInter` 别名）+ **IntConicConic::Perform(Lin,Lin) 1:1**（_1.cxx L1381-2233 + DomainIntersection/LineLineGeometricIntersection/FindPositionLL/getDomainParametrs/computeIntPoint(gka 0022833)/CheckLLCoincidence/SegmentToPoint 七助手 + 文件局部 TOLERANCE_ANGULAIRE=1e-15 照抄）。④ HLRBRep_EdgeData（EMaskFlags 位旗标全套 + Set 的内核边界拆分）+ **HLRBRep_Intersector**（hxx L26-98 + cxx L86-646：双边 perform 的 decalage walk-away 循环（aDecalagea1/b1/a2/b2 四系数 ×2 加倍与 -1 失效）、SimulateOnePoint（determine_transition_in_out）、Load/perform_line CS 分支（quadric 直连 / polyhedron 八角参数窗 pmin/pmax 剔除 + ThePolygonOfInterCSurf(L,pmin,pmax,3)；aMinNbHLRSamples=4 注入 CInter）。锚点 6 个：CInter Line×Line 中点交叉 + Other×Other 多边形精化链（经 CPolyTool/ProjPCur/ExactIntersectionPoint）、CurveTool 投影线、Intersector 双边交叉/SimulateOnePoint/线-平面 CS w=1.5。**3b 已知遗留**：IntConicConic 其余 12 个 conic×conic 重载 unimplemented（HLR 圆/圆弧边对消费前按 E-E imp-par 先例补）；Contap quadric-exact 死分支待 3d；EdgeData::Set 的 TopoDS→CurveView 桥在 3f。基线：algo lib **192** / kernel 648 / builder 分阶段 76+1 / pavefiller 26）
- [x] **Stage 3b** HLRBRep 相交层 + The* 实例化（2026-09-05 session 5 全关）
- [x] **Stage 3c** BRepApprox_Approx 链路（**3c-1 ✅ `d1dc1ad1`**：AppParCurves_BSpFunction/BSpGradient.gxx → app_par_curves_bsp.rs + Approx_BSplComputeLine.gxx → bspl_compute_line.rs，AppDef_BSplineCompute 实例化，OCCT ComputeCurve 死声明省略。锚点 2 个：四分之一圆 17 点逼近达容差、两点 Interpol 度 1。**3c-2 ✅ 2026-09-05 session 8**：① 引擎接缝——`impl ApproxIntWLine for ApproxLine` + `impl ApproxIntMultiLine for TheMultiLineOfApprox`（brep_approx/ 目录化：mod.rs 1912 + seams.rs 425，锚点 4 个）；② 函数族——`brep_approx_prm_prm.rs`（ZerParFunc/Int2S/PrmPrmSvSurfaces + choix_ref/compute_tangence，锚点 6）+ `brep_approx_imp_prm.rs`（ZerImpFunc/ImpPrmSvSurfaces 含 Singular/NonSingular/FillInitialVectorOfSolution 全链，锚点 3）+ `TheMultiLineOfApprox<'a>` 生命周期槽（OCCT void* 语义）；③ 本体——`brep_approx_approx.rs`（perform_wline/perform_prm_prm quadric switch/perform_implicit/SetParameters/accessors，747 行，锚点 4）+ 引擎 `perform_impl(cut)` 零行为重构。怪癖照抄：P2DOnFirst=ApproxU1V1、MakeMLBetween 空线失败语义、gxx L611 SetCoord(tu2,tu2) 笔误、quadric switch default 臂 Quad 默认构造。基线 algo lib **240**）
- [x] **Stage 3d** HLRTopoBRep 包（✅ 2026-09-05 session 9：`geomalgo/hatch/hatcher.rs`（Hatcher 主引擎 2385 行总量 1:1，锚点 5）+ `topalgo/brep_top_adaptor/tool.rs`（BRepTopAdaptor_Tool）+ `hlr/topo_brep/`（Data/FaceData/VData/DSFiller/FaceIsoLiner/OutLiner，两波五代理；锚点 24：Data 8 + DSFiller 8 + FaceIsoLiner 5 + OutLiner 3；kernel topods.rs +3 BRepBuilder 方法）。全链 OutLiner::Fill → DSFiller::Insert → Data 就绪。基线 algo lib 299 / kernel 649）
- [x] **Stage 3e** HLRBRep 干涉数据结构（✅ 2026-09-05 session 8 并行完成：`hlr/brep/` 八件——BiPoint/BiPnt2D（hxx-only）/AreaLimit/EdgeIList/VertexList/FaceIterator/EdgeInterferenceTool/EdgeBuilder（521 全 1:1）；锚点 11 个全绿。遗留：EdgeInterferenceTool 的 `Data` trait 待 3f HLRBRep_Data 实现）
- [x] **Stage 3f** HLRBRep Data（拆子模块）（✅ 2026-09-05 session 10：`hlr/brep/data/` = mod.rs（Data\<'a\> 70 字段契约 + 文件静态常量/counters）+ update.rs（ctor/Write/Update/InitBoundSort/InitEdge/exploration + lxx，8 锚点）+ classify.rs（2416 行：NextInterference…IsBadFace 16 方法 + 文件静态，7 锚点）+ tableau_rejection.rs（cxx L77-492，5 锚点）；连带 face_data.rs（FaceData 方法+lxx 13 对位访问器，3 锚点）与 edge_face_tool.rs（UVPoint/CurvatureValue + ExtPF 本地翻译，6 锚点）；uv_point 7 参契约裁决（Data.my_brep kernel 上下文 + my_e_map 槽位）；hlr:: 树 125/125 绿）
- [x] **Stage 3g** Hider→InternalAlgo→Algo→HLRToShape（✅ session 10 六件全落 + **session 11 烟囱闭环**：`hlr/tests.rs` 4 个端到端 smoke——**盒体与 OCCT 完全一致**（VCompound 9 边 146.96939 = OCCT lprops、HCompound 3 边 48.98980、12 条投影线段逐一解析对拍 1e-6、可见/隐藏归属一致）+ 圆柱部分锚定（φ=45° 轮廓/椭圆弧分段/seam 像）+ 稳定性 3 连跑。锚点 26+4+2 转绿 = hlr:: 151/0/0）
- [x] **Stage 4a** gen_hlr_ref.py（✅ session 11：4/6 JSON + tools/occt-hlr-runner ref_output.txt）
- [◐] **Stage 4b** hlr/commands.rs（无头验收走 occt-test-gen hlr_translate 的 HLR_HELPERS + rcad smoke 路径，HLRTest Tcl 层未单独移植；如后续需要交互式 DRAW 等价再单列）
- [x] **Stage 4c** occt-test-gen hlr 接入（✅ session 12：tools/occt-test-gen/src/hlr_translate.rs（try_translate_hlr_script + generate_hlr_batch + HLR_HELPERS 镜像 ViewerTest VComputeHLR L3299-3362 组成）+ main.rs 4 接入点（mod 声明 / batch 分支 exact_hlr|poly_hlr / translate_draw_script_inner 分发 / 生成文件 generated_occt_boolean_hlr_{grid}.rs）；80 例 locate_data_file 永久排除、Plate 因无 BRepOffsetAPI_MakeFilling 跳过、poly_hlr 88 例按 Stage 5 跳过）
- [◐] **Stage 4d** exact_hlr 断言闭环（session 12：hlr/acceptance.rs 三测试——box 全绿（消费 ref_output.txt，9/3 边 + 三 mass 全中）；bug25813_1 与 ptorus 两测试就位带完整 OCCT 真值分解注释、#[ignore] 待最后两缺口（见文末 session 12 追加段的机制级诊断）。阶段 2 已再修 6 处形式偏差，combined mass 266.5264 与 OCCT 精确一致）
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

### 本 session 追加 3（2026-09-05，session 4：3a 完成 + 3b InterCSurf）

| commit | 内容 |
|---|---|
| `55549b30` | 3a 修复三件：Curve2d 枚举委托 / UpdateEdge 原地身份（COW 分裂根因）/ locations 表对齐（un-ignore 边界 On） |
| `ccde4456` | 3a part 1：BRepAdaptorSurface + HLRBRep_Surface/Curve + B{Surface,Curve}Tool（1840 行） |
| `8eb52a47` | 3a part 2：ClPropsBase 泛型引擎 + LPropStatus OCCT 判别序 + HLRBRep_CLProps + Curve::tangent + LineTool |
| `e226a8f9` | 3b part 1：SlPropsBase 泛型引擎（基本形式+DirectPolynomialRoots+脐点捷径）+ HLRBRep_SLProps |
| `d92659f9` | 3b part 2：HLRBRep_LineTool/SurfaceTool 工具绑定 + HLRBRep_InterCSurf 实例化（线/平面、线/圆柱锚点） |
| `4833acd8` | §7 文档（3a 全关 + 3b 两项落地） |

**本 session 新确立的翻译事实（沿例勿改）：**

1. **InterCSurf 的 TheCurve 是 gp_Lin**（InterCSurf.hxx L44-53 别名表）——HLR 的线/面相交输入是视线/轮廓直线，不是任意曲线适配器。因此该实例化走「landed HInter 引擎直通 + 固定类型参数」零新引擎代码；与 CInter 不同（TheCurve=HLRBRep_CurvePtr 必须泛型化 IntCurveCurveGen 引擎本体）。
2. **HLRBRep_LineTool 的怪癖照抄清单**：Intervals 空函数体；D2/D3 = ElCLib::D1 + 零向量；DN 仅 N=1 返回方向；Circle/Ellipse/Hypr/Parab 返回默认构造（rcad 用 panic 标记不可达）；Bezier/BSpline 返回 null handle；**NbSamples=3**（非 Trait 默认 10）。
3. **HLRBRep_SurfaceTool 全成员一行转发**到 HLRBRep_Surface 适配器；UTrim/VTrim = BRepAdaptor_Surface::UTrim（窗口收窄同一曲面，GeomAdaptor 语义）；BasisCurve/BasisSurface/Direction/OffsetValue/Bezier/BSpline = Standard_NoSuchObject 路径（BRepAdaptor_Surface 对非偏移/非拉伸面的行为）；NbSamplesU/V 分派与 Adaptor3d_HSurfaceTool 相同（用 trait 默认即可）。
4. **HLRBRep_Surface 的裸指针**：`my_proj: *const Projector` 按 OCCT 同形保留（`unsafe impl Send`；Projector 由 Data 拥有、Surface 只读引用——与 OCCT HLRBRep_Data 对象图相同）。
5. **HLRBRep_Curve 的空函数体**：`d3_2d`（cxx L406）与 `dn`（cxx L410-413 返回 gp_Vec2d()）是 OCCT 字面行为；CLProps 的 D3 评估经 ATool 调到空体后 deriv[2] 保持零——切线搜索自然落到更高阶失败，与 OCCT 一致。
6. **坑 (4)**：kernel 泛型文件 `use glam::Vec3` 会引入 f32 类型——必须 `use crate::geom::Vec3`（=DVec3 别名）。
7. **坑 (5)**：枚举判别序被隐式数值比较依赖时（LPropStatus `>= Defined`、`== Undefined`），重排为 OCCT 序后必须全库 grep 使用点；旧实现的 rcad-only 状态（Zero）追加尾部而非删除，保旧语义零回归。
8. **BRepBuilder COW 陷阱（坑 (6)，2a-4 边界 On 的第二层根因）**：`Arc::make_mut` 式 `edge_mut` 在其它句柄（wire）已引用 TShape 时克隆分身——BRepBuilder 的 UpdateEdge 系方法（补 pcurve/3D 曲线/公差/退化标记）必须用 `edge_mut_inplace`；语义 = OCCT BRep_Builder::UpdateEdge 原地改 TShape。

### 下一 session 入口（按序）

1. **Stage 3c-2**：approx_int.rs（ApproxInt_Approx.gxx 的 GeomInt 实例）泛型化 + BRepApprox_Approx/ApproxLine + BRepApprox_SurfaceTool + PrmPrmSvSurfaces/Int2S 函数族绑定（BRepApprox_TheComputeLineOfApprox = 已落地的 BSplineCompute 换 BRepApprox MultiLine/MyLineTool 绑定）。
2. 顺手项（HLR 圆/圆弧边对消费前必须补）：IntConicConic 其余 11 个 conic×conic 重载——Lin×Ells（_1.cxx L2861）/Circ×Circ（L807）照 _1.cxx 闭式移植（Lin×Circ 已落地 `2f1ac08e`），其余按 E-E imp-par 先例评估委托。
3. 之后 3d（HLRTopoBRep：DSFiller/Data/FaceIsoLiner/OutLiner + Geom2dHatch_Hatcher）→ 3e（干涉数据结构 + EdgeData::Set 的 TopoDS→CurveView 桥）→ 3f（Data 2683 拆子模块，持 Arc\<BRep\>）→ 3g（Hider→InternalAlgo→Algo→HLRToShape + 烟囱测试）。

### 本 session 追加 4（2026-09-05，session 5：3b 全关 + Intersector）

| commit | 内容 |
|---|---|
| `fe89ac91` | **ContapDomain for BRepTopAdaptor_TopolTool**：Curve2dAdaptor::as_any（DownCast）+ HVertexBehavior/HVertexHandle 多态顶点句柄 + BRep 侧四类转 Arc\<BRep\> 所有权 + ContapDomain::edge()→Option\<Shape\>（13 文件，contap 链零回归） |
| `729f0366` | **IntCurveCurveGen 泛型化**：`IntCurveCurveGen<C, PT, IC, IPP>` 引擎 + IntConicCurveMember/IntPCurvePCurveMember 成员接口 + SetMinNbSamples（gxx L1023）+ GInter 类型别名；Extrema 定位机制泛型化（LocatorCurveTool + 共享 FindParameter 体） |
| `2f9b2d3e` | **CInter 家族 + Intersector**：HLRBRep_CurveTool / TheProjPCurOfCInter / TheIntersectorOfCInter / TheIntConicCurveOfCInter / IntConicCurveOfCInter / CPolyTool+TheIntPCurvePCurveOfCInter / CInter + HLRBRep_EdgeData + HLRBRep_Intersector（99/858）+ **IntConicConic::Perform(Lin,Lin) 1:1**（_1.cxx L1381-2233 + 七助手）；锚点 6 个 |

**本 session 新确立的翻译事实（沿例勿改）：**

1. **DownCast 的 rcad 形**：OCCT `occ::down_cast<T>(handle)` → `Curve2dAdaptor::as_any() -> &dyn Any`（必需方法，6 个实现者各一行）；BRepTopolTool::Initialize(C)/Orientation(C) 经 `a.as_any().downcast_ref::<BRepCurve2d>()` 照抄，失败 panic = Standard_ConstructionError。
2. **多态顶点句柄**：`handle<Adaptor3d_HVertex>` → `HVertexHandle = Arc<dyn HVertexBehavior>`（value/parameter/resolution/orientation/is_same/topo_vertex）；BRepTopAdaptor_HVertex::IsSame 的静态下行转换 → `other.topo_vertex()` None 臂（混合类型单域内不出现，返回 false 注释说明）；base Identical = `v1.is_same(v2.as_ref())` 多态分发，两路语义自然一致。**坑 (6)**：`&Arc<dyn T>` → `&dyn T` 不自动 deref 强转（E0277），调用点必须 `.as_ref()`。
3. **BRep 侧所有权**：BRepCurve2d/BRepHVertex/FClass2dTopol/BRepTopolTool 全部从 `&'a BRep` 转 `Arc<BRep>`（HLRBRep_Data 3f 的所有权设计）——借用版删除不保留；BRepHVertex 持 owned curve（OCCT ctor 的 handle 拷贝语义）；ContapDomain::arc() 返回 `Arc::new(curve.clone())`（廉价：Arc clone + 小字段）。
4. **实例形别名先例扩展**：GInter/CInter 均为 `pub type X = IntCurveCurveGen<...>` 类型别名级实例化；子引擎成员经 `IntConicCurveMember`/`IntPCurvePCurveMember` trait（base()/base_mut() + perform 系）注入——IC 需要 `Default`（TheIntConicCurveOfCInter 的 Default = bare()）；成员是泛型类型时 Default 按 `(C, PT, JT)` 具体实例 impl。
5. **定位机制单源**：Extrema_GCurveLocator/GFuncExtPC/GenLocateExtPC/TheProjPCur::FindParameter 是 Geom2dInt 与 HLRBRep 两份逐字相同的实例化体 → 泛型 `LocatorCurveTool<C>`（first/last/value/d0/d1/d2/dn/get_type/nb_samples/eps_x）+ `proj_cur_find_parameter_{bounded,unbounded}` 共享体；Geom2d 侧委托保持原签名（`'a` 生命周期弹性）。**坑**：`impl LocatorCurveTool<dyn Curve2dAdaptor>` 不带 'a 会在委托点报 "requires that 'a must outlive 'static"——impl 必须写 `impl<'a> LocatorCurveTool<dyn Curve2dAdaptor + 'a>`。
6. **IntConicConic::Perform(Lin,Lin)**（_1.cxx L1381-2233）：文件局部 `#define TOLERANCE_ANGULAIRE 1.e-15`（≠ IntImpParGen 的 1e-8，NIZHNY-MKK 2000 修改）照抄为模块私有常量；OCCT 未初始化 out 参数（Pos1a/Pos2a/Pos1b/Pos2b）在 rcad 以 Position::Middle 初始化（domain_intersection 必写、早退分支不读，语义不变）；FindPositionLL 的 `double& Param` 归一 inout → `&mut f64`；computeIntPoint/CheckLLCoincidence/SegmentToPoint/DomainIntersection/getDomainParametrs 为文件静态自由函数。
7. **Intersector 的借用形状**：OCCT `HLRBRep_EdgeData*` 可变指针在 rcad 全部走 immutable（`geometry()` + `status_ref()`）——引擎只读边数据；`my_surface: Option<*const Surface<'a>>` 保留 OCCT `HLRBRep_Surface*` 裸指针对（load 后解引用 unsafe 限一处）。HLRBRep_CurveTool::NbSamples(C,U1,U2) 的 BSpline 臂用单 C0 区间约定（2*Degree，钳 [2,50]），C1 区间机制落地时同步。
8. **坑 (7)**：heredoc 大文件追加第 6 次截断（本次截在 compute_int_point 中段）——一律 Write 落盘 + `cat >>` 拼接；截断后先 grep 定位节标记行号、`sed -i 'N,Md'` 清尾再重拼，不要盲目 `git checkout`（会连带回滚已成功的 Edit）。
9. **回归必含 builder 分阶段**：`cargo test -p rcad-algo --test builder_stage_tests`（76）+ `builder_stage_smoke`（1）+ `pavefiller_stage_tests`（26）——builder 分阶段测试基线与 TKHLR 锚点一并跑。

### 下一 session 入口（2026-09-05 session 1 留，已被 session 2 取代，见文末）


### 本 session 追加 2（2026-09-05，Stage 2b Contap 包 1:1 完成）

| commit | 内容 |
|---|---|
| （本轮代码提交） | **2b 全包**：`hlr/contap/`（18 文件）+ `geomalgo/extrema_gen_ext_pc2d.rs` + `FunctionSetRoot`/`IWalking` 泛型化。锚点 24 个全绿；基线 algo lib 174+1i / kernel 645 / tkhelix 16 / pavefiller 26 / tkgeom_algo 134+1 |

**本 session 新确立的翻译事实（沿例勿改）：**

1. **实例形 adaptor trait（`SurfaceAdapter`）**：OCCT 把 `handle<Adaptor3d_Surface>` 在 Contap 里到处传，rcad 对应为 object-safe 的 instance trait（方法带 `&self`），`SurfaceHandle = Arc<dyn SurfaceAdapter>` 就是 handle。这与 `int_curve_surface::HSurfaceTool`（static 工具形，带关联类型、非 object-safe）是两个面：Contap 引擎直接调 instance trait（每个 OCCT `Adaptor3d_HSurfaceTool::X(S,...)` = `s.x(...)` 一次 delegation），`GeomTool`（= HSurfaceTool over GeomSurfaceAdapter）只用于喂 Adaptor3d TopolTool。后续 HLRBRep_Surface 实现 SurfaceAdapter 即接入 Contap。
2. **无内核曲面的引擎实例化桥**：rcad IntWalk 引擎是 Surface3 实例（perform 收 `caro: &Surface3` + domain 窗口）。Contap_TheIWalking 的组装 = `IWFunction` trait impl 中的 `set_surface(&Surface3)` 忽略入参、对**已存 adapter** 重跑 `Set`（OCCT `Func.Set(Caro)` 在 Contour L152 已做过、gxx L247 重做是幂等的——语义等价），`surface_value` 委托 adapter；Contour 调 engine 时经 `kernel_surface()` 取底层 Surface3。**偏差已注释**：engine 的 tolerance 用 `u_resolution(domain)` 近似，OCCT 是 `UResolution(Caro, Confusion())`——只影响步进容差。
3. **死分支照 OCCT 字面守卫保留**：`IntStart_SearchOnBoundaries.gxx` 的 quadric-exact IntCS 分支以 `TypeQuad != OtherSurface` 为门——Contap 实例的 `ArcFunction.myQuad` 从不赋值（`IntSurf_Quadric()` 默认构造即 OtherSurface，cxx L38-46），故该分支在 OCCT 里就是死代码。rcad 保留守卫与结构、分支体记 documented-unimplemented（3b 随 Adaptor3d_CurveOnSurface 落地）。`TreatLC`/`IsRegularity` 同理：base TopolTool 的 `edge()` 为 None → 字面早退。
4. **jag 940620 怪癖（锚点测试记录）**：`ComputeTransitionOnLine` 的 `det < RealEpsilon() → Undecided` 对圆柱第二条侧影线（u=3π/2 处 v1=+1、det=−1）成立——OCCT 字面即返回 Undecided，不是 bug。集成锚点断言 `[Out, Undecided]`。
5. **数据类内部可变性裁剪**：OCCT `Contap_Line` 与 IWLine 共享 `IntSurf_LineOn2S` handle——rcad 侧 IWLine 在 SetLineOn2S 后不再被读、后续变更全部经 Contap_Line，因此 owned clone 语义等价（文件头注明）；顶点序列同理用 `Vec<Point>` + `&mut` 访问，无需 RefCell。
6. **算子方向坑（第三次出现）**：OCCT `gp_Vec2d(A, B)` = B − A。`Extrema_GFuncExtPC` 的 `PPc(myP, myPc)` = `pc − p`（写反会翻转 IsMin 分类；零点集合不变，所以 LocateExtPC 类只找零点的消费者永远发现不了）。已修 `extrema_gen_ext_pc2d.rs`；geom2d_int 私有 GFuncExtPC 的同款符号问题已记头部快照（暂不动，FindParameter 只用零点）。
7. **kernel 缺口记录**：`Curve2d` 枚举级 `Curve2dEval::is_closed/is_periodic` 未按 variant 委托（trait 默认 false）——3a 顺手修；本次锚点改用 adaptor 层 `period()`。
8. **heredoc 截断第 5 次**——大段代码一律 Write 落盘 + `head`/`cat` 拼接，不在 bash heredoc 里内联大文件。

### 下一 session 入口（按序）

1. **Stage 3a**：HLRBRep adaptor/tool 层——Curve(215/632/169)/Surface(197/335/271)/CurveTool/SurfaceTool/BCurveTool/BSurfaceTool/LineTool/CLProps/SLProps；HLRBRep_Surface 实现 `SurfaceAdapter`、HLRBRep_Curve 接 `Curve2dAdaptor`，即同时接入 Contap 与 HInter。
2. 顺手项 (a)：`Curve2d` 枚举级 is_closed/is_periodic 委托修复；(b) FaceShapeSource/BRepTool locations 表约定对齐，un-ignore `fclass2d_topol_perform_boundary_on`；(c) `BRepTopAdaptor_TopolTool` 补 `ContapDomain` impl（3b 用）。
3. Stage 3b：HLRBRep 相交层 + The*OfInterCSurf 实例化（HInter 同套 trait）+ `Adaptor3d_CurveOnSurface`（补 TheSearch 的 quadric-exact 分支）。

### 本 session 追加（2026-09-05，session 8：Stage 3c-2 全关 + Geom2dHatch 叶子 + Stage 3e）

| commit | 内容 |
|---|---|
| （本轮代码提交 1） | **3c-2 全关**：brep_approx 目录化（mod.rs 数据层 + seams.rs 接缝 = `ApproxIntWLine for ApproxLine` + `ApproxIntMultiLine for TheMultiLineOfApprox`，P2DOnFirst=ApproxU1V1 字面照抄）+ `brep_approx_prm_prm.rs`（ZerParFunc/Int2S/PrmPrmSvSurfaces + choix_ref/compute_tangence，2000 行）+ `brep_approx_imp_prm.rs`（ZerImpFunc/ImpPrmSvSurfaces，1786 行，gxx L611 SetCoord(tu2,tu2) 笔误照抄）+ `brep_approx_approx.rs`（BRepApproxApprox 壳：perform_wline/quadric switch/perform_implicit，747 行）+ 引擎 `perform_impl(cut)` 零行为重构（Perform(WLine) Init(true) vs PrmPrm/ImpPrm Init(cut) 的 gxx 字面差异）；`TheMultiLineOfApprox<'a>` 生命周期槽（OCCT void* PtrOnmySvSurfaces 语义）。锚点 17 个（4 seams + 6 prm_prm + 3 imp_prm + 4 approx） |
| （本轮代码提交 2） | **并行加做**：`geomalgo/hatch/` 目录化（9 文件 2545 行：HatchGen 六类 + Elements/Intersector/FClass2dOfClassifier/Classifier/Hatching；OCC12627 magic 拷贝构造、OtherSegment 探测常数与放行语义照抄；锚点 8 个）+ `hlr/brep/` 3e 八件（~2400 行：EdgeBuilder 19 方法全 1:1、HasArea right==myLimits 约定、Orientation() WNT 警告怪 return 照抄；锚点 11 个）+ 本文档更新 |

**本 session 新确立的翻译事实（沿例勿改）：**

1. **void* 槽的生命周期形**：OCCT `void* PtrOnmySvSurfaces` → `Option<Arc<dyn SvSurfaces + 'a>>`（TheMultiLineOfApprox<'a>）——比裸指针干净且 msvc 兼容；`BRepAdaptorSurface<'a>` 的 PhantomData 签名不动，函数族（ZerParFunc/Int2S/PrmPrmSvSurfaces/ZerImpFunc/ImpPrmSvSurfaces）持 `&'a BRepAdaptorSurface` 借用。
2. **具体类型的方法解析优先命中固有方法**——经 trait 缝调用必须显式限定 `<TheMultiLineOfApprox as ApproxIntMultiLine>::make_ml_between(...)`（固有版存在同名方法时）。
3. **引擎 seam 的 value 分路径**：WLineApprox Bezier 路径结果在 `my_bez_to_bspl.value()`，非 Bezier 才在 `my_value` 字段——消费/测试一律走 `apx.value()` 访问器。
4. **`Arc<dyn Trait>` 尾随生命周期**：`Arc<dyn SvSurfaces>` 默认 'static；带借用实现者的 slot 必须写 `Arc<dyn SvSurfaces + 'a>` 并让容器类型带 'a（协变可收 'static 测试 stub）。
5. **Geom2dHatch 的 Elements 拷贝构造是 OCC12627 magic**（不复制 map、迭代器归零）——手写 Clone 非 derive；NCollection_DataMap 的 Bind 覆盖返回 false 语义 → `Vec<(usize, T)>` 保序实现。
6. **TopClass 两个 .pxx 是 C++ 模板头**（Classifier2d/FaceClassifier）——按 BRepClass 先例把模板体内联进 Geom2dHatch 具体类。
7. **EdgeBuilder 的句柄相等含 null==null=true**（OCCT 句柄比较语义）→ `Arc::ptr_eq` 前置 null 检查；析构函数 → `Drop{destroy()}`。
8. **多代理集成时序（session 8 实证）**：主代理预创建骨架文件 + 预注册 mod.rs → 上游叶子代理（A/B）先行 → 有编译依赖的壳代理（C）后发，brief 里定死跨代理类型契约 → 等待期主代理做 diff 审计 + 怪癖抽查 + 下一 Stage 勘察 + 文档草稿。
9. **heredoc 截断第 7 次**（seams.rs 测试段）——继续执行 Write 落盘 + `cat >>` 拼接纪律。

### 下一 session 入口（按序）

1. **Stage 3d**：Geom2dHatch_Hatcher 主引擎（1919 cxx + 209 hxx + 257 lxx——叶子层已全落，五入口 AddElement/AddHatching/Trim/ComputeDomains/Domain；注意 myElements 的 `change_find` 改元素方向、`myHatchings` 整数键 map）→ HLRTopoBRep 包（DSFiller 758 → Data 303 拆子模块 → FaceData/VData → FaceIsoLiner 488 → OutLiner 344；DSFiller 直接消费 BRepApprox_Approx + Contap + BRepTopAdaptor，全就绪）。
2. 之后 3f（HLRBRep Data 2683 拆子模块，持 Arc<BRep>；落 EdgeInterferenceTool 的 Data trait）→ 3g（Hider→InternalAlgo→Algo→HLRToShape + 烟囱测试）→ Stage 4 验收闭环。

### 本 session 追加（2026-09-05，session 9：Stage 3d 全关）

| commit | 内容 |
|---|---|
| （本轮代码提交） | **3d 全关**：`geomalgo/hatch/hatcher.rs`（Hatcher 主引擎 1995 行：trim 三重载/global_transition/compute_domains + lxx 全内联；cxx L934 Param 覆盖死代码与 L545 无条件 cout 照抄；锚点 5）+ `topalgo/brep_top_adaptor/tool.rs`（BRepTopAdaptorTool）+ `hlr/topo_brep/`（Data/FaceData/VData + DSFiller 1422 行 + FaceIsoLiner ~990 行 + OutLiner ~1050 行；锚点 24）+ kernel topods.rs +3 BRepBuilder 方法（add_to_edge/remove_from_edge/update_vertex_on_edge）。两波五代理（F/G → H/I/J）契约化并行，H↔I 契约一次对接成功 |

**本 session 新确立的翻译事实（沿例勿改）：**

1. **map 键 = `Shape::ptr_id()`**：`TopTools_ShapeMapHasher`（TShape 指针身份，IsSame 忽略 location）→ `Vec<(Shape, V)>` 插入序 + ptr_id 查找（data.rs 的 `map_is_bound`/`map_pos` 先例；`expect` = `Standard_NoSuchObject`）。DataMap 迭代序（OCCT 未定义）→ 确定性插入序。
2. **NCollection_List::InsertBefore 的 PInsertBefore 语义**：插入后迭代器仍指原当前节点——Vec 编码下 insert 后索引 +1（data_edge_vertex_iteration 锚定）。
3. **OCCT DataMap::Bind 是值拷贝** → MST 存的 Tool 类型需 `#[derive(Clone)]`（tool/topol_tool_brep/fclass2d_topol/class2d 四件）。
4. **DSFiller 的 Walking 臂三分支**：`ipL-ipF<1` 失败 / `<5` 度1 BSpline（knots=1..nbp、mults 首末=2）/ else BRepApprox_Approx 链（`Approx.Value(1)/Value(2)` 的 Index 参数被 gxx L487-497 忽略——rcad 调 `value()` 一次，curve(1) 取 3d、curve2d(2) 取 2d 槽）；`L447-450 if (Maxz > Maxx) { Maxx = Maxy; }` OCCT 笔误逐字照抄。
5. **FaceIsoLiner 的怪癖**：MakeEdge 将 INTERNAL 方向重置为 FORWARD（TopoDS_Builder.cxx L28-33）；`ComputeDomains(IndH)` 消费点每轮恰一个 bound hatching → 无参 `compute_domains()` 数据流等价；V1U/V2U 跨域残留未重置（OCCT 字面）。
6. **OutLiner 的怪癖**：cxx L106-123 注释掉的 splitted 预扫描 → `NF = F.EmptyCopied()` 无条件照抄；SameEdge 判据字面常量 0.34/0.66/1e-14；TopAbs_ShapeEnum 复杂度排序与 rcad ShapeType 相反 → `occt_shape_rank` 换算。
7. **rcad 数据模型显式化（备案）**：DSFiller::insert 多 `&mut BRep` arena 参数（OCCT 全局 TShape 图）；kernel 缺口授权模式 = "仅当 OCCT 对应物存在而 rcad 缺，逐条报告"（H 补 3 个 BRepBuilder 方法；I 用文件内私有 stand-in 走通，未动 kernel）。
8. **两波契约化并行的时序模板**：主代理预读消费者源文（DSFiller.cxx 758 行）把对接点写进 brief → 上游叶子代理先行 → 消费者代理按契约并行 → 集成时主代理 diff 审计 + 怪癖抽查。
9. heredoc 大段追加继续绕行（Write 落盘 + `cat >>`）。

### 下一 session 入口（按序）

1. **Stage 3f**：HLRBRep_Data（2683 行拆 `brep/data/` 子模块，持 Arc\<BRep\>；落 3e 的 EdgeInterferenceTool `Data` trait；EdgeData::Set 的 TopoDS→CurveView 桥接；HLRBRep_Surface/Curve 与 Projector 的所有权收口）。
2. **Stage 3g**：Hider(852) → ShapeBounds → InternalAlgo(1020) → Algo → HLRToShape + 盒体/圆柱烟囱测试。
3. **Stage 4**：4a gen_hlr_ref.py（无头 vcomputehlr spike）→ 4b commands.rs → 4c occt-test-gen 接入 → 4d exact_hlr 3 用例闭环 + module-map 更新。

### 本 session 追加（2026-09-05，session 10：Stage 3f + 3g 全关，HLR 精确管线全链贯通）

| commit | 内容 |
|---|---|
| （本轮代码提交） | **3f+3g 全关**：`hlr/brep/data/`（Data 70 字段契约 + update/classify 2416 行/tableau_rejection 473 行）+ face_data + edge_face_tool + hider 1343 行（Hide 主循环全量 + EIT Data trait impl 接通）+ algo（组合+Deref）+ internal_algo 1654 行 + hlr_to_shape 1146 行（含 HLRBRep.cxx 包静态）+ shape_to_hlr + shape_bounds。两波六代理（K/L/M→N/O1/O2）+ uv_point 契约裁决；锚点 50+（data 20、face_data 3、edge_face_tool 6、hider 6、algo 4、internal_algo 4、hlr_to_shape 4、shape 6）；algo lib **349**（2 ignore = hider 双锚，阻塞在 intersector 层） |

**本 session 新确立的翻译事实（沿例勿改）：**

1. **struct 字段契约由主代理预落**（data/mod.rs 70 字段按 hxx L134-286 顺序）——三个代理并发写 impl 零冲突；子模块可访问父模块私有字段（seams.rs 先例）。
2. **跨代理契约冲突的裁决流**：冲突方 SendMessage 主代理 → 主代理改共享契约（Data.my_brep）→ SendMessage 双方接线责任——一轮消息收敛，不重发任务。
3. **`Arc<dyn Any + Send + Sync>` 作为 OCCT `handle<Standard_Transient>` 的不透明语义**（ShapeBounds.SDataHandle）——消费侧以 `down_cast` 取回。
4. **Data<'static> 的 leak 路线**（shape_to_hlr）：`brep.clone()` 与 curve adaptor 各 leak 一个 `&'static`（init_edge 的 Arc::as_ptr 先例同族）；HLR 会话对象数量有限，泄漏有界。
5. **OCCT 继承 → Rust 组合 + Deref/DerefMut**（Algo { internal: InternalAlgo }）——HLRToShape 的 `myAlgo->X()` 直接 deref 映射。
6. **catch_unwind(AssertUnwindSafe(...)) 对应 OCCT try/catch**（Hider 的 try 体抽私有方法，catch 分支忽略 Result）。
7. **怪癖照抄清单（session 10）**：InitEdge 的 `myHideCount` 双自增；Vertical 判据只查 0..6 维；`if (Maxz>Maxx){Maxx=Maxy;}`；`& 0x80008000` i32 位测（const MASK + wrapping_sub）；NextEdge 非测试路径不写 HideCount；ExploreShape 的 flag[0] 落槽；REVERSED `visible = fd.Back()`；heapsort `k <<= 1` 逐行；注释掉的 TriOk 冒泡排序保留为注释。
8. **kernel-context 显式化**：uv_point(brep, ..., edge, ...) 的 brep/edge 是 OCCT 全局 TShape 图的 rcad 显式参数（DSFiller insert 先例的推广）。

### 下一 session 入口（按序）

1. **intersector 交点语义对齐**（hider 双锚 un-ignore 的前置）：`Intersector::perform` 对 2D 直线横穿返回段而非带 In/Out transition 的点——对照 OCCT IntRes2d 的 $OCCT_SRC 行为逐分支核查。
2. **3g 收尾**：盒体/圆柱端到端烟囱测试（HLRBRep_Algo Add→Projector→Update→Hide→HLRToShape 全链 + nbshapes 结构断言）。
3. **Stage 4a**：gen_hlr_ref.py（勘察报告由代理 P 交付：VComputeHLR 调用链/无头可行性/6 自包含用例显式参数）→ 4b commands.rs → 4c occt-test-gen 接入 → 4d exact_hlr 3 用例闭环 + module-map 更新。

### 本 session 追加 2（2026-09-05，session 11：3g 烟囱闭环 + 4a 基础）

| commit | 内容 |
|---|---|
| （本轮代码提交 2） | **3g 烟囱闭环**：`hlr/tests.rs` 4 个端到端 smoke（盒体=OCCT 真值全对拍 + 圆柱部分锚定 + 稳定性）+ Q 的 hider 修复（双锚转绿，0 ignore）+ R 的 fclass2d IndexedVertexMap + kernel brep_lib/MakeEdge2d（14 锚点）+ S 的 MakeEdge2d 消费端接线 + `tools/gen-occt-ref/gen_hlr_ref.py`（4/6 JSON，bug25813_1 自检精确 PASS）+ 文档 |

### 下一 session 入口（按序）

1. **烟囱暴露的 5 缺口对齐**（Stage 4d 主攻前置）：③ Contap 圆柱第二轮廓（φ=225°）→ ④ seam 的 rg1_line（edges_to_faces 双面记录）→ ⑤ 共享边 Used 去重 → ① InternalAlgo arena 接线（产线化测试侧注入）→ ② fclass2d_topol CurveOnPlane 回退。
2. **4a 尾巴**：Plate/bug25813_3 的参考生成（DRAWEXE lprops 对 HLR 结果 ACCESS VIOLATION——C++ runner 直链 TKHLR+TKGProp 替代 Tcl，occt-bool-runner 方式）。
3. **4d**：exact_hlr 用例的 rcad 管线断言（消费 occt_hlr_*.json；盒体烟囱已证管线对 box 类真值精确）。

### 本 session 追加 3（2026-09-06，session 11 续：四缺口对齐 + 4a C++ runner 收官）

| commit | 内容 |
|---|---|
| （本轮代码提交 3） | **缺口①②③对齐 + 4a runner**：g_inter.rs DomainIntersection 语义（域外解丢弃/域端钳位 Head-End，+39 行）→ 圆柱双轮廓 + 每圆 3 段弧；kernel curve_on_plane（BRep_Tool.cxx L367-450 回退，~70 行）→ 无 pcurve 平面脸正确分类；InternalAlgo 会话 BRep（my_brep 捕获 + update clone，smoke 走自然路径断言不变）；tools/occt-hlr-runner（661 行 C++，box/bug25813_1 逐位吻合，ptorus viewer-path 302.685 定论 tuple 崩溃 = OCCT 8.0.0 exact-algo bug） |

**新确立的翻译事实（沿例勿改）：**

1. **域外解处理属求交器语义**：IntCurve_IntConicConic_1.cxx L641-714 的 DomainIntersection（丢弃/钳位/Head-End 赋位）是求交器（g_inter）的职责，不是消费者的——分类器射线上的域外交点必须在此层处理。
2. **回退语句要落在 OCCT 所在的层**：FClass2d init 无回退；BRep_Tool::CurveOnSurface 的 L367-372 "Try projection on plane" 才是回退点——修在 kernel 的 curve_on_surface 尾部，而不是给 FClass2d 加层。
3. **顺带发现并修正与 OCCT 相悖的旧断言**要给出二进制级证据（W/T 用安装版 OCCT 直接求证），不许凭推理改期望。

### 下一 session 入口（按序）

1. **烟囱剩余 3 缺口**：④ seam rg1_line（kernel 需存 edge-face regularity；锚点 ShapeToHLR.cxx L118-131 BRep_Tool::Continuity(F1,F2)→reg1/regn）→ ⑤ 共享边 Used 去重（HLRToShape DrawFace 的 Used 流转）→ Hider 层的底圆远半弧隐藏（hider/HidingStartLevel 对圆柱 case）。
2. **4d**：exact_hlr 消费 occt_hlr_*.json 的 rcad 断言（bug25813_1 204.19 + box 烟囱真值已备；Plate/ptorus 用 viewer-path 302.685/406.283 或 C++ runner 数值）+ module-map 更新。

### 本 session 追加 4（2026-09-06，session 11 续 2：缺口④⑤ + Hider 圆柱远半弧全关）

| commit | 内容 |
|---|---|
| （本轮代码提交 4） | **缺口④⑤ + Hider 远半弧全关**：kernel GeomAbsShape 真实枚举序 + CurveOn2Surfaces regularity 存储 + BRepBuilder::continuity + brep_tool_continuity 真实现（seam=CN → RgNLineVCompound）；Y 核实缺口⑤已在 3f/3g 关闭并加去重 pin 断言；kernel Lin::transform 方向修复 + Curve::First/LastParameter 的 Parameter2d 包装 → 圆柱 smoke 与 OCCT 真值逐位一致（V=5 弧 75.671 / H=1 弧 25.2237 / seam 在 RgN）；tktopalgo_gtests Trsf 字面量补齐（36/36） |

### 下一 session 入口（按序）

1. **4d**：exact_hlr 消费 occt_hlr_*.json 的 rcad 断言（`tools/gen-occt-ref/gen_hlr_ref.py` 已产 bug25813_1+poly×3；box/ptorus/Plate 用 occt-hlr-runner 的 ref_output 值）+ occt-test-gen hlr 接入（4c）+ module-map 更新。
2. **残余记录项**：④ 附带——`brep_tool_is_closed_edge_face`（ExploreFace 的 Dbl/IsReallyClosed）按 face ptr_id 匹配而非曲面值匹配，对 OutLinedShape 面副本会把 seam 的 Dbl 判 false（X 报告，不影响 typ 过滤）；烟囱 fixture 的 pcurve 注册注释已过时（U 报告）。

### 本 session 追加 5（2026-09-06，session 12：Stage 4c 全关 + 4d 断言就位 + 阶段2 再修 6 处）

| commit | 内容 |
|---|---|
| （本轮代码提交） | **4c 全关 + 4d 断言就位**：`tools/occt-test-gen/src/hlr_translate.rs`（根仓库）+ main.rs 4 接入点；`hlr/acceptance.rs`（box 全绿消费 ref_output.txt；bug25813_1/ptorus 就位 #[ignore]）；阶段 2 修复——Contap destination 上界（Contap_Contour.cxx L1645 `Array1(1,NbPointRst+1)`）+ vectg.Reverse（L856-858）+ kernel BRep_Tool::CurveOnSurface 闭合缝 REVERSED 取 PCurve2（BRep_Tool.cxx L339/353-357）+ fclass2d 有向 occurrence 走查 + TopoDS EmptyCopy 旗标重置（TopoDS_TShape.hxx L169）+ surf_function d2d 归一化候选修复（cxx L271 gp_Dir2d）。基线 algo lib **364+2i** / builder **76+1** / pavefiller **26** / tktopalgo **36** |

**剩余两缺口的机制级诊断（下一 session 直接续打，勿重推导）：**

1. **bug25813_1 顶缘弧隐藏区间丢失（差 16.9334）**——诊断链已闭合到最后一步：
   - 已验证正确：orient_out_line 改写 outline 边 orientation（Data::Update L843）→ next_interference 可见它们；求交参数 (1.4289, 3.2835) **精确正确**（x'=±8 穿越点，Parameter3d 转换在 rejected_point L2270 已有）；Builder 产出 mask=3 的完整区间 [1.4289, 3.2835]（其 2D 弧长 ≈16.93 与 OCCT 隐藏块吻合）；hiding_start_level=-2 两干涉都通过 Level 过滤。
   - 卡点：`classify(2.356)`（区间中点）被 k=7（dim14/15=视图深度 z）编码盒门拒绝 → aTestState=OUT → `ES.Hide` 被跳过（OCCT Hider.cxx L609-628 同构）。中点 3D (-7.07,7.07,30) 深度 9.16 确在 face4 盒 [10.79,36.56] 外。
   - **矛盾（下一 session 首查）**：按射线几何（视线 P0+u·(1,-1,1)，u∈[1.415,12.7] 穿过小圆柱径向脚印）中点**应被遮挡**，但深度比较显示小圆柱面在视线上的深度 (20.1) **大于**点的深度 (15.86)——即点在遮挡面之前，几何上不该被隐藏。而 OCCT 真值 204.19 里确有 16.93 被隐藏。⇒ ①用安装版 OCCT 直接验证被隐藏的 16.93 到底在哪条边（此前"OCCT 真值分解"是代理从 DRAWEXE 提取的间接结论，与射线分析矛盾，二者必有一错）；②核对 rcad `Project(P,x,y,z)` 的 z 语义与 OCCT HLRAlgo_Projector::Project 逐行（含 dp/ perspective 分支）。
2. **ptorus IWalking 0 线**——Contap 侧已通（8 restriction + 20 interior + 8 departure）；surf_function d2d 归一化候选修复已落（cxx L271）但未验证；RCAD_IWALK_DEBUG 探针已埋（i_walking/function_set_root）。首查走行方向/步进在圆环参数化下的符号与归一化。
3. **小圆柱 seam 落 Iso 而非 RgN**（生成 harness 计入 Iso 多 17.97）：nbIso=0 时 OCCT IsoLineVCompound 恒空；疑似 bfuse 丢失 seam 的 CN regularity（bop 层 CurveOn2Surfaces 传递）或 Data 的 Iso 旗标判定。cylinder smoke 的 seam 归 RgN 正确，对照融合前后 regularity 存储即可定位。

**临时探针（下 session 修完即删，已随本提交入库存档）**：`RCAD_HLR_TRACE`（classify.rs/update.rs/internal_algo.rs 的 trace_enabled 块 + hider.rs builder 探针 + classify REJECT 探针）、`RCAD_IWALK_DEBUG`（i_walking/function_set_root 的 TEMP-DEBUG dump）。

### 下一 session 入口（按序）

1. 按"剩余两缺口的机制级诊断"续打：先用安装版 OCCT（tools/occt-hlr-runner 加 Hider/Data 层探针或 debug 构建重编 TKHLR）裁决 16.93 的真实归属，再修 classify/Project 语义或接受几何结论；ptorus 验证 d2d 候选修复 + IWalking 方向链；seam regularity 经 bfuse 的传递。
2. 三缺口全绿后：删全部临时探针 → un-ignore 两个 acceptance 测试 → `cargo test -p occt-generated-tests --test generated_occt_boolean_hlr_exact_hlr` 全绿 → module-map 基线数字终版。
3. 残余记录项（非阻塞，同 session 11 清单）不变。

### 本 session 追加 6（2026-09-06，session 12 续：ptorus 连破三关，断点收敛到 DS 边提取）

| 修复 | 内容 |
|---|---|
| d2d 归一化（已验证生效） | `contap/surf_function.rs`：`d2d = (-fpv/d, fpu/d)`（OCCT cxx L271 `d2d` 是 gp_Dir2d 构造即归一化；旧代码存未归一化向量）→ IWalking 从 0 线变 52 线 |
| init_edge MST 过期 arena（已验证生效） | `brep/data/update.rs`：MST 命中分支**不再覆盖** `Data::my_brep`——工具快照拍摄于 ds_filler::insert 入口，早于 FaceIsoLiner 轮廓边入 arena（len 7 vs 索引 12 越界 panic 的根因）；kernel 上下文恒为 Data::update 设置的完整会话 clone（OCCT 语义：MST 读全局 TShape 图，永远最新） |

**ptorus 现状**：IWalking `nb_lines=52`（含重复片段，起点去重疑似失效——OCCT 用 etat 取负标记 crossing point 防重建，rcad 已对齐该机制但 52>1 说明仍有偏差）→ Contour `nb_lines=51`（采纳正常）→ **断点 = FaceIsoLiner→DS 轮廓边→HLRToShape 提取段**（visible mass 仍 0；注意隐藏阶段已经能引用到轮廓边——说明 DS 边存在，疑点收窄到边的 3D 曲线/range 建造或 OutLine 旗标）。探针已埋好：`RCAD_IWALK_DEBUG`（i_walking domain/wd1/wd2/perform end + contour perform.rs 的 CWDBG after line_constructor）。

### 本 session 追加 7（session 12 续 2：方法论回归——源码逐行审查替代运行时探针）

按 AGENTS.md 阶段 2 纪律（不 println、从拓扑差异回溯到函数、逐行核对 OCCT 源码）完成 ptorus 断链的静态审查，结论：

**已逐行核对为 1:1 的链路（从 Contap 到提取的 8 个环节）**：
1. `contap/contour/perform.rs`（fo.perform 内部：52 线采纳 51）
2. `topo_brep/ds_filler.rs` insert（fo.perform 调用点 L175，done=true nb_lines=51/26）
3. `ds_filler.rs` insert_face（51 条全部 IType::Walking；Walking 分支造边 + `add_int_l(f).push(e)` L751）
4. `topo_brep/data.rs` is_int_l_face_edge / is_spl_e_edge_edge（cxx L116-131 形式一致）
5. `out_liner.rs` process_face（IntL 块：INTERNAL 化 → SameEdge 检查 → B.Add(W,E) → B.Add(NF,W)，与 cxx L157-257 一致）
6. `out_liner.rs` build_shape（cxx L307-340 一致）
7. `brep/shape_to_hlr.rs` explore_shape（**OriginalShape 的 shell 探索**，传原始面 F 给 explore_face——cxx L254-310 一致）
8. `shape_to_hlr.rs` explore_face（DS 查询用原始面 F、拓扑遍历用 FM(i)=NF——cxx L187-243 一致；Int 旗标 → SetWEdge）

**RunTime 实测（CDBG，生成 harness 临时探针）**：RgNV=2（缝，来自 Load 的 regn=CN ✓），V/OutV/Rg1V/IsoV/H 全 0。typ=2（OutLine）提取条件 = 面局部 `Itf.Internal()` 旗标（HLRToShape.cxx L176-186）→ 该旗标来自 SetWEdge 的 `Int = IsIntLFaceEdge(原始F, E)`。

**剩余唯一疑点（下 session 首查，读源码可决）**：IntL 的键与内容在 UpdateEdgeData 时刻的同一性——
(a) `ds.add_int_l(f)` 的 f（insert 探索路径的原始面）与 explore_face 的 f（explore_shape 探索路径的原始面）ptr_id 是否同一 Arc；
(b) insert_face 压入 IntL 的边 TShape 与 process_face 写入 NF 的边 TShape 是否同一（make_edge 是否被调用两次各造一份）；
(c) 另发现真缺口：**rcad 缺 `HLRToShape::OutLineHCompound` 访问器**（hxx L122，ViewerTest L3290 消费）——补齐属形式完成项。

### 本 session 追加 8（session 12 终：ptorus 根因实锤——TShape Arc 身份分裂）

探针实证（SFDBG/CWDBG）：IntL 非空（51/26 条，键 ptr 正确）、建边链全 in-place（add_edge/add_to_edge/update_edge_pcurve 均 edge_mut_inplace）——但 **explore_face 探索到的 NF 边 Arc ≠ IntL 边 Arc**（e_same_any=false 全体成立）。OCCT 语义：B.Add 存同一 TShape 句柄、BRep_Builder 原地变更、身份永不分裂；rcad 的 `wire_mut`/`face_mut`/`shell_mut` 仍是 `Arc::make_mut`（topods.rs L1652-1670）——当 wire/face 的 Arc 被多处引用（arena + 局部句柄 + 父容器）时 make_mut 克隆分身，导致 NF 实际存储的 wire/edge 与 IntL/DS 记账的句柄分裂 → Int 旗标永不命中 → 轮廓边从所有 compound 消失（仅 RgN 缝存活的准确解释）。

**修复方向（下 session 首任务）**：容器类（wire/face/shell/compound）的变异器与 builder_add 家族改为身份保持语义（仿 edge_mut_inplace 的裸指针原地变异，或建后不再变异的构造式用法）；审计清单 = topods.rs L3215-3320 全部 make_mut 变异器 + builder_add/make_wire/make_shell 调用序。修完 SFDBG 的 e_same_any 应转 true → RgNV=2 之外应出现轮廓边 → 对照 OCCT 302.685。

### 本 session 追加 9（session 13：容器变异器已改 in-place；断点再收窄一环）

1. **已完成**：`topods.rs` 的 wire_mut/face_mut/shell_mut/solid_mut 全部从 `Arc::make_mut` 改为身份保持的原地变异（同 edge_mut_inplace 契约）。
2. **实测排除**：修后 e_same 仍 false → 分裂点不在容器变异器。进一步实证：被探索的 NF 边是**晚于 IntL 的新 TShape**（ptr 段位证明）——它们是 `process_edges` 的 **SplE 分段**（`empty_copy` 新建，OCCT cxx L713-758 忠实翻译，本身正确）。
3. **真正的机制（已还原）**：OCCT 的轮廓边查找走 `IsSplEEdgeEdge(IntLEdge, E)`——IntL 存原始边，NF 存 SplE 分段，通过 `EdgeHasSplE(IntLEdge) → SplE 列表含 E` 匹配。rcad 的 is_spl_e_edge_edge 已对齐。断点 = **Walking 分支的内部点没有走 insert_vertex 注册**（→ SplE 表按边注册依赖 init_vertex，若无注册则 process_edges 不为这些边建 SplE / 或建了但查询回退原边比较必 false）。对照点：OCCT DSFiller cxx L161-171（walking 分支的 `if (P.IsInternal()) InsertVertex(P, tol, E, DS)`）与 rcad insert_face 非 Restriction 分支的对应语句。
4. 下 session 首查：rcad insert_face 非 Restriction 分支的 insert_vertex 调用与 OCCT L146-192 逐行对照（重点 L161-171 内部点注册），补齐后 SFDBG 的 int 应转 true。

### 本 session 追加 10（session 13 续：SplE 机制还原完毕，收口点钉死）

静态还原完成 OCCT 的完整匹配机制：insert_face 把 IntL **原始边**推入列表；process_edges（insert 尾）为 DS 注册表中**有内部点注册**的边造 SplE 分段（empty_copy 新 TShape）；process_face 的 IntL 块对 EdgeHasSplE 的边**改加 SplE 分段**；提取查询 IsIntLFaceEdge 走 IsSplEEdgeEdge（SplE 列表含被探索边 → true；无 SplE → 回退原边 IsSame）。

单次运行实测：IntL=51（原始边 ptr 1b08fe6-ec 段）、被探索边 = 1b08fef-ff 段（= SplE 分段，证明 process_face 的 EdgeHasSplE 分支**生效**）→ 但 is_int_l_face_edge 仍 false ⇒ **SplE 映射的键（DS my_edges_vertices 注册表的边句柄）与 IntL 边句柄不一致**，或 SplE 分段与被探索分段非同一 Arc。收口动作（纯读源码可决，无需探针）：(1) rcad data.rs my_edges_vertices 的注册键来源——insert_face Walking 分支谁调 init_vertex（对照 OCCT：walking 分支无 InsertVertex，注册表应仅有 restriction/iso 边 → 则 SplE 空 → process_face 应走 else 加原始边 → 与实测"探索到 SplE 分段"矛盾——即 rcad 的注册表被某种自播种污染，process_edges 对未注册边造了 SplE——对 OCCT L713-758 的 InitEdge/MoreEdge 语义逐行核对即可定位）；(2) 修复后验证 = SFDBG int=true → OutLineV 非空 → mass≈302.685。

### 本 session 追加 11（session 13 终 2：终极对照数据 + 下一动作）

终极探针数据（单次运行）：`has_int_l=true int_l_len=51 has_spl=false`（IntL 边的 SplE 表**为空**——即 process_face 走的是 else 分支、把 **IntL 原始边** 加进了 NF）→ 按此，探索到的边必须就是 IntL 原始 Arc、IsSame 回退必命中。但实测被探索边 ptr（2bca3f2xxxx，连续段位）≠ IntL ptr（2bca3ea-ee 段）——**"process_face 存入的句柄"与"探索器产出的句柄"不一致**。

下一动作（精确到函数，读源码即可决）：(1) rcad `TopExpExplorer`（brep 遍历器）——`current()` 返回的是父 TShape 里存的 child 句柄，还是按 index 从 arena 重解析的句柄；(2) `insert()` 的 `top_exp_faces(s,...)` 与 `load()` 的面探索是否产出同一 Arc 的面句柄（两处探索路径的 Arc 一致性）；(3) 71 条被探索边实例 vs 预期 53（2 缝 + 51 轮廓）——多出的 18 条从哪个 wire 来。三者都是纯读码核对，无歧义空间。

### 本 session 追加 12（session 14：Int 旗标已全部正确——断点移至 Hider 自隐藏）

容器变异器 in-place 修复 + 全链对齐后，SFDBG 实测：**51 条轮廓边全部 int=true / e_same_any=true**（e_kind=B BSpline，idx 9+，IntL 注册链完全打通）；被探索的缝边 SplE 分段（idx 150+，Circle 类）为 process_edges 正常产物。**补齐真缺口 `HLRToShape::OutLineHCompound/OfShape`（lxx L140-150）。**

**当前断点（最后一段）**：轮廓边旗标全对，但 visible mass 仍 0 → **Hider 对 torus 自隐藏（单面 solid 的轮廓被自己的面隐藏）把全部轮廓判成了 hidden**（OCCT 应为近半 visible 302.685 / 远半 hidden）。下一动作：对照 Hider.cxx 自隐藏路径的 Classify 深度语义（dz 符号 / myProj.Project 的 z 约定 / RejectedPoint 的 dz≥TolZ 判据 L2349-2361）与 rcad classify.rs/hider.rs 对应段；重点核对 z 是否"沿视线方向深度"及 near/fare 判定符号。

### 本 session 追加 13（session 15：旗标全通、缝分割生效；唯一断点 = Hider 自隐藏 CS 深度）

容器变异器修复的实效已验证：RgN 可见缝从 2 → 5（缝被 SplE 正确分段后提取），RgNH+IsoH=1（隐藏侧也出现了）。轮廓边 51 条全部 int=true 注册。DrawFace/DrawEdge/InternalCompound（含 Used 重置、HideCount 延迟、leftover 遍历）逐行对照均一致。

**唯一剩余断点**：torus 轮廓边的 EdgeStatus 全部 hidden（visible 区间空）。机制上 = 自隐藏 classify 的射线-面 CS 深度判定把近半轮廓也判了 IN。剩余对照 = `my_proj.shoot` 射线方向/参数化 + `elclib::line_parameter` + `perform_line` 多面体求交的 w 值符号（OCCT Data.cxx L2192-2262 已确认；rcad classify.rs L2059-2091 结构一致）→ 需打印一次 near-half 轮廓片段的 (w, w_lim) 实测值对拍即可定位（classify hits 探针已埋，trace_cl 开关）。修复后验证：OutV 非空 → mass≈302.685。
