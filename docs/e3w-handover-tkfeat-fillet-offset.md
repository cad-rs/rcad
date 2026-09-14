# E3-W 交接：三域（TKFeat / TKFillet / TKOffset）翻译推进 —— 2026-09-12

> **一句话（追加 39 收尾态，2026-09-14 —— 最新，读这一行即可开工）**：**工作口径仍是"先译后调"**。本轮单线（epsilon 同族全库收敛 + canonical 第三轮修正 + 队列 3 对拍）：**① canonical `epsilon_of` 修正为 OCCT 精确 nextafter 体**（位递增体在 x==±RealLast ⇒ OCCT 给 **0**、±inf ⇒ **-inf**、-0.0 ⇒ **+5e-324** 三类边缘偏离；`hlr/intrv::Interval::new()` 恰好用 `epsilon(±RealFirst/Last)` 作默认容差，先收敛后修 canonical 会把 HLR 区间容差从 0 污染成 +inf ⇒ **先修 canonical 再收敛是硬顺序**）· **② 16 份本地 Epsilon 副本全库收敛到 kernel canonical**（队列点名的 4 份 + 名无关终检兜底再抓 6 份漏网 + 普查扩大 6 份；每站点回源核验 OCCT 调用方；其中 3 份是 `f64::EPSILON*v` 相对式**公式真偏离**）· **③ `bspl_lib.rs` 6 处 `Epsilon(1.)*x` 公式偏离修复**（`Epsilon(1.)*x` ≠ `Epsilon(x)`，仅 x 为 2 的幂时相等；两处还多了非 OCCT 的 `.abs()`）· **④ 队列 3 对拍收口**（16 域网格逐格同基线；`feat_featlf` a3 失败层 = `tool_rehost.rs:1117` "Courbes non jointives"，不在 bean_face 爆炸半径链上；`loc_ope_split_drafts` 无直测——记档）。
> **门槛（当前基线）**：**450/0/0 · 737/0 · 36/36 · 26/26 · 76/76 · 1/1**（kernel 736 → **737** = 1 个 canonical 边缘判别测试）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。
> **★★ 本轮最值钱的发现：canonical 的边缘语义必须对照 C 标准库逐输入推演。** "有限输入上等价"不等于"函数等价"——位递增体在全部正常输入（含追加 38 的 gtest 与判别性测试）上与 OCCT **逐位相同**，却在 x==y/±inf/-0.0 四类边缘上偏离，而**唯一的现役消费者恰好把这些边缘值当输入**（`Intrv_Interval` 默认容差 = `(float)Epsilon(±RealFirst/RealLast)` = OCCT 的 0.0）。另：**"测试失败先怀疑测试"再 +1**——canonical 边缘测试第一版期望值 `epsilon_of(+inf) == -RealLast` 写错，真值是 `-inf`（`RealLast() - inf` 的 IEEE 有限减无穷）。
> **★ 第二条：名无关终检必须做两遍。** 先按实现拼写 grep（`nextafter`）普查得 17 文件；收敛完按函数名模式 `fn (epsilon|standard_epsilon|...)\(` 终检**又抓到 6 份漏网**（`next_up()`/位递增拼写，文本里没有 "next_after"）——坑 28/追加 33 第四次应验。共收敛 **16 份**，其中 3 份是 `f64::EPSILON*v` 相对式**公式真偏离**（自称 OCCT Epsilon、实为别的公式，阈值差 25%~57%）。
> **★ 第三条（零可见翻转第十二轮）**：9 处真实阈值修正（最高 ~1.57 倍）在八网格与 16 个域网格**全部零可见**——这些比较点在既有用例里都不在阈值边缘；canonical 修正是**潜在险情排除**（险情路径在 HLR 域、两套网格都照不到）⇒ "判别性单测是唯一验收手段"第十三次重申。
> （历史：追加 21–38 收尾态见下方存档块；追加 38 = OCCT `AdjustPeriodic` gtest 移植 + `epsilon_of` 复制族首收敛 + `bean_face_intersector` 误标翻译；追加 37 = `Ax3` 平行分支 + `Init` 的 location + `AdjustPeriodic` 收敛；追加 36 = `gp_Trsf2d` 全类 + 判别性单测；追加 35 = BSpline iso 精确化 + `CurveOnPlane`。）

> **（存档）一句话（追加 38 收尾态，2026-09-14 —— 已被追加 39 取代）**：**工作口径仍是"先译后调"**。本轮三线并行（测试移植 + 收敛 + 误标修复）：**① 移植 OCCT 自带的 `AdjustPeriodic` gtest**（canonical 体的唯一 ground truth，此前 rcad 对它零测试）· **② `epsilon_of` 复制族收敛**（**前提被推翻**）· **③ `bean_face_intersector` 的 `AdjustPeriodic` 误标**（确认为真，按 OCCT 1:1 翻译）+ 两处返回 `MIN_POSITIVE` 的 `epsilon` 修正。
> **门槛（当前基线）**：**450/0/0 · 736/0 · 36/36 · 26/26 · 76/76 · 1/1**（algo 451 → **450** = 删掉的自失效测试；kernel 733 → **736** = 3 个新测试）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。
> **★★ 本轮最值钱的发现：OCCT 自带的测试并不覆盖它自己的守卫。** 把 `Epsilon(ULast)` 那行守卫**整个禁用**，OCCT 那五个用例**照样全过** —— 而该守卫正是 `epsilon_of`（追加 37 刚修的那个公式）的**唯一消费者**，此前**零覆盖**。补的判别性测试用**零长度区间**（无守卫 ⇒ 除以 0 ⇒ **NaN**）；我第一版用"周期很小但不为零" **不判别**（wrap 恰好落回同样的值，守卫禁用下仍通过）。⇒ **移植 OCCT 自带测试 ≠ 覆盖完整，必须自己扰动一遍。**
> **★★ 第二条：前提推翻是"好消息型"。** 记档说"同一个错误公式被复制到多处"，实测 **HEAD 里每一份副本的公式本来就是对的** —— 错误公式**只存在于 canonical 体本身**（追加 37 刚修的那个）。⇒ **"多处复制了同一错误"必须先逐份核验再动手**，否则会去"修"一批本来正确的代码。
> **★ 真实的翻译缺口（批次 3）**：`bean_face_intersector` 里那个 `AdjustPeriodic` 是自创重写（`1e-12*period` 容差 + `ceil` 夹取 + 自创 `ok` 门），OCCT 是"无条件 `solutionIsValid = true` 并忽略返回值"⇒ 已 1:1 翻译。**三处实质变化**：修剪柱面（UV 窗口窄于周期）上原本被丢弃的交点现在保留；守卫 eps 6.28e-12 → 8.88e-16；平移公式 `ceil` → `modf/trunc`。
> **★ 方法学（第十一轮"零可见翻转"）**：批次 3 改变了接受集合与下游 UV 值（爆炸半径 `LocOpe_SplitDrafts`）、批次 4 改了两个容差 helper 的零行为 —— 八网格与 16 个域网格**全部零可见**。
> （历史：追加 21–37 收尾态见下方存档块；追加 37 = `Ax3` 平行分支 + `Init` 的 location + `AdjustPeriodic` 收敛与 `epsilon_of` 修正；追加 36 = `gp_Trsf2d` 全类 + 判别性单测；追加 35 = BSpline iso 精确化 + `CurveOnPlane`。）

> **（存档）一句话（追加 37 收尾态，2026-09-14 —— 已被追加 38 取代）**：**工作口径仍是"先译后调"**。本轮三线并行（内核翻译补全 + 收敛 + 勘误）：**① `Ax3` 补回 OCCT 的两个平行分支**，并把 `direct()` 从存储字段改回**推导式**（删掉陈旧 `sense` 模型）· **② `Init` 的 `my_location` 补齐**（追加 36 报出的更大缺口，记录为真）· **③ `AdjustPeriodic` 收敛**（**记档位置写错**：那份在 `fillet/**` 不在 `bop/**`，且 OCCT 有**四个**同名函数）+ **canonical `epsilon_of` 公式 1:1 修正**。
> **门槛（当前基线）**：**451/0/0 · 733/0 · 36/36 · 26/26 · 76/76 · 1/1**（kernel 728 → **733**）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。
> **★★ 本轮最值钱的产出仍是"三次回源核验救了命"**：① `AdjustPeriodic` 的记档**位置**写错（在 fillet 不在 bop）；② `bop/**` 那两处**根本不是同一个 OCCT 函数**（OCCT 有四个同名 `AdjustPeriodic`，判"是否同一个"必须逐站点判）；③ **canonical `epsilon_of` 的公式与锚点都是错的** —— 它自称 `Precision.hxx` 的 `Epsilon`、公式 `Max(|v|, RealSmall()) * RealEpsilon()`，而 `Precision.hxx` **没有** `Epsilon`，OCCT 唯一的 `Epsilon(double)` 是 **nextafter 间隙**（`Standard_Real.hxx` L242-246）。⇒ 与追加 35/36 合起来，**"canonical / 记档"必须逐条回源核验**已是本域的**常态动作**。
> **★ 真正的翻译缺口（本轮批次 1）**：`Ax3::set_direction` / `set_x_direction` **各缺一个平行分支** —— 给定方向与主方向（或当前 X）(反)平行时 `V ^ (X ^ V)` 退化为零，OCCT 改为**重贴轴标签**，rcad 旧版却代入**零向量** ⇒ 产出**退化框架**；`direct()` 也从存储字段改回 OCCT 的**推导式**（`(vxdir ^ vydir)·N > 0`）。
> **★ 方法学（第十轮"零可见翻转"）**：内核形态学（`direct()` 由存储改推导）、`Init` 的 location、`AdjustPeriodic` 收敛、`epsilon_of` 公式 —— 八网格与 16 个域网格**全部零可见**。
> （历史：追加 21–36 收尾态见下方存档块；追加 36 = `gp_Trsf2d` 全类 + `gp.rs` 判别性单测 + 拆文件与 location 勘误；追加 35 = `gp_Trsf2d` 译完 + BSpline iso 精确化 + `CurveOnPlane`；追加 34 = 判别性单测 + `Trsf2d` 两缺陷 + `Trimmed`/`dn_at`。）

> **（存档）一句话（追加 36 收尾态，2026-09-14 —— 已被追加 37 取代）**：**工作口径仍是"先译后调"**。本轮三线并行（翻译 + 合规 + 勘误）：**① `gp_Trsf2d` 全类译完**（余下 10 函数 + `Power`，载体抽为独立模块 `trsf2d.rs`）· **② `gp.rs` 判别性单测补齐**（12 → 17，含 3 次扰动实证）· **③ `brep_fill_sweep.rs` 拆回 2000 行内 + location 比对勘误**（记档为真，已修）。
> **门槛（当前基线）**：**451/0/0 · 728/0 · 36/36 · 26/26 · 76/76 · 1/1**（algo 438 → **451**、kernel 723 → **728**）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。**改动文件全部 ≤ 2000 行**。
> **★★ 本轮最有价值的产出是三条"OCCT 自身怪癖"的如实记档**（都容易被人"顺手修正"成与 OCCT 不一致）：**（a）`gp_Trsf2d::Invert()` 对 `gp_PntMirror` 不是数学逆** —— 只反 offset，中心对称绕 `c` 变成绕 `-c`（`t_inv(t(p)) == p - 4c`），**3D 的 `gp_Trsf::Invert` 完全一样**（`gp_Trsf.cxx` L406-409）⇒ 这是 OCCT 的怪癖，**测试要钉 OCCT 行为、不要钉"往返等于原值"**；**（b）`SetValues` 的 `M.Divide(s)` 经公开 API 不可观测**（后续 `Orthogonalize` 归一化把正因子精确抵消）⇒ **如实记为不可观测，不假装覆盖**；**（c）`SetTranslationPart` 从 `gp_Rotation` 出发得到 `CompoundTrsf`**（cxx L113-116 的 catch-all）。
> **★ 新形态的问题：本轮出现三处"测试写错"**（我两处 + 子代理一处），**全部是"代码对、期望值错"**（`value(3,3)` 误算成 0；`SetTranslationPart` 的标签；`Invert` 对 PntMirror 的非逆）。⇒ **"测试失败先怀疑测试"升为常规动作**，且**期望值必须来自 OCCT 公式或独立几何定义，不能来自"我以为"**（已写入 §0）。
> **★ 两条卡的结局**：`gp.rs` 的 `Trsf` 助手**已全部有判别性测试 ⇒ 关闭**（下一步把同一手法用到其他内核助手）；**`mySn` 卡证伪 ⇒ 关闭**（`mySn` 只在 `LocOpe_Gluer` 且 rcad 侧已完整；`feat_featrf` 的真实前墙是 `loc_ope_glued_shape.rs:193` 的 `Standard_ConstructionError`，属**调试/状态类**而非翻译缺口）。
> **★ 方法学（第九轮"零可见翻转"）**：`Trsf2d` 全类、location 比对、守卫勘误在八网格与 16 个域网格上**全部零可见**。
> （历史：追加 21–35 收尾态见下方存档块；追加 35 = `gp_Trsf2d` 译完 + BSpline iso 收敛到精确体 + `CurveOnPlane`；追加 34 = 判别性单测 + `Trsf2d` 两缺陷 + `Trimmed`/`dn_at`；追加 33 = `UpdateEdge` 全库收敛 + 两个内核缺陷。）

> **（存档）一句话（追加 35 收尾态，2026-09-14 —— 已被追加 36 取代）**：**工作口径仍是"先译后调"**。本轮按你的指令**三线并行且全取翻译/收敛类**：**① `gp_Trsf2d` 载体译完**（`TrsfForm` 补真实 2D 形态 + `Multiply` **全 15 分支** 1:1）· **② BSpline `UIso/VIso` 的三份 de Boor 近似收敛到精确 `BSplSLib::Iso`** · **③ `CurveOnPlane` 回退补齐**（记档为真）。
> **★★ 本轮最有价值的产出是"两次推翻前提"**：**（a）** 队列说 `TrsfForm` 缺"四个"变体 —— 实际 `gp_Trsf2d` **没有 `Ax2Mirror`**（点对称/线对称才是 2D 的），只缺**三个**；**（b）** 批次 2 的任务前提"canonical 体已是精确翻译、两份副本是近似"**只有一半为真** —— 所谓 canonical 体**自己就是第三份 de Boor 近似**。子代理因任务书里"禁止收敛到另一份近似"这一条而先读了 OCCT，**否则就会把两份近似合并成一份近似：看起来收敛了、实际什么都没修**。⇒ **动手前必须验证记档与任务前提**（本轮升为与"判别性必须实证"并列的硬规矩）。
> **门槛（当前基线）**：**438/0/0 · 723/0 · 36/36 · 26/26 · 76/76 · 1/1**（algo 430 → **438**）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。
> **★ 本轮量出的真实行为差异**（此前一直错着）：BSpline iso 在**参数域外**旧体夹到边界极、精确体按请求**外推**（实例极差 **0.894**）；**权重全为 2.0 的均匀有理面**旧判据不是 OCCT 的 `Rational()`（权重留 2.0 vs 写 1.0）；**周期面**旧体完全没有周期处理且硬编码 `is_periodic: false`（当前潜伏）。
> **★ 方法学（第八轮"零可见翻转"，连"改了分支走向"也不可见）**：批次 3 把一条分支从"返回 `None`"改成"投影返回"，批次 2 换掉整条求值公式 —— 八网格与 16 个域网格**依然全部零变化**。
> **★ 可复用手法（本轮新增一条）**：**代数自洽 ≠ 几何正确** —— 组合律测试（`A.multiplied(B) == A(B(P))`）**抓不到**矩阵转置（它对同一个错矩阵自洽），只有断言独立几何量的测试能抓；**两类断言都要留**（另有"测试失败先怀疑测试"一例：我写的 `value` 断言误以为 `set_scale` 保留矩阵，实际 OCCT 会 `matrix.SetIdentity()`）。
> （历史：追加 21–34 收尾态见下方存档块；追加 34 = 判别性单测 + `Trsf2d` 两缺陷 + `Trimmed`/`dn_at`；追加 33 = `UpdateEdge` 全库收敛 + 两个内核缺陷；追加 32 = `set_rotation` 转置修复。）

> **（存档）一句话（追加 34 收尾态，2026-09-14 —— 已被追加 35 取代）**：**工作口径仍是"先译后调"**（先 1:1 译全 body、把重复载体收敛到真身，代码基本译完再回头调试）。
> 本轮（追加 34）**三线并行，做队列第 1、2、4 项**：**① `gp.rs` 判别性单测补齐（4 → 12 个，含 5 次扰动实证）** · **② `composite_surface.rs` 的本地 `Trsf2d` 审计 —— 再抓 2 个同族缺陷（其中一个可达）** · **③ `Trimmed` 的 `dn` 臂 + `GeomAdaptor_Surface::dn_at` 全阶（记档的"无 OCCT 对应"被推翻，panic 清零）**。另把缺失的 `Trsf::SetTransformation(FromA1,ToA2)` 重载译全。
> **★ 本轮最值钱的发现（动摇了一次"测试已覆盖"的判断）**：`transform_vec`（`VectorialPart`，含 scale）与 `transform_dir`（`HVectorialPart`，原始 matrix）这对配对的**第一版测试在扰动下仍然通过 —— 根本不具判别性**。原因是 `transform_dir` 末尾的**归一化抵消了正缩放**，只有**负缩放**经 `scale<0` 的反转才暴露（双重取反）。改成同时跑 `scale = ±3` 后才失败。⇒ **不写扰动实验就会留下"看似覆盖、实则恒真"的假回归测试**。
> **门槛（当前基线）**：**430/0/0 · 723/0 · 36/36 · 26/26 · 76/76 · 1/1**（algo 427 → **430**、kernel 714 → **723**，全是新测试）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与追加 32 基线逐项相同**。
> **★ 方法学（第七轮"零可见翻转"，本轮形态最极端）**：`dn_at` 从**一个 panic** 变成**真的算值**、`set_scale` 从**一次 abort** 变成**返回退化变换** —— 八网格与 16 个域网格**依然全部零变化**。前几轮还能解释成"改了但没走到"，本轮是**把不可达路径变成可达**仍然零翻转 ⇒ **两套网格合起来的"覆盖率"对本域进展几乎无分辨力，"判别性单测"应当作唯一验收手段**（第七次重申，且口径已从"补用例"改为"补判别性单测 + 扰动实证"）。
> **★ 可复用手法（本轮新增两条）**：**① 判别性必须用扰动实证**（改回旧形 → 确认目标测试 FAILED → 还原；"我推导过它会失败"不算数 —— 本轮同一轮里手工推导与"看起来该失败"各错过一次）；**② "注释说缺 ≠ 真缺"第三次应验**（`Trimmed` 的 `EvalDN` 一直在 OCCT 里，记档却写了"无对应"）。
> （历史：追加 21–33 收尾态见下方存档块；追加 33 = `UpdateEdge` 全库收敛 + `set_displacement`/`multiplied` 两个内核缺陷；追加 32 = `set_rotation` 转置修复；追加 25 = `blend_simple_a1` 通过；追加 28 = `LocOpe_SplitShape` 全类；追加 30 = `Geom_OsculatingSurface`。）

> **（存档）一句话（追加 33 收尾态，2026-09-14 —— 已被追加 34 取代）**：**工作口径仍是"先译后调"**（先 1:1 译全 body、把重复载体收敛到真身，代码基本译完再回头调试）。
> 本轮（追加 33）**把追加 32 队列的第 2、3 项一次做完**：`UpdateEdge` 家族**全库收敛完毕**（**9 个站点**改走唯一真身 `topods::update_curves_range`，含**审计清单外的 `fillet/hbuilder.rs`**）+ `gp.rs` 矩阵约定**全量审计**（813 行全读、**32 项判 ALIGNED**）。
> **★ 本轮最重要的产出：再挖出并修掉 2 个同类内核缺陷** —— **`Trsf::set_displacement` 的两处转置乘**（`gp_Trsf.cxx` L218-240）与 **`Trsf::multiplied` 漏乘 `scale`**（`gp_Trsf.cxx` L544-559）。两处都**潜伏**（唯一调用点恰好落在退化配置上），已补**判别性**回归测试。
> **门槛（当前基线）**：**427/0/0 · 714/0 · 36/36 · 26/26 · 76/76 · 1/1**（kernel 712 → **714** = 两个新测试）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格的真实断言数与追加 32 基线逐项相同**。
> **★ 方法学（第六轮同一结论，本轮证据最硬）**：本轮修的是**真实内核缺陷**，却在八网格与 16 个域网格上**全部零可见** —— 因为**唯一调用点恰好都落在缺陷的退化配置上**（`Ax3::new()` = 零位置 + 单位姿态；`sweep_section_generator` = scale 1）。
> ⇒ **"触发用例缺口"卡应改写为"补判别性单测"**：这类错误的暴露条件可精确刻画（**轴位置、角度是否对称、缩放值、FromA1 姿态**四个变量），一条 20 行单测即可永久钉死，比新造几何用例便宜得多。
> **★ 本轮最值钱的手法**：**判别性必须实证，不能靠推导自证** —— 修复后把公式**临时改回旧形**跑一遍，确认新测试真的 FAILED 再还原（本轮两个新测试双双 FAILED，其余四个仍 ok）。
> **★ 第二条手法**：**清点同族站点时不要按变量名 grep** —— `fillet/hbuilder.rs` 用 `edp.range` 做变量名，直接漏过追加 32 的清点；改用名无关模式 `is_infinite_value\([A-Za-z_0-9]*\.range\[`（**与坑 28 同族，同一处栽了第二次**）。
> （历史：追加 21–32 收尾态见下方存档块；追加 32 = `set_rotation` 转置修复 + `UpdateEdge` 首次收敛；追加 25 = `blend_simple_a1` 通过；追加 28 = `LocOpe_SplitShape` 全类；追加 30 = `Geom_OsculatingSurface`。）

> **（存档）一句话（追加 32 收尾态，2026-09-13 —— 已被追加 33 取代）**：**工作口径仍是"先译后调"**（先 1:1 译全 body、把重复载体收敛到真身，代码基本译完再回头调试）。
> 本轮（追加 29–32）**四批、纯翻译为主**：`Approx_SameParameter`（992 行）· `UIso/VIso` 重复收敛（删 9 处 + 修 canonical 真身的基本面型臂缺陷）· `GeomLib_CheckCurveOnSurface`（774 行）· `Geom_OsculatingSurface`（839 行）+ `Geom_OffsetSurfaceUtils` ⇒ **canonical `UIso/VIso` 的 `Offset` 臂打通** · 曲面导数层补全（`BSplSLib::RationalDerivative`/`D0–DN`/Bezier `EvalD0–DN`，修掉 `bspline_surface_dn` 越界与 `BezierSurface::point_at` 有理求值缺陷）· `UpdateEdge` 家族收敛成唯一真身 · `Surface3::dn` 余留五面型补全。
> **★ 最重要的产出是修掉一个真正的 kernel 缺陷**：**`gp_Trsf::SetRotation` 的 `loc` 用了转置矩阵**（OCCT `gp_XYZ::Multiply(const gp_Mat&)` = `<me> = theMatrix * <me>`，`gp_XYZ.hxx` L308-309）⇒ **任何「轴不过原点」的旋转此前都错**；由 `bcut_simple_g6` 回归（面积 41187.4 → 61953.63，拓扑全过）逼出，见 §0.0h 的定位链。
> **门槛（当前基线）**：**427/0/0 · 712/0 · 36/36 · 26/26 · 76/76 · 1/1**；**八网格 8/8**（375/378/379/373/12/102/83/110）；**15 个域网格的真实断言数与基线逐项相同**。
> **★ 方法学（连续四轮同一结论，已升为最高优先队列项）**：追加 29–32 的**全部**修复（含内核几何修正）在 15 个域网格上**零可见翻转** ⇒ **域网格对内核几何/容差/区间类修正近乎不敏感，唯一的探测器是八网格**；本轮 `g6` 的覆盖是**回归逼出来的**、不是用例设计出来的。
> **★ 本轮最值钱的手法**：**数值不等价时先看"自洽对"** —— 差分导数与 `point_at` 按构造自洽，故「旧 `point_at`+差分 = OCCT 参考」而「新 `point_at`+差分 ≠ 参考」可**立刻**判定新体有误（无需读 358 行代码）。
> （历史：追加 21 收尾态见下方存档块；追加 25 = `blend_simple_a1` 通过；追加 28 = `LocOpe_SplitShape` 全类；追加 30 = `Geom_OsculatingSurface`。）

> **（存档）一句话（追加 21 收尾态）**：**工作口径仍是"先译后调"**——先把缺失的 body 逐一 1:1 译全、把重复/过期载体收敛到真身，等代码基本译完再回头调试与修测试（用户明确指令）。
> 本轮（追加 19–21）共 **9 批**：geomplate ProjLib 三处接线 · 池外 `BRep_Tool::Curve` · `GeomLib::BuildCurve3d` 家族 · 池外读取续链 · blend 侧两处字面翻译 ·
> **`BRepTools_Modifier` 家族翻译 + feat/offset 两域消费者接线** · **`Geom_Surface::UIso/VIso` 唯一真身**。
> **★ 已拿到域内第一个真实断言通过**：`draft_angle` **49/49 → 50 passed / 48 failed**（`draft_angle_b3` 通过，库内 panic 28 → 27）；`feat_featlf` a3 深入一层。
> **off-gate 效果汇总**：`offset_shape_type_i` **池外 panic 清零**（a3/a4/d2/d3 → 测试断言）· `blend_simple` q4/q7 深入两层并与 a3/a4 汇合 · `offset_shape_type_a` a4 推进（⚠ 该例失败点**不稳定**）· `blend_simple` a2/p8/p9 离开 GAP panic。
> 六门槛与八网格**全程零回归**（八网格始终 8/8，重编 exe 后实测）。
> **★ 三条方法学硬教训**（坑 28/31/33）：**同一 helper 行号会掩盖层推进（必须看 backtrace）**、**不要从截断清单数数**、**失败地图默认不稳定（HashMap 哈希序）——连跑 ≥3 次报集合**；
> 另有两条"旧笔记会过期"（坑 29/30）：**架构难点要回查 OCCT 基类**、**`GetType()` 常量返回决定分支**。
> （前三轮：**追加 18** = 翻译补全轮；**追加 17** = D3 结案 + `featrf_a1` init 首次通过；**追加 16** = 0a 收尾 / a1 拓扑全等 / TKOffset 定界。）

## 0. 新 session 一句话提示词（直接粘贴 —— 追加 39 收尾态，2026-09-14）

> 读 `rcad/docs/e3w-handover-tkfeat-fillet-offset.md`（本交接：门槛实测值 / 提交链 / 三域队列 / 配方 / 坑清单；
> **先读顶部"追加 39 收尾态"一行 + §0.0o**）与 `rcad/docs/tkfeat-fillet-offset-port-plan.md` §E3-W 追加 11–39
> （权威脉络，**追加 39 是当前状态**）；
> 先 `cd rcad` 跑 6 条门槛确认 **450/0/0 · 737/0 · 36/36 · 26/26 · 76/76 · 1/1**，
> 再 `cd /c/Users/lilu/works/rcad-pro && cargo test --no-run -p occt-generated-tests`（**重编 exe，否则八网格会拿旧产物误判**）
> 后 `bash output/run_eight_grids.sh` 确认**八网格 8/8**（375/378/379/373/12/102/83/110）；
> 然后**按 §0.0o 的队列继续"先译后调"**——**先只做 1:1 翻译/接线**（逐行语句对照 + OCCT 行号锚点 +
> 禁载体/禁等价替换/禁凑结果 + 架构差异先消灭再对齐），**代码基本译完前不针对性修测试数值**；
> **涉及布尔层的代码直接复用已对齐的 `bop/**`（TKBO）实现、不要重写第二份**（`BRepAlgoAPI_*` 包装层按 OCCT 形式接线到既有真身）；
> **立卡前先按 OCCT 函数名 grep 全库**（`panic!("GAP…")`/`unimplemented!` 的文案会过期，多例真身其实早已在库）；
> **失败层深度自检必须看 backtrace 而不是失败地图行号**（追加 19 实测：同一 `topods.rs:1799` 背后换了调用帧，行号完全看不出推进 —— 见坑 28）；
> **★ 改内核几何/容差/区间类代码必须跑八网格**：追加 29–39 证明 16 个域网格对这类修正**近乎不敏感**，
> **唯一的探测器是八网格**（追加 32 的 `bcut_simple_g6` 就是靠它抓到的；回归时先 `git stash -- <可疑目录>` 二分）；
> ⚠ 但追加 33 给出了**边界**：**唯一调用点落在退化配置上**时，真实内核缺陷连八网格也抓不到（第六轮"零可见翻转"即此）；
> ⚠ 追加 39 再给一条：**阈值比较类修正若既有用例都不在阈值边缘，八网格也照不到**（9 处真实阈值修正零可见）——
> 唯一兜底是**判别性单测 + 扰动实证**。
> **★ 数值不等价时先看"自洽对"**：差分导数与 `point_at` 按构造自洽 ⇒ 「旧 `point_at`+差分 = OCCT 参考」而
> 「新 `point_at`+差分 ≠ 参考」即可**立刻**判定**新体**有误，无需通读代码（追加 32 定位内核缺陷的关键手法）；
> **★ 手写矩阵约定必查左右乘**：根因是 `gp_XYZ::Multiply(gp_Mat)` 的语义（`<me> = theMatrix * <me>`）被写成转置；
> 这类错误在"轴过原点"/"对称角"/"scale = 1"/"FromA1 单位姿态"下**完全不可见**，写/改任何 `Trsf`/`Mat` 助手时逐句对 OCCT
> （追加 33 又抓到**两处**同族：`set_displacement` 两处转置乘、`multiplied` 漏乘 `scale`）；
> **★ 判别性单测必须实证**（追加 33 新立，追加 34 升为硬规矩）：修完把公式**临时改回旧形**跑一遍，确认新测试真的 FAILED 再还原 ——
> 「我推导过它一定会失败」不算数：追加 33 手工推导先错过一次，**追加 34 更狠 —— `transform_vec`/`transform_dir` 的第一版测试
> 在扰动下依然通过**（归一化抵消了正缩放，只有 `scale = ±3` 里的**负值**经末尾取反才暴露），**不扰动就得不到真判别性**；
> **★ 别把"注释说缺"当成"真缺"**（追加 34 第三次应验）：`Trimmed` 的 `EvalDN` 一直在 OCCT 里（`Geom_RectangularTrimmedSurface.cxx` L419-429），
> 记档却写了"无对应" —— 动手前先 grep OCCT 真身，见坑 7 / 追加 32 的 `EmptyCopy`；
> **★★ 动手前先验证"记档与任务前提"本身**（追加 35 升为硬规矩，**与"判别性必须实证"并列**；追加 36/37/39 各再证一次）：追加 37 三次翻车 ——
> `AdjustPeriodic` 的记档**位置**写错（在 `fillet/**` 不在 `bop/**`）；`bop/**` 那两处**根本不是同一个 OCCT 函数**（OCCT 有**四个**同名者，
> 判"是否同一个"必须逐站点判）；**canonical `epsilon_of` 的公式与锚点都错**（自称 `Precision.hxx::Epsilon`，而该头文件没有它，
> OCCT 的 `Epsilon(double)` 是 **nextafter 间隙**，`Standard_Real.hxx` L242-246）。**"canonical"两个字不等于"已对齐"** —— 必须回源；
> **追加 39 第三次应验**：修好的 canonical 仍是位递增体，在 x==±RealLast（OCCT 给 0）/±inf（给 -inf）/-0.0（给 +5e-324）边缘偏离，
> 而 `Intrv_Interval` 默认容差恰好吃这些值（OCCT 给 0.0，位递增给 +inf）——**先修 canonical 再收敛副本，顺序错了会污染消费者**；
> 派子代理时务必把"若前提不成立就停下来报告"写进任务书；
> **★ canonical/数学函数的边缘语义必须对照 C 标准库逐输入推演**（追加 39 新立）：x==y / ±inf / ±0.0 / NaN 四类逐个过
> （C nextafter：x==y 返回 y；`nextafter(inf, RealLast)` 下行到 RealLast），"有限输入上等价"不等于"函数等价"；
> **★ 代数自洽 ≠ 几何正确**（追加 35 新立）：组合律断言（`A.multiplied(B) == A(B(P))`）**抓不到矩阵转置**（对同一个错矩阵自洽），
> 只有断言**独立几何量**（Rodrigues / Householder / 绕心旋转 / 投影）的测试能抓 —— **两类断言都要留，不可互相替代**；
> **反过来，测试失败先怀疑测试**：追加 35 一例的期望值算错（误以为 `set_scale` 保留矩阵，OCCT 实际 `matrix.SetIdentity()`），
> **追加 36 一轮内出现三处**，**追加 39 边缘测试期望值又写错一次**（`epsilon_of(+inf)` 真值是 `-inf`，不是 `-RealLast`——
> IEEE 有限减无穷）⇒ 这条**升为常规动作**，且**期望值只能来自 OCCT 公式或独立几何定义逐项推演，不能来自"我以为"**；
> **★ 遇到 OCCT 自身的怪癖，钉 OCCT 行为而不是钉"数学上对的性质"**（追加 36 三例，最容易帮倒忙）：`gp_Trsf2d::Invert()` 对
> `gp_PntMirror` **不是**数学逆（只反 offset，3D `gp_Trsf::Invert` 一样，`gp_Trsf.cxx` L406-409）；`SetValues` 的 `M.Divide(s)`
> 经公开 API **不可观测**（后续 `Orthogonalize` 抵消）⇒ 如实记"不可观测"、不假装覆盖；`SetTranslationPart` 从 `gp_Rotation`
> 出发得到 `CompoundTrsf`（cxx L113-116 的 catch-all）。**若照"数学直觉"去写断言或去"修正"实现，反而会与 OCCT 分叉**；
> **★ 移植了 OCCT 自带测试 ≠ 覆盖完整**（追加 38 新立）：OCCT 自己的 `TEST(ElclibTests, AdjustPeriodic)` **五个用例都不覆盖它自己的
> `Epsilon(ULast)` 守卫**（禁掉守卫五个用例照样全过），而该守卫是那条公式的**唯一消费者**。移植值得做，**移植完必须自己扰动一遍**；
> 判别性输入的挑选同样要实证 —— 追加 38 第一版用"周期很小但不为零"**不判别**（wrap 恰好落回同值），只有**零长度区间**（除以 0 ⇒ NaN）才判别；
> **★ "多处复制了同一错误"要先逐份核验**（追加 38）：记档说错误公式被复制到多处，实测**每一份副本都是对的**、错误只在 canonical 体一处 ——
> 若不复核就会去"修"一批本来就对的代码（**好消息型推翻**，同样省下了白工）。**追加 39 的镜像教训**：另一些站点则**真是**复制病
> （3 份 `f64::EPSILON*v` 相对式 + bspl_lib 6 处 `Epsilon(1.)*x`）——**先逐份核验、再决定修不修**，两个方向都会翻车；
> **★ 清点同族站点不要按变量名 grep**（追加 33 新立）：用名无关模式
> `is_infinite_value\([A-Za-z_0-9]*\.range\[` 才不漏 —— `fillet/hbuilder.rs` 的 `edp.range` 就漏过了一整轮清点（与坑 28 同族）。
> **追加 39 第四次应验 + 升级为"终检兜底"**：按实现拼写 grep（`nextafter`）普查后，收敛完按函数名模式
> `fn (epsilon|standard_epsilon|standard_real_epsilon|occt_epsilon|epsilon_of)\(` 终检**又抓到 6 份漏网**（`next_up()` 拼写）——
> **普查和终检要用两个不同的名无关模式各跑一遍**；
> **★ 收敛重复体到 canonical 之前，先证明 canonical 本身逐输入等价于 OCCT**（追加 39 新立）：否则收敛会把 canonical 的缺陷
> 扇形污染到全部消费者（`Intrv_Interval` 默认容差 0→+inf 即险情实例）；
> 每 Edit 后 `cargo check -p rcad-algo`；**跑非门槛网格一律加 `timeout`**（§6 坑 16）；探针即用即清
> （提交前 `git diff | grep -c "+.*eprintln"` = 0，OCCT 侧探针用后 `git checkout --` 还原 + 重建 DLL）；
> 每批做完跑**六门槛 + 八网格 + 该域网格**，按坑 21 的**失败层深度**（不是通过数）自检，更新 port-plan §E3-W 追加，
> 并提交**两仓库**（rcad + 根仓库指针，rcad 推得上就推）。

### 0.0o 追加 39 收尾态（2026-09-14；**最新** —— 单线：epsilon 同族全库收敛（16 份）+ canonical 第三轮修正（OCCT 精确 nextafter 体）+ `bspl_lib` 6 处公式偏离 + 队列 3 对拍收口）

> **提交链**：规则不变 —— **rcad 顶尖 = 本交接文件所在提交**（`cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 追加 39 的成对 hash 按 §0.0h 同格式追加。（开工前两个 hash **必须成对**各自确认一遍。）

- **门槛与网格（全部在树实测）**：六门槛 **450/0/0 · 737/0 · 36/36 · 26/26 · 76/76 · 1/1**（kernel 736 → **737** = 1 个 canonical 边缘判别测试）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。合并后一次实测；探针 = 0；扰动实验已还原；净删 114 行（18 文件 +159/−273）。
- **★ 批次 1（canonical 第三轮修正）**：追加 37/38 修好的 `epsilon_of`（位递增）在四类边缘输入上偏离 OCCT（`Standard_Real.hxx` L242-248）：`x==±RealFirst/RealLast` ⇒ OCCT **0**（C nextafter 的 x==y→返回 y）、位递增给 ±inf；`x=±inf` ⇒ OCCT **-inf**（`nextafter(inf,RealLast)=RealLast`，`RealLast-inf=-inf`）、位递增给 NaN；`x=-0.0` ⇒ OCCT **+5e-324**（三元式取 `>=0.0` 臂）、位递增给 -5e-324；NaN 双方一致。**负载性发现**：`hlr/intrv::Interval::new()` 用 `epsilon(±RealFirst/Last)` 作默认容差（`Intrv_Interval.cxx` L39-40，OCCT 给 `(float)0.0`）——先收敛后修会把 HLR 区间容差从 0 污染成 +inf ⇒ **先修 canonical 再收敛是硬顺序**。落地：OCCT 原文三元式 + kernel `pub(crate) next_after`（C 语义全分支）；唯一消费者 `elclib_adjust_periodic` 守卫不受影响（inf 被其前 `is_infinite_value` 早退）。判别性测试 `epsilon_of_matches_nextafter_at_the_range_edges`（6 断言）+ 扰动实证（改回位递增 ⇒ 只有它 FAILED）。**期望值写错过一次**（`epsilon_of(+inf)` 真值 `-inf`，非 `-RealLast`）——"测试失败先怀疑测试"再 +1。
- **★ 批次 2（16 份本地体全库收敛）**：队列点名 4 份（`hlr/intrv`、`fillet/chfi3d_perform_elspine`、`geomalgo/geom_int_line_constructor`、`kernel/geom/mod.rs`——后两份即队列 2 的判定：**两份 `standard_epsilon` 就是同一个 OCCT 表达式**（`GeomInt_LineTool.cxx` L332/386 的 `Epsilon(firstp/lastp)`；rcad 侧 `_included=true` 是"IntPatchLine 存闭区间"架构注记、分支静态死））+ 名无关终检兜底抓到 6 份漏网（`next_up`/位递增拼写）+ 普查扩大 6 份。全部 `use ...::epsilon_of as <原名>` re-export 模式（调用点零改动）。**每站点回源核验**：`LocOpe_WiresOnShape.cxx` L1147/1234+ · `math_FunctionSetRoot.cxx` L873 · `ShapeAnalysis_TransferParametersProj.cxx` L333/356 · `Draft_Modification_1.cxx` L2198/2206 · `BRepLib_FindSurface.cxx` L499/548。**3 份公式真偏离**（自称 OCCT Epsilon、实为 `f64::EPSILON*v` 相对式）：`loc_ope_wires_on_shape_b`（10+ 活跃调用点）、`transfer_parameters_proj`（~25% 阈值差）、`draft_modification_1_b`（`Epsilon(2π)` 处 ~1.57 倍）。判无关：`color.rs::epsilon()`、`int_conic_conic_circ_circ` 的 `next_after`（取值用途）、`fclass2d_topol::safe_increment`（OCCT safeIncrement）。
- **★ 批次 3（`bspl_lib.rs` 6 处 `Epsilon(1.)*x` 公式偏离）**：无参 `epsilon()` 助手（= `f64::EPSILON`）被当相对式用，**仅 x 为 2 的幂时与 OCCT `Epsilon(x)` 相等**。逐站点修复：L271→OCCT L263 `Epsilon(min(|K|,|U|))`；L1743→L2125 与 L2509→L1911 `max(Tolerance, Epsilon(u/au))`（并删非 OCCT 的 `.abs()`）；L2719/2731→L614/628 三项和；L2863→L789 `Epsilon(|Knots|)` + nextafter 拼写统一 kernel `next_after`。助手删除（Rule 4）。
- **★ 批次 4（队列 3 对拍收口）**：16 域网格逐格同基线（`feat_featrevol` 46=1真+45占位等逐格对上）；`feat_featlf` a3 失败层 = `brep_sweep/tool_rehost.rs:1117` `Standard_Failure: Courbes non jointives`（忠实 OCCT 异常 = 上游几何分歧），**不在**追加 38 批次 3 的爆炸半径链上；`loc_ope_split_drafts*.rs` **无直测**（0 `#[test]`）——记档。
- **★ 方法学（第十二轮"零可见翻转"）**：9 处真实阈值修正（最高 ~1.57 倍）在两套网格全部零可见（比较点都不在阈值边缘）；canonical 修正是**潜在险情排除**（险情路径在 HLR 域、网格照不到）⇒ 判别性单测是唯一验收手段（第十三次重申）。**两条新规矩**：① **普查与终检用两个不同名无关模式各跑一遍**（实现拼写 grep + 函数名 grep 兜底）；② **收敛前先证明 canonical 逐输入等价 OCCT**（x==y/±inf/±0.0/NaN 四类推演），否则收敛会把 canonical 缺陷扇形污染到全部消费者。
- **⏭ 下一轮队列（按序）**：
  1. **canonical `epsilon_of` 消费者阈值抽查**（追加 38 队列 5 顺延）：`Intf_InterferencePolygon2d`、`GeomInt_IntSS_1` 等分支判据的比较方向/操作数逐站点核验。
  2. **`loc_ope_split_drafts` 直测**（新立，便宜）：0 个 `#[test]`，链上 `adjust_periodic_pair`（eps 1e-7→0.0）与 bean_face 的行为变化目前只有域网格计数兜底。
  3. **`tkgeom_algo_gtests::...::endless_loop_prevention` 既有 panic**（追加 38 队列 4 顺延）：`GCPnts_TangentialDeflection is not available in rcad-kernel`（`base/extrema_ext_pc.rs:119`）。
  4. **`intrv::Interval::new()` 消费点判别性测试**（新立，2 行）：默认容差 `Epsilon(±RealFirst/Last)==0` 已由 canonical 边缘测试钉住，消费点本身无断言。
  5. **`REAL_FIRST/REAL_LAST` 常量收敛**（新立，低优）：`hlr/intrv/mod.rs` L64-67 与 `kernel/precision.rs` L92-98 重复（值相同）。
  6. 追加 38 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · 池外读取链复核 · blend 剩余十例 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · 过期锚点勘误）。
- **复测脚本**：`rcad/temp/run_gates.sh` · `rcad/temp/run_domain_grids.sh` · `rcad/temp/validate_round5.sh`。

### 0.0n 追加 38 收尾态（2026-09-14；**已被追加 39 取代 —— 见 §0.0o** —— 三线并行：OCCT gtest 移植 + `epsilon_of` 复制族收敛（前提推翻）+ `bean_face_intersector` 误标翻译）

> **提交链**：规则不变 —— **rcad 顶尖 = 本交接文件所在提交**（`cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 追加 38 的成对 hash 按 §0.0h 同格式追加。（开工前两个 hash **必须成对**各自确认一遍。）

- **门槛与网格（全部在树实测）**：六门槛 **450/0/0 · 736/0 · 36/36 · 26/26 · 76/76 · 1/1**（algo 451 → **450**、kernel 733 → **736**）；**八网格 8/8**；**16 个域网格与基线逐项相同**。合并后一次实测；探针 = 0；扰动实验全部还原。
- **★ 批次 1（主代理）：移植 OCCT 自己的 `AdjustPeriodic` 测试，并发现它的覆盖缺口**（`rcad-kernel/src/math/el.rs`）。
  OCCT `GTests/ElCLib_Test.cxx` **L83-124** 五个用例**逐字移植**（同实参同期望值，`Precision::Confusion()` ↔ `precision::CONFUSION`）+ 无限区间一例。此前 rcad 对 canonical 体**零测试**。
  **★ 缺口**：这五个用例**完全不覆盖 `Epsilon(ULast)` 守卫** —— 把守卫**整个禁用**，五个用例**照样全过**；而该守卫是 `epsilon_of`（追加 37 刚修的公式）的**唯一消费者** ⇒ 此前**零覆盖**。
  **补判别性测试**：**零长度区间**（`u_first == u_last`）是唯一能暴露守卫的输入 —— 无守卫时 `a_period == 0.0` 作除数 ⇒ **NaN**。
  **★ 并实证了一个"不判别"的反例**：第一版用"周期很小但不为零"（跨 binade）**在守卫禁用下仍通过**（wrap 恰好落回同值）⇒ 只有零长度判别。
- **★ 批次 2（子代理）：`epsilon_of` 复制族收敛 —— 前提被推翻**。
  **记档"同一错误公式被复制到多处"是错的**：HEAD 里**每份副本的公式本来都对**（四种写法数值相同）。错误公式**只在 canonical 体一处**（即追加 37 修掉的那个）⇒ 收敛是**严格数值无操作**。
  **仍落地的结构收敛**：canonical 扩为 `pub fn epsilon_of` · **删 6 份本地体**、调用点指向 kernel · `direct_polynomial_roots.rs` 的 `epsilon` 保留为 **re-export**（消费者在 `topalgo/**`，超出该子代理范围）。
  **★ canonical 修正的实际影响**：旧公式在 `10.0` 给 `2.22e-15`、OCCT 给 `1.78e-15`，而 `10.0` 正是 `adjust_u_periodic` 的比较值 ⇒ **一条分支判据上曾有 25% 的阈值误差**，现已消除。
  **另报 4 份清单外副本（未改）**：两份 `standard_epsilon`（公式正确）；两份 `epsilon` **偏离 OCCT**（见批次 4）。
- **★ 批次 3（子代理）：`bean_face_intersector` 的 `AdjustPeriodic` 误标 —— 确认为真，1:1 翻译**。
  OCCT `IntTools_BeanFaceIntersector.cxx` **L599-649**：L612 `aEps = Epsilon(aUPeriod)`、L614-620 调 #2、**L621-623 无条件** `solutionIsValid = true; bUCorrected = true; U = aNewU;`，**忽略返回值**（V 块 L638 同理）。旧 rcad 是自创重写（`1e-12*period` + `ceil` 夹取 + 自创 `ok` 门）。
  **三处实质行为变化**：① **移除 rcad 独有的拒绝** —— UV 窗口**窄于周期**的 U 周期面（**半圆柱这类修剪柱面**）旧码**丢掉**窗口外的交点、OCCT 保留（带平移 U）；② 守卫 eps `6.28e-12` → `Epsilon(period)` ≈ `8.88e-16`；③ 平移 `ceil` → `modf/trunc`（恰为周期整数倍且越界时**差一个周期**）。
  **Job 2**：`adjust_periodic_pair` 判为 **#3** 的 1:1 形态 ⇒ 定为 canonical（`pub(crate)`），**四处 eps `PCONFUSION`(1e-7) → `0.0`**（OCCT 用声明默认）—— **真实行为偏离**：越界 ≤1e-7 的参数此前**不平移**、现在平移一整周期（同点，UV 差 2π）。`face_make_curve.rs` 的 `adjust_periodic_uv` 收成三行 shim。
  **爆炸半径**：`GeomIntLineConstructor` → `GeomIntIntSS` → `feat/loc_ope_split_drafts*.rs`，**不在布尔网格路径上**。algo 451 → 450（删的就是它那个主体已消失的测试）。
- **★ 批次 4（主代理）：两处返回 `MIN_POSITIVE` 的 `epsilon` 修正**。
  `hlr/intrv/mod.rs:45` 与 `fillet/chfi3d_perform_elspine.rs:65/79`：OCCT 的 `nextafter(0,+Inf)` 是**最小正次正规数 5e-324**，二者却给 `MIN_POSITIVE`(2.2e-308) —— **差 4.5e16 倍**。已改为 `f64::from_bits(1)`（负向取负）。两个 helper 的文档都自称 OCCT `Epsilon`，属"自称对齐、实际偏离"同类。
- **★ 方法学（第十一轮"零可见翻转"）**：批次 3 **改变了接受集合**与下游 UV 值、批次 4 改了两个容差 helper 的零行为 —— 八网格与 16 个域网格**全部零可见**。
  **★ 两条新教训**：① **移植 OCCT 自带测试 ≠ 覆盖完整**（必须自己扰动）；② **"多处复制了同一错误"要先逐份核验**（本轮是**好消息型**推翻，没有代码被污染）。
- **⏭ 下一轮队列（按序）**：
  1. **两处 `epsilon` 收敛到 kernel canonical**（新立）：`hlr/intrv/mod.rs`、`fillet/chfi3d_perform_elspine.rs` —— 公式已对齐，但仍是两份本地实现。
  2. **两份 `standard_epsilon`**（`geomalgo/geom_int_line_constructor.rs:1165`、`kernel/src/geom/mod.rs:743`）：公式正确，**是否同一个 OCCT 表达式**待判。
  3. **`bean_face_intersector` 的行为变化需域网格之外的对拍**：爆炸半径在 `LocOpe_SplitDrafts` ⇒ 核对该链失败层深度（`feat_featlf`/`feat_featrf` + `loc_ope_split_drafts` 直测）。
  4. **`tkgeom_algo_gtests::geom_fill_corrected_frenet_tests::endless_loop_prevention` 既有 panic**（子代理报出，与本轮无关）：`GCPnts_TangentialDeflection is not available in rcad-kernel`（`base/extrema_ext_pc.rs:119`）。
  5. **canonical `epsilon_of` 的消费者阈值抽查**（可选）：`Intf_InterferencePolygon2d`、`GeomInt_IntSS_1` 等分支判据。
  6. 追加 37 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · 池外读取链复核 · blend 剩余十例 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · 过期锚点勘误）。
- **复测脚本**：`rcad/temp/run_gates.sh` · `rcad/temp/run_domain_grids.sh` · `rcad/temp/validate_round5.sh`。

### 0.0m 追加 37 收尾态（2026-09-14；**已被追加 38 取代 —— 见 §0.0n** —— 三线并行：`Ax3` 补回平行分支 + `Init` 的 `my_location` + `AdjustPeriodic` 收敛与 `epsilon_of` 修正）

> **提交链**：规则不变 —— **rcad 顶尖 = 本交接文件所在提交**（`cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 追加 37 的成对 hash 按 §0.0h 同格式追加。（开工前两个 hash **必须成对**各自确认一遍。）

- **门槛与网格（全部在树实测）**：六门槛 **451/0/0 · 733/0 · 36/36 · 26/26 · 76/76 · 1/1**（kernel 728 → **733**）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。三线 + 批次 4 **合并后一次实测**；探针 = 0；扰动实验全部还原；改动文件全部 ≤ 2000 行。
- **★ 批次 1（主代理）：`Ax3` 形态学补全 —— 本轮真正的翻译缺口**（`rcad-kernel/src/math/gp.rs`）。
  **缺口**：`Ax3::set_direction`（`gp_Ax3.hxx` **L506-534**）与 `set_x_direction`（**L540-570**）**各缺一个平行分支**。给定方向与主方向（或当前 X）**(反)平行**时 `theV ^ (vxdir ^ theV)` 退化为零，OCCT 改为**重贴轴标签**：
  `set_direction`：`aDot>0` ⇒ `vxdir ← old vydir, vydir ← old axis`（L511-515）；`aDot<0` ⇒ 仅 `vxdir ← old axis`（L516-519）。
  `set_x_direction`：`aDot>0` ⇒ `axis ← old vxdir, vydir = -vydir`（L545-549）；`aDot<0` ⇒ 仅 `axis ← old vxdir`（L550-553）。
  **rcad 旧版代入零向量 ⇒ 产出退化框架**（单位性/正交性全破坏）。判据 `1-|aDot| <= Precision::Angular()` ↔ rcad `precision::ANGULAR`。
  **同时**把 `direct()` 从存储字段改为 OCCT 的**推导式**（`gp_Ax3.hxx` **L208**：`(vxdir ^ vydir)·axis.Direction() > 0`），删除只被构造函数置 `true` 的 `sense` 字段（追加 33 审计点名的"陈旧模型"）。
  **4 个新测试（gp 17 → 21）+ 3 次扰动实证**（关平行分支 ⇒ 两个平行测试失败；`direct()` 恒真 ⇒ 间接框架测试失败），**每个扰动只杀自己的目标**。
- **★ 批次 2（子代理）：`Init` 的 `my_location`（追加 36 报出）—— 记档为真，已修**。
  OCCT `BRepLib_FindSurface.cxx` **L289** 把 `myLocation` 传作 index 重载**出参**，`BRep_Tool.cxx` **L528** 赋 `L = E.Location() * GC->Location()`，消费于 L296/L337/L620；miss 是 `L.Identity()`（L534-537）。rcad 丢弃出参 ⇒ `my_location` 恒 0，L673 的 `loc_j == self.my_location` 拿 pcurve 哈希比一个**从不更新**的值。
  **落地**：私有 index 重载返回 OCCT 组合后的 location，用**边包装 location** 表达（rcad 把 pcurve 存在边的局部框架）——**与孪生再宿主 `build_curves3d.rs:677-678` 的既有替换一致**，未新造规则；miss 返回 0。旁证：`fs.location()` 写进面的曲面放置（`make_face.rs:180`、`brep_feat_rib_slot_b.rs:146`），对应 `BRepLib_MakeFace.cxx:206` ⇒ 必须是**绝对 location**。
- **★ 批次 3（子代理）：`AdjustPeriodic` 收敛 —— 记档"位置"写错**。
  **记档说残留两份在 `bop/**`；实际在 `fillet/chfi2d_ana_fillet_algo.rs` 与 `fillet/chfi2d_builder.rs`**（追加 29/30 正文本来就写着 fillet，是后续引用把位置串错）；而 `bop/**` 那两处**不是同一个 OCCT 函数**。
  **OCCT 有四个同名 `AdjustPeriodic`**：**#1** `ElCLib::AdjustPeriodic`（`ElCLib.cxx` **L115-149**）· **#2** `GeomInt::AdjustPeriodic`（`GeomInt.cxx` **L21-48**，`theEps` 默认 **0.0**，返回 bool）· **#3** `GeomInt_LineConstructor.cxx` 文件静态包装（L737-816）· **#4** `GeomInt_LineTool.cxx:128`。判"是否同一个"**必须逐站点判**。
  **收敛**：#1 → kernel `elclib_adjust_periodic`（**删 4 份 / 改 7 个调用点**）；#2 → `geomalgo::geom_int_line_constructor::geom_int_adjust_periodic`（**删 3 份 / 改 7 个调用点**），并按各 OCCT caller 传的 eps 设定，其中 `face_make_curve.rs` **1e-9 → 0.0**（1e-9 不是任何 OCCT caller 传的值）。
  **四份 #1 的守卫差异（如实记档）**：kernel 用 `epsilon_of(u_last)`（OCCT 的参数 / 非 OCCT 的公式）· bop 用 `standard_epsilon(a_period)`（OCCT 的公式 / **参数错** ⇒ 溢出守卫对小周期不触发）· fillet ×2 与 `tool_rehost` 用 `|ulast|*DBL_EPSILON`（≈kernel）⇒ **收敛到 kernel 体**（守卫行为与 OCCT 一致），bop 体因参数错被否决。
  **未改并记档**：`bean_face_intersector.rs:1919` **标错**（自创 ceil 夹取 + `ok` 门，收敛会**改变管线控制流** ⇒ 需独立卡）；`adjust_periodic_uv` 是 #3，又被 `geomalgo` 的 `adjust_periodic_pair` 重复（后者传 `PCONFUSION`，OCCT 传 0.0）。
  **★ 并推翻一条文档断言**：kernel `el.rs` 称 canonical 体"含 `Epsilon(theULast)`"，而 `epsilon_of` 自称 `Precision.hxx::Epsilon`、公式 `Max(|v|,RealSmall())*RealEpsilon()` —— **两半都错**。
- **★ 批次 4（主代理）：canonical `epsilon_of` 1:1 修正**（`rcad-kernel/src/base/extrema_ext_elc.rs`）。
  按 OCCT 的 **nextafter 间隙**重写（`f64::from_bits(bits ± 1)`），锚点/文档一并改正；**1 个判别性测试**（`epsilon_of(1.5)` 必须是**一 ulp** = 2.220446049250313e-16，旧公式给 3.33e-16 即判别点）+ 扰动实证。
  **数值影响**：同一 binade 内旧式与 OCCT 最多差 **2 倍**；`1.0` 等 2 的幂处**完全相同**；只对**次正规**输入真正咬合。
- **★ 方法学（第十轮"零可见翻转"）**：内核形态学（`direct()` 存储→推导）、`Init` 的 location、`AdjustPeriodic` 收敛、`epsilon_of` 公式 —— 八网格与 16 个域网格**全部零可见**。
  **★ 本轮三次"回源核验救命"**（位置错 / 不同函数 / canonical 公式错）⇒ 与追加 35/36 合并，**"canonical / 记档必须逐条回源核验"已是常态动作**。
- **⏭ 下一轮队列（按序）**：
  1. **`epsilon_of` 复制族收敛**（新立）：仓库里还有 **~7 份本地 `epsilon_of`**（`brep_fill_sweep.rs:743`、`geom_int_int_ss_1.rs:348`、`intf_interference_polygon2d.rs:24`、`int_conic_conic_lin_ells.rs:30`、`law_bspline.rs:29`、`geom2d_convert/bspline_curve.rs:53` 等），**每份都可能复制了同一错误公式** ⇒ 逐份对 OCCT 锚点核验后再合并到 kernel 真身。
  2. **`bean_face_intersector.rs:1919` 的 `AdjustPeriodic` 误标**（新立）：收敛会改变管线控制流（当前 rcad 拒绝 OCCT 接受的解）⇒ 需独立评估。
  3. **`adjust_periodic_uv`(#3) 与其重复 `adjust_periodic_pair`**（新立，便宜）。
  4. **移植 OCCT `ElCLib_Test.cxx:83` 的 `TEST(ElclibTests, AdjustPeriodic)`（5 断言）**守护 canonical 体。
  5. **其余内核助手的判别性测试**：`gp.rs` 的 `Trsf` 助手卡已关（21 测试）；同一手法继续扫 `Lin::*`、`Ax2::from_direction` 的分支树。
  6. 追加 36 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · 池外读取链复核 · blend 剩余十例 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · 过期锚点勘误）。
- **复测脚本**：`rcad/temp/run_gates.sh` · `rcad/temp/run_domain_grids.sh` · `rcad/temp/validate_round5.sh`。

### 0.0l 追加 36 收尾态（2026-09-14；**已被追加 37 取代 —— 见 §0.0m** —— 三线并行：`gp_Trsf2d` 全类译完（抽模块）+ `gp.rs` 判别性单测 + 超长文件拆分与 location 勘误）

> **提交链**：规则不变 —— **rcad 顶尖 = 本交接文件所在提交**（`cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 追加 36 的成对 hash 按 §0.0h 同格式追加。（开工前两个 hash **必须成对**各自确认一遍。）

- **门槛与网格（全部在树实测）**：六门槛 **451/0/0 · 728/0 · 36/36 · 26/26 · 76/76 · 1/1**（algo 438 → **451**；kernel 723 → **728**）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。三线**合并后一次实测**；探针 = 0；扰动实验全部还原；**改动文件全部 ≤ 2000 行**。
- **★ 批次 1（主代理）：`gp.rs` 判别性单测补齐 + `eval_b.rs` 守卫政策勘误**。
  **5 个新测试（12 → 17）**：`vectorial_part_folds_in_the_scale` · `value_folds_in_the_scale_and_reads_the_location` · `is_negative_reads_the_scale_not_the_form` · `set_translation_part_follows_the_occt_form_transitions` · `set_scale_factor_follows_the_occt_form_transitions`（12 例表驱动）。
  **3 次扰动实证**：`is_negative` 改读 `form` ⇒ 只有它失败；`vectorial_part` 忘记乘 scale ⇒ 它 + 消费者 `transform_vec` 一起失败（正确归因）；`set_translation_part` 的 `_` 臂改成无条件 `CompoundTrsf` ⇒ 只有它失败 —— 判别性正落在「**Rotation + 零位移应保持 `gp_Rotation`**」（OCCT L264-280 的 `if (!loc_null)`）。
  **`eval_b.rs` 勘误**：该模块头注称 `Geom_UndefinedDerivative` 守卫"无 rcad 对应、故不复现"——**假的**。`BSplSLib::DN`（本模块再宿主的**叶子**，`BSplSLib.cxx` L1519）确实无阶数检查，但守卫属 `Geom_*` **包装层**，rcad **确实复现了**（`eval.rs` Trimmed 臂 / `extrusion_utils.rs` / `revolution_utils.rs` / `eval_c.rs` / `curve_dn.rs`）。**并记下**：OCCT 把守卫放在 eval-rep 短路的**两侧**（`Geom_BSplineSurface::EvalDN` L279-285 在**前**、`Geom_BezierSurface::EvalDN` L1674-1677 在**后**）⇒ **不能把共享守卫提到两者之上**。
- **★ 批次 2（子代理）：`gp_Trsf2d` 全类译完 + 抽模块**。
  译全余下 **10 函数 + `Power`**（全部逐语句 + OCCT 行号）：`SetTransformation(Ax2d,Ax2d)`（cxx L48-70）· `SetTransformation(Ax2d)`（L72-84）· `SetTranslationPart`（L86-117）· `SetScaleFactor`（L120-196，3 外层 × 6 内层标签迁移）· `VectorialPart`（L198-214，含对角分支）· `RotationPart`（L216-219）· `Invert`（L221-251，4 分支含 2 raise）· `Power`（L384-548）· `PreMultiply`（L550-671，15 分支）· `SetValues`（L676-710）· `Orthogonalize`（L721-749）；新助手 `mat2_transpose`（`gp_Mat2d.hxx` L394-399）、`gp_xy_normalize`（`gp_XY.hxx` L410-417）。
  因加测试会越 2000 行，**整个载体字节级抽到新模块** `shhealing/shape_extend/trsf2d.rs`（1631 行），`composite_surface.rs` 1866 → **746**；`mod.rs` 声明并 re-export，**旧路径仍可解析**。
  **13 新测试（`trsf2d` 模块 22）+ 12 次判别性扰动实验**（注入错形 ⇒ 目标 FAILED ⇒ 还原，文件 md5 复原）。
  **★ 三条如实记档的 OCCT 怪癖**：**（a）`Invert()` 对 `gp_PntMirror` 不是数学逆**（cxx L231-234 只反 offset ⇒ `t_inv(t(p)) == p - 4c`；**3D `gp_Trsf::Invert` 完全相同**，`gp_Trsf.cxx` L406-409）⇒ 翻译与之对齐、**测试钉 OCCT 行为**；**（b）`SetValues` 的 `M.Divide(s)`（cxx L701）经公开 API 不可观测**（`s = sqrt(|det|)>0`，随后 `Orthogonalize` 归一化精确抵消）⇒ **如实记"不可观测"，不假装覆盖**；**（c）`SetTranslationPart` 从 `gp_Rotation` 出发得到 `CompoundTrsf`**（cxx L113-116 catch-all）。
  另记：`VectorialPart` 对角专用分支与 `Orthogonalize` 行 pass 经公开 setter 不可达/不可单独观测（已注明）。
- **★ 批次 3（子代理）：超长文件拆分 + location 比对勘误（记档为真）**。
  **`brep_fill_sweep.rs` 2072 → 1822 行**：把 `bspline_surface_iso_tests`（追加 35 新增）与相邻 `offset_surface_iso_tests` 整体移到兄弟文件（71 + 179 行），按既有 `#[path]` 约定声明；**测试名与全路径逐字节不变**、**测试数 438 → 438 不变**、生产代码零改动。
  **location 勘误为真**：OCCT `BRep_Tool.cxx` **L350** `const TopLoc_Location loc = L.Predivided(E.Location());`，存进表示的也是同一 predivided 值（`BRep_Builder.cxx` L645/L692），`IsCurveOnSurface` 按**相等**比（`BRep_CurveOnSurface.cxx` L59-63）⇒ 调用方须先把**原始 `L`** 除以边 location。rcad 的 `_stored` 却拿 `key.1`（predivided 哈希）比**原始** `the_location`（**两个量不同**）；同族 `shhealing/shape_analysis/edge.rs:131` 才正确。已对齐（位置亦照 OCCT L350 提到 rep 循环之前）。**行为影响**：identity location 边不变；**loc 边旧代码永远匹配不上**、改后能匹配 ⇒ 严格更接近 OCCT。
  **⚠ 另报一处更大的相邻缺口（未改，见队列第 1 项）**：同文件 `Init`（L651）**丢弃** index 重载的 location 出参、从不赋 `self.my_location`（恒 0），而 OCCT 设 `myLocation = E.Location() * GC->Location()`（`BRep_Tool.cxx` L528）。
- **★ 方法学（第九轮"零可见翻转"）**：本轮全部改动在八网格与 16 个域网格上**零可见**。
  **★ 新形态问题：一轮内三处"测试写错"**（我两处 + 子代理一处），**全部"代码对、期望错"** ⇒ **"测试失败先怀疑测试"升为常规动作**，且**期望值只能来自 OCCT 公式或独立几何定义**。
- **⏭ 下一轮队列（按序）**：
  1. **`brep_lib_find_surface.rs` 的 `Init` 丢弃 location 出参 + `my_location` 恒 0**（新立，**比已修的那处更大**）：OCCT `BRep_Tool.cxx` L528 / `BRepLib_FindSurface.cxx` L289/L337。
  2. **`gp.rs` 的 `Trsf` 助手卡关闭**（15 助手 × 17 测试已全覆盖）⇒ 把同一手法（写测试 → 扰动 → 确认 FAILED → 还原 + 独立几何量断言）用到**其他无测试的内核助手**：`Ax2::from_direction` 的分支树、`Lin::*`、`Ax3::set_x_direction` 等。
  3. **`Trsf2d` 的三处"不可观测语句"如实留档**（`SetValues` 的 `M.Divide(s)`、`VectorialPart` 对角分支、`Orthogonalize` 行 pass）：**不要假装覆盖**。
  4. **`mySn` 卡证伪 ⇒ 关闭**：`mySn` 只存在于 `LocOpe_Gluer`（全 TKFeat grep 唯一）且 rcad 侧已完整；下一墙在 `feat/loc_ope_glued_shape.rs:193`（OCCT `LocOpe_GluedShape.cxx` L106-110 的"边面祖先数 != 2"），属**调试/状态类**。
  5. **`geom/eval_b.rs` 勘误卡关闭**（本轮已改完）。
  6. 追加 33 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · 过期锚点勘误）。
- **复测脚本**：`rcad/temp/run_gates.sh` · `rcad/temp/run_domain_grids.sh` · `rcad/temp/validate_round5.sh`。

### 0.0k 追加 35 收尾态（2026-09-14；**已被追加 36 取代 —— 见 §0.0l** —— 三线并行、全翻译类：`gp_Trsf2d` 译完 + BSpline iso 收敛到精确体 + `CurveOnPlane` 回退）

> **提交链**：规则不变 —— **rcad 顶尖 = 本交接文件所在提交**（`cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 追加 35 的成对 hash 按 §0.0h 同格式追加。（开工前两个 hash **必须成对**各自确认一遍。）

- **门槛与网格（全部在树实测）**：六门槛 **438/0/0 · 723/0 · 36/36 · 26/26 · 76/76 · 1/1**（rcad-algo 430 → **438**；kernel 不变）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与基线逐项相同**。三线**合并后一次实测**；探针 = 0。
- **★ 批次 1（主代理）：`gp_Trsf2d` 载体译完**（`shhealing/shape_extend/composite_surface.rs`）。
  **① 记档纠错**：队列写"缺 `Rotation/PntMirror/Ax1Mirror/Ax2Mirror` **四个**"—— **错**。`gp_Trsf2d` 的形态集无 `Ax2Mirror`（双面对称只在 3D；2D 是点对称 `PntMirror` 与线对称 `Ax1Mirror`），**只缺三个**。
  **② 新译**（逐语句 + OCCT 行号入注释）：`set_rotation(P,Ang)`（hxx **L246-257**）、`set_mirror_pnt(P)`（hxx **L259-266**）、`set_mirror_ax2d`（cxx **L31-46**）、`value(r,c)`（hxx **L301-314**）。
  **③ `multiplied` 重写为 `gp_Trsf2d::Multiply` 全 15 分支 1:1**（cxx **L253-384**）：旧体只实现 4 条、其余塌成 `CompoundTrsf`；现在分支顺序、每支的 `Tloc`/`scale`/`matrix` 组合与 **form 标签**全跟随 C++（L282 双点对称 ⇒ `Translation`、L288 双线对称 ⇒ `Rotation`、L323 `shape = T.shape`）。两个矩阵助手：`mat2_mul_xy` = `gp_XY::Multiply(gp_Mat2d)` = **矩阵左乘**（`gp_XY.hxx` **L401-406**）；`mat2_mul` = `gp_Mat2d::Multiply` = `this*theOther`（`gp_Mat2d.hxx` **L326-334**）。
  **④ 6 个新测试 + 2 次扰动实证**（扰动 ⇒ 只有目标测试失败）。**组合律测试**（`A.multiplied(B)` == `A(B(P))`，5 形态 × 25 组合 × 2 点）一条断言覆盖 15 分支的代数。
  **★ 测试侧两条教训**：**（a）测试失败先怀疑测试** —— 我的 `value` 期望值算错（误以为 `set_scale` 保留矩阵，OCCT hxx L270-277 实际 `matrix.SetIdentity()`），**代码是对的**；**（b）代数自洽 ≠ 几何正确** —— 转置 `set_rotation` 矩阵后组合律测试**不失败**（对同一错矩阵自洽），只有断言独立几何量的测试失败 ⇒ **两类断言都要留**。
- **★★ 批次 2（子代理）：BSpline `UIso/VIso` 收敛 —— 推翻任务前提**。
  **前提只有一半为真**：任务书假设"canonical 体精确、两份副本近似"，实测**canonical 体自己就是第三份 de Boor 近似**。子代理因任务书禁止"收敛到另一份近似"而**先读 OCCT 再动手**，否则会把两份近似合并成一份、**看起来收敛实际没修**。
  **落地**：canonical 臂改接**精确** `bspl_slib_iso`（`BSplSLib::Iso` 1:1，`BSplSLib.cxx` **L1617-1740**），逐语句对齐 `Geom_BSplineSurface::UIso`（`Geom_BSplineSurface_1.cxx` **L598-630**）/`VIso`（**L775-807**），含 `Rational()` 选权、`NoMults()` 扁平节点、**逆向**的节点/次数/周期传递；删掉两份副本（`brep_fill_nsections.rs:1027`、`geomfill/nsections.rs:280`），5 个调用点改接，三文件 `de_boor_homo` **grep = 0**。
  **★ 量出的真实行为差异（此前一直错着）**：① **域外参数**（V=1.4 于 [0,1] 面）旧体夹到边界极 `(0,1,2)`、精确体**外推** `(0,1.4,2.8)`，**极差 0.894**；② **权重全 2.0 的均匀有理面**旧判据 `any weight != 1.0` **不是** OCCT 的 `Rational()`（权重留 2.0 vs 写 1.0）；③ **周期面**旧体无任何周期处理且硬编码 `is_periodic: false`（潜伏：树里无周期面）。域内非周期两者数学等价（三次 max 极差 3.1e-16）。
- **★ 批次 3（子代理）：`CurveOnPlane` 回退 —— 记档为真，已补齐**。
  **字形容易误读**：`CurveOnPlane` 是**函数**（`BRep_Tool::CurveOnPlane`，`BRep_Tool.cxx` **L379-450**），**不是**表示类 —— `BRep_CurveOnPlane` 类**在本 OCCT 树里不存在**。它是 `BRep_Tool::CurveOnSurface(E,S,L,First,Last,theIsStored)`（L327-373）的尾巴，由 **L367-372 `// Curve is not found. Try projection on plane`** 进入（投影得到、**不存**）。
  **缺口在 topalgo 再宿主**：`topalgo/brep_lib_find_surface.rs:118` 未命中返回 `None`，OCCT 则投影返回；唯一调用者 `is_2d_connected`（镜像 `BRepLib_FindSurface::Is2DConnected`，L76-99）把 `None` 变 `false` ⇒ **与 OCCT 走不同分支**。已补回退 + 新 `brep_tool_curve_on_plane`（1:1 `BRep_Tool.cxx` L379-450：`Geom_RectangularTrimmedSurface` 基面解包、3D 曲线存在性、`aCurveLocation.Predivided(L)` + 缩放重算参数、复用既有投影）。
  **已核对"不是缺口"**：kernel `BRepTool::curve_on_surface`（`topods.rs:2408`）**早有**回退（`topods.rs:2509-2511`）；`build_curves3d.rs:636` 的 `..._index` **正确地没有**回退（OCCT index 重载也没有）。
  **⚠ 另报相邻不一致（未改，见队列第 2 项）**：topalgo `_stored` 拿**原始** `the_location` 比对存哈希，而 OCCT（`BRep_Tool.cxx:345`）与同族 `shhealing/shape_analysis/edge.rs:124` 比 `L.Predivided(E.Location())`。
- **★ 方法学（第八轮"零可见翻转"，本轮连"分支走向改变"也不可见）**：批次 3 把一条分支从"返回 `None`"改成"投影返回"、批次 2 换掉整条求值公式 —— 八网格与 16 个域网格**依然全部零变化**。
- **⏭ 下一轮队列（按序）**：
  1. **`brep_fill_sweep.rs` 拆出测试模块**（新立，合规类）：本轮 +75 行测试后该文件 **2072 行**，越过 2000 行红线。
  2. **`brep_lib_find_surface.rs` 的 `_stored` location 口径勘误**（新立，便宜）：与 OCCT 及同族不一致。
  3. **`gp.rs` 判别性单测继续**（最高优先，第八次）：仍零测试的 `vectorial_part` / `value` / `is_negative` / `set_translation_part` / `set_scale_factor`；**标准动作 = 写测试 → 扰动 → 确认 FAILED → 还原**，且**同时保留"代数自洽"与"独立几何量"两类断言**。
  4. **`Trsf2d` 余下的 2D 批**（按需）：`Invert` / `Power` / `PreMultiply` / `Orthogonalize` / `SetValues` / `SetTransformation(Ax2d)` / `SetScaleFactor` / `SetTranslationPart`（`gp_Trsf2d.cxx` 还剩 L86-252、L384-749）。
  5. **`geom/eval_b.rs` 守卫政策勘误**（追加 34 立卡，仍在）。
  6. **`mySn` 的构造**（TKFeat，`feat_featrf` 直接前墙）。
  7. 追加 33 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · 过期锚点勘误）。
- **复测脚本**：`rcad/temp/run_gates.sh` · `rcad/temp/run_domain_grids.sh` · `rcad/temp/validate_round5.sh`。

### 0.0j 追加 34 收尾态（2026-09-14；**已被追加 35 取代 —— 见 §0.0k** —— 三线并行：`gp.rs` 判别性单测（4→12）+ `Trsf2d` 审计（再抓 2 缺陷）+ `Trimmed`/`dn_at` 补全）

> **提交链**：规则不变 —— **rcad 顶尖 = 本交接文件所在提交**（`cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 追加 34 的成对 hash 按 §0.0h 同格式追加。（开工前两个 hash **必须成对**各自确认一遍。）

- **门槛与网格（全部在树实测）**：六门槛 **430/0/0 · 723/0 · 36/36 · 26/26 · 76/76 · 1/1**（rcad-algo 427 → **430**；rcad-kernel 714 → **723**；增量全是新测试）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与追加 32 基线逐项相同**。三线**合并后一次实测**；探针 = 0，扰动实验全部还原。
- **★ 批次 1（主代理）：`gp.rs` 判别性单测补齐（队列第 1 项）+ 新译一个缺失重载（队列第 3 项）**。
  **测试 4 → 12**。新增：`transform_vec_folds_in_scale_while_transform_dir_does_not` · `apply_scales_the_matrix_term_only` · `invert_round_trips_scale_rotation_and_translation` · `set_transformation_ax3_maps_its_frame_to_the_standard_frame` · `to_daffine3_agrees_with_apply` · `set_transformation_from_to_changes_between_frames`（另有追加 33 的两个内核缺陷测试）。
  **★ 5 次扰动实验实证判别性**（改回旧形 ⇒ 目标测试 FAILED ⇒ 还原；每次**只有**目标测试失败）。
  **★★ 本轮最值钱的意外（P1）**：`transform_vec`（用 `VectorialPart()`，含 scale）与 `transform_dir`（用 `HVectorialPart()`，原始 matrix）这对**第一版测试在扰动下仍然通过 —— 根本不具判别性**。原因：`transform_dir` 末尾的**归一化抵消了正缩放**，只有**负缩放**经 `scale < 0` 的反转才暴露（错体先被 `3R` 缩放、归一化后又被末尾取反 = 双重取反抵消）。测试改为**同时跑 `scale = ±3`** 后 P1 才失败。⇒ **不做扰动实验就会留下"看似覆盖配对、实则恒真"的假回归测试。**
  **新译 `Trsf::set_transformation_from_to`**（OCCT `gp_Trsf::SetTransformation(FromA1,ToA2)`，`gp_Trsf.cxx` **L172-194**）—— **与 `set_displacement` 是两条不同路径**：前者**换坐标**，后者**搬框架**；OCCT 明说两者 vectorial part **互为逆**（`gp_Trsf.hxx` **L136-137**）。
  **三处逐字差异（写反就错）**：① 这里 `gp_Mat MA1` **不 Transpose**（`SetDisplacement` 里有）；② `MA1loc.Multiply(matrix)` 乘的是 **`matrix` 而不是 `MA1`**；③ `gp_XYZ::Multiply(const gp_Mat&)` = `<me> = theMatrix * <me>`（`gp_XYZ.hxx` **L308-309**），两处都是**矩阵左乘**。测试用**独立推导**（同一世界点投影到两个框架）断言 + `FromA1 == ToA2 ⇒ 恒等` 的退化检查；扰动 P5（改用 `SetDisplacement` 的转置约定）⇒ 目标测试 FAILED。
- **★ 批次 2（子代理）：`shhealing/shape_extend/composite_surface.rs` 本地 `Trsf2d` 审计（队列第 2 项）**。
  对象是**本地** `Trsf2d`（L758-891，`gp_Trsf2d` 的简化 GAP 载体，**不是** `rcad-kernel::math::gp`）；唯一生产/消费点是 `ShapeExtend_CompositeSurface::GlobalToLocalTransformation`（`ShapeExtend_CompositeSurface.cxx` L369-393，构造 `Trsf = Scale * Shift`）。
  **9 项审计：7 ALIGNED / 2 DEFECT**：
  1. **`multiplied` 的 catch-all**（`gp_Trsf2d.cxx` **L365-381**）：旧体 `tloc = t.loc * res.scale` 既**漏了矩阵因子**（OCCT `Tloc = matrix * T.loc`，左乘，`gp_XY.hxx` L401-406）又**整条 `matrix.Multiply(T.matrix)` 没做**。反例：`self={matrix=[[1,2],[3,5]],loc=0,scale=2}`、`t={matrix=[[2,0],[0,3]],loc=(1,1),scale=1}` ⇒ OCCT `(10,28)`，旧 rcad `(4,8)`。**与追加 33 的两个 gp 缺陷同族。**
  2. **`set_scale` 里一个 OCCT 没有的 raise**：`assert!(scale.abs() > f64::MIN_POSITIVE)` 出自 **3D** 的 `gp_Trsf.cxx` L164-165，**2D 类不 raise**（`gp_Trsf2d.hxx` L270-277），语句顺序也不同。**★ 而且可达**：唯一调用点的 `scalev = (v2-v1)/(…)` 在 V 跨度为零的 patch 上就是 **0** ⇒ OCCT 返回退化变换，rcad 直接 abort。**"自创的守卫"与"缺失的步骤"同为缺陷**（对准禁忌"不自创方法"）。
  **判 UNCERTAIN 且未改（如实记档）**：`Translation∘Translation`/`Translation∘Scale` 的 form 标签与 OCCT 特例不同（rcad 恒 `CompoundTrsf`），**无数值影响**；另记 `TrsfForm` 缺 `Rotation/PntMirror/Ax1Mirror/Ax2Mirror` 四个变体 ⇒ catch-all 保留 GAP 注释。3 个新测试 + 3 次扰动实验。
- **★ 批次 3（子代理）：`Trimmed` 的 `dn` 臂 + `dn_at` 全阶（队列第 4 项）**。
  **★ 先推翻记档**：`Surface3::dn` 的注释断称「`Trimmed` 包装**没有** OCCT `Geom_Surface::EvalDN` 对应」—— **假的**。OCCT `Geom_RectangularTrimmedSurface::EvalDN` 确实存在（**L419-429**，`final`），体就是「范围守卫 + `basisSurf->EvalDN(U,V,Nu,Nv)`」。**"注释说缺 ≠ 真缺"第三次应验。**
  已加 `Trimmed` 臂（含 OCCT 自己的 `Nu+Nv<1 || Nu<0 || Nv<0` throw）。**注意**：该守卫属各 `Geom_Surface::EvalDN`，**`GeomAdaptor_Surface::EvalDN` 本身没有** ⇒ `dn_at` 不加守卫（逐字对齐的结果，不是遗漏）。
  `dn_at`（`geom_adaptor_surface.rs:435`）从「只支持 `(1,0)`/`(0,1)`、其余 **panic**」扩到 **OCCT 全分支**（`GeomAdaptor_Surface.cxx` **L1697-1814**）：边界 snap → `mySurfaceType` switch（BSpline / extrusion / revolution / offset / 五个二次曲面 → `ElSLib::DN`；Bezier+Other 的 `break` → `mySurface->EvalDN` 尾）。**panic 清零**，只剩 rcad 独有的 Ruled/Coons/Pipe/TriBezier 保留 GAP（OCCT 无对应类）。**叶子全部复用既有真身**，未新写数学。
  `geom_adaptor_transformed_surface.rs:437` **判为无需改动**（OCCT 对应 `GeomAdaptor_TransformedSurface::EvalDN` L306-316 就是 identity 分支 + `transform_vec`，rcad 已逐行等价，自动继承扩宽覆盖）。3 个新测试，期望值由 OCCT 公式独立给出；柱面用例专门打**此前必然 panic 的阶数**。
  **⚠ 跨文件不一致（另报，未改，见队列第 2 项）**：`geom/eval_b.rs` L32-35 把「`Geom_UndefinedDerivative` 守卫无 rcad 对应」写成**政策**，但 `Surface3::dn` 的兄弟臂**已经都在 assert** ⇒ 政策文与现状矛盾。
- **★ 方法学（第七轮"零可见翻转"，本轮形态最极端）**：`dn_at` 从**一个 panic** 变成**真的算值**、`set_scale` 从**一次 abort** 变成**返回退化变换** —— 八网格与 16 个域网格**依然全部零变化**。前几轮尚可解释为"改了但没走到"，本轮是**把不可达路径变成可达**仍零翻转 ⇒ **两套网格合起来的覆盖率对本域进展几乎无分辨力，"判别性单测"是唯一验收手段**（第七次重申）。
- **⏭ 下一轮队列（按序）**：
  1. **★ 判别性单测继续**（最高优先，第七次）：`gp.rs` 仍零测试的 `vectorial_part` / `value` / `is_negative` / `set_translation_part` / `set_scale_factor` 及 `multiplied` 的特例分支。**标准动作 = 写测试 → 扰动实现 → 确认目标测试 FAILED → 还原。**
  2. **`geom/eval_b.rs` 守卫政策勘误**（新立，便宜）：与 `Surface3::dn` 兄弟臂现状矛盾，择一统一 + 两处互相引用。
  3. **`TrsfForm` 补四个缺失变体**（新立）：`Trsf2d` 的 catch-all 因此只有 General 分支；要覆盖 `Rotation/PntMirror/Ax1Mirror/Ax2Mirror` 需先补枚举。
  4. **`Trsf2d` 补 `matrix` 公开 setter / `SetValues`**：子代理的测试不得不直接构造结构体字面量 ⇒ 2D 载体尚未译完（TKMath 2D 批）。
  5. **`mySn` 的构造**（TKFeat，`feat_featrf` 直接前墙）。
  6. 追加 33 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · BSpline VIso 两处内嵌副本 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · `brep_tool_curve_on_surface` 缺 `CurveOnPlane` 回退 · 过期锚点勘误）。
- **复测脚本**：`rcad/temp/run_gates.sh` · `rcad/temp/run_domain_grids.sh` · `rcad/temp/validate_round5.sh`（门槛 + 重编 + 八网格 + 域网格一次跑完）。

### 0.0i 追加 33 收尾态（2026-09-14；**已被追加 34 取代 —— 见 §0.0j** —— 队列表末两项：`UpdateEdge` 全库收敛完毕 + `gp.rs` 矩阵约定全量审计（再挖出 2 个内核缺陷））

> **提交链**：规则不变 —— **rcad 顶尖 = 本交接文件所在提交**（`cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 追加 33 的成对 hash 按 §0.0h 同格式追加到该列表。（开工前两个 hash **必须成对**各自确认一遍。）

- **门槛与网格（全部在树实测）**：六门槛 **427/0/0 · 714/0 · 36/36 · 26/26 · 76/76 · 1/1**（kernel 712 → **714** = 两个判别性回归测试）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**16 个域网格与追加 32 基线逐项相同**。
- **★ 批次 1（主代理）：`UpdateEdge` 家族**全库收敛完毕**（追加 32 队列第 3 项）**。
  **9 个站点**改走唯一真身 `rcad_kernel::topods::update_curves_range`（`topods.rs:3929`）：`fillet/hbuilder_face/classify.rs::bb_update_edge_pcurve`（本域最后一处手搓体）· `fillet/hbuilder_face.rs` · **`fillet/hbuilder.rs:464`（★ 审计清单外的第 9 处）** · `bop/algo/pave_filler.rs` ×2 · `bop/algo/pave_filler_make_blocks.rs` ×2 · `bop/ds/mod.rs` ×2。
  **漏网原因**：追加 32 的清点是按 `ed.range` / `the_edge.range` 这类**字面变量名** grep 的，而该站点用的是 **`edp.range`** ⇒ 直接漏过。**改用名无关模式** `is_infinite_value\([A-Za-z_0-9]*\.range\[` 后，全库只剩 canonical 体内两行。**这与坑 28（"不要从截断清单数数"）同族 —— 同一处栽了第二次。**
  三个文件里 `is_infinite_value` 随之成为**死导入并删除**（`hbuilder_face.rs`/`hbuilder.rs` 保留同行的 `CONFUSION`）。6 文件 **+54/−111（净 −57 行）**，探针 = 0。
- **★ 批次 2（子代理，只读）：`gp.rs` 矩阵约定全量审计（813 行全读）**。
  **32 项判 ALIGNED**，含最易写反的配对：`transform_vec` 必须用 **`VectorialPart()`（含 scale）**、`transform_dir` 必须用 **`HVectorialPart()`（原始 matrix）** —— rcad **两处都取对**（互换即双重取反缺陷）。`invert` 里的 `M^T·loc` **合法**（OCCT 本身先就地转置），与 `set_rotation` 的误用被审计**明确区分**。
  **★ 查出 2 个同类 DEFECT**：
  1. **`Trsf::set_displacement`（`gp.rs:365-433`）两处转置乘** —— `MA1loc.Multiply(MA1)` 与 `MA1loc.Multiply(matrix)` 语义都是**矩阵左乘**（`gp_Trsf.cxx` L218-240），rcad 按行向量右乘写（= `MA1ᵀ L`、`Mᵀ MA1loc`）⇒ 平移项错。**FromA1 为零位置+单位姿态时两项错值都被乘 0**，而唯一调用点正是 `Ax3::new()`。
  2. **`Trsf::multiplied`（`gp.rs:620-644`）漏乘 `scale`** —— OCCT general 分支（`gp_Trsf.cxx` **L544-559**）为 `Tloc = matrix * T.loc; if (scale != 1.0) Tloc *= scale;`，rcad 少了 `scale` ⇒ `A(B(p))` 在 `A` 含非单位缩放时少缩放一次平移。唯一调用点（`sweep_section_generator.rs:478,482`）**scale 全为 1**。
  **根因侧观察**：本轮之前 `gp.rs` **只有 4 个测试**（3 × `set_rotation` + 1 × `set_translation`），`set_displacement`/`multiplied`/`transform_*`/`invert` **零测试** —— 与追加 32 同构：**不是没对齐，是没断言**。
  另记**结构性缺口（非缺陷）**：`PreMultiply`/`Power`/`Orthogonalize`/`SetTransformation(FromA1,ToA2)`/`SetValues`/`SetScale`/`SetMirror` 无 rcad 对应。
- **★ 批次 3（主代理）：两处修复 + 判别性测试（实证判别性）**。
  修复照 OCCT：`set_displacement` 两处改**矩阵左乘**；`multiplied` 补 `if (scale != 1.0) Tloc *= scale`。
  新增测试：`gp::tests::set_displacement_carries_from_frame_onto_to_frame`（**FromA1 既移位又倾斜 + ToA2 绕 Z 转 90°**；断言定义性不变量 `F(FromA1.Location()) == ToA2.Location()` 与三条轴向映射）· `gp::tests::multiplied_scales_the_left_translation`（**scale 3 × 平移**；断言 `A(B(p)) = 3p + (3,6,9)`）。
  **★ 判别性已实证**：把公式**临时改回旧形**再跑 ⇒ **两个新测试双双 FAILED、其余四个仍 ok**；还原后 **6/6 通过**。故二者是回归测试而非恒真断言。
- **★ 方法学（第六轮"零可见翻转"，本轮证据最硬）**：本轮修的是**真实内核缺陷**，但在八网格与 16 个域网格上**全部零可见**，因为**唯一调用点恰好都落在缺陷的退化配置上**。
  ⇒ **"触发用例缺口"卡改写为"补判别性单测"**：暴露条件可精确刻画（**轴位置 / 角度是否对称 / 缩放值 / FromA1 姿态**），一条 20 行单测即可永久钉死，比新造几何用例便宜得多。**优先覆盖 `gp.rs` 里仍零测试的 `apply`/`invert`/`transform_vec`/`transform_dir`/`set_transformation_ax3`。**
- **⏭ 下一轮队列（按序）**：
  1. **★ 补判别性单测**（最高优先，第六次重申；口径已由"新用例"改为"判别性单测"）：范式 = 上文四个变量；先扫 `gp.rs` 零测试的助手。
  2. **`shhealing/shape_extend/composite_surface.rs` 的本地 `Trsf2d` 同类审计**（新立，未审）：其 `multiplied`（L410）同样把 scale 与 translation 组合 ⇒ 潜在同族。**注意它是本地类型，不是 `rcad-kernel::math::gp`。**
  3. **`gp.rs` 结构性缺口补译**（按需）：`SetTransformation(FromA1,ToA2)`（`gp_Trsf.cxx` L172-194，与 `set_displacement` 是**两条不同路径**）· `PreMultiply`（L713）· `Power` · `Orthogonalize` · `SetMirror`/`SetScale`/`SetValues`。
  4. **`Surface3::dn` 的 `Trimmed` 臂**（便宜）+ `geom_adaptor_surface.rs::dn_at` 全阶。
  5. **`mySn` 的构造**（TKFeat，`feat_featrf` 前墙）。
  6. 追加 32 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · BSpline VIso 两处内嵌副本 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · `brep_tool_curve_on_surface` 缺 `CurveOnPlane` 回退 · 过期锚点勘误）。
- **复测脚本（本轮新增，可复用）**：`rcad/temp/run_gates.sh`（六门槛）· `rcad/temp/run_domain_grids.sh`（16 域网格）· `rcad/temp/validate_round5.sh`（门槛 + 重编 + 八网格 + 域网格一次跑完）。

### 0.0h 追加 32 收尾态（2026-09-13；**已被追加 33 取代 —— 见 §0.0i** —— 翻译优先轮第四批 + 一次回归的从根修复）

> **提交链**：**rcad 顶尖 = 本交接文件所在提交**（用 `cd rcad && git log --oneline -1` 即得）；根仓库指针 = 本文件所在的根提交。
> 逐轮（rcad × 根）：`783b0713 × 4a4fa29`（追加 32 代码）· `b651f36b × 102d130`（追加 31）· `c01540da × 34a225d`（追加 30）·
> `ebc7a30f × 09c11bf`（追加 29）· `df607f8c × 9aa9858`（追加 28）· `3fa1fb7d × 7f1fcb8`（追加 27）。
> 开工前用 `cd rcad && git log --oneline -1` 与 `cd .. && git ls-tree HEAD rcad` 各自确认一遍（**两个 hash 必须成对**）。

- **门槛与网格（全部在树实测）**：六门槛 **427/0/0 · 712/0 · 36/36 · 26/26 · 76/76 · 1/1**；**八网格 8/8**（375/378/379/373/12/102/83/110）；**15 个域网格与追加 31 基线逐项相同**。
- **★ 本轮最重要的产出：修掉一个真正的 kernel 缺陷 —— `gp_Trsf::SetRotation` 的 `loc` 用了转置矩阵。**
  OCCT `gp_Trsf::SetRotation`（`gp_Trsf.cxx` L90-99）里 `loc.Multiply(matrix)` 的语义由 `gp_XYZ::Multiply(const gp_Mat&)` 定义为 **`<me> = theMatrix * <me>`**（`gp_XYZ.hxx` **L308-309**，**矩阵左乘**）。rcad 写成了行向量右乘（= `M^T l`），且注释自称"row-vector times matrix" ⇒ **任何「轴不过原点」的旋转平移项都是错的**（轴过原点时 `loc = 0`，缺陷不可见）。
  **发现路径**：子代理接入解析 `Geom_SurfaceOfRevolution::EvalD0`（走 `Trsf`）后，八网格 `bcut_simple` 由 110/0 **回归**为 109/1（`g6` 面积 41187.4 → 61953.63，**拓扑全过只有面积错**）。`git stash -- libs/rcad-kernel/src/geom/` ⇒ g6 通过 ⇒ 定位到该批次；组内二分后**关键判读**是「差分导数与 `point_at` 恒自洽 ⇒ 旧 `point_at`+差分 = 41187.4 正确、新 `point_at`+差分 = 72764.7 错误 ⇒ **新解析体不等价**」；再用临时单测把 `revolution_eval_d0` 与 Rodrigues 式直接对拍，得 `diff = (38.94, 0, 0)`（**纯平移项偏移**）⇒ 追到 `set_rotation`。修正后临时测试 **`diff = 0`**，g6 通过且**批次 2 的解析体全部保留**。
  **为何长期潜伏**：旧 `RevolutionSurface::point_at` 是手搓 Rodrigues、**绕过了 `Trsf`**，掩盖了内核缺陷；既有的 `set_rotation_keeps_axis_location_invariant` 用 **π 旋转**，而 `M(π)` 对称 ⇒ `M^T L == M L`，**该测试原理上抓不到**。已补判别性测试 `set_rotation_about_off_origin_axis_matches_rodrigues`。
- **本轮另外三批**：
  1. **`UpdateEdge` 家族收敛成唯一真身** `rcad_kernel::topods::update_curves_range`（`topods.rs:3929`；缺陷前身 `pc_parameter_range` 删除）：两个 kernel 方法再宿主（`update_edge_pcurve_closed` **删掉 OCCT 本就没有的形参**）、**16 个调用点**随之改、另收敛 **7 处**散手副本（含审计清单外的 `brep_sweep/brep_sweep_builder.rs`）；`brep_algo/tool.rs::builder_range_edge` **修好**（OCCT 会把区间写到**每一个**表示）；`brep_fill_filling.rs` 修掉一个真缺陷（`BRep_TEdge::EmptyCopy` **会**拷 curve 表示，故有限 3D 覆盖必须生效）。
  2. **`Surface3::dn` 余留五面型补全**（Offset / LinearExtrusion / Ellipsoid / Helicoid / Revolution，全部解析体）：新增 `geom/{eval_c,revolution_utils,curve_dn,offset_surface_utils_b}.rs`。**Ellipsoid 的参数化差异（OCCT 纬度 vs rcad 余纬）已显式转换并双锚点记档**，未用等价替换糊过去。Ruled/Coons/Pipe/TriBezier **确认无 OCCT 对应**（枚举了全部 `Geom_Surface` 子类与 `GeomEval_*`），保留 GAP。
  3. **两处接线（主代理）**：`EvalAndUpdateTol` 真调 `GeomLibCheckCurveOnSurface`；`CurveOnSurface::ShallowCopy` 归位为 trait override。
- **★ 方法学（本批给出了与前三批不同的证据）**：本轮的**内核修复是真实行为修正**（所有轴不过原点的旋转），但在域网格上仍**零可见翻转**。
  ⇒ 八网格里**能感知内核几何修正的用例极少**，而这次覆盖是**回归逼出来的**、不是用例设计出来的。**"触发用例缺口"卡的优先级应再提一档**（第五次重申）。
- **⏭ 下一轮队列（按序）**：
  1. **★ 触发用例缺口 + 与 OCCT 参考对拍**（最高优先，第五次）。
  2. **★ `gp.rs` 同类「转置/左右乘」审计**（新立，高优先）：本轮证明这类手写矩阵约定极易写反且长期不可见 ⇒ 逐个核对 `Trsf`/`Mat`/`Ax2`/`Ax3` 里所有矩阵助手（`multiply`/`transformed`/`apply`/`transform_vec`/`set_displacement`/`set_rotation` …）的 OCCT 语义。
  3. **`UpdateEdge` 家族最后收敛**：`fillet/hbuilder_face/classify.rs:1344` + `fillet/hbuilder_face.rs:1022` 重定向到 `update_curves_range`；`bop/**` 三处内联体归并；`fillet` 的 `bb_update_edge_pcurve` 归并。
  4. **`Surface3::dn` 的 `Trimmed` 臂**（便宜）+ `geom_adaptor_surface.rs::dn_at` 全阶。
  5. **`mySn` 的构造**（TKFeat，`feat_featrf` 前墙）。
  6. 追加 31 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · BSpline VIso 两处内嵌副本 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · `brep_tool_curve_on_surface` 缺 `CurveOnPlane` 回退 · 过期锚点勘误）。

### 0.0g 追加 31 收尾态（2026-09-13；**已被追加 32 取代 —— 见 §0.0h** —— 翻译优先轮第三批）

- **门槛与网格（全部在树实测）**：六门槛 **427/0/0 · 701/0 · 36/36 · 26/26 · 76/76 · 1/1**（⚠ `rcad-kernel --lib` 691 → **701**）**；八网格 8/8**（375/378/379/373/12/102/83/110）—— **本批改了 `bop/**` 布尔核心，八网格仍零回归**；**15 个域网格与追加 30 基线逐项相同**。
- **本轮三批**：
  1. **曲面导数层补全（两个真缺陷）**：`bspline_surface_dn` 的**索引越界**（任何 BSpline `Surface3::dn(1,0)` panic）**不是放大缓冲而是换引擎** —— 改为 `Geom_BSplineSurface::EvalDN` 的 1:1（`BSplSLib::DN` L1519-1605）；`BSplSLib::RationalDerivative`（2D，L87-300）落地；新建 `geom/eval_b.rs`（388 行）= `BSplSLib::D0/D1/D2/DN`；Bezier `EvalD0–DN` 全部走解析体。**★ 并修掉 `BezierSurface::point_at` 的有理求值缺陷**（原按「逐 V 列有理 u 求值 + 单位权 V 组合」，不是张量有理求值；权重除 `w[1][1]=2` 全 1 时改前 `(0.75,1,0)`、改后 `(1,1,0)`）。**★ 从 OCCT 捞回的关键细节**：`D1/D2/D3` 省略末实参 ⇒ `All` 默认 **true**（有理表偏移 +6/+3、+9/+3/+18/+6/+12），只有 `DN` 传 `false`。+10 单测（解析精确值 + 独立商法则 oracle）。
  2. **`UpdateCurves` 两步规则 10+ 站点对齐**（`brep_fill_sweep_b.rs` 两处优先站点 · `brep_algo/tool.rs` · `bop/ds/mod.rs` 两处**删掉 OCCT 本就没有的形参** · `bop/algo/builder.rs` 调用点 · `pave_filler.rs` · `pave_filler_make_blocks.rs` ×2 · `brep_fill_filling.rs` ×2）。一处审计后判正确未改（`brep_fill_evolved_d.rs:179`）。
     **★ 新立最高价值翻译批**：**kernel 侧 `BRepBuilder::update_edge_pcurve`（`topods.rs:3342`）+ `update_edge_pcurve_closed`（`:3383`）仍未对齐，是 55 个调用点的漏斗**（`brep_sweep/**`、`feat/loc_ope_generator.rs`、`feat/loc_ope_split_drafts*`、`shhealing/**` 8 文件、`topalgo/brep_lib/make_face.rs` 8 处、`hlr/topo_brep/ds_filler.rs`）。正确体已散在 **6 处**（子代理拒绝加第 7 份）⇒ **下一轮应把它收敛成一份 canonical 并把两个 kernel 方法再宿主上去**。
  3. **两处接线（主代理）**：`offset/draft_modification.rs::brep_tools_eval_and_update_tol` 由硬编码 `ct_is_done=false / error_status=2` 改为真调 `GeomLibCheckCurveOnSurface`（按 `BRepTools.cxx` L1276-1333）；`Adaptor3d_CurveOnSurface::ShallowCopy` 由临时自由函数**归位**为 `base/proj_lib/adaptor.rs` 的 trait override。
- **★ 方法学（第四次同一结论，且本批给出了更强证据）**：本批改了 `bop/**` 布尔核心（pcurve 区间口径）、换了导数引擎、修了 `point_at` 缺陷、把 `EvalAndUpdateTol` 从假值变真值 —— **15 个域网格 + 八网格全部零可见翻转**。
  ⇒ 八网格的**拓扑断言对 pcurve 区间口径不敏感** ⇒ 最高优先队列项现在是：**为 pcurve 区间补定向断言**（不是加拓扑用例），并与 OCCT 参考值对拍。现有定向单测只验证自洽性。
- **⏭ 下一轮队列（按序）**：
  1. **★ pcurve 区间的定向断言 + 与 OCCT 参考对拍**（最高优先，第四次）。
  2. **★ kernel `UpdateEdge` 再宿主对齐 + 6 处正确体收敛成一份**（最高价值翻译批，55 调用点）。
  3. **`Surface3::dn` 余留面型**：Offset（**体已在** `offset_surface_eval_d1`，只差接进 `SurfaceEval for OffsetSurface` —— 最便宜）· LinearExtrusion · Ellipsoid · Helicoid · Revolution（需 `Geom_RevolutionUtils`）。
  4. **`brep_algo/tool.rs::builder_range_edge`** 未传播到 `ed.pcurves`。
  5. **`mySn` 的构造**（TKFeat，`feat_featrf` 前墙）。
  6. 追加 30 余项不变（`extrema_gen_ext_cs.rs` 过期 PSO 栈 · BSpline VIso 两处内嵌副本 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · `brep_tool_curve_on_surface` 缺 `CurveOnPlane` 回退 · 过期锚点勘误）。

### 0.0f 追加 30 收尾态（2026-09-13；**已被追加 32 取代 —— 见 §0.0h** —— 翻译优先轮第二批）

- **门槛与网格（全部在树实测）**：六门槛 **427/0/0 · 691/0 · 36/36 · 26/26 · 76/76 · 1/1**（⚠ 基线：rcad-algo 422 → **427**，rcad-kernel 688 → **691**）；**八网格 8/8**（375/378/379/373/12/102/83/110）；**15 个域网格与追加 29 基线逐项相同**（零回归）。
- **本轮两批（本域迄今最大单轮翻译量，子代理合计新增约 4,600 行）**：
  1. **`GeomLib_CheckCurveOnSurface` 真身**（OCCT 774 行 → 新建 `rcad-kernel/src/base/geom_lib/check_curve_on_surface.rs` 1561 行）：`TargetFunc`(Hessian) / `Local` / `Perform` / `FillSubIntervals` / `PSO_Perform` / `MinComputing` 全部逐句。**落位按 module-map 判为 kernel `base/`**（TKGeomBase），**不是** `geomalgo/`（TKGeomAlgo）。删掉两份旧实现（algo 的 GAP 载体 + kernel 的 **257 采样替身**），消费点重定向。
     **★ 关键正确性发现**：OCCT 的 `Geom_BSplineCurve::Knots()` 返回**去重**结点而 `KnotSequence()` 才是扁平向量；rcad 存扁平 ⇒ 必须做拆分再宿主，**否则每个重复结点都会多出一个边界**。
  2. **`Geom_OsculatingSurface`（839 行）+ `Geom_OffsetSurfaceUtils` + `Geom_ExtrusionUtils` 落地 ⇒ canonical `UIso/VIso` 的 `Offset` 臂打通**（`brep_fill/brep_fill_sweep.rs`，含两个 evaluator 类 L505-597 与 AdvApprox 尾部）。
     **★ AdvApprox 接线现在忠实**：`derive=1` 到达真解析 `EvaluateD1`（`ElSLib::DN`/`BSplSLib::DN` + 奇异点的真 `Geom_OsculatingSurface`），**不再喂差分导数** —— 上一轮拒绝接线的理由（近似 D1）已解除。
     **仍 GAP**：`BSplSLib::RationalDerivative`（有理 BSpline 底）· `Geom_Surface::EvalD1..DN` 对 **Bezier**/Revolution/Offset/Ruled/Coons/Pipe/Ellipsoid/Helicoid/TriBezier 底（Bezier 尤其值钱：`Geom_BezierSurface::D1` 无翻译，且 `is_q_punctual` 需要它）。
- **★ 方法学（连续第三轮同一结论，已升为最高优先队列项）**：本轮（+上轮）的**全部**修复 —— `Approx_SameParameter`、`CheckCurveOnSurface`、`Offset` iso 臂、基本面型 iso 臂 —— 在 **15 个域网格上零可见翻转**。
  ⇒ **域网格覆盖不到这些路径**。**"触发用例缺口"卡升为最高优先**：先用探针统计现有域网格里哪些例真走到这些路径，再补**与 OCCT 参考值对拍**的定向用例（现有定向单测只验证自洽性，不验证与 OCCT 一致）。
- **⏭ 下一轮队列（按序）**：
  1. **★ 触发用例缺口 + 与 OCCT 参考对拍**（最高优先）。
  2. **`geom/eval.rs::bspline_surface_dn` 索引越界缺陷**（L3390-3459：导数缓冲按单极点分配而写 `n+1` 个 ⇒ **任何 BSpline `Surface3::dn(1,0)` 越界**；且不除以权重和）。
  3. **`bop/int_tools/extrema_gen_ext_cs.rs` 的过期 PSO 栈**（RNG 种子与 OCCT 8.0.0 不同 ⇒ `Extrema_GenExtCS` 随机序列不同）。
  4. **`offset/draft_modification.rs:276-289` 的 `EvalAndUpdateTol` 接线**（前置已满足）。
  5. **`BSplSLib::RationalDerivative` + `Geom_BezierSurface::D1`**（`Offset` evaluator 的下一步）。
  6. **`Adaptor3d_CurveOnSurface::ShallowCopy` 移入 `base/proj_lib/adaptor.rs`**。
  7. **`mySn` 的构造**（TKFeat，`feat_featrf` 直接前墙）。
  8. 追加 29 余项不变（BSpline VIso 两处内嵌副本 · TKOffset 10 处 `UpdateCurves` 同族缺陷 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · `brep_tool_curve_on_surface` 缺 `CurveOnPlane` 回退 · 过期锚点勘误）。

### 0.0e 追加 29 收尾态（2026-09-13；**已被追加 32 取代 —— 见 §0.0h** —— 翻译优先轮）

- **门槛与网格（全部在树实测）**：六门槛 **422/0/0**（⚠ 基线 418 → **422**）**· 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；**八网格 8/8**（375/378/379/373/12/102/83/110）；**15 个域网格与追加 28 基线逐项相同**（零回归）。
- **本轮三批**：
  1. **`Approx_SameParameter` 全类译全**（992 行 cxx → 1465 行；原来是 `new` 即 GAP 的骨架）。`chfi3d_same_parameter` 的修正步换成 OCCT L1612-1621 的字面语句序，GAP 消除。+3 单测。
  2. **`UIso/VIso` 重复收敛**：删 **9** 处重复（含清单外 grep 出的 4 处，其中 `brep_fill/generator.rs::surface_uiiso` **忽略 `u`** 是真缺陷），消费点全部重定向到 canonical 真身。判非重复而保留的 2 处已立卡。
  3. **★★ 收敛暴露并修掉 canonical 真身的**基本面型臂**缺陷**：`brep_fill/brep_fill_sweep.rs::{surface_uiso,surface_viso}` 的 **Cylinder U/V-iso 互换**（U-iso 应是母线 **Line**、V-iso 应是 **Circle**）、**Sphere UIso 帧错（经线差 90°）且缺 `Trimmed(−π/2, π/2)` 包裹**、**Sphere VIso 的 sin/cos 互换**、**Cone VIso 漏 `Radius +`**。修法 = 委派到 kernel 已 1:1 的 `elslib_iso.rs`（OCCT 这些函数本就是 `ElSLib` 的薄包装）。
     **为何长期潜伏**：该文件单测**只覆盖 Bezier**。已补 `elementary_surface_isos_follow_the_occt_wrappers`（5 面型 × U/V-iso 共 10 条，9 点几何对拍 + 类型断言）。
     ⚠ **连带承认**：追加 27 把 fillet 的 iso 分支导流到这两个函数，因此**圆柱 iso-v 在追加 27 之后拿到错曲线**（旧代码用的 `cylinder_v_iso` 是对的）；八网格/域网格当时未捕获。本轮修 canonical 即一并修好，且受益面是**所有** canonical 消费者（brep_fill sweep / fillet `SplitSurf`+`ComputeArete` / offset `EnlargeGeometry`）。
- **★ 本轮最重要的方法学发现**：**连续两轮"潜伏缺口清除"在 15 个域网格上零可见翻转**（本轮改了基本面型 iso 的真实行为、删了 9 处重复、让 `chfi3d_same_parameter` 真跑，全都没翻）。
  ⇒ **结论：域网格的覆盖不足以验证这类修复**。**新立卡（队列第 1 项）**：为基本面型 iso / `Approx_SameParameter` / `SameParameter` 找或造**能区分**的用例（先用探针统计现有域网格里到底哪些例真走到这些路径）。
  ⇒ 也再次印证：**"通过数不变"不能当作"没进展"**，反之"零回归"也不能当作"改动被验证过"。
- **⏭ 下一轮队列（按序）**：
  1. **★ 触发用例缺口**（本轮新立，高优先）：为基本面型 iso 臂 / `Approx_SameParameter` / `chfi3d_same_parameter` 补能区分的用例（定向单测或探针定位）。
  2. **`mySn` 的构造**（TKFeat，`feat_featrf` 直接前墙，已定界到 `LocOpe_Revol` / `BRepFeat_MakeRevolutionForm`）。
  3. **`UIso/VIso` 的 `Offset` 臂**：前置 `Geom_OffsetSurfaceUtils::EvaluateD1` + `Geom_OsculatingSurface`（839 行）+ `GeomEval_RepSurfaceDesc`。
  4. **BSpline VIso 两处内嵌副本**（`brep_fill_nsections.rs:1027` / `geomfill/nsections.rs:280`；`de_boor_homo` 近似 vs 真身精确 `BSplSLib::Iso`）。
  5. 追加 28 余项不变（TKOffset 10 处 `UpdateCurves` 同族缺陷优先 · 池外读取链复核 · blend 剩余十例 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 · `BRepFill_Pipe` 收敛 · `BRepExtrema*`/`GeomIntIntSS` 重复 · `brep_tool_curve_on_surface` 缺 `CurveOnPlane` 回退 · 过期锚点勘误）。

### 0.0d 追加 28 收尾态（2026-09-13；**已被追加 32 取代 —— 见 §0.0h**）

- **门槛与网格（全部在树实测）**：六门槛 **418/0/0**（⚠ **基线已由 415 升到 418** = 本批新增 3 个 `brep_fill_sweep` 单测，含转正的那个）**· 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；**八网格 8/8**（375/378/379/373/12/102/83/110，重编 exe 后）；**15 个域网格与追加 27 基线逐项相同**（零回归，`blend_simple` a1 保持通过）。
- **本轮四批**：
  1. **★ TKFeat：`LocOpe_SplitShape` 全类译全**（4 个新文件 ~3473 行：`feat/loc_ope_split_shape{,_b}.rs` + **两个此前缺失的 OCCT 类** `topalgo/brep_tools_wire_explorer.rs`(885) 与 `topalgo/brep_lib_make_wire.rs`(512)）；`loc_ope_spliter.rs` 的 stub 整体删除。
     **效果 = 失败点换层**：`feat_featrf_a1` 由"测试断言 L866 面积 0"变为 **`feat/loc_ope_glued_shape.rs:193` 的 `Standard_ConstructionError`** ⇒ spliter 报 done、generator 跑起来、gluer 走到 OCCT 的 `OrientedFaces()`。
     **新墙已定性**：该 panic 是 **OCCT 自己的 throw 的 1:1**（`LocOpe_GluedShape.cxx` **L106-108**）⇒ **`mySn` 非闭合**（`myGShape` 的一张 4 边粘合面有一条边只有 1 个祖先面）⇒ **下一手 = `mySn` 构造（`LocOpe_Revol` / `BRepFeat_MakeRevolutionForm`）**。
     ⚠ **勿误读**：坑 16 的 `panic < assert` 是按**发现成本**排序，不是质量排序。本轮是"静默空结果 → OCCT 自己的显式错误"，方向正确（禁静默错几何），但**不是通过数提升**。
  2. **`Geom_BezierSurface::UIso/VIso` 两臂完成**（OCCT `Geom_BezierSurface.cxx` L1769-1810 / L1821-1863；含 `NCollection_Array2` 的 `RowLength()==NbColumns()` 反直觉命名这个坑）。**`Offset` 臂有意保留 GAP**，缺件已定位：`Geom_OffsetSurfaceUtils::EvaluateD1` + `Geom_OsculatingSurface`（rcad 无）+ `GeomEval_RepSurfaceDesc`。`AdvApprox` 虽在库，但其 `SimpleApprox::Perform` 用 `derive=1` 喂**近似 D1** ⇒ 明确不接。
  3. **★ kernel 缺陷（跨域，对拍 OCCT 实证后修复）**：`bspl_slib_iso`（`rcad-kernel/src/math/bspl_lib.rs`）的 `weight_at` 在 `is_u == false` 支读**转置**权重。OCCT `BSplSLib.cxx` **L1678-1681** 为 `P = IsU ? Poles(index,j) : Poles(j,index)` / `w = IsU ? (*Weights)(index,j) : (*Weights)(j,index)` ⇒ **极点与权重同一对索引**；两个调用点已排好实参，故访问器两向都应是朴素 `Weights(i,j)`。已修。**影响面 = 全部有理 BSpline/Bezier 的 VIso**；子代理因它 `#[ignore]` 的新测试已转正。
  4. **TKFillet：`ChFi3d_CheckSameParameter` 1:1**（OCCT `ChFi3d_Builder_0.cxx` L1565-1596 的 45 采样循环全译）+ **`ChFi3d_SameParameter` 真跑检查**（旧实现是 `*tolreached = tol3d` 的空体）。**修正步 `Approx_SameParameter` 保留为显式 GAP** —— 其真身仍是骨架（`geomalgo/approx_same_parameter.rs::new` = GAP panic，992 行 cxx 未译）。签名改回 OCCT 的 `-> bool`。
     **实测：新 panic 路径在 15 个域网格上未被触发** ⇒ 属**潜伏口径修正**，非行为翻转（旧空体给的 `tolreached` 口径不对但下游分支恰好未受影响）。
- **⏭ 下一轮队列（按序）**：
  1. **`mySn` 的构造**（TKFeat，`feat_featrf` 的直接前墙，已定界到 `LocOpe_Revol` / `BRepFeat_MakeRevolutionForm`）。
  2. **`Approx_SameParameter` 全类**（TKGeomBase/Approx，992 行）—— `chfi3d_same_parameter` 修正步的真缺口。
  3. **`UIso/VIso` 重复收敛**（E3-Q 用户定规）：库中至少 **5 份**，唯一真身 = `brep_fill/brep_fill_sweep.rs::{surface_uiso,surface_viso}`；收敛后 offset/fillet 各域**白得 Bezier 臂**。
  4. **`Geom_Surface::UIso/VIso` 的 `Offset` 臂**（前置 = `Geom_OffsetSurfaceUtils::EvaluateD1` + `Geom_OsculatingSurface`）。
  5. 追加 27 队列余项不变（TKOffset 10 处 `UpdateCurves` 同族缺陷优先 · TKOffset 池外读取链复核 · blend 剩余十例 · `split_edge.rs:1392` 另一半 · `hbuilder.rs` pcurve 键 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 · `builder_set_degenerated` 的 fork 风险 …）。

### 0.0c 追加 27 收尾态（2026-09-13；**读这一节即可开工**，权威脉络见 port-plan 追加 27）

- **门槛与网格（全部在树实测，含本轮三批）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；**八网格 8/8**（375/378/379/373/12/102/83/110，重编 exe 后）；**15 个域网格的真实断言数与追加 25/26 基线逐项相同**（`blend_simple` **1 过/10 败**，`a1` 保持通过；`draft_angle` 1/48；`feat_featrevol` 1/44；`fillet2d_fillet2d` 5/0；`mkface_after_offset` 2/0；`mkface_after_extsurf_and_offset` 16/0；其余 0/N——后四者的余数是 `*_geometry_loads` 恒过占位）。
- **本轮三批（全部"翻译/接线"类 ⇒ 通过数不变是预期，关键指标是零回归）**：
  1. **TKOffset 入池（原队列第 2 项 = 根）已落地**：`BRepTools_Quilt` + `BRepAlgo_FaceRestrictor` 产物**入池**（13 文件 / 9 调用点随之改签名；`BRepBuilder::remove_from_compound` 新增）。设计判定：quilt **改收前导 `brep: &mut BRep`**、不自持池（自持池会让产物按原 index 别名到消费者池里 = 坑 2，比现状更糟）。`shell_registry` 镜像整体删除（池槽位**就是** OCCT 共享 TShape）。
  2. **TKFillet `ChFi3d_ComputeArete` 的四个 stand-in 分支全部 1:1 译全**：iso u/v（原只支持圆柱 + 用**位置互换**代替 `ReversedParameter`+`Reverse`，**正是坑 15 禁的手写互换**）→ 接线 `surface_uiso/surface_viso` 真身；`IFlag==0` 非 iso（原丢 `pardeb/parfin` 与 `C3d`）→ `BRepAdaptor_Curve::D1` + `ChFi3d_BuildPCurve` + Bnd_Box2d 事后检查 + `GeomLib::BuildCurve3d`；`else` → 接线**早已在库**的 `chfi3d_project_pcurv` + UV1 对齐；并补上 OCCT L2000 的 `tolreached = tol3d;`。顺带把 `ElCLib::AdjustPeriodic` 收敛成 **`rcad_kernel::math::el::elclib_adjust_periodic` 唯一真身**（fillet 那份是 `while` 循环简化版、kernel 那份是 `pub(crate)`，双双删除）。
  3. **TKFeat 队列第 0 项被实测推翻并改判**（见下）。
- **★ TKFeat 第 0 项改判（最重要的结论，覆盖 §4.6 TKFeat 0 的旧表述）**：探针实测失败点是 **`LocOpe_Gluer::Perform` 的第二个 early-out `!the_split.is_done()`**（`loc_ope_gluer.rs:545`），**`LocOpe_Generator::Perform` 根本没被进入**——旧表述"查 `LocOpe_Generator::Perform` 为何不达"是错的。真缺口 = **`LocOpe_SplitShape` 是 stub 类**（`feat/loc_ope_spliter.rs:287-384`；OCCT `LocOpe_SplitShape.cxx` 1776 行，约 **1,340 行未译**，grep 全库 0 命中 ⇒ **真身不在库**）。**分类 = 功能缺失，且属"翻译/接线"** ⇒ 是下一轮首选。
- **⏭ 下一轮队列（按序，取代 §4.6 里被覆盖的项）**：
  1. **翻译 `LocOpe_SplitShape` 全类**（TKFeat；`feat_featrf_a1` 的直接前墙）。
  2. **补 `Geom_Surface::UIso/VIso` 的 Bezier / Offset 臂**（本轮因批次 2 提级；`bspl_slib_iso` 已是 `BSplSLib::Iso` 的 1:1，Bezier 只缺隐式 knot/mult 描述）。
  3. **`ChFi3d_SameParameter` 仍是无体 stand-in**（`chfi3d_builder_0.rs:358-366`）—— 新立卡，`chfi3d_compute_arete` 之后同族的下一个翻译点。
  4. TKOffset 追加 26 第 1 项的 10 处 `UpdateCurves` 两步规则同族缺陷（`brep_fill_sweep_b.rs:234/:259` 优先）。
  5. TKOffset 池外读取链**复核**（入池后 `curve_pool_free`/`edge_data_pool_free` 等守卫应已无活路径；勿重复立卡）。
  6. 其余既有队列不变（blend 剩余十例 · `split_edge.rs:1392` 另一半 · `hbuilder.rs` pcurve 键 · `elclib_adjust_periodic` 残留两份 · `builder.rs`/`pave_filler.rs` 拆分 …）。
- **⚠ 本轮引入、已量化、实测未触发的新风险**：批次 2 把**所有**面型导流到 `surface_uiso/surface_viso`，其对 **Bezier / Offset / 及 rcad 独有变体（Ellipsoid/Helicoid/Pipe/Ruled/Coons/TriBezier）** 仍是 GAP `panic!`；改动前这些面型是**静默留空**。15 个域网格无新增失败 ⇒ 本批保留，但补臂列为队列第 2 项（**保留 OCCT 失败路径的 GAP 优于静默错几何**）。

### 0.0 当前状态速览（追加 25 收尾，2026-09-13；**rcad 顶尖 = 本交接文件所在提交**（用 `cd rcad && git log -1 --oneline` 即得；写就时基线为 `0de47019`）/ 根仓库指针 = 本文件所在提交，**均已推送**）

- **门槛与网格**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（375·378·379·373·12·102·83·110，**重编 exe 后**实测）。域网格**真实断言通过数**：`draft_angle` **1/49（`b3`，★ 本轮域内首个真实通过）** · `feat_featrevol` **1/45（`a5`）** · `blend_simple` **1/11（`a1`，★ 追加 25 的修复）** · `fillet2d_fillet2d` 10/10 · `fillet2d_chamfer2d` 2/2 · `mkface_after_offset` 4/4 · `mkface_after_extsurf_and_offset` 32/32 · 其余 0（`feat_featlf` 0/15 · `feat_featprism` 0/6 · `feat_featrf` 0/5 · `blend_complex` 0/2 · `offset_shape_type_a` 0/1 · `offset_shape_type_i` 0/12 · `offset_faces_type_i` 0/8 · `thrusection_specific` 0/26）。
- **★ 工作模式（用户指令，追加 18 起生效）**：**先完成代码的等价实现（近乎 1:1 的翻译），代码基本译完再开始调试/修测试**。⇒ 队列里**先取"翻译/接线"项**，取"调试/定界"项前先确认没有未译的 body 挡在前面。
- **本轮已落地（12 批，见 §0.7）**：geomplate 的 ProjLib 三处接线 · 池外 `BRep_Tool::Curve` · `GeomLib::BuildCurve3d` 家族 · 池外读取续链 · blend 侧两处字面翻译 · `BRepTools_Modifier` 家族（翻译 + feat/offset 两域接线）· `Geom_Surface::UIso/VIso` 唯一真身 · **TKBool/TKBO 复用三批**（brep_fill 布尔 API 接线 / `BRepAlgo_Loop` 归位 / `BRepOffset_SimpleOffset` 1:1 + MakerVolume 接线）。
- **★ 域网格实测失败地图（下一轮的对照基线；口径 = 逐例 file:line）**：
  - `offset_shape_type_i`：a1/a2 → `brep_offset_inter2d.rs:1060`（`EdgeInter: E2 carries no pcurve`，OCCT 同处也 raise ⇒ **状态**）；**a3/a4/d2/d3 → 已离开库内，落到测试断言（310/422/534/646）** —— 批次 2 + 补记 1 的直接收益，**该网格池外 panic 已清零**；e1/e2/e3/e4/e6/e7 → 测试断言（758/870/948/1060/1171/1283）。
  - `offset_shape_type_a`：a4 → `brep_algo/image.rs:159`（本轮由 `brep_offset_make_offset_c.rs:52` 推到这里）。
  - `offset_faces_type_i`：a1/a2 → `brep_offset_inter2d.rs:1060`；a5/a6/g1/g2/g6/g7 → 测试断言。
  - `blend_simple`（11 例）：a1 → 测试断言 **L113**（面积，门 = 无界 pcurve）；a2/p8/p9 → `proj_lib_h_comp_projected_curve.rs:449`（`D0` 的 `Standard_DomainError` ⇒ **状态**类）；**a3/a4/q4/q7 → 四例汇合在 `chfi3d_builder_2b.rs:583`**（q4 与 q7 由批次 5 深入两层汇合到这里 ⇒ **该处是本域下一堵最值得优先对齐的墙**）；q1 → 测试断言 L963；q2 → `chfi3d_builder_c1.rs:1162`；x1 → `base/convert/mod.rs:2008`。
  - `offset_shape_type_a`：a4 → `brep_offset_make_offset_c.rs:52`（**GAP**：`BRepAdaptor_Curve(E)` 的 curve-on-surface 3D 回退未译）。⚠ **待归因**：批次 3 曾把它推到 `brep_algo/image.rs:159`，批次 4/5 之后**退回**到这里（见追加 19 补记 2；`image.rs:159` 本身**不是缺陷**——OCCT 自己的 `BRepAlgo_Image::Root` L168 就抛那句复制粘贴的 `FirstImageFrom`）。
  - `feat_featrf`：a1 → 测试断言 **L866**（面积 0）；a4/a5/a7/a9 → 测试断言（**未劣化**）。
  - `feat_featlf`（15 例）：全部为测试断言（本轮实测**库内 panic 已清零**，与追加 18 记的 3 例库内 panic 相比又前进一层）。
  - `feat_featrevol`：44 例全为测试断言 + a5 通过（1/45）。
  - `offset_shape_type_a`：a4 → **不稳定**：`brep_offset_make_offset_c.rs:52`（**GAP**：`BRepAdaptor_Curve(E)` 的 curve-on-surface 3D 回退未译）或 `brep_algo/image.rs:159`，**同 exe 连跑 5 次为 3:2 分裂**（见坑 33）。注：`image.rs:159` **不是缺陷**——OCCT 自己的 `BRepAlgo_Image::Root` L168 就抛那句复制粘贴的 `FirstImageFrom`。
  - `draft_angle`（49 例）：库内 panic **28 = 17 × `draft_modification_1_b.rs` + 6 × `_1_c.rs` + 2 × `make_revol.rs` + 2 × `approx_int.rs` + 1 × `brep_offset_api_draft_angle.rs`**（其余 21 例测试断言）。**★ 在 `33ad8514`（追加 18 收尾树）上实测同样是 28 且构成逐项相同** ⇒ **本 session 五个批次对 `draft_angle` 零影响**；交接档旧记的"25"是**过期货**（见坑 7、补记 3）。
  - ⚠ **失败地图本身不稳定（★ 本轮最重要的方法学发现，见坑 33）**：同一棵树、同一份 exe，**连跑 5 次可以给出不同的"第一个失败点"**（实测 `offset_shape_type_a` a4 为 `image.rs:159` ×3 与 `brep_offset_make_offset_c.rs:52` ×2）。
    根源是 rcad 大量遍历 `std::collections::HashMap`/`HashSet`（**每进程随机哈希种子**）；OCCT 的 `NCollection_Sequence`/`IndexedDataMap` 是插入序确定的。
    ⇒ **报地图要连跑 ≥3 次、报集合与分裂比例**；**归因必须换树复测**；并**新立卡**：核对 offset/几何管线的容器选用与遍历序（这是**潜在的行为差异**，不只是报数问题）。
- **三域下一步（详见 §4.6）**：**TKOffset** = ① 池外读取**续链**（`check_same_range` / `gcurve_range` / `brep_tool_curve_on_surface_index` / `brep_tool_range_on_surface` / `brep_tool_degenerated` / `brep_tool_tolerance` 收敛到受守卫的 edge-data 读取；**写回侧池外无池可变，需另法，勿硬凑**）→ ② `BRepTools_Quilt`/`FaceRestrictor` **产物入池**（根；入池后整条池读链一次解开，且拓扑计数随之正确；**牵涉拓扑计数 ⇒ 全套复测**）；**TKFillet** = 接力项 A 的"无界 pcurve 附着点"（a1 面积 −2e100 的门，本轮未触及）+ blend a2/p8/p9 的新墙（状态类）；**TKFeat** = `LocOpe_Generator::Perform` 的 `IsDone`（a1 粘合路径下一墙）。

### 0.7 追加 19–25 本轮落地（**翻译/接线 + 首例测试修复**，2026-09-13；二十批，rcad 提交链见 §3）

> **★ 追加 25 收尾态 = 当前状态**（**首例测试断言修复** + kernel 容差补齐）：
> 批次 O **`blend_simple_a1` 通过**（`f8e9728d`）—— 无界 pcurve 结案：`fillet/hbuilder.rs::build_faces` 的 pcurve 写入只做了 `BRep_Builder::UpdateCurves` 的**前半**（二维曲线区间播种），**漏了"有限 3D 区间覆盖"**（`BRep_Builder.cxx` **L157-165**）；
> 批次 N **kernel 侧容差真身**（`0de47019`：新 `kernel/src/topo/brep_tool.rs`，删 `brep_adaptor.rs` 的无下限副本，`BRepTool::{tolerance,vertex_tolerance}` 改委托）。
> **★ 域网格首次出现真实断言翻转**：`blend_simple` **11/11 → 12 passed / 10 failed**，**失败地图 diff 恰好一行**（a1 消失）。
> **本轮最值钱的通用教训**：OCCT 的 `BRep_Builder::UpdateCurves` 是**两步区间规则**（二维播种 + **有限 3D 覆盖**），rcad 只落了前半；
> 且**同族第二处**（`hbuilder_face.rs:1019-1029`）**早已修过**，`hbuilder.rs` 是被漏掉的那一处
> ⇒ **凡"写 pcurve 区间"的位置都要按这条规则逐点核对**（已列为队列第 1 项）。

> **（上一状态，追加 23）给真 TKBO 体补 OCCT 公共 facade 轮**：
> 批次 H `bop/**` **补 OCCT 公共 facade**（`82383294`：`BOPAlgo_Builder::{add_argument, perform_with_filler, perform_internal1, build_bop, build_bop_states, clear}` + **history 读面** + `BRepAlgoAPI_BuilderAlgo` 带 filler 形态 + `BRepTools_History::from_algorithm/merge_algorithm` + **MakerVolume 的 images 不再被丢弃**）·
> 批次 I **`BRep_Tool::Tolerance` 收敛出带下限的唯一真身**（`8ad5f981`，删 5 份重复；★ **另发现尚存 24 份重复**）·
> 批次 J **brep_fill 的 #7 偏差消灭**（`90375e64`，新 facade 的第一个消费者）。
> **本轮三条硬信息**：① facade 补齐时**顺手修掉两处翻译 bug**（`TopExp_Explorer(S, SOLID)` 对 solid 参数必须产出 `S`；MakerVolume 组合 Builder 未继承 `myFillHistory`）；② **#8/#9 本轮未消灭，障碍已精确记档**（#8 = `PaveFiller::set_glue` 表达不了 `GlueShift`；#9 = DS 深拷贝使 history merge 会**静默为空**，须先经 `ds.argument_remap`）；③ `bop/algo/builder.rs` 现 **8012 行**（基线即超 2000 行规范）⇒ 已立卡。
> 验收：六门槛全绿 · 八网格 **8/8** · **十个域网格逐例失败地图与追加 22 基线逐字节相同**（零回退）。

| 批次 | 提交 | 内容 | 实测 |
|------|------|------|------|
| **M（追加 24）** | `7d725b00` | **`bop/**`+`brep_fill`**：架构差异 **#8/#9 消灭** —— 枚举版 `SetGlue` 补在 `BOPalgo_PaveFiller`/`BOPalgo_Builder`（**不在** `BOPalgo_Options`），旧的 bool `set_glue`（无 OCCT 重载、零调用者）删除；**Builder 趟改为只跑一次**（BuildBOP 复用 images，此前每个都重跑并重分割 DS = 行为纠正）；#9 经 `ds.argument_remap` 后**真的执行且不再静默为空** | 六门槛 + 八网格基线；⚠ **`BRepFill_Draft::Fuse` 从测试不可达** ⇒ 两处消灭**无网格暴露** |
| **L（追加 24）** | `a2116505` | **容差切片 A**：feat/fillet/hlr **八份** re-host 收敛并入唯一样本；**并修掉一处潜伏缺陷** —— `feat/loc_ope_pipe.rs` 的 reader **只有 Vertex 臂**却接 **Face**（`LocOpe_Pipe.cxx` L237）⇒ **静默返回 0.0** | 四域网格 3×3 一致 |
| **K（追加 24）** | `cbba7f55` | **容差切片 B**：shhealing/topalgo **13 份**收敛（27 文件，+68/−171），含 `brep_check_result.rs` 三个 per-kind helper 合并 ⇒ **repo 级 `brep_tool_tolerance` 定义由 ~30 份降到 3 份** | 六门槛 + 八网格基线 |
| **O（追加 25）** | `f8e9728d` | **fillet**：**`blend_simple_a1` 通过** —— `hbuilder.rs::build_faces` 的 pcurve 写入补上 `BRep_Builder::UpdateCurves` 的**第二步**（有限 3D 区间覆盖，`BRep_Builder.cxx` L157-165）；同族第二处（`hbuilder_face.rs:1019-1029`）此前已修，`hbuilder.rs` 是被漏掉的那处 | **`blend_simple` 11/11 → 12 passed / 10 failed**；面积 **59527.876** vs OCCT 59527.9；`step-topo-diff` fully matches；**失败地图 diff 恰好一行** |
| **N（追加 25）** | `0de47019` | **kernel**：新 `kernel/src/topo/brep_tool.rs::brep_tool_tolerance(&TShape)`（三臂带 `Confusion` 下限），删 `brep_adaptor.rs` 的无下限副本，`BRepTool::{tolerance, vertex_tolerance}` 改委托；校正两处错锚（`BRepAdaptor_*::Tolerance` 实在 L92-95/L146-149） | 六门槛 + 八网格基线；⚠ **设计点**：池读者必须读**池槽**（首版改读 `Shape::data` 当场打破 kernel 单测 688/1，改回即 689/0） |
| **J（追加 23）** | `90375e64` | **brep_fill**：`BRepAlgoAPI_Section(Sol1, Sol2, aPF)` 改用 `SectionOp::from_shapes_with_filler` + `build_with_filler` ⇒ **架构差异 #7 消灭** | 十域地图**逐字节相同**（零回退） |
| **I（追加 23）** | `8ad5f981` | **容差**：`brep_algo/tool.rs::brep_tool_tolerance` 成为**三臂带 `Precision::Confusion` 下限**的唯一样本，删 5 份重复并改 import 行 | 六门槛 + 八网格基线；★ **尚存 24 份重复**（队列第 2 项） |
| **H（追加 23）** | `82383294` | **`bop/**` 补 OCCT 公共 facade**：Builder 的 `add_argument`/`perform_with_filler`/`perform_internal1`/`build_bop`/`build_bop_states`/`clear` + `BOPAlgo_BuilderShape` 的 history 读面 + `BRepAlgoAPI_BuilderAlgo` 带 filler 形态 + `BRepTools_History` 模板构造（trait 承载）+ MakerVolume 持有并暴露 images | 八网格 8/8 零回归；新路径用**一次性集成测试**端到端验过（用后已删） |
| **G（追加 22）** | `b8711d80` | **offset**：`BRepOffset_SimpleOffset` 1:1（`BRepOffset_SimpleOffset.cxx` L1-427，六个 override 的 GAP panic 清零）+ `BOPAlgo_MakerVolume` 桩 → `bop/algo/maker_volume.rs` | 域网格不变；⚠ **该 mapper 无测试覆盖**（唯一调用者无调用者）⇒ 形式对齐、运行时未验证 |
| **F（追加 22）** | `507ca574` | **feat**：`feat/loc_ope_generator_b.rs` 的本地 `BRepAlgoLoop`（头注"NOT YET PORTED"）→ 真身 `brep_algo/loop.rs` | feat 四网格失败层**完全不动**（无案例走到该路径）⇒ 属潜伏 `unimplemented!()` 清除 |
| **E（追加 22）** | `f3556fea` | **brep_fill**：删三个布尔桩 + 两个偏译枚举 + **一个遮蔽真身的 `BRepToolsHistory` 重复**；七处调用点接真 TKBO 体（arch diff #7/#8/#9 已记档） | 八网格 8/8 零回归 |
| **D（追加 21）** | `a374a286` | **offset 域接线**：删两个本地 `BRepToolsModifier` 载体（draft_angle / make_simple_offset），改用真身；新增 `impl BRepToolsModification for DraftModification`（六个 override **委托**既有 inherent 方法 ⇒ 一个引擎） | **`draft_angle` 49/49 → 50/48（`b3` 通过）**；库内 panic **28 → 27**；其余域网格不变 |
| **C（追加 21）** | `ffed3791` | **feat 域接线**：删 `loc_ope_prism.rs` 的本地载体 + 本地 `BRepToolsTrsfModification`，四个 `int_perf` 改用真身；`brep_fill_evolved_c.rs` 改用 OCCT **两参构造** `BRepTools_Modifier(S,M)` | `feat_featlf` **a3 深入一层**；featprism/featrevol/featrf 不变 |
| **A（追加 20）** | `a1490588` | **`BRepTools_Modifier` + `BRepTools_Modification` 1:1 翻译**（新 `topalgo/brep_tools_modification.rs` + `brep_tools_modifier.rs`，共 ~2,170 行；**不接线消费者**） | 无在役消费者 ⇒ 无可观测变化属预期 |
| **B（追加 20）** | `745a63aa` | **`Geom_Surface::UIso/VIso` 八臂并集唯一真身**（新 `geomalgo/geom_surface_iso.rs`）+ 关掉 `approx_curve_on_surface.rs` 的 `AdvApprox_PrefAndRec` / 1D 子空间 GAP | 域网格**逐例失败地图与基线逐字节相同**（无回退）；顺带修一处保真缺陷（Sphere UIso 缺 `Geom_TrimmedCurve` 包装） |
| 1（追加 19） | `d7beea3c` | **`ProjLib_HCompProjectedCurve` 在 `GeomPlate_BuildPlateSurface` 三处接线**（metrics 比较 cxx L1746-1802 / ProjectCurve L254-303 / ProjectedCurve L307-349）——真身与 `Adaptor3d_CurveOnSurface` 真身早已在库 | `blend_simple` **a2/p8/p9** 由 GAP panic 推进到 `proj_lib_h_comp_projected_curve.rs:449` |
| 2 | `830eb2c2` | **池外 `BRep_Tool::Curve`**：`topods.rs` 新增 `curve_pool_free`（OCCT **BRep_Tool.cxx L172-196**）+ `shape_is_in_pool`；`build_curves3d.rs::brep_tool_curve` 与 `topexp.rs::{brep_tool_curve_loc,brep_tool_range}` 按守卫分流 | `offset_shape_type_i` **a3/a4/d2/d3** 由 `brep_tool_curve` 内推进到 `BRepLib::check_same_range`（backtrace 取证） |
| 3 | `e841e805` | **`GeomLib::BuildCurve3d` 家族 1:1**（新增 `ApproxAFunction::with_cut_tool` + 子空间存储、`kernel_curve2d` 桥、`GeomLib_CurveOnSurfaceEvaluator`、`build_curve3d`、`isIsoLine`/`buildC3dOnIsoLine`、`GeomLib_MakeCurvefromApprox`）+ **四处 GAP 载体收敛** | `offset_shape_type_a` **a4** 推进到 `brep_algo/image.rs:159`（⚠ 该例失败点不稳定，见坑 33）；`draft_angle` **零影响**（旧记"25 → 20"已撤销，正确数字 = **28**，见补记 3） |
| 4 | `fe297616` | **池外读取续链**：kernel `edge_data_pool_free` + `build_curves3d.rs::edge_data` 守卫访问器，该链上**十处**池索引读（`check_same_range`/`same_range`/`gcurve_range` + 四个 `brep_tool_*` re-host）全部收敛 | `offset_shape_type_i` **a3/a4/d2/d3 离开库内**（→ 测试断言 310/422/534/646）⇒ **该网格池外 panic 清零** |
| 5 | `d7587fdb` | **blend 侧两处字面翻译缺陷**：`bsplclib_resolution` 的有符号 `ii - Deg1` clamp（OCCT `BSplCLib.cxx` L4481-4490，原 usize 下溢）+ `add_singular_point` 的 **1-based `jalons.Value(jj)`** 读取（OCCT `BRepBlend_Walking.cxx` L160，原越界） | `blend_simple` **q4/q7 双双深入两层**并与 a3/a4 **汇合到 `chfi3d_builder_2b.rs:583`** |

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
| blend_simple | **1/11** | **追加 25 起 `a1` 通过**（无界 pcurve 修复，面积 59527.876 vs OCCT 59527.9）⇒ 该域**首次出现真实断言翻转**；其余十例见 §1 的逐例地图 |
| blend_complex | **0/2** | 全败 |
| feat_featlf | **0/15** | 全败；`is_done` 门**11/15 已越过**（追加 14），**kernel panic 清零**，剩余全部为下游墙（§4.1） |
| feat_featprism / featrf | **0/6** / **0/5** | 全败 |
| feat_featrevol | **1/45**（a5 **已通过**，2026-09-12 尾实测；其余 44 败） | 追加 15 尾修正：旧记档"0/45"**已过期**——`feat_featrevol_a5` 的真实断言**现在通过**（作者已用 path-scoped `git stash` + 重编证明**不是**批次 B 带来的；最可能来自本 session 早期的**扫掠子形状/面构造**修复，未逐条归因）。**未逐格复核，仅此一条已确认** |
| offset_shape_type_a / _i / _i_c | **0/1** / **0/12** / **0/19** | 全败（**逐例失败地图见 §0.0**；追加 19 后 `_i` 的 a3/a4/d2/d3 池外帧前进一层，`_a` 的 a4 推进到 `brep_algo/image.rs:159`） |
| offset_faces_type_i | **0/8** | 全败 |
| draft_angle | **1/49**（★ 追加 21 起不再是 0/49） | `draft_angle_b3`（零角度 draft）**已通过**（走完 `is_done`+拓扑+面积断言）；当前库内 panic **27** = 17 × `draft_modification_1_b.rs` + 6 × `_1_c.rs` + 2 × `make_revol.rs` + 2 × `approx_int.rs`（其余测试断言）。**该数字的归属见追加 19 补记 3：`28` 在追加 18 收尾树上实测相同 ⇒ 早期"25/20"均为过期货** |
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

## 3. 提交链与落地内容（**十三轮**：追加 13…23 / 24 / **追加 25 = 本轮**）

**追加 25（本轮，2026-09-13；rcad `main` 顶尖 `0de47019`，根 `main` 顶尖 = 本文件所在指针 sync，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`0de47019`（**核心 N**：kernel 侧容差真身 `kernel/src/topo/brep_tool.rs`，删 `brep_adaptor.rs` 的无下限副本，`BRepTool::{tolerance,vertex_tolerance}` 改委托）
← `f8e9728d`（**核心 O**：**`blend_simple_a1` 通过** —— `hbuilder.rs::build_faces` 补上 `BRep_Builder::UpdateCurves` 的**有限 3D 区间覆盖**步；面积 **59527.876** vs OCCT 59527.9）
← `fb8449c8`（= 追加 24 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `0de47019`）← `6054bfb`。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
**十个域网格与追加 24 基线对拍，唯一差异 = `blend_simple_a1` 从失败地图消失**（`blend_simple` 11/11 → **12/10**），其余逐字节相同；探针 = 0。
**★ 本轮是"先译后调"路线在 TKFillet 域的第一次完整兑现**：形式对齐（补上 OCCT 的第二步区间规则）+ 同族第二处的既有先例 ⇒ 测试断言自然达标，且**无需任何数值调参**。

**追加 24（上一轮，2026-09-13；rcad `main` 顶尖 `fb8449c8`，根 `main` 顶尖 `6054bfb`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`7d725b00`（**核心 M**：架构差异 **#8/#9 消灭** —— 枚举版 `SetGlue` 补在 `BOPAlgo_PaveFiller`(hxx L152-153/cxx L107-110)/`BOPAlgo_Builder`(hxx L122-126)，旧 bool `set_glue`（无 OCCT 重载、零调用者）删除；**Builder 趟只跑一次**（BuildBOP 复用 images，此前每个都重跑并重分割 DS）；#9 经 `ds.argument_remap` 后真的执行）
← `cbba7f55`（**核心 K**：shhealing/topalgo **13 份**容差 re-host 收敛，27 文件 +68/−171；repo 级定义 **~30 → 3 份**）
← `a2116505`（**核心 L**：feat/fillet/hlr **八份**收敛 + **修掉 `loc_ope_pipe.rs` 的 face/vertex 混用**（静默返回 0.0））
← `8a60b90b`（= 追加 23 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `7d725b00`）← `21197ed` ← `410bfca`。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
**十个域网格逐例失败地图与追加 23 基线逐字节相同**（唯一差异是已知不稳定的 `offset_shape_type_a` a4）；探针 = 0；42 文件、**+247/−411（净 −164 行）**。
**本轮新立卡 3 条**（§4.6 第 1-3 项）：`rcad-kernel` 两份无下限容差读者的**内核侧决策** · **过期锚点批量勘误**（5 处，纯注释）· `builder.rs` 8026 行与 `pave_filler.rs` 6256 行的**拆分**。

**追加 23（上一轮，2026-09-13；rcad `main` 顶尖 `8a60b90b`，根 `main` 顶尖 `21197ed`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`90375e64`（**核心 J**：brep_fill —— 新 facade 的第一个消费者，`BRepAlgoAPI_Section(Sol1,Sol2,aPF)` 改用 `SectionOp::from_shapes_with_filler` + `build_with_filler` ⇒ **架构差异 #7 消灭**）
← `8ad5f981`（**核心 I**：`BRep_Tool::Tolerance` 收敛为**三臂带下限的唯一样本**，删 5 份重复；★ 另发现尚存 **24 份**）
← `82383294`（**核心 H**：`bop/**` 补 OCCT 公共 facade —— Builder 的 `add_argument`/`perform_with_filler`/`perform_internal1`/`build_bop`/`build_bop_states`/`clear` + history 读面 + `BRepAlgoAPI_BuilderAlgo` 带 filler 形态 + `BRepTools_History` 模板构造 + MakerVolume images；顺手修两处翻译 bug）
← `54d6eb00`（= 追加 22 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `90375e64`）← `410bfca` ← `6dd53db`。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
**十个域网格逐例失败地图与追加 22 基线逐字节相同**（`draft_angle` 仍 50/48，`b3` 未回归）⇒ 三批**零回退**；探针 = 0。
**本轮新立卡 4 条**（§4.6 第 1-4 项）：#8/#9 收尾（前置 = `PaveFiller::SetGlue(BOPAlgo_GlueEnum)`）· **24 份容差重复收敛** · `builder.rs` 8012 行拆分 · `brep_check_result.rs` 过期锚点。

**追加 22（上一轮，2026-09-13；rcad `main` 顶尖 `54d6eb00`，根 `main` 顶尖 `410bfca`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`b8711d80`（**核心 G**：offset —— `BRepOffset_SimpleOffset` 1:1（cxx L1-427，六个 override 的 GAP panic 清零）+ `BOPAlgo_MakerVolume` 桩 → 真 TKBO 体）
← `507ca574`（**核心 F**：feat —— 本地 `BRepAlgoLoop`（头注"NOT YET PORTED"**已过期**）→ 真身 `brep_algo/loop.rs`）
← `f3556fea`（**核心 E**：brep_fill —— 三个布尔桩 + 两个偏译枚举 + **一个遮蔽真身的 `BRepToolsHistory` 重复**全部删除，七处调用点接真 TKBO 体；arch diff #7/#8/#9 记档）
← `517eab4a`（= 追加 21 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `b8711d80`）← `6dd53db` ← `24ec0c3`。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
**域网格与追加 21 基线逐例对拍：九个网格逐字节相同**（含 `draft_angle` 仍 50/48、`b3` 未回归），唯一差异是已知不稳定的 `offset_shape_type_a` a4；探针 = 0。
**本轮性质**：**潜伏缺口清除轮**（无通过数变化——这些载体在测量网格上不可达），产出**五条证据化的新立卡**（§4.6 第 1-5 项）。

**追加 21（上一轮，2026-09-13；rcad `main` 顶尖 `517eab4a`，根 `main` 顶尖 `6dd53db`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`a374a286`（**核心 D**：offset 域接线 —— 删两个本地 `BRepToolsModifier` 载体 + `impl BRepToolsModification for DraftModification`；**`draft_angle_b3` 通过**，库内 panic 28 → 27）
← `ffed3791`（**核心 C**：feat 域接线 —— 删 `loc_ope_prism.rs` 本地载体 + 本地 `BRepToolsTrsfModification`，四个 `int_perf` 与 `brep_fill_evolved_c.rs` 改用真身；`feat_featlf` a3 深入一层）
← `050c1323`（= 追加 20 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `a374a286`）← `24ec0c3` ← `9190c73`（追加 19/20）。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
**域网格与追加 20 基线逐例对拍，变化只有三处且全为正向**（`draft_angle_b3` 通过、`feat_featlf` a3 深入、`offset_shape_type_a` a4 在已知不稳定点间翻转）；探针 = 0。

**追加 20（上一轮，2026-09-13；rcad `main` 顶尖 `050c1323`，根 `main` 顶尖 `24ec0c3`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`745a63aa`（**核心 B**：`Geom_Surface::UIso/VIso` 八臂并集唯一真身 `geomalgo/geom_surface_iso.rs` + 关掉 `approx_curve_on_surface.rs` 的 `AdvApprox_PrefAndRec` / 1D 子空间 GAP）
← `a1490588`（**核心 A**：`BRepTools_Modifier` + `BRepTools_Modification` 1:1 翻译，~2,170 行，**消费者未接线**）
← `11027474`（= 追加 19 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `745a63aa`）← `9190c73` ← `e61aafd`（追加 19）。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
**域网格逐例失败地图与追加 19 基线逐字节相同**（唯一差异是那 1 例已知不稳定）⇒ 两批**无失败层回退**；探针 = 0。
**★ 本轮最重要产出**：把"哈希序"**量化**（172 失败例 ×3 次连跑中仅 1 例失败点不稳定；48 个全绿在域用例 ×5 次零抖动）⇒ **降级为潜在隐患，不做重构**（见坑 33 与本轮正文）。

**追加 19（上一轮，2026-09-13；rcad `main` 顶尖 `11027474`，根 `main` 顶尖 `9190c73`，均已推送）**
— rcad `main`（自上而下 = 新到旧）：
`d7587fdb`（**核心 5**：blend 侧两处字面翻译缺陷 —— `bsplclib_resolution` 的有符号 `ii - Deg1` clamp（OCCT `BSplCLib.cxx` L4481-4490）+ `add_singular_point` 的 1-based `jalons.Value(jj)` 读取（OCCT `BRepBlend_Walking.cxx` L160））
← `ab0f7c9f`（docs：追加 19 补记 1）
← `fe297616`（**核心 4**：池外读取续链 —— kernel `edge_data_pool_free` + `build_curves3d.rs::edge_data` 守卫，十处池索引读收敛；`offset_shape_type_i` 池外 panic 清零）
← `e841e805`（**核心 3**：`GeomLib::BuildCurve3d` 家族 1:1 + 四处 GAP 载体收敛 —— `AdvApprox_PrefAndRec`、`ApproxAFunction::with_cut_tool` + 1D/2D 子空间存储、`kernel_curve2d` 桥、`GeomLib_CurveOnSurfaceEvaluator`、`GeomLib::isIsoLine`/`buildC3dOnIsoLine`、`GeomLib_MakeCurvefromApprox`）
← `830eb2c2`（**核心 2**：池外 `BRep_Tool::Curve` —— `topods.rs::curve_pool_free` + `shape_is_in_pool`，`brep_tool_curve`/`brep_tool_curve_loc`/`brep_tool_range` 按守卫分流）
← `d7beea3c`（**核心 1**：`ProjLib_HCompProjectedCurve` 在 `GeomPlate_BuildPlateSurface` 三处接线）
← `33ad8514`（= 追加 18 链尾）。
根 `main`：本轮 pointer sync（`rcad` 指针 → `e841e805`）← `7528964` ← `030caf6`（追加 18）。
**本轮验证（全部在树实测）**：六门槛 **415/0/0 · 689/0 · 36/36 · 26/26 · 76/76 · 1/1**；八网格 **8/8**（**重编 exe 后**）；
域网格逐格通过数不变（失败层深度见 §0.7 / 追加 19）；探针 = 0。
**本轮新坑 28/29/30/31/32/33** 见 §6 尾（同一 helper 行号掩盖层推进 / 旧笔记的"架构难点"会过期 / `GetType()` 常量返回决定分支 /
**截断清单数数** / 换树复测归因 / **★ 失败地图默认不稳定（HashMap 哈希序）**）。
**★ 本轮的两次自纠**：追加 19 正文的"`draft_angle` 25 → 20"与补记 2 的"两处待归因"**都已撤销**——
前者是截断清单数数（坑 31），后者经**换树复测 + 同树 5 次复跑**证明是**测量假象**（坑 32/33）。**以追加 19 补记 3 为准。**

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
5. **调试类（按用户口径排在翻译/接线项之后）**：blend a2/p8/p9 的 `D0` 域错误；`offset_shape_type_i` a1/a2 的 `EdgeInter: E2 carries no pcurve`；`draft_angle` 的 28 例库内 panic（**在追加 18 收尾树上同样存在**，同 OCCT 亦 raise 的状态类）。
6. **★ 已完成（追加 21）：接线 `BRepTools_Modifier` 消费者** —— `feat/**`（`ffed3791`）与 `offset/**`（`a374a286`）两域都已收敛，`draft_angle_b3` 因此通过。
   **残留两处（已记档，非本批失误）**：① `BRepOffset_MakeSimpleOffset::Perform` 的 `aBB.Degenerated(anEdge,true)`（OCCT L197-198）走既有 `builder_set_degenerated` 载体，而它用 **`Arc::make_mut`** ⇒ 若 map 仍别名该边 TShape，标志会落在 **fork** 上（坑 17 家族；就地形态私有于 `topalgo/brep_tools_modifier.rs`）——**建议单独立卡**；② `BRepOffsetSimpleOffset` 类体（`BRepOffset_SimpleOffset.cxx` L1-427）**仍未译**（现只有接口 + 同 GAP panic）。
7. **补 `Geom_Surface::UIso/VIso` 的 `Bezier` / `Offset` / `LinearExtrusion` 三臂**（OCCT 确有其 override ⇒ 真缺口；三臂之外的面型是 rcad 独有变体、OCCT 无对应）。另：`geom_convert_curve_to_bspline` 的两份重复应收敛（`GeomConvert::CurveToBSplineCurve` 只有一个 OCCT 实现）。
8. **哈希序（已量化降级，见坑 33）**：**不是**高优先级行为差异；只按热点核对，不搞全局容器替换。

> **★ 追加 27 更新（2026-09-13，勿按本节的 0 项取活）**：下面**第 0 项已改判**——实测失败点是 `LocOpe_Gluer::Perform` 的 `!the_split.is_done()`（`loc_ope_gluer.rs:545`），`LocOpe_Generator::Perform` **未被进入**；
> 真缺口 = **`LocOpe_SplitShape` 是 stub 类**（`feat/loc_ope_spliter.rs:287-384` ↔ OCCT `LocOpe_SplitShape.cxx` 1776 行，约 1,340 行未译，真身不在库）⇒ 下一手 = **翻译该全类**。详见 §0.0c 与 port-plan 追加 27。

**TKFeat（`feat/**` + `topalgo/int_curves_face_intersector.rs`）**
0. **★ 第 0 项（最急，a1 的直接下一墙，已定界到函数）**：**`LocOpe_Generator::Perform` 的 `IsDone`**（`feat/loc_ope_generator.rs` / `loc_ope_generator_b.rs`）。
   `geom_proj_lib::curve2d` 的 `Trimmed` 解包**已落地**（追加 17 补记 1）：`BRepFeat::IsInside` 恢复正常分类 ⇒ `LFPerform` 的粘合判定 `collage=true ope=Fuse`（与 OCCT 一致）⇒ 进入 `theOpe=1` 的粘合路径。**新墙**：`theGlue.Perform()` 后 rcad **`IsDone=false`** vs OCCT **`IsDone=1` / `resultFaces=9`**（面积 109.511）。OCCT 的 `myDone` 来自 `myDone = theGen.IsDone()`（`LocOpe_Gluer.cxx` **L230**，`theGen` = `LocOpe_Generator`）⇒ **下一手 = 查 `LocOpe_Generator::Perform` 为何不达**（对拍通道现成：`RCAD_WS_PROBE` + `output/build_tkfeat.bat`；OCCT 侧插桩点在 `LocOpe_Generator.cxx` 的 `Perform` / `myDone` 赋值处）。
1. **`Propagate` 的 `BOPAlgo_Section` 共面 FF**（`FalseSide ×5` 根因，追加 15 已定界）：`build_section` 逐字取 `PaveBlocksSc`，而该集合的唯一生产者是 **FF 曲线**；b3 实测 `sc_pb=0 sc_v=0` + DS 的 FF 记录 `curves=0 points=0` ⇒ 墙在**共面/同曲面面片对的 FF 分支**。**下一手必须是 OCCT 侧对拍**（在 `BOPAlgo_PaveFiller` 的 FF 段 / `IntTools_FaceFace` 同曲面分支插桩，跑 featlf b3），先拿"OCCT 在该对上产出了什么"再改 rcad。
2. **featprism 结果根**（`NoExtFace ×3` 的真因，④ 已证伪"欠计数"假说）：OCCT 侧 `nbshapes r` = **11 面 SOLID**（V12/E20/W12/F11/SHELL1/SOLID1），rcad 侧 `root_shape`（生成测试 helper）拿到 **3 面 Shell**。**先判定两条完全不同的修法**：是**结果池缺根**（`feat/brep_feat_make_prism.rs` / `loc_ope_*.rs` / `brep_feat_form*.rs` / `brep_feat_rib_slot*.rs`）还是**生成器 helper 取错**（`tools/occt-test-gen`）——判定依据 = 直接从结果 BRep 数出顶层 Solid 的面数是否等于 11。
3. **b6/b7**：`PtOnEdgeVertex`（`brep_feat_rib_slot_b.rs:1110-1115`）同一个 `my_sbase` 问题；D2 **已修**（追加 17 又补了 `basis_curve_of`），复测。
4. 其余既有队列（`NoFaceProf ×3`、`BRepTools_Modifier::Perform` GAP）见 §4.1；feat 侧实测失败地图见 §0.0。

> **★ 追加 27 更新（2026-09-13）**：`ChFi3d_ComputeArete`（OCCT `ChFi3d_Builder_0.cxx` L1984-2180）的**四个 stand-in 分支已全部 1:1 译全**——iso u/v 接线 `surface_uiso/surface_viso` 真身（并消灭了手写位置互换，坑 15）；`IFlag==0` 非 iso 补全 `ChFi3d_BuildPCurve` + Bnd_Box2d 事后检查 + `GeomLib::BuildCurve3d`（原实现丢 `pardeb/parfin` 与 `C3d`）；`else` 接线早已在库的 `chfi3d_project_pcurv`；补 OCCT L2000 的 `tolreached = tol3d;`。
> `ElCLib::AdjustPeriodic` 收敛为 `rcad_kernel::math::el::elclib_adjust_periodic` **唯一真身**（旧 fillet 简化版 + kernel `pub(crate)` 版双删）。**接力项 = ① `ChFi3d_SameParameter`（仍是无体 stand-in，`chfi3d_builder_0.rs:358-366`）② 补 `UIso/VIso` 的 Bezier/Offset 臂 ③ 既有 blend 剩余失败点。**

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


> **★ 追加 27 更新（2026-09-13）**：本节**第 2 项（入池，根）已落地**——`BRepTools_Quilt` + `BRepAlgo_FaceRestrictor` 产物入池（13 文件 / 9 调用点改签名；新增 `BRepBuilder::remove_from_compound`；`shell_registry` 镜像删除，池槽位即 OCCT 共享 TShape）。
> 设计判定：quilt **改收前导 `brep: &mut BRep`**、不自持池（自持池 ⇒ 产物按原 index 别名到消费者池 = 坑 2，比现状更糟）。
> **⇒ 接力项改为**：① 池外读取链**复核**（入池后 `curve_pool_free`/`edge_data_pool_free`/`build_curves3d.rs::edge_data` 等守卫应已无活路径，**勿重复立卡**）；② 第 1 项的"写回侧"随"池外无池可变"一并消失；③ 追加 26 第 1 项的 10 处 `UpdateCurves` 同族缺陷（`brep_fill_sweep_b.rs:234/:259` 优先）。

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

31. **（追加 19，最便宜也最容易自欺）不要从**被截断/分页**的失败地图清单里数数**：
    追加 19 正文的"`draft_angle` 库内 panic 25 → 20"就是**用一个 `| head -80` 截断的清单**数出来的（列表尾部被切掉，少算了 8 行），
    实测是 **28**，且本 session 五个批次**三次实测构成完全一致** ⇒ 那个"推进 5 例"的结论**是假的**。
    **做法**：一律 `sed -n '/^--- <grid>/,/^--- <next>/p'` 取**完整**分区后再 `grep -c` / `sort | uniq -c`；
    报数前把**分区总行数**与"通过数 + 失败数"对账（`draft_angle`：28 库内 + 21 断言 = 49）。
    另：**同一批次的"树内实测"要做两次以上**（本轮批次 1 后、3 后、5 后各一次）——只测一次容易把**基线漂移**当成自己的功劳。
    与坑 7 / 坑 8 同族，但这条是**操作纪律**而非环境问题。

32. **（追加 19）怀疑"是不是我改坏的"时，**先跑同树 5 次**，再决定要不要换树复测**（原写法已由坑 33 取代，见下）。
    本轮两处"疑似回归"最后都是**测量假象**：一处是记档漂移，一处是**同一份 exe 的运行间分裂**。
    **决定性的两个廉价手段**：① **同树复跑**（同一 exe 连跑 5 次，成本≈0）；② **换树复测**（`git checkout -b tmp <上一轮 tip>` + 根仓库重编 + `git checkout main` + 删临时分支）。
    **在拿到那两样之前，文档里只能写"待归因"，不能写结论**——本轮正文先误报"推进"，再用两次勘误修（坑 31/33），就是没有守住这条。

33. **（追加 19，最贵的方法学坑）rcad 的逐例失败地图**默认不稳定**，根源是 `std::collections::HashMap`/`HashSet` 的**每进程随机哈希种子**：
    **同一棵树、同一份已编译的 exe，连跑 5 次可以给出不同的"第一个失败点"**（本 session 实测 `offset_shape_type_a` a4：
    `brep_algo/image.rs:159` ×3 vs `brep_offset_make_offset_c.rs:52` ×2）。**这不止影响报数**——凡"多解/多候选取第一个"的
    路径，rcad 的行为本身就随进程漂移，而 OCCT 的 `NCollection_Sequence`/`IndexedDataMap` 是**插入序确定**的。
    **做法**：① 失败地图**连跑 ≥3 次**，出现过多个点就报**集合** + 分裂比例，不要报单点；
    ②**归因必须换树复测**（`git checkout -b tmp <上一轮 tip>` → 根仓库重编 → 跑同一张图 → `git checkout main` + 删分支），
    不能拿两次单跑的差异当因果；③ 同树复跑成本≈0，**没有任何理由省这一步**；
    ④ **★ 已量化（追加 20 收尾，据此降级）**：对 **11 个域网格 × 3 次连跑 = 172 个失败用例**统计"逐例失败点集合"，**只有 1 例**
    （`offset_shape_type_a` 的 a4）出现过两个失败点，**其余 171 例三次完全一致**；且 **4 个在域全绿网格的 48 例 × 5 次连跑零抖动**、
    八网格 ~8 次全跑计数一致。⇒ **哈希序今天只表现为"某一例先报哪个失败点"，没有任何证据表明它改变过**结果**。**
    **它是潜在隐患（OCCT 用插入序容器，长期值得对齐），但不是"高优先级行为差异"，更不要为此做大规模容器替换重构**——
    先按 OCCT 的容器类型在**具体热点**上核对，不搞全局换代。

## 7. 资产位置（本轮新增/更新）

| 资产 | 路径 |
|------|------|
| **追加 19 新增**：`ProjLib_HCompProjectedCurve` 三处接线（metrics 比较 / ProjectCurve / ProjectedCurve） | `libs/rcad-algo/src/geomalgo/geomplate/build_plate_surface.rs`（metrics / 非线性回退）+ `build_plate_surface_b.rs::{project_curve, projected_curve}`；真身 `geomalgo/proj_lib_h_comp_projected_curve.rs` + kernel `base/proj_lib/adaptor.rs::CurveOnSurface` |
| **追加 19 新增**：池外 `BRep_Tool::Curve`（`curve_pool_free` + `shape_is_in_pool`）与**续链** `edge_data_pool_free` | `libs/rcad-kernel/src/topo/topods.rs`（紧随 `curve_on_surface_pool_free`）；消费点 `topalgo/brep_lib/build_curves3d.rs::{brep_tool_curve, edge_data}`、`shhealing/shape_upgrade/unify_same_domain/topexp.rs::{brep_tool_curve_loc, brep_tool_range}` |
| **追加 20 新增**：`BRepTools_Modifier` + `BRepTools_Modification` 真身（~2,170 行；**消费者未接线**） | `libs/rcad-algo/src/topalgo/brep_tools_modifier.rs`（驱动体）+ `brep_tools_modification.rs`（接口 + `TrsfModification`） |
| **追加 20 新增**：`Geom_Surface::UIso/VIso` **唯一真身**（八臂并集） | `libs/rcad-algo/src/geomalgo/geom_surface_iso.rs`；消费者 `approx_curve_on_surface.rs` / `geom_lib_iso_line.rs` 改为委托（本地副本已删） |
| **追加 20 新增**：确定性测量脚本（可复用） | `output/determinism_check.sh`（域网格 ×3 次连跑 → 逐例失败点集合，报告多成员用例）；原始地图在 `output/determinism_maps/` |
| **追加 19 新增**：blend 链两处 1:1 修复 | `libs/rcad-kernel/src/geom/bspline_ops.rs::bsplclib_resolution`（有符号 clamp）、`libs/rcad-algo/src/fillet/brep_blend_walking.rs::add_singular_point`（1-based 读取） |
| **★ 追加 19 立卡 → 追加 20 已完成**：`Geom_Surface::UIso/VIso` 原为 rcad **两份**偏译 | **已合并为唯一真身** `geomalgo/geom_surface_iso.rs`（八臂并集），两个消费者改为委托（原卡：OCCT 里 `UIso/VIso` 是**单一虚函数**，不是 `isIsoLine`/`buildC3dOnIsoLine` 那种"OCCT 自己也有两份"的静态）。仍缺 `Bezier`/`Offset`/`LinearExtrusion` 三臂 |
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
