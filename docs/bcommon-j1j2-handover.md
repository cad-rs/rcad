# bcommon J1/J2 修复交接（2026-09-07 session 19 终）

> **一句话任务**：修复 `bcommon_simple J1/J2` 布尔失败——严格 1:1 对齐 OCCT（静态逐行审查，禁动态调试探针），验收 = 生成测试 j1/j2 转绿 + 布尔网格全绿 + 既有全套回归不破。

## 失败事实（全部已实证）

### J1
```tcl
pcylinder c1 20 100
pcylinder c2 20 100
ttranslate c2 0 0 50
bcommon result c1 c2
checkprops result -s 8796.46
```
- 配置：**同径同轴圆柱沿轴平移 50**——两侧面是同一圆柱面（same-domain 面）。
- OCCT 期望：common = c1 裁到 z∈[50,100]（侧膜 6283.19 + 双平面帽 2×1256.64 = 8796.46）。
- rcad 实际：**15079.644737231005 = 完整 c1 的面积**——侧面未被 z=50 圆分割。

### J1 拓扑 diff（step-topo-diff，已实测）
```
OCCT ref:  3 faces, 3 edges, 2 vertices; MANIFOLD_SOLID_BREP=1, CLOSED_SHELL=1
           EDGE_CURVE=3; 平面帽各 1 条边（整圆）
rcad:      3 faces, 2 edges, 2 vertices; OPEN_SHELL=3, MANIFOLD_SOLID_BREP=0
           平面帽 0 条边环（无边界！）；缺 1 条边（侧面 seam）
           圆柱面 4e 邻接签名不同
```
结构解读：**面存在但边界/环数据丢失**——平面帽没有 wire 环，侧面没有 seam、未被裁剪，壳不闭合。面积 15079 = 未裁侧膜全面积 + 双平面（无环面按曲面界算）。

### J2
```tcl
box b 10 10 10
copy b c
pcylinder s 2 4
ttranslate s 5 5 -2
bcut rr c s
explode rr so
bcommon result rr_1 c
checkprops result -s 625.133
```
- 配置：**rr_1（带孔盒）⊂ c（原盒）**，5 个盒面同域 + copy 形状。
- OCCT 真值（bool-runner 实测）：Area 625.133，V=60 E=30 W=9 F=8 SOLID=1，CYLINDER(1)+PLANE(7)。
- rcad 失败：**VERTEX 12 vs 参考 JSON 10**（多 2 顶点；参考 JSON = tests/occt/step_reference/occt_boolean_bcommon_simple_j2.json，V10/E?/F8）。
- 注意：j2 测试在拓扑断言处 panic，STEP 导出（maybe_export_step）在其之后 → 未导出；若需 STEP，把拓扑断言块临时让路或用 RCAD_STEP_ONLY（不行——拓扑断言不受 STEP_ONLY 保护）。

### 影响面
bcommon_simple 网格 166 测试中仅 j1/j2 两个 `draw_script_rcad_equivalent` 失败（j3..j9 及其余全绿）。Volumemaker/其他网格未复测（生成文件已按 --merge-groups 重生成在 tests/occt/tests/，gitignored）。

## 二分结论（勿再二分）

J1 在 **4f638749（08-31，最早可用锚点）到 HEAD 的每一个提交上都逐位同样失败**（15079.644737231005）。ef78ceb3/070c8ee8/98ee4daf/3839b351/3a1be962/65c43d63/013af8fd/686f7a55 全测过，全坏。

⇒ **不是近期回归**。09-01 基线"bcommon 83/83"之所以绿，是当时生成的 harness 断言集与现在不同（旧生成文件已随 gitignore 丢失，无法直接比对；当前生成器对 bcommon 发出 `-s` 面积断言 + nbshapes 拓扑断言）。即：**这是长期存在的缺陷，被重新生成的 harness 首次暴露**。

## 静态审计入口（新 session 直接开工）

失败结构 = DS→BRep **面/环装配**阶段：面有、环无。按 OCCT 1:1 对齐顺序审计：

1. **rcad** `libs/rcad-algo/src/bop/algo/builder.rs`（BooleanBuilder：BuildFace/wire 装配/面重建）＋子模块——对照 **OCCT `BOPAlgo_BuilderFace.cxx`**（Perform/PerformLoops/PerformAreas/PerformSDFaces——同域面用对方边分割）与 `BOPAlgo_Builder.cxx`（BuildShape/BuildSolid）。
2. **同域面链**：J1 侧面同域 → PaveFiller FF 阶段（`bop/pavefiller*` / `geomalgo/int_patch/imp_prm_intersection.rs`）的 same-domain 检测与 SD 连接 → `bop/ds/common_block.rs`（BOPDS_CommonBlock 已 1:1）→ BuilderFace 消费 CommonBlock pave blocks 分割。
3. **J2 特有**：`copy`（同一 box 的两个形状）+ explode so——检查 copy 的 TShape/位置语义与 rr_1/c 的同域面匹配。
4. 关键症状映射：平面帽 0 边环 = 帽面的 wire 未从 pave blocks/Loops 装配；缺 seam = 侧面 wire 未闭合；OPEN_SHELL = 面间共享边未共 TShape（extract/assembly 阶段 identity 断裂——注意 `98ee4daf` 提交"extract_solids shares TShapes"与此相关，但其时间点已排除）。

## 复现命令（根仓库）

```bash
# 失败复现（j1 面积断言）
cargo test -p occt-generated-tests --test generated_occt_boolean_bcommon_simple j1_draw_script
# STEP 导出（j2 拓扑断言在导出前 panic，需注释断言或先行处理）
RCAD_STEP_DIR="tests/occt/step_output" RCAD_STEP_ONLY=1 cargo test -p occt-generated-tests --test generated_occt_boolean_bcommon_simple j1_draw_script
# 拓扑对比
PATH="/c/tools/opencascade-8.0.0-vc14-64/win64/vc14/bin:$PATH" \
STEP_TOPOLOGY_DUMP="tools/step-topo-dump/build/Release/step_topology_dump.exe" \
tools/step-topo-diff/target/release/step-topo-diff.exe \
  tests/occt/step_output/occt_boolean_bcommon_simple_j1.step \
  tests/occt/step_output/ref/bcommon_simple_J1.step
# OCCT 真值（J2 已在册；J1 未注册可仿 cases/*.hpp 添加）
PATH="/c/tools/occt-debug/win64/vc14/bind:$PATH" tools/occt-bool-runner/build/Debug/occt_bool_runner.exe bcommon_simple J2
```

## 纪律（沿袭）

- **静态逐行对比修复**，禁运行时探针/println；每改一处 `cargo check -p rcad-algo`；不简化、不自创、等价替换=没对齐。
- 布尔网格需 `--merge-groups` 生成：`cargo run -p occt-test-gen -- --batch-boolean --merge-groups`（生成文件 gitignored）。
- 布尔网格基线作废重录：09-01 的"83/83"等数字按当前 harness 首次复测后重建。
- 冻结文件勿动：`algo_ext/mod.rs`（TEMP-EXCLUDED 注释版）、`fillet/topopebrepbuild.rs`（WIP）、`tests/tkgeom_algo_gtests.rs`——当前工作树已恢复为冻结态（未提交），**不要 git add/checkout 它们**。
- rcad-render 的 rcad_algorithms 导入损坏为预存问题（bc0bfa85 删库），不在范围。
- TKHLR exact 线路已于本 session 收口（box/bug25813_1/ptorus 三验收全绿，提交 ce92a09d 及之前），勿回退。

## 结案记录（2026-09-07 后续 session）

**j1/j2 已全绿（bcommon_simple 83/83）。根因不在布尔管线——本档此前的定性
（"DS→BRep 环装配缺陷"）是错的；缺陷在测试生成器。**

- **真根因**：根仓库 `tools/occt-test-gen/src/main.rs` 的 draw-script
  `ttranslate` 圆柱折叠分支先把平移写进 `*base`，但重发的
  `make_cylinder_brep` 行用了**旧 base 值**——平移被静默丢弃。j1 变成
  `common(c1, c2@原点)`（两圆柱完全重合），结果=完整 c1，面积
  15079.644737231005 = 12566.37 + 2×1256.64 逐位吻合；j2 的孔打在盒角而非
  (5,5,-2) 中心（V12 vs V10）。这也解释了"每个提交逐位同败"——生成器对
  所有提交产出同样的坏测试，二分结论在此意义上成立但指向了错误层。
- **修复**：发射折叠后的新 base（tools/occt-test-gen 根仓库提交）。
- **取证更正**：交接时引用的 step-topo-diff"平面帽 0 边环/OPEN_SHELL=3"
  来自过期 STEP 产物；当前代码导出的 j1 是结构完整的 c1（2V/3E/3F，
  计数断言恰好通过，仅面积断言暴露）。
- **volumemaker A1 为另一独立预存失败**（嵌套球 mkvolume→cut 返回未裁剪
  外球）：用旧生成器 stash 实测逐位相同失败，与本修复无关，待查。
- TKHLR(session16-19)对布尔管线的改动（Line3 清扫/FunctionSetRoot/
  Ellipse2d）经逐 commit 排查均与本问题无关；另发现 e65d5ee1 提交信息
  声称的 elclib 去归一化不在该提交的树中（elclib.rs 仍双重 normalize），
  属声称未落地，ulp 级，不影响本问题。
