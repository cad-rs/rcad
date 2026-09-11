# E3-X 交接：g6 已通过 · 八网格 8/8 全绿（2026-09-12）

> **一句话**：bcut_simple g6 于 2026-09-11/12 通过（拓扑 10F/27E/18V 与 OCCT 全等 + SA
> 41187.4），八网格首次全绿。本轮（追加 9+10）共 4 处 1:1 修复、1 个 ~1900 行移植单元。
> 下一项 = E3-V 队列第 3 项（feat/offset/blend 各域）。

## 0. 当前门槛（全部实测，2026-09-12 凌晨）

| 基线 | 值 |
|------|-----|
| rcad-algo lib | **412/0/0** |
| rcad-kernel lib | **678/0** |
| pavefiller_stage_tests | **26/26** |
| builder_stage_tests | **76/76** |
| builder_stage_smoke | **1/1** |
| bopfuse / bopcommon / bopcut / boptuc | 375 / 378 / 379 / 373 |
| splitter / bfuse_simple / bcommon_simple | 12 / 102 / 83 |
| **bcut_simple** | **110/110（g6 已通过）** |

## 1. 提交链（本次会话，均未推送）

- 根仓库 `c96a0e3`（sync: rcad 指针 a5fced7c）← `b90eead`（occt-bool-runner 新增
  bcut_simple **G6 用例**——脱离 DRAWEXE 跑真 g6 几何，debug OCCT 无 Draw 模块时的
  唯一通道，**保留复用**）
- rcad `a5fced7c`（**追加 10 核心**：rst_int.rs 移植 + imp_prm 接线 + PerformEF
  懒初始化 + 回转圆柱 frame Y）← `2513a6dd`（**追加 9**：ProjLib_Plane::Project
  in-plane BSpline）← `e7472c51`（追加 8 文档）
- OCCT 源探针已全部 `git checkout` 还原（TKBO WireSplitter、TKPrim Rotation 均 0 残留）

## 2. g6 三层真因链与修复（细节见 tkfeat-fillet-offset-port-plan.md §E3-W 追加 9/10）

1. **`IntPatch_RstInt::PutVertexOnLine` 缺失** → 新文件
   `libs/rcad-algo/src/geomalgo/int_patch/rst_int.rs`（~1900 行 1:1，
   IntPatch_RstInt.cxx L433-1355 + PolyLine/PolyArc/CSFunction/矩形域 TopolTool +
   全部静态函数），接线于 `imp_prm_intersection.rs` 的 L1768-1774 循环。
   - 唯一 GAP：`FindParameter` 的 Restriction 分支（需 `IntPatch_HInterTool::Project`），
     返回 false = OCCT 失败路径，待出现含 RLine 的 ImpPrm 用例再补。
   - `IntPatchVertex` 新增 `vtx_on_s1/vtx_on_s2: Option<(u16,u16)>`（域顶点身份）。
2. **PerformEF 只读 PB 访问器** → `pave_filler.rs` 改 `change_pave_blocks(n_e)`
   （OCCT ChangePaveBlocks 懒初始化语义，PaveFiller_5.cxx L245）。
3. **回转圆柱 frame 丢 Y** → `tool_rehost.rs::cylinder()` 改 `y_dir: Some(axe_rev.3)`
   （OCCT Load 的 YReverse 后 Ax3 为左手系，u 与被扫圆同相；`CylindricalSurface.y_dir`
   字段本为此场景而设）。
4. （追加 9）**`build_projected_pcurve` 缺 ProjLib_Plane::Project 分支** →
   `face_make_curve.rs` 对"BSpline 落在平面内"走精确 pole 映射
   （kernel `PlaneProjector::project_bspline`），消除 ComputeTolerance 23.33 毒化。

## 3. 下一步队列（按优先级，来自 E3-W 追加 10）

1. **E3-V 队列第 3 项照旧**：feat 域 BRepBuilderAPI_Transform/Trsf::set_rotation/
   BRepCheck_Analyzer/Geom2dAPI_ExtCC2d；offset 域 GeomAPI_ProjectPointOnCurve 真身 +
   ExtentEdge 面盒延伸 + GeomInt_IntSS 专批（~1,860 行）；blend a1 上游
   ChFiKPart_ComputeData::Compute；a1/TKOffset TrimEdges。
2. rst_int.rs 的 FindParameter Restriction 分支（IntPatch_HInterTool::Project）。
3. StepWriter 周期面 seam 保真度（rev2 STEP 9 vs 8 EDGE_CURVE，追加 5 记录）。
4. 清理照旧：`rcad-kernel/src/math/gprop/` 死目录、`face_classifier.rs` 的 `[FC]`
   eprintln、`CurveToolHandle::for_curve3` resolution=1.0 占位、
   `normal_projection.rs` BRepLibMakeWire 载体切真身、
   `tools/occt-test-gen/tests/` 历史遗留跟踪文件、prism 构造迁移至
   `rcad-modeling/prim/primapi/`。

## 4. 操作配方（复用，勿重摸）

```bash
# 分阶段基线（每次改动后）
cd /c/Users/lilu/works/rcad-pro/rcad
cargo test -p rcad-algo --lib                 # 412/0/0
cargo test -p rcad-kernel --lib               # 678/0
cargo test -p rcad-algo --test pavefiller_stage_tests   # 26/26
cargo test -p rcad-algo --test builder_stage_tests      # 76/76
cargo test -p rcad-algo --test builder_stage_smoke      # 1/1

# 八网格（exe 必须取根仓库 target/debug/deps 最新产物；tests/occt/target 是旧产物陷阱）
cd /c/Users/lilu/works/rcad-pro && bash output/run_eight_grids.sh

# g6 单例（含 STEP 导出 + 拓扑对比）
cargo test --no-run -p occt-generated-tests --test generated_occt_boolean_bcut_simple
EXE=$(ls -t target/debug/deps/generated_occt_boolean_bcut_simple-*.exe | head -1)
RCAD_STEP_DIR=tests/occt/step_output RCAD_STEP_ONLY=1 "$EXE" \
  --exact occt_boolean_bcut_simple_g6::occt_boolean_bcut_simple_g6_draw_script_rcad_equivalent
STEP_TOPOLOGY_DUMP=tools/step-topo-dump/build/Release/step_topology_dump.exe \
PATH="/c/tools/opencascade-8.0.0-vc14-64/win64/vc14/bin:$PATH" \
  tools/step-topo-diff/target/debug/step-topo-diff.exe \
  tests/occt/step_output/occt_boolean_bcut_simple_g6.step \
  tests/occt/step_output/ref/bcut_simple_g6.step

# OCCT 侧插桩对拍（最速定界手段）
# 1) 改 C:/Users/lilu/works/OCCT/src/...（必须用 Edit 工具，Bash 写 OCCT 源被钩子拦）
# 2) cmd.exe //c "C:\Users\lilu\works\rcad-pro\output\build_tkprim.bat"   # TKBO 同理换 target
# 3) cp -f C:/tools/occt-debug/win64/vc14/bind/<DLL>.dll tools/occt-bool-runner/build/Debug/
#    grep -ac <探针串> tools/occt-bool-runner/build/Debug/<DLL>.dll       # 验证 >0
# 4) cd tools/occt-bool-runner/build/Debug && RCAD_WS_PROBE=1 ./occt_bool_runner.exe bcut_simple G6
# 5) 用后 git checkout -- 还原 OCCT 源（探针即用即清）
```

## 5. 本轮固化坑（勿重复）

- **gp_Ax3 的 Y 方向陷阱**：`gp_Ax3(P,N,Vx)` 在 N 被反转后其 YDirection ≠ N×Vx——
  GeomAdaptor Load 的 `YReverse()` 把它翻回"按未反转轴计算"。凡 rcad 以
  `y_dir: None`（重算 `axis×ref_dir`）消费这类 frame 会引入 π 相位差。
  `CylindricalSurface.y_dir: Some(...)` 是既定通路；**cone/sphere/torus 结构体尚无
  此字段**，遇翻转轴用例需同样扩展。
- **OCCT `Change*` 访问器带懒初始化副作用**：换成只读访问器会让"从未被触碰的形状"
  （孤立闭合圆等）静默跳过整个阶段。
- **探针必须 `--nocapture`** 才可见；提交前 `git diff | grep "+.*eprintln"` 应为 0。
- 后台 cargo 编译期间**不要并行改源码**（本轮 staged tests 被探针清理竞态打断过一次，
  重跑即好）。
- `temp/` 未跟踪、`output/` gitignored——探针输出/临时 bat 一律放这两处。
