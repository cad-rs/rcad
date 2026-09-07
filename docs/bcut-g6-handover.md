# bcut_simple G6 修复交接(2026-09-07 session 20)

> **一句话任务**:修复 bcut_simple G6(profile rev + revol 360° + bcut,checkprops -s 41187.4)
> 超时。方法论:严格 1:1 翻译 OCCT,不自创。已完成大半(PaveFiller 124s→<1s,
> 全部回归绿),剩余两个翻译单元在文末,按清单机械执行即可。

## session 21 增补(同日):A/B 两单元已完成

- **A 完成**:bean_face_intersector.rs 的 `IntCurveSurfaceHInter` 自创
  Phase1/Phase2/refine_crossing 采样投影全删,结构体改为持有
  `geomalgo::int_curve_surface::Intersection`(PARAMEQUAL 去重基类),
  `perform` 经 `hinter_adaptor.rs`(BRepAdaptorCurveTool/BRepAdaptorSurfaceTool
  两个 tool trait + IntCurveSurfaceHInter 的 HInterHost 实现)直通已 1:1 的
  inter_impl 引擎——Line×Revolution 走 PerformConicSurfLine L560-745 非解析
  分支(面窗口 polyhedron nbsu/nbsv>=20 + IntfTool::lin_box + 每段
  2 点 polygon + internal_perform 干涉 + ZerCSParFunc/IntCS 精化)。
  bean 的 distance_simple/distance_with_uv 同时改为 OCCT L397/L465 形式:
  `myContext->ProjPS(face)` 缓存项目器(context.rs 新增
  `ProjectOnSurface::new_init`;每 bean init 一次 10x10 种子网格,
  查询只做 Newton),不再每次调用重建 32x32 全网格。
- **B 完成**:kernel function_set_root.rs 补齐 OCCT `State = F.GetStateNumber()`
  全部调用点(trait 新增带默认实现的 `get_state_number`);新翻译
  `Extrema_FuncExtCS`(Extrema_FuncExtCS.cxx 全文,D1/D2 函数集 +
  GetStateNumber 去重记录解);ExtremaGenExtCS::perform 按
  Extrema_GenExtCS.cxx L299-423 重写(aNbIntC 闭合/周期分段循环、每段
  GlobMin* + math_FunctionSetRoot(myF, Tol)、filter_false_extrema 假极值
  过滤);自创 `cs_newton_refine` 数值 Hessian Newton 已删。
- **回归**:rcad lib 369/369、kernel 664/664、builder_stage 76+smoke 1、
  pavefiller_stage 26/26 全绿;布尔网格(子测试级,含 geometry_loads+
  draw_script 两条/case):bcommon 166/166、bcut 201 过+17i(g6 两条
  filtered)、bfuse 196/196+8i、bopcommon 755/755+1i、bopcut 729/729+29i、
  bopfuse 748/748+2i;boptuc 741 过/4 失败(ze7/ze8/ze9/zf1)、
  splitter 16 过/2 失败(a2/b2)——经 stash 反向验证,这 6 个失败在
  不含本 session 改动的工作树上同样失败,属冻结 fillet 线/已提交代码的
  既有失败,非本次引入。
- **遗留(下一单元)**:g6 仍超时,但慢点已从 EF 采样(本单元已删)移到
  主 PF 的 **EE 递归**(Line×Circle/BSpline 近交对)。诊断结论
  (RCAD_HANGDBG 临时插桩实测,已还原):FindSolutionsRec 对近切线
  对的 3 叉树达 ~10^5 节点 × ~0.7ms/节点(解析曲线点值+逐节点
  BndBuildBox/find_parameters walk),单对 ~100s。OCCT 同算法每节点
  仅 ~10μs(C++),故 OCCT 全程秒级。下一单元:1:1 审查
  IntTools_EdgeEdge 的 FindParameters 自适应步长(k 增长/aMaxDt 上限)、
  BndLib_Add3dCurve 弧盒与 Resolution 语义,找形式差异使树收敛。
  另:tkdata_remaining_gtests 编译失败(minor_dir)系 2850e5e7
  (Ellipse2d frame 修改)未同步测试所致,与本线无关,待单独修。

## G6 是什么

```
box b 100 100 40
profile rev S b_4 F 50 20 Y 50 C 10 180 Y -50 C 10 180
revol rev2 rev 0 0 50 0 1 0 360
bcut result b rev2        # checkprops result -s 41187.4
```
rev2 = 18 段轮廓(2 直线 + 2×180° 圆弧折线)绕 Y 轴 360° 旋转体,
每段 profile 直线扫成 `Surface3::Revolution`(线 profile 旋转面,
对应 OCCT `Geom_SurfaceOfRevolution`,自然 V 域无界)。

## 已完成(已提交,rcad 子仓库 e38 后续、根仓库指针已 bump)

1. **`surface_type` 映射**(bean_face_intersector.rs):Revolution→
   `GeomAbs_SurfaceOfRevolution`、LinearExtrusion→`SurfaceOfExtrusion`。
   原来折叠进 OtherSurface,错误启用 ComputeLocalized
   (IntTools_BeanFaceIntersector.cxx L327-338 只允许 Bezier/Other/BSpline)。
   这是 PaveFiller 124s 的根因。
2. **`BRepAdaptorSurface::with_uv_bounds`**:adaptor 参数域=面 UV rect
   (IntTools_Context::SurfaceAdaptor = BRepAdaptor_Surface(theFace),
   BRepTools::UVBounds);EdgeFace::new 与 context::compute_ef 按
   `ds.face_uv_boundary()` 初始化。
3. **`Extrema_GenExtCS` 翻译**(新文件 bop/int_tools/extrema_gen_ext_cs.rs):
   math_PSO/math_PSOParticlesPool/math_BullardGenerator(math_PSO.cxx
   L125-268)、GlobOptFuncCS value/gradient(GlobOptFuncCS.cxx L30-56)、
   GlobOptFuncConicS、GlobMinConicS/GlobMinGenCS/GlobMinCQuadric
   (Extrema_GenExtCS.cxx L427-823);ExtremaExtCS::perform 按
   Extrema_ExtCS.cxx 分发(二次曲面走原路径;其余 NbT=12,NbU=NbV=10,周期 13)。
4. **投影面域限制**:`DS::face_restricted_surface`、`ProjectOnSurface::init`
   按 GeomAPI Init(aS,U1,U2,V1,V2,Tol) 限制搜索域 + GeomGridEval 种子网格、
   Distance/distance_with_uv 用受限 adaptor。

验证:lib 369 / pavefiller 26 / builder 76+1 / kernel 664 全绿;
17 布尔网格中 16 个与基线一致(volumemaker A1 为另一预存失败);
bcut 109 过 + g6 超时(外部行为与基线一致)。

## 剩余翻译单元(按序执行)

### A. HInter 非解析路径(IntCurveSurfaceHInter::perform 的采样回退=自创,必须替换)

OCCT 真身:IntCurveSurface_Inter.pxx PerformConicSurfLine(L560-745)
Line×SurfaceOfRevolution → 精确分支全部 fall through →:

1. `nbsu = NbSamplesU(S,U1,U2); nbsv = NbSamplesV(S,V1,V2)`(≥1);
   UV 无界时 Revolution 走 `EstLimForInfRevl`(IntCurveSurface_InterUtils.pxx)
   估有限界(rcad 用面窗口 adaptor 时 UV 已有限,此步自然跳过);
   `nbsu = max(nbsu,20); nbsv = max(nbsv,20)`。
2. `PolyhedronType polyhedron(S, nbsu, nbsv, U1, V1, U2, V2)`
   — IntCurveSurface_ThePolyhedronOfHInter.cxx:UV 网格点 + 每格两三角形,
   顶点存 (point, u, v);`Bounding()` 总盒。
3. `Intf_Tool bndTool; bndTool.LinBox(line, poly.Bounding(), boxLine)`
   — 直线∩总盒 → [pinf, psup] 段(slab 法,<1e-10 时扩 1e-10)。
4. 每段:`polygon(curve, pinf, psup, )`(2 点)→
   `IntCurveSurface_TheInterferenceOfHInter`(Intf_InterferencePolygonPolyhedron.gxx):
   段×三角形盒过滤 → 段-三角形求交 → Intf_SectionPoint(t, u, v)。
5. 每个 SectionPoint 经 `IntCurveSurface_TheExactHInter`
   (IntImp_IntCS.gxx L26-153)精化:
   - 函数 IntImp_ZerCSParFunc.gxx:F = Psurf(u,v) - Pcurv(w)(3 方程),
     雅可比 [Du, Dv, -Dw],Root()=|F|²。
   - `math_FunctionSetRoot Rsnld(F); Rsnld.Perform(F, (U,V,W), BornInf, BornSup)`,
     Tolerance=(UResolution, VResolution, CurveResolution)(Confusion);
     失败重试 w0/w1(≤3 次);接受条件 |Root| ≤ TolTangency²(≥SquareConfusion)。
   - **需要 N 变量 math_FunctionSetRoot**:rcad 现有 function_set_root.rs 是
     2 变量特化(IntPatch 用)。按 math_FunctionSetRoot.cxx(OCCT
     FoundationClasses/TKMath/math,~600 行)把同一算法泛化为 N 变量
     (向量用 Vec<[f64;N]> 或小栈数组),或并列写 3 变量实例
     `function_set_root3.rs`(同一 .cxx 逐行,N=3)。
6. 段的两个入/出点构成 IntersectionSegment(线在面内的穿越段)。

替换范围:bean_face_intersector.rs `IntCurveSurfaceHInter::perform` 的
`_ => Vec::new()`(非解析)分支之后的全部采样代码(Phase1/Phase2/
refine_crossing)。二次曲面精确分支(Line×Plane/Cylinder/Sphere/Cone/Torus)
保持不动(已 1:1)。

### B. PostTreatFF 嵌套子 PaveFiller 的代价(g6 挂死的表现层)

实证(pave_filler_make_blocks.rs 加临时步进打印实测,已还原):
- 主 PF 全部阶段秒级;挂死在 post_treat_ff 内的 `a_pf.perform(&a_ps)`
  (子 PaveFiller)——OCCT 同样有该子 PF(PostTreatFF L1389-1397),
  机制忠实;挂是因为子 PF 的 EF 又走 18ms/对 的采样路径(A 完成后自然消失)。
- 另 VF=4.7s:compute_vf 的 proj_ps 已限域+网格缓存,若仍慢,检查
  `put_paves_on_curve`/VertexFace 是否绕过 ProjectOnSurface。
- 遗留自创(必须随 A 一起删):HInter 的 Phase1/Phase2/refine_crossing
  采样;`cs_newton_refine`(被 GenExtCS 使用,A5 完成后换 N 变量 FunctionSetRoot)。

### C. 复现与验收

```bash
# 挂死复现(根仓库)
timeout 300 cargo test -p occt-generated-tests --test generated_occt_boolean_bcut_simple g6 -- --nocapture
# 逐例验收(全绿目标 110/110)
EXE=$(ls -t target/debug/deps/generated_occt_boolean_bcut_simple-*.exe | head -1)
"$EXE" --list | grep draw_script | ... # 逐例 --exact + timeout 20
# 回归: lib 369 / pavefiller 26 / builder 76+1 / kernel 664 / 17 网格
```

### D. 纪律

- 全部 OCCT 源:TKGeomAlgo/IntCurveSurface(IntCurveSurface_Inter.pxx
  L560-745、PolyhedronUtils.pxx、ThePolyhedronOfHInter.cxx)、TKGeomAlgo/IntImp
  (IntImp_IntCS.gxx、IntImp_ZerCSParFunc.gxx)、TKMath/math
  (math_FunctionSetRoot.cxx、Intf 目录:Intf_InterferencePolygonPolyhedron.gxx
  等)。逐行 1:1,每函数标注行号,每处 cargo check。
- 冻结:fillet 线(topopebrepbuild.rs、algo_ext/mod.rs 的 TEMP-EXCLUDED、
  tkgeom_algo_gtests.rs)留待后续推进,勿动。
