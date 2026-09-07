# bcut_simple G6 修复交接(2026-09-07 session 20)

> **一句话任务**:修复 bcut_simple G6(profile rev + revol 360° + bcut,checkprops -s 41187.4)
> 超时。方法论:严格 1:1 翻译 OCCT,不自创。已完成大半(PaveFiller 124s→<1s,
> 全部回归绿),剩余两个翻译单元在文末,按清单机械执行即可。

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
