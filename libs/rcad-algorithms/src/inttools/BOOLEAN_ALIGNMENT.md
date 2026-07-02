# OCCT Boolean Pipeline → rcad 对齐状态

每条记录：OCCT函数（源文件:行号）→ rcad函数（文件:行号）· 状态标记

## 第 0 层：顶层调度

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `BRepAlgoAPI_Fuse::Perform` | `boolean_op_with_retry` (lib.rs) | ✅ | 架构等价 |
| `BOPAlgo_BOP::Perform` | `bop_occt_union::fuse` (bop_occt_union.rs) | ⏳ | rcad 跳过了部分 Glue/非 Destructive 分支 |
| `BOPAlgo_BOP::PerformInternal1` | `BooleanBuilder::build_with_history` (builder.rs:795) | ⏳ | 步骤结构相同，但 rcad 合并了部分 BuildResult 步骤 |

## 第 1 层：PaveFiller — 交线计算与 DS 填充

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `BOPAlgo_PaveFiller::Perform` | `PaveFiller::perform` (mod.rs) | ⏳ | 整体流程一致，部分检查/容差计算不同 |
| `BOPAlgo_PaveFiller::PerformEE` | EE 交线计算 (edge_edge.rs) | ⏳ | 结构等价，容差常数不同 |
| `BOPAlgo_PaveFiller::PerformEF` | EF 交线 (edge_face.rs) | ⏳ | 部分对齐 |
| `BOPAlgo_PaveFiller::PerformVF` | VF 交线 | ❌ | rcad 有 unique 逻辑 |
| `BOPAlgo_PaveFiller::PerformVV` | VV → PairIterator | ✅ | 架构 A3 已修复 |
| `BOPAlgo_PaveFiller::MakeBlocks` | `make_blocks` (make_blocks.rs:15) | ⏳ | 循环+变量结构与 OCCT L725-1107 大致对齐，但细节差异多 |
| `IntTools_Context::IsVertexOnLine` | `is_vertex_on_line` (paves.rs:449) | ✅ | **刚对齐** |

### PutPaveOnCurve 系列

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `PutPaveOnCurve` (单vertex) | `put_pave_on_curve` (paves.rs:244) | ⏳ | OCCT L2950+ 有 extended tolerance 逻辑；rcad 省略了部分分支 |
| `PutPavesOnCurve` (批量) | `put_paves_on_curve` (paves.rs:314) | ⏳ | 结构类似 |
| `PutBoundPaveOnCurve` | 内联在 make_blocks.rs:290-379 | ⏳ | 逻辑分散，没有独立函数 |
| `PutStickPavesOnCurve` | `put_stick_paves_on_curve` (paves.rs) | ⏳ | 未逐行对齐 |
| `FilterPavesOnCurves` | `filter_paves_on_curves` (paves.rs) | ⏳ | 部分对齐 |

### Face-Face 交线

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `IntTools_FaceFace::Perform` | `ff_intersect::perform` (ff_intersect.rs) | ⏳ | 整体流程一致 |
| `IntPatch_Intersection` (球-平面) | `intersect_plane_sphere_faces` (analytic_plane.rs:411) | ⏳ | 结果处理不同 |
| `IntPatch_Intersection` (柱-平面) | `intersect_plane_cylinder_faces` (analytic_plane.rs) | ⏳ |  |
| `IntPatch_Intersection` (球-球) | `analytic_sphere.rs` | ⏳ |  |
| `ProjLib::MakePCurveOfType` (平面) | `circle_pcurve_on_plane` (pcurve_derive.rs:31) | ⏳ | BSpline 分支 knot rescale 已修；Circle2d 分支对齐 ✅ |
| `ProjLib::MakePCurveOfType` (球) | `circle_pcurve_on_sphere` (pcurve_derive.rs:130) | ❌→✅ | **刚修**：knot rescale + isIsoU/isIsoV 标注（Line2d 路径因缺 TrimmedCurve2 裁剪暂未启用） |
| `ProjLib::MakePCurveOfType` (柱) | `circle_pcurve_on_cylinder` (pcurve_derive.rs:184) | ⏳ | knot rescale 已修 |
| `ProjLib_Sphere::Project(gp_Lin)` | 无 | ❌ | 线-球投影未实现 |
| `ProjLib_Sphere::Project(gp_Elips)` | 无 | ❌ | 椭圆-球投影 fallthrough 到通用 |

### PaveFiller 辅助函数

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `SubShapesOnIn` | 内联在 make_blocks.rs:124-166 | ⏳ | 结构不同 |
| `SharedEdges` | 内联在 make_blocks.rs:168-176 | ✅ |  |
| `GetStickVertices` | 内联在 make_blocks.rs:181-237 | ⏳ | 过滤条件差异 |
| `IsValidBlockForFaces` | 内联在 make_blocks.rs:438-471 | ⏳ | pcurve 求值/容差不同 |
| `FindValidRange` | `find_valid_range` (make_blocks.rs) | ⏳ |  |
| `PutClosingPaveOnCurve` | `put_closing_pave_on_curve` (make_blocks.rs) | ⏳ |  |
| `CorrectTolerances` | `correct_tolerances` (boptools/extra.rs) | ⏳ |  |

## 第 2 层：Builder — 结果构建

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `BOPAlgo_Builder::PerformInternal1` | `BooleanBuilder::build_with_history` (builder.rs:795) | ⏳ | 步骤顺序一致，部分跳过(TreatDim, Refine) |
| `FillImagesVertices` | `fill_images_vertices` (filler.rs:11) | ✅ | 简短的 SD 登记，已对齐 |
| `FillImagesEdges` | `fill_images_edges` (filler.rs:49) | ⏳ | 逻辑相同，DS 索引方案差异 |
| `FillImagesFaces` | `fill_images_faces` (filler.rs:181) | ⏳ | **待检查** |
| `FillImagesContainers(WIRE)` | `fill_images_containers_wires` (filler.rs:96) | ⏳ | 逻辑等价，wire 表示不同 |
| `FillImagesContainers(SHELL)` | `fill_images_containers_shells` (filler.rs) | ⏳ |  |
| `FillImagesSolids` | `fill_images_solids` (filler.rs) | ⏳ |  |
| `BuildResult(VERTEX)` | `build_result_occt(Vertex)` (result_build.rs:1142) | ⏳ | |
| `BuildResult(EDGE)` | `build_result_occt(Edge)` | ⏳ | |
| `BuildResult(FACE)` | `build_result_occt(Face)` | ⏳ | rcad 在 fill_images_faces 就直接创建 TShape，不是 batch |
| `BuildResult(SHELL)` | `build_result_occt(Shell)` | ⏳ | |
| `BuildResult(SOLID)` | `build_result_occt(Solid)` | ⏳ | |
| `BuildShape` | `build_shape` (result_build.rs:1159) | ⏳ | |
| `PostTreat` | `post_treat` (result_build.rs:1169) | ⏳ | 只做 tolerance correction |
| `MakePCurve` (face 边) | `build_face_reps` / edge_builders.rs | ⏳ | **待检查** |

## 第 3 层：BOPDS — 数据结构

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `BOPDS_DS` | `DS` (bopds/ds/types.rs) | ⏳ | 架构不同(Rust style)但功能等价 |
| `BOPDS_Curve` | `IntersectionCurve` (types.rs:360) | ✅ | 字段对齐 |
| `BOPDS_PaveBlock` | `PaveBlock` (pave.rs) | ✅ | |
| `BOPDS_Pave` | `Pave` (pave.rs) | ✅ | |
| `BOPDS_ShapeInfo` | `ShapeInfo` (types.rs) | ✅ | 标记系统已补(A4) |
| `BOPDS_Iterator` | `PairIterator`/Bvh | ⏳ | BVH 替代 O(n²) |
| `BOPDS_SubIterator` | 类似 | ⏳ | |
| `BRep_Builder` 增量构建 | `ResultBuilder` + `build_edges()` | ✅ | A1 已修 |

## 第 4 层：辅助工具

| OCCT | rcad | 状态 | 差异 |
|------|------|------|------|
| `Geom2d_TrimmedCurve` | `TrimmedCurve2` (geom.rs:740) | ✅ | 结构对齐 |
| `Geom2d_BSplineCurve` | `BSplineCurve2` (geom.rs:702) | ✅ | |
| `GeomLib::SameRange` | `same_range_2d` (boptools/extra.rs:262) | ✅ | 逐分支对齐：Line(translate), Circle(rotate frame), Trimmed(recurse), BSpline(reparam knots) |
| `gp_Circ2d` (坐标方向轴) | `Circle2d { x_dir, y_dir }` (geom.rs:628) | ✅ | **刚补齐** — 新增 x_dir/y_dir + `rotate_center()`，`point_at` 使用定向参数化 |
| `BOPTools_AlgoTools::CorrectTolerances` | `correct_tolerances` (tolerance.rs) | ⏳ | |
| `IntTools_Tools::ComputeTolerance` | `estimate_pcurve_deviation` (boptools/extra.rs:336) | ✅ | 25点均匀采样 + 1.00001 margin，对齐 OCCT L737-779 |
| `IntTools_Context::ProjectPointOnEdge` | `closest_point_on_curve` (projection.rs) | ⏳ | 功能等价但实现不同 |
| `IntTools_Context::ProjPT` (缓存投影器) | 无 | ❌ | rcad 每次新建，不缓存 |
| `Extrema_LocateExtPC` | `closest_point_on_curve` 局部模式 (64 iter Newton) | ⏳ | 功能等价但形式不对 |
| `Extrema_ExtPC` | 无独立函数 | ❌ | rcad 用单一 closest_point_on_curve 替代两级投影 |

## 关键未对齐点总结

### 可立即对齐（有 OCCT 源码，rcad 有对应函数）

| 优先级 | OCCT 函数 | rcad 文件 | 问题 |
|--------|-----------|-----------|------|
| **P0** | `FillImagesFaces` (BOPAlgo_Builder_1.cxx:128-170) | filler.rs:181 | 面分裂+分类，直接影响 F=7→final BRep |
| **P0** | `PutPaveOnCurve` (BOPAlgo_PaveFiller_6.cxx:2940-3010) | paves.rs:244 | 顶点在曲线上投影+登记，影响 V=15 |
| **P1** | `MakeBlocks` 主循环 (L725-1107) | make_blocks.rs:15 | 多重循环+变量，影响全部 |
| **P1** | `IsValidBlockForFaces` | make_blocks.rs:438 | pcurve 求值方式不同 |

### rcad 独有 / OCCT 没有（需删除）

| rcad 概念 | 文件 | 说明 |
|-----------|------|------|
| SubFace / split_planar_face | builder/ | OCCT 无对应。AGENTS.md 明确说不要对齐 |
| merge_close_vertices | brep_repair/ | OCCT 无此步骤，顶点合并在 DS 内完成 |
| SubShape classification | classify.rs | 与 OCCT 不同 |

### OCCT 有 / rcad 没有（需补齐）

| OCCT 函数 | 用途 | 影响 |
|-----------|------|------|
| `DoSplitSEAMOnFace` | 在处理周期面(polar)时分割 seam 边 | 球面 seam 处理可能缺失 |
| `ComputeTolReached3d` | 3D 公差传播 | 精度差异 |
| `GeomAPI_ProjectPointOnCurve` 缓存 | 投影器缓存提升性能 | 不影响正确性 |
| `ProjLib_Sphere::Project(gp_Lin)` | 线在球面的投影 | 影响线-球相交场景 |

## 对齐优先级建议

1. **P0**: `FillImagesFaces` — 面分裂，直接影响拓扑对齐
2. **P0**: `PutPaveOnCurve` — 顶点管理，直接影响 V=15
3. **P1**: `MakeBlocks` 主循环 — 综合影响
4. **P1**: `IsValidBlockForFaces` — pcurve 求值
