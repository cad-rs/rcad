---
name: pavefiller_ve_alignment_handover
description: Handover for continuing after_PerformVE alignment in a new session
metadata:
  type: project
  project: boolean-alignment
---

# PaveFiller VE/EE alignment handover

## Session result

`after_PerformVV` 全面对齐通过。所有 25 个分阶段测试的 VV 阶段已与 OCCT 参考一致。

## 当前状态

| 阶段 | 通过率 |
|------|--------|
| after_Init | 25/25 ✅ |
| after_Prepare | 25/25 ✅ |
| after_PerformVV | 25/25 ✅ |
| after_PerformVE | 15 fail at nPB/VE, 10 pass |
| after_PerformEE | 10 fail at nV/nPB |

## after_PerformVE 失败模式

所有失败都是 rcad 生成的 PaveBlock 数（nPB）或 VE 干涉数（VE）低于 OCCT 参考。典型：

- `nPB=1/3` — 只有 1 个 PaveBlock，OCCT 有 3 个
- `VE=0/1` — 0 个 VE 干涉，OCCT 有 1 个
- `VE=0/2` — 0 个 VE 干涉，OCCT 有 2 个

失败全部在 `check_pf_stages` 的断言中触发（桩文件 rcad/libs/rcad-algo/tests/pavefiller_stage_tests.rs）。

## 关键文件

- **pave_filler.rs**: `perform_ve`, `fill_shrunk_data`, `intersect_ve`, `init_pave_blocks`
- **iterator.rs**: `BOPDS_Iterator::intersect` — VE pair 选择（已用 BVH 对齐）
- **algo_tools.rs**: `compute_vv`、`is_on_pave_1` 等
- **box_tree.rs**: `BoxTree` + `PairSelector`（已用 LBVH 对齐）
- **bvh_tree.rs**: `BvhTreeBase`（SOA 格式，已对齐）

## 已对齐的部分

compute_vv（sum+Confusion）、BoxTree（LBVH+Morton+PairSelector）、BvhTreeBase（SOA node storage）、BOPDS_Iterator::intersect（BVH-based pair selection）、is_on_pave_1、fill_map、make_blocks、init_pave_blocks/change_pave_blocks 分离、BOPAlgo_BPC、IsBasedOnPlane 等。

## 已知问题

- `rcad-modeling/Cargo.toml` 移除了 `rcad-algorithms` dev-dep（源文件未使用），`tests/primitives.rs` 改用 `rcad_kernel::surface_area`/`volume`
- `BRep::edge_mut/face_mut/shell_mut` 改用 `Arc::make_mut` 而非 `get_mut`（支持 COW 语义）

## 一句话提示词，在新 session 中执行

> 继续对齐 after_PerformVE：检查 pave_filler.rs 的 perform_ve 和 fill_shrunk_data，对比 OCCT BOPAlgo_PaveFiller_2.cxx L171-238，找出 PaveBlock 和 VE 干涉数量差异的原因。
