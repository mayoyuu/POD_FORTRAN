# STT 传播器允许 1 阶选择 — 设计规格

**日期**: 2026-06-16
**目标**: 放宽 STT 传播器的阶数下限，从 `>= 2` 改为 `>= 1`，使 order=1（纯 STM 线性传播）可用。

---

## 1. 动机

当前 `uq_stt_propagator` 的阶数验证要求 `order >= 2`，无法选择 1 阶（STM-only 模式）。但 `test_stt_crtbp` 测试需要 order=1 来验证 STM 退化（与 DA 线性传播、MC 参考解对比）。底层代码（`compute_stt_moments`、`stt_eval_deviates`、`compute_stt_aug_derivatives`）已正确处理 order=1，只需放宽验证边界。

## 2. 变更范围

**唯一修改文件**: `src/lib/uncertainty/propagation/pod_uq_prop_stt_module.f90`

| 位置 | 行 | 当前 | 改为 |
|------|----|------|------|
| 模块注释 | 5 | `STT 阶数 (2~6)` | `STT 阶数 (1~6)` |
| `set_stt_order` 验证 | 87 | `order >= 2` | `order >= 1` |
| `set_stt_order` 错误消息 | 90 | `order must be 2..6` | `order must be 1..6` |
| `stt_propagate` 验证 | 216 | `stt_order < 2` | `stt_order < 1` |
| `stt_propagate_deviates` 验证 | 569 | `stt_order < 2` | `stt_order < 1` |

## 3. 无需变更的部分

以下代码已正确处理 order=1，无需修改：

- `compute_stt_moments`: STM 线性项无条件计算，高阶修正通过 `if (order >= N)` 守卫
- `stt_eval_deviates`: order=1 时仅执行 STM 矩阵乘法，高阶循环 `do p=2,order` 为空
- `compute_stt_aug_derivatives`: 同上，高阶 RHS 循环 `do p=2,order` 为空
- `stt_store%init`: 当 order=1 时 `max_sz = max(0, 1) = 1`，正确分配
- `stt_flat_offset`: p=1 时 `do pp=1,0` 为空，off 计算正确

## 4. 行为说明

order=1 时：
- 增广系统维度: 6 (标称) + 36 (STM) = 42
- 仅积分标称轨道 + STM，不计算高阶张量
- 矩计算: 均值 = x*(标称终点)，协方差 = STM · P₀ · STMᵀ
- 与 DA 的 order=1 结果应完全一致（均为纯线性传播）
