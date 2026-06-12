# STT 精度发散修复 — 设计规格

**日期**: 2026-06-12
**目标**: 修复 STT 计算中的两个缺陷：组合数学张量收缩遗漏（精度发散）和 DIAGNOSTIC 慢路径（性能断崖）

---

## 1. 问题诊断

### 1.1 症状

teststt 测试结果 (order=2, vars=6, deviates mode):

```
FAIL: mean mismatch, max diff = 2.1452E-04
STT mean: 8.498253E-01 -3.818007E-05 ... 4.790272E-01 ...
DA  mean: 8.498769E-01 -5.724562E-05 ... 4.788127E-01 ...

STT cov diag: ... 3.036965E-07 ...   ← 分量 5 与 DA (5.192975E-07) 偏差 ~42%
DA  cov diag: ... 5.192975E-07 ...
```

均值最大差异 ~2.1e-4，STT 协方差对角元分量 5 严重偏低。

### 1.2 根因

**Bug 1 — 组合数学张量收缩遗漏 (精度发散)**：
`compute_stt_rhs_order` 中非均匀项累加时只除以等大小块排列因子 `eq_fac`，未乘以赋值分布的多重性系数。导致部分 contractions 被低计。

**Bug 2 — DIAGNOSTIC 慢路径 (性能断崖)**：
`compute_stt_aug_derivatives` 使用三重嵌套循环的逐分量旧版 `compute_stt_rhs`，而非批量优化的 `compute_stt_rhs_order`。每次调用需要 pack/unpack stt_store 对象，且每个 (i_comp, cidx) 对都重新计算块指配。

---

## 2. 修复方案

### 2.1 Fix 1: 添加多重量到 gen_block_assignments_v2

**文件**: `src/lib/math/pod_stt_tensor_module.f90`

**变更 1**: `gen_block_assignments_v2` 新增输出参数 `mults_out(:)`

```fortran
subroutine gen_block_assignments_v2(tup, p, part_sizes, k, blocks_out, n_assign, mults_out)
    ...
    integer, intent(out) :: mults_out(:)  ! 新增
```

在 `enum_dist_v2` 中找到有效赋值后紧跟在 `bpos = bpos + 1` 之后计算多重性：

```fortran
bpos = bpos + 1
mult = 1
do i = 1, uniq_n
    fac_i = factorial(cnts(i))
    do blk = 1, kk
        fac_i = fac_i / factorial(gba_dist(i, blk))
    end do
    mult = mult * fac_i
end do
mults_out(bpos) = mult
```

**变更 2**: `compute_stt_rhs_order` 分配并使用 `blk_mults`

- 新增 `blk_mults(pl%n_parts, MAX_ASSIGN_PER_PART)` 与 `blk_assigns` 并列分配
- 传入 `gen_block_assignments_v2` 作为额外实参
- 非均匀项累加改为: `dphi = dphi + term * real(blk_mults(ip, ia), DP) / real(eq_fac, DP)`
- 在结束前正确释放 `blk_mults`

### 2.2 Fix 2: 移除 DIAGNOSTIC 路径

**文件**: `src/lib/uncertainty/propagation/pod_uq_prop_stt_module.f90`

**变更**:

1. 删除局部变量 `stt_tmp`, `pp`, `np2`, `cidx2`, `tup2`
2. 删除整个 DIAGNOSTIC 块 (lines ~169-195):
   - `stt_tmp%init(order)` + pack 循环
   - 逐分量 `compute_stt_rhs` 调用
   - `stt_tmp%destroy()`
3. 替换为优化循环:

```fortran
! ---- 5. 高阶 STT RHS (order 2..max_order) ----
do p = 2, order
    call compute_stt_rhs_order(p, y(7:), f_tensors, order, dydt, pos)
    pos = pos + 6 * stt_sizes(p)
end do
```

---

## 3. 预期结果

修复后的 teststt:

- **均值匹配**: max diff < 1e-6 (当前 2.1e-4)
- **协方差 STT-vs-DA**: rel err 保持或改善 (当前 9.2e-4)
- **协方差 DA-vs-MC**: 不受影响 (当前 2.5e-3)
- **性能**: 积分器每步开销显著降低 (消除 O(n_parts × n_stt × 6^k) pack 成本)

---

## 4. 文件变更清单

| 文件 | 变更类型 |
|------|---------|
| `src/lib/math/pod_stt_tensor_module.f90` | 修改 — `gen_block_assignments_v2` 添加 mults_out, `compute_stt_rhs_order` 使用多重性 |
| `src/lib/uncertainty/propagation/pod_uq_prop_stt_module.f90` | 修改 — `compute_stt_aug_derivatives` 替换 DIAGNOSTIC 路径为优化路径 |

---

## 5. 不纳入范围

- 不改动 `compute_stt_rhs` (旧版逐分量 RHS) — 仅优化版 `compute_stt_rhs_order` 受影响
- 不修改 `stt_eval_deviates` 或矩计算公式
- 不涉及公共接口变化
