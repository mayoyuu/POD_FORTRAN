# STT 误差传播 — 设计规格

**日期**: 2026-06-05
**目标**: 在 POD_FORTRAN 中实现基于 State Transition Tensor (STT) 的非线性误差传播方法，先做 CRTBP 下的实现，再分析高精度力模型的挑战。

---

## 1. 数学背景

STT 方法将轨道偏差 δx(t) 展开为初始偏差 δx° 的 m 阶泰勒级数：

δx_i(t) = Σ_{p=1}^{m} (1/p!) Φ_{i, k₁...kₚ} δx°_{k₁} ... δx°_{kₚ}

其中 Φ_{i, k₁...kₚ} 是状态转移张量（STT），满足微分方程（1~4 阶见公式 17-20，高阶由 Faà di Bruno 公式递推）：

Φ̇_{i, k₁...kₚ} = f*_{i, α} Φ_{α, k₁...kₚ} + [涉及低阶 STT 与力场导数 f* 的组合项]

初始条件：Φ_{i,a}(t°) = δ_{i,a}（STM 为单位阵），高阶张量初值为零。

STT 求解后，利用非线性矩公式计算传播后的均值与协方差。

---

## 2. 架构设计

### 2.1 新增模块

```
src/lib/math/pod_stt_tensor_module.f90          — 对称张量存储、索引、Faà di Bruno 组合
src/lib/forcemodel/pod_crtbp_derivatives_module.f90 — CRTBP 力场任意阶解析导数
src/lib/uncertainty/propagation/pod_uq_prop_stt_module.f90 — STT 传播器（继承 uq_propagator_base）
```

**依赖关系**：

```
pod_uq_prop_stt_module
 ├── pod_stt_tensor_module          (张量组合/索引)
 ├── pod_crtbp_derivatives_module   (力场解析导数)
 ├── pod_uq_base_module             (基类)
 ├── pod_integrator_module          (复制并改造的 RKF78，支持可变维度)
 └── pod_global / pod_config
```

STT 传播器**不依赖** DACE 引擎（pod_dace_classes），完全基于实数运算。

### 2.2 与现有架构的集成

STT 传播器通过 UQ 传播体系接入，与现有 `uq_da_propagator`、`uq_propagation` 等模块同级。用户调用方式：

```fortran
type(uq_stt_propagator) :: propagator
call propagator%set_stt_order(4)           ! 选择 STT 阶数
call propagator%set_integrator(METHOD_RKF78)
call propagator%propagate(t_start, t_end, input_state, output_state)
```

---

## 3. 子模块设计

### 3.1 pod_stt_tensor_module — 对称张量基础库

**职责**：提供 STT 计算所需的对称张量数据结构和组合算法。

**核心功能**：

1. **对称索引压缩映射**
   - 利用 FT 中 n 选 p 的可重复组合，将对称指标组 (k₁≤k₂≤...≤kₚ) 映射为词典序整数索引
   - 预计算双向查找表（order ≤ max_stt_order, 维度固定 6）
   - 提供 `sym_index_to_tuple(p, idx)` 和 `tuple_to_sym_index(tuple)` 转换

2. **STT 存储**
   - 在传播器生命周期内分配 `sttp(6, n_p)` 数组，其中 n_p = C(6+p-1, p)
   - 支持根据压缩索引读写对称张量元素

3. **Faà di Bruno 组合引擎**
   - 预计算 1~m 阶的 integer partition 表
   - 对每个阶数 p，穷举 partition λ ⊢ p，生成对应的指标置换组合
   - 对每个 partition λ = (λ₁,...,λₖ)，按 STT 方程组装 Σ f*_{i, α₁...αₖ} · Φ_{α₁, part1} · ... · Φ_{αₖ, partk}
   - 结果写入 Φ̇_{i, k₁...kₚ} 的对应压缩索引位置

4. **矩计算**
   - 根据 STT 展开公式计算非线性均值修正和协方差修正（2~6 阶截断）
   - 输入：标称终点状态 x*(6)，STT(1~m)，初始协方差 P(6,6)
   - 输出：修正后的均值 μ(6)，协方差 Σ(6,6)

**接口**：

```fortran
public :: init_stt_indexing, stt_order_p_size
public :: sym_index_to_tuple, tuple_to_sym_index
public :: compute_stt_rhs
public :: compute_stt_moments
```

### 3.2 pod_crtbp_derivatives_module — CRTBP 解析导数

**职责**：计算 CRTBP 力场沿标称轨道的 1~m 阶偏导数 f*_{i, k₁...kₚ}。

**核心递推**：

- 定义 U_n = 1/r^(2n+1)，递推公式 ∂/∂q U_n = -(2n+1)·(q-q₀)·U_{n+1}
- 对 r₁（初始天体）和 r₂（目标天体）分别递推
- 离心项 (x²+y²)/2 的导数平凡（一阶和二阶非零）
- 速度分量的偏导数：仅 f₄ 含 +2vy，f₅ 含 -2vx，二阶及以上速度导数为零
- β 混合偏导数：由于 CRTBP 中位置和速度的可分离性，所有交叉偏导数在二阶以上为零

**接口**：

```fortran
public :: crtbp_force_derivatives
! 输入: x(6), mu, max_order
! 输出: f_tensors — 压缩对称格式的力场导数张量（1~max_order 阶）
```

**关键实现点**：利用 U_n 递推一次遍历即可产生所有阶数的位置偏导数，避免逐阶重复计算。力场导数同样使用压缩对称格式存储，与张量模块保持一致。

### 3.3 pod_uq_prop_stt_module — STT 传播器

**职责**：将标称轨道 + STT 联合积分，继承 `uq_propagator_base`。

**自定义数据类型**：

```fortran
type, extends(uq_propagator_base) :: uq_stt_propagator
    integer :: stt_order = 2           ! STT 阶数 (2~6+)
    integer :: total_stt_size = 0      ! 增广 ODE 总维度
    type(stt_index_tables) :: tables   ! 预计算的索引表
contains
    procedure :: propagate => stt_propagate
    procedure :: set_stt_order
    procedure :: get_method_name
end type
```

**传播流程**：

```
propagate():
  1. 从 input_state 获取初始均值 μ₀(6) 和协方差 P₀(6,6)
  2. 预计算组合/索引表，分配增广状态 (6 + total_stt_size)
  3. 初始条件：x(1:6) = μ₀，STM = I₆，高阶 STT = 0
  4. 调用改造版 RKF78 自适应积分器
  5. 从终点增广状态解包标称轨道 x* 和所有 STT
  6. 调用 compute_stt_moments 计算终点均值与协方差
  7. 填充 output_state
```

**联合积分 RHS**：

每一步计算 `compute_stt_derivatives(aug_state, time)`：
1. 从 aug_state 解包 x(6) 和 STT
2. 调用 `crtbp_force_derivatives(x, mu, stt_order)` 获取 f* 张量
3. 对每个阶数 p=1..stt_order 调用 `compute_stt_rhs` 填充 STT 导数
4. 打包返回 augment_derivatives

**积分器改造**：

从 `pod_integrator_module` 复制 RKF78 核心逻辑到 STT 模块内部，将固定维度 6 替换为可配置维度。保留 WRMS 步长控制、安全因子等机制不变。

**强制初始化检查**：propagate 入口处验证 `total_stt_size > 0`，否则报错提示先调用 `set_stt_order`。

### 3.4 配置文件扩展

在 `pod_config_module` 中新增：

```fortran
integer :: default_stt_order = 2       ! STT 默认阶数
real(DP) :: stt_integration_abs_tol = 1.0d-12  ! STT 积分绝对容差
```

---

## 4. 数据流

```
Input: 初始均值 μ₀(6), 协方差 P₀(6,6)

1. 初始化增广状态
   aug = [μ₀(1:6) | I₆ | 0_{126} | 0_{336} | ...]  (取决于阶数)

2. 联合积分循环 (自适应 RKF78)
   每步:
     x = aug(1:6)
     f_tensors = crtbp_force_derivatives(x, mu, order)
     STT_rhs = compute_stt_rhs(f_tensors, STT, order)
     aug_derivs = [crtbp_derivatives_real(x) | I·f* | STT_rhs(2) | STT_rhs(3) | ... ]

3. 终点: x*(6), STT(1..order)

4. 矩计算
     μ = x* + 非线性修正(Φ, P₀)
     Σ = Φ·P₀·Φᵀ + 高阶修正(Φ, P₀)

Output: μ(6), Σ(6,6)
```

---

## 5. 验证策略

### 5.1 与 DA 交叉验证

CRTBP 下，DA 和 STT 应产生一致结果。验证方式：
- 相同初始条件（均值、协方差）
- DA 阶数 = STT 阶数
- 比较终点均值、协方差的 Frobenius 范数差异
- 预期差异在 DA/STT 的截断误差量级（10⁻¹⁰ ~ 10⁻¹² 量级）

### 5.2 与 MC 基准对比

- 用大量（10⁵）MC 粒子的统计矩作为参考
- 比较 DA、STT、MC 三者的均值/协方差

### 5.3 阶数收敛性

- 固定传播时间，m=1..6 阶
- 验证高阶修正项递减，矩值趋于稳定

### 5.4 STM 退化测试

- STT order=1 退化为标准 STM 传播
- 与理论线性传播结果对比（均值不变，协方差 = STM·P₀·STMᵀ）

---

## 6. 高精度力模型的挑战

### 6.1 力项分类与导数可行性

| 力项 | 数学模型 | 解析递推 |
|------|---------|:---:|
| N 体引力 (×11) | GM·r/|r|³ | 可（同 CRTBP 的 1/r 递推，乘以 11 个天体） |
| 地球非球形引力场 (70×70) | Σ P_{nm}(sin φ)·(R/r)^n·三角项 | 极其复杂 |
| 月球非球形引力场 | 同上 | 极其复杂 |
| 太阳光压 (SRP) | 柱/锥阴影模型 | 不可（不光滑） |
| 相对论修正 | v²/c² 多项式 | 可（多项式） |
| SPICE 坐标旋转 | 时间函数 | 平凡（∂R/∂x=0） |

### 6.2 核心问题

1. **球谐展开的组合爆炸**：5000+ 项 × 6 阶混合偏导数 × 笛卡尔-球坐标链式法则。每项连带勒让德多项式 P_{nm} 虽可递推，但笛卡尔偏导数需要展开到底的所有交叉项

2. **阴影函数不光滑**：SRP 的阴影边界处导数不连续甚至不存在，STT 需要逐点精确导数值，无法像 DA 那样通过多项式"平滑跨越"

3. **每步开销**：若退而用 DACE 辅助求 f*（每步创建/销毁 DA + 求值 + 提取系数），每步开销 = 13 阶段 × DA 力模型（N 体+球谐+SRP+相对论），约为纯 DA 传播的若干倍，失去 STT 相对于 DA 的效率意义

4. **稀疏性丧失**：CRTBP 中大量导数恒为零（速度分量、混合偏导数），高精度力模型几乎所有分量都耦合（球谐项对所有位置坐标产生非零导数），Faà di Bruno 组合项数全量级激活

### 6.3 可能的缓解路径

- **只对部分力项用 STT**：N 体 + 相对论用解析导数，球谐部分忽略高阶耦合或用低阶展开近似
- **DACE 辅助 + 大步长**：接受每步的 DACE 开销，但利用 STT 的线性齐次结构使用比 DA 积分器更大的步长
- **预处理球谐导数**：用符号推导为常用阶数（如 ≤ 70 阶）预生成球谐导数的查表，避免运行时解析递推

### 6.4 结论

CRTBP 下 STT 完全可行且高效（解析导数 + 纯实数积分，无 DACE 开销）。高精度力模型下，纯解析路径不现实；DACE 辅助路径引入了"不如直接 DA 传播"的困境。建议 CRTBP 部分完成后，根据性能对比结果和实际需求决定是否推进高精度 STT。

---

## 7. 文件变更清单

### 新增文件

| 文件 | 描述 |
|------|------|
| `src/lib/math/pod_stt_tensor_module.f90` | 对称张量存储、索引、Faà di Bruno、矩计算 |
| `src/lib/forcemodel/pod_crtbp_derivatives_module.f90` | CRTBP 解析力场导数递推 |
| `src/lib/uncertainty/propagation/pod_uq_prop_stt_module.f90` | STT 传播器 |

### 修改文件

| 文件 | 变更 |
|------|------|
| `src/lib/system/pod_config_module.f90` | 新增 `default_stt_order`、`stt_integration_abs_tol` |

### 测试文件

| 文件 | 描述 |
|------|------|
| `test/test_stt_crtbp.f90` | STT vs DA vs MC 对比测试、阶数收敛性、STM 退化测试 |

---

## 8. 不纳入范围

- 高阶 STT (m > 6) 的优化（Faà di Bruno 项数指数增长，6 阶以上不保证实时性）
- 高精度力模型的 STT 实现（仅限分析和讨论）
- STT 的 GPU 加速
- 定轨滤波器（EMDAC/UKF）中的 STT 应用
