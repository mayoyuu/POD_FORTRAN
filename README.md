# POD System: Precision Orbit Determination & Propagation Framework

![Fortran](https://img.shields.io/badge/Fortran-90%2F2008-734f96?logo=fortran)
![Platform](https://img.shields.io/badge/platform-Linux-blue)
![License](https://img.shields.io/badge/license-MIT-green)

## 📌 Overview (项目简介)
POD System 是一个基于 Fortran 开发的高精度轨道传播与不确定性量化（UQ）框架。本项目不仅支持传统的高精度数值积分（如 RKF45/RKF78），还深度集成了**微分代数（Differential Algebra, DA）**技术，能够高效处理复杂引力场下的高阶状态转移矩阵（STM）计算以及非线性误差传播分析。

适用于深空探测、地月空间（Cislunar）轨道动力学演化及蒙特卡洛（MC）协方差分析。

## ✨ Key Features (核心特性)
* **多模式轨道传播：** 支持真实物理量积分与无量纲化动力学方程积分。
* **微分代数（DA）引擎集成：** 借助 DACE 库，支持任意阶泰勒多项式展开，一键获取高阶非线性动力学特征。
* **高精度引力模型：** 内置日地月多体引力及高阶地球/月球非球形引力场（支持 SPICE 星历调用）。
* **统一的 UQ 接口：** 极简的 API 设计，只需切换 `METHOD_MC` 或 `METHOD_DA` 即可无缝对比蒙特卡洛与微分代数的误差传播结果。

## 🛠️ Dependencies (依赖库)
本项目依赖以下外部库：
* **SPICE Toolkit (`libspicelib.a`)**: 用于获取高精度行星星历与时间系统转换。
* **DACE (Differential Algebra Core Engine)**: 提供底层微分代数数据结构与 C/C++ 接口 (`libdace_s.a`)。
* **BLAS/LAPACK**: 用于高效的线性代数运算。
* **external/dace_build**: 包含 DACE 的 Fortran 包装器`libdace_wrapper.a`，简化 Fortran 与 DACE 的交互。
* **Fortran Package Manager (fpm)** 或 **CMake**: 用于项目构建。

## 🚀 Installation (安装指南)

1. **克隆仓库及子模块（获取 DACE 依赖）：**
   ```bash
   git clone --recursive [https://github.com/mayoyuu/POD_FORTRAN.git](https://github.com/mayoyuu/POD_FORTRAN.git)
   cd POD_System

2. **编译外部依赖：**

   Bash

   ```
   git clone https://github.com/dacelib/dace.git
   bash setup_env.sh
   ```

1. **配置物理环境：** 修改 `config/pod_config.txt` 中的引力模型与积分器容差设置。

## 📖 Quick Start (快速开始)

项目提供了顶层封装接口，可以极其方便地进行轨道传播与误差分析。

### 1. 标称轨道传播 (RKF78)

Fortran

```
use pod_engine_module, only: pod_engine_init
use pod_orbit_propagation, only: orbit_state, propagation_result, propagate_orbit
use pod_da_integrator_module, only: METHOD_RKF78

! 初始化物理环境
call pod_engine_init('config/pod_config.txt')

! 装载初始状态
initial_state%state = [100000.0_DP, 50000.0_DP, 20000.0_DP, 1.5_DP, 2.5_DP, 0.5_DP]

! 传播 1 天 (86400 秒)
call propagate_orbit(initial_state, 86400.0_DP, METHOD_RKF78, result)
```

### 2. 微分代数 (DA) 误差传播分析

仅需设定 DA 阶数并调用 `run_uq_propagation`：

Fortran

```
use pod_uq_propagation, only: run_uq_propagation, METHOD_DA

! 设置展开阶数等初始条件...
call run_uq_propagation( &
    nominal_state = nominal_orbit, &
    initial_cov   = initial_covariance, &
    epoch0        = epoch_start, &
    t_start       = 0.0_DP, &
    t_end         = 86400.0_DP, &
    method_switch = METHOD_DA, &      ! 核心开关：切换至 DA 算法
    n_particles   = 10000, &
    initial_state_out = initial_dist, &
    final_state_out   = final_dist_da &
)
```

## 📂 Repository Structure (项目结构)

*为了保持整洁，这里只展示核心源码结构。*

Plaintext

```
.
├── app/                  # 包含主程序 demo (如 pod_demo.f90)
├── config/               # 物理引擎配置文件
├── src/                  
│   ├── functions/        # 顶层功能封装 (传播、定轨API)
│   └── lib/              # 底层核心库 (力模型、积分器、DA核心、参考系等)
├── test/                 # 单元测试与集成测试用例
└── setup_env.sh          # 环境配置脚本
```

## 📝 License & Citation

This project is licensed under the MIT License - see the LICENSE file for details. *(如果你有相关的学术论文正在撰写或准备发表，可以在这里提示：“If you use this code in your research, please cite: [TBD]”*