# POD System: Precision Orbit Determination & Propagation Framework

## 📌 Overview

POD System 是一个基于 Fortran 开发的高精度轨道传播、定轨（Orbit Determination）与不确定性量化（UQ）框架。

本项目不仅支持传统的高精度数值积分与**无迹卡尔曼滤波 (UKF)**，还深度集成了**微分代数 (Differential Algebra, DA)** 技术。通过集成 **EMDAC-n ** 滤波器，能够高效处理高阶状态转移矩阵（STM）计算、非线性误差传播分析以及复杂动力学环境下的高精度定轨任务。

## ⚠️ 编译器要求与严重内存泄漏风险

> **必须使用 GNU Fortran（gfortran）13.3.0 或更高版本。低于 13.3.0 的版本不受支持，请勿用于本项目的 DA / ADS / DAMC 长时或批量传播。**

DA 运算符会创建持有 DACE 原生句柄的临时对象，其释放依赖 Fortran 自动析构（finalization）。已在 **gfortran 11.4.0** 上复现函数返回临时对象未执行析构：100 次表达式后仍遗留 100 个模拟句柄；同一检查在 **13.3.0** 的 `-O0` / `-O3` 下均正常释放。此类问题在真实 DA 计算中可能造成持续、严重的内存泄漏，即使数值结果正常也不能认为内存安全。13.3.0 是本项目当前支持下限，并非声称所有更低版本均已逐一测试。

多版本可以共存，无需修改系统默认编译器。先确认所选编译器的实际版本，再让 fpm 明确使用它（`gfortran-13` 的名称本身不保证版本已达到 13.3.0）：

```bash
gfortran-13 --version                  # 确认版本 >= 13.3.0
export FPM_FC="$(command -v gfortran-13)"  # 或填写新编译器的完整路径
export FC="$FPM_FC"                   # 同时供 CMake 等工具使用
"$FPM_FC" --version
```

在项目根目录加载原有依赖环境后，使用选定编译器构建：

```bash
source ./setup_env.sh
fpm build --compiler "$FPM_FC" --profile release --flag "-ffree-line-length-none"
```

**升级后必须重新编译全部项目 Fortran 源码（包括 DA 模块）**，不要混用旧编译器的 `.mod`、`.o`、项目静态库或可执行文件；旧构建目录可先改名备份。仅重编 C++ `dace_wrapper`、更换运行库或设置环境变量不会修复旧可执行文件中的析构调用。新终端及后台批量任务也需选择相同编译器，并先进行短时真实 DA 资源释放检查，再启动长时传播。满足版本要求不等于已经排除所有内存泄漏。

## ✨ Key Features

- **多模式轨道传播** — 支持真实物理量积分与无量纲化动力学方程积分
- **先进定轨算法** — 内置 **UKF** 以及基于微分代数的 **EMDAC-n** 高阶非线性滤波器
- **微分代数 (DA) 引擎** — 借助 DACE 库，支持任意阶泰勒多项式展开，一键获取高阶非线性动力学特征
- **高精度引力模型** — 支持日地月多体引力及高阶非球形引力场（支持 SPICE 星历调用）
- **统一的 UQ 接口** — 极简 API 设计，支持 `METHOD_MC` 或 `METHOD_DA` 两种不确定性传播方法切换

## 🛠️ Dependencies

在 Linux 服务器（无 `sudo` 权限）环境下，请按照以下步骤配置依赖：

### 1. DACE (微分代数库)

```bash
git clone https://github.com/dacelib/dace.git dace
cd dace && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/local_dace ..
make && make install
```

### 2. SPICE Toolkit (星历工具)

```bash
mkdir -p ~/spice_fortran && cd ~/spice_fortran
wget https://naif.jpl.nasa.gov/pub/naif/toolkit//FORTRAN/PC_Linux_GFORTRAN_64bit/packages/toolkit.tar.Z
tar -xzvf toolkit.tar.Z
```

### 3. OpenBLAS (线性代数库)

```bash
git clone https://github.com/OpenMathLib/OpenBLAS.git
cd OpenBLAS
make -j4
make PREFIX=$HOME/local_openblas install
```

### 4. 项目本地静态库 (`lib/`)

`fpm.toml` 通过 `link = ["spicelib", "m", "dace_wrapper", "dace_s", ...]` 链接静态库，
而 `setup_env.sh` 只把 `$PROJECT_ROOT/lib` 加入 `LIBRARY_PATH`，因此需要在该目录准备三份库文件：

```bash
# 1) 编译 DACE 的 C++ 包装库（首次或 wrapper 改动后执行）
cd "$PROJECT_ROOT/external/dace_build"
mkdir -p build && cd build && cmake .. && make dace_wrapper && cd "$PROJECT_ROOT"

# 2) 建立链接（链接真实文件即可，便于后续更新）
mkdir -p "$PROJECT_ROOT/lib"
ln -sfn "$PROJECT_ROOT/external/dace_build/build/libdace_wrapper.a" "$PROJECT_ROOT/lib/libdace_wrapper.a"
ln -sfn "$HOME/local_dace/lib/libdace_s.a"                          "$PROJECT_ROOT/lib/libdace_s.a"
ln -sfn "$HOME/spice_fortran/toolkit/lib/spicelib.a"                "$PROJECT_ROOT/lib/libspicelib.a"
```

## 🚀 Installation

### 环境配置

每次编译前，在项目根目录执行以下命令激活环境（请先根据实际安装路径修改 `setup_env.sh` 中的路径）：

```bash
source ./setup_env.sh
```

### 编译构建

```bash
git clone --recursive https://github.com/mayoyuu/POD_FORTRAN.git
cd POD_FORTRAN
source ./setup_env.sh          # 激活环境
mkdir build && cd build
cmake -DCMAKE_Fortran_COMPILER="$FPM_FC" ..
make
```

## 📖 Quick Start

### 1. EMDAC-N 定轨

```fortran
use pod_emdac_runner_module, only: run_emdac_orbit_determination

call run_emdac_orbit_determination( &
    obs_file          = 'obs_data.txt', &       ! 观测数据文件
    site_json_file    = 'station_config.json', & ! 测站信息
    gmm_in_switch     = .true., &               ! 是否启用 GMM 初始化
    initial_json_file = 'initial_guess.json', & ! 初始名义轨道与协方差
    output_file_name  = 'od_result.out', &      ! 结果输出路径
    opt_particles     = 10000, &                ! 粒子数
    max_da_order      = 5, &                    ! DA 阶数
    opt_em_max_iter   = 20, &                   ! EM 算法最大迭代次数
    opt_em_tol        = 1.0d-6, &               ! 收敛容差
    n_components      = 3)                      ! 高斯混合分量数
```

### 2. UT 定轨

```fortran
use pod_ut_runner_module, only: run_ut_orbit_determination

call run_ut_orbit_determination( &
    obs_file          = 'obs_data.txt', &
    site_json_file    = 'station_config.json', &
    initial_json_file = 'initial_guess.json', &
    output_file_name  = 'od_result.out')
```

### 3. DA 轨道传播

```fortran
use pod_da_orbit_propagation, only: da_orbit_state, propagate_da_orbit, &
                                    extract_stm_from_result

type(da_orbit_state) :: state
type(da_propagation_result) :: result
real(DP) :: stm(6, 6)

state%nominal_state = [x, y, z, vx, vy, vz]
state%epoch         = tdb_epoch
state%da_order      = 3

call propagate_da_orbit(state, 86400.0_DP * 4, 2, result)
call extract_stm_from_result(result, stm)
```

### 4. UQ 不确定性传播

```fortran
use pod_uq_propagation, only: run_uq_propagation, METHOD_DA, METHOD_MC

! DA 方法（高效，一次传播获得完整统计信息）
call run_uq_propagation(nominal_state = orbit, initial_cov = cov, &
    method_switch = METHOD_DA, final_state_out = result_da)

! MC 方法（参考基准，需大量粒子采样）
call run_uq_propagation(nominal_state = orbit, initial_cov = cov, &
    method_switch = METHOD_MC, n_particles = 10000, &
    final_state_out = result_mc)
```

## 📘 Example Programs

`example/` 目录下提供了四个可直接运行的示例程序：

| 示例 | 文件 | 说明 |
|------|------|------|
| EMDAC-N 定轨 | `example_emdac_orbit_determination.f90` | 基于 DA 的高阶非线性定轨 |
| UT 定轨 | `example_ut_orbit_determination.f90` | 基于无迹变换的定轨 |
| DA 轨道传播 | `example_da_orbit_propagation.f90` | DA 轨道传播与 STM 提取 |
| UQ 传播 | `example_uq_propagation.f90` | DA vs MC 不确定性传播对比 |

```bash
# 使用 fpm 运行
fpm run --example example_emdac_orbit_determination
fpm run --example example_ut_orbit_determination
fpm run --example example_da_orbit_propagation
fpm run --example example_uq_propagation

# 或直接运行编译后的可执行文件
./build/example/example_emdac_orbit_determination
```

> **提示：** 运行前请确保已正确配置 SPICE 内核路径 (`src/lib/system/pod_spice.f90`) 且 DACE 库已正确链接。示例中的文件路径可能需要根据实际数据位置进行调整。

## 📂 Repository Structure

```
.
├── app/                    # 主程序入口 (pod_demo.f90)
├── example/                # 示例程序（可直接运行）
├── src/
│   ├── functions/
│   │   ├── orbitimprove/   # 定轨 API (EMDAC-N, UT 滤波器)
│   │   └── orbitprop/      # 传播 API (实数/DA 传播, UQ)
│   └── lib/
│       ├── api/            # 引擎初始化接口
│       ├── data/           # 数据格式处理
│       ├── forcemodel/     # 力模型（实数/DA 版）
│       ├── frame/          # 坐标框架
│       ├── integrator/     # 积分器（实数/DA 版）
│       ├── math/           # 数学工具
│       ├── measurement/    # 观测模型
│       ├── statistics/     # 统计与 GMM
│       ├── system/         # 系统模块 (SPICE, DACE, 全局参数)
│       ├── timesystem/     # 时间系统
│       ├── uncertainty/    # UQ 数据结构、传播、变换
│       └── utils/          # 通用工具
├── test/                   # 单元测试与集成测试
├── config/                 # 配置文件 (SPICE 内核路径等)
├── input/                  # 示例输入数据 (OPM, OBS)
├── external/               # 外部依赖包装器 (dace_wrapper)
├── fpm.toml                # Fortran Package Manager 配置
└── setup_env.sh            # 自动化环境配置脚本
```

## 🔐 Git 贡献指南

若在服务器上遇到 GitHub 推送权限问题，请使用 SSH 方式：

```bash
# 生成密钥
ssh-keygen -t ed25519 -f ~/my_ssh_keys/github_key

# 配置远程仓库
git remote set-url origin git@github.com:mayoyuu/POD_FORTRAN.git
git config core.sshCommand "ssh -i ~/my_ssh_keys/github_key -o StrictHostKeyChecking=no"

# 创建新分支开发
git checkout -b feature-new-algo
```

---

## 📝 License & Citation

本项目采用 MIT 开源协议。引用本项目请注明作者及 GitHub 地址。
