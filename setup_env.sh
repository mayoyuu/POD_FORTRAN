#!/bin/bash
# POD_Fortran 开发环境配置脚本
# 使用方法: source ./setup_env.sh
#
# 请根据实际安装路径修改以下路径

# 获取脚本所在目录的绝对路径
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# DACE 微分代数库
export DACE_ROOT="${DACE_ROOT:-$HOME/local_dace}"
export CPLUS_INCLUDE_PATH="$DACE_ROOT/include:$CPLUS_INCLUDE_PATH"
export LIBRARY_PATH="$DACE_ROOT/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$DACE_ROOT/lib:$LD_LIBRARY_PATH"

# SPICE Toolkit (Fortran) —— 提供 libspicelib.a
export SPICE_ROOT="${SPICE_ROOT:-$HOME/spice_fortran/toolkit}"
export LIBRARY_PATH="$SPICE_ROOT/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$SPICE_ROOT/lib:$LD_LIBRARY_PATH"

# OpenBLAS 线性代数库
export LIBRARY_PATH="$HOME/local_openblas/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$HOME/local_openblas/lib:$LD_LIBRARY_PATH"

# 项目本地库目录
export LIBRARY_PATH="$PROJECT_ROOT/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$PROJECT_ROOT/lib:$LD_LIBRARY_PATH"

echo "✅ POD_Fortran 开发环境已激活"
echo "   DACE:      $DACE_ROOT"
echo "   OpenBLAS:  $HOME/local_openblas"
echo "   SPICE:     $SPICE_ROOT"
echo "   Project:   $PROJECT_ROOT"
