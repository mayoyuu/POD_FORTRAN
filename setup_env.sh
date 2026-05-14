#!/bin/bash
# POD_Fortran 开发环境配置脚本
# 使用方法: source ./setup_env.sh
#
# 请根据实际安装路径修改以下路径

# 获取脚本所在目录的绝对路径
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# DACE 微分代数库
export CPLUS_INCLUDE_PATH="$HOME/local_dace/include:$CPLUS_INCLUDE_PATH"
export LIBRARY_PATH="$HOME/local_dace/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$HOME/local_dace/lib:$LD_LIBRARY_PATH"

# OpenBLAS 线性代数库
export LIBRARY_PATH="$HOME/local_openblas/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$HOME/local_openblas/lib:$LD_LIBRARY_PATH"

# 项目本地库目录
export LIBRARY_PATH="$PROJECT_ROOT/lib:$LIBRARY_PATH"
export LD_LIBRARY_PATH="$PROJECT_ROOT/lib:$LD_LIBRARY_PATH"

echo "✅ POD_Fortran 开发环境已激活"
echo "   DACE:      $HOME/local_dace"
echo "   OpenBLAS:  $HOME/local_openblas"
echo "   Project:   $PROJECT_ROOT"
