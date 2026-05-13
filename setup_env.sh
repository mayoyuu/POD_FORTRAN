#!/bin/bash
# 获取脚本所在目录的绝对路径
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

export LIBRARY_PATH="$PROJECT_ROOT/lib:$LIBRARY_PATH"
# # 添加项目的 lib 目录，以及刚安装的 OpenBLAS 目录
# export LIBRARY_PATH="$PROJECT_ROOT/lib:/homes/wxsong/local_openblas/lib:$LIBRARY_PATH"
# # 如果是动态链接，顺便把 LD_LIBRARY_PATH 也加上防患于未然
# export LD_LIBRARY_PATH="/homes/wxsong/local_openblas/lib:$LD_LIBRARY_PATH"

echo "✅ POD_Fortran 开发环境已激活"


# 脚本运行需要是source ./sh