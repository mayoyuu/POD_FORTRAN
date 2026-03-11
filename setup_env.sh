#!/bin/bash
# 获取脚本所在目录的绝对路径
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 将 lib 目录加入库搜索路径
export LIBRARY_PATH="$PROJECT_ROOT/lib:$LIBRARY_PATH"

echo "✅ POD_Fortran 开发环境已激活"


# 脚本运行需要是source ./sh