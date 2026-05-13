cmake_minimum_required(VERSION 3.10)
project(DaceWrapper CXX)

# 设置 C++ 标准 (DACE 通常需要 C++11 或更高)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 告诉 CMake 你的 DACE 根目录在哪里
set(DACE_ROOT "/home/songyu/dace")

# ==========================================
# 核心：将所有 DACE 相关的头文件路径一网打尽
# ==========================================
include_directories(
    ${DACE_ROOT}/interfaces/cxx/include  # 解决 <dace/dace.h>
    ${DACE_ROOT}/core/include            # 解决 <dace/dacecore.h>
    ${DACE_ROOT}/include                 # 备用核心路径
    ${DACE_ROOT}/build/core/include           # 解决 <dace/config.h> (如果你之前编译在 build 下)
)

# 生成我们要的静态库：libdace_wrapper.a
add_library(dace_wrapper STATIC dace_wrapper.cpp)

message(STATUS "DaceWrapper 配置完成，准备编译！")