cmake_minimum_required(VERSION 3.10)
project(DaceWrapper CXX)

# 设置 C++ 标准 (DACE 通常需要 C++11 或更高)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 未指定构建类型时默认 RelWithDebInfo（与 DACE 导出的 CMake 配置一致）
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# ==========================================
# DACE 根目录：默认 $HOME/local_dace，可用 -DDACE_ROOT=... 覆盖
# 可通过 setup_env.sh 中的 DACE_ROOT 环境变量指定
# ==========================================
if(NOT DACE_ROOT)
    if(DEFINED ENV{DACE_ROOT})
        set(DACE_ROOT "$ENV{DACE_ROOT}")
    else()
        set(DACE_ROOT "$ENV{HOME}/local_dace")
    endif()
endif()
set(DACE_ROOT "${DACE_ROOT}" CACHE PATH "DACE 安装根目录")

if(NOT EXISTS "${DACE_ROOT}/include/dace/dace.h")
    message(FATAL_ERROR
        "在 ${DACE_ROOT} 下找不到 include/dace/dace.h。\n"
        "请用 -DDACE_ROOT=/path/to/dace 指定正确的 DACE 安装路径。")
endif()

# 生成我们要的静态库：libdace_wrapper.a
add_library(dace_wrapper STATIC dace_wrapper.cpp)

# 优先使用 DACE 安装时导出的 CMake 包（daceConfig.cmake -> dace::dace / dace::dace_s）
list(APPEND CMAKE_PREFIX_PATH "${DACE_ROOT}")
find_package(dace CONFIG QUIET)

if(dace_FOUND)
    # 与 fpm.toml 中链接的 dace_s（静态库）保持一致
    target_link_libraries(dace_wrapper PUBLIC dace::dace_s)
    message(STATUS "DaceWrapper: 使用 DACE CMake 包 (${DACE_ROOT})")
else()
    # 回退：手工指定头文件与库
    target_include_directories(dace_wrapper PUBLIC "${DACE_ROOT}/include")
    target_link_libraries(dace_wrapper PUBLIC "${DACE_ROOT}/lib/libdace_s.a")
    message(STATUS "DaceWrapper: 未找到 DACE CMake 包，回退到手工路径 (${DACE_ROOT})")
endif()

message(STATUS "DaceWrapper 配置完成，准备编译！")