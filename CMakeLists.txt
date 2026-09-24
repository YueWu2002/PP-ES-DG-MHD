cmake_minimum_required(VERSION 3.15)

project(PPES-DG-MHD
    VERSION 1.0.0
    LANGUAGES Fortran
    DESCRIPTION "positivity-preserving entropy stable discontinuous Galerkin method for ideal magnetohydrodynamics equations"
)
# 自动搜索Fortran编译器并设置CMAKE_Fortran_COMPILER变量（也即编译器的绝对路径）

# 默认构建为Release模式，如果用户没有指定CMAKE_BUILD_TYPE
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    message(STATUS "No build type selected, defaulting to Release")
endif()

# 强制检查 CMAKE_BUILD_TYPE
set(ALLOWED_TYPES Debug Release)
if(NOT CMAKE_BUILD_TYPE IN_LIST ALLOWED_TYPES)
    message(FATAL_ERROR "Invalid CMAKE_BUILD_TYPE: '${CMAKE_BUILD_TYPE}'. Please use -DCMAKE_BUILD_TYPE=Debug or Release")
endif()

# 设置Fortran .mod文件输出目录
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mod)

# 匹配当前目录下的所有.f90文件，并将它们存储在SRC_LIST变量中（请勿使用GLOB_RECURE，否则会递归搜索）
file(GLOB SRC_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")

# 添加可执行文件以及源文件列表
add_executable(main ${SRC_LIST})


# 1. 设置默认值
option(USE_OMP "Enable OpenMP support" ON)
# 2. 严格校验输入 (处理用户通过 -D 传入非法字符串的情况)
set(ALLOWED_VALUES ON OFF TRUE FALSE 0 1) # 容错处理：包含常见的布尔表示
if(NOT USE_OMP IN_LIST ALLOWED_VALUES)
    message(FATAL_ERROR "Invalid USE_OMP: '${USE_OMP}'. Please use ON OFF TRUE FALSE 0 1.")
endif()

# 链接OpenMP库
if (USE_OMP)
    message(STATUS "OpenMP support enabled")
    find_package(OpenMP REQUIRED)
    target_link_libraries(main PRIVATE OpenMP::OpenMP_Fortran)
else()
    message(STATUS "OpenMP support disabled")
endif()

# 链接MKL或者Lapack库
# MKL config paras can be found in the file MKLConfig.cmake
set(MKL_LINK dynamic)
set(MKL_INTERFACE lp64)
set(MKL_THREADING sequential)
find_package(MKL CONFIG) # use CONFIG mode to find the modern CMake targets provided by MKL
if(MKL_FOUND)
    # message(STATUS "${MKL_IMPORTED_TARGETS}") # 显示找到的库的绝对路径
    target_link_libraries(main PRIVATE MKL::MKL)
else()
    message(STATUS "MKL not found. Look for Lapack instead.")
    find_package(LAPACK REQUIRED) # 强制必须找到
    target_link_libraries(main PRIVATE LAPACK::LAPACK)
endif()

# 开启链接期优化
set_property(TARGET main PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)


if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    # ifx + release/debug
    target_compile_options(main PRIVATE
        -xHost
        -r8
        -fpp
        -g
        -traceback
        -no-wrap-margin
        -qopt-report=3
        -warn declarations
        $<$<CONFIG:Release>:-O3>
        # $<$<CONFIG:Release>:-fpe0># catch FP exceptions for floating-point exceptions (invalid, zero divide, overflow) (costly)
        $<$<CONFIG:Debug>:-O0>
        $<$<CONFIG:Debug>:-check all>
        $<$<CONFIG:Debug>:-init=arrays>
        $<$<CONFIG:Debug>:-init=snan>
        $<$<CONFIG:Debug>:-fpe0>
    )

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # gfortran + release/debug
    target_compile_options(main PRIVATE
        -g
        -fdefault-real-8
        -mtune=native
        -march=native
        -m64
        -fbacktrace
        -ffree-line-length-none
        -cpp
        -fimplicit-none
        $<$<CONFIG:Release>:-O3>
        # $<$<CONFIG:Release>:-ffpe-trap=zero,invalid,overflow> # catch FP exceptions for floating-point exceptions (invalid, zero divide, overflow) (costly)
        $<$<CONFIG:Debug>:-Og>
        $<$<CONFIG:Debug>:-Wimplicit-interface>
        $<$<CONFIG:Debug>:-Wimplicit-procedure>
        $<$<CONFIG:Debug>:-finit-real=nan>
        $<$<CONFIG:Debug>:-fcheck=all>
        $<$<CONFIG:Debug>:-ffpe-trap=zero,invalid,overflow>
    )
    
else()
    message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# 输出目录与文件名
set_target_properties(main PROPERTIES
    OUTPUT_NAME main.exe
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
)

# 输出检查
message(STATUS "CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}")
message(STATUS "CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# 将完整的编译命令输出到compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)