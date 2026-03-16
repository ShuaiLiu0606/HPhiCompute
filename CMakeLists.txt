cmake_minimum_required(VERSION 3.16.3 FATAL_ERROR)

# 设置项目名称
project(build_my_pro)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 启用调试符号
set(CMAKE_BUILD_TYPE Debug)

# 头文件的搜索路径
include_directories(/opt/intel/oneapi/mkl/latest/include)


include_directories(${PROJECT_SOURCE_DIR}/include)

# 库的搜索路径
link_directories(/opt/intel/oneapi/mkl/latest/lib/intel64)

# 查找 OpenMP
find_package(OpenMP REQUIRED)

# 设置库输出目录
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# 设置可执行文件输出目录
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/output)

# 搜索源文件
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)

# 生成可执行文件
add_executable(my_pro ${SRC_LIST})

target_link_libraries(my_pro PRIVATE mkl_rt pthread m dl)

# 添加 stdc++fs 库的链接
#target_link_libraries(my_pro PUBLIC libmkl_rt.so pthread libm.so libdl.so stdc++fs)


# 链接 OpenMP 库
if(OpenMP_CXX_FOUND)
    target_link_libraries(my_pro PRIVATE OpenMP::OpenMP_CXX)
endif()

