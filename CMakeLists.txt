cmake_minimum_required(VERSION 3.30)
project(filtersnp)

set(CMAKE_CXX_STANDARD 26)

find_package(spdlog REQUIRED)

# 添加 src目录
file(GLOB SRC_FILES src/*.cpp)

# 添加可执行文件
add_executable(filtersnp main.cpp ${SRC_FILES})

# 添加include目录，让编译器能找到 创建的模块
target_include_directories(filtersnp PRIVATE include)

# 引入hts库
target_link_libraries(filtersnp PRIVATE spdlog::spdlog hts)