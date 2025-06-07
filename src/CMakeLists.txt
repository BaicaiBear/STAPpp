cmake_minimum_required(VERSION 3.15)
project(stap++ LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# option(STAPPP_DEBUG "Output detailed results for debugging." OFF)
# if(STAPPP_DEBUG)
#     add_definitions(-D_DEBUG_)
# endif()

file(GLOB_RECURSE HEADER_FILES CONFIGURE_DEPENDS h/*.h)
file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS cpp/*.cpp)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/h" PREFIX "Header Files" FILES ${HEADER_FILES})
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/cpp" PREFIX "Source Files" FILES ${SOURCE_FILES})

add_executable(stap++ ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(stap++ PRIVATE h ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/eigen-3.4.0)
