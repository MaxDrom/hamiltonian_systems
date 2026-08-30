add_library(parallel INTERFACE)
find_package(OpenMP REQUIRED)
target_link_libraries(parallel INTERFACE OpenMP::OpenMP_CXX)
target_include_directories(parallel INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})