project(newton)
add_library(newton INTERFACE)
find_package(OpenMP REQUIRED)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED)

target_link_libraries(newton INTERFACE precision Boost::headers)
target_include_directories(newton INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
