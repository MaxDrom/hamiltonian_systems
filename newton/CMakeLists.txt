project(newton)
add_library(newton STATIC src/newton.cpp)
find_package(OpenMP REQUIRED)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED)

target_link_libraries(newton PRIVATE precision Boost::headers OpenMP::OpenMP_CXX)
target_include_directories(newton PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
