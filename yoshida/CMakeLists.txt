add_library(yoshida INTERFACE src/yoshida.cpp)
find_package(OpenMP REQUIRED)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED)

target_link_libraries(yoshida INTERFACE runge_kutta precision Boost::headers OpenMP::OpenMP_CXX)
target_include_directories(yoshida INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
