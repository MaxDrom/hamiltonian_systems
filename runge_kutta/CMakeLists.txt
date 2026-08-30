project(runge_kutta)
add_library(runge_kutta STATIC src/rk4.cpp src/brk.cpp src/utilities.cpp)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(Boost REQUIRED)
target_link_libraries(runge_kutta PRIVATE precision Boost::headers newton)
target_include_directories(runge_kutta PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
