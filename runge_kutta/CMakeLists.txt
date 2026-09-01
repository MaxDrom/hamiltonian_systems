add_library(runge_kutta INTERFACE)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

find_package(Boost REQUIRED)
target_link_libraries(runge_kutta INTERFACE precision Boost::headers newton)
target_include_directories(runge_kutta INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
