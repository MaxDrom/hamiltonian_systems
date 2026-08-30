find_package(doctest REQUIRED)
include(doctest)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED)

function(make_test test_name)
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE gmp mpfr Boost::headers doctest::doctest auto_diff precision newton runge_kutta parallel)    
    doctest_discover_tests(${test_name})
endfunction()
add_subdirectory(auto_diff)
