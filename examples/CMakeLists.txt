if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED)

function(make_example example_name)
    add_executable(${example_name} ${example_name}.cpp)
    target_link_libraries(${example_name} PRIVATE Boost::headers precision newton parallel runge_kutta auto_diff mpfr gmp)
    add_compile_options(-Wall -O3 -march=native)
endfunction()

function(copy_file example_name file_name)
    add_custom_command(TARGET ${example_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/${file_name}
        ${CMAKE_CURRENT_BINARY_DIR}/${file_name})
endfunction()

add_subdirectory(sitnikov)
add_subdirectory(pendulum)
add_subdirectory(three_body)
