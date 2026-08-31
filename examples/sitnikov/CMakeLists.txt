make_example(sitnikov)
add_custom_command(TARGET sitnikov POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/plot_scr
    ${CMAKE_CURRENT_BINARY_DIR}/plot_scr)
