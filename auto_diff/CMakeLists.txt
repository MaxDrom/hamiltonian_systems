add_library(auto_diff INTERFACE)
target_link_libraries(auto_diff INTERFACE precision)
target_include_directories(auto_diff INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
