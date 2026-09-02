add_library(precision INTERFACE)
if(DEFINED USE_MULTIPRECISION)
    target_compile_definitions(precision INTERFACE USE_MULTIPRECISION)
endif()

target_include_directories(precision INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
