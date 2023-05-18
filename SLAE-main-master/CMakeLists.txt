cmake_minimum_required(VERSION 3.0.0)
project(slae)
#[[//add_compile_definitions(-Wall -Werror)]]
set(CMAKE_CXX_STANDART 17)

add_subdirectory(src)
add_subdirectory(gtest)

if (WITH_TESTS)
    add_compile_definitions(COMPILE_WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
    message("TESTING: ON")
else ()
    message("TESTING: OFF")
endif ()

#[[
//target_link_libraries(cmake_demo fmt::fmt)]]
