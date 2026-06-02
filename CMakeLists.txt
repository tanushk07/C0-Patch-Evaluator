cmake_minimum_required(VERSION 3.20)

project(c0Evaluator)

set(C0EVAL_SOURCES
        Core/Main.cpp
        Core/Evaluator.cpp
        Core/Evaluator.h
        Tests/Tests.cpp
        Utils/Vector3.h
)

add_executable(C0Eval ${C0EVAL_SOURCES})

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${C0EVAL_SOURCES})

target_compile_features(C0Eval PRIVATE cxx_std_17)

target_include_directories(C0Eval PRIVATE
        Core
        Utils
)
