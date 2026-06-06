cmake_minimum_required(VERSION 3.20)

project(c0Evaluator)

set(C0EVAL_SOURCES
        Core/Main.cpp
        Core/Evaluator.cpp
        Core/Evaluator.h
        Tests/Tests.h
        Utils/Vector3.h
        Utils/MakeSphereMesh.h
        Utils/MakeTorusMesh.h
        Utils/Mesh.h
        Utils/ErrorStats.h
)

set(C0EVAL_DOCS
        README.md
        Docs/comparison.png
        Docs/Notes.md
        Docs/edge_curve.svg
        Docs/param_triangle.svg
        Docs/torus_metric.svg
)

add_executable(C0Eval ${C0EVAL_SOURCES})

add_custom_target(C0EvalDocs SOURCES ${C0EVAL_DOCS})

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES
        ${C0EVAL_SOURCES}
        ${C0EVAL_DOCS}
)

target_compile_features(C0Eval PRIVATE cxx_std_17)

target_include_directories(C0Eval PRIVATE
        Core
        Utils
        Docs
        Tests
)
