# c0Evaluator

C++/CMake project for evaluating flat triangle patches against Nagata patches on a sphere.

## Build

```powershell
cmake -S . -B cmake-build-debug
cmake --build cmake-build-debug --config Debug
```

## Run

```powershell
.\cmake-build-debug\Debug\C0Eval.exe
```

The current entry point is `Core/Main.cpp`, which runs `Evaluator::RunSimplificationExperiment()`.

## Project Layout

- `Core/` - evaluator implementation and entry point
- `Tests/` - test helper code
- `Utils/` - vector math utilities

Generated build files and mesh outputs are ignored by Git.
