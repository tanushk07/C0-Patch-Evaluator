# C0 Nagata Patch Evaluator

A C++ implementation and evaluation of **C0 Nagata patches** (Nagata, 2005),
a quadratic interpolation scheme that reconstructs a smooth curved surface
from a triangle's three vertices and their three normals.

This project implements the patch evaluator, verifies it against analytic
surfaces (sphere and torus), and quantifies how much more accurately Nagata
patches approximate a curved surface compared to flat triangles at the same
mesh resolution.

## What it does

- **Evaluates** a C0 Nagata patch: given 3 vertices + 3 normals, computes the
  patch surface point and surface normal at any parameter `(eta, zeta)`.
- **Verifies** the evaluator with four tests (curvature orthogonality, vertex
  recovery, normal recovery, and accuracy against a known surface).
- **Compares** flat-triangle vs Nagata error across mesh resolutions, measured
  as distance from the true analytic surface.
- **Exports** OBJ files with per-vertex error coloring (green = on surface,
  red = far from it) and Nagata surface normals, for visual inspection.

## Build

```powershell
cmake -S . -B cmake-build-debug
cmake --build cmake-build-debug --config Debug
```

## Run

```powershell
.\cmake-build-debug\Debug\C0Eval.exe
```

Running the program:
1. Runs the four correctness tests on a single patch and a sphere mesh.
2. Prints flat-vs-Nagata max error across resolutions for a sphere and a torus.
3. Writes four error-colored OBJ files (sphere + torus, flat + Nagata) to the
   working directory.

## Output files

| File | Surface | Shows |
|------|---------|-------|
| `sphere_flat_error.obj`   | Sphere | Flat-triangle error (per-vertex color) |
| `sphere_nagata_error.obj` | Sphere | Nagata-patch error + surface normals |
| `torus_flat_error.obj`    | Torus  | Flat-triangle error |
| `torus_nagata_error.obj`  | Torus  | Nagata-patch error + surface normals |

Open the colored OBJs in a viewer that reads per-vertex color (e.g.
[3dviewer.net](https://3dviewer.net), or Blender with Solid shading →
Color → Attribute). For the smooth-normal view, use the Nagata sphere with
Shade Smooth.

## Project layout

- `Core/` — evaluator (`Evaluator.cpp/.h`) and entry point (`Main.cpp`)
- `Tests/` — the four correctness tests
- `Utils/` — vector math, mesh structs, sphere/torus mesh generators, error stats

## Scope

This implements the **C0** patch (position-continuous across edges). The G1
extension (tangent-plane continuity) is not implemented; see the notes for why
C0 leaves a visible crease at shared edges and what G1 would add.

## Reference

T. Nagata, "Simple local interpolation of surfaces using normal vectors,"
*Computer Aided Geometric Design*, 2005.