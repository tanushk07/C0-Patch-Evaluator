# C0 Nagata Patch Evaluator

Reconstructing smooth curved surfaces from a triangle mesh using **C0 Nagata patches** (Nagata, 2005): each triangle's three vertices and their three normals define a quadratic patch that curves to follow the surface. This project implements the patch evaluator, verifies it against analytic surfaces (sphere and torus), and quantifies how much more accurately Nagata patches approximate a curved surface than flat triangles at the same mesh resolution.

## Documentation

**[Docs/Notes.md](Docs/Notes.md)** is the full write-up behind this code: the motivation, the derivation of the patch coefficients from first principles, the verification tests, the error metrics for the sphere and torus, the resolution experiment, and the findings and known limitations.

## Headline result

Error against the true surface, measured as the maximum distance from the analytic shape, across mesh resolutions. Both methods use the same control mesh and the same sampling; only flat-vs-curved differs.

**Sphere (radius 1):**

| Triangles | Flat max error | Nagata max error |
|----------:|---------------:|-----------------:|
| 32        | 0.3210         | 0.06066          |
| 128       | 0.09121        | 0.004534         |
| 512       | 0.02375        | 0.0002890        |
| 2048      | 0.006002       | 0.00001812       |
| 8192      | 0.001505       | 0.000001192      |

Flat-triangle error falls about 4x per resolution doubling (order h²); Nagata error falls about 16x (order h⁴). **Nagata at 128 triangles (0.0045) is already more accurate than flat triangles at 2048 (0.0060): roughly 16x fewer control-mesh triangles for the same accuracy, and the gap widens as the mesh refines.**

**Torus (centerline radius 5, tube radius 3):**

| Triangles | Flat max error | Nagata max error |
|----------:|---------------:|-----------------:|
| 32        | 2.460          | 0.4853           |
| 128       | 0.7652         | 0.05102          |
| 512       | 0.2063         | 0.005941         |
| 2048      | 0.05264        | 0.0007186        |
| 8192      | 0.01323        | 0.00008845       |

Same pattern as the sphere. The absolute numbers are larger only because the torus is a bigger shape (radii 5 and 3, not a unit surface); the convergence behavior is what matters.

## Visual comparison

![Flat vs Nagata error maps](Docs/comparison.png)

Green means on the true surface, red means far from it, scaled so that the flat mesh's worst error maps to full red. Flat triangles (the faceted, red-blotched shapes) sag away from the surface between vertices; Nagata patches (uniformly green) stay on it. The spiky renders show the per-vertex surface normals the patches reconstruct.

## Build and run

```powershell
cmake -S . -B cmake-build-debug
cmake --build cmake-build-debug --config Debug
.\cmake-build-debug\Debug\C0Eval.exe
```

Running the program:
1. Runs four correctness tests on a single patch and a sphere mesh.
2. Prints flat-vs-Nagata max error across resolutions for a sphere and a torus.
3. Writes four error-colored OBJ files (sphere and torus, flat and Nagata) to the working directory.

Open the colored OBJ files in a viewer that reads per-vertex color (for example [3dviewer.net](https://3dviewer.net),https://threejs.org/editor/, or Blender with Solid shading then Color set to Attribute). For the smooth-normal view, use the Nagata sphere with Shade Smooth.

## What is implemented

- C0 Nagata patch evaluation: surface point and surface normal at any parameter `(eta, zeta)`.
- Four correctness tests: curvature orthogonality, vertex recovery, normal recovery, and accuracy against a known analytic surface.
- A shape-independent error metric (passed in as a function), so the same code measures both the sphere and the torus.
- Flat-vs-Nagata error comparison across mesh resolutions.
- Error-colored OBJ export, with reconstructed normals on the Nagata meshes.

## Verification

All four tests pass on the symmetric sphere octant: curvature orthogonality, vertex recovery, normal recovery, and accuracy against the analytic sphere.

<details>
<summary>Full program output</summary>

```
 Test 1: Orthogonality check
  v0->v1:  |n_s.(d-c)|=0  |n_e.(d+c)|=0  OK
  v1->v2:  |n_s.(d-c)|=0  |n_e.(d+c)|=0  OK
  v0->v2:  |n_s.(d-c)|=0  |n_e.(d+c)|=0  OK

Test 2: Vertex recovery
  v0 (0,0):  error=0  PASS
  v1 (0,1):  error=0  PASS
  v2 (1,1):  error=0  PASS

Test 3: Normal recovery at vertices
  n0 (0,0):  |1-|dot||=0  PASS
  n1 (0,1):  |1-|dot||=0  PASS
  n2 (1,1):  |1-|dot||=0  PASS

Test 4: Mesh accuracy (2048 triangles, 20 samples across the parameter grid)
Samples:           458304
Max flat error:    0.00600249
Max Nagata error:  1.81198e-05
Avg flat error:    0.00233273
Avg Nagata error:  2.47684e-06
Improvement ratio: 331.266x
  Nagata < flat?     PASS

Results: 4/4 tests passed

Sphere:
res     triangles       maxFlat         maxNagata
4       32              0.320982        0.0606602
8       128             0.0912061       0.00453424
16      512             0.0237502       0.000288963
32      2048            0.00600249      1.81198e-05
64      8192            0.00150472      1.19209e-06

Torus:
res     triangles       maxFlat         maxNagata
4       32              2.46012         0.485282
8       128             0.765181        0.0510187
16      512             0.206297        0.00594068
32      2048            0.0526421       0.000718594
64      8192            0.0132301       8.84533e-05
```

The 331x figure in Test 4 is the max-error ratio at a single resolution (2048 triangles). The more meaningful comparison is triangle count for equal accuracy, which is about 16x.

</details>

## Project layout

- `Core/` — evaluator (`Evaluator.cpp`, `Evaluator.h`) and entry point (`Main.cpp`)
- `Tests/` — the four correctness tests
- `Utils/` — vector math, mesh structs, sphere and torus mesh generators, error-stats struct
- `Docs/`: [Notes.md](Docs/Notes.md) (full derivation and findings), plus the figures used in this README

## Scope and notes

This implements the **C0** patch only (position-continuous across shared edges). The G1 extension (tangent-plane continuity) is not implemented; C0 was the requested scope and is the foundation G1 builds on. See `NOTES.md` for the coefficient derivation, the C0-vs-G1 distinction, a sign discrepancy found between the source papers, and known limitations.

## Reference

T. Nagata, "Simple local interpolation of surfaces using normal vectors,"
*Computer Aided Geometric Design*, 2005.

Y. Nishidate et al., "Ray-tracing method for isotropic inhomogeneous
refractive-index media from arbitrary discrete input." (This paper's form of the
$\mathbf{c}_{11}$ coefficient matches the sign derived here.)
