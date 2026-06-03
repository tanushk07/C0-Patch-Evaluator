#include "Evaluator.h"
#include "Tests.h"
#include "MakeSphereMesh.h"
#include "MakeTorusMesh.h"
#include "Vector3.h"
#include <iostream>

using namespace std;

int main()
{
    //Phase 1: Correctness tests on a single C0 Nagata patch. A small symmetric triangle on a sphere of radius R is enough to check the curvature parameter, vertex recovery, and normals.
    const float R = 1.0f;

    vector<Vector3> verts = { Vector3(1,0,0), Vector3(0,1,0), Vector3(0,0,1) };
    vector<Vector3> norms = verts;                 // on a unit sphere, normal == position
    coefficients c = Evaluator::MakeCoefficients(verts, norms);

    ErrorMetric sphereError = [R](Vector3 p) {
        return std::abs(p.length() - R);
    };

    Mesh sphere;
    MakeSphereMesh(R, 32, 32, sphere);

    int passed = 0, total = 4;
    if (Tests::Test1_CheckOrthogonality(verts, norms)) passed++;
    if (Tests::Test2_VertexRecovery(c, verts))         passed++;
    if (Tests::Test3_NormalRecovery(c, norms))         passed++;
    if (Tests::Test4_Accuracy(sphere, sphereError))    passed++;
    cout << "\nResults: " << passed << "/" << total << " tests passed\n";

    // Phase 2: Quantitative comparison across mesh resolutions.
    // Prints, for sphere and torus, how flat-triangle error and
    // Nagata-patch error shrink as the control mesh is refined.
    Evaluator::RunSimplificationExperiment();

    // Phase 3: Visual artifacts (OBJ with per-vertex error color).
    // Flat and Nagata share ONE error scale so the colors are
    // directly comparable: green = on the surface, red = far off.
    const int displayN = 6;                        // samples per patch edge, display only

    // Sphere maps, scaled by the flat mesh's worst error
    ErrorStats sphereStats = Evaluator::MeasureError(sphere, sphereError, 20);
    Evaluator::WriteFlatErrorObj  (sphere, "sphere_flat_error.obj",   displayN, sphereError, sphereStats.maxFlat);
    Evaluator::WriteNagataErrorObj(sphere, "sphere_nagata_error.obj", displayN, sphereError, sphereStats.maxFlat);

    // Torus maps, same idea with the torus metric
    const float Rc = 5.0f, Rt = 3.0f;
    ErrorMetric torusError = [Rc, Rt](Vector3 p) {
        float q = std::sqrt(p.x*p.x + p.y*p.y) - Rc;   // distance from centerline ring (in xy-plane)
        float d = std::sqrt(q*q + p.z*p.z);            // distance from the tube's center circle
        return std::abs(d - Rt);                       // how far off the tube surface
    };

    Mesh torus;
    MakeTorusMesh(Rc, Rt, 20, 10, torus);
    ErrorStats torusStats = Evaluator::MeasureError(torus, torusError, 20);
    Evaluator::WriteFlatErrorObj  (torus, "torus_flat_error.obj",   displayN, torusError, torusStats.maxFlat);
    Evaluator::WriteNagataErrorObj(torus, "torus_nagata_error.obj", displayN, torusError, torusStats.maxFlat);

    cout << "\nWrote 4 error-mapped OBJ files (sphere + torus, flat + Nagata).\n";

    return (passed == total) ? 0 : 1;
}