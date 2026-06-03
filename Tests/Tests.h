#pragma once

#include "Evaluator.h"
#include <algorithm>
#include "Vector3.h"
#include <iostream>
#include <vector>

#include "ErrorStats.h"
#include "MakeSphereMesh.h"

using namespace std;
class Tests
{
public:
    static bool Test1_CheckOrthogonality(vector<Vector3>& verts, vector<Vector3>& norms)
    {
        //Tangent of the plane at vertices should be perpendicular to its normal. This checks the correctness of Curvature Parameter.
        cout << "\n Test 1: Orthogonality check\n";
        int edges[][2] = {{0,1}, {1,2}, {0,2}};
        const char* names[] = {"v0->v1", "v1->v2", "v0->v2"};
        bool ok = true;

        for (int e = 0; e < 3; e++)
        {
            int s = edges[e][0], t = edges[e][1];
            Vector3 d = verts[t] - verts[s];
            Vector3 c = Evaluator::Curvature(d, norms[s], norms[t]);

            float e1 = std::abs(norms[s].dot(d - c));
            float e2 = std::abs(norms[t].dot(d + c));
            bool pass = (e1 < TOL && e2 < TOL);
            ok &= pass;

            cout << "  " << names[e]
                 << ":  |n_s.(d-c)|=" << e1
                 << "  |n_e.(d+c)|=" << e2
                 << "  " << (pass ? "OK" : "FAIL") << "\n";
        }
        return ok;
    }

    static bool Test2_VertexRecovery(const coefficients& c, vector<Vector3>& verts)
    {
        //Path Eval should throw the same vertice at respective boundaries of eta and zeta
        
        cout << "\nTest 2: Vertex recovery\n";
        float params[][2] = {{0,0}, {0,1}, {1,1}};
        const char* names[] = {"v0 (0,0)", "v1 (0,1)", "v2 (1,1)"};
        bool ok = true;

        for (int i = 0; i < 3; i++)
        {
            Vector3 got = Evaluator::EvalPatch(c, params[i][0], params[i][1]);
            float err = (got - verts[i]).length();
            bool pass = err < TOL;
            ok &= pass;
            cout << "  " << names[i] << ":  error=" << err
                 << "  " << (pass ? "PASS" : "FAIL") << "\n";
        }
        return ok;
    }

    static bool Test3_NormalRecovery(const coefficients& c, vector<Vector3>& norms)
    {
        //Perpendicular resultant of delta eta and delta zeta should be same in direction as the x vertice normal at (eta,zeta) = x
        
        cout << "\nTest 3: Normal recovery at vertices\n";
        float params[][2] = {{0,0}, {0,1}, {1,1}};
        const char* names[] = {"n0 (0,0)", "n1 (0,1)", "n2 (1,1)"};
        bool ok = true;

        for (int i = 0; i < 3; i++)
        {
            Vector3 got = Evaluator::EvalNormal(c, params[i][0], params[i][1]);
            Vector3 exp = norms[i].normalize();
            float err = 1.0f - std::abs(got.dot(exp));  // 0 when parallel. This test only checks that normal lies on the same line as the true one, not that it points the same way.
            bool pass = err < TOL;
            ok &= pass;
            cout << "  " << names[i] << ":  |1-|dot||=" << err
                 << "  " << (pass ? "PASS" : "FAIL") << "\n";
        }
        return ok;
    }

    static bool Test4_Accuracy(const Mesh& mesh, const ErrorMetric &Error, int sampleN = 20)
    {
        cout << "\nTest 4: Mesh accuracy  "
             << ", " << mesh.triangles.size() << " triangles, "
             << sampleN << " samples/edge per triangle)\n";

        ErrorStats stats = Evaluator::MeasureError(mesh,Error,sampleN);

        cout << "Samples:           " << stats.samples << endl;
        cout << "Max flat error:    " << stats.maxFlat << endl;
        cout << "Max Nagata error:  " << stats.maxNagata << endl;
        cout << "Avg flat error:    " << stats.avgFlat << endl;
        cout << "Avg Nagata error:  " << stats.avgNagata << endl;
        cout << "Improvement ratio: " << stats.ratio << "x  (scale-invariant)\n";

        bool pass = (stats.maxNagata < stats.maxFlat);
        cout << "  Nagata < flat?     " << (pass ? "PASS" : "FAIL") << "\n";
        return pass;
    }
};
