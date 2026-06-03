#pragma once

#include "Evaluator.h"
#include <algorithm>
#include "Vector3.h"
#include <iostream>
#include <vector>
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

    static bool Test4_Accuracy(const Mesh& mesh, float distToSurface, int sampleN = 20)
    {
        cout << "\nTest 4: Mesh accuracy (distToSurface " << distToSurface
             << ", " << mesh.triangles.size() << " triangles, "
             << sampleN << " samples/edge per triangle)\n";

        float  maxFlat = 0, maxNagata = 0;
        double sumFlat = 0, sumNagata = 0;   // double: many samples accumulate, float drifts
        long long count = 0;

        for (const Tri& tri : mesh.triangles)
        {
            vector<Vector3> verts = {mesh.vertices[tri.a], mesh.vertices[tri.b], mesh.vertices[tri.c]};
            vector<Vector3> norms = {mesh.normals[tri.a],  mesh.normals[tri.b],  mesh.normals[tri.c]};

            // skip zero-area (pole) triangles so we measure exactly what gets rendered
            float area2 = ((verts[1] - verts[0]).cross(verts[2] - verts[0])).length();
            if (area2 < 1e-6f) continue;

            coefficients c = Evaluator::MakeCoefficients(verts, norms);

            for (int j = 0; j <= sampleN; j++)
            {
                float zeta = static_cast<float>(j) / sampleN;
                for (int i = 0; i <= j; i++)
                {
                    float eta = static_cast<float>(i) / sampleN;

                    float errN = std::abs(Evaluator::EvalPatch(c, eta, zeta).length() - distToSurface);
                    float errF = std::abs(Evaluator::EvalFlat(verts, eta, zeta).length() - distToSurface);

                    maxNagata = std::max(maxNagata, errN);
                    maxFlat   = std::max(maxFlat,   errF);
                    sumNagata += errN;
                    sumFlat   += errF;
                    count++;
                }
            }
        }

        float ratio = (maxNagata > 1e-10f) ? maxFlat / maxNagata : 0;

        cout << "Samples:           " << count << endl;
        cout << "Max flat error:    " << maxFlat << endl;
        cout << "Max Nagata error:  " << maxNagata << endl;
        cout << "Avg flat error:    " << (count ? sumFlat / count : 0) << endl;
        cout << "Avg Nagata error:  " << (count ? sumNagata / count : 0) << endl;
        cout << "Improvement ratio: " << ratio << "x  (scale-invariant)\n";

        bool pass = (maxNagata < maxFlat);
        cout << "  Nagata < flat?     " << (pass ? "PASS" : "FAIL") << "\n";
        return pass;
    }
};
