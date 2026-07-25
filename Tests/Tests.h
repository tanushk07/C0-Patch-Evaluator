#pragma once

#include "Evaluator.h"
#include <algorithm>
#include "Vector3.h"
#include <iostream>
#include <vector>
#include <map>
#include "ErrorStats.h"
#include "MakeSphereMesh.h"
#include <utility>
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
            double len = d.length();
            if (len <= 0.0) { ok = false; continue; }
            
            double e1 = std::abs(norms[s].normalize().dot(d - c)) / len;
            double e2 = std::abs(norms[t].normalize().dot(d + c)) / len;
            bool pass = (e1 < TOL && e2 < TOL);
            ok &= pass;

            cout << "  " << names[e]
                 << ":  |n_s.(d-c)|/|d|=" << e1
                 << "  |n_e.(d+c)|/|d|=" << e2
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
             << sampleN << " samples across the whole triangle's parameter grid)\n";

        ErrorStats stats = Evaluator::MeasureError(mesh,Error,sampleN);

        cout << "Samples:           " << stats.samples << endl;
        cout << "Max flat error:    " << stats.maxFlat << endl;
        cout << "Max Nagata error:  " << stats.maxNagata << endl;
        cout << "Avg flat error:    " << stats.avgFlat << endl;
        cout << "Avg Nagata error:  " << stats.avgNagata << endl;
        cout << "Ratio at this resolution: " << stats.ratio << "x\n";

        bool pass = (stats.maxNagata < stats.maxFlat);
        cout << "  Nagata < flat?     " << (pass ? "PASS" : "FAIL") << "\n";
        return pass;
    }
    
    static bool Test5_ResidualSweep(const char* name, const Mesh& mesh)
    {
        cout << "\nTest 5: Orthogonality residual sweep  (" << name << ")\n";
 
        double worstE0 = 0.0f, worstE1 = 0.0f;
        long long checked = 0, failed = 0, parallel = 0, antiparallel = 0;
 
        for (const auto& t : mesh.triangles)
        {
            int v[3] = {t.a, t.b, t.c};
            for (int e = 0; e < 3; ++e)
            {
                int s = v[e], f = v[(e + 1) % 3];
                Vector3 d = mesh.vertices[f] - mesh.vertices[s];
 
                CurvatureStatus st = CurvatureStatus::Ok;
                Vector3 c = Evaluator::Curvature(d, mesh.normals[s], mesh.normals[f], &st);
 
                if (st == CurvatureStatus::Parallel)     { parallel++;     continue; }
                if (st == CurvatureStatus::Antiparallel) { antiparallel++; continue; }
                double len = d.length();
                if (len <= 0.0) continue;
                double e0 = std::abs(mesh.normals[s].normalize().dot(d - c)) / len;
                double e1 = std::abs(mesh.normals[f].normalize().dot(d + c)) / len;
 
                worstE0 = std::max(worstE0, e0);
                worstE1 = std::max(worstE1, e1);
                checked++;
                if (e0 > TOL || e1 > TOL) failed++;
            }
        }
 
        cout << "  edges checked        " << checked << "\n";
        cout << "  worst |n0.(d-c)|/|d| " << worstE0 << "\n";
        cout << "  worst |n1.(d+c)|/|d| " << worstE1 << "\n";
        cout << "  flat edges (skipped) " << parallel << "\n";
        cout << "  opposed normals      " << antiparallel
             << (antiparallel ? "   <-- investigate" : "") << "\n";
        cout << "  edges over tolerance " << failed
             << "  " << (failed == 0 ? "PASS" : "FAIL") << "\n";
        return failed == 0;
    }
 
    // Point 5: closure. Edges are counted from the triangle list, never assumed.
    static bool Test6_Topology(const char* name, const Mesh& mesh, int expectedChi)
    {
        cout << "\nTest 6: Topology  (" << name << ")\n";
 
        map<pair<int,int>, int> edgeCount;
        auto addEdge = [&](int a, int b) {
            if (a > b) std::swap(a, b);
            edgeCount[{a, b}]++;
        };
 
        long long degenerate = 0, flipped = 0;
        for (const auto& t : mesh.triangles)
        {
            if (t.a == t.b || t.b == t.c || t.a == t.c) degenerate++;
            else
            {
                Vector3 e1 = mesh.vertices[t.b] - mesh.vertices[t.a];
                Vector3 e2 = mesh.vertices[t.c] - mesh.vertices[t.a];
                Vector3 gn = e1.cross(e2);
                if (0.5f * gn.length() < 1e-6f) degenerate++;
                else
                {
                    // winding: does the face normal agree with the vertex normals?
                    Vector3 avg = mesh.normals[t.a] + mesh.normals[t.b] + mesh.normals[t.c];
                    if (gn.dot(avg) < 0.0f) flipped++;
                }
            }
            addEdge(t.a, t.b); addEdge(t.b, t.c); addEdge(t.c, t.a);
        }
 
        long long boundary = 0, nonManifold = 0;
        for (const auto& e : edgeCount)
        {
            if (e.second == 1) boundary++;
            else if (e.second > 2) nonManifold++;
        }
 
        long long V = mesh.vertices.size(), E = edgeCount.size(), F = mesh.triangles.size();
        long long chi = V - E + F;
 
        cout << "  V=" << V << "  E=" << E << "  F=" << F << "\n";
        cout << "  Euler characteristic " << chi << "  (expected " << expectedChi << ")  "
             << (chi == expectedChi ? "PASS" : "FAIL") << "\n";
        cout << "  boundary edges       " << boundary << "  "
             << (boundary == 0 ? "PASS" : "FAIL") << "\n";
        cout << "  non-manifold edges   " << nonManifold << "  "
             << (nonManifold == 0 ? "PASS" : "FAIL") << "\n";
        cout << "  degenerate triangles " << degenerate << "  "
             << (degenerate == 0 ? "PASS" : "FAIL") << "\n";
        cout << "  inward-wound faces   " << flipped << "  "
             << (flipped == 0 ? "PASS" : "FAIL") << "\n";
 
        return chi == expectedChi && boundary == 0 && nonManifold == 0
               && degenerate == 0 && flipped == 0;
    }
};
