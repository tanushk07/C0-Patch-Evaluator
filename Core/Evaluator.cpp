#include "Evaluator.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <fstream>
using namespace std;

Vector3 Evaluator::Curvature(Vector3 d, Vector3 nStart, Vector3 nEnd)
{
    float cs = nStart.dot(nEnd);
    if (std::abs(cs) > 1.0f - 1e-6f) return Vector3(0, 0, 0);

    float a = nStart.dot(d);
    float b = nEnd.dot(d);

    return ((a + cs * b) * nStart - (cs * a + b) * nEnd) / (1 - cs * cs);
}

coefficients Evaluator::MakeCoefficients(const vector<Vector3>& vertices, const vector<Vector3>& normals)
{
    Vector3 curvA = Curvature(vertices[1] - vertices[0], normals[0], normals[1]);
    Vector3 curvB = Curvature(vertices[2] - vertices[1], normals[1], normals[2]);
    Vector3 curvC = Curvature(vertices[2] - vertices[0], normals[0], normals[2]);

    coefficients co;
    co.c00 = vertices[0];
    co.c01 = (vertices[1] - vertices[0]) - curvA;
    co.c02 = curvA;
    co.c10 = (vertices[2] - vertices[1]) + curvA - curvC;
    co.c20 = curvB;
    co.c11 = curvC - curvA - curvB;

    return co;
}

Vector3 Evaluator::EvalPatch(const coefficients& c, float eta, float zeta)
{
    return c.c00
         + c.c10 * eta  + c.c01 * zeta
         + c.c11 * eta * zeta
         + c.c20 * eta * eta + c.c02 * zeta * zeta;
}

Vector3 Evaluator::EvalNormal(const coefficients& c, float eta, float zeta)
{
    Vector3 dEta  = c.c10 + c.c11 * zeta + 2.0f * c.c20 * eta;
    Vector3 dZeta = c.c01 + c.c11 * eta  + 2.0f * c.c02 * zeta;
    return (dEta.cross(dZeta)).normalize();
}

Vector3 Evaluator::EvalFlat(vector<Vector3>& verts, float eta, float zeta)
{
    return verts[0] + (verts[1] - verts[0]) * zeta + (verts[2] - verts[1]) * eta;
}

void Evaluator::WriteFlatObj(const Mesh &sphereMesh, const char* filename)
{
    std::ofstream file(filename);
    
    for (auto v: sphereMesh.vertices)
    {
        file <<"v "<< v.x << " " << v.y << " " << v.z << "\n";
    }
    for (auto n: sphereMesh.triangles)
    {
        file <<"f "<< n.a + 1 << " " << n.b +1  << " " << n.c +1 << "\n";
    }
}

pair<float,float> Evaluator::MeasureError(const Mesh &Mesh, float radius, int N)
{
    float maxFlat =0, maxNagata = 0;
    for (auto T: Mesh.triangles)
    {
        vector<Vector3> verts = {Mesh.vertices[T.a], Mesh.vertices[T.b], Mesh.vertices[T.c]};
        vector<Vector3> norms = {Mesh.normals[T.a], Mesh.normals[T.b], Mesh.normals[T.c]};
        coefficients c = MakeCoefficients(verts, norms);
        for (int j = 0; j <= N; j++)
        {
            float zeta = (float)j / N;
            for (int i = 0; i <= j; i++)
            {
                float eta = (float)i / N;

                Vector3 pN = EvalPatch(c, eta, zeta);
                Vector3 pF = EvalFlat(verts, eta, zeta);
                
                maxNagata  = max(maxNagata, abs(pN.length() - radius));
                maxFlat = max(maxFlat, abs(pF.length() - radius));
            }
        }
    }
    return {maxFlat,maxNagata};
}
void Evaluator::RunSimplificationExperiment()
{
    float radius = 1.0f;
    cout << "res\ttriangles\tmaxFlat\t\tmaxNagata\n";
    for (int res : {4, 8, 16, 32, 64})
    {
        Mesh mesh = MakeSphereMesh(radius, res, res);
        pair<float,float> p = MeasureError(mesh, radius, 20);   // sampleN = 20 for accurate measurement
        cout << res << "\t" << mesh.triangles.size()
             << "\t\t" << p.first << "\t" << p.second << "\n";
    }
}
void Evaluator::WriteNagataObj(Mesh& sphereMesh, const char * filename, int N)
{
	ofstream file(filename);
    long long vertexcounter =0;
    
    for (auto p : sphereMesh.triangles)
    {
        vector<Vector3> verts = {sphereMesh.vertices[p.a], sphereMesh.vertices[p.b], sphereMesh.vertices[p.c]};
        vector<Vector3> norms = {sphereMesh.normals[p.a], sphereMesh.normals[p.b], sphereMesh.normals[p.c]};
        coefficients c = MakeCoefficients(verts, norms);
        
        vector<vector<int>> localInd(N+1);
        
        for (int i=0;i<=N;++i)
        {
            float zeta = static_cast<float>(i) / N;
            localInd[i].resize(i+1);
            for (int j = 0; j <=i; ++j)
            {
                float eta = static_cast<float>(j) / N;
                Vector3 patchVertex = EvalPatch(c, eta, zeta);
                file << "v " << patchVertex.x << " " << patchVertex.y << " " << patchVertex.z << "\n";
                localInd[i][j]=vertexcounter++;
            }
        }
        
        for (int j = 1; j <= N; ++j)
        {
            for (int i = 0; i < j; ++i)
            {
                int A = localInd[j][i];
                int B = localInd[j - 1][i];
                int C = localInd[j][i + 1];
                file << "f " << (A + 1) << " " << (B + 1) << " " << (C + 1) << "\n";

                if (i < j - 1)              
                {
                    int D = localInd[j - 1][i + 1];
                    file << "f " << (B + 1) << " " << (D + 1) << " " << (C + 1) << "\n";
                }
            }
        }
    }

}


Mesh Evaluator::MakeSphereMesh(float radius, int rings, int sectors)
{
    Mesh sphereMesh;

    for (int i = 0; i <= rings; ++i)
    {
        for (int j = 0; j <= sectors; ++j)
        {
            float theta = pi * i / rings;
            float phi = 2.0f * pi * j / sectors;
            float x = radius * std::sin(theta) * std::cos(phi);
            float y = radius * std::sin(theta) * std::sin(phi);
            float z = radius * std::cos(theta);

            Vector3 p = {x, y, z};
            sphereMesh.vertices.push_back(p);
            sphereMesh.normals.push_back(p.normalize());
        }
    }

    for (int i = 0; i < rings; ++i)
    {
        for (int j = 0; j < sectors; ++j)
        {
            int row1 = i * (sectors + 1);
            int row2 = (i + 1) * (sectors + 1);

            int i0 = row1 + j;
            int i1 = row1 + j + 1;
            int i2 = row2 + j;
            int i3 = row2 + j + 1;

            sphereMesh.triangles.push_back({i0, i2, i1});
            sphereMesh.triangles.push_back({i1, i2, i3});
        }
    }

    return sphereMesh;
}
