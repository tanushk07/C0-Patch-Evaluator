#include "Evaluator.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <fstream>
#include "MakeSphereMesh.h"
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




void Evaluator::WriteNagataObj( const Mesh& Mesh, const char * filename, int N)
{
	ofstream file(filename);
    long long vertexcounter =0;
    
    for (auto p : Mesh.triangles)
    {
        vector<Vector3> verts = {Mesh.vertices[p.a], Mesh.vertices[p.b], Mesh.vertices[p.c]};
        vector<Vector3> norms = {Mesh.normals[p.a], Mesh.normals[p.b], Mesh.normals[p.c]};
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
                Vector3 patchNormal = EvalNormal(c, eta, zeta);
                
                Vector3 ref = norms[0] * (1.0f - zeta)
                            + norms[1] * (zeta - eta)
                            + norms[2] * eta;
                if (patchNormal.dot(ref) < 0.0f) patchNormal = -patchNormal;
                file << "v " << patchVertex.x << " " << patchVertex.y << " " << patchVertex.z << "\n";
                file << "vn " << patchNormal.x << " " << patchNormal.y << " " << patchNormal.z << "\n";
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
                file << "f " << (A+1) << "//" << (A+1) << " " << (C+1) << "//" << (C+1) << " " << (B+1) << "//" << (B+1) << "\n";


                if (i < j - 1)              
                {
                    int D = localInd[j - 1][i + 1];
                    file << "f " << (B+1) << "//" << (B+1) << " " << (C+1) << "//" << (C+1) << " " << (D+1) << "//" << (D+1) << "\n";
                }
            }
        }
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
            float zeta = static_cast<float>(j) / N;
            for (int i = 0; i <= j; i++)
            {
                float eta = static_cast<float>(i) / N;

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
        Mesh mesh;
        MakeSphereMesh(radius, res, res,mesh);
        pair<float,float> p = MeasureError(mesh, radius, 20);   // sampleN = 20 for accurate measurement
        cout << res << "\t" << mesh.triangles.size()
             << "\t\t" << p.first << "\t" << p.second << "\n";
    }
}

