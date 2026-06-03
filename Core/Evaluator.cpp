#include "Evaluator.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include "Vector3.h"
#include <utility>
#include <fstream>
#include <MakeTorusMesh.h>

#include "ErrorStats.h"
#include "MakeSphereMesh.h"
using namespace std;

Vector3 Evaluator::Curvature(Vector3 d, Vector3 nStart, Vector3 nEnd)
{
    float cs = nStart.dot(nEnd);
    if (std::abs(cs) > 1.0f - 1e-6f) return {0, 0, 0};

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

Vector3 Evaluator::EvalFlat(const vector<Vector3>& verts, float eta, float zeta)
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
void Evaluator::WriteFlatErrorObj(const Mesh &Mesh, const char* filename, int N,const ErrorMetric& error, float maxError)
{
    std::ofstream file(filename);
    long long vertexcounter = 0;
    for (auto v: Mesh.triangles)
    {
        vector<Vector3> verts = {Mesh.vertices[v.a], Mesh.vertices[v.b], Mesh.vertices[v.c]};
        vector<vector<int>> localInd(N+1);
        for (int i = 0; i <= N; ++i)
        {
            float zeta = static_cast<float>(i) / N;
            localInd[i].resize(i+1);
            for (int j = 0; j <= i; ++j)
            {
                float eta = static_cast<float>(j) / N;
                Vector3 p = EvalFlat(verts,eta,zeta);
                float err = error(p);
                float t = min(max(err/maxError,0.f),1.f);
                file << "v " << p.x << " " << p.y << " " << p.z << " "<< t <<" "<< (1.f - t) << " 0\n";
                localInd[i][j]=vertexcounter++;
            }
        }
        
        for (int i = 1; i <= N; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                int A = localInd[i][j], B = localInd[i-1][j], C = localInd[i][j+1];
                file << "f "<< A+1 << " " << C+1 << " " << B+1 << "\n";
                if (j < i-1) 
                { 
                    int D = localInd[i-1][j+1];
                    file << "f " << B+1 << " " << C+1 << " " << D+1 << "\n";
                }
            }
            
        }
    }
}


void Evaluator::WriteNagataErrorObj( const Mesh& Mesh, const char * filename, int N, const ErrorMetric& error, float maxError)
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
                
                float err = error(patchVertex);
                float t = min(max(err/maxError,0.f),1.f);
                float r = t, g=1.f-t, b=0.f;
                file << "v " << patchVertex.x << " " << patchVertex.y << " " << patchVertex.z 
                << " " << r <<" "<<g<<" "<<b << "\n";
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


ErrorStats Evaluator::MeasureError(const Mesh &Mesh,const ErrorMetric &Error , int N)
{
    float  maxFlat = 0, maxNagata = 0;
    double sumFlat = 0, sumNagata = 0;   // double: many samples accumulate, float drifts
    long long count = 0;
    
    for (auto T: Mesh.triangles)
    {
        vector<Vector3> verts = {Mesh.vertices[T.a], Mesh.vertices[T.b], Mesh.vertices[T.c]};
        vector<Vector3> norms = {Mesh.normals[T.a], Mesh.normals[T.b], Mesh.normals[T.c]};
        
        // skip zero-area (pole) triangles so we measure exactly what gets rendered
        float area2 = ((verts[1] - verts[0]).cross(verts[2] - verts[0])).length();
        if (area2 < 1e-6f) continue;
        
        coefficients c = MakeCoefficients(verts, norms);
        for (int j = 0; j <= N; j++)
        {
            float zeta = static_cast<float>(j) / N;
            for (int i = 0; i <= j; i++)
            {
                float eta = static_cast<float>(i) / N;

                Vector3 pN = EvalPatch(c, eta, zeta);
                Vector3 pF = EvalFlat(verts, eta, zeta);
                
                float errF = Error(pF);
                float errN = Error(pN);
                count++;
                sumFlat += errF;
                sumNagata += errN;
                maxNagata  = max(maxNagata, abs(errN));
                maxFlat = max(maxFlat, abs(errF));
            }
        }
    }
    
    float ratio = (maxNagata > 1e-10f) ? maxFlat / maxNagata : 0;
    float avgFlat = count ? sumFlat / count : 0;
    float avgNagata = count ? sumNagata / count : 0;
    
    ErrorStats stats;
    stats.ratio = ratio;
    stats.maxFlat = maxFlat;
    stats.maxNagata = maxNagata;
    stats.avgFlat = avgFlat;
    stats.avgNagata = avgNagata;
    stats.samples = count;
    
    return stats;
}
void Evaluator::RunSimplificationExperiment()
{
    float radius = 1.0f;
    cout<<"Sphere:\n";
    cout << "res\ttriangles\tmaxFlat\t\tmaxNagata\n";
    for (int res : {4, 8, 16, 32, 64})
    {
        Mesh mesh;
        MakeSphereMesh(radius, res, res,mesh);
        ErrorMetric sphereError = [radius](Vector3 p) {
            return std::abs(p.length() - radius);
        };
        ErrorStats p = MeasureError(mesh, sphereError, 20);   // sampleN = 20 for accurate measurement
        cout << res << "\t" << mesh.triangles.size()
             << "\t\t" << p.maxFlat << "\t" << p.maxNagata << "\n";
    }
    
    cout<<"Torus:\n";
    
    cout << "res\ttriangles\tmaxFlat\t\tmaxNagata\n";
    for (int res : {4, 8, 16, 32, 64})
    {
        Mesh mesh;
        
        float Rc = 5.0f, Rt = 3.0f;
        ErrorMetric torusError = [Rc, Rt](Vector3 p) {
            float q = std::sqrt(p.x*p.x + p.y*p.y) - Rc; // it is the distance from centerline ring, in xy-plane
            float d = std::sqrt(q*q + p.z*p.z);          // distance from the tube's center circle
            return std::abs(d - Rt);                     // how far off the tube surface
        };
        MakeTorusMesh(Rc, Rt, res,res,mesh);
        ErrorStats p = MeasureError(mesh, torusError, 20);   // sampleN = 20 for accurate measurement
        cout << res << "\t" << mesh.triangles.size()
             << "\t\t" << p.maxFlat << "\t" << p.maxNagata << "\n";
    }
}

