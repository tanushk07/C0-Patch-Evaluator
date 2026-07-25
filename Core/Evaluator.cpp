#include "Evaluator.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include "Vector3.h"
#include <utility>
#include <fstream>
#include "MakeTorusMesh.h"
#include "ErrorStats.h"
#include "MakeSphereMesh.h"
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

namespace {
    // Inner width of each column (padding of one space is added on each side).
    constexpr int kColW[]   = {5, 8, 12, 12, 8, 8, 9};
    constexpr int kNumCols  = 7;

    // Format a double into a fixed-precision string so we can right-align it
    // inside a cell (streaming directly into cout makes width control fiddly).
    string Fmt(double v, int prec, bool sci)
    {
        ostringstream os;
        os << (sci ? scientific : fixed) << setprecision(prec) << v;
        return os.str();
    }

    // "+------+--------+ ... +" border line matching the column widths.
    string Separator()
    {
        string s = "+";
        for (int i = 0; i < kNumCols; ++i)
            s += string(kColW[i] + 2, '-') + "+";
        return s;
    }

    // Print one "| a | b | ... |" row from already-formatted cell strings.
    void PrintRow(const string cells[])
    {
        cout << "|";
        for (int i = 0; i < kNumCols; ++i)
            cout << " " << right << setw(kColW[i]) << cells[i] << " |";
        cout << "\n";
    }

    void PrintSweepHeader(const char* title)
    {
        const string head[kNumCols] =
            {"res", "tris", "maxFlat", "maxNagata", "flat/x", "nag/x", "ratio"};
        cout << "\n" << title << "\n"
             << Separator() << "\n";
        PrintRow(head);
        cout << Separator() << "\n";
    }

    void PrintSweepRow(int res, size_t tris, double maxFlat, double maxNagata,
                       double prevF, double prevN)
    {
        string cells[kNumCols];
        cells[0] = to_string(res);
        cells[1] = to_string(tris);
        cells[2] = Fmt(maxFlat,   3, true);
        cells[3] = Fmt(maxNagata, 3, true);
        if (prevF > 0)
        {
            cells[4] = Fmt(prevF / maxFlat, 2, false);
            cells[5] = Fmt(prevN / maxNagata, 2, false);
        }
        else
        {
            cells[4] = "-";
            cells[5] = "-";
        }
        cells[6] = Fmt(maxFlat / maxNagata, 1, false) + "x";
        PrintRow(cells);
    }

    void PrintSweepFooter()
    {
        cout << Separator() << "\n";
    }
}

// A triangle is degenerate if two corners are the same vertex, or if its area is
// negligible compared with its own size. The test is RELATIVE (area against the
// longest edge squared) so it is scale independent: a mesh scaled up by 1000
// classifies exactly the same triangles as the original.
bool Evaluator::IsDegenerate(const Mesh& mesh, const Tri& t)
{
    if (t.a == t.b || t.b == t.c || t.a == t.c) return true;

    const Vector3& v0 = mesh.vertices[t.a];
    const Vector3& v1 = mesh.vertices[t.b];
    const Vector3& v2 = mesh.vertices[t.c];

    double twiceArea = ((v1 - v0).cross(v2 - v0)).length();
    double e0 = (v1 - v0).length();
    double e1 = (v2 - v1).length();
    double e2 = (v0 - v2).length();
    double longest = max(e0, max(e1, e2));

    if (longest <= 0.0) return true;
    return twiceArea < 1e-6 * longest * longest;
}

Vector3 Evaluator::Curvature(Vector3 d, Vector3 nStart, Vector3 nEnd, CurvatureStatus* status)
{
    nStart = nStart.normalize();
    nEnd   = nEnd.normalize();
    
    double cs = nStart.dot(nEnd);
    cs = std::max(-1.0, std::min(1.0, cs));
    if (cs > 1.0f - 1e-6f)
    {
        if (status) *status = CurvatureStatus::Parallel;
        return {0, 0, 0};
    }
    
    if (cs < -1.0f + 1e-6f)
    {
        if (status) *status = CurvatureStatus::Antiparallel;
        return {0, 0, 0};
    }
    if (status) *status = CurvatureStatus::Ok;

    
    double a = nStart.dot(d);
    double b = nEnd.dot(d);

    return ((a + cs * b) * nStart - (cs * a + b) * nEnd) / (1 - cs * cs);
}

coefficients Evaluator::MakeCoefficients(const vector<Vector3>& vertices, const vector<Vector3>& normals)
{
    coefficients co;
    Vector3 curvA = Curvature(vertices[1] - vertices[0], normals[0], normals[1], &co.status[0]);
    Vector3 curvB = Curvature(vertices[2] - vertices[1], normals[1], normals[2], &co.status[1]);
    Vector3 curvC = Curvature(vertices[2] - vertices[0], normals[0], normals[2], &co.status[2]);

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

void Evaluator::WriteFlatErrorObj(const Mesh &Mesh, const char* filename, int N,const ErrorMetric& error, float maxError)
{
    std::ofstream file(filename);
    long long vertexcounter = 0;
    for (auto v: Mesh.triangles)
    {
        if (IsDegenerate(Mesh, v)) continue;
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
        if (IsDegenerate(Mesh, p)) continue;
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
    double  maxFlat = 0, maxNagata = 0;
    double sumFlat = 0, sumNagata = 0;
    double count = 0;
    
    for (auto T: Mesh.triangles)
    {
        vector<Vector3> verts = {Mesh.vertices[T.a], Mesh.vertices[T.b], Mesh.vertices[T.c]};
        vector<Vector3> norms = {Mesh.normals[T.a], Mesh.normals[T.b], Mesh.normals[T.c]};
        
        if (IsDegenerate(Mesh, T)) continue;
        
        coefficients c = MakeCoefficients(verts, norms);
        for (int j = 0; j <= N; j++)
        {
            float zeta = static_cast<float>(j) / N;
            for (int i = 0; i <= j; i++)
            {
                float eta = static_cast<float>(i) / N;

                Vector3 pN = EvalPatch(c, eta, zeta);
                Vector3 pF = EvalFlat(verts, eta, zeta);
                
                double errF = Error(pF);
                double errN = Error(pN);
                count++;
                sumFlat += errF;
                sumNagata += errN;
                maxNagata  = max(maxNagata, abs(errN));
                maxFlat = max(maxFlat, abs(errF));
            }
        }
    }
    
    double ratio = (maxNagata > 1e-10f) ? maxFlat / maxNagata : 0;
    double avgFlat = count ? sumFlat / count : 0;
    double avgNagata = count ? sumNagata / count : 0;
    
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
    // The drop columns are the real result. Flat triangles converge at order h^2
    // (error falls ~4x per refinement), Nagata patches at order h^4 (~16x). The
    // flat/Nagata ratio therefore GROWS ~4x each refinement, so no single ratio
    // describes the method: quote the convergence orders instead.
    const float radius = 1.0f;
    double prevF = 0, prevN = 0;

    PrintSweepHeader("Sphere:");
    for (int res : {4, 8, 16, 32, 64})
    {
        Mesh mesh;
        MakeSphereMesh(radius, res, res, mesh);
        ErrorMetric sphereError = [radius](Vector3 p) {
            return std::abs(p.length() - radius);
        };
        ErrorStats p = MeasureError(mesh, sphereError, 20);
        PrintSweepRow(res, mesh.triangles.size(), p.maxFlat, p.maxNagata, prevF, prevN);
        prevF = p.maxFlat;
        prevN = p.maxNagata;
    }
    PrintSweepFooter();

    prevF = 0; prevN = 0;

    PrintSweepHeader("Torus:");
    for (int res : {4, 8, 16, 32, 64})
    {
        Mesh mesh;
        float Rc = 5.0f, Rt = 3.0f;
        ErrorMetric torusError = [Rc, Rt](Vector3 p) {
            double q = std::sqrt(p.x*p.x + p.y*p.y) - Rc;
            double d = std::sqrt(q*q + p.z*p.z);
            return std::abs(d - Rt);
        };
        MakeTorusMesh(Rc, Rt, res, res, mesh);
        ErrorStats p = MeasureError(mesh, torusError, 20);
        PrintSweepRow(res, mesh.triangles.size(), p.maxFlat, p.maxNagata, prevF, prevN);
        prevF = p.maxFlat;
        prevN = p.maxNagata;
    }
    PrintSweepFooter();
}