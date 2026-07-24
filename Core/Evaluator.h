#pragma once

#include <functional>
#include <vector>
#include <utility>
#include "Vector3.h"
using namespace std;

inline constexpr float TOL = 1e-5f;
inline constexpr float pi = 3.14159265358979323846f;
using ErrorMetric = function<float(Vector3)>;
struct ErrorStats;
struct Mesh;

enum class CurvatureStatus { Ok, Parallel, Antiparallel };
struct coefficients
{
    Vector3 c00,c11,c10,c01,c20,c02;
    CurvatureStatus status[3];
};
    
class Evaluator
{
    public:
    static Vector3 Curvature(Vector3 d, Vector3 nStart, Vector3 nEnd, CurvatureStatus *status = nullptr);
    static coefficients MakeCoefficients( const vector<Vector3>& vertices, const vector<Vector3>& normals);
    static Vector3 EvalPatch(const coefficients& c, float eta, float zeta);
    static Vector3 EvalNormal(const coefficients& c, float eta, float zeta);
    static Vector3 EvalFlat(const std::vector<Vector3>& verts, float eta, float zeta);
    static void WriteFlatObj(const Mesh& sphereMesh, const char* filename);
    static ErrorStats MeasureError(const Mesh& Mesh, const ErrorMetric &Error, int N);
    static void RunSimplificationExperiment();
    static void WriteFlatErrorObj(const Mesh& Mesh, const char* filename, int N,const ErrorMetric& error, float maxError);
    static void WriteNagataErrorObj( const Mesh& Mesh, const char * filename, int N, const ErrorMetric& error, float maxError);
};

