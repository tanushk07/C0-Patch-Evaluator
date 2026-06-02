#pragma once

#include <vector>
#include <utility>
#include "Vector3.h"
using namespace std;

inline constexpr float TOL = 1e-5f;
inline constexpr float pi = 3.14159265358979323846f;

struct coefficients
{
    Vector3 c00,c11,c10,c01,c20,c02;
};
struct Tri{int a,b,c;};
struct Mesh
{
    vector<Vector3> vertices, normals;
    vector<Tri> triangles;
};
class Evaluator
{
    public:
    static Vector3 Curvature(Vector3 d, Vector3 nStart, Vector3 nEnd);
    static coefficients MakeCoefficients( const vector<Vector3>& vertices, const vector<Vector3>& normals);
    static Vector3 EvalPatch(const coefficients& c, float eta, float zeta);
    static Vector3 EvalNormal(const coefficients& c, float eta, float zeta);
    static Vector3 EvalFlat(std::vector<Vector3>& verts, float eta, float zeta);

    static Mesh MakeSphereMesh(float radius, int rings, int sectors);

    static void WriteFlatObj(const Mesh& sphereMesh, const char* filename);
    void WriteNagataObj(Mesh& sphereMesh, const char* filename, int N);
    pair<float, float> MeasureError(const Mesh& Mesh, float radius, int N);
    void RunSimplificationExperiment();
};

