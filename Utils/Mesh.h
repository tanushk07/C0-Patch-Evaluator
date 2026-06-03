#pragma once
#include <vector>
#include "Vector3.h"


struct Tri{int a,b,c;};
struct Mesh
{
    std::vector<Vector3> vertices, normals;
    std::vector<Tri> triangles;
};
