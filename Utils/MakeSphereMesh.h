#pragma once

#include "Vector3.h"
#include "Mesh.h"

using namespace std;

class MakeSphereMesh
{
public:
    
    MakeSphereMesh(float radius, int rings, int sectors,Mesh &sphereMesh)
    {
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
    }

};
