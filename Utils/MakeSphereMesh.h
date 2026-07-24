#pragma once

#include "Vector3.h"
#include "Mesh.h"
#include "Evaluator.h"
using namespace std;

class MakeSphereMesh
{
public:
    
    MakeSphereMesh(float radius, int rings, int sectors,Mesh &sphereMesh)
    {
        if (rings < 2 || sectors < 3) return;
        
        auto idx = [sectors](int r, int j) {
            return 1 + (r - 1) * sectors + ((j % sectors) + sectors) % sectors;
        };
        const int northPole = 0;
        const int southPole = 1 + (rings - 1) * sectors;
        sphereMesh.vertices.push_back({0.0f, 0.0f, radius});
        sphereMesh.normals.push_back({0.0f, 0.0f, 1.0f});
        
        for (int r = 1; r <= rings-1; ++r)
        {
            float theta = pi * static_cast<float>(r) / rings;
            float st = std::sin(theta), ct = std::cos(theta);
            for (int j = 0; j < sectors; ++j)
            {
                float phi = 2.0f * pi * j / sectors;
                float x = radius * st * std::cos(phi);
                float y = radius * st * std::sin(phi);
                float z = radius * ct;

                Vector3 p = {x, y, z};
                sphereMesh.vertices.push_back(p);
                sphereMesh.normals.push_back(p.normalize());
            }
        }
        sphereMesh.vertices.push_back({0.0f, 0.0f, -radius});
        sphereMesh.normals.push_back({0.0f, 0.0f, -1.0f});
        
        for (int j = 0; j < sectors; ++j)
            sphereMesh.triangles.push_back({northPole, idx(1, j), idx(1, j + 1)});
        
        for (int r = 1; r <= rings - 2; ++r)
        {
            for (int j = 0; j < sectors; ++j)
            {
                int a = idx(r,     j);
                int b = idx(r,     j + 1);
                int c = idx(r + 1, j);
                int d = idx(r + 1, j + 1);
                sphereMesh.triangles.push_back({a, c, b});
                sphereMesh.triangles.push_back({b, c, d});
            }
        }
        
        for (int j = 0; j < sectors; ++j)
            sphereMesh.triangles.push_back({southPole, idx(rings - 1, j + 1), idx(rings - 1, j)});
    }

};
