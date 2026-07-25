#pragma once
#include "Evaluator.h"
#include "Mesh.h"

class MakeTorusMesh
{
public:
    static int id(int i, int j,int Nu,int Nv) {
        return ((i + Nu) % Nu) * Nv + ((j + Nv) % Nv);
    }
    MakeTorusMesh(float CenterLineRadius, float TubeRadius, int majorSegments, int minorSegments, Mesh &TorusMesh)
    {
		if (majorSegments < 3 || minorSegments < 3) return;
        
        const long long n = majorSegments * minorSegments;
        TorusMesh.vertices.resize(n);
        TorusMesh.normals.resize(n);
        for (int i = 0; i < majorSegments; ++i)
        {
            float theta = 2.0f * pi * static_cast<float>(i) / majorSegments;
            
            float ct = std::cos (theta);
            float st = std::sin (theta);
            for (int j = 0; j < minorSegments; ++j)
            {
                float phi = 2.0f * pi * static_cast<float>(j) / minorSegments;
                float cp = std::cos (phi);
                float sp = std::sin (phi);
                float rho = CenterLineRadius + TubeRadius*cp;
                Vector3 p = {rho* ct,rho*st,TubeRadius*sp};
                Vector3 n = {cp*ct,cp*st,sp};
                int ids = id(i,j,majorSegments,minorSegments);
                TorusMesh.vertices[ids] = p;
                TorusMesh.normals[ids] = n;
            }
        }
        
        for (int j = 0; j < majorSegments; ++j)
        {
            for (int i = 0; i < minorSegments; ++i)
            {                
                int i0 = id(j,i,majorSegments,minorSegments);
                int i1 = id(j+1,i,majorSegments,minorSegments);
                int i2 = id(j,i+1,majorSegments,minorSegments);
                int i3 = id(j+1,i+1,majorSegments,minorSegments);
                
                TorusMesh.triangles.push_back({i0,i1,i2});
                TorusMesh.triangles.push_back({i1,i3,i2});
            }
        }
    }
};
