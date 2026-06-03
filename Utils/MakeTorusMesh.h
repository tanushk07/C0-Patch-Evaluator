#pragma once
#include "Evaluator.h"
#include "MakeSphereMesh.h"

class MakeTorusMesh
{
public:
    MakeTorusMesh(float CenterLineRadius, float TubeRadius, int majorSegments, int minorSegments, Mesh &TorusMesh)
    {
		if (majorSegments < 3 || minorSegments < 3) return;

        for (int i = 0; i <= majorSegments; ++i)
        {
            float theta = 2.0f * pi * static_cast<float>(i) / majorSegments;
            
            float ct = std::cos (theta);
            float st = std::sin (theta);
            for (int j = 0; j <= minorSegments; ++j)
            {
                float phi = 2.0f * pi * static_cast<float>(j) / minorSegments;
                float cp = std::cos (phi);
                float sp = std::sin (phi);
                float rho = CenterLineRadius + TubeRadius*cp;
                Vector3 p = {rho* ct,rho*st,TubeRadius*sp};
                Vector3 n = {cp*ct,cp*st,sp};
                TorusMesh.vertices.push_back(p);
                TorusMesh.normals.push_back(n);
            }
        }
        
        for (int j = 0; j < majorSegments; ++j)
        {
            for (int i = 0; i < minorSegments; ++i)
            {
                int row1 = j* (minorSegments+1);
                int row2 = (j+1)* (minorSegments+1);
                
                int i0 = row1+i;
                int i1 = row1+i+1;
                int i2 = row2+i;
                int i3 = row2+i+1;
                
                TorusMesh.triangles.push_back({i0,i2,i1});
                TorusMesh.triangles.push_back({i1,i2,i3});
            }
        }
    }
};
