#pragma once
#include <cmath>

class Vector3
{
public:
    double x, y, z;
    Vector3(double x,double y,double z):x(x),y(y),z(z){};
    Vector3():x(0),y(0),z(0){};
    Vector3 operator+(Vector3 other) const
    {
        return Vector3(x+other.x,y+other.y, z+other.z);
    }
    Vector3 operator-(Vector3 other) const
    {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }
    inline Vector3 operator*(float other) const
    {
        return Vector3(x*other , y*other , z*other);
    }
    
    Vector3 operator/(float other) const
    {
        return Vector3(x/other , y/other , z/other);
    }
    double length() const
    {
        return (sqrt(x*x + y*y + z*z));
    }
    
    double dot(Vector3 other) const
    {
        return (x*other.x + y*other.y + z*other.z);
    }
    
    Vector3 operator-() const
    {
        return (Vector3(-x, -y, -z));
    }
    
    Vector3 cross(Vector3 other) const
    {
        return Vector3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }
    
    Vector3 normalize() const
    {
        float len = length();
        if (len > 0)
            return Vector3(x / len, y / len, z / len);
        return Vector3(0, 0, 0);
    }
    friend inline Vector3 operator*(float s, Vector3 v) { return v * s; }
};
