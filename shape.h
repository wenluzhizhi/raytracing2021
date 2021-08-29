#ifndef SHAPE_H
#define SHAPE_H
#include <iostream>
#include <glm/glm.hpp>
using namespace std;
using namespace glm;

// 光线
typedef struct Ray
{
    vec3 startPoint = vec3(0, 0, 0); // 起点
    vec3 direction = vec3(0, 0, 0);  // 方向
} Ray;

typedef struct Material
{
    bool isEmissive = false;        // 是否发光
    vec3 normal = vec3(0, 0, 0);    // 法向量
    vec3 color = vec3(0, 0, 0);     // 颜色
    double specularRate = 0.0f;     // 反射光占比
    double roughness = 1.0f;        // 粗糙程度
    double refractRate = 0.0f;      // 折射光占比
    double refractAngle = 1.0f;     // 折射率
    double refractRoughness = 0.0f; // 折射粗糙度
} Material;

typedef struct HitResult
{
    bool isHit = false;            // 是否命中
    double distance = 0.0f;        // 与交点的距离
    vec3 hitPoint = vec3(0, 0, 0); // 光线命中点
    Material material;             // 命中点的表面材质
} HitResult;

class Shape
{
public:
    Shape() {}
    virtual HitResult intersect(Ray ray)
    {
        return HitResult();
    }
};


// 三角形
class Triangle : public Shape
{
public:
    Triangle(){}
    Triangle(vec3 P1, vec3 P2, vec3 P3, vec3 C) 
    { 
        p1 = P1, p2 = P2, p3 = P3; 
        material.normal = normalize(cross(p2 - p1, p3 - p1)); material.color = C;
    }
    vec3 p1, p2, p3;    // 三顶点
    Material material;  // 材质

    // 与光线求交
    HitResult intersect(Ray ray) 
    { 
        HitResult res;

        vec3 S = ray.startPoint;        // 射线起点
        vec3 d = ray.direction;         // 射线方向
        vec3 N = material.normal;       // 法向量
        if (dot(N, d) > 0.0f) N = -N;   // 获取正确的法向量

        // 如果视线和三角形平行
        if (fabs(dot(N, d)) < 0.00001f) return res;

        // 距离
        float t = (dot(N, p1) - dot(S, N)) / dot(d, N);
        if (t < 0.0005f) return res;    // 如果三角形在相机背面

        // 交点计算
        vec3 P = S + d * t;

        // 判断交点是否在三角形中
        vec3 c1 = cross(p2 - p1, P - p1);
        vec3 c2 = cross(p3 - p2, P - p2);
        vec3 c3 = cross(p1 - p3, P - p3);
        vec3 n = material.normal;   // 需要 "原生法向量" 来判断
        if (dot(c1, n) < 0 || dot(c2, n) < 0 || dot(c3, n) < 0) return res;

        // 装填返回结果
        res.isHit = true;
        res.distance = t;
        res.hitPoint = P;
        res.material = material;
        res.material.normal = N;    // 要返回正确的法向
        return res; 
    };
};


class Sphere : public Shape
{
public:
    Sphere(){}
    Sphere(vec3 o, double r, vec3 c) { O = o; R = r; material.color = c; }
    vec3 O;             // 圆心
    double R;           // 半径
    Material material;  // 材质

    // 与光线求交
    HitResult intersect(Ray ray)
    {
        HitResult res;

        vec3 S = ray.startPoint;        // 射线起点
        vec3 d = ray.direction;         // 射线方向

        float OS = length(O - S);
        float SH = dot(O - S, d);
        float OH = sqrt(pow(OS, 2) - pow(SH, 2));

        if (OH > R) return res; // OH大于半径则不相交

        float PH = sqrt(pow(R, 2) - pow(OH, 2));

        float t1 = length(SH) - PH;
        float t2 = length(SH) + PH;
        float t = (t1 < 0) ? (t2) : (t1);   // 最近距离
        vec3 P = S + t * d;     // 交点

        // 防止自己交自己
        if (fabs(t1) < 0.0005f || fabs(t2) < 0.0005f) return res;

        // 装填返回结果
        res.isHit = true;
        res.distance = t;
        res.hitPoint = P;
        res.material = material;
        res.material.normal = normalize(P - O); // 要返回正确的法向
        return res;
    }
};


#endif      