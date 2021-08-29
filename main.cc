#include <iostream>
#include <vector>
#include <random>
#include <stdlib.h>
#include <glm/glm.hpp> // 数学库支持
#include "svpng.inc"   // png输出 ref: https://github.com/miloyip/svpng
#include <omp.h>       // openmp多线程加速
#include <string.h>
#include <memory.h>

#ifdef __unix
#define fopen_s(pFile, filename, mode) ((*(pFile)) = fopen((filename), (mode))) == NULL
#endif

using namespace glm;
using namespace std;

const vec3 RED(1, 0.5, 0.5);
const vec3 GREEN(0.5, 1, 0.5);
const vec3 BLUE(0.5, 0.5, 1);
const vec3 YELLOW(1.0, 1.0, 0.1);
const vec3 CYAN(0.1, 1.0, 1.0);
const vec3 MAGENTA(1.0, 0.1, 1.0);
const vec3 GRAY(0.5, 0.5, 0.5);
const vec3 WHITE(1, 1, 1);


std::uniform_real_distribution<> dis(0.0, 1.0);
random_device rd;
mt19937 gen(rd());
double randf()
{
  return dis(gen);
}

vec3 randomVec3()
{
  vec3 d;
  do
  {
    d = 2.0f * vec3(randf(), randf(), randf()) - vec3(1, 1, 1);
  } while (dot(d, d) > 1.0);
  return normalize(d);
}

vec3 randomDirection(vec3 n)
{
  // 法向半球
  vec3 d;
  do
  {
    d = randomVec3();
  } while (dot(d, n) < 0.0f);
  return d;
}


typedef struct Ray
{
  vec3 startPoint = vec3(0, 0, 0);
  vec3 direction = vec3(0, 0, 0);
} Ray;

typedef struct Material
{
  bool isEmissive = false;
  vec3 normal = vec3(0, 0, 0);
  vec3 color = vec3(0, 0, 0);
  double specularRate = 0.0;
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
  virtual HitResult intersect(Ray ray) { return HitResult(); }
};

class Sphere : public Shape
{
public:
  Sphere() {}
  Sphere(vec3 o, double r, vec3 c)
  {
    O = o;
    R = r;
    material.color = c;
  }
  vec3 O;            // 圆心
  double R;          // 半径
  Material material; // 材质
  HitResult intersect(Ray ray)
  {
    HitResult res;

    vec3 S = ray.startPoint; // 射线起点
    vec3 d = ray.direction;  // 射线方向

    float OS = length(O - S);
    float SH = dot(O - S, d);
    float OH = sqrt(pow(OS, 2) - pow(SH, 2));

    if (OH > R)
      return res; // OH大于半径则不相交

    float PH = sqrt(pow(R, 2) - pow(OH, 2));

    float t1 = length(SH) - PH;
    float t2 = length(SH) + PH;
    float t = (t1 < 0) ? (t2) : (t1); // 最近距离
    vec3 P = S + t * d;               // 交点

    // 防止自己交自己
    if (fabs(t1) < 0.0005f || fabs(t2) < 0.0005f)
      return res;

    // 装填返回结果
    res.isHit = true;
    res.distance = t;
    res.hitPoint = P;
    res.material = material;
    res.material.normal = normalize(P - O); // 要返回正确的法向
    return res;
  }
};

class Triangle : public Shape
{
public:
  Triangle() {}
  Triangle(vec3 P1, vec3 P2, vec3 P3, vec3 C)
  {
    p1 = P1, p2 = P2, p3 = P3;
    material.normal = normalize(cross(p2 - p1, p3 - p1));
    material.color = C;
  }
  vec3 p1, p2, p3;   // 三顶点
  Material material; // 材质

  // 与光线求交
  // 与光线求交
  HitResult intersect(Ray ray)
  {
    HitResult res;

    vec3 S = ray.startPoint;  // 射线起点
    vec3 d = ray.direction;   // 射线方向
    vec3 N = material.normal; // 法向量
    if (dot(N, d) > 0.0f)
      N = -N; // 获取正确的法向量

    // 如果视线和三角形平行
    if (fabs(dot(N, d)) < 0.00001f)
      return res;

    // 距离
    float t = (dot(N, p1) - dot(S, N)) / dot(d, N);
    if (t < 0.0005f)
      return res; // 如果三角形在相机背面

    // 交点计算
    vec3 P = S + d * t;

    // 判断交点是否在三角形中
    vec3 c1 = cross(p2 - p1, P - p1);
    vec3 c2 = cross(p3 - p2, P - p2);
    vec3 c3 = cross(p1 - p3, P - p3);
    vec3 n = material.normal; // 需要 "原生法向量" 来判断
    if (dot(c1, n) < 0 || dot(c2, n) < 0 || dot(c3, n) < 0)
      return res;

    // 装填返回结果
    res.isHit = true;
    res.distance = t;
    res.hitPoint = P;
    res.material = material;
    res.material.normal = N; // 要返回正确的法向
    return res;
  };
};
const int scale = 1;
const int WIDTH = 256 * scale;
const int HEIGHT = 256 * scale;

// 相机参数
const double SCREEN_Z = 1.1;      // 视平面 z 坐标
const vec3 EYE = vec3(0, 0, 4.0); // 相机位置
HitResult shoot(vector<Shape *> &shapes, Ray ray)
{
  HitResult res, r;
  res.distance = 1145141919.810f; // inf

  // 遍历所有图形，求最近交点
  for (auto &shape : shapes)
  {
    r = shape->intersect(ray);
    if (r.isHit && r.distance < res.distance)
      res = r; // 记录距离最近的求交结果
  }

  return res;
}
vec3 pathTracing(vector<Shape *> &shapes, Ray ray, int depth)
{
  if (depth > 8)
    return vec3(0);
  HitResult res = shoot(shapes, ray);

  if (!res.isHit)
    return vec3(0);

  if (res.material.isEmissive)
    return res.material.color;

  double r = randf();
  float p = 0.8;
  if (r > p)
    return vec3(0);

  Ray randomRay;
  randomRay.startPoint = res.hitPoint;
  randomRay.direction = randomDirection(res.material.normal);

  vec3 color = vec3(0);
  float cosine = fabs(dot(-ray.direction, res.material.normal));
  
  r = randf();
  if (r < res.material.specularRate)  // 镜面反射
    {
        randomRay.direction = normalize(reflect(ray.direction, res.material.normal));
        color = pathTracing(shapes, randomRay, depth + 1) * cosine;
    }
    else    // 漫反射
    {
        vec3 srcColor = res.material.color;
        vec3 ptColor = pathTracing(shapes, randomRay, depth+1) * cosine;
        color = ptColor * srcColor;    // 和原颜色混合
    }
    return color/p;

}

void imshow(double *SRC)
{

  unsigned char *image = new unsigned char[WIDTH * HEIGHT * 3]; // 图像buffer
  unsigned char *p = image;
  double *S = SRC; // 源数据

  FILE *fp;
  fopen_s(&fp, "image.png", "wb");

  for (int i = 0; i < HEIGHT; i++)
  {
    for (int j = 0; j < WIDTH; j++)
    {
      *p++ = (unsigned char)clamp((*S++) * 255, 0.0, 255.0); // R 通道
      *p++ = (unsigned char)clamp((*S++) * 255, 0.0, 255.0); // G 通道
      *p++ = (unsigned char)clamp((*S++) * 255, 0.0, 255.0); // B 通道
    }
  }

  svpng(fp, WIDTH, HEIGHT, image, 0);
}


void testLoadPic()
{

  const int SAMPLE = 2000;
  const double BRIGHTNESS = (2.0f * 3.1415926f) * (1.0f / double(SAMPLE));
  vector<Shape *> shapes; // 几何物体的集合
  // 三角形
  shapes.push_back(new Triangle(vec3(-0.5, -0.5, -0.5), vec3(0.5, -0.5, -0.5), vec3(0, -0.5, 0.5), CYAN));
  // 底部平面
  shapes.push_back(new Triangle(vec3(10, -1, 10), vec3(-10, -1, -10), vec3(-10, -1, 10), WHITE));
  shapes.push_back(new Triangle(vec3(10, -1, 10), vec3(10, -1, -10), vec3(-10, -1, -10), WHITE));
  // 光源
  Triangle l1 = Triangle(vec3(0.6, 0.99, 0.4), vec3(-0.2, 0.99, -0.4), vec3(-0.2, 0.99, 0.4), WHITE);
  Triangle l2 = Triangle(vec3(0.6, 0.99, 0.4), vec3(0.6, 0.99, -0.4), vec3(-0.2, 0.99, -0.4), WHITE);
  l1.material.isEmissive = true;
  l2.material.isEmissive = true;
  shapes.push_back(&l1);
  shapes.push_back(&l2);

  Sphere sp1 = Sphere(vec3(-0.6, -0.8, 0.6), 0.2, WHITE);
  sp1.material.specularRate = 0.5;

  shapes.push_back(&sp1);
  shapes.push_back(new Sphere(vec3(-0.1, -0.7, 0.2), 0.3, WHITE));
  shapes.push_back(new Sphere(vec3(0.5, -0.6, -0.5), 0.4, WHITE));

  double *image = new double[WIDTH * HEIGHT * 3];
  memset(image, 0.0, sizeof(double) * WIDTH * HEIGHT * 3);

  // omp_set_num_threads(50); // 线程个数
  //#pragma omp parallel for
  for (int k = 0; k < SAMPLE; k++)
  {
    std::cout << k << std::endl;
    double *p = image;
    for (int i = 0; i < HEIGHT; i++)
    {
      for (int j = 0; j < WIDTH; j++)
      {
        // 像素坐标转投影平面坐标
        double x = 2.0 * double(j) / double(WIDTH) - 1.0;
        double y = 2.0 * double(HEIGHT - i) / double(HEIGHT) - 1.0;

        vec3 coord = vec3(x, y, SCREEN_Z);       // 计算投影平面坐标
        vec3 direction = normalize(coord - EYE); // 计算光线投射方向

        // 生成光线
        Ray ray;
        ray.startPoint = coord;
        ray.direction = direction;

        // 与场景的交点
        HitResult res = shoot(shapes, ray);
        vec3 color = vec3(0, 0, 0);

        if (res.isHit)
        {
          // 命中光源直接返回光源颜色
          if (res.material.isEmissive)
          {
            color = res.material.color;
          }
          // 命中实体则选择一个随机方向重新发射光线并且进行路径追踪
          else
          {
            // 根据交点处法向量生成交点处反射的随机半球向量
            Ray randomRay;
            randomRay.startPoint = res.hitPoint;
            randomRay.direction = randomDirection(res.material.normal);

            // 颜色积累
            vec3 srcColor = res.material.color;
            vec3 ptColor = pathTracing(shapes, randomRay, 0);
            color = ptColor * srcColor; // 和原颜色混合
            color *= BRIGHTNESS;
          }
        }

        *p += color.x;
        p++; // R 通道
        *p += color.y;
        p++; // G 通道
        *p += color.z;
        p++; // B 通道
      }
    }
  }

  imshow(image);
}

int main()
{

  testLoadPic();
  return 0;
}