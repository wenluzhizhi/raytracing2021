#include <iostream>
#include <vector>
#include <random>
#include <stdlib.h>
#include <glm/glm.hpp> // 数学库支持
#include "svpng.inc"   // png输出 ref: https://github.com/miloyip/svpng
#include <omp.h>       // openmp多线程加速
#include "shape.h"
#include <cstring>

using namespace glm;
using namespace std;

#ifdef __unix
#define fopen_s(pFile, filename, mode) ((*(pFile)) = fopen((filename), (mode))) == NULL
#endif

std::uniform_real_distribution<> dis(0.0, 1.0);
random_device rd;
mt19937 gen(rd());
double randf()
{
  return dis(gen);
}

const int WIDTH = 256;
const int HEIGHT = 255;
const int SAMPLE = 1000;
const double BRIGHTNESS = (2.0f * 3.1415926f) * (1.0f / double(SAMPLE));
const double SCREEN_Z = 1.1;
const vec3 EYE = vec3(0, 0, 4.0); // 相机位置

const vec3 RED(1, 0.5, 0.5);
const vec3 GREEN(0.5, 1, 0.5);
const vec3 BLUE(0.5, 0.5, 1);
const vec3 YELLOW(1.0, 1.0, 0.1);
const vec3 CYAN(0.1, 1.0, 1.0);
const vec3 MAGENTA(1.0, 0.1, 1.0);
const vec3 GRAY(0.5, 0.5, 0.5);
const vec3 WHITE(1, 1, 1);

void imgShow(double *SRC)
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
      *p++ = (unsigned char)clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0); // R 通道
      *p++ = (unsigned char)clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0); // G 通道
      *p++ = (unsigned char)clamp(pow(*S++, 1.0f / 2.2f) * 255, 0.0, 255.0); // B 通道
    }
  }
  svpng(fp, WIDTH, HEIGHT, image, 0);
};

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
  return normalize(randomVec3() + n);
}

HitResult shoot(vector<Shape *> &shapes, Ray ray)
{
  HitResult res, r;
  res.distance = 1145141919.810f;
  for (auto &shape : shapes)
  {
    r = shape->intersect(ray);
    if (r.isHit && r.distance < res.distance)
    {
      res = r;
    }
  }
  return res;
}

vec3 pathTracing(vector<Shape *> shapes, Ray ray, int depth)
{
  if (depth > 8)
    return vec3(0);
  HitResult res = shoot(shapes, ray);
  if (!res.isHit)
    return vec3(0); // 未命中
  if (res.material.isEmissive)
    return res.material.color;

  double r = randf();
  float P = 0.8;
  if (r > P)
    return vec3(0);
  Ray randomRay;
  randomRay.startPoint = res.hitPoint;
  randomRay.direction = randomDirection(res.material.normal);

  vec3 color = vec3(0);
  float cosine = fabs(dot(-ray.direction, res.material.normal));

  r = randf();
  if (r < res.material.specularRate) // 镜面反射
  {
    vec3 ref = normalize(reflect(ray.direction, res.material.normal));
    randomRay.direction = mix(ref, randomRay.direction, res.material.roughness);
    color = pathTracing(shapes, randomRay, depth + 1) * cosine;
  }
  else if (res.material.specularRate <= r && r <= res.material.refractRate) // 折射
  {
    vec3 ref = normalize(refract(ray.direction, res.material.normal, float(res.material.refractAngle)));
    randomRay.direction = mix(ref, -randomRay.direction, res.material.refractRoughness);
    color = pathTracing(shapes, randomRay, depth + 1) * cosine;
  }
  else // 漫反射
  {
    vec3 srcColor = res.material.color;
    vec3 ptColor = pathTracing(shapes, randomRay, depth + 1) * cosine;
    color = ptColor * srcColor; // 和原颜色混合
  }

  return color / P;
}

int main()
{
  cout << "main function running" << endl;
  vector<Shape *> shapes;

  // 背景盒子
  // bottom
  shapes.push_back(new Triangle(vec3(1, -1, 1), vec3(-1, -1, -1), vec3(-1, -1, 1), RED));
  shapes.push_back(new Triangle(vec3(1, -1, 1), vec3(1, -1, -1), vec3(-1, -1, -1), RED));
  // top
  Triangle t1 = Triangle(vec3(1, 1, 1), vec3(-1, 1, 1), vec3(-1, 1, -1), WHITE);
  t1.material.isEmissive = true;
  shapes.push_back(&t1);
  shapes.push_back(new Triangle(vec3(1, 1, 1), vec3(-1, 1, -1), vec3(1, 1, -1), WHITE));
  // back
  shapes.push_back(new Triangle(vec3(1, -1, -1), vec3(-1, 1, -1), vec3(-1, -1, -1), CYAN));
  shapes.push_back(new Triangle(vec3(1, -1, -1), vec3(1, 1, -1), vec3(-1, 1, -1), CYAN));
  // left
  shapes.push_back(new Triangle(vec3(-1, -1, -1), vec3(-1, 1, 1), vec3(-1, -1, 1), BLUE));
  shapes.push_back(new Triangle(vec3(-1, -1, -1), vec3(-1, 1, -1), vec3(-1, 1, 1), BLUE));
  // right
  shapes.push_back(new Triangle(vec3(1, 1, 1), vec3(1, -1, -1), vec3(1, -1, 1), RED));
  shapes.push_back(new Triangle(vec3(1, -1, -1), vec3(1, 1, 1), vec3(1, 1, -1), RED));

  Sphere s1 = Sphere(vec3(-0.65, -0.7, 0.0), 0.3, GREEN);
  Sphere s2 = Sphere(vec3(0.0, -0.3, 0.0), 0.4, WHITE);
  Sphere s3 = Sphere(vec3(0.65, 0.1, 0.0), 0.3, BLUE);
  s1.material.specularRate = 0.3;
  s1.material.roughness = 0.1;

  s2.material.specularRate = 0.3;
  s2.material.refractRate = 0.95;
  s2.material.refractAngle = 0.1;
  //s2.material.refractRoughness = 0.05;

  s3.material.specularRate = 0.3;

  shapes.push_back(&s1);
  shapes.push_back(&s2);
  shapes.push_back(&s3);

  double *image = new double[WIDTH * HEIGHT * 3];
  memset(image, 0.0, sizeof(double) * WIDTH * HEIGHT * 3);

  cout << randf() << endl;
  cout << randf() << endl;

  for (int k = 0; k < SAMPLE; k++)
  {
    cout << k << endl;
    double *p = image;
    for (int i = 0; i < HEIGHT; i++)
    {
      for (int j = 0; j < WIDTH; j++)
      {
        double x = 2.0 * double(j) / double(WIDTH) - 1.0;
        double y = 2.0 * double(HEIGHT - i) / double(HEIGHT) - 1.0;

        x += (randf() - 0.5f) / double(WIDTH);
        y += (randf() - 0.5f) / double(HEIGHT);

        vec3 coord = vec3(x, y, SCREEN_Z);       // 计算投影平面坐标
        vec3 direction = normalize(coord - EYE); // 计算光线投射方向

        // 生成光线
        Ray ray;
        ray.startPoint = coord;
        ray.direction = direction;

        HitResult res = shoot(shapes, ray);
        vec3 color = vec3(0, 0, 0);

        if (res.isHit)
        {
          if (res.material.isEmissive)
          {
            color = res.material.color;
          }
          else
          {
            Ray randomRay;
            randomRay.startPoint = res.hitPoint;
            randomRay.direction = randomDirection(res.material.normal);
            double r = randf();
            if (r < res.material.specularRate) // 镜面反射
            {
              vec3 ref = normalize(reflect(ray.direction, res.material.normal));
              randomRay.direction = mix(ref, randomRay.direction, res.material.roughness);
              color = pathTracing(shapes, randomRay, 0);
            }
            else if (res.material.specularRate <= r && r <= res.material.refractRate) // 折射
            {
              vec3 ref = normalize(refract(ray.direction, res.material.normal, float(res.material.refractAngle)));
              randomRay.direction = mix(ref, -randomRay.direction, res.material.refractRoughness);
              color = pathTracing(shapes, randomRay, 0);
            }
            else // 漫反射
            {
              vec3 srcColor = res.material.color;
              vec3 ptColor = pathTracing(shapes, randomRay, 0);
              color = ptColor * srcColor; // 和原颜色混合
            }

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
  imgShow(image);
  return 0;
}