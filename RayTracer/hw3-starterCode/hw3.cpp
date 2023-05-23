/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Angelina Liang
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#include <cmath>
#include <limits>
#include <iostream>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

#define PI 3.1415926535

// intersection types
#define NONE_INRCT -1
#define SPHERE_INRCT 0
#define TRI_INRCT 1

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

struct Vector3
{
    double x, y, z;
    Vector3(double mx = 0.0, double my = 0.0, double mz = 0.0) {
        x = mx;
        y = my;
        z = mz;
    }
    Vector3 addToSelf(Vector3 other) { // this + other
        return Vector3(x + other.x, y + other.y, z + other.z);
    }
    Vector3 minusToSelf(Vector3 other) { // this - other
        return Vector3(x - other.x, y - other.y, z - other.z);
    }
    Vector3 Multiply(double scaler) {
        return Vector3(x * scaler, y * scaler, z * scaler);
    }
    Vector3 Negate() {
        return Vector3(-x, -y, -z);
    }
    double Magnitude() {
        return sqrt(x * x + y * y + z * z);
    }
    Vector3 Normalize() {
        double temp = this->Magnitude();
        return this->Multiply(1.0 / temp);
    }
    double dotToSelf(Vector3 other) { // this dot other
        return x * other.x + y * other.y + z * other.z;
    }
    Vector3 crossToSelf(Vector3 other) { // this cross other
        double xx = y * other.z - z * other.y;
        double yy = z * other.x - x * other.z;
        double zz = x * other.y - y * other.x;
        return Vector3(xx, yy, zz);

    }
    void Set(Vector3 a) { // this = a
        x = a.x;
        y = a.y;
        z = a.z;
    }
    void Clamp() {
        if (x > 1.0) x = 1.0;
        if (y > 1.0) y = 1.0;
        if (z > 1.0) z = 1.0;
    }
};

Vector3 Add(Vector3 a, Vector3 b) { // a+b
    return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vector3 Minus(Vector3 a, Vector3 b) { // a-b
    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}
double Dot(Vector3 a, Vector3 b) { // a dot b
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vector3 Cross(Vector3 a, Vector3 b) { // a cross b
    double xx = a.y * b.z - a.z * b.y;
    double yy = a.z * b.x - a.x * b.z;
    double zz = a.x * b.y - a.y * b.x;
    return Vector3(xx, yy, zz);
}


struct Ray
{
    Vector3 start;
    Vector3 dir;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

/*
* Send rays from camera location
*/
Ray GenerateRay(double x, double y) {
    Ray res;
    double aspect_ratio = (static_cast<double>(WIDTH)) / (static_cast<double>(HEIGHT));
    double angle = fov / 2.0 * PI / 180.0;

    double x_min = -aspect_ratio * tan(angle);
    double x_max = aspect_ratio * tan(angle);
    double y_min = -tan(angle);
    double y_max = tan(angle);

    // set ray data
    res.start.x = 0.0;
    res.start.y = 0.0;
    res.start.z = 0.0;

    res.dir.x = x_min + (x / (static_cast<double>(WIDTH))) * 2.0 * x_max;
    res.dir.y = y_min + (y / (static_cast<double>(HEIGHT))) * 2.0 * y_max;
    res.dir.z = -1.0;
    res.dir = res.dir.Normalize();
    //std::cout << res.dir.x << " " << res.dir.y << " " << res.dir.z << std::endl;

    return res;
}

/*
* Intersection check - Sphere
*/
bool isSphereIntersect(Ray ray, double& min_t, int& sph_idx) {
    // loop over all spheres, return true if intersects, 
    // saves the min t among all intersections 
    // and which sphere it intersects with
    sph_idx = -1;

    for (int i = 0; i < num_spheres; ++i) {
        Vector3 sphere_pos(
            spheres[i].position[0], 
            spheres[i].position[1], 
            spheres[i].position[2]);

        double b = 2.0 * (ray.dir.x * (ray.start.x - sphere_pos.x)
            + ray.dir.y * (ray.start.y - sphere_pos.y) 
            + ray.dir.z * (ray.start.z - sphere_pos.z));

        double c = (ray.start.x - sphere_pos.x) * (ray.start.x - sphere_pos.x)
            + (ray.start.y - sphere_pos.y) * (ray.start.y - sphere_pos.y)
            + (ray.start.z - sphere_pos.z) * (ray.start.z - sphere_pos.z)
            - spheres[i].radius * spheres[i].radius;

        double checker = b * b - 4.0 * c;

        // abort if negate
        if (checker >= 0.0){
            // get t0, t1
            double t0 = (-b - sqrt(checker)) / 2.0;
            double t1 = (-b + sqrt(checker)) / 2.0;

            // check if t0,1 > 0, and get min(t0,t1)
            double min = DBL_MAX;
            if (t0 > 0.0 && t1 > 0.0) {
                min = min(t0, t1);
            }
            else if (t0 > 0.0 && t1 < 0.0) {
                min = t0;
            }
            else if (t1 > 0.0 && t0 < 0.0) {
                min = t1;
            }
            else
            {
                // now, both t0 and t1 are negative
                continue;
            }

            // update total min
            if (min < min_t) {
                min_t = min;
                sph_idx = i;
            }
        }
    }

    return (sph_idx != -1);
}

/*
* Intersection check - Triangle
*/
bool isIntersectInTriangle(Ray ray, Vector3 normal, double t, Vector3 a, Vector3 b, Vector3 c, 
    double& alpha, double& beta) {
    // check if the ray is intersecting in the range of triangle
    // record alpha and beta for phong calculation

    // intersection point
    Vector3 p0 = Add(ray.start, ray.dir.Multiply(t));

    //// tri area
    Vector3 ab = Minus(b, a); Vector3 ac = Minus(c, a);
    Vector3 ABC = Cross(ab, ac);
    double area = ABC.Magnitude() / 2.0;

    // edge checks
    Vector3 ap0 = Minus(p0, a);
    Vector3 ab_ap0 = Cross(ab, ap0);
    double ab_checker = Dot(ab_ap0, normal);
    if (ab_checker < 0.0) return false;

    Vector3 bc = Minus(c, b);
    Vector3 bp0 = Minus(p0, b);
    Vector3 bc_bp0 = Cross(bc, bp0);
    double bc_checker = Dot(bc_bp0, normal);
    if (bc_checker / area < 0.0) return false;

    Vector3 ca = Minus(a, c);
    Vector3 cp0 = Minus(p0, c);
    Vector3 ca_cp0 = Cross(ca, cp0);
    double ca_checker = Dot(ca_cp0, normal);
    if (ca_checker / area < 0.0) return false;

    // record alpha, beta
    double bcp0 = Cross(bc, bp0).Magnitude() / 2.0;
    double cap0 = Cross(ca, cp0).Magnitude() / 2.0;
    alpha = bcp0 / area;
    beta = cap0 / area;

    //std::cout << alpha << " " << beta << " " << (1-alpha-beta) << std::endl;
    return true;
}

bool isTriangleIntersect(Ray ray, double& min_t, int& tri_idx, double& alpha, double& beta) {
    // for all triangles, check intersection
    // records the min t and corresponding idx of plane
    // records alpha and beta that is used for phong shading calculation
    tri_idx = -1;

    // for all triangles check for plane intersection
    for (int i = 0; i < num_triangles; ++i) {
        Vector3 a(
            triangles[i].v[0].position[0],
            triangles[i].v[0].position[1],
            triangles[i].v[0].position[2]);
        Vector3 b(
            triangles[i].v[1].position[0],
            triangles[i].v[1].position[1],
            triangles[i].v[1].position[2]);
        Vector3 c(
            triangles[i].v[2].position[0],
            triangles[i].v[2].position[1],
            triangles[i].v[2].position[2]);

        // get normal of plane
        Vector3 ab = Minus(b,a);
        Vector3 ac = Minus(c,a);
        Vector3 normal = Cross(ab,ac);
        normal = normal.Normalize();

        // check n dot d
        double nd = Dot(normal, ray.dir);
        if (std::abs(nd) > 1e-8) {
            double t = (Dot(normal, Minus(a, ray.start))) / nd;
            if (t > 0.0 && t < min_t
                && isIntersectInTriangle(ray, normal, t, a, b, c, alpha, beta)) { // validity check
                min_t = t;
                tri_idx = i;
            }
        }

    }

    return (tri_idx != -1);
}

/*
*  Generates per pixel color
*/
int GetClosestIntersectionType(double t_sph, double t_tri, int idx_sph, int idx_tri, 
    double& t, int& idx) {
    // return intersection type for further decision on shading method
    // decide on final t and idx values
    if (idx_sph == -1 && idx_tri == -1) {
        return NONE_INRCT;
    }
    else if (idx_sph == -1 && idx_tri != -1) { // only intersects with tri
        t = t_tri;
        idx = idx_tri;
        return TRI_INRCT;
    }
    else if (idx_sph != -1 && idx_tri == -1) { // only intersects with sphere
        t = t_sph;
        idx = idx_sph;
        return SPHERE_INRCT;
    }
    else {
        if (t_tri > t_sph) { // sphere is closer
            t = t_sph;
            idx = idx_sph;
            return SPHERE_INRCT;
        }
        if (t_tri < t_sph) { // tri is closer
            t = t_tri;
            idx = idx_tri;
            return TRI_INRCT;
        }
    }
    return NONE_INRCT;
}

bool isBlockedBySphere(Ray ray, Vector3 light) {
    // test for shadow ray sphere intersection
    for (int i = 0; i < num_spheres; ++i) {
        Vector3 sphere_pos(
            spheres[i].position[0],
            spheres[i].position[1],
            spheres[i].position[2]);

        double b = 2.0 * (ray.dir.x * (ray.start.x - sphere_pos.x)
            + ray.dir.y * (ray.start.y - sphere_pos.y)
            + ray.dir.z * (ray.start.z - sphere_pos.z));

        double c = (ray.start.x - sphere_pos.x) * (ray.start.x - sphere_pos.x)
            + (ray.start.y - sphere_pos.y) * (ray.start.y - sphere_pos.y)
            + (ray.start.z - sphere_pos.z) * (ray.start.z - sphere_pos.z)
            - spheres[i].radius * spheres[i].radius;

        double checker = b * b - 4.0 * c;

        if (checker >= 1e-8) {
            double t0 = (-b - sqrt(checker)) / 2.0;
            double t1 = (-b + sqrt(checker)) / 2.0;

            t0 = t0 < 0 ? 0 : t0;
            t1 = t1 < 0 ? 0 : t1;
            double t = t0 < t1 ? t0 : t1;

            double dis_to_light = Minus(light, ray.start).Magnitude();

            if (t > 1e-8 && t < dis_to_light) return true;

        }
    }
    return false;
}

bool isBlockedByTri(Ray ray, Vector3 light) {
    // test for shadow ray triangle intersection
    for (int i = 0; i < num_triangles; ++i) {
        Vector3 a(
            triangles[i].v[0].position[0],
            triangles[i].v[0].position[1],
            triangles[i].v[0].position[2]);
        Vector3 b(
            triangles[i].v[1].position[0],
            triangles[i].v[1].position[1],
            triangles[i].v[1].position[2]);
        Vector3 c(
            triangles[i].v[2].position[0],
            triangles[i].v[2].position[1],
            triangles[i].v[2].position[2]);
        // get normal of plane
        Vector3 ab = Minus(b, a);
        Vector3 ac = Minus(c, a);
        Vector3 normal = Cross(ab, ac);
        normal = normal.Normalize();

        double nd = Dot(normal, ray.dir);
        if (nd > 1e-8) {
            double t = (Dot(normal, Minus(a, ray.start))) / nd;
            double checker = ray.dir.Multiply(t).Magnitude();
            double dis_to_light = Minus(light, ray.start).Magnitude();
            double dummy1, dummy2;
            if (t >= 1e-8
                && isIntersectInTriangle(ray, normal, t, a, b, c, dummy1, dummy2)
                && checker < dis_to_light) {
                return true;
            }
        }
    }
    return false;
}

Vector3 Reflect(Vector3 l, Vector3 n) { // reflection helper
    // reflect l about n
    // 2*(l dot n)*n-l;
    double L_dot_N = Dot(l, n);
    Vector3 right = n.Multiply(L_dot_N);
    right = right.Multiply(2.0); // 2*(l dot n)*n
    return Minus(right, l).Normalize();
}

Vector3 AddAmbient(Vector3 color) { // add ambient light to color
    Vector3 res;
    // get ambient light
    Vector3 amb(
        ambient_light[0],
        ambient_light[1],
        ambient_light[2]);

    // add ambient
    res = color.addToSelf(amb);
    res.Clamp();

    return res;
}

Vector3 GenerateSpherePhongColor(Ray ray, Ray shadow, int light_idx, double t, int sph_idx) {
    // I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
    // generate light color for one light

    Vector3 res;

    // get shpere info
    Vector3 sphere_pos(
        spheres[sph_idx].position[0],
        spheres[sph_idx].position[1],
        spheres[sph_idx].position[2]);

    // get L, N
    Vector3 light(
        lights[light_idx].position[0],
        lights[light_idx].position[1],
        lights[light_idx].position[2]);
    Vector3 intersection = ray.dir.Multiply(t);
    Vector3 L = Minus(light, intersection);
    L = L.Normalize();
    Vector3 N = Minus(intersection, sphere_pos);
    N = N.Normalize();

    double L_N = Dot(L, N);
    L_N = L_N > 0.0 ? L_N : 0.0;


    // get diffuse kd
    Vector3 sphere_diffuse(
        spheres[sph_idx].color_diffuse[0],
        spheres[sph_idx].color_diffuse[1],
        spheres[sph_idx].color_diffuse[2]);

    // kd * (L dot N)
    sphere_diffuse = sphere_diffuse.Multiply(L_N);

    res = res.addToSelf(sphere_diffuse);

    // get R, V
    Vector3 V = intersection.Negate().Normalize();
    Vector3 R = Reflect(L, N);
    R = R.Normalize();

    double R_V = Dot(R, V);
    R_V = R_V > 0.0 ? R_V : 0.0;

    // get shine
    double sh = spheres[sph_idx].shininess;

    // (R dot V) ^ sh
    double R_V_sh = pow(R_V, sh);

    // get speculer ks
    Vector3 sphere_spec(
        spheres[sph_idx].color_specular[0],
        spheres[sph_idx].color_specular[1],
        spheres[sph_idx].color_specular[2]);

    // ks * (R dot V) ^ sh
    sphere_spec = sphere_spec.Multiply(R_V_sh);

    res = res.addToSelf(sphere_spec);

    // lightColor* (kd * (L dot N) + ks * (R dot V) ^ sh)
    // for all channels
    res.x *= lights[light_idx].color[0];
    res.y *= lights[light_idx].color[1];
    res.z *= lights[light_idx].color[2];

    return res;
}

Vector3 GenerateTriPhongColor(Ray ray, Ray shadow, int light_idx, double t, int tri_idx,
    double alpha, double beta) {
    // I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
    // ganrentees alpha != -1 && beta != -1
    Vector3 res;

    // get L
    Vector3 light(
        lights[light_idx].position[0],
        lights[light_idx].position[1],
        lights[light_idx].position[2]);
    Vector3 intersection = ray.dir.Multiply(t);
    Vector3 L = Minus(light, intersection);
    L = L.Normalize();

    // get N
    double gamma = 1.0 - alpha - beta;
    Vector3 na(
        triangles[tri_idx].v[0].normal[0],
        triangles[tri_idx].v[0].normal[1],
        triangles[tri_idx].v[0].normal[2]);
    Vector3 nb(
        triangles[tri_idx].v[1].normal[0],
        triangles[tri_idx].v[1].normal[1],
        triangles[tri_idx].v[1].normal[2]);
    Vector3 nc(
        triangles[tri_idx].v[2].normal[0],
        triangles[tri_idx].v[2].normal[1],
        triangles[tri_idx].v[2].normal[2]);

    Vector3 N = Add(Add(na.Multiply(alpha),nb.Multiply(beta)), nc.Multiply(gamma));
    N = N.Normalize();

    double L_N = Dot(L, N);
    L_N = L_N > 0.0 ? L_N : 0.0;
    
    // get diffuse
    Vector3 da(
        triangles[tri_idx].v[0].color_diffuse[0],
        triangles[tri_idx].v[0].color_diffuse[1],
        triangles[tri_idx].v[0].color_diffuse[2]);
    Vector3 db(
        triangles[tri_idx].v[1].color_diffuse[0],
        triangles[tri_idx].v[1].color_diffuse[1],
        triangles[tri_idx].v[1].color_diffuse[2]);
    Vector3 dc(
        triangles[tri_idx].v[2].color_diffuse[0],
        triangles[tri_idx].v[2].color_diffuse[1],
        triangles[tri_idx].v[2].color_diffuse[2]);
    // kd = alpha*da + beta*db + gamma*dc
    Vector3 kd = Add(Add(da.Multiply(alpha), db.Multiply(beta)), dc.Multiply(gamma));

    // kd * (L dot N)
    kd = kd.Multiply(L_N);

    res = res.addToSelf(kd);

    // get R, V
    Vector3 V = ray.dir.Negate();
    Vector3 R = Reflect(L, N);
    R = R.Normalize();

    double R_V = Dot(R, V);
    R_V = R_V > 0.0 ? R_V : 0.0;

    // get shine
    double sh = triangles[tri_idx].v[0].shininess * alpha
        + triangles[tri_idx].v[1].shininess * beta 
        + triangles[tri_idx].v[2].shininess * gamma;

    // (R dot V) ^ sh
    double R_V_sh = pow(R_V, sh);

    // get ks
    Vector3 sa(
        triangles[tri_idx].v[0].color_specular[0],
        triangles[tri_idx].v[0].color_specular[1],
        triangles[tri_idx].v[0].color_specular[2]);
    Vector3 sb(
        triangles[tri_idx].v[1].color_specular[0],
        triangles[tri_idx].v[1].color_specular[1],
        triangles[tri_idx].v[1].color_specular[2]);
    Vector3 sc(
        triangles[tri_idx].v[2].color_specular[0],
        triangles[tri_idx].v[2].color_specular[1],
        triangles[tri_idx].v[2].color_specular[2]);

    // ks = alpha*sa + beta*sb + gamma*sc
    Vector3 ks = Add(Add(sa.Multiply(alpha), sb.Multiply(beta)), sc.Multiply(gamma));

    // ks * (R dot V) ^ sh
    ks = ks.Multiply(R_V_sh);

    res = res.addToSelf(ks);

    // lightColor* (kd * (L dot N) + ks * (R dot V) ^ sh)
    // for all channels
    res.x *= lights[light_idx].color[0];
    res.y *= lights[light_idx].color[1];
    res.z *= lights[light_idx].color[2];

    return res;
}


Vector3 GenerateIllumination(Ray ray, int type, double t, int idx, double alpha, double beta) {
    // for all lights, first draw a shadow ray and check intersection
    // if intersects, there is a shadow
    // else create coloring
    // t is final t value, idx is the intersection idx of object
    Vector3 res;

    for (int i = 0; i < num_lights; ++i) {
        // create shadow ray
        Ray shadow;
        Vector3 intersection = ray.dir.Multiply(t);
        shadow.start.Set(intersection);

        Vector3 light(
            lights[i].position[0],
            lights[i].position[1], 
            lights[i].position[2]);
        shadow.dir.Set(Minus(light, shadow.start));
        shadow.dir = shadow.dir.Normalize();

        // check if in shadow
        bool shadow_meets_sphere = isBlockedBySphere(shadow, light);
        bool shadow_meets_tri = isBlockedByTri(shadow, light);
        if (shadow_meets_sphere || shadow_meets_tri) {
            res = res.addToSelf(Vector3(0.0, 0.0, 0.0));
        }
        else {
            // get corresponding color
            if (type == SPHERE_INRCT) { // sphere intersection
                res = res.addToSelf(GenerateSpherePhongColor(ray, shadow, i, t, idx));
            }
            else if (type == TRI_INRCT) { // tri intersection
                res = res.addToSelf(GenerateTriPhongColor(ray, shadow, i, t, idx, alpha, beta));
            }
        }
    }

    // add ambient
    res = AddAmbient(res).Multiply(255.0);

    return res;
}

Vector3 GenerateColor(Ray ray) {
    // Genral per pixel color generator
    // checks for intersection, then draw color based on intersection type
    Vector3 res;

    // get corresponding t and idx values from intersections
    double t_sph = DBL_MAX;
    int idx_sph = -1;
    bool has_sphere_intersection = isSphereIntersect(ray, t_sph, idx_sph);

    double alpha = -1; double beta = -1;
    double t_tri = DBL_MAX;
    int idx_tri = -1;
    bool has_tri_intersection = isTriangleIntersect(ray, t_tri, idx_tri, alpha , beta);

    // decide on final t and idx values, generate color
    double t; int idx;
    int type = GetClosestIntersectionType(t_sph, t_tri, idx_sph, idx_tri, t, idx);
    if (type == NONE_INRCT) {
        res = Vector3(255.0, 255.0, 255.0); // background
    }
    else {
        res.Set(GenerateIllumination(ray, type, t, idx, alpha, beta));
    }

    return res;
}


//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      Ray ray = GenerateRay(static_cast<double>(x), static_cast<double>(y));
      Vector3 color = GenerateColor(ray);
      plot_pixel(x, y, color.x, color.y, color.z);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

