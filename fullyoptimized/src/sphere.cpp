#include "utils.h"
#include "material.h"
#include "argparser.h"
#include "sphere.h"
#include "vertex.h"
#include "mesh.h"
#include "ray.h"
#include "hit.h"

// ====================================================================
// ====================================================================
double getMin(double t1, double t2)
{
	if (t1<t2)
		return t1;
	return t2;
}
bool Sphere::intersect(const Ray &r, Hit &h) const {



  // ==========================================
  // ASSIGNMENT:  IMPLEMENT SPHERE INTERSECTION
  // ==========================================

  // plug the explicit ray equation into the implict sphere equation and solve
  Vec3f temp = r.getOrigin()-center;
  double b = 2 * r.getDirection().Dot3(temp);
  double c = (temp).Dot3(temp) - radius*radius;
  double d = b*b-4*c;
  if (d<0)
	  return false;
  double t1 = (-b+sqrt(d))/2;
  double t2 = (-b-sqrt(d))/2;
  double t;
  if (t1<=EPSILON)
	  t=t2;
  else if (t2<=EPSILON)
	  t=t1;
  else
	  t=getMin(t1,t2);
  if (t<=EPSILON)
	  return false;
  if (t<h.getT())
  {
	  Vec3f normal = r.getOrigin()+t*r.getDirection()-center;
  	  h.set(t,material,normal*(1/normal.Length()));
  }
  else
	  return false;
  // return true if the sphere was intersected, and update
  // the hit data structure to contain the value of t for the ray at
  // the intersection point, the material, and the normal
  return true;
} 

// ====================================================================
// ====================================================================

// helper function to place a grid of points on the sphere
Vec3f ComputeSpherePoint(double s, double t, const Vec3f center, double radius) {
  double angle = 2*M_PI*s;
  double y = -cos(M_PI*t);
  double factor = sqrt(1-y*y);
  double x = factor*cos(angle);
  double z = factor*-sin(angle);
  Vec3f answer = Vec3f(x,y,z);
  answer *= radius;
  answer += center;
  return answer;
}

void Sphere::addRasterizedFaces(Mesh *m, ArgParser *args) {
  
  // and convert it into quad patches for radiosity
  int h = args->sphere_horiz;
  int v = args->sphere_vert;
  assert (h % 2 == 0);
  int i,j;
  int va,vb,vc,vd;
  Vertex *a,*b,*c,*d;
  int offset = m->numVertices(); //vertices.size();

  // place vertices
  m->addVertex(center+radius*Vec3f(0,-1,0));  // bottom
  for (j = 1; j < v; j++) {  // middle
    for (i = 0; i < h; i++) {
      double s = i / double(h);
      double t = j / double(v);
      m->addVertex(ComputeSpherePoint(s,t,center,radius));
    }
  }
  m->addVertex(center+radius*Vec3f(0,1,0));  // top

  // the middle patches
  for (j = 1; j < v-1; j++) {
    for (i = 0; i < h; i++) {
      va = 1 +  i      + h*(j-1);
      vb = 1 + (i+1)%h + h*(j-1);
      vc = 1 +  i      + h*(j);
      vd = 1 + (i+1)%h + h*(j);
      a = m->getVertex(offset + va);
      b = m->getVertex(offset + vb);
      c = m->getVertex(offset + vc);
      d = m->getVertex(offset + vd);
      m->addRasterizedPrimitiveFace(a,b,d,c,material);
    }
  }

  for (i = 0; i < h; i+=2) {
    // the bottom patches
    va = 0;
    vb = 1 +  i;
    vc = 1 + (i+1)%h;
    vd = 1 + (i+2)%h;
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(d,c,b,a,material);
    // the top patches
    va = 1 + h*(v-1);
    vb = 1 +  i      + h*(v-2);
    vc = 1 + (i+1)%h + h*(v-2);
    vd = 1 + (i+2)%h + h*(v-2);
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(b,c,d,a,material);
  }
}
