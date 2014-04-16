#define _USE_MATH_DEFINES
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

bool Sphere::intersect(const Ray &r, Hit &h) const {

	// ==========================================
	// ASSIGNMENT:  IMPLEMENT SPHERE INTERSECTION
	// ==========================================
	// plug the explicit ray equation into the implict sphere equation and solve
	// return true if the sphere was intersected, and update the hit
	// data structure to contain the value of t for the ray at the
	// intersection point, the material, and the normal

	//WORKINJG COPY OF CODE
	
	//Variables for solving the quadratic equation from sphere/ray intersection
	glm::vec3 o = r.getOrigin() - this->center;
	glm::vec3 d = r.getDirection();
	glm::normalize(d);
	float A = glm::dot(d, d);
	float B = 2 * (glm::dot(d, o));
	float C = glm::dot(o, o) - this->radius*this->radius;
	float t, t0, t1;
	glm::vec3 hit;

	//discriminant check
	float disc = B*B - 4*A*C;
	if (disc < 0)
	{
		return false;
	}
	else if (disc = 0)
	{
		t = (-B) / (2 * A);
	}
	else
	{
		float discSqrt = sqrtf(disc);;

		t0 = (-B - discSqrt) / (2 * A);
		t1 = (-B + discSqrt) / (2 * A);
		if (abs(t0) < .001)
			t0 = 0;
		if (abs(t1) < .001)
			t1 = 0;

		if (t0 > t1)
		{
			float temp = t0;
			t0 = t1;
			t1 = temp;
		}
		if (t1 < 0)
		{
			return false;
		}
		if (t0 < 0)
		{
			t = t1;
		}
		else
		{
			t = t0;
		}
	}
	
	float fd = glm::distance(r.getOrigin() + r.getDirection()*t, this->center);
	if (0 <= (this->radius*this->radius) - (fd*fd))
	{
		float sqrtfd = sqrt((this->radius*this->radius) - (fd*fd));
		float adjust = abs(sqrtfd);
		t -= adjust;
	}
	hit = r.getOrigin() + r.getDirection()*t;
	glm::vec3 normal = glm::normalize(hit - this->center);

	if (t < h.getT())
	{
		h.set(t, this->material, normal);
	}
	return true;
	
	/*
	//Variables for solving the quadratic equation from sphere/ray intersection
	glm::vec3 o = r.getOrigin() - this->center;
	glm::vec3 d = r.getDirection();
	glm::normalize(d);
	float A = glm::dot(d, d);
	float B = 2 * (glm::dot(d, o));
	float C = glm::dot(o, o) - this->radius*this->radius;
	float t;
	glm::vec3 hit;

	//discriminant check
	float disc = B * B - 4 * A * C;
	if (disc < 0)
	{
		return false;
	}
	else if (disc = 0)
	{
		t = (-B) / (2 * A);
	}
	else
	{
		float det = sqrtf(disc);
		if (-B - det < 0 && -B + det < 0)
		{
			return false;
		}
		if ()
		t0 = (-B - discSqrt) / (2 * A);
		t1 = (-B + discSqrt) / (2 * A);
		if (abs(t0) < .001)
			t0 = 0;
		if (abs(t1) < .001)
			t1 = 0;

		if (t0 > t1)
		{
			float temp = t0;
			t0 = t1;
			t1 = temp;
		}
		if (t1 < 0)
		{
			return false;
		}
		if (t0 < 0)
		{
			t = t1;
		}
		else
		{
			t = t0;
		}
	}

	float fd = glm::distance(r.getOrigin() + r.getDirection()*t, this->center);
	if (0 <= (this->radius*this->radius) - (fd*fd))
	{
		float sqrtfd = sqrt((this->radius*this->radius) - (fd*fd));
		float adjust = abs(sqrtfd);
		t -= adjust;
	}
	hit = r.getOrigin() + r.getDirection()*t;
	glm::vec3 normal = glm::normalize(hit - this->center);

	if (t < h.getT())
	{
		h.set(t, this->material, normal);
	}
	return true;
	*/
} 

// ====================================================================
// ====================================================================

// helper function to place a grid of points on the sphere
glm::vec3 ComputeSpherePoint(float s, float t, const glm::vec3 center, float radius) {
  float angle = 2*M_PI*s;
  float y = -cos(M_PI*t);
  float factor = sqrt(1-y*y);
  float x = factor*cos(angle);
  float z = factor*-sin(angle);
  glm::vec3 answer = glm::vec3(x,y,z);
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
  m->addVertex(center+radius*glm::vec3(0,-1,0));  // bottom
  for (j = 1; j < v; j++) {  // middle
    for (i = 0; i < h; i++) {
      float s = i / float(h);
      float t = j / float(v);
      m->addVertex(ComputeSpherePoint(s,t,center,radius));
    }
  }
  m->addVertex(center+radius*glm::vec3(0,1,0));  // top

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
