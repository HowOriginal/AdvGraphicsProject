#define _USE_MATH_DEFINES
#include "glCanvas.h"
#include <windows.h>
#include "raytracer.h"
#include "material.h"
#include "argparser.h"
#include "raytree.h"
#include "utils.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "photon_mapping.h"


// ===========================================================================
// casts a single ray through the scene geometry and finds the closest hit
bool RayTracer::CastRay(const Ray &ray, Hit &h, bool use_rasterized_patches) const {
  bool answer = false;

  // intersect each of the quads
  for (int i = 0; i < mesh->numOriginalQuads(); i++) {
    Face *f = mesh->getOriginalQuad(i);
    if (f->intersect(ray,h,args->intersect_backfacing)) answer = true;
  }

  // intersect each of the primitives (either the patches, or the original primitives)
  if (use_rasterized_patches) {
    for (int i = 0; i < mesh->numRasterizedPrimitiveFaces(); i++) {
      Face *f = mesh->getRasterizedPrimitiveFace(i);
      if (f->intersect(ray,h,args->intersect_backfacing)) answer = true;
    }
  } else {
    int num_primitives = mesh->numPrimitives();
    for (int i = 0; i < num_primitives; i++) {
      if (mesh->getPrimitive(i)->intersect(ray,h)) answer = true;
    }
  }
  return answer;
}

// ===========================================================================
// does the recursive (shadow rays & recursive rays) work
glm::vec3 RayTracer::TraceRay(Ray &ray, Hit &hit, int bounce_count) const {
	//Move ray forward a distance of epsilon to avoid speckling
	ray = Ray(ray.getOrigin() + ray.getDirection() * .01f, ray.getDirection());

	// First cast a ray and see if we hit anything.
	hit = Hit();
	bool intersect = CastRay(ray, hit, false);

	// if there is no intersection, simply return the background color
	if (intersect == false) {
		return glm::vec3(srgb_to_linear(mesh->background_color.r),
			srgb_to_linear(mesh->background_color.g),
			srgb_to_linear(mesh->background_color.b));
	}

	// otherwise decide what to do based on the material
	Material *m = hit.getMaterial();
	assert(m != NULL);
	// rays coming from the light source are set to white, don't bother to ray trace further.
	if (glm::length(m->getEmittedColor()) > 0.001) {
		return m->getEmittedColor();
	}


	glm::vec3 normal = hit.getNormal();
	glm::vec3 point = ray.pointAtParameter(hit.getT());
	glm::vec3 answer;

	//  start with the indirect light (ambient light)
	glm::vec3 diffuse_color = m->getDiffuseColor(hit.get_s(), hit.get_t());
	if (args->gather_indirect) {
		// photon mapping for more accurate indirect light
		answer = diffuse_color * (photon_mapping->GatherIndirect(point, normal, ray.getDirection()) + args->ambient_light);
	}
	else {
		// the usual ray tracing hack for indirect light
		answer = diffuse_color * args->ambient_light;
	}

	// ----------------------------------------------
	// add contributions from each light that is not in shadow
	int num_lights = mesh->getLights().size();
	for (int i = 0; i < num_lights; i++)
	{
		Face *f = mesh->getLights()[i];
		glm::vec3 lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
		glm::vec3 myLightColor;
		glm::vec3 lightCentroid = f->computeCentroid();
		glm::vec3 dirToLightCentroid = glm::normalize(lightCentroid - point);

		float distToLightCentroid = glm::length(lightCentroid - point);
		myLightColor = lightColor / float(M_PI*distToLightCentroid*distToLightCentroid);

		//For each shadow sample, cast a light ray
		if (args->num_shadow_samples <= 1)
		{
			args->num_shadow_samples = 1;
		}

		//Adaptive samples
		//std::vector<glm::vec3> shading;
		glm::vec3 shadeanswer(0, 0, 0);
		std::vector<glm::vec3> shading;
		int shadecounter = 0;
		for (int i = 0; i < 4; i++)
		{
			glm::vec3 lightCorner;
			lightCorner = (*f)[i]->get();
			dirToLightCentroid = glm::normalize(lightCorner - point);

			//Cast a ray towards the light source
			Ray r(point, dirToLightCentroid);
			Hit Lhit = Hit();
			bool Lintersect = CastRay(r, Lhit, false);
			// If there is no intersection, cast a ray to the light source
			// The ray will hit an object or the light source, if it is the light the distance of the ray Hit will be ~distToLightCentroid
			distToLightCentroid = glm::length(lightCorner - point);
			if (Lhit.getT() > distToLightCentroid - (float)EPSILON)
			{
				RayTree::AddShadowSegment(r, 0, distToLightCentroid);
				glm::vec3 part = m->Shade(ray, hit, dirToLightCentroid, myLightColor, args);
				shading.push_back(part);
				part /= 4;
				shadeanswer += part;
				shadecounter++;
			}
		}
		bool withinBounds = true;
		//No obstructions
		if (shadecounter == 4)
		{
			answer += shadeanswer;
		}
		//Some obstructions
		else
		{
			// shading is a vector for the pixel corners
			glm::vec3 moreshade(0, 0, 0);
			std::vector<glm::vec3> aslight;
			aslight.push_back(((*f)[0]->get() + (*f)[1]->get()) / (float)2);
			aslight.push_back(((*f)[0]->get() + (*f)[3]->get()) / (float)2);
			aslight.push_back(((*f)[2]->get() + (*f)[1]->get()) / (float)2);
			aslight.push_back(((*f)[2]->get() + (*f)[3]->get()) / (float)2);

			aslight.push_back(f->computeCentroid());

			for (int i = 0; i < aslight.size(); i++)
			{
				lightCentroid = aslight[i];

				dirToLightCentroid = glm::normalize(lightCentroid - point);

				//Cast a ray towards the light source
				Ray r(point, dirToLightCentroid);
				Hit Lhit = Hit();
				bool Lintersect = CastRay(r, Lhit, false);
				// If there is no intersection, cast a ray to the light source
				// The ray will hit an object or the light source, if it is the light the distance of the ray Hit will be ~distToLightCentroid
				distToLightCentroid = glm::length(lightCentroid - point);
				if (Lhit.getT() > distToLightCentroid - (float)EPSILON)
				{
					RayTree::AddShadowSegment(r, 0, distToLightCentroid);
					glm::vec3 part = m->Shade(ray, hit, dirToLightCentroid, myLightColor, args);
					part /= aslight.size();
					moreshade += part;
				}
			}
			//shadeanswer *= 4;
			shadeanswer += moreshade;
			shadeanswer /= 2;
			answer += moreshade;
		}
	}

	// ----------------------------------------------
	// add contribution from reflection, if the surface is shiny
	glm::vec3 reflectiveColor = m->getReflectiveColor();
	if (glm::length(reflectiveColor) <= .001)
	{
		return answer;
	}

	//Reflective recursive rays
	if (bounce_count > 0)
	{
		//Error limit
		const float ERR = 1e-12;
		//R = V - 2N(V.N)  /or/ R = V - N(V.N) - N(V.N)
		glm::vec3 reflection = ray.getDirection() - normal*(glm::dot(ray.getDirection(), normal)) - normal*(glm::dot(ray.getDirection(), normal));
		reflection += normal * ERR;
		Ray reflect = Ray(point, reflection);
		Hit Rhit = Hit();
		glm::vec3 reflectanswer = TraceRay(reflect, Rhit, bounce_count - 1);
		if (Rhit.getMaterial() != NULL)
		{
			RayTree::AddReflectedSegment(reflect, 0, glm::length(reflect.pointAtParameter(Rhit.getT()) - point));
		}
		else
		{
			RayTree::AddReflectedSegment(reflect, 0, 1);
		}
		answer += reflectanswer * reflectiveColor;

	}
	return answer;
}


void RayTracer::initializeVBOs() {
  glGenBuffers(1, &pixels_a_VBO);
  glGenBuffers(1, &pixels_b_VBO);
  glGenBuffers(1, &pixels_indices_a_VBO);
  glGenBuffers(1, &pixels_indices_b_VBO);
  render_to_a = true;
}


void RayTracer::resetVBOs() {

  pixels_a.clear();
  pixels_b.clear();

  pixels_indices_a.clear();
  pixels_indices_b.clear();

  render_to_a = true;
}

void RayTracer::setupVBOs() {
	glBindBuffer(GL_ARRAY_BUFFER, pixels_a_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(VBOPosNormalColor)*pixels_a.size(), pixels_a.size() > 0 ? &pixels_a[0] : NULL, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, pixels_b_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(VBOPosNormalColor)*pixels_b.size(), pixels_b.size() > 0 ? &pixels_b[0] : NULL, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pixels_indices_a_VBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		sizeof(VBOIndexedTri)* pixels_indices_a.size(),
		pixels_indices_a.size() > 0 ? &pixels_indices_a[0] : NULL, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pixels_indices_b_VBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		sizeof(VBOIndexedTri)* pixels_indices_b.size(),
		NULL, GL_STATIC_DRAW);
}

void RayTracer::drawVBOs() {
  // turn off lighting
  glUniform1i(GLCanvas::colormodeID, 0);
  // turn off depth buffer
  glDisable(GL_DEPTH_TEST);

  if (render_to_a) {
    drawVBOs_b();
    drawVBOs_a();
  } else {
    drawVBOs_a();
    drawVBOs_b();
  }

  glEnable(GL_DEPTH_TEST);
}

void RayTracer::drawVBOs_a() {
  if (pixels_a.size() == 0) return;
  glBindBuffer(GL_ARRAY_BUFFER, pixels_a_VBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,pixels_indices_a_VBO); 
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
  glDrawElements(GL_TRIANGLES,
                 pixels_indices_a.size()*3,
                 GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(3);
}

void RayTracer::drawVBOs_b() {
  if (pixels_b.size() == 0) return;
  glBindBuffer(GL_ARRAY_BUFFER, pixels_b_VBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,pixels_indices_b_VBO); 
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor),(void*)sizeof(glm::vec3) );
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT,GL_FALSE,sizeof(VBOPosNormalColor), (void*)(sizeof(glm::vec3)*2 + sizeof(glm::vec4)));
  glDrawElements(GL_TRIANGLES,
                 pixels_indices_b.size()*3,
                 GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(3);
}


void RayTracer::cleanupVBOs() {
  glDeleteBuffers(1, &pixels_a_VBO);
  glDeleteBuffers(1, &pixels_b_VBO);
}

