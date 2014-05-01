#include "raytracer.h"
#include "material.h"
#include "vectors.h"
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
Vec3f RayTracer::TraceRay(Ray &ray, Hit &hit, int bounce_count,int count) const {

  // First cast a ray and see if we hit anything.
  hit = Hit();
  bool intersect = CastRay(ray,hit,false);

  // if there is no intersection, simply return the background color
  if (intersect == false) {
    return Vec3f(srgb_to_linear(mesh->background_color.r()),
		 srgb_to_linear(mesh->background_color.g()),
		 srgb_to_linear(mesh->background_color.b()));
    
  }
  // otherwise decide what to do based on the material
  Material *m = hit.getMaterial();
  assert (m != NULL);

  // rays coming from the light source are set to white, don't bother to ray trace further.
  if (m->getEmittedColor().Length() > 0.001) {
    return m->getEmittedColor();
  } 
 
  
  Vec3f normal = hit.getNormal();
  Vec3f point = ray.pointAtParameter(hit.getT());
  Vec3f answer;

  // ----------------------------------------------
  //  start with the indirect light (ambient light)
  Vec3f diffuse_color = m->getDiffuseColor(hit.get_s(),hit.get_t());
  if (args->gather_indirect) {
    // photon mapping for more accurate indirect light
    answer = diffuse_color * (photon_mapping->GatherIndirect(point, normal, ray.getDirection()) + args->ambient_light);
  } else {
    // the usual ray tracing hack for indirect light
    answer = diffuse_color * args->ambient_light;
  }      
  // ----------------------------------------------
  // add contributions from each light that is not in shadow
	int num_lights = mesh->getLights().size();
	for (int i = 0; i < num_lights; i++) {

		Face *f = mesh->getLights()[i];
		Vec3f lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
		Vec3f myLightColor;
		Vec3f lightCentroid = f->computeCentroid();
		Vec3f dirToLightCentroid = lightCentroid - point;
		dirToLightCentroid.Normalize();

		float distToLightCentroid = (lightCentroid - point).Length();
		myLightColor = lightColor / float(M_PI*distToLightCentroid*distToLightCentroid);

		//For each shadow sample, cast a light ray
		if (args->num_shadow_samples <= 1)
		{
			args->num_shadow_samples = 1;
		}
		//Adaptive samples
		//std::vector<glm::vec3> shading;
		Vec3f shadeanswer(0, 0, 0);
		std::vector<Vec3f> shading;
		int shadecounter = 0;
		for (int i = 0; i < 4; i++)
		{
			Vec3f lightCorner;
			lightCorner = (*f)[i]->get();
			dirToLightCentroid = (lightCorner - point);
			dirToLightCentroid.Normalize();

			//Cast a ray towards the light source
			Ray r(point, dirToLightCentroid);
			Hit Lhit = Hit();
			bool Lintersect = CastRay(r, Lhit, false);
			// If there is no intersection, cast a ray to the light source
			// The ray will hit an object or the light source, if it is the light the distance of the ray Hit will be ~distToLightCentroid
			distToLightCentroid = (lightCorner - point).Length();
			if (Lhit.getT() > distToLightCentroid - (float)EPSILON)
			{
				RayTree::AddShadowSegment(r, 0, distToLightCentroid);
				Vec3f part = m->Shade(ray, hit, dirToLightCentroid, myLightColor, args);
				shading.push_back(part);
				part *= .25;
				shadeanswer += part;
				shadecounter++;
			}
		}

		bool withinBounds = true;
		//No obstructions
		if (shadecounter == 4)
		{
			//return shadeanswer;
			answer += shadeanswer;
		}
		//All obstructions
		else if (shadecounter==0)
		{
			answer += shadeanswer;
		}
		//Some obstructions
		else
		{
			shadeanswer = Vec3f(0, 0, 0);
			//return Vec3f(1,0,0);

			for (int i = 0; i < args->num_shadow_samples; i++)
			{
				lightCentroid = f->RandomPoint();

				dirToLightCentroid = (lightCentroid - point);
				dirToLightCentroid.Normalize();

				//Cast a ray towards the light source
				Ray r(point, dirToLightCentroid);
				Hit Lhit = Hit();
				bool Lintersect = CastRay(r, Lhit, false);
				// If there is no intersection, cast a ray to the light source
				// The ray will hit an object or the light source, if it is the light the distance of the ray Hit will be ~distToLightCentroid
				distToLightCentroid = (lightCentroid - point).Length();
				if (Lhit.getT() > distToLightCentroid - (float)EPSILON)
				{
					RayTree::AddShadowSegment(r, 0, distToLightCentroid);
					Vec3f part = m->Shade(ray, hit, dirToLightCentroid, myLightColor, args);

					shadeanswer += part;
				}
			}
			shadeanswer *= 1.0/args->num_shadow_samples;
			answer += shadeanswer;
		}
	}
  // ----------------------------------------------
  // add contribution from reflection, if the surface is shiny
  Vec3f reflectiveColor = m->getReflectiveColor();



  // =================================
  // ASSIGNMENT:  ADD REFLECTIVE LOGIC
  // =================================
  if (count<bounce_count)
  {
	  Ray r2( point,ray.getDirection()-2*(ray.getDirection().Dot3(normal))*normal);
	  Hit h2;
	  answer+=reflectiveColor*TraceRay(r2,h2,bounce_count,count+1);
	  RayTree::AddReflectedSegment(r2,0,h2.getT());
  }


  return answer; 
}

