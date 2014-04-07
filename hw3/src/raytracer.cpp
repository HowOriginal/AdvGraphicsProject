#define _USE_MATH_DEFINES
#include "glCanvas.h"

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
		return glm::vec3(1, 1, 1);
	}


	glm::vec3 normal = hit.getNormal();
	glm::vec3 point = ray.pointAtParameter(hit.getT());
	glm::vec3 answer;

	bool from_outside;
	float origin_check = glm::dot(ray.getDirection(), normal) / (glm::length(normal) * glm::length(ray.getDirection()));
	if (origin_check > 0)
	{
		from_outside = false;
	}
	else
	{
		from_outside = true;
	}
	//std::cout << from_outside << "\n";
	// ----------------------------------------------
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
	for (int i = 0; i < num_lights; i++) {

		Face *f = mesh->getLights()[i];
		glm::vec3 lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
		glm::vec3 myLightColor;
		glm::vec3 lightCentroid = f->computeCentroid();
		glm::vec3 dirToLightCentroid = glm::normalize(lightCentroid - point);



		// ===========================================
		// ASSIGNMENT:  ADD SHADOW & SOFT SHADOW LOGIC
		// ===========================================

		float distToLightCentroid = glm::length(lightCentroid - point);
		myLightColor = lightColor / float(M_PI*distToLightCentroid*distToLightCentroid);

		//For each shadow sample, cast a light ray
		if (args->num_shadow_samples <= 1)
		{
			//Cast a ray towards the light source
			Ray r(point, dirToLightCentroid);
			Hit Lhit = Hit();
			bool Lintersect = CastRay(r, Lhit, false);
			// If there is no intersection, cast a ray to the light source
			// The ray will hit an object or the light source, if it is the light the distance of the ray Hit will be ~distToLightCentroid
			if (Lhit.getT() > distToLightCentroid - 1)
			{
				RayTree::AddShadowSegment(r, 0, distToLightCentroid);
				answer += m->Shade(ray, hit, dirToLightCentroid, myLightColor, args);
			}
		}
		else
		{
			//Stratified sampling
			int omega = ceil(sqrt(args->num_shadow_samples));

			//Find possible open stratum and choose a number of them to use for random points
			std::vector< std::pair<int, int> > grid;
			std::vector< std::pair<int, int> > stratum;
			
			for (int l = 0; l < omega; l++)
			{
				for (int p = 0; p < omega; p++)
				{
					std::pair<int,int> temp(l,p);
					grid.push_back(temp);
				}
			}
			for (int l = 0; l < args->num_shadow_samples; l++)
			{
				int pick = rand() % grid.size();
				stratum.push_back(grid[pick]);
				grid.erase(grid.begin() + pick);
			}

			glm::vec3 A = (*f)[0]->get();                       //A
			glm::vec3 colDir = (*f)[3]->get() - (*f)[0]->get(); //AD
			glm::vec3 rowDir = (*f)[1]->get() - (*f)[0]->get(); //AB
			colDir /= omega;
			rowDir /= omega;
			//End stratified code
			for (int i = 0; i < args->num_shadow_samples; i++)
			{
				//Stratified code
				
				glm::vec3 col = colDir;
				col *= stratum[i].first;
				glm::vec3 row = rowDir;
				row *= stratum[i].second;
				glm::vec3 stratumpoint = A + col + row +colDir * (float)GLOBAL_MTRAND.rand() + rowDir * (float)GLOBAL_MTRAND.rand();

				lightCentroid = stratumpoint;
				//End stratified code

				//Random sampling
				//lightCentroid = f->RandomPoint();

				dirToLightCentroid = glm::normalize(lightCentroid - point);

				//Cast a ray towards the light source
				Ray r(point, dirToLightCentroid);
				Hit Lhit = Hit();
				bool Lintersect = CastRay(r, Lhit, false);
				// If there is no intersection, cast a ray to the light source
				// The ray will hit an object or the light source, if it is the light the distance of the ray Hit will be ~distToLightCentroid
				if (Lhit.getT() > distToLightCentroid - 1)
				{
					RayTree::AddShadowSegment(r, 0, distToLightCentroid);
					glm::vec3 part = m->Shade(ray, hit, dirToLightCentroid, myLightColor, args);
					part /= args->num_shadow_samples;
					answer += part;
				}
			}
		}
		/*
		//Cast a ray towards the light source
		Ray r(point, dirToLightCentroid);
		Hit Lhit = Hit();
		bool Lintersect = CastRay(r, Lhit, false);
		// If there is no intersection, cast a ray to the light source
		// The ray will hit an object or the light source, if it is the light the distance of the ray Hit will be ~distToLightCentroid
		if (Lhit.getT() > distToLightCentroid - (1e-4))
		{
			RayTree::AddShadowSegment(r, 0, distToLightCentroid);
			answer += m->Shade(ray, hit, dirToLightCentroid, myLightColor, args);
		}
		*/

		//Normal for debugging
		//Ray n(point, normal);
		//RayTree::AddMainSegment(n, 0, .25);



		// add the lighting contribution from this particular light at this point
		// (fix this to check for blockers between the light & this surface)
		//answer += m->Shade(ray,hit,dirToLightCentroid,myLightColor,args);
	}

	// ----------------------------------------------
	// add contribution from reflection, if the surface is shiny
	glm::vec3 reflectiveColor = m->getReflectiveColor();
	//std::cout << glm::length(reflectiveColor) << "\n";
	if (glm::length(reflectiveColor) > .001)
	{
		//answer += reflectiveColor;
		//answer *= m->getReflectiveColor();
	}
	else
	{
		return answer;
	}

  // =================================
  // ASSIGNMENT:  ADD REFLECTIVE LOGIC
  // =================================

  //Reflective recursive rays
  if (bounce_count > 0)
  {
	  if (args->num_glossy_samples > 1)
	  {
		  const float ERR = 1e-12;
		  glm::vec3 reflection = ray.getDirection() - normal*(glm::dot(ray.getDirection(), normal)) - normal*(glm::dot(ray.getDirection(), normal));
		  //Generate additional reflection rays
		  std::vector<glm::vec3> glossyrays;
		  glossyrays.push_back(reflection + (normal * ERR));
		  for (int n = 0; n < args->num_glossy_samples; n++)
		  {
			  //Higher = more refined
			  float randrestraint = 3;
			  glm::vec3 randomvec((GLOBAL_MTRAND.rand() - .5) / randrestraint, (GLOBAL_MTRAND.rand() - .5) / randrestraint, (GLOBAL_MTRAND.rand() - .5) / randrestraint);
			  glm::vec3 newreflection = reflection + randomvec;
			  glm::normalize(newreflection);
			  glossyrays.push_back(newreflection);
		  }
		  for (int n = 0; n < glossyrays.size(); n++)
		  {
			  Ray reflect = Ray(point, glossyrays[n]);
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
			  reflectanswer /= glossyrays.size();
			  answer += reflectanswer * reflectiveColor;
		  }
	  }
	  else
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
			  if (glm::length(Rhit.getMaterial()->getReflectiveColor()) > .001)
			  {
				  //answer += reflectanswer;
			  }
			  //answer += reflectanswer * m->getReflectiveColor();
			  //std::cout << Rhit.getMaterial()->getReflectiveColor().x << " " << Rhit.getMaterial()->getReflectiveColor().y << " " << Rhit.getMaterial()->getReflectiveColor().z << "\n";
		  }
		  else
		  {
			  RayTree::AddReflectedSegment(reflect, 0, 1);
			  //answer += reflectanswer;
		  }
		  //if (Rhit.getMaterial()->getDiffuseColor())
		  //answer += reflectanswer * Rhit.getMaterial()->getReflectiveColor();
		  answer += reflectanswer * reflectiveColor;
	  }
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

