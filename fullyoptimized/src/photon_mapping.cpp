#include "glCanvas.h"

#include <iostream>
#include <algorithm>

#include "argparser.h"
#include "photon_mapping.h"
#include "mesh.h"
#include "matrix.h"
#include "face.h"
#include "primitive.h"
#include "kdtree.h"
#include "utils.h"
#include "raytracer.h"

// ==========
// DESTRUCTOR
PhotonMapping::~PhotonMapping() {
  // cleanup all the photons
  delete kdtree;
}

// ========================================================================
// Recursively trace a single photon

void PhotonMapping::TracePhoton(const Vec3f &position, const Vec3f &direction, 
				const Vec3f &energy, int iter) {


  // ==============================================
  // ASSIGNMENT: IMPLEMENT RECURSIVE PHOTON TRACING
  // ==============================================

  // Trace the photon through the scene.  At each diffuse or
  // reflective bounce, store the photon in the kd tree.

  // One optimization is to *not* store the first bounce, since that
  // direct light can be efficiently computed using classic ray
  // tracing.

  //do ray cast
  Ray r(position,direction*(1/direction.Length()));
  Hit h;
  raytracer->CastRay(r,h,true);
  if (h.getT()>1000)
	  return;
  MTRand mtrand;
  Vec3f refl = h.getMaterial()->getReflectiveColor();
  Vec3f diff = h.getMaterial()->getDiffuseColor();
  double ran=mtrand.rand();
  if (iter==0)
	  ran= mtrand.rand(refl.Length()+diff.Length());
  //std::cout<<iter<<" "<<h.getT()<<" "<<refl.Length()+diff.Length()<<std::endl;
  //send reflective photon
  if (iter<args->num_bounces&&ran<=refl.Length())
	  TracePhoton(r.pointAtParameter(h.getT()),r.getDirection()-2*(r.getDirection().Dot3(h.getNormal()))*h.getNormal(),energy,iter+1);
  else if (iter<args->num_bounces&&ran<=refl.Length()+diff.Length())
	  TracePhoton(r.pointAtParameter(h.getT()),RandomDiffuseDirection(h.getNormal()),energy,iter+1);
  else
  {
	  Photon p(position,direction,energy,iter);
	  kdtree->AddPhoton(p);
  }


}


// ========================================================================
// Trace the specified number of photons through the scene

void PhotonMapping::TracePhotons() {
  std::cout << "trace photons" << std::endl;

  // first, throw away any existing photons
  delete kdtree;

  // consruct a kdtree to store the photons
  BoundingBox *bb = mesh->getBoundingBox();
  Vec3f min = bb->getMin();
  Vec3f max = bb->getMax();
  Vec3f diff = max-min;
  min -= 0.001*diff;
  max += 0.001*diff;
  kdtree = new KDTree(BoundingBox(min,max));

  // photons emanate from the light sources
  const std::vector<Face*>& lights = mesh->getLights();

  // compute the total area of the lights
  double total_lights_area = 0;
  for (unsigned int i = 0; i < lights.size(); i++) {
    total_lights_area += lights[i]->getArea();
  }

  // shoot a constant number of photons per unit area of light source
  // (alternatively, this could be based on the total energy of each light)
  for (unsigned int i = 0; i < lights.size(); i++) {  
    double my_area = lights[i]->getArea();
    int num = args->num_photons_to_shoot * my_area / total_lights_area;
    // the initial energy for this photon
    Vec3f energy = my_area/double(num) * lights[i]->getMaterial()->getEmittedColor();
    Vec3f normal = lights[i]->computeNormal();
    for (int j = 0; j < num; j++) {
      Vec3f start = lights[i]->RandomPoint();
      // the initial direction for this photon (for diffuse light sources)
      Vec3f direction = RandomDiffuseDirection(normal);
      TracePhoton(start,direction,energy,0);
    }
  }
}


// ======================================================================
// During ray tracing, when a diffuse (or partially diffuse) object is
// hit, gather the nearby photons to approximate indirect illumination
bool comparePhotons(const Photon& p1,const Photon& p2)
{
	double dis1 = (p1.getPosition()-p1.getPoint()).Length();
	double dis2 = (p2.getPosition()-p2.getPoint()).Length();
	return dis1<dis2;
}
Vec3f PhotonMapping::GatherIndirect(const Vec3f &point, const Vec3f &normal, const Vec3f &direction_from) const {


  if (kdtree == NULL) { 
    std::cout << "WARNING: Photons have not been traced throughout the scene." << std::endl;
    return Vec3f(0,0,0); 
  }


  // ================================================================
  // ASSIGNMENT: GATHER THE INDIRECT ILLUMINATION FROM THE PHOTON MAP
  // ================================================================


  // collect the closest args->num_photons_to_collect photons
  // determine the radius that was necessary to collect that many photons
  // average the energy of those photons over that radius
  double radius=1;

  std::vector<Photon> potentials;
  while (potentials.size()<args->num_photons_to_collect)
  {
	  BoundingBox b(point-Vec3f(radius,radius,radius),point+Vec3f(radius,radius,radius));

	  kdtree->CollectPhotonsInBox(b,potentials);
	  for (int i=0;i<potentials.size();i++)
	  {
		  if (normal.Dot3(potentials[i].getDirectionFrom())>0)
		  {
			  potentials.erase(potentials.begin()+i);
			  i--;
		  }
		  else if ((potentials[i].getPosition()-point).Length()>radius)
		  {
			  potentials.erase(potentials.begin()+i);
			  i--;
		  }
	  }

	  radius*=2;


  }
  assert(potentials.size()>=args->num_photons_to_collect);
  for (unsigned int i=0;i<potentials.size();i++)
	  potentials[i].setPoint(point);
  std::sort(potentials.begin(),potentials.end(),comparePhotons);

  radius = (potentials[args->num_photons_to_collect-1].getPosition()-point).Length()/2;
  Vec3f color;
  for (int i=0;i<args->num_photons_to_collect;i++)
	  color+=potentials[i].getEnergy();
  color*=1/(radius*radius);
  //std::cout<<radius<<"\n";
  // return the color
  return color;
}



// ======================================================================
// PHOTON VISUALIZATION FOR DEBUGGING
// ======================================================================


void PhotonMapping::initializeVBOs() {
  glGenBuffers(1, &photon_verts_VBO);
  glGenBuffers(1, &photon_direction_indices_VBO);
  glGenBuffers(1, &kdtree_verts_VBO);
  glGenBuffers(1, &kdtree_edge_indices_VBO);
}

void PhotonMapping::setupVBOs() {
  photon_verts.clear();
  photon_direction_indices.clear();
  kdtree_verts.clear();
  kdtree_edge_indices.clear();

  // initialize the data
  int dir_count = 0;
  int edge_count = 0;
  BoundingBox *bb = mesh->getBoundingBox();
  double max_dim = bb->maxDim();

  if (kdtree == NULL) return;
  std::vector<const KDTree*> todo;  
  todo.push_back(kdtree);
  while (!todo.empty()) {
    const KDTree *node = todo.back();
    todo.pop_back(); 
    if (node->isLeaf()) {
      const std::vector<Photon> &photons = node->getPhotons();
      int num_photons = photons.size();
      for (int i = 0; i < num_photons; i++) {
	const Photon &p = photons[i];
	Vec3f energy = p.getEnergy()*args->num_photons_to_shoot;
	const Vec3f &position = p.getPosition();
	Vec3f other = position + p.getDirectionFrom()*0.02*max_dim;
	photon_verts.push_back(VBOPosColor(position,energy));
	photon_verts.push_back(VBOPosColor(other,energy));
	photon_direction_indices.push_back(VBOIndexedEdge(dir_count,dir_count+1)); dir_count+=2;
      }

      // initialize kdtree vbo
      const Vec3f& min = node->getMin();
      const Vec3f& max = node->getMax();
      kdtree_verts.push_back(VBOPos(Vec3f(min.x(),min.y(),min.z())));
      kdtree_verts.push_back(VBOPos(Vec3f(min.x(),min.y(),max.z())));
      kdtree_verts.push_back(VBOPos(Vec3f(min.x(),max.y(),min.z())));
      kdtree_verts.push_back(VBOPos(Vec3f(min.x(),max.y(),max.z())));
      kdtree_verts.push_back(VBOPos(Vec3f(max.x(),min.y(),min.z())));
      kdtree_verts.push_back(VBOPos(Vec3f(max.x(),min.y(),max.z())));
      kdtree_verts.push_back(VBOPos(Vec3f(max.x(),max.y(),min.z())));
      kdtree_verts.push_back(VBOPos(Vec3f(max.x(),max.y(),max.z())));

      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count  ,edge_count+1)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+1,edge_count+3)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+3,edge_count+2)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+2,edge_count  )); 

      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+4,edge_count+5)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+5,edge_count+7)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+7,edge_count+6)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+6,edge_count+4)); 

      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count  ,edge_count+4)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+1,edge_count+5)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+2,edge_count+6)); 
      kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+3,edge_count+7)); 


      edge_count += 8;

    } else {
      todo.push_back(node->getChild1());
      todo.push_back(node->getChild2());
    } 
  }
  assert (2*photon_direction_indices.size() == photon_verts.size());
  int num_directions = photon_direction_indices.size();
  int num_edges = kdtree_edge_indices.size();

  // cleanup old buffer data (if any)
  cleanupVBOs();

  // copy the data to each VBO
  if (num_directions > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,photon_verts_VBO); 
    glBufferData(GL_ARRAY_BUFFER,
		 sizeof(VBOPosColor) * num_directions * 2,
		 &photon_verts[0],
		 GL_STATIC_DRAW); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,photon_direction_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * num_directions,
		 &photon_direction_indices[0], GL_STATIC_DRAW);
  } 

  if (num_edges > 0) {
    glBindBuffer(GL_ARRAY_BUFFER,kdtree_verts_VBO); 
    glBufferData(GL_ARRAY_BUFFER,
		 sizeof(VBOPos) * num_edges * 2,
		 &kdtree_verts[0],
		 GL_STATIC_DRAW); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,kdtree_edge_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * num_edges,
		 &kdtree_edge_indices[0], GL_STATIC_DRAW);
  } 
}

void PhotonMapping::drawVBOs() {

  glDisable(GL_LIGHTING);
  if (args->render_photons) {
    int num_directions = photon_direction_indices.size();
    if (num_directions > 0) {
      // render directions
      glLineWidth(1);
      glBindBuffer(GL_ARRAY_BUFFER, photon_verts_VBO);
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, sizeof(VBOPosColor), BUFFER_OFFSET(0));
      glEnableClientState(GL_COLOR_ARRAY);
      glColorPointer(3, GL_FLOAT, sizeof(VBOPosColor), BUFFER_OFFSET(12));
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, photon_direction_indices_VBO);
      glDrawElements(GL_LINES, num_directions*2, GL_UNSIGNED_INT, BUFFER_OFFSET(0));
      glDisableClientState(GL_COLOR_ARRAY);
      glDisableClientState(GL_VERTEX_ARRAY);
      // render hit points
      glPointSize(3);
      glBindBuffer(GL_ARRAY_BUFFER, photon_verts_VBO);
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, sizeof(VBOPosColor)*2, BUFFER_OFFSET(0));
      glEnableClientState(GL_COLOR_ARRAY);
      glColorPointer(3, GL_FLOAT, sizeof(VBOPosColor)*2, BUFFER_OFFSET(12));
      glDrawArrays(GL_POINTS, 0, num_directions);
      glDisableClientState(GL_COLOR_ARRAY);
      glDisableClientState(GL_VERTEX_ARRAY);
    }
  }
  if (args->render_kdtree) {
    int num_edges = kdtree_edge_indices.size();
    if (num_edges > 0) {
      glDisable(GL_LIGHTING);
      // render edges
      glLineWidth(1);
      glColor3f(0,1,1);
      glBindBuffer(GL_ARRAY_BUFFER, kdtree_verts_VBO);
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, sizeof(VBOPos), BUFFER_OFFSET(0));
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, kdtree_edge_indices_VBO);
      glDrawElements(GL_LINES, num_edges*2, GL_UNSIGNED_INT, BUFFER_OFFSET(0));
      glDisableClientState(GL_VERTEX_ARRAY);
    }
  }
}

void PhotonMapping::cleanupVBOs() {
  glDeleteBuffers(1, &photon_verts_VBO);
  glDeleteBuffers(1, &photon_direction_indices_VBO);
  glDeleteBuffers(1, &kdtree_verts_VBO);
  glDeleteBuffers(1, &kdtree_edge_indices_VBO);
}

