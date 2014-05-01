#include "glCanvas.h"
#include "radiosity.h"
#include "mesh.h"
#include "face.h"
#include "glCanvas.h"
#include "sphere.h"
#include "raytree.h"
#include "raytracer.h"
#include "utils.h"

// ================================================================
// CONSTRUCTOR & DESTRUCTOR
// ================================================================
Radiosity::Radiosity(Mesh *m, ArgParser *a) {
  mesh = m;
  args = a;
  num_faces = -1;  
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
  Reset();
}

Radiosity::~Radiosity() {
  Cleanup();
  cleanupVBOs();
}

void Radiosity::Cleanup() {
  delete [] formfactors;
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;
  num_faces = -1;
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
}

void Radiosity::Reset() {
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;

  // create and fill the data structures
  num_faces = mesh->numFaces();
  area = new double[num_faces];
  undistributed = new Vec3f[num_faces];
  absorbed = new Vec3f[num_faces];
  radiance = new Vec3f[num_faces];
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    f->setRadiosityPatchIndex(i);
    setArea(i,f->getArea());
    Vec3f emit = f->getMaterial()->getEmittedColor();
    setUndistributed(i,emit);
    setAbsorbed(i,Vec3f(0,0,0));
    setRadiance(i,emit);
  }

  // find the patch with the most undistributed energy
  findMaxUndistributed();
}


// =======================================================================================
// =======================================================================================

void Radiosity::findMaxUndistributed() {
  // find the patch with the most undistributed energy 
  // don't forget that the patches may have different sizes!
  max_undistributed_patch = -1;
  total_undistributed = 0;
  total_area = 0;
  double max = -1;
  for (int i = 0; i < num_faces; i++) {
    double m = getUndistributed(i).Length() * getArea(i);
    total_undistributed += m;
    total_area += getArea(i);
    if (max < m) {
      max = m;
      max_undistributed_patch = i;
    }
  }
  assert (max_undistributed_patch >= 0 && max_undistributed_patch < num_faces);
}


void Radiosity::ComputeFormFactors() {
  assert (formfactors == NULL);
  assert (num_faces > 0);
  formfactors = new double[num_faces*num_faces];


  int isPrint=false;
  // =====================================
  // ASSIGNMENT:  COMPUTE THE FORM FACTORS
  // =====================================
  for (int i=0;i<num_faces;i++)
  {
	  Face* f1 =mesh->getFace(i);
	  double sum=0;
	  for (int j=0;j<num_faces;j++)
	  {
		  if (i==j)
		  {
			  formfactors[i*num_faces+j]=0;
			  continue;
		  }
		  double eps=.3;
		  Face*f2 = mesh->getFace(j);
		  formfactors[i*num_faces+j]=0;
		  for (int k=0;k<args->num_form_factor_samples;k++)
  		  {
			  if (k==0)
			  {
				  Vec3f dir = (f2->computeCentroid()-f1->computeCentroid());
				  dir.Normalize();
				  Ray r(f1->computeCentroid(),dir);
				  Hit h;
				  raytracer->CastRay(r,h,true);
				  float occlusion=1;
				  for (int l=0;l<args->num_shadow_samples;l++)
				  {
					  Ray r2(f2->computeCentroid(),-dir);
					  Hit h2;
					  raytracer->CastRay(r2,h2,true);
					  if (fabs(h2.getT()-h.getT())>eps)
					  {
						  if (isPrint)
							  std::cout<<h2.getT()<<" "<<h.getT()<<"\n";
						  occlusion-=1.0/args->num_shadow_samples;
					  }
				  }
				  double cosi = fabs(f1->computeNormal().Dot3(r.getDirection())/f1->computeNormal().Length());
				  double cosj = fabs(f2->computeNormal().Dot3(r.getDirection())/f2->computeNormal().Length());
				  double rad = (f2->computeCentroid()-f1->computeCentroid()).Length();
				  formfactors[i*num_faces+j]+=occlusion*(1/f1->getArea())*cosi*cosj/(3.14159265358979*rad*rad);
			  }
			  else
			  {
				  Vec3f p1 =f1->RandomPoint();
				  Vec3f p2 =f2->RandomPoint();
				  Vec3f dir = (p2-p1);
				  dir.Normalize();
				  Ray r(p1,dir);
				  Hit h;
				  raytracer->CastRay(r,h,true);
				  float occlusion=1;
				  for (int l=0;l<args->num_shadow_samples;l++)
				  {
					  Ray r2(p2,-dir);
					  Hit h2;
					  raytracer->CastRay(r2,h2,true);
					  if (fabs(h2.getT()-h.getT())>EPSILON)
						  occlusion-=1.0/args->num_shadow_samples;
				  }
				  double cosi = fabs(f1->computeNormal().Dot3(r.getDirection())/f1->computeNormal().Length());
				  double cosj = fabs(f2->computeNormal().Dot3(r.getDirection())/f2->computeNormal().Length());
				  double rad = (p2-p1).Length();
				  formfactors[i*num_faces+j]+=occlusion*(1/f1->getArea())*cosi*cosj/(3.14159265358979*rad*rad);
			  }
		  }
		  formfactors[i*num_faces+j]/=args->num_form_factor_samples;
		  sum+=formfactors[i*num_faces+j];
	  }
	  sum=1/sum;
	  for (int j=0;j<num_faces;j++)
	  {
		  formfactors[i*num_faces+j]*=sum;
	  }
	  isPrint=false;
  }


}

// ================================================================
// ================================================================

double Radiosity::Iterate() {
  if (formfactors == NULL) 
	  ComputeFormFactors();
  assert (formfactors != NULL);


  // ==========================================
  // ASSIGNMENT:  IMPLEMENT RADIOSITY ALGORITHM
  // ==========================================
  double total=0;
  double maxUn=0;
  int maxind=0;
  for (int i=0;i<num_faces;i++)
  {

	  radiance[i]+=undistributed[max_undistributed_patch]*mesh->getFace(i)->getMaterial()->getDiffuseColor()*formfactors[max_undistributed_patch+i*num_faces];
	  absorbed[i]+=undistributed[max_undistributed_patch]*(Vec3f(1,1,1)-mesh->getFace(i)->getMaterial()->getDiffuseColor())*formfactors[max_undistributed_patch+i*num_faces];
	  undistributed[i]+=undistributed[max_undistributed_patch]*mesh->getFace(i)->getMaterial()->getDiffuseColor()*formfactors[max_undistributed_patch+i*num_faces];
	  total+=undistributed[i].Length();
	  if (undistributed[i].Length()>maxUn)
	  {
		  maxind = i;
		  maxUn=undistributed[i].Length();
	  }
  }
  undistributed[max_undistributed_patch].Scale(0);
  max_undistributed_patch=maxind;
  // return the total light yet undistributed
  // (so we can decide when the solution has sufficiently converged)
  return total;
}


// =======================================================================================
// VBO & DISPLAY FUNCTIONS
// =======================================================================================

// for interpolation
void CollectFacesWithVertex(Vertex *have, Face *f, std::vector<Face*> &faces) {
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == f) return;
  }
  if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2] && have != (*f)[3]) return;
  faces.push_back(f);
  for (int i = 0; i < 4; i++) {
    Edge *ea = f->getEdge()->getOpposite();
    Edge *eb = f->getEdge()->getNext()->getOpposite();
    Edge *ec = f->getEdge()->getNext()->getNext()->getOpposite();
    Edge *ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
    if (ea != NULL) CollectFacesWithVertex(have,ea->getFace(),faces);
    if (eb != NULL) CollectFacesWithVertex(have,eb->getFace(),faces);
    if (ec != NULL) CollectFacesWithVertex(have,ec->getFace(),faces);
    if (ed != NULL) CollectFacesWithVertex(have,ed->getFace(),faces);
  }
}

// different visualization modes
Vec3f Radiosity::setupHelperForColor(Face *f, int i, int j) {
  assert (mesh->getFace(i) == f);
  assert (j >= 0 && j < 4);
  if (args->render_mode == RENDER_MATERIALS) {
    return f->getMaterial()->getDiffuseColor();
  } else if (args->render_mode == RENDER_RADIANCE && args->interpolate == true) {
    std::vector<Face*> faces;
    CollectFacesWithVertex((*f)[j],f,faces);
    double total = 0;
    Vec3f color = Vec3f(0,0,0);
    Vec3f normal = f->computeNormal();
    for (unsigned int i = 0; i < faces.size(); i++) {
      Vec3f normal2 = faces[i]->computeNormal();
      double area = faces[i]->getArea();
      if (normal.Dot3(normal2) < 0.5) continue;
      assert (area > 0);
      total += area;
      color += area * getRadiance(faces[i]->getRadiosityPatchIndex());
    }
    assert (total > 0);
    color /= total;
    return color;
  } else if (args->render_mode == RENDER_LIGHTS) {
    return f->getMaterial()->getEmittedColor();
  } else if (args->render_mode == RENDER_UNDISTRIBUTED) { 
    return getUndistributed(i);
  } else if (args->render_mode == RENDER_ABSORBED) {
    return getAbsorbed(i);
  } else if (args->render_mode == RENDER_RADIANCE) {
    return getRadiance(i);
  } else if (args->render_mode == RENDER_FORM_FACTORS) {
    if (formfactors == NULL) ComputeFormFactors();
    double scale = 0.2 * total_area/getArea(i);
    double factor = scale * getFormFactor(max_undistributed_patch,i);
    return Vec3f(factor,factor,factor);
  } else {
    assert(0);
  }
  exit(0);
}


void Radiosity::initializeVBOs() {
  // create a pointer for the vertex & index VBOs
  glGenBuffers(1, &mesh_quad_verts_VBO);
  glGenBuffers(1, &mesh_quad_indices_VBO);
  glGenBuffers(1, &mesh_textured_quad_indices_VBO);
  glGenBuffers(1, &mesh_interior_edge_indices_VBO);
  glGenBuffers(1, &mesh_border_edge_indices_VBO);
}


void Radiosity::setupVBOs() {
  mesh_quad_verts.clear();
  mesh_quad_indices.clear();
  mesh_textured_quad_indices.clear();
  mesh_border_edge_indices.clear();
  mesh_interior_edge_indices.clear();

  // initialize the data in each vector
  int num_faces = mesh->numFaces();
  assert (num_faces > 0);
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    Edge *e = f->getEdge();
    for (int j = 0; j < 4; j++) {
      Vec3f pos = ((*f)[j])->get();
      Vec3f normal = f->computeNormal();
      Vec3f color = setupHelperForColor(f,i,j);
      color = Vec3f(linear_to_srgb(color.r()),
		    linear_to_srgb(color.g()),
		    linear_to_srgb(color.b()));
      mesh_quad_verts.push_back(VBOPosNormalColorTexture(pos,normal,color,(*f)[j]->get_s(),(*f)[j]->get_t()));
      if (e->getOpposite() == NULL) { 
	mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+j,i*4+(j+1)%4));
      } else if (e->getStartVertex()->getIndex() < e->getEndVertex()->getIndex()) {
	mesh_interior_edge_indices.push_back(VBOIndexedEdge(i*4+j,i*4+(j+1)%4));
      }
      e = e->getNext();
    }
    if (f->getMaterial()->hasTextureMap()) {
      mesh_textured_quad_indices.push_back(VBOIndexedQuad(i*4,i*4+1,i*4+2,i*4+3));
    } else {
      mesh_quad_indices.push_back(VBOIndexedQuad(i*4,i*4+1,i*4+2,i*4+3));
    }
    // also outline the max_undistributed patch
    if (args->render_mode == RENDER_FORM_FACTORS && i == max_undistributed_patch) {
      mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+0,i*4+1));
      mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+1,i*4+2));
      mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+2,i*4+3));
      mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+3,i*4+0));
    }
  }
  assert ((int)mesh_quad_verts.size() == num_faces*4);
  assert ((int)mesh_quad_indices.size() + (int)mesh_textured_quad_indices.size() == num_faces);

  // cleanup old buffer data (if any)
  cleanupVBOs();

  // copy the data to each VBO
  glBindBuffer(GL_ARRAY_BUFFER,mesh_quad_verts_VBO); 
  glBufferData(GL_ARRAY_BUFFER,
	       sizeof(VBOPosNormalColorTexture) * num_faces * 4,
	       &mesh_quad_verts[0],
	       GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_quad_indices_VBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
	       sizeof(VBOIndexedQuad) * mesh_quad_indices.size(),
	       &mesh_quad_indices[0], GL_STATIC_DRAW);
  if (mesh_textured_quad_indices.size() > 0) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_textured_quad_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedQuad) * mesh_textured_quad_indices.size(),
		 &mesh_textured_quad_indices[0], GL_STATIC_DRAW);
  }
  if (mesh_interior_edge_indices.size() > 0) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_interior_edge_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * mesh_interior_edge_indices.size(),
		 &mesh_interior_edge_indices[0], GL_STATIC_DRAW);
  }
  if (mesh_border_edge_indices.size() > 0) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_border_edge_indices_VBO); 
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * mesh_border_edge_indices.size(),
		 &mesh_border_edge_indices[0], GL_STATIC_DRAW);
  }

  // WARNING: this naive VBO implementation only allows a single texture
  int num_textured_materials = 0;
  for (unsigned int mat = 0; mat < mesh->materials.size(); mat++) {
    Material *m = mesh->materials[mat];
    if (m->hasTextureMap()) {
      glBindTexture(GL_TEXTURE_2D,m->getTextureID());
      num_textured_materials++;
    }
  }
  assert (num_textured_materials <= 1);
}


void Radiosity::drawVBOs() {
  // =====================
  // DRAW ALL THE POLYGONS
  if (args->render_mode == RENDER_MATERIALS) {
    glEnable(GL_LIGHTING);
  } else {
    glDisable(GL_LIGHTING);
  }
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.1,4.0);
  int num_faces = mesh->numFaces();
  assert ((int)mesh_quad_indices.size() + (int)mesh_textured_quad_indices.size() == num_faces);

  glBindBuffer(GL_ARRAY_BUFFER, mesh_quad_verts_VBO);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(0));
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(12));
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(24));

  // draw non textured faces
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_quad_indices_VBO);
  glDrawElements(GL_QUADS, 
		 mesh_quad_indices.size()*4,
		 GL_UNSIGNED_INT,
		 BUFFER_OFFSET(0));

  // draw textured faces
  if (args->render_mode == RENDER_MATERIALS) {
    glEnable(GL_TEXTURE_2D);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glTexCoordPointer( 2, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(36));
  }
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_textured_quad_indices_VBO);
  glDrawElements(GL_QUADS, 
		 mesh_textured_quad_indices.size()*4,
		 GL_UNSIGNED_INT,
		 BUFFER_OFFSET(0));
  if (args->render_mode == RENDER_MATERIALS) {
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_2D);
  }

  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  glDisable(GL_POLYGON_OFFSET_FILL);  


  // =====================
  // DRAW WIREFRAME
  if (args->wireframe) {
    glDisable(GL_LIGHTING);
    if (mesh_interior_edge_indices.size() > 0) {
      glLineWidth(1);
      glColor3f(0,0,0);
      glBindBuffer(GL_ARRAY_BUFFER, mesh_quad_verts_VBO);
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(0));
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_interior_edge_indices_VBO);
      glDrawElements(GL_LINES, mesh_interior_edge_indices.size()*2, GL_UNSIGNED_INT, 0);
      glDisableClientState(GL_VERTEX_ARRAY);
    }
    if (mesh_border_edge_indices.size() > 0) {
      glLineWidth(3);
      glColor3f(1,0,0);
      glBindBuffer(GL_ARRAY_BUFFER, mesh_quad_verts_VBO);
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(0));
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_border_edge_indices_VBO);
      glDrawElements(GL_LINES, mesh_border_edge_indices.size()*2, GL_UNSIGNED_INT, 0);
      glDisableClientState(GL_VERTEX_ARRAY);
    }
  }
  HandleGLError(); 
}


void Radiosity::cleanupVBOs() {
  glDeleteBuffers(1, &mesh_quad_verts_VBO);
  glDeleteBuffers(1, &mesh_quad_indices_VBO);
  glDeleteBuffers(1, &mesh_textured_quad_indices_VBO);
  glDeleteBuffers(1, &mesh_interior_edge_indices_VBO);
  glDeleteBuffers(1, &mesh_border_edge_indices_VBO);
}

