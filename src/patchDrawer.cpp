/*
 * Patch Drawer.
 * Written by Bryce Summers.
 */

#include "patchDrawer.h"

namespace CMU462
{

  void PatchDrawer::drawCatmullClarkQuadPatch(FaceIter & face)
  {
    BezierPatch patch;
    patch.loadControlPoints(face);
    drawBicubicPatch(patch);
  }

  void PatchDrawer::drawTangentPatches(FaceIter & face)
  {
    BezierPatch patch;
    patch.loadControlPoints(face);
    drawUVTangents(patch);
  }

  /*
   * Draws an evaluatable patch as an array of quadrilaterals by evaluating the
   * the 'push forward' of a u-v space grid as a surface in space.
   */
  void PatchDrawer::drawBicubicPatch(BezierPatch & patch)
  {

        // The prescision by which we are drawing the Bicubic Patch.
        int n = 4;

	// Rasterize the Bicubic patch via the parametric function.
	for(int i = 0; i < n; i++)
	for(int j = 0; j < n; j++)
	{
	  double u1 = i*1.0/n;
	  double v1 = j*1.0/n;

  	  double u2 = (i+1)*1.0/n;
	  double v2 = (j+1)*1.0/n;

	  // Compute the Geometry points.
	  Vector3D p1 = patch.evaluateGeometryPatch(u1, v1);
	  Vector3D p2 = patch.evaluateGeometryPatch(u1, v2);
	  Vector3D p3 = patch.evaluateGeometryPatch(u2, v2);
	  Vector3D p4 = patch.evaluateGeometryPatch(u2, v1);

	  // Compute the normals to the surface using the tangent patches.
	  Vector3D n1 = patch.evaluateNormal(u1, v1);
	  Vector3D n2 = patch.evaluateNormal(u1, v2);
	  Vector3D n3 = patch.evaluateNormal(u2, v2);
	  Vector3D n4 = patch.evaluateNormal(u2, v1);
  

	  // FIXME: directly evaluate the normal form the patch.
	  Vector3D normal = cross(p2 - p1, p1 - p4).unit();

	  // Per Polygon normal.
	  //glNormal3dv( &normal.x );
	  //glNormal3dv( &n1.x );

	  // Draw a quadrilateral.
	  // Using the positions from the geometry patch.
	  // and using the tangents
	  
	  glBegin(GL_POLYGON);

	  glNormal3dv( &n1.x);
	  glVertex3dv( &p1.x );

	  glNormal3dv( &n2.x);
	  glVertex3dv( &p2.x );

	  glNormal3dv( &n3.x);
	  glVertex3dv( &p3.x );

	  glNormal3dv( &n4.x);
	  glVertex3dv( &p4.x );

	  glEnd();

	}
  }

  void PatchDrawer::drawUVTangents(BezierPatch & patch)
  {
    std::vector<Vector3D> g;
    patch.ejectGeometryPatchControlPoints(g);

    double N = 4;
    
    for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
    {
      if(i == 1 || i == 2 || j == 1 || j == 2)
      {
	continue;
      }
      
      Vector3D pos = g[j*4 + i];

      double u = i / (N - 1.0);
      double v = j / (N - 1.0);

      //      cout << u << ", " << v << endl;
      
      Vector3D du = patch.evaluateTangentUPatch(u, v, 0, 0).unit();
      Vector3D dv = patch.evaluateTangentVPatch(u, v, 0, 0).unit();
      Vector3D posU = pos + du/10.0;
      Vector3D posV = pos + dv/10.0;

      // Draw U tangent vectors.
      /*
      glBegin(GL_LINES);
      glVertex3dv( &pos.x  );
      glVertex3dv( &posU.x );
      glEnd();

      // Draw V tangent vectors.
      glBegin(GL_LINES);
      glVertex3dv( &pos.x  );
      glVertex3dv( &posV.x );
      glEnd();
      */

      // Draw Tangent plane triangles.

      Vector3D normal = cross(du, dv).unit();
      glNormal3dv( &normal.x );

      glBegin(GL_POLYGON);

      glVertex3dv( &pos.x  );
      glVertex3dv( &posU.x );
      glVertex3dv( &posV.x );

      glEnd();
    }
  }

  // Rasterizes the control mesh's face with a simple draw polygon call.
  void PatchDrawer::drawControlFace(FaceIter & face)
  {
	// Start specifying the polygon.
	glBegin(GL_POLYGON);

	// Set the normal of this face.
	Vector3D normal = face->normal();
	glNormal3dv( &normal.x );

	// iterate over this polygon's vertices
	HalfedgeIter h = face->halfedge();
	do
	{
	  // Draw this vertex.
	  Vector3D position = h->vertex()->position;
	  glVertex3dv( &position.x );

	  // go to the next vertex in this polygon
	  h = h->next();

	  // end of iteration over polygon vertices
	} while( h != face->halfedge() );
	// Finish drawing the polygon.
	glEnd();
  }
}
