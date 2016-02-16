/*
 * Patch Drawer.
 * Written by Bryce Summers.
 */

#include "patchDrawer.h"

namespace CMU462
{


  /* In this function, we compute a set of 16 control points P cooresponding to
   * a bicubic spline from a quadrilateral face in a Catmull-Clark control mesh.
   *
   * We compute the points in P according to the paper:
   *   'Approximating Catmull-Clark Subdivision Surfaces with Bicubic Patches.'
   *
   * Thus far, we are only computing geometry patches and are not using an
   * sophisticated schemes to make the tangents smooth up to and including
   * edge cases.
   *
   * 1. We need to find and label every vertex on the face.
   *    We will need to deduce relevant vertices in the one neighborhood of
   *    the face on the fly to account for potential extraordinary
   *    vertices of degree != 4.
   *    one neighborhood of the face.
   *
   * 2. Compute the 16 bicubic patch control points.
   *    - 3 types: 4 interior points, 8 edge points, and 4 corner points.
   *
   * Note: We are computing 16 local control points using a the region of
   *       16 control mesh points.
   * In this way we are deriving one set of 16 points in terms of a larger
   * region of 16 points. This is only the case for quad faces without an
   * extraordinary vertex though.
   *
   * 3. We then hand these points down to the bicubic patch drawing
   *    function.
   */
  void PatchDrawer::drawCatmullClarkQuadPatch(FaceIter & face)
  {
    	std::vector<Vector3D> cpts;
	computeControlPoints(face, cpts);
	drawBicubicPatch(cpts);
  }

  void PatchDrawer::computeControlPoints(FaceIter & face,
					 std::vector<Vector3D> & cpts)
  {
	// We assume that the face has 4 vertices.

	if(face->degree() != 4)
	{
	  cerr << "ERROR: PatchDrawer::drawCatmullClarkQuadPatch: face is not a quadrilateral.";
	  //exit(-7);
	  return;
	}

	HalfedgeIter edge = face->halfedge();

	// First label the 4 points on the face and the 4 halfedges.
	VertexIter v1, v2, v3, v4;
	HalfedgeIter e1, e2, e3, e4;

	e1 = edge;
	e2 = e1->next();
	e3 = e2->next();
	e4 = e3->next();

	v1 = e1->vertex();
	v2 = e2->vertex();
	v3 = e3->vertex();
	v4 = e4->vertex();

	Vector3D & p1 = v1->position;
	Vector3D & p2 = v2->position;
	Vector3D & p3 = v3->position;
	Vector3D & p4 = v4->position;

	int degree_1 = v1->degree();
	int degree_2 = v2->degree();
	int degree_3 = v3->degree();
	int degree_4 = v4->degree();

	// -- Declare the 16 control points.
	Vector3D b11, b12, b21, b22, // Interior Points.
                 b10, b20, b01, b02, // Edge Points.
                 b13, b23, b31, b32,
                 b00, b03, b30, b33; // Corner Points.

	// -- Interior Points.
	// First off we will compute the interior points,
	// since we don't have to find anything.

	// 2 --- 1  Numbers represent weights in a normalized mask.
	// |     |
	// | x   |  <-- 4 Interior points, in for orientations, biased towards the
	// n --- 2      4 corners.


	// Note all cases should be cyclic permutations of all indices modulo 4.
	b11 = (p1*degree_1 + p2*2 + p3*1 + p4*2)/(degree_1 + 5);
	b12 = (p2*degree_2 + p3*2 + p4*1 + p1*2)/(degree_2 + 5);
	b22 = (p3*degree_3 + p4*2 + p1*1 + p2*2)/(degree_3 + 5);
	b21 = (p4*degree_4 + p1*2 + p2*1 + p3*2)/(degree_4 + 5);

	// -- Now we will find all of the edge points.

	// 2  --- 1  Numbers represent weights in a normalized mask.
	// |      |
	// 2n x-y 4 <-- Edge points x and y along a given edge.
	// |      |
	// 2  --- 1

	Vector3D neighbor1, neighbor2;

	// Half Edges on the twin faces.
	// These will be used later when deriving the corner points.
	HalfedgeIter e5, e6, e7, e8;

	// Procceed in 4 orientations, which are cyclic permutations
	// of each other.

	// Edge points along edge 4.
	e8 = getEdgeNeighborPositions(e4, neighbor1, neighbor2);

	b10 = 2*p2 + 1*p3 +
	      2*degree_1*p1 + 4*p4 +
	      2*neighbor1 + 1*neighbor2;
	b10 /= (10 + 2*degree_1);

	b20 = 1*p2 +        2*p3 +
	      4*p1 +        2*degree_4*p4 +
	      1*neighbor1 + 2*neighbor2;
	b20 /= (10 + 2*degree_4);

	// Edge points along edge 1.
	e5 = getEdgeNeighborPositions(e1, neighbor1, neighbor2);

	b02 = 2*p3 + 1*p4 +
	      2*degree_2*p2 + 4*p1 +
	      2*neighbor1 + 1*neighbor2;
	b02 /= (10 + 2*degree_2);

	b01 = 1*p3 +        2*p4 +
	      4*p2 +        2*degree_1*p1 +
	      1*neighbor1 + 2*neighbor2;
	b01 /= (10 + 2*degree_1);


	// Edge points along edge 2.
	e6 = getEdgeNeighborPositions(e2, neighbor1, neighbor2);

	b23 = 2*p4 + 1*p1 +
	      2*degree_3*p3 + 4*p2 +
	      2*neighbor1 + 1*neighbor2;
	b23 /= (10 + 2*degree_3);

	b13 = 1*p4 +        2*p1 +
	      4*p3 +        2*degree_2*p2 +
	      1*neighbor1 + 2*neighbor2;
	b13 /= (10 + 2*degree_2);

	// Edge points along edge 3.
	e7 = getEdgeNeighborPositions(e3, neighbor1, neighbor2);

	b31 = 2*p1 + 1*p2 +
	      2*degree_4*p4 + 4*p3 +
	      2*neighbor1 + 1*neighbor2;
	b31 /= (10 + 2*degree_4);

	b32 = 1*p1 +        2*p2 +
	      4*p4 +        2*degree_3*p3 +
	      1*neighbor1 + 2*neighbor2;
	b32 /= (10 + 2*degree_3);



	// -- Now we need to compute the corner points and we will then be all done.

	// First we will find all of the far halfedges on the 4 twin quadrilateral faces.

	//             v2    v3
	//       1 --- 4 --- 1  The number represent weights in a normalized mask.
	//    e5 |     |     |
	//       4 --\ |x    | <-- Corner point for each of 4 oriantations.
	//             n^2-- 4 v4  (v1 is on the inside with the n^2)
	// neig1 4 --/ |     |
	//       |     |     |
	// neig2 1 --- 4 --- 1
	//                e8
	//
	// A corner is defined as a sort of centroid of its neighboring vertices.

	// Lets start with the v1 corner point (b00) and then procceed via cyclic
	// permutations of the code.

	addCornerPointNeighborhood(b00, e4);
	b00 += degree_1*degree_1*p1;
	b00 /= (degree_1*(5 + degree_1));

	// Now lets move on to the v2 corner point (b03);
	addCornerPointNeighborhood(b03, e1);
	b03 += degree_2*degree_2*p2;
	b03 /= (degree_2*(5 + degree_2));

	// Now lets move on to the v3 corner point (b33);
	addCornerPointNeighborhood(b33, e2);
	b33 += degree_3*degree_3*p3;
	b33 /= (degree_3*(5 + degree_3));

	// Now lets move on to the v4 corner point (b30);
	addCornerPointNeighborhood(b30, e3);
	b30 += degree_4*degree_4*p4;
	b30 /= (degree_4*(5 + degree_4));

	// Step 3: Now that we have all of the control points,
	// we package them up and pass them to the Bicubic Patch drawer.

	cpts.push_back(b00);
	cpts.push_back(b01);
	cpts.push_back(b02);
	cpts.push_back(b03);

	cpts.push_back(b10);
	cpts.push_back(b11);
	cpts.push_back(b12);
	cpts.push_back(b13);

	cpts.push_back(b20);
	cpts.push_back(b21);
	cpts.push_back(b22);
	cpts.push_back(b23);

	cpts.push_back(b30);
	cpts.push_back(b31);
	cpts.push_back(b32);
	cpts.push_back(b33);
  }

  void PatchDrawer::addCornerPointNeighborhood(Vector3D & cpt, HalfedgeIter edge)
  {
	HalfedgeIter h0 = edge -> twin();
	HalfedgeIter h  = h0;

	do
	{
	  HalfedgeIter temp = h -> next();
	  cpt += 4*temp->vertex() -> position + temp -> next() -> vertex() -> position;
	  h = h -> twin() -> next();
	}while(h != h0);
  }

  HalfedgeIter PatchDrawer::getEdgeNeighborPositions(HalfedgeIter & edge,
						     Vector3D & e1,
						     Vector3D & e2)
  {
	// First go to the twin facing quadrilateral and walk to the other side.
	HalfedgeIter e = edge->twin()->next()->next();
	e2 = e->vertex() -> position;
	e1 = e->next() -> vertex() -> position;

	return e;
  }

  /*
   * Uses combinations of the control points and the bezier basic functions to parametrically describe a bicubic patch. We then draw it as quads.
   *
   */
  void PatchDrawer::drawBicubicPatch(std::vector<Vector3D> & control_points)
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

	  Vector3D p1 = evaluatePatch(control_points, u1, v1);
	  Vector3D p2 = evaluatePatch(control_points, u1, v2);
	  Vector3D p3 = evaluatePatch(control_points, u2, v2);
	  Vector3D p4 = evaluatePatch(control_points, u2, v1);

	  Vector3D normal = cross(p2 - p1, p1 - p4).unit();
	  glNormal3dv( &normal.x );

	  // Draw a quadrilateral.
	  glBegin(GL_POLYGON);

	  glVertex3dv( &p1.x );
	  glVertex3dv( &p2.x );
	  glVertex3dv( &p3.x );
	  glVertex3dv( &p4.x );

	  glEnd();

	}

  }

  // A parametric function that describes the surface.
  Vector3D evaluatePatch(
	   std::vector<Vector3D> & control_points,
	   double u, double v,
	   int partial_u, int partial_v)
  {
	Vector3D sum(0,0,0);

	for(int i = 0; i < 4; i++)
	for(int j = 0; j < 4; j++)
	{
	  sum += control_points[i + 4*j] *
	         Bernstein(u, i, partial_u) *
	         Bernstein(v, j, partial_v);
	}

	return sum;
  }

  double Bernstein(double u, int index, int derivative)
  {
    switch(derivative)
    {
      // Normal, undifferentiared value.
      case 0:
	switch(index)
	{
	  case 0: return pow(1.0 - u, 3);
	  case 1: return 3*u*pow(1.0 - u, 2);
	  case 2: return 3*u*u*(1.0 - u);
	  case 3: return pow(u, 3);
	}

      // First derivatives.
      case 1:
	switch(index)
	{
	  case 0: return -3*pow(1.0 - u, 2);
	  case 1: return (u - 1.0)*(9.0*u - 3.0);
	  case 2: return (6.0 - 9.0*u)*u;
	  case 3: return 3.0*u*u;
	}
	
      case 2:
	switch(index)
	{
	  case 0: return 6.0 - 6.0*u;
	  case 1: return 18.0*u - 12.0;
	  case 2: return 6.0 - 18.0*u;
	  case 3: return 6.0*u;
	}
	
      case 3:
	switch(index)
	{
	  case 0: return -6.0;
	  case 1: return 18.0;
	  case 2: return -18.0;
	  case 3: return 6.0;
	}
    }

    if(derivative > 3)
    {
      return 0.0;
    }

    if(index < 0 || index > 3)
    {
	cerr << "ERROR: PatchDrawer::Bernstein index value not in range";
	return 0.0;
    }

    if(derivative < 0)
    {
	cerr << "ERROR: PatchDrawer::Bernstein derivative value should not be negative.";
	return 0.0;
    }

    cerr << "ERROR: PatchDrawer::Bernstein Something went wrong.";
    return 0.0;
    
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