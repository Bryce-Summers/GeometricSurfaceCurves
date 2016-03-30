/*
 * Patch Drawer.
 * Written by Bryce Summers.
 */

#include "bezierPatch.h"
#include <cmath>
#include <math.h>

#define PI M_PI

namespace CMU462
{

  // Direct reduction.
  void BezierPatch::loadControlPoints(FaceIter & face)
  {
    loadControlPoints(face->halfedge());
  }

  void BezierPatch::loadControlPoints(HalfedgeIter & edge)
  {
    // Compute the geompetry position patch.
    computeGeometryControlPoints(edge, geometry);

    // As in the Loop-Shaffer paper, the orientation for ev should be
    // along the same edge as the geometry control points.
    // eu is oriented to the previous edge of the face.
    /*      + <------
     *      |       .
     * eu   |      /|\
     *     \|/      |
     *      . ------>
     *      edge = ev
     */
    HalfedgeIter ev = edge;
    HalfedgeIter eu = edge -> next() -> next() -> next();
    
    // Compute the partial u tangent field patch.
    computeTangentControlPoints(ev, tangents_v, false);

    // Compute the partial v tangent field patch.
    computeTangentControlPoints(eu, tangents_u, true);
  }

  // This is a helper function used index the geometry patch.
  inline Vector3D b(int col, int row, std::vector<Vector3D> & cpts)
  {
    return cpts[row*4 + col];
  }

  // Standard unmodified u direction vectors.
  inline Vector3D u(int col, int row, std::vector<Vector3D> & cpts)
  {
    return 3*(b(col + 1, row, cpts) - b(col, row, cpts));
  }

  // Standard unmodified v direction vectors.
  inline Vector3D v(int col, int row, std::vector<Vector3D> & cpts)
  {
    return 3*(b(col, row + 1, cpts) - b(col, row, cpts));
  }

  void BezierPatch::computeTangentControlPoints(HalfedgeIter & edge,
						std::vector<Vector3D> &cpts,
						bool transpose)
  {

    if(edge -> face() -> degree() != 4)
    {
	cerr << "ERROR: bezierPatch::computeTangentControlPoints: face is not a quadrilateral.";
	//exit(-7);
	return;
    }

    // First we need to compute the control points for the input
    // half edge orientation.
    std::vector<Vector3D> g;
    computeGeometryControlPoints(edge, g);

    /**********************
     *   V02  e2    V32
     *    +<---------+
     *    |          .
     *    |         /|\
     *e3  |          |  e1
     *   \|/         |
     *    .          |
     *    +--------->+
     *   V00  e0    V30
     *********************/

    // Compute the 4 edges of this face.
    HalfedgeIter e0 = edge;
    HalfedgeIter e1 = e0 -> next();
    HalfedgeIter e2 = e1 -> next();
    HalfedgeIter e3 = e2 -> next();
    
    // Because of symmetry we will assume that cpts is the partial_v patch.

    // -- Let us list out the 12 control points that we wish to compute.
    Vector3D V02, V12, V22, V32;
    Vector3D V01, V11, V21, V31;
    Vector3D V00, V10, V20, V30;

    // Set the middle row tangent vectors, because they don't need any modification.
    // Vij = 4(G(i, j + 1) - G(i, j)).
    V01 = v(0, 1, g);
    V11 = v(1, 1, g);
    V21 = v(2, 1, g);
    V31 = v(3, 1, g);

    // I am assuming that we want the input points to point straight up,
    // instead of right. (+v, instead of +u).
    
    /*
    V30 = computeCornerTangentPoint(e1);
    V32 = computeCornerTangentPoint(e2 -> twin() -> next());
    V02 = computeCornerTangentPoint(e2 -> twin() -> next() -> next()
				                 -> next() -> twin());
    V00 = computeCornerTangentPoint(e3 -> twin());
    */

    V30 = computeCornerTangentPoint(e1);
    V32 = computeCornerTangentPoint(e2);
    V02 = computeCornerTangentPoint(e3);
    V00 = computeCornerTangentPoint(e0);

    // Right facing, like in the paper.
    /*
    V30 = computeCornerTangentPoint(e1 -> twin() -> next());
    V32 = computeCornerTangentPoint(e2 -> twin() -> next() ->
				    twin() -> next());
    V02 = computeCornerTangentPoint(e2 -> twin());
    V00 = computeCornerTangentPoint(e0);
    */

    /*
    V02 *= -1;
    V30 *= -1;
    */
    

    // The paper says that I need to flip two of these corner points,
    // but I think that I have properly oriented them.


    // Now all that remains is for us to fix V10, V20, V12, and V22.

    // The valence of the lower left corner.
    const int n00 = e0 -> vertex() -> degree();

    // The valence of the lower right corner.
    const int n01 = e1 -> vertex() -> degree();

    // The valence of the upper left corner.
    const int n10 = e3 -> vertex() -> degree();

    // The valence of the upper right corner.
    const int n11 = e2 -> vertex() -> degree();

    const int c00 = cos(2*PI/n00);
    const int c01 = cos(2*PI/n01);
    const int c10 = cos(2*PI/n10);
    const int c11 = cos(2*PI/n11);

    /**********************
     *   V02  e2    V32
     *    +<---------+
     *    |  v12 v22 .
     *    |         /|\
     *e3  |          |  e1
     *   \|/ v10 v20 |
     *    .          |
     *    +--------->+
     *   V00  e0    V30
     *********************/
    
    V10 = 1.0/3*(2*c00*u(1, 0, g) - c01*u(0,0, g)) +
          3*(b(1, 1, g) - b(1, 0, g));
    V20 = 1.0/3*(  c00*u(2, 0, g) - 2*c01*u(1,0, g)) +
          3*(b(2, 1, g) - b(2, 0, g));

    // V12 and V22 are implemented in an analagous way.
    V12 = 1.0/3*(2*c10*u(1, 3, g) -   c11*u(0,3, g)) + 3*(b(1, 3, g) - b(1, 2, g));
    V22 = 1.0/3*(  c10*u(2, 3, g) - 2*c11*u(1,3, g)) + 3*(b(2, 3, g) - b(2, 2, g));

    // Use default geometry patch tangent patches, but they are not continuous.
    
    V00 = v(0, 0, g);
    V01 = v(0, 1, g);
    V02 = v(0, 2, g);
    V10 = v(1, 0, g);
    V11 = v(1, 1, g);
    V12 = v(1, 2, g);
    V20 = v(2, 0, g);
    V21 = v(2, 1, g);
    V22 = v(2, 2, g);
    V30 = v(3, 0, g);
    V31 = v(3, 1, g);
    V32 = v(3, 2, g);
    

    if(!transpose)
    {
      // -- Let us list out the 12 control points that we wish to compute.
      // V02, V12, V22, V32;
      // V01, V11, V21, V31;
      // V00, V10, V20, V30;

      // Row major order.
      cpts.push_back(V00);
      cpts.push_back(V10);
      cpts.push_back(V20);
      cpts.push_back(V30);
      
      cpts.push_back(V01);
      cpts.push_back(V11);
      cpts.push_back(V21);
      cpts.push_back(V31);

      cpts.push_back(V02);
      cpts.push_back(V12);
      cpts.push_back(V22);
      cpts.push_back(V32);


    }
    else // Transposed.
    {
      // V02, V12, V22, V32;
      // V01, V11, V21, V31;
      // V00, V10, V20, V30;

      // Column major order.
      cpts.push_back(V00);
      cpts.push_back(V01);
      cpts.push_back(V02);

      cpts.push_back(V10);
      cpts.push_back(V11);
      cpts.push_back(V12);

      cpts.push_back(V20);
      cpts.push_back(V21);
      cpts.push_back(V22);

      cpts.push_back(V30);
      cpts.push_back(V31);
      cpts.push_back(V32);
    }

    return;
  }

  // Returns the Tangent Control point associated with corner point at the beginning
  // of the given half edge, where the half edge should be pointing in the positive
  // v direction.
  // NOTE: if something is wrong, try rotating the input edge which is used.
  Vector3D BezierPatch::computeCornerTangentPoint(HalfedgeIter & up_half_edge)
  {
    int n = up_half_edge -> vertex() -> degree();

    HalfedgeIter edge = up_half_edge -> next();

    // The output will be the weighted sum of vectors along a mask.
    Vector3D out(0, 0, 0);

    // We will precompute some constant terms in the weighting expression.
    double cos_pi_n = cos(PI/n);
    double radical = sqrt(4 + cos_pi_n*cos_pi_n);
    const double beta_constant  = 1.0/(n*radical);
    const double alpha_constant = 1.0/n + cos_pi_n*beta_constant;

    for(int i = 0; i < n; i++)
    {
      // ASSUMPTION: edge should be originating from the alpha_i location
      //             and going to the beta_i location.

      // Compute the mask coeficients.
      double alpha_i = alpha_constant*cos(2*PI*i/n);
      double beta_i  = beta_constant *cos((2*i + 1)*PI/n);
      out += alpha_i*edge -> vertex() -> position;

      // Transition to the other half edge in the face which doesn't touch the
      // corner point, which is the center of the one-ring.
      // After this our half edge is originating from beta_i and is going towards
      // alpha_i + 1.
      edge = edge -> next();

      out += beta_i*edge -> vertex() -> position;

      // Now transition to the next face ready for the next loop iteration!
      edge = edge -> next() -> twin() -> next();
    }

    return out;
  }

  /* In this function, we compute a set of 16 control points P cooresponding to
   * a bicubic spline from a quadrilateral face in a Catmull-Clark control mesh.
   *
   * We compute the points in P according to the paper:
   *   'Approximating Catmull-Clark Subdivision Surfaces with Bicubic Patches.'
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
  void BezierPatch::computeGeometryControlPoints(HalfedgeIter & edge,
					 std::vector<Vector3D> &cpts)
  {
    // We assume that the face is a quadrilateral with 4 vertices.

    if(edge -> face() -> degree() != 4)
    {
	cerr << "ERROR: bezierPatch::computeControlPoints: face is not a quadrilateral.";
	//exit(-7);
	return;
    }

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

    // -- Declare the 16 geometry patch control points.
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

    // First of all, ensure the geometry patch is clear.
    cpts.clear();
    
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

  void BezierPatch::addCornerPointNeighborhood(Vector3D & cpt, HalfedgeIter edge)
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

  HalfedgeIter BezierPatch::getEdgeNeighborPositions(HalfedgeIter & edge,
						     Vector3D & e1,
						     Vector3D & e2)
  {
	// First go to the twin facing quadrilateral and walk to the other side.
	HalfedgeIter e = edge->twin()->next()->next();
	e2 = e->vertex() -> position;
	e1 = e->next() -> vertex() -> position;

	return e;
  }

  
  Vector3D BezierPatch::evaluateGeometryPatch(double u, double v,
				 int partial_u, int partial_v)
  {
    Vector3D sum(0, 0, 0);

    // Control points are interpolated.
    for(int i = 0; i < 4; i++)// row
    for(int j = 0; j < 4; j++)// column
    {
      sum += geometry[i + 4*j] * // 4 is num columns.
	     Bernstein(u, i, partial_u) *
	     Bernstein(v, j, partial_v);
    }

    return sum;
  }

  Vector3D BezierPatch::evaluateNormal(double u, double v,
				       int partial_u, int partial_v)
  {
    Vector3D du = evaluateTangentUPatch(u, v, 0, 0);
    Vector3D dv = evaluateTangentVPatch(u, v, 0, 0);

    // Cross the two partials to get the normal.
    return cross(du, dv).unit();
  }

  // Note: The partial results in the geometry differentiated by u 1 + partial
  // since this is a tangent patch.
  Vector3D BezierPatch::evaluateTangentUPatch(double u, double v,
				 int partial_u, int partial_v)
  {
    // Evaluate the du patch.
    Vector3D du(0, 0, 0);

    // 4 by 3 (Row major)
    for(int j = 0; j < 4; j++)// row
    for(int i = 0; i < 3; i++)// column
    {
      // Since there are only 3 horizontal control points,
      // we use the quadratic bernstein basis functions in the u direction.
      du += tangents_u[i + 3*j]  * // 3 = num columns.
            Bernstein_2(u, i, partial_u) *
            Bernstein  (v, j, partial_v);
    }

    return du;
  }
  
  // Note: The partial results in the geometry differentiated by u 1 + partial
  // since this is a tangent patch.
  Vector3D BezierPatch::evaluateTangentVPatch(double u, double v,
				 int partial_u, int partial_v)
  {
    // Evaluate the dv patch.
    Vector3D dv(0, 0, 0);

    // 3 by 4 (Row major)
    for(int j = 0; j < 3; j++)// row
    for(int i = 0; i < 4; i++)// column
    {
      dv += tangents_v[i + 4*j]  * // 4 = num columns
            Bernstein  (u, i, partial_u) *
            Bernstein_2(v, j, partial_v);
    }

    return dv;
  }

  void BezierPatch::ejectGeometryPatchControlPoints(std::vector<Vector3D> & out)
  {
    for(auto iter = geometry.begin(); iter != geometry.end(); ++iter)
    {
      out.push_back(*iter);
    }
  }
  

}// CMU462 namespace.
