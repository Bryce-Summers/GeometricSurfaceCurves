#ifndef PATCH_DRAWER_H
#define PATCH_DRAWER_H

#include "halfEdgeMesh.h" // Connectivity data structure.
#include "PatchFunctions.h" // Math: Berstein polynomials, control point evaluation.

/*
 * Patch Drawer,
 * Written by Bryce Summers on 11/29/2015.
 *
 * This Class handles the drawing and representation of Bicubic patches.
 * The bicubic patches will be an approximation of Catmull-Clark subdivision
 * surfaces on quadrilateral meshes.
 *
 * - HalfEdge quad mesh to control points.
 * - Use the control points to derive a parametric function f(x, y)
 * - Use the parametric function to draw the mesh.
 *
 *
 * The mathematical definition for the evaluation of a bernstein basis function.
 * and the evaluation of a bicubic patch from a standard set of 16 control points
 * are defined outside of a patch drawer for use anywhere, whereas the particular
 * methods for converting a half edge mesh into set of canonical control points
 * is kept within a class, because we may want to use different schemes for other
 * types of input in the future, such as non quadrilateral mesh inputs.
 */

using namespace std;

namespace CMU462
{

  class PatchDrawer
  {

      public:

         PatchDrawer (){};
         ~PatchDrawer(){}

	 /* The Face
	  * - should be a quadrilateral.
	  * - represented in a half edge structure representing a quad mesh.
	  * - The quad mesh serves as a control mesh that cooresponds to the
	  *   limit surface that we will approximate using a bicubic patch.
	  */
	 void drawCatmullClarkQuadPatch(FaceIter & face);

	 // Derives the 16 control points associated with the given
	 // quadrilateral face.
	 // pushes them onto the given control_point vector.
	 // OUTPUT: a list of control points in row major order.
	 // Every 4th point starts a new row.
	 // B00, B01, B02, B03,
	 // B10, B11, B12, B13,
	 // B20, B21, B22, B23,
	 // B30, B31, B32, B33,
	 //
	 // Half Edge and u/v coordinate orientation information.
	 // The input face's canonical half edge cooresponds to upper edge in this
	 // diagram between u,v coordinates (0, 0) --> (1, 0).
	 // (0,0) ----> (1,0) (u, v)
	 //   .           |
	 //  /|\          |
	 //   |           |
	 //   |          \|/
	 //   |           .
	 // (0,1) <---- (1, 1)
	 // ENSURES: clears the input array of control points.
	 void computeControlPoints(FaceIter & face,
				   std::vector<Vector3D> & control_point);

	 // Works the same as the face version, except it orients the control points
	 // to the given edge, instead of the canonical halfedge for a face.
	 // ENSURES: clears the input array beforehand.
	 //          There will be exactly 16 points afterwards.
	 void computeControlPoints(HalfedgeIter & edge,
				   std::vector<Vector3D> & control_point);
	 
	 // Takes 16 control points and rasterizes a bicubic patch from a set
	 // of evenly spaced quadrilaterals.
	 void drawBicubicPatch(std::vector<Vector3D> & control_points);

	 // Simply draws the face as a polygon.
	 // This is equivilant to drawing the control mesh.
	 void drawControlFace(FaceIter & face);

      private:

	 // Returns the positions of both vertices on the twin face of the
	 // given edge that do not touch the given edge.
	 // These are needed for deriving edge points.
	 // Returns an iter to the half edge pointing from vertex2 to vertex1.
	 HalfedgeIter getEdgeNeighborPositions(HalfedgeIter & edge,
					       Vector3D & v1,
					       Vector3D & v2);

	 
	 void addCornerPointNeighborhood(Vector3D & cpt, HalfedgeIter h0);

   };
}

#endif // STUDENT_CODE_H
