#ifndef PATCH_DRAWER_H
#define PATCH_DRAWER_H

#include "halfEdgeMesh.h"   // Connectivity data structure.
#include "PatchFunctions.h" // Math: Berstein polynomials.
#include "bezierPatch.h" // Bezier Surface representation.

/*
 * Patch Drawer,
 * Written by Bryce Summers on 11/29/2015.
 *
 * This Class handles the drawing Bicubic patches.
 * The bicubic patches will be an approximation of Catmull-Clark subdivision
 * surfaces on quadrilateral meshes.
 *
 * - HalfEdge quad mesh to control points.
 * - Use the control points to derive a parametric function f(x, y)
 * - Use the parametric function to draw the mesh.
 *
 * We keep the representation of the patches encapsulated inside of the bezier
 * patch class so that this class may focus on just drawing patches to the
 * screen from evaluations of the positions and normals of the bezier patch.
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

	 void drawTangentPatches(FaceIter & face);
	 
	 // Takes 16 control points and rasterizes a bicubic patch from a set
	 // of evenly spaced quadrilaterals.
	 void drawBicubicPatch(BezierPatch & patch);

	 // Draws lines from the geometry control points in the direction
	 // of the u and v tangent vectors at those points.
	 void drawUVTangents(BezierPatch & patch);

	 // Simply draws the face as a polygon.
	 // This is equivilant to drawing the control mesh.
	 void drawControlFace(FaceIter & face);

   };
}

#endif // STUDENT_CODE_H
