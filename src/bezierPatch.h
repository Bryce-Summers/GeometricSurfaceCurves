#ifndef BEZIER_PATCH_H
#define BEZIER_PATCH_H

#include "halfEdgeMesh.h"   // Connectivity data structure.
#include "PatchFunctions.h" // Math: Bernstein polynomials.
#include "CMU462/CMU462.h"  // Standard 462 Vectors, etc.

/*
 * The Bezier Patch class controls the representation and evaluation of bezier
 * surfaces that approximate the positions and tangents of points on a limit
 * Catmull - Clark subdivision surface.
 *
 * This class represents the surface associated with a particular quadrilateral
 * face in a halfedgemesh by a geometry patch and a pair of tangent patches,
 * as described in the paper "Approximating Catmull-Clark Subdivision Surfaces
 * with Bicubic Patches" by Charels Loop and Scott Schaefer.
 *
 * The geometry patches have (4 by 4) control points,
 * The tangent patches du is (4 by 3) and dv is (3 by 4).
 *
 */

using namespace std;

namespace CMU462
{

  class BezierPatch
  {
  private:

    // Here is the labeling scheme for the 16 control points for the geometry
    // patch.
    // B00, B01, B02, B03,
    // B10, B11, B12, B13,
    // B20, B21, B22, B23,
    // B30, B31, B32, B33,

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

    // All of the control points arrays are packed in row major order.
    std::vector<Vector3D> geometry; // 4 by 4.

    // -> u03, -> u13, -> u23
    // -> u02, -> u12, -> u22
    // -> u01, -> u11, -> u21
    // -> u00, -> u10, -> u20
    std::vector<Vector3D> tangents_u; // 4 by 3 (Row major)

    //  .    .    .    .
    // /|\  /|\  /|\  /|\
    //  |    |    |    |
    // v02, v12, v22, v32
    //  .    .    .    .
    // /|\  /|\  /|\  /|\
    //  |    |    |    |
    // v01, v11, v21, v31
    //  .    .    .    .
    // /|\  /|\  /|\  /|\
    //  |    |    |    |
    // v00, v10, v20, v30
    std::vector<Vector3D> tangents_v; // 3 by 4 (Row major order.)
    
  public:
    
    BezierPatch (){};
    ~BezierPatch(){}

    // This function derives the 16 control points for the
    // geometry patch associated with the given quadrilateral face.
    // This function also computes the 12 control points for both
    // the u and v tangent patches.
    // ENSURES: clears the input array of control points.
    void loadControlPoints(FaceIter & face);

    // This function works the same as the face version,
    // except that it orients the control points to the given edge,
    // instead of the canonical halfedge for a face.
    // ENSURES: clears the input array beforehand.
    //          There will be exactly 16 points afterwards.
    void loadControlPoints(HalfedgeIter & edge);

    // Evaluates a Bicubic patch based on a standard set of 16 control points.
    // The (0, 0) u,v location is oriented with the first control point,
    // the collumns vary over u and the rows vary over v.
    //
    // Half Edge and u/v coordinate orientation information.
    // (0,0) ----> (1,0) (u, v)
    //   .     0     |
    //  /|\          |
    //   |  3     1  |
    //   |          \|/
    //   |     2     .
    // (0,1) <---- (1, 1)
    // Evaluates the 'push forward' of the geometry patch
    // at the given u, v coordinates.
    // if partial_u and are used
    // REQUIRES: The geometry patch and the tangent patches need to
    //           already be previously loaded.
    Vector3D evaluateGeometryPatch(double u, double v,
				   int partial_u = 0, int partial_v = 0);

    // Note: The partial results in the geometry differentiated by u 1 + partial
    // since this is a tangent patch.
    Vector3D evaluateTangentUPatch(double u, double v,
				   int partial_u = 0, int partial_v = 0);

    // Note: The partial results in the geometry differentiated by u 1 + partial
    // since this is a tangent patch.
    Vector3D evaluateTangentVPatch(double u, double v,
				   int partial_u = 0, int partial_v = 0);
    
    // Evaluates the 'push forward' of the tangent patches.
    // The output vector is normalized.
    Vector3D evaluateNormal(double u, double v,
			    int partial_u = 0, int partial_v = 0);

    // Copies the geometry patch's contents into the given out vector.
    // This function does not clear the out vector.
    void ejectGeometryPatchControlPoints(std::vector<Vector3D> & out);

  private:

    // Populate the input array with the control points associated with the
    // given half edge.
    // This function may be used to directly determine the geometry patch and
    // may also be used for the twin edges needed for the tangent patches.
    // REQUIRES: The HE mesh is entirely quadrilateral,
    //        extraordinary vertices are handled correctly.
    // ENSURES: the given array is cleared internally and will only contain the 16 relevant points in row major order.
    // Evaluates a Bicubic patch based on a standard set of 16 control points.
    // The (0, 0) u,v location is oriented with the first control point,
    // the collumns vary over u and the rows vary over v.
    //
    // Half Edge and u/v coordinate orientation information.
    // B00   edge.
    // (0,0) ----> (1,0) (u, v)
    //   .     0     |
    //  /|\          |
    //   |  3     1  |
    //   |          \|/
    //   |     2     .
    // (0,1) <---- (1, 1)
    void computeGeometryControlPoints(HalfedgeIter & edge,
				      std::vector<Vector3D> &cpts);

    // populates the given tangent patch with appropiate control vectors.
    // Because of symmetry, this function merely populates the Dv (3 by 4) patch.
    // If you desire to populate the Du patch, call this function on the next
    // halfedge and then transpose the output.
    // ASSUMPTION: all faces should be quadrilaterals.
    void computeTangentControlPoints(HalfedgeIter & edge,
				     std::vector<Vector3D> &cpts,
				     bool transpose);

    // -- These are some helper function used in the computation of control
    //    point positions.

    // Returns the positions of both vertices on the twin face of the
    // given edge that do not touch the given edge.
    // These are needed for deriving edge points.
    // Returns an iter to the half edge pointing from vertex2 to vertex1.
    HalfedgeIter getEdgeNeighborPositions(HalfedgeIter & edge,
					  Vector3D & v1,
					  Vector3D & v2);

    void addCornerPointNeighborhood(Vector3D & cpt, HalfedgeIter h0);


    // Tangent helper functions.
    Vector3D computeCornerTangentPoint(HalfedgeIter & up_half_edge);
  };
}

#endif // BEZIER_PATCH_H
