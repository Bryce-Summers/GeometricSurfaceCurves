#ifndef PATCH_FUNCTIONS
#define PATCH_FUNCTIONS

/* Factored out by Bryce Summers around Midnight on Wednesday, Febuary 24.
 *
 * This class includes foundational mathematical functions for specifying bicubic
 * patches.
 *
 * This class also stores some useful functions for transversing a set of
 * bicubic patches,
 * whose connectivity is represented by a half edge mesh.
 *
 *
 */

#include "critical_points.h"
#include "halfEdgeMesh.h"
#include "CMU462/CMU462.h" // Standard 462 Vectors, etc.

using namespace std;

namespace CMU462
{
  class Bernstein_3
  {
  public:
    static bool initialised;// Static bool should be initialized to 0 = FALSE.
    static Polynomial B[4][4];
    static void init();
  };

  // The definition of Berstein basis functions at namespace level,
  // so that they can be accessed everywhere. Not just in the patch drawer.
  // Returns the interpolated value for the indexed
  // berstein basis function for a bicubic patch.
  // Also can return the evaluation of the derivatives.
  // Returns B_[index,3] (input)
  // REQUIRES: input in range [0, 1].
  // REQUIRES: index in range [0, 3].
  // REQUIRES: derivative >= 0. 0 --> non diferentiated function evaluation.
  
  double Bernstein(double input, int index, int derivative = 0);

  // Returns the polynomial representation, instead of the numeric evaluation
  // for the given Bernstein polynomial.
  Polynomial Bernstein_poly(int index, int derivative = 0);

  // Evaluates a Bicubic patch based on a standard set of 16 control points.
  // The (0, 0) u,v location is on oriented with the first control poin,
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
  //
  Vector3D evaluatePatch(std::vector<Vector3D> & control_points,
			 double u, double v,
			 int partial_u = 0, int partial_v = 0);
  
  // Checks the given u and v coordinates for being out of the [0, 1) x [0, 1)
  // If they are then it transitions to the correct new face and
  // localizes the u and v value to the new face.
  // RETURNS true iff a transition occured.
  // The face iter is updated as necessary.
  // This function works hand and hand with face to patch control point
  // translation function. They need to be consistent.
  // ENSURES: Applies up to 1 transition. This function should be called repeatedly
  // until it no longer provides a valid transition. For instance this function should
  // be called twice when u < 0 and v < 0.
  bool transitionIfNeccessary(FaceIter & current_face, double & u, double & v,
			      bool & stop);

  // The input edge is associated with a Face.
  // The face is associated with a canonical and arbitrary starting half edge.
  // this method returns how many times the next() function needs to be called
  // to return the input edge.
  // (0,0) ----> (1,0) (u, v)
  //   .     0     |
  //  /|\          |
  //   |  3     1  |
  //   |          \|/
  //   |     2     .
  // (0,1) <---- (1, 1)
  // The return value is the edge associated number.
  int determineOrientation(HalfedgeIter edge);

  // Given a half edge, populates the data values with the u and v location of the
  // origination vertex associated with the edge and direction values that the half
  // edge is pointing in in terms coordinates on the bicubic patch.
  // ENSURES: naturally u and v will end up with values of 0 or 1 and
  // du and dv will take values of -1, 0, or 1.
  // either du or dv will be 0 and the other one will be -1 or 1.
  // (0,0) ----> (1,0) (u, v)
  //   .     0     |
  //  /|\          |
  //   |  3     1  |
  //   |          \|/
  //   |     2     .
  // (0,1) <---- (1, 1)
  void determineLocationAndHeading(HalfedgeIter edge, // INPUT.
				   double & u,   // OUT
				   double & v,   // OUT
				   double & du,  // OUT
				   double & dv); // OUT

  // REQUIRES: x_in specifies the location along the given edge
  // as a 1 dimensional coordinate from [0, 1],
  // where 0 cooresponds to the origin point for the edge,
  // and 1 cooresponds to the destination or ending point for the edge.
  // returns u and v in the canonical orientation assocated with edge->face();
  void orient_coordinate_along_edge(HalfedgeIter edge,
				    double x_in,
				    double & u,
				    double & v);
  
  // WARNING: These cubic and quadratic roots codes have neither been tested,
  //          nor used.
  
  // Finds the 3 roots of the cubic polynomial with the given coeficients.
  // p(x) = ax^3 + bx^2 + cx + d.
  // RETURNS true iff the roots are real only.
  // REQUIRES: a != 0.
  bool findCubicRoots(double a, double b, double c, double d, // IN
		      double & r1,  // OUT
		      double & r2,  // OUT
		      double & r3); // OUT
    
  bool findCubicRoots(double a, double b, double c, double d, // IN
		      Complex & r1,  // OUT
		      Complex & r2,  // OUT
		      Complex & r3); // OUT

  // Finds the 2 roots of quadratic polynomial with the given coeficients.
  // p(x) = ax^2 + bx^2 + c
  // Returns only Real parts.
  // RETURNS true iff the results are real only.
  bool findQuadraticRoots(double a, double b, double c,
			  double & r1, double & r2);
  
  // Returns full complex roots.
  // RETURNS truee iff the results are real only.
  bool findQuadraticRoots(double a, double b, double c,
			  Complex & r1, Complex & r2);
}


#endif // PATCH_FUNCTIONS
