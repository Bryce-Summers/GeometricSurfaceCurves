#ifndef CURVE_TRACER_H
#define CURVE_TRACER_H

// Input Control Mesh faces.
#include "halfEdgeMesh.h"// FIXME: Determine whether we actually need this.

// Face to Bicubic Patch transformation.
#include "patchDrawer.h" // FIXME: Determine whether we actually need this.

// Spline Representation of Traced Curves.
#include "pointCurve.h"

#include "critical_points.h"

#include "UnionFind/UF_Serial.h"

/*
 * Contour Tracer,
 * Written by Bryce Summers on 1/29/2016. (Started)
 *
 * This class contains functionality for specifying implicitly defined curves.
 *  as well as their 1st and second order partial derivatives.
 *
 * While the patch drawer class specifies the continuous function f defining the
 * mesh surface S, which is parameterized with u and v values,
 * this class specifies particular equations that define points on curves on
 * this surface.
 */

using namespace std;

namespace CMU462
{

  // FIXME : Evenetually we are going to abstract this to a subclass.
  // We should make a base class that has all of the high level algorithms for tracing gradients and paths,
  // given a specification for the function and its 1st and 2nd order derivatives in u and v.
  class Curve_Silhouette
  {
    public:

    // The critical points will be accurate to within a manhattan metric distance
    // of tolerance*2;
    const double TOLERANCE          = .00000001;
    const double CURVE_TRACING_STEP = .000001; // This needs to be larger than the Tolerance.
    // This is the coeficient for tracing along the gradient of the
    // magnitude of the gradient..
    const int GRAD_SQR_COEF = 10000;

    // This is the coeficient for tracing along the gradient.
    const int GRAD_COEF = 10000;
    
    // The viewing eye ray.
    Vector3D E;
    
    Curve_Silhouette(Vector3D e){this->E = e;};
    ~Curve_Silhouette(){}


    // Finds exactly one point on each silhouette curve.
    // num_critical_points -> IN
    //  - specifies the number of critical points on the mesh that generated
    //    these critical points, All critical points should have unique indices
    //    between [0 and num_critical_points) stored within themselves.
    //
    // edges -> IN a set of point curves which connect and associate critical points
    //          and contain the locations of origination points.
    void find_unique_silhouette_points(
			       int num_critical_points,
			       std::vector<PointCurve> & edges,
			       std::vector<Critical_Point> & silhouette_points);
    
    // Constructs the Morse Smale Complex
    // critical points and edges will be assumed to be empty.
    // mesh   -> IN  control cage with connectivity information.
    // cp_out -> OUT All critical points found for the silhouette function.
    // edges  -> OUT All of the morse smale complex lines, represented by
    //               gradient following integral lines.
    void morse_smale_complex(HalfedgeMesh& mesh,
			     std::vector<Critical_Point> & cp_out,
			     std::vector<PointCurve> & edges);

    // -- Critical Point finding functionality.
    // Populates the given list with critical points found on the silhouette curve
    // function defined on the given mesh.
    // Pass NULL values if you don't want either type of critical point.
    // Populates mesh faces with the location of any critical point located on them.
    void findCriticalPoints (HalfedgeMesh& mesh,
			     std::vector<Critical_Point> * saddle_points,
			     std::vector<Critical_Point> * all_points
			     );

    // Finds all critical points.
    void findCriticalPoints(HalfedgeMesh& mesh,
			    std::vector<Critical_Point> & points)
    {
      findCriticalPoints(mesh, NULL, &points);
    }
    

  private:

    // Checks the given u and v coordinates for being out of the [0, 1) x [0, 1)
    // If they are then it transitions to the correct new face and
    // localizes the u and v value to the new face.
    // Returns true iff a transition occured.
    // The face iter is updated as necessary.
    // This function works hand and hand with face to patch control point
    // translation function. They need to be consistent.
    // ENSURES: Applies up to 1 transition. This function should be called repeatedly
    // until it no longer provides a valid transition. For instance this function should
    // be called twice when u < 0 and v < 0.
    bool transitionIfNeccessary(FaceIter & current_face, double & u, double & v);

    inline int determineOrientation(HalfedgeIter edge);
    
    // Traces a curve from the given critical point in the given direction to a MIN or MAX critical point.
    // If the initial directional gradient is positive, then this function will trace a postive  going gradient.
    // If the initial directional gradient is negative, then this function will trace a negative going gradient.
    // For creating a MS complex, these should originate form saddle points.
    // returns false if an ending critical point was not found, for instance in meshes with boundaries.
    // ASSUMPTION: if we choose 2 perpendicular directions, then they will go in opposite gradient directions.
    // REQUIRES: curve should be empty coming in.
    bool trace_gradient(Critical_Point cp,    // Starting point.
			double du, double dv, // Direction.
			PointCurve & curve    // OUT: Output point curve that is populated.
			);
    
    // Follow the gradient of f^2 down to a critical point.
    // returns true if a critical point was found in [0,1] x [0,1],
    // otherwise returns false if the search goes out of bounds.
    bool search_for_critical_point(std::vector<Vector3D> & control_points,
				   double u, double v,
				   Critical_Point & p);

    // Classifies critical points based on the signs of the leading principle minors,
    // of the hessian.
    Critical_Point_Type classify_point(std::vector<Vector3D> & control_points,
				       double u, double v);
    // -- Function specifications.
    
    // position on surface or a partial derivative in the given combination of
    // u and v partials.
    // P(c,u,v, 0, 0) 
    Vector3D P(std::vector<Vector3D> & control_points,
	       double u, double v,
	       int partial_u = 0, int partial_v = 0);

    // Normal vector to surface at point P(u,v);
    Vector3D N(std::vector<Vector3D> & control_points, double u, double v);

    // F = N dot E, this defines the silhouette when it equals 0.
    double F(std::vector<Vector3D> & control_points, double u, double v);

    // -- 1st order partial derivatives.
    double F_u(std::vector<Vector3D> & control_points, double u, double v);
    double F_v(std::vector<Vector3D> & control_points, double u, double v);

    // -- 2nd order partial derivatives.
    double F_uu(std::vector<Vector3D> & control_points, double u, double v);
    // Note: F_uv = F_vu, because of 2nd partial calculus.
    double F_uv(std::vector<Vector3D> & control_points, double u, double v);
    double F_vv(std::vector<Vector3D> & control_points, double u, double v);

    // Returns the gradient of f at the given uv coordinates.
    Vector2D grad_f(std::vector<Vector3D> & control_points, double u, double v);

    // The gradient of the magnitude of the gradient of F.
    // Using this result as a descent direction should yield all critical points.
    Vector2D grad_f_sqr(std::vector<Vector3D> & control_points,
			double u, double v);
  };
  
}

#endif // CURVE_TRACER_H
