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

#include "PatchFunctions.h"

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
  
  // Curve data encapsulates the generation, drawing, and exportation of
  // point curves related to visual features of input meshes.
  class CurveTracer
  {
    // ----------------
    // -- Data Storage.
    // ----------------

  private:
    std::vector<Critical_Point> points;
    std::vector<PointCurve> curves;

    PatchDrawer patcher; // Used for computing the control points.

    // The number of points traced per patch for u/v aligned parameter curves.
    int STEPS_PER_PATCH = 10;
    
  public:

    // ----------------
    // -- Data Generation.
    // ----------------

    
    // -- Traces all of the silhouette curves for the current viewpoint.
    void trace_all_silhouettes(HalfedgeMesh& mesh, Vector3D eye_direction);

    // -- Traces the loop including the curve cooresponding to the
    //    given input edge from the control mesh.
    void trace_parametric_curve(HalfedgeIter edge, Vector3D eye_direction);


    // ------------------------------
    // -- Data Drawing to the screen.
    // ------------------------------

    void drawCriticalPoints();
    void drawCurves();

    // Clears all of the data from this object.
    void clearData();

    // Exports all of the data as paths and points in an svg file.
    void export_to_svg(size_t width, size_t height);

  };


  // FIXME : Eventually we are going to abstract this to a subclass.
  // We should make a base class that has all of the high level algorithms for tracing gradients and paths,
  // given a specification for the function and its 1st and 2nd order derivatives in u and v.
  class Curve_Silhouette
  {
    
  private:
    PatchDrawer patcher; // Used for computing the control points.
    
    public:

    // The critical points will be accurate to within a manhattan metric distance
    // of tolerance*2;
    const double TOLERANCE          = .00000001;
    const double LEVEL_TOLERANCE = .00001;
    const double CURVE_TRACING_STEP = .000001; // This needs to be larger than the Tolerance.
    // This is the coeficient for tracing along the gradient of the
    // magnitude of the gradient..
    const int GRAD_SQR_COEF = 100;

    // This is the coeficient for tracing along the gradient.
    const int GRAD_COEF = 100;

    // How much to go in the perpendicular gradient direction.
    const double GRAD_PERP_DIST = .1;
    
    // The viewing eye ray.
    Vector3D E;
    
    Curve_Silhouette(Vector3D e){this->E = e;};
    ~Curve_Silhouette(){}

    // returns true iff the eye can see this location.
    // reduces to a sign check of f.
    bool visible(std::vector<Vector3D> & control_points, double u, double v);

    // Compiles a list of all silhouette points on patch boundaries.
    // GIVEN: A halfedge mesh that specifies the set of patches.
    //
    // RETURNS: pts_out contains a set of points with useful metadata attached.
    //
    void find_silhouette_points(HalfedgeMesh& mesh,
				std::vector<Critical_Point> & pts_out);
    
    // Given a list of silhouette point locations,
    // Creates a Point curve associated with each one of the silhouette points.
    // to edges.
    // REQUIRES: Assumes for now that the mesh is closed and therefore all 0
    //           level sets are closed.
    // REQUIREs: the critical points must be of the origination type.
    //           They also ought to be close to a 0 location for relevance.
    //          -one_patch specifies whether the tracing should span multiple patches
    //           if it is false, then the curve will only be traced until it exits
    //           the patch.
    void trace_zero_level_sets(std::vector<Critical_Point> & silhouette_points,
			       std::vector<PointCurve> & edges);
 
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

    // Returns false if it encounters a boundary.
    // Follows the 0 set perpendicular to the gradient in a loop.
    bool trace_zero_curve(Critical_Point start, PointCurve & curve);

    // -- Moves the given tracing states down to the nearest point at the
    //    Desired Level.
    void moveOntoLevel(double & u, double & v,
		       FaceIter & current_face,
		       std::vector<Vector3D> & control_points,
		       double level = 0.0);

    // Moves the given tracing states perpendicular to the gradient
    // in a consistent direction, (always lefthand or righthand rule).
    // Returns the direction moved in.
    Vector3D movePerpGrad(double & u, double & v,
			  FaceIter & current_face,
			  std::vector<Vector3D> & control_points);

    // Performs as many transitions as are necessary and updates all of the
    // input values according to the ending location on the geometric surface.
    // Marks any critical points encountered as visited.
    void performTransitions(FaceIter & current_face,
			    double & u,
			    double & v,
			    std::vector<Vector3D> & control_points);



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
    Vector2D grad_mag_f(std::vector<Vector3D> & control_points,
			double u, double v);

    // The gradient of the sqr of f, minnimums coorespond to 0 points.
    Vector2D grad_f_2(std::vector<Vector3D> & control_points,
		      double u, double v);

    // Appends all points along the given edge which lie along a silhouette curve.
    // Ignores roots outside of the patch's bounds [0, 1) x [0, 1)
    //
    // edge is the location of intersects.
    // edge -> halfedge() is the canonical oriented side that the roots will be
    //    calculated from along the boundary u in [0, 1).
    // edge -> halfedge() -> face() is the face that defines the canonical patch
    //    coordinates, which will be stored within the Critical points
    //    themselves for location equivalence checks.
    //    
    void findRoots_F(EdgeIter edge,
		     std::vector<Critical_Point> & roots);
      
  };
  
}

#endif // CURVE_TRACER_H
