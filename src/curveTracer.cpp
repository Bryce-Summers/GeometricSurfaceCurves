/*
 * CurveTracer.
 * Written by Bryce Summers.
 *
 * Contains the specification of the calculus for Silhouette curves.
 *
 * Also provides routines for calculating the morse smale complex of the silhouette curve
 * defining function.
 */

#include "curveTracer.h"
#include "style.h"
#include "svg_exporter.h"
#include <sstream>

namespace CMU462
{

  void Curve_Silhouette::find_silhouette_points(HalfedgeMesh& mesh,
			      std::vector<Critical_Point *> & pts_out)
  {
    // For every edge,
    // 1. Compute the coeficients of the silhouette function polynomial along that
    //    edge.
    // 2. Points are at locations along the edge cooresponding to roots of this
    // polynomial.

    for( EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++ )
    {
      findRoots_F(e, pts_out);
    }
  }
  
  void Curve_Silhouette::trace_zero_level_sets(
       std::vector<Critical_Point *> & silhouette_points,
       std::vector<PointCurve * > & edges)
  {
    for(auto iter = silhouette_points.begin();
	iter != silhouette_points.end();
	++iter)
    {
      Critical_Point * point = *iter;
      
      // Don't revisit critical points.
      // -- IGNORE Critical points that are not of the origination type.
      if(point -> visited || point -> type != ORIGINATION)
      {
	cout << "Ignoring an origination point!\n";
	continue;
      }

      point -> visited = true;

      PointCurve * curve = new PointCurve();
      edges.push_back(curve);

      // Create a reference to the allocated curve object.
      trace_zero_curve(point, curve);
      //return;// REMOVE me!

    }

  }

  // This is very similar to the gradient tracer.
  // It traces the 0 level set perpendicular to the gradient.
  bool Curve_Silhouette::trace_zero_curve(Critical_Point * start,
					  PointCurve * curve)
  {
    // Compute the canonical control points for the initial face.
    FaceIter current_face = start -> face;
    patch.loadControlPoints(current_face);

    // Compute a scalar 'dir' that represents whether we wish to go up or down.
    // The gradient will be multiplied by -1 if we wish to head towards a minnimum.
    //Vector2D grad = grad_f(control_points, start.u, start.v);

    // The curve starts and ends at the given start point.
    curve -> p1 = start;
    curve -> p2 = start;

    // We will now trace the gradient in a circle.
    double u = start -> u;
    double v = start -> v;

    // Make sure we start on the silhouette curve.
    //    moveOntoLevel(u, v, current_face, control_points);

    // Store the locations of the starting point.
    double u_original = u;
    double v_original = v;
    bool stop = false;
    
    // NOTE: every point that is added, needs to also have a tangent added.
    Vector3D p0 = P(patch, u, v);
    curve -> addPoint(p0);
    curve -> addTangent(movePerpGrad (u, v, current_face, patch, stop));
    moveOntoLevel(u, v, current_face, patch);

    int count = 0;
    const int min_count = 10;
    
    while(count < min_count || (!stop && count < 10000))
	   // || // FIXME: Improve these bounds.
	    //(P(control_points, u, v) - p0).norm() > GRAD_PERP_DIST)))
      // I would like to use the STOP variables for this.
    {
      count++;
      Vector3D point = P(patch, u, v);
      curve -> addPoint(point);

      Vector3D tangent = movePerpGrad (u, v, current_face, patch, stop);
      curve -> addTangent(tangent);

      moveOntoLevel(u, v, current_face, patch);
    }

    // Add one last point at the end to elliminate the gap.
    curve -> addPoint(P(patch, u, v));
    Vector3D tangent = movePerpGrad (u, v, current_face, patch, stop);
    curve -> addTangent(tangent);
    
    cout << count << endl;

    return true;
  }


  // CURRENTLY ONLY SUPPORTS MOVING to 0 level.
  void Curve_Silhouette::moveOntoLevel(double & u, double & v,
				       FaceIter & current_face,
				       BezierPatch & patch,
				       double level)
  {
    int count = 0;

    double val = abs(F(patch, u, v));
    double val_old = val + 1.0;
    
    while(count < 1000 && val > LEVEL_TOLERANCE && val < val_old)
    {
      count++;
      val_old = val;
      
      // Downhill direction line search.
      Vector2D neg_grad = -grad_f_2(patch, u, v);

      moveOntoLevel_step(u, v, // start
			 patch, // geometry.
			 neg_grad // direction
			 );

      performTransitions(current_face, u, v, patch);

      val = abs(F(patch, u, v));
    }

    //cout << "Exit loop." << endl;
  }

  // Implementation From Keenan Crane's notes on line search using the
  // Armijo and Wolfe conditions.
  // Note: it is possible that we will need to be wary of patch transions in
  // these calculations.
  // The candidate location will be returned in the out_location vector.
  void Curve_Silhouette::moveOntoLevel_step(
            double & u, double & v, // Start.
	    BezierPatch & patch,
	    Vector2D dir // Line search direction.
			  )
  {
    
    // The function we are evaluating is f_2.

    double alpha = 0.0;
    double beta  =  -1;
    double t = 1.0; // Step size.

    // 0 < c1 < c2 < 1
    const double c1 = .8;
    const double c2 = .9;

    // The directional derivative in the gradient's direction.
    const double Du_x0 = dir.norm();
    const Vector2D x0 = Vector2D(u, v);
    const double f_x0 = F_sqr(patch, x0.x, x0.y);

    // Loop until we have found a viable candidate.
    for(int i = 0; i < 100; i++)
    {
      Vector2D x1 = x0 + t*dir;

      // Evaluate the objective.
      double f_x1 = F_sqr(patch, x1.x, x1.y);

      // The Armijo condition states that we wish to find a point that
      // strictly decreases the objective.
      bool Armijo = f_x1 <= f_x0 + c1*t*Du_x0;

      // The Wolfe condition states that we want to decrease the magnitude of
      // the gradient.
      bool Wolfe = false; // We only compute the Wolfe condition as necessary.

      double Du_x1;
      
      // Evaluate the gradient.
      if(Armijo)
      {
	//dot(grad_f_2(patch, x1.x, x1.y), dir);
	Du_x1 = grad_f_2(patch, x1.x, x1.y).norm();
	Wolfe = abs(Du_x1) <= c2*abs(Du_x0);
      }

      if(!Armijo)
      {
	beta = t;
      }
      else if(!Wolfe)
      {
	alpha = t;
      }
      else // We have found a candidate point.
      {
	u = x1.x; // Update the incoming locations.
	v = x1.y;

	/*
	cout << "f^2(x1) = " << f_x1  << endl;
	cout << "Du(x1) = " << Du_x1 << endl;
	cout << "f(x1) = " << F(control_points, u, v) << endl;
	*/

	return;
      }

      // While the Armijo condition is still met, we expand the search interval.
      if(beta < 0.0)
      {
	t = 2*alpha;
      }
      else // Otherwise we essentially perform binary search on the interval.
	   // (Alpha, beta)
      {
	t = (alpha + beta)/2.0;
      }
    }

    //    return;
    std::cerr << "Curve Tracer Error, we should never come here.\n";
  }

  bool Curve_Silhouette::performTransitions(
       FaceIter & current_face,
       double & u,
       double & v,
       BezierPatch & patch)
  {
    bool transition = false;
    
    // NOTE: may mutate current_face, u, and v.
    bool stop = false;
    bool previously_visited_edge = false;
    while(transitionIfNeccessary(current_face, u, v, stop))
    {
      transition = true;
      previously_visited_edge |= stop;
    }

    // If a transition occured, then we need to recalculate the control points.
    // FIXME : It might be a good idea to store critical points on bicubic patches.
    if(transition)
    {
      patch.loadControlPoints(current_face);
    }
    
    return previously_visited_edge;
  }

  Vector3D Curve_Silhouette::movePerpGrad(double & u, double & v,
				      FaceIter & current_face,
				      BezierPatch & patch,
				      bool & stop
					  )
  {

    Vector2D grad = grad_f(patch, u, v);

    grad = grad.unit();
    
    u += -GRAD_PERP_DIST*grad.y;
    v +=  GRAD_PERP_DIST*grad.x;
          
    stop = performTransitions(current_face, u, v, patch);

    // Returns the tangent vector in the direction of the perpendicular
    // gradient movement.
    Vector3D world_direction = P(patch, u, v, 1, 0)*(-grad.y) +
                               P(patch, u, v, 1, 0)*( grad.x);
    return world_direction;
  }

  void Curve_Silhouette::find_unique_silhouette_points(int num_critical_points,
					  std::vector<PointCurve*> & edges,
			  std::vector<Critical_Point *> & silhouette_points)
  {
    // Now use Union find on the list of curves to derive one unique
    // origination point for every level set curve that we are interested in.
    UF_Serial UF(num_critical_points);

    // Union all edges that do not contain an origination point.
    // This decomposes the MS complex into connected components fully above the level set
    // and those that are completely below the level set.
    // The level set is the set of points satisfying the implicit equation of the final curves that we wish to extract.
    for(auto iter = edges.begin(); iter != edges.end(); ++iter)
    {
      if((*iter) -> has_crossing_point == false)
      {
	Critical_Point * cp1 = (*iter) -> p1;
	Critical_Point * cp2 = (*iter) -> p2;

	// Union these critical point sets.
	// The entire set will be above or below the level set.
	UF.op_union(cp1 -> index, cp2 -> index);
      }
    }

    // FIXME : If we want we could store those sets with crossing points, then
    // simply iterate over them and perhaps store a count of how many distinct
    // connected components remain and stop once there is only 1.
    
    // For all MS edges,
    // Add an origination point if the edge's two sets are not yet connected.
    for(auto iter = edges.begin(); iter != edges.end(); ++iter)
    {
      if((*iter) -> has_crossing_point == true)
      {
	Critical_Point * cp1 = (*iter) -> p1;
	Critical_Point * cp2 = (*iter) -> p2;

	if(!UF.connected(cp1 -> index, cp2 -> index))
	{
	  // Push
	  silhouette_points.push_back((*iter) -> level_set_crossing_point);

	  // Union these critical point sets.
	  UF.op_union(cp1 -> index, cp2 -> index);
	}
      }
    }
    
    return;
  }

  void Curve_Silhouette::morse_smale_complex(HalfedgeMesh& mesh,
			   std::vector<Critical_Point *> & cp_all, // This goes out.
			   std::vector<PointCurve*> & edges
					     )
  {
    // To build the MS complex, we will trace 4 gradient curves from
    // each saddle point to two min's and max's each.
    std::vector<Critical_Point *> saddle_points;

    // Find all of the critical points,
    // populate the saddle point array.
    // populate the output critical points array for whomever out there might want them,
    // for instance for visualization purposes.
    findCriticalPoints(mesh, &saddle_points, &cp_all);


    // There will be 4*|saddle_points| edges in the MS complex,
    // because every saddle point extends 4 edges.
    
    // We will trace 4 curves for each saddle point.
    // Therefore we tell the edges vector how much space we will need.
    edges.reserve(saddle_points.size()*4);
    
    // For every saddle point, we will trace 4 gradient curves, which are MS complex edges.
    for(auto iter = saddle_points.begin(); iter != saddle_points.end(); ++iter)
    {
      Critical_Point * saddle = *iter;
      double amount = CURVE_TRACING_STEP; // FIXME, choose a more appropiate value.
      cout << "Tracing a Gradient" << endl;
      
      // ASSUMPTION: if we choose 2 perpendicular directions,
      //             then they will go in opposite gradient directions.
      // Here we are tracing gradients in 4 directions out from the saddle point.
      for(int sx = -1; sx <= 1; sx += 2)
      for(int sy = -1; sy <= 1; sy += 2)
      {

	PointCurve * curve = new PointCurve();
	// These act as edges in the MS Complex.
	edges.push_back(curve);

	// Create a reference to the allocated curve object.
 	trace_gradient(saddle, amount*sx, amount*sy, curve);
      }
    }

    return;
  }

  // REQUIRES: curve should be empty coming in.
  // Returns false if no ending critical point was found.
  bool Curve_Silhouette::trace_gradient(
		Critical_Point * cp,    // Starting point.
		double du, double dv, // Direction.
		PointCurve * curve    // IN: point curve that is populated.
					)
  {

    // Compute the canonical control points for the initial face.
    FaceIter current_face = cp -> face;
    patch.loadControlPoints(current_face);

    
    // Compute the a scalar 'dir' that represents whether we wish to go up or down.
    // The gradient will be multiplied by -1 if we wish to head towards a minnimum.
    Vector2D grad = grad_f(patch, cp -> u + du, cp -> v + dv);
    double   dot  = grad.x*du + grad.y*dv; // Positive dot product then <du,dv>
    double   dir  = dot > 0 ? 1.0 : -1.0;  // heads uphill.   

    // The curve starts at the original given critical point.
    curve -> p1 = cp;
    curve -> addPoint(cp -> location);

    // Now follow the gradient until we find a max/min.
    // Add these points to the curves as we go along.
    // We transition from one bicubic patch to another
    // when we go out of [0,1] x [0,1] bounds.
    double u = cp -> u;
    double v = cp -> v;
    double u_new = u + du;
    double v_new = v + dv;

    // Silhouette curve function value
    double f = F(patch, u, v);
    // Curve heading in the direction of a Silhouette.
    // aka heading towards 0 level set.
    bool check_crossing = f*dir < 0;
    
    while( abs(u_new - u) > TOLERANCE || abs(v_new - v) > TOLERANCE )
    {      
      u = u_new;
      v = v_new;

      Vector3D location = P(patch, u, v);
      curve -> addPoint(location);

      // --  Handle the logic of switching between bicubic patches.
      performTransitions(current_face, u, v, patch);

      // Add the silhouette crossing point to the curve if found.
      // crossing found if the curve is no longer headed towards the 0 level set.
      if(check_crossing && F(patch, u, v)*dir > 0)
      {
	cout << "Adding a crossing Point" << f << endl;
	
	Critical_Point * crossing = new Critical_Point();
	crossing -> location = location;
	crossing -> u = u;
	crossing -> v = v;
	crossing -> type = ORIGINATION;
	crossing -> index = -1;// FIXME Use something useful if needed.
	crossing -> face = current_face;
	crossing -> f_val = 0; // Pretend, close to 0.
	
	curve -> level_set_crossing_point = crossing;
	curve -> has_crossing_point = true;

	// Don't check anymore, because all integral curves are monotonic
	// increasing / decreasing. There will only be 1.
	check_crossing = false;
      }

      
      grad = dir*grad_f(patch, u, v);
      
      u_new = u + GRAD_COEF*grad.x;
      v_new = v + GRAD_COEF*grad.y;
    }

    if(current_face -> critical_point == NULL)
    {
      cerr << "ERROR: CurveTracer: Function converged," <<
	"but no critical point was previously stored on this face." << endl;
      return false;
    }

    curve -> p2 = current_face -> critical_point;
    curve -> addPoint(curve -> p2 -> location);
    return true;
  }

  // Populates the given list with critical points found on the silhouette curve
  // function defined on the given mesh.
  void Curve_Silhouette::findCriticalPoints (HalfedgeMesh& mesh,
					     std::vector<Critical_Point*> * saddle_points,
					     std::vector<Critical_Point*> * all_points)
  {
    int index = 0;
    
    for(FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++)
    {
      patch.loadControlPoints(f);
    
      Critical_Point * cp = new Critical_Point();

      // FIXME : Our strategy for finding points still needs some thought.
      // IDEA1 : Perhaps critical points can only occur when curvature is present
      // in the control mesh.
      // IDEA2: Maybe search during subdivison.
      // IDEAL : analytically find the answer.
      // Maybe solve 1D problems and then trace along those routes.

      // Right now we just start at each corner of the patch and see if they find a
      // critical point inside of the patch.
      bool found_point = false;
      for(int u = 0; u <= 1; u++)
      {
      for(int v = 0; v <= 1; v++)
      {
	// Find Critical Points..
	if(search_for_critical_point(patch, u, v, cp))
	{
	  cp -> type = classify_point(patch, u, v);
	  cp -> face = f;// The cp needs to know its face for direct access
	              // to the HalfEdge Mesh connectivity structure.
	  cp -> index = index; // Every cp gets a unique index.
	  index++;
	  cp -> f_val = F(patch, u, v); // Store the function evaluation at this point.	  
	  
	  // Keep special track of saddle points,
	  // because we will use them in the building of a morse smale complex.
	  if(saddle_points != NULL && cp -> type == SADDLE)
	  {
	    saddle_points -> push_back(cp);
	  }
	  
	  // Record all points in the larger output array.
	  if(all_points != NULL)
	  {
	    all_points -> push_back(cp);
	  }

	  // Tell the face what its critical point is.
	  f -> critical_point     = cp;

	  // FIXME : We are currently only looking for one critical point per patch.
	  // This may not be an assumption that holds water.
	  found_point = true;
	  break;
	}

	// Mark this face as not having a valid critical point thus far.
	f -> critical_point = NULL;
      }
      // Double break;
      if(found_point){break;}
      }      
    }
  }
  
  // Follow the gradient of down to a critical point.
  // returns true if a critical point was found in [0,1] x [0,1],
  // otherwise returns false if the tracing goes out of bounds.
  bool Curve_Silhouette::search_for_critical_point(
			      BezierPatch & patch,
			      double u, double v,
			      Critical_Point * p)
  {
    Vector2D grad;
    do
    {
      grad = GRAD_SQR_COEF*grad_mag_f(patch, u, v);

      u -= grad.x;
      v -= grad.y;

      if(u < 0 || u > 1 || v < 0 || v > 1)
      {
	return false;
      }
    }while(abs(grad.x) > TOLERANCE || abs(grad.y) > TOLERANCE);

    // Populate the found critical point.
    p -> u = u;
    p -> v = v;
    p -> location = P(patch, u, v);

    return true;
  }

  Critical_Point_Type Curve_Silhouette::classify_point(
		      BezierPatch & patch,
		      double u, double v)
  {
    // Compute entires in the 2 by 2 symmetric hessian matrix.
    double A = F_uu(patch, u, v);
    double B = F_uv(patch, u, v);
    double D = F_vv(patch, u, v);

    // Compute the coefficients of the characteristic quadratic
    // formula Q for the hessian matrix.
    double a = 1, b = -(A + D), c = A*D - B*B;

    // Compute the Discriminant of Q.
    double disc = sqrt(b*b - 4*a*c);
    double two_a = 2*a;
        
    // Compute the eigenvalues, which are the roots of Q.
    // Using the quaratic formula.
    double val_1 = (-b + disc)/two_a;
    double val_2 = (-b - disc)/two_a;
    
    // First principle minors, correctly characterize positive definitness,
    // but not others.
    /*
    double pm_1 = A;
    double pm_2 = A*D - B*B;
    */

    //    cout << "Eigen Values: " << val_1 << ", " << val_2 << endl;
    
    /*
    if(abs(dot(N(control_points, u, v).unit(), E)) < .9)
    {
	return SADDLE;1
    }
    */

    if(abs(val_1) < TOLERANCE || abs(val_2) < TOLERANCE)
    {
      return SADDLE;
    }
    
    if(val_1 > 0 && val_2 > 0)
    {
      return MIN;
    }

    if(val_1 < 0 && val_2 < 0)
    {
      return MAX;
    }

    return SADDLE;
  };
    
  // position on surface.
  Vector3D Curve_Silhouette::P(BezierPatch & patch,
			       double u, double v, int partial_u, int partial_v)
  {
    /*
	Vector3D sum(0,0,0);

	for(int i = 0; i < 4; i++)
	for(int j = 0; j < 4; j++)
	{
	  sum += control_points[i + 4*j] * Bernstein(u, i) * Bernstein(v, j);
	}

	return sum;
    */
    // Evaluates a (u,v)th order partial derivative of the patch position.

    /*
    if(partial_u >= 1 && partial_v == 0)
    {
      return patch.evaluateTangentUPatch(u, v, partial_u - 1, 0);
    }
    else if(partial_u == 0 && partial_v >= 1)
    {
      return patch.evaluateTangentUPatch(u, v, 0, partial_v - 1);
    }
    */
    
    return patch.evaluateGeometryPatch(u, v, partial_u, partial_v);
  }

  // Normal vector to surface at point P(u,v);
  Vector3D Curve_Silhouette::N(BezierPatch & patch,
			       double u, double v)
  {
    // FIXME: I should instead use the tangent patch.

    //return patch.evaluateNormal(u, v);
    
    return cross(P(patch, u, v, 1, 0), // P_u
		 P(patch, u, v, 0, 1));// P_v
    
  }

  // F = N dot E, this defines the silhouette when it equals 0.
  double Curve_Silhouette::F(BezierPatch & patch,
			     double u, double v)
  {
    return dot(N(patch, u, v), E);
  }

  // Returns the square of F.
  double Curve_Silhouette::F_sqr(BezierPatch & patch,
				 double u, double v)
  {
    double f    = F(patch, u, v);
    return 1000*f*f;
  }

  // -- 1st order partial derivatives.
  double Curve_Silhouette::F_u(BezierPatch & patch,
			       double u, double v)
  {
    return dot(E, cross(P(patch, u, v, 2, 0),
			P(patch, u, v, 0, 1)) +
	          cross(P(patch, u, v, 1, 0),
			P(patch, u, v, 1, 1)));
  }
  
  double Curve_Silhouette::F_v(BezierPatch & patch,
			       double u, double v)
  {
    return dot(E, cross(P(patch, u, v, 1, 1),
			P(patch, u, v, 0, 1)) +
	          cross(P(patch, u, v, 1, 0),
			P(patch, u, v, 0, 2)));
  }

  // -- 2nd order partial derivatives.
  double Curve_Silhouette::F_uu(BezierPatch & patch,
				double u, double v)
  {
    return dot(E, cross(P(patch, u, v, 3, 0),
			P(patch, u, v, 0, 1)) +
	       
	        2*cross(P(patch, u, v, 2, 0),
			P(patch, u, v, 1, 1)) +

	          cross(P(patch, u, v, 1, 0),
			P(patch, u, v, 2, 1))
	       );
  }    
  
  // Note: F_uv = F_vu, because of 2nd partial calculus.
  double Curve_Silhouette::F_uv(BezierPatch & patch,
				double u, double v)
  {
    return dot(E, cross(P(patch, u, v, 2, 1),
			P(patch, u, v, 0, 1)) +
	       
	          cross(P(patch, u, v, 2, 0),
			P(patch, u, v, 0, 2)) +

	          cross(P(patch, u, v, 1, 0),
			P(patch, u, v, 1, 2))
	       );
  }
  
  double Curve_Silhouette::F_vv(BezierPatch & patch,
				double u, double v)
  {
    return dot(E, cross(P(patch, u, v, 1, 2),
			P(patch, u, v, 0, 1)) +
	       
	        2*cross(P(patch, u, v, 1, 1),
			P(patch, u, v, 0, 2)) +

	        2*cross(P(patch, u, v, 1, 0),
			P(patch, u, v, 0, 3))
	       );
  }


  // Returns the gradient of f at the given uv coordinates.
  Vector2D Curve_Silhouette::grad_f(BezierPatch & patch,
				    double u, double v)
  {
    return Vector2D(
		    F_u(patch, u, v),
		    F_v(patch, u, v)
		   );
  }

  // The gradient of the magnitude of the gradient of f.
  Vector2D Curve_Silhouette::grad_mag_f(BezierPatch & patch,
					double u, double v)
  {
    double f_u  = F_u(patch, u, v);
    double f_v  = F_v(patch, u, v);
    double f_uu = F_uu(patch, u, v);
    double f_uv = F_uv(patch, u, v);
    double f_vv = F_vv(patch, u, v);

    // <
    // d/du (f_u^2 + f_v^2),
    // d/dv ...
    // >
    
    return Vector2D(
		    2*f_u*f_uu + 2*f_v*f_uv,
		    2*f_u*f_uv + 2*f_v*f_vv
		   );
  }
  
  Vector2D Curve_Silhouette::grad_f_2(BezierPatch & patch,
				      double u, double v)
  {
    double f    = F(patch, u, v);
    double f_u  = F_u(patch, u, v);
    double f_v  = F_v(patch, u, v);

    return Vector2D(
		    2*f*f_u,
		    2*f*f_v
		    );
  }

  bool Curve_Silhouette::visible(BezierPatch & patch,
				 double u,
				 double v)
  {
    return F(patch, u, v) > 0;
  }
  
  void Curve_Silhouette::findRoots_F(EdgeIter edge,
				     std::vector<Critical_Point *> & out_points)
  {

    // Compute the canonical halfedge associated with this edge.
    HalfedgeIter half_edge = edge -> halfedge();
    
    // We will solve in 'edge' oriented space when v = 0.
    // In other words we will compute control patches oriented with edge at
    // the top and then apply some math to compute the polynomial
    // along this border.
    //    v = 0
    //u=0       u=1
    //    v = 1
    // We will then compute the locations of roots in edge relative space and
    // convert them to canonical u,v point locations.

    // We will compute the polynomial as follows, assuming v = 0.0:
    // f(u) = Eye dot (P_u x P_v), evaluated as polynomials, not numbers.

    // Compute the edge oriented control points, NOT CANONICAL.
    std::vector<Vector3D> control_points;
    patch.loadControlPoints(half_edge);
    patch.ejectGeometryPatchControlPoints(control_points);
    
    PolynomialVector3D P_u; // Starts at 0.0;

    for(int i = 0; i <= 3; i++)
    {
      P_u += Bernstein_poly(i, 1)*control_points[i];
    }

    PolynomialVector3D P_v;

    for(int j = 0; j <= 3; j++)
    {
      P_v += 3*Bernstein_poly(j, 0)*(control_points[4 + j] - control_points[j]);
    }

    // Combine the eye vector direction and the 2 tangent directions in u and v.
    // using the triple product E dot (P_u x P_v).
    // to compute the explicit polynomial for the evaluation of the silhouette
    // function at any point along 'edge'.
    Polynomial silhouette_func = dot(E, cross(P_u, P_v));

    std::vector<double> roots;
    silhouette_func.computeRealRoots(roots, 0.0, 1.0, TOLERANCE);

    // Now that we have the edge local positions, we need to convert them into
    // useful critical points in canonical surface u,v space.
    (edge -> intersects).clear();
    for(auto iter = roots.begin(); iter != roots.end(); ++iter)
    {
      double u, v;
      double root = *iter;
      orient_coordinate_along_edge(half_edge, root, u, v);

      // A new output root point.
      Critical_Point * out = new Critical_Point();
      out_points.push_back(out);
      int index = out_points.size() - 1;

      out -> location = P(patch, root, 0); // Used in point drawing.
      out -> face = half_edge -> face();
      out -> u = u;// These are in canonical patch orientation.
      out -> v = v;
      out -> type = ORIGINATION;
      out -> index = index;// FIXME: We might need more useful index value.
      out -> f_val = 0.0; // By definition silhouette points have f values of 0.0;
      out -> visited = false;
      
      (edge -> intersects).push_back(out);
    }
    
  }

  // -- Critical Point and curve drawing functions.
  void CurveTracer::drawCriticalPoints()
  {
    //glDisable(GL_DEPTH_TEST);
    glPointSize(20.0);

    // Choose a draw color based on what type of point we are drawing.
    // Then draw the point.
    for(std::vector<Critical_Point*>::iterator iter = points.begin();
	iter != points.end(); ++iter)
    {
        Critical_Point * point = *iter;
	Vector3D p = point -> location;
	//	 cout << p << endl;
	switch(point -> type)
	{
	  case MIN: glColor3f(0.0, 0.0, 1.0);
	    break;
	  case MAX: glColor3f(1.0, 0.0, 0.0);
	    break;
	  case SADDLE: glColor3f(0.0, 1.0, 0.0);//green
	    break;
	  case ORIGINATION: glColor3f(1.0, 1.0, 1.0);
	    break;
	  case TERMINATION:
	    break;
	  default: break;
	}

	glBegin( GL_POINTS );
	glVertex3d( p.x, p.y, p.z );
	glEnd();
    }

    //glEnable( GL_DEPTH_TEST );
  }

  void CurveTracer::drawCurves()
  {
    // -- Draw all each one of the morse smale edges.
    for(auto iter = curves.begin(); iter != curves.end(); ++iter)
    {
      PointCurve * curve = *iter;
      curve -> draw();
    }
  }


  void CurveTracer::trace_all_silhouettes(
       HalfedgeMesh& mesh, Vector3D eye_direction)
  {
    Curve_Silhouette tracer(eye_direction);
    tracer.find_silhouette_points(mesh, points);
    tracer.trace_zero_level_sets (points, curves);
  }

  
  void CurveTracer::trace_parametric_curve(HalfedgeIter edge, Vector3D eye_direction)
  {
    
    // Setup the tracing parameters.
    std::vector<Vector3D> control_points;
    Curve_Silhouette curve_tracer(eye_direction); // Used for visibility checks.
    double u, v, du, dv;
    bool visible = true;

    // Allocate a new Curve Object to store the results in.
    PointCurve * curve = new PointCurve();
    curves.push_back(curve);
    curve -> visible = visible;

    
    HalfedgeIter e0 = edge;
    
    do
    {
      determineLocationAndHeading(edge, u, v, du, dv);
      du /= STEPS_PER_PATCH;
      dv /= STEPS_PER_PATCH;
      
      // Compute the current local Bicubic Patch.
      BezierPatch patch;
      patch.loadControlPoints(edge -> face());

      // Add regularly spaced points on the halfedge defined curve.
      for(int i = 0; i < STEPS_PER_PATCH; i++)
      {
	u += du;
	v += dv;

	// If visibility has changed, split the curve.
	if(curve_tracer.visible(patch, u, v) != visible)
	{
	  visible = !visible;
	  curve = new PointCurve();
	  curve -> visible = visible;
	}
	
	curve -> addPoint(patch.evaluateGeometryPatch(u, v));
      }

      edge = edge -> next();

      // Check for curve termination conditions.
      VertexIter vertex = edge->vertex();
      if(vertex -> isBoundary()) {break;}
      if(vertex -> degree() != 4){break;}
      edge = edge -> twin() -> next();

      //cout << &edge << endl;

    }while(edge != e0); 
  }

  void CurveTracer::clearData()
  {
    for(auto iter = points.begin(); iter != points.end(); ++iter)
    {
      delete (*iter);
    }

    for(auto iter = curves.begin(); iter != curves.end(); ++iter)
    {
      delete (*iter);
    }
    
    points.clear();
    curves.clear();
  }

  // Exports all of the data as paths and points in an svg file.
  void CurveTracer::export_to_svg(size_t screen_w, size_t screen_h)
  {
    
    cout << "CurveTracer: Exporting SVG file" << endl;
    
    SVG_Exporter out;

    SVG_Style style;
    style.stroke       = "green";
    style.fill         = "none";
    style.stroke_width = "5";

    std::stringstream w_string;
    std::stringstream h_string;
    w_string << screen_w;
    h_string << screen_h;

    out.beginSVG("curve_svg.svg",
		 w_string.str(), // width.
		 h_string.str(), // height.
		 "0 0 " + w_string.str() + " " + h_string.str(),
		 // Canonical coordinates.
		 style
		 );

    out.addTitle("Silhouette Curves.");
    out.addDescription("This is a cool svg exported from Bryce's curve extraction program!");

    out.loadProjectionMatrix(screen_w, screen_h);
    
    out.beginGroup(style);

    // Put all of the curves in the svg.
    for(auto iter = curves.begin(); iter != curves.end(); ++iter)
    {
      out.g_curve(**iter);
    }
    
    out.endGroup();
    
    out.endSVG();

    cout << "CurveTracer: Done Exporting SVG" << endl;
  }

}
