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

namespace CMU462
{

  void Curve_Silhouette::trace_zero_level_sets(
       std::vector<Critical_Point> & silhouette_points,
       std::vector<PointCurve> & edges)
  {
    for(auto iter = silhouette_points.begin();
	iter != silhouette_points.end();
	++iter)
    {
      // -- IGNORE Critical points that are not of the origination type.
      if(iter->type != ORIGINATION)
      {
	  continue;
      }
      edges.push_back(PointCurve());

      // Create a reference to the allocated curve object.
      PointCurve & curve = edges[edges.size() - 1];
      
      trace_zero_curve(*iter, curve);
    }
    
  }

  // This is very similar to the gradient tracer.
  // It traces the 0 level set perpendicular to the gradient.
  bool Curve_Silhouette::trace_zero_curve(Critical_Point start, PointCurve & curve)
  {
    // Compute the canonical control points for the initial face.
    FaceIter current_face = start.face;
    std::vector<Vector3D> control_points;
    patcher.computeControlPoints(current_face, control_points);
    
    // Compute the a scalar 'dir' that represents whether we wish to go up or down.
    // The gradient will be multiplied by -1 if we wish to head towards a minnimum.
    Vector2D grad = grad_f(control_points, start.u, start.v);
    

    // The curve starts and ends at the given start point.
    curve.p1 = start;
    curve.p2 = start;

    // We will now trace the gradient in a circle.
    double u = start.u;
    double v = start.v;

    // Make sure we start on the silhouette curve.
    moveOntoLevel(u, v, current_face, control_points);

    // Store the locations of the starting point.
    double u_original = u;
    double v_original = v;

    curve.addPoint(P(control_points, u, v));
    
    movePerpGrad (u, v, current_face, control_points);
    moveOntoLevel(u, v, current_face, control_points);
    
    while(abs(u - u_original) + abs(v - v_original) > GRAD_PERP_DIST)
    {
      curve.addPoint(P(control_points, u, v));

      movePerpGrad (u, v, current_face, control_points);
      moveOntoLevel(u, v, current_face, control_points);
    }

    return true;
  }

  void Curve_Silhouette::moveOntoLevel(double & u, double & v,
				       FaceIter & current_face,
				       std::vector<Vector3D> & control_points,
				       double level)
  {

    while(F(control_points, u, v) > TOLERANCE)
    {
      Vector2D grad = grad_f_2(control_points, u, v);

      // Walk downhill. (negative uphill.)
      u -= GRAD_SQR_COEF*grad.x;
      v -= GRAD_SQR_COEF*grad.y;

      performTransitions(current_face, u, v, control_points);
    }

  }

  void Curve_Silhouette::performTransitions(
       FaceIter & current_face,
       double & u,
       double & v,
       std::vector<Vector3D> & control_points)
  {
    bool transition = false;
    
    // NOTE: may mutate current_face, u, and v.
    
    while(transitionIfNeccessary(current_face, u, v))
    {
      transition = true;
    }

    // If a transition occured, then we need to recalculate the control points.
    // FIXME : It might be a good idea to store critical points on bicubic patches.
    if(transition)
    {
      control_points.clear();
      patcher.computeControlPoints(current_face, control_points);
    }
  }

  void Curve_Silhouette::movePerpGrad(double & u, double & v,
				      FaceIter & current_face,
				      std::vector<Vector3D> & control_points)
  {

    Vector2D grad = grad_f(control_points, u, v);

    grad = grad.unit();
    
    u += -GRAD_PERP_DIST*grad.y;
    v +=  GRAD_PERP_DIST*grad.x;
    
    performTransitions(current_face, u, v, control_points);
  }

  void Curve_Silhouette::find_unique_silhouette_points(int num_critical_points,
					  std::vector<PointCurve> & edges,
			  std::vector<Critical_Point> & silhouette_points)
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
      if(iter -> has_crossing_point == false)
      {
	Critical_Point & cp1 = iter -> p1;
	Critical_Point & cp2 = iter -> p2;

	// Union these critical point sets.
	// The entire set will be above or below the level set.
	UF.op_union(cp1.index, cp2.index);
      }
    }

    // FIXME : If we want we could store those sets with crossing points, then
    // simply iterate over them and perhaps store a count of how many distinct
    // connected components remain and stop once there is only 1.
    
    // For all MS edges,
    // Add an origination point if the edge's two sets are not yet connected.
    for(auto iter = edges.begin(); iter != edges.end(); ++iter)
    {
      if(iter -> has_crossing_point == true)
      {
	Critical_Point & cp1 = iter -> p1;
	Critical_Point & cp2 = iter -> p2;

	if(!UF.connected(cp1.index, cp2.index))
	{
	  // Push
	  silhouette_points.push_back(iter-> level_set_crossing_point);

	  // Union these critical point sets.
	  UF.op_union(cp1.index, cp2.index);
	}
      }
    }
    
    return;
  }

  void Curve_Silhouette::morse_smale_complex(HalfedgeMesh& mesh,
			   std::vector<Critical_Point> & cp_all, // This goes out.
			   std::vector<PointCurve> & edges
					     )
  {
    // To build the MS complex, we will trace 4 gradient curves from
    // each saddle point to two min's and max's each.
    std::vector<Critical_Point> saddle_points;

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
      Critical_Point & saddle = *iter;
      double amount = CURVE_TRACING_STEP; // FIXME, choose a more appropiate value.
      cout << "Tracing a Gradient" << endl;
      
      // ASSUMPTION: if we choose 2 perpendicular directions,
      //             then they will go in opposite gradient directions.
      // Here we are tracing gradients in 4 directions out from the saddle point.
      for(int sx = -1; sx <= 1; sx += 2)
      for(int sy = -1; sy <= 1; sy += 2)
      {
	// These act as edges in the MS Complex.
	edges.push_back(PointCurve());

	// Create a reference to the allocated curve object.
	PointCurve & curve = edges[edges.size() - 1];
 	trace_gradient(saddle, amount*sx, amount*sy, curve);
      }
    }

    return;
  }

  // REQUIRES: curve should be empty coming in.
  // Returns false if no ending critical point was found.
  bool Curve_Silhouette::trace_gradient(
		Critical_Point cp,    // Starting point.
		double du, double dv, // Direction.
		PointCurve & curve    // OUT: Output point curve that is populated.
					)
  {

    // Compute the canonical control points for the initial face.
    FaceIter current_face = cp.face;
    std::vector<Vector3D> control_points;
    patcher.computeControlPoints(current_face, control_points);

    
    // Compute the a scalar 'dir' that represents whether we wish to go up or down.
    // The gradient will be multiplied by -1 if we wish to head towards a minnimum.
    Vector2D grad = grad_f(control_points, cp.u + du, cp.v + dv);
    double   dot  = grad.x*du + grad.y*dv; // Positive dot product then <du,dv>
    double   dir  = dot > 0 ? 1.0 : -1.0;  // heads uphill.   

    // The curve starts at the original given critical point.
    curve.p1 = cp;
    curve.addPoint(cp.location);

    // Now follow the gradient until we find a max/min.
    // Add these points to the curves as we go along.
    // We transition from one bicubic patch to another
    // when we go out of [0,1] x [0,1] bounds.
    double u = cp.u;
    double v = cp.v;
    double u_new = u + du;
    double v_new = v + dv;

    // Silhouette curve function value
    double f = F(control_points, u, v);
    // Curve heading in the direction of a Silhouette.
    // aka heading towards 0 level set.
    bool check_crossing = f*dir < 0;
    
    while( abs(u_new - u) > TOLERANCE || abs(v_new - v) > TOLERANCE )
    {      
      u = u_new;
      v = v_new;

      Vector3D location = P(control_points, u, v);
      curve.addPoint(location);

      // Perform transitions between bicubic patches.
      // There may be up to two transitions if we have run off the corner.
      bool transition = false;

      // --  Handle the logic of switching between bicubic patches.
      performTransitions(current_face, u, v, control_points);
      

      // Add the silhouette crossing point to the curve if found.
      // crossing found if the curve is no longer headed towards the 0 level set.
      if(check_crossing && F(control_points, u, v)*dir > 0)
      {
	cout << "Adding a crossing Point" << f << endl;
	
	Critical_Point crossing;
	crossing.location = location;
	crossing.u = u;
	crossing.v = v;
	crossing.type = ORIGINATION;
	crossing.index = -1;// FIXME Use something useful if needed.
	crossing.face = current_face;
	crossing.f_val = 0; // Pretend, close to 0.
	
	curve.level_set_crossing_point = crossing;
	curve.has_crossing_point = true;

	// Don't check anymore, because all integral curves are monotonic
	// increasing / decreasing. There will only be 1.
	check_crossing = false;
      }

      
      grad = dir*grad_f(control_points, u, v);
      
      u_new = u + GRAD_COEF*grad.x;
      v_new = v + GRAD_COEF*grad.y;
    }

    if(!current_face->has_critical_point)
    {
      cerr << "ERROR: CurveTracer: Function converged," <<
	"but no critical point was previously stored on this face." << endl;
      return false;
    }

    curve.p2 = current_face->critical_point;
    curve.addPoint(curve.p2.location);
    return true;
  }

  // Populates the given list with critical points found on the silhouette curve
  // function defined on the given mesh.
  void Curve_Silhouette::findCriticalPoints (HalfedgeMesh& mesh,
					     std::vector<Critical_Point> * saddle_points,
					     std::vector<Critical_Point> * all_points)
  {
    int index = 0;
    
    for(FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++)
    {
      std::vector<Vector3D> control_points;
      patcher.computeControlPoints(f, control_points);
    
      Critical_Point cp;

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
	if(search_for_critical_point(control_points, u, v, cp))
	{
	  cp.type = classify_point(control_points, u, v);
	  cp.face = f;// The cp needs to know its face for direct access
	              // to the HalfEdge Mesh connectivity structure.
	  cp.index = index; // Every cp gets a unique index.
	  index++;
	  cp.f_val = F(control_points, u, v); // Store the function evaluation at this point.	  
	  
	  // Keep special track of saddle points,
	  // because we will use them in the building of a morse smale complex.
	  if(saddle_points != NULL && cp.type == SADDLE)
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
	  f -> has_critical_point = true;

	  // FIXME : We are currently only looking for one critical point per patch.
	  // This may not be an assumption that holds water.
	  found_point = true;
	  break;
	}

	// Mark this face as not having a valid critical point thus far.
	f -> has_critical_point = false;
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
			      std::vector<Vector3D> & control_points,
			      double u, double v,
			      Critical_Point & p)
  {
    Vector2D grad;
    do
    {
      grad = GRAD_SQR_COEF*grad_mag_f(control_points, u, v);

      u -= grad.x;
      v -= grad.y;

      if(u < 0 || u > 1 || v < 0 || v > 1)
      {
	return false;
      }
    }while(abs(grad.x) > TOLERANCE || abs(grad.y) > TOLERANCE);

    // Populate the found critical point.
    p.u = u;
    p.v = v;
    p.location = P(control_points, u, v);

    return true;
  }

  Critical_Point_Type Curve_Silhouette::classify_point(
		      std::vector<Vector3D> & control_points,
		      double u, double v)
  {
    // Compute entires in the 2 by 2 symmetric hessian matrix.
    double A = F_uu(control_points, u, v);
    double B = F_uv(control_points, u, v);
    double D = F_vv(control_points, u, v);

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
  Vector3D Curve_Silhouette::P(std::vector<Vector3D> & control_points,
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
    return evaluatePatch(control_points, u, v, partial_u, partial_v);
  }

  // Normal vector to surface at point P(u,v);
  Vector3D Curve_Silhouette::N(std::vector<Vector3D> & control_points,
			       double u, double v)
  {
    return cross(P(control_points, u, v, 1, 0), // P_u
		 P(control_points, u, v, 0, 1));// P_v
  }

  // F = N dot E, this defines the silhouette when it equals 0.
  double Curve_Silhouette::F(std::vector<Vector3D> & control_points,
			     double u, double v)
  {
    return dot(N(control_points, u, v), E);
  }

  // -- 1st order partial derivatives.
  double Curve_Silhouette::F_u(std::vector<Vector3D> & control_points,
			       double u, double v)
  {
    return dot(E, cross(P(control_points, u, v, 2, 0),
			P(control_points, u, v, 0, 1)) +
	          cross(P(control_points, u, v, 1, 0),
			P(control_points, u, v, 1, 1)));
  }
  
  double Curve_Silhouette::F_v(std::vector<Vector3D> & control_points,
			       double u, double v)
  {
    return dot(E, cross(P(control_points, u, v, 1, 1),
			P(control_points, u, v, 0, 1)) +
	          cross(P(control_points, u, v, 1, 0),
			P(control_points, u, v, 0, 2)));
  }

  // -- 2nd order partial derivatives.
  double Curve_Silhouette::F_uu(std::vector<Vector3D> & control_points,
				double u, double v)
  {
    return dot(E, cross(P(control_points, u, v, 3, 0),
			P(control_points, u, v, 0, 1)) +
	       
	        2*cross(P(control_points, u, v, 2, 0),
			P(control_points, u, v, 1, 1)) +

	          cross(P(control_points, u, v, 1, 0),
			P(control_points, u, v, 2, 1))
	       );
  }    
  
  // Note: F_uv = F_vu, because of 2nd partial calculus.
  double Curve_Silhouette::F_uv(std::vector<Vector3D> & control_points,
				double u, double v)
  {
    return dot(E, cross(P(control_points, u, v, 2, 1),
			P(control_points, u, v, 0, 1)) +
	       
	          cross(P(control_points, u, v, 2, 0),
			P(control_points, u, v, 0, 2)) +

	          cross(P(control_points, u, v, 1, 0),
			P(control_points, u, v, 1, 2))
	       );
  }
  
  double Curve_Silhouette::F_vv(std::vector<Vector3D> & control_points,
				double u, double v)
  {
    return dot(E, cross(P(control_points, u, v, 1, 2),
			P(control_points, u, v, 0, 1)) +
	       
	        2*cross(P(control_points, u, v, 1, 1),
			P(control_points, u, v, 0, 2)) +

	        2*cross(P(control_points, u, v, 1, 0),
			P(control_points, u, v, 0, 3))
	       );
  }


  // Returns the gradient of f at the given uv coordinates.
  Vector2D Curve_Silhouette::grad_f(std::vector<Vector3D> & control_points,
				    double u, double v)
  {
    return Vector2D(
		    F_u(control_points, u, v),
		    F_v(control_points, u, v)
		   );
  }

  // The gradient of the magnitude of the gradient of f.
  Vector2D Curve_Silhouette::grad_mag_f(std::vector<Vector3D> & control_points,
					double u, double v)
  {
    double f_u  = F_u(control_points, u, v);
    double f_v  = F_v(control_points, u, v);
    double f_uu = F_uu(control_points, u, v);
    double f_uv = F_uv(control_points, u, v);
    double f_vv = F_vv(control_points, u, v);

    // <
    // d/du (f_u^2 + f_v^2),
    // d/dv ...
    // >
    
    return Vector2D(
		    2*f_u*f_uu + 2*f_v*f_uv,
		    2*f_u*f_uv + 2*f_v*f_vv
		   );
  }

  Vector2D Curve_Silhouette::grad_f_2(std::vector<Vector3D> & control_points,
				      double u, double v)
  {
    double f    = F_u(control_points, u, v);
    double f_u  = F_u(control_points, u, v);
    double f_v  = F_v(control_points, u, v);

    return Vector2D(
		    2*f*f_u,
		    2*f*f_v
		    );
  }

  bool Curve_Silhouette::visible(std::vector<Vector3D> & control_points,
				 double u,
				 double v)
  {
    return F(control_points, u, v) > 0;
  }

  // Finds the points along a 1D parameter aligned line in a bicubic patch.
  // where F (The silhouette defining function) is 0.
  // returns the values for the parameter along the degree of freedom where 0's occur.
  // roots that are outside of the range 0 to 1 will be disregarded for bicubic patch
  // defined meshes.
  void Curve_Silhouette::findRoots_F(std::vector<Vector3D> & control_points,
				     HalfedgeIter edge,
				     std::vector<double> & roots)
  {
    int orientation = determineOrientation(edge);

    switch(orientation)
    {
      case 0:
	findRoots_F_v(control_points, 0, roots);
	return;
      case 1:
	findRoots_F_u(control_points, 1, roots);
	return;
      case 2:
	findRoots_F_v(control_points, 1, roots);
	return;
      case 3:
	findRoots_F_u(control_points, 0, roots);
	return;
    }
  }
    
  void Curve_Silhouette::findRoots_F_u(std::vector<Vector3D> & control_points,
		     double u,
		     std::vector<double> & roots
		     )
  {
    
  }
  
  void Curve_Silhouette::findRoots_F_v(std::vector<Vector3D> & control_points,
		     double v,
		     std::vector<double> & roots
		     )
  {

  }


  // -- Critical Point and curve drawing functions.

  void CurveTracer::drawCriticalPoints()
  {
    //glDisable(GL_DEPTH_TEST);
    glPointSize(20.0);


    for(std::vector<Critical_Point>::iterator iter = critical_points.begin();
	iter != critical_points.end(); ++iter)
    {
	Vector3D p = iter -> location;
	//	 cout << p << endl;
	switch(iter -> type)
	{
	  case MIN: glColor3f(0.0, 0.0, 1.0);
	    break;
	  case MAX: glColor3f(1.0, 0.0, 0.0);
	    break;
	  case SADDLE: glColor3f(0.0, 1.0, 0.0);//green
	    break;
	  case ORIGINATION: glColor3f(1.0, 1.0, 1.0);
	    break;
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
	iter -> draw();
    }
  }


  void CurveTracer::trace_all_silhouettes(
       HalfedgeMesh& mesh, Vector3D eye_direction)
  {
    Curve_Silhouette curve_tracer(eye_direction);
    critical_points.clear();
    curves.clear();

    //curve_tracer.findCriticalPoints(*mesh, critical_points);

    // Find the Entire Morse Smale Complex.
    curve_tracer.morse_smale_complex(mesh, critical_points, curves);
    curve_tracer.find_unique_silhouette_points(
					       critical_points.size(),
					       curves,
					       critical_points);

    curve_tracer.trace_zero_level_sets(critical_points, curves);
      
    //morse_smale_edges.clear();
  }

  
  void CurveTracer::trace_parametric_curve(HalfedgeIter edge, Vector3D eye_direction)
  {
    
    // Setup the tracing parameters.
    std::vector<Vector3D> control_points;
    Curve_Silhouette curve_tracer(eye_direction); // Used for visibility checks.
    double u, v, du, dv;
    bool visible = true;

    // Allocate a new Curve Object to store the results in.
    curves.push_back(PointCurve());
    PointCurve * curve = &curves[curves.size() - 1];
    curve -> visible = visible;

    
    HalfedgeIter e0 = edge;
    
    do
    {
      determineLocationAndHeading(edge, u, v, du, dv);
      du /= STEPS_PER_PATCH;
      dv /= STEPS_PER_PATCH;
      
      // Compute the current local Bicubic Patch.
      patcher.computeControlPoints(edge->face(), control_points);

      // Add regularly spaced points on the halfedge defined curve.
      for(int i = 0; i < STEPS_PER_PATCH; i++)
      {
	u += du;
	v += dv;

	// If visibility has changed, split the curve.
	if(curve_tracer.visible(control_points, u, v) != visible)
	{
	  visible = !visible;
	  curves.push_back(PointCurve());
	  curve = &curves[curves.size() - 1];
	  curve -> visible = visible;
	}
	
	curve -> addPoint(evaluatePatch(control_points, u, v));
      }

      edge = edge->next();

      // Check for curve termination conditions.
      VertexIter vertex = edge->vertex();
      if(vertex -> isBoundary()) {break;}
      if(vertex -> degree() != 4){break;}
      edge = edge -> twin() -> next();

      cout << &edge << endl;

    }while(edge != e0); 
  }

  void CurveTracer::clearData()
  {
    critical_points.clear();
    curves.clear();
  }

  // Exports all of the data as paths and points in an svg file.
  void CurveTracer::export_to_svg()
  {
    cout << "IMPLEMENT Export_to_svg in curveTracer.cpp" << endl;
  }

}
