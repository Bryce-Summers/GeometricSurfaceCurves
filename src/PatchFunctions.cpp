/* Patch Functions for transversal and specification. */

#include "PatchFunctions.h"

namespace CMU462
{

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
  
  bool transitionIfNeccessary(FaceIter & current_face,
			      double & u,
			      double & v)
  {
    /*       0                 0      
     *  .--------->       .---------> 
     * /|\   .    |      /|\   .    | 
     *  |   /|\   |       |   /|\   | 
     *3 |    |    | 1   3 |    |    | 1
     *  |    |    | o1  o2|    |    |
     *  |        \|/      |        \|/
     *  <---------.       <---------.
     *       2                 2
     *
     * Edge 1 in the starting face goes to edge 3 in the twin face.
     * orientation1 (o1) is therefore 1, and orientation2 (o2) is therefore 3.
     *
     * We can determine how many clockwise turns are needed to make the coordinates 
     * consistent after the transformation as follows:
     *
     * o2 == (o1 + 2) % 4 --> 0 clockwise turns.
     * o2 == (o1 + 3) % 4 --> 1 clockwise turns.
     * o2 == (o1 + 4) % 4 --> 2 clockwise turns.
     * o2 == (o1 + 5) % 4 --> 3 clockwise turns.
     *
     * Here are the transformations for u and v after a number of turns.
     * (u,v) --> (u - .5, v - .5) [Translation centers coordinate system around origin.]
     *
     * [cos(90) -sin(90) ]   x number of clockwise turns.
     * [sin(90)  cos(90) ]
     *
     * [   0       -1      ] [u] = [-v ]    u' = -v;
     * [   1        0      ] [v] = [ u ]    v' =  u;
     *
     * (u, v) --> (u + .5, v + .5) [Translate the cordinate back into [0, 1) x [0, 1)
     */

    // NOTE: We are making a major assumption that we will only be
    //       performing 1 transformation at a time.
    
    HalfedgeIter e0, e1, e2, e3;
    e0 = current_face -> halfedge(); // Top.
    e1 = e0 -> next(); // Right.
    e2 = e1 -> next(); // Bottom.
    e3 = e2 -> next(); // Left.

    int o1; // The index of the original face's edge.
    int o2; // The index of the twin face's edge.
    HalfedgeIter twin; // The twin transition edge with index o2.

    // Compute the transition variables and apply local transformation.
    if(v < 0.0)
    {
      twin = e0 -> twin();
      o1 = 0;
      v = v + 1.0;
    }
    else if(u > 1.0)
    {
      twin = e1 -> twin();
      o1   = 1;
      u    = u - 1.0;
    }
    else if(v > 1.0)
    {
      twin = e2 -> twin();
      o1   = 2;
      v    = v - 1.0;
    }
    else if(u < 0.0)
    {
      twin = e3 -> twin();
      o1   = 3;
      u    = u + 1.0;
    }
    else
    {
      return false;
    }

    o2 = determineOrientation(twin);
    current_face = twin -> face();

    u -= .5;
    v -= .5;
    /*
     * o2 == (o1 + 2) % 4 --> 0 clockwise turns.
     * o2 == (o1 + 3) % 4 --> 1 clockwise turns.
     * o2 == (o1 + 4) % 4 --> 2 clockwise turns.
     * o2 == (o1 + 5) % 4 --> 3 clockwise turns.
     */
    int rotations = (o2 + 4 - o1 + 2) % 4;

    /*
     * u' = -v;
     * v' =  u;
     */

    for(int i = 0; i < rotations; i++)
    {
      double u_new = -v;
      double v_new = u;

      u = u_new;
      v = v_new;
    }

    u += .5;
    v += .5;

    return true;
  }

  // Returns how deep the given edge is in its on edge's cycle.
  inline int determineOrientation(HalfedgeIter edge)
  {
    HalfedgeIter e0 = edge->face()->halfedge();
    int output = 0;
    while(e0 != edge)
    {
      output++;
      e0 = e0 -> next();
    }

    return output;
  }

  void determineLocationAndHeading(HalfedgeIter edge,
				   double & u,
				   double & v,
				   double & du,
				   double & dv)
  {
    int orientation = determineOrientation(edge);

    switch(orientation)
    {
      case 0:
	 u = 0.0;  v = 0.0;
	du = 1.0; dv = 0.0;
	return;
      case 1:
	 u = 1.0;  v = 0.0;
	du = 0.0; dv = 1.0;
	return;
      case 2:
	 u = 1.0;   v = 1.0;
	du = -1.0; dv = 0.0;
	return;
      case 3:
	 u =  0.0;  v =  1.0;
	du =  0.0; dv = -1.0;
	return;
    }

    cerr << "DetermineLocationAndHeading: Something has gone terribly wrong." << endl;
  }

  bool findCubicRoots(double a, double b, double c, double d, // IN
		      double & r1,  // OUT
		      double & r2,  // OUT
		      double & r3)  // OUT
  {

  }
  
  bool findCubicRoots(double a, double b, double c, double d, // IN
		      Complex r1,  // OUT
		      Complex r2,  // OUT
		      Complex r3)  // OUT
  {

    // Algorithm from https://en.wikipedia.org/wiki/Cubic_function.
    double disc = 18*a*b*c*d - 4*b*b*b*d + b*b*c*c - 4*a*c*c*c - 27*a*a*d*d;

    if(disc > 0)
    {
      // Equation has 3 distinct real roots.
    }
    else if(disc == 0.0)
    {
      // The equation has a multiple root and all its roots are real.
    }
    else // dist < 0.
    {
      // The equation has one real root and two nonreal complex conjugate roots.
    }

    double delta0, delta1, delta2;
    delta0 = b*b - 3*a*c;
    delta1 = 2*b*b*b - 9*a*b*c + 27*a*a*d;
  
    // delta1^2 - 4*delta^3
    double diff = 27*a*a*disc;

    Complex bigC(diff, 0);
    bigC = ((delta1 + bigC.pow(.5))/2);
    bigC = bigC.pow(.33);
    
    // 3rd order Roots of Unity.
    Complex u1(1, 0);
    Complex u2(-1/2, sqrt(3)/2);
    Complex u3(-1/2, -sqrt(3)/2);

    double a_constant = -1.0/(3*a);
    
    r1 = a_constant*(b + u1*bigC + delta0/(u1*bigC));
    r2 = a_constant*(b + u2*bigC + delta0/(u2*bigC));
    r3 = a_constant*(b + u3*bigC + delta0/(u3*bigC));
    
  }

  bool findQuadraticRoots(double a, double b, double c,
			  double & r1, double & r2)
  {
    Complex r1_temp;
    Complex r2_temp;

    bool out = findQuadraticRoots(a, b, c, r1_temp, r2_temp);

    r1 = r1_temp.Re();
    r2 = r2_temp.Re();
    return  out;
  }
  
  bool findQuadraticRoots(double a, double b, double c,
			  Complex & r1,
			  Complex & r2)
  {
      // Quadratic formula.
      // (-b +- sqrt(b^2 - 4ac))/(2a)
      
      // Discirminant.
      double disc = b*b - 4*a*c;
      double a2 = a*2;

      // Handle Complex return values.
      if(disc < 0)
      {
	r1 = Complex(-b/a2, sqrt(-disc)/a2);
	r2 = Complex(r1.Re(), -r1.Im());
	return false;
      }

      r1 = Complex((-b + disc)/a2, 0.0);
      r2 = Complex((-b - disc)/a2, 0.0);

      return true;
      
  }

}