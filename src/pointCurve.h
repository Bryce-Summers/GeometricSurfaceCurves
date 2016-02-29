#ifndef POINT_CURVE_H
#define POINT_CURVE_H

/*
 * Point Wise curve representation.
 * Written by Bryce Summers on 2/8/2016
 *
 * Specifies a list of points that define a curve through interpolation over splines.
 * I will probably start by using something simple like ye old -A/16 + 9B/16 + 9C/16 - D/16 formula.
 * // FIXME: We might want to try using some other type of spline.
 * 
 * This class will also provide functionality for converting the representation to a 
 * 2D spline through a perspective projection of the points.
 */

#include "critical_points.h"

using namespace std;

namespace CMU462
{

  class PointCurve
  {    
   
  private:
    // Stores a list of points on a curve.
    std::vector<Vector3D> points;

    // Store the list of points as points in 2D space after a viewing projection.
    // We probably won't need this until we want to output these curves as svg
    // elements to a file.
    
    //--//std::vector<Vector2D> projected_points.

  public:

    // -- Constructor.
    PointCurve(){};
    ~PointCurve(){}

    void addPoint(Vector3D p){ points.push_back(p); }

    // Draws this point curve to the screen
    // as a set of lines.
    // with the current gl settings, e.g colors, thickness,
    // perspective projection.
    void draw();

    bool visible = true;
    
    // Extra information that can be used with these curves.
    Critical_Point p1;
    Critical_Point p2;

    // Not really a critical point, but the Critical point structure encodes all of
    // the information that we need to use this point in the future.
    // The level set crossing point is only guranteed to be close to the level set.
    Critical_Point level_set_crossing_point;
    bool has_crossing_point = false; // Stores whether we actually have a level set crossing point.

    // Add spline subdivision methods.
    // Write functionality to convert this into 2D points based on the current Opengl settings.
  };
  
}


#endif // POINT_CURVE_H
