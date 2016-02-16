#ifndef CRITICAL_POINT_H
#define CRITICAL_POINT_H

/*
 * Critical Point specification.
 * Refactored to an independant file on 2/15/2016.
 */

#include "halfEdgeMesh.h"

extern FaceIter;

using namespace std;


namespace CMU462
{
  enum Critical_Point_Type{MIN, MAX, SADDLE, ORIGINATION};
  
  class Critical_Point
  {
  public:
    Vector3D location; // This Critical Point's location in 3D space.
    FaceIter face;     // The face defining the bicubic patch that this
                       // critical point resides on.
    double u, v;       // the coordinates of this cp on the bicubic patch.
    Critical_Point_Type type; // MIN, MAX, or SADDLE.
    int index; // Every critical point is given a unique index between [0, |critical points|]
    double f_val; // The value of the function at this critical point.
  };

}

#endif // CRITICAL_POINT_H
