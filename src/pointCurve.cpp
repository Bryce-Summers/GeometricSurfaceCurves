/*
 * Point Curves
 * Written by Bryce Summers on 2/16/2016.
 */

#include "pointCurve.h"

namespace CMU462
{

  // Draw the line.
  void PointCurve::draw()
  {
    glBegin(GL_LINE_STRIP);

    for(auto iter = points.begin(); iter != points.end(); ++iter)
    {
      glVertex3dv( &(iter -> x) );
    }
    
    glEnd();
  }

}
