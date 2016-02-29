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
    if(visible)
    {
      glColor3f(0.0, 1.0, 0.0);//green
    }
    else
    {
      glColor3f(0.0, 0.0, 1.0);//blue
    }
    
    glBegin(GL_LINE_STRIP);

    for(auto iter = points.begin(); iter != points.end(); ++iter)
    {
      glVertex3dv( &(iter -> x) );
    }
    
    glEnd();
  }

}
