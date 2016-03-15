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

  inline Vector2D w2s(Vector3D & vec_in, Matrix4x4 & Proj, size_t width, size_t height)
  {
    Vector4D vec = Vector4D(vec_in);
    vec.w = 1.0;

    // Project into 3d Homogeneous coordinates.
    vec  = Proj*vec;

    // Convert from 3d Homgeneous to [-1, -1, -1], [1, 1, 1] cube.
    vec /= vec.w;
    
    Vector2D out;
    out.x = (vec.x + 1)*width/2.0;
    out.y = (-vec.y + 1)*height/2.0;

    return out;   
  }
  
  // Exports this point Curve as a 2D svg path.
  void PointCurve::export_svg_path(std::ostream& file, Matrix4x4 &P,
				   size_t w, size_t h)
  {
    /*
    <path d="M 100 100 L 300 100 L 200 300 z"
        fill="red" stroke="blue" stroke-width="3" />
    */
    
    // Paths are given as a list of commands (M, L, B)
    // and coordinate locations. (100, 100)
    // M means move, with pen up.
    // L means straight line to, pen down.
    // H x+ means horizontal line.
    // V y+ means vertical line.
    // B
    // Z means 'close the path'
    // capitlization matters:
    // -- CAPITAL means absolute coordinates,
    // -- lowercase means relative coordinates.

    // Start of the path element.

    int len = points.size();
    if(len < 2)
    {
      cout << "PointCurve: " <<
	      "PATH's must contain at least  is not long enough." << endl;
      return;
    }

    file << "<path d=\"";

    // -- Start the curve.
    Vector2D p0 = w2s(points[0], P, w, h);
    Vector2D p1 = w2s(points[1], P, w, h);
    // Positive tangent offset control point from p0.
    Vector2D c0 = p0;// + Vector4D(tangents[0]);
    // Negative tangent offset control point from p1.
    Vector2D c1 = p1 ;//- Vector4D(tangents[1]);

    file << "M" << p0.x <<","<< p0.y << " ";
    file << "C" << c0.x <<","<< c0.y << " ";
    file <<        c1.x <<","<< c1.y << " ";
    file <<        p1.x <<","<< p1.y << " ";

    for(int i = 2; i < len; i++)
    {
      Vector2D pos  = w2s(points[i], P, w, h);
      Vector2D ctrl = pos;// - Vector4D(tangents[i]);

      file << "S" << ctrl.x <<","<< ctrl.y << " ";
      file << " " << pos.x  <<","<< pos.y  << " ";
    }

    // Draw the connecting Bezier curve.
    if(closed)
    {
      Vector2D pos  = w2s(points[0], P, w, h);
      Vector2D ctrl = pos;// - Vector4D(tangents[0]);

      file << "S" << ctrl.x <<","<< ctrl.y << " ";
      file << " " << pos.x  <<","<< pos.y  << " ";
      
      // file << " z"; // Close with line.
    }
    
    // End of the path element.
    file << "\"/>" << endl;
   }

  
}
