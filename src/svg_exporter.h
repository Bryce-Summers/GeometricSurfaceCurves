#ifndef SVG_EXPORTER_H
#define SVG_EXPORTER_H

#include <string>   // SVG creation requires quite a bit of string creation.
#include <iostream> // '<<' output stream output syntax.
#include <fstream>  // file output stream.
#include <stack>

#include "style.h" // Exportation Stylization.
#include "pointCurve.h"

using namespace std;

/*
 * SVG Exportation class.
 * Generates and exports a complete svg file.
 *
 * Written by Bryce Summers on 3/14/2016.
 *
 * The user of this class is reponsible for making sure the coordinates line up
 * and that 
 */

/* Here is a simple example svg from https://www.w3.org/TR/SVG/paths.html
 *
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="4cm" height="4cm" viewBox="0 0 400 400"
     xmlns="http://www.w3.org/2000/svg" version="1.1">
  <title>Example triangle01- simple example of a 'path'</title>
  <desc>A path that draws a triangle</desc>
  <rect x="1" y="1" width="398" height="398"
        fill="none" stroke="blue" />
  <path d="M 100 100 L 300 100 L 200 300 z"
        fill="red" stroke="blue" stroke-width="3" />
</svg>
*
*/

using namespace std;

namespace CMU462
{

  class SVG_Exporter
  {
    private:
    ofstream file; // Store the output file stream.

    Matrix4x4 Projection;

    size_t screen_w;
    size_t screen_h;

    public:
    // -- Constructor.
    SVG_Exporter(){};
    ~SVG_Exporter(){}

    // Opens
    // REQUIRES: file_name should be an svg.
    void beginSVG(string file_name,
		  string width,
		  string height,
		  string viewBox,
		  SVG_Style & default_style);
    void endSVG();

    // Information meta components.
    void addTitle(string title);
    void addDescription(string description);

    // Start a grouped collection of elements.
    // They will all be drawn in a particular style.
    // PUSHES the previous input style onto the style stack.
    void beginGroup(SVG_Style & style_in);
    // Signals the end of a group.
    // POPS the previous style from the stack.
    // Puts an end group tag in the file.
    void endGroup();

    // BE SURE TO LOAD THIS MATRIX before calling geometry functions.
    // It caches the opengl projection Matrix that converts
    // coordinates from model space to screen space.
    // Also loads in the screen width and screen height.
    void loadProjectionMatrix(size_t w, size_t h);
    
    // Supported Geometries.
    // Draws the given point curve as an SVG path in the current style.
    // Uses the 2D representation with its tangent vectors for the B-splines.
    void g_curve(PointCurve & curve);

    // void g_line();
  };

}

#endif // SVG_EXPORTER_H
