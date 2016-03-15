#ifndef STYLE_H
#define STYLE_H

#include "CMU462/CMU462.h"
#include <string>

/*
 * Draw Stylization Information Structures.
 * Refactored to its own file by Bryce Summers on 3/14/2016.
 *
 * Uses:
 * 1. Used to represent the style for drawing the models to the screen.
 *
 * 2. Used to represent the styles for SVG representatins of geometric figures.
 */

using namespace std;

namespace CMU462 {

  // This structure is used to store a color pallete for this application.
  struct DrawStyle
  {
     Color halfedgeColor;
     Color vertexColor;
     Color edgeColor;
     Color faceColor;

     float strokeWidth;
     float vertexRadius;
  };

  // Specifies the desired behavior for drawing svg elements.
  class SVG_Style
  {
  public:
    string stroke;
    string fill;
    string stroke_width;
  };

  std::ostream& operator<<( std::ostream& os, const SVG_Style& s);
}

#endif // STYLE_H
