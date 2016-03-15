#include "svg_exporter.h"

using namespace std;

// -- Constructor.

void SVG_Exporter::beginSVG(string file_name,
	      string width,
	      string height,
	      string viewBox,
	      SVG_Style & default_style)
{
  // Initialize the style attributes and style stack.

  // FIXME: Initialize svg wide styles.

  // Open the desired output file for writing.
  file.open(file_name);

  /* Write the following header to the file:
   *
   *
    <?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
    <svg width="4cm" height="4cm" viewBox="0 0 400 400"
    xmlns="http://www.w3.org/2000/svg" version="1.1">
   *
   */
 
  file << "<?xml version=\"1.0\" standalone=\"no\"?>" << endl;
  file << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << endl;
  file << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl;
  file << "<svg width="
       << "\"" << width  << "\" "
       << "height="
       << "\"" << height << "\" "
       << "viewBox="
       << "\"" << viewBox << "\"" << endl
       << "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << endl;
}
void SVG_Exporter::endSVG()
{
  file << "</svg>" << endl;
  file.close();
}

// Information meta components.
void SVG_Exporter::addTitle(string title)
{
  file << "<title>" << title << "</title>" << endl;
}

void SVG_Exporter::addDescription(string description)
{
  file << "<desc>" << description << "</desc>" << endl;
}

// Start a grouped collection of elements.
// They will all be drawn in a particular style.
// PUSHES the previous input style onto the style stack.
void SVG_Exporter::beginGroup(SVG_Style & style_in)
{
  file << "<g " << style_in << ">" << endl;
}

// Signals the end of a group.
// POPS the previous style from the stack.
// Puts an end group tag in the file.
void SVG_Exporter::endGroup()
{
  file << "</g>" << endl;
}

// Supported Geometries.
// Draws the given point curve as an SVG path in the current style.
// Uses the 2D representation with its tangent vectors for the B-splines.
void SVG_Exporter::g_curve(PointCurve & curve)
{
  // Appends an svg path element to the file representing the given
  // curve's projection onto the current opengl projection screen.
  curve.export_svg_path(file, Projection, screen_w, screen_h);
}

void SVG_Exporter::loadProjectionMatrix(size_t w, size_t h)
{
  // Populate SVG instance level variables.
  screen_w = w;
  screen_h = h;
  
  GLdouble projMatrix [16];
  GLdouble modelMatrix[16];

  for(int i = 0; i < 16; i++)
  {
    projMatrix[i]  = 0.0;
    modelMatrix[i] = 0.0;
  }

  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetDoublev(GL_MODELVIEW_MATRIX,  modelMatrix);


  Matrix4x4 P;
  Matrix4x4 M;

  for(int r = 0; r < 4; r++)
  for(int c = 0; c < 4; c++)
  {
    P(r, c) = projMatrix [4*c + r];
    M(r, c) = modelMatrix[4*c + r];
  }

  Projection = P*M;
}
