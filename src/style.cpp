#include "style.h"

namespace CMU462 {

  // Converts an SVG Style into xml properties that may be embedded in
  // svg tags, such as groups.
  std::ostream& operator<<( std::ostream& os, const SVG_Style& s)
  {
    os << " stroke=\"" << s.stroke << "\" fill=\"" << s.fill
       << "\" stroke-width=\"" << s.stroke_width << "\"";

    return os;
  }


}
