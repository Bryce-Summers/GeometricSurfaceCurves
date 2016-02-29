// This file provides the implementations for 3D Vectors
// Containing Polynomials.

#include "polynomialVector3D.h"

namespace CMU462 {

  std::ostream& operator<<( std::ostream& os, const PolynomialVector3D & v ) {
    os << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
    return os;
  }

} // namespace CMU462
