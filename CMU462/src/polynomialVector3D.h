#ifndef CMU462_POLYNOMIAL_VECTOR3D_H
#define CMU462_POLYNOMIAL_VECTOR3D_H

#include <ostream>
#include <cmath>
#include "polynomial.h"
#include "vector3D.h"

namespace CMU462 {

/**
 * Defines 3D vectors.
 */
class PolynomialVector3D {
 public:

  // components
  Polynomial x;
  Polynomial y;
  Polynomial z;

  /**
   * Constructor.
   * Initializes tp vector (0,0,0).
   */
  PolynomialVector3D()
  {
    // Initialize all polynomials to 0.
    x = Polynomial();
    y = Polynomial();
    z = Polynomial();
  }

  /**
   * Constructor.
   * Initializes to vector (x,y,z).
   * Makes copies of each polynomial, so as not to disturb the
   * input polynomials.
   */
  PolynomialVector3D(Polynomial x_in, Polynomial y_in, Polynomial z_in)
  {
    x = x_in;
    y = y_in;
    z = z_in;
  }

  /**
   * Constructor.
   * Initializes to vector (c,c,c)
   */
  PolynomialVector3D( Polynomial c)
  {
    x = c;
    y = c;
    z = c;
  }

  /**
   * Constructor.
   * Initializes from existing vector
   */
  PolynomialVector3D(const PolynomialVector3D &v)
  {
    x = v.x.copy();
    y = v.y.copy();
    z = v.z.copy();
  }

  // returns reference to the specified component (0-based indexing: x, y, z)
  inline Polynomial& operator[] ( const int& index ) {
    return ( &x )[ index ];
  }

  // returns const reference to the specified component (0-based indexing: x, y, z)
  inline const Polynomial& operator[] ( const int& index ) const {
    return ( &x )[ index ];
  }

  // negation
  inline PolynomialVector3D operator-( void ){
    Polynomial x_new = -x;
    Polynomial y_new = -y;
    Polynomial z_new = -z;
    PolynomialVector3D result = PolynomialVector3D( x_new, y_new, z_new);
    return result;
  }

  // addition
  inline PolynomialVector3D operator+( const PolynomialVector3D& v ) const {
    return PolynomialVector3D( x + v.x, y + v.y, z + v.z );
  }

  // subtraction
  inline PolynomialVector3D operator-( const PolynomialVector3D& v ) const {
    return PolynomialVector3D( x - v.x, y - v.y, z - v.z );
  }

  // right scalar multiplication
  inline PolynomialVector3D operator*( const double& c ) const {
    return PolynomialVector3D( x * c, y * c, z * c );
  }

  // scalar division
  inline PolynomialVector3D operator/( const double& c ) const {
    const double rc = 1.0/c;
    return PolynomialVector3D( rc * x, rc * y, rc * z );
  }

  // addition / assignment
  inline void operator+=( const PolynomialVector3D& v ) {
    x += v.x; y += v.y; z += v.z;
  }

  // subtraction / assignment
  inline void operator-=( const PolynomialVector3D& v ) {
    x -= v.x; y -= v.y; z -= v.z;
  }

  // scalar multiplication / assignment
  inline void operator*=( const double& c ) {
    x *= c; y *= c; z *= c;
  }

  // scalar division / assignment
  inline void operator/=( const double& c ) {
    (*this) *= ( 1./c );
  }

  /**
   * Returns Euclidean length.
   */
  /*
  inline Polynomial norm( void ) const {
    return sqrt( x*x + y*y + z*z );
  }
  */

  /**
   * Returns Euclidean length squared.
   */
  inline Polynomial norm2( void ) const {
    return x*x + y*y + z*z;
  }

}; // class Vector3D

// left scalar multiplication
inline PolynomialVector3D operator* ( const double& c, const PolynomialVector3D& v ) {
  return PolynomialVector3D( c * v.x, c * v.y, c * v.z );
}

// dot product (a.k.a. inner or scalar product)
inline Polynomial dot( const PolynomialVector3D& u, const PolynomialVector3D& v ){
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

// cross product
inline PolynomialVector3D cross( const PolynomialVector3D& u, const PolynomialVector3D& v ) {
  return PolynomialVector3D( u.y*v.z - u.z*v.y,
			     u.z*v.x - u.x*v.z,
			     u.x*v.y - u.y*v.x );
}

// -- Mixed Type operators.
 
inline PolynomialVector3D operator* ( const Polynomial p, const Vector3D& v ) {
  return PolynomialVector3D( p*v.x, p*v.y, p*v.z);
}

inline PolynomialVector3D operator* (const Vector3D& v,  const Polynomial p) {
  return p*v;
}

// dot product (a.k.a. inner or scalar product)
inline Polynomial dot(const PolynomialVector3D& u, const Vector3D& v ){
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline Polynomial dot(const Vector3D& v, const PolynomialVector3D& u){
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

// prints components
std::ostream& operator<<( std::ostream& os, const PolynomialVector3D& v );

} // namespace CMU462

#endif // CMU462_POLYNOMIAL_VECTOR3D_H
