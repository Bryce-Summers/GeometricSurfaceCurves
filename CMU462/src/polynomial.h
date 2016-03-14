#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

/*
 * A Polynomial Class.
 * Single variable, real (double) coeficients.
 * Written by Bryce Summers on 2/25/2016.
 *
 * USAGE: Polynomials are only the size of a pointer to a deque.
 *        Therefore you are free to stack allocate them to your heart's
 *        content.
 */

#include <deque>
#include <vector>
#include <cmath>
#include <stack>
#include <ostream>

namespace CMU462
{

  class Polynomial;
  
  // Sturm Chains are sequences of polynomials.
  typedef std::vector<Polynomial> SturmChain;
  
  class Polynomial
  {
  private:

    // A list of the coeficients of this polynomial where
    // data[i] represents the coeficient of the x^i term in the polynomial.
    // We will use a deque to make in place operations such as differentiation efficient.
    std::deque<double> data;

  public:
    // -- Constructor.
    // Takes in a vector of coeficients, doesn't mutate the input.
    Polynomial(std::vector<double> & coefs);
    Polynomial()
    {
      //data = new std::deque<double>();
    } // Empty polynomial of degree -1. (0.0)
    ~Polynomial()
    {
      //delete data;
    }

    // Evaluate the polynomial at the given location.
    double eval(double x) const;
    
    // Appends a leading coeficient term of degree degree() + 1.
    void addLeadingTerm(double coef);
    void removeLeadingTerm();
    void addTrailingTerm(double coef);
    void removeTrailingTerm();

    void addTrailingTerms(double c1);
    void addTrailingTerms(double c1, double c2);
    void addTrailingTerms(double c1, double c2, double c3);
    void addTrailingTerms(double c1, double c2, double c3, double c4);
    
    // REQUIRES: 0 <= power.
    // RETURNS the coeficient if it exists and returns 0.0 if power > degree().
    double getCoef(int power) const;

    // Accessing operators.
    inline double& operator[] (int index) {
      return data[index];
    }

    // Returns the exponent of the highest coeficient in this polynomial.
    int degree() const;

    // Calculus operations.
    Polynomial differentiate() const;
    void differentiate_in_place();
    Polynomial integrate(double constant = 0.0) const;

    // Assignment Operator.
    void operator=(const Polynomial &other);
    
    // Arithmetic operations in both assignment and input preserving varieties.
    void         operator+=(const Polynomial & other) ;
    Polynomial   operator+ (const Polynomial & other) const;
    
    void         operator-=(const Polynomial & other) ;
    Polynomial   operator- (const Polynomial & other) const;
    
    void         operator*=(const Polynomial & other) ;
    Polynomial   operator* (const Polynomial & other) const;
    void         operator/=(const Polynomial & other) ;
    Polynomial   operator/ (const Polynomial & other) const;

    // Returns the remainder of this / other.
    void         operator%=(const Polynomial & other) ;
    Polynomial   operator% (const Polynomial & other) const;

    void         operator*=(double scalar) ;
    Polynomial   operator* (double scalar) const;
    Polynomial   operator- () const;

    void         operator/=(double scalar) ;
    Polynomial   operator/ (double scalar) const;

    // This does all of the leg work for computing quotients and remainders.
    void division(const Polynomial & other, // INPUT.
		  Polynomial & quotient, Polynomial & remainder) const; // OUTPUT
  
    // Returns all of the real roots for this polynomial.
    // Uses the Sturm chain to numerically find the roots
    // within the interval lower to upper.
    // All roots will be computer within the given tolerance of their actual
    // locations. TOLERANCE therefore specifies the prescision of the results.
    void computeRealRoots(std::vector<double> & roots,
			  double lower, double upper,
			  double TOLERANCE) const;

    Polynomial copy () const;

    // Returns the Sturm Chain associated with this polynomial.
    void computeSturmChain(SturmChain & chain) const;

    // RETURNS true iff this polynomial is equal to 0, i.e degree() < 0.
    bool isZero() const;

    // Makes this polynomial identical to the other.
    void load(const Polynomial & other);

    // degree -1 trivial polynomial.
    inline static Polynomial ZERO()
    {
      return Polynomial();
    };
    
  private:

    // Private constructor, initializes a polynomial of the given degree.
    // with all 0.0 coeficients.
    Polynomial(int degree);
    
    // REQUIRES: power <= degree() of this polynomial.
    void setCoef(int power, double coef);

    // Returns the number of times the sign flips in the sturm chain.
    int sturm_sign_flips(SturmChain & sturm_chain, double input) const;

    // Removes leading 0's.
    // simplifies this polynomial to canonical form.
    void reduce();

    class Interval
    {
    public:
      double lower, upper;
      int flips_lower, flips_upper;
    };
	
  };
  
  inline void operator*=(double s, Polynomial &a)
  {
    a *= s;
  }
  
  inline Polynomial operator*(double s, const Polynomial &a)
  {
    return a * s;
  }
  
  std::ostream& operator<<( std::ostream& os, const Polynomial &poly);
}

#endif // POLYNOMIAL_H
