#include "polynomial.h"
#include <stdexcept>
#include <iostream>

namespace CMU462
{

  Polynomial::Polynomial(std::vector<double> & coefs)
  {
    //data = new std::deque<double>();
    // -- Simply copy the incoming coefs.
    // ASSUMES they are in order of lowest to highest degree.
    // coef of x^i is represented by coefs[i];
    for(auto iter = coefs.begin(); iter != coefs.end(); ++iter)
    {
      data.push_back(*iter);
    }
  }

  Polynomial::Polynomial(int degree)
  {
    //  data = new std::deque<double>();
    for(int i = 0; i <= degree; i++)
    {
      data.push_back(0.0);
    }
  }

  double Polynomial::eval(double x) const
  {

    // Uses Horner's method.
    int deg = degree();
    double out = 0.0;

    for(int i = deg; i >= 0; i--)
    {
      out = out*x + getCoef(i);
    }

    return out;
  }

  // 0.0 has a degree of < 0. (Theoretically it is -infinity.)
  int Polynomial::degree() const
  {
    return data.size() - 1;
  }

  void Polynomial::addLeadingTerm(double coef)
  {
    data.push_back(coef);
  }

  void Polynomial::removeLeadingTerm()
  {
    data.pop_back();
  }

  // Increases the degree of all terms and sets the coeficient of x^0 to 0.0.
  void Polynomial::addTrailingTerm(double coef)
  {
    data.push_front(coef);
  }
  
  void Polynomial::addTrailingTerms (double c1)
  {
    data.push_front(c1);
  }
  
  void Polynomial::addTrailingTerms(double c1, double c2)
  {
    data.push_front(c1);
    data.push_front(c2);
  }

  void Polynomial::addTrailingTerms(double c1, double c2, double c3)
  {
    data.push_front(c1);
    data.push_front(c2);
    data.push_front(c3);
  }
  
  void Polynomial::addTrailingTerms(double c1, double c2, double c3, double c4)
  {
    data.push_front(c1);
    data.push_front(c2);
    data.push_front(c3);
    data.push_front(c4);
  }

  // Removes the lowest degree term and lowers the degree of all terms by 1.
  void Polynomial::removeTrailingTerm()
  {
    data.pop_front();
  }
  
  double Polynomial::getCoef(int power) const
  {
   
    return power > degree()
           ? 0.0
           : data[power];
  }

  void Polynomial::setCoef(int power, double coef)
  {
    data[power] = coef;
  }
  
  Polynomial Polynomial::differentiate() const
  {    
    int deg = degree();
    Polynomial out;

    // Copy over coeficients at a 1 shift and multiplied by the exponents.
    for(int i = 1; i <= deg; i++)
    {
      out.addLeadingTerm(data[i]*i);
    }
    
    return out;
  }
  
  void Polynomial::differentiate_in_place()
  {
    int deg = degree();

    for(int i = 1; i <= deg; i++)
    {
      data[i] *= i;
    }

    // Remove the previous degree 0 term.
    removeLeadingTerm();
  }

  Polynomial Polynomial::operator+(const Polynomial & other) const
  {
    int deg = std::max(degree(), other.degree());
    Polynomial out;

    // Copy over coeficients at a 1 shift and multiplied by the exponents.
    for(int i = 0; i <= deg; i++)
    {
      out.addLeadingTerm(getCoef(i) + other.getCoef(i));
    }

    out.reduce();
    return out;
  }
  
  void Polynomial::operator+=(const Polynomial & other)
  {
    int deg = degree();

    // Add all of the similar terms.
    for(int i = 0; i <= deg; i++)
    {
      data[i] += other.getCoef(i);
    }

    // Append all of the terms that other has and this does not.
    int deg_other = other.degree();
    for(int i = deg + 1; i <= deg_other; i++)
    {
      addLeadingTerm(other.getCoef(i));
    }

    reduce();
  }
  
  Polynomial Polynomial::operator-(const Polynomial &other) const
  {
    int deg = std::max(degree(), other.degree());
    Polynomial out;

    // Copy over coeficients at a 1 shift and multiplied by the exponents.
    for(int i = 0; i <= deg; i++)
    {
      out.addLeadingTerm(getCoef(i) - other.getCoef(i));
    }

    out.reduce();
    return out;
  }
  
  void Polynomial::operator-=(const Polynomial & other)
  {
    int deg = degree();

    // Subtract out all of the terms shared by the polynomials.
    for(int i = 0; i <= deg; i++)
    {
      data[i] -= other.getCoef(i);
    }

    // Subtract out all of the terms contained in other, but not in this.
    int deg_2 = other.degree();
    for(int i = deg + 1; i <= deg_2; i++)
    {
      addLeadingTerm(-(other.getCoef(i)));
    }

    reduce();
  }  

  Polynomial Polynomial::operator*(const Polynomial & other) const
  {
    // Compute the size of the output Polynomial.
    int end_degree = degree() + other.degree();

    // Compute which polynomial we will use as the first operand,
    // because we can conserve additions in this way.
    const Polynomial *p1;
    const Polynomial *p2;
    int min_degree;
		
    if(degree() < other.degree())
    {
      min_degree = degree();
      p1 = this;
      p2 = &other;
    }
    else
    {
      min_degree = other.degree();
      p1 = &other;
      p2 = this;
    }
		
    // Start with an empty polynomial.
    Polynomial output = ZERO();

    // Now perform polynomial multiplication.
    for(int r = 0; r <= min_degree; r++)
    {
      double coef1 = p1 -> getCoef(r);

      Polynomial temp(end_degree + 1);
			
      for(int c = 0; c <= p2 -> degree(); c++)
      {
	double coef2 = p2 -> getCoef(c);
	
	// Multiply powers of x and set them properly.
	temp[c + r] = coef1*coef2;
      }

      temp.reduce();

      // Add up the component sums.
      output += temp;
    }

    output.reduce();
    return output;
  }

  // Negation operation.
  Polynomial Polynomial::operator-() const
  {
    return (*this)*(-1.0);
  }
  
  Polynomial Polynomial::operator*(double scalar) const
  {
    // Multiplying by 0 produces the zero polynomial.
    if(scalar == 0.0)
    {
      return ZERO();
    }

    Polynomial out;
    int deg = degree();

    // Copy over coeficients at a 1 shift and multiplied by the exponents.
    for(int i = 0; i <= deg; i++)
    {
      out.addLeadingTerm(getCoef(i)*scalar);
    }

    out.reduce();
    return out;
  }

  void Polynomial::operator*=(double scalar)
  {
    // Multiplying by 0 produces the zero polynomial.
    if(scalar == 0.0)
    {
      data.clear();
      return;
    }
    
    int deg = degree();

    // Copy over coeficients at a 1 shift and multiplied by the exponents.
    for(int i = 0; i <= deg; i++)
    {
      data[i] *= scalar;
    }
  }

  void Polynomial::division(const Polynomial & other,// IN
			    Polynomial & quotient_out,
			    Polynomial & remainder_out)// OUT
    const
  {
    Polynomial dividend = this -> copy();
    Polynomial divisor  = other.copy();
    
    int deg_1 = this -> degree();
    int deg_2 = divisor.degree();

    // 0.0 --> ERROR.
    if(deg_2 < 0)
    {
      throw new std::runtime_error("Polynomial Division by 0 Error!");
      return;
    }

    int deg_quotient = deg_1 - deg_2;

    if(deg_quotient < 0)
    {
      quotient_out = ZERO();
      remainder_out = copy();
      return;
    }

    // Initialize the output quotient.

    // a / b, a is the dividend, b is the divisor.
    // a - b, a is the minuend,  b is the subtrahend.

    Polynomial quotient = Polynomial(deg_quotient);
    for(int i = 0; i < deg_quotient; i++)
    {
      divisor.addTrailingTerm(0.0);
    }

    for(int i = deg_quotient; i >= 0; i--)
    {
      int j = deg_1 - deg_quotient + i;
      double coef = dividend[j] / divisor[j];
      quotient[i] = coef;
      Polynomial subtrahend = divisor*coef;
      
      // Remove terms that should cancel out to 0 anyways.
      // This prevents numerical cancellation errors.
      dividend  .removeLeadingTerm();
      subtrahend.removeLeadingTerm();
      dividend -= subtrahend;
      divisor.removeTrailingTerm();
    }

    dividend.reduce();
    quotient.reduce();
  
    remainder_out = dividend;
    quotient_out  = quotient;
   
    return;
  }

  Polynomial Polynomial::operator/ (const Polynomial & other) const
  {
    Polynomial quotient;
    Polynomial remainder;

    division(other, quotient, remainder);
   
    return quotient;
  }

  // Returns the remainder of this / other.
  Polynomial Polynomial::operator% (const Polynomial & other) const
  {
    Polynomial quotient;
    Polynomial remainder;

    division(other, quotient, remainder);

    return remainder;
  }  

  void Polynomial::operator/=(double scalar)
  {
    (*this) *= (1.0 / scalar);
  }
  
  Polynomial Polynomial::operator/ (double scalar) const
  {
    return (*this) * (1.0 / scalar);
  }

  // Reduces this polnomial to canonical form by removing leading 0's.
  void Polynomial::reduce()
  {
    while(degree() >= 0 && getCoef(degree()) == 0.0)
    {
      removeLeadingTerm();
    }
  }

  // -- Copying and copy (load) based operations.
  
  Polynomial Polynomial::copy() const
  {
    Polynomial out;

    for(auto iter = data.begin(); iter != data.end(); ++iter)
    {
      out.addLeadingTerm(*iter);
    }

    return out;
  }

  // Copies the given polynomial's coeficients into this polynomial.
  void Polynomial::load(const Polynomial & other)
  {
    const std::deque<double> & data_in = other.data;

    data.clear();

    for(auto iter = data_in.begin(); iter != data_in.end(); ++iter)
    {
      data.push_back(*iter);
    }
  }
  
  void Polynomial::operator/=(const Polynomial & other)
  {
    Polynomial result = (*this) / other;
    load(result);
  }

  void Polynomial::operator*=(const Polynomial & other)
  {
    Polynomial result = (*this) * other;
    load(result);
  }

  void Polynomial::operator%=(const Polynomial & other)
  {
    Polynomial result = (*this) % other;
    load(result);
  }

  bool Polynomial::isZero() const
  {
    return data.size() == 0;
  }

  // Returns all of the real roots for this polynomial.
  // Uses the Sturm chain to numerically find the roots.
  // populates the given roots array with all roots found within the
  // interval (lower, upper]
  void Polynomial::computeRealRoots(std::vector<double> & roots, // OUT.
				    double lower_in, double upper_in,
				    double TOLERANCE) // IN.
    const
  {
    SturmChain sturm_chain;
    computeSturmChain(sturm_chain);

    Interval start_interval;
    start_interval.lower = lower_in;
    start_interval.upper = upper_in;
    start_interval.flips_lower = sturm_sign_flips(sturm_chain, lower_in);
    start_interval.flips_upper = sturm_sign_flips(sturm_chain, upper_in);

    // Terminate if we don't detect any real roots in the input interval.
    if(start_interval.flips_lower == start_interval.flips_upper)
    {
      return; // No roots in interval.
    }

    // Keep a stack of intervals;
    // In this implementation we will push intervals lower, then upper.
    // When popping upper will come out first, then lower.
    std::stack<Interval> S;
    S.push(start_interval);

    // Binary search for all roots.
    // THEOREM: the number of distinct roots in the interval (a, b]
    //          is sign_flips(a) - sign_flips(b).
    while(!S.empty())
    {
      Interval I = S.top();
      S.pop();

      // Zoom in.
      while(I.upper - I.lower > TOLERANCE)
      {
	double diff = I.upper - I.lower;
	double mid  = I.lower + diff/2;

	int flips_mid = sturm_sign_flips(sturm_chain, mid);

	int roots_left  = I.flips_lower - flips_mid;
	int roots_right = flips_mid - I.flips_upper;

	// CASE 1: roots only in left interval.
	if(roots_left > 0 && roots_right == 0)
	{
	  I.upper = mid;
	  I.flips_upper = flips_mid;
	  continue;
	}
	// CASE 2: roots only in the right interval.
	else if(roots_right > 0 && roots_left == 0)
	{
	  I.lower = mid;
	  I.flips_lower = flips_mid;
	  continue;
	}
	else if(roots_left == 0 && roots_right == 0)
	{
	  //throw new std::runtime_error("ERROR: We've lost the roots!!");
	  std::cout << "We have lost a roots..." << std::endl;
	  break;
	}

	// CASE 3: roots on both sides of the midpoint.
	// Split the intervals, store the right interval for later.

	// We construct a new interval and store it in the stack for
	// later processing.
	Interval other;
	other.lower = mid;
	other.upper = I.upper;
	other.flips_lower  = flips_mid;
	other.flips_upper = I.flips_upper;
	S.push(other);

	// The local I interval will carry on computations for the left half.
	// NOTE: we perform these modifications after the data in I needed
	// for the right interval computations has already been used.
	I.upper = mid;
	I.flips_upper = flips_mid;
      }

      // Found a root.
      roots.push_back(I.lower);
    }
    
  }

  void Polynomial::operator=(const Polynomial &other)
  {
    Polynomial copy = other.copy();
    load(copy);
  }

  // Computes the entire Sturm Chain for use in root finding.
  void Polynomial::computeSturmChain(SturmChain & chain) const
  {
    Polynomial p1 = copy();
    chain.push_back(p1);

    if(p1.isZero())
    {     
      return;
    }

    Polynomial p2 = p1. differentiate();
    chain.push_back(p2);

    while(!p2.isZero())
    {
      Polynomial remainder_3;


      remainder_3 = -(p1 % p2);
      
      chain.push_back(remainder_3);

      p1 = p2;
      p2 = remainder_3;
    }
  }
  
  // Returns the number of times the sign flips in the sturm chain.
  int Polynomial::sturm_sign_flips(SturmChain & sturm_chain,
				   double input) const
  {
    int flips = 0;
    int sign = 0;

    for(auto iter = sturm_chain.begin();
	iter != sturm_chain.end();
	++iter)
    {
      double val = iter -> eval(input);
      int val_sign = (val > 0.0) - (val < 0.0);

      if(val_sign == 0) continue;
      if(sign == 0)
      {
	sign = val_sign;
	continue;
      }

      // If signs are flipped add a flip count.
      if(val_sign * sign < 0)
      {
	flips++;
	sign = val_sign;
      }
    }

    return flips;    
  }

  std::ostream& operator<<( std::ostream& os, const Polynomial &poly)
  {
    int deg = poly.degree();

    if(deg < 0)
    {
      os << "0";
      return os;
    }

    for(int i = deg; i > 0; i--)
    {
      os << poly.getCoef(i);
      os << "x^" << i << " + ";
    }

    os << poly.getCoef(0);
    return os;
  }

}
