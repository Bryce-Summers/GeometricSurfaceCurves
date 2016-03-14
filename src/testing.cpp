#include "testing.h"

namespace CMU462
{
  using namespace std;

  void Tester::test()
  {
    test_polynomials();
  }
  
  void Tester::test_polynomials()
  {
    Polynomial ZERO = Polynomial::ZERO();
    cout << ZERO << endl;

    Polynomial p = ZERO.copy();
    
    p.addTrailingTerm(3);
    p.addTrailingTerm(-5);
    p.addTrailingTerm(10);
    p.addTrailingTerm(-3);

    Polynomial p2 = ZERO;//ZERO.copy();//new Polynomial();
    p2.addTrailingTerm(3);
    p2.addTrailingTerm(1);
    //    p2 -> addTrailingTerm(10);

    cout << "a:   " << p  << endl;
    cout << "b:   " << p2 << endl;
    
    Polynomial p12 = p + p2;
    cout << "a+b: " << p12 << endl;

    Polynomial sub = p - p2;
    cout << "a-b: " << sub << endl;

    Polynomial scale = p *.5;
    cout << "a*.5: " << scale << endl;

    Polynomial mult = p*p2;
    cout << "a*b: " << mult << endl;

    // Division and Remainder.
    Polynomial div = p/p2;
    cout << "a/b: " << p/p2  << endl;
    //    cout << "a%b: " << p%p2 << endl;

    // WARNING: Memory Explosion if not handled with care.
    
    std::vector<double> roots;

    /*
    p -> computeRealRoots(roots,
			  -1000, 1000,
			  .000001);
    */

    Polynomial p3 = ZERO;
    // This polynomial has roots 1, 2, 3, 4, 5. They are all found accuratly.
    
    p3.addTrailingTerm(1);
    p3.addTrailingTerm(-15);
    p3.addTrailingTerm(85);
    p3.addTrailingTerm(-225);
    p3.addTrailingTerm(274);
    p3.addTrailingTerm(-120);
    

    // This polynomial has roots:
    // -1.32, 1.32, 4.56, 4.58, 4.6.
    // Thus far we only find -1.32, 1.32, and 4.53658.
    // We don't find the close last 3 values accuratly.
    /*
    p3.addTrailingTerm(1);
    p3.addTrailingTerm(-13.74);
    p3.addTrailingTerm(61.1864);
    p3.addTrailingTerm(-72.1295);
    p3.addTrailingTerm(-109.647);
    p3.addTrailingTerm(167.393);
    */

    std::cout << "Quintic Polynomial: " << p3 << std::endl;
    
    p3.computeRealRoots(roots,
			-1000, 1000,
			.0000001);

    for(auto iter = roots.begin(); iter != roots.end(); ++iter)
    {
      cout << "Root: " << *iter << " -> " << p3.eval(*iter)<< endl;
    }

    Polynomial x, y, z;

    x.addTrailingTerm(1);
    x.addTrailingTerm(1);

    y.addTrailingTerm(2);
    y.addTrailingTerm(2);

    z.addTrailingTerm(3);
    z.addTrailingTerm(3);

    PolynomialVector3D vec(x, y, z);

    cout << "vec - vec = " << vec - vec << endl;
    cout << "vec + vec = " << vec + vec << endl;
    
    cout << vec << endl;
    vec = dot(vec, vec);
    cout << vec << endl;
  }
}
