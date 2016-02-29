#ifndef TESTING_H
#define TESTING_H

/*
 * Testing Class for My classes.
 * Written by Bryce Summers on 2/26/2016.
 *
 * Tests:
 *  - Polynomials.
 */

#include "CMU462/CMU462.h"
#include <vector>
#include <iostream>

namespace CMU462
{

  class Tester
  {

  public:
    Tester(){};
    ~Tester(){};

    void test();
    void test_polynomials();
  };
}

#endif // TESTING_H
