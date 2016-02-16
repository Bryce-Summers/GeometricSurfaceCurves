#ifndef UF_SERIAL_H
#define UF_SERIAL_H

#include "UF_ADT.h"
/*
 * Written By Bryce Summers on 4/12/2015.
 *
 * Purpose : This is a standard serial implementation of a Union Find Data structure.
 *
 * Updated: 4/19/2015 : Converted find from a recursive to an iterative implementation.
 */

class UF_Serial : public UF_ADT
{
    public:
        UF_Serial(int size);
        ~UF_Serial();

	// void op_union(EdgeList edgeList);
        bool op_union(int v1, int v2);
        int op_find(int vertex);
        bool connected(int v1, int v2);

    protected:
    private:

        // Links two root elements.
	// REQUIRES : v1 and v2 should be root nodes.
	// Returns true iff the two root nodes are not equal and have been linked together.
        bool link(int v1, int v2);

    // The Data.
    int * parents;
    int * ranks;
    int size;
};

#endif // UF_SERIAL_H
