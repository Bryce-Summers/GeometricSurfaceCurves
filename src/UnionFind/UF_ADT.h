#ifndef UF_ADT_H
#define UF_ADT_H

#include <stdlib.h> // Needed for mallocing.

/*
 * Written by Bryce Summers on 4/12/2015.
 *
 * This class defines the abstract data type of a Union Find structure.
 *
 * These structures will support operations on fixed contiguous integer valued indices.
 * The interpretation of the results of these union and find operations will be left to the code utilizing the structure.
 */

class UF_ADT
{
    public:

        // Should only be called once!
        UF_ADT();

	virtual ~UF_ADT() = 0;

        // This function should be overriden if a parallel implementation
        // wishes to parallelize the insertion of edges.
        // Returns true iff the two vertices were not connected initially and are now connected.
        virtual bool op_union(int v1, int v2) = 0;

        // Guaranteed to be valid for sequential implementations.
        // For concurrent usage using comparison checks, please use connected instead.
        virtual int op_find(int vertex) = 0;


        /* ENSURES:
         * If the two vertices are connected at the beginning --> true.
         * disconnected at beginning, but connected mid operation --> undefined. //FIXME : IDEAL : true.
         * disconnected at beginning, disconnected at end --> false.
         */
        virtual bool connected (int v1, int v2) = 0;

	//virtual int getSize();

    protected:
    private:
};

#endif // UF_STRUCTURE_H
