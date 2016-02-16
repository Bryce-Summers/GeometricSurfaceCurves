#include "UF_ADT.h"


UF_ADT::UF_ADT()
{

}

UF_ADT::~UF_ADT()
{

}


/* Resizes the F structure to the given size.
 * and disassociates all of the sets.
 */
 /*
void UF_ADT::makeSize(int size)
{

    free(this->parents);
    free(this->ranks);

    this->parents = (int*) malloc(sizeof(int)*size);
    this->ranks   = (int*) malloc(sizeof(int)*size);
    this->size = size;

    for(int i = 0; i < size; i++)
    {
        parents[i] = i;
        ranks[i]   = 0;
    }

}
*/

