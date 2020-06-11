/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:34:42
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-11 18:35:29
*/

#include "environment.h"

int worldsize, worldrank, nodesize, noderank;
MPI_Comm nodecom;

void checkEnvironment()
{
    if (NUM_DIM==1 and worldsize!=1) throw inconsistentEnvironmentException();
}
