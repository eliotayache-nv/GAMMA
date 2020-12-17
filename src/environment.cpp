/*
* @Author: Eliot Ayache
* @Date:   2020-06-11 18:34:42
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-12-17 10:44:10
*/

#include "environment.h"

int worldsize, worldrank, nodesize, noderank;
MPI_Comm nodecom;

void radiation_InitialChecks();

void checkEnvironment()
{
  if (NUM_D==1 and worldsize!=1) throw inconsistentEnvironmentException();

  #if LOCAL_SYNCHROTRON_ == ENABLED_
    radiation_InitialChecks();
  #endif
}

void radiation_InitialChecks(){

  printf("Local Synchrotron calculation ENABLED\n");
  if (SHOCK_DETECTION_ != ENABLED_){
    printf("Shock detection needs to be ENABLED!\n");
    exit(30);
  }
  
}

