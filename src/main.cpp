/*
* @Author: eliotayache
* @Date:   1020-05-05 10:06:26
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-09-28 18:59:27
*/

#include "main.h"
#include "environment.h"
#include "simu.h"
#include "err.h"
#include "mpisetup.h"


static void chbindir(const char binPath[]);


int main(int argc, char *argv[]){

  int status;
  Simu simu;
  Flags flags;

  flags.load(argc, argv);
  flags.checkApplicable();

  chbindir(argv[0]);
  mpi_init(&argc, &argv);
  checkEnvironment();

  if (flags.resume){
    DIR* dir = opendir("../results/Last");
    simu.reinitialise(dir);    
  }
  else
    simu.initialise();

  simu.run();

  mpi_finalise();
  return 0;

}


void Flags::load(int argc, char *argv[]){

  if (argc == 1) return; // no command-line arguments passed
  for (int a = 0; a < argc; ++a){
    if (!strncmp(argv[a],"-r",2))
      resume = true;
    if (!strncmp(argv[a],"-w",2))
      overwrite = true;
  }

}


// Check that the combination of flags is allowed
void Flags::checkApplicable(){

  if (!overwrite){
    DIR* dir = opendir("../results/Last");
    if (dir) throw OverwriteOutputException();
  }

  if (overwrite and resume){
    throw WrongCombinatinoOfFlagsException();
  }

}



static void chbindir(const char binPath[]){

  char pathCopy[300];
  int status = 0;

  strcpy(pathCopy, binPath); 
  status = chdir(dirname(pathCopy));
  if (status == -1) throw StartupFileManagementException();
  return;
  
}