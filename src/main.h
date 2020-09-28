#ifndef MAIN_H_
#define MAIN_H_

#include "environment.h"

class Flags
{

public:
    
  Flags(){

    overwrite = false;
    resume = false;

  }

  ~Flags(){}

  bool overwrite;
  bool resume;

  void load(int argc, char *argv[]);
  void checkApplicable();

};


#endif