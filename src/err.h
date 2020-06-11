#ifndef ERR_H_
#define ERR_H_ 

#include <iostream>
#include <exception>

// Global variables
// Debug variables:
extern bool     g_printErr;         // flag to print or not
extern bool     g_lockPrintErr;     // lock printErr if verbose is 0 or 2

// Global functions
void setPrintErr(bool b);
void printErr(const char* format,...); 


using namespace std;

struct MPITooManyTasksException : public exception
{
  const char * what () const throw ()
  {
    return "Too many MPI tasks, Number of tasks must be equal to number of nodes.";
  }
};

struct NegativePressureException : public exception
{
  const char * what () const throw ()
  {
    return "Negative Pressure Exception!";
  }
};

struct PrimitiveReconstructionDomainException : public exception
{
  const char * what () const throw ()
  {
    return "Pressure ran out of definition domain in primitive reconstruction!";
  }
};

struct PrimitiveReconstructionFailedException : public exception
{
  const char * what () const throw ()
  {
    return "Timeout on primitive reconstruction!";
  }
};

struct StartupFileManagementException : public exception
{
  const char * what () const throw ()
  {
    return "Error when handling files on startup";
  }
};


#endif