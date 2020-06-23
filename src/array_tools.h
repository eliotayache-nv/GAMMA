//////////////////////////////////////////////////////////////////////////////////
// arraytools.h
//
// Created pre 28 July, 2010 by Weiqun Zhang
// Last modified August 2, 2010 by HJvE
//
// When creating multidimensional arrays by manually creating a list of pointers
// the array will not occupy one continuous chunk of memory. The routines here
// provide a workaround. 
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ARRAY_TOOLS_
#define ARRAY_TOOLS_

#include <stdlib.h>
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// prior declarations, to prevent spurious compiler warnings

template <class T> T* array_1d(const int ni);
template <class T> T** array_2d(const int ni, const int nj);
template <class T> T*** array_3d(const int ni, const int nj, const int nk);
template <class T> void delete_array_1d(T* p);
template <class T> void delete_array_2d(T** pp);
template <class T> void delete_array_3d(T*** ppp);
template <class T> void arrcpy(T* src, T* dest, int len);

////////////////////////////////////////////////////////////////////////////////

template <class T> void arrcpy(T* src, T* dest, int len){
  
  for (int i = 0; i < len; ++i){
    dest[i] = src[i];
  }

}

////////////////////////////////////////////////////////////////////////////////

template <class T> T* array_1d(const int ni)
{
  T* p;

  try {
    p = new T[ni];
  }
  catch (bad_alloc&) {
    cout<<"arr_1d: memory allocation error, p = new T["<<ni<<"]"<<endl;
    exit(1);
  }
  return p;
}

////////////////////////////////////////////////////////////////////////////////

template <class T> T** array_2d(const int ni, const int nj)
{
  T** pp;

  try {
    pp = new T*[ni];
  }
  catch (bad_alloc&) {
    cout<<"arr_2d: memory allocation error, pp = new T*["<<ni<<"]"<<endl;
    exit(1);
  }

  try {
    pp[0] = new T[ni*nj];
  }
  catch (bad_alloc&) {
    cout<<"arr_2d: memory allocation error, pp[0] = new T["<<ni*nj<<"]"<<endl;
    exit(1);
  }

  for (int i=1; i<ni; ++i) { pp[i] = pp[0] + i * nj; }

  return pp;
}

template <class T> T** array_2d_nogst(T** in, int nj, const int ngst){

  T** pp;

  int nj_new = nj-2*ngst;
  try {
    pp = new T*[nj_new];
  }
  catch (bad_alloc&) {
    cout<<"arr_2d: memory allocation error, pp = new T*["<<nj_new<<"]"<<endl;
    exit(1);
  }

  for (int j=0; j<nj_new; ++j) { pp[j] = in[ngst+j] + ngst; }

  return pp;

}

template <class T> T** array_2d_dyn(const int ni, const int nj, const int nmax)
{
  T** pp;

  try {
    pp = new T*[ni];
  }
  catch (bad_alloc&) {
    cout<<"arr_2d: memory allocation error, pp = new T*["<<ni<<"]"<<endl;
    exit(1);
  }

  try {
    pp[0] = new T[nmax];
  }
  catch (bad_alloc&) {
    cout<<"arr_2d: memory allocation error, pp[0] = new T["<<ni*nj<<"]"<<endl;
    exit(1);
  }

  for (int i=1; i<ni; ++i) { pp[i] = pp[0] + i * nj; }

  return pp;
}

////////////////////////////////////////////////////////////////////////////////

template <class T> T*** array_3d(const int ni, const int nj, const int nk)
{
  T*** ppp;

  try {
    ppp = new T**[ni];
  }
  catch (bad_alloc&) {
    cout<<"arr_3d: memory allocation error, ppp = new T**["<<ni<<"]"<<endl;
    exit(1);
  }

  try {
    ppp[0] = new T*[ni*nj];
  }
  catch (bad_alloc&) {
    cout<<"arr_3d: memory allocation error, ppp[0] = new T*["<<ni*nj<<"]"<<endl;
    exit(1);
  }
  ppp[0][0] = new T[ni*nj*nk];
  if (ppp[0][0] == 0) {
    cout<<"arr_3d: memory allocation error, ppp[0][0] = new T["
      <<ni*nj*nk<<"]"<<endl;
    exit(1);
  }

  for (int i=1; i < ni; ++i) { ppp[i] = ppp[0] + i * nj; }
  
  // loop over 2D entries to specify pointers. We need to avoid the 0,0 entry  
  for (int i=1; i < ni; ++i) 
    for (int j=0; j < nj; ++j)
      ppp[i][j] = ppp[0][0] + i * nj * nk + j * nk ;
   
  for (int j=1; j < nj; j++)
    ppp[0][j] =  ppp[0][0] + j * nk;
 
  return ppp;
}

////////////////////////////////////////////////////////////////////////////////

template <class T> void delete_array_1d(T* p)
{
  delete [] p;
}

template <class T> void delete_array_2d(T** pp)
{
  delete [] pp[0];
  delete [] pp;
}

template <class T> void delete_array_3d(T*** ppp)
{
  delete [] ppp[0][0];
  delete [] ppp[0];
  delete [] ppp;
}

#endif // ARRAY_TOOLS_
