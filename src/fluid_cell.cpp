#include "fluid_cell.h"
#include "environment.h"

c_fluid_cell::c_fluid_cell()
{
    normalized = false;
    #if ndim_ == 1
        G.u[x_] = 0.;
    #endif
    merged = false;
    split = false;
    isShocked = false;
    FshockCount = 0;
    RshockCount = 0;

    deriv_v  = 0.;
    deriv2_v = 0.;
}

c_fluid_cell::~c_fluid_cell()
{

}


void c_fluid_cell :: cons2rsqrd()
{
    S.D *= G.pos[x_] * G.pos[x_];
    S.m *= G.pos[x_] * G.pos[x_];
    S.E *= G.pos[x_] * G.pos[x_];

    S.rho *= G.pos[x_] * G.pos[x_];
    S.v *= G.pos[x_] * G.pos[x_];
    S.p *= G.pos[x_] * G.pos[x_];
    S.u *= G.pos[x_] * G.pos[x_];
}


void c_fluid_cell :: rsqrd2cons()
{
    S.D /= G.pos[x_] * G.pos[x_];
    S.m /= G.pos[x_] * G.pos[x_];
    S.E /= G.pos[x_] * G.pos[x_];

    S.rho /= G.pos[x_] * G.pos[x_];
    S.v /= G.pos[x_] * G.pos[x_];
    S.p /= G.pos[x_] * G.pos[x_];
    S.u /= G.pos[x_] * G.pos[x_];
}


void c_fluid_cell::assignInterfaces(
    c_cell_interface *pIL,
    c_cell_interface *pIR)
{
    G.pIL = pIL;
    G.pIR = pIR;
}














