/*
* @Author: Eliot Ayache
* @Date:   2020-06-24 14:51:12
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-24 15:38:05
*/

#include <stdio.h>


class State
{
public:
    State(){}
    ~State(){}
    
    int m;
    void statefunc();

};

void State::statefunc(){

    m*=10;

}


class Cell{

public:
    Cell(){}
    ~Cell(){}
    
    int n;
    State S;

    void cellfunc();

};

void Cell::cellfunc(){

    n+=1;

}


class Grid
{
public:
    Grid(){}
    ~Grid(){}

    Cell c;

    template <class T> void apply(void (T::*func)());
    template <> void apply<State>(void (State::*func)()); 

};


template <class T> void Grid::apply(void (T::*func)()){
    
    (c.*func)();

}

template <> void Grid::apply<State>(void (State::*func)()){

    (c.S.*func)();
    
}


int main(int argc, char const *argv[])
{
    Grid g;
    g.c.n = 1;
    g.c.S.m = 2;

    g.apply(&Cell::cellfunc);
    g.apply(&State::statefunc);

    printf("%d %d\n",g.c.n, g.c.S.m);

    return 0;
}













