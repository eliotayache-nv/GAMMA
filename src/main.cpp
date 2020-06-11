/*
* @Author: eliotayache
* @Date:   1020-05-05 10:06:26
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2020-06-11 15:34:22
*/

#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>
#include <unistd.h>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#include "environment.h"
#include "simu.h"
#include "grid.h"
#include "err.h"
#include "mpisetup.h"
#include "array_tools.h"

void chbindir(const char binPath[]);

int main(int argc, char *argv[])
{
    c_simu simu;

    chbindir(argv[0]);

    mpi_init(argc, argv);
    simu.initialise();

    // hydro.loadConfig();
    // grid.create();
    // grid.initialState();

    // int** arr;
    // arr = array_2d_dyn<int>(10,10,120);

    // for (int i = 0; i < 10; ++i)
    // {
    //     for (int j = 0; j < 10; ++j)
    //     {
    //         arr[i][j] = 10 * i + j;
    //     }
    // }

    // for (int i = 0; i < 10; ++i)
    // {
    //     for (int j = 0; j < 10; ++j)
    //     {
    //         printf("%2d ", arr[i][j]);
    //     }
    //     printf("\n");
    // }

    return 0;
}

void chbindir(const char binPath[])
{
    char pathCopy[300];
    int status = 0;

    strcpy(pathCopy, binPath); 
    status = chdir(dirname(pathCopy));
    if (status == -1) throw StartupFileManagementException();
    return;
}