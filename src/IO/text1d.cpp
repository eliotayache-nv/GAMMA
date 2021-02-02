/*
* @Author: Eliot Ayache
* @Date:   2020-09-28 16:57:12
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-02-02 20:54:20
*/

#include "../simu.h"
#include "../mpisetup.h"
#include "../environment.h"
#include "../constants.h"
#include "../array_tools.h"

class Data
{
public:

  Data(char *line): line(line){

    #if LOCAL_SYNCHROTRON_ == DISABLED_
      sscanf(line, 
        "%le %d %d %le %le %le %le %le %le %le %le %le %le\n",
        &t, &nact, &i, 
        &x, &dx, &dlx, &rho, &vx, &p, &D, &sx, &tau, &trac
        );
    #else
      sscanf(line, 
        "%le %d %d %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
        &t, &nact,  &i, 
        &x, &dx, &dlx, &rho, &vx, &p, &D, &sx, &tau, &trac,
        &Sd, &gmin, &gmax
        );
    #endif

  }

  ~Data(){}
  
  char *line;
  int i, nact;
  double t, x, dx, dlx, rho, vx, p, D, sx, tau, trac, Sd, gmin, gmax;

  void toCell(Grid *g){

    Cell *c = &(g->Cinit[i]);

    c->G.x[x_] = x;
    c->G.dx[x_] = dx;
    c->computeAllGeom();

    c->S.prim[RHO] = rho;
    c->S.prim[VV1] = vx;
    c->S.prim[PPP] = p;
    c->S.prim[TR1] = trac;

    #if LOCAL_SYNCHROTRON_ == ENABLED_
      c->S.prim[GMN] = pow(rho,1./3.) / gmin;
      c->S.prim[GMX] = pow(rho,1./3.) / gmax;
      // if (gmax > 10) printf("blah %d %d %le %le\n", j, i, gmax, c->S.prim[GMX]);
    #endif

  }

}; 


static void openLastSnapshot(DIR* dir, vector<Data> *data, long int *it, double *t){

  struct dirent  *dirp;       
  vector<string>  files;     
  char strfile[512];
  char strFilePath[512] = "../results/Last/";
  char *addr;

  while ((dirp = readdir(dir)) != NULL) {
    std::string fname = dirp->d_name;
    if(fname.find("phys") != std::string::npos)
      files.push_back(fname);
  }
  std::sort(files.begin(), files.end());
  strcpy(strfile,files[files.size()-1].c_str());
  sscanf(strfile, "phys%ld.h5", it);

  addr = strcat(strFilePath, strfile);
  if (worldrank == 0) printf("resuming setup from file: %s | ", addr);

  FILE *snap = fopen(addr, "r");
  char line[512];
  fgets(line, sizeof(line), snap);
  while (fgets(line, sizeof(line), snap)) {
    Data datapoint(line);
    data->push_back(datapoint);
  }
  *t = data->at(0).t;
  if (worldrank == 0) printf("tstart = %le\n", *t);
  fclose(snap);
  closedir(dir);

}


static void reloadFromData(Grid* g, vector<Data> *data);
void Simu::reinitialise(DIR* dir){

  std::vector<Data> data;
  openLastSnapshot(dir, &data, &it, &t);
  loadParams(&par);
  grid.initialise(par);   // this is unchanged from IC startup
  reloadFromData(&grid, &data);
  mpi_distribute(&grid);
  grid.prepForRun();

  data.clear();
  std::vector<Data>().swap(data); // freeing memory  

}


void reloadFromData(Grid *g, vector<Data> *data){

  int ngst = g->ngst;
  for (int c = 0; c < (int) data->size(); ++c){

    (*data)[c].toCell(g);

    int nact = (*data)[c].nact;
    g->nact   = nact;
    g->ntrack = nact + 2*ngst;
    g->iRbnd  = nact + ngst;
  }

}


void Grid::printCols(int it, double t){

  #if SHOCK_DETECTION_ == ENABLED_
      for (int i = 0; i < ntrack; ++i){
        Ctot[i].S.prim[TR1+1] = (double) Ctot[i].isShocked;
      }
  #endif

  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Datatype cell_mpi = {0}; 
  // generate_mpi_cell(&cell_mpi);

  int size = nde_nax[MV];;
  int jsize = 1;
  Cell    *Cdump = array_1d<Cell>(nax[MV]);
  s_cell *SCdump = array_1d<s_cell>(nax[MV]);
  for (int i = 0; i < nde_nax[MV]; ++i){
    toStruct(Ctot[i], &SCdump[i]);
  }

  std::copy_n(&Ctot[0], size, &Cdump[0]);

  // for (int j = 1; j < worldsize; ++j){
  //   int o[NUM_D]; // origin
  //   MPI_Recv(      &sizes[j],        1,  MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   MPI_Recv(              o,    NUM_D,  MPI_INT, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   int i1 = o[F1];
  //   int i2 = o[MV];
  //   MPI_Recv(&SCdump[i1][i2], sizes[j], cell_mpi, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   MPI_Recv(      &jsize[j],        1,  MPI_INT, j, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   MPI_Recv(   &ntrackd[i1+ngst], jsize[j],  MPI_INT, j, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // }

  std::stringstream ss;
  ss << std::setw(10) << std::setfill('0') << it;
  std::string s = "../results/Last/phys" + ss.str() + ".out";
  const char* strfout = s.c_str();
  FILE* fout = fopen(strfout, "w");

  #if LOCAL_SYNCHROTRON_ == ENABLED_
  fprintf(fout, "t nact i x dx dlx rho vx p D sx tau trac Sd gmin gmax\n");
  #else
  fprintf(fout, "t nact i x dx dlx rho vx p D sx tau trac Sd\n");
  #endif
  for (int i = ngst; i < ntrack-ngst; ++i) {
    toClass(SCdump[i], &Cdump[i]);
    double lfac = Cdump[i].S.lfac();
    int nactd = ntrack-2*ngst;
    #if LOCAL_SYNCHROTRON_ == ENABLED_
      Cdump[i].radiation_apply_trac2gammae();
      fprintf(fout, "%1.15le %d %d %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le\n", 
        t,
        nactd,
        i-ngst,
        Cdump[i].G.x[x_],
        Cdump[i].G.dx[x_],
        Cdump[i].G.dl[x_],
        Cdump[i].S.prim[RHO],
        Cdump[i].S.prim[UU1]/lfac,
        Cdump[i].S.prim[PPP],
        Cdump[i].S.cons[DEN],
        Cdump[i].S.cons[SS1],
        Cdump[i].S.cons[TAU],
        Cdump[i].S.prim[TR1],
        Cdump[i].S.prim[TR1+1],
        Cdump[i].S.prim[GMN],
        Cdump[i].S.prim[GMX]);
    #else
      fprintf(fout, "%1.15le %d %d %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le\n", 
        t,
        nactd,
        i-ngst,
        Cdump[i].G.x[x_],
        Cdump[i].G.dx[x_],
        Cdump[i].G.dl[x_],
        Cdump[i].S.prim[RHO],
        Cdump[i].S.prim[UU1]/lfac,
        Cdump[i].S.prim[PPP],
        Cdump[i].S.cons[DEN],
        Cdump[i].S.cons[SS1],
        Cdump[i].S.cons[TAU],
        Cdump[i].S.prim[TR1],
        Cdump[i].S.prim[TR1+1]);
    #endif
  }
  fclose(fout);

  delete_array_1d<Cell>(Cdump);
  delete_array_1d<s_cell>(SCdump);

  // }else{
  //   int size  = nde_nax[F1] * nde_nax[MV];  // size includes MV ghost cells
  //   int jsize = nde_nax[F1]-2*ngst;
  //   MPI_Send( &size,     1,  MPI_INT, 0, 0, MPI_COMM_WORLD);
  //   MPI_Send(origin, NUM_D,  MPI_INT, 0, 1, MPI_COMM_WORLD);
  //   s_cell **SC = array_2d<s_cell>(nde_nax[F1],nde_nax[MV]);
  //   for (int j = 0; j < nde_nax[F1]; ++j){
  //     for (int i = 0; i < nde_nax[MV]; ++i){
  //       toStruct(Ctot[j][i], &SC[j][i]);
  //     }
  //   }
  //   MPI_Send(&SC[0][0],     size, cell_mpi, 0, 2, MPI_COMM_WORLD);
  //   MPI_Send(&jsize,     1,  MPI_INT, 0, 3, MPI_COMM_WORLD);
  //   MPI_Send(&ntrack[ngst], jsize, MPI_INT, 0, 4, MPI_COMM_WORLD);
  //   delete_array_2d<s_cell>(SC);
  // }

  // MPI_Barrier(MPI_COMM_WORLD);

}
