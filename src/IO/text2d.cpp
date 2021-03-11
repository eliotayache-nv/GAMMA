/*
* @Author: Eliot Ayache
* @Date:   2020-09-28 16:57:12
* @Last Modified by:   Eliot Ayache
* @Last Modified time: 2021-03-11 15:08:09
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
        "%le %d %d %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
        &t, &nact, &j, &i, 
        &x, &y, &dx, &dy, &dlx, &dly, &rho, &vx, &vy, &p, &D, &sx, &sy, &tau, &trac
        );
    #else
      sscanf(line, 
        "%le %d %d %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
        &t, &nact, &j, &i, 
        &x, &y, &dx, &dy, &dlx, &dly, &rho, &vx, &vy, &p, &D, &sx, &sy, &tau, &trac,
        &Sd, &gmin, &gmax, &gmid, &psyn
        );
    #endif

  }

  ~Data(){}
  
  char *line;
  int j, i, nact;
  double t, x, y, dx, dy, dlx, dly, rho, vx, vy, p, D, sx, sy, tau, trac, Sd, gmin, gmax, gmid, psyn;

  void toCell(Grid *g){

    Cell *c = &(g->Cinit[j][i]);

    c->G.x[x_] = x;
    c->G.x[y_] = y;
    c->G.dx[x_] = dx;
    c->G.dx[y_] = dy;
    c->computeAllGeom();

    c->S.prim[RHO] = rho;
    c->S.prim[VV1] = vx;
    c->S.prim[VV2] = vy;
    c->S.prim[PPP] = p;
    c->S.prim[TR1] = trac;

    #if LOCAL_SYNCHROTRON_ == ENABLED_
      c->S.prim[GMN] = pow(rho,1./3.) / gmin;
      c->S.prim[GMX] = pow(rho,1./3.) / gmax;
      c->S.prim[GMD] = pow(rho,1./3.) / gmid;
      c->S.prim[GMX] = psyn;
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
  int o = g->origin[F1];
  int e = o + g->nde_ncell[F1]-1;
  for (int c = 0; c < (int) data->size(); ++c){
    int j = (*data)[c].j;

    (*data)[c].toCell(g);

    int nact = (*data)[c].nact;
    int jtrgt = j-o+ngst;
    if (j>=o and j<=e){
      g->nact[jtrgt]   = nact;
      g->ntrack[jtrgt] = nact + 2*ngst;
      g->iRbnd[jtrgt]  = nact + ngst;
    }
  }

}


void Grid::printCols(int it, double t){

  #if SHOCK_DETECTION_ == ENABLED_
    for (int j = 0; j < nde_nax[F1]; ++j){
      for (int i = 0; i < ntrack[j]; ++i){
        Ctot[j][i].S.prim[TR1+1] = (double) Ctot[j][i].Sd;
      }
    }
  #endif

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Datatype cell_mpi = {0}; 
  generate_mpi_cell(&cell_mpi);

  if (worldrank == 0){

    int sizes[worldsize];
    int jsize[worldsize];
    int ntrackd[ncell[F1]+2*ngst];
    Cell    **Cdump = array_2d<Cell>(nax[F1], nax[MV]);   // nax = ncell+2ngst (includes ghosts)
    s_cell **SCdump = array_2d<s_cell>(nax[F1], nax[MV]);
    for (int j = 0; j < nde_nax[F1]; ++j){
      for (int i = 0; i < nde_nax[MV]; ++i){
        toStruct(Ctot[j][i], &SCdump[j][i]);
      }
    }

    // bluntly copying Ctot into lower part of Cdump
    sizes[0] = nde_nax[F1] * nde_nax[MV];
    jsize[0] = nde_nax[F1];
    std::copy_n(&Ctot[0][0], sizes[0], &Cdump[0][0]);
    std::copy_n(&ntrack[0], nde_nax[F1], &ntrackd[0]);  

    for (int j = 1; j < worldsize; ++j){
      int o[NUM_D]; // origin
      MPI_Recv(      &sizes[j],        1,  MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(              o,    NUM_D,  MPI_INT, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int i1 = o[F1]+ngst;
      int i2 = o[MV];
      MPI_Recv(&SCdump[i1][i2], sizes[j], cell_mpi, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(      &jsize[j],        1,  MPI_INT, j, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(   &ntrackd[i1], jsize[j],  MPI_INT, j, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    std::stringstream ss;
    ss << std::setw(10) << std::setfill('0') << it;
    std::string s = "../results/Last/phys" + ss.str() + ".out";
    const char* strfout = s.c_str();
    FILE* fout = fopen(strfout, "w");

    #if LOCAL_SYNCHROTRON_ == ENABLED_
      fprintf(fout, "t nact j i x y dx dy dlx dly rho vx vy p D sx sy tau trac Sd gmin gmax gmid psyn\n");
    #else
      fprintf(fout, "t nact j i x y dx dy dlx dly rho vx vy p D sx sy tau trac Sd\n");
    #endif
    for (int j = ngst; j < ncell[F1]+ngst; ++j){
      for (int i = ngst; i < ntrackd[j]-ngst; ++i) {
        toClass(SCdump[j][i], &Cdump[j][i]);
        double lfac = Cdump[j][i].S.lfac();
        int nactd = ntrackd[j]-2*ngst;
        #if LOCAL_SYNCHROTRON_ == ENABLED_
          Cdump[j][i].radiation_apply_trac2gammae();
          fprintf(fout, "%1.15le %d %d %d %1.15le %1.15le %1.15le %1.15le %1.15le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n", 
            t,
            nactd,
            j-ngst,
            i-ngst,
            Cdump[j][i].G.x[x_],
            Cdump[j][i].G.x[y_],
            Cdump[j][i].G.dx[x_],
            Cdump[j][i].G.dx[y_],
            Cdump[j][i].G.dl[x_],
            Cdump[j][i].G.dl[y_],
            Cdump[j][i].S.prim[RHO],
            Cdump[j][i].S.prim[UU1]/lfac,
            Cdump[j][i].S.prim[UU2]/lfac,
            Cdump[j][i].S.prim[PPP],
            Cdump[j][i].S.cons[DEN],
            Cdump[j][i].S.cons[SS1],
            Cdump[j][i].S.cons[SS2],
            Cdump[j][i].S.cons[TAU],
            Cdump[j][i].S.prim[TR1],
            Cdump[j][i].S.prim[TR1+1],
            Cdump[j][i].S.prim[GMN],
            Cdump[j][i].S.prim[GMX],
            Cdump[j][i].S.prim[GMD],
            Cdump[j][i].S.prim[PSN]);
        #else
          fprintf(fout, "%1.15le %d %d %d %1.15le %1.15le %1.15le %1.15le %1.15le %1.15le %le %le %le %le %le %le %le %le %le %le\n", 
            t,
            nactd,
            j-ngst,
            i-ngst,
            Cdump[j][i].G.x[x_],
            Cdump[j][i].G.x[y_],
            Cdump[j][i].G.dx[x_],
            Cdump[j][i].G.dx[y_],
            Cdump[j][i].G.dl[x_],
            Cdump[j][i].G.dl[y_],
            Cdump[j][i].S.prim[RHO],
            Cdump[j][i].S.prim[UU1]/lfac,
            Cdump[j][i].S.prim[UU2]/lfac,
            Cdump[j][i].S.prim[PPP],
            Cdump[j][i].S.cons[DEN],
            Cdump[j][i].S.cons[SS1],
            Cdump[j][i].S.cons[SS2],
            Cdump[j][i].S.cons[TAU],
            Cdump[j][i].S.prim[TR1],
            Cdump[j][i].S.prim[TR1+1]
            );
        #endif
      }
    }
    fclose(fout);

    delete_array_2d<Cell>(Cdump);
    delete_array_2d<s_cell>(SCdump);

  }else{
    int size  = (nde_nax[F1]-2*ngst) * nde_nax[MV];  // size includes MV ghost cells
    int jsize = nde_nax[F1]-2*ngst;
    MPI_Send( &size,     1,  MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(origin, NUM_D,  MPI_INT, 0, 1, MPI_COMM_WORLD);
    s_cell **SC = array_2d<s_cell>((nde_nax[F1]-2*ngst),nde_nax[MV]);
    for (int j = 0; j < nde_nax[F1]-2*ngst; ++j){
      for (int i = 0; i < nde_nax[MV]; ++i){
        toStruct(Ctot[j+ngst][i], &SC[j][i]);
      }
    }
    MPI_Send(&SC[0][0],      size, cell_mpi, 0, 2, MPI_COMM_WORLD);
    MPI_Send(&jsize,            1,  MPI_INT, 0, 3, MPI_COMM_WORLD);
    MPI_Send(&ntrack[ngst], jsize,  MPI_INT, 0, 4, MPI_COMM_WORLD);
    delete_array_2d<s_cell>(SC);
  }

  MPI_Barrier(MPI_COMM_WORLD);

}
