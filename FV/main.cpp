/* 2D Navier Stokes solver for a structured quadrilateral grid
 *
 * T. Sailor Koeplinger - Feb2024
 */
// C/C++ libraries
#include <iostream>
#include <cmath>
#include <unistd.h>

// External Libraries
#include <mpi.h>
#include <petscts.h>

// Project Headers
#include "FileIO.h"
#include "Indexing.h"
#include "SpatialDiscretization.h"
#include "MeshModule.h"
#include "StateVariables.h"
#include "LUtools.h"
#include "Jacobian.h"
#include "Thermo.h"

  // List of context needed: air,nx,ny,CFL,Elemvar, geofa
typedef struct {
    Thermo          air;
    PetscInt        nx;
    PetscInt        ny;     //////   	MAY HAVE ISSUES WITH CONVERSION TO REGULAR INT  
    PetscScalar     CFL;
    *State          ElemVar;
    *double         geofa;
} AppCtx	


double find_dt(Thermo& air, int nx, int ny, double CFL, const double* uRef, State& var, double* geofa){
    double vmax, dt, mindx;

    //Max info speed = v+c
    vmax = sqrt(var.v2) + var.a;

    //Min grid spacing
    mindx = 1.0;
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            mindx = fmin(mindx, geofa[IJK(i,j,0,nx,6)]);
            mindx = fmin(mindx, geofa[IJK(i,j,3,nx,6)]);
        }
    }

    dt = CFL * mindx / vmax;
    ASSERT(!__isnan(dt), "NAN dt")
    return dt;
}

void calculate_residual(int nx, int ny, double* res, double* ressum){
    double res2[NVAR] = {0.0};

    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            int iu = IJK(i,j,0,nx-1,NVAR);

            for (int k=0; k<NVAR; k++){
                ASSERT(!__isnan(res[iu+k]),"res NaN")
                res2[k] += res[iu+k]*res[iu+k];
            }
        }
    }

    ressum[0] = sqrt(res2[0]);
    ressum[1] = sqrt(res2[1]);
    ressum[2] = sqrt(res2[2]);
    ressum[3] = sqrt(res2[3]);
}
// ======================================================================================
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~==========
// ======================================================================================
PetscErrorCode pre_step_routine(TS ts) {
  AppCtx            *app_ctx;
  Vec               U;
  PetscReal         dt, norm;
  MPI_Comm          comm;
  //
  PetscFunctionBeginUser;
  comm = PETSC_COMM_WORLD;
  //
  int world_size,bnum;
  double dt, dt_glob;
  //
  PetscCallMPI(MPI_Comm_size(comm, &world_size));
  PetscCallMPI(MPI_Comm_rank(comm, &bnum));
  //
  PetscCall(TSGetApplicationContext(ts, &app_ctx));
  //
  // List of context needed: air,nx,ny,CFL,Elemvar, geofa
  //Get needed values from ts context
  PetscCall(TSGetTimeStep(ts, &dt));
  PetscCall(TSGetSolution(ts, &U));
  // get local solution vector
  PetscCall(VecGetArray(U, &unk));
  //
  // ==========================================================================
  // Compute Timestep
  // 
  // Find global timestep based off of CFl condition
  dt = find_dt(app_ctx->air, app_ctx->nx, app_ctx->ny ,app_ctx->CFL, 
	       unk,          app_ctx->ElemVar[0],     app_ctx->geofa);
  if(DEBUG) {printf("::%3d::Calculated Timestep..... \n", bnum);}
  // 
  // Get the smalest dt out of all processes
  PetscCall(MPI_Allreduce(&dt,&dt_glob, 1, MPI_DOUBLE, MPI_MIN, comm));
  PetscCall(TSSetTimeStep(ts, dt_glob));
  if(DEBUG) {printf("::%3d::Communicated Timestep..... \n", bnum);}
  //
  // ==========================================================================
  //
  // Cleanup, cleanup, everybody everywhere....
  PetscCall(VecRestoreArray(U, &unk));
  PetscFunctionReturn(PETSC_SUCCESS);
  //
}
//
// ======================================================================================
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~==========
// ======================================================================================
//
PetscErrorCode calc_rhs(TS ts, PetscReal t, Vec u, Vec Fv, void *ac) {
    // INPUS
    //       ts  :: timestepping context struct
    //       t   :: solution time
    //       u   :: solution vector (global)
    //       F   :: RHS vector (global)
    //       ctx :: user specified context (big boy)
    Vec               Fu
    PetscScalar       dt, norm;
    PetscScalar       *unk, *res, *dv;
    PetscInt          world_size,bnum;
    MPI_Comm          comm;
    //
    PetscFunctionBeginUser;
    comm = PETSC_COMM_WORLD;
    //
    double ressum;
    //
    PetscCallMPI(MPI_Comm_size(comm, &world_size));
    PetscCallMPI(MPI_Comm_rank(comm, &bnum));
    //
    PetscCall(TSGetApplicationContext(ts, &ac));
    // Need for AppCtx: ivisc, accur, iaxi, mxangle, bbounds, bids, nx, ny, air, elemvar
    //                  uFSm ibound, geoel, geofa, yfa, xfa, unk,
    //                  ux(tmp nullptr), uy(same), resx(same), resy(same)
    //
    PetscCall(VecDuplicate(Fv,&Fu));
    //
    // Step1 calculate dudt
    //.... Separate out input (unk), and output (res)   
    PetscCall(VecGetArray(U,  &unk));
    PetscCall(VecGetArray(Fu, &res));
    PetscCall(VecGetArray(Fv, &dv ));
    //
    calc_dudt(ac->ivisc, ac->accur,  ac->iaxi,  ac->mxangle, ac->bbounds,
              ac->bids,  ac->nx,     ac->ny,    ac->air,     ac-> ElemVar,
	      ac->uFS,   ac->ibound, ac->geoel, ac->geofa,   ac->yfa,  ac->xfa,
	      /**/unk,   ac->ux,     ac->uy,    /**/res,     ac->resx, ac->resy);
    if(DEBUG) {printf("::%3d::Calculated dudt..... \n", bnum);}
    //calculate the right hand side residual term (change of conserved quantities)
    calculate_residual(nx, ny, res, ressum);
    /*
        if ( accur == 1 ) {
            calculate_residual(nx, ny, resx, ressumx);
            calculate_residual(nx, ny, resy, ressumy);
        }
    */
    //
    //========== Solve linear system on each element (turns chg in conservatives to change in solution variables)
    // I'm sure there's a better way to do this...
    auto D = (double**)malloc((NVAR) * sizeof(double*));
    for (int k = 0; k < NVAR; k++)
        D[k] = (double*)malloc( (NVAR) * sizeof(double));
    //
    //PetscCall(TSGetTimeStep(ts, &dt));
    //
    for (int i=0; i<(ac->nx-1); i++) {
        for (int j=0; j<(ac->ny-1); j++) {
            double *unkij = &(unk[IJK(i, j, 0, nx-1, NVAR)]);
            double LUtol = 1e-16;
            int iel = IJ(i,j,nx-1);
            int N = NVAR;
            int P[NVAR+1]{}; //permutation vector for pivoting

            //get the rhs block needed
            double *b = &(res[IJK(i, j, 0, nx - 1, NVAR)]);
            double *xLU = &(dv[ IJK(i, j, 0, nx - 1, NVAR)]);

            //Evaluate the jacobian / Implicit matrix
            BuildJacobian(1.0, unkij, ac->ElemVar[iel], D);
            LUPDecompose(D, N, LUtol, P);
            LUPSolve(D, P, b, N, xLU);

            //axi modification
            double ycc = 1.0;
            if (iaxi==1) {
                ycc = ac->geoel[IJK(i, j, 2, nx - 1, 3)];
                for (int k = 0; k < NVAR; k++) {
                    xLU[k] *= (1.0 / ycc);
                }
            }
    /*   Additions for 2nd order
            if (accur == 1) {
                ////  CURRENTLY USING THE CELL CENTERED JACOBIAN VAL
                double *bx = &(resx[IJK(i, j, 0, nx - 1, NVAR)]);
                double *by = &(resy[IJK(i, j, 0, nx - 1, NVAR)]);
                double *xLUx = &(dvx[IJK(i, j, 0, nx - 1, NVAR)]);
                double *xLUy = &(dvy[IJK(i, j, 0, nx - 1, NVAR)]);
                LUPSolve(D, P, bx, N, xLUx);
                LUPSolve(D, P, by, N, xLUy);

                //axi modification
                if (iaxi==1) {
                    for (int k = 0; k < NVAR; k++) {
                        xLUx[k] *= (1.0 / ycc);
                        xLUy[k] *= (1.0 / ycc);
                    }
                }
            }
    */
        }
    }
    //
    // dt in the jacobian is replaced with 1.0, which will turn this from a 
    // dimensional delta_V to a dv/dt
    // Rest is code for accumilating the residual for solution monitoring
    // .... sounds like TS may already do this so not implementing for now
    /*
    if (accur ==1){
        for (int i=0; i<NVAR; i++) {
            ressum[i] += ressumx[i] + ressumy[i];
        }
    }
    if (iter==0) {
        for (int i=0; i<NVAR; i++){
            res0[i] = ressum[i];
        }
    }
    restotal = 0.0;
    double ressumsum{0.0}, res0sum{0.0};
    for (int i=0; i<NVAR; i++){
        ASSERT(ressum[i] >= 0.0, "Nonpositive Residual")
        if (res0[i] < 1e-16) res0[i] = fmax(ressum[i], 1e-16);
        ressumsum += ressum[i];
        res0sum += res0[i];
    }

    double rss_gather[2] = {ressumsum, res0sum};
    PetscBarrier(PETSC_NULLPTR);
    for (int iblk=1; iblk<world_size; iblk++) {
        double buffer[2] = {rss_gather[0], rss_gather[1]};
        MPI_Status status;

        if (bnum==0) {
            MPI_Recv(&buffer, 2, MPI_DOUBLE, iblk, 0, PETSC_COMM_WORLD, &status);
            rss_gather[0] += buffer[0];
            rss_gather[1] += buffer[1];
        } else if (iblk==bnum) {
            MPI_Send(&buffer, 2, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD);
        }
    }
    */
    //
    //
    //
    // 
    PetscCall(VecRestoreArray(U, &unk));
    PetscCall(VecRestoreArray(Fu, &res));
    PetscCall(VecRestoreArray(Fv, &dv ));
    PetscCall(VecDestroy(&Fu));
    PetscFunctionreturn(PETSC_SUCCESS);
} // calc_rhs
//
// ======================================================================================
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~==========
// ======================================================================================
//
int solve_with_petsc(int argc, char **argv,
	             int mxiter. int nelem) {
    // Input param
    PetscInt n_solvec_local = nelem*NVAR;
    PetscInt n_solvec_globl; //Will be calculated later
    // Function Local Stuff
    TS       ts; // Time Stepping context
    TSAdapt  adapt; //Adaptive timestep context
    Vec      U;
    MPI_Comm comm;
    AppCtx   app_ctx;
    PetscScalar *petsc_array;
    
    // Load Information onto Application Context
    app_ctx -> geofa = geofa;
    app_ctx -> geoel = geoel;
    //.... and so on

    // Initial Stuff
    comm = PETSC_COMM_WORLD;
    PetscCall(TSCreate(comm, &ts));
    PetscCall(TSSetProblemType(ts, TS_NONLINEAR));
    PetscCall(TSSetApplicationContext(ts, &app_ctx));
    
    //  Initialize Solution Vector
    //....  Get total vector length
    MPI_Allreduce(&n_solvec_local, &n_solvec_global, 1, MPI_INT, MPI_SUM, comm);
    //....  Create global vector
    PetscCall(VecCreateMPI(comm, n_solvec_local, n_solvec_globl, &U));
    //....  Copy in the current process's portion of vector
    PetscCall(VecGetArray(U, &petsc_array));
    PetscCall(PetscMemcpy(petsc_array, unk, n_solvec_local*sizeof(PetscScalar)));
    PetscCall(VecRestoreArray(U, &petsc_array));
    //....  Hand off to TS object
    PetscCall(TSSetSolution(ts, U));
    
    // Time Stepping Method Setup
    TSSetType(ts, TSEULER); 		// Time Stepping Method
    TSSetTime(ts, 0.0);     		// Initial Time
    TSSetTimeStep(ts,1.0);  		// Initial timestep (overwritten in prestep)
    TSSetMaxSteps(ts, mxiter); 		// Maximum number of time steps

    // Set Up Function evaluation Stuff
    PetscCall(TSSetRHSFunction(ts, NULL, calc_rhs, app_ctx);
    TSSetPreStep(ts, pre_step_routine); // Function called at the beginning of each time step
					// It's used to calculate the timestep based on a CFl condition
    //TSSetPostStage(ts, .....) 		// Function called after each stage (can be used for filtering)
    //
    PetscCall(TSSetFromOptions(ts));
    PetscCall(TSMonitorSet(ts, TSMonitorStdio, NULL, NULL));
    //
    PetscCall(TSSolve(ts, U));
    //
    PetscCall(TSDestroy(ts));
    PetscCall(VecDestroy(U));
    //
    return 0;
}
// ======================================================================================
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~==========
// ======================================================================================
int main(int argc, char **argv) {
    // ========== Input Parameters	
    double p0, u0, tol, CFL, T0, v0, rho0, gam, damp, duscale,mxangle;
    int accur, iaxi, ivisc, niter, printiter, saveiter;
    // ========== Mesh Variables	
    int nx, ny, npoin, nb, nelem, nface;
    int* ibound;
    double *geoel, *geofa, *yfa, *xfa, *x, *y;
    // ========== Solution Variables	
    //Vec     unk_gl, res_gl, dv_gl; // Global
    double *unk, *res, *dv;
    double *ux = nullptr, *uy = nullptr, *resx=nullptr, *resy=nullptr,  
	   *dvx=nullptr,  *dvy=nullptr;

     /*
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int bnum;
    MPI_Comm_rank(MPI_COMM_WORLD, &bnum);
    */
    
    //PetscCall(PetscInitialize(&argc, &argv, NULL, help));
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    PetscMPIInt bnum,world_size;
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &world_size));
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &bnum));

    if(bnum==0) {printf("PAR2D running on %3d processes.\n",world_size);}

    //=========================================================================
    //=========  Read input file  =============================================
    //=========================================================================
    Thermo air = Thermo();
    // I know this is awful, fight me about it
    std::string line;
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "P0 " << line << std::endl;}
        p0 = stod(line);    
    std::getline(std::cin, line);
    if(bnum==0) { std::cout << "u0 " << line << std::endl;}
        u0 = stod(line);    
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "v0 " << line << std::endl;}
        v0 = stod(line);    
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "T0 " << line << std::endl;}
        T0 = stod(line);    
        rho0 = p0 / (air.Rs[0]*T0);
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "gam " << line << std::endl;}
        gam   = stod(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "accur " << line << std::endl;}
        accur = stoi(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "iaxi " << line << std::endl;}
        iaxi  = stoi(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "ivisc " << line << std::endl;}
        ivisc = stoi(line);
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "niter " << line << std::endl;}
        niter     = stoi(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "cfl " << line << std::endl;}
        CFL       = stod(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "tol " << line << std::endl;}
        tol       = stod(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "printiter " << line << std::endl;}
        printiter = stoi(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "saveiter " << line << std::endl;}
        saveiter  = stoi(line);
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "damp " << line << std::endl;}
        damp    = stod(line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "duscale " << line << std::endl;}
        duscale = stod(line);
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    std::getline(std::cin, line);
    if(bnum==0) {std::cout << "mxangle (not imp)" << line << std::endl;}
        mxangle = stod(line);
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    int mxiter = niter; //maximum number of iteration before stopping
    double time = 0.0;

    if (bnum==0) printf("==================== Loading Mesh ====================\n");
    //=========================================================================
    //==========  Load Mesh  ==================================================
    //=========================================================================

    //read in mesh file
    int nblock_file,ncell_glob;
    FILE* fconn = nullptr;
    int bbounds[4];
    int bids[4]; // bottom, right, top left
    for (int iblk=0; iblk<world_size; iblk++) {
        if (iblk == bnum) {
            printf("::%3d:: opening mesh file \n", bnum);
            read_mesh(bnum, &nx, &ny, &ibound, &x, &y);

            printf("::%3d:: opening conn file \n", bnum);
            fconn = fopen("./Grid/grid.conn", "r");
            if (fconn == nullptr) {
                printf("::%3d:: Failed to open conn file.\n", bnum);
                exit(1);
            }

            fscanf(fconn, "%d,%d", &nblock_file, &ncell_glob);
            if (nblock_file != world_size) {
                printf("Discrepancy between #blocks and #procs!");
                exit(1);
            }
            for (int iblk = 0; iblk <= bnum; iblk++) {
                if (iblk == bnum) {
                    fscanf(fconn, "%*d: (%d,%d) (%d,%d) (%d,%d) (%d,%d)",
                           &bbounds[3], &bids[3], &bbounds[1], &bids[1], &bbounds[0], &bids[0],
                           &bbounds[2], &bids[2]);
                } else {
                    fscanf(fconn, "%*d: (%*d,%*d) (%*d,%*d) (%*d,%*d) (%*d,%*d)");
                }
            }
            //bids[0]--;
            //bids[1]--;
            //bids[2]--;
            //bids[3]--;
            fclose(fconn);
        }
	PetscBarrier(PETSC_NULLPTR);
    }

    // Calculate nums based off of grid file
    npoin = nx*ny;
    nb = 2*nx + 2*ny;
    nelem = (nx-1)*(ny-1);
    nface = 2*nelem + nx + ny;
    //Find elem volume and centroid, face len and normals
    for (int iblk=0; iblk<world_size; iblk++) {
        if (iblk==bnum) {
            printf("::%3d::Calculating Grid Metrics..... \n", bnum);
            calc_geoel_geofa(nx, ny, x, y, &geoel, &geofa, &yfa, &xfa);
        }
        PetscBarrier(PETSC_NULLPTR);
    }
    //
    //
    PetscBarrier(PETSC_NULLPTR);
    if (bnum==0) printf("==================== Initializing ====================\n");
    //=========================================================================
    //==========  Setup Solution Variables  ===================================
    //=========================================================================
    //
    // Allocate Vectors with PETSc
    int n_solvec_local = NVAR*nelem;
    //int n_solvec_globl = NVAR*ncell_glob;
    // Global
    //PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n_solvec_local, n_solvec_globl, &unk_gl));
    //PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n_solvec_local, n_solvec_globl, &res_gl));
    //PetscCall(VecCreateMPI(PETSC_COMM_WORLD, n_solvec_local, n_solvec_globl, &dv_gl));
    // Local
    auto* unk  = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* res  = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* dv   = (double*)malloc(NVAR*nelem*sizeof(double));
    //PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_solvec_local, &unk);
    //PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_solvec_local, &res);
    //PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_solvec_local, &dv);
    //
    if (accur == 1) {
        printf("PETSc suport of higher order hasn't been implimented yet.\n"); exit(1);
        ux = (double *) malloc(NVAR * nelem * sizeof(double));    //xi  derivative
        uy = (double *) malloc(NVAR * nelem * sizeof(double));    //eta derivative
        resx = (double*)malloc(NVAR*nelem*sizeof(double));
        resy = (double*)malloc(NVAR*nelem*sizeof(double));
        dvx  = (double*)malloc(NVAR*nelem*sizeof(double));
        dvy  = (double*)malloc(NVAR*nelem*sizeof(double));
    }

    //initialize solution on mesh (zero aoa)
    double uFS[4], uBP[4];
    uFS[0] = rho0;
    uFS[1] = u0; //0.0
    uFS[2] = v0;
    uFS[3] = T0;
    
    uBP[0] = rho0;
    uBP[1] = u0;
    uBP[2] = v0;
    uBP[3] = T0;

    //Plenum State
    ///back pressure here

    for (int ielem=0; ielem<nelem; ielem++){
        unk[NVAR*ielem]   = uFS[0]; // *0.2
        unk[NVAR*ielem+1] = uFS[1]; // 600.0
        unk[NVAR*ielem+2] = uFS[2];
        unk[NVAR*ielem+3] = uFS[3];

        if (accur==1) {
            ux[NVAR * ielem] = 0.0;
            ux[NVAR * ielem + 1] = 0.0;
            ux[NVAR * ielem + 2] = 0.0;
            ux[NVAR * ielem + 3] = 0.0;

            uy[NVAR * ielem] = 0.0;
            uy[NVAR * ielem + 1] = 0.0;
            uy[NVAR * ielem + 2] = 0.0;
            uy[NVAR * ielem + 3] = 0.0;
        }
    }

    PetscBarrier(PETSC_NULLPTR);
    if (bnum==0) printf("===== Generating Mesh and Initial State Tecplot Files ====\n");
    print_elem_stats("MeshVolumeStats", nx, ny, geoel);
    if (accur==1){
        print_state_DGP1(time, 0, bnum,"Initial State", nx, ny, air, x, y, unk, ux, uy, geoel);
    } else {
        print_state(0, bnum,"Initial State", nx, ny, air, x, y, unk, geoel);
    }

    //Find timestep based off of CFL limit for initial condition (dt = CFL dx / c )
    double dt;
    State* ElemVar;//
    ElemVar = (State*)malloc(nelem*sizeof(State));// [nelem];
    //Set up structures for calculating/containing non-state variables on each element
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++) {
            int ie = IJ(i,j,nx-1);
            double* unkel = &(unk[IJK(i,j,0,nx-1,NVAR)]);

            //if (i > 95){
            //    unkel[1] = 1000.0;
            //}

            ElemVar[ie].Initialize(unkel);
            ElemVar[ie].UpdateState(air);
        }
    }
    //Same memory to be used for each local matrix (chg this if making parallel [jk lol])
    auto D = (double**)malloc((NVAR) * sizeof(double*));
    for (int k = 0; k < NVAR; k++)
        D[k] = (double*)malloc( (NVAR) * sizeof(double));

    //save residual history
    FILE *fres;
    if (bnum==0) {
        fres = fopen("./Outputs/res.tec", "w");
        if (fres == nullptr) {
            printf("~~~~~~~~~ Failed to save residual file, error:%d\n", errno);
        } else {
            fprintf(fres, "Residual history\n");
//        fprintf(fres, "%d,\t%le,\n",0,1.0);
        }
    }

    PetscBarrier(PETSC_NULLPTR);
    if (bnum==0) printf("==================== Starting Solver ====================\n");


    double res0[NVAR]{};
    double ressum[NVAR], ressumx[NVAR], ressumy[NVAR], restotal;
    int iter;
    for (iter=0; iter<mxiter; iter++){
        //Explicit Euler Time Integration

        //Find global timestep based off of CFl condition
        dt = find_dt(air, nx, ny, CFL, unk, ElemVar[0], geofa);
        if(DEBUG) {printf("::%3d::Calculated Timestep..... \n", bnum);}
        // sync timestep across threads
        PetscBarrier(PETSC_NULLPTR);
        for (int iblk = 1; iblk < world_size; iblk++) {
            double buffer = dt;
            MPI_Status status;
            if (bnum == 0) {
                MPI_Recv(&buffer, 1, MPI_DOUBLE, iblk, 0, PETSC_COMM_WORLD, &status);
                dt = fmax(dt, buffer);
            } else if (iblk == bnum) {
                MPI_Send(&buffer, 1, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD);
            }
        }
        if(DEBUG) {printf("::%3d::Gathered Timestep..... \n", bnum);}
	fflush(stdout);
        PetscBarrier(PETSC_NULLPTR);
        for (int iblk = 1; iblk < world_size; iblk++) {
            double buffer = dt;
            MPI_Status status;
            if (bnum == 0) {
                MPI_Send(&buffer, 1, MPI_DOUBLE, iblk, 0, PETSC_COMM_WORLD);
            } else if (iblk == bnum) {
                MPI_Recv(&buffer, 1, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD, &status);
                dt = buffer;
            }
        }
        if(DEBUG) {printf("::%3d::Scattered Timestep..... \n", bnum);}
        time += dt;
        PetscBarrier(PETSC_NULLPTR);

        //calculate the right hand side residual term (change of conserved quantities)
        calc_dudt(ivisc, accur, iaxi, mxangle, bbounds, bids, nx, ny, air, ElemVar, uFS,
                  ibound, geoel, geofa, yfa, xfa, unk, ux, uy, res, resx, resy);
        calculate_residual(nx, ny, res, ressum);
        if(DEBUG) {printf("::%3d::Calculated dudt..... \n", bnum);}

        if ( accur == 1 ) {
            calculate_residual(nx, ny, resx, ressumx);
            calculate_residual(nx, ny, resy, ressumy);
        }

        //========== Solve linear system on each element (turns chg in conservatives to change in solution variables)
        int flg = 0;
        for (int i=0; i<nx-1; i++) {
            for (int j=0; j<ny-1; j++) {
                double *unkij = &(unk[IJK(i, j, 0, nx-1, NVAR)]);
                double LUtol = 1e-16;
                int iel = IJ(i,j,nx-1);
                int N = NVAR;
                int P[NVAR+1]{}; //permutation vector for pivoting

                //get the rhs block needed
                double *b = &(res[IJK(i, j, 0, nx - 1, NVAR)]);
                double *xLU = &(dv[ IJK(i, j, 0, nx - 1, NVAR)]);

                //Evaluate the jacobian / Implicit matrix
                BuildJacobian(dt, unkij, ElemVar[iel], D);
                LUPDecompose(D, N, LUtol, P);
                LUPSolve(D, P, b, N, xLU);

                //axi modification
                double ycc = 1.0;
                if (iaxi==1) {
                    ycc = geoel[IJK(i, j, 2, nx - 1, 3)];
                    for (int k = 0; k < NVAR; k++) {
                        xLU[k] *= (1.0 / ycc);
                    }
                }

                if (accur == 1) {
                    ////  CURRENTLY USING THE CELL CENTERED JACOBIAN VALUE FOR THE 1ST ORDER MODES
                    double *bx = &(resx[IJK(i, j, 0, nx - 1, NVAR)]);
                    double *by = &(resy[IJK(i, j, 0, nx - 1, NVAR)]);
                    double *xLUx = &(dvx[IJK(i, j, 0, nx - 1, NVAR)]);
                    double *xLUy = &(dvy[IJK(i, j, 0, nx - 1, NVAR)]);
                    LUPSolve(D, P, bx, N, xLUx);
                    LUPSolve(D, P, by, N, xLUy);

                    //axi modification
                    if (iaxi==1) {
                        for (int k = 0; k < NVAR; k++) {
                            xLUx[k] *= (1.0 / ycc);
                            xLUy[k] *= (1.0 / ycc);
                        }
                    }
                }
            }
        }

        //perform iteration

        //Mechanism for switching between 1st and 2nd order
        //if (iter<= 1.5*(3000.0*(0.3/CFL)*(101.0/101.0))) damp = 0.0;
        //if (iter >= 10000) damp = 0.95;


        if (accur ==1){
            for (int i=0; i<NVAR; i++) {
                ressum[i] += ressumx[i] + ressumy[i];
            }
        }
        if (iter==0) {
            for (int i=0; i<NVAR; i++){
                res0[i] = ressum[i];
            }
        }
        restotal = 0.0;
        double ressumsum{0.0}, res0sum{0.0};
        for (int i=0; i<NVAR; i++){
            ASSERT(ressum[i] >= 0.0, "Nonpositive Residual")
            if (res0[i] < 1e-16) res0[i] = fmax(ressum[i], 1e-16);
            ressumsum += ressum[i];
            res0sum += res0[i];
        }

        double rss_gather[2] = {ressumsum, res0sum};
        PetscBarrier(PETSC_NULLPTR);
        for (int iblk=1; iblk<world_size; iblk++) {
            double buffer[2] = {rss_gather[0], rss_gather[1]};
            MPI_Status status;

            if (bnum==0) {
                MPI_Recv(&buffer, 2, MPI_DOUBLE, iblk, 0, PETSC_COMM_WORLD, &status);
                rss_gather[0] += buffer[0];
                rss_gather[1] += buffer[1];
            } else if (iblk==bnum) {
                MPI_Send(&buffer, 2, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD);
            }
        }

        for (int ielem=0; ielem<nelem; ielem++){
            int iu = NVAR*ielem;
            unk[iu  ] += dv[iu  ];
            unk[iu+1] += dv[iu + 1];
            unk[iu+2] += dv[iu + 2];
            unk[iu+3] += dv[iu + 3];

            ASSERT(unk[iu] > 0.0, "Nonpositive density")

            //unk[iu+3] = fmax(unk[iu+3], 201.0); //limit temperature
            ElemVar[ielem].UpdateState(air);

            if (accur ==1) {
                ux[iu  ] += dvx[iu  ];
                ux[iu+1] += dvx[iu + 1];
                ux[iu+2] += dvx[iu + 2];
                ux[iu+3] += dvx[iu + 3];

                uy[iu  ] += dvy[iu  ];
                uy[iu+1] += dvy[iu + 1];
                uy[iu+2] += dvy[iu + 2];
                uy[iu+3] += dvy[iu + 3];

                ux[iu  ] *= damp; 
                ux[iu+1] *= damp;
                ux[iu+2] *= damp;
                ux[iu+3] *= damp;

                uy[iu  ] *= damp;
                uy[iu+1] *= damp;
                uy[iu+2] *= damp;
                uy[iu+3] *= damp;

                //Slope limiting
                int iuim, iuip, iujm, iujp;
                iuim = iu - IJK(1,0,0,nx-1,NVAR);
                iuip = iu + IJK(1,0,0,nx-1,NVAR);
                iujm = iu - IJK(0,1,0,nx-1,NVAR);
                iujp = iu + IJK(0,1,0,nx-1,NVAR);
                for (int kvar=0;kvar<NVAR;kvar++){
                    double du;
                    int nu = nelem*NVAR;
                    if (iuip < nu-1  and iuim >= 0.0) {
                       du = duscale*((unk[iuip + kvar] - unk[iuim + kvar]) / 2.0);
                       ux[iu+kvar] = sign(ux[iu+kvar])*fmin(fabs(ux[iu+kvar]), fabs(du));
                    } else if (iuim < 0.0) {
                        //bottom boundary
                        du = duscale*((unk[iuip + kvar] - unk[iu + kvar]) / 1.0);
                        ux[iu+kvar] = sign(ux[iu+kvar])*fmin(fabs(ux[iu+kvar]), fabs(du));
                    } else {
                        //top boundary
                        du = duscale*((unk[iu + kvar] - unk[iuim + kvar]) / 1.0);
                        ux[iu+kvar] = sign(ux[iu+kvar])*fmin(fabs(ux[iu+kvar]), fabs(du));
                    }

                    if (iujp <= nu-1 and iujm >= 0) {
                        du = duscale*((unk[iujp + kvar] - unk[iujm + kvar]) / 2.0);
                        uy[iu+kvar] = sign(uy[iu+kvar])*fmin(fabs(uy[iu+kvar]), fabs(du));
                    } else if (iujm < 0.0) {
                        //bottom boundary
                        du = duscale*((unk[iujp + kvar] - unk[iu + kvar]) / 1.0);
                        uy[iu+kvar] = sign(uy[iu+kvar])*fmin(fabs(uy[iu+kvar]), fabs(du));
                    } else {
                        //top boundary
                        du = duscale*((unk[iu + kvar] - unk[iujm + kvar]) / 1.0);
                        uy[iu+kvar] = sign(uy[iu+kvar])*fmin(fabs(uy[iu+kvar]), fabs(du));
                    }
                }

            }
        }
        /*
        if (accur ==1){
            for (int i=0; i<NVAR; i++) {
                ressum[i] += ressumx[i] + ressumy[i];
            }
        }

        if (iter==0) {
            for (int i=0; i<NVAR; i++){
                res0[i] = ressum[i];
            }
        }
        restotal = 0.0;
        for (int i=0; i<NVAR; i++){
            ASSERT(ressum[i] >= 0.0 && !__isnan(ressum[i]), "Invalid Residual")
            if (res0[i] < 1e-16) res0[i] = fmax(ressum[i], 1e-16);
            restotal += ressum[i] / res0[i];
        }
        */
        if (iter > 0 and iter % saveiter == 0) {
            //printf("Saving current Solution\n");
            if (accur == 1) {
                print_state_DGP1(time, iter, bnum, "Final State", nx, ny, air, x, y, unk, ux, uy, geoel);
            } else {
                print_state(iter, bnum, "Final State", nx, ny, air, x, y, unk, geoel);
            }
        }
        PetscBarrier(PETSC_NULLPTR);
        if (bnum==0) {
            double relres_gather = rss_gather[0] / rss_gather[1];
            fprintf(fres, "%d,\t%le\n", iter, relres_gather);


            if (iter % printiter == 0) {
                printf("Iter:%7d\tdt:%7.4e \t\t RelativeTotalResisual:  %8.5e ", \
                        iter, dt, relres_gather);
                if (accur == 1) {
                    printf("\t damp:%f, duscale:%f", damp, duscale);
                }
                printf("\n");
            }
            if (world_size > 1) {
                if (relres_gather < tol) {
                    int buf = 1;
                    MPI_Send(&buf, 1, MPI_INT, 1, 0, PETSC_COMM_WORLD);
                    printf("Thread %3d breaking main loop.\n", bnum);
                    break;  // and damp >= 1) break;
                } else {
                    int buf = 0;
                    MPI_Send(&buf, 1, MPI_INT, 1, 0, PETSC_COMM_WORLD);
                }
            } else { // running single core
                if (relres_gather < tol) {
                    break; 
                }
            }
        } else { // not the head node
            int buf{0};
            if (bnum < world_size-1) {
                MPI_Status status;
                MPI_Recv(&buf, 1, MPI_INT, bnum - 1, 0, PETSC_COMM_WORLD, &status);
                MPI_Send(&buf, 1, MPI_INT, bnum+1, 0, PETSC_COMM_WORLD);
            } else {
                MPI_Status status;
                MPI_Recv(&buf, 1, MPI_INT, bnum - 1, 0, PETSC_COMM_WORLD, &status);
            }
            if (buf ==1) {
                printf("Thread %3d breaking main loop.\n", bnum);
                break;
            }
        }
    }

    PetscBarrier(PETSC_NULLPTR);
    sleep(1);
    if (bnum==0) {
        printf("==================== Calculation Finished ====================\n");
    }
    printf("::%3d:: Saving Solution File..... \n",bnum);

    if (accur==1){
        print_state_DGP1(time,iter,bnum,"Final State", nx, ny, air, x, y, unk, ux, uy, geoel);
    } else {
        print_state(iter,bnum,"Final State", nx, ny, air, x, y, unk, geoel);
    }
    //print_state_axi("Final State", nx, ny, air, x, y, unk, geoel);

    free(ElemVar);
    free(unk);
    free(ux);
    free(uy);
    free(res);
    free(resx);
    free(resy);
    free(dv);
    printf("%3d Complete.\n",bnum);
    //MPI_Finalize();
    PetscCall(PetscFinalize());
}
