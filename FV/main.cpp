/* 2D Navier Stokes solver for a structured quadrilateral grid
 *
 * T. Sailor Koeplinger - Feb2024
 */
#include <iostream>
#include <cmath>
#include <mpi.h>
#include <unistd.h>
#include "FileIO.h"
#include "Indexing.h"
#include "SpatialDiscretization.h"
#include "MeshModule.h"
#include "StateVariables.h"
#include "LUtools.h"
#include "Jacobian.h"
#include "Thermo.h"


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

void vec_copy(double n, double* a, double* b){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

int main() {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get the rank of the process
    int bnum;
    MPI_Comm_rank(MPI_COMM_WORLD, &bnum);

    //Read in setup file
    double p0, u0, tol, CFL, T0, v0, rho0;
    tol = 1e-6;
    CFL = 0.1;
    Thermo air = Thermo();

    p0 = 1000.0;
    u0 = 1736.0;
    T0 = 300.0;
    rho0 = p0 / (air.Rs[0]*T0);
    v0 = 0.0;
    int mxiter = NITER; //maximum number of iteration before stopping
    int printiter = 1;
    int saveiter = 10;
    double damp = 0.95;///fmin(iter / (3000.0*0.3/CFL),1.0);  //coarse mesh 3000, fine 7500
    double duscale = 0.7;

    if (bnum==0) printf("==================== Loading Mesh ====================\n");
    //==================== Load Mesh ====================
    int nx, ny, npoin, nb, nelem, nface;
    double *x, *y;
    int* ibound;
    double* geoel;
    double* geofa;
    double* yfa;
    double* xfa;

    //read in mesh file
    int nblock_file;
    FILE* fconn = nullptr;
    int bbounds[4];
    int bids[4]; // bottom, right, top left
    //for (int iblk=0; iblk<world_size; iblk++) {
    //    if (iblk == bnum) {
            printf("Reading Mesh File..... \n");
            read_mesh(bnum, &nx, &ny, &ibound, &x, &y);

            printf("Reading Block Connectivity File..... \n");
            fconn = fopen("../Case/grid.conn", "r");
            if (fconn == nullptr) {
                printf("::%3d:: Failed to open conn file.\n", bnum);
                exit(1);
            }
            printf("::%3d:: opened conn file \n", bnum);

            fscanf(fconn, "%d", &nblock_file);
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
            bids[0]--;
            bids[1]--;
            bids[2]--;
            bids[3]--;
            fclose(fconn);
    //    }
    //}

    // Calculate nums based off of grid file
    npoin = nx*ny;
    nb = 2*nx + 2*ny;
    nelem = (nx-1)*(ny-1);
    nface = 2*nelem + nx + ny;
    //Find elem volume and centroid, face len and normals
    //for (int iblk=0; iblk<world_size; iblk++) {
    //    if (iblk==bnum) {
            printf("::%3d::Calculating Grid Metrics..... \n", bnum);
            calc_geoel_geofa(nx, ny, x, y, &geoel, &geofa, &yfa, &xfa);
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD);
    //}

    MPI_Barrier(MPI_COMM_WORLD);
    if (bnum==0) printf("==================== Initializing ====================\n");
    //==================== Setup for Sim ====================
    auto* unk  = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* res  = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* dv   = (double*)malloc(NVAR*nelem*sizeof(double));

    double *ux = nullptr, *uy = nullptr;

    if (ACCUR == 1) {
        ux = (double *) malloc(NVAR * nelem * sizeof(double));    //xi  derivative
        uy = (double *) malloc(NVAR * nelem * sizeof(double));    //eta derivative
    }
    auto* resx = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* resy = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* dvx  = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* dvy  = (double*)malloc(NVAR*nelem*sizeof(double));

    //initialize solution on mesh (zero aoa)
    double uFS[4], uBP[4];
    uFS[0] = rho0;
    uFS[1] = u0;
    uFS[2] = v0;
    uFS[3] = T0;

    //Plenum State
    ///back pressure here

    for (int ielem=0; ielem<nelem; ielem++){
        unk[NVAR*ielem]   = uFS[0];
        unk[NVAR*ielem+1] = uFS[1];
        unk[NVAR*ielem+2] = uFS[2];
        unk[NVAR*ielem+3] = uFS[3];

        if (ACCUR==1) {
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

    MPI_Barrier(MPI_COMM_WORLD);
    if (bnum==0) printf("===== Generating Mesh and Initial State Tecplot Files ====\n");
    print_elem_stats("MeshVolumeStats", nx, ny, geoel);
    if (ACCUR==1){
        print_state_DGP1(bnum,"Initial State", nx, ny, air, x, y, unk, ux, uy, geoel);
    } else {
        print_state(bnum,"Initial State", nx, ny, air, x, y, unk, geoel);
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
        fres = fopen("../Outputs/res.tec", "w");
        if (fres == nullptr) {
            printf("~~~~~~~~~ Failed to save residual file, error:%d\n", errno);
        } else {
            fprintf(fres, "Residual history\n");
//        fprintf(fres, "%d,\t%le,\n",0,1.0);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (bnum==0) printf("==================== Starting Solver ====================\n");


    double res0[NVAR]{};
    double ressum[NVAR], ressumx[NVAR], ressumy[NVAR], restotal;
    int iter;
    for (iter=0; iter<mxiter; iter++){
        //Explicit Euler Time Integration

        //Find global timestep based off of CFl condition
        dt = find_dt(air, nx, ny, CFL, unk, ElemVar[0], geofa);
        // sync timestep across threads
        MPI_Barrier(MPI_COMM_WORLD);
            for (int iblk = 1; iblk < world_size; iblk++) {
                double buffer = dt;
                MPI_Status status;
                if (bnum == 0) {
                    MPI_Recv(&buffer, 1, MPI_DOUBLE, iblk, 0, MPI_COMM_WORLD, &status);
                    dt = fmax(dt, buffer);
                } else if (iblk == bnum) {
                    MPI_Send(&buffer, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            for (int iblk = 1; iblk < world_size; iblk++) {
                double buffer = dt;
                MPI_Status status;
                if (bnum == 0) {
                    MPI_Send(&buffer, 1, MPI_DOUBLE, iblk, 0, MPI_COMM_WORLD);
                } else if (iblk == bnum) {
                    MPI_Recv(&buffer, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
                    dt = buffer;
                }
            }

        //calculate the right hand side residual term (change of conserved quantities)
        calc_dudt(bbounds, bids, nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, yfa, xfa, unk, ux, uy, res, resx, resy);
        calculate_residual(nx, ny, res, ressum);

        if ( ACCUR == 1 ) {
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
                if (IAXI==1) {
                    ycc = geoel[IJK(i, j, 2, nx - 1, 3)];
                    for (int k = 0; k < NVAR; k++) {
                        xLU[k] *= (1.0 / ycc);
                    }
                }

                if (ACCUR == 1) {
                    ////  CURRENTLY USING THE CELL CENTERED JACOBIAN VALUE FOR THE 1ST ORDER MODES
                    double *bx = &(resx[IJK(i, j, 0, nx - 1, NVAR)]);
                    double *by = &(resy[IJK(i, j, 0, nx - 1, NVAR)]);
                    double *xLUx = &(dvx[IJK(i, j, 0, nx - 1, NVAR)]);
                    double *xLUy = &(dvy[IJK(i, j, 0, nx - 1, NVAR)]);
                    LUPSolve(D, P, bx, N, xLUx);
                    LUPSolve(D, P, by, N, xLUy);

                    //axi modification
                    if (IAXI==1) {
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
        //if (iter<= 400) damp = 0.0;


        if (ACCUR ==1){
            for (int i=0; i<NVAR; i++) {
                ressum[i] += damp*ressumx[i] + damp*ressumy[i];
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
        MPI_Barrier(MPI_COMM_WORLD);
        for (int iblk=1; iblk<world_size; iblk++) {
            double buffer[2] = {rss_gather[0], rss_gather[1]};
            MPI_Status status;

            if (bnum==0) {
                MPI_Recv(&buffer, 2, MPI_DOUBLE, iblk, 0, MPI_COMM_WORLD, &status);
                rss_gather[0] += buffer[0];
                rss_gather[1] += buffer[1];
            } else if (iblk==bnum) {
                MPI_Send(&buffer, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
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

            if (ACCUR ==1) {
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
        if (ACCUR ==1){
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
            if (ACCUR == 1) {
                print_state_DGP1(bnum, "Final State", nx, ny, air, x, y, unk, ux, uy, geoel);
            } else {
                print_state(bnum, "Final State", nx, ny, air, x, y, unk, geoel);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (bnum==0) {
            double relres_gather = rss_gather[0] / rss_gather[1];
            fprintf(fres, "%d,\t%le\n", iter, relres_gather);


            if (iter % printiter == 0) {
                printf("Iter:%7d\tdt:%7.4e \t\t RelativeTotalResisual:  %8.5e ", \
                        iter, dt, relres_gather);
                if (ACCUR == 1) {
                    printf("\t damp:%f, duscale:%f", damp, duscale);
                }
                printf("\n");
            }
            if (world_size > 1) {
                if (relres_gather < tol) {
                    int buf = 1;
                    MPI_Send(&buf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                    printf("Thread %3d breaking main loop.\n", bnum);
                    break;  // and damp >= 1) break;
                } else {
                    int buf = 0;
                    MPI_Send(&buf, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                }
            }
        } else {
            int buf{0};
            if (bnum < world_size-1) {
                MPI_Status status;
                MPI_Recv(&buf, 1, MPI_INT, bnum - 1, 0, MPI_COMM_WORLD, &status);
                MPI_Send(&buf, 1, MPI_INT, bnum+1, 0, MPI_COMM_WORLD);
            } else {
                MPI_Status status;
                MPI_Recv(&buf, 1, MPI_INT, bnum - 1, 0, MPI_COMM_WORLD, &status);
            }
            if (buf ==1) {
                printf("Thread %3d breaking main loop.\n", bnum);
                break;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    sleep(1);
    if (bnum==0) {
        printf("==================== Calculation Finished ====================\n");
    }
    printf("::%3d:: Saving Solution File..... \n",bnum);

    if (ACCUR==1){
        print_state_DGP1(bnum,"Final State", nx, ny, air, x, y, unk, ux, uy, geoel);
    } else {
        print_state(bnum,"Final State", nx, ny, air, x, y, unk, geoel);
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
    MPI_Finalize();
}
