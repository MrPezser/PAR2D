//
// Created by Tsail on 1/28/2024.
//

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

#include "SpatialDiscretization.h"
#include "Indexing.h"
#include "BoundaryConditions.h"
#include "EulerFlux.h"
#include "StateVariables.h"
#include "DGP1Tools.h"

void generate_ghost_cells(int nx, int ny, double* unk, double* ux, double* uy, State* ElemVar, Thermo air, int* ibound,
                            double* geofa, double* uFS, double* uGBot, double* uGTop, double* uGLeft, double* uGRight,
                            State* BotVar, State* TopVar, State* LeftVar, State* RightVar){
    double unkelij[NVAR];
    State varij = State();

    //bottom side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[i];
        int iint = IJK(i,0,0, nx-1, 4);
        int iel = IJ(i, 0, nx-1);

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(i, 0, 1,nx,6)];
        normy = geofa[IJK(i, 0, 2,nx,6)];
        //==========Ghost State
        BotVar[i].Initialize(&(uGBot[IJ(0,i,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS,&unk[iint], ElemVar[iel],
                       &(uGBot[IJ(0,i,NVAR)]));
        BotVar[i].UpdateState(air);
    }
    //right side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[j+ nx-1];
        int iint = IJK(nx-2,j,0, nx-1, NVAR);
        int iel = IJ(nx-2, j, nx-1);

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(0, j, 4,nx,6)];
        normy = geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        RightVar[j].Initialize(&(uGRight[IJ(0,j,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS,&unk[iint], ElemVar[iel],
                       &(uGRight[IJ(0,j,NVAR)]));
        RightVar[j].UpdateState(air);
    }

    //top side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        int ib = 2*(nx-1) + ny - 1 - i - 1;
        btype = ibound[ib];
        int iint = IJK(i,ny-2,0, nx-1, NVAR);
        int iel = IJ(i, ny-2, nx-1);

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(i, ny-1, 1,nx,6)];
        normy = -geofa[IJK(i, ny-1, 2,nx,6)];
        //==========Ghost State
        TopVar[i].Initialize(&(uGTop[IJ(0,i,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS, &unk[iint], ElemVar[iel],
                       &(uGTop[IJ(0,i,NVAR)]));
        TopVar[i].UpdateState(air);
    }

    //left side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        int jb = 2*(ny-1)+2*(nx-1)-j-1;
        btype = ibound[jb];
        int iint = IJK(0,j,0, nx-1, 4);
        int iel = IJ(0, j, nx-1);

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(0, j, 4,nx,6)];
        normy = -geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        LeftVar[j].Initialize(&(uGLeft[IJ(0,j,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS, &unk[iint], ElemVar[iel],
                       &(uGLeft[IJ(0,j,NVAR)]));
        LeftVar[j].UpdateState(air);
    }
}

void viscous(int nx, double normy, double normx, double* uLeft, State& varL, double* uRight, State varR, double* dc, double* visc_contrib){

    if (IVISC==0) {return;}

    // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
    ///Need to make axisymmatric modification and add extra termsnn/JE
    //printf("Need to update viscous fluxes for axisymmetric and DP higher order\n");
    //exit(0);
    double tau[4], Sij[4], st[4], dc2i, trace, vflux[2];
    double uL, vL, uR, vR;

    double mu = 0.5 * (varL.mu + varR.mu);

    dc2i = 1.0 / (dc[0] * dc[0] + dc[1] * dc[1]);

    uL = uLeft[1];
    vL = uLeft[2];
    uR = uRight[1];
    vR = uRight[2];
    //stress tensor
    st[0] = (uR - uL) * dc[0] * dc2i;//  du/dx
    st[1] = (uR - uL) * dc[1] * dc2i;//  du/dy
    st[2] = (vR - vL) * dc[0] * dc2i;//  dv/dx
    st[3] = (vR - vL) * dc[1] * dc2i;//  dv/dy
    trace = (1.0 / 3.0) * (st[0] + st[3]);

    // Viscous strain rate
    Sij[0] = st[0] - trace;                //S_1,1
    Sij[1] = 0.5 * (st[1] + st[2]);         //S_1,2
    Sij[2] = Sij[1];                        //S_2,1
    Sij[3] = st[1] - trace;                //S_2,2

    //Viscous Stress  (Stoke's Law for Monoatoms)
    tau[0] = 2 * mu * Sij[0];
    tau[1] = 2 * mu * Sij[1];
    tau[2] = 2 * mu * Sij[2];
    tau[3] = 2 * mu * Sij[3];

    // Viscous Contributions to Flux
    vflux[0] = tau[0] * normx + tau[2] * normy;
    vflux[1] = tau[1] * normx + tau[3] * normy;

    visc_contrib[0] = vflux[0];
    visc_contrib[1] = vflux[1];
    visc_contrib[2] = ((uL * tau[0] + vL * tau[2]) * normx + (uL * tau[1] + vL * tau[3]) * normy);

    visc_contrib[3] = vflux[0];
    visc_contrib[4] = vflux[1];
    visc_contrib[5] = ((uR * tau[0] + vR * tau[2]) * normx + (uR * tau[1] + vR * tau[3]) * normy);

    for (int iv = 0; iv < 6; iv++) {
        if (__isnan(visc_contrib[iv])) {
            printf("asdfaed\n");
        }
        ASSERT(!__isnan(visc_contrib[iv]), "Visc Flux NaN");
    }
}

void calc_dudt(int* bbounds, int* bids, int nx, int ny, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
               double* geofa, double* yfa, double* xfa, double* unk, double* ux, double* uy, double* dudt, double* duxdt, double* duydt) {
    int nelem = (nx-1)*(ny-1);
    double *rhsel, *rhselx = nullptr, *rhsely = nullptr, parr;
    rhsel  = (double*)malloc(NVAR*nelem*sizeof(double));

    if (ACCUR==1) {
        rhselx = (double *) malloc(NVAR * nelem * sizeof(double));
        rhsely = (double *) malloc(NVAR * nelem * sizeof(double));
    }

    for(int i=0; i<NVAR*nelem; i++) {
        rhsel[i] = 0.0;
        if (ACCUR==1) {
            rhselx[i] = 0.0;
            rhsely[i] = 0.0;
        }
    }

    //Calculate boundary cell state (ghost state)
    double *uGBot, *uGRight, *uGTop, *uGLeft;
    State *BotVar, *TopVar, *RightVar, *LeftVar;
    if (ACCUR == 0) {
        uGBot   = (double*)malloc((NVAR * (nx - 1))*sizeof(double));
        uGRight = (double*)malloc((NVAR * (ny - 1))*sizeof(double));
        uGTop   = (double*)malloc((NVAR * (nx - 1))*sizeof(double));
        uGLeft  = (double*)malloc((NVAR * (ny - 1))*sizeof(double));
        BotVar   = (State*)malloc(((nx - 1))*sizeof(State));
        RightVar = (State*)malloc(((ny - 1))*sizeof(State));
        TopVar   = (State*)malloc(((nx - 1))*sizeof(State));
        LeftVar  = (State*)malloc(((ny - 1))*sizeof(State));
        generate_ghost_cells(nx, ny, unk, ux, uy, ElemVar, air, ibound, geofa, uFS, uGBot, uGTop, uGLeft, uGRight,
                             BotVar, TopVar, LeftVar, RightVar);
    } else if (ACCUR==1){
        uGBot   = (double*)malloc(2*(NVAR * (nx - 1))*sizeof(double));
        uGRight = (double*)malloc(2*(NVAR * (ny - 1))*sizeof(double));
        uGTop   = (double*)malloc(2*(NVAR * (nx - 1))*sizeof(double));
        uGLeft  = (double*)malloc(2*(NVAR * (ny - 1))*sizeof(double));
        BotVar   = (State*)malloc(2*((nx - 1))*sizeof(State));
        RightVar = (State*)malloc(2*((ny - 1))*sizeof(State));
        TopVar   = (State*)malloc(2*((nx - 1))*sizeof(State));
        LeftVar  = (State*)malloc(2*((ny - 1))*sizeof(State));
        DGP1_ghost_cell_generator(nx, ny, unk, ux, uy, ElemVar, air, ibound, geofa, uFS, uGBot, uGTop, uGLeft, uGRight,
                BotVar, TopVar, LeftVar, RightVar);
    } else {
        printf("ACCUR must be 1 or 0.");
        exit(1);
    }

    //Get interior boundary information from other processes
    MPI_Barrier(MPI_COMM_WORLD);
    int bnum;
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &bnum);
    int blockflag[4] = {0,0,0,0};


    for (int iblk=0; iblk<world_size; iblk++) {
        if (world_size==1) continue;
        // Loop through all the blocks
        int iblk2; //tranfer target
        double* usend;
        double* urecv;
        MPI_Status status;

        // If the block is a neighbor initiate send/recieve with said neighbor
        // If the block is the current block, initiate send/recieve
        if ((iblk == bids[3] || (iblk==bnum && bids[3]>bnum)) && blockflag[0]==0 ) {
            // ~~~~~~~~~~~~~~~~~ LEFT BOUNDARY
            iblk2 = bids[3];
            int ntrans = NVAR * (ny - 1);
            ntrans *= (ACCUR + 1);
            // printf("left  comm, proc: %3d, tgt: %3d, iblk: %3d\n", bnum, iblk2, iblk);
            usend = (double *) malloc(ntrans * sizeof(double));
            urecv = (double *) malloc(ntrans * sizeof(double));

            int indjmp = IJ(NVAR - 1, ny - 2, NVAR);
            for (int j = 0; j < ny - 1; j++) {
                if (ACCUR == 0) {
                    for (int k = 0; k < NVAR; k++) {
                        int iu = IJK(0, j, k, nx - 1, NVAR);
                        usend[IJ(k, j, NVAR)] = unk[iu];
                    }
                } else {
                    int iu = IJK(0, j, 0, nx - 1, NVAR);
                    int ie = IJ(0, j, nx - 1);
                    double xsi  = -1.0;
                    double eta1 =  1.0 / sqrt(3.0);
                    double eta2 = -1.0 / sqrt(3.0);
                    double uPoint1[NVAR], uPoint2[NVAR];
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi, eta1, uPoint1);
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi, eta2, uPoint2);
                    for (int k = 0; k < NVAR; k++) {
                        usend[IJ(k, j, 2*NVAR)] = uPoint1[k];
                        usend[IJ(k+NVAR, j, 2*NVAR)] = uPoint2[k];
                    }
                }
            }

            MPI_Sendrecv(usend, ntrans, MPI_DOUBLE, iblk2,0,
                         urecv, ntrans, MPI_DOUBLE, iblk2,0,
                         MPI_COMM_WORLD, &status);

            if (ACCUR==0) {
                for (int j = 0; j < ntrans; j++) {
                    uGLeft[j] = urecv[j];
                }
                for (int j = 0; j < ny - 1; j++) {
                    LeftVar[j].UpdateState(air);
                }
            } else {
                for (int j = 0; j < ny - 1; j++) {
                    int iuEx1 = IJK(0,j,0,2,NVAR);
                    int ieEx1 =  IJ(0,j,2);
                    int iuEx2 = IJK(1,j,0,2,NVAR);
                    int ieEx2 =  IJ(1,j,2);

                    for (int k=0; k<NVAR; k++) {
                        uGLeft[iuEx1+k] = urecv[IJ(k, j, 2*NVAR)];
                        uGLeft[iuEx2+k] = urecv[IJ(k+NVAR, j, 2*NVAR)];
                    }
                    LeftVar[ieEx1].UpdateState(air);
                    LeftVar[ieEx2].UpdateState(air);
                }
            }

            blockflag[0] = 1;
            free(usend);
            free(urecv);
        }
        if ((iblk == bids[1] || (iblk==bnum && bids[1]>bnum)) && blockflag[1]==0 ) {
            // ~~~~~~~~~~~~~~~~~~ Right Boundary
            iblk2 = bids[1];
            int ntrans = NVAR*(ny-1);
            ntrans *= (ACCUR+1);
            //printf("right comm, proc: %3d, tgt: %3d, iblk: %3d\n", bnum, iblk2, iblk);
            usend = (double*)malloc(ntrans*sizeof(double));
            urecv = (double*)malloc(ntrans*sizeof(double));

            for (int j=0; j<ny-1; j++) {
                if (ACCUR == 0) {
                    for (int k = 0; k < NVAR; k++) {
                        usend[IJ(k,j,NVAR)] = unk[IJK(nx-2,j,k,nx-1,NVAR)];
                    }
                } else {
                    int iu = IJK(nx-2,j,0,nx-1,NVAR);
                    int ie = IJ(nx-2, j, nx - 1);
                    double xsi  =  1.0;
                    double eta1 =  1.0 / sqrt(3.0);
                    double eta2 = -1.0 / sqrt(3.0);
                    double uPoint1[NVAR], uPoint2[NVAR];
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi, eta1, uPoint1);
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi, eta2, uPoint2);
                    for (int k = 0; k < NVAR; k++) {
                        usend[IJ(k, j, 2*NVAR)] = uPoint1[k];
                        usend[IJ(k+NVAR, j, 2*NVAR)] = uPoint2[k];
                    }
                }
            }

            MPI_Sendrecv(usend, ntrans, MPI_DOUBLE, iblk2,0,
                         urecv, ntrans, MPI_DOUBLE, iblk2,0,
                         MPI_COMM_WORLD, &status);

            if (ACCUR==0) {
                for (int j = 0; j < ntrans; j++) {
                    uGRight[j] = urecv[j];
                }
                for (int j = 0; j < ny - 1; j++) {
                    RightVar[j].UpdateState(air);
                }
            } else {
                for (int j = 0; j < ny - 1; j++) {
                    int iuEx1 = IJK(0,j,0,2,NVAR);
                    int ieEx1 =  IJ(0,j,2);
                    int iuEx2 = IJK(1,j,0,2,NVAR);
                    int ieEx2 =  IJ(1,j,2);
                    for (int k = 0; k < NVAR; k++) {
                        uGRight[iuEx1+k] = urecv[IJ(k, j, 2 * NVAR)];
                        uGRight[iuEx2+k] = urecv[IJ(k + NVAR, j, 2 * NVAR)];
                    }
                    RightVar[ieEx1].UpdateState(air);
                    RightVar[ieEx2].UpdateState(air);
                }
            }

            blockflag[1] = 1;
            free(usend);
            free(urecv);
        }
        if ((iblk == bids[0] || (iblk==bnum && bids[0]>bnum)) && blockflag[2]==0 ) {
            // ~~~~~~~~~~~~~~~~~~~~ Bottom Boundary
            iblk2 = bids[0];
            int ntrans = NVAR*(nx-1);
            ntrans *= (ACCUR+1);
            //printf("bot   comm, proc: %3d, tgt: %3d, iblk: %3d\n", bnum, iblk2, iblk);
            usend = (double*)malloc(ntrans*sizeof(double));
            urecv = (double*)malloc(ntrans*sizeof(double));

            for (int i=0; i<nx-1; i++) {
                if (ACCUR == 0) {
                    for (int k = 0; k < NVAR; k++) {
                        int iu = IJK(i, 0, k, nx - 1, NVAR);
                        usend[IJ(k, i, NVAR)] = unk[iu];
                    }
                } else {
                    int iu = IJK(i, 0, 0, nx - 1, NVAR);
                    int ie = IJ(i, 0, nx - 1);
                    double eta  = -1.0;
                    double xsi1 =  1.0 / sqrt(3.0);
                    double xsi2 = -1.0 / sqrt(3.0);
                    double uPoint1[NVAR], uPoint2[NVAR];
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi1, eta, uPoint1);
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi2, eta, uPoint2);
                    for (int k = 0; k < NVAR; k++) {
                        usend[IJ(k, i, 2 * NVAR)] = uPoint1[k];
                        usend[IJ(k + NVAR, i, 2 * NVAR)] = uPoint2[k];
                    }
                }
            }
            MPI_Sendrecv(usend, ntrans, MPI_DOUBLE, iblk2,0,
                         urecv, ntrans, MPI_DOUBLE, iblk2,0,
                         MPI_COMM_WORLD, &status);

            if (ACCUR ==0) {
                for (int i = 0; i < ntrans; i++) {
                    uGBot[i] = urecv[i];
                }
                for (int i = 0; i < nx - 1; i++) {
                    BotVar[i].UpdateState(air);
                }
            } else {
                for (int i = 0; i < nx - 1; i++) {
                    int iuEx1 = IJK(0,i,0,2,NVAR);
                    int ieEx1 =  IJ(0,i,2);
                    int iuEx2 = IJK(1,i,0,2,NVAR);
                    int ieEx2 =  IJ(1,i,2);

                    for (int k=0; k<NVAR; k++) {
                        uGBot[iuEx1+k] = urecv[IJ(k, i, 2*NVAR)];
                        uGBot[iuEx2+k] = urecv[IJ(k+NVAR, i, 2*NVAR)];
                    }
                    BotVar[ieEx1].UpdateState(air);
                    BotVar[ieEx2].UpdateState(air);
                }
            }

            blockflag[2] = 1;
            free(usend);
            free(urecv);
        }
        if ((iblk == bids[2] || (iblk==bnum && bids[2]>bnum)) && blockflag[3]==0 ) {
            // Top Boundary
            iblk2 = bids[2];
            int ntrans = NVAR*(nx-1);
            ntrans *= (ACCUR+1);
            //printf("top   comm, proc: %3d, tgt: %3d, iblk: %3d\n", bnum, iblk2, iblk);
            usend = (double*)malloc(ntrans*sizeof(double));
            urecv = (double*)malloc(ntrans*sizeof(double));

            for (int i=0; i<nx-1; i++) {
                if (ACCUR == 0) {
                    for (int k = 0; k < NVAR; k++) {
                        int iu = IJK(i, ny - 2, k, nx - 1, NVAR);
                        usend[IJ(k, i, NVAR)] = unk[iu];
                    }
                } else {
                    int iu = IJK(i, ny-2, 0, nx - 1, NVAR);
                    int ie = IJ(i, ny-2, nx - 1);
                    double eta  =  1.0;
                    double xsi1 =  1.0 / sqrt(3.0);
                    double xsi2 = -1.0 / sqrt(3.0);
                    double uPoint1[NVAR], uPoint2[NVAR];
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi1, eta, uPoint1);
                    get_u_val(&(unk[iu]), ElemVar[ie], air, &(ux[iu]), &(uy[iu]), xsi2, eta, uPoint2);
                    for (int k = 0; k < NVAR; k++) {
                        usend[IJ(k, i, 2*NVAR)] = uPoint1[k];
                        usend[IJ(k+NVAR, i, 2*NVAR)] = uPoint2[k];
                    }
                }
            }
            MPI_Sendrecv(usend, ntrans, MPI_DOUBLE, iblk2,0,
                         urecv, ntrans, MPI_DOUBLE, iblk2,0,
                         MPI_COMM_WORLD, &status);

            if (ACCUR==0) {
                for (int i = 0; i < ntrans; i++) {
                    uGTop[i] = urecv[i];
                }
                for (int i = 0; i < nx - 1; i++) {
                    TopVar[i].UpdateState(air);
                }
            } else {
                for (int i = 0; i < nx - 1; i++) {
                    int iuEx1 = IJK(0,i,0,2,NVAR);
                    int ieEx1 =  IJ(0,i,2);
                    int iuEx2 = IJK(1,i,0,2,NVAR);
                    int ieEx2 =  IJ(1,i,2);
                    for (int k = 0; k < NVAR; k++) {
                        uGTop[iuEx1+k] = urecv[IJ(k, i, 2 * NVAR)];
                        uGTop[iuEx2+k] = urecv[IJ(k + NVAR, i, 2 * NVAR)];
                    }
                    TopVar[ieEx1].UpdateState(air);
                    TopVar[ieEx2].UpdateState(air);
                }
            }

            blockflag[3] = 1;
            free(usend);
            free(urecv);
        }
    }// iblk loop


    //====================Evaluate Flux Contributions====================
    //dudt = sum(flux_in * face_length) / volume
    //==========Fully Interior Faces


    // I|xi Fluxes
    for (int i=1; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            // ~~~~~~~~~~ Inviscid fluxes ~~~~~~~~~~
            //face above point i,j
            double len, fNormal[2], fflux[NVAR], rFace, *uLeft, *uRight, yCenter[2];
            State varL = State();
            State varR = State();
            len = geofa[IJK(i,j,3,nx,6)];
            fNormal[0] = geofa[IJK(i,j,4,nx,6)];
            fNormal[1] = geofa[IJK(i,j,5,nx,6)];
            if (IAXI==1) {
                rFace = yfa[IJK(i, j, 1, nx, 2)];
                yCenter[0] = geoel[IJK(i-1, j, 2, nx-1, 3)];
                yCenter[1] = geoel[IJK(i,   j, 2, nx-1, 3)];
            } else {
                rFace = 0.0;
                yCenter[0] = 0.0;
                yCenter[1] = 0.0;
            }

            int iuL = IJK(i-1,j,0,nx-1,NVAR);
            int iuR = IJK(i  ,j,0,nx-1,NVAR);
            int ieL = IJ(i-1, j, nx-1);
            int ieR = IJ(i,   j, nx-1);



            if (ACCUR==1) {
                DGP1_xsi_face_integral(ieL, ieR, iuL, iuR, unk, ElemVar, ux, uy, yCenter, air,
                                       rFace, fNormal, len, rhsel, rhselx, rhsely);
            } else {
                uLeft = &unk[iuL];
                uRight = &unk[iuR];
                varL = ElemVar[ieL];
                varR = ElemVar[ieR];

                LDFSS(fNormal[0], fNormal[1], len, rFace, uLeft, varL, uRight, varR, fflux, &parr);

                //Add flux contribution to elements
                for (int kvar=0; kvar<NVAR; kvar++) {
                    if (kvar == NSP + 1) { //Axisymmetric pressure correction
                        rhsel[iuL + kvar] -= (fflux[kvar] + (yCenter[0] * parr));
                        rhsel[iuR + kvar] += (fflux[kvar] + (yCenter[1] * parr));
                    } else {
                        rhsel[iuL + kvar]  -= fflux[kvar];
                        rhsel[iuR + kvar]  += fflux[kvar];
                    }
                }
            }

            if (IVISC==1) {
                // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
                double vflux[6], dc[2];
                //mirror the next interior cell to the boundary for that flux contrib
                dc[0] = geoel[IJK(i, j, 1, nx - 1, 3)] - geoel[IJK(i - 1, j, 1, nx - 1, 3)];
                dc[1] = geoel[IJK(i, j, 2, nx - 1, 3)] - geoel[IJK(i - 1, j, 2, nx - 1, 3)];
                varL = ElemVar[ieL];
                varR = ElemVar[ieR];
                viscous(nx, fNormal[1], fNormal[0], &(unk[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

                //Add flux contribution to elements
                rhsel[iuL + 1] += len * vflux[0];
                rhsel[iuL + 2] += len * vflux[1];
                rhsel[iuL + 3] += len * vflux[2];

                rhsel[iuR + 1] -= len * vflux[3];
                rhsel[iuR + 2] -= len * vflux[4];
                rhsel[iuR + 3] -= len * vflux[5];
            }

        }
    }


    // J|eta fluxes
    for (int i=0; i<nx-1; i++){
        for (int j=1; j<ny-1; j++){
            //faces to the left/above point i,j
            double len, fNormal[2], fflux[NVAR], rFace, *uLeft, *uRight, yCenter[2];
            State varL = State();
            State varR = State();
            len = geofa[IJK(i,j,0,nx,6)];
            fNormal[0] = geofa[IJK(i,j,1,nx,6)];
            fNormal[1] = geofa[IJK(i,j,2,nx,6)];
            if (IAXI) {
                rFace = yfa[IJK(i, j, 0, nx, 2)];
                yCenter[0] = geoel[IJK(i, j,   2, nx-1, 3)];
                yCenter[1] = geoel[IJK(i, j-1, 2, nx-1, 3)];
            } else {
                rFace = 1.0;
                yCenter[0] = 0.0;
                yCenter[1] = 0.0;
            }

            int iuL = IJK(i,j  ,0,nx-1,NVAR);
            int iuR = IJK(i,j-1,0,nx-1,NVAR);
            int ieL = IJ(i, j  , nx-1);
            int ieR = IJ(i, j-1, nx-1);



            if (ACCUR==1) {
                DGP1_eta_face_integral(ieL, ieR, iuL, iuR, unk, ElemVar, ux, uy, yCenter, air,
                                       rFace, fNormal, len, rhsel, rhselx, rhsely);
            } else {
                uLeft = &unk[iuL];
                uRight = &unk[iuR];
                varL = ElemVar[ieL];
                varR = ElemVar[ieR];

                LDFSS(fNormal[0], fNormal[1], len, rFace, uLeft, varL, uRight, varR, fflux, &parr);

                //Add flux contribution to elements
                for (int kvar=0; kvar<NVAR; kvar++) {
                    if (kvar == NSP + 1) { //Axisymmetric pressure correction
                        rhsel[iuL + kvar] -= (fflux[kvar] + (yCenter[0] * parr));
                        rhsel[iuR + kvar] += (fflux[kvar] + (yCenter[1] * parr));
                    } else {
                        rhsel[iuL + kvar]  -= fflux[kvar];
                        rhsel[iuR + kvar]  += fflux[kvar];
                    }
                }
            }

            if (IVISC==1) {
                // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
                double vflux[6], dc[2];
                //mirror the next interior cell to the boundary for that flux contrib
                dc[0] = geoel[IJK(i, j - 1, 1, nx - 1, 3)] - geoel[IJK(i, j, 1, nx - 1, 3)];
                dc[1] = geoel[IJK(i, j - 1, 2, nx - 1, 3)] - geoel[IJK(i, j, 2, nx - 1, 3)];
                varL = ElemVar[ieL];
                varR = ElemVar[ieR];
                viscous(nx, fNormal[0], fNormal[1], &(unk[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

                //Add flux contribution to elements
                rhsel[iuL + 1] += len * vflux[0];
                rhsel[iuL + 2] += len * vflux[1];
                rhsel[iuL + 3] += len * vflux[2];

                rhsel[iuR + 1] -= len * vflux[3];
                rhsel[iuR + 2] -= len * vflux[4];
                rhsel[iuR + 3] -= len * vflux[5];
            }
        }
    }



    for (int iu=0; iu<nelem*NVAR; iu++){
        if (__isnan(rhsel[iu])){
            printf("NAN RHSEL\n");
        }
    }


    //========== Boundary faces ==========

    double len, normx, normy, fflux[NVAR], yface, *uRight, *uLeft;
    State varL = State();
    State varR = State();

    // I|xi fluxes
    for (int j=0; j<ny-1; j++){
        //--------------------  LEFT BOUNDARY
        len   = geofa[IJK(0, j, 3, nx, 6)];
        normx = geofa[IJK(0, j, 4, nx, 6)];
        normy = geofa[IJK(0, j, 5, nx, 6)];
        yface = yfa[IJK(0, j, 1, nx, 2)];
        int iuL = IJ(0, j, NVAR);
        int iuR = IJK(0, j, 0, nx - 1, NVAR);
        int ieL = j;
        int ieR = IJ(0, j, nx - 1);
        //double ycL = 2.0 * geoel[IJK(0, j, 2, nx-1, 3)] - geoel[IJK(0, j+1, 2, nx-1, 3)]; //extrapolate to get gost y coord
        double ycR = geoel[IJK(0, j, 2, nx - 1, 3)];


        int iFaceType, ieEx, iuEx;
        iFaceType = 4; //number CCW starting from bot
        ieEx = IJ(0,j,2);
        iuEx = IJK(0,j,0,2,NVAR);
        double yCenter, fNormal[2]{ normx, normy},
                        fNormalL[2]{normx, normy},
                        fNormalR[2]{normx, normy};
        double rFace;
        if (IAXI==1) {
            yCenter = ycR;
            rFace = yface;
        } else {
            yCenter = 0.0;
            rFace = 0.0;
        }

        //find neighboring normal vectors
        if (j != 0) {
            fNormalL[0] = geofa[IJK(0, j - 1, 4, nx, 6)];
            fNormalL[1] = geofa[IJK(0, j - 1, 5, nx, 6)];
        }
        if (j!= ny-2){
            fNormalR[0] = geofa[IJK(0, j + 1, 4, nx, 6)];
            fNormalR[1] = geofa[IJK(0, j + 1, 5, nx, 6)];
        }

        if (ACCUR==1) {
            DGP1_boundary_face_integral(ieR, ieEx, iuR, iuEx, unk, ElemVar, ux, uy, iFaceType, uGLeft, LeftVar,
                                        yCenter, air, rFace, fNormal, fNormalL, fNormalR, len, rhsel, rhselx, rhsely);
        } else {
            uLeft = &uGLeft[iuL];
            uRight = &unk[iuR];
            varL = LeftVar[ieL];
            varR = ElemVar[ieR];

            LDFSS(normx, normy, len, rFace, uLeft, varL, uRight, varR, fflux, &parr);

            //Add flux contribution to elements
            for (int kvar=0; kvar<NVAR; kvar++) {
                if (kvar == NSP + 1) { //Axisymmetric pressure correction
                    rhsel[iuR + kvar] += (fflux[kvar] + (yCenter * parr));
                } else {
                    rhsel[iuR + kvar]  += fflux[kvar];
                }
            }
        }

        if(IVISC==1) {
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6]{};
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(1, j, 1, nx - 1, 3)] - geoel[IJK(0, j, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(1, j, 2, nx - 1, 3)] - geoel[IJK(0, j, 2, nx - 1, 3)];
            varL = LeftVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(uGLeft[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuR + 1] -= len * vflux[3];
            rhsel[iuR + 2] -= len * vflux[4];
            rhsel[iuR + 3] -= len * vflux[5];
        }




        //--------------------  RIGHT BOUNDARY
        len = geofa[IJK(nx - 1, j, 3, nx, 6)];
        normx = geofa[IJK(nx - 1, j, 4, nx, 6)];
        normy = geofa[IJK(nx - 1, j, 5, nx, 6)];
        yface = yfa[IJK(nx - 1, j, 1, nx, 2)];
        iuL = IJK(nx - 2, j, 0, nx - 1, NVAR);
        iuR = IJ(0, j, NVAR);
        ieL = IJ(nx - 2, j, nx - 1);
        ieR = j;
        double ycL = geoel[IJK(nx - 2, j, 2, nx - 1, 3)];
        //int iFaceType, ieEx, iuEx;
        iFaceType = 2; //number CCW starting from bot = 1
        ieEx = IJ(0,j,2);
        iuEx = IJK(0,j,0,2,NVAR);
        //double yCenter, fNormal[2]{normx, normy}, fNormalL[2]{normx, normy}, fNormalR[2]{normx, normy};
        fNormal[0] = normx;
        fNormal[1] = normy;

        //find neighboring normal vectors
        if (j != 0) {
            fNormalL[0] = geofa[IJK(nx - 1, j-1, 4, nx, 6)];
            fNormalL[1] = geofa[IJK(nx - 1, j-1, 5, nx, 6)];
        }
        if (j!= ny-2){
            fNormalR[0] = geofa[IJK(nx - 1, j+1, 4, nx, 6)];
            fNormalR[1] = geofa[IJK(nx - 1, j+1, 5, nx, 6)];
        }

        //double rFace;
        if (IAXI) {
            rFace = yface;
            yCenter = ycL;
        } else {
            rFace = 0.0;
        }

        if (ACCUR == 1) {
            DGP1_boundary_face_integral(ieL, ieEx, iuL, iuEx, unk, ElemVar, ux, uy, iFaceType, uGRight, RightVar,
                                        yCenter, air, rFace, fNormal, fNormalL, fNormalR, len, rhsel, rhselx, rhsely);
        } else {
            uLeft = &unk[iuL];
            uRight = &uGRight[iuR];
            varL = ElemVar[ieL];
            varR = RightVar[ieR];

            LDFSS(normx, normy, len, rFace, uLeft, varL, uRight, varR, fflux, &parr);

            //Add flux contribution to elements
            for (int kvar=0; kvar<NVAR; kvar++) {
                if (kvar == NSP + 1) { //Axisymmetric pressure correction
                    rhsel[iuL + kvar] -= (fflux[kvar] + (yCenter * parr));
                } else {
                    rhsel[iuL + kvar]  -= fflux[kvar];
                }
            }
        }

        if(IVISC==1) {
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            //mirror the next interior cell to the boundary for that flux contrib
            double dc[2], vflux[6]{};
            dc[0] = geoel[IJK(nx - 2, j, 1, nx - 1, 3)] - geoel[IJK(nx - 3, j, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(nx - 2, j, 2, nx - 1, 3)] - geoel[IJK(nx - 3, j, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = RightVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(uGRight[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL + 1] += len * vflux[0];
            rhsel[iuL + 2] += len * vflux[1];
            rhsel[iuL + 3] += len * vflux[2];
        }

    }


    // J|eta fluxes
    for (int i=0; i<nx-1; i++) {
        //BOTTOM BOUNDARY
        len   = geofa[IJK(i,0,0,nx,6)];
        normx = geofa[IJK(i,0,1,nx,6)];
        normy = geofa[IJK(i,0,2,nx,6)];
        yface = yfa[IJK(i,0,0,nx,2)];
        int iuL = IJK(i, 0, 0, nx - 1, NVAR);
        int iuR = IJ(0,i,NVAR);
        int ieL = IJ(i,0,nx-1);
        int ieR = i;
        double ycL = geoel[IJK(i, 0, 2, nx - 1, 3)];
        int iFaceType, ieEx, iuEx;
        iFaceType = 1; //number CCW starting from bot = 1
        ieEx = IJ(0,i,2);
        iuEx = IJK(0,i,0,2,NVAR);
        double yCenter, fNormal[2]{normx, normy}, fNormalL[2]{normx, normy}, fNormalR[2]{normx, normy};


        //find neighboring normal vectors
        if (i != 0) {
            fNormalL[0] = geofa[IJK(i-1,0,1,nx,6)];
            fNormalL[1] = geofa[IJK(i-1,0,2,nx,6)];
        }
        if (i!= nx-2){
            fNormalR[0] = geofa[IJK(i+1,0,1,nx,6)];
            fNormalR[1] = geofa[IJK(i+1,0,2,nx,6)];
        }

        double rFace;
        if (IAXI) {
            rFace = yface;
            yCenter = ycL;
        } else {
            rFace = 0.0;
            yCenter = 0.0;
        }

        if (ACCUR==1) {
            DGP1_boundary_face_integral(ieL, ieEx, iuL, iuEx, unk, ElemVar, ux, uy, iFaceType, uGBot, BotVar,
                                        yCenter, air, rFace, fNormal, fNormalL, fNormalR, len, rhsel, rhselx, rhsely);
            //printf("(i,j,) rhx,x,y: (%2d,%2d) %f,%f,%f\n\n",
            //       i,0,rhsel[iuL],rhselx[iuL],rhsely[iuL]);
            //printf("\n");
        } else {
            uLeft = &unk[iuL];
            uRight = &uGBot[iuR];
            varL = ElemVar[ieL];
            varR = BotVar[ieR];

            LDFSS(normx, normy, len, rFace, uLeft, varL, uRight, varR, fflux, &parr);

            //Add flux contribution to elements
            for (int kvar=0; kvar<NVAR; kvar++) {
                if (kvar == NSP + 1) { //Axisymmetric pressure correction
                    rhsel[iuL + kvar] -= (fflux[kvar] + (yCenter * parr));
                } else {
                    rhsel[iuL + kvar]  -= fflux[kvar];
                }
            }
        }
        if(IVISC==1) {
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i, 0, 1, nx - 1, 3)] - geoel[IJK(i, 1, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(i, 0, 2, nx - 1, 3)] - geoel[IJK(i, 1, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = BotVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(uGBot[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL + 1] += len * vflux[0];
            rhsel[iuL + 2] += len * vflux[1];
            rhsel[iuL + 3] += len * vflux[2];
        }

        //TOP BOUNDARY
        len = geofa[IJK(i, ny - 1, 0, nx, 6)];
        normx = geofa[IJK(i, ny - 1, 1, nx, 6)];
        normy = geofa[IJK(i, ny - 1, 2, nx, 6)];
        yface = yfa[IJK(i, ny - 1, 0, nx, 2)];

        iuL = IJ(0, i, NVAR);
        iuR = IJK(i, ny - 2, 0, nx - 1, NVAR);
        ieL = i;
        ieR = IJ(i, ny - 2, nx - 1);

        double ycR = geoel[IJK(i, ny - 2, 2, nx - 1, 3)];


        //int iFaceType, ieEx, iuEx;
        iFaceType = 3; //number CCW starting from bot = 1
        ieEx = IJ(0,i,2);
        iuEx = IJK(0,i,0,2,NVAR);
        //double yCenter, fNormal[2]{normx, normy}, fNormalL[2]{normx, normy}, fNormalR[2]{normx, normy};


        //find neighboring normal vectors
        fNormal[0] = normx;
        fNormal[1] = normy;
        if (i != 0) {
            fNormalL[0] = geofa[IJK(i-1, ny - 1, 1, nx, 6)];
            fNormalL[1] = geofa[IJK(i-1, ny - 1, 2, nx, 6)];
        }
        if (i!= nx-2){
            fNormalR[0] = geofa[IJK(i+1, ny - 1, 1, nx, 6)];
            fNormalR[1] = geofa[IJK(i+1, ny - 1, 2, nx, 6)];
        }

        //double rFace;
        if (IAXI) {
            rFace = yface;
            yCenter = ycR;
        } else {
            rFace = 0.0;
            yCenter = 0.0;
        }

        if (ACCUR==1) {
            DGP1_boundary_face_integral(ieR, ieEx, iuR, iuEx, unk, ElemVar, ux, uy, iFaceType, uGTop, TopVar,
                                        yCenter, air, rFace, fNormal, fNormalL, fNormalR, len, rhsel, rhselx, rhsely);
        } else {
            uLeft = &uGTop[iuL];
            uRight = &unk[iuR];
            varL = TopVar[ieL];
            varR = ElemVar[ieR];

            LDFSS(normx, normy, len, rFace, uLeft, varL, uRight, varR, fflux, &parr);

            //Add flux contribution to elements
            for (int kvar=0; kvar<NVAR; kvar++) {
                if (kvar == NSP + 1) { //Axisymmetric pressure correction
                    rhsel[iuR + kvar] += (fflux[kvar] + (yCenter * parr));
                } else {
                    rhsel[iuR + kvar]  += fflux[kvar];
                }
            }
        }

        /*
        //DG extension (eta flux on 'horizontal' face)
        double xsiR, etaR;
        xsiR = 0.0;
        etaR = 1.0;
        get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRight);
        varR.Initialize(uRight);
        varR.UpdateState(air);

        //Find interface flux
        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, len, yface, &(uGTop[iuL]), TopVar[i], uRight, varR, fflux, &parr);


        //Add flux contribution to elements
        rhsel[iuR] += fflux[0];
        rhsel[iuR + 1] += fflux[1];
        rhsel[iuR + 2] += (fflux[2] + (ycR * parr));
        rhsel[iuR + 3] += fflux[3];
         */


        if(IVISC==1) {
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            //mirror the next interior cell to the boundary for that flux contrib
            double dc[2], vflux[6]{};
            dc[0] = geoel[IJK(i, ny - 3, 1, nx - 1, 3)] - geoel[IJK(i, ny - 2, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(i, ny - 3, 2, nx - 1, 3)] - geoel[IJK(i, ny - 2, 2, nx - 1, 3)];
            varL = TopVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(uGTop[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuR + 1] -= len * vflux[3];
            rhsel[iuR + 2] -= len * vflux[4];
            rhsel[iuR + 3] -= len * vflux[5];
        }

        }

    //====================Combine to Find du/dt====================
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            int iu = IJK(i,j,0,nx-1,NVAR);
            double vol = geoel[IJK(i,j,0,nx-1,3)];
            dudt[iu  ] = rhsel[iu  ] / vol;
            dudt[iu+1] = rhsel[iu+1] / vol;
            dudt[iu+2] = rhsel[iu+2] / vol;
            dudt[iu+3] = rhsel[iu+3] / vol;

            //printf("(i,j,) rhx,x,y: (%2d,%2d) %f,%f,%f\n",
            //       i,j,rhsel[iu],rhselx[iu],rhsely[iu]);

            if (ACCUR==1) {
                duxdt[iu]     = 3.0 * rhselx[iu]     / vol;
                duxdt[iu + 1] = 3.0 * rhselx[iu + 1] / vol;
                duxdt[iu + 2] = 3.0 * rhselx[iu + 2] / vol;
                duxdt[iu + 3] = 3.0 * rhselx[iu + 3] / vol;

                duydt[iu]     = 3.0 * rhsely[iu]     / vol;
                duydt[iu + 1] = 3.0 * rhsely[iu + 1] / vol;
                duydt[iu + 2] = 3.0 * rhsely[iu + 2] / vol;
                duydt[iu + 3] = 3.0 * rhsely[iu + 3] / vol;
            }

            if (__isnan(dudt[iu]) or __isnan(dudt[iu+1]) or __isnan(dudt[iu+2]) or __isnan(dudt[iu+3])) {
                printf("dug\n");
            }

            ASSERT(!__isnan(dudt[iu]), "NaN Encounterd in Flux formulation")

        }
    }
    if (ACCUR==1) {
        // SINCE M IS DIAGONAL MATRIX, THE BELOW ALREADY INCLUDES THE MULTIPLE 3/VOL FROM ITS INVERSION
        DGP1_volume_integral(nx, ny, 1.0, xfa, yfa, geoel, unk, ElemVar, duxdt, duydt);
    }

    for (int i=0; i<nx-1; i++) {
        for (int j = 0; j < ny - 1; j++) {
            int iu = IJK(i, j, 0, nx - 1, NVAR);
            double vol = geoel[IJK(i, j, 0, nx - 1, 3)];

            //printf("(i,j,) dudt:dx,dy: (%2d,%2d) %f,%f\n",
            //       i, j, duxdt[iu], duydt[iu]);
        }
    }

    /*
    for (int iu=0; iu<nelem*NVAR ; iu++){
        if (dudt[iu] > 1e-6){
            printf("%le\n",dudt[iu]);
        }
    }
     */

    free(rhsel);
    if (ACCUR==1) {
        free(rhselx);
        free(rhsely);
    }

    free(BotVar);
    free(TopVar);
    free(RightVar);
    free(LeftVar);

    free(uGBot);
    free(uGRight);
    free(uGTop);
    free(uGLeft);
}
