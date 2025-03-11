//
// Functions to aid in the addition of DG for higher order extension
// Created by Tsail on 5/29/2024.
//

#include "DGP1Tools.h"
#include "EulerFlux.h"
#include "BoundaryConditions.h"

void get_u_val(const double* unk, State var, Thermo air,const double* ux, const double* uy, double xsi, double eta, double* uout){
    if (ACCUR == 1) {
        for (int k = 0; k < NVAR; k++) {
            if (k==0) {

                //find pressure and derivatives
                double p = var.p;
                double dpdxsi, dpdeta, pout;
                dpdxsi = ux[0] * air.Rs[0] * unk[3] +
                         unk[0] * air.Rs[0] * ux[3];  // = (dp/drho)(drho/dxsi) + (dp/dT)(dT/dxsi)
                dpdeta = uy[0] * air.Rs[0] * unk[3] +
                         unk[0] * air.Rs[0] * uy[3];  // = (dp/drho)(drho/deta) + (dp/dT)(dT/deta)

                //extrapolate with pressure instead of density
                pout = p + (xsi*dpdxsi) + (eta*dpdeta);

                uout[k] = pout / (air.Rs[0]*(unk[3] + (xsi * ux[3]) + (eta * uy[3])));

            } else{
                uout[k] = unk[k] + (xsi * ux[k]) + (eta * uy[k]);
            }
        }
    } else {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k];
        }
    }
}

void get_u_val_standardrecon(const double* unk,const double* ux, const double* uy, double xsi, double eta, double* uout){
    if (ACCUR == 1) {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k] + (xsi * ux[k]) + (eta * uy[k]);
        }
    } else {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k];
        }
    }
}


void DGP1_volume_integral(int nx, int ny, double vol, double* xfa, double* yfa, double* geoel, double* unk, State* ElemVar,
                          double* duxdt, double* duydt){
    // NOTE: Includes axisymetric flux modification
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++) {
            int iu = IJK(i,j,0,nx-1,NVAR);
            int iel = IJ(i,j,nx-1);

            //Calculate cell centered flux
            double Fx[NVAR], Fy[NVAR];
            double* unkij = &(unk[iu]);
            State var = ElemVar[iel];

            //Find Flux Vector
            double ycc;
            if (IAXI==1) {
                ycc = geoel[IJK(i,j,2,nx-1, 3)];
            } else {
                ycc = 1.0;
            }
            double rho = unkij[0];
            Fx[0] = ycc *  rho * unkij[1];
            Fx[1] = ycc * (rho * unkij[1] * unkij[1] + var.p);
            Fx[2] = ycc *  rho * unkij[1] * unkij[2];
            Fx[3] = ycc *  rho * unkij[1] * var.h0;//unkij[1] * (rho*var.e + var.p);

            Fy[0] = ycc *  rho * unkij[2];
            Fy[1] = ycc *  rho * unkij[1] * unkij[2];
            Fy[2] = ycc * (rho * unkij[2] * unkij[2] + var.p);
            Fy[3] = ycc *  rho * unkij[2] * var.h0;//unkij[2] * (rho*var.e + var.p);

            // Derivatives of Basis Functions
            double xL,xR,xD,xU,yL,yR,yD,yU;
            double db1dx, db1dy, db2dx,db2dy;
            xL = xfa[IJK(i,j,1,nx,2)];
            xR = xfa[IJK(i+1,j,1,nx,2)];
            xD = xfa[IJK(i,j,0,nx,2)];
            xU = xfa[IJK(i,j+1,0,nx,2)];

            yL = yfa[IJK(i,j,1,nx,2)];
            yR = yfa[IJK(i+1,j,1,nx,2)];
            yD = yfa[IJK(i,j,0,nx,2)];
            yU = yfa[IJK(i,j+1,0,nx,2)];

            //d(xsi)/dx
            //d(xsi)/dy
                //xsi plane normal
            double nxsi[3];
            nxsi[0] = -0.5*(yU-yD);
            nxsi[1] =  0.5*(xU-xD);
            nxsi[2] =  -0.25 * ((xU-xD) * (yL-yR) - (yU-yD) * (xL-xR));
            db1dx = nxsi[0] / nxsi[2];
            db1dy = nxsi[1] / nxsi[2];

            //d(eta)/dx
            //d(eta)/dy
                //eta plane normal
            double neta[3];
            neta[0] = -0.5*(yL-yR);
            neta[1] =  0.5*(xL-xR);
            neta[2] =  -0.25 * ((xL-xR) * (yD-yU) - (yL-yR) * (xD-xU));
            db2dx = neta[0] / neta[2];
            db2dy = neta[1] / neta[2];

            //Assemble volume contribution
            for (int kvar=0; kvar<NVAR; kvar++){
                duxdt[iu+kvar] += 3.0 * (Fx[kvar]*db1dx + Fy[kvar]*db1dy);
                duydt[iu+kvar] += 3.0 * (Fx[kvar]*db2dx + Fy[kvar]*db2dy);
            }

        }
    }
}


void DGP1_xsi_face_integral(int ieL, int ieR, int iuL, int iuR,double* unk, State* ElemVar, double* ux, double* uy,
                        const double* yCenter, Thermo air, double rFace, double* fNormal, double len,
                        double* rhsel, double* rhselx, double* rhsely){

    //Input left and right variable/state information
    //Output addition of flux contribution to respective elemets
    double fflux[NVAR], uLFace[NVAR], uRFace[NVAR], parr;
    State varL = State();
    State varR = State();
    varL.Initialize(uLFace);
    varR.Initialize(uRFace);

    //Cell centers used for adding in pressure term
    double ycL, ycR;

    if (IAXI == 1) {
        ycL = yCenter[0];
        ycR = yCenter[1];
    } else {
        ycL = 0.0;
        ycR = 0.0;
        rFace = 0.0;
    }

    //DG extension (xsi flux on vertical face)
    double weight = 0.5;
    double xsiL, etaL, xsiR, etaR;

    //1st point
    xsiL = 1.0;
    etaL = 1.0/sqrt(3.0);
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = -1.0;
    etaR = 1.0/sqrt(3.0);
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!__isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

    //2nd point
    xsiL = 1.0;
    etaL = -1.0/sqrt(3.0);
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = -1.0;
    etaR = -1.0/sqrt(3.0);
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!__isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

}

void DGP1_eta_face_integral(int ieL, int ieR, int iuL, int iuR,double* unk, State* ElemVar, double* ux, double* uy,
                            const double* yCenter, Thermo air, double rFace, double* fNormal, double len,
                            double* rhsel, double* rhselx, double* rhsely){
    //Input left and right variable/state information
    //Output addition of flux contribution to respective elemets
    double fflux[NVAR], uLFace[NVAR], uRFace[NVAR], parr;
    State varL = State();
    State varR = State();
    varL.Initialize(uLFace);
    varR.Initialize(uRFace);

    //Cell centers used for adding in pressure term
    double ycL, ycR;
    if (IAXI) {
        ycL = yCenter[0];
        ycR = yCenter[1];
    } else {
        ycL = 0.0;
        ycR = 0.0;
        rFace = 0.0;
    }

    //DG extension (xsi flux on vertical face)
    double weight = 0.5;
    double xsiL, etaL, xsiR, etaR;

    //1st point
    xsiL = 1.0/sqrt(3.0);
    etaL = -1.0;
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = 1.0/sqrt(3.0);
    etaR = 1.0;
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!__isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

    //2nd point
    xsiL = -1.0/sqrt(3.0);
    etaL = -1.0;
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = -1.0/sqrt(3.0);
    etaR = 1.0;
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!__isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

}

void DGP1_boundary_face_integral(int ieIn, int ieEx, int iuIn, int iuEx,double* unk, State* ElemVar, double* ux, double* uy,
                            int iFaceType, double* unkExt, State* EVExt, double yCenter, Thermo air, double rFace,
                            double* fNormal, double* fNormalL, double* fNormalR, double len,
                            double* rhsel, double* rhselx, double* rhsely){
    //Input left and right variable/state information
    //Output addition of flux contribution to respective elemets
    double fflux[NVAR], uInFace[NVAR], *uExFace, parr;
    double normal[2];
    uExFace = &(unkExt[iuEx]);
    State varIn = State();
    State varEx = EVExt[ieEx];
    varIn.Initialize(uInFace);
    //varEx.Initialize(uExFace);
    //varEx.UpdateState(air);

    ///////////////////////////////////////////////////////
    //fNormalL = fNormal;
    //fNormalR = fNormal;
    ///////////////////////////////////////////////////////

    //Cell centers used for adding in pressure term
    double ycIn;
    if (IAXI==1) {
        ycIn = yCenter;
    } else {
        ycIn = 0.0;
    }

    //Stuff for the type of face
    double sideMult{0.0}; //multiplier to correct for is the face is a right face or a left face according to convention.
    switch (iFaceType) {
        case 1 : {//horizontal face - bottom boundary
            sideMult = -1.0;
            break;
        }
        case 2 : {//vertical face - right boundary
            sideMult = -1.0;
            break;
        }
        case 3 : {//horizontal face - top boundary
            sideMult = 1.0;
            break;
        }
        case 4 : {//vertical face - left boundary
            sideMult = 1.0;
            break;
        }
        default: {
            ASSERT(fabs(sideMult) > 0.0, "Invalid Face Type")
        }
    }

    //DG extension (xsi flux on vertical face)
    double weight = 0.5;
    double xsi{}, eta{};

    //1st point
    switch (iFaceType) {
        case 1 : {//horizontal face - bottom boundary
            xsi = 1.0 / sqrt(3.0);
            eta = -1.0;
            double temp = fabs(xsi);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalR[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalR[1];
            break;
        }
        case 2 : {//vertical face - right boundary
            xsi = 1.0;
            eta = 1.0 / sqrt(3.0);
            double temp = fabs(eta);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalR[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalR[1];
            break;
        }
        case 3 : {//horizontal face - top boundary
            xsi = 1.0 / sqrt(3.0);
            eta = 1.0;
            double temp = fabs(xsi);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalR[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalR[1];
            break;
        }
        case 4 : {//vertical face - left boundary
            xsi = -1.0;
            eta = 1.0 / sqrt(3.0);
            double temp = fabs(eta);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalR[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalR[1];
            break;
        }
        default: {
            ASSERT(iFaceType >= 1 and iFaceType <= 4, "Invalid Face Type")
        }
    }

    //Get the interior values at the quadrature point
    get_u_val(&(unk[iuIn]), ElemVar[ieIn], air, &(ux[iuIn]), &(uy[iuIn]), xsi, eta, uInFace);
    varIn.UpdateState(air);

    //Find interface flux
    if (iFaceType == 1 or iFaceType == 2) {
        LDFSS(normal[0], normal[1], len, rFace, uInFace, varIn,
              uExFace, varEx, fflux, &parr);
    } else {
        ASSERT(iFaceType == 3 or iFaceType == 4, "Invalid iFaceType")
        LDFSS(normal[0], normal[1], len, rFace, uExFace, varEx,
              uInFace, varIn, fflux, &parr);
    }

    if (iFaceType == 10) {
        printf("(in f1) rhx,x,y: %f,%f,%f\n",
               rhsel[iuIn], rhselx[iuIn], rhsely[iuIn]);
    }
    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuIn + kvar]  += sideMult * weight * (fflux[kvar] + (ycIn*parr));
            rhselx[iuIn + kvar] += sideMult * weight * (fflux[kvar] + (ycIn*parr)) * xsi;
            rhsely[iuIn + kvar] += sideMult * weight * (fflux[kvar] + (ycIn*parr)) * eta;
        } else {
            rhsel[iuIn + kvar]  += sideMult * weight * fflux[kvar];
            rhselx[iuIn + kvar] += sideMult * weight * (fflux[kvar]) * xsi;
            rhsely[iuIn + kvar] += sideMult * weight * (fflux[kvar]) * eta;
        }
    }

    ASSERT(!__isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

    //2nd point
    ieEx += IJ(1,0,2);
    iuEx += IJK(1,0,0,2,NVAR);

    uExFace = &(unkExt[iuEx]);
    varEx = EVExt[ieEx];
    //varEx.Initialize(uExFace);
    //varEx.UpdateState(air);

    switch (iFaceType) {
        case 1 : {//horizontal face - bottom boundary
            xsi = -1.0 / sqrt(3.0);
            eta = -1.0;

            double temp = fabs(xsi);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalL[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalL[1];
            break;
        }
        case 2 : {//vertical face - right boundary
            xsi = 1.0;
            eta = -1.0 / sqrt(3.0);

            double temp = fabs(eta);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalL[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalL[1];
            break;
        }
        case 3 : {//horizontal face - top boundary
            xsi = -1.0 / sqrt(3.0);
            eta = 1.0;

            double temp = fabs(xsi);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalL[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalL[1];
            break;
        }
        case 4 : { //vertical face - left boundary
            xsi = -1.0;
            eta = -1.0 / sqrt(3.0);

            double temp = fabs(eta);
            normal[0] = (1.0 - 0.5*temp)*fNormal[0] + 0.5*temp*fNormalL[0];
            normal[1] = (1.0 - 0.5*temp)*fNormal[1] + 0.5*temp*fNormalL[1];
            break;
        }
        default: {
            ASSERT(iFaceType >= 1 and iFaceType <= 4, "Invalid Face Type")
        }
    }

    //Get the interior values at the quadrature point
    get_u_val(&(unk[iuIn]), ElemVar[ieIn], air, &(ux[iuIn]), &(uy[iuIn]), xsi, eta, uInFace);
    varIn.UpdateState(air);

    //Find interface flux
    if (iFaceType == 1 or iFaceType == 2) {
        LDFSS(normal[0], normal[1], len, rFace, uInFace, varIn,
              uExFace, varEx, fflux, &parr);
    } else {
        ASSERT(iFaceType == 3 or iFaceType == 4, "Invalid iFaceType")
        LDFSS(normal[0], normal[1], len, rFace, uExFace, varEx,
              uInFace, varIn, fflux, &parr);
    }

    //Add flux contribution to elements
    if (iFaceType == 10) {
        printf("(in f1) rhx,x,y: %f,%f,%f\n",
               rhsel[iuIn], rhselx[iuIn], rhsely[iuIn]);
    }
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuIn + kvar]  += sideMult * weight * (fflux[kvar] + (ycIn*parr));
            rhselx[iuIn + kvar] += sideMult * weight * (fflux[kvar] + (ycIn*parr)) * xsi;
            rhsely[iuIn + kvar] += sideMult * weight * (fflux[kvar] + (ycIn*parr)) * eta;
        } else {
            rhsel[iuIn + kvar]  += sideMult * weight * fflux[kvar];
            rhselx[iuIn + kvar] += sideMult * weight * (fflux[kvar]) * xsi;
            rhsely[iuIn + kvar] += sideMult * weight * (fflux[kvar]) * eta;
        }
    }
    if (iFaceType == 10) {
        printf("(in f2) rhx,x,y: %f,%f,%f\n",
               rhsel[iuIn], rhselx[iuIn], rhsely[iuIn]);
    }
    ASSERT(!__isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")
}

void get_boundary_point(int btype, double normx, double normy, double* uFS, double* unkiint, State EViel, Thermo& air,
                        double* uxiint, double* uyiint, double xsi, double eta, State& ElemVarieex, double* uGiuex){
    double unkelij[NVAR];
    State varij = State();

    get_u_val(unkiint, EViel, air, uxiint, uyiint, xsi, eta, unkelij);
    varij.Initialize(unkelij);
    varij.UpdateState(air);
    //==========Ghost State
    boundary_state(btype,air,normx,normy,uFS,unkelij, varij,
                   uGiuex);
    ElemVarieex.Initialize(uGiuex);
    ElemVarieex.UpdateState(air);
}

void DGP1_ghost_cell_generator(int nx, int ny, double* unk, double* ux, double* uy, State* ElemVar, Thermo air, int* ibound,
                               double* geofa, double* uFS, double* uGBot, double* uGTop, double* uGLeft, double* uGRight,
                               State* BotVar, State* TopVar, State* LeftVar, State* RightVar){
    double unkelij[NVAR];
    State varij = State();

    //bottom side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[i];
        int iuint = IJK(i, 0, 0, nx - 1, NVAR);
        int iel = IJ(i, 0, nx-1);
        int iuEx = IJK(0,i,0,2,NVAR);
        int ieEx = IJ(0,i,2);

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(i, 0, 1,nx,6)];
        normy = geofa[IJK(i, 0, 2,nx,6)];

        //DG extension
        double xsi, eta;
        //point 1
        xsi = 1.0 / sqrt(3.0);
        eta = -1.0;

        double normR[2]{normx, normy};
        if (i != nx-2){
            normR[0] = geofa[IJK(i+1, 0, 1,nx,6)];
            normR[1] = geofa[IJK(i+1, 0, 2,nx,6)];
        }
        normx = normx*(1.0 - 0.5*xsi) + 0.5*xsi*normR[0];
        normy = normy*(1.0 - 0.5*xsi) + 0.5*xsi*normR[1];


        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, BotVar[ieEx], &(uGBot[iuEx]));


        //point 2
        iuEx = IJK(1,i,0,2,NVAR);
        ieEx =  IJ(1,i,2);
        xsi = -1.0 / sqrt(3.0);
        eta = -1.0;

        normx = geofa[IJK(i, 0, 1,nx,6)];
        normy = geofa[IJK(i, 0, 2,nx,6)];
        double normL[2]{normx, normy};
        if (i != 0){
            normL[0] = geofa[IJK(i-1, 0, 1,nx,6)];
            normL[1] = geofa[IJK(i-1, 0, 2,nx,6)];
        }
        normx = normx*(1.0 + 0.5*xsi) - 0.5*xsi*normL[0];
        normy = normy*(1.0 + 0.5*xsi) - 0.5*xsi*normL[1];


        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, BotVar[ieEx], &(uGBot[iuEx]));

    }
    //right side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[j+ nx-1];
        int iuint = IJK(nx - 2, j, 0, nx - 1, NVAR);
        int iel = IJ(nx-2, j, nx-1);
        int iuEx = IJK(0,j,0,2,NVAR);
        int ieEx = IJ(0,j,2);

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(nx-1, j, 4,nx,6)];
        normy = geofa[IJK(nx-1, j, 5,nx,6)];

        //DG extension
        double xsi, eta;
        //point 1
        xsi = 1.0;
        eta = 1.0 / sqrt(3.0);

        double normR[2]{normx, normy};
        if (j != ny-2){
            normR[0] = geofa[IJK(nx-1, j+1, 4,nx,6)];
            normR[1] = geofa[IJK(nx-1, j+1, 5,nx,6)];
        }
        normx = normx*(1.0 - 0.5*eta) + 0.5*eta*normR[0];
        normy = normy*(1.0 - 0.5*eta) + 0.5*eta*normR[1];


        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, RightVar[ieEx], &(uGRight[iuEx]));
        //point 2
        iuEx = IJK(1,j,0,2,NVAR);
        ieEx = IJ(1,j,2);
        xsi = 1.0;
        eta = -1.0 / sqrt(3.0);

        normx = geofa[IJK(nx-1, j, 4,nx,6)];
        normy = geofa[IJK(nx-1, j, 5,nx,6)];
        double normL[2]{normx, normy};
        if (j != 0){
            normL[0] = geofa[IJK(nx-1, j-1, 4,nx,6)];
            normL[1] = geofa[IJK(nx-1, j-1, 5,nx,6)];
        }
        normx = normx*(1.0 + 0.5*eta) - 0.5*eta*normL[0];
        normy = normy*(1.0 + 0.5*eta) - 0.5*eta*normL[1];


        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, RightVar[ieEx], &(uGRight[iuEx]));
    }

    //top side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        int ib = nx-2-i;
        btype = ibound[ib+nx+ny-2];
        int iuint = IJK(i, ny - 2, 0, nx - 1, NVAR);
        int iel = IJ(i, ny-2, nx-1);
        int iuEx = IJK(0,i,0,2,NVAR);
        int ieEx = IJ(0,i,2);

        double normx, normy;
        normx = -geofa[IJK(i, ny-1, 1,nx,6)];
        normy = -geofa[IJK(i, ny-1, 2,nx,6)];

        //DG extension
        double xsi, eta;
        //point 1
        xsi = 1.0 / sqrt(3.0);
        eta = 1.0;

        double normR[2]{normx, normy};
        if (i != nx-2){
            normR[0] = -geofa[IJK(i+1, ny-1, 1,nx,6)];
            normR[1] = -geofa[IJK(i+1, ny-1, 2,nx,6)];
        }
        normx = normx*(1.0 - 0.5*xsi) + 0.5*xsi*normR[0];
        normy = normy*(1.0 - 0.5*xsi) + 0.5*xsi*normR[1];

        //printf("p1 nx,ny = %f,%f\t\t", normx, normy);


        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, TopVar[ieEx], &(uGTop[iuEx]));
        //point 2
        iuEx = IJK(1,i,0,2,NVAR);
        ieEx =  IJ(1,i,2);
        xsi = -1.0 / sqrt(3.0);
        eta = 1.0;


        normx = -geofa[IJK(i, ny-1, 1,nx,6)];
        normy = -geofa[IJK(i, ny-1, 2,nx,6)];
        double normL[2]{normx, normy};
        if (i != 0){
            normL[0] = -geofa[IJK(i-1, ny-1, 1,nx,6)];
            normL[1] = -geofa[IJK(i-1, ny-1, 2,nx,6)];
        }
        normx = normx*(1.0 + 0.5*xsi) - 0.5*xsi*normL[0];
        normy = normy*(1.0 + 0.5*xsi) - 0.5*xsi*normL[1];

        //printf("p2 nx,ny = %f,%f\n", normx, normy);

        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, TopVar[ieEx], &(uGTop[iuEx]));
    }

    //left side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        int jb = (ny-2)-j;
        btype = ibound[jb+(2*nx)+ny-3];
        int iuint = IJK(0, j, 0, nx - 1, NVAR);
        int iel = IJ(0, j, nx-1);
        int iuEx = IJK(0,j,0,2,NVAR);
        int ieEx = IJ(0,j,2);

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(0, j, 4,nx,6)];
        normy = -geofa[IJK(0, j, 5,nx,6)];

        //DG extension
        double xsi, eta;
        //point 1
        xsi = -1.0;
        eta = 1.0 / sqrt(3.0);

        double normR[2]{normx, normy};
        if (j != ny-2){
            normR[0] = -geofa[IJK(0, j+1, 4,nx,6)];
            normR[1] = -geofa[IJK(0, j+1, 5,nx,6)];
        }
        normx = normx*(1.0 - 0.5*eta) + 0.5*eta*normR[0];
        normy = normy*(1.0 - 0.5*eta) + 0.5*eta*normR[1];

        double uFSangle[4];
        double angle = (MXANGLE * M_PI / 180.0) * (double)j / (double)(ny-2);
        uFSangle[0] = uFS[0];
        uFSangle[1] = uFS[1]*cos(angle);
        uFSangle[2] = uFS[1]*sin(angle);
        uFSangle[3] = uFS[3];

        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, LeftVar[ieEx], &(uGLeft[iuEx]));
        //point 2
        iuEx = IJK(1,j,0,2,NVAR);
        ieEx = IJ(1,j,2);
        xsi = -1.0;
        eta = -1.0 / sqrt(3.0);

        normx = -geofa[IJK(0, j, 4,nx,6)];
        normy = -geofa[IJK(0, j, 5,nx,6)];
        double normL[2]{normx, normy};
        if (j != 0){
            normL[0] = -geofa[IJK(0, j-1, 4,nx,6)];
            normL[1] = -geofa[IJK(0, j-1, 5,nx,6)];
        }
        normx = normx*(1.0 + 0.5*eta) - 0.5*eta*normL[0];
        normy = normy*(1.0 + 0.5*eta) - 0.5*eta*normL[1];

        angle = (MXANGLE * M_PI / 180.0) * (double)j / (double)(ny-2);
        uFSangle[0] = uFS[0];
        uFSangle[1] = uFS[1]*cos(angle);
        uFSangle[2] = uFS[1]*sin(angle);
        uFSangle[3] = uFS[3];

        get_boundary_point(btype, normx, normy, uFS, &(unk[iuint]), ElemVar[iel], air, &(ux[iuint]),
                           &(uy[iuint]), xsi, eta, LeftVar[ieEx], &(uGLeft[iuEx]));
    }
}


void DGP1_DDG_viscous(double* uL, double* uR, State VarR, State VarL, double* uxL, double* uyL, double* uxR, double* uyR,
                            const double* yCenter, Thermo air, double rFace, double* fNormal, double len, double dc
                            double* vflx){
    // Code for computing the viscous interface flux according to the direct discontinuous galerkin method (Liu&Yan) 

    double mu = 0.5 * (varL.mu + varR.mu);

    dci = 1.0 / sqrt(dc[0] * dc[0] + dc[1] * dc[1]);

    //  ========== Convert primative variables at point to conserved vars
    double cvl[NVAR], cvr[NVAR];
   
   // Species density fractions remain 
    for (int isp=0; isp<NSP; isp++) {
       cvl[isp] = uL[isp];
       cvr[isp] = ur[isp];
    }
    //
    // velocities ==>> momentums
    cvl[NSP] = uL[NSP] * VarL.rho; 
    cvr[NSP] = ur[NSP] * VarR.rho; 
    cvl[NSP+1] = uL[NSP+1] * VarL.rho; 
    cvr[NSP+1] = ur[NSP+1] * VarR.rho; 
    //
    // temp ==>> energy
    cvl[NSP+2] = varL.e;
    cvr[NSP+2] = varR.e;
    //
    //
    //  ========== Compute the interface 1st derivatives DDG(CV)
    // dcv_dx_int = B0*(jump(u))/delta + (average dx) + B1*delta*(jump(ddudxdx))
    double cvi[NVAR];
    

}




















