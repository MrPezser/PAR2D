//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_FILEIO_H
#define FVEULER_FILEIO_H

#include "Thermo.h"

void read_mesh(int bnum, int* nx, int* ny, int** ibound, double** x, double** y);
void print_elem_stats(const char *title, int nx, int ny, const double* geoel);
void print_state(int iter, int bnum, const char *title, int nx, int ny, Thermo& air, double* x, double* y, double* unk, double* geoel );
void print_state_axi(const char *title, int nx, int ny, Thermo& air, double* x, double* y, double* unk, double* geoel );
void print_state_DGP1(double t, int iter, int bnum, const char *title, int nx, int ny, Thermo& air, double* x, double* y, double* unk, double* ux, double* uy, double* geoel );

#endif //FVEULER_FILEIO_H
