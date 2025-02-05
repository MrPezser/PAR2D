//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_INDEXING_H
#define FVEULER_INDEXING_H
#include <cstdio>
#include <cstdlib>

#define IJ(i, j, ni)  (((j)*(ni)) + (i))
#define IJK(i, j, k, ni, nk)  ((((j)*(ni)) + (i))*(nk) + (k))
#define NVAR 4
#define NDEGR 1
#define NSP 1
#define IVISC 1
#define IGAM (1.4) //(-1.28)
#define ACCUR (0)
#define IAXI (0) ///2D WORKING SUPER GREAT
#define MXANGLE (6.0)
#define NITER (1e4)

#define sign(x)  ((std::signbit(x) ?  -1 : 1))
#define ASSERT(cond, msg) if(!(cond)){printf("Failed Assert: %s:%u %s\n %s\n", __FILE__, __LINE__, #cond, msg); exit(0);}
#define CHECKD(cond, msg, val) if(!(cond)){printf("Failed Assert: %s:%u %s\n %s\n Val:%lf\n", __FILE__, __LINE__, #cond, msg, val); exit(0);}
#define FORCEVALD(cond, msg, var, val) if(!(cond)){printf("Limiting Value: %s se to %lf\n",msg,val); var = val;}

#endif //FVEULER_INDEXING_H
