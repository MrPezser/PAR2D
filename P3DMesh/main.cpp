#include <stdio.h>


/// File for transforming Plot3d file format into structured coordinate list format

int main(int argc, char *argv[]) {
	
	if (argc <= 1) {printf("Please provide file name.\n"); return 1;}

	FILE* fin = fopen(argv[1]);

	int nblock;
	fscanf(fin,"%d",&nblock);

	if (nblock > 1) {printf("Code in multiblock support you lazy fool.\n"); return 1;}

	int ni, nj, nk;
	fscanf(fin,"%d %d %d",&ni, &nj, &nk);

	int npoin = ni * nj * nk;
	int nline = ceil(npoin * 0.25);
	auto* x = (double*)malloc(npoin*sizeof(double));
	auto* y = (double*)malloc(npoin*sizeof(double));
	auto* z = (double*)malloc(npoin*sizeof(double));
	
	int ip=0;
	for (int iline=0; iline<nline; iline++) {

		if (iline == nline-1){
	        int p1,p2,p3,p4;
		while (ip < npoin){
		

		}
		

		fscanf(fin,"%d %d %d %d",x[ip],x[ip+1],x[ip+2],x[ip+3]);
		ip += 4;
	}

	for (int k=0,k<ni,k++) {
	for (int j=0,j<nj,k++) {
	for (int i=0,i<nk,k++) {

		

	}
	}
	}

	return 0;
}
