#include <stdio.h>
#include <cmath>

/// File for transforming Plot3d file format into structured coordinate list format

int main(int argc, char *argv[]) {
	
	if (argc <= 1) {printf("Please provide input and output file name.\n"); return 1;}

    printf("Reading in file: %s\n",argv[1]);

	FILE* fin = fopen(argv[1],"r");

	int nblock;
	fscanf(fin,"%d",&nblock);
    printf("Found %d blocks.\n",nblock);

	if (nblock > 1) {printf("Code in multiblock support you lazy fool.\n"); return 1;}

	int ni, nj, nk;
	fscanf(fin,"%d %d %d",&ni, &nj, &nk);
    printf("Block dim = (%d, %d, %d)\nProceeding to read in coordinate data.\n", ni, nj, nk);
    if (nk > 1) {printf("\nWARNING\n!!!!!!!!!!!!!!!!!!!!\nThis is a 3D mesh!!!\n!!!!!!!!!!!!!!!!!!\n\n");}

	int npoin = ni * nj * nk;
	int nline = ceil(npoin * 0.25);
	auto* xyz = (double*)malloc(3*npoin*sizeof(double));
	double *x, *y, *z;
	
    double row[4]{-999.0};
    int ind = 0;
	while (true) {
        int asdf = (fscanf(fin,"%e %e %e %e",&(row[0]),&(row[1]),&(row[2]),&(row[3])) >= 1);
        if (asdf < 1) {break;}
        //printf("%d\n",asdf);
        //
        for (int i=0; i<4; i++) {
        //
            if (ind+i < 3*npoin) {
                xyz[ind+i] = row[i];
            } else {
                break;
            }
        }
        //
    }

    printf("Finished Reading in from file\n");
    
    x = xyz;
    y = &xyz[npoin];
    //z = &xyz[2*npoin];

    printf("Writing to File: %s\n", argv[2]);
    FILE* fout = fopen(argv[2],"w");
    fprintf(fout, "%d %d\n",ni,nj);

    for (int j=0; j<nj; j++) {
    for (int i=0; i<ni; i++) {
        int ind = nj*i + j; ///////////////////////////////////////////////////////////////////////////////////////
        fprintf(fout,"%.15le %.15le\n", x[ind], y[ind]);
    }
    }
        
    fclose(fin);
    fclose(fout);
	return 0;
}






