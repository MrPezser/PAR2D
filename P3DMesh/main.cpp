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
	fscanf(fin,"%d %d %d",&nj, &ni, &nk);
    printf("Block dim = (%d, %d, %d)\nProceeding to read in coordinate data.\n", ni, nj, nk);
    if (nk > 1) {printf("\nWARNING\n!!!!!!!!!!!!!!!!!!!!\nThis is a 3D mesh!!!\n!!!!!!!!!!!!!!!!!!\n\n");}

	int npoin = ni * nj * nk;
	int nline = ceil(npoin * 0.25);
	auto* xyz = (double*)malloc(3*npoin*sizeof(double));
	double *x, *y, *z;
	
    double row0[1],row1[1],row2[1],row3[1];
    int ind = 0;
	while (true) {
        //
        int asdf = fscanf(fin,"%le %le %le %le%*[^\n]",row0,row1,row2,row3);
        if (asdf < 1) {break;}
        printf("%d\t%e\n",asdf,row1[0]);
        //
        if (ind < 3*npoin) {
            xyz[ind] = row0[0];
        } else {
            break;
        }
        //
        if (ind+1 < 3*npoin) {
            xyz[ind+1] = row1[0];
        } else {
            break;
        }
        //
        if (ind+2 < 3*npoin) {
            xyz[ind+2] = row2[0];
        } else {
            break;
        }
        //
        if (ind+3 < 3*npoin) {
            xyz[ind+3] = row3[0];
        } else {
            break;
        }
        //
        ind += 4;
    }

    printf("Finished Reading in from file\n");
    
    //x = xyz;
    //y = &xyz[npoin];
    //z = &xyz[2*npoin];

    printf("Writing to File: %s\n", argv[2]);
    FILE* fout = fopen(argv[2],"w");
    fprintf(fout, "%d %d\n",ni,nj);
    
    double xx[ni][nj], yy[ni][nj];

    for (int i=0; i<ni; i++) {
    for (int j=0; j<nj; j++) {
        //int ind = ni*j + i; ///////////////////////////////////////////////////////////////////////////////////////
        int ind = nj*i + j;
        //fprintf(fout,"%.15le %.15le\n", xyz[ind], xyz[ind+npoin]);
        xx[i][j] = xyz[ind];
        yy[i][j] = xyz[ind+npoin];
    }
    }
    for (int j=0; j<nj; j++) {
    for (int i=0; i<ni; i++) {
        fprintf(fout,"%.15le %.15le\n", xx[i][j], yy[i][j]);
    }
    }
    
    // Boundary Faces
    fprintf(fout,"BC\n");

    int nbfac = 2*(ni-1) + 2*(nj-1);
    int ipjm = ni-1;        //Last index of bottom face
    int ipjp = ipjm + nj-1; // " " right face
    int imjp = ipjp + ni-1; // top
    int imjm = imjp + nj-1; // left

    for (int ib=0; ib<nbfac; ib++) {
        int bid = -1;
        
        if (ib < ipjm){
            // Bottom Surface
            bid = 4;
        } else if (ib < ipjp){
            // Right Boundary
            bid = 3;
        } else if (ib < imjp){
            // Top Boundary
            bid = 4;
        } else if (ib < imjm){
            // Left Boundary
            bid = 1;
        } else {
            printf("Whoops!\n");
            exit(0);
        }


        fprintf(fout,"%d\n",bid);
    }
        
    fclose(fin);
    fclose(fout);
	return 0;
}






