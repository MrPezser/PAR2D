#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int main(int argc, char** argv) {
    // read in a gridpro style multiblock file and partition into a number of smaller domains
    // gridin.dat = source file
    // gridin.conn = number, connectivity, and BC of original block file
    //
    // grid.conn = output block connectivity and BC
    // subgrids/grid.block###.dat = output grid blocks
    int nblock;

    // Read input file for grid metadata
    FILE* fgconn;
    fgconn = fopen("./gridinconn.dat","r");

    printf("========== READING GRID ==========\n");
    FILE* fgrid;
    fgrid = fopen("./grid_in.dat","r");

    int ni, nj;
    int npoin = ni * nj;
    int nb = 2*(ni-1) + 2*(nj-1);
    double x[nj][ni][nblk], y[nj][ni][nblk];
    int ibound[nb];

    for (int b=0; b<nblock; b++){ 
    fscanf(fgrid, "%d %d", &ni, &nj);
    for (int j=0; j<nj; j++){
    for (int i=0; i<ni; i++){
        fscanf(fgrid, "%le %le", &x[j][i], &y[j][i]);
    }
    }
    }
    
}

int block_decomp(int npart, int ipart, int jpart) {
    //Partition the given grid into an input number of partitions with given padding


    //read BC assignments
    fscanf(fgrid,"%*s");
    for(int ib=0; ib<nb; ib++){
        fscanf(fgrid, "%d", &(ibound[ib]));
    }

    fclose(fgrid);


   int istep, jstep;
   istep = ceil((double)ni/ipart - 1e-8);
   jstep = ceil((double)nj/jpart - 1e-8);
   printf("istep, jstep : %d, %d\n", istep, jstep);    

   int block_num = 0;
   for (int jblock=0; jblock<jpart; jblock++) {
   for (int iblock=0; iblock<ipart; iblock++) {
        // Get boundary of partition
        int istart = istep*iblock;
        int jstart = jstep*jblock;
        int iend = MIN(ni-1, istep*(iblock+1));
        int jend = MIN(nj-1, jstep*(jblock+1));

        //Add on padding
        //istart = MAX(istart-1, 0); 
        //jstart = MAX(jstart-1, 0);
        block_num++;

        printf("========== WRITING BLOCK %d ==========\n", block_num );
        
        char fsgname[50];
        sprintf(fsgname,"./subgrids/grid.block%d.dat",block_num);
        printf("%s\n",fsgname);
        
        FILE* fsubgrid = fopen(fsgname, "w");
        if (fsubgrid == NULL) {printf("Failed to open file.\n"); exit(0);}

        fprintf(fgrid, "npart: %d %d \n", ipart, jpart);
        fprintf(fgrid, "ipart: %d %d \n", iblock,jblock);
        fprintf(fgrid, "ngrid: %d %d \n", iend-istart+1, jend-jstart+1);
        for (int j=jstart; j<=jend; j++){
        for (int i=istart; i<=iend; i++){
            fprintf(fgrid, "%le %le \n", x[j][i], y[j][i]);
        }
        }
        // Print out Boundary information for the cell and grid connectivity
        fprintf(fsubgrid,"BC\n");

        int ileft,iright,idown,iup;
        int bleft,bright,bdown,bup;
        ileft  = iblock == 0        ? 0 : -1;
        iright = iblock == ipart-1  ? 1 : -1;
        idown  = jblock == 0        ? 2 : -1;
        iup    = jblock == jpart-1  ? 3 : -1;
        //i<dir> = boundary type

        bleft  = ileft  >= 0 ? -1 : block_num-1;
        bright = iright >= 0 ? -1 : block_num+1;
        bdown  = idown  >= 0 ? -1 : block_num-ipart;
        bup    = iup    >= 0 ? -1 : block_num+ipart;
        //b<dir> = block number

        //Loop over bot
        for (int i=istart; i<iend; i++){
            if (idown==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
               int ibc = i;
               fprintf(fsubgrid,"%d\n",ibound[ibc]); 
            }
        }
        //Loop over right
        for (int j=jstart; j<jend; j++){
            if (iright==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
               int ibc = j + ni - 1;
               fprintf(fsubgrid,"%d\n",ibound[ibc]); 
            }
        }

        int ibend = (ni-1) - istart;
        int ibstart = (ni-1) - iend + 1;
        //Loop over top
        for (int i=ibstart-1; i<ibend; i++){
            if (iup==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
               int ibc = ni + nj - 2 + i;
               fprintf(fsubgrid,"%d\n",ibound[ibc]); 
            }
        }
        int jbend = (nj-1) - jstart;
        int jbstart = (nj-1) - jend + 1;
        //Loop over left
        for (int j=jbstart-1; j<jbend; j++){
            if (iblock == 0 && jblock == 0) {printf("j = %d\n",j);}

            if (ileft==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
               int ibc = nj + 2*ni - 3 + j;
               fprintf(fsubgrid,"%d\n",ibound[ibc]); 
            }
        }
        
        //print to grid conn file 
        fprintf(fgconn,"%d: (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n", 
                        block_num, ileft, bleft, iright,bright, idown ,bdown, iup, bup); 
        
        fclose(fsubgrid);
   }
   }
    
   //idk why this is here
   // for (int iblk=0; iblk<npart; iblk++){
   // }
}
