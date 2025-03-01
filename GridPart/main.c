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
    int nblock, nbtot;
    //
    nbtot = 0;
    //
    // Read input file for grid metadata
    printf("========== READING GCONN ==========\n");
    FILE* fgconnin;
    //
    fgconnin = fopen("./gridin.conn","r");
    fscanf(fgconnin,"%d",&nblock);
    printf("nblock = %d\n",nblock);
    //
    int ixm[nblock], ixp[nblock], iym[nblock], iyp[nblock];
    int ipart[nblock], jpart[nblock];
    for (int ib=0; ib<nblock; ib++){
        // Boundary info 
        // Neighbor ib >= 0
        // External boundary ib = boundary
        fscanf(fgconnin,"%d %d %d %d", &ixm[ib], &ixp[ib], &iym[ib], &iyp[ib]); 
        fscanf(fgconnin,"%d %d", &ipart[ib], &jpart[ib]); 
        nbtot += ipart[ib]*jpart[ib];
    }
    fclose(fgconnin);

    printf("========== READING GRIDS ==========\n");
    FILE* fgrid;
    fgrid = fopen("./gridin.dat","r");

    int ni[nblock], nj[nblock];
    double **x[nbtot], **y[nbtot];
    //int npoin = ni * nj;
    //int nb = 2*(ni-1) + 2*(nj-1);
    //int ibound[nb];

    for (int b=0; b<nblock; b++){ 
    printf("========== GRID: %d\n",b);
    fscanf(fgrid, "%d %d", &ni[b], &nj[b]);
    printf("(ni,nj) : (%4d,%4d)\n",ni[b],nj[b]);

    x[b] = (double**)malloc(nj[b]*sizeof(double*));
    y[b] = (double**)malloc(nj[b]*sizeof(double*));
        for (int j=0; j<nj[b]; j++){
            x[b][j] = (double*)malloc(ni[b]*sizeof(double));
            y[b][j] = (double*)malloc(ni[b]*sizeof(double));
        }
        // read block in
        for (int j=0; j<nj[b]; j++){
        for (int i=0; i<ni[b]; i++){
            fscanf(fgrid, "%le %le %*e", &x[b][j][i], &y[b][j][i]);
        }
        }
    }
    fclose(fgrid);
    
    printf("========== Decomposing Grids ==========\n");
    //Create arrays containing the boundary block numbers of each block;
    int mxpart = 0, iblk = -1;
    for (int b=0;b<nblock;b++){mxpart = MAX(mxpart, ipart[b]); mxpart = MAX(mxpart, jpart[b]);}
    //
    int blkind[nblock][mxpart][mxpart];       
    for (int b=0;b<nblock;b++){
        for (int j=0;j<jpart[b];j++){
        for (int i=0;i<ipart[b];i++){
            iblk++;
            blkind[b][i][j] = iblk;
        }    
        }    
    }
    printf("========== GRIDS INDEXED\n");

    // ~~~~~~~~~~ Transform from input blocks to total blocks
    //Check that partitions and grid sized on on neighboring blocks match up
    FILE* fgconn;
    fgconn = fopen("./grid.conn","w");
    fprintf(fgconn,"%d\n",nbtot);
    //
    for (int b=0; b<nblock; b++){
        printf("========== PARTITIONING GRID %d \n",b);
        //Check for partition or grid errors in the input grid
        if(ixm[b] >= 0){
            if(nj[b]    != nj[ixm[b]])     {printf("Grid(%d) Number Mismatch XM\t%d~%d\n",b, nj[b],nj[ixm[b]]);exit(100);}
            if(jpart[b] != jpart[ixm[b]])  {printf("Grid(%d) Partition Mismatch XM\n",b);exit(100);}
        }
        if(ixp[b] >= 0){
            if(nj[b]    != nj[ixp[b]])     {printf("Grid(%d) Number Mismatch XP\t%d~%d\n",b, nj[b],nj[ixp[b]]);exit(100);}
            if(jpart[b] != jpart[ixp[b]])  {printf("Grid(%d) Partition Mismatch XP\n",b);exit(100);}
        }
        if(iym[b] >= 0){
            if(ni[b]    != ni[iym[b]])     {printf("Grid(%d) Number Mismatch YM\t%d~%d\n",b, ni[b],ni[iym[b]]);exit(100);}
            if(ipart[b] != ipart[iym[b]])  {printf("Grid(%d) Partition Mismatch YM\n",b);exit(100);}
        }
        if(iyp[b] >= 0){
            if(ni[b]    != ni[iyp[b]])     {printf("Grid(%d) Number Mismatch YP\t%d~%d\n",b, ni[b],ni[iyp[b]]);exit(100);}
            if(ipart[b] != ipart[iyp[b]])  {printf("Grid(%d) Partition Mismatch YP\n",b);exit(100);}
        }
        //Decompose current block
        int istep, jstep;
        istep = ceil((double)ni[b]/ipart[b] - 1e-8);
        jstep = ceil((double)nj[b]/jpart[b] - 1e-8);
        printf("B(%d) istep, jstep : %d, %d\n", b, istep, jstep);    
        
        // Loop through subblock
        for (int jblock=0; jblock<jpart[b]; jblock++) {
        for (int iblock=0; iblock<ipart[b]; iblock++) {
            printf("========== WRITING SUBGRID %d:(%2d, %2d) \n",b,iblock,jblock);
            // Get boundary of partition
            int istart = istep*iblock;
            int jstart = jstep*jblock;
            int iend = MIN(ni[b]-1, istep*(iblock+1));
            int jend = MIN(nj[b]-1, jstep*(jblock+1));
            //Add on padding
            //istart = MAX(istart-1, 0); 
            //jstart = MAX(jstart-1, 0);
            int block_num = blkind[b][jblock][iblock];
            //  
            printf("========== WRITING BLOCK %d ==========\n", block_num );
            //
            char fsgname[50];
            sprintf(fsgname,"./subgrids/grid.block%d.dat",block_num);
            printf("%s\n",fsgname);
            //
            FILE* fsubgrid = fopen(fsgname, "w");
            if (fsubgrid == NULL) {printf("Failed to open file: %s.\n",fsgname); exit(0);}
            //
            fprintf(fsubgrid, "npartb: %d %d \n", ipart[b], jpart[b]);
            fprintf(fsubgrid, "ipartb: %d %d \n", iblock,jblock);
            fprintf(fsubgrid, "ngrid: %d %d \n", iend-istart+1, jend-jstart+1);
            for (int j=jstart; j<=jend; j++){
            for (int i=istart; i<=iend; i++){
                fprintf(fsubgrid, "%le %le \n", x[b][j][i], y[b][j][i]);
            }
            }
            //
            // Print out Boundary information for the cell and grid connectivity
            fprintf(fsubgrid,"BC\n");
            //
            int ileft,iright,idown,iup;
            int bleft,bright,bdown,bup;
            //i<dir> = boundary type
            //ileft  = iblock == 0        ? 0 : -1;
            //iright = iblock == ipart-1  ? 1 : -1;
            //idown  = jblock == 0        ? 2 : -1;
            //iup    = jblock == jpart-1  ? 3 : -1;
            //b<dir> = block number
            //bleft  = ileft  >= 0 ? -1 : block_num-1;
            //bright = iright >= 0 ? -1 : block_num+1;
            //bdown  = idown  >= 0 ? -1 : block_num-ipart;
            //bup    = iup    >= 0 ? -1 : block_num+ipart;
            //Corrections to recognize original block structure
            //LEFT
            if (ixm[b]>= 0){
                int ibn = ixm[b];
                int iind = ipart[ibn] - 1;
                bleft = blkind[ibn][jblock][iind];
                ileft = -1; 
            } else {
                bleft = -1;
                ileft = -ixm[b];
            }
            //RIGHT
            if (ixp[b]>= 0){
                int ibn = ixp[b];
                int iind = 0;
                bright = blkind[ibn][jblock][iind]; 
                iright = -1; 
            } else {
                bright = -1;
                iright = -ixp[b];
            }
            //BOT
            if (iym[b]>= 0){
                int ibn = iym[b];
                int jind = jpart[ibn] - 1;
                bdown = blkind[ibn][jind][iblock];
                idown = -1; 
            } else {
                bdown = -1;
                idown = -iym[b];
            }
            //TOP
            if (iyp[b]>= 0){
                int ibn = iyp[b];
                int jind = 0;
                bup = blkind[ibn][jind][iblock]; 
                iup = -1; 
            } else {
                bup = -1;
                iup = -iyp[b];
            }
            //Loop over bot
            for (int i=istart; i<iend; i++){
                if (idown==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
                   //int ibc = i;
                   fprintf(fsubgrid,"%d\n",idown); 
                }
            }
            //Loop over right
            for (int j=jstart; j<jend; j++){
                if (iright==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
                   //int ibc = j + ni - 1;
                   fprintf(fsubgrid,"%d\n",iright); 
                }
            }

            int ibend = (ni[b]-1) - istart;
            int ibstart = (ni[b]-1) - iend + 1;
            //Loop over top
            for (int i=ibstart-1; i<ibend; i++){
                  if (iup==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
                   //int ibc = ni + nj - 2 + i;
                   fprintf(fsubgrid,"%d\n",iup); 
                 }
            }
            int jbend = (nj[b]-1) - jstart;
            int jbstart = (nj[b]-1) - jend + 1;
            //Loop over left
            for (int j=jbstart-1; j<jbend; j++){
//                if (iblock == 0 && jblock == 0) {printf("j = %d\n",j);}

                if (ileft==-1) {fprintf(fsubgrid,"%d\n", -1);} else {
                    //int ibc = nj + 2*ni - 3 + j;
                    fprintf(fsubgrid,"%d\n",ileft); 
                }
            }
           
            //print to grid conn file 
            fprintf(fgconn,"%d: (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n", 
                           block_num, ileft, bleft, iright,bright, idown ,bdown, iup, bup); 
           
            fclose(fsubgrid);
       
        }
        }
    }
    fclose(fgconn);
    free(x);
    free(y);

}
/*
int block_decomp(int ipart, int jpart) {
    //Partition the given grid into an input number of partitions with given padding
    int npart = ipart * jpart;
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
*/
