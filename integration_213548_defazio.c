
#include <string.h>
#include <stdio.h>
#include "mpi.h"

int main(int argc,char* argv[]){
     int myRank;//rank del processo
    int nProc;//num of process
    int MASTER=0;

    double start,end; //time start,end

    int integrationMethod=0;
    MPI_Init(&argc,&argv);              
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);    

    
     if(myRank == MASTER)
    for(int i=1;i<argc;i++)
    {
        printf("%s\n",argv[i]);
        
        
        if(strcmp(argv[i], "trap") == 0) 
            integrationMethod=0;
        else
            integrationMethod=1;
     
    }
   

    

    ///can chose integration method 
    // possible trapezoids or simpson method
    




      MPI_Barrier(MPI_COMM_WORLD);
      if(myRank == MASTER){
         start = MPI_Wtime();
        
    }
   

    printf("processo : %d \n",myRank);


   




        MPI_Barrier(MPI_COMM_WORLD);
     if(myRank == MASTER){
         end = MPI_Wtime();

        printf("time final : %f \n",end-start);    
      }

    MPI_Finalize();




    return 0;

}