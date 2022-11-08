
#include <string.h>
#include <stdio.h>
#include "mpi.h"

int main(int argc,char* argv[]){

    //variable declaration
    int myRank;//rank del processo
    int nProc;//num of process
    int MASTER=0;
    int messageTag=0;
    double h;  //base width of segmet or trapezoid
    double start,end; //time start,end
    int integrationMethod=0; // 0 -> trapezoid rule / 1-> simpson rule



    //MPI INIZITIALIZATION
    MPI_Init(&argc,&argv);              
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);    


    //same value for both algorithem
    h = (b-a)/n;  
	double local_n= n/nProc; //this is or number of trap or number of subintervals for each process
    
    if(myRank == MASTER)
        for(int i=1;i<argc;i++)
        {
            printf("%s\n",argv[i]);
            
            
            if(strcmp(argv[i], "-s") == 0) 
                integrationMethod=1;
            else
                integrationMethod=0;
        
        }



  
    if( integrationMethod == 0){  //TRAPEZOIDS RULE
      //here set the start and end of integration interval
      
        double intervallo=local_n *h;
        double local_a= a + myRank * intervallo;
        double local_b=local_a + intervallo;

        Trapezoids(local_a,local_b,local_n,h)
        //here call the integration method 

    }else{ //SIMPSON RULE
        local_a = a+my_rank*(b-a)/p;	
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