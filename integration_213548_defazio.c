
#include <string.h>
#include <stdio.h>
#include "mpi.h"
int flagFunction=0;
double Trapezoids(double,double,double,double);
int main(int argc,char* argv[]){

    //variable declaration
    int myRank;//rank del processo
    int nProc;//num of process
    int MASTER=0;
    int messageTag=0;
    double h;  //base width of segmet or trapezoid
    double start,end; //time start,end
    int integrationMethod=0; // 0 -> trapezoid rule / 1-> simpson rule
    double integral;
    


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

        integral = Trapezoids(local_a,local_b,local_n,h)
        //here call the integration method 
        printf(integral);
    }else{ //SIMPSON RULE
        local_a = a+my_rank*(b-a)/p;	
    }

    

    




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
double Trapezoids(double local_a,double local_b,double local_n,double h ){
    double integral;
    double x;
    int i;

    integral = (function(local_a) + function(local_b))/2.0;
    x =local_a;
    for(i=1;i<=local_n-1;i++){
        x=x+h;
        integral=integral +function(x);
    }

    integral=integral *h ;
    return integral;
}

double function(double x){
	double return_val;
	// simply choose a function based on flag
	// 0 => f(x)=cos(x)
	// 1 => f(x)=sin(x)
	// 2 => f(x)=tan(x)
	// 3 => f(x)=1/x
	switch (flagFunction) 
	{
    case 0:
      return_val = sin(x);
      break;
    case 1:
      return_val = cos(x);
      break;   
    case 2:
      return_val = tan(x);
      break;
    case 3:
      return_val = 1/x;
      break;   
    default:
      return_val = cos(x);
      break;
  }
	/*Add your functions here, should be able to switch functions depending on the user's flag passed in*/
	return return_val;
}