
#include <string.h>
#include <stdio.h>
#include "mpi.h"
#include "stdlib.h"
#include "math.h"
int flagFunction=0;
double Trapezoids(double,double,double,double);
double Simpson(double ,double  ,double );
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

    double total=0;
    
    double a,b,n;

    //MPI INIZITIALIZATION
    MPI_Init(&argc,&argv);              
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);    

    MPI_Barrier(MPI_COMM_WORLD);
    if(myRank == MASTER){
        start = MPI_Wtime();    
    }

    	if(argc == 6){
		n = atoi(argv[1]);
		a = atof(argv[2]);
		b = atof(argv[3]);
		flagFunction = atoi(argv[4]); /*Note: you should include some functionality to switch between different functions*/
    integrationMethod=atof(argv[5]);
    }else{
        n=1024;
        a=0;
        b=3;
        flagFunction=4;
        integrationMethod=0;
    }
    
    //same value for both algorithem
    h = (b-a)/n; 
   
	double local_n= n/nProc; //this is or number of trap or number of subintervals for each process

    if(myRank == MASTER){
        
        printf("n: %f a:%f b:%f method : ",n,a,b);

            if(integrationMethod==0){
                printf("Trapezoid Methods ");
            }
            else{
                  printf("SIMPSON Methods ");
            }
            
            switch (flagFunction)
            {
            case 0:
                printf(" ,Function: SIN\n");
                 break;
            case 1:
                  printf(", Function: COS\n");
            break;   
            case 2:
                    printf(", Function: TAN\n");
            break;
            case 3:
                   printf(", Function: 1/X\n");
            break;   
            case 4:
                  printf(" ,Function: X^2\n");
            
            break;
            case 5:
                  printf(" , Function:  1+X/(X+2.5)");
                  break;

             case 6:
                  printf(" , Function:  15 * cos(x)*sin(x) * sqrt(2*x) * 3.1459 * 7;");
                  break;

             case 7:
                  printf(" , Function:  pow( abs(x),(15 * cos(x)* sqrt(2*x) * 3.1459)) * sqrt(2*x);");
                  break;
            default:
                   printf(" ,Function: COS\n");
            break;
            }
            
           
                
           }

    double local_a;

  
    if( integrationMethod == 0){  //TRAPEZOIDS RULE
      //here set the start and end of integration interval
      
        double intervallo=local_n *h;
        local_a= a + myRank * intervallo;
        double local_b=local_a + intervallo;

        integral = Trapezoids(local_a,local_b,local_n,h);
        //here call the integration method 
     
    }else{ //SIMPSON RULE
        local_a = a+myRank*(b-a)/nProc;	
       integral = Simpson(local_a,local_n,h);
    }

    

    
    MPI_Reduce(&integral, &total, 1, MPI_DOUBLE, MPI_SUM, MASTER , MPI_COMM_WORLD);



    
   

  

       
   




    MPI_Barrier(MPI_COMM_WORLD);
    if(myRank == MASTER){
        end = MPI_Wtime();
        printf("\ntotal: %f \n",total); 
        printf("time final : %f \n",end-start);    
    }

    MPI_Finalize();




    return 0;

}
double chosenFunction(double x){
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
    case 4:
      return_val = x*x;
      
      break;
     case 5:
           return_val = (1+x)/(x+2.5);
      break;
        case 6:
           return_val =   15 * cos(x)*sin(x) * sqrt(2*x) * 3.1459 * 7 * 20 * tan(x);
      break;
    
    
    default:
      return_val = cos(x);
      break;
  }
	/*Add your functions here, should be able to switch functions depending on the user's flag passed in*/
	return return_val;
}


double Simpson(double local_a, double local_n, double h){
	/*Write your code here to do the integration*/
	double result = 0;
	for(int i =0; i<=local_n; i++)
	{
		// if i is 0 or n, only add f(x_i)
		if((i==0)||(i==local_n))
			result += chosenFunction(local_a + i* h);
		else if(i%2==0)
			// then if i is even, add 2*f(x_i)
			result += 2* chosenFunction(local_a + i*h);
		else
			// else if i is odd, add 4*f(x_i)
			result += 4* chosenFunction(local_a + i*h);	
	}
	result *= h/3;  // multiply the above series by h/3
	return(result);
}

double Trapezoids(double local_a,double local_b,double local_n,double h ){
    double integral;

    double x;
    int i;
    
    integral = (chosenFunction(local_a) + chosenFunction(local_b))/2.0;
    x =local_a;

    for(i=1;i<=local_n-1;i++){
        x=x+h;
        integral=integral + chosenFunction(x);
    }

    integral=integral *h ;
   
    return integral;
}

