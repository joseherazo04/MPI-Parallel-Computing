/*
PARALELIZACION DEL ALGORITMO DE GAUSS SEIDEL
*/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define MAX_ITER 100   //maximun number of iteration
#define MAX 100        //maximum value of the matrix element
#define TOL 0.000001   //tolerance

int np, myrank, communication,cnt_iter = 0;        //number of processes and rank
float *vector;

double tstartdiffdone, tfinishdiffdone, TotalTimediffdone;
double tstartrow, tfinishrow, TotalTimerow;


//**********************GENERATES RAMDOMS*************************//

float rand_float(int max){        // Generate a random float number with the maximum value of max
  return ((float)rand()/(float)(RAND_MAX)) * max;
}

//!!!---NUEVA VERSION INITIALIZES MATRIX FOR VECTOR
void allocate_init_vector(int n, int m){
  int i, j;

  vector = (float *) malloc(n * m * sizeof(float));

  for(i = 0; i < n; i++) {

    for (j = 0; j < m; j++)
      vector[i*m+j] = rand_float(MAX);
  }
}


//******************RESERVA ESPACIO EN MEM***********************//

void allocate_mem(int n, int m){

  vector = (float *) malloc(n * m * sizeof(float));

}


//******************   IMPRIME  VECTOR   ***********************//
void printVector(int n, int m){
    
    int i,j;

    for (i = 0; i < n; i++) {
      printf("\n");
      for (j = 0; j < m; j++){
	printf(" %.2f",vector[i*m+j]);
      }
    }
    printf("\n");
}

//***********************SOLVE MATRIX****************************//

void solver(int n, int m){
  float diff = 0, temp, difftemp=0;
  int done = 0,  i, j;

  while (!done && (cnt_iter < MAX_ITER)){


    diff = 0;
    
    //------ACTUALIZAR ELEMENTOS----
    for (i = 1; i < n - 1; i++) {

      for (j = 1; j < m - 1; j++) {

	temp = vector[i*m+j];
	vector[i*m+j] = 0.2 * (vector[i*m + j] + vector[i*m + j-1] + vector[(i-1)*m + j] + vector[i*m + j+1] + vector[(i+1)*m + j]);
	diff += abs(vector[i*m + j] - temp); //Calcula mi diff

      }
    }

    //----COMUNICACIÓN DE LAS FILAS VECINAS--

    tstartrow= MPI_Wtime(); //Medicion del tiempo de comunicacion de filas vecinas

    if(myrank%2==0) { //Ranks pares envían y luego reciben

      if(myrank!=np-1){

        MPI_Send(&vector[(n-2)*m],m,MPI_FLOAT,myrank+1,1,MPI_COMM_WORLD); //envía penúltima
        MPI_Recv(&vector[(n-1)*m],m,MPI_FLOAT,myrank+1,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);  //recibe en última
      }

      if(myrank!=0){

        MPI_Send(&vector[m],m,MPI_FLOAT,myrank-1,1,MPI_COMM_WORLD); //envía segunda
        MPI_Recv(&vector[0],m,MPI_FLOAT,myrank-1,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);  //recibe en primera
      }

    }

    if(myrank%2==1) { //Ranks impares reciben luego envían

      MPI_Recv(&vector[0],m,MPI_FLOAT,myrank-1,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);  //recibe en primera
      MPI_Send(&vector[m],m,MPI_FLOAT,myrank-1,1,MPI_COMM_WORLD); //envía segunda

      if(myrank!=np-1){

        MPI_Recv(&vector[(n-1)*m],m,MPI_FLOAT,myrank+1,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);  //recibe en última
        MPI_Send(&vector[(n-2)*m],m,MPI_FLOAT,myrank+1,1,MPI_COMM_WORLD); //envía penúltima
      }

    }

    tfinishrow = MPI_Wtime();
    TotalTimerow = TotalTimerow + tfinishrow - tstartrow;
    
    
    //-----COMUNICACION DIFF Y DONE----

    tstartdiffdone= MPI_Wtime(); //Medicion del tiempo de comunicacion de diff y done

    switch(communication){
      case 0: { //LLAMADAS P2P

        //----ENVIAR DIFF y Esperar DONE de 0----
        if(myrank!=0) {

          MPI_Send(&diff,1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
          MPI_Recv(&done, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       }

       //----SI SOY 0 RECIBO DE TODOS SU DIFF y digo si done
       if(myrank==0) {
         int l; 

         for(l=1;l<np;l++) {
           MPI_Recv(&difftemp, 1, MPI_FLOAT, l, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
           diff=diff+difftemp;
         }

         //Calcular si terminamos
         if (diff/n/n < TOL)  done=1;

         //Enviar done a todos
         for(l=1;l<np;l++) 
           MPI_Send(&done,1,MPI_INT,l,1,MPI_COMM_WORLD); 
       
       }
     
        break;
      }
   
      case 1: { //LLAMADAS COLECTIVAS
        MPI_Reduce(&diff, &difftemp, 1, MPI_FLOAT, MPI_SUM, 0,MPI_COMM_WORLD);

        //Calcular si terminamos
        if (difftemp/n/n < TOL)  done=1;

        MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
        break;
      }

    }

    tfinishdiffdone = MPI_Wtime();
    TotalTimediffdone = TotalTimediffdone + tfinishdiffdone - tstartdiffdone;
   
    //CUENTA ITERACIÓN    
    cnt_iter ++;   

  }
}


//*********************** M A I N ****************************//

int main(int argc, char *argv[]) {
  
  int n,tam;
  float **a;
  double tstart, tfinish, TotalTime;
  double tstartscatter, tfinishscatter, TotalTimescatter;
  double tstartgather, tfinishgather, TotalTimegather;
  double tstartsolver, tfinishsolver, TotalTimesolver;
  double tfinishMPI, TotalTimeMPI; //utilizado para medir el tiempo que se tarda en cargar la matriz y alojar memoria

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //id de proceso
  MPI_Comm_size(MPI_COMM_WORLD, &np);     //num procs para el comunicador COMM_WORLD

  tstart = MPI_Wtime();

  //-----CHECKS IF EVERYTHING IT'S OK!-----
  if (argc < 3 || np < 2) {

    if (myrank == 0) {
      printf("\nRemember args or create more than ONE process please!\n\n");
    }

    MPI_Finalize();
    exit(1);
  }

  //-----RECEIVES PROGRAM ARGS!-----
  n = atoi(argv[1]);              //num de elementos en el arreglo
  communication =  atoi(argv[2]); //tipo de comunication


  //-----PROCESO 0 CREA LA MATRIZ -------

  if(myrank==0)//inicializa matriz aplanada
    allocate_init_vector(n,n);


  //-----PROCESSES EXCEPT 0 RESERVE MEM--

  tam=n/np+2; //filas que tendrá cada proceso 
              //every process will have its own matrix or block

  if(myrank!=0&&myrank!=np-1)
    allocate_mem(tam,n);//here everyone has its space for receiving data
                            
  if(myrank==np-1)
    allocate_mem(tam-1,n);//de acuerdo con el análisis el último proceso no necesita tantas filas (le vale tam-1) 

  tfinishMPI = MPI_Wtime();
  TotalTimeMPI = tfinishMPI - tstart;

  //------------SCATTER SWITCH-----------

  tstartscatter= MPI_Wtime(); //Medicion del tiempo de scatter

  switch (communication){
    case 0: { // P2P communication for scattering
      
      //PROCESS 0 SENDS MATRIX P2P
      if(myrank==0) {

        int e=tam-2, proc=1, i; //total de elementos a enviar por proceso (conteo de filas) y num de proceso

        for(i=e;i<n;i+=e) { //PROCESO DE ENVIO MUY EXTRAÑO PERO PERMITE QUE EL PROGRAMA SE AJUSTE AL NÚM DE NODOS ELEGIDOS DEPENDIENDO DEL TAM DE LA MATRIZ
        
          if(proc!=np-1)
            MPI_Send(&vector[i*n-n],n*tam,MPI_FLOAT,proc,1,MPI_COMM_WORLD);  
          else          
            MPI_Send(&vector[i*n-n],n*tam-n,MPI_FLOAT,proc,1,MPI_COMM_WORLD);
         
          proc++; //siguiente proceso	
        }
      }
      
      if(myrank!=0&&myrank!=np-1) //LOS DEMÁS MENOS EL ÚLTIMO RECIBEN
          MPI_Recv(&vector[0], n*tam, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
      if(myrank==np-1) //ÚLTIMO RECIBE MENOS
          MPI_Recv(&vector[0], n*tam-n, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
          

    break;
    }

    case 1: { // Collective communication for scattering

      int sendcount[np],displ[np];
      int e=tam-2, proc,i;

      //proceso 0 no recibe 
      sendcount[0]=0; displ[0]=0;

      for(i=e,proc=1;i<n;i+=e,proc++){ //Carga los desplazamientos
          displ[proc]=(i-1)*n;          
      }

      for(i=1;i<np-1;i++) //carga todos los sizes, una fila menos para el last one
        sendcount[i]=n*tam;
      sendcount[np-1]=n*tam-n;

      MPI_Scatterv(vector,sendcount,displ,MPI_FLOAT,vector,tam*n,MPI_FLOAT,0,MPI_COMM_WORLD);

    break;
    }
  }

  tfinishscatter = MPI_Wtime();
  TotalTimescatter = tfinishscatter - tstartscatter;

 
  //------------SOLVER---------------

  tstartsolver= MPI_Wtime(); //Medicion del tiempo de scatter

  if(myrank==0) solver(tam-1, n);
  if(myrank!=0&&myrank!=np-1) solver(tam, n);
  if(myrank==np-1) solver(tam-1,n);

  tfinishsolver = MPI_Wtime();
  TotalTimesolver = tfinishsolver - tstartsolver;

  
  //-----GATHERING THE RESULTS-------

  tstartgather= MPI_Wtime(); //Medicion del tiempo de gather

  switch (communication){
    case 0: { // P2P communication for gathering

       //-----TODOS ENVÍAN SUS PEDAZOS A 0
      if(myrank!=0&&myrank!=np-1) //LOS DEMÁS MENOS EL ÚLTIMO ENVÍAN
         MPI_Send(&vector[n], n*(tam-2), MPI_FLOAT, 0, 1, MPI_COMM_WORLD); 
              
      if(myrank==np-1)
         MPI_Send(&vector[n], n*(tam-3), MPI_FLOAT, 0, 1, MPI_COMM_WORLD); 

     //PROCESS 0 RECEIVES MATRIX P2P
      if(myrank==0) {
       
        int e=tam-2, proc=1, i; //total de elementos a recibir por proceso (conteo de filas) y num de proceso

        for(i=e;i<n;i+=e) {//PROCESO DE ENVIO MUY EXTRAÑO PERO PERMITE QUE EL PROGRAMA SE AJUSTE AL NÚM DE NODOS ELEGIDOS DEPENDIENDO DEL TAM DE LA MATRIZ
 
          if(proc!=np-1)
            MPI_Recv(&vector[i*n],n*(tam-2),MPI_FLOAT,proc,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
          else          
            MPI_Recv(&vector[i*n],n*(tam-3),MPI_FLOAT,proc,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          proc++; //siguiente proceso	
        }
      }
    
    break;
    }

    case 1: { // Collective communication for gathering

      int recvcount[np],displ[np];
      int e=tam-2, proc,i;

      //proceso 0 no recibe 
      displ[0]=0;

      for(i=e,proc=1;i<n;i+=e,proc++){ //Carga los desplazamientos
          displ[proc]=i*n;          
      }

      for(i=1;i<np;i++) //carga todos los sizes, una fila menos para el last one
        recvcount[i]=n*(tam-2);
      recvcount[0]=0;

      MPI_Gatherv(vector+n,n*(tam-2),MPI_FLOAT,vector,recvcount,displ,MPI_FLOAT,0,MPI_COMM_WORLD);

    
    break;
    }
  }

  tfinishgather = MPI_Wtime();
  TotalTimegather = tfinishgather - tstartgather;

  //-----------LET'S FINISH----------
  tfinish = MPI_Wtime();
  TotalTime = tfinish - tstart;

  if(myrank==0) {

   printf("\n\n * * * * * * * * * * * REGISTRO * * * * * * * * * * * * * \n");
   printf("\n - - - - -Matrix size                      : %d", n); 
   printf("\n - - - - -Communication                    : %d", communication);
   printf("\n - - - - -# procesos                       : %d", np);
   printf("\n - - - - -TIEMPO SOLUCION                  : %lf s", TotalTime);
   printf("\n - - - - -tiempo mem. and MPI              : %lf s", TotalTimeMPI);
   printf("\n - - - - -tiempo scatter matriz            : %lf s", TotalTimescatter);
   printf("\n - - - - -tiempo gather  matriz            : %lf s", TotalTimegather);
   printf("\n - - - - -tiempo en solver                 : %lf s", TotalTimesolver);
   printf("\n - - - - -tiempo total comu. diff y done   : %lf s", TotalTimediffdone);
   printf("\n - - - - -tiempo total comu. filas vecinas : %lf s", TotalTimerow);
   printf("\n - - - - -iteraciones                      : %i ", cnt_iter);
   printf("\n - - - - -t. media por itera.              : %lf s", TotalTimesolver/cnt_iter);
   printf("\n - - - - -t. media comu. diff y done       : %lf s", TotalTimediffdone/cnt_iter);
   printf("\n - - - - -t. media comu. filas vecinas     : %lf s", TotalTimerow/cnt_iter);
   printf("\n\n * * * * * * * * * *  FIN REGISTRO * * * * * * * * * * * * \n");
   
   }

  MPI_Finalize();
  return 0;
}
