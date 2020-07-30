# Parallel Computing with MPI

Parallel implementation of Gauss-Seidel algorithm used to solve linear system of equations using **MPI (Message Passing Interface)**. Aditional information about the method on *“Parallel Computer Architecture. A Hardware/Software Approach. David Culler, Jaswinder Pal Singh, Anoop Gupta”*.

# Description

The program automatically creates a nxn matrix and initialices it with random numbers then solves this matrix using a parallel implementaion of Gauss-Seidel algorithm.

```Note```: some comments and parts of the code are written in Spanish.

# Usage

**Inputs** 

- ```n```: positive integer number representing the matrix size to be create (different of the ```-n``` input of ```mpiexec``` command)
- ```communication```: either 0 or 1 for point-to-point communication or collective communication

**Example**

To create a 1026x1026 matrix and solve it using collective communication with 4 processes the command will look like this:

```mpiexec –n 4 –f machine_file ./gs_mpi 1026 1```

**Output**

At the end of the execution some interesting information is printed. The information covers the following data:
- ```Matrix size```: matrix size                 
- ```Communication```: communication type         
- ```# procesos```: total number of processes                      
- ```TIEMPO SOLUCION```: time to converge                 
- ```tiempo mem. and MPI```: time loading matrix and MPI stuffs              
- ```tiempo scatter matriz```: time consumed for scattering matrix            
- ```tiempo gather  matriz```: time consumed for gathering matrix           
- ```tiempo en solver```: time spent in solving                  
- ```tiempo total comu. diff y done```: time consumed to communicate some variables   
- ```tiempo total comu. filas vecinas```: time consumed communicating rows between processes 
- ```iteraciones```: total number of iterations                    
- ```t. media por itera.```: mean time per iteration              
- ```t. media comu. diff y done```: mean time communicating some varibles       
- ```t. media comu. filas vecinas```: mean time spent for communicating rows between processes  

# End notes

This code was part of my studies and I hope it can be useful for you in your studies on Parallel Computing. 
This piece of code can be used for testing performance on different computers and other some interesting stuff you can imagine.




