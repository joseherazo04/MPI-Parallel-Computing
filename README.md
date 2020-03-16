# Parallel-Computing

Parallel implementation of Gauss-Seidel algorithm used to solve linear system of equations using MPI. Aditional information about the method on “Parallel Computer Architecture. A Hardware/Software Approach. David Culler, Jaswinder Pal Singh, Anoop Gupta”.

# Description

The program automatically creates a nxn matrix and initialices it with random numbers then solves this matrix using a parallel implementaion of Gauss-Seidel algorithm.

# Usage

**Inputs** 

- ```n```: positive integer number representing the matrix size to be create (different of the ```-n``` input of ```mpiexec``` command)
- ```communication```: either 0 or 1 for point-to-point communication or collective communication

**Example**

To create a 1026x1026 matrix and solve it using collective communication with 4 processes the command will look like this:

```mpiexec –n 4 –f machine_file ./gs_mpi 1026 1```




