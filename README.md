# Unnormalized Spectral Clustering Algorithm

This project implements a version of the unnormalized spectral clustering algorithm. 
The algorithm is implemented in Python and also utilizes a C module called spkmeans.c for efficient computation. The project provides an API through a Python file called spkmeans.py, which can be executed in the command line with the following arguments: K, goal, and file name.
K - number of iteration
goal - spk,wam,ddg,ggl,jacobi.
file name - datapoint in txt file

The project aslo provides an API through a C file called spkmeans.c, which can be executed in the command line with the following arguments: goal, and file name
goal - spk,wam,ddg,ggl.

goal specifies the task to be performed:

# project logics

spk: Perform a full spectral k-means algorithm on the data points in the given file. 
calculate as described below:
1. Form the weighted adjacency matrix
2. compute the graph laplacian L
3. compute the first K eigen vectors of L
4. let U be the matrix contain the eigen vectors of L as columns
5. treating each rwo of U as datapoint in R^k, cluster them using K - means algorithm.


wam: Calculate and output the weighted adjacency matrix for the data points in the given file.
ddg: Calculate and output the diagonal degree matrix for the data points in the given file.
gl: Calculate the graph Laplacian for the data points in the given file.
Jacobi: Calculate and output the eigenvalues and eigenvectors of the data points in the given file.

# working with the project
The Python API interacts with the C module spkmeansmodule.c, which provides efficient implementations of certain calculations. 
The C module can be used directly by calling spkmeans.c with the arguments: goal, file_name.
The C module can perform spk,wam,ddg,ggl.





# Project Structure

The project consists of the following files:

1. spkmeans.py: The Python file responsible for the API of the project. It can be executed in the command line to perform various tasks related to spectral clustering.
2. spkmeansmodule.c: The C API module that provides the interface between the Python code and the C implementation in the spkmeans.c.
3. spkmeans.h: The header file that defines all function prototypes used in spkmeansmodule.c and implemented in spkmeans.c.
4. spkmeans.c: The C module that contains the efficient implementations of the calculations needed for spectral clustering. can be used as API for directly working with the C algorithm.
5. kmeans.c: The c module that contain the efficient implementations for kmeans algorithm.
6. kmeansmodule.c: The Python C API module that provides the interface between the Python code and the C implementation in the kmeans.c.
7. setup.py: The build script used to create the necessary files that allow spkmeans.py to import mykmeansssp.
8. Makefile: A script for building the spkmeans executable.
