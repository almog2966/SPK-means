#ifndef KMEANS_H_
#define KMEANS_H_


void printWam(double** A,int n, int d);
double **wam(double **A,int n, int d);

double **ddg(double **A,int n , int d);
void printDdg(double **A, int n, int d);

void printJacobi(double **A,int n, int d);
double **jacobi(double **A,int n ,int d);

void printGl(double **A, int n,int d);
double **gl(double **A, int n,int d);

void freeMatrix(double **temp,int n);
void printMatrix(double** A, int n, int d);

double **spk(double **A,int n, int d,int *k);


#endif
