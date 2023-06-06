#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"



char *goal,*fileName;
char c;
FILE *ifp;
double **dataPoints;

/*matrix creation*/
double** createZeroMatrix(int n, int d);
double **createOneMatrix(int n);


/*matrix opeartion*/
void subsMatrix(double **target,double** A, double** B, int n, int d);
double** transposeMatrix(double** A, int n, int d);
double** mulMatrix(double **A, double **B,int n ,int m,int d);
double *matrixDiag(double **A,int n);
void printVector(double *vector,int d);

void copyMatrix(double** A, double ** target, int n, int d);
void copyVector(double *A, double *B, int n);

double **sortMatrix(double **A,double *vector, int n,int d);


/*jacobi auxiliary functions*/
int* calcPivotIndices(double **A,int n);
double* obtainCandS(double Aii,double Aij,double Ajj);
void pivotMoveJacobi(double **A, double **mulPmatrix, int n);
double calcOffConvergence(double **A,int n);
double **createMatrixP(int n, int i, int j, double c, double s);


/*numeric opretaions*/
int sign(double n);
double calcAbs(double num);




int main(int argc, char* argv[])
{
    int d,n,first,i,j;
    if(argc == 3)
    {
        goal=argv[1];
        fileName = argv[2];
    } 
    else
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    ifp = fopen(fileName, "r");
    if (ifp == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    d = 1;
    n = 0;
    first = 1;
    while ((c = fgetc(ifp)) != EOF)
    {
        if (c == ',' && first)
            d++;
        if(c == '\n')
        {
            n++;
            first = 0;
        }
    }
    rewind(ifp);
    dataPoints = (double **)malloc(n*sizeof(double*));
    if(dataPoints == NULL)
    {
        fclose(ifp);
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i = 0; i<n; i++)
    {
        dataPoints[i] = (double*)malloc(d*sizeof(double));
        if(dataPoints[i] == NULL)
            freeMatrix(dataPoints,n);
        for(j=0; j<d;j++)
        {
            if(fscanf(ifp,"%lf", &dataPoints[i][j])!=EOF)
                fgetc(ifp); /* skip ',' */
        }
    }
    fclose(ifp);
    if(!strcmp("wam",goal))
        printWam(dataPoints,n,d);
    else if(!strcmp("ddg",goal))
        printDdg(dataPoints,n,d);
    else if(!strcmp("gl",goal))
        printGl(dataPoints,n,d);
    else if(!strcmp("jacobi",goal))
        if(n==d)
            printJacobi(dataPoints,n,d);
        else/*if matrix is not symetric*/
        {
            printf("An Error Has Occurred\n");
            freeMatrix(dataPoints,n);
            exit(1);
        }   
    else{
        printf("An Error Has Occurred\n");
        freeMatrix(dataPoints,n);
        exit(1);
    }
    freeMatrix(dataPoints,n);
    return 0;
    
}

/* create new Matrix size n x d */
double** createZeroMatrix(int n, int d)
{
    int i,j;
    double** zeroMatrix = (double **)malloc(n * sizeof(double*));
    if(zeroMatrix == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i = 0; i<n; i++)
    {
        zeroMatrix[i] = (double*)malloc(d * sizeof(double));
        if(zeroMatrix[i] == NULL)
        {
            freeMatrix(zeroMatrix,n);
            printf("An Error Has Occurred\n");
            exit(1);
        }
            
        for(j=0; j<d;j++)
        {
            zeroMatrix[i][j] = 0;
        }
    }
    return zeroMatrix;
}

/*create I matrix size n x n*/
double **createOneMatrix(int n)
{
    int i,j;
    double** oneMat = (double **)malloc(n * sizeof(double*));
    if(oneMat == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i = 0; i<n; i++)
    {
        oneMat[i] = (double*)malloc(n * sizeof(double));
        if(oneMat[i] == NULL)
            freeMatrix(oneMat,n);
        for(j=0; j<n;j++)
        {
            if(i==j)
                oneMat[i][j] = 1;
            else
                oneMat[i][j] = 0;
        }
    }
    return oneMat;
}

/* res = A - B */
void subsMatrix(double **target,double **A, double **B, int n, int d)
{
    int i = 0;
    int j = 0;
    for(i=0;  i<n; i++){
        for (j=0; j< d; j++){
            target[i][j]= A[i][j]- B[i][j];
        }
    }
}

/*matrix multiplaction return new matrix: AxB*/
double** mulMatrix(double **A, double **B,int n ,int m, int d)
{
    int i,j,k;
    double sum;
    double **result = createZeroMatrix(n,d);
    for(i=0;  i<n; i++){
        for (j=0; j< d; j++){
            sum=0;
            for (k=0; k< m; k++){
                sum+= A[i][k] * B[k][j];
            }
            result[i][j]=sum;
        }
    }
    return result;
}

/* return A transpose return new Matrix */
double** transposeMatrix(double **A, int n, int d)
{
    double **transpose;
    int i = 0;
    int j = 0;
    transpose = createZeroMatrix(n, d);    
    for (i=0; i< n;i++){
        for (j=0; j<d; j++){
            transpose[i][j]= A[j][i];
        }
    }
    return transpose;
}

/*return array of diag elemnts*/
double *matrixDiag(double **A,int n)
{
    int i = 0;
    double *diagArray = (double*)malloc(sizeof(double)*n);
    if(diagArray == NULL)
    {
        printf("An Error Has Occurred\n");
        freeMatrix(A,n);
        exit(1);
    }
    for(i=0; i<n; i++)
        diagArray[i]=A[i][i];
    return diagArray; 
}



/*copy Matrix A to Matrix targer*/
void copyMatrix(double** A, double ** target, int n, int d)
{
    int i,j;
    for (i=0; i<n; i++){
        for(j=0; j<d; j++){
            target[i][j]= A[i][j];
        }
    }
}

/*copy A to B*/
void copyVector(double *A, double *B, int n)
{
    int i;
    for(i = 0; i < n; i++){
        B[i] = A[i];
    }
}

/*free matrix*/
void freeMatrix(double **temp,int n)
{
    int i = 0;
    if(temp == NULL)
        return;
    for(i = 0; i<n; i++)
    {
        if(temp[i] != NULL)
            free(temp[i]);
    }
    free(temp);
}

/*print matrix, each Aij with 4 digits after point*/
void printMatrix(double **A, int n, int d)
{
    int i = 0;
    int j = 0;
    for (i=0; i< n; i++){
        for (j=0; j< d-1;j++){
            printf("%.4f,", A[i][j]);
        }
        printf("%.4f\n", A[i][d-1]);
    }
}

/*print values of array with 4 digits after point*/
void printVector(double *vector,int d)
{
    int i;
    for(i=0; i<d; i++)
    {
        if(i==d-1)
            printf("%.4f\n",vector[i]);
        else
            printf("%.4f,",vector[i]);
    }
}
/*return the wam matrix of A matrix*/
double **wam(double **A,int n, int d)
{
    double **wamMatrix;
    int i,j,l;
    double norm = 0;
    wamMatrix = createZeroMatrix(n,n);
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
        {
            if(i==j)
                wamMatrix[i][j] = 0;
            else
            {
                norm = 0;
                for(l=0; l<d; l++)
                {
                    norm += pow(A[i][l] - A[j][l],2);
                }
                wamMatrix[i][j] = exp(-1*norm/2);
            }
        }

    return wamMatrix;
}
/*print the wam matrifx of A*/
void printWam(double** A,int n, int d)
{
    double **resultWam;
    resultWam = wam(A,n,d);
    printMatrix(resultWam,n,n);
    freeMatrix(resultWam,n);
}

/* create Diagonal from A matirx */
double **ddg(double **A,int n , int d)
{
    int i,j;
    double **resultDdg,**wamMatrix;
    wamMatrix = wam(A,n,d);
    
    resultDdg = createZeroMatrix(n,n);
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            resultDdg[i][i] += wamMatrix[i][j];
        }
    }
    freeMatrix(wamMatrix,n);
    return resultDdg;
}

/*print the diagonal matrifx of wam of A*/
void printDdg(double **A, int n, int d)
{
    double **ddgMatrix;
    ddgMatrix = ddg(A,n,d);
    printMatrix(ddgMatrix,n,n);
    freeMatrix(ddgMatrix,n);

}

/*return the Laplacian Graph from A matrix*/
double** gl(double **A,int n, int d)
{
    double** wamMatrix, **glMatrix, **ddgMatrix;
    glMatrix = createZeroMatrix(n,n);
    wamMatrix = wam(A,n,d);
    ddgMatrix = ddg(A,n,d);
    subsMatrix(glMatrix,ddgMatrix,wamMatrix,n,n);
    freeMatrix(ddgMatrix,n);
    freeMatrix(wamMatrix,n);
    return(glMatrix);

}

/*print the Laplacian Graph from A matrix*/
void printGl(double **A,int n, int d)
{
    double **glMatrix;
    glMatrix = gl(A,n,d);
    printMatrix(glMatrix,n,n);
    freeMatrix(glMatrix,n);
}


/*
    jacobi auxiliary functions:
*/



/*create rotation P matrix with given indices and c,s*/
double **createMatrixP(int n, int i, int j, double c, double s)
{
    double **matrixP = createOneMatrix(n);

    matrixP[i][i] = c;
    matrixP[j][j] = c;
    matrixP[i][j] = s;
    matrixP[j][i] = -1 * s;
    return matrixP;
}

/*return abs value of num*/
double calcAbs(double num)
{
    if(num>0)
        return num;
    return -1*num;
}

/*calc Aij: the off diagonal element with largest abs value*/
int* calcPivotIndices(double **A,int n)
{
    int i,j;
    double max = -1;
    double curr = 0;
    int *indices = (int*)malloc(2*sizeof(int));
    if(indices == NULL)
    {
        printf("An Error Has Occurred\n");
        freeMatrix(A,n);
        exit(1);
    }
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            if(i!=j)
            {
                curr = calcAbs(A[i][j]);
                if(curr > max)
                {
                    indices[0] = i;
                    indices[1] = j;
                    max = curr;
                }
            }
                
    return indices;
}

/*return 1 if n>=0 else -1*/
int sign(double n)
{
    if(n >= 0 || n == -0)
        return 1;
    return -1;
}

/*calculating c and s for roatation P matrix*/
double* obtainCandS(double Aii,double Aij,double Ajj)
{
    double c,t;
    double theta;
    double *cAndS = (double*)malloc(sizeof(double)*2);
    if(cAndS == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    theta = (Ajj - Aii)/(2*Aij);
    t =  sign(theta)/(fabs(theta) + sqrt(pow(theta,2)+1));
    c = 1/sqrt(t*t+1);
    cAndS[0] = c;
    cAndS[1] = t*c; /*t*c==s*/
    return cAndS;
}


/*calculating A' = P^t * A * P, V = P1 * P2 * P3 .... all in place */
void pivotMoveJacobi(double **A, double **mulPmatrix, int n)
{
    double c,s;
    int *indices;
    int l,i_M,j_M,k;
    double **matrixP,**resultMatrix,**AtagMatrix;
    double *cAnds;
    AtagMatrix = createZeroMatrix(n,n);
    indices = calcPivotIndices(A,n);
    i_M = indices[0];
    j_M = indices[1];
    cAnds = obtainCandS(A[i_M][i_M],A[i_M][j_M],A[j_M][j_M]);
    c = cAnds[0];
    s = cAnds[1];
    matrixP = createMatrixP(n,i_M,j_M,c,s);
    /* claculating V = P1 * P2 * .... inplace */
    resultMatrix = mulMatrix(mulPmatrix,matrixP,n,n,n); /*res = mulP * P*/
    copyMatrix(resultMatrix,mulPmatrix,n,n);/*from res to mulP*/
    /*calc A' in place*/
    for (k = 0; k < n; k++ )
        for (l = 0; l < n; l++)
            if (k != i_M && k != j_M) 
            {
                if (l==i_M)
                    {
                    AtagMatrix[i_M][k]=c*A[i_M][k] - s*A[j_M][k];
                    AtagMatrix[k][i_M]=c*A[i_M][k] - s*A[j_M][k];
                } 
                else if (l==j_M)
                    {
                    AtagMatrix[j_M][k]=s*A[i_M][k] + c*A[j_M][k];
                    AtagMatrix[k][j_M]=s*A[i_M][k] + c*A[j_M][k];
                } 
                else 
                    AtagMatrix[k][l] = A[k][l];
            } 
    AtagMatrix[i_M][i_M] = c*c*A[i_M][i_M] + s*s*A[j_M][j_M] - 2 * s * c * A[i_M][j_M];
    AtagMatrix[j_M][j_M] = s*s*A[i_M][i_M] + c*c*A[j_M][j_M] + 2 * s * c * A[i_M][j_M];
    AtagMatrix[i_M][j_M] = (c*c - s*s)*(A[i_M][j_M]) + s*c*(A[i_M][i_M]-A[j_M][j_M]);
    AtagMatrix[j_M][i_M] = (c*c - s*s)*(A[i_M][j_M]) + s*c*(A[i_M][i_M]-A[j_M][j_M]);
    copyMatrix(AtagMatrix,A,n,n);
    freeMatrix(resultMatrix,n);
    freeMatrix(AtagMatrix,n);
    freeMatrix(matrixP,n);
    free(cAnds);
    free(indices);
}

/*return off(A)*/
double calcOffConvergence(double **A,int n)
{
    double offA=0;
    int i,j;
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            if(i!=j)
                offA += pow(A[i][j],2);
    return offA;
}

/*calculate and return n*n+1 matrix, last row is eignValues, other is eigenVectors*/
double **jacobi(double **A,int n ,int d)
{
    int i;
    double **vMatrix,**AtagMatrix;
    double epsilon = 1;
    int iter = 1;
    double offA,offAtag;
    double *eigenValues;
    if(d!=n)
    {
        printf("An Error Has Occurred\n");
        freeMatrix(A,n);
        exit(1);
    }
    vMatrix = createOneMatrix(n);
    AtagMatrix = createZeroMatrix(n,n);
    copyMatrix(A,AtagMatrix,n,n);
    while(epsilon > 0.00001 && iter <= 100)
    {
        
        offA = calcOffConvergence(AtagMatrix,n);
        pivotMoveJacobi(AtagMatrix,vMatrix,n);
        offAtag = calcOffConvergence(AtagMatrix,n);
        epsilon = offA - offAtag;
        iter++;
    }
    eigenValues = matrixDiag(AtagMatrix,n);
    freeMatrix(AtagMatrix,n);
    AtagMatrix = createZeroMatrix(n+1,n);
    copyMatrix(vMatrix,AtagMatrix,n,n);
    freeMatrix(vMatrix,n);
    for(i=0; i<n; i++)
        AtagMatrix[n][i] = eigenValues[i];
    free(eigenValues);
    return AtagMatrix;
}

/*return the min value of array*/
double minList(double* vector, int n)
{
    int i;
    double min;
    min = vector[0];
    for(i = 1; i < n; i++){
        if(vector[i] < min){
            min = vector[i];
        }
    }
    return min;
}

double minList(double* vector, int n);
double argMinList(double* vector, int n);


/*return the index of the min element int array*/
double argMinList(double* vector, int n)
{
    int i,k;
    double min;
    min = vector[0];
    k = 0;
    for(i = 1; i < n; i++){
        if(vector[i] < min){
            min = vector[i];
            k = i;
        }
    }
    return k;
}

/*return K number of clusetr as described in the eigengam heruistic
sorting the eigenValues(and vecotr) ascending and returning the argmax(delta) 
*/
int calcK(double* eigenValues,int n)
{
    int i,maxIdx;
    double maxDelta;
    double *deltaArray = (double*)malloc(sizeof(double) * (n-1));
    if(deltaArray == NULL)
    {
        printf("An Error Has Occurred\n");
        free(eigenValues);
        exit(1);
    }
    for(i=0; i<n-1; i++)
        deltaArray[i] = calcAbs(eigenValues[i] - eigenValues[i+1]);
    maxDelta = deltaArray[0];
    maxIdx = 0;
    for(i=0;i<n/2;i++)
        if(deltaArray[i]>maxDelta)
        {
            maxDelta = deltaArray[i];
            maxIdx = i;
        }
    free(deltaArray);
    return maxIdx+1;
}

/*sort matrix in accordance to eigen values in vecotr in ascending order*/
double **sortMatrix(double **A,double *vector, int n,int d)
{
    double tmp;
    double *tmpVector;
    double **Atranspoes = transposeMatrix(A,n,d);
    int i,j;
    for(i=0; i<d; i++)
        for(j=i+1; j<n;j++)
            if(vector[i]>vector[j])
            {
                tmp = vector[i];
                tmpVector = Atranspoes[i];
                vector[i] = vector[j];
                Atranspoes[i] = Atranspoes[j];
                vector[j] = tmp;
                Atranspoes[j] = tmpVector;

            }
    Atranspoes = transposeMatrix(Atranspoes,d,n);
    return Atranspoes;

}

/*full spk as described, return new Matrix*/
double **spk(double **A,int n, int d,int *k)
{

    double **Umatrix,**jacobiMatrix,**glMatrix,**sortedMatrix;
    glMatrix = gl(A,n,d);
    jacobiMatrix = jacobi(glMatrix,n,n);   
    sortedMatrix = sortMatrix(jacobiMatrix,jacobiMatrix[n],n,n);
    if( *k == -1)
        *k = calcK(jacobiMatrix[n],n);
    Umatrix = createZeroMatrix(n,*k);
    copyMatrix(sortedMatrix,Umatrix,n,*k);
    freeMatrix(jacobiMatrix,n+1);
    freeMatrix(sortedMatrix,n);
    freeMatrix(glMatrix,n);
    return Umatrix;
}

/*print jacboi, first row is eigen values, after its the eigen vecotr matrix as column*/
void printJacobi(double **A,int n, int d)
{
    int i;
    double **jacobiRes;
    jacobiRes = jacobi(A,n,d);
    
    for(i=0; i <n; i++)
    {
        if(i!=n-1)
            printf("%.4f,", jacobiRes[n][i]);
        else
            printf("%.4f\n", jacobiRes[n][i]);
    }
    printMatrix(jacobiRes,n,n);
    freeMatrix(jacobiRes,n+1);
}



