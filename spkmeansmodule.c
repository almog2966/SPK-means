#define PY_SSIZE_T_CLEAN 
#include <Python.h>
#include "kmeans.h"
#include "spkmeans.h"



int i,j,n,d,N,*k,kInput;
int isjacobi;
int isSpk;

/*translate the Pyobject to 2D matrix, insert n,d sizes*/
double **getDataPoints(PyObject *self, PyObject *args)
{
    double **matrix;
    double *row;
    PyObject *dataRow, *dataPoints, *point;
    if(isSpk)
    {
        if(!PyArg_ParseTuple(args, "iO",&kInput,&dataPoints)) 
            return NULL;
    }
    else
    {
        kInput = -1;
        if(!PyArg_ParseTuple(args, "O",&dataPoints)) 
            return NULL;
    }
    k = &kInput;
    if (!PyList_Check(dataPoints))
        return NULL;
    n = PyList_Size(dataPoints);
    dataRow = PyList_GetItem(dataPoints,0);
    d = PyList_Size(dataRow);
    matrix = (double**)malloc(sizeof(double*)*n);
    if(matrix == NULL)
    {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i=0; i<n; i++)
    {
        row = (double*)malloc(sizeof(double)*d);
        if(matrix == NULL)
        {   
            freeMatrix(matrix,n);
            printf("An Error Has Occurred\n");
            exit(1);
        }
        matrix[i] = row;
    }
    for(i=0; i<n; i++)
    {
        dataRow = PyList_GetItem(dataPoints,i);
        for(j=0; j<d; j++)
        {
            point = PyList_GetItem(dataRow,j);
            matrix[i][j] = PyFloat_AsDouble(point);
        }
    }
        
    return matrix;
}

/*translate 2D matrix to PyObject insert n,d sizes*/
static PyObject* retDataPoints(double **matrix)
{
    PyObject *retMatrix;
    int vectorSize;
    if(isSpk)
        vectorSize = *k;
    else
        vectorSize = d;
    if(isjacobi == 1)
        n++;
    retMatrix = PyList_New(n);
    for(i=0; i<n; i++)
    {
        PyObject *row = PyList_New(vectorSize);
        for(j=0; j<vectorSize; j++)
        {
            PyObject *point = Py_BuildValue("d",matrix[i][j]);
            PyList_SetItem(row,j,point);

        }
        PyList_SetItem(retMatrix,i,row);
    }
    return Py_BuildValue("O",retMatrix);

}

/*called form python, calculate and return PyObject of wam matrix*/
static PyObject *wamPython(PyObject *self, PyObject *args){
    PyObject *returnData;
    double **matrix,**wamMatrix;
    isjacobi = 0;
    isSpk = 0;
    matrix = getDataPoints(self, args);
    wamMatrix = wam(matrix, n, d);
    d=n;
    returnData = retDataPoints(wamMatrix);
    freeMatrix(wamMatrix,n);
    freeMatrix(matrix, n);   
    return returnData;
}

/*called form python, calculate and return PyObject of spk matrix*/
static PyObject *spkPython(PyObject *self, PyObject *args){
    PyObject *returnData;
    double **matrix,**centoids;
    isjacobi = 0;
    isSpk = 1;
    matrix = getDataPoints(self, args);
    centoids = spk(matrix, n, d,k);
    returnData = retDataPoints(centoids);
    freeMatrix(centoids,n);
    freeMatrix(matrix, n);   
    return returnData;
}

/*called form python, calculate and return PyObject of ddg matrix*/
static PyObject *ddgPython(PyObject *self, PyObject *args){
    PyObject *returnData;
    double **matrix,**ddgMatrix;
    isjacobi = 0;
    isSpk = 0;
    matrix = getDataPoints(self, args);
    ddgMatrix = ddg(matrix, n, d);
    d=n;
    returnData = retDataPoints(ddgMatrix);
    freeMatrix(ddgMatrix,n);
    freeMatrix(matrix, n);   
    return returnData;
}

/*called form python, calculate and return PyObject of gl matrix*/
static PyObject *glPython(PyObject *self, PyObject *args){
    PyObject *returnData;
    double **matrix,**glMatrix;
    isjacobi = 0;
    isSpk = 0;
    matrix = getDataPoints(self, args);
    glMatrix = gl(matrix, n, d);
    d=n;
    returnData = retDataPoints(glMatrix);
    freeMatrix(matrix, n);   
    return returnData;
}

/*called form python, calculate and return PyObject of jacobi matrix*/
static PyObject *jacobiPython(PyObject *self, PyObject *args){
    isjacobi = 1;
    isSpk = 0;
    PyObject *returnData;
    double **matrix,**jacobiMatrix;
    matrix = getDataPoints(self, args);
    jacobiMatrix = jacobi(matrix, n, d);
    returnData = retDataPoints(jacobiMatrix);
    freeMatrix(jacobiMatrix,n);
    freeMatrix(matrix, n-1);  
    return returnData;
}



static PyMethodDef spkmeansMethods[] = {
    {"spk",(PyCFunction)spkPython,METH_VARARGS,PyDoc_STR("spk function")},
    {"wam",(PyCFunction)wamPython,METH_VARARGS,PyDoc_STR("wam function")},
    {"ddg",(PyCFunction)ddgPython,METH_VARARGS,PyDoc_STR("ddg function")},
    {"gl",(PyCFunction)glPython,METH_VARARGS,PyDoc_STR("gl function")},
    {"jacobi",(PyCFunction)jacobiPython,METH_VARARGS,PyDoc_STR("jacobi function")},
      
    {NULL, NULL, 0, NULL}                                     
};

static struct PyModuleDef spkmeanmoudle = {
    PyModuleDef_HEAD_INIT,
    "myspkmeans", 
    NULL,      
    -1,         
    spkmeansMethods  
};

PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&spkmeanmoudle);
    if(!m)
    {
        return NULL;
    }
    return m;
}

