#define PY_SSIZE_T_CLEAN 
#include <Python.h>
#include "kmeans.h"

double *datapoints_c,*centroids_c,*ret_centroids;
int i,j;
static double* getDataPoints(PyObject *datapoints_py,int n_rows, int n_col)
{
    
    datapoints_c = (double *)malloc(n_col*n_rows*sizeof(double));
    PyObject *dataRow, *dataPoint;

    if(datapoints_c == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    for(i=0; i<n_rows; i++)
    {
        dataRow = PyList_GetItem(datapoints_py, i);
        for(j=0; j<n_col; j++)
        {
            dataPoint = PyList_GetItem(dataRow,j);
            datapoints_c[i*n_col + j] = PyFloat_AsDouble(dataPoint);
        }
    }
    return datapoints_c;
}
static double* getCentroids(PyObject *centroids_py,int k, int n_col){
    
    centroids_c = (double *)malloc(n_col*k*sizeof(double));
    PyObject *dataRow, *dataPoint;

    if(centroids_c == NULL){
        printf("An Error Has Occurred\n");
        free(datapoints_c);
        exit(1);
    }
    
    for(i=0; i<k; i++)
    {
        dataRow = PyList_GetItem(centroids_py, i);
        for(j=0; j<n_col; j++)
        {
            dataPoint = PyList_GetItem(dataRow,j);
            centroids_c[i*n_col + j] = PyFloat_AsDouble(dataPoint);
        }
    }
    return centroids_c;
}

static PyObject* kmeans_func(PyObject *self, PyObject *args)
{
    int k,n_rows,n_col,max_iter;
    double eps;
    PyObject *datapoints_py, *centroids_py;
    if(!PyArg_ParseTuple(args, "iiiidOO", &k, &n_rows,&n_col,&max_iter,&eps,&datapoints_py,&centroids_py)) {
        return NULL; 
    }
    if (!PyList_Check(datapoints_py)){
        printf("111\n");
        return NULL;
    }
    if(!PyList_Check(centroids_py)){
        printf("222\n");
        return NULL;
    }
    /*transfer pyton data to c data*/
    datapoints_c = getDataPoints(datapoints_py,n_rows,n_col);
    centroids_c = getCentroids(centroids_py,k,n_col);

    /*calling the kmeans func*/
    ret_centroids = kmeans_calc(k,n_rows,n_col,max_iter,eps,datapoints_c,centroids_c);

    PyObject* ret_pythonList = PyList_New(k);
    /*transfering the data back to python object*/
    for(i = 0; i< k; i++)
    {
        PyObject* row_list = PyList_New(n_col);
        for(j = 0; j< n_col; j++)
        {
            PyObject* point = Py_BuildValue("d", ret_centroids[i*n_col + j]);
            PyList_SetItem(row_list,j,point);
        }
        PyList_SetItem(ret_pythonList,i,row_list);
    }
    if(centroids_c == NULL)
    {
        free(centroids_c);
    }
    if(datapoints_c == NULL)
    {
        free(datapoints_c);
    }

    return Py_BuildValue("O", ret_pythonList); 
}

static PyMethodDef kmeansMethods[] = {
    {"fit",                                                    
     (PyCFunction)kmeans_func,                                      
     METH_VARARGS,                                             
     PyDoc_STR("kmeans algorithm for given datapoints and info")}, 
    {NULL, NULL, 0, NULL}                                      
};

static struct PyModuleDef kmeanmoudle = {
    PyModuleDef_HEAD_INIT,
    "kmeans_c", 
    NULL,      
    -1,         
    kmeansMethods  
};

PyMODINIT_FUNC
PyInit_kmeans_c(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeanmoudle);
    if(!m)
    {
        return NULL;
    }
    return m;
}
