#include <Python.h>
#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

static size_t 
find_pivot_row(double **A, size_t column, size_t num_rows, size_t row)
{
    size_t i;
    for (i = row; i < num_rows; i++){
        if (A[i][column] != 0.0f){
            return i;
        }
    }
    return num_rows;
}

static void 
row_echelon_transformation(double **A, size_t rows, size_t cols)
{
    size_t row_cnt = 0, col_cnt = 0, pivot_row;
    double pivot, front, *swap;

    while((row_cnt < (rows - 1)) && (col_cnt < (cols - 1))) {

        pivot_row = find_pivot_row(A, col_cnt, rows, row_cnt);

        if (row_cnt < pivot_row && pivot_row < rows){
            swap = A[row_cnt];
            A[row_cnt] = A[pivot_row];
            A[pivot_row] = swap;
        }
        else if (pivot_row == rows){
            col_cnt++;
            continue;
        }

        pivot = A[row_cnt][col_cnt];
        size_t i, j;
        for (i = row_cnt + 1; i < rows; i++){
            front = A[i][col_cnt];
            for (j = 0; j < cols; j++){
                A[i][j] = A[i][j]  - (front/pivot)*(A[row_cnt][j]);
            }
        }
        row_cnt++; col_cnt++;
    }
}

static PyObject *
row_echelonize(PyObject *self, PyObject *args)
{
    PyObject *matrix = NULL;
    PyObject *row = NULL;

    if (!PyArg_ParseTuple(args, "O", &matrix)){
        return Py_None;
    }
    Py_INCREF(matrix);

    if (PyList_Check(matrix)) {
        row = PyList_GetItem(matrix, 0);
        Py_INCREF(row);
        if (!PyList_Check(row)){
            Py_DECREF(row);
            return matrix; /* matrix is no good */
        }
    }
    else {
        return matrix; /* matrix is no good */
    }

    size_t rows;
    size_t cols;


    assert(row != NULL);
    rows = (size_t)PyList_GET_SIZE(matrix);
    cols = (size_t)PyList_GET_SIZE(row);

    double **A = (double **)malloc(sizeof(double *)*rows);
    assert(A != NULL);
    size_t i;
    for (i = 0; i < cols; i++){
        A[i] = (double *)malloc(sizeof(double)*cols);
        assert(A[i] != NULL);
    }

    PyObject *item;
    size_t ii, jj;
    for (ii = 0; ii < rows; ii++){
        row = PyList_GetItem(matrix, ii);
        for (jj = 0; jj < cols; jj++){
            item = PyList_GetItem(row, jj);
            if (PyFloat_Check(item)){
                A[ii][jj] = PyFloat_AsDouble(item);
            }
            else if (PyLong_Check(item)){
                A[ii][jj] = PyLong_AsDouble(item);
            }
        }
    }
    Py_DECREF(matrix);

    row_echelon_transformation(A, rows, cols);

    PyObject *value;
    matrix = PyList_New(rows);
    for (ii = 0; ii < rows; ii++){
        row = PyList_New(cols);
        for (jj = 0; jj < cols; jj++){
            value = PyFloat_FromDouble(A[ii][jj]);
            PyList_SetItem(row, jj, value);
        }
        PyList_SetItem(matrix, ii, row);
    }
    return matrix;
}

static PyMethodDef GaussianMethods[] = {
    {"c_row_echelonize", row_echelonize, METH_VARARGS, "transforms a matrix into row-echelon form"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef gaussianmodule = {
    PyModuleDef_HEAD_INIT,
    "c_gaussian",
    NULL, 
    -1,
    GaussianMethods
};

PyMODINIT_FUNC
PyInit_c_gaussian(void)
{
    return PyModule_Create(&gaussianmodule);
}