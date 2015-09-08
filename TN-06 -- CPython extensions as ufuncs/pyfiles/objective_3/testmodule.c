// 1. Include Python API and required libraries
#include "Python.h"
#include "math.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/halffloat.h"

/*
 * testmodule.c
 * This is the C code for creating your own
 * Numpy ufunc for a expsin function.
 *
 * Each function of the form type_expsin defines the
 * expsin function for a different numpy dtype. Each
 * of these functions must be modified when you
 * create your own ufunc. The computations that must
 * be replaced to create a ufunc for
 * a different funciton are marked with BEGIN
 * and END.
 *
 * Details explaining the Python-C API can be found under
 * 'Extending and Embedding' and 'Python/C API' at
 * docs.python.org .
 *
 */

// 3. Methods table
static PyMethodDef TestMethods[] = {
    {NULL, NULL, 0, NULL}
};


/* The loop definition must precede the PyMODINIT_FUNC. */
static void long_double_expsin(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char * in = args[0];
    char * out = args[1];
    npy_intp in_step = steps[0]; 
    npy_intp out_step = steps[1];
    
    long double tmp;

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        tmp = *(long double *)in;
	*((long double *)out)  = exp(2*tmp) + sin(tmp);
        /*END main ufunc computation*/

        in += in_step; 
        out += out_step;
    }
}

static void double_expsin(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char * in = args[0];
    char * out = args[1];
    npy_intp in_step = steps[0]; 
    npy_intp out_step = steps[1];
    
    double tmp;

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        tmp = *(double *)in;
	*((double *)out)  = exp(2*tmp) + sin(tmp);
        /*END main ufunc computation*/

        in += in_step; 
        out += out_step;
    }
}

static void float_expsin(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char * in = args[0];
    char * out = args[1];
    npy_intp in_step = steps[0]; 
    npy_intp out_step = steps[1];
    
    float tmp;

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        tmp = *(float *)in;
	*((float *)out)  = exp(2*tmp) + sin(tmp);
        /*END main ufunc computation*/

        in += in_step; 
        out += out_step;
    }
}

static void half_float_expsin(char **args, npy_intp *dimensions,
                            npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char * in = args[0];
    char * out = args[1];
    npy_intp in_step = steps[0]; 
    npy_intp out_step = steps[1];
    
    float tmp;

    for (i = 0; i < n; i++) {
        /*BEGIN main ufunc computation*/
        tmp = *(npy_half *)in;
        tmp = npy_half_to_float(tmp);
	tmp  = exp(2*tmp) + sin(tmp);
        *((npy_half *)out) = npy_float_to_half(tmp);
        /*END main ufunc computation*/

        in += in_step; 
        out += out_step;
    }
}


/*This a pointer to the above function*/
PyUFuncGenericFunction funcs[4] = {&long_double_expsin,
				     &double_expsin,
				     &float_expsin,
				     &half_float_expsin};

/* These are the input and return dtypes of expsin.*/
static char types[8] = {NPY_LONGDOUBLE, NPY_LONGDOUBLE,
			  NPY_DOUBLE, NPY_DOUBLE,
			  NPY_FLOAT, NPY_FLOAT,
			  NPY_HALF, NPY_HALF};

/* This is extra parameters, if any, needed to compute the function eval. */
static void *data[4] = {NULL,
			  NULL,
			  NULL,
			  NULL};

PyMODINIT_FUNC initnpufunc(void)
{
    PyObject *m, *expsin, *d;


    m = Py_InitModule("npufunc", TestMethods);
    if (m == NULL) {
        return;
    }

    import_array();
    import_umath();

    expsin = PyUFunc_FromFuncAndData(funcs, data, types, 4, 1, 1,
                                    PyUFunc_None, "expsin",
                                    "expsin_docstring", 0);

// ntypes = 4 varieties of this ufunc (4 dtypes)

    d = PyModule_GetDict(m);

    PyDict_SetItemString(d, "expsin", expsin);
    Py_DECREF(expsin);
}
