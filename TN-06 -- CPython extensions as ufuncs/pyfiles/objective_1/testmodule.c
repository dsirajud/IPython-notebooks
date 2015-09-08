// 1. Include Python API and required libraries
#include <Python.h>
#include <math.h>

/*
 * expsinmodule.c
 * This is the C code for a non-numpy Python extension to
 * define the expsin function, where expsin(p) = exp(2*x) + sin(x).
 * This function will not work on numpy arrays automatically.
 * numpy.vectorize must be called in python to generate
 * a numpy-friendly function.
 *
 * Details explaining the Python-C API can be found under
 * 'Extending and Embedding' and 'Python/C API' at
 * docs.python.org .
 */


// 2. Function declaration
static PyObject * test_expsin (PyObject * self, PyObject * args);

// 3. Methods table
static PyMethodDef TestMethods[] = {
  {"expsin",
   test_expsin,
   METH_VARARGS,
   "returns the the value expsin(p) = exp(2*x) + sin(x)"},
  {NULL, NULL, 0, NULL}

};

// 4. Function definition
  static PyObject * test_expsin (PyObject * self, PyObject * args)
  {

    double p;

    if (!PyArg_ParseTuple(args, "d", &p))
	return NULL;

    p = exp(2*p) + sin(p);

    return Py_BuildValue ("d", p);
	
  }


// 5. Initialization function 
PyMODINIT_FUNC
inittest(void)
{
  PyObject * m;
  m = Py_InitModule("test", TestMethods);
  if (m == NULL) 
    return;
}
