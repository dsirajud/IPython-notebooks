// 1. Include Python API and required libraries
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

// 2. function declaration
static PyObject *evolve(PyObject *self, PyObject *args);

// 3. methods table
static PyMethodDef Evolve1Methods[] = {
  { "evolve", 
    evolve,
    METH_VARARGS,
    "Doc string."},
  { NULL, NULL, 0, NULL }
};

// 5. function initialization

PyMODINIT_FUNC initevolve(void) {
  (void) Py_InitModule("evolve", methods);
  import_array();
}

/*******************************************************************************/
/* 1.5 define convenient macros for array indexing 
 *     (must be named py_m for 1D array, or py_r for 2D)
 */

#define m(x0) (*(npy_float64*)((PyArray_DATA(py_m) +                \
                                (x0) * PyArray_STRIDES(py_m)[0])))

#define m_shape(i) (py_m->dimensions[(i)])

#define r(x0, x1) (*(npy_float64*)((PyArray_DATA(py_r) +                \
                                    (x0) * PyArray_STRIDES(py_r)[0] +   \
                                    (x1) * PyArray_STRIDES(py_r)[1])))

#define r_shape(i) (py_r->dimensions[(i)])

#define F(x0, x1) (*(npy_float64*)((PyArray_DATA(py_F) +                \
                                    (x0) * PyArray_STRIDES(py_F)[0] +   \
                                    (x1) * PyArray_STRIDES(py_F)[1])))

#define F_shape(i) (py_F->dimensions[(i)])

#define v(x0, x1) (*(npy_float64*)((PyArray_DATA(py_v) +                \
                                    (x0) * PyArray_STRIDES(py_v)[0] +   \
                                    (x1) * PyArray_STRIDES(py_v)[1])))

#define v_shape(i) (py_v->dimensions[(i)])

/*  example usage of above macros:
 *
 *  can loop through a 2D numpy array and intialize to zero per:
 *
 * for (i = 0; i < r_shape(0); ++i) 
 * {
 *     for (j = 0; j < r_shape(1); ++j)
 *        r(i,j)= 0; // zero all elements;
 * }
 *
 */

/*******************************************************************************/







// 4. function definitions 

/* evolve1 communicates with the Python interpreter
 * by taking Python inputs, parsing to C objects
 * computation is done in C, and evolve1 updates
 * translates the C objects to PyObjects
 */

static inline void compute_F(npy_int64 N,
			     PyArrayObject *py_m,
			     PyArrayObject *py_r,
			     PyArrayObject *py_F)
{

  npy_int64 i, j;
  npy_float64 rx, ry, Fx, Fy, r3, tmp;

  // set all forces to zero
  for (i = 0; i < N; ++i)
    {
      F(i,0) = F(i,1) = 0;
    }


  for (i = 0; i < N; ++i)
    {
      for (j = i+1; j < N; ++j)
	{
	  rx = r(i,0) - r(j,0);
	  ry = r(i,1) - r(j,1);

	  r3 = sqrt(rx * rx + ry * ry);
	  r3 = r3 * r3 * r3; // |r| ** 3

          tmp = m(i)*m(j) / r3;
          Fx += tmp * rx;
          Fy += tmp * ry;

	  F(i,0) += Fx;
	  F(i,1) += Fy;

	  // Newton's 3rd: F(i,j) = -F(j,i)
	  F(j,0) -= Fx;
	  F(j,1) -= Fy;

	}
    }

static PyObject * evolve (PyObject * self, PyObject * args)
{

  // Declare variables
  npy_int64 N, steps, step, i;
  npy_float64 dt;
  PyArrayObject *py_m, *py_r, *py_v, *py_F;

  // Parse argments
  if (!PyArg_ParseTuple(args, "ldllO!O!O!O!",
			&dt,
			&steps,
			&N,
			&PyArray_Type, &py_m,
			&PyArray_Type, &py_r,
			&PyArray_Type, &py_v,
			&PyArray_Type, &py_F))
    {
      return NULL;
    }


  // Evolve world instance
  for (step = 0; step < steps; ++steps)
    {
      compute_F(N, py_m, py_r, py_F);

      for (i = 0; i < N; ++i)
	{
	  v(i,0) += F(i,0) * dt / m(i);
	  v(i,1) += F(i,1) * dt / m(i);

	  r(i,0) += v(i,0) * dt;
	  r(i,1) += v(i,1) * dt;
	}
    }
  Py_RETURN_NONE;
}


