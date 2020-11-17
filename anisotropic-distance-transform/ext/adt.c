/*
 * License: MIT
 *
 * Ideas taken from:
 * http://www.cs.brown.edu/~pff/dt/
 * http://www.cs.auckland.ac.nz/~rklette/TeachAuckland.html/mm/MI30slides.pdf
 *
 * Original implementation: Jan Hosang, https://github.com/hosang/adt
 * Anisotropic extension by Sebastian Kosch
 */

#include "Python.h"
#include "numpy/arrayobject.h"

#include <math.h>

#define SQ(s)s*s

#ifndef INFINITY
#define INFINITY 1.0/0.0;
#endif

// declaration up here; actual definition is at the bottom
static PyObject *adt(PyObject *self, PyObject *args);

static PyMethodDef adt_methods[] = {
  { "adt",adt, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef mod_def =
  {
   PyModuleDef_HEAD_INIT,
   "adt", /* name of module */
   "",      /* module documentation, may be NULL */
   -1,      /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   adt_methods
  };

PyMODINIT_FUNC PyInit_adt(void) {
    PyObject *m;
    m = PyModule_Create(&mod_def);
    if (!m) return NULL;
    import_array();
    return m;
}

/*
 * Function: dt1d_isotropic
 * ----------------------------
 * Performs a squared distance transform on a 1D array of values f, writing the result to out.
 * This version is *isotropic*, i.e. the calculation is the same in all directions. This function
 * is called if alpha == beta. See below for the *anisotropic* version.
 *
 *   f: A 1D array of initial distances.
 *   out: A 1D array that the result will be written to.
 *   len: The number of (valid) values in f.
 *   alpha: a multiplier stretching each parabola vertically
 *   _beta: assumed equal to alpha, thus irrelevant
 *   v: a pre-allocated buffer to keep track of the vertices of the parabolas that form the envelope.
 *   z: a pre-allocated buffer to keep track of the intersections between then parabolas v[k].
 */
static void dt1d_isotropic(float *f, float *out, int len, float alpha, float _beta, int *v, float *z) {
    int k = 0; // the number of parabolas forming the lower envelope
    v[0] = 0;
    z[0] = -INFINITY; // v[0] is valid from z[0] to z[1], etc.
    z[1] = INFINITY;
    int q;
    float s;
    for (q = 1; q < len; q++) { // go through all points, left to right (or top to bottom)
        while (1) {
            s = ((f[q]/alpha + q * q) - (f[v[k]]/alpha + v[k] * v[k])) / (2 * (q - v[k]));

            // s could be to the left of the last "valid from"
            // if that is so,
            if (s <= z[k]) {
                k--; // if wraps the last parabola, get rid of that last one
                if (k < 0) {
                    k = 0;
                    v[k] = q;
                }
                continue;
            } else {
                k++;
                v[k] = q;
                z[k] = s;
                z[k+1] = INFINITY;
                break;
            }
        }
    }

    k = 0; // start at the left
    for (q = 0; q < len; q++) {
        while (z[k+1] < q) {
            k++;
        }
        out[q] = alpha * (q-v[k])*(q-v[k]) + f[v[k]];
    }
}


/*
 * Function: dt1d_anisotropic
 * ----------------------------
 * Performs a squared distance transform on a 1D array of values f, writing the result to out.
 * This version is *anisotropic*, i.e. forward distances are weighed differently than backward distances.
 * This is called if alpha != beta.
 *
 *   f: A 1D array of initial distances.
 *   out: A 1D array that the result will be written to.
 *   len: The number of (valid) values in f.
 *   alpha: A factor applied to the *left* arm of every parabola.
 *   beta: A factor applied to the *right* arm of every parabola.
 *   v: a pre-allocated buffer to keep track of the vertices of the parabolas that form the envelope.
 *   z: a pre-allocated buffer to keep track of the intersections between then parabolas v[k].
 */
static void dt1d_anisotropic(float *f, float *out, int len, float alpha, float beta, int *v, float *z) {
    int k = 0; // the number of parabolas forming the lower envelope
    v[0] = 0;
    z[0] = -INFINITY; // v[0] is valid from z[0] to z[1], etc.
    z[1] = INFINITY;
    int q;
    float s;

    for (q = 1; q < len; q++) { // go through all points, left to right (or top to bottom)
        while (1) {
            // Find the intersection between this new parabola q and the last (relevant) one v[k].
            // There are three situations we need to consider:
            //   - right arm of v[k] intersects with right arm of q
            //   - left arm of v[k] intersects with left arm of q
            //   - right arm of v[k] intersects with left arm of q

            // We can pre-compute a few values that will be reused:
            float qmv = (q - v[k]);
            float qmvt2 = 2*(q - v[k]);
            float qmvsq = qmv * qmv;

            float fqda = f[q] / alpha;
            float fqdb = f[q] / beta;
            float fvda = f[v[k]] / alpha;
            float fvdb = f[v[k]] / beta;

            if (qmvsq + fvdb < fqdb) {
                // right arm of v[k] intersects with right arm of q
                // (if f[q] is high enough, relative to f[v[k]], that its vertex is inside v[k])
                s = ((fqdb + q * q) - (fvdb + v[k]*v[k])) / qmvt2;

            } else if (qmvsq + fqda < fvda) {
                // left arm of v[k] intersects with left arm of q
                // (if f[v[k]] is high enough, relative to f[q], that its vertex is inside q)
                s = ((fqda + q * q) - (fvda + v[k]*v[k])) / qmvt2;

            } else {
                // right arm of v[k] intersects with left arm of q
                // (if neither vertex is inside the other)
                float sqrt_comp = sqrt(alpha*(beta*qmvsq + f[v[k]] - f[q]) + beta*(f[q] - f[v[k]]));
                float fac_comp = alpha * q - beta * v[k];
                float amb = alpha - beta;

                s = -(sqrt_comp - fac_comp) / amb;
            }

            // s could be to the left of the last "valid from"
            // if that is so,
            if (s <= z[k]) {
                k--; // if wraps the last parabola, get rid of that last one
                if (k < 0) {
                    k = 0;
                    v[k] = q;
                }
                continue;
            } else {
                k++;
                v[k] = q;
                z[k] = s;
                z[k+1] = INFINITY;
                break;
            }
        }
    }

    k = 0; // start at the left
    for (q = 0; q < len; q++) {
        while (z[k+1] < q) {
            k++;
        }
        if (q < v[k]) {
            out[q] = alpha * (q-v[k])*(q-v[k]) + f[v[k]];
        } else {
            out[q] = beta * (q-v[k])*(q-v[k]) + f[v[k]];
        }
    }
}


/*
 * Functions: strided_load, strided_store
 * ----------------------------
 * Copies data into a buffer column-wise.
 *
 *   dest: A 1D buffer to copy into.
 *   src: A 1D buffer to copy from.
 *   num: The number of values to copy.
 *   stride: The width (or height) of the array
 */

static inline void strided_load(float *dest, float *src, size_t num, size_t stride) {
    int i;
    for (i = 0; i < num; i++) {
        dest[i] = src[stride * i];
    }
}

static inline void strided_store(float *dest, float *src, size_t num, size_t stride) {
    int i;
    for (i = 0; i < num; i++) {
        dest[stride * i] = src[i];
    }
}

/*
 * Function: dt2d
 * ----------------------------
 * Goes through all rows and columns to apply the dt1d function to each.
 *
 *   f: A 2D buffer (height x width)
 *   height: The number of rows
 *   width: The number of columns
 *   row_alpha: square of left_factor (see below)
 *   row_beta: square of right_factor (see below)
 *   col_alpha: square of top_factor (see below)
 *   col_beta: square of bottom_factor (see below)
 */

static void dt2d(float *f, int height, int width,
                 float row_alpha, float row_beta, float col_alpha, float col_beta) {

    int longest_dim = (width > height) ? width : height;

    int *vstack = malloc(longest_dim * sizeof(int)); // pre-allocated buffer to store v[k]
    float *zstack = malloc(longest_dim * sizeof(float)); // pre-allocated buffer to store z[k]

    float *line_in = malloc(longest_dim * sizeof(float)); // line_in is as long as the wider dimension.
    float *line_out = malloc(longest_dim * sizeof(float));

    int y, x; // for loop indices

    // transform rows
    void (*row_dt1d_func)(float*, float*, int, float, float, int*, float*) =
        (row_alpha == row_beta ? &dt1d_isotropic : &dt1d_anisotropic);

    for (y = 0; y < height; y++) {
        // move one row of the image into line_in
        memcpy(line_in, f + (y * width), width * sizeof(float));
        // input, output, len, ints, floats
        (*row_dt1d_func)(line_in, f + (y * width), width, row_alpha, row_beta, vstack, zstack);
    }

    // transform cols
    void (*col_dt1d_func)(float*, float*, int, float, float, int*, float*) =
        (col_alpha == col_beta ? &dt1d_isotropic : &dt1d_anisotropic);

    for (x = 0; x < width; x++) {
        strided_load(line_in, f + x, height, width);  // copy one column into line_in
        memcpy(line_out, line_in, height * sizeof(float)); // copy the column into line_out also
        (*col_dt1d_func)(line_in, line_out, height, col_alpha, col_beta, vstack, zstack); // apply the transform to line_in, and overwrite line_out where needed
        strided_store(f + x, line_out, height, width); // write the result back to f.
    }

    // Outputs are squared distances, so we take the square root of everything
    int index = 0;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            index = (y * width) + x;
            f[index] = sqrt(f[index]);
        }
    }

    free(line_in);
    free(line_out);
    free(vstack);
    free(zstack);
}

/*
 * Function: not_floatmatrix
 * ----------------------------
 * Checks whether the input is a 2D image of floats
 *
 *   mat: A 2D buffer (height x width)
 */
int not_floatmatrix(PyArrayObject *mat) { 
    if (mat->descr->type_num != NPY_FLOAT32 || mat->nd != 2) {
        PyErr_SetString(PyExc_ValueError,
                        "In not_floatmatrix: array must be of type Float, and 2-dimensional (n x m).");
        return 1;
    }
    return 0;
}

/*
 * Function: adt
 * ----------------------------
 * The main function – this is called from Python to perform the
 * distance transform on the Numpy array mat.
 * This function expects five PyObjects in *args:
 *
 *   mat: A 2D buffer (height x width)
 *   left_factor: A multiplier that squishes the distance to the left
 *   right_factor: A multiplier that squishes the distance to the right
 *   top_factor: A multiplier that squishes the distance to the top
 *   bottom_factor: A multiplier that squishes the distance to the bottom
 *
 * The output is written back to the same Numpy array.
 */
static PyObject *adt(PyObject *self, PyObject *args) {
    PyArrayObject *mat;
    float left_factor, right_factor, top_factor, bottom_factor;
    if (!PyArg_ParseTuple(args, "O!ffff", &PyArray_Type, &mat,
                          &left_factor, &right_factor,
                          &top_factor, &bottom_factor))
        return NULL;

    if (mat == NULL) return NULL;
    if (not_floatmatrix(mat)) return NULL;

    int h = mat->dimensions[0];
    int w = mat->dimensions[1];

    float *cout = (float *) mat->data;
    dt2d(cout, h, w, SQ(left_factor), SQ(right_factor), SQ(top_factor), SQ(bottom_factor));
    Py_INCREF(Py_None);
    return Py_None;
}
