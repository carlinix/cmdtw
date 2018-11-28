/**
Author: rcarlini@gmail.com
 */
#include "cmdtw.h"
#include <math.h>
#include <stdlib.h>


/* Global Variables Used to work with some distance functions,
  the idea is to preserve the signature of functions */
double *dist_mapped_aux = NULL;  /* When using mapped distances */
int *dist_mapped_len_aux = NULL; /* the lengh of dist_mapped squared matrix -
                                    or the p value for */


/*path patterns for step functions*/
/*
 'Sakoe-Chiba 1/2 sym'
 [5,-1,-3,-1,-2,-1,-1,-2,-1,-3,-1] <- 5 end points (-1-,3),(-1,-2),...
 [3, 0,-1, 0,-2,-1,-3, 0, 0, 0, 0] <- subpath to the first end point
 [2, 0,-1,-1,-2, 0, 0, 0, 0, 0, 0],
 [1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0],
 [2,-1, 0,-2,-1, 0, 0, 0, 0, 0, 0],
 [3,-1, 0,-2, 0,-3,-1, 0, 0, 0, 0],
 */

const int dp1_path_pattern[6][11] = {
   { 2, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0},
    {1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0},
    {1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

const int dp2_path_pattern[6][11] = {
    { 3, 0,-1,-1, 0,-1,-1, 0, 0, 0, 0},
    {1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0},
    {1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

const int p1div2_path_pattern[6][11] = {
    {5,-1,-3,-1,-2,-1,-1,-2,-1,-3,-1},
    {3, 0,-1, 0,-2,-1,-3, 0, 0, 0, 0},
    {2, 0,-1,-1,-2, 0, 0, 0, 0, 0, 0},
    {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
    {2,-1, 0,-2,-1, 0, 0, 0, 0, 0, 0},
    {3,-1, 0,-2, 0,-3,-1, 0, 0, 0, 0}};

const int p1_path_pattern[6][11] =  {
    {3,-1,-2,-1,-1,-2,-1, 0, 0, 0, 0},
    {2, 0,-1,-1,-2, 0, 0, 0, 0, 0 ,0},
    {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
    {2,-1, 0,-2,-1, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

const int p2_path_pattern[6][11] = {
    {3,-2,-3,-1,-1,-3,-2, 0, 0, 0, 0},
    {3, 0,-1,-1,-2,-2,-3, 0, 0, 0 ,0},
    {1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0},
    {3,-1, 0,-2,-1,-3,-2, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};


/*******************************************************************************
 * Helper functions
 ******************************************************************************/
/**
 * 
 * @param matrix
 * @param M
 * @param N
 */
void _print_matrix(double *matrix, int M, int N) {
    register int i, j ;
    printf("Allocated Cost Matrix C[%d][%d]:\n",M,N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++)
            printf("%.2lf\t\t",matrix[idx(i, j, N)]);
        printf("\n");       
    }    
}

/**
 * @brief  Finds the minimum of 2 doubles ant its position
 * @param a a double
 * @param b b double
 * @return min(a,b) and its position
 *
 */
 struct t_item min2idx(double a, double b) {
    struct t_item item;    
    if (a < b) { 
        item.val = a; 
        item.idx = 0;
    }   
    else { 
        item.val = b; 
        item.idx = 1  ;
    }
    return item;
}

/**
 * @brief  Finds the minimum of 3 doubles 
 * @param x  x 
 * @param y  y
 * @param z z
 * @return min(x,y,z)
 *
 */
 double min3(double x, double y, double z) {
    if (x < y && x < z)
        return x;
    else if (y < x && y < z)
        return y;
    else
        return z;
}

/**
 * @brief  Finds the minimum of 3 doubles and its position (0, 1 or 2)
 * @param x
 * @param y
 * @param z
 * @return struct t_item, t_item = { min(x,y,z), position }
 */
 struct t_item min3idx(double x, double y, double z) {/*0, 1, 2*/
    struct t_item item;
    if (x < y && x < z) {
        item.val = x;
        item.idx = 0;
    } else if (y < x && y < z) {
        item.val = y;
        item.idx = 1;
    } else {
        item.val = z;
        item.idx = 2;
    }
    return item;
}

/**
 * @brief Finds the minimum of n doubles
 * @param arr array
 * @param n length of the array
 * @return  min(arr)
 */
double min_n(double* arr, int n) {

    double min = arr[0];
    register int i = 1;
    for(; i < n; i++) {
        if(min > arr[i])
            min = arr[i];
    }
    return min;
}

/**
 * @brief Finds the minimum of n doubles and its position
 * @param arr array
 * @param n length of the array
 * @return  min(arr) and its position
 */
struct t_item min_nidx(double* arr, int n) {
    struct t_item item;
    item.val = arr[0];;
    item.idx = 0;
    register int i = 1;
    for(; i < n; i++) {
        if(item.val > arr[i]) {
            item.val = arr[i];
            item.idx = i;
        }
    }
    return item;
}

/**
 * @bref Rounds a number (inline function)
 * @param number
 * @return number + 1 or number - 1
 */
 extern double _round(double number);

/*******************************************************************************
 * Distance functions
 ******************************************************************************/

/**
 * @brief  Multidimensional CityBlock distance
 * @param a multidimensional point
 * @param b  multidimensional point
 * @param dim number of dimensions
 * 
 * @return  Manhattan distance between two multidimensional points
 *
  """
    Computes the City Block (Manhattan) distance.

    Computes the Manhattan distance between two 1-D arrays `u` and `v`,
    which is defined as

    .. math::

       \\sum_i {\\left| u_i - v_i \\right|}.

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    cityblock : double
        The City Block (Manhattan) distance between vectors `u` and `v`.

    """
 */
double cityblock(double* a, double* b, int dim) {
    register int i = 0;
    double sum = 0.0;
    for (; i < dim; i++)
        sum += fabs(a[i] - b[i]);
    return sum;
}

/**
 * @brief  Multidimensional Euclidian distance
 * @param a multidimensional point
 * @param b  multidimensional point
 * @param dim number of dimensions
 * 
 * @return  Euclidian distance between two multidimensional points

   Computes the Euclidean distance between two 1-D arrays.

    The Euclidean distance between 1-D arrays `u` and `v`, is defined as

    .. math::

       {||u-v||}_2

    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    euclidean : double
        The Euclidean distance between vectors `u` and `v`.

 * 
 */
double euclidean(double* a, double* b, int dim) {
    register int i = 0;
    double sum = 0.0;
    for (; i < dim; i++)
        sum += pow((a[i] - b[i]), 2);
    return sqrt(sum);
}


/**
 * @brief  Multidimensional Euclidian distance
 * @param a multidimensional point
 * @param b  multidimensional point
 * @param dim number of dimensions
 *
 * @return  Euclidian distance between two multidimensional points
 *

   """
    Computes the squared Euclidean distance between two 1-D arrays.

    The squared Euclidean distance between `u` and `v` is defined as

    .. math::

       {||u-v||}_2^2.


    Parameters
    ----------
    u : (N,) array_like
        Input array.
    v : (N,) array_like
        Input array.

    Returns
    -------
    sqeuclidean : double
        The squared Euclidean distance between vectors `u` and `v`.

 */
double sqeuclidean(double* a, double* b, int dim) {
    register int i = 0;
    double sum = 0.0;
    for (; i < dim; i++)
        sum += pow((a[i] - b[i]), 2);
    return sum;
}

/**
 * @brief
 * @param a
 * @param b
 * @param dim
 * @return 
 */
double indexdist(double *a, double *b, int dim) {

        return (!dist_mapped_len_aux || *dist_mapped_len_aux < 0) ? INFINITY :
        dist_mapped_aux[(int)idx(_round(a[0]), _round(b[0]), *dist_mapped_len_aux)];
}

/**
 * @brief Chooses right distance functio
 * @param dist_type dtw_settings.dist_type
 * 
 * @return pointer to a distance function
 * 
 */
dist_fptr choose_dist(struct t_dtw_settings t_s) {
    switch (t_s.dist_type) {
        case(_EUCLID):
            return &euclidean;
        case(_EUCLID_SQUARED):
            return &sqeuclidean;
        case(_CITYBLOCK):
            return &cityblock;
        case(_INDEXMATRIX):
            dist_mapped_aux = t_s.dist_mapped;
            dist_mapped_len_aux = t_s.dist_mapped_len;
            return &indexdist;
        default:
            return &cityblock;
    }
}


/**
 * @brief Execute a function distance between two matrices
 * @param a - pointer to matrix A
 * @param b - pointer to matrix B
 * @param r  - pointer to matrix C (dest)
 * @param ra - number of rows A
 * @param rb - number of rows B
 * @param dim - number of columns in A and B
 * @param type - type of distance
 * @return 0 ok
 *
 */
int m_dist(double *a, double *b, double *r,  int ra, int rb, int dim, int type) {
    register int i = 0;
    register int j = 0;

    struct t_dtw_settings t_s;
    t_s.dist_type = type;
    dist_fptr dist = NULL;

    if (type == _INDEXMATRIX) {
        return 2;
    }

    dist = choose_dist(t_s);
    if (dim < 1)
        return 1;

    for (i = 0; i < ra; i++)
        for(j = 0; j < rb; j++)
            r[idx(i,j,rb)] = dist((a+(i*dim)), (b+(j*dim)), dim);

    return 0;
}


/*******************************************************************************
 * Window functions
 ******************************************************************************/

/**
 * @brief Don't use a global constraint (window)
 * @param i not used
 * @param j not used
 * @param k not used
 * @param I not used
 * @param J not used
 * 
 * @return always return true for every input value
 * 
 */
bool nowindow(int i, int j, double k, double I, double J) {
    return true;
}

/**
 * @brief Itakura global constraints
 * @param i row index
 * @param j column index
 * @param k k = 1/s = I/J = len_ref/len_query
 * @param I length of reference sequence, len_ref
 * @param J length of query sequence, len_query
 * @return boolean value, if cost[i][j] is legal or not
 */
bool itakura(int i, int j, double k, double I, double J) {
    return  j <  2*i                &&  
            i <= 2*j                &&  
            i >= I-1 -2*(J-j)       &&  
            j >  J-1 -2*(I-i);
}

/**
 * @brief Sakoe-Chiba band global constraint
 * @param i  row index
 * @param j column index
 * @param p width of the window
 * @param I not used
 * @param J not used
 * @return boolean value, if cost[i][j] is legal or not
 */
bool scband(int i, int j, double p, double I, double J) {
    int r = 0;
    //r = (p < 1) ?  max2(I,J) * 0.2 : p;
    r = max2(p,(fabs(I - J) + 2));
   // printf("R %f -> %d\n",p,r);
    return abs(i - j) < r;
}

/**
 * @brief Palival global constraint, it is similar to scband, but adapts to the
 * length of sequences
 * @param i row index
 * @param j column index
 * @param r width of the window
 * @param k k = 1/s = I/J = len_ref/len_query
 * @param I length of reference sequence, len_ref
 * @param J length of query sequence, len_query
 * @return boolean value, if cost[i][j] is legal or not
 */
bool palival(int i, int j, double r, double I, double J) { 
    return fabs(i*J/I - j) < r; 
} 

/**
 * 
 * @brief Palival global constraint, it is similar to scband, but adapts to the
 * length of sequences
 * @param i row index
 * @param j column index
 * @param r width of the window
 * @param k k = 1/s = I/J = len_ref/len_query
 * @param I length of reference sequence, len_ref
 * @param J length of query sequence, len_query
 * @return boolean value, if cost[i][j] is legal or not
 */
bool palival_mod(int i, int j, double p, double I, double J) {
    double k = I/J;
    double r = p*I;
    return i > k*j - r && i < k*j + r;
}

/**
 *   FIXME: NOTE: function is partly redundant at the moment, it always returns  itakura
 *  Chooses right window function
 *  int dp_type: dtw_settings.win , step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 * 
 * @param dtw_settings
 * @return 
 */
window_fptr choose_window(struct t_dtw_settings *dtw_settings) {
    switch (dtw_settings->window_type) {
        case(_SCBAND):
            return &scband;
        case(_PALIVAL):
            return &palival;
        case(_ITAKURA):
            return &itakura;
        case(_PALIVAL_MOD):
            return &palival_mod;
        default:
            return &nowindow;
    }
}

/**
 * @bref Computes extra size for the cost matrix (dtw_settings.offset)
 * @param dp_type dtw_settings.dp_type, step pattern type
 * @return int (1,2 or 3)
 * 
 */
int extra_size(int dp_type) {
    switch (dp_type) {
        case(_DP1):
        case(_DP2):
        case(_DP3):
        case(_SCP0SYM):
        case(_SCP0ASYM):
            return 1;
        case(_SCP1SYM):
        case(_SCP1ASYM):
            return 2;
        default:
            return 3;
    }
}

/**
 * @brief Computes and return parameter for the window function
 * @param dtw_settings structure with dtw settings
 * @param len_ref length of the reference sequence
 * @param len_query length of the query sequence
 * 
 * @return  window parameter
 * 
 */
double choose_window_param(struct t_dtw_settings *dtw_settings, int len_ref,
        int len_query) {
    double p = 0;
    if (dtw_settings->window_type == _ITAKURA)
        p = len_ref / (double) len_query;
    else
        p = dtw_settings->window_param;
    return p;
}

/*******************************************************************************
 * Step functions
 ******************************************************************************/

/**
 * NOTE: ALL STEP FUNCTIONS uses the same args defined in macro  _DP_ARGS. 
 * Step functions like step_pattern_type and step_pattern_typedir are pretty 
 * similar, step_pattern_type are used in computing dtw without path (without
 * traceback). step_pattern_typedir are used in computing dtw with path 
 * (traceback) 
 */

/* --------------------- Sakoe-Chiba classifications --------------------- */
/**
 * @brief Sakoe-Chiba classification p = 0, symmetric step pattern
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s 
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double
 * 
 * @return  double, value to assign cost_matrix[i][j] 
 * 
 * This function is a alias for dp1
 *  Step pattern dp1 - :
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   2*d(r[i],q[j])
 *     ) 
 */
double p0_sym(_DP_ARGS) {
    return dp1(ref, query, dim, cost_matrix, i, j, t_s, size2, dist);
}

/**
 * @brief Sakoe-Chiba classification p = 0, symmetric step pattern
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s 
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double
 * 
 * @return  t_item, value to assign cost_matrix[i][j] and its position
 * 
 * This function is a alias for dp1dir
 *  Step pattern dp1dir - :
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   2*d(r[i],q[j])
 *     ) 
 */
struct t_item p0_symdir(_DP_ARGS) {
    return dp1dir(ref, query, dim, cost_matrix, i, j, t_s, size2, dist);
}

/**
 * @brief  Sakoe-Chiba classification p = 0, asymmetric step pattern:
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s 
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double
 * @return value to assign cost_matrix[i][j] 
 * 
 *  Step pattern:
 *  min(      
 *      cost_matrix[i][j-1]   +  0    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   d(r[i],q[j])
 *     ) 
 */
double p0_asym(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3(   cost_matrix[idx(i, j-1, size2)],         //0
                   cost_matrix[idx(i-1, j, size2)]   + d,   //1
                   cost_matrix[idx(i-1, j-1, size2)] + d);  //2
}

/**
 * @brief  Sakoe-Chiba classification p = 0, asymmetric step pattern:
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s 
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double
 * @return t_item, value to assign cost_matrix[i][j] and its position
 * 
 *  Step pattern:
 *  min(      
 *      cost_matrix[i][j-1]   +  0    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   d(r[i],q[j])
 *     ) 
 */
struct t_item p0_asymdir(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3idx(cost_matrix[idx(i, j-1, size2)],         //0
                   cost_matrix[idx(i-1, j, size2)]   + d,   //1
                   cost_matrix[idx(i-1, j-1, size2)] + d);  //2
}

/**
 * @brief Sakoe-Chiba classification p = 1, symmetric step pattern
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j] 
 * 
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 *  min(
 *      cost_matrix[i-1][j-2] + 2d(r[i],q[j-1]) + d(r[i],q[j]),
 *      cost_matrix[i-1][j-1] + 2d(r[i],q[j]),
 *      cost_matrix[i-2][j-1] + 2d(r[i-1],q[j]) + d(r[i],q[j]),
       )
 */
double p1_sym(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0 , dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);

    return min3(cost_matrix[idx(i-1, j-2, size2)] + 2*d01 + d00,
                cost_matrix[idx(i-1, j-1, size2)] + 2*d00,
                cost_matrix[idx(i-2, j-1, size2)] + 2*d10 + d00);
}

/**
 * @brief Sakoe-Chiba classification p = 1, symmetric step pattern
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return t_item, value to assign cost_matrix[i][j] and its position
 * 
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 *  min(
 *      cost_matrix[i-1][j-2] + 2d(r[i],q[j-1]) + d(r[i],q[j]),
 *      cost_matrix[i-1][j-1] + 2d(r[i],q[j]),
 *      cost_matrix[i-2][j-1] + 2d(r[i-1],q[j]) + d(r[i],q[j]),
       )
 */
struct t_item p1_symdir(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0 ,dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset,0,dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset), 0), 0, dim)],
                               &query[idx(j - t_s->offset, 0, dim)], dim);

    return min3idx(cost_matrix[idx(i-1, j-2, size2)] + 2*d01 + d00,
                cost_matrix[idx(i-1, j-1, size2)] + 2*d00,
                cost_matrix[idx(i-2, j-1, size2)] + 2*d10 + d00);
}

/**
 * @brief Sakoe-Chiba classification p = 1, symmetric step pattern
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return t_item, value to assign cost_matrix[i][j] and its position
 *
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 *  min(
 *      cost_matrix[i-1][j-2] + (d(r[i],q[j-1]) + d(r[i],q[j]))/2,
 *      cost_matrix[i-1][j-1] + d(r[i],q[j]),
 *      cost_matrix[i-2][j-1] + d(r[i-1],q[j]) + d(r[i],q[j]),
 )
 */
struct t_item p1_asymdir(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    
    return min3idx(cost_matrix[idx(i-1, j-2, size2)] + (d01 + d00)/2,
                cost_matrix[idx(i-1, j-1, size2)] + d00,
                cost_matrix[idx(i-2, j-1, size2)] + d10 + d00);
}

/**
 * @brief Sakoe-Chiba classification p = 1, asymmetric step pattern
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j]
 *
 * Sakoe-Chiba classification p = 1, symmetric step pattern:
 *  min(
 *      cost_matrix[i-1][j-2] + (d(r[i],q[j-1]) + d(r[i],q[j]))/2,
 *      cost_matrix[i-1][j-1] + d(r[i],q[j]),
 *      cost_matrix[i-2][j-1] + d(r[i-1],q[j]) + d(r[i],q[j]),
 )
 */
double p1_asym(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0 , dim)], dim);
    
    return min3(cost_matrix[idx(i-1, j-2, size2)] + (d01 + d00)/2,
                cost_matrix[idx(i-1, j-1, size2)] + d00,
                cost_matrix[idx(i-2, j-1, size2)] + d10 + d00);
}

/**
 * @brief Sakoe-Chiba classification p = 2, symmetric step pattern:
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j]
 *
 * min(
 *     cost_matrix[i-2][j-3] + 2d(r[i-1],q[j-2]) +
 *                             2d(r[i],q[j-1])   +
 *                             d(r[i],q[j]),
 *
 *     cost_matrix[i-1][j-1] + 2d(r[i],q[j]),
 *
 *     cost_matrix[i-3][j-2] + 2d(r[i-2],q[j-1]) +
 *                             2d(r[i-1],q[j])   +
 *                             d(r[i],q[j])
 *     )
 */
 struct t_item p2_symdir(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    
    double d12 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    
    double d21 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    
    return min3idx(cost_matrix[idx(i-2, j-3, size2)] + 2*d12 + 2*d01 + d00,
                cost_matrix[idx(i-1, j-1, size2)] + 2*d00,
                cost_matrix[idx(i-3, j-2, size2)] + 2*d21 + 2*d10 + d00);
}
 
/**
 * @brief Sakoe-Chiba classification p = 2, symmetric step pattern:
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j]
 *
 * min(
 *     cost_matrix[i-2][j-3] + 2d(r[i-1],q[j-2]) +
 *                             2d(r[i],q[j-1])   +
 *                             d(r[i],q[j]),
 *
 *     cost_matrix[i-1][j-1] + 2d(r[i],q[j]),
 *
 *     cost_matrix[i-3][j-2] + 2d(r[i-2],q[j-1]) +
 *                             2d(r[i-1],q[j])   +
 *                             d(r[i],q[j])
 *     )
 */
 double p2_sym(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    
    double d12 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    
    double d21 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    
    return min3(cost_matrix[idx(i-2, j-3, size2)] + 2*d12 + 2*d01 + d00,
                cost_matrix[idx(i-1, j-1, size2)] + 2*d00,
                cost_matrix[idx(i-3, j-2, size2)] + 2*d21 + 2*d10 + d00);
}
 
/**
 * @brief Sakoe-Chiba classification p = 2, asymmetric step pattern:
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j]
 *
 *  min(     
 *     cost_matrix[i-2][j-3] + 2( d(r[i-1],q[j-2]) +
 *                                d(r[i],q[j-1])   +
 *                                d(r[i],q[j]) ),
 *
 *     cost_matrix[i-1][j-1] + d(r[i],q[j]),
 *
 *     cost_matrix[i-3][j-2] + d(r[i-2],q[j-1]) +
 *                             d(r[i-1],q[j])   +
 *                             d(r[i],q[j])
 *     )

 */
 struct t_item p2_asymdir(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0 , dim)],
                      &query[idx(j - t_s->offset, 0 , dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d12 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    double d21 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    return min3idx(cost_matrix[idx(i-2, j-3, size2)] + 2.0*(d12 + d01 + d00)/3.0,
                cost_matrix[idx(i-1, j-1, size2)] + d00,
                cost_matrix[idx(i-3, j-2, size2)] + d21 + d10 + d00);
}
 
 /**
 * @brief Sakoe-Chiba classification p = 2, asymmetric step pattern:
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j]
 *
 *  min(     
 *     cost_matrix[i-2][j-3] + 2( d(r[i-1],q[j-2]) +
 *                                d(r[i],q[j-1])   +
 *                                d(r[i],q[j]) ),
 *
 *     cost_matrix[i-1][j-1] + d(r[i],q[j]),
 *
 *     cost_matrix[i-3][j-2] + d(r[i-2],q[j-1]) +
 *                             d(r[i-1],q[j])   +
 *                             d(r[i],q[j])
 *     )
 */
 double p2_asym(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    
    double d12 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    
    double d21 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    

    return min3(cost_matrix[idx(i-2, j-3, size2)] + 2.0*(d12 + d01 + d00)/3.0,
                cost_matrix[idx(i-1, j-1, size2)] + d00,
                cost_matrix[idx(i-3, j-2, size2)] + d21 + d10 + d00);
}
 
 /**
 * @brief Sakoe-Chiba classification p = 0.5, symmetric step pattern
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j]
 */
struct t_item p1div2_symdir(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d02 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d20 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double arr[5] = {
        cost_matrix[idx(i-1, j-3, size2)] + 2.0*d02 + d01 + d00,
        cost_matrix[idx(i-1, j-2, size2)] + 2.0*d01 + d00,
        cost_matrix[idx(i-1, j-1, size2)] + 2.0*d00,
        cost_matrix[idx(i-2, j-1, size2)] + 2.0*d10 + d00,
        cost_matrix[idx(i-3, j-1, size2)] + 2.0*d20 + d10 + d00};
    return min_nidx(arr,5);
}
 
/**
 * @brief Sakoe-Chiba classification p = 0.5, symmetric step pattern
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return value to assign cost_matrix[i][j]
 */
double p1div2_sym(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d02 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d20 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    
    double arr[5] = {
        cost_matrix[idx(i-1, j-3, size2)] + 2.0*d02 + d01 + d00,
        cost_matrix[idx(i-1, j-2, size2)] + 2.0*d01 + d00,
        cost_matrix[idx(i-1, j-1, size2)] + 2.0*d00,
        cost_matrix[idx(i-2, j-1, size2)] + 2.0*d10 + d00,
        cost_matrix[idx(i-3, j-1, size2)] + 2.0*d20 + d10 + d00 };
    return min_n(arr,5);
}

/**
 * @brief Sakoe-Chiba classification p = 0.5, asymmetric step pattern:
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return 
 */
struct t_item p1div2_asymdir(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d02 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d20 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
   
    double arr[5] = {cost_matrix[idx(i-1, j-3, size2)] + (d02 + d01 + d00)/3.0,
                     cost_matrix[idx(i-1, j-2, size2)] + (d01 + d00)/2,
                     cost_matrix[idx(i-1, j-1, size2)] + d00,
                     cost_matrix[idx(i-2, j-1, size2)] + d10 + d00,
                     cost_matrix[idx(i-3, j-1, size2)] + d20 + d10 + d00};
    return min_nidx(arr,5);
}

/**
 * 
 * @param ref
 * @param query
 * @param dim
 * @param cost_matrix
 * @param i
 * @param j
 * @param t_s
 * @param size2
 * @param dist
 * @return 
 */
double p1div2_asym(_DP_ARGS) {
    double d00 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d01 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-1 - t_s->offset),0), 0, dim)], dim);
    double d02 = dist(&ref[idx(i - t_s->offset, 0, dim)],
                      &query[idx(max2((j-2 - t_s->offset),0), 0, dim)], dim);
    double d10 = dist(&ref[idx(max2((i-1 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
    double d20 = dist(&ref[idx(max2((i-2 - t_s->offset),0), 0, dim)],
                      &query[idx(j - t_s->offset, 0, dim)], dim);
   
    double arr[5] = {cost_matrix[idx(i-1, j-3, size2)] + (d02 + d01 + d00)/3.0,
                     cost_matrix[idx(i-1, j-2, size2)] + (d01 + d00)/2,
                     cost_matrix[idx(i-1, j-1, size2)] + d00,
                     cost_matrix[idx(i-2, j-1, size2)] + d10 + d00,
                     cost_matrix[idx(i-3, j-1, size2)] + d20 + d10 + d00};
    return min_n(arr,5);
}

/**
 * @brief Step patern dpw - weights
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s 
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return t_item, value to assign cost_matrix[i][j] and pos
 * 
 *  min(      
 *      cost_matrix[i][j-1]   +   a*d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   b*d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   c*d(r[i],q[j])
 *     )
 * where a,b,c are weights
 */
struct t_item dpwdir(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + t_s->weights.a * d,
            cost_matrix[idx(i - 1, j, size2)] + t_s->weights.b * d,
            cost_matrix[idx(i - 1, j - 1, size2)] + t_s->weights.c * d);
}

/**
 * 
 * @brief Step patern dpw - weights
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return double, value to assign cost_matrix[i][j]
 * 
 *  min(      
 *      cost_matrix[i][j-1]   +   a*d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   b*d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   c*d(r[i],q[j])
 *     )
 * where a,b,c are weights
 */
double dpw(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3(cost_matrix[idx(i, j - 1, size2)] + t_s->weights.a * d,
            cost_matrix[idx(i - 1, j, size2)] + t_s->weights.b * d,
            cost_matrix[idx(i - 1, j - 1, size2)] + t_s->weights.c * d);
}

/**
 * @brief Step pattern dp1dir - 
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s not used
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return t_item, value to assign cost_matrix[i][j] and pos
 * 
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   2*d(r[i],q[j])
 *     )
 */
struct t_item dp1dir(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + d,
            cost_matrix[idx(i - 1, j, size2)] + d,
            cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d);
}

/**
 * 
 * @brief Step pattern dp1 - 
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return double, value to assign cost_matrix[i][j]
 * 
 *  Step pattern dp1 - :
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   2*d(r[i],q[j])
 *     ) 
 */
double dp1(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3(cost_matrix[idx(i, j - 1, size2)] + d,
            cost_matrix[idx(i - 1, j, size2)] + d,
            cost_matrix[idx(i - 1, j - 1, size2)] + 2 * d);
}

/**
 * @brief Step pattern dp2dir - 
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s not used
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return t_item, value to assign cost_matrix[i][j] and pos
 * 
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   d(r[i],q[j])
 *     )
 */
struct t_item dp2dir(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3idx(cost_matrix[idx(i, j - 1, size2)] + d,
            cost_matrix[idx(i - 1, j, size2)] + d,
            cost_matrix[idx(i - 1, j - 1, size2)] + d);
}

/**
 * 
 * @brief Step pattern dp2 - 
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return double, value to assign cost_matrix[i][j]
 * 
 *  Step pattern dp2 - :
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j]),    
 *      cost_matrix[i-1][j-1] +   d(r[i],q[j])
 *     )
 */
double dp2(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min3(cost_matrix[idx(i, j - 1, size2)] + d,
            cost_matrix[idx(i - 1, j, size2)] + d,
            cost_matrix[idx(i - 1, j - 1, size2)] + d);
}

/**
 * @brief Step pattern dp3dir - 
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s not used
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return t_item, value to assign cost_matrix[i][j] and pos
 * 
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j])
 *     )
 */
struct t_item dp3dir(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min2idx(cost_matrix[idx(i, j - 1, size2)] + d,
            cost_matrix[idx(i - 1, j, size2)] + d);
}

/**
 * 
 * @brief Step pattern dp3 - 
 * @param ref reference sequence
 * @param query query sequence
 * @param dim dimension of ref/query points
 * @param cost_matrix cost matrix
 * @param i index of the cost matrix
 * @param j column index of the cost matrix
 * @param t_s
 * @param size2 cost matrix columns count 
 * @param dist pointer to distance function double

 * @return double, value to assign cost_matrix[i][j]
 * 
 *  Step pattern dp3 - :
 *  min(      
 *      cost_matrix[i][j-1]   +   d(r[i],q[j])    
 *      cost_matrix[i-1][j]   +   d(r[i],q[j])
 *     )
 */
double dp3(_DP_ARGS) {
    double d = dist(&ref[idx(i - t_s->offset, 0, dim)],
                    &query[idx(j - t_s->offset, 0, dim)], dim);
    return min2(cost_matrix[idx(i, j - 1, size2)] + d,
            cost_matrix[idx(i - 1, j, size2)] + d);
}

/**
 * @brief Chooses right step function(without traceback)
 * @param dp_type dtw_settings.dp_type, step function type
 * 
 * @return pointer to a step function without traceback
 */
dp_fptr choose_dp(int dp_type) {
    switch (dp_type) {
        case(_DPW):
            return &dpw;
        case(_DP1):
            return &dp1;
        case(_DP2):
            return &dp2;
        case(_DP3):
            return &dp3;
        case(_SCP0SYM):
            return &p0_sym;
        case(_SCP0ASYM):
            return &p0_asym;
        case(_SCP1DIV2SYM):
            return &p1div2_sym;
        case(_SCP1DIV2ASYM):
            return &p1div2_asym;
        case(_SCP1SYM):
            return &p1_sym;
        case(_SCP1ASYM):
            return &p1_asym;
        case(_SCP2SYM):
            return &p2_sym;
        case(_SCP2ASYM):
            return &p2_asym;
        default:
            return &dp2;
    }
}

/**
 * @brief Chooses right step function(with traceback)
 * @param dp_type dtw_settings.dp_type, step function type
 * 
 * @return pointer to a step function without traceback
 *
 */
dpdir_fptr choose_dpdir(int dp_type) {
    switch (dp_type) {
        case(_DPW):
            return &dpwdir;
        case(_DP1):
            return &dp1dir;
        case(_DP2):
            return &dp2dir;
        case(_DP3):
            return &dp3dir;
        case(_SCP0SYM):
            return &p0_symdir;
        case(_SCP0ASYM):
            return &p0_asymdir;
        case(_SCP1DIV2SYM):
            return &p1div2_symdir;
        case(_SCP1DIV2ASYM):
            return &p1div2_asymdir;
        case(_SCP1SYM):
            return &p1_symdir;
        case(_SCP1ASYM):
            return &p1_asymdir;
        case(_SCP2SYM):
            return &p2_symdir;
        case(_SCP2ASYM):
            return &p2_asymdir;
        default:
            return &dp2dir;
    }
}

/*******************************************************************************
 * Path functions
 ******************************************************************************/
 
/**
 * @brief  Compute full path from path_points and path_pattern.
 * @param pattern path pattern path pattern
 * @param len_ref length of the reference sequence
 * @param len_query length of the query sequence
 * @param path_points path points as array of directions, where direction is 
 * idx returned by dp_dir functions
 * @param path_points_count length of path points = len_ref + len_query
 * @param path  path array
 * 
 * @return length of the constructed path
 * 
 *  Compute full path from path_points and path_pattern. This function is 
 *  necessary for the multi-step step functions (like p1div2_sym), and also 
 *  to transform "directions" to coordinates(matrix indices).
 */
int create_path_from_pattern(const int pattern[6][11], int len_ref, 
        int len_query, int* path_points, int path_points_count,
        struct t_path_element* path) {
    
    int path_idx = 1;
    int i = 0;
    int j = 0;
    path[0].i = len_ref - 1;
    path[0].j = len_query - 1;

    for (; i < path_points_count; i++) {
        int path_idx_tmp = path_idx;
        for (j = 1; j < 2 * pattern[path_points[i] + 1][0] + 1; j += 2) {
            path[path_idx].i = path[path_idx_tmp - 1].i +
                    pattern[path_points[i] + 1][j];
            path[path_idx].j = path[path_idx_tmp - 1].j +
                    pattern[path_points[i] + 1][j + 1];
            path_idx++;
        }
    }
    return path_idx;
}

/**
 * @brief Chooses right path_patterns(2d array) based on current step function
 * @param dtw_settings structure with dtw settings
 *
 * @return path_pattern
 * 
 */
const int (*choose_path_pattern(struct t_dtw_settings dtw_settings))[11] {
    switch (dtw_settings.dp_type) {
        case(_DP1):
        case(_DP2):
        case(_SCP0ASYM):
            return dp2_path_pattern;
        case(_DP3):
            return dp1_path_pattern;
        case(_SCP1DIV2SYM):
        case(_SCP1DIV2ASYM):
            return p1div2_path_pattern;
        case(_SCP1SYM):
        case(_SCP1ASYM):
            return p1_path_pattern;
        case(_SCP2SYM):
        case(_SCP2ASYM):
            return p2_path_pattern;
        default:
            return dp2_path_pattern;
    }
}

/**
 *
 * @brief Computes path(direction array) from directions matrix
 * @param dir_matrix directions matrix
 * @param path_points  direction array(path points)
 * @param len_ref length of the reference sequence
 * @param len_query length of the query sequence
 * @param pattern :path pattern
 * 
 * @return path_points length
 * 
 */
int direct_matrix_to_path_points(int* dir_matrix, int *path_points, int len_ref,
        int len_query, const int pattern[6][11]) {

    long tin_idx = len_ref * len_query - 1;
    int tout_idx = 0;

    register int i = 0;
    register int j = 0;

    path_points[tout_idx] = dir_matrix[tin_idx];

    for (; tin_idx >= 1;) {
        i = pattern[0][2 * dir_matrix[tin_idx] + 1];
        j = pattern[0][2 * dir_matrix[tin_idx] + 2];
        tin_idx += i * len_query + j;
        tout_idx++;
        if (tin_idx < 1)
            break;
        path_points[tout_idx] = dir_matrix[tin_idx];
    }
    return tout_idx;
}

/**
 * @brief Prepares cost matrix. Set extra rows and columns to INFINITY.
 * @param matrix pointer to cost matrix
 * @param len_ref length of the query sequence
 * @param len_query length of the reference sequence
 * @param dtw_settings structure with dtw settings
 * 
 * Prepares cost matrix. Set extra rows and columns to INFINITY. Set all matrix
 * to INFINITY, if there is a window
 */
void fill_matrix(double *matrix, int len_ref, int len_query,
        struct t_dtw_settings dtw_settings) {
    register int M = len_ref + dtw_settings.offset;
    register int N = len_query + dtw_settings.offset;
    register int i = 0;
    register int j = 0;

    /* if there is a window or complicated step pattern */
    if (dtw_settings.window_type != _NOWINDOW ||
            (dtw_settings.dp_type != _DP1 && dtw_settings.dp_type != _DP2 &&
            dtw_settings.dp_type != _DP3 && dtw_settings.dp_type != _SCP0SYM)) {
                
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                matrix[idx(i, j, N)] = INFINITY;
    } 
    else {
        for (i = 0; i < M; i++)
            for (j = 0; j < dtw_settings.offset; j++)
                matrix[idx(i, j, N)] = INFINITY;

        for (i = 0; i < dtw_settings.offset; i++)
            for (j = 0; j < N; j++)
                matrix[idx(i, j, N)] = INFINITY;
    }
    matrix[0] = 0.0;
#ifdef DEBUG
   _print_matrix(matrix, M, N);
#endif
}

/**
 * @brief Dynamic Time Warping algorithm(with traceback)
 * @param ref reference sequence
 * @param query query sequence
 * @param dim number of dimensions of the sequences
 * @param len_ref length of the reference sequence
 * @param len_query length of the query sequence
 * @param dist pointer to a distance function
 * @param dp_dir pointer to a step function
 * @param window  pointer to a window function
 * @param p window parameter
 * @param cost_matrix  cost matrix
 * @param dir_matrix direction matrix
 * @param dtw_settings structure with dtw settings  
 * @return distance between reference and query sequences
 */
double cdtwpath(double* ref, double* query, int dim, int len_ref,
        int len_query, dist_fptr dist, dpdir_fptr dp_dir, window_fptr window,
        double p, double *cost_matrix, int *dir_matrix, 
        struct t_dtw_settings dtw_settings) {

    int off = dtw_settings.offset;
    struct t_item item = {0, 0};
    /* extending matrix */
    int M = len_ref + off;
    int N = len_query + off;
    int i = 0;
    int j = 0;
    double w = 0;
    double s = 0;
    
    bool fast_glob = (dtw_settings.window_type == _PALIVAL ||
            dtw_settings.window_type == _PALIVAL_MOD);
    /* no window or fast window case */
    if (fast_glob || dtw_settings.window_type == _NOWINDOW) {
        if (fast_glob) {
            if (dtw_settings.window_type == _PALIVAL_MOD)
                w = dtw_settings.window_param * (double)len_ref;
            else
                w = dtw_settings.window_param;
            s = len_query / (double)len_ref;
        } else {
            w = INFINITY;
            s = 1;
        }
        // off - off
        cost_matrix[idx(off, off, N)] = dist(&ref[idx(off-off, 0, dim)],
                                             &query[idx(off-off, 0, dim)],
                                             dim);
        for (j = max2(off + 1, _round(s * (off - w))); 
                j < min2(N, _round(s * (off + w) + 1)); j++) {
            item = dp_dir(ref, query, dim, cost_matrix, off, j, &dtw_settings,
                    N, dist);
            cost_matrix[idx(off, j, N)] = item.val;
            dir_matrix[idx(0, j - off, N - off)] = item.idx;
        }

        for (i = off + 1; i < M; i++) {
            for (j = max2(off, _round(s * (i - w))); 
                    j < min2(N, _round(s * (i + w) + 1)); j++) {
                item = dp_dir(ref, query, dim, cost_matrix, i, j, &dtw_settings,
                        N, dist);
                cost_matrix[idx(i, j, N)] = item.val;
                dir_matrix[idx(i - off, j - off, N - off)] = item.idx;
            }
        }
    } /* slow window case */
    else {
        // off - off
        cost_matrix[idx(off, off, N)] = dist(&ref[idx(off-off, 0, dim)],
                                             &query[idx(off-off, 0, dim)],
                                             dim);
        for (j = off + 1; j < N; j++) {
            if (window(off, j, p, len_ref, len_query)) {
                item = dp_dir(ref, query, dim, cost_matrix, off, j, 
                        &dtw_settings, N, dist);
                cost_matrix[idx(off, j, N)] = item.val;
                dir_matrix[idx(0, j - off, N - off)] = item.idx;
            }
        }
        for (i = off + 1; i < M; i++) {
            for (j = off; j < N; j++) {
                if (window(i, j, p, len_ref, len_query)) {
                    item = dp_dir(ref, query, dim, cost_matrix, i, j, 
                            &dtw_settings, N, dist);
                    cost_matrix[idx(i, j, N)] = item.val;
                    dir_matrix[idx(i - off, j - off, N - off)] = item.idx;
                }
            }
        }
    }
    return cost_matrix[idx(M - 1, N - 1, N)];
}

/**
 * 
 * @brief Dynamic Time Warping algorithm(without traceback)
 * @param ref reference sequence
 * @param query query sequence
 * @param dim number of dimensions of the sequences
 * @param len_ref length of the reference sequence
 * @param len_query length of the query sequence
 * @param dist pointer to a distance function
 * @param dp pointer to a step function
 * @param window pointer to a window function
 * @param p  window parameter
 * @param cost_matrix cost matrix
 * @param dtw_settings structure with dtw settings
 * @return distance between reference and query sequences 
 * 
 */
double cdtwnopath(double* ref, double* query, int dim, int len_ref,
        int len_query, dist_fptr dist, dp_fptr dp, window_fptr window,
        double p, double *cost_matrix, struct t_dtw_settings dtw_settings) {

    int off = dtw_settings.offset;
    /* memory was already allocated */
    /* extending matrix */
    int M = len_ref + off;
    int N = len_query + off;
    register int i = 0;
    register int j = 0;
    double w = 0;
    double s = 0;
    bool fast_glob = (dtw_settings.window_type == _PALIVAL ||
            dtw_settings.window_type == _PALIVAL_MOD);
    
    /* no window or fast window case */
    if (fast_glob || (dtw_settings.window_type ==_NOWINDOW)) {
        if (fast_glob) {
            if (dtw_settings.window_type == _PALIVAL_MOD)
                w = dtw_settings.window_param * (double)len_ref;
            else
                w = dtw_settings.window_param;
            s = len_query / (double)len_ref;
        } 
        else {
            w = INFINITY;
            s = 1;
        }
        // primeiro elemento nao é o do offset... caralho
        cost_matrix[idx(off, off, N)] = dist(&ref[idx(off-off, 0, dim)],
                                             &query[idx(off-off, 0, dim)],
                                             dim);
        for (j = max2(off + 1, _round(s * (off - w))); 
                j < min2(N, _round(s * (off + w) + 1)); j++) {
            cost_matrix[idx(off, j, N)] = dp(ref, query, dim, cost_matrix, off,
                    j, &dtw_settings, N, dist);
        }

        for (i = off + 1; i < M; i++) {
            for (j = max2(off, _round(s * (i - w))); 
                    j < min2(N, _round(s * (i + w) + 1)); j++)
                cost_matrix[idx(i, j, N)] = dp(ref, query, dim, cost_matrix, i,
                        j, &dtw_settings, N, dist);
        }
    }/* slow window case */
    else {
        cost_matrix[idx(off, off, N)] = dist(&ref[idx(off-off, 0, dim)],
                                             &query[idx(off-off, 0, dim)],
                                             dim);
        for (j = off + 1; j < N; j++) {
            if (window(off, j, p, len_ref, len_query))
                cost_matrix[idx(off, j, N)] = dp(ref, query, dim, cost_matrix, 
                        off, j, &dtw_settings, N, dist);
        }
        for (i = off + 1; i < M; i++) {
            for (j = off; j < N; j++) {
                if (window(i, j, p, len_ref, len_query))
                    cost_matrix[idx(i, j, N)] = dp(ref, query, dim, cost_matrix,
                            i, j, &dtw_settings, N, dist);
            }
        }
    }
    return cost_matrix[idx(M - 1, N - 1, N)];
}

/**
 * @brief Dynamic Time Warping -  main entry for the cmdtw 
 * @param ref reference sequence
 * @param query query sequence
 * @param dim number of dimensions for each point
 * @param len_ref length of the reference sequence
 * @param len_query length of the query sequence
 * @param cost_matrix cost matrix, pointer to an allocated memory for the 2 
 * dimensional matrix. (with extra size)
 * @param path warping path, pointer to an allocated array, size is maximum 
 * possible path length
 * @param true_path_len length of the computed warping pat
 * @param dtw_settings structure with dtw settings
 * 
 * @return distance between reference and query sequences 
 */
double cmdtw(double* ref_flat, double* query_flat, int dim, int len_ref,
        int len_query, double* cost_matrix,
        struct t_path_element* path, int *true_path_len,
        struct t_dtw_settings dtw_settings) {

    double distance = 0; /* init distance */
    double p = 0; /* init window function param */

    //register int k = 0;
    //register int o = 0;

   // double** ref = malloc(len_ref * sizeof (double *));
   // double** query = malloc(len_query * sizeof (double *));

    
  //  for (k = 0; k < len_ref; k++) {
  //      ref[k] = malloc(dim * sizeof (double));
  //  }
  //  for (k = 0; k < len_query; k++) {
  //      query[k] = malloc(dim * sizeof (double));
  //  }

    dp_fptr dp = NULL; /* distance function used */
    dpdir_fptr dp_dir = NULL; /* init pointer to step function*/

    /* pointer to window function */
    window_fptr window = choose_window(&dtw_settings);
    
    /* pointer to distance function */
    //dist_fptr dist = choose_dist(dtw_settings.dist_type);
    dist_fptr dist = choose_dist(dtw_settings);

    int *dir_matrix = NULL;
    int *path_points = NULL;
    int path_points_count = 0;
    
    /* I've made all things to use double **, but numpy can use only double*... 
     * so, It's necessary do this POG */
//     #define DEBUG 1
#ifdef DEBUG
//    printf("DIM %d REF_LEN %d QUERY_LEN %d offst  %d\n", dim, len_ref, len_query,
//           dtw_settings.offset);
//    printf("Ref\n");
//    for(o = 0; o < len_ref; o++) {
//        for(k = 0; k < dim; k++)
//            printf("%.2f ",ref_flat[idx(o, k, dim)]);
//        printf("\n");
//    }
#endif

//    for(o = 0; o < len_ref; o++)
//        for(k = 0; k < dim; k++)
//            ref[o][k] = ref_flat[idx(o, k, dim)];

#ifdef DEBUG
 //   printf("Query\n");
//    for(o = 0; o < len_query; o++) {
//        for(k = 0; k < dim; k++)
//            printf("%.2f ",query_flat[idx(o, k, dim)]);
//        printf("\n");
//    }
#endif

//    for(o = 0; o < len_query; o++)
//        for(k = 0; k < dim; k++)
//            query[o][k] = query_flat[idx(o, k, dim)];
    
    /* assign step function */
    if (dtw_settings.compute_path)
        dp_dir = choose_dpdir(dtw_settings.dp_type);
    else
        dp = choose_dp(dtw_settings.dp_type);

    /* assign window parameter */
    p = choose_window_param(&dtw_settings, len_ref, len_query);


    /* Allocate Memory for cost_matrix */
    fill_matrix(cost_matrix, len_ref, len_query, dtw_settings);

    /* dtw without traceback case(only cost_matrix and distance) */
    if (dtw_settings.compute_path == false) {
        distance = cdtwnopath(ref_flat, query_flat, dim, len_ref, len_query, dist, dp,
                window, p, cost_matrix, dtw_settings);
    }
    /* dtw with traceback case */
    else {
        /* allocate direction matrix */
        dir_matrix = (int*) calloc((len_ref * len_query), sizeof(int));
        /* call cdtwpath, computes distance, cost matrix and direction matrix */
        distance = cdtwpath(ref_flat, query_flat, dim, len_ref, len_query, dist, dp_dir,
                window, p, cost_matrix, dir_matrix, dtw_settings);

        /* if distance is INFINITY there is not any path */
        if (distance == INFINITY) {
            *true_path_len = 0;
            return INFINITY;
        }
        /* allocate path points(direction array) */
        path_points = (int*) calloc((len_ref + len_query), sizeof(int));

        /* compute path(directions array) from direction matrix */
        path_points_count = direct_matrix_to_path_points(dir_matrix, 
                path_points, len_ref, len_query, 
                choose_path_pattern(dtw_settings));

        /* cleaning */
        free(dir_matrix);
        /* compute warping path(finally, as array of the path elements */
        *true_path_len = create_path_from_pattern(
                choose_path_pattern(dtw_settings), len_ref, len_query, 
                path_points, path_points_count, path);
        /* cleaning */
        free(path_points);
    }
    
    /* cleaning */
//    for (k = 0; k < len_ref; k++) {
//        free(ref[k]);
//    }
//    for (k = 0; k < len_query; k++) {
//        free(query[k]);
//    }
//    free(ref);
//    free(query);
    
    /* euclidian metric case */
  //  if(dtw_settings.dist_type == _EUCLID)
  //      distance = sqrt(distance);
#ifdef DEBUG
    printf("Distance C: %lf\n",distance);
    _print_matrix(cost_matrix, len_ref+dtw_settings.offset, len_query+dtw_settings.offset);

    for (int i = 0; i< *true_path_len; i++)
        printf("(%d,%d)\n", path[i].i, path[i].j);
#endif
    /* normalization case */
    if(dtw_settings.norm)
        return distance/(double)(len_ref+len_query);  //TOTO: ADD NORM FACTORS
    /* there is no normalization */
    else
        return distance;
}

