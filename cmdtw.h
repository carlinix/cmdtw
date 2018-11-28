/**
 * @file cmdtw.h
 * @Author Ricardo Carlini Sperandio
 * @date 25/03/2016
 * @brief Brief description of file
 *
 * Detailed description of file
 */

#ifndef cmdtw_h
#define cmdtw_h

#include <stdio.h>
#include <math.h>

/* distance types */
#define _EUCLID 11          /*euclidean distance*/
#define _EUCLID_SQUARED 12  /*squared euclidean distance*/
#define _MANHATTAN 13       /*manhattan distance*/
#define _CITYBLOCK _MANHATTAN
#define _INDEXMATRIX 14

/*step pattern types, read struct dtw_settings for more info*/
#define _DPW            20 /* With weight */
#define _DP1            21 /* a = b = 1 ; c=2 */
#define _DP2            22 /* a = b = c */
#define _DP3            23 /* without diagonal */
#define _SCP0SYM        24
#define _SCP0ASYM       25
#define _SCP1DIV2SYM    26
#define _SCP1DIV2ASYM   27
#define _SCP1SYM        28
#define _SCP1ASYM       29
#define _SCP2SYM        210
#define _SCP2ASYM       211

#define _NOWINDOW       00
#define _SCBAND         31
#define _PALIVAL        32
#define _ITAKURA        33
#define _PALIVAL_MOD    35


/* macros to min/max */
#define min2(a,b) ((a) < (b) ? a : b)
#define max2(a,b) ((a) > (b) ? a : b)

/* matrix mapping */
#define idx(i,j,J) ( (i)*(J) + (j) )


/* boolean helpper */
typedef int bool;
#define true 1
#define false 0



/*******************************************************************************
 * Structures
 ******************************************************************************/

/*weights for step pattern*/
struct t_weights {
    double a;
    double b;
    double c;
};

struct t_path_element {
    int i;
    int j;
};

struct t_dtw_settings {
    int compute_path;
    int dist_type;
    int dp_type;
    double *dist_mapped;
    int *dist_mapped_len;
    int window_type;
    double window_param;
    int norm;
    int offset;
    struct t_weights weights;
};

struct t_item {
    double val;
    int idx;
};

/*******************************************************************************
 * Helper functions
 ******************************************************************************/

/**
 * @brief  Finds the minimum of 2 doubles ant its position
 * @param a a double
 * @param b b double
 * @return min(a,b) and its position
 *
 */
inline struct t_item min2idx(double a, double b);

/**
 * @brief  Finds the minimum of 3 doubles
 * @param x  x
 * @param y  y
 * @param z z
 * @return min(x,y,z)
 *
 */
inline double min3(double x, double y, double z);

/**
 * @brief  Finds the minimum of 3 doubles and its position (0, 1 or 2)
 * @param x
 * @param y
 * @param z
 * @return struct t_item, t_item = { min(x,y,z), position }
 */
inline struct t_item min3idx(double x, double y, double z);

/**
 * @brief Finds the minimum of n doubles
 * @param arr array
 * @param n length of the array
 * @return  min(arr)
 */
double min_n(double* arr, int n);

/**
 * @brief Finds the minimum of n doubles and its position
 * @param arr array
 * @param n length of the array
 * @return  min(arr) and its position
 */
struct t_item min_nidx(double* arr, int n);

///**
// * @bref Rounds a number 
// * @param number
// * @return number + 1 or number - 1
// */
//inline double _round(double number);


/**
 * @bref Rounds a number 
 * @param number
 * @return number + 1 or number - 1
 */
 inline double _round(double number) {
    return floor(number + 0.5);
}

/*******************************************************************************
 * Distance functions
 ******************************************************************************/
/*step functions args*/
#define _DP_ARGS    double* ref,\
double* query,\
int dim,\
double* cost_matrix,\
int i,\
int j,\
struct t_dtw_settings* t_s,\
int size2, \
double (*dist)(double *a, double *b, int dim)


/* pointer to distance function */
typedef double (*dist_fptr)(double *a, double *b, int dim);

/* pointer to step function */
typedef double (*dp_fptr)(_DP_ARGS);

/* pointer to step function with traceback */
typedef struct t_item(*dpdir_fptr)(_DP_ARGS);


/**
 * @brief  Multidimensional CityBlock distance
 * @param a multidimensional point
 * @param b  multidimensional point
 * @param dim number of dimensions
 * 
 * @return  Manhattan distance between two multidimensional points
 * 
 */
double manhattan(double* a, double* b, int dim);

/**
 * @brief  Multidimensional Euclidian distance
 * @param a multidimensional point
 * @param b  multidimensional point
 * @param dim number of dimensions
 * 
 * @return  Euclidian distance between two multidimensional points
 * 
 */
double euclid(double* a, double* b, int dim);

//FIXME: how can I do this ?
/**
 * @brief
 * @param a
 * @param b
 * @param dim
 * @return 
 */
double indexdist(double *a, double *b, int dim);

/**
 * @brief Chooses right distance functio
 * @param dist_type dtw_settings.dist_type
 * 
 * @return pointer to a distance function
 * 
 */
//dist_fptr choose_dist(int dist_type);
dist_fptr choose_dist(struct t_dtw_settings);
/*******************************************************************************
 * Window functions
 ******************************************************************************/
/* pointer to window function */
typedef bool (*window_fptr)(int i, int j, double r, double I, double J);

/**
 * @bref Computes extra size for the cost matrix (dtw_settings.offset)
 * @param dp_type dtw_settings.dp_type, step pattern type
 * @return int (1,2 or 3)
 * 
 */
int extra_size(int dp_type);

/**
 * @brief This function is for testing only
 * @param i not used
 * @param j not used
 * @param k not used
 * @param I not used
 * @param J not used
 * 
 * @return always return true
 * 
 */
bool nowindow(int i, int j, double k, double I, double J);

/**
 * @brief Sakoe-Chiba band global constraint
 * @param i  row index
 * @param j column index
 * @param r width of the window
 * @param I not used
 * @param J not used
 * @return boolean value, if cost[i][j] is legal or not
 */
bool scband(int i, int j, double r, double I, double J);

/**
 * @brief Itakura global constraints
 * @param i row index
 * @param j column index
 * @param k k = 1/s = I/J = len_ref/len_query
 * @param I length of reference sequence, len_ref
 * @param J length of query sequence, len_query
 * @return boolean value, if cost[i][j] is legal or not
 */
bool itakura(int i, int j, double k, double I, double J);

/**
 *   FIXME: NOTE: function is partly redundant at the moment, it always returns  itakura
 *  Chooses right window function
 *  int dp_type: dtw_settings.win , step function type
 *  returns: dpdir_fptr, pointer to a step function with traceback
 * 
 * @param dtw_settings
 * @return 
 */
window_fptr choose_window(struct t_dtw_settings *dtw_settings);

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
        int len_query);

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
struct t_item dpwdir(_DP_ARGS);
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
double dpw(_DP_ARGS);

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
struct t_item dp1dir(_DP_ARGS);

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
double dp1(_DP_ARGS);
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
struct t_item dp2dir(_DP_ARGS);

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
double dp2(_DP_ARGS);

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
struct t_item dp3dir(_DP_ARGS);

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
double dp3(_DP_ARGS);

/**
 * @brief Chooses right step function(without traceback)
 * @param dp_type dtw_settings.dp_type, step function type
 * 
 * @return pointer to a step function without traceback
 */
dp_fptr choose_dp(int dp_type);

/**
 * @brief Chooses right step function(with traceback)
 * @param dp_type dtw_settings.dp_type, step function type
 * 
 * @return pointer to a step function without traceback
 *
 */
dpdir_fptr choose_dpdir(int dp_type);

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
        struct t_path_element* path);

/**
 * @brief Chooses right path_patterns(2d array) based on current step function
 * @param dtw_settings structure with dtw settings
 *
 * @return path_pattern
 * 
 */
const int (*choose_path_pattern(struct t_dtw_settings dtw_settings))[11];
    
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
        int len_query, const int pattern[6][11]);

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
        struct t_dtw_settings dtw_settings);


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
        struct t_dtw_settings dtw_settings) ;

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
        double p, double *cost_matrix, struct t_dtw_settings dtw_settings) ;

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
double cmdtw(double* ref, double* query, int dim, int len_ref,
        int len_query, double* cost_matrix,
        struct t_path_element* path, int *true_path_len,
        struct t_dtw_settings dtw_settings);


/*******************************************************************************
 * Distance functions
 ******************************************************************************/



/*******************************************************************************
 * Window functions
 ******************************************************************************/


#endif /* cmdtw_h */
