/**
 * Author: rcarlini@gmail.com
 *  
 */
#include <stdio.h>
#include <stdlib.h>
#include "cmdtw.h"


void print_matrix(double* matrix, int size1, int size2) {
    register int i = 0;
    register int j = 0;
    printf("\n");
    for(i = 0; i<size1; i++) {
        for(j = 0; j<size2; j++)
            printf("%.2f\t  ", (double)matrix[idx(i, j, size2)]);
        printf("\n");
    }
}

void print_intarr(int* arr, int n) {
    register int i = 0;
    for(i = 0; i<n; i++)
        printf("%4i",arr[i]);
    printf("\n");
}

int main(int argc, const char *argv[]) {
    //insert code here...
    printf("Hello, World!\n");


    struct t_dtw_settings dtw_settings;

    dtw_settings.compute_path = false;
    dtw_settings.dist_type = _EUCLID;
    dtw_settings.dp_type = _DP2;
    dtw_settings.window_type = _SCBAND;
    dtw_settings.window_param = 2.0;
    dtw_settings.norm = false;
    dtw_settings.offset = extra_size(dtw_settings.dp_type);
    dtw_settings.weights.a = 1.0;
    dtw_settings.weights.b = 1.0;
    dtw_settings.weights.c = 1.0;

    int dim = 1;
    int len_ref = 5;
    //sizeof(r) / sizeof(r[0]);
    int len_query = 5;
    //sizeof(q) / sizeof(q[0]);
    ///*
    
#if 0
    double **r = NULL;
    double **q = NULL;

    r = malloc(len_ref * sizeof (double *));
    q = malloc(len_query * sizeof (double *));
    register int k, j;
    for (j = 0; j < len_ref; j++) {
        for (k = 0; k < dim; k++) {
            r[j] = malloc(dim * sizeof (double));
            q[j] = malloc(dim * sizeof (double));
        }
    }
    //r[0][0] = 0;
    r[0][0] = 5;
    r[1][0] = 1;
    r[2][0] = 1;
    r[3][0] = 1;
    r[4][0] = 1;
    printf("%lf\n",*(*(r+0)+2));
    exit(0);

    //q[0][0] = 0;
    q[0][0] = 10;
    q[1][0] = 1;
    q[2][0] = 1;
    q[3][0] = 1;
    q[4][0] = 2;
    //*/

#endif
    double * r = (double *) calloc((len_ref) * dim, sizeof(double));
                                      
       double * q = (double *) calloc((len_query) * dim, sizeof(double));
         r[idx(0, 0, dim)] = 5;
         r[idx(1, 0, dim)] = 1;
         r[idx(2, 0, dim)] = 1;
         r[idx(3, 0, dim)] = 1;
         r[idx(4, 0, dim)] = 1;
         
                
       
         q[idx(0, 0, dim)] = 10;
         q[idx(1, 0, dim)] = 1;
         q[idx(2, 0, dim)] = 1;
         q[idx(3, 0, dim)] = 1;
         q[idx(4, 0, dim)] = 2;
         
    
    double x;

   // double *cost = (double *) malloc(sizeof(double) * ((len_ref +
   //         dtw_settings.offset) * (len_query + dtw_settings.offset)));
    double * cost = (double *) calloc(((len_ref + dtw_settings.offset) *
                                       (len_query + dtw_settings.offset)),
                                      sizeof(double));


    struct t_path_element *path = (struct t_path_element *) calloc((len_ref + len_query), sizeof(struct t_path_element));
    int path_len = 0;
    x = cmdtw(r, q, dim, len_ref, len_query, cost, path, &path_len, dtw_settings);
    printf("CUSTO: %lf",x);
    print_matrix(cost, len_ref+dtw_settings.offset, len_query+dtw_settings.offset);
    for (int i = 0; i<path_len; i++)
        printf("(%d,%d)\n", path[i].i, path[i].j);
    return 0;
}
