#include "2dbat.sort.h"

// 2DBAT user defined functions
// SORT related

// Sort arrays
int array1D_comp(const void *a, const void *b);
int array2D_comp(const void *a, const void *b);
int array3D_comp(const void *a, const void *b);
int comparisonFunctionInt(const void *a, const void *b);
int comparisonFunctionFloat(const void *a, const void *b);
int comparisonFunctionDouble(const void *a, const void *b);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* QSORT COMPARE FUNCTION - SORTS BY COLUMN 1
    aescending AND THEN BY COLUMN 1 ASCENDING */
int array1D_comp(const void *a, const void *b)
{
    double *arr1 = *((double **) a);
    double *arr2 = *((double **) b);

    double diff1 = arr1[0] - arr2[0];
    if (diff1) return diff1;
    return arr2[0] - arr1[0];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* QSORT COMPARE FUNCTION - SORTS BY COLUMN 2
    aescending AND THEN BY COLUMN 2 ASCENDING */
int array2D_comp(const void *a, const void *b)
{
    double *arr1 = *((double **) a);
    double *arr2 = *((double **) b);

    double diff1 = arr1[1] - arr2[1];
    if (diff1) return diff1;
    return arr2[1] - arr1[1];
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* QSORT COMPARE FUNCTION - SORTS BY COLUMN 3
    aescending AND THEN BY COLUMN 3 ASCENDING */
int array3D_comp(const void *a, const void *b)
{
    double *arr1 = *((double **) a);
    double *arr2 = *((double **) b);

    double diff1 = arr1[2] - arr2[2];

    if (diff1) return diff1;

    return arr2[2] - arr1[2];
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* comparison function inside qsort() for getMedian */
int comparisonFunctionDouble(const void *a, const void *b)
{
    if (*(double*)a <  *(double*)b) return -1;
    if (*(double*)a == *(double*)b) return  0;

    return 1;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* comparison function inside qsort() for */
int comparisonFunctionFloat(const void *a, const void *b)
{
    if (*(float*)a <  *(float*)b) return -1;
    if (*(float*)a == *(float*)b) return  0;

    return 1;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* comparison function inside qsort() for integer array */
int comparisonFunctionInt(const void *a, const void *b)
{
    if (*(int*)a <  *(int*)b) return -1;
    if (*(int*)a == *(int*)b) return  0;

    return 1;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int sshist(double *x, int n_data)
{
    int i, j, idx;
    int N_MIN = 4, N_MAX = 2000;

    //double hist_x[9999], hist_y[9999];
    double *hist_x, *hist_y;
    double x_max, x_min, cmin;
    double N[N_MAX-N_MIN];
    double C[N_MAX-N_MIN], C_t[N_MAX-N_MIN], D[N_MAX-N_MIN];
    double l, u;
    double sum, v, k;
    gsl_histogram *h_input;


    // dynamic x_t temp array
    double *x_t = (double*)malloc(sizeof(double) * n_data); // n_lines x 1

    if(N_MAX > n_data) N_MAX = n_data;

    for(i=0; i<n_data; i++)
    {
        x_t[i] = x[i];
    }
    qsort(x_t, n_data, sizeof(x_t[0]), comparisonFunctionDouble);

    x_min = x_t[0]; // minimum
    x_max = x_t[n_data-1]; // maximum

    for(i=0; i<(N_MAX-N_MIN); i++)
    {
        N[i] = (double)(i+N_MIN);
        D[i] = (x_max-x_min) / N[i];
    }

    for(i=0; i<N_MAX-N_MIN; i++)
    {
        h_input = gsl_histogram_alloc(N[i]);
        gsl_histogram_set_ranges_uniform(h_input, x_min, x_max);
    
        // assign histograms
        for (j=0; j<n_data; j++)
        {
            gsl_histogram_increment(h_input, x[j]);
        }

        hist_x = (double*)malloc(sizeof(double) * h_input->n); // h_input->n x 1
        hist_y = (double*)malloc(sizeof(double) * h_input->n); // h_input->n x 1

        for(j=0; j<h_input->n; j++)
        {
            gsl_histogram_get_range(h_input, j, &l, &u);
            hist_x[j] = (l+u)/2.;
            hist_y[j] = gsl_histogram_get(h_input, j);
        }

        k = gsl_stats_mean(hist_y, 1, h_input->n);
        sum = 0;
        for(j=0; j<h_input->n; j++)
        {
            sum += (hist_y[j]-k)*(hist_y[j]-k);
        }

        v = sum / N[i]; // variance of event count
        C[i] = (2*k-v)/(D[i]*D[i]);
        C_t[i] = (2*k-v)/(D[i]*D[i]);

        free(hist_x);
        free(hist_y);
    }

    // Optimal bin size selecction
    qsort(C_t, N_MAX-N_MIN, sizeof(C_t[0]), comparisonFunctionDouble);
    cmin = C_t[0];

    for(i=0; i<N_MAX-N_MIN; i++)
    {
        if( (C[i] - cmin) < 0.0001)
        {
            idx = i;
            break;
        }
    }

    free(x_t);
    gsl_histogram_free(h_input);

    return idx; // return optimal bin size
}


// --- End of line

