/*
 * =============================================================
 * graphcutsweep.c - Compute the cut values for a sweep across
 * an ordering of vertices in a graph.
 *
 * David Gleich
 * =============================================================
 */

/*
 * 19 February 2007
 * Updated to use Matlab 2006b sparse matrix interface
 */

#define NUM_SWEEPS 7

#define char16_t UINT16_T

#include "mex.h"
#include "math.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif // MX_API_VER

#ifndef _WIN32
double min(double a, double b)
{
    return ((a)*((a) < (b)) + (b)*((a) >= (b)));
}

double max(double a, double b)
{
    return ((a)*((a) > (b)) + (b)*((a) <= (b)));
}
#endif 

/*
 * Compute the cut sweeps.
 * The first set of paramters is the matrix, the second set is the ordering
 * the final paramter is the output matrix of cut values.
 *
 * cv[NUM_SWEEPS*i+SWEEP_INDEX] = cut value for SWEEP_INDEX with the
 * split at vertex i.
 *
 * @param norder the number of elements in sorder
 */
void cut_sweeps(int n, double nnz, 
                mwIndex *A_row, mwIndex *A_col, double *A_val, 
                int *sorder, mwSize norder, double *cv)
{
    mwSize i, j;
    
    int *vindex;

    /* nnz is 2*sum(Aij) */
    double cut = 0;
    double assoc_a = cut;
    double assoc_b = nnz - cut;
    double assoc_aa = 0;
    double assoc_bb = nnz;
    
    vindex = mxMalloc(n*sizeof(int));
    for (i = 0; i < n; i++) { vindex[i] = 0; }
    
    /* mexPrintf("n=%i; nnz=%i\n", n, nnz); */
    
    /* only through n-1 because at the end we get back to all
     * vertices on one side. */
    for (i = 0; i < norder; i++) {
        /* add perm[i] to set A and recompute */
		int vertex = sorder[i];
        
        double edges_a = 0.0;
        double edges_b = 0.0;
        
        double num_edges = 0.0;
        
        /* A_col[sorder[i]+1] - A_col[sorder[i]]; */
        
        for (j = A_col[vertex]; j < A_col[vertex+1]; ++j)
        {
            num_edges += A_val[j];
            /* mexPrintf("%i -> %i\n", vertex, A_row[j]); */
            if (vindex[A_row[j]])
			{
				/* this neighbor is in A */
				edges_a+=A_val[j];
			}
			else
			{
				/* this neighbor is in B */
				edges_b+=A_val[j];
			}
        }
        
        
        
        /* update the cuts and associations */
		cut -= (double)edges_a;
		cut += (double)edges_b;
		assoc_a += (double)(num_edges);
		assoc_b -= (double)(num_edges);
        assoc_aa += 2.0*(double)edges_a;
        assoc_bb -= 2.0*(double)edges_b;
        
        /*mexPrintf("%f\t%f\t%f\t%f\t%f\n", cut, assoc_a, assoc_b, 
            assoc_aa, assoc_bb); */
        
        vindex[vertex] = 1;
        
        cv[NUM_SWEEPS*i + 0] = cut/assoc_a + cut/assoc_b;
        cv[NUM_SWEEPS*i + 1] = assoc_aa/(double)(i+1) + assoc_bb/(double)(n-i-1);
        cv[NUM_SWEEPS*i + 2] = cut/(double)(i+1) + cut/(double)(n-i-1);
        cv[NUM_SWEEPS*i + 3] = cut/min(assoc_a, assoc_b);
        cv[NUM_SWEEPS*i + 4] = cut/(double)min(i+1, n-i-1);
        cv[NUM_SWEEPS*i + 5] = cut;
        cv[NUM_SWEEPS*i + 6] = assoc_a;
    }
        
}

/**
 * Calculate the various cut values for a graph based on a potential split.
 *
 * The formulae for the cut values are
 *
 *              cut(A,B)     cut(A,B)
 * NCut(A,B) = ---------- + ----------
 *             assoc(A,V)   assoc(B,V)
 *
 *
 *                 assoc(A,A)   assoc(B,B)
 * AvgAssoc(A,B) = ---------- + ----------
 *                     |A|         |B|
 *
 *               cut(A,B)   cut(A,B)
 * AvgCut(A,B) = -------- + --------
 *                 |A|        |B|
 *
 *                             cut(A,B)  
 * QuotientCut(A,B) = --------------------------- 
 *                    min(assoc(A,V), assoc(B,V))
 *                  
 *
 * where
 *
 * cut(A,B) = sum of edges between A and B
 * assoc(A,B) = sum of edges between A and B
 * V is used to denote the set of vertices in the graph
 *
 * We compute this function incrementally over all possible cut values.
 *
 * Observe that if we have a vertex i that is moving from set B to set
 * A, then we can partition the edges of that vertex into a set 
 * of edges pointing to A and pointing to B.
 *
 * Since cut(A,B) counts the number of edges between the sets,
 * cut(A+v_i,B-v_i) = cut(A,B) - edges_to_A + edges_to_B
 * since the edges between vertex i and A no longer count, but
 * we get to add the edges that are between vertex i and B.
 *
 * Likewise, we increase the number of edges connecting set A
 * and the rest of the graph so 
 * assoc(A+v_i,V) = assoc(A,V) + edges_to_B
 * assoc(A+v_i,V) = assoc(B,V) - edges_to_A
 *
 * Similarily,
 * assoc(A+v_i,A) = assoc(A,A) + edges_to_A
 * assoc(B+v_i,B) = assoc(B,B) - edges_to_B
 *
 * This is the gateway routine into the function.
 * We take in a sparse matrix and an ordering of the vertices.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize i, j, k;
    mwSize n, mrows, ncols;
    
    /* sparse matrix */
    mwSize A_nz;
    mwIndex *A_row, *A_col;
    double *A_val;
    
    mwSize nperm;
    
    double sum_A;
    
    /* ordering */
    int *sorder;
    
    /* cut values */
    double *cv;
    
    /* argument check */
    
    if (nrhs != 2) 
    {
        mexErrMsgTxt("Two inputs required.");
    } 
    else if (nlhs > 2) 
    {
        mexErrMsgTxt("Too many output arguments");
    }
    
    /* The input must be a noncomplex sparse matrix.*/
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if (mrows != ncols ||
        !mxIsSparse(prhs[0]) ||
        !mxIsDouble(prhs[0]) || 
        mxIsComplex(prhs[0])) 
    {
        mexErrMsgTxt("Input must be a noncomplex square sparse matrix.");
    }
    
    /* okay, the first input passes */
    n = mrows;
    
    /* The second input must be a permutation. */
    if (mxGetNumberOfElements(prhs[1]) > n ||
        mxIsSparse(prhs[1]))
    {
        mexErrMsgTxt("Invalid permutation.");
    }
    
    /* Get the sparse matrix */
    A_nz = mxGetNzmax(prhs[0]);
    A_val = mxGetPr(prhs[0]);
    A_row = mxGetIr(prhs[0]);
    A_col = mxGetJc(prhs[0]);
    
    nperm = mxGetNumberOfElements(prhs[1]);
    
    /* Get the permutation */
    {
        double *p = mxGetPr(prhs[1]);
        sorder = mxMalloc(nperm*sizeof(int));
        for (k = 0; k < nperm; k++)
        {
            sorder[k] = p[k] - 1; /* we are 0 based here */
        }
    }
      
    /* Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(NUM_SWEEPS,nperm,mxREAL);
    
    /* Create a C pointer to a copy of the output matrix. */
    cv = mxGetPr(plhs[0]);
    
    sum_A = 0.0;
    for (k = 0; k < A_col[n]; k++)
    {
        sum_A += A_val[k];
    }
    
    cut_sweeps(n, sum_A, A_row, A_col, A_val, sorder, nperm, cv);
    
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[1]) = sum_A;
    
    mxFree(sorder);
    
  
}

