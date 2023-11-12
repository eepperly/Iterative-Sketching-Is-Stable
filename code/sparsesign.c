#include <math.h> 
#include "mex.h"

#define RADEMACHER() (2 * (rand() % 2) -1)

int log_base_d(long a, int d) {
    int output = 0;
    while (a > 0) {
        a /= d;
        output += 1;
    }
    return output-1;
}

void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* If the mex file was built using interleaved complex flag, display
     * an error and exit.
     */
    #if MX_HAS_INTERLEAVED_COMPLEX
        mexErrMsgIdAndTxt( "MATLAB:sparsesign:invalidBuild",
            "This example supports separate complex representation only."
            "\nPlease do not use interleaved complex flag with mex command."
            "\n\nPlease use following command to compile using separate complex:"
            "\nmex -R2017b fulltosparse.c");
    #else
        /* Check for proper number of input and output arguments */
        if (nrhs != 3) {
            mexErrMsgIdAndTxt( "MATLAB:sparsesign:invalidNumInputs",
                    "Three input arguments required.");
        }
        if(nlhs > 1){
            mexErrMsgIdAndTxt( "MATLAB:sparsesign:maxlhs",
                    "Too many output arguments.");
        }

        /* Check data type of input argument  */
	for (int i = 0; i < 3; ++i) {
	    if( !mxIsNumeric(prhs[0]) || 
		mxIsComplex(prhs[0]) ||
		mxGetNumberOfElements(prhs[0])!=1 ) {
		mexErrMsgIdAndTxt("MATLAB:sparsesign:notScalar","Inputs must be scalars.");
	    }
	}

	mwSize d = (mwSize) mxGetScalar(prhs[0]);
	mwSize m = (mwSize) mxGetScalar(prhs[1]);
	mwSize zeta = (mwSize) mxGetScalar(prhs[2]);
	mwSize nnz = m*zeta;

    int idx_per_rand = log_base_d((long) RAND_MAX + 1, d);
    int bit_per_rand = log_base_d((long) RAND_MAX + 1, 2);

	if (zeta > d) zeta = d;

	double lowval = -1/sqrt((double) zeta);
    double increment = -2*lowval;

    plhs[0] = mxCreateSparse(d,m,nnz,false);
    double *vals  = mxGetPr(plhs[0]);
    mwIndex *rows = mxGetIr(plhs[0]);
    mwIndex *colstarts = mxGetJc(plhs[0]);

	// Set values
    int myrand = rand();
	for (int i = 0; i+bit_per_rand < nnz; i += bit_per_rand) {
        for (int j = i; j < i+bit_per_rand; ++j) {
	        vals[j] = (myrand % 2) * increment + lowval;
            myrand = myrand >> 1;
        }
        myrand = rand();
	}
    for (int i = bit_per_rand*(nnz/bit_per_rand); i < nnz; ++i) {
        vals[i] = (myrand % 2) * increment + lowval;
        myrand = myrand >> 1;
    }

	// Set column starts
	for (int i = 1; i < m+1; ++i){
	    colstarts[i] = i*zeta;
	}

	// Set row indices
    myrand = rand();
    int ir = 0;
    for (int i = 0; i < m*zeta; i += zeta) {
	    int idx = 0;
	    while (idx < zeta) {
		    rows[i+idx] = myrand % d;
            ir++;
            if (ir == idx_per_rand) {
                ir = 0;
                myrand = rand();
            } else {
                myrand /= d;
            }
		    int j = 0;
		    for (; j < idx; ++j) {
		        if (rows[i+idx] == rows[i+j]) break;
		    }
		    idx += (int) (j == idx);
	    }
	}
    #endif
}