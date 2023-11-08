#include <math.h> 
#include "mex.h"

#define RADEMACHER() (2 * (rand() % 2) -1)

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

	if (zeta > d) zeta = d;

	double val = 1/sqrt((double) zeta);

        plhs[0] = mxCreateSparse(d,m,nnz,false);
        double *vals  = mxGetPr(plhs[0]);
        mwIndex *rows = mxGetIr(plhs[0]);
        mwIndex *colstarts = mxGetJc(plhs[0]);

	// Set values
	for (int i = 0; i < nnz; ++i) {
	    vals[i] = RADEMACHER() * val;
	}

	// Set column starts
	for (int i = 1; i < m+1; ++i){
	    colstarts[i] = i*zeta;
	}

	// Set row indices
        for (int i = 0; i < m*zeta; i += zeta) {
	    int idx = 1;
	    rows[i] = rand() % d;
	    while (idx < zeta) {
		rows[i+idx] = rand() % d;
		int j = 0;
		for (; j < idx; ++j) {
		    if (rows[i+idx] == rows[i+j]) break;
		}
		idx += (int) (j == idx);
	    }
	}
    #endif
}
