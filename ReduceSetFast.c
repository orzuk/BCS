/* Cluster database sequences - unite sequences with few mismatches
 * =============================================================
 * ReduceSet=ReduceSetFast (numOfCreatures,seqLen,seqInt,diff)
 * numofCreatures - the number of sequences in the database
 * seqLen - the length of each sequence
 * seqInt - the original sequence database matrix
 * diff - the maximal difference for removal of a sequence
 * return a vector [reduceSet] of indices of the reduced sequence database
 * =============================================================
 */

#include "stdio.h"
#include "mex.h"
#include "matrix.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int numOfCreatures,seqLen;
    unsigned char *seqInt;
    int diff;

    int cCreature,dCreature,cNuc,ccNuc;
    double *outDat;
    int cdiff;
    
    numOfCreatures = (int) mxGetScalar(prhs[0]);
    seqLen = (int) mxGetScalar(prhs[1]);
    seqInt = (unsigned char *)  mxGetPr(prhs[2]);
    diff = mxGetScalar(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix(1,numOfCreatures, mxREAL);
    outDat=mxGetPr(plhs[0]);
    
    for (cCreature=0;cCreature<numOfCreatures;cCreature++) {
        if (outDat[cCreature]==0)
            for (dCreature=cCreature+1;dCreature<numOfCreatures;dCreature++) {
                if (outDat[dCreature]==0) {
                    cdiff=0;
                    ccNuc=0;
                    for (cNuc=0;cNuc<seqLen;cNuc++) {
                        if (seqInt[cCreature+ccNuc]!=seqInt[dCreature+ccNuc]) {
                            cdiff++;
                            if (cdiff>diff) {
                                cNuc=seqLen;
                            }
                        }
                        ccNuc+= numOfCreatures;
                    }
                    if (cdiff<=diff) {
                        outDat[dCreature]=1;
                    }
                }
            }
    }
}
