#include "mex.h"
#include <stdlib.h>

/* Semi-collapsed Gibbs sampler - without integrating out Phi
 * plhs[0] = *z
 * plhs[1] = **Nwt
 * plhs[2] = **Ndt
 * 
 * prhs[0] = *z
 * prhs[1] = **Nwt
 * prhs[2] = **Ndt
 * prhs[3] = *w
 * prhs[4] = *d
 * prhs[5] = *alpha
 * prhs[6] = **PHIwt
 *
 */

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  int N;
  int W;
  int D;
  int T;

  double alpha;

  double *w;
  double *d;
  double *z;
  double *Nwt;
  double *Ndt;
  double *PHIwt;
 
  int i, t, wi, di;
  double Z, U, cumprob;
  double *probs;
 
  double *z_;
  double *Nwt_;
  double *Ndt_;
  
  N = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  W = mxGetM(prhs[6]);
  T = mxGetN(prhs[6]);
  D = mxGetM(prhs[2]);
 
  z     = mxGetPr(prhs[0]);
  Nwt   = mxGetPr(prhs[1]);
  Ndt   = mxGetPr(prhs[2]);
  w     = mxGetPr(prhs[3]);
  d     = mxGetPr(prhs[4]);
  alpha = *mxGetPr(prhs[5]);  	
  PHIwt = mxGetPr(prhs[6]);
   
  probs = mxMalloc(T * sizeof(double));
  
  plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(W,T,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(D,T,mxREAL);
  
  z_   = mxGetPr(plhs[0]);
  Nwt_ = mxGetPr(plhs[1]);
  Ndt_ = mxGetPr(plhs[2]);
  
  /* copy input arrays to output arrays */
  for (i = 0; i < N; i++) {
    w[i]--;
    d[i]--;
    z_[i] = z[i]-1;
  }

  for (i = 0; i < W*T; i++) Nwt_[i] = Nwt[i];
  for (i = 0; i < D*T; i++) Ndt_[i] = Ndt[i];

  /*printf("aplha=%.2f\n", alpha);*/
   
  /******************************************************/
  for (i = 0; i < N; i++) {

    wi = (int)(w[i]);
    di = (int)(d[i]);
    
    t = z_[i];
    Nwt_[t*W + wi]--; 
    Ndt_[t*D + di]--;

    Z = 0;
    for (t = 0; t < T; t++) {
      probs[t] = ( PHIwt[t*W + wi] )  * (Ndt_[t*D + di]+alpha);
      Z += probs[t];
    }
    
    U = Z * drand48();
    cumprob = probs[0];
    t = 0;
    while (U > cumprob) {
      t++;
      cumprob += probs[t];
    }
   
    z_[i] = t;
    Nwt_[t*W + wi]++;
    Ndt_[t*D + di]++;

  }

  for (i = 0; i < N; i++) {
    z_[i]++;
    w[i]++;
    d[i]++;
  }
 
  mxFree(probs);

} 
