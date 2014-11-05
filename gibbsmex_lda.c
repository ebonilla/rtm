#include "mex.h"
#include <stdlib.h>

/*
 * plhs[0] = *z
 * plhs[1] = **Nwt
 * plhs[2] = **Ndt
 * 
 * prhs[0] = *z
 * prhs[1] = **Nwt
 * prhs[2] = **Ndt
 * prhs[3] = *Nt
 * prhs[4] = *w
 * prhs[5] = *d
 *
 */

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  int N;
  int W;
  int D;
  int T;

  double alpha = 999999999.99; /* alpha is set below */

  double *w;
  double *d;
  double *z;
  double *Nwt;
  double *Ndt;
  double *Nt;

  int i, t, wi, di;
  double Z, U, cumprob;
  double *probs;

  double *z_;
  double *Nwt_;
  double *Ndt_;
  
  N = mxGetM(prhs[0]) * mxGetN(prhs[0]);
  W = mxGetM(prhs[1]);
  T = mxGetN(prhs[1]);
  D = mxGetM(prhs[2]);

  alpha = 0.05 * N / (D*T);   /*mexPrintf("alpha = %f\n", alpha);*/

  z    = mxGetPr(prhs[0]);
  Nwt  = mxGetPr(prhs[1]);
  Ndt  = mxGetPr(prhs[2]);
  Nt   = mxGetPr(prhs[3]);
  w    = mxGetPr(prhs[4]);
  d    = mxGetPr(prhs[5]);
  
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

  
  /******************************************************/
  for (i = 0; i < N; i++) {

    wi = (int)(w[i]);
    di = (int)(d[i]);
    
    t = z_[i];
    Nt[t]--;     
    Nwt_[t*W + wi]--;
    Ndt_[t*D + di]--;

    Z = 0;
    for (t = 0; t < T; t++) {
      
      probs[t]=(Nwt_[t*W+wi]+0.01)*(Ndt_[t*D+di]+alpha)/(Nt[t]+0.01*W);
      
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
    Nt[t]++;     
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
