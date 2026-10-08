/* Additional native harness; original controller/operator/ODE are unchanged.
   Noise is generated and seeded in Python. No real-time execution is claimed. */
#include "sil_core.c"
API void sil_advance(double t,double *state,const int *levels,double h,double Tc,int aged,double slope);

/* Synchronous acquisition counterpart of sil_validation.run, with no input
   or gate-delay intervention. Preserve the four-substep integration grid. */
API int sil_mc_run(int n,int bits,int aged,double Tc,const double *noise,
                  double *current,double *sensed,int *states,double *junction){
  double state[99]={0},sense[3],capture[3]={0};int old[3]={0},q[3],freeq[3],k,j,m,h;
  double lsb=bits?64.0/pow(2.0,bits):0,limit=bits?pow(2.0,bits-1):0;
  for(j=0;j<24;j++)junction[j]=Tc;
  for(k=0;k<=n;k++){
    for(j=0;j<3;j++){
      sense[j]=capture[j]+noise[3*k+j];
      if(bits){double a=sense[j]/lsb,lo=floor(a),frac=a-lo;
        /* Explicit nearest-even rounding, matching numpy.rint. */
        double rounded=frac<.5?lo:frac>.5?lo+1:(fmod(fabs(lo),2)==0?lo:lo+1);
        if(rounded < -limit)rounded=-limit;if(rounded > limit-1)rounded=limit-1;
        sense[j]=rounded*lsb;}
      sensed[3*k+j]=sense[j];
    }
    if(k==n)break;
    sil_control(k*50e-6,sense,old,q);
    for(j=0;j<3;j++){states[3*k+j]=q[j];freeq[j]=state[j]>=0?(old[j]<q[j]?old[j]:q[j]):(old[j]>q[j]?old[j]:q[j]);}
    for(h=0;h<4;h++)sil_advance(k*50e-6+h*.25e-6,state,freeq,.25e-6,Tc,aged,1);
    for(h=0;h<4;h++)sil_advance(k*50e-6+1e-6+h*12.25e-6,state,q,12.25e-6,Tc,aged,1);
    for(j=0;j<3;j++){current[3*(k+1)+j]=capture[j]=state[j];if(fabs(state[j])>=31.99)return -5;old[j]=q[j];}
    for(j=0;j<24;j++){double T=Tc;for(m=0;m<4;m++)T+=state[3+4*j+m];junction[24*(k+1)+j]=T;if(T<25-1e-9||T>175+1e-9)return -6;}
  }
  return 0;
}

/* Assemble the same nominal midpoint operator using oracle junction values.
   Oracle compensation is declared and does not emulate a calibrated sensor. */
API void sil_mc_observe(int n,const double *sensed,const int *states,
                       const double *junction,double *X,double *y,int *keep){
  int k,j,zero[3]={0};double v[24],r[24];
  for(k=0;k<n;k++){
    for(j=0;j<24;j++){double T=(junction[24*k+j]+junction[24*(k+1)+j])/2;
      r[j]=j%2?.015:.020;
      v[j]=(j%2?(T<=125?1.45-.0005*(T-25):1.40):1.65+.002*(T-25))-(j%2?20:40)*r[j];}
    keep[k]=sil_observe(k*50e-6,sensed+3*k,sensed+3*(k+1),states+3*k,k?states+3*(k-1):zero,v,r,0,X+100*k,y+2*k);
  }
}

/* Additional QR fit with the full 2x2 target covariance. sigma_weight is a
   declared reference weight, not a fitted feedback/quantization likelihood. */
API int sil_fit_reporting(int n,const double *X,const double *y,const int *keep,
                         double sigma_weight,double *target,double *cov,double *diag){
  int k,h,j,b,d,m=0,row,prev=-2;double *a,*rhs,last=0;
  double diagonal=(40.125*40.125+39.875*39.875)*sigma_weight*sigma_weight,off=40.125*(-39.875)*sigma_weight*sigma_weight;
  double prior[2][NB]={{0}},prior_y[2]={0},scales[NB]={0},beta[NB],largest=0,smallest=1e300,z[2][NB]={{0}};
  for(k=0;k<n;k++)if(keep[k])m+=2;
  diag[0]=m/2;if(m<NB||sigma_weight<=0)return -1;
  a=(double*)calloc((size_t)m*NB,sizeof(double));rhs=(double*)calloc(m,sizeof(double));if(!a||!rhs){free(a);free(rhs);return -2;}
  row=0;
  for(k=0;k<n;k++)if(keep[k]){
    double sub=k==prev+1?off/last:0;last=sqrt(diagonal-sub*sub);
    for(h=0;h<2;h++){
      for(b=0;b<NB;b++){double value=0;for(j=0;j<NP;j++)value+=X[(k*2+h)*NP+j]*basis[j][b];value=(value-sub*prior[h][b])/last;prior[h][b]=value;a[row*NB+b]=value;scales[b]+=value*value;}
      rhs[row]=(y[k*2+h]-sub*prior_y[h])/last;prior_y[h]=rhs[row];row++;
    }prev=k;
  }
  for(b=0;b<NB;b++){scales[b]=sqrt(scales[b]);if(scales[b]<1e-14){free(a);free(rhs);return -3;}for(row=0;row<m;row++)a[row*NB+b]/=scales[b];}
  for(b=0;b<NB;b++){
    double norm=0,alpha,v0,vnorm=1;for(row=b;row<m;row++)norm=hypot(norm,a[row*NB+b]);
    if(norm<1e-10){free(a);free(rhs);return -4;}alpha=a[b*NB+b]>=0?-norm:norm;v0=a[b*NB+b]-alpha;
    for(row=b+1;row<m;row++)a[row*NB+b]/=v0;
    for(row=b+1;row<m;row++)vnorm+=a[row*NB+b]*a[row*NB+b];
    for(d=b+1;d<NB;d++){double dot=a[b*NB+d];for(row=b+1;row<m;row++)dot+=a[row*NB+b]*a[row*NB+d];dot*=2/vnorm;a[b*NB+d]-=dot;for(row=b+1;row<m;row++)a[row*NB+d]-=dot*a[row*NB+b];}
    {double dot=rhs[b];for(row=b+1;row<m;row++)dot+=a[row*NB+b]*rhs[row];dot*=2/vnorm;rhs[b]-=dot;for(row=b+1;row<m;row++)rhs[row]-=dot*a[row*NB+b];}
    a[b*NB+b]=alpha;largest=fmax(largest,norm);smallest=fmin(smallest,norm);
  }
  for(b=NB-1;b>=0;b--){double value=rhs[b];for(j=b+1;j<NB;j++)value-=a[b*NB+j]*beta[j];beta[b]=value/a[b*NB+b];}
  target[0]=target[1]=0;
  for(b=0;b<NB;b++){double w0=(basis[0][b]-basis[7][b])/scales[b],w1=(basis[24][b]-basis[30][b])/scales[b];target[0]+=w0*beta[b];target[1]+=w1*beta[b];
    for(j=0;j<b;j++){w0-=a[j*NB+b]*z[0][j];w1-=a[j*NB+b]*z[1][j];}z[0][b]=w0/a[b*NB+b];z[1][b]=w1/a[b*NB+b];}
  for(h=0;h<2;h++)for(j=0;j<2;j++){cov[2*h+j]=0;for(b=0;b<NB;b++)cov[2*h+j]+=z[h][b]*z[j][b];}
  diag[1]=NB;diag[2]=smallest;diag[3]=largest;free(a);free(rhs);return 0;
}
