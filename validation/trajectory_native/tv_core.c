/* Portable double precision controller, residual operator and batch GLS.
   Plant integration lives in Python and is not linked into this DLL. */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include "sil_basis.h"
#include "tv_config.h"
#define API __declspec(dllexport)
#define NP 50
#define NB 36
static const double c[2][3]={{.7071067811865475244,-.7071067811865475244,0},{.4082482904638630164,.4082482904638630164,-.8164965809277260327}};
static const double shift[3]={0,-2.0943951023931954923,2.0943951023931954923};
static const int paths[6][2]={{0,-1},{2,5},{7,-1},{1,-1},{4,3},{6,-1}};
static int id(int q,double i){return 1-q+(i<0?3:0);}
static double signum(double x){return (x>0)-(x<0);}
static double forward(int path,double i){
  double v=0;int j,d;
  for(j=0;j<2;j++){d=paths[path][j];if(d>=0)v+=(d%2?1.15:.85)+fabs(i)*(d%2?.015:.020);}
  return signum(i)*v;
}
API void sil_control(double t,const double *sense,const int *old,int *chosen){
  int a,b,d,j,q[3];double best=1e300,amp=tv_amp*(1+tv_mod*(5*sin(2*3.141592653589793*7*(t+50e-6))+3*sin(2*3.141592653589793*13*(t+50e-6)))/16);
  for(a=-1;a<=1;a++)for(b=-1;b<=1;b++)for(d=-1;d<=1;d++){
    double u[3],mean=0,cost=0;q[0]=a;q[1]=b;q[2]=d;
    if(abs(a-old[0])>1||abs(b-old[1])>1||abs(d-old[2])>1)continue;
    for(j=0;j<3;j++){u[j]=150*q[j]-forward(id(q[j],sense[j]),sense[j])-tv_E*sin(2*3.141592653589793*25*t+shift[j]);mean+=u[j]/3;}
    for(j=0;j<3;j++){
      double pred=sense[j]+(50e-6/tv_L)*(u[j]-mean-tv_R*sense[j]);
      double ref=amp*sin(2*3.141592653589793*25*(t+50e-6)+shift[j]+tv_phase);
      cost+=(pred-ref)*(pred-ref)+.005*(q[j]-old[j])*(q[j]-old[j]);
    }
    if(cost<best){best=cost;memcpy(chosen,q,3*sizeof(int));}
  }
}
/* Per-period operator. v and r are the declared 24-device nominal drops,
   optionally conditioned on model junction temperature; never aging labels. */
API int sil_observe(double t,const double *lo,const double *hi,const int *q,
                    const int *old,const double *v,const double *r,
                    double known_delay,double *X,double *y){
  double mid[3],B[3][24]={{0}},nom[3],n[3],z[3];int p,h,j,k,lev[3];
  double weight[3]={known_delay/50e-6,.02,1-.02-known_delay/50e-6};
  memset(X,0,2*NP*sizeof(double));
  for(p=0;p<3;p++){
    mid[p]=(lo[p]+hi[p])/2;z[p]=signum(mid[p]);
    if(fabs(lo[p])<=.5||fabs(hi[p])<=.5||lo[p]*hi[p]<=0)return 0;
    lev[0]=old[p];lev[1]=mid[p]>=0?(old[p]<q[p]?old[p]:q[p]):(old[p]>q[p]?old[p]:q[p]);lev[2]=q[p];
    n[p]=abs(q[p]-lev[1]);nom[p]=-tv_E*sin(2*3.141592653589793*25*(t+25e-6)+shift[p])-tv_R*mid[p];
    for(h=0;h<3;h++){
      int path=id(lev[h],mid[p]);nom[p]+=weight[h]*150*lev[h];
      for(j=0;j<2;j++){int device=paths[path][j];if(device>=0)B[p][8*p+device]+=weight[h];}
    }
    for(j=0;j<24;j++)nom[p]-=B[p][j]*(z[p]*v[j]+mid[p]*r[j]);
  }
  for(k=0;k<2;k++){
    y[k]=0;
    for(p=0;p<3;p++){
      y[k]+=c[k][p]*((tv_L/50e-6)*(hi[p]-lo[p])-nom[p]);
      for(j=0;j<24;j++){X[k*NP+j]-=c[k][p]*z[p]*B[p][j];X[k*NP+24+j]-=c[k][p]*mid[p]*B[p][j];}
      X[k*NP+48]-=c[k][p]*n[p]*z[p];X[k*NP+49]-=c[k][p]*mid[p];
    }
  }
  return 1;
}
/* Fixed structural quotient basis from an independent separating fixture.
   Runtime whitening, column scaling and Householder QR use only this record. */
API int sil_fit(int n,const double *X,const double *y,const int *keep,double *target,double *diag){
  int k,h,j,b,d,m=0,row,prev=-2;double *a,*rhs,last=0,diagonal=(40.125*40.125+39.875*39.875)*.002*.002,off=40.125*(-39.875)*.002*.002;
  double prior[2][NB]={{0}},prior_y[2]={0},scales[NB]={0},beta[NB],largest=0,smallest=1e300;
  for(k=0;k<n;k++)if(keep[k])m+=2;
  diag[0]=m/2;if(m<NB)return -1;
  a=(double*)calloc((size_t)m*NB,sizeof(double));rhs=(double*)calloc(m,sizeof(double));if(!a||!rhs){free(a);free(rhs);return -2;}
  row=0;
  for(k=0;k<n;k++)if(keep[k]){
    double sub=k==prev+1?off/last:0;last=sqrt(diagonal-sub*sub);
    for(h=0;h<2;h++){
      for(b=0;b<NB;b++){
        double value=0;for(j=0;j<NP;j++)value+=X[(k*2+h)*NP+j]*basis[j][b];
        value=(value-sub*prior[h][b])/last;prior[h][b]=value;a[row*NB+b]=value;scales[b]+=value*value;
      }
      rhs[row]=(y[k*2+h]-sub*prior_y[h])/last;prior_y[h]=rhs[row];row++;
    }
    prev=k;
  }
  for(b=0;b<NB;b++){scales[b]=sqrt(scales[b]);if(scales[b]<1e-14){free(a);free(rhs);return -3;}for(row=0;row<m;row++)a[row*NB+b]/=scales[b];}
  for(b=0;b<NB;b++){
    double norm=0,alpha,v0,vnorm=0;
    for(row=b;row<m;row++)norm=hypot(norm,a[row*NB+b]);
    if(norm<1e-10){free(a);free(rhs);return -4;}
    alpha=a[b*NB+b]>=0?-norm:norm;v0=a[b*NB+b]-alpha;
    for(row=b+1;row<m;row++)a[row*NB+b]/=v0;
    vnorm=1;for(row=b+1;row<m;row++)vnorm+=a[row*NB+b]*a[row*NB+b];
    for(d=b+1;d<NB;d++){
      double dot=a[b*NB+d];for(row=b+1;row<m;row++)dot+=a[row*NB+b]*a[row*NB+d];dot*=2/vnorm;
      a[b*NB+d]-=dot;for(row=b+1;row<m;row++)a[row*NB+d]-=dot*a[row*NB+b];
    }
    {double dot=rhs[b];for(row=b+1;row<m;row++)dot+=a[row*NB+b]*rhs[row];dot*=2/vnorm;rhs[b]-=dot;for(row=b+1;row<m;row++)rhs[row]-=dot*a[row*NB+b];}
    a[b*NB+b]=alpha;largest=fmax(largest,norm);smallest=fmin(smallest,norm);
  }
  for(b=NB-1;b>=0;b--){double value=rhs[b];for(j=b+1;j<NB;j++)value-=a[b*NB+j]*beta[j];beta[b]=value/a[b*NB+b];}
  target[0]=target[1]=0;
  for(b=0;b<NB;b++){target[0]+=(basis[0][b]-basis[7][b])*beta[b]/scales[b];target[1]+=(basis[24][b]-basis[30][b])*beta[b]/scales[b];}
  diag[1]=NB;diag[2]=smallest;diag[3]=largest;free(a);free(rhs);return 0;
}
API void sil_timing(int repeats,double *times){
  LARGE_INTEGER freq,start,end;double sense[3]={12,-17,5};int old[3]={0,0,0},q[3],k;
  QueryPerformanceFrequency(&freq);
  for(k=0;k<repeats;k++){QueryPerformanceCounter(&start);sil_control(k*50e-6,sense,old,q);QueryPerformanceCounter(&end);times[k]=(end.QuadPart-start.QuadPart)*1e6/freq.QuadPart;memcpy(old,q,3*sizeof(int));}
}
