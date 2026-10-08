/* Independently integrated circuit ODE. This unit knows neither the estimator
   incidence matrix nor the quotient basis and never constructs y=X theta. */
#include <math.h>
#include <string.h>
#include "tv_config.h"
#define API __declspec(dllexport)
static const double rt[4]={.08245484,.144197,.2151774,.1581708};
static const double tt[4]={7.3e-5,7e-4,.01235548,.08020881};
static const double rd[4]={.6701584,.775759,.3540826,0};
static const double td[4]={3.4e-4,4.7e-3,.04680901,1};
static void ode(double t,const double *x,const int *q,double Tc,int aged,double slope,double *out){
  double current[3],pole[3],power[24]={0},mean=0;int phase,j,m;
  for(phase=0;phase<3;phase++){
    int route[2],count=1;double drop=0,i=x[phase];current[phase]=i;
    /* Explicit directed physical routes, independent of regression assembly. */
    if(q[phase]>0){route[0]=i>=0?0:1;}
    else if(q[phase]<0){route[0]=i>=0?7:6;}
    else if(i>=0){route[0]=2;route[1]=5;count=2;}
    else{route[0]=4;route[1]=3;count=2;}
    for(j=0;j<count;j++){
      int device=8*phase+route[j],diode=device%2;double T=Tc,voltage,resistance;
      for(m=0;m<4;m++)T+=x[3+4*device+m];
      resistance=(diode?.015:.020)*slope;
      /* Piecewise-linear anchor interpolation in 25--175 C. */
      if(diode)voltage=T<=125?1.45-.0005*(T-25):1.40;
      else voltage=1.65+.002*(T-25);
      voltage-=resistance*(diode?20:40);voltage+=fabs(i)*resistance;
      if(aged&&device==0)voltage+=.05+.0002*fabs(i);
      drop+=voltage;power[device]=fabs(i)*voltage;
    }
    pole[phase]=150*q[phase]-(i>0?1:i<0?-1:0)*drop-tv_E*sin(157.07963267948966*t-(phase==0?0:phase==1?2.09439510239319549:-2.09439510239319549));
    mean+=pole[phase]/3;
  }
  for(phase=0;phase<3;phase++)out[phase]=(pole[phase]-mean-tv_R*current[phase])/tv_L;
  for(j=0;j<24;j++)for(m=0;m<4;m++)out[3+4*j+m]=((j%2?rd[m]:rt[m])*power[j]-x[3+4*j+m])/(j%2?td[m]:tt[m]);
}
static void integrate(double t,double *x,const int *q,double h,double Tc,int aged,double slope,int depth){
  double k1[99],k2[99],k3[99],k4[99],s2[99],s3[99],s4[99],end[99];int j,cross=0;
  ode(t,x,q,Tc,aged,slope,k1);for(j=0;j<99;j++)s2[j]=x[j]+h*k1[j]/2;
  ode(t+h/2,s2,q,Tc,aged,slope,k2);for(j=0;j<99;j++)s3[j]=x[j]+h*k2[j]/2;
  ode(t+h/2,s3,q,Tc,aged,slope,k3);for(j=0;j<99;j++)s4[j]=x[j]+h*k3[j];
  ode(t+h,s4,q,Tc,aged,slope,k4);for(j=0;j<99;j++)end[j]=x[j]+h*(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;
  for(j=0;j<3;j++){
    int lo=x[j]<0? -1:x[j]>0,hi=lo,values[4]={s2[j]<0?-1:s2[j]>0,s3[j]<0?-1:s3[j]>0,s4[j]<0?-1:s4[j]>0,end[j]<0?-1:end[j]>0},m;
    for(m=0;m<4;m++){if(values[m]<lo)lo=values[m];if(values[m]>hi)hi=values[m];}if(lo!=hi)cross=1;
  }
  if(cross&&depth<8){integrate(t,x,q,h/2,Tc,aged,slope,depth+1);integrate(t+h/2,x,q,h/2,Tc,aged,slope,depth+1);}
  else memcpy(x,end,sizeof(end));
}
API void sil_advance(double t,double *state,const int *levels,double h,double Tc,int aged,double slope){integrate(t,state,levels,h,Tc,aged,slope,0);}

API void tv_ode(double t,const double *x,const int *q,double Tc,double *out){ode(t,x,q,Tc,0,1,out);}
