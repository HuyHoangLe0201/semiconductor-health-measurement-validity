/* Additional native harness; original controller/operator/ODE are unchanged.
   Noise is generated and seeded in Python. No real-time execution is claimed. */
#include "tv_core.c"
double tv_L=.002,tv_R=.25,tv_E=20,tv_amp=16,tv_mod=1,tv_phase=-.45;
API void tv_set(double L,double R,double E,double A,double m,double phase){tv_L=L;tv_R=R;tv_E=E;tv_amp=A;tv_mod=m;tv_phase=phase;}
API void sil_advance(double t,double *state,const int *levels,double h,double Tc,int aged,double slope);

/* Synchronous acquisition counterpart of sil_validation.run, with no input
   or gate-delay intervention. Preserve the four-substep integration grid. */
static int tv_refine=4,tv_replay=0;
API void tv_set_replay(int r){tv_replay=r;}
API void tv_set_refine(int r){tv_refine=r;}
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
    if(tv_replay)memcpy(q,states+3*k,3*sizeof(int));else sil_control(k*50e-6,sense,old,q);
    for(j=0;j<3;j++){states[3*k+j]=q[j];freeq[j]=state[j]>=0?(old[j]<q[j]?old[j]:q[j]):(old[j]>q[j]?old[j]:q[j]);}
    for(h=0;h<tv_refine;h++)sil_advance(k*50e-6+h*1e-6/tv_refine,state,freeq,1e-6/tv_refine,Tc,aged,1);
    for(h=0;h<tv_refine;h++)sil_advance(k*50e-6+1e-6+h*49e-6/tv_refine,state,q,49e-6/tv_refine,Tc,aged,1);
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

