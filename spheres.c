/* here's spheres.c; i'm sending it from raul's machine since i can't log */
/* on over the ocean from here.  it takes on the command line */

/* spheres Omega Omega_vac h n Omega_bh2 Radius z_i delta_z */

/* where Radius is the radius of the spheres in Mpc, and the rest is */
/* self-explanatory.  for example, the command line that generated the */
/* data in the figure that follows was: */

/* spheres 1 0 0.5 1 0.025 20 10 3 0.001 */

/* the output is in the file patchy.out.  the first column is ell.  the */
/* second is our calculation; the third is our calculation but doing the */
/* w integral only approximately instead of numerically; the fourth and */
/* fifth column's are the same for the GH calculation. */

/* the program still has a graceful-exit problem; it ends with a "Too */
/* many steps in routine QROMB" error, but it does so only after its */
/* gotten past the point of interest. */

/*    calculate P_perp for reionization with uncorrelated spheres.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vish.h"
#include "RECIPES/nrutil.c"
#include "RECIPES/qromb.c"
#include "RECIPES/qgaus.c"
#include "RECIPES/trapzd.c"
#include "RECIPES/polint.c"
#include "RECIPES/splint.c"
#include "RECIPES/spline.c"

#define sqr(x) ((x)*(x))
#define PI (3.141592653589793)
#define TINY (0.00000000001)
#define cube(x) ((x)*(x)*(x))
#define SMALLNUM (1.e-10)   /* for ==0.0 checks */


#define H0INV 3000.0         /* h^-1 Mpc */
#define Mpc_to_cm 3.086e24   /* cm / Mpc */
#define NHAT 1.123e-5        /* cm^-3; rho_Baryon=NHAT Omega_Bh2 */
#define SIGMA_T 6.6e-25      /* cm^2;  thomson cross-section */
#define S8_CLUS 0.60         /* sigma_8 Omega_0^0.53 from clusters */
#define T0_muK  2.73e6       /* T_0 / muK */
#define DUMMYK 600.

static double Omega,Omega_vac,hhh,nnn,Omega_bh2;
static double dummyk,Radius,ell;
static double *wz_array,*w_array,*w_2,*D_array,*D_2,*z_2;
static double *k_array,*Sk_array,*Sk_2;
              /* these arrays will be required for the prefactor
		 and possibly for future more general visibility fns */
static double D_0;

static double w_low,w_high;    /* w limits of patchiness */

static double xe_bar(double w);
static double w_integrand(double w);
static double ggg(double w);
static double Cl_integrand(double w);
static double Cl_GHintegrand(double w);
static double S_GH(double k);
static double P_integrand(double u);
static double P_integral;
static double sigma8_integrand(double u);
static double v2_integrand(double u);

main(int argc,char *argv[])
{
  double kph,k,lognumber,x,y,z=0.,z_i,delta_z,w;
  double delta_H,prefactor,normalization,w_r,Delta_w,tau_r,Cl,logell;
  double Cl_approx;
  double z_r,D_r,a_r,addot_over_adot,adot_over_a,Ddot_over_D,ClGH;
  double sigma8,v2,dummykph,u,w_i,GH,numerical;
  double aDdot_over_D;
  
  /*
    double S_interp(),S(),integrand(),P(),E(),www(),zzz(),
    wintegrand(),qromb(),Dintegrand(),D(),sigma8_integrand(),uintegrand(); */
  int n;
  FILE *zoutfile,*outfile,*Soutfile,*xoutfile,*Poutfile,*youtfile,*patchyfile;

  Omega=atof(argv[1]);
  Omega_vac=atof(argv[2]);
  hhh=atof(argv[3]);
  nnn=atof(argv[4]);
  Omega_bh2=atof(argv[5]);
  Radius=atof(argv[6]);       /* comoving Radius of sphere in Mpc  */
  z_i=atof(argv[7]);
  delta_z=atof(argv[8]);

  wz_array=vector(1,1100);
  w_array=vector(1,1100);
  w_2=vector(1,1100);
  z_2=vector(1,1100);
  D_array=vector(1,1100);
  D_2=vector(1,1100);
  k_array=vector(1,281);
  Sk_array=vector(1,281);
  Sk_2=vector(1,281);

  /*  diagnostic to test P_integrand
  dummyk=DUMMYK;
  for(u= -1.-log10(DUMMYK);u<=7.-log10(DUMMYK);u+=0.01) {
    printf("%e  %e\n",u,P_integrand(u));
  }
  */

  dummyk=DUMMYK;
  P_integral=log(10.)*qromb(P_integrand,-1.-log10(DUMMYK),7.-log10(DUMMYK));

  /*  set up lookup table to evaluate www(z) and zzz(w) */
  wz_array[1]=0.;  w_array[1]=0.;
  for(n=2;n<=1100;n++) {
    z+=1.;
    wz_array[n]=z;
    w_array[n]=qromb(wintegrand,0.,z);
  }
  spline(wz_array,w_array,111,1.e40,1.e40,w_2);
  spline(w_array,wz_array,111,1.e40,1.e40,z_2);


  /*  set up lookup table to evaluate D(z)  */
  for(n=1;n<=1100;n++) {
    z=wz_array[n];
    D_array[n]=5.*Omega/2.*E(z)*qromb(Dintegrand,1.e-5,1./(1.+z));
  }
  spline(wz_array,D_array,1100,1.e40,1.e40,D_2);

  D_0=D(0.);

 /*  print file of z, w(z), and D(z) as a diagnostic  */
  zoutfile=fopen("z.out","w");
  for(z=0.;z<=1100.;z+=1.) {
    fprintf(zoutfile,"%f  %f  %f\n",z,www(z),D(z)/D_0);
  }
  fclose(zoutfile);

  w_low=www(z_i-0.5*delta_z);  w_high=www(z_i+0.5*delta_z);
  w_i=www(z_i);

  /*   debugging routine to check the power spectrum  */
  Poutfile=fopen("P.out","w");
  for(lognumber= -3.;lognumber <= 0.; lognumber += 0.1) {
    kph=pow(10.,lognumber);
    fprintf(Poutfile,"%e  %e\n",kph,P(kph*hhh));
  }
  fclose(Poutfile);


  /*  debugging routine to check the y integrand
  youtfile=fopen("y.out","w");
  dummyk=dummykph*6000.;
  for(u= -2.-log10(dummyk);u<=6.-log10(dummyk);u+=0.01) {
    fprintf(youtfile,"%e  %e  %e\n",u,dummykph*u,uintegrand(u));
  }
  fclose(youtfile);  */

  /* set up lookup table to evaluate S(k)  */
  Soutfile=fopen("S.out","w");
  n=0;
  for(lognumber=-5.; lognumber<=2.; lognumber += 0.025) {
    n++;
    kph = pow(10.,lognumber);
    k=6000.*kph;
    k_array[n]=k;   Sk_array[n]=S(k);
  }
  spline(k_array,Sk_array,281,1.e40,1.e40,Sk_2);
  
  /*  print out S(k) and S_interp(k) as diagnostic  */
  for(lognumber= -4.83; lognumber<=1.96; lognumber += 0.01) {
    kph = pow(10.,lognumber);
    k=6000.*kph;
    fprintf(Soutfile,"%e  %e  %e  %e\n",
	    hhh*kph*Radius,k*S(k),k*S_interp(k),k*S_GH(k));
  }
  fclose(Soutfile);

  /*  here we include the normalization of the power spectrum  */
  delta_H=1.94e-5 * pow(Omega,-0.785-0.05*log(Omega))
               * exp((nnn-1.)+1.97*sqr(nnn-1.));
  /*  this is the cobe power-spectrum normalization from bunn&white */
  
  normalization = 2.*PI*PI*sqr(delta_H)/8.;


  addot_over_adot = (hhh/H0INV)*(Omega_vac-Omega/2.);
  adot_over_a = (hhh/H0INV);
  Ddot_over_D = addot_over_adot - adot_over_a 
    + 5.*Omega/2.*adot_over_a  / D_0;
  aDdot_over_D=2.*(H0INV/hhh)*Ddot_over_D;

  
  /*  check v2_integrand
  for(u=log10(6000*0.001);u<=log10(6000.*100.);u+=0.1) {
    printf("%e  %e\n",u,v2_integrand(u));
  }
  exit(1);
  */

  /*  calculate sigma_8  and v2
  sigma8=9./2./PI/PI * log(10.)*qromb(sigma8_integrand,log10(6000*0.001),
				      log10(6000.*100.));
  sigma8 *= 2.*PI*PI*sqr(delta_H)/8.;

  dummyk=DUMMYK;
  v2= 1./2./PI/PI * log(10.) * DUMMYK
   * qromb(P_integrand,-1.-log10(DUMMYK),7.-log10(DUMMYK));
  v2 *= 2.*PI*PI*sqr(delta_H)/8. * sqr(aDdot_over_D);

  v2= 1./2./PI/PI * log(10.) * qromb(v2_integrand,log10(6000*0.001),
				      log10(6000.*100.));
  v2 *= 2.*PI*PI*sqr(delta_H)/8. * sqr(aDdot_over_D);

  printf("sigma8=%e   <v2>^{1/2}=%e\n",sqrt(sigma8),sqrt(v2));
  exit(1);  */

  /*  print Cl_GHintegrand as diagnostic
  for(w=w_low;w<=w_high;w+=(w_high-w_low)/100.) {
    printf("%e  %e\n",w,Cl_GHintegrand(w));
  }
  */

  /*  Cl=0.5*cube(Radius)/sqrt(PI)*normalization
       *qromb(Cl_integrand,w_low,w_high);
  printf("%e  %e\n",ell,sqr(ell)*Cl);  */
  
  /*  ell=1000.;
  for(w=w_low;w<w_high;w+=(w_high-w_low)/100.) {
    printf("%e  %e\n",w,w_integrand(w)/(w_high-w_low));
  }
  exit(1);  */

  /*   now we do the calculation of the C_l's  */
  patchyfile=fopen("patchy.out","w");
  for(logell=1.;logell<=4.;logell+=0.1) {
    ell=pow(10.,logell);
    Cl=cube(Radius*hhh/2./H0INV)*pow(2.*PI,1.5)/4./sqr(PI)*normalization
      *qromb(Cl_integrand,w_low,w_high);
    Cl_approx=cube(Radius*hhh/2./H0INV)*pow(2.*PI,1.5)/4./sqr(PI)*normalization
	*Cl_integrand(w_i)*(w_high-w_low);
    ClGH=cube(Radius*hhh/2./H0INV)*pow(2.*PI,1.5)/sqr(PI)*normalization
      *P_integral*DUMMYK*qromb(Cl_GHintegrand,w_low,w_high);
    GH=pow(1.-w_i,-6.)/w_i/w_i/6.
      *exp(-sqr(ell*Radius*hhh/w/6000.)/2.)*(w_high-w_low);
    numerical=qromb(w_integrand,w_low,w_high);
    fprintf(patchyfile,"%e  %e  %e  %e  %e\n",
	    ell,sqr(ell)*Cl/2./PI,sqr(ell)*Cl_approx/2./PI,sqr(ell)*ClGH/2./PI,
	    sqr(ell)*ClGH/2./PI*GH/numerical);
    /*  printf("%e  %e  %e\n",ell,numerical,GH);  */
  }
  fclose(patchyfile);
}


double S(double k)
{
  dummyk=k;

  if(k<30000./Radius/hhh) {
    return k*log(10.)*qromb(uintegrand,-2.-log10(k),6.-log10(k));
    /* to be conservative, the upper limit here is probably a bit larger
       than it needs to be.  its probably lower by an order of magnitude 
       (for EdS) for k=1.e6, and lower by a factor of 100 for k=10.  */
  /* we integrate over u=log10(y) instead of y since the integrand
     is smoother in u than in y  */
  }
  else {
    return 0.;
  }
}


double S_GH(double k)
{
  double kp=hhh*k/2.0/3000.;
  
  return (4./3.)*exp(-sqr(kp*Radius)/2.)*P_integral*(DUMMYK);
}



double uintegrand(double u)
{
  double kp,alpha,y,thing;

  y=pow(10.,u);

  kp=hhh*dummyk/2.0/3000.;

  alpha=sqr(kp*Radius)*y;

  thing=sqr(kp*Radius)/2.;

  /*  return y * P(kp*y) * exp(-alpha*y-sqr(kp*Radius)) * 
    (alpha*cosh(alpha)-sinh(alpha))/alpha/alpha/alpha;   */
  /*  printf("%e  %e  %e\n",y,alpha-thing,alpha*y);  */
  /*  return y * P(kp*y) * 
    (alpha*0.5*(exp(alpha-thing-alpha*y)+exp(-alpha-thing-alpha*y))
     -0.5*(exp(alpha-thing-alpha*y)
	   -exp(-alpha-thing-alpha*y)))/alpha/alpha/alpha; */

  /*  return y*P(kp*y)*exp(-sqr(kp*Radius)/2.)*4./3.;  for GH case */

  if(alpha>0.01) {
   return 2.* y * P(kp*y) *(exp(alpha-thing-thing*y*y)*(alpha-1.)
      +exp(-alpha-thing-thing*y*y)*(alpha+1.))/alpha/alpha/alpha;
  }
  else {
    return y * P(kp*y) * exp(-thing-thing*y*y)*4./3.;
  }
}


double P(double kp)   /* current power spectrum (from bardeen et al.)  */
{
  double thing,transfer,k,q,aa,ab,al,q2,t2;
  
  k=2.0*kp*3000/hhh;
  /*  kp *= 1./Omega/sqr(hhh); */    /*  now scaled by Omega h^2 */
  /*  thing=pow( 6.4*kp + 3.0*kp*sqrt(3.0*kp) + sqr(1.7*kp) , 1.13);
  transfer=pow(1.+ thing, -1./1.13);
  */
  

  aa=pow(46.9*Omega*sqr(hhh),0.67)*(1+pow(32.1*Omega*sqr(hhh),-0.532));
  ab=pow(12.0*Omega*sqr(hhh),0.424)*(1+pow(45.*Omega*sqr(hhh),-0.582));
  al=pow(aa,-Omega_bh2/sqr(hhh)/Omega)*pow(ab,-cube(-Omega_bh2/sqr(hhh)/Omega));
    q2=kp/Omega/sqr(hhh)*exp(Omega_bh2/sqr(hhh)*(1.+1./Omega));

  q=kp/Omega/sqr(hhh)/sqrt(al)*sqr(1.0104)/sqrt(1-Omega_bh2/sqr(hhh)/Omega);
  transfer=log(1.+2.34*q)/(2.34*q)
    *pow(1.+3.89*q+sqr(16.1*q)+cube(5.46*q)+sqr(sqr(6.71*q)),-0.25);
  /*  printf("\n al=%f   q=%f  q2=%f a1=%f a2=%f\n",al,q/kp,q2/kp,aa,ab);*/

  return pow(k/2.,nnn)*sqr(transfer);
}


double wintegrand(double z)
{
  return 0.5/E(z);
}

double www(double z)    /*  calculate w(z) by cubic interpolation  */
{
  double y;
  
  splint(wz_array,w_array,w_2,1100,z,&y);
  return y;
}


double zzz(double w)   /* calculate z(w) by cubic interpolation  */
{
  double y;
  
  splint(w_array,wz_array,z_2,1100,w,&y);
  return y;
}


double E(double z)
{
  return sqrt(Omega*(1.+z)*(1.+z)*(1.+z)+Omega_vac 
	      +(1.-Omega-Omega_vac)*sqr(1.+z));
}


double D(double z)   /* calculate D(z) by cubic interpolation  */
{
  double y;
  
  splint(wz_array,D_array,D_2,1100,z,&y);
  return y;
}


double Dintegrand(double alpha)
{
  double EEE;

  EEE=E(1./alpha-1.);

  return 1./cube(alpha)/cube(EEE);
}
  
double sigma8_integrand(double u)
{
  double k=pow(10.,u),window;
  double kph=k/6000.;

  window=j1(kph*8.)/(kph*8.);
  return k*k*k*sqr(window)*P(kph*hhh);
}


double v2_integrand(double u)
{
  double k=pow(10.,u);
  double kph=k/6000.;

  return k*P(kph*hhh);
}


double j1(double x)
{
  return sin(x)/x/x - cos(x)/x;
}


double S_interp(double k)     /*  evaluate S(k) by interpolation  */
{
  double y;
  

  if(k<30000./Radius/hhh) {
    splint(k_array,Sk_array,Sk_2,281,k,&y);
    return y;
  }
  else {
    return 0.;
  }
}


double Cl_integrand(double w)
{
  double g=ggg(w),xe=xe_bar(w),z=zzz(w),RS_chi;

  double addot_over_adot,adot_over_a,Ddot_over_D,aDdot_over_D;
  
  addot_over_adot = (hhh/H0INV)*(Omega_vac-Omega*cube(1.+z)/2.)
    /E(z);
  adot_over_a = (hhh/H0INV)*E(z);
  Ddot_over_D = addot_over_adot - adot_over_a 
    + 5.*Omega/2.*adot_over_a * sqr((1.+z)/E(z)) / D(z);
  aDdot_over_D=2.*(H0INV/hhh)/(1.+z)*Ddot_over_D*D(z)/D_0;   /* changed */

  if(fabs(1.-Omega-Omega_vac)<SMALLNUM) {
    RS_chi=w;
  }
  else {
    if(1.-Omega-Omega_vac>0.) {
      RS_chi=0.5/sqrt(1.-Omega-Omega_vac) 
	       * sinh(2.0*sqrt(1.-Omega-Omega_vac)*w);
    }
    else {
      RS_chi=0.5/sqrt(-1.+Omega+Omega_vac) 
	       * sin(2.0*sqrt(-1.+Omega+Omega_vac)*w);
    }
  }

/*  return sqr(g/RS_chi)*(1-xe)*xe*sqr(aDdot_over_D)
    *(ell/RS_chi)*S_interp(ell/RS_chi);  */
  return sqr(g/RS_chi)*(1-xe)*xe*sqr(aDdot_over_D)
    *S_interp(ell/RS_chi);
}


double Cl_GHintegrand(double w)
{
  double g=ggg(w),xe=xe_bar(w),z=zzz(w),RS_chi;

  double addot_over_adot,adot_over_a,Ddot_over_D,aDdot_over_D;

  addot_over_adot = (hhh/H0INV)*(Omega_vac-Omega*cube(1.+z)/2.)
    /E(z);
  adot_over_a = (hhh/H0INV)*E(z);
  Ddot_over_D = addot_over_adot - adot_over_a 
    + 5.*Omega/2.*adot_over_a * sqr((1.+z)/E(z)) / D(z);
  aDdot_over_D=2.*(H0INV/hhh)/(1.+z)*Ddot_over_D*D(z)/D_0;  /* changed */

  if(fabs(1.-Omega-Omega_vac)<SMALLNUM) {
    RS_chi=w;
  }
  else {
    if(1.-Omega-Omega_vac>0.) {
      RS_chi=0.5/sqrt(1.-Omega-Omega_vac) 
	       * sinh(2.0*sqrt(1.-Omega-Omega_vac)*w);
    }
    else {
      RS_chi=0.5/sqrt(-1.+Omega+Omega_vac) 
	       * sin(2.0*sqrt(-1.+Omega+Omega_vac)*w);
    }
  }

  return (1./6.)*sqr(g/RS_chi)*(1-xe)*xe*sqr(aDdot_over_D)
    *exp(-sqr(ell*Radius*hhh/w/6000.)/2.);
}


double ggg(double w)
{ 
  double tau_prefactor=(NHAT*Omega_bh2*SIGMA_T*H0INV*Mpc_to_cm/hhh);

  double z=zzz(w);
  return 2.*tau_prefactor*sqr(1.+z);
}


double xe_bar(double w)
{
  if(w<w_low || w>w_high) {
    return 0.;
  }
  else {
    return 1.-(w-w_low)/(w_high-w_low);
  }
}
    

double P_integrand(double u)
{
  double kp,alpha,y,thing;

  y=pow(10.,u);

  kp=hhh*dummyk/2.0/3000.;

  return y * P(kp*y);
}


double w_integrand(double w)
{
  return pow(1.-w,-6.)/w/w*xe_bar(w)*(1.-xe_bar(w))
    *exp(-sqr(ell*Radius*hhh/w/6000.)/2.);
}


