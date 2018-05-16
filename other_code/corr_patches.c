/* From: Marc Kamionkowski <kamion@cuphy3.phys.columbia.edu> */
/* To: jaffe@physics.Berkeley.EDU */
/* Subject: corr_patches.c */
/* Date: Mon, 1 Jun 1998 17:40:04 -0400 */

/* here's the code that calculates that term.  i looked through and  */
/* found no egregious normalization errors.  do you want to take a */
/* look?  it operates exactly like vish.c but changes that quantity  */
/* in the integrand. */

/*    calculate Cl's from correlated patches.
      take Omega_0 hhh nnn Omegabh2 Omega_vac z_i delta_z on the command line

      assume reionization occurs at z_i with width delta_z and
      xe_bar goes from 0 to 1 in that interval.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vish.h"
#include "/home/kamion/recipes_c/nrutil.c"
#include "/home/kamion/recipes_c/qromb.c"
#include "/home/kamion/recipes_c/qgaus.c"
#include "/home/kamion/recipes_c/gauleg.c"
#include "/home/kamion/recipes_c/trapzd.c"
#include "/home/kamion/recipes_c/polint.c"
#include "/home/kamion/recipes_c/splint.c"
#include "/home/kamion/recipes_c/spline.c"

#define ZMAX 1100        /* surface of last scattering, more or less */
#define DELZT 0.25       /* tau(z), g(z) calculated with this interval */
#define DELZ  2.0        /* w(z) calculated with this interval */
#define NARR  ((int)floor(ZMAX/DELZ))   /* size of arrays */
#define NARRtau ((int)ceil(ZMAX/DELZT)) /* >(NARR*DELZ/DELZT)
                                           big so we don't fall off the end */
#define NARRPk 402

#define H0INV 3000.0         /* h^-1 Mpc */
#define Mpc_to_cm 3.086e24   /* cm / Mpc */
#define NHAT 1.123e-5        /* cm^-3; rho_Baryon=NHAT Omega_Bh2 */
#define SIGMA_T 6.6e-25      /* cm^2;  thomson cross-section */
#define S8_CLUS 0.60         /* sigma_8 Omega_0^0.53 from clusters */
#define T0_muK  2.73e6       /* T_0 / muK */

/* switches! */
/* #define OMEGA_EQ_1 */
/* #define USE_S_APPROX */
/* #define SQRT_SING */   /* deal with sqrt singularity of xintegrand */

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define PI M_PI             /* (3.141592653589793) */
#define LOG10 M_LN10        /* (log(10.0)) */
#define SMALLNUM (1.e-10)   /* for ==0.0 checks */

static double patch_xintegrand(double w);  /* xintegrand for corr patches */
static double patch_Pintegrand(double w);  /* P integrand for corr patches */
static double patch_uintegrand(double w);  /* u integrand for corr patches */
static double patch_S(double k);    /*  S(k) for corr patches */
static double patch_S_interp(double k);    /*  S(k) for corr patches */
static double patch_p_proj(double k);    /*  p_proj for corr patches */
static double ggg(double w);             /*  visibility function  */
static double xe_bar(double w);    /*  ionization fraction as fn of time  */
static double Omega,hhh,nnn,Omega_Bh2,Omega_k, Omega_vac;
static double dummyk,dummyy, kappa_g,w_low,w_high,w_i;
static double *wz_array,*w_array,*w_2,*D_array,*D_2,*z_2, *logT2_array, *Pk_2;
static double *g_array, *g_2, *tau_array, *tau_2, *tauz_array, *logk_array;
static double *logpatchS_array, *logpatchS_logk_array, *logpatchS_2;
static double D_0, ETA0;
static double alpha_g,  zmax_g,  x0_g,  xbar_g, tau_r;
static int narrpk, narrtau, do_read_Pk;
static double tau_prefactor;

int main(int argc, char *argv[])
{
    double kph,k,lognumber,x;  /* y*/
    double z, zt, w, RS_chi;
    double delta_H,prefactor,normalization,w_r,Delta_w,Cl,ell;
    double z_r,D_r,a_r,addot_over_adot,adot_over_a,Ddot_over_D;
    double sigma82, gg, z_i, delta_z;
    double z_m, x_0, xbar, alpha, g_prefactor;
    double DT_T2, dlog10ell, llCl, llCl_prev=0.0, deriv, deriv_prev=0.0, gsum=0.0;
    int n, i, do_Gauss_vis, findzm=0;
    char Pkfilename[120];
    FILE *zoutfile,*outfile,*Soutfile,*xoutfile,*Poutfile, *Pkfile;/* *youtfile*/
    
    Omega=atof(argv[1]);
    hhh=atof(argv[2]);
    nnn=atof(argv[3]);
    Omega_Bh2=atof(argv[4]);
    Omega_vac=atof(argv[5]);
    z_i=atof(argv[6]);
    delta_z=atof(argv[7]);
    Omega_k=1.0-Omega-Omega_vac;

    
    printf("Omega=%f, Omega_k=%f, Omega_vac=%f\n", Omega, Omega_k, Omega_vac);
    printf("hhh=%f\n", hhh);
    printf("nnn=%f\n", nnn);
    printf("Omega_Bh2=%f\n", Omega_Bh2);
    printf("z_i=%f   delta_z=%f\n",z_i,delta_z);

    /* get vars the easy way */
    printf("\nRead T(k) [>0 = yes]? ");
    scanf("%d", &do_read_Pk);
    if (do_read_Pk>0) printf("\nP(k) filename: "); else printf("\nDummy: ");
    scanf("%s", Pkfilename);
    if (do_read_Pk>0) Pkfile=fopen(Pkfilename, "r");

    
    logk_array=vector(1,NARRPk);
    logT2_array=vector(1,NARRPk);
    Pk_2=vector(1,NARRPk);
    logpatchS_logk_array=vector(1,NARRPk);
    logpatchS_array=vector(1,NARRPk);
    logpatchS_2=vector(1,NARRPk);
    
    g_array=vector(1,NARRtau);
    g_2=vector(1,NARRtau);
    tau_array=vector(1,NARRtau);
    tau_2=vector(1,NARRtau);
    tauz_array=vector(1,NARRtau);
    
    wz_array=vector(1,NARR);
    w_array=vector(1,NARR);
    w_2=vector(1,NARR);
    z_2=vector(1,NARR);
    D_array=vector(1,NARR);
    D_2=vector(1,NARR);
    
    /* get vars the easy way */

    tau_prefactor=(NHAT*Omega_Bh2*SIGMA_T*H0INV*Mpc_to_cm/hhh);

    /* read in P(k) & set up lookup table to evaluate */
    if (do_read_Pk>0) {
        narrpk=read_trans(Pkfile, logk_array, logT2_array, NARRPk);
        printf("narrpk=%d\n", narrpk);
        printf("min logk P(k)=%f %f\n", logk_array[1], logT2_array[1]);
        printf("max logk P(k)=%f %f\n", logk_array[narrpk], logT2_array[narrpk]);
        /* known slopes */
        spline(logk_array, logT2_array, narrpk, 0.0, -4.0, Pk_2); 
    }

    /*  set up lookup table to evaluate www(z), zzz(w) & tau(z) */
    /*  nb. tabulate www(z) at dz=10; tau(z) at dz=0.25 */
    printf("Setting www, zzz... "); fflush(stdout);
    z=0.0;
    wz_array[1] = w_array[1] = 0.;
    for(n=2;n<=NARR;n++) {
        z+=DELZ;
        wz_array[n]=z;
        w_array[n]=qromb(wintegrand,0.,z);
    }
    spline(wz_array,w_array,NARR,1.e40,1.e40,w_2);
    spline(w_array,wz_array,NARR,1.e40,1.e40,z_2);
    
    printf("logpatchS(logk)... "); fflush(stdout);
    for(i=1; i<=NARRPk; i++) {
      lognumber=-3.0+(i-1)*(7.0/(NARRPk-1));
      kph=pow(10.0, lognumber);
      k=2.0*H0INV*kph;
      logpatchS_logk_array[i]=log10(k);
      logpatchS_array[i]=log10(patch_S(k));
    } /* known slopes */
    spline(logpatchS_logk_array, logpatchS_array, NARRPk, 0.0, -3.0, 
	   logpatchS_2);
    
    /*  set up lookup table to evaluate D(z)  */
    printf("D(z)\n"); fflush(stdout);
    for(n=1;n<=NARR;n++) {
        z=wz_array[n];
        D_array[n]=5.*Omega/2.*E(z)*qromb(Dintegrand,1.e-5,1./(1.+z));
    }
    spline(wz_array,D_array,NARR,1.e40,1.e40,D_2);

    w_low=www(z_i-0.5*delta_z);  w_high=www(z_i+0.5*delta_z);
    w_i=www(z_i);

    D_0=D(0.);
    ETA0=www(wz_array[NARR]);
    printf("D_0=%f; eta_0=%f\n", D_0, ETA0); 
    fflush(stdout);

    /*  print file of z, w(z), and D(z) as a diagnostic
    printf("writing z.out... "); fflush(stdout);
    zoutfile=fopen("z.out","w");
    for(z=0.;z<=1100.;z+=1.) {
        fprintf(zoutfile,"%f  %f  %f  %f\n",z,www(z),zzz(www(z)),D(z)/D_0);
    }
    fclose(zoutfile);
    */

    /*   debugging routine to check the power spectrum
    printf("P.out... "); fflush(stdout); 
    Poutfile=fopen("P.out","w"); 
    for(lognumber= -3.;lognumber <= 0.; lognumber += 0.1) { 
        kph=pow(10.,lognumber); 
        fprintf(Poutfile,"%f %f %g\n",kph, kph*hhh, P(kph*hhh)); 
    } 
    fclose(Poutfile);   */


     /*  debugging routine to check the x integrand
    printf("x.out... "); fflush(stdout);  
    xoutfile=fopen("x.out","w"); 
    for(x=-1.;x<=1.;x+=0.001) { 
        dummyk = 6.e3;  dummyy=1.; 
        fprintf(xoutfile,"%f  %g\n",x,patch_xintegrand(x)); 
    }
    fclose(xoutfile);  */

    /*  debugging routine to check the y integrand; not currently operational
        youtfile=fopen("y.out","w");
        dummyk=2.0*H0INV*10.;
        for(lognumber=6.-log10(dummy)-6.;lognumbery=1.e-20; y<= 1.e6/dummyk; 
        y+=1.e6/dummyk/1000.) {
        fprintf(youtfile,"%e  %e\n",y,yintegrand(y));
        }
        fclose(youtfile);
        */

    /* calculate patch_S(k) as a diagnostic */
    printf("patch_S.out\n"); fflush(stdout); 
    Soutfile=fopen("patch_S.out","w");
    for(lognumber=-6.; lognumber<=6.; lognumber += 0.1) {
        kph = pow(10.,lognumber);
        k=2.0*H0INV*kph;
        fprintf(Soutfile,"%e %e %e %e\n",kph,patch_S(k),patch_S_interp(k));
    }
    fclose(Soutfile);

    /*  here we include the normalization of the power spectrum  */
    if (fabs(Omega_k)<SMALLNUM) {
        delta_H=1.94e-5 * pow(Omega,-0.785-0.05*log(Omega))
            * exp((nnn-1.)+1.97*sqr(nnn-1.));
    } else if (fabs(Omega_vac)<SMALLNUM) {
        delta_H=1.95e-5* pow(Omega,-0.35-0.19*log(Omega)-0.17*(nnn-1.0))
            * exp(-(nnn-1.0)-0.14*sqr(nnn-1.0));
    } else if (nnn==1.0) {
        delta_H=2.422 - 1.166*exp(Omega) + 0.800*exp(Omega_vac)
            + 3.780*Omega - 2.267*Omega*exp(Omega_vac) + 0.487*sqr(Omega)
            + 0.561*Omega_vac + 3.392*Omega_vac * exp(Omega) 
            - 8.568*Omega * Omega_vac + 1.080 * sqr(Omega_vac);
    } else {
        printf("Can't COBE-normalize this spectrum!\n");
    }
        /*  this is the cobe power-spectrum normalization from bunn&white */
    normalization = sqr(2.*PI*PI*sqr(delta_H)/8.);

    /*  print results to vishniac.out  */
    printf("writing corr_patches.out...\n"); fflush(stdout); 
    outfile=fopen("corr_patches.out","w");


    DT_T2=0.0;
    dlog10ell=0.1;
    for(lognumber=1.;lognumber<=6;lognumber+=dlog10ell) {
        ell=pow(10.,lognumber);
	Cl=g_prefactor*normalization*patch_p_proj(ell);
        /* find the maximum by keeping track of the sign of the derivative */
        llCl=ell*ell*Cl;
        deriv=llCl-llCl_prev;
        if (deriv<0.0 && deriv_prev>0.0) {
            printf("deriv ~ 0 at ell=%f\n", ell);
        }
        deriv_prev=deriv; llCl_prev=llCl;
        DT_T2 += llCl;
        
        fprintf(outfile,"%e  %e\n",ell,llCl);
    }
    fclose(outfile);
    /* calculate <DT/T^2> = \sum_l (l C_l)/(2pi) */
    DT_T2 *= dlog10ell*LOG10/(2.0*PI);

    /*  calculate sigma_8  */
    sigma82=9./2./PI/PI *LOG10*qromb(sigma8_integrand,log10(2.0*H0INV*0.00001),
                                       log10(2.0*H0INV*10000.));
    sigma82 *= 2.*PI*PI*sqr(delta_H)/8.;

    printf("   <DT/T^2>=%g  rms DT/T=%g, clus. norm. rms DT/T=%g\n",
           DT_T2, sqrt(DT_T2), 
           sqr(S8_CLUS)*pow(Omega,-1.06)*sqrt(DT_T2)/sigma82);

    printf("muK: <DT^2>=%g  rms   DT=%g, clus. norm. rms  DT =%g\n",
           DT_T2*sqr(T0_muK), sqrt(DT_T2)*T0_muK, 
           sqr(S8_CLUS)*pow(Omega,-1.06)*sqrt(DT_T2)*T0_muK/sigma82);
    printf("sigma8=%e\n",sqrt(sigma82));

    return 0;
}
  

#define YYY 0.5
double patch_p_proj(double kappa)
{
    static double prefac_loc=1./(16.0*PI*PI);
    kappa_g=kappa;
    return prefac_loc*qgaus_ahj(patch_Pintegrand, w_low, w_high);
}


double patch_Pintegrand(double ww) /* y=eta/eta_0 */
{
    double RS_chi, gg, SS;
    double zz=zzz(ww);  
    double DD=D(zz);
    double aa = (2.*H0INV/hhh)/(1.+zz);
    double addot_over_adot = (hhh/H0INV)*(Omega_vac-Omega*cube(1.+zz)/2.)
        /E(zz);
    double adot_over_a = (hhh/H0INV)*E(zz);
    double Ddot_over_D = addot_over_adot - adot_over_a 
        + 5.*Omega/2.*adot_over_a * sqr((1.+zz)/E(zz)) / DD;
    /* printf("y,w,z,D,a=%f %f %f %f %f\n",y,ww,zz,DD,aa); */
    if (fabs(Omega_k)<SMALLNUM) RS_chi=ww;
    else RS_chi=0.5/sqrt(Omega_k) * sinh(2.0*sqrt(Omega_k)*ww);
    SS=patch_S_interp(kappa_g/RS_chi);
    gg=ggg(ww)*xe_bar(ww);
    return 2.*(1.-xe_bar(ww))*sqr(aa * Ddot_over_D * (DD/D_0)*gg/RS_chi)*SS;
}


double ggg(double w)
{
  double tau_prefactor=(NHAT*Omega_Bh2*SIGMA_T*H0INV*Mpc_to_cm/hhh);

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



double patch_S_interp(double k)   /* current power spectrum   */
     /*  kp = physical wave number [Mpc^-1] */
{
    double logk, SS;
    
    logk=log10(k); 
    if (logk <= logpatchS_logk_array[1]) 
        SS=logpatchS_array[1];
    else if (logk >= logpatchS_logk_array[NARRPk]) {
      /* extrapolate with k^-3 */
        SS=logpatchS_array[NARRPk]-3.0*(logk-logpatchS_logk_array[NARRPk]);
    }
    else {
        splint(logpatchS_logk_array,logpatchS_array,logpatchS_2,NARRPk,
	       logk,&SS);
    }
    return pow(10.0,SS);
}

double patch_S(double k)
{
    double lk=log10(k), dn=-3.0-lk, up=6.0-lk;
    dummyk=k;
    /* return k*LOG10*qromb(patch_uintegrand,-2.-lk,6.-lk); */
    if (dn<0.0 && up>0.0) /* integrable singularity at y=1? */
        return k*LOG10*(qromb(patch_uintegrand, dn, 0.0) + 
                        qromb(patch_uintegrand, 0.0, up));
    else 
        return  k*LOG10*(qromb(patch_uintegrand, dn, up));
    /* to be conservative, the upper limit here is probably a bit larger 
       than it needs to be.  its probably lower by an order of magnitude 
       (for EdS) for k=1.e6, and lower by a factor of 100 for k=10.  */

    /* we integrate over u=log10(y) instead of y since the integrand
       is smoother in u than in y  */
}

double patch_xintegrand(double x)
{
    double kp;

    kp=hhh*dummyk/2.0/H0INV;
  
    return P(kp*sqrt(1.+sqr(dummyy)-2.*x*dummyy)) 
      * (1.-x*x);
}

double patch_uintegrand(double u)
{
    double kp,y;

    y=pow(10.,u);
    kp=hhh*dummyk/2.0/H0INV;
    dummyy=y;

    return y * P(kp*y) * qgaus_ahj(patch_xintegrand,-1.,1.);
}


double wintegrand(double z)
{
    return 0.5/E(z);
}

double www(double z)    /*  calculate w(z) by cubic interpolation  */
{
    double y;
  
    splint(wz_array,w_array,w_2,NARR,z,&y);
    return y;
}


double zzz(double w)   /* calculate z(w) by cubic interpolation  */
{
    double y;
  
    splint(w_array,wz_array,z_2,NARR,w,&y);
    return y;
}


double E(double z)
{
    return sqrt(Omega*cube(1.+z) + Omega_k*sqr(1.0+z) + Omega_vac);
}


double D(double z)   /* calculate D(z) by cubic interpolation  */
{
    double y;
  
    splint(wz_array,D_array,D_2,NARR,z,&y);
    return y;
}


double Dintegrand(double alpha)
{
    double EEE;

    EEE=E(1./alpha-1.);
    return 1./cube(alpha)/cube(EEE);
}
  
double sigma8_integrand(u)
     double u;
{
    double k=pow(10.,u),window;
    double kph=k/(2.0*H0INV);

    window=jj1(kph*8.)/(kph*8.);
    return k*k*k*sqr(window)*P(kph*hhh);
}

double jj1(double x)  /* careful; math.h has j1()=J_1 */
{
    return sin(x)/x/x - cos(x)/x;  
    /* should check for x=0? [j1(x)=x/3 + O(x^3)] */
}

double P_cdm(double kp)   /* current power spectrum (from bardeen et al.)  */
     /* kp= physical wave number  */
{
    double transfer,k,q;
    static double aa, ab, al;
    static int first=1;
    /*  double thing,t2,q2 */
    if (first==1) {
        aa=pow(46.9*Omega*sqr(hhh),0.67)*(1+pow(32.1*Omega*sqr(hhh),-0.532));
        ab=pow(12.0*Omega*sqr(hhh),0.424)*(1+pow(45.*Omega*sqr(hhh),-0.582));
        al=pow(aa,-Omega_Bh2/sqr(hhh)/Omega)*
            pow(ab,-cube(-Omega_Bh2/sqr(hhh)/Omega));
        first=0;
    }

    k=2.0*kp*H0INV/hhh;
    /*  kp *= 1./Omega/sqr(hhh); */    /*  now scaled by Omega h^2 */
    /*  thing=pow( 6.4*kp + 3.0*kp*sqrt(3.0*kp) + sqr(1.7*kp) , 1.13);
        transfer=pow(1.+ thing, -1./1.13);
        */
    /* q2=kp/Omega/sqr(hhh)*exp(Omega_Bh2/sqr(hhh)*(1.+1./Omega)); */

    q=kp/Omega/sqr(hhh)/sqrt(al)*sqr(1.0104)/sqrt(1.0-Omega_Bh2/sqr(hhh)/Omega);
    transfer=log(1.+2.34*q)/(2.34*q)
        *pow(1.+3.89*q+sqr(16.1*q)+cube(5.46*q)+sqr(sqr(6.71*q)),-0.25);
    /*  printf("\n al=%f   q=%f  q2=%f a1=%f a2=%f\n",al,q/kp,q2/kp,aa,ab);*/
  
    return pow(k/2.,nnn)*sqr(transfer);
}

/* power spectrum */
double P(double kp)
{
    double P_cdm(double), P_trans(double);
    if (do_read_Pk>0) return P_trans(kp);
    else return P_cdm(kp);
}

double P_trans(double kp)   /* current power spectrum   */
     /*  kp = physical wave number [Mpc^-1] */
{
    double logk, T2;
    
    logk=log10(kp); /* Mpc^-1 */
    if (logk <= logk_array[1]) 
        T2=1.0;
    else if (logk >= logk_array[narrpk]) { /* extrapolate with k^4 */
        T2=logT2_array[narrpk]-4.0*(logk-logk_array[narrpk]);
        T2=pow(10.0,T2);
    }
    else {
        splint(logk_array,logT2_array,Pk_2,narrpk,logk,&T2);
        T2=pow(10.0,T2);   
    }
    return pow(kp*H0INV/hhh,nnn)*T2;
}

/* read a transfer function */
int read_trans(FILE* Pkfile, double *logk_array, double *logT2_array,
               int narrin) /* narrin not used yet */
{
    double ns=nnn;
    double k, trans;
    int i, narr;
    FILE *Poutfile;
    
    i=1;
    Poutfile=fopen("Pk.out","w"); 
    while (!feof(Pkfile)) {
        fscanf(Pkfile, "%lf %lf", &k, &trans);
        logk_array[i]=log10(k);
        logT2_array[i]=2.0*log10(trans);
        fprintf(Poutfile,"%f %f %f %f %f\n", 
                k, trans, logk_array[i], logT2_array[i], 
                pow(k,ns)*pow(10.0,logT2_array[i]));
        i++;
    }
    narr=i-2;
    fclose(Pkfile); fclose(Poutfile);
    return narr;
}

/* integrate func from a->b with NGAU-pt gaussian quadrature */
/* only calculate the weights and abcissas once, and scale (a,b)->(-1,1) */
#define NGAU 500
double qgaus_ahj(double (*func)(double), double a, double b)
{
    static int first=1;
    double sum=0.0;
    static double *x, *w;
    int i;

    if (first==1) {
        x=vector(1,NGAU);
        w=vector(1,NGAU);
        gauleg(-1.0, 1.0, x, w, NGAU);
        first=0;
    }
    for (i=1; i<=NGAU; i++) sum += w[i]*(*func)(0.5*((b+a)+(b-a)*x[i]));
    return 0.5*(b-a)*sum;
}

