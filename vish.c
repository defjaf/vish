/*    calculate Cl's for vishniac effect in flat Universe for given
      Omega (NR matter density), hhh (hubble constant in 100 km/sec/Mpc
      and nnn (spectral index).  this version assumes a gaussian
      differential optical depth with comoving-distance distribution
      centered at w_0 of width Deltaw and width tau_r

      AHJ: modified 08/00 for
      restrict calculation to (z_max,z_min)

      AHJ: modified 07/00 for

      1. arbitrary ionization fraction xe(z) read in from a file
         [do_xe_arr]

      2. allow reading of a "window function" to multiply the built-in
         transfer function.

      3. allow a window function such that W(k) is just a low or high-k cutoff
         (this seems to give problems with convergence of the integrals at
         high k)

      4. allow a window function acting on S(k) rather than P(k)
      ****************** CHECK: Pk_window(kp) vs (k)  *************************

      AHJ: modified 11/97 for
      1. visibility functions due to an ionization fraction
            |xbar                                           eta > eta_1
      x_e = |x0 (1+z)^b = x0 (a/a0)^-b = x0 (eta/eta0)^-2b  eta_m<eta<eta_1
            |0                                              eta < eta_m
      (need params x0, xbar, xm, b==-alpha)

      nb. program only needs one of x_e, tau_r, z_r (=zzz(eta0-eta_m))

      2. arbitrary P(k), read in from a file [phys k units]
      check extrapolation for isocurv models

      THIS PROGRAM USES UNITS SUCH THAT (a0 H0 = 2)!!!!!!

 */

/* ahj 05/26:
   modified to calculate C_l for patchy reionization:
   need W(kR, w) - window function defining smoothing of density -> ionization
   Xe(w) - mean ionization fraction evolution
   b(w) - "bias" for clustering of high peaks responsible for ionization

   2 terms:   <(v xe)_1 (v xe)_2> + <(v xe)_1 (v d)_2>
globals:       do_XeXe_corr           do_Xed_corr

for this model, we use
    Xe(z) = x0[1/2 - (z-z_i)/dz_i] for z_i-dz_i/2<z<z_i+dz_i
    b(z) = const. = 1
    W(k) = j1(kR)/kR
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vish.h"
// #include "nrutil.c"
// #include "qromb.c"
// #include "qgaus.c"
// #include "gauleg.c"
// #include "trapzd.c"
// #include "polint.c"
// #include "splint.c"
// #include "spline.c"

#define ZMAX 1100        /* surface of last scattering, more or less */
#define DELZT 0.25       /* tau(z), g(z) calculated with this interval */
#define DELZ  2.0        /* w(z) calculated with this interval */
#define NARR  ((int)floor(ZMAX/DELZ))   /* size of arrays */
#define NARRtau ((int)ceil(ZMAX/DELZT)) /* >(NARR*DELZ/DELZT)
                                           big so we don't fall off the end */
#define NARRPk 402
/* #define NARRPk 102*/

#define NARRXE 1400

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
/* #define USE_LNG */     /* calculate and store ln g(z) rather than g(z) */
/* #define YSCALE */      /* scale P_integrand to Pintegrand(YYY)*/

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define PI M_PI             /* (3.141592653589793) */
#define LOG10 M_LN10        /* (log(10.0)) */
#define SMALLNUM (1.e-15)   /* for ==0.0 checks */

static double Omega,hhh,nnn,Omega_Bh2,Omega_k, Omega_vac;
static double dummyk,dummyy, kappa_g;
static double *wz_array,*w_array,*w_2,*D_array,*D_2,*z_2, *logT2_array, *Pk_2;
static double *g_array, *g_2, *tau_array, *tau_2, *tauz_array, *logk_array;
static double *xe_array, *xe_2, *xe_z_array;
static double *logS_array, *logS_logk_array, *logS_2;
static double D_0, ETA0;
static double alpha_g,  zmax_g,  x0_g,  xbar_g, tau_r;
static int narrpk, narrxe, narrtau, do_read_Pk, do_window_Sk=0, do_xe_arr;
static int do_XeXe_corr, do_Xed_corr, do_linear_turnon;
static double tau_prefactor;
static double Radius, z_i, delta_zi;
static double WkR0, kphys_cut;
static double z_max, z_min;
static int do_zslice=0;

int main(int argc, char *argv[])
{
    double kph,k,lognumber,x;  /* y*/
    double z, zt, w, RS_chi;
    double delta_H,prefactor,normalization,w_r,Delta_w,Cl,ell;
    double z_r,D_r,a_r,addot_over_adot,adot_over_a,Ddot_over_D;
    double sigma82, gg;
    double z_m, x_0, xbar, alpha, g_prefactor, SSk;
    double DT_T2, dlog10ell, llCl, llCl_prev=0.0, deriv, deriv_prev=0.0, gsum=0.0;
    double logkmin, logkmax;
    int n, i, do_Gauss_vis, findzm=0;
    char Pkfilename[120], xefilename[120];
    FILE *zoutfile,*outfile,*Soutfile,*xoutfile,*Poutfile, *Pkfile, *xefile;
    /* *youtfile */
    int iahj;
    printf("argc=%d\n", argc);
    for (iahj=0; iahj<5; iahj++) printf("%d: %s\n", iahj, argv[iahj]);

    Omega=strtod(argv[1], (char **)NULL);
    hhh=strtod(argv[2], (char **)NULL);
    nnn=strtod(argv[3], (char **)NULL);
    Omega_Bh2=strtod(argv[4], (char **)NULL);
    Omega_vac=strtod(argv[5], (char **)NULL);
    Omega_k=1.0-Omega-Omega_vac;
    /*   tau_r=atof(argv[5]);
         w_r=atof(argv[6]);
         Delta_w=atof(argv[7]); */

    printf("Omega=%f, Omega_k=%f, Omega_vac=%f\n", Omega, Omega_k, Omega_vac);
    printf("hhh=%f\n", hhh);
    printf("nnn=%f\n", nnn);
    printf("Omega_Bh2=%f\n", Omega_Bh2);
    /* printf("tau_r=%f\n", tau_r);
       printf("w_r=%f\n", w_r);
       printf("Delta_w=%f\n", Delta_w); */

    logk_array=vector(1,NARRPk);
    logT2_array=vector(1,NARRPk);
    Pk_2=vector(1,NARRPk);
    logS_logk_array=vector(1,NARRPk);
    logS_array=vector(1,NARRPk);
    logS_2=vector(1,NARRPk);

    xe_z_array=vector(1,NARRXE);
    xe_array=vector(1,NARRXE);
    xe_2=vector(1,NARRXE);

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
    printf("\nRead T(k) [1 = yes, 2 = read window fn W(k), 3=cut >k, 4=cut<k]: ");
    printf("\n          [         +10 = apply to S(k) rather than P(k)      ]: ");
    scanf("%d", &do_read_Pk);
    if (do_read_Pk>10) { do_read_Pk -= 10; do_window_Sk=1; }
    if (do_read_Pk==1) {
        printf("\nP(k) filename: ");
        scanf("%s", Pkfilename);
    } else if (do_read_Pk==2) {
        printf("\n R_0 [Mpc/h] W^2(k R0) filename: ");
        scanf("%lf %s", &WkR0, Pkfilename);
    } else if (do_read_Pk==3 || do_read_Pk==4) {
        printf("\n K cutoff (h/Mpc)");
        scanf("%lf", &kphys_cut);
        kphys_cut *= hhh;
        printf("do_read_Pk=%d; kcut=%f h/Mpc = %f / Mpc\n", do_read_Pk, kphys_cut/hhh, kphys_cut);
    }
    else {
        printf("\nDummy: ");
        scanf("%s", Pkfilename);
    }
    if (do_read_Pk==1 || do_read_Pk==2) Pkfile=fopen(Pkfilename, "r");

    printf("\n1=Gaussian visibility function, 2=restrict to redshift slice: ");
    scanf("%d", &do_Gauss_vis);
    if (do_Gauss_vis==0 || do_Gauss_vis==2) {
        printf("  [only set one of x_0 or tau_r != 0]");
        printf("\nModel Ionization history, z_m, x_0, tau_r, xbar, alpha: ");
        scanf("%lf %lf %lf %lf %lf", &z_m, &x_0, &tau_r, &xbar, &alpha);
        if (do_Gauss_vis==0) {
            printf("\nDummy [3 doubles]: ");
            scanf("%lf %lf %lf", &tau_r, &w_r, &Delta_w);
        } else if (do_Gauss_vis==2) {
            printf("\nz_min, z_max, dummy (<0: ignore)");
            scanf("%lf %lf %lf", &z_min, &z_max, &Delta_w);
            printf("z_min=%f, z_max=%f\n", z_min, z_max);
            do_zslice=1; do_Gauss_vis=0;
        }
    } else {
        printf("\nDummy [5 doubles]: ");
        scanf("%lf %lf %lf %lf %lf", &z_m, &x_0, &tau_r, &xbar, &alpha);
        printf("Gaussian ionization history, tau_r, w_r, Delta_w: ");
        scanf("%lf %lf %lf", &tau_r, &w_r, &Delta_w);
    }
    printf("\nz_m=%f; x_0=%f; xbar=%f; alpha=%f\n", z_m, x_0, xbar, alpha);
    printf("tau_r=%f, w_r=%f, Delta_w=%f\n", tau_r, w_r, Delta_w);
    zmax_g=z_m; alpha_g=alpha; x0_g=x_0; xbar_g=xbar;

    /* ahj 6/98: read in parameters for patchy reionization */
    printf("\ndo <Xe Xe>, <Xe d> [>0 = yes]? ");
    scanf("%d %d", &do_XeXe_corr, &do_Xed_corr);
    if (do_XeXe_corr || do_Xed_corr) printf("\nRadius, z_i, delta z: ");
    else printf("\nDummy [3 doubles]: ");
    scanf("%lf %lf %lf", &Radius, &z_i, &delta_zi);

    /* ahj 07/00: read in parameters for ionization fraction file */
    printf("\n Read ionization fraction from file? [>0 = yes]? ");
    scanf("%d", &do_xe_arr);
    if (do_xe_arr) printf("\nFilename: "); else printf("\nDummy [string]: ");
    scanf("%s", xefilename);
    if (do_xe_arr) xefile=fopen(xefilename, "r");

    if (do_XeXe_corr || do_Xed_corr) {
        printf("R=%f, z_i=%f, delta zi=%f\n", Radius, z_i, delta_zi);
        z_m=zmax_g=z_i+delta_zi/2.0;
    }
    if (do_XeXe_corr < 0 && do_Xed_corr < 0) {
        do_linear_turnon=1;
        do_XeXe_corr=do_Xed_corr=0;
    }

#if defined(PRINT_PARAMS)
    printf("Parameters read. \n");
    printf("do_read_Pk=%d\n", do_read_Pk);
    printf("do_linear_turnon=%d\n", do_linear_turnon);
    printf("do_xe_arr=%d\n", do_xe_arr);
    printf("do_zslice=%d\n", do_zslice);
    printf("do_window_SK=%d\n", do_window_Sk);
    printf("do_XeXe_corr=%d\n", do_XeXe_corr);
    printf("do_Xed_corr=%d\n", do_Xed_corr);
    printf("do_Gauss_vis=%d\n", do_Gauss_vis);
#endif

    tau_prefactor=(NHAT*Omega_Bh2*SIGMA_T*H0INV*Mpc_to_cm/hhh);

    /* read in P(k) & set up lookup table to evaluate */
    if (do_read_Pk==1 || do_read_Pk==2) {
        narrpk=read_trans(Pkfile, logk_array, logT2_array, NARRPk);
        printf("narrpk=%d\n", narrpk);
        printf("min logk P(k)=%f %f\n", logk_array[1], logT2_array[1]);
        printf("max logk P(k)=%f %f\n",
               logk_array[narrpk], logT2_array[narrpk]);

        if (do_read_Pk==1) /* known slopes */
            spline(logk_array, logT2_array, narrpk, 0.0, -4.0, Pk_2);
        else if (do_read_Pk==2)
            spline(logk_array, logT2_array, narrpk, 0.0,  0.0, Pk_2);
    }

    /* read in xe(z) & set up lookup table to evaluate */
    if (do_xe_arr>0) {
        narrxe=read_xe(xefile, xe_z_array, xe_array, NARRXE);
        printf("narrxe=%d\n", narrxe);
        printf("min z xe(z)=%f %f\n", xe_z_array[1], xe_array[1]);
        printf("max z xe(z))=%f %f\n", xe_z_array[narrxe], xe_array[narrxe]);
        /* known slopes */
        spline(xe_z_array, xe_array, narrxe, 0.0, 0.0, xe_2);
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

    if (do_Gauss_vis<1) {
        if (fabs(z_m)<SMALLNUM) {
            zmax_g=z_m=ZMAX-1.0; findzm=1;
        } else if (fabs(x_0)<SMALLNUM) { /* need x_0 as fn. of tau_r */
            /* use tau_integrand with unnormalized ionization history */
            x0_g=1.0;
            x_0 = tau_r/(tau_prefactor*qromb(tau_integrand, 0.0, z_m));
            x0_g=x_0; /* set normalization */
        } else {
            tau_r = tau_prefactor*qromb(tau_integrand, 0.0, z_m);
        }
        printf("tau "); fflush(stdout);
        zt=0.0;
        tauz_array[1] = tau_array[1]=0.0;
        narrtau=0;
        for(n=2; n<=NARRtau; n++) {
            zt+=DELZT;
            tauz_array[n]=zt;
            if (zt <= z_m) {
                tau_array[n]=qromb(tau_integrand, 0.0, zt);
                if (findzm==1 && tau_prefactor*tau_array[n] >= tau_r) {
                    /* linear interp. to find z_r */
                    z_m=zt+DELZT*(tau_r/tau_prefactor-tau_array[n])/
                        (tau_array[n]-tau_array[n-1]);
                    zmax_g=z_m;
                    tau_array[n]=tau_r;
                    if (narrtau==0) narrtau=n-1;
                    findzm=0;
                    printf ("[z_m=%g]... ", z_m); fflush(stdout);
                }
            } else {
                tau_array[n]=tau_r;
                if (narrtau==0) narrtau=n-1;
            }
        }
        printf ("[narrtau=%d]... ", narrtau);
        spline(tauz_array,tau_array,narrtau,1.e40,1.e40,tau_2);

        printf("vis. fn... "); fflush(stdout);
        for(n=1; n<=NARRtau; n++) {
#if defined(USE_LNG)
            g_array[n]=lngfun(tauz_array[n], z_m, x_0, xbar, alpha);
#else
            g_array[n]=gfun(tauz_array[n], z_m, x_0, xbar, alpha);
#endif
        }
        spline(tauz_array,g_array,NARRtau,1.e40,1.e40,g_2);
        g_prefactor=1.0;

        printf("logS(logk)... "); fflush(stdout);
        for(i=1; i<=NARRPk; i++) {
	  if (do_XeXe_corr || do_Xed_corr) {
	    logkmin=-6.0; logkmax=6.0;
	  } else {
	    logkmin=-3.0; logkmax=4.0;
	  }
          lognumber=logkmin+(i-1)*((logkmax-logkmin)/(NARRPk-1));
          kph=pow(10.0, lognumber);
          k=2.0*H0INV*kph;
          if (i%100 == 1)  printf(" [i=%d, k=%6g]", i, k); fflush(stdout);
          logS_logk_array[i]=log10(k);
          logS_array[i]=log10((SSk=S(k))<=SMALLNUM?SMALLNUM:SSk);
          /*            printf("%d: logk, logS=%f %f\n", i,
                          logS_logk_array[i],logS_array[i]); fflush(stdout);  */
        } /* known slopes */
        spline(logS_logk_array, logS_array, NARRPk, 0.0, -3.0, logS_2);
    }

    /*  set up lookup table to evaluate D(z)  */
    printf("D(z)\n"); fflush(stdout);
    for(n=1;n<=NARR;n++) {
        z=wz_array[n];
        D_array[n]=5.*Omega/2.*E(z)*qromb(Dintegrand,1.e-5,1./(1.+z));
    }
    spline(wz_array,D_array,NARR,1.e40,1.e40,D_2);

    D_0=D(0.);
    ETA0=www(wz_array[NARR]);
    printf("D_0=%f; eta_0=%f\n", D_0, ETA0);
    if (do_Gauss_vis<1) {
        printf("ionization history: x_0=%f, tau_r=%f, z_m=%g\n", x_0, tau_r, z_m);
    }
    fflush(stdout);

    /*  print file of z, w(z), and D(z) as a diagnostic  */
    printf("writing z.out... "); fflush(stdout);
    zoutfile=fopen("z.out","w");
    for(z=0.;z<=1100.;z+=1.) {
        fprintf(zoutfile,"%f  %f  %f  %f\n",z,www(z),zzz(www(z)),D(z)/D_0);
    }
    fclose(zoutfile);

    /*   debugging routine to check the power spectrum  */
    printf("P.out... "); fflush(stdout);
    Poutfile=fopen("P.out","w");
    for(lognumber= -3.;lognumber <= 0.; lognumber += 0.1) {
        kph=pow(10.,lognumber);
        fprintf(Poutfile,"%f %f %g\n",kph, kph*hhh, P(kph*hhh));
    }
    fclose(Poutfile);

    if (do_Gauss_vis<1) {
        /*  print file of z, tau(z), g(z) as a diagnostic  */
        printf("tau.out... "); fflush(stdout);
        zoutfile=fopen("tau.out","w");
        for(z=0.;z<=1100.;z+=1.) {
#if defined(USE_LNG)
            fprintf(zoutfile,"%f %f %f %fq %g\n",z,www(z),x_e(z),tau_prefactor*tau(z),
                    gg=g_prefactor*exp(lngfun(z, zmax_g, x0_g, xbar_g, alpha_g));
#else
            fprintf(zoutfile,"%f %f %f %f  %g\n",z,www(z), x_e(z), tau_prefactor*tau(z),
                    gg=g_prefactor*gfun(z, zmax_g, x0_g, xbar_g, alpha_g));
#endif
            gsum+=gg/E(z);
        }
        gsum *= 0.5; /* 1/(a0 H0) */
        fclose(zoutfile);
        zoutfile=fopen("proj.out","w");
        printf("proj.out [k=%4.4f]... ", kappa_g);
        for(kappa_g=1e-6; kappa_g<=1.1e7; kappa_g*=100.0) {
            for(w=0.;w<=1.001;w+=0.01) {
                fprintf(zoutfile,"%e  %f  %g\n",kappa_g, w,Pintegrand(w));
            }
        }
        fclose(zoutfile);
    }

    /*  debugging routine to check the x integrand  */
    printf("x.out... "); fflush(stdout);
    xoutfile=fopen("x.out","w");
    for(x=-1.;x<=1.;x+=0.001) {
        dummyk = 6.e3;  dummyy=1.;
        fprintf(xoutfile,"%f  %g\n",x,xintegrand_old(x));
    }
    fclose(xoutfile);

    /*  debugging routine to check the y integrand; not currently operational
        youtfile=fopen("y.out","w");
        dummyk=2.0*H0INV*10.;
        for(lognumber=6.-log10(dummy)-6.;lognumbery=1.e-20; y<= 1.e6/dummyk;
        y+=1.e6/dummyk/1000.) {
        fprintf(youtfile,"%e  %e\n",y,yintegrand(y));
        }
        fclose(youtfile);
        */

    /* calculate S(k) as a diagnostic */
    printf("S.out\n"); fflush(stdout);
    Soutfile=fopen("S.out","w");
    for(lognumber=-6.; lognumber<=6.; lognumber += 0.1) {
        kph = pow(10.,lognumber);
        k=2.0*H0INV*kph;
        fprintf(Soutfile,"%e %e %e %e\n",kph,S(k),S_interp(k), Sapprox(k));
    }
    fclose(Soutfile);

    if (do_Gauss_vis>=1) {
        /* all of these assume a gaussian g(w).... */
        z_r=zzz(w_r);  D_r=D(z_r);   a_r = (2.*H0INV/hhh)/(1.+z_r);
        addot_over_adot = (hhh/H0INV)*(Omega_vac-Omega*cube(1.+z_r)/2.)
            /E(z_r);
        adot_over_a = (hhh/H0INV)*E(z_r);
        Ddot_over_D = addot_over_adot - adot_over_a
            + 5.*Omega/2.*adot_over_a * sqr((1.+z_r)/E(z_r)) / D_r;

        /*  this is the prefactor multiplying the C_l's  */
        prefactor = sqr(tau_r)/16./pow(PI,2.5)/Delta_w / sqr(w_r)
            * sqr( a_r * Ddot_over_D * sqr(D_r/D_0));
        prefactor /= 2.0; /* factor of 1/2 from
                             <qperpi qperpj> = 1/2<qq>(dij + khati khat j) */
        /* through here .... */
    } else
        printf("Check: -ln(1-\\int dw g(w))=%f ?=? tau_r =%f\n",
               -log(1.0-gsum),tau_r);


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
    printf("writing vishniac.out...\n"); fflush(stdout);
    outfile=fopen("vishniac.out","w");
    if (do_Gauss_vis>=1) {
        if (fabs(Omega_k)<SMALLNUM) RS_chi=w_r;
        else RS_chi=0.5/sqrt(Omega_k) * sinh(2.0*sqrt(Omega_k)*w_r);
    }

    DT_T2=0.0;
    dlog10ell=0.1;
    for(lognumber=1.;lognumber<=6;lognumber+=dlog10ell) {
        ell=pow(10.,lognumber);
        if (do_Gauss_vis>=1) {
            Cl=prefactor*normalization*S(ell/RS_chi);
        } else {
            Cl=g_prefactor*normalization*p_proj(ell);
        }
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

double tau_integrand(double z)
{
    return x_e(z) * sqr(1.0+z)/E(z);
}



/* ionization fraction */
double x_e(double z)
{
    double xx;
    /* globals: double alpha_g, double zmax_g, double x0_g, double xbar_g */
    /* globals: Radius, z_i, delta_zi [patchy] */
    if (do_linear_turnon || do_XeXe_corr || do_Xed_corr) { /* patchy */
        if (z<z_i-delta_zi/2.0) return x0_g;
        else if (z>z_i+delta_zi/2.0) return 0.0;
        else return x0_g*(0.5-(z-z_i)/delta_zi);
    } else {  /* standard calc. */
        if (z > zmax_g) return 0.0;
        else  {
            if (do_xe_arr) splint(xe_z_array, xe_array, xe_2, narrxe, z, &xx);
            else xx=x0_g*pow(1.0+z, alpha_g);
            return (xx>xbar_g) ? xbar_g : xx;
        }
    }
}

double tau(double z)
{
    double y;

    if (z < zmax_g) {
        splint(tauz_array,tau_array,tau_2,narrtau,z,&y);
        return y;
    } else /* kludge since it gets multiplied later...*/
        return tau_r/tau_prefactor;
}

/* the visibility function, g(eta) = tau' e^-tau
   given the model
 *    x_e = x_0 (1.0+z)^alpha = x_0 (eta/eta0)^-2alpha;   eta_m < eta < eta_1
 */
#if defined(OMEGA_EQ_1)   /* this is WRONG!!! */
double lngfun(double eta, double z_m, double x0, double xbar, double alpha)
{
    /* // this code only works for Omega==1!! // */
    static double eta_1_0=pow(x0/xbar, 1.0/2.0/alpha);  /* rel. to eta_0 */
    static double eta_m_0=pow(1.0 + z_m, -0.5);
    static tau_prefactor=(1.0/3.0)*a0*ETA0*SIGMA_T*N0*xbar;
    static tau_prime_prefactor=x0*N0*SIGMA_T*a0;
    double tau_prime, tau;
    double y=eta/ETA0;

    if (eta < eta_m) return 0.0;
    else if (y < eta_1_0) {
        tau_prime=sqr(y)*pow(y, -4.0-2.0*alpha);
        tau=(1.0/cube(eta_1_0)-1.0)+
            (x0/xbar)*3.0/(3.0+2.0*alpha)*(pow(y, -3.0-2.0*alpha)-1);
    }
    else {
        tau_prime=sqr(y)*pow(y, -4.0)*xbar/x0;
        tau=(1.0/cube(eta_1_0)-1.0)+
            (x0/xbar)*3.0/(3.0+2.0*alpha)*(pow(eta_m_0, -3.0-2.0*alpha)-1);
    }
    return tau_prime_prefactor*tau_prime * exp(-tau_prefactor*tau);
}
#else  /* but this is OK ! */
double lngfun(double z, double z_m, double x0, double xbar, double alpha)
{
    static double lng_prefactor;
    double xe=x_e(z);

    lng_prefactor=log(NHAT*Omega_Bh2*SIGMA_T*2.0*H0INV*Mpc_to_cm/hhh);
    alpha_g=alpha; zmax_g=z_m; x0_g=x0; xbar_g=xbar;

    /*    printf("lngfun: tau_prefactor, lng_prefactor, z, tau(z)=%f
          %f %f %f\n", tau_prefactor, lng_prefactor, z, tau(z)); */
    if (xe==0.0) return log(xe);
    else
        return lng_prefactor+log(x_e(z)*sqr(1.0+z))-tau(z)*tau_prefactor;
    /* NB. this returns -Inf for x_e=0 which is OK! */
}

double gfun(double z, double z_m, double x0, double xbar, double alpha)
{
    static double g_prefac_loc;
    double xe=x_e(z);

    g_prefac_loc=NHAT*Omega_Bh2*SIGMA_T*2.0*H0INV*Mpc_to_cm/hhh;
    alpha_g=alpha; zmax_g=z_m; x0_g=x0; xbar_g=xbar;

    /*    printf("lngfun: tau_prefactor, lng_prefactor, z, tau(z)=%f
          %f %f %f\n", tau_prefactor, lng_prefactor, z, tau(z)); */
    if (xe==0.0)
        return 0.0;
    else
        return g_prefac_loc*xe*sqr(1.0+z)*exp(-tau(z)*tau_prefactor);
}
#endif

#define YYY 0.5
double p_proj(double kappa)
{
    static double prefac_loc;
    static int first=1;
    static double ymin=0.0, ymax=1.0;
    if (first && (do_XeXe_corr || do_Xed_corr)) {
        ymax=1.-www(z_i-delta_zi/2.0)/ETA0;
        ymin=1.-www(z_i+delta_zi/2.0)/ETA0;
	printf("[ymin, ymax=%f %f]", ymin, ymax); fflush(stdout);
    } else if (first && do_zslice) {
        ymax = 1.-www(z_min)/ETA0;
        ymin = z_max > 0.0 ? 1.-www(z_max)/ETA0 : 0.0;
        printf("[ymin, ymax=%f %f]", ymin, ymax); fflush(stdout);
    }
    first=0;
    kappa_g=kappa;
#if defined(YSCALE)
    /*    printf("kappa, pi(yyy)=%g %g\n", kappa,prefac_loc*(16.0*PI*PI));*/
    prefac_loc=Pintegrand_1(YYY)/(16.0*PI*PI);
    return prefac_loc*qgaus_ahj(Pintegrand, ymin, ymax);
#else
    prefac_loc=1./(16.0*PI*PI);
    return prefac_loc*qgaus_ahj(Pintegrand_1, ymin, ymax);
#endif     /* with 8Pi^2 prefactor and 1/2 from <qi.qj>/2 */
}

/* scaled integrand */
double Pintegrand(double y)
{
#if defined(YSCALE)
    double prefac_loc=1.0/Pintegrand_1(YYY);
#else   /* actually, not used in this case */
    static double prefac_loc=1.;
#endif
    return Pintegrand_1(y)*prefac_loc;
}

double Pintegrand_1(double y) /* y=eta/eta_0 */
{
    double RS_chi, gg, SS, xz;
    double ww=ETA0*(1.0-y);
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
#if defined(USE_S_APPROX)
    SS=Sapprox(kappa_g/RS_chi);
#else
    SS=S_interp(kappa_g/RS_chi);
#endif
#define CLIPFAC
#if defined(CLIPFAC)
    if (do_XeXe_corr || do_Xed_corr) {
        xz=x_e(zz);
        SS*=sqr(1.0-xz);  /* is this right? or  (1-xz)/xz */
    }
#endif
#if defined(USE_LNG)
    gg=lngfun(zz,zmax_g,x0_g,xbar_g,alpha_g);
    /*     return exp(2.0*(log(aa * Ddot_over_D * sqr(DD/D_0)/ww)+gg))*SS; */
    return exp(2.0*(log(aa * Ddot_over_D * sqr(DD/D_0)/RS_chi)+gg))*SS;
#else
    gg=gfun(zz,zmax_g,x0_g,xbar_g,alpha_g);
    /*    return sqr(aa * Ddot_over_D * sqr(DD/D_0)*gg/ww)*SS; */
    return sqr(aa * Ddot_over_D * sqr(DD/D_0)*gg/RS_chi)*SS;
#endif
}

double S_interp(double k)   /* current power spectrum   */
     /*  kp = physical wave number [Mpc^-1] */
     /* nb. the extrapolations outside of the logS_logk_array are
        actually only correct for n_s=1, but it doesn't really matter! */
{
    double logk, SS;

    /**** CHECK!!! ****/
    double kp = hhh*k/2.0/H0INV;

    logk=log10(k);
    if (logk <= logS_logk_array[1])
        SS=logS_array[1];
    else if (logk >= logS_logk_array[NARRPk]) {  /* extrapolate with k^-3 */
        SS=logS_array[NARRPk]-3.0*(logk-logS_logk_array[NARRPk]);
    } else if (logk >= logS_logk_array[NARRPk]) {  /* extrapolate with k^-3 */
        SS=logS_array[NARRPk]-3.0*(logk-logS_logk_array[NARRPk]);
    }
    else {
        splint(logS_logk_array,logS_array,logS_2,NARRPk,logk,&SS);
    }

    if (do_window_Sk) {
        if (do_read_Pk==2) return  Pkwindow(kp)*pow(10.0,SS);
        else if (do_read_Pk==3)
            return (kp>=kphys_cut) ? pow(10.0,SS) : 0.0;
        else if (do_read_Pk==4)
            return (kp<=kphys_cut) ? pow(10.0,SS) : 0.0;
        else return -1;
    } else return pow(10.0,SS);
}

double S(double k)
{
    double lk=log10(k), dn=-3.0-lk, up=6.0-lk;
    dummyk=k;
    /* return k*LOG10*qromb(uintegrand,-2.-lk,6.-lk); */
    if (dn<0.0 && up>0.0) /* integrable singularity at y=1? */
#if defined(DO_GAUSS)
        return k*LOG10*(qgaus_ahj(uintegrand, dn, 0.0) +
                        qgaus_ahj(uintegrand, 0.0, up));
    else
        return  k*LOG10*(qgaus_ahj(uintegrand, dn, up));
#else
        return k*LOG10*(qromb(uintegrand, dn, 0.0) +
                        qromb(uintegrand, 0.0, up));
    else
        return  k*LOG10*(qromb(uintegrand, dn, up));
#endif
    /* to be conservative, the upper limit here is probably a bit larger
       than it needs to be.  its probably lower by an order of magnitude
       (for EdS) for k=1.e6, and lower by a factor of 100 for k=10.  */

    /* we integrate over u=log10(y) instead of y since the integrand
       is smoother in u than in y  */
}

double Sapprox(double k)
{
    static int firstS=1;
    static double kp_cut, kp;
    double aa,ab,al, I_n(double);

    if (firstS==1) {
        aa=pow(46.9*Omega*sqr(hhh),0.67)*(1.0+pow(32.1*Omega*sqr(hhh),-0.532));
        ab=pow(12.0*Omega*sqr(hhh),0.424)*(1.0+pow(45.*Omega*sqr(hhh),-0.582));
        al=pow(aa,-Omega_Bh2/sqr(hhh)/Omega)*
            pow(ab,-cube(-Omega_Bh2/sqr(hhh)/Omega));
        kp_cut=Omega*sqr(hhh)*sqrt(al)/sqr(1.0104)
            *sqrt(1.-Omega_Bh2/sqr(hhh)/Omega);
        firstS=0;
    }

    kp=k*hhh/(2.0*H0INV);
    return k * pow(k/2.,2.*nnn) * I_n(kp/kp_cut);

    /*  the normalization prefactor is NOT included, so this Sapprox()
        should be directly comparable to S_interp() or S() */
}

double I_n(double q)
{
  double h,log10q0,sigma1,sigma2,squarething,logh;

  logh = -4.716 -0.636*(nnn-1.)+0.832*sqr(nnn-1.);
  log10q0=0.477+1.36*(nnn-1.0)+0.975*sqr(nnn-1.);
  sigma1=0.6+0.33*(nnn-1.)+0.34*sqr(nnn-1.);
  sigma2=0.87+1.105*(nnn-1.)+1.19*sqr(nnn-1.);

  h=pow(10.,logh);

  if(log10(q)<log10q0){
    squarething = (log10(q)-log10q0)/sigma1;
  }
  else {
    squarething = (log10(q)-log10q0)/sigma2;
  }
  return h*pow(2.,2.*nnn+1.)/pow(q,2.0*nnn+3.)*exp( -0.5*sqr(squarething));
}

double I_n_old(double q)
{
    double a,b,c,squarething;

    a=1.e-4*(1.25 - 4.04*(nnn-1.0) +5.85*sqr(nnn-1.0) );
    b= -0.216+0.697*(nnn-1.)+0.54*sqr(nnn-1.0);
    c= 0.65 + 0.23*(nnn-1.0)+0.14*sqr(nnn-1.0);

    squarething = log10(q)-b;
    return a/pow(q,2.0*(nnn+1.0))*exp( -sqr(squarething)/c);
}
/* take into account the inverse-square-root singularity at x=y=1 */
double xintegrand(double t)
{
    return 2.0*t*xintegrand_old(1.0-t*t);
}

double xintegrand_old(double x)
{
    double fac;
    double alph2=1.+sqr(dummyy)-2.*x*dummyy;
    double kp=hhh*dummyk/2.0/H0INV;

    /* select which correlation function to evaluate based on global variables
       do_*_corr==1 */
    if (do_XeXe_corr && do_Xed_corr) {
        fac=alph2*wind(kp*sqrt(alph2))-sqr(dummyy)*wind(kp*dummyy);
        fac*=fac+1.-2.*x*dummyy;
    } else if (do_XeXe_corr) {
        fac=alph2*wind(kp*sqrt(alph2))-sqr(dummyy)*wind(kp*dummyy);
        fac*=fac;
    } else if (do_Xed_corr) {
        fac=alph2*wind(kp*sqrt(alph2))-sqr(dummyy)*wind(kp*dummyy);
        fac*=1.-2.*x*dummyy;
    } else {
        fac=1.-2.*x*dummyy;
        fac*=fac;
    }

    return P(kp*sqrt(alph2))*(1.-x*x)*fac/sqr(alph2);
}

/* window function for density -> ionization fraction  */
double wind(double kph)  /* wavenumber in Mpc */
{
    /* global: Radius */#define GAUSSIAN_WINDOW
#if defined(GAUSSIAN_WINDOW)
    return exp(-0.5 * kph*kph * Radius*Radius);
#else
    return jj1(kph*Radius)/(kph*Radius); /* top hat */
#endif
}

double uintegrand(double u)
{
    double kp,y;

    y=pow(10.,u);
    kp=hhh*dummyk/2.0/H0INV;
    dummyy=y;

#if defined(SQRT_SING)
    return y * P(kp*y) * qgaus_ahj(xintegrand, 0.0, sqrt(2.0));
#else
    return y * P(kp*y) * qgaus_ahj(xintegrand_old,-1.,1.);
#endif
}


double yintegrand(double y)
{
    double kp;

    kp=hhh*dummyk/2.0/H0INV;
    dummyy=y;

    return P(kp*y) * qgaus_ahj(xintegrand,-1.,1.);
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


double g(double z) /* calculate visibility fn by interpolation */
{
    double y;

    splint(tauz_array,g_array,g_2,NARRtau,z,&y);
    return y;
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
        aa=pow(46.9*Omega*sqr(hhh),0.67)*(1.0+pow(32.1*Omega*sqr(hhh),-0.532));
        ab=pow(12.0*Omega*sqr(hhh),0.424)*(1.0+pow(45.*Omega*sqr(hhh),-0.582));
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

    if (do_read_Pk==1) return P_trans(kp);
    else if (!do_window_Sk && do_read_Pk) {
        if (do_read_Pk==2) return Pkwindow(kp)*P_cdm(kp);
        else if (do_read_Pk==3) return (kp>=kphys_cut) ? P_cdm(kp) : 0.0;
        else if (do_read_Pk==4) return (kp<=kphys_cut) ? P_cdm(kp) : 0.0;
        else return -1;
    }
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

double Pkwindow(double kp)   /* power spectrum window fn  */
     /*  kp = physical wave number [Mpc^-1] */
     /* W tabulated at W(x); x=k * R0/(Mpc/h) */
{
    double logk, T2;
    double x=WkR0/hhh * kp;

    logk=log10(x); /* Mpc^-1 */
    if (logk<=logk_array[1]) T2=logT2_array[1];
    else if (logk>=logk_array[narrpk]) T2=logT2_array[narrpk];
    else splint(logk_array,logT2_array,Pk_2,narrpk,logk,&T2);
    return (T2>0.0) ? 1.0 : pow(10.0,T2);
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
        if (k==0.0) k=1.0e-30;
        if (trans==0.0) trans=1.0e-30;
        logk_array[i]=log10(k);
        if (do_read_Pk==1) {
            logT2_array[i]=2.0*log10(trans);
            fprintf(Poutfile,"%f %f %f %f %f\n",
                    k, trans, logk_array[i], logT2_array[i],
                    pow(k,ns)*pow(10.0,logT2_array[i]));
        } else if (do_read_Pk==2) logT2_array[i]=log10(trans);
        i++;
    }
    narr=i-2;
    fclose(Pkfile); fclose(Poutfile);
    return narr;
}

/* integrate func from a->b with NGAU-pt gaussian quadrature */
/* only calculate the weights and abcissas once, and scale (a,b)->(-1,1) */
#define NGAU 300
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

int read_xe(FILE *xefile, double *xe_z_array, double *xe_2_array,
               int narrin) /* narrin not used yet */
{
    double z, xe;
    int i, narr;

    i=1;
    while (!feof(xefile)) {
        fscanf(xefile, "%*lf %lf %lf", &z, &xe);
        xe_array[i]=xe;
        xe_z_array[i]=z;
        i++;
    }
    narr=i-2;
    fclose(xefile);
    return narr;
}
