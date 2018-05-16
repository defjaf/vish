double tau_integrand(double z);
double x_e(double z);
double tau(double z);
double lngfun(double z, double z_m, double x0, double xbar, double alpha); 
double gfun(double z, double z_m, double x0, double xbar, double alpha); 
double p_proj(double kappa);
double Pintegrand(double y);   /* y=eta/eta_0 */
double Pintegrand_1(double y); /* y=eta/eta_0 */
double S_interp(double);
double S(double k);
double Sapprox(double k);
double xintegrand(double x);
double xintegrand_old(double x);
double wind(double k);
double uintegrand(double u);
double yintegrand(double y);
double P_cdm(double kp);      /* current power spectrum (from bardeen et al.);  */
double wintegrand(double  z);
double www(double z);         /*  calculate w(z); by cubic interpolation  */
double zzz(double w);         /* calculate z(w); by cubic interpolation   */
double E(double z);
double g(double z);           /* calculate visibility fn by interpolation */
double D(double z);           /* calculate D(z); by cubic interpolation   */
double Dintegrand(double alpha);
double sigma8_integrand(double u);
double jj1(double x);
double P(double kp);   /* current power spectrum   */
int read_trans(FILE* Pkfile, double *logk_array, double *logT2_array,
                 int narrin); /* narrin not used yet */
double qgaus_ahj(double (*func)(double), double a, double b);
double qromb(double (*func)(double), double a, double b);
double qgaus(double (*func)(double), double a, double b);
void splint(double *, double *, double *, int, double, double *);
void spline(double *, double *, int, double, double, double *);
void gauleg(double, double, double *, double *, int);
int read_xe(FILE *xefile, double *xe_z_array, double *xe_2_array,
               int narrin);
double Pkwindow(double kp);
