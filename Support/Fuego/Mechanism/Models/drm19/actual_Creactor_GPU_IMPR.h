/*save inv molecular weights into array */
__device__ void imolecularWeight_d(double * iwt);
/* Returns R, Rc, Patm */
__device__ void ckrp_d( double * ru, double * ruc, double * pa);
/*Compute P = rhoRT/W(x) */
__device__ void ckpx_d(double * rho, double * T, double * x, double * P);
/*Compute P = rhoRT/W(y) */
__device__ void ckpy_d(double * rho, double * T, double * y_wk, double * P);
/*Compute rho = P*W(y)/RT */
__device__ void ckrhoy_d(double * P, double * T, double * y_wk, double * rho);
/*convert y[species] (mass fracs) to c[species] (molar conc) */
__device__ void ckytcr_d(double * rho, double * T, double * y_wk, double * c);
/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
__device__ void ckcvms_d(double * T, double * cvms);
/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
__device__ void ckcpms_d(double * T, double * cpms);
/*Returns internal energy in mass units (Eq 30.) */
__device__ void ckums_d(double * T, double * ums);
/*Returns enthalpy in mass units (Eq 27.) */
__device__ void ckhms_d(double * T, double * hms);
/*Returns the mean specific heat at CP (Eq. 34) */
__device__ void ckcpbs_d(double * T, double * y_wk, double * cpbs);
/*Returns the mean specific heat at CV (Eq. 36) */
__device__ void ckcvbs_d(double * T, double * y_wk, double * cvbs);
/*Returns mean enthalpy of mixture in mass units */
__device__ void ckhbms_d(double * T, double * y_wk, double * hbms);
/*get mean internal energy in mass units */
__device__ void ckubms_d(double * T, double * y_wk, double * ubms);
/*compute the production rate for each species */
__device__ void ckwc_d(double * T, double * C, double * wdot);
/*compute the production rate for each species */
__device__ void productionRate_d(double * wdot, double * sc, double T);
__device__ void comp_qfqr_d(double *  qf, double * qr, double * sc, double * tc, double invT);
/*compute the production rate for each species */
__device__ void dwdot_d(double * J, double * sc, double * Tp, int * consP);
/*compute the reaction Jacobian */
__device__ void ajacobian_d(double * J, double * sc, double T, int consP);
/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void gibbs_d(double * species, double *  tc);
/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cv_R_d(double * species, double *  tc);
/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cp_R_d(double * species, double *  tc);
/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void dcvpRdT_d(double * species, double *  tc);
/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesInternalEnergy_d(double * species, double *  tc);
/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesEnthalpy_d(double * species, double *  tc);
/*save molecular weights into array */
__device__ void molecularWeight_d(double * wt);
/* get temperature given internal energy in mass units and mass fracs */
__device__ void get_t_given_ey_d_(double * e, double * y_wk, double * t, int * ierr);
/* get temperature given enthalpy in mass units and mass fracs */
__device__ void get_t_given_hy_d_(double * h, double * y_wk, double * t, int * ierr);
