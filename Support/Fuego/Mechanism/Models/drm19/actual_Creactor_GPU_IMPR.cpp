#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


/*save inv molecular weights into array */
__device__ void imolecularWeight_d(double * iwt)
{
    iwt[0] = 0.496047; /*H2 */
    iwt[1] = 0.992093; /*H */
    iwt[2] = 0.062502; /*O */
    iwt[3] = 0.031251; /*O2 */
    iwt[4] = 0.058798; /*OH */
    iwt[5] = 0.055508; /*H2O */
    iwt[6] = 0.030297; /*HO2 */
    iwt[7] = 0.071291; /*CH2 */
    iwt[8] = 0.071291; /*CH2(S) */
    iwt[9] = 0.066511; /*CH3 */
    iwt[10] = 0.062332; /*CH4 */
    iwt[11] = 0.035701; /*CO */
    iwt[12] = 0.022722; /*CO2 */
    iwt[13] = 0.034461; /*HCO */
    iwt[14] = 0.033304; /*CH2O */
    iwt[15] = 0.032222; /*CH3O */
    iwt[16] = 0.035645; /*C2H4 */
    iwt[17] = 0.034409; /*C2H5 */
    iwt[18] = 0.033256; /*C2H6 */
    iwt[19] = 0.035697; /*N2 */
    iwt[20] = 0.025033; /*AR */

    return;
}


/* Returns R, Rc, Patm */
__device__ void ckrp_d( double * ru, double * ruc, double * pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
__device__ void ckpx_d(double * rho, double * T, double * x, double * P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
__device__ void ckpy_d(double * rho, double * T, double * y_wk, double * P)
{
    double imw[21];/* inv molecular weight array */
    imolecularWeight_d(imw);
    double YOW = 0;/* for computing mean MW */
    YOW += y_wk[0]*imw[0]; /*H2 */
    YOW += y_wk[1]*imw[1]; /*H */
    YOW += y_wk[2]*imw[2]; /*O */
    YOW += y_wk[3]*imw[3]; /*O2 */
    YOW += y_wk[4]*imw[4]; /*OH */
    YOW += y_wk[5]*imw[5]; /*H2O */
    YOW += y_wk[6]*imw[6]; /*HO2 */
    YOW += y_wk[7]*imw[7]; /*CH2 */
    YOW += y_wk[8]*imw[8]; /*CH2(S) */
    YOW += y_wk[9]*imw[9]; /*CH3 */
    YOW += y_wk[10]*imw[10]; /*CH4 */
    YOW += y_wk[11]*imw[11]; /*CO */
    YOW += y_wk[12]*imw[12]; /*CO2 */
    YOW += y_wk[13]*imw[13]; /*HCO */
    YOW += y_wk[14]*imw[14]; /*CH2O */
    YOW += y_wk[15]*imw[15]; /*CH3O */
    YOW += y_wk[16]*imw[16]; /*C2H4 */
    YOW += y_wk[17]*imw[17]; /*C2H5 */
    YOW += y_wk[18]*imw[18]; /*C2H6 */
    YOW += y_wk[19]*imw[19]; /*N2 */
    YOW += y_wk[20]*imw[20]; /*AR */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute rho = P*W(y)/RT */
__device__ void ckrhoy_d(double * P, double * T, double * y_wk, double * rho)
{
    double YOW = 0;
    double imw[21];/* inv molecular weight array */
    imolecularWeight_d(imw);
    double tmp[21];

    for (int i = 0; i < 21; i++)
    {
        tmp[i] = y_wk[i]*imw[i];
    }
    for (int i = 0; i < 21; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
__device__ void ckytcr_d(double * rho, double * T, double * y_wk, double * c)
{
    double imw[21];/* inv molecular weight array */
    imolecularWeight_d(imw);
    for (int i = 0; i < 21; i++)
    {
        c[i] = (*rho)  * y_wk[i] * imw[i];
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
__device__ void ckcvms_d(double * T, double * cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R_d(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124383662212169e+07; /*H2 */
    cvms[1] *= 8.248767324424338e+07; /*H */
    cvms[2] *= 5.196763628636074e+06; /*O */
    cvms[3] *= 2.598381814318037e+06; /*O2 */
    cvms[4] *= 4.888768810227566e+06; /*OH */
    cvms[5] *= 4.615239012974499e+06; /*H2O */
    cvms[6] *= 2.519031701678171e+06; /*HO2 */
    cvms[7] *= 5.927466067445207e+06; /*CH2 */
    cvms[8] *= 5.927466067445207e+06; /*CH2(S) */
    cvms[9] *= 5.530081023953346e+06; /*CH3 */
    cvms[10] *= 5.182630712527496e+06; /*CH4 */
    cvms[11] *= 2.968349425484326e+06; /*CO */
    cvms[12] *= 1.889234139098090e+06; /*CO2 */
    cvms[13] *= 2.865242610581105e+06; /*HCO */
    cvms[14] *= 2.769058254894261e+06; /*CH2O */
    cvms[15] *= 2.679121853578248e+06; /*CH3O */
    cvms[16] *= 2.963733033722604e+06; /*C2H4 */
    cvms[17] *= 2.860941121011349e+06; /*C2H5 */
    cvms[18] *= 2.765040511976673e+06; /*C2H6 */
    cvms[19] *= 2.968047434442088e+06; /*N2 */
    cvms[20] *= 2.081333233203164e+06; /*AR */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
__device__ void ckcpms_d(double * T, double * cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R_d(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124383662212169e+07; /*H2 */
    cpms[1] *= 8.248767324424338e+07; /*H */
    cpms[2] *= 5.196763628636074e+06; /*O */
    cpms[3] *= 2.598381814318037e+06; /*O2 */
    cpms[4] *= 4.888768810227566e+06; /*OH */
    cpms[5] *= 4.615239012974499e+06; /*H2O */
    cpms[6] *= 2.519031701678171e+06; /*HO2 */
    cpms[7] *= 5.927466067445207e+06; /*CH2 */
    cpms[8] *= 5.927466067445207e+06; /*CH2(S) */
    cpms[9] *= 5.530081023953346e+06; /*CH3 */
    cpms[10] *= 5.182630712527496e+06; /*CH4 */
    cpms[11] *= 2.968349425484326e+06; /*CO */
    cpms[12] *= 1.889234139098090e+06; /*CO2 */
    cpms[13] *= 2.865242610581105e+06; /*HCO */
    cpms[14] *= 2.769058254894261e+06; /*CH2O */
    cpms[15] *= 2.679121853578248e+06; /*CH3O */
    cpms[16] *= 2.963733033722604e+06; /*C2H4 */
    cpms[17] *= 2.860941121011349e+06; /*C2H5 */
    cpms[18] *= 2.765040511976673e+06; /*C2H6 */
    cpms[19] *= 2.968047434442088e+06; /*N2 */
    cpms[20] *= 2.081333233203164e+06; /*AR */
}


/*Returns internal energy in mass units (Eq 30.) */
__device__ void ckums_d(double * T, double * ums)
{
    double imw[21];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    imolecularWeight_d(imw);
    speciesInternalEnergy_d(ums, tc);
    for (int i = 0; i < 21; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
__device__ void ckhms_d(double * T, double * hms)
{
    double imw[21];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    imolecularWeight_d(imw);
    speciesEnthalpy_d(hms, tc);
    for (int i = 0; i < 21; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns the mean specific heat at CP (Eq. 34) */
__device__ void ckcpbs_d(double * T, double * y_wk, double * cpbs)
{
    double result = 0; 
    double imw[21];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[21], tresult[21]; /* temporary storage */
    imolecularWeight_d(imw);
    cp_R_d(cpor, tc);
    for (int i = 0; i < 21; i++)
    {
        tresult[i] = cpor[i]*y_wk[i]*imw[i];

    }
    for (int i = 0; i < 21; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
__device__ void ckcvbs_d(double * T, double * y_wk, double * cvbs)
{
    double result = 0; 
    double imw[21];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[21]; /* temporary storage */
    imolecularWeight_d(imw);
    cv_R_d(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y_wk[0]*imw[0]; /*H2 */
    result += cvor[1]*y_wk[1]*imw[1]; /*H */
    result += cvor[2]*y_wk[2]*imw[2]; /*O */
    result += cvor[3]*y_wk[3]*imw[3]; /*O2 */
    result += cvor[4]*y_wk[4]*imw[4]; /*OH */
    result += cvor[5]*y_wk[5]*imw[5]; /*H2O */
    result += cvor[6]*y_wk[6]*imw[6]; /*HO2 */
    result += cvor[7]*y_wk[7]*imw[7]; /*CH2 */
    result += cvor[8]*y_wk[8]*imw[8]; /*CH2(S) */
    result += cvor[9]*y_wk[9]*imw[9]; /*CH3 */
    result += cvor[10]*y_wk[10]*imw[10]; /*CH4 */
    result += cvor[11]*y_wk[11]*imw[11]; /*CO */
    result += cvor[12]*y_wk[12]*imw[12]; /*CO2 */
    result += cvor[13]*y_wk[13]*imw[13]; /*HCO */
    result += cvor[14]*y_wk[14]*imw[14]; /*CH2O */
    result += cvor[15]*y_wk[15]*imw[15]; /*CH3O */
    result += cvor[16]*y_wk[16]*imw[16]; /*C2H4 */
    result += cvor[17]*y_wk[17]*imw[17]; /*C2H5 */
    result += cvor[18]*y_wk[18]*imw[18]; /*C2H6 */
    result += cvor[19]*y_wk[19]*imw[19]; /*N2 */
    result += cvor[20]*y_wk[20]*imw[20]; /*AR */

    *cvbs = result * 8.31451e+07;
}


/*Returns mean enthalpy of mixture in mass units */
__device__ void ckhbms_d(double * T, double * y_wk, double * hbms)
{
    double result = 0;
    double imw[21];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[21], tmp[21]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy_d(hml, tc);
    imolecularWeight_d(imw);
    int id;
    for (id = 0; id < 21; ++id) {
        tmp[id] = y_wk[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 21; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in mass units */
__device__ void ckubms_d(double * T, double * y_wk, double * ubms)
{
    double result = 0;
    double imw[21];/* inv molecular weight array */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[21]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy_d(ums, tc);
    imolecularWeight_d(imw);
    /*perform dot product + scaling by wt */
    result += y_wk[0]*ums[0]*imw[0]; /*H2 */
    result += y_wk[1]*ums[1]*imw[1]; /*H */
    result += y_wk[2]*ums[2]*imw[2]; /*O */
    result += y_wk[3]*ums[3]*imw[3]; /*O2 */
    result += y_wk[4]*ums[4]*imw[4]; /*OH */
    result += y_wk[5]*ums[5]*imw[5]; /*H2O */
    result += y_wk[6]*ums[6]*imw[6]; /*HO2 */
    result += y_wk[7]*ums[7]*imw[7]; /*CH2 */
    result += y_wk[8]*ums[8]*imw[8]; /*CH2(S) */
    result += y_wk[9]*ums[9]*imw[9]; /*CH3 */
    result += y_wk[10]*ums[10]*imw[10]; /*CH4 */
    result += y_wk[11]*ums[11]*imw[11]; /*CO */
    result += y_wk[12]*ums[12]*imw[12]; /*CO2 */
    result += y_wk[13]*ums[13]*imw[13]; /*HCO */
    result += y_wk[14]*ums[14]*imw[14]; /*CH2O */
    result += y_wk[15]*ums[15]*imw[15]; /*CH3O */
    result += y_wk[16]*ums[16]*imw[16]; /*C2H4 */
    result += y_wk[17]*ums[17]*imw[17]; /*C2H5 */
    result += y_wk[18]*ums[18]*imw[18]; /*C2H6 */
    result += y_wk[19]*ums[19]*imw[19]; /*N2 */
    result += y_wk[20]*ums[20]*imw[20]; /*AR */

    *ubms = result * RT;
}


/*compute the production rate for each species */
__device__ void ckwc_d(double * T, double * C, double * wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate_d(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}

/*compute the production rate for each species */
__device__ void productionRate_d(double * wdot, double * sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    double qdot, q_f[84], q_r[84];
    comp_qfqr_d(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 21; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[16] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[9] -= 2 * qdot;
    wdot[18] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[8] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[22]-q_r[22];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[13] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[27]-q_r[27];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[1] -= qdot;
    wdot[3] -= 2 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[30]-q_r[30];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[19] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[31]-q_r[31];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[0] -= qdot;
    wdot[0] += 2 * qdot;
    wdot[1] -= 2 * qdot;

    qdot = q_f[34]-q_r[34];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[0] += qdot;
    wdot[1] -= 2 * qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[1] -= qdot;
    wdot[4] += 2 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[40]-q_r[40];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[41]-q_r[41];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[44]-q_r[44];
    wdot[2] += qdot;
    wdot[4] -= 2 * qdot;
    wdot[5] += qdot;

    qdot = q_f[45]-q_r[45];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[48]-q_r[48];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[50]-q_r[50];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[52]-q_r[52];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[54]-q_r[54];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[55]-q_r[55];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[56]-q_r[56];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[59]-q_r[59];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[60]-q_r[60];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[61]-q_r[61];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[7] -= qdot;
    wdot[9] += 2 * qdot;
    wdot[10] -= qdot;

    qdot = q_f[63]-q_r[63];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[19] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[64]-q_r[64];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[67]-q_r[67];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[68]-q_r[68];
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[69]-q_r[69];
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[70]-q_r[70];
    wdot[8] -= qdot;
    wdot[9] += 2 * qdot;
    wdot[10] -= qdot;

    qdot = q_f[71]-q_r[71];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[73]-q_r[73];
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[74]-q_r[74];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[75]-q_r[75];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[76]-q_r[76];
    wdot[1] += qdot;
    wdot[9] -= 2 * qdot;
    wdot[17] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[78]-q_r[78];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[79]-q_r[79];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[80]-q_r[80];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[81]-q_r[81];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    return;
}

__device__ void comp_qfqr_d(double *  qf, double * qr, double * sc, double * tc, double invT)
{

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    qf[0] = sc[1]*sc[7];
    qr[0] = sc[9];

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    qf[1] = sc[1]*sc[9];
    qr[1] = sc[10];

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    qf[2] = sc[1]*sc[13];
    qr[2] = sc[14];

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    qf[3] = sc[1]*sc[14];
    qr[3] = sc[15];

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    qf[4] = sc[1]*sc[16];
    qr[4] = sc[17];

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    qf[5] = sc[1]*sc[17];
    qr[5] = sc[18];

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    qf[6] = sc[0]*sc[11];
    qr[6] = sc[14];

    /*reaction 8: 2 CH3 (+M) <=> C2H6 (+M) */
    qf[7] = sc[9]*sc[9];
    qr[7] = sc[18];

    /*reaction 9: O + H + M <=> OH + M */
    qf[8] = sc[1]*sc[2];
    qr[8] = sc[4];

    /*reaction 10: O + CO + M <=> CO2 + M */
    qf[9] = sc[2]*sc[11];
    qr[9] = sc[12];

    /*reaction 11: H + O2 + M <=> HO2 + M */
    qf[10] = sc[1]*sc[3];
    qr[10] = sc[6];

    /*reaction 12: 2 H + M <=> H2 + M */
    qf[11] = sc[1]*sc[1];
    qr[11] = sc[0];

    /*reaction 13: H + OH + M <=> H2O + M */
    qf[12] = sc[1]*sc[4];
    qr[12] = sc[5];

    /*reaction 14: HCO + M <=> H + CO + M */
    qf[13] = sc[13];
    qr[13] = sc[1]*sc[11];

    /*reaction 15: O + H2 <=> H + OH */
    qf[14] = sc[0]*sc[2];
    qr[14] = sc[1]*sc[4];

    /*reaction 16: O + HO2 <=> OH + O2 */
    qf[15] = sc[2]*sc[6];
    qr[15] = sc[3]*sc[4];

    /*reaction 17: O + CH2 <=> H + HCO */
    qf[16] = sc[2]*sc[7];
    qr[16] = sc[1]*sc[13];

    /*reaction 18: O + CH2(S) <=> H + HCO */
    qf[17] = sc[2]*sc[8];
    qr[17] = sc[1]*sc[13];

    /*reaction 19: O + CH3 <=> H + CH2O */
    qf[18] = sc[2]*sc[9];
    qr[18] = sc[1]*sc[14];

    /*reaction 20: O + CH4 <=> OH + CH3 */
    qf[19] = sc[2]*sc[10];
    qr[19] = sc[4]*sc[9];

    /*reaction 21: O + HCO <=> OH + CO */
    qf[20] = sc[2]*sc[13];
    qr[20] = sc[4]*sc[11];

    /*reaction 22: O + HCO <=> H + CO2 */
    qf[21] = sc[2]*sc[13];
    qr[21] = sc[1]*sc[12];

    /*reaction 23: O + CH2O <=> OH + HCO */
    qf[22] = sc[2]*sc[14];
    qr[22] = sc[4]*sc[13];

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    qf[23] = sc[2]*sc[16];
    qr[23] = sc[9]*sc[13];

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    qf[24] = sc[2]*sc[17];
    qr[24] = sc[9]*sc[14];

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    qf[25] = sc[2]*sc[18];
    qr[25] = sc[4]*sc[17];

    /*reaction 27: O2 + CO <=> O + CO2 */
    qf[26] = sc[3]*sc[11];
    qr[26] = sc[2]*sc[12];

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    qf[27] = sc[3]*sc[14];
    qr[27] = sc[6]*sc[13];

    /*reaction 29: H + 2 O2 <=> HO2 + O2 */
    qf[28] = sc[1]*sc[3]*sc[3];
    qr[28] = sc[3]*sc[6];

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    qf[29] = sc[1]*sc[3]*sc[5];
    qr[29] = sc[5]*sc[6];

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    qf[30] = sc[1]*sc[3]*sc[19];
    qr[30] = sc[6]*sc[19];

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    qf[31] = sc[1]*sc[3]*sc[20];
    qr[31] = sc[6]*sc[20];

    /*reaction 33: H + O2 <=> O + OH */
    qf[32] = sc[1]*sc[3];
    qr[32] = sc[2]*sc[4];

    /*reaction 34: 2 H + H2 <=> 2 H2 */
    qf[33] = sc[0]*sc[1]*sc[1];
    qr[33] = sc[0]*sc[0];

    /*reaction 35: 2 H + H2O <=> H2 + H2O */
    qf[34] = sc[1]*sc[1]*sc[5];
    qr[34] = sc[0]*sc[5];

    /*reaction 36: 2 H + CO2 <=> H2 + CO2 */
    qf[35] = sc[1]*sc[1]*sc[12];
    qr[35] = sc[0]*sc[12];

    /*reaction 37: H + HO2 <=> O2 + H2 */
    qf[36] = sc[1]*sc[6];
    qr[36] = sc[0]*sc[3];

    /*reaction 38: H + HO2 <=> 2 OH */
    qf[37] = sc[1]*sc[6];
    qr[37] = sc[4]*sc[4];

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    qf[38] = sc[1]*sc[10];
    qr[38] = sc[0]*sc[9];

    /*reaction 40: H + HCO <=> H2 + CO */
    qf[39] = sc[1]*sc[13];
    qr[39] = sc[0]*sc[11];

    /*reaction 41: H + CH2O <=> HCO + H2 */
    qf[40] = sc[1]*sc[14];
    qr[40] = sc[0]*sc[13];

    /*reaction 42: H + CH3O <=> OH + CH3 */
    qf[41] = sc[1]*sc[15];
    qr[41] = sc[4]*sc[9];

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    qf[42] = sc[1]*sc[18];
    qr[42] = sc[0]*sc[17];

    /*reaction 44: OH + H2 <=> H + H2O */
    qf[43] = sc[0]*sc[4];
    qr[43] = sc[1]*sc[5];

    /*reaction 45: 2 OH <=> O + H2O */
    qf[44] = sc[4]*sc[4];
    qr[44] = sc[2]*sc[5];

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    qf[45] = sc[4]*sc[6];
    qr[45] = sc[3]*sc[5];

    /*reaction 47: OH + CH2 <=> H + CH2O */
    qf[46] = sc[4]*sc[7];
    qr[46] = sc[1]*sc[14];

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    qf[47] = sc[4]*sc[8];
    qr[47] = sc[1]*sc[14];

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    qf[48] = sc[4]*sc[9];
    qr[48] = sc[5]*sc[7];

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    qf[49] = sc[4]*sc[9];
    qr[49] = sc[5]*sc[8];

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    qf[50] = sc[4]*sc[10];
    qr[50] = sc[5]*sc[9];

    /*reaction 52: OH + CO <=> H + CO2 */
    qf[51] = sc[4]*sc[11];
    qr[51] = sc[1]*sc[12];

    /*reaction 53: OH + HCO <=> H2O + CO */
    qf[52] = sc[4]*sc[13];
    qr[52] = sc[5]*sc[11];

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    qf[53] = sc[4]*sc[14];
    qr[53] = sc[5]*sc[13];

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    qf[54] = sc[4]*sc[18];
    qr[54] = sc[5]*sc[17];

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    qf[55] = sc[6]*sc[7];
    qr[55] = sc[4]*sc[14];

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    qf[56] = sc[6]*sc[9];
    qr[56] = sc[3]*sc[10];

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    qf[57] = sc[6]*sc[9];
    qr[57] = sc[4]*sc[15];

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    qf[58] = sc[6]*sc[11];
    qr[58] = sc[4]*sc[12];

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    qf[59] = sc[3]*sc[7];
    qr[59] = sc[4]*sc[13];

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    qf[60] = sc[0]*sc[7];
    qr[60] = sc[1]*sc[9];

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    qf[61] = sc[7]*sc[9];
    qr[61] = sc[1]*sc[16];

    /*reaction 63: CH2 + CH4 <=> 2 CH3 */
    qf[62] = sc[7]*sc[10];
    qr[62] = sc[9]*sc[9];

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    qf[63] = sc[8]*sc[19];
    qr[63] = sc[7]*sc[19];

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    qf[64] = sc[8]*sc[20];
    qr[64] = sc[7]*sc[20];

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    qf[65] = sc[3]*sc[8];
    qr[65] = sc[1]*sc[4]*sc[11];

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    qf[66] = sc[3]*sc[8];
    qr[66] = sc[5]*sc[11];

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    qf[67] = sc[0]*sc[8];
    qr[67] = sc[1]*sc[9];

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    qf[68] = sc[5]*sc[8];
    qr[68] = sc[5]*sc[7];

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    qf[69] = sc[8]*sc[9];
    qr[69] = sc[1]*sc[16];

    /*reaction 71: CH2(S) + CH4 <=> 2 CH3 */
    qf[70] = sc[8]*sc[10];
    qr[70] = sc[9]*sc[9];

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    qf[71] = sc[8]*sc[11];
    qr[71] = sc[7]*sc[11];

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    qf[72] = sc[8]*sc[12];
    qr[72] = sc[7]*sc[12];

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    qf[73] = sc[8]*sc[12];
    qr[73] = sc[11]*sc[14];

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    qf[74] = sc[3]*sc[9];
    qr[74] = sc[2]*sc[15];

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    qf[75] = sc[3]*sc[9];
    qr[75] = sc[4]*sc[14];

    /*reaction 77: 2 CH3 <=> H + C2H5 */
    qf[76] = sc[9]*sc[9];
    qr[76] = sc[1]*sc[17];

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    qf[77] = sc[9]*sc[13];
    qr[77] = sc[10]*sc[11];

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    qf[78] = sc[9]*sc[14];
    qr[78] = sc[10]*sc[13];

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    qf[79] = sc[9]*sc[18];
    qr[79] = sc[10]*sc[17];

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    qf[80] = sc[5]*sc[13];
    qr[80] = sc[1]*sc[5]*sc[11];

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    qf[81] = sc[3]*sc[13];
    qr[81] = sc[6]*sc[11];

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    qf[82] = sc[3]*sc[15];
    qr[82] = sc[6]*sc[14];

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    qf[83] = sc[3]*sc[17];
    qr[83] = sc[6]*sc[16];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 21; ++i) {
        mixture += sc[i];
    }

    /* Evaluate the kfs */
    double k_f, Corr;
    double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  H + CH2 (+M) <=> CH3 (+M)
    k_f = 1.0000000000000002e-06 * 25000000000000000 
               * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 3.2000000000000002e+27 
               * exp(-3.1400000000000001  * tc[0] - 0.50321666580471969  * 1230 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(78) > 1.e-100 ? (1.-0.68000000000000005)*exp(-tc[1] / 78) : 0.) 
        + (fabs(1995) > 1.e-100 ? 0.68000000000000005 * exp(-tc[1]/1995) : 0.) 
        + (4 == 4 ? exp(-5590 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    // (1):  H + CH3 (+M) <=> CH4 (+M)
    k_f = 1.0000000000000002e-06 * 12700000000000000 
               * exp(-0.63 * tc[0] - 0.50321666580471969 * 383 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 2.4769999999999999e+33 
               * exp(-4.7599999999999998  * tc[0] - 0.50321666580471969  * 2440 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(74) > 1.e-100 ? (1.-0.78300000000000003)*exp(-tc[1] / 74) : 0.) 
        + (fabs(2941) > 1.e-100 ? 0.78300000000000003 * exp(-tc[1]/2941) : 0.) 
        + (4 == 4 ? exp(-6964 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[1] *= Corr * k_f;
    // (2):  H + HCO (+M) <=> CH2O (+M)
    k_f = 1.0000000000000002e-06 * 1090000000000 
               * exp(0.47999999999999998 * tc[0] - 0.50321666580471969 * -260 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.35e+24 
               * exp(-2.5699999999999998  * tc[0] - 0.50321666580471969  * 1425 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(271) > 1.e-100 ? (1.-0.78239999999999998)*exp(-tc[1] / 271) : 0.) 
        + (fabs(2755) > 1.e-100 ? 0.78239999999999998 * exp(-tc[1]/2755) : 0.) 
        + (4 == 4 ? exp(-6570 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[2] *= Corr * k_f;
    // (3):  H + CH2O (+M) <=> CH3O (+M)
    k_f = 1.0000000000000002e-06 * 540000000000 
               * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * 2600 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    redP = Corr / k_f * 1e-12 * 2.2e+30 
               * exp(-4.7999999999999998  * tc[0] - 0.50321666580471969  * 5560 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(94) > 1.e-100 ? (1.-0.75800000000000001)*exp(-tc[1] / 94) : 0.) 
        + (fabs(1555) > 1.e-100 ? 0.75800000000000001 * exp(-tc[1]/1555) : 0.) 
        + (4 == 4 ? exp(-4200 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[3] *= Corr * k_f;
    // (4):  H + C2H4 (+M) <=> C2H5 (+M)
    k_f = 1.0000000000000002e-06 * 1080000000000 
               * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * 1820 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.1999999999999999e+42 
               * exp(-7.6200000000000001  * tc[0] - 0.50321666580471969  * 6970 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(210) > 1.e-100 ? (1.-0.97529999999999994)*exp(-tc[1] / 210) : 0.) 
        + (fabs(984) > 1.e-100 ? 0.97529999999999994 * exp(-tc[1]/984) : 0.) 
        + (4 == 4 ? exp(-4374 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[4] *= Corr * k_f;
    // (5):  H + C2H5 (+M) <=> C2H6 (+M)
    k_f = 1.0000000000000002e-06 * 5.21e+17 
               * exp(-0.98999999999999999 * tc[0] - 0.50321666580471969 * 1580 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.9900000000000001e+41 
               * exp(-7.0800000000000001  * tc[0] - 0.50321666580471969  * 6685 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(125) > 1.e-100 ? (1.-0.84219999999999995)*exp(-tc[1] / 125) : 0.) 
        + (fabs(2219) > 1.e-100 ? 0.84219999999999995 * exp(-tc[1]/2219) : 0.) 
        + (4 == 4 ? exp(-6882 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[5] *= Corr * k_f;
    // (6):  H2 + CO (+M) <=> CH2O (+M)
    k_f = 1.0000000000000002e-06 * 43000000 
               * exp(1.5 * tc[0] - 0.50321666580471969 * 79600 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 5.0699999999999998e+27 
               * exp(-3.4199999999999999  * tc[0] - 0.50321666580471969  * 84350 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(197) > 1.e-100 ? (1.-0.93200000000000005)*exp(-tc[1] / 197) : 0.) 
        + (fabs(1540) > 1.e-100 ? 0.93200000000000005 * exp(-tc[1]/1540) : 0.) 
        + (4 == 4 ? exp(-10300 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[6] *= Corr * k_f;
    // (7):  2 CH3 (+M) <=> C2H6 (+M)
    k_f = 1.0000000000000002e-06 * 21200000000000000 
               * exp(-0.96999999999999997 * tc[0] - 0.50321666580471969 * 620 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    redP = Corr / k_f * 1e-12 * 1.7700000000000001e+50 
               * exp(-9.6699999999999999  * tc[0] - 0.50321666580471969  * 6220 *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (fabs(151) > 1.e-100 ? (1.-0.53249999999999997)*exp(-tc[1] / 151) : 0.) 
        + (fabs(1038) > 1.e-100 ? 0.53249999999999997 * exp(-tc[1]/1038) : 0.) 
        + (4 == 4 ? exp(-4970 * invT) : 0.) );
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[7] *= Corr * k_f;
    // (8):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 5e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    qf[8] *= Corr * k_f;
    // (9):  O + CO + M <=> CO2 + M
    k_f = 1.0000000000000002e-12 * 602000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 3000 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[3] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 3.5 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.5 - 1)*sc[20];
    qf[9] *= Corr * k_f;
    // (10):  H + O2 + M <=> HO2 + M
    k_f = 1.0000000000000002e-12 * 2.8e+18 
               * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0.75 - 1)*sc[11] + ( 1.5 - 1)*sc[12] + ( 1.5 - 1)*sc[18] + ( 0 - 1)*sc[19] + ( 0 - 1)*sc[20];
    qf[10] *= Corr * k_f;
    // (11):  2 H + M <=> H2 + M
    k_f = 1.0000000000000002e-12 * 1e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 0 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.63 - 1)*sc[20];
    qf[11] *= Corr * k_f;
    // (12):  H + OH + M <=> H2O + M
    k_f = 1.0000000000000002e-12 * 2.2e+22 
               * exp(-2 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 3 - 1)*sc[18] + ( 0.38 - 1)*sc[20];
    qf[12] *= Corr * k_f;
    // (13):  HCO + M <=> H + CO + M
    k_f = 1.0000000000000002e-06 * 1.87e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * 17000 * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    qf[13] *= Corr * k_f;
    // (14):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 50000 
               * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * 6290 * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    // (15):  O + HO2 <=> OH + O2
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    // (16):  O + CH2 <=> H + HCO
    k_f = 1.0000000000000002e-06 * 80000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    // (17):  O + CH2(S) <=> H + HCO
    k_f = 1.0000000000000002e-06 * 15000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    // (18):  O + CH3 <=> H + CH2O
    k_f = 1.0000000000000002e-06 * 84300000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    // (19):  O + CH4 <=> OH + CH3
    k_f = 1.0000000000000002e-06 * 1020000000 
               * exp(1.5 * tc[0] - 0.50321666580471969 * 8600 * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    // (20):  O + HCO <=> OH + CO
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    // (21):  O + HCO <=> H + CO2
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    // (22):  O + CH2O <=> OH + HCO
    k_f = 1.0000000000000002e-06 * 39000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 3540 * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    // (23):  O + C2H4 <=> CH3 + HCO
    k_f = 1.0000000000000002e-06 * 19200000 
               * exp(1.8300000000000001 * tc[0] - 0.50321666580471969 * 220 * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    // (24):  O + C2H5 <=> CH3 + CH2O
    k_f = 1.0000000000000002e-06 * 132000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    // (25):  O + C2H6 <=> OH + C2H5
    k_f = 1.0000000000000002e-06 * 89800000 
               * exp(1.9199999999999999 * tc[0] - 0.50321666580471969 * 5690 * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    // (26):  O2 + CO <=> O + CO2
    k_f = 1.0000000000000002e-06 * 2500000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 47800 * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    // (27):  O2 + CH2O <=> HO2 + HCO
    k_f = 1.0000000000000002e-06 * 100000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 40000 * invT);
    Corr  = 1.0;
    qf[27] *= Corr * k_f;
    // (28):  H + 2 O2 <=> HO2 + O2
    k_f = 1.0000000000000002e-12 * 3e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[28] *= Corr * k_f;
    // (29):  H + O2 + H2O <=> HO2 + H2O
    k_f = 1.0000000000000002e-12 * 9.38e+18 
               * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[29] *= Corr * k_f;
    // (30):  H + O2 + N2 <=> HO2 + N2
    k_f = 1.0000000000000002e-12 * 3.75e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[30] *= Corr * k_f;
    // (31):  H + O2 + AR <=> HO2 + AR
    k_f = 1.0000000000000002e-12 * 7e+17 
               * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[31] *= Corr * k_f;
    // (32):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 83000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 14413 * invT);
    Corr  = 1.0;
    qf[32] *= Corr * k_f;
    // (33):  2 H + H2 <=> 2 H2
    k_f = 1.0000000000000002e-12 * 90000000000000000 
               * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[33] *= Corr * k_f;
    // (34):  2 H + H2O <=> H2 + H2O
    k_f = 1.0000000000000002e-12 * 6e+19 
               * exp(-1.25 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[34] *= Corr * k_f;
    // (35):  2 H + CO2 <=> H2 + CO2
    k_f = 1.0000000000000002e-12 * 5.5e+20 
               * exp(-2 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[35] *= Corr * k_f;
    // (36):  H + HO2 <=> O2 + H2
    k_f = 1.0000000000000002e-06 * 28000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 1068 * invT);
    Corr  = 1.0;
    qf[36] *= Corr * k_f;
    // (37):  H + HO2 <=> 2 OH
    k_f = 1.0000000000000002e-06 * 134000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 635 * invT);
    Corr  = 1.0;
    qf[37] *= Corr * k_f;
    // (38):  H + CH4 <=> CH3 + H2
    k_f = 1.0000000000000002e-06 * 660000000 
               * exp(1.6200000000000001 * tc[0] - 0.50321666580471969 * 10840 * invT);
    Corr  = 1.0;
    qf[38] *= Corr * k_f;
    // (39):  H + HCO <=> H2 + CO
    k_f = 1.0000000000000002e-06 * 73400000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[39] *= Corr * k_f;
    // (40):  H + CH2O <=> HCO + H2
    k_f = 1.0000000000000002e-06 * 23000000000 
               * exp(1.05 * tc[0] - 0.50321666580471969 * 3275 * invT);
    Corr  = 1.0;
    qf[40] *= Corr * k_f;
    // (41):  H + CH3O <=> OH + CH3
    k_f = 1.0000000000000002e-06 * 32000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[41] *= Corr * k_f;
    // (42):  H + C2H6 <=> C2H5 + H2
    k_f = 1.0000000000000002e-06 * 115000000 
               * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * 7530 * invT);
    Corr  = 1.0;
    qf[42] *= Corr * k_f;
    // (43):  OH + H2 <=> H + H2O
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * 3430 * invT);
    Corr  = 1.0;
    qf[43] *= Corr * k_f;
    // (44):  2 OH <=> O + H2O
    k_f = 1.0000000000000002e-06 * 35700 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * -2110 * invT);
    Corr  = 1.0;
    qf[44] *= Corr * k_f;
    // (45):  OH + HO2 <=> O2 + H2O
    k_f = 1.0000000000000002e-06 * 29000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * -500 * invT);
    Corr  = 1.0;
    qf[45] *= Corr * k_f;
    // (46):  OH + CH2 <=> H + CH2O
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[46] *= Corr * k_f;
    // (47):  OH + CH2(S) <=> H + CH2O
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[47] *= Corr * k_f;
    // (48):  OH + CH3 <=> CH2 + H2O
    k_f = 1.0000000000000002e-06 * 56000000 
               * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * 5420 * invT);
    Corr  = 1.0;
    qf[48] *= Corr * k_f;
    // (49):  OH + CH3 <=> CH2(S) + H2O
    k_f = 1.0000000000000002e-06 * 25010000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[49] *= Corr * k_f;
    // (50):  OH + CH4 <=> CH3 + H2O
    k_f = 1.0000000000000002e-06 * 100000000 
               * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * 3120 * invT);
    Corr  = 1.0;
    qf[50] *= Corr * k_f;
    // (51):  OH + CO <=> H + CO2
    k_f = 1.0000000000000002e-06 * 47600000 
               * exp(1.228 * tc[0] - 0.50321666580471969 * 70 * invT);
    Corr  = 1.0;
    qf[51] *= Corr * k_f;
    // (52):  OH + HCO <=> H2O + CO
    k_f = 1.0000000000000002e-06 * 50000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[52] *= Corr * k_f;
    // (53):  OH + CH2O <=> HCO + H2O
    k_f = 1.0000000000000002e-06 * 3430000000 
               * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * -447 * invT);
    Corr  = 1.0;
    qf[53] *= Corr * k_f;
    // (54):  OH + C2H6 <=> C2H5 + H2O
    k_f = 1.0000000000000002e-06 * 3540000 
               * exp(2.1200000000000001 * tc[0] - 0.50321666580471969 * 870 * invT);
    Corr  = 1.0;
    qf[54] *= Corr * k_f;
    // (55):  HO2 + CH2 <=> OH + CH2O
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[55] *= Corr * k_f;
    // (56):  HO2 + CH3 <=> O2 + CH4
    k_f = 1.0000000000000002e-06 * 1000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[56] *= Corr * k_f;
    // (57):  HO2 + CH3 <=> OH + CH3O
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[57] *= Corr * k_f;
    // (58):  HO2 + CO <=> OH + CO2
    k_f = 1.0000000000000002e-06 * 150000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 23600 * invT);
    Corr  = 1.0;
    qf[58] *= Corr * k_f;
    // (59):  CH2 + O2 <=> OH + HCO
    k_f = 1.0000000000000002e-06 * 13200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 1500 * invT);
    Corr  = 1.0;
    qf[59] *= Corr * k_f;
    // (60):  CH2 + H2 <=> H + CH3
    k_f = 1.0000000000000002e-06 * 500000 
               * exp(2 * tc[0] - 0.50321666580471969 * 7230 * invT);
    Corr  = 1.0;
    qf[60] *= Corr * k_f;
    // (61):  CH2 + CH3 <=> H + C2H4
    k_f = 1.0000000000000002e-06 * 40000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[61] *= Corr * k_f;
    // (62):  CH2 + CH4 <=> 2 CH3
    k_f = 1.0000000000000002e-06 * 2460000 
               * exp(2 * tc[0] - 0.50321666580471969 * 8270 * invT);
    Corr  = 1.0;
    qf[62] *= Corr * k_f;
    // (63):  CH2(S) + N2 <=> CH2 + N2
    k_f = 1.0000000000000002e-06 * 15000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 600 * invT);
    Corr  = 1.0;
    qf[63] *= Corr * k_f;
    // (64):  CH2(S) + AR <=> CH2 + AR
    k_f = 1.0000000000000002e-06 * 9000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 600 * invT);
    Corr  = 1.0;
    qf[64] *= Corr * k_f;
    // (65):  CH2(S) + O2 <=> H + OH + CO
    k_f = 1.0000000000000002e-06 * 28000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[65] *= Corr * k_f;
    // (66):  CH2(S) + O2 <=> CO + H2O
    k_f = 1.0000000000000002e-06 * 12000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[66] *= Corr * k_f;
    // (67):  CH2(S) + H2 <=> CH3 + H
    k_f = 1.0000000000000002e-06 * 70000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[67] *= Corr * k_f;
    // (68):  CH2(S) + H2O <=> CH2 + H2O
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[68] *= Corr * k_f;
    // (69):  CH2(S) + CH3 <=> H + C2H4
    k_f = 1.0000000000000002e-06 * 12000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * -570 * invT);
    Corr  = 1.0;
    qf[69] *= Corr * k_f;
    // (70):  CH2(S) + CH4 <=> 2 CH3
    k_f = 1.0000000000000002e-06 * 16000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * -570 * invT);
    Corr  = 1.0;
    qf[70] *= Corr * k_f;
    // (71):  CH2(S) + CO <=> CH2 + CO
    k_f = 1.0000000000000002e-06 * 9000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[71] *= Corr * k_f;
    // (72):  CH2(S) + CO2 <=> CH2 + CO2
    k_f = 1.0000000000000002e-06 * 7000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[72] *= Corr * k_f;
    // (73):  CH2(S) + CO2 <=> CO + CH2O
    k_f = 1.0000000000000002e-06 * 14000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[73] *= Corr * k_f;
    // (74):  CH3 + O2 <=> O + CH3O
    k_f = 1.0000000000000002e-06 * 26750000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 28800 * invT);
    Corr  = 1.0;
    qf[74] *= Corr * k_f;
    // (75):  CH3 + O2 <=> OH + CH2O
    k_f = 1.0000000000000002e-06 * 36000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 8940 * invT);
    Corr  = 1.0;
    qf[75] *= Corr * k_f;
    // (76):  2 CH3 <=> H + C2H5
    k_f = 1.0000000000000002e-06 * 4990000000000 
               * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * 10600 * invT);
    Corr  = 1.0;
    qf[76] *= Corr * k_f;
    // (77):  CH3 + HCO <=> CH4 + CO
    k_f = 1.0000000000000002e-06 * 26480000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    Corr  = 1.0;
    qf[77] *= Corr * k_f;
    // (78):  CH3 + CH2O <=> HCO + CH4
    k_f = 1.0000000000000002e-06 * 3320 
               * exp(2.8100000000000001 * tc[0] - 0.50321666580471969 * 5860 * invT);
    Corr  = 1.0;
    qf[78] *= Corr * k_f;
    // (79):  CH3 + C2H6 <=> C2H5 + CH4
    k_f = 1.0000000000000002e-06 * 6140000 
               * exp(1.74 * tc[0] - 0.50321666580471969 * 10450 * invT);
    Corr  = 1.0;
    qf[79] *= Corr * k_f;
    // (80):  HCO + H2O <=> H + CO + H2O
    k_f = 1.0000000000000002e-06 * 2.244e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * 17000 * invT);
    Corr  = 1.0;
    qf[80] *= Corr * k_f;
    // (81):  HCO + O2 <=> HO2 + CO
    k_f = 1.0000000000000002e-06 * 7600000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 400 * invT);
    Corr  = 1.0;
    qf[81] *= Corr * k_f;
    // (82):  CH3O + O2 <=> HO2 + CH2O
    k_f = 1.0000000000000002e-06 * 4.2799999999999999e-13 
               * exp(7.5999999999999996 * tc[0] - 0.50321666580471969 * -3530 * invT);
    Corr  = 1.0;
    qf[82] *= Corr * k_f;
    // (83):  C2H5 + O2 <=> HO2 + C2H4
    k_f = 1.0000000000000002e-06 * 840000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * 3875 * invT);
    Corr  = 1.0;
    qf[83] *= Corr * k_f;

    /*compute the Gibbs free energy */
    double g_RT[21];
    gibbs_d(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    qr[0] *= qf[0] / (exp(g_RT[1] + g_RT[7] - g_RT[9]) * refCinv);
    qr[1] *= qf[1] / (exp(g_RT[1] + g_RT[9] - g_RT[10]) * refCinv);
    qr[2] *= qf[2] / (exp(g_RT[1] + g_RT[13] - g_RT[14]) * refCinv);
    qr[3] *= qf[3] / (exp(g_RT[1] + g_RT[14] - g_RT[15]) * refCinv);
    qr[4] *= qf[4] / (exp(g_RT[1] + g_RT[16] - g_RT[17]) * refCinv);
    qr[5] *= qf[5] / (exp(g_RT[1] + g_RT[17] - g_RT[18]) * refCinv);
    qr[6] *= qf[6] / (exp(g_RT[0] + g_RT[11] - g_RT[14]) * refCinv);
    qr[7] *= qf[7] / (exp(2*g_RT[9] - g_RT[18]) * refCinv);
    qr[8] *= qf[8] / (exp(g_RT[1] + g_RT[2] - g_RT[4]) * refCinv);
    qr[9] *= qf[9] / (exp(g_RT[2] + g_RT[11] - g_RT[12]) * refCinv);
    qr[10] *= qf[10] / (exp(g_RT[1] + g_RT[3] - g_RT[6]) * refCinv);
    qr[11] *= qf[11] / (exp(-g_RT[0] + 2*g_RT[1]) * refCinv);
    qr[12] *= qf[12] / (exp(g_RT[1] + g_RT[4] - g_RT[5]) * refCinv);
    qr[13] *= qf[13] / (exp(-g_RT[1] - g_RT[11] + g_RT[13]) * refC);
    qr[14] *= qf[14] / exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    qr[15] *= qf[15] / exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    qr[16] *= qf[16] / exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13]);
    qr[17] *= qf[17] / exp(-g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13]);
    qr[18] *= qf[18] / exp(-g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14]);
    qr[19] *= qf[19] / exp(g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10]);
    qr[20] *= qf[20] / exp(g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13]);
    qr[21] *= qf[21] / exp(-g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13]);
    qr[22] *= qf[22] / exp(g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14]);
    qr[23] *= qf[23] / exp(g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16]);
    qr[24] *= qf[24] / exp(g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17]);
    qr[25] *= qf[25] / exp(g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18]);
    qr[26] *= qf[26] / exp(-g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12]);
    qr[27] *= qf[27] / exp(g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14]);
    qr[28] *= qf[28] / (exp(g_RT[1] + 2*g_RT[3] - g_RT[3] - g_RT[6]) * refCinv);
    qr[29] *= qf[29] / (exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]) * refCinv);
    qr[30] *= qf[30] / (exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19]) * refCinv);
    qr[31] *= qf[31] / (exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20]) * refCinv);
    qr[32] *= qf[32] / exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    qr[33] *= qf[33] / (exp(g_RT[0] - 2*g_RT[0] + 2*g_RT[1]) * refCinv);
    qr[34] *= qf[34] / (exp(-g_RT[0] + 2*g_RT[1] + g_RT[5] - g_RT[5]) * refCinv);
    qr[35] *= qf[35] / (exp(-g_RT[0] + 2*g_RT[1] + g_RT[12] - g_RT[12]) * refCinv);
    qr[36] *= qf[36] / exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    qr[37] *= qf[37] / exp(g_RT[1] - 2*g_RT[4] + g_RT[6]);
    qr[38] *= qf[38] / exp(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10]);
    qr[39] *= qf[39] / exp(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13]);
    qr[40] *= qf[40] / exp(-g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14]);
    qr[41] *= qf[41] / exp(g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15]);
    qr[42] *= qf[42] / exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18]);
    qr[43] *= qf[43] / exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    qr[44] *= qf[44] / exp(-g_RT[2] + 2*g_RT[4] - g_RT[5]);
    qr[45] *= qf[45] / exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    qr[46] *= qf[46] / exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14]);
    qr[47] *= qf[47] / exp(-g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14]);
    qr[48] *= qf[48] / exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9]);
    qr[49] *= qf[49] / exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9]);
    qr[50] *= qf[50] / exp(g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10]);
    qr[51] *= qf[51] / exp(-g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12]);
    qr[52] *= qf[52] / exp(g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13]);
    qr[53] *= qf[53] / exp(g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14]);
    qr[54] *= qf[54] / exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18]);
    qr[55] *= qf[55] / exp(-g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14]);
    qr[56] *= qf[56] / exp(-g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10]);
    qr[57] *= qf[57] / exp(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15]);
    qr[58] *= qf[58] / exp(-g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12]);
    qr[59] *= qf[59] / exp(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13]);
    qr[60] *= qf[60] / exp(g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9]);
    qr[61] *= qf[61] / exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16]);
    qr[62] *= qf[62] / exp(g_RT[7] - 2*g_RT[9] + g_RT[10]);
    qr[63] *= qf[63] / exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19]);
    qr[64] *= qf[64] / exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20]);
    qr[65] *= qf[65] / (exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11]) * refC);
    qr[66] *= qf[66] / exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11]);
    qr[67] *= qf[67] / exp(g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9]);
    qr[68] *= qf[68] / exp(g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8]);
    qr[69] *= qf[69] / exp(-g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16]);
    qr[70] *= qf[70] / exp(g_RT[8] - 2*g_RT[9] + g_RT[10]);
    qr[71] *= qf[71] / exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11]);
    qr[72] *= qf[72] / exp(-g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12]);
    qr[73] *= qf[73] / exp(g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14]);
    qr[74] *= qf[74] / exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15]);
    qr[75] *= qf[75] / exp(g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14]);
    qr[76] *= qf[76] / exp(-g_RT[1] + 2*g_RT[9] - g_RT[17]);
    qr[77] *= qf[77] / exp(g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13]);
    qr[78] *= qf[78] / exp(g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14]);
    qr[79] *= qf[79] / exp(g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18]);
    qr[80] *= qf[80] / (exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13]) * refC);
    qr[81] *= qf[81] / exp(g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13]);
    qr[82] *= qf[82] / exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15]);
    qr[83] *= qf[83] / exp(g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17]);


    return;
}

/*compute the reaction Jacobian */
__device__ void dwdot_d(double * J, double * sc, double * Tp, int * consP)
{
    double c[21];

    for (int k=0; k<21; k++) {
        c[k] = 1.e6 * sc[k];
    }

    ajacobian_d(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<21; k++) {
        J[462+k] *= 1.e-6;
        J[k*22+21] *= 1.e6;
    }

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void gibbs_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.615080560000000e+03 * invT
            +4.095940888000000e+00
            -3.992015430000000e+00 * tc[0]
            +1.200658760000000e-03 * tc[1]
            -7.696564016666666e-07 * tc[2]
            +3.234277775000000e-10 * tc[3]
            -6.820573500000000e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.600404010000000e+04 * invT
            +2.200146820000000e+00
            -3.762678670000000e+00 * tc[0]
            -4.844360715000000e-04 * tc[1]
            -4.658164016666667e-07 * tc[2]
            +3.209092941666667e-10 * tc[3]
            -8.437085950000000e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.049681630000000e+04 * invT
            +4.967723077000000e+00
            -4.198604110000000e+00 * tc[0]
            +1.183307095000000e-03 * tc[1]
            -1.372160366666667e-06 * tc[2]
            +5.573466508333334e-10 * tc[3]
            -9.715736850000000e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +9.786011000000000e+02 * invT
            -1.104597300000000e+01
            -2.106204000000000e+00 * tc[0]
            -3.608297500000000e-03 * tc[1]
            -8.897453333333333e-07 * tc[2]
            +6.148030000000000e-10 * tc[3]
            -1.037805000000000e-13 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.284162650000000e+04 * invT
            -4.007435600000004e-01
            -4.306465680000000e+00 * tc[0]
            +2.093294460000000e-03 * tc[1]
            -8.285713450000000e-06 * tc[2]
            +4.992721716666666e-09 * tc[3]
            -1.152545020000000e-12 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 19: N2 */
        species[19] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.858657000000000e+03 * invT
            -1.383808430000000e+00
            -3.092887670000000e+00 * tc[0]
            -2.742148580000000e-04 * tc[1]
            -2.108420466666667e-08 * tc[2]
            +7.328846300000000e-12 * tc[3]
            -5.870618800000000e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +4.626360400000000e+04 * invT
            -3.297092110000000e+00
            -2.874101130000000e+00 * tc[0]
            -1.828196460000000e-03 * tc[1]
            +2.348243283333333e-07 * tc[2]
            -2.168162908333333e-11 * tc[3]
            +9.386378350000000e-16 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +5.092599970000000e+04 * invT
            -6.334463270000000e+00
            -2.292038420000000e+00 * tc[0]
            -2.327943185000000e-03 * tc[1]
            +3.353199116666667e-07 * tc[2]
            -3.482550000000000e-11 * tc[3]
            +1.698581825000000e-15 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 14: CH2O */
        species[14] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.278325200000000e+02 * invT
            +8.412240000000000e-01
            -3.770799000000000e+00 * tc[0]
            -3.935748500000000e-03 * tc[1]
            +4.427306666666667e-07 * tc[2]
            -3.287025833333333e-11 * tc[3]
            +1.056308000000000e-15 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.285752000000000e+04 * invT
            -1.150777788000000e+01
            -1.954656420000000e+00 * tc[0]
            -8.698636100000001e-03 * tc[1]
            +1.330344446666667e-06 * tc[2]
            -1.460147408333333e-10 * tc[3]
            +7.482078800000000e-15 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 19: N2 */
        species[19] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 20: AR */
        species[20] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
    }
    return;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void dcvpRdT_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
        /*species 1: H */
        species[1] =
            +7.05332819e-13
            -3.99183928e-15 * tc[1]
            +6.90244896e-18 * tc[2]
            -3.71092933e-21 * tc[3];
        /*species 2: O */
        species[2] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 3: O2 */
        species[3] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 4: OH */
        species[4] =
            -2.40131752e-03
            +9.23587682e-06 * tc[1]
            -1.16434000e-08 * tc[2]
            +5.45645880e-12 * tc[3];
        /*species 5: H2O */
        species[5] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: CH2 */
        species[7] =
            +9.68872143e-04
            +5.58979682e-06 * tc[1]
            -1.15527346e-08 * tc[2]
            +6.74966876e-12 * tc[3];
        /*species 8: CH2(S) */
        species[8] =
            -2.36661419e-03
            +1.64659244e-05 * tc[1]
            -2.00644794e-08 * tc[2]
            +7.77258948e-12 * tc[3];
        /*species 9: CH3 */
        species[9] =
            +2.01095175e-03
            +1.14604371e-05 * tc[1]
            -2.06135228e-08 * tc[2]
            +1.01754294e-11 * tc[3];
        /*species 10: CH4 */
        species[10] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 11: CO */
        species[11] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 13: HCO */
        species[13] =
            -3.24392532e-03
            +2.75598892e-05 * tc[1]
            -3.99432279e-08 * tc[2]
            +1.73507546e-11 * tc[3];
        /*species 14: CH2O */
        species[14] =
            -9.90833369e-03
            +7.46440016e-05 * tc[1]
            -1.13785578e-07 * tc[2]
            +5.27090608e-11 * tc[3];
        /*species 15: CH3O */
        species[15] =
            +7.21659500e-03
            +1.06769440e-05 * tc[1]
            -2.21329080e-08 * tc[2]
            +8.30244000e-12 * tc[3];
        /*species 16: C2H4 */
        species[16] =
            -7.57052247e-03
            +1.14198058e-04 * tc[1]
            -2.07476626e-07 * tc[2]
            +1.07953749e-10 * tc[3];
        /*species 17: C2H5 */
        species[17] =
            -4.18658892e-03
            +9.94285614e-05 * tc[1]
            -1.79737982e-07 * tc[2]
            +9.22036016e-11 * tc[3];
        /*species 18: C2H6 */
        species[18] =
            -5.50154270e-03
            +1.19887658e-04 * tc[1]
            -2.12539886e-07 * tc[2]
            +1.07474308e-10 * tc[3];
        /*species 19: N2 */
        species[19] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
        /*species 20: AR */
        species[20] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
        /*species 1: H */
        species[1] =
            -2.30842973e-11
            +3.23123896e-14 * tc[1]
            -1.42054571e-17 * tc[2]
            +1.99278943e-21 * tc[3];
        /*species 2: O */
        species[2] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 3: O2 */
        species[3] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 4: OH */
        species[4] =
            +5.48429716e-04
            +2.53010456e-07 * tc[1]
            -2.63838467e-10 * tc[2]
            +4.69649504e-14 * tc[3];
        /*species 5: H2O */
        species[5] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: CH2 */
        species[7] =
            +3.65639292e-03
            -2.81789194e-06 * tc[1]
            +7.80538647e-10 * tc[2]
            -7.50910268e-14 * tc[3];
        /*species 8: CH2(S) */
        species[8] =
            +4.65588637e-03
            -4.02383894e-06 * tc[1]
            +1.25371800e-09 * tc[2]
            -1.35886546e-13 * tc[3];
        /*species 9: CH3 */
        species[9] =
            +7.23990037e-03
            -5.97428696e-06 * tc[1]
            +1.78705393e-09 * tc[2]
            -1.86861758e-13 * tc[3];
        /*species 10: CH4 */
        species[10] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 11: CO */
        species[11] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 13: HCO */
        species[13] =
            +4.95695526e-03
            -4.96891226e-06 * tc[1]
            +1.76748533e-09 * tc[2]
            -2.13403484e-13 * tc[3];
        /*species 14: CH2O */
        species[14] =
            +9.20000082e-03
            -8.84517626e-06 * tc[1]
            +3.01923636e-09 * tc[2]
            -3.53542256e-13 * tc[3];
        /*species 15: CH3O */
        species[15] =
            +7.87149700e-03
            -5.31276800e-06 * tc[1]
            +1.18332930e-09 * tc[2]
            -8.45046400e-14 * tc[3];
        /*species 16: C2H4 */
        species[16] =
            +1.46454151e-02
            -1.34215583e-05 * tc[1]
            +4.41668769e-09 * tc[2]
            -5.02824244e-13 * tc[3];
        /*species 17: C2H5 */
        species[17] =
            +1.73972722e-02
            -1.59641334e-05 * tc[1]
            +5.25653067e-09 * tc[2]
            -5.98566304e-13 * tc[3];
        /*species 18: C2H6 */
        species[18] =
            +2.16852677e-02
            -2.00512134e-05 * tc[1]
            +6.64236003e-09 * tc[2]
            -7.60011560e-13 * tc[3];
        /*species 19: N2 */
        species[19] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 20: AR */
        species[20] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cv_R_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +2.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +3.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +1.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +3.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 19: N2 */
        species[19] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +1.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +1.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 10: CH4 */
        species[10] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +2.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +9.54656420e-01
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 19: N2 */
        species[19] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cp_R_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 11: CO */
        species[11] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00
            +7.21659500e-03 * tc[1]
            +5.33847200e-06 * tc[2]
            -7.37763600e-09 * tc[3]
            +2.07561000e-12 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00
            +7.87149700e-03 * tc[1]
            -2.65638400e-06 * tc[2]
            +3.94443100e-10 * tc[3]
            -2.11261600e-14 * tc[4];
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesInternalEnergy_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +2.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +3.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 11: CO */
        species[11] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +1.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +3.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +1.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +1.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 11: CO */
        species[11] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +2.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +9.54656420e-01
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 20: AR */
        species[20] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesEnthalpy_d(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 11: CO */
        species[11] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +2.10620400e+00
            +3.60829750e-03 * tc[1]
            +1.77949067e-06 * tc[2]
            -1.84440900e-09 * tc[3]
            +4.15122000e-13 * tc[4]
            +9.78601100e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: CH2 */
        species[7] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 * invT;
        /*species 8: CH2(S) */
        species[8] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 * invT;
        /*species 9: CH3 */
        species[9] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 10: CH4 */
        species[10] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 11: CO */
        species[11] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 13: HCO */
        species[13] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 14: CH2O */
        species[14] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 15: CH3O */
        species[15] =
            +3.77079900e+00
            +3.93574850e-03 * tc[1]
            -8.85461333e-07 * tc[2]
            +9.86107750e-11 * tc[3]
            -4.22523200e-15 * tc[4]
            +1.27832520e+02 * invT;
        /*species 16: C2H4 */
        species[16] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 17: C2H5 */
        species[17] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 * invT;
        /*species 18: C2H6 */
        species[18] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 19: N2 */
        species[19] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 20: AR */
        species[20] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
    }
    return;
}


/*save molecular weights into array */
__device__ void molecularWeight_d(double * wt)
{
    wt[0] = 2.015940; /*H2 */
    wt[1] = 1.007970; /*H */
    wt[2] = 15.999400; /*O */
    wt[3] = 31.998800; /*O2 */
    wt[4] = 17.007370; /*OH */
    wt[5] = 18.015340; /*H2O */
    wt[6] = 33.006770; /*HO2 */
    wt[7] = 14.027090; /*CH2 */
    wt[8] = 14.027090; /*CH2(S) */
    wt[9] = 15.035060; /*CH3 */
    wt[10] = 16.043030; /*CH4 */
    wt[11] = 28.010550; /*CO */
    wt[12] = 44.009950; /*CO2 */
    wt[13] = 29.018520; /*HCO */
    wt[14] = 30.026490; /*CH2O */
    wt[15] = 31.034460; /*CH3O */
    wt[16] = 28.054180; /*C2H4 */
    wt[17] = 29.062150; /*C2H5 */
    wt[18] = 30.070120; /*C2H6 */
    wt[19] = 28.013400; /*N2 */
    wt[20] = 39.948000; /*AR */

    return;
}

/* get temperature given internal energy in mass units and mass fracs */
__device__ void get_t_given_ey_d_(double * e, double * y_wk, double * t, int * ierr)
{
    const int maxiter = 200;
    const double tol  = 1.e-6;
    double ein  = *e;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    ckubms_d(&tmin, y_wk, &emin);
    ckubms_d(&tmax, y_wk, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        ckcvbs_d(&tmin, y_wk, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        ckcvbs_d(&tmax, y_wk,&cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        ckubms_d(&t1,y_wk,&e1);
        ckcvbs_d(&t1,y_wk,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* get temperature given enthalpy in mass units and mass fracs */
__device__ void get_t_given_hy_d_(double * h, double * y_wk, double * t, int * ierr)
{
    const int maxiter = 200;
    const double tol  = 1.e-6;
    double hin  = *h;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    ckhbms_d(&tmin, y_wk,  &hmin);
    ckhbms_d(&tmax, y_wk,  &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        ckcpbs_d(&tmin, y_wk, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        ckcpbs_d(&tmax, y_wk, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        ckhbms_d(&t1,y_wk,&h1);
        ckcpbs_d(&t1,y_wk,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}


/*compute the reaction Jacobian */
__device__ void ajacobian_d(double * J, double * sc, double T, int consP)
{


    for (int i=0; i<484; i++) {
        J[i] = 0.0;
    }

    double wdot[21];
    for (int k=0; k<21; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 21; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[21];
    gibbs_d(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[21];
    speciesEnthalpy_d(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[21];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 25000000000000000
                * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -0.80000000000000004 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 3.2000000000000002e+27 * exp(-3.1400000000000001 * tc[0] - 0.50321666580471969 * 3.2000000000000002e+27 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -3.1400000000000001 * invT + 0.50321666580471969 * 1230 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.68000000000000005)*exp(-T/78);
    Fcent2 = 0.68000000000000005 * exp(-T/1995);
    Fcent3 = exp(-5590 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/78
        -Fcent2/1995
        + Fcent3*5590*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[9]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
        J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[7];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
        J[31] += dqdci;               /* dwdot[CH3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[117] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[119] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[1];
        J[155] -= dqdci;              /* dwdot[H]/d[CH2] */
        J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
        J[163] += dqdci;              /* dwdot[CH3]/d[CH2] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
        J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[249] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[251] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[271] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[273] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[403] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[405] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[447] -= dqdci;              /* dwdot[CH2]/d[AR] */
        J[449] += dqdci;              /* dwdot[CH3]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[7];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f*sc[1];
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac - k_r;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+7] -= dqdc[k];
            J[22*k+9] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[469] -= dqdT; /* dwdot[CH2]/dT */
    J[471] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 12700000000000000
                * exp(-0.63 * tc[0] - 0.50321666580471969 * 383 * invT);
    dlnkfdT = -0.63 * invT + 0.50321666580471969 *  383  * invT2;
    /* pressure-fall-off */
    k_0 = 2.4769999999999999e+33 * exp(-4.7599999999999998 * tc[0] - 0.50321666580471969 * 2.4769999999999999e+33 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.7599999999999998 * invT + 0.50321666580471969 * 2440 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.78300000000000003)*exp(-T/74);
    Fcent2 = 0.78300000000000003 * exp(-T/2941);
    Fcent3 = exp(-6964 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/74
        -Fcent2/2941
        + Fcent3*6964*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[1] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[10]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[9] -= dqdci;                /* dwdot[CH3]/d[H2] */
        J[10] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[9];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[32] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[120] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[1];
        J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac - k_r;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[273] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[274] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[449] -= dqdci;              /* dwdot[CH3]/d[AR] */
        J[450] += dqdci;              /* dwdot[CH4]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[9];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*sc[1];
        dqdc[10] = 2*dcdc_fac - k_r;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+9] -= dqdc[k];
            J[22*k+10] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[471] -= dqdT; /* dwdot[CH3]/dT */
    J[472] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 1090000000000
                * exp(0.47999999999999998 * tc[0] - 0.50321666580471969 * -260 * invT);
    dlnkfdT = 0.47999999999999998 * invT + 0.50321666580471969 *  -260  * invT2;
    /* pressure-fall-off */
    k_0 = 1.35e+24 * exp(-2.5699999999999998 * tc[0] - 0.50321666580471969 * 1.35e+24 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.5699999999999998 * invT + 0.50321666580471969 * 1425 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.78239999999999998)*exp(-T/271);
    Fcent2 = 0.78239999999999998 * exp(-T/2755);
    Fcent3 = exp(-6570 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/271
        -Fcent2/2755
        + Fcent3*6570*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[1] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[14]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[13] -= q; /* HCO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[13];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f*sc[1];
        J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[300] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[321] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[453] -= dqdci;              /* dwdot[HCO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[13];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac + k_f*sc[1];
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+13] -= dqdc[k];
            J[22*k+14] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = 1.0000000000000002e-06 * 540000000000
                * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * 2600 * invT);
    dlnkfdT = 0.45400000000000001 * invT + 0.50321666580471969 *  2600  * invT2;
    /* pressure-fall-off */
    k_0 = 2.2e+30 * exp(-4.7999999999999998 * tc[0] - 0.50321666580471969 * 2.2e+30 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -4.7999999999999998 * invT + 0.50321666580471969 * 5560 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.75800000000000001)*exp(-T/94);
    Fcent2 = 0.75800000000000001 * exp(-T/1555);
    Fcent3 = exp(-4200 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/94
        -Fcent2/1555
        + Fcent3*4200*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[1] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[15]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[14] -= q; /* CH2O */
    wdot[15] += q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        J[15] += dqdci;               /* dwdot[CH3O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[14];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
        J[37] += dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[125] += dqdci;              /* dwdot[CH3O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[235] += dqdci;              /* dwdot[CH3O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[256] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        J[257] += dqdci;              /* dwdot[CH3O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[278] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[279] += dqdci;              /* dwdot[CH3O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  + k_f*sc[1];
        J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[323] += dqdci;              /* dwdot[CH3O]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  - k_r;
        J[331] -= dqdci;              /* dwdot[H]/d[CH3O] */
        J[344] -= dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[410] -= dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[411] += dqdci;              /* dwdot[CH3O]/d[C2H6] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[14];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac + k_f*sc[1];
        dqdc[15] = dcdc_fac - k_r;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+14] -= dqdc[k];
            J[22*k+15] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[476] -= dqdT; /* dwdot[CH2O]/dT */
    J[477] += dqdT; /* dwdot[CH3O]/dT */

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = 1.0000000000000002e-06 * 1080000000000
                * exp(0.45400000000000001 * tc[0] - 0.50321666580471969 * 1820 * invT);
    dlnkfdT = 0.45400000000000001 * invT + 0.50321666580471969 *  1820  * invT2;
    /* pressure-fall-off */
    k_0 = 1.1999999999999999e+42 * exp(-7.6200000000000001 * tc[0] - 0.50321666580471969 * 1.1999999999999999e+42 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -7.6200000000000001 * invT + 0.50321666580471969 * 6970 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.97529999999999994)*exp(-T/210);
    Fcent2 = 0.97529999999999994 * exp(-T/984);
    Fcent3 = exp(-4374 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/210
        -Fcent2/984
        + Fcent3*4374*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[17];
    Kc = refCinv * exp(g_RT[1] + g_RT[16] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[17]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[16] -= q; /* C2H4 */
    wdot[17] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[16];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[38] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[126] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[236] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[258] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[259] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[280] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[281] += dqdci;              /* dwdot[C2H5]/d[CO2] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[1];
        J[353] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[368] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[369] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[390] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[456] -= dqdci;              /* dwdot[C2H4]/d[AR] */
        J[457] += dqdci;              /* dwdot[C2H5]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[16];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f*sc[1];
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+16] -= dqdc[k];
            J[22*k+17] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[478] -= dqdT; /* dwdot[C2H4]/dT */
    J[479] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = 1.0000000000000002e-06 * 5.21e+17
                * exp(-0.98999999999999999 * tc[0] - 0.50321666580471969 * 1580 * invT);
    dlnkfdT = -0.98999999999999999 * invT + 0.50321666580471969 *  1580  * invT2;
    /* pressure-fall-off */
    k_0 = 1.9900000000000001e+41 * exp(-7.0800000000000001 * tc[0] - 0.50321666580471969 * 1.9900000000000001e+41 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -7.0800000000000001 * invT + 0.50321666580471969 * 6685 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.84219999999999995)*exp(-T/125);
    Fcent2 = 0.84219999999999995 * exp(-T/2219);
    Fcent3 = exp(-6882 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/125
        -Fcent2/2219
        + Fcent3*6882*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[18];
    Kc = refCinv * exp(g_RT[1] + g_RT[17] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[18]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[17] -= q; /* C2H5 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[39] -= dqdci;               /* dwdot[C2H5]/d[H] */
        J[40] += dqdci;               /* dwdot[C2H6]/d[H] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[127] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[237] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[259] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[281] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H5] */
        dqdci =  + k_f*sc[1];
        J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[392] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac - k_r;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[413] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[457] -= dqdci;              /* dwdot[C2H5]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[17];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[1];
        dqdc[18] = 3*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+17] -= dqdc[k];
            J[22*k+18] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[479] -= dqdT; /* dwdot[C2H5]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = 1.0000000000000002e-06 * 43000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * 79600 * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  79600  * invT2;
    /* pressure-fall-off */
    k_0 = 5.0699999999999998e+27 * exp(-3.4199999999999999 * tc[0] - 0.50321666580471969 * 5.0699999999999998e+27 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -3.4199999999999999 * invT + 0.50321666580471969 * 84350 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.93200000000000005)*exp(-T/197);
    Fcent2 = 0.93200000000000005 * exp(-T/1540);
    Fcent3 = exp(-10300 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/197
        -Fcent2/1540
        + Fcent3*10300*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[0] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[14]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[11] -= q; /* CO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac + k_f*sc[11];
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[220] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac + k_f*sc[0];
        J[242] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[264] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[308] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[319] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac;
        J[396] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H2]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac + k_f*sc[11];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac + k_f*sc[0];
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+0] -= dqdc[k];
            J[22*k+11] -= dqdc[k];
            J[22*k+14] += dqdc[k];
        }
    }
    J[462] -= dqdT; /* dwdot[H2]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 8: 2 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[9]*sc[9];
    k_f = 1.0000000000000002e-06 * 21200000000000000
                * exp(-0.96999999999999997 * tc[0] - 0.50321666580471969 * 620 * invT);
    dlnkfdT = -0.96999999999999997 * invT + 0.50321666580471969 *  620  * invT2;
    /* pressure-fall-off */
    k_0 = 1.7700000000000001e+50 * exp(-9.6699999999999999 * tc[0] - 0.50321666580471969 * 1.7700000000000001e+50 * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -9.6699999999999999 * invT + 0.50321666580471969 * 6220 * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.53249999999999997)*exp(-T/151);
    Fcent2 = 0.53249999999999997 * exp(-T/1038);
    Fcent3 = exp(-4970 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/151
        -Fcent2/1038
        + Fcent3*4970*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[18];
    Kc = refCinv * exp(2*g_RT[9] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[9]) + (h_RT[18]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[9] += -2 * dqdci;           /* dwdot[CH3]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[119] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2*sc[9];
        J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[216] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*dcdc_fac;
        J[229] += -2 * dqdci;         /* dwdot[CH3]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[251] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[273] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*dcdc_fac - k_r;
        J[405] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[449] += -2 * dqdci;         /* dwdot[CH3]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*2*sc[9];
        dqdc[10] = 2*dcdc_fac;
        dqdc[11] = 1.5*dcdc_fac;
        dqdc[12] = 2*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = 3*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = 0.69999999999999996*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+9] += -2 * dqdc[k];
            J[22*k+18] += dqdc[k];
        }
    }
    J[471] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 9: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.69999999999999996 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 5e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[4]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[24] -= dqdci;               /* dwdot[O]/d[H] */
        J[26] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[45] -= dqdci;               /* dwdot[H]/d[O] */
        J[46] -= dqdci;               /* dwdot[O]/d[O] */
        J[48] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[89] -= dqdci;               /* dwdot[H]/d[OH] */
        J[90] -= dqdci;               /* dwdot[O]/d[OH] */
        J[92] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[114] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[246] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[444] += dqdci;              /* dwdot[OH]/d[AR] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = 6*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = 1.5*q_nocor;
        dqdc[12] = 2*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.69999999999999996*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+2] -= dqdc[k];
            J[22*k+4] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[466] += dqdT; /* dwdot[OH]/dT */

    /*reaction 10: O + CO + M <=> CO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[3] + ( 6 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 3.5 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.5 - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = 1.0000000000000002e-12 * 602000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 3000 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3000  * invT2;
    /* reverse */
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[2] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[12]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[12] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[11];
        J[46] -= dqdci;               /* dwdot[O]/d[O] */
        J[55] -= dqdci;               /* dwdot[CO]/d[O] */
        J[56] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[O2] */
        dqdci = (6 - 1)*q_nocor;
        J[68] -= dqdci;               /* dwdot[O]/d[O2] */
        J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
        J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*q_nocor;
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[122] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[232] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*q_nocor + k_f*sc[2];
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.5 - 1)*q_nocor - k_r;
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[408] += dqdci;              /* dwdot[CO2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.5 - 1)*q_nocor;
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[452] += dqdci;              /* dwdot[CO2]/d[AR] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[11];
        dqdc[3] = 6*q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = 6*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = 1.5*q_nocor + k_f*sc[2];
        dqdc[12] = 3.5*q_nocor - k_r;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.5*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+11] -= dqdc[k];
            J[22*k+12] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[474] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 11: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0.75 - 1)*sc[11] + ( 1.5 - 1)*sc[12] + ( 1.5 - 1)*sc[18] + ( 0 - 1)*sc[19] + ( 0 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.8e+18
                * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -0.85999999999999999 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[25] -= dqdci;               /* dwdot[O2]/d[H] */
        J[28] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (0 - 1)*q_nocor + k_f*sc[1];
        J[67] -= dqdci;               /* dwdot[H]/d[O2] */
        J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (0.75 - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (1.5 - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[270] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (1.5 - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[399] -= dqdci;              /* dwdot[O2]/d[C2H6] */
        J[402] += dqdci;              /* dwdot[HO2]/d[C2H6] */
        /* d()/d[N2] */
        dqdci = (0 - 1)*q_nocor;
        J[419] -= dqdci;              /* dwdot[H]/d[N2] */
        J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
        /* d()/d[AR] */
        dqdci = (0 - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] =  + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = 0.75*q_nocor;
        dqdc[12] = 1.5*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 1.5*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+3] -= dqdc[k];
            J[22*k+6] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[465] -= dqdT; /* dwdot[O2]/dT */
    J[468] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 12: 2 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 0 - 1)*sc[12] + ( 3 - 1)*sc[18] + ( 0.63 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[1];
    k_f = 1.0000000000000002e-12 * 1e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + 2*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1]) + (h_RT[0]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (0 - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2*sc[1];
        J[22] += dqdci;               /* dwdot[H2]/d[H] */
        J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[221] += -2 * dqdci;         /* dwdot[H]/d[CH4] */
        /* d()/d[CO2] */
        dqdci = (0 - 1)*q_nocor;
        J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[397] += -2 * dqdci;         /* dwdot[H]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.63 - 1)*q_nocor;
        J[440] += dqdci;              /* dwdot[H2]/d[AR] */
        J[441] += -2 * dqdci;         /* dwdot[H]/d[AR] */
    }
    else {
        dqdc[0] =  - k_r;
        dqdc[1] = q_nocor + k_f*2*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.63*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+0] += dqdc[k];
            J[22*k+1] += -2 * dqdc[k];
        }
    }
    J[462] += dqdT; /* dwdot[H2]/dT */
    J[463] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 13: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 3 - 1)*sc[18] + ( 0.38 - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-12 * 2.2e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (0.72999999999999998 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[26] -= dqdci;               /* dwdot[OH]/d[H] */
        J[27] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[89] -= dqdci;               /* dwdot[H]/d[OH] */
        J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (3.6499999999999999 - 1)*q_nocor - k_r;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (0.38 - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[444] -= dqdci;              /* dwdot[OH]/d[AR] */
        J[445] += dqdci;              /* dwdot[H2O]/d[AR] */
    }
    else {
        dqdc[0] = 0.72999999999999998*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = 3.6499999999999999*q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = 0.38*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+4] -= dqdc[k];
            J[22*k+5] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[466] -= dqdT; /* dwdot[OH]/dT */
    J[467] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 14: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[5] + ( 2 - 1)*sc[10] + ( 1.5 - 1)*sc[11] + ( 2 - 1)*sc[12] + ( 3 - 1)*sc[18];
    /* forward */
    phi_f = sc[13];
    k_f = 1.0000000000000002e-06 * 1.87e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * 17000 * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  17000  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = refC * exp(-g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13]) + (h_RT[1] + h_RT[11]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[11] += dqdci;               /* dwdot[CO]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[11];
        J[23] += dqdci;               /* dwdot[H]/d[H] */
        J[33] += dqdci;               /* dwdot[CO]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[111] += dqdci;              /* dwdot[H]/d[H2O] */
        J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (2 - 1)*q_nocor;
        J[221] += dqdci;              /* dwdot[H]/d[CH4] */
        J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (1.5 - 1)*q_nocor - k_r*sc[1];
        J[243] += dqdci;              /* dwdot[H]/d[CO] */
        J[253] += dqdci;              /* dwdot[CO]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (2 - 1)*q_nocor;
        J[265] += dqdci;              /* dwdot[H]/d[CO2] */
        J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[287] += dqdci;              /* dwdot[H]/d[HCO] */
        J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[C2H6] */
        dqdci = (3 - 1)*q_nocor;
        J[397] += dqdci;              /* dwdot[H]/d[C2H6] */
        J[407] += dqdci;              /* dwdot[CO]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor - k_r*sc[11];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 2*q_nocor;
        dqdc[11] = 1.5*q_nocor - k_r*sc[1];
        dqdc[12] = 2*q_nocor;
        dqdc[13] = q_nocor + k_f;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = 3*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] += dqdc[k];
            J[22*k+11] += dqdc[k];
            J[22*k+13] -= dqdc[k];
        }
    }
    J[463] += dqdT; /* dwdot[H]/dT */
    J[473] += dqdT; /* dwdot[CO]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 15: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-06 * 50000
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * 6290 * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  6290  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[4] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[44] -= dqdci;               /* dwdot[H2]/d[O] */
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 16: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[47] += dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[50] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[68] -= dqdci;               /* dwdot[O]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[134] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 17: O + CH2 <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 80000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[51] -= dqdci;               /* dwdot[CH2]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[2];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[156] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 18: O + CH2(S) <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = 1.0000000000000002e-06 * 15000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[8] -= q; /* CH2(S) */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[52] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[178] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[189] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[294] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 19: O + CH3 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[9];
    k_f = 1.0000000000000002e-06 * 84300000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[9]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 20: O + CH4 <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = 1.0000000000000002e-06 * 1020000000
                * exp(1.5 * tc[0] - 0.50321666580471969 * 8600 * invT);
    dlnkfdT = 1.5 * invT + 0.50321666580471969 *  8600  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[54] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 21: O + HCO <=> OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[55] += dqdci;               /* dwdot[CO]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[244] -= dqdci;              /* dwdot[O]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 22: O + HCO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[12] += q; /* CO2 */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[298] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 23: O + CH2O <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = 1.0000000000000002e-06 * 39000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 3540 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3540  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[58] -= dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = 1.0000000000000002e-06 * 19200000
                * exp(1.8300000000000001 * tc[0] - 0.50321666580471969 * 220 * invT);
    dlnkfdT = 1.8300000000000001 * invT + 0.50321666580471969 *  220  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[13];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[9] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[13] += q; /* HCO */
    wdot[16] -= q; /* C2H4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[60] -= dqdci;               /* dwdot[C2H4]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[13];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[214] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[9];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[295] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[302] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[354] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[361] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[365] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[368] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = 1.0000000000000002e-06 * 132000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[14];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[9] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[14] += q; /* CH2O */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[61] -= dqdci;               /* dwdot[C2H5]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[14];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[215] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[9];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[325] -= dqdci;              /* dwdot[C2H5]/d[CH2O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[383] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[388] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = 1.0000000000000002e-06 * 89800000
                * exp(1.9199999999999999 * tc[0] - 0.50321666580471969 * 5690 * invT);
    dlnkfdT = 1.9199999999999999 * invT + 0.50321666580471969 *  5690  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[17];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[4] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[61] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[62] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[4];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[378] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 27: O2 + CO <=> O + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = 1.0000000000000002e-06 * 2500000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 47800 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  47800  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[12];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[2] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[12];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[55] -= dqdci;               /* dwdot[CO]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[244] += dqdci;              /* dwdot[O]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[266] += dqdci;              /* dwdot[O]/d[CO2] */
    J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = 1.0000000000000002e-06 * 100000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 40000 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  40000  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[80] -= dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[13];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[145] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[146] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 29: H + 2 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[3];
    k_f = 1.0000000000000002e-12 * 3e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[1] + 2*g_RT[3] - g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + 2*h_RT[3]) + (h_RT[3] + h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[3];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2*sc[3] - k_r*sc[6];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1.0000000000000002e-12 * 9.38e+18
                * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -0.76000000000000001 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[5]) + (h_RT[5] + h_RT[6]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[5];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = 1.0000000000000002e-12 * 3.75e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[19];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[19]) + (h_RT[6] + h_RT[19]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[19];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[19];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[19];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[419] -= dqdci;              /* dwdot[H]/d[N2] */
    J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
    J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = 1.0000000000000002e-12 * 7e+17
                * exp(-0.80000000000000004 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -0.80000000000000004 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[20];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[20]) + (h_RT[6] + h_RT[20]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[20];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[20];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[20];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[441] -= dqdci;              /* dwdot[H]/d[AR] */
    J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
    J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 33: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 83000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 14413 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  14413  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[24] += dqdci;               /* dwdot[O]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[45] -= dqdci;               /* dwdot[H]/d[O] */
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 34: 2 H + H2 <=> 2 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[1]*sc[1];
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[0];
    Kc = refCinv * exp(g_RT[0] - 2*g_RT[0] + 2*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + 2*h_RT[1]) + (2*h_RT[0]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*2*sc[0];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0]*2*sc[1];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 35: 2 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[1]*sc[5];
    k_f = 1.0000000000000002e-12 * 6e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = refCinv * exp(-g_RT[0] + 2*g_RT[1] + g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[5]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[1]*sc[5];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*sc[0];
    J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 36: 2 H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[1]*sc[12];
    k_f = 1.0000000000000002e-12 * 5.5e+20
                * exp(-2 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = refCinv * exp(-g_RT[0] + 2*g_RT[1] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[1] + h_RT[12]) + (h_RT[0] + h_RT[12]) + 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2*sc[1]*sc[12];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[1]*sc[1] - k_r*sc[0];
    J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 37: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 1068 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  1068  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[3] += dqdci;                /* dwdot[O2]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] += dqdci;               /* dwdot[O2]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[66] += dqdci;               /* dwdot[H2]/d[O2] */
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[132] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 38: H + HO2 <=> 2 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 134000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 635 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  635  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[4];
    Kc = exp(g_RT[1] - 2*g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (2*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2*sc[4];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[136] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 660000000
                * exp(1.6200000000000001 * tc[0] - 0.50321666580471969 * 10840 * invT);
    dlnkfdT = 1.6200000000000001 * invT + 0.50321666580471969 *  10840  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[0] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    J[10] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[32] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[198] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 40: H + HCO <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = 1.0000000000000002e-06 * 73400000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[11];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[11];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[11] += dqdci;               /* dwdot[CO]/d[H2] */
    J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[242] += dqdci;              /* dwdot[H2]/d[CO] */
    J[243] -= dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 41: H + CH2O <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = 1.0000000000000002e-06 * 23000000000
                * exp(1.05 * tc[0] - 0.50321666580471969 * 3275 * invT);
    dlnkfdT = 1.05 * invT + 0.50321666580471969 *  3275  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[0] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[13] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[308] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 42: H + CH3O <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = 1.0000000000000002e-06 * 32000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[15] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[37] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[331] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[339] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = 1.0000000000000002e-06 * 115000000
                * exp(1.8999999999999999 * tc[0] - 0.50321666580471969 * 7530 * invT);
    dlnkfdT = 1.8999999999999999 * invT + 0.50321666580471969 *  7530  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[0] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
    J[18] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[40] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[0];
    J[374] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 44: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * 3430 * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[27] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 45: 2 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[4];
    k_f = 1.0000000000000002e-06 * 35700
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * -2110 * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  -2110  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + 2*g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[4]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[48] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[49] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2*sc[4];
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[92] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[112] += dqdci;              /* dwdot[O]/d[H2O] */
    J[114] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[466] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = 1.0000000000000002e-06 * 29000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * -500 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -500  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[113] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[116] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[137] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 47: OH + CH2 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 56000000
                * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * 5420 * invT);
    dlnkfdT = 1.6000000000000001 * invT + 0.50321666580471969 *  5420  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] += q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[95] += dqdci;               /* dwdot[CH2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[7];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[159] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[205] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = 1.0000000000000002e-06 * 25010000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[8] += q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[96] += dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[206] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = 1.0000000000000002e-06 * 100000000
                * exp(1.6000000000000001 * tc[0] - 0.50321666580471969 * 3120 * invT);
    dlnkfdT = 1.6000000000000001 * invT + 0.50321666580471969 *  3120  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[119] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[120] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[4];
    J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 52: OH + CO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = 1.0000000000000002e-06 * 47600000
                * exp(1.228 * tc[0] - 0.50321666580471969 * 70 * invT);
    dlnkfdT = 1.228 * invT + 0.50321666580471969 *  70  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[33] -= dqdci;               /* dwdot[CO]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[268] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 53: OH + HCO <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = 1.0000000000000002e-06 * 50000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = 1.0000000000000002e-06 * 3430000000
                * exp(1.1799999999999999 * tc[0] - 0.50321666580471969 * -447 * invT);
    dlnkfdT = 1.1799999999999999 * invT + 0.50321666580471969 *  -447  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[13];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[123] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[5];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[313] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = 1.0000000000000002e-06 * 3540000
                * exp(2.1200000000000001 * tc[0] - 0.50321666580471969 * 870 * invT);
    dlnkfdT = 2.1200000000000001 * invT + 0.50321666580471969 *  870  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[18]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    J[128] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[5];
    J[378] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[379] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[7];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[7];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[139] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[160] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[314] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[10];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[3] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[76] += dqdci;               /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[142] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[201] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[3];
    J[223] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[226] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] += dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[147] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[336] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[11];
    k_f = 1.0000000000000002e-06 * 150000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 23600 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  23600  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[11]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[144] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[270] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 13200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 1500 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  1500  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[7];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[73] -= dqdci;               /* dwdot[CH2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[3];
    J[157] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 500000
                * exp(2 * tc[0] - 0.50321666580471969 * 7230 * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  7230  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[0];
    J[154] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[9];
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[9];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    J[170] += dqdci;              /* dwdot[C2H4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[359] -= dqdci;              /* dwdot[CH2]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 63: CH2 + CH4 <=> 2 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[10];
    k_f = 1.0000000000000002e-06 * 2460000
                * exp(2 * tc[0] - 0.50321666580471969 * 8270 * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  8270  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[9];
    Kc = exp(g_RT[7] - 2*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[10]) + (2*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* CH2 */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[10];
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[164] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2*sc[9];
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[7];
    J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[19];
    k_f = 1.0000000000000002e-06 * 15000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 600 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  600  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[19];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[19]) + (h_RT[7] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[19];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[425] += dqdci;              /* dwdot[CH2]/d[N2] */
    J[426] -= dqdci;              /* dwdot[CH2(S)]/d[N2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[20];
    k_f = 1.0000000000000002e-06 * 9000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 600 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  600  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[20];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[20]) + (h_RT[7] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[20];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[20];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[447] += dqdci;              /* dwdot[CH2]/d[AR] */
    J[448] -= dqdci;              /* dwdot[CH2(S)]/d[AR] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[1] + h_RT[4] + h_RT[11]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[67] += dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[180] += dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = 1.0000000000000002e-06 * 12000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[5] += q; /* H2O */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[8];
    k_f = 1.0000000000000002e-06 * 70000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[8] -= dqdci;                /* dwdot[CH2(S)]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[0];
    J[176] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[9];
    k_f = 1.0000000000000002e-06 * 12000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * -570 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -570  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[9];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[192] += dqdci;              /* dwdot[C2H4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[360] -= dqdci;              /* dwdot[CH2(S)]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 71: CH2(S) + CH4 <=> 2 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = 1.0000000000000002e-06 * 16000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * -570 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -570  * invT2;
    /* reverse */
    phi_r = sc[9]*sc[9];
    Kc = exp(g_RT[8] - 2*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (2*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[10];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[186] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2*sc[9];
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[8];
    J[228] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = 1.0000000000000002e-06 * 9000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[11];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[7] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[11];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[11];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[249] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 7000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[7]*sc[12];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[7] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[12];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[271] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = 1.0000000000000002e-06 * 14000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[11]*sc[14];
    Kc = exp(g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[11] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    wdot[12] -= q; /* CO2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[188] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[14];
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[254] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8];
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[11];
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[319] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[320] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[474] -= dqdT;               /* dwdot[CO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 26750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 28800 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  28800  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[15];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[59] += dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[81] += dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[200] += dqdci;              /* dwdot[O]/d[CH3] */
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[2];
    J[332] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = 1.0000000000000002e-06 * 36000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 8940 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  8940  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 77: 2 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[9];
    k_f = 1.0000000000000002e-06 * 4990000000000
                * exp(0.10000000000000001 * tc[0] - 0.50321666580471969 * 10600 * invT);
    dlnkfdT = 0.10000000000000001 * invT + 0.50321666580471969 *  10600  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + 2*g_RT[9] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2*h_RT[9]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[17] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[31] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2*sc[9];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[375] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[383] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[471] += -2 * dqdT;          /* dwdot[CH3]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[13];
    k_f = 1.0000000000000002e-06 * 26480000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 0 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[11];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[13]) + (h_RT[10] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[13];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[209] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[211] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[11];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[9];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[14];
    k_f = 1.0000000000000002e-06 * 3320
                * exp(2.8100000000000001 * tc[0] - 0.50321666580471969 * 5860 * invT);
    dlnkfdT = 2.8100000000000001 * invT + 0.50321666580471969 *  5860  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[13];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[14]) + (h_RT[10] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[212] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[13];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[233] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[10];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[9];
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[318] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[18];
    k_f = 1.0000000000000002e-06 * 6140000
                * exp(1.74 * tc[0] - 0.50321666580471969 * 10450 * invT);
    dlnkfdT = 1.74 * invT + 0.50321666580471969 *  10450  * invT2;
    /* reverse */
    phi_r = sc[10]*sc[17];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[18]) + (h_RT[10] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[18];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[216] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[17];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    J[238] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[10];
    J[383] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[384] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[9];
    J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = 1.0000000000000002e-06 * 2.244e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * 17000 * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  17000  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[1] + h_RT[5] + h_RT[11]) - 1);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[13] - k_r*sc[1]*sc[11];
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[5];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = 1.0000000000000002e-06 * 7600000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 400 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  400  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    J[79] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[11];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[145] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = 1.0000000000000002e-06 * 4.2799999999999999e-13
                * exp(7.5999999999999996 * tc[0] - 0.50321666580471969 * -3530 * invT);
    dlnkfdT = 7.5999999999999996 * invT + 0.50321666580471969 *  -3530  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[14] += q; /* CH2O */
    wdot[15] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[81] -= dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[147] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[323] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[336] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[344] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = 1.0000000000000002e-06 * 840000000000
                * exp(0 * tc[0] - 0.50321666580471969 * 3875 * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3875  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[16];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[6] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[16] += q; /* C2H4 */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[82] += dqdci;               /* dwdot[C2H4]/d[O2] */
    J[83] -= dqdci;               /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[16];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[148] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[149] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[6];
    J[355] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[358] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[369] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[3];
    J[377] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[380] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[390] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

    double c_R[21], dcRdT[21], e_RT[21];
    double * eh_RT;
    if (consP) {
        cp_R_d(c_R, tc);
        dcvpRdT_d(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R_d(c_R, tc);
        dcvpRdT_d(dcRdT, tc);
        speciesInternalEnergy_d(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 21; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[462+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 21; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 21; ++m) {
            dehmixdc += eh_RT[m]*J[k*22+m];
        }
        J[k*22+21] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[483] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}

/* End of file  */
