#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <AMReX_Gpu.H>

extern "C"
{
AMREX_GPU_HOST_DEVICE void get_imw(amrex_real imw_new[]);
AMREX_GPU_HOST_DEVICE void get_mw(amrex_real mw_new[]);
void egtransetEPS(amrex_real *  EPS);
void egtransetSIG(amrex_real* SIG);
void atomicWeight(amrex_real *  awt);
void molecularWeight(amrex_real *  wt);
AMREX_GPU_HOST_DEVICE void gibbs(amrex_real *  species, amrex_real *  tc);
AMREX_GPU_HOST_DEVICE void helmholtz(amrex_real *  species, amrex_real *  tc);
AMREX_GPU_HOST_DEVICE void speciesInternalEnergy(amrex_real *  species, amrex_real *  tc);
AMREX_GPU_HOST_DEVICE void speciesEnthalpy(amrex_real *  species, amrex_real *  tc);
AMREX_GPU_HOST_DEVICE void speciesEntropy(amrex_real *  species, amrex_real *  tc);
AMREX_GPU_HOST_DEVICE void cp_R(amrex_real *  species, amrex_real *  tc);
AMREX_GPU_HOST_DEVICE void cv_R(amrex_real *  species, amrex_real *  tc);
void equilibriumConstants(amrex_real *  kc, amrex_real *  g_RT, amrex_real T);
AMREX_GPU_HOST_DEVICE void productionRate(amrex_real *  wdot, amrex_real *  sc, amrex_real T);
AMREX_GPU_HOST_DEVICE void comp_qfqr(amrex_real *  q_f, amrex_real *  q_r, amrex_real *  sc, amrex_real *  tc, amrex_real invT);
#ifndef AMREX_USE_CUDA
void comp_k_f(amrex_real *  tc, amrex_real invT, amrex_real *  k_f);
void comp_Kc(amrex_real *  tc, amrex_real invT, amrex_real *  Kc);
#else
AMREX_GPU_DEVICE amrex_real W_spec_d(amrex_real Q[84], int ispec);
AMREX_GPU_DEVICE amrex_real Q_reac_d(amrex_real rho, amrex_real T, amrex_real Y[21], int reacID);
#endif
AMREX_GPU_HOST_DEVICE void progressRate(amrex_real *  qdot, amrex_real *  speciesConc, amrex_real T);
AMREX_GPU_HOST_DEVICE void progressRateFR(amrex_real *  q_f, amrex_real *  q_r, amrex_real *  speciesConc, amrex_real T);
void CKINIT();
void CKFINALIZE();
void GET_REACTION_MAP(int *  rmap);
void SetAllDefaults();
#ifdef AMREX_USE_CUDA
void AllocateOnDevice();
void DeallocateOnDevice();
#endif
void CKINDX(int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, amrex_real *  rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, amrex_real *  rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(amrex_real *  ru, amrex_real *  ruc, amrex_real *  pa);
void CKPX(amrex_real *  rho, amrex_real *  T, amrex_real *  x, amrex_real *  P);
AMREX_GPU_HOST_DEVICE void CKPY(amrex_real *  rho, amrex_real *  T, amrex_real *  y, amrex_real *  P);
void CKPC(amrex_real *  rho, amrex_real *  T, amrex_real *  c, amrex_real *  P);
void CKRHOX(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  rho);
void CKRHOY(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  rho);
void CKRHOC(amrex_real *  P, amrex_real *  T, amrex_real *  c, amrex_real *  rho);
void CKWT(amrex_real *  wt);
void CKAWT(amrex_real *  awt);
void CKMMWY(amrex_real *  y, amrex_real *  wtm);
void CKMMWX(amrex_real *  x, amrex_real *  wtm);
void CKMMWC(amrex_real *  c, amrex_real *  wtm);
AMREX_GPU_HOST_DEVICE void CKYTX(amrex_real *  y, amrex_real *  x);
void CKYTCP(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  c);
AMREX_GPU_HOST_DEVICE void CKYTCR(amrex_real *  rho, amrex_real *  T, amrex_real *  y, amrex_real *  c);
void CKXTY(amrex_real *  x, amrex_real *  y);
void CKXTCP(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  c);
void CKXTCR(amrex_real *  rho, amrex_real *  T, amrex_real *  x, amrex_real *  c);
void CKCTX(amrex_real *  c, amrex_real *  x);
void CKCTY(amrex_real *  c, amrex_real *  y);
void CKCPOR(amrex_real *  T, amrex_real *  cpor);
void CKHORT(amrex_real *  T, amrex_real *  hort);
void CKSOR(amrex_real *  T, amrex_real *  sor);
void CKCVML(amrex_real *  T, amrex_real *  cvml);
void CKCPML(amrex_real *  T, amrex_real *  cvml);
void CKUML(amrex_real *  T, amrex_real *  uml);
void CKHML(amrex_real *  T, amrex_real *  uml);
void CKGML(amrex_real *  T, amrex_real *  gml);
void CKAML(amrex_real *  T, amrex_real *  aml);
void CKSML(amrex_real *  T, amrex_real *  sml);
AMREX_GPU_HOST_DEVICE void CKCVMS(amrex_real *  T, amrex_real *  cvms);
AMREX_GPU_HOST_DEVICE void CKCPMS(amrex_real *  T, amrex_real *  cvms);
AMREX_GPU_HOST_DEVICE void CKUMS(amrex_real *  T, amrex_real *  ums);
AMREX_GPU_HOST_DEVICE void CKHMS(amrex_real *  T, amrex_real *  ums);
void CKGMS(amrex_real *  T, amrex_real *  gms);
void CKAMS(amrex_real *  T, amrex_real *  ams);
void CKSMS(amrex_real *  T, amrex_real *  sms);
void CKCPBL(amrex_real *  T, amrex_real *  x, amrex_real *  cpbl);
void CKCPBS(amrex_real *  T, amrex_real *  y, amrex_real *  cpbs);
void CKCVBL(amrex_real *  T, amrex_real *  x, amrex_real *  cpbl);
AMREX_GPU_HOST_DEVICE void CKCVBS(amrex_real *  T, amrex_real *  y, amrex_real *  cpbs);
void CKHBML(amrex_real *  T, amrex_real *  x, amrex_real *  hbml);
void CKHBMS(amrex_real *  T, amrex_real *  y, amrex_real *  hbms);
void CKUBML(amrex_real *  T, amrex_real *  x, amrex_real *  ubml);
AMREX_GPU_HOST_DEVICE void CKUBMS(amrex_real *  T, amrex_real *  y, amrex_real *  ubms);
void CKSBML(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  sbml);
void CKSBMS(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  sbms);
void CKGBML(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  gbml);
void CKGBMS(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  gbms);
void CKABML(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  abml);
void CKABMS(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  abms);
AMREX_GPU_HOST_DEVICE void CKWC(amrex_real *  T, amrex_real *  C, amrex_real *  wdot);
void CKWYP(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  wdot);
void CKWXP(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  wdot);
AMREX_GPU_HOST_DEVICE void CKWYR(amrex_real *  rho, amrex_real *  T, amrex_real *  y, amrex_real *  wdot);
void CKWXR(amrex_real *  rho, amrex_real *  T, amrex_real *  x, amrex_real *  wdot);
void CKQC(amrex_real *  T, amrex_real *  C, amrex_real *  qdot);
void CKKFKR(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  q_f, amrex_real *  q_r);
void CKQYP(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  qdot);
void CKQXP(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  qdot);
void CKQYR(amrex_real *  rho, amrex_real *  T, amrex_real *  y, amrex_real *  qdot);
void CKQXR(amrex_real *  rho, amrex_real *  T, amrex_real *  x, amrex_real *  qdot);
void CKNU(int * kdim, int * nuki);
void CKNCF(int * mdim, int * ncf);
void CKABE(amrex_real *  a, amrex_real *  b, amrex_real *  e );
void CKEQC(amrex_real *  T, amrex_real *  C , amrex_real *  eqcon );
void CKEQYP(amrex_real *  P, amrex_real *  T, amrex_real *  y, amrex_real *  eqcon);
void CKEQXP(amrex_real *  P, amrex_real *  T, amrex_real *  x, amrex_real *  eqcon);
void CKEQYR(amrex_real *  rho, amrex_real *  T, amrex_real *  y, amrex_real *  eqcon);
void CKEQXR(amrex_real *  rho, amrex_real *  T, amrex_real *  x, amrex_real *  eqcon);
AMREX_GPU_HOST_DEVICE void DWDOT(amrex_real *  J, amrex_real *  sc, amrex_real *  T, int * consP);
AMREX_GPU_HOST_DEVICE void DWDOT_PRECOND(amrex_real *  J, amrex_real *  sc, amrex_real *  Tp, int * HP);
AMREX_GPU_HOST_DEVICE void SLJ_PRECOND_CSC(amrex_real *  Jsps, int * indx, int * len, amrex_real * sc, amrex_real * Tp, int * HP, amrex_real * gamma);
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_PRECOND(int * nJdata, int * consP);
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
#ifndef AMREX_USE_CUDA
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * indx, int * consP);
#else
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_PRECOND(int * rowPtr, int * colIndx, int * consP);
#endif
AMREX_GPU_HOST_DEVICE void aJacobian(amrex_real *  J, amrex_real *  sc, amrex_real T, int consP);
AMREX_GPU_HOST_DEVICE void aJacobian_precond(amrex_real *  J, amrex_real *  sc, amrex_real T, int HP);
AMREX_GPU_HOST_DEVICE void dcvpRdT(amrex_real *  species, amrex_real *  tc);
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(amrex_real *  e, amrex_real *  y, amrex_real *  t, int *ierr);
void GET_T_GIVEN_HY(amrex_real *  h, amrex_real *  y, amrex_real *  t, int *ierr);
void GET_CRITPARAMS(amrex_real *  Tci, amrex_real *  ai, amrex_real *  bi, amrex_real *  acentric_i);
/*vector version */
void VCKYTX(int *  np, amrex_real *  y, amrex_real *  x);
void VCKHMS(int *  np, amrex_real *  T, amrex_real *  ums);
void VCKWYR(int *  np, amrex_real *  rho, amrex_real *  T,
            amrex_real *  y,
            amrex_real *  wdot);
#ifndef AMREX_USE_CUDA
void vproductionRate(int npt, amrex_real *  wdot, amrex_real *  c, amrex_real *  T);
void VCKPY(int *  np, amrex_real *  rho, amrex_real *  T, amrex_real *  y, amrex_real *  P);
void vcomp_k_f(int npt, amrex_real *  k_f_s, amrex_real *  tc, amrex_real *  invT);
void vcomp_gibbs(int npt, amrex_real *  g_RT, amrex_real *  tc);
void vcomp_Kc(int npt, amrex_real *  Kc_s, amrex_real *  g_RT, amrex_real *  invT);
void vcomp_wdot_1_50(int npt, amrex_real *  wdot, amrex_real *  mixture, amrex_real *  sc,
                amrex_real *  k_f_s, amrex_real *  Kc_s,
                amrex_real *  tc, amrex_real *  invT, amrex_real *  T);
void vcomp_wdot_51_84(int npt, amrex_real *  wdot, amrex_real *  mixture, amrex_real *  sc,
                amrex_real *  k_f_s, amrex_real *  Kc_s,
                amrex_real *  tc, amrex_real *  invT, amrex_real *  T);
#endif
/*Transport function declarations */
void egtransetLENIMC(int* LENIMC);
void egtransetLENRMC(int* LENRMC);
void egtransetNO(int* NO);
void egtransetKK(int* KK);
void egtransetNLITE(int* NLITE);
void egtransetPATM(amrex_real* PATM);
void egtransetWT(amrex_real* WT);
void egtransetEPS(amrex_real* EPS);
void egtransetSIG(amrex_real* SIG);
void egtransetDIP(amrex_real* DIP);
void egtransetPOL(amrex_real* POL);
void egtransetZROT(amrex_real* ZROT);
void egtransetNLIN(int* NLIN);
void egtransetCOFETA(amrex_real* COFETA);
void egtransetCOFLAM(amrex_real* COFLAM);
void egtransetCOFD(amrex_real* COFD);
void egtransetKTDIF(int* KTDIF);
}

namespace thermo
{

    extern amrex_real fwd_A[84], fwd_beta[84], fwd_Ea[84];
    extern amrex_real low_A[84], low_beta[84], low_Ea[84];
    extern amrex_real rev_A[84], rev_beta[84], rev_Ea[84];
    extern amrex_real troe_a[84],troe_Ts[84], troe_Tss[84], troe_Tsss[84];
    extern amrex_real sri_a[84], sri_b[84], sri_c[84], sri_d[84], sri_e[84];
    extern amrex_real activation_units[84], prefactor_units[84], phase_units[84];
    extern int is_PD[84], troe_len[84], sri_len[84], nTB[84], *TBid[84];
    extern int NuVals[1764];
    extern amrex_real *TB[84];

    extern amrex_real fwd_A_DEF[84], fwd_beta_DEF[84], fwd_Ea_DEF[84];
    extern amrex_real low_A_DEF[84], low_beta_DEF[84], low_Ea_DEF[84];
    extern amrex_real rev_A_DEF[84], rev_beta_DEF[84], rev_Ea_DEF[84];
    extern amrex_real troe_a_DEF[84],troe_Ts_DEF[84], troe_Tss_DEF[84], troe_Tsss_DEF[84];
    extern amrex_real sri_a_DEF[84], sri_b_DEF[84], sri_c_DEF[84], sri_d_DEF[84], sri_e_DEF[84];
    extern amrex_real activation_units_DEF[84], prefactor_units_DEF[84], phase_units_DEF[84];
    extern int is_PD_DEF[84], troe_len_DEF[84], sri_len_DEF[84], nTB_DEF[84], *TBid_DEF[84];
    extern amrex_real *TB_DEF[84];

    extern std::vector<int> rxn_map;

}
