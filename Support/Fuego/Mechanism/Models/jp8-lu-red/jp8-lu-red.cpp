#include "chemistry_file.H"

#ifndef AMREX_USE_CUDA
namespace thermo
{
    double fwd_A[0], fwd_beta[0], fwd_Ea[0];
    double low_A[0], low_beta[0], low_Ea[0];
    double rev_A[0], rev_beta[0], rev_Ea[0];
    double troe_a[0],troe_Ts[0], troe_Tss[0], troe_Tsss[0];
    double sri_a[0], sri_b[0], sri_c[0], sri_d[0], sri_e[0];
    double activation_units[0], prefactor_units[0], phase_units[0];
    int is_PD[0], troe_len[0], sri_len[0], nTB[0], *TBid[0];
    double *TB[0];
    std::vector<std::vector<double>> kiv(0); 
    std::vector<std::vector<double>> nuv(0); 

    double fwd_A_DEF[0], fwd_beta_DEF[0], fwd_Ea_DEF[0];
    double low_A_DEF[0], low_beta_DEF[0], low_Ea_DEF[0];
    double rev_A_DEF[0], rev_beta_DEF[0], rev_Ea_DEF[0];
    double troe_a_DEF[0],troe_Ts_DEF[0], troe_Tss_DEF[0], troe_Tsss_DEF[0];
    double sri_a_DEF[0], sri_b_DEF[0], sri_c_DEF[0], sri_d_DEF[0], sri_e_DEF[0];
    double activation_units_DEF[0], prefactor_units_DEF[0], phase_units_DEF[0];
    int is_PD_DEF[0], troe_len_DEF[0], sri_len_DEF[0], nTB_DEF[0], *TBid_DEF[0];
    double *TB_DEF[0];
    std::vector<int> rxn_map;
};

using namespace thermo;
#endif

/* Inverse molecular weights */
/* TODO: check necessity on CPU */
static AMREX_GPU_DEVICE_MANAGED double imw[29] = {
    1.0 / 154.297990,  /*POSF10325 */
    1.0 / 28.054180,  /*C2H4 */
    1.0 / 16.043030,  /*CH4 */
    1.0 / 2.015940,  /*H2 */
    1.0 / 42.081270,  /*C3H6 */
    1.0 / 56.108360,  /*C4H81 */
    1.0 / 56.108360,  /*iC4H8 */
    1.0 / 78.114720,  /*C6H6 */
    1.0 / 92.141810,  /*C6H5CH3 */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 15.035060,  /*CH3 */
    1.0 / 30.026490,  /*CH2O */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 26.038240,  /*C2H2 */
    1.0 / 30.070120,  /*C2H6 */
    1.0 / 42.037640,  /*CH2CO */
    1.0 / 41.073300,  /*aC3H5 */
    1.0 / 77.106750,  /*C6H5 */
    1.0 / 91.133840,  /*C6H5CH2 */
    1.0 / 108.097580,  /*C6H4O2 */
    1.0 / 106.125270,  /*C6H5CHO */
    1.0 / 28.013400};  /*N2 */

/* Inverse molecular weights */
/* TODO: check necessity because redundant with molecularWeight */
static AMREX_GPU_DEVICE_MANAGED double molecular_weights[29] = {
    154.297990,  /*POSF10325 */
    28.054180,  /*C2H4 */
    16.043030,  /*CH4 */
    2.015940,  /*H2 */
    42.081270,  /*C3H6 */
    56.108360,  /*C4H81 */
    56.108360,  /*iC4H8 */
    78.114720,  /*C6H6 */
    92.141810,  /*C6H5CH3 */
    1.007970,  /*H */
    15.999400,  /*O */
    17.007370,  /*OH */
    33.006770,  /*HO2 */
    18.015340,  /*H2O */
    34.014740,  /*H2O2 */
    31.998800,  /*O2 */
    15.035060,  /*CH3 */
    30.026490,  /*CH2O */
    28.010550,  /*CO */
    44.009950,  /*CO2 */
    26.038240,  /*C2H2 */
    30.070120,  /*C2H6 */
    42.037640,  /*CH2CO */
    41.073300,  /*aC3H5 */
    77.106750,  /*C6H5 */
    91.133840,  /*C6H5CH2 */
    108.097580,  /*C6H4O2 */
    106.125270,  /*C6H5CHO */
    28.013400};  /*N2 */

AMREX_GPU_HOST_DEVICE
void get_imw(double imw_new[]){
    for(int i = 0; i<29; ++i) imw_new[i] = imw[i];
}

/* TODO: check necessity because redundant with CKWT */
AMREX_GPU_HOST_DEVICE
void get_mw(double mw_new[]){
    for(int i = 0; i<29; ++i) mw_new[i] = molecular_weights[i];
}


#ifndef AMREX_USE_CUDA
/* Initializes parameter database */
void CKINIT()
{

    rxn_map = {};

    SetAllDefaults();
}

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<0; ++i) {
        rmap[i] = rxn_map[i] + 1;
    }
}

#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  printf("No reactions in this model");
  abort();
  return 0;
}

void ResetAllParametersToDefault()
{
    for (int i=0; i<0; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }

        fwd_A[i]    = fwd_A_DEF[i];
        fwd_beta[i] = fwd_beta_DEF[i];
        fwd_Ea[i]   = fwd_Ea_DEF[i];

        low_A[i]    = low_A_DEF[i];
        low_beta[i] = low_beta_DEF[i];
        low_Ea[i]   = low_Ea_DEF[i];

        rev_A[i]    = rev_A_DEF[i];
        rev_beta[i] = rev_beta_DEF[i];
        rev_Ea[i]   = rev_Ea_DEF[i];

        troe_a[i]    = troe_a_DEF[i];
        troe_Ts[i]   = troe_Ts_DEF[i];
        troe_Tss[i]  = troe_Tss_DEF[i];
        troe_Tsss[i] = troe_Tsss_DEF[i];

        sri_a[i] = sri_a_DEF[i];
        sri_b[i] = sri_b_DEF[i];
        sri_c[i] = sri_c_DEF[i];
        sri_d[i] = sri_d_DEF[i];
        sri_e[i] = sri_e_DEF[i];

        is_PD[i]    = is_PD_DEF[i];
        troe_len[i] = troe_len_DEF[i];
        sri_len[i]  = sri_len_DEF[i];

        activation_units[i] = activation_units_DEF[i];
        prefactor_units[i]  = prefactor_units_DEF[i];
        phase_units[i]      = phase_units_DEF[i];

        nTB[i]  = nTB_DEF[i];
        if (nTB[i] != 0) {
           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);
           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);
           for (int j=0; j<nTB[i]; j++) {
             TB[i][j] = TB_DEF[i][j];
             TBid[i][j] = TBid_DEF[i][j];
           }
        }
    }
}

void SetAllDefaults()
{
    for (int i=0; i<0; i++) {
        if (nTB_DEF[i] != 0) {
            nTB_DEF[i] = 0;
            free(TB_DEF[i]);
            free(TBid_DEF[i]);
        }

        fwd_A_DEF[i]    = fwd_A[i];
        fwd_beta_DEF[i] = fwd_beta[i];
        fwd_Ea_DEF[i]   = fwd_Ea[i];

        low_A_DEF[i]    = low_A[i];
        low_beta_DEF[i] = low_beta[i];
        low_Ea_DEF[i]   = low_Ea[i];

        rev_A_DEF[i]    = rev_A[i];
        rev_beta_DEF[i] = rev_beta[i];
        rev_Ea_DEF[i]   = rev_Ea[i];

        troe_a_DEF[i]    = troe_a[i];
        troe_Ts_DEF[i]   = troe_Ts[i];
        troe_Tss_DEF[i]  = troe_Tss[i];
        troe_Tsss_DEF[i] = troe_Tsss[i];

        sri_a_DEF[i] = sri_a[i];
        sri_b_DEF[i] = sri_b[i];
        sri_c_DEF[i] = sri_c[i];
        sri_d_DEF[i] = sri_d[i];
        sri_e_DEF[i] = sri_e[i];

        is_PD_DEF[i]    = is_PD[i];
        troe_len_DEF[i] = troe_len[i];
        sri_len_DEF[i]  = sri_len[i];

        activation_units_DEF[i] = activation_units[i];
        prefactor_units_DEF[i]  = prefactor_units[i];
        phase_units_DEF[i]      = phase_units[i];

        nTB_DEF[i]  = nTB[i];
        if (nTB_DEF[i] != 0) {
           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);
           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);
           for (int j=0; j<nTB_DEF[i]; j++) {
             TB_DEF[i][j] = TB[i][j];
             TBid_DEF[i][j] = TBid[i][j];
           }
        }
    }
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<0; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]); TB_DEF[i] = 0; 
    free(TBid_DEF[i]); TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

#else
/* TODO: Remove on GPU, right now needed by chemistry_module on FORTRAN */
AMREX_GPU_HOST_DEVICE void CKINIT()
{
}

AMREX_GPU_HOST_DEVICE void CKFINALIZE()
{
}

#endif


/*A few mechanism parameters */
void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 5;
    *kk = 29;
    *ii = 0;
    *nfit = -1; /*Why do you need this anyway ?  */
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char cstr[1000];
    char *saveptr;
    char *p; /*String Tokens */
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            break;
        }
        cstr[i] = line[i];
    }
    cstr[i] = '\0';

    p = strtok_r(cstr," ", &saveptr);
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok_r(NULL, " ", &saveptr);
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*5; i++) {
        kname[i] = ' ';
    }

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 4*lenkname + 0 ] = 'A';
    kname[ 4*lenkname + 1 ] = 'R';
    kname[ 4*lenkname + 2 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*29; i++) {
        kname[i] = ' ';
    }

    /* POSF10325  */
    kname[ 0*lenkname + 0 ] = 'P';
    kname[ 0*lenkname + 1 ] = 'O';
    kname[ 0*lenkname + 2 ] = 'S';
    kname[ 0*lenkname + 3 ] = 'F';
    kname[ 0*lenkname + 4 ] = '1';
    kname[ 0*lenkname + 5 ] = '0';
    kname[ 0*lenkname + 6 ] = '3';
    kname[ 0*lenkname + 7 ] = '2';
    kname[ 0*lenkname + 8 ] = '5';
    kname[ 0*lenkname + 9 ] = ' ';

    /* C2H4  */
    kname[ 1*lenkname + 0 ] = 'C';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = 'H';
    kname[ 1*lenkname + 3 ] = '4';
    kname[ 1*lenkname + 4 ] = ' ';

    /* CH4  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = 'H';
    kname[ 2*lenkname + 2 ] = '4';
    kname[ 2*lenkname + 3 ] = ' ';

    /* H2  */
    kname[ 3*lenkname + 0 ] = 'H';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* C3H6  */
    kname[ 4*lenkname + 0 ] = 'C';
    kname[ 4*lenkname + 1 ] = '3';
    kname[ 4*lenkname + 2 ] = 'H';
    kname[ 4*lenkname + 3 ] = '6';
    kname[ 4*lenkname + 4 ] = ' ';

    /* C4H81  */
    kname[ 5*lenkname + 0 ] = 'C';
    kname[ 5*lenkname + 1 ] = '4';
    kname[ 5*lenkname + 2 ] = 'H';
    kname[ 5*lenkname + 3 ] = '8';
    kname[ 5*lenkname + 4 ] = '1';
    kname[ 5*lenkname + 5 ] = ' ';

    /* iC4H8  */
    kname[ 6*lenkname + 0 ] = 'I';
    kname[ 6*lenkname + 1 ] = 'C';
    kname[ 6*lenkname + 2 ] = '4';
    kname[ 6*lenkname + 3 ] = 'H';
    kname[ 6*lenkname + 4 ] = '8';
    kname[ 6*lenkname + 5 ] = ' ';

    /* C6H6  */
    kname[ 7*lenkname + 0 ] = 'C';
    kname[ 7*lenkname + 1 ] = '6';
    kname[ 7*lenkname + 2 ] = 'H';
    kname[ 7*lenkname + 3 ] = '6';
    kname[ 7*lenkname + 4 ] = ' ';

    /* C6H5CH3  */
    kname[ 8*lenkname + 0 ] = 'C';
    kname[ 8*lenkname + 1 ] = '6';
    kname[ 8*lenkname + 2 ] = 'H';
    kname[ 8*lenkname + 3 ] = '5';
    kname[ 8*lenkname + 4 ] = 'C';
    kname[ 8*lenkname + 5 ] = 'H';
    kname[ 8*lenkname + 6 ] = '3';
    kname[ 8*lenkname + 7 ] = ' ';

    /* H  */
    kname[ 9*lenkname + 0 ] = 'H';
    kname[ 9*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 10*lenkname + 0 ] = 'O';
    kname[ 10*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 11*lenkname + 0 ] = 'O';
    kname[ 11*lenkname + 1 ] = 'H';
    kname[ 11*lenkname + 2 ] = ' ';

    /* HO2  */
    kname[ 12*lenkname + 0 ] = 'H';
    kname[ 12*lenkname + 1 ] = 'O';
    kname[ 12*lenkname + 2 ] = '2';
    kname[ 12*lenkname + 3 ] = ' ';

    /* H2O  */
    kname[ 13*lenkname + 0 ] = 'H';
    kname[ 13*lenkname + 1 ] = '2';
    kname[ 13*lenkname + 2 ] = 'O';
    kname[ 13*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 14*lenkname + 0 ] = 'H';
    kname[ 14*lenkname + 1 ] = '2';
    kname[ 14*lenkname + 2 ] = 'O';
    kname[ 14*lenkname + 3 ] = '2';
    kname[ 14*lenkname + 4 ] = ' ';

    /* O2  */
    kname[ 15*lenkname + 0 ] = 'O';
    kname[ 15*lenkname + 1 ] = '2';
    kname[ 15*lenkname + 2 ] = ' ';

    /* CH3  */
    kname[ 16*lenkname + 0 ] = 'C';
    kname[ 16*lenkname + 1 ] = 'H';
    kname[ 16*lenkname + 2 ] = '3';
    kname[ 16*lenkname + 3 ] = ' ';

    /* CH2O  */
    kname[ 17*lenkname + 0 ] = 'C';
    kname[ 17*lenkname + 1 ] = 'H';
    kname[ 17*lenkname + 2 ] = '2';
    kname[ 17*lenkname + 3 ] = 'O';
    kname[ 17*lenkname + 4 ] = ' ';

    /* CO  */
    kname[ 18*lenkname + 0 ] = 'C';
    kname[ 18*lenkname + 1 ] = 'O';
    kname[ 18*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 19*lenkname + 0 ] = 'C';
    kname[ 19*lenkname + 1 ] = 'O';
    kname[ 19*lenkname + 2 ] = '2';
    kname[ 19*lenkname + 3 ] = ' ';

    /* C2H2  */
    kname[ 20*lenkname + 0 ] = 'C';
    kname[ 20*lenkname + 1 ] = '2';
    kname[ 20*lenkname + 2 ] = 'H';
    kname[ 20*lenkname + 3 ] = '2';
    kname[ 20*lenkname + 4 ] = ' ';

    /* C2H6  */
    kname[ 21*lenkname + 0 ] = 'C';
    kname[ 21*lenkname + 1 ] = '2';
    kname[ 21*lenkname + 2 ] = 'H';
    kname[ 21*lenkname + 3 ] = '6';
    kname[ 21*lenkname + 4 ] = ' ';

    /* CH2CO  */
    kname[ 22*lenkname + 0 ] = 'C';
    kname[ 22*lenkname + 1 ] = 'H';
    kname[ 22*lenkname + 2 ] = '2';
    kname[ 22*lenkname + 3 ] = 'C';
    kname[ 22*lenkname + 4 ] = 'O';
    kname[ 22*lenkname + 5 ] = ' ';

    /* aC3H5  */
    kname[ 23*lenkname + 0 ] = 'A';
    kname[ 23*lenkname + 1 ] = 'C';
    kname[ 23*lenkname + 2 ] = '3';
    kname[ 23*lenkname + 3 ] = 'H';
    kname[ 23*lenkname + 4 ] = '5';
    kname[ 23*lenkname + 5 ] = ' ';

    /* C6H5  */
    kname[ 24*lenkname + 0 ] = 'C';
    kname[ 24*lenkname + 1 ] = '6';
    kname[ 24*lenkname + 2 ] = 'H';
    kname[ 24*lenkname + 3 ] = '5';
    kname[ 24*lenkname + 4 ] = ' ';

    /* C6H5CH2  */
    kname[ 25*lenkname + 0 ] = 'C';
    kname[ 25*lenkname + 1 ] = '6';
    kname[ 25*lenkname + 2 ] = 'H';
    kname[ 25*lenkname + 3 ] = '5';
    kname[ 25*lenkname + 4 ] = 'C';
    kname[ 25*lenkname + 5 ] = 'H';
    kname[ 25*lenkname + 6 ] = '2';
    kname[ 25*lenkname + 7 ] = ' ';

    /* C6H4O2  */
    kname[ 26*lenkname + 0 ] = 'C';
    kname[ 26*lenkname + 1 ] = '6';
    kname[ 26*lenkname + 2 ] = 'H';
    kname[ 26*lenkname + 3 ] = '4';
    kname[ 26*lenkname + 4 ] = 'O';
    kname[ 26*lenkname + 5 ] = '2';
    kname[ 26*lenkname + 6 ] = ' ';

    /* C6H5CHO  */
    kname[ 27*lenkname + 0 ] = 'C';
    kname[ 27*lenkname + 1 ] = '6';
    kname[ 27*lenkname + 2 ] = 'H';
    kname[ 27*lenkname + 3 ] = '5';
    kname[ 27*lenkname + 4 ] = 'C';
    kname[ 27*lenkname + 5 ] = 'H';
    kname[ 27*lenkname + 6 ] = 'O';
    kname[ 27*lenkname + 7 ] = ' ';

    /* N2  */
    kname[ 28*lenkname + 0 ] = 'N';
    kname[ 28*lenkname + 1 ] = '2';
    kname[ 28*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(double *  ru, double *  ruc, double *  pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double *  rho, double *  T, double *  x, double *  P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*154.297990; /*POSF10325 */
    XW += x[1]*28.054180; /*C2H4 */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*42.081270; /*C3H6 */
    XW += x[5]*56.108360; /*C4H81 */
    XW += x[6]*56.108360; /*iC4H8 */
    XW += x[7]*78.114720; /*C6H6 */
    XW += x[8]*92.141810; /*C6H5CH3 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*17.007370; /*OH */
    XW += x[12]*33.006770; /*HO2 */
    XW += x[13]*18.015340; /*H2O */
    XW += x[14]*34.014740; /*H2O2 */
    XW += x[15]*31.998800; /*O2 */
    XW += x[16]*15.035060; /*CH3 */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*26.038240; /*C2H2 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*42.037640; /*CH2CO */
    XW += x[23]*41.073300; /*aC3H5 */
    XW += x[24]*77.106750; /*C6H5 */
    XW += x[25]*91.133840; /*C6H5CH2 */
    XW += x[26]*108.097580; /*C6H4O2 */
    XW += x[27]*106.125270; /*C6H5CHO */
    XW += x[28]*28.013400; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE void CKPY(double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*POSF10325 */
    YOW += y[1]*imw[1]; /*C2H4 */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*C3H6 */
    YOW += y[5]*imw[5]; /*C4H81 */
    YOW += y[6]*imw[6]; /*iC4H8 */
    YOW += y[7]*imw[7]; /*C6H6 */
    YOW += y[8]*imw[8]; /*C6H5CH3 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*OH */
    YOW += y[12]*imw[12]; /*HO2 */
    YOW += y[13]*imw[13]; /*H2O */
    YOW += y[14]*imw[14]; /*H2O2 */
    YOW += y[15]*imw[15]; /*O2 */
    YOW += y[16]*imw[16]; /*CH3 */
    YOW += y[17]*imw[17]; /*CH2O */
    YOW += y[18]*imw[18]; /*CO */
    YOW += y[19]*imw[19]; /*CO2 */
    YOW += y[20]*imw[20]; /*C2H2 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*CH2CO */
    YOW += y[23]*imw[23]; /*aC3H5 */
    YOW += y[24]*imw[24]; /*C6H5 */
    YOW += y[25]*imw[25]; /*C6H5CH2 */
    YOW += y[26]*imw[26]; /*C6H4O2 */
    YOW += y[27]*imw[27]; /*C6H5CHO */
    YOW += y[28]*imw[28]; /*N2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


#ifndef AMREX_USE_CUDA
/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}
#endif


/*Compute P = rhoRT/W(c) */
void CKPC(double *  rho, double *  T, double *  c,  double *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*154.297990; /*POSF10325 */
    W += c[1]*28.054180; /*C2H4 */
    W += c[2]*16.043030; /*CH4 */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*42.081270; /*C3H6 */
    W += c[5]*56.108360; /*C4H81 */
    W += c[6]*56.108360; /*iC4H8 */
    W += c[7]*78.114720; /*C6H6 */
    W += c[8]*92.141810; /*C6H5CH3 */
    W += c[9]*1.007970; /*H */
    W += c[10]*15.999400; /*O */
    W += c[11]*17.007370; /*OH */
    W += c[12]*33.006770; /*HO2 */
    W += c[13]*18.015340; /*H2O */
    W += c[14]*34.014740; /*H2O2 */
    W += c[15]*31.998800; /*O2 */
    W += c[16]*15.035060; /*CH3 */
    W += c[17]*30.026490; /*CH2O */
    W += c[18]*28.010550; /*CO */
    W += c[19]*44.009950; /*CO2 */
    W += c[20]*26.038240; /*C2H2 */
    W += c[21]*30.070120; /*C2H6 */
    W += c[22]*42.037640; /*CH2CO */
    W += c[23]*41.073300; /*aC3H5 */
    W += c[24]*77.106750; /*C6H5 */
    W += c[25]*91.133840; /*C6H5CH2 */
    W += c[26]*108.097580; /*C6H4O2 */
    W += c[27]*106.125270; /*C6H5CHO */
    W += c[28]*28.013400; /*N2 */

    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double *  P, double *  T, double *  x,  double *  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*154.297990; /*POSF10325 */
    XW += x[1]*28.054180; /*C2H4 */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*42.081270; /*C3H6 */
    XW += x[5]*56.108360; /*C4H81 */
    XW += x[6]*56.108360; /*iC4H8 */
    XW += x[7]*78.114720; /*C6H6 */
    XW += x[8]*92.141810; /*C6H5CH3 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*17.007370; /*OH */
    XW += x[12]*33.006770; /*HO2 */
    XW += x[13]*18.015340; /*H2O */
    XW += x[14]*34.014740; /*H2O2 */
    XW += x[15]*31.998800; /*O2 */
    XW += x[16]*15.035060; /*CH3 */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*26.038240; /*C2H2 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*42.037640; /*CH2CO */
    XW += x[23]*41.073300; /*aC3H5 */
    XW += x[24]*77.106750; /*C6H5 */
    XW += x[25]*91.133840; /*C6H5CH2 */
    XW += x[26]*108.097580; /*C6H4O2 */
    XW += x[27]*106.125270; /*C6H5CHO */
    XW += x[28]*28.013400; /*N2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double *  P, double *  T, double *  y,  double *  rho)
{
    double YOW = 0;
    double tmp[29];

    for (int i = 0; i < 29; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double *  P, double *  T, double *  c,  double *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*154.297990; /*POSF10325 */
    W += c[1]*28.054180; /*C2H4 */
    W += c[2]*16.043030; /*CH4 */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*42.081270; /*C3H6 */
    W += c[5]*56.108360; /*C4H81 */
    W += c[6]*56.108360; /*iC4H8 */
    W += c[7]*78.114720; /*C6H6 */
    W += c[8]*92.141810; /*C6H5CH3 */
    W += c[9]*1.007970; /*H */
    W += c[10]*15.999400; /*O */
    W += c[11]*17.007370; /*OH */
    W += c[12]*33.006770; /*HO2 */
    W += c[13]*18.015340; /*H2O */
    W += c[14]*34.014740; /*H2O2 */
    W += c[15]*31.998800; /*O2 */
    W += c[16]*15.035060; /*CH3 */
    W += c[17]*30.026490; /*CH2O */
    W += c[18]*28.010550; /*CO */
    W += c[19]*44.009950; /*CO2 */
    W += c[20]*26.038240; /*C2H2 */
    W += c[21]*30.070120; /*C2H6 */
    W += c[22]*42.037640; /*CH2CO */
    W += c[23]*41.073300; /*aC3H5 */
    W += c[24]*77.106750; /*C6H5 */
    W += c[25]*91.133840; /*C6H5CH2 */
    W += c[26]*108.097580; /*C6H4O2 */
    W += c[27]*106.125270; /*C6H5CHO */
    W += c[28]*28.013400; /*N2 */

    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT( double *  wt)
{
    get_mw(wt);
}


/*get atomic weight for all elements */
void CKAWT( double *  awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double *  y,  double *  wtm)
{
    double YOW = 0;
    double tmp[29];

    for (int i = 0; i < 29; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *  x,  double *  wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*154.297990; /*POSF10325 */
    XW += x[1]*28.054180; /*C2H4 */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*42.081270; /*C3H6 */
    XW += x[5]*56.108360; /*C4H81 */
    XW += x[6]*56.108360; /*iC4H8 */
    XW += x[7]*78.114720; /*C6H6 */
    XW += x[8]*92.141810; /*C6H5CH3 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*17.007370; /*OH */
    XW += x[12]*33.006770; /*HO2 */
    XW += x[13]*18.015340; /*H2O */
    XW += x[14]*34.014740; /*H2O2 */
    XW += x[15]*31.998800; /*O2 */
    XW += x[16]*15.035060; /*CH3 */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*26.038240; /*C2H2 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*42.037640; /*CH2CO */
    XW += x[23]*41.073300; /*aC3H5 */
    XW += x[24]*77.106750; /*C6H5 */
    XW += x[25]*91.133840; /*C6H5CH2 */
    XW += x[26]*108.097580; /*C6H4O2 */
    XW += x[27]*106.125270; /*C6H5CHO */
    XW += x[28]*28.013400; /*N2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double *  c,  double *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*154.297990; /*POSF10325 */
    W += c[1]*28.054180; /*C2H4 */
    W += c[2]*16.043030; /*CH4 */
    W += c[3]*2.015940; /*H2 */
    W += c[4]*42.081270; /*C3H6 */
    W += c[5]*56.108360; /*C4H81 */
    W += c[6]*56.108360; /*iC4H8 */
    W += c[7]*78.114720; /*C6H6 */
    W += c[8]*92.141810; /*C6H5CH3 */
    W += c[9]*1.007970; /*H */
    W += c[10]*15.999400; /*O */
    W += c[11]*17.007370; /*OH */
    W += c[12]*33.006770; /*HO2 */
    W += c[13]*18.015340; /*H2O */
    W += c[14]*34.014740; /*H2O2 */
    W += c[15]*31.998800; /*O2 */
    W += c[16]*15.035060; /*CH3 */
    W += c[17]*30.026490; /*CH2O */
    W += c[18]*28.010550; /*CO */
    W += c[19]*44.009950; /*CO2 */
    W += c[20]*26.038240; /*C2H2 */
    W += c[21]*30.070120; /*C2H6 */
    W += c[22]*42.037640; /*CH2CO */
    W += c[23]*41.073300; /*aC3H5 */
    W += c[24]*77.106750; /*C6H5 */
    W += c[25]*91.133840; /*C6H5CH2 */
    W += c[26]*108.097580; /*C6H4O2 */
    W += c[27]*106.125270; /*C6H5CHO */
    W += c[28]*28.013400; /*N2 */

    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE void CKYTX(double *  y,  double *  x)
{
    double YOW = 0;
    double tmp[29];

    for (int i = 0; i < 29; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 29; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


#ifndef AMREX_USE_CUDA
/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, double *  y,  double *  x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}
#else
/*TODO: remove this on GPU */
void VCKYTX(int *  np, double *  y,  double *  x)
{
}
#endif


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double *  P, double *  T, double *  y,  double *  c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 29; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 29; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 29; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE void CKYTCR(double *  rho, double *  T, double *  y,  double *  c)
{
    for (int i = 0; i < 29; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double *  x,  double *  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*154.297990; /*POSF10325 */
    XW += x[1]*28.054180; /*C2H4 */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*42.081270; /*C3H6 */
    XW += x[5]*56.108360; /*C4H81 */
    XW += x[6]*56.108360; /*iC4H8 */
    XW += x[7]*78.114720; /*C6H6 */
    XW += x[8]*92.141810; /*C6H5CH3 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*17.007370; /*OH */
    XW += x[12]*33.006770; /*HO2 */
    XW += x[13]*18.015340; /*H2O */
    XW += x[14]*34.014740; /*H2O2 */
    XW += x[15]*31.998800; /*O2 */
    XW += x[16]*15.035060; /*CH3 */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*26.038240; /*C2H2 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*42.037640; /*CH2CO */
    XW += x[23]*41.073300; /*aC3H5 */
    XW += x[24]*77.106750; /*C6H5 */
    XW += x[25]*91.133840; /*C6H5CH2 */
    XW += x[26]*108.097580; /*C6H4O2 */
    XW += x[27]*106.125270; /*C6H5CHO */
    XW += x[28]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*154.297990*XWinv; 
    y[1] = x[1]*28.054180*XWinv; 
    y[2] = x[2]*16.043030*XWinv; 
    y[3] = x[3]*2.015940*XWinv; 
    y[4] = x[4]*42.081270*XWinv; 
    y[5] = x[5]*56.108360*XWinv; 
    y[6] = x[6]*56.108360*XWinv; 
    y[7] = x[7]*78.114720*XWinv; 
    y[8] = x[8]*92.141810*XWinv; 
    y[9] = x[9]*1.007970*XWinv; 
    y[10] = x[10]*15.999400*XWinv; 
    y[11] = x[11]*17.007370*XWinv; 
    y[12] = x[12]*33.006770*XWinv; 
    y[13] = x[13]*18.015340*XWinv; 
    y[14] = x[14]*34.014740*XWinv; 
    y[15] = x[15]*31.998800*XWinv; 
    y[16] = x[16]*15.035060*XWinv; 
    y[17] = x[17]*30.026490*XWinv; 
    y[18] = x[18]*28.010550*XWinv; 
    y[19] = x[19]*44.009950*XWinv; 
    y[20] = x[20]*26.038240*XWinv; 
    y[21] = x[21]*30.070120*XWinv; 
    y[22] = x[22]*42.037640*XWinv; 
    y[23] = x[23]*41.073300*XWinv; 
    y[24] = x[24]*77.106750*XWinv; 
    y[25] = x[25]*91.133840*XWinv; 
    y[26] = x[26]*108.097580*XWinv; 
    y[27] = x[27]*106.125270*XWinv; 
    y[28] = x[28]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double *  P, double *  T, double *  x,  double *  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double *  rho, double *  T, double *  x, double *  c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*154.297990; /*POSF10325 */
    XW += x[1]*28.054180; /*C2H4 */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*42.081270; /*C3H6 */
    XW += x[5]*56.108360; /*C4H81 */
    XW += x[6]*56.108360; /*iC4H8 */
    XW += x[7]*78.114720; /*C6H6 */
    XW += x[8]*92.141810; /*C6H5CH3 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*17.007370; /*OH */
    XW += x[12]*33.006770; /*HO2 */
    XW += x[13]*18.015340; /*H2O */
    XW += x[14]*34.014740; /*H2O2 */
    XW += x[15]*31.998800; /*O2 */
    XW += x[16]*15.035060; /*CH3 */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*26.038240; /*C2H2 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*42.037640; /*CH2CO */
    XW += x[23]*41.073300; /*aC3H5 */
    XW += x[24]*77.106750; /*C6H5 */
    XW += x[25]*91.133840; /*C6H5CH2 */
    XW += x[26]*108.097580; /*C6H4O2 */
    XW += x[27]*106.125270; /*C6H5CHO */
    XW += x[28]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double *  c, double *  x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 29; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 29; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double *  c, double *  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*154.297990; /*POSF10325 */
    CW += c[1]*28.054180; /*C2H4 */
    CW += c[2]*16.043030; /*CH4 */
    CW += c[3]*2.015940; /*H2 */
    CW += c[4]*42.081270; /*C3H6 */
    CW += c[5]*56.108360; /*C4H81 */
    CW += c[6]*56.108360; /*iC4H8 */
    CW += c[7]*78.114720; /*C6H6 */
    CW += c[8]*92.141810; /*C6H5CH3 */
    CW += c[9]*1.007970; /*H */
    CW += c[10]*15.999400; /*O */
    CW += c[11]*17.007370; /*OH */
    CW += c[12]*33.006770; /*HO2 */
    CW += c[13]*18.015340; /*H2O */
    CW += c[14]*34.014740; /*H2O2 */
    CW += c[15]*31.998800; /*O2 */
    CW += c[16]*15.035060; /*CH3 */
    CW += c[17]*30.026490; /*CH2O */
    CW += c[18]*28.010550; /*CO */
    CW += c[19]*44.009950; /*CO2 */
    CW += c[20]*26.038240; /*C2H2 */
    CW += c[21]*30.070120; /*C2H6 */
    CW += c[22]*42.037640; /*CH2CO */
    CW += c[23]*41.073300; /*aC3H5 */
    CW += c[24]*77.106750; /*C6H5 */
    CW += c[25]*91.133840; /*C6H5CH2 */
    CW += c[26]*108.097580; /*C6H4O2 */
    CW += c[27]*106.125270; /*C6H5CHO */
    CW += c[28]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*154.297990*CWinv; 
    y[1] = c[1]*28.054180*CWinv; 
    y[2] = c[2]*16.043030*CWinv; 
    y[3] = c[3]*2.015940*CWinv; 
    y[4] = c[4]*42.081270*CWinv; 
    y[5] = c[5]*56.108360*CWinv; 
    y[6] = c[6]*56.108360*CWinv; 
    y[7] = c[7]*78.114720*CWinv; 
    y[8] = c[8]*92.141810*CWinv; 
    y[9] = c[9]*1.007970*CWinv; 
    y[10] = c[10]*15.999400*CWinv; 
    y[11] = c[11]*17.007370*CWinv; 
    y[12] = c[12]*33.006770*CWinv; 
    y[13] = c[13]*18.015340*CWinv; 
    y[14] = c[14]*34.014740*CWinv; 
    y[15] = c[15]*31.998800*CWinv; 
    y[16] = c[16]*15.035060*CWinv; 
    y[17] = c[17]*30.026490*CWinv; 
    y[18] = c[18]*28.010550*CWinv; 
    y[19] = c[19]*44.009950*CWinv; 
    y[20] = c[20]*26.038240*CWinv; 
    y[21] = c[21]*30.070120*CWinv; 
    y[22] = c[22]*42.037640*CWinv; 
    y[23] = c[23]*41.073300*CWinv; 
    y[24] = c[24]*77.106750*CWinv; 
    y[25] = c[25]*91.133840*CWinv; 
    y[26] = c[26]*108.097580*CWinv; 
    y[27] = c[27]*106.125270*CWinv; 
    y[28] = c[28]*28.013400*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double *  T, double *  cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double *  T, double *  hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double *  T, double *  sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double *  T,  double *  cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double *  T,  double *  cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *  T,  double *  uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double *  T,  double *  hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double *  T,  double *  gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double *  T,  double *  aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double *  T,  double *  sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_HOST_DEVICE void CKCVMS(double *  T,  double *  cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 5.388605515859278e+05; /*POSF10325 */
    cvms[1] *= 2.963733033722604e+06; /*C2H4 */
    cvms[2] *= 5.182630712527496e+06; /*CH4 */
    cvms[3] *= 4.124383662212169e+07; /*H2 */
    cvms[4] *= 1.975822022481736e+06; /*C3H6 */
    cvms[5] *= 1.481866516861302e+06; /*C4H81 */
    cvms[6] *= 1.481866516861302e+06; /*iC4H8 */
    cvms[7] *= 1.064397337659278e+06; /*C6H6 */
    cvms[8] *= 9.023601772094556e+05; /*C6H5CH3 */
    cvms[9] *= 8.248767324424338e+07; /*H */
    cvms[10] *= 5.196763628636074e+06; /*O */
    cvms[11] *= 4.888768810227566e+06; /*OH */
    cvms[12] *= 2.519031701678171e+06; /*HO2 */
    cvms[13] *= 4.615239012974499e+06; /*H2O */
    cvms[14] *= 2.444384405113783e+06; /*H2O2 */
    cvms[15] *= 2.598381814318037e+06; /*O2 */
    cvms[16] *= 5.530081023953346e+06; /*CH3 */
    cvms[17] *= 2.769058254894261e+06; /*CH2O */
    cvms[18] *= 2.968349425484326e+06; /*CO */
    cvms[19] *= 1.889234139098090e+06; /*CO2 */
    cvms[20] *= 3.193192012977835e+06; /*C2H2 */
    cvms[21] *= 2.765040511976673e+06; /*C2H6 */
    cvms[22] *= 1.977872687429646e+06; /*CH2CO */
    cvms[23] *= 2.024310196648431e+06; /*aC3H5 */
    cvms[24] *= 1.078311561568864e+06; /*C6H5 */
    cvms[25] *= 9.123405751365244e+05; /*C6H5CH2 */
    cvms[26] *= 7.691670803361185e+05; /*C6H4O2 */
    cvms[27] *= 7.834618465517213e+05; /*C6H5CHO */
    cvms[28] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE void CKCPMS(double *  T,  double *  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 5.388605515859278e+05; /*POSF10325 */
    cpms[1] *= 2.963733033722604e+06; /*C2H4 */
    cpms[2] *= 5.182630712527496e+06; /*CH4 */
    cpms[3] *= 4.124383662212169e+07; /*H2 */
    cpms[4] *= 1.975822022481736e+06; /*C3H6 */
    cpms[5] *= 1.481866516861302e+06; /*C4H81 */
    cpms[6] *= 1.481866516861302e+06; /*iC4H8 */
    cpms[7] *= 1.064397337659278e+06; /*C6H6 */
    cpms[8] *= 9.023601772094556e+05; /*C6H5CH3 */
    cpms[9] *= 8.248767324424338e+07; /*H */
    cpms[10] *= 5.196763628636074e+06; /*O */
    cpms[11] *= 4.888768810227566e+06; /*OH */
    cpms[12] *= 2.519031701678171e+06; /*HO2 */
    cpms[13] *= 4.615239012974499e+06; /*H2O */
    cpms[14] *= 2.444384405113783e+06; /*H2O2 */
    cpms[15] *= 2.598381814318037e+06; /*O2 */
    cpms[16] *= 5.530081023953346e+06; /*CH3 */
    cpms[17] *= 2.769058254894261e+06; /*CH2O */
    cpms[18] *= 2.968349425484326e+06; /*CO */
    cpms[19] *= 1.889234139098090e+06; /*CO2 */
    cpms[20] *= 3.193192012977835e+06; /*C2H2 */
    cpms[21] *= 2.765040511976673e+06; /*C2H6 */
    cpms[22] *= 1.977872687429646e+06; /*CH2CO */
    cpms[23] *= 2.024310196648431e+06; /*aC3H5 */
    cpms[24] *= 1.078311561568864e+06; /*C6H5 */
    cpms[25] *= 9.123405751365244e+05; /*C6H5CH2 */
    cpms[26] *= 7.691670803361185e+05; /*C6H4O2 */
    cpms[27] *= 7.834618465517213e+05; /*C6H5CHO */
    cpms[28] *= 2.968047434442088e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE void CKUMS(double *  T,  double *  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 29; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE void CKHMS(double *  T,  double *  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 29; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


#ifndef AMREX_USE_CUDA
/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[29];

    for (int i=0; i<(*np); i++) {
        tc[0] = 0.0;
        tc[1] = T[i];
        tc[2] = T[i]*T[i];
        tc[3] = T[i]*T[i]*T[i];
        tc[4] = T[i]*T[i]*T[i]*T[i];

        speciesEnthalpy(h, tc);

        hms[0*(*np)+i] = h[0];
        hms[1*(*np)+i] = h[1];
        hms[2*(*np)+i] = h[2];
        hms[3*(*np)+i] = h[3];
        hms[4*(*np)+i] = h[4];
        hms[5*(*np)+i] = h[5];
        hms[6*(*np)+i] = h[6];
        hms[7*(*np)+i] = h[7];
        hms[8*(*np)+i] = h[8];
        hms[9*(*np)+i] = h[9];
        hms[10*(*np)+i] = h[10];
        hms[11*(*np)+i] = h[11];
        hms[12*(*np)+i] = h[12];
        hms[13*(*np)+i] = h[13];
        hms[14*(*np)+i] = h[14];
        hms[15*(*np)+i] = h[15];
        hms[16*(*np)+i] = h[16];
        hms[17*(*np)+i] = h[17];
        hms[18*(*np)+i] = h[18];
        hms[19*(*np)+i] = h[19];
        hms[20*(*np)+i] = h[20];
        hms[21*(*np)+i] = h[21];
        hms[22*(*np)+i] = h[22];
        hms[23*(*np)+i] = h[23];
        hms[24*(*np)+i] = h[24];
        hms[25*(*np)+i] = h[25];
        hms[26*(*np)+i] = h[26];
        hms[27*(*np)+i] = h[27];
        hms[28*(*np)+i] = h[28];
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}
#else
/*TODO: remove this on GPU */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
}
#endif


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *  T,  double *  gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 29; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *  T,  double *  ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 29; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *  T,  double *  sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 5.388605515859278e+05; /*POSF10325 */
    sms[1] *= 2.963733033722604e+06; /*C2H4 */
    sms[2] *= 5.182630712527496e+06; /*CH4 */
    sms[3] *= 4.124383662212169e+07; /*H2 */
    sms[4] *= 1.975822022481736e+06; /*C3H6 */
    sms[5] *= 1.481866516861302e+06; /*C4H81 */
    sms[6] *= 1.481866516861302e+06; /*iC4H8 */
    sms[7] *= 1.064397337659278e+06; /*C6H6 */
    sms[8] *= 9.023601772094556e+05; /*C6H5CH3 */
    sms[9] *= 8.248767324424338e+07; /*H */
    sms[10] *= 5.196763628636074e+06; /*O */
    sms[11] *= 4.888768810227566e+06; /*OH */
    sms[12] *= 2.519031701678171e+06; /*HO2 */
    sms[13] *= 4.615239012974499e+06; /*H2O */
    sms[14] *= 2.444384405113783e+06; /*H2O2 */
    sms[15] *= 2.598381814318037e+06; /*O2 */
    sms[16] *= 5.530081023953346e+06; /*CH3 */
    sms[17] *= 2.769058254894261e+06; /*CH2O */
    sms[18] *= 2.968349425484326e+06; /*CO */
    sms[19] *= 1.889234139098090e+06; /*CO2 */
    sms[20] *= 3.193192012977835e+06; /*C2H2 */
    sms[21] *= 2.765040511976673e+06; /*C2H6 */
    sms[22] *= 1.977872687429646e+06; /*CH2CO */
    sms[23] *= 2.024310196648431e+06; /*aC3H5 */
    sms[24] *= 1.078311561568864e+06; /*C6H5 */
    sms[25] *= 9.123405751365244e+05; /*C6H5CH2 */
    sms[26] *= 7.691670803361185e+05; /*C6H4O2 */
    sms[27] *= 7.834618465517213e+05; /*C6H5CHO */
    sms[28] *= 2.968047434442088e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *  T, double *  x,  double *  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[29]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE void CKCPBS(double *  T, double *  y,  double *  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[29], tresult[29]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 29; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 29; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *  T, double *  x,  double *  cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[29]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE void CKCVBS(double *  T, double *  y,  double *  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[29]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*POSF10325 */
    result += cvor[1]*y[1]*imw[1]; /*C2H4 */
    result += cvor[2]*y[2]*imw[2]; /*CH4 */
    result += cvor[3]*y[3]*imw[3]; /*H2 */
    result += cvor[4]*y[4]*imw[4]; /*C3H6 */
    result += cvor[5]*y[5]*imw[5]; /*C4H81 */
    result += cvor[6]*y[6]*imw[6]; /*iC4H8 */
    result += cvor[7]*y[7]*imw[7]; /*C6H6 */
    result += cvor[8]*y[8]*imw[8]; /*C6H5CH3 */
    result += cvor[9]*y[9]*imw[9]; /*H */
    result += cvor[10]*y[10]*imw[10]; /*O */
    result += cvor[11]*y[11]*imw[11]; /*OH */
    result += cvor[12]*y[12]*imw[12]; /*HO2 */
    result += cvor[13]*y[13]*imw[13]; /*H2O */
    result += cvor[14]*y[14]*imw[14]; /*H2O2 */
    result += cvor[15]*y[15]*imw[15]; /*O2 */
    result += cvor[16]*y[16]*imw[16]; /*CH3 */
    result += cvor[17]*y[17]*imw[17]; /*CH2O */
    result += cvor[18]*y[18]*imw[18]; /*CO */
    result += cvor[19]*y[19]*imw[19]; /*CO2 */
    result += cvor[20]*y[20]*imw[20]; /*C2H2 */
    result += cvor[21]*y[21]*imw[21]; /*C2H6 */
    result += cvor[22]*y[22]*imw[22]; /*CH2CO */
    result += cvor[23]*y[23]*imw[23]; /*aC3H5 */
    result += cvor[24]*y[24]*imw[24]; /*C6H5 */
    result += cvor[25]*y[25]*imw[25]; /*C6H5CH2 */
    result += cvor[26]*y[26]*imw[26]; /*C6H4O2 */
    result += cvor[27]*y[27]*imw[27]; /*C6H5CHO */
    result += cvor[28]*y[28]*imw[28]; /*N2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *  T, double *  x,  double *  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[29]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
AMREX_GPU_HOST_DEVICE void CKHBMS(double *  T, double *  y,  double *  hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[29], tmp[29]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 29; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 29; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *  T, double *  x,  double *  ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[29]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 29; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
AMREX_GPU_HOST_DEVICE void CKUBMS(double *  T, double *  y,  double *  ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[29]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*POSF10325 */
    result += y[1]*ums[1]*imw[1]; /*C2H4 */
    result += y[2]*ums[2]*imw[2]; /*CH4 */
    result += y[3]*ums[3]*imw[3]; /*H2 */
    result += y[4]*ums[4]*imw[4]; /*C3H6 */
    result += y[5]*ums[5]*imw[5]; /*C4H81 */
    result += y[6]*ums[6]*imw[6]; /*iC4H8 */
    result += y[7]*ums[7]*imw[7]; /*C6H6 */
    result += y[8]*ums[8]*imw[8]; /*C6H5CH3 */
    result += y[9]*ums[9]*imw[9]; /*H */
    result += y[10]*ums[10]*imw[10]; /*O */
    result += y[11]*ums[11]*imw[11]; /*OH */
    result += y[12]*ums[12]*imw[12]; /*HO2 */
    result += y[13]*ums[13]*imw[13]; /*H2O */
    result += y[14]*ums[14]*imw[14]; /*H2O2 */
    result += y[15]*ums[15]*imw[15]; /*O2 */
    result += y[16]*ums[16]*imw[16]; /*CH3 */
    result += y[17]*ums[17]*imw[17]; /*CH2O */
    result += y[18]*ums[18]*imw[18]; /*CO */
    result += y[19]*ums[19]*imw[19]; /*CO2 */
    result += y[20]*ums[20]*imw[20]; /*C2H2 */
    result += y[21]*ums[21]*imw[21]; /*C2H6 */
    result += y[22]*ums[22]*imw[22]; /*CH2CO */
    result += y[23]*ums[23]*imw[23]; /*aC3H5 */
    result += y[24]*ums[24]*imw[24]; /*C6H5 */
    result += y[25]*ums[25]*imw[25]; /*C6H5CH2 */
    result += y[26]*ums[26]*imw[26]; /*C6H4O2 */
    result += y[27]*ums[27]*imw[27]; /*C6H5CHO */
    result += y[28]*ums[28]*imw[28]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double *  P, double *  T, double *  x,  double *  sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[29]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 29; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *  P, double *  T, double *  y,  double *  sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[29]; /* temporary storage */
    double x[29]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*POSF10325 */
    YOW += y[1]*imw[1]; /*C2H4 */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*C3H6 */
    YOW += y[5]*imw[5]; /*C4H81 */
    YOW += y[6]*imw[6]; /*iC4H8 */
    YOW += y[7]*imw[7]; /*C6H6 */
    YOW += y[8]*imw[8]; /*C6H5CH3 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*OH */
    YOW += y[12]*imw[12]; /*HO2 */
    YOW += y[13]*imw[13]; /*H2O */
    YOW += y[14]*imw[14]; /*H2O2 */
    YOW += y[15]*imw[15]; /*O2 */
    YOW += y[16]*imw[16]; /*CH3 */
    YOW += y[17]*imw[17]; /*CH2O */
    YOW += y[18]*imw[18]; /*CO */
    YOW += y[19]*imw[19]; /*CO2 */
    YOW += y[20]*imw[20]; /*C2H2 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*CH2CO */
    YOW += y[23]*imw[23]; /*aC3H5 */
    YOW += y[24]*imw[24]; /*C6H5 */
    YOW += y[25]*imw[25]; /*C6H5CH2 */
    YOW += y[26]*imw[26]; /*C6H4O2 */
    YOW += y[27]*imw[27]; /*C6H5CHO */
    YOW += y[28]*imw[28]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(154.297990*YOW); 
    x[1] = y[1]/(28.054180*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(42.081270*YOW); 
    x[5] = y[5]/(56.108360*YOW); 
    x[6] = y[6]/(56.108360*YOW); 
    x[7] = y[7]/(78.114720*YOW); 
    x[8] = y[8]/(92.141810*YOW); 
    x[9] = y[9]/(1.007970*YOW); 
    x[10] = y[10]/(15.999400*YOW); 
    x[11] = y[11]/(17.007370*YOW); 
    x[12] = y[12]/(33.006770*YOW); 
    x[13] = y[13]/(18.015340*YOW); 
    x[14] = y[14]/(34.014740*YOW); 
    x[15] = y[15]/(31.998800*YOW); 
    x[16] = y[16]/(15.035060*YOW); 
    x[17] = y[17]/(30.026490*YOW); 
    x[18] = y[18]/(28.010550*YOW); 
    x[19] = y[19]/(44.009950*YOW); 
    x[20] = y[20]/(26.038240*YOW); 
    x[21] = y[21]/(30.070120*YOW); 
    x[22] = y[22]/(42.037640*YOW); 
    x[23] = y[23]/(41.073300*YOW); 
    x[24] = y[24]/(77.106750*YOW); 
    x[25] = y[25]/(91.133840*YOW); 
    x[26] = y[26]/(108.097580*YOW); 
    x[27] = y[27]/(106.125270*YOW); 
    x[28] = y[28]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    result += x[9]*(sor[9]-log((x[9]+1e-100))-logPratio);
    result += x[10]*(sor[10]-log((x[10]+1e-100))-logPratio);
    result += x[11]*(sor[11]-log((x[11]+1e-100))-logPratio);
    result += x[12]*(sor[12]-log((x[12]+1e-100))-logPratio);
    result += x[13]*(sor[13]-log((x[13]+1e-100))-logPratio);
    result += x[14]*(sor[14]-log((x[14]+1e-100))-logPratio);
    result += x[15]*(sor[15]-log((x[15]+1e-100))-logPratio);
    result += x[16]*(sor[16]-log((x[16]+1e-100))-logPratio);
    result += x[17]*(sor[17]-log((x[17]+1e-100))-logPratio);
    result += x[18]*(sor[18]-log((x[18]+1e-100))-logPratio);
    result += x[19]*(sor[19]-log((x[19]+1e-100))-logPratio);
    result += x[20]*(sor[20]-log((x[20]+1e-100))-logPratio);
    result += x[21]*(sor[21]-log((x[21]+1e-100))-logPratio);
    result += x[22]*(sor[22]-log((x[22]+1e-100))-logPratio);
    result += x[23]*(sor[23]-log((x[23]+1e-100))-logPratio);
    result += x[24]*(sor[24]-log((x[24]+1e-100))-logPratio);
    result += x[25]*(sor[25]-log((x[25]+1e-100))-logPratio);
    result += x[26]*(sor[26]-log((x[26]+1e-100))-logPratio);
    result += x[27]*(sor[27]-log((x[27]+1e-100))-logPratio);
    result += x[28]*(sor[28]-log((x[28]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double *  P, double *  T, double *  x,  double *  gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[29]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 29; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double *  P, double *  T, double *  y,  double *  gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[29]; /* temporary storage */
    double x[29]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*POSF10325 */
    YOW += y[1]*imw[1]; /*C2H4 */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*C3H6 */
    YOW += y[5]*imw[5]; /*C4H81 */
    YOW += y[6]*imw[6]; /*iC4H8 */
    YOW += y[7]*imw[7]; /*C6H6 */
    YOW += y[8]*imw[8]; /*C6H5CH3 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*OH */
    YOW += y[12]*imw[12]; /*HO2 */
    YOW += y[13]*imw[13]; /*H2O */
    YOW += y[14]*imw[14]; /*H2O2 */
    YOW += y[15]*imw[15]; /*O2 */
    YOW += y[16]*imw[16]; /*CH3 */
    YOW += y[17]*imw[17]; /*CH2O */
    YOW += y[18]*imw[18]; /*CO */
    YOW += y[19]*imw[19]; /*CO2 */
    YOW += y[20]*imw[20]; /*C2H2 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*CH2CO */
    YOW += y[23]*imw[23]; /*aC3H5 */
    YOW += y[24]*imw[24]; /*C6H5 */
    YOW += y[25]*imw[25]; /*C6H5CH2 */
    YOW += y[26]*imw[26]; /*C6H4O2 */
    YOW += y[27]*imw[27]; /*C6H5CHO */
    YOW += y[28]*imw[28]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(154.297990*YOW); 
    x[1] = y[1]/(28.054180*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(42.081270*YOW); 
    x[5] = y[5]/(56.108360*YOW); 
    x[6] = y[6]/(56.108360*YOW); 
    x[7] = y[7]/(78.114720*YOW); 
    x[8] = y[8]/(92.141810*YOW); 
    x[9] = y[9]/(1.007970*YOW); 
    x[10] = y[10]/(15.999400*YOW); 
    x[11] = y[11]/(17.007370*YOW); 
    x[12] = y[12]/(33.006770*YOW); 
    x[13] = y[13]/(18.015340*YOW); 
    x[14] = y[14]/(34.014740*YOW); 
    x[15] = y[15]/(31.998800*YOW); 
    x[16] = y[16]/(15.035060*YOW); 
    x[17] = y[17]/(30.026490*YOW); 
    x[18] = y[18]/(28.010550*YOW); 
    x[19] = y[19]/(44.009950*YOW); 
    x[20] = y[20]/(26.038240*YOW); 
    x[21] = y[21]/(30.070120*YOW); 
    x[22] = y[22]/(42.037640*YOW); 
    x[23] = y[23]/(41.073300*YOW); 
    x[24] = y[24]/(77.106750*YOW); 
    x[25] = y[25]/(91.133840*YOW); 
    x[26] = y[26]/(108.097580*YOW); 
    x[27] = y[27]/(106.125270*YOW); 
    x[28] = y[28]/(28.013400*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(gort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(gort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(gort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(gort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(gort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(gort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(gort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(gort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(gort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(gort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(gort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(gort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(gort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(gort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(gort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(gort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(gort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(gort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(gort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(gort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(gort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(gort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(gort[28]+log((x[28]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double *  P, double *  T, double *  x,  double *  abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[29]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 29; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double *  P, double *  T, double *  y,  double *  abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[29]; /* temporary storage */
    double x[29]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*POSF10325 */
    YOW += y[1]*imw[1]; /*C2H4 */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*C3H6 */
    YOW += y[5]*imw[5]; /*C4H81 */
    YOW += y[6]*imw[6]; /*iC4H8 */
    YOW += y[7]*imw[7]; /*C6H6 */
    YOW += y[8]*imw[8]; /*C6H5CH3 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*OH */
    YOW += y[12]*imw[12]; /*HO2 */
    YOW += y[13]*imw[13]; /*H2O */
    YOW += y[14]*imw[14]; /*H2O2 */
    YOW += y[15]*imw[15]; /*O2 */
    YOW += y[16]*imw[16]; /*CH3 */
    YOW += y[17]*imw[17]; /*CH2O */
    YOW += y[18]*imw[18]; /*CO */
    YOW += y[19]*imw[19]; /*CO2 */
    YOW += y[20]*imw[20]; /*C2H2 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*CH2CO */
    YOW += y[23]*imw[23]; /*aC3H5 */
    YOW += y[24]*imw[24]; /*C6H5 */
    YOW += y[25]*imw[25]; /*C6H5CH2 */
    YOW += y[26]*imw[26]; /*C6H4O2 */
    YOW += y[27]*imw[27]; /*C6H5CHO */
    YOW += y[28]*imw[28]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(154.297990*YOW); 
    x[1] = y[1]/(28.054180*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(2.015940*YOW); 
    x[4] = y[4]/(42.081270*YOW); 
    x[5] = y[5]/(56.108360*YOW); 
    x[6] = y[6]/(56.108360*YOW); 
    x[7] = y[7]/(78.114720*YOW); 
    x[8] = y[8]/(92.141810*YOW); 
    x[9] = y[9]/(1.007970*YOW); 
    x[10] = y[10]/(15.999400*YOW); 
    x[11] = y[11]/(17.007370*YOW); 
    x[12] = y[12]/(33.006770*YOW); 
    x[13] = y[13]/(18.015340*YOW); 
    x[14] = y[14]/(34.014740*YOW); 
    x[15] = y[15]/(31.998800*YOW); 
    x[16] = y[16]/(15.035060*YOW); 
    x[17] = y[17]/(30.026490*YOW); 
    x[18] = y[18]/(28.010550*YOW); 
    x[19] = y[19]/(44.009950*YOW); 
    x[20] = y[20]/(26.038240*YOW); 
    x[21] = y[21]/(30.070120*YOW); 
    x[22] = y[22]/(42.037640*YOW); 
    x[23] = y[23]/(41.073300*YOW); 
    x[24] = y[24]/(77.106750*YOW); 
    x[25] = y[25]/(91.133840*YOW); 
    x[26] = y[26]/(108.097580*YOW); 
    x[27] = y[27]/(106.125270*YOW); 
    x[28] = y[28]/(28.013400*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(aort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(aort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(aort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(aort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(aort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(aort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(aort[12]+log((x[12]+1e-100))+logPratio);
    result += x[13]*(aort[13]+log((x[13]+1e-100))+logPratio);
    result += x[14]*(aort[14]+log((x[14]+1e-100))+logPratio);
    result += x[15]*(aort[15]+log((x[15]+1e-100))+logPratio);
    result += x[16]*(aort[16]+log((x[16]+1e-100))+logPratio);
    result += x[17]*(aort[17]+log((x[17]+1e-100))+logPratio);
    result += x[18]*(aort[18]+log((x[18]+1e-100))+logPratio);
    result += x[19]*(aort[19]+log((x[19]+1e-100))+logPratio);
    result += x[20]*(aort[20]+log((x[20]+1e-100))+logPratio);
    result += x[21]*(aort[21]+log((x[21]+1e-100))+logPratio);
    result += x[22]*(aort[22]+log((x[22]+1e-100))+logPratio);
    result += x[23]*(aort[23]+log((x[23]+1e-100))+logPratio);
    result += x[24]*(aort[24]+log((x[24]+1e-100))+logPratio);
    result += x[25]*(aort[25]+log((x[25]+1e-100))+logPratio);
    result += x[26]*(aort[26]+log((x[26]+1e-100))+logPratio);
    result += x[27]*(aort[27]+log((x[27]+1e-100))+logPratio);
    result += x[28]*(aort[28]+log((x[28]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE void CKWC(double *  T, double *  C,  double *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double *  P, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*POSF10325 */
    YOW += y[1]*imw[1]; /*C2H4 */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*C3H6 */
    YOW += y[5]*imw[5]; /*C4H81 */
    YOW += y[6]*imw[6]; /*iC4H8 */
    YOW += y[7]*imw[7]; /*C6H6 */
    YOW += y[8]*imw[8]; /*C6H5CH3 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*OH */
    YOW += y[12]*imw[12]; /*HO2 */
    YOW += y[13]*imw[13]; /*H2O */
    YOW += y[14]*imw[14]; /*H2O2 */
    YOW += y[15]*imw[15]; /*O2 */
    YOW += y[16]*imw[16]; /*CH3 */
    YOW += y[17]*imw[17]; /*CH2O */
    YOW += y[18]*imw[18]; /*CO */
    YOW += y[19]*imw[19]; /*CO2 */
    YOW += y[20]*imw[20]; /*C2H2 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*CH2CO */
    YOW += y[23]*imw[23]; /*aC3H5 */
    YOW += y[24]*imw[24]; /*C6H5 */
    YOW += y[25]*imw[25]; /*C6H5CH2 */
    YOW += y[26]*imw[26]; /*C6H4O2 */
    YOW += y[27]*imw[27]; /*C6H5CHO */
    YOW += y[28]*imw[28]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double *  P, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE void CKWYR(double *  rho, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, double *  rho, double *  T,
	    double *  y,
	    double *  wdot)
{
#ifndef AMREX_USE_CUDA
    double c[29*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<29*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
#endif
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double *  rho, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*154.297990; /*POSF10325 */
    XW += x[1]*28.054180; /*C2H4 */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*42.081270; /*C3H6 */
    XW += x[5]*56.108360; /*C4H81 */
    XW += x[6]*56.108360; /*iC4H8 */
    XW += x[7]*78.114720; /*C6H6 */
    XW += x[8]*92.141810; /*C6H5CH3 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*17.007370; /*OH */
    XW += x[12]*33.006770; /*HO2 */
    XW += x[13]*18.015340; /*H2O */
    XW += x[14]*34.014740; /*H2O2 */
    XW += x[15]*31.998800; /*O2 */
    XW += x[16]*15.035060; /*CH3 */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*26.038240; /*C2H2 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*42.037640; /*CH2CO */
    XW += x[23]*41.073300; /*aC3H5 */
    XW += x[24]*77.106750; /*C6H5 */
    XW += x[25]*91.133840; /*C6H5CH2 */
    XW += x[26]*108.097580; /*C6H4O2 */
    XW += x[27]*106.125270; /*C6H5CHO */
    XW += x[28]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*POSF10325 */
    YOW += y[1]*imw[1]; /*C2H4 */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*H2 */
    YOW += y[4]*imw[4]; /*C3H6 */
    YOW += y[5]*imw[5]; /*C4H81 */
    YOW += y[6]*imw[6]; /*iC4H8 */
    YOW += y[7]*imw[7]; /*C6H6 */
    YOW += y[8]*imw[8]; /*C6H5CH3 */
    YOW += y[9]*imw[9]; /*H */
    YOW += y[10]*imw[10]; /*O */
    YOW += y[11]*imw[11]; /*OH */
    YOW += y[12]*imw[12]; /*HO2 */
    YOW += y[13]*imw[13]; /*H2O */
    YOW += y[14]*imw[14]; /*H2O2 */
    YOW += y[15]*imw[15]; /*O2 */
    YOW += y[16]*imw[16]; /*CH3 */
    YOW += y[17]*imw[17]; /*CH2O */
    YOW += y[18]*imw[18]; /*CO */
    YOW += y[19]*imw[19]; /*CO2 */
    YOW += y[20]*imw[20]; /*C2H2 */
    YOW += y[21]*imw[21]; /*C2H6 */
    YOW += y[22]*imw[22]; /*CH2CO */
    YOW += y[23]*imw[23]; /*aC3H5 */
    YOW += y[24]*imw[24]; /*C6H5 */
    YOW += y[25]*imw[25]; /*C6H5CH2 */
    YOW += y[26]*imw[26]; /*C6H4O2 */
    YOW += y[27]*imw[27]; /*C6H5CHO */
    YOW += y[28]*imw[28]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 
    c[9] = PWORT * y[9]*imw[9]; 
    c[10] = PWORT * y[10]*imw[10]; 
    c[11] = PWORT * y[11]*imw[11]; 
    c[12] = PWORT * y[12]*imw[12]; 
    c[13] = PWORT * y[13]*imw[13]; 
    c[14] = PWORT * y[14]*imw[14]; 
    c[15] = PWORT * y[15]*imw[15]; 
    c[16] = PWORT * y[16]*imw[16]; 
    c[17] = PWORT * y[17]*imw[17]; 
    c[18] = PWORT * y[18]*imw[18]; 
    c[19] = PWORT * y[19]*imw[19]; 
    c[20] = PWORT * y[20]*imw[20]; 
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double *  P, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 
    c[9] = 1e6 * (*rho) * y[9]*imw[9]; 
    c[10] = 1e6 * (*rho) * y[10]*imw[10]; 
    c[11] = 1e6 * (*rho) * y[11]*imw[11]; 
    c[12] = 1e6 * (*rho) * y[12]*imw[12]; 
    c[13] = 1e6 * (*rho) * y[13]*imw[13]; 
    c[14] = 1e6 * (*rho) * y[14]*imw[14]; 
    c[15] = 1e6 * (*rho) * y[15]*imw[15]; 
    c[16] = 1e6 * (*rho) * y[16]*imw[16]; 
    c[17] = 1e6 * (*rho) * y[17]*imw[17]; 
    c[18] = 1e6 * (*rho) * y[18]*imw[18]; 
    c[19] = 1e6 * (*rho) * y[19]*imw[19]; 
    c[20] = 1e6 * (*rho) * y[20]*imw[20]; 
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[29]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*154.297990; /*POSF10325 */
    XW += x[1]*28.054180; /*C2H4 */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*2.015940; /*H2 */
    XW += x[4]*42.081270; /*C3H6 */
    XW += x[5]*56.108360; /*C4H81 */
    XW += x[6]*56.108360; /*iC4H8 */
    XW += x[7]*78.114720; /*C6H6 */
    XW += x[8]*92.141810; /*C6H5CH3 */
    XW += x[9]*1.007970; /*H */
    XW += x[10]*15.999400; /*O */
    XW += x[11]*17.007370; /*OH */
    XW += x[12]*33.006770; /*HO2 */
    XW += x[13]*18.015340; /*H2O */
    XW += x[14]*34.014740; /*H2O2 */
    XW += x[15]*31.998800; /*O2 */
    XW += x[16]*15.035060; /*CH3 */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*28.010550; /*CO */
    XW += x[19]*44.009950; /*CO2 */
    XW += x[20]*26.038240; /*C2H2 */
    XW += x[21]*30.070120; /*C2H6 */
    XW += x[22]*42.037640; /*CH2CO */
    XW += x[23]*41.073300; /*aC3H5 */
    XW += x[24]*77.106750; /*C6H5 */
    XW += x[25]*91.133840; /*C6H5CH2 */
    XW += x[26]*108.097580; /*C6H4O2 */
    XW += x[27]*106.125270; /*C6H5CHO */
    XW += x[28]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 0; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim,  int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 29 * kd; ++ id) {
         nuki[id] = 0; 
    }
}


#ifndef AMREX_USE_CUDA
/*Returns a count of species in a reaction, and their indices */
/*and stoichiometric coefficients. (Eq 50) */
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
    if (*i < 1) {
        /*Return max num species per reaction */
        *nspec = 0;
    } else {
        if (*i > 0) {
            *nspec = -1;
        } else {
            *nspec = kiv[*i-1].size();
            for (int j=0; j<*nspec; ++j) {
                ki[j] = kiv[*i-1][j] + 1;
                nu[j] = nuv[*i-1][j];
            }
        }
    }
}
#endif


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim,  int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 29; ++ id) {
         ncf[id] = 0; 
    }

    /*POSF10325 */
    ncf[ 0 * kd + 2 ] = 11; /*C */
    ncf[ 0 * kd + 1 ] = 22; /*H */

    /*C2H4 */
    ncf[ 1 * kd + 2 ] = 2; /*C */
    ncf[ 1 * kd + 1 ] = 4; /*H */

    /*CH4 */
    ncf[ 2 * kd + 2 ] = 1; /*C */
    ncf[ 2 * kd + 1 ] = 4; /*H */

    /*H2 */
    ncf[ 3 * kd + 1 ] = 2; /*H */

    /*C3H6 */
    ncf[ 4 * kd + 2 ] = 3; /*C */
    ncf[ 4 * kd + 1 ] = 6; /*H */

    /*C4H81 */
    ncf[ 5 * kd + 2 ] = 4; /*C */
    ncf[ 5 * kd + 1 ] = 8; /*H */

    /*iC4H8 */
    ncf[ 6 * kd + 1 ] = 8; /*H */
    ncf[ 6 * kd + 2 ] = 4; /*C */

    /*C6H6 */
    ncf[ 7 * kd + 2 ] = 6; /*C */
    ncf[ 7 * kd + 1 ] = 6; /*H */

    /*C6H5CH3 */
    ncf[ 8 * kd + 2 ] = 7; /*C */
    ncf[ 8 * kd + 1 ] = 8; /*H */

    /*H */
    ncf[ 9 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 10 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 11 * kd + 0 ] = 1; /*O */
    ncf[ 11 * kd + 1 ] = 1; /*H */

    /*HO2 */
    ncf[ 12 * kd + 1 ] = 1; /*H */
    ncf[ 12 * kd + 0 ] = 2; /*O */

    /*H2O */
    ncf[ 13 * kd + 1 ] = 2; /*H */
    ncf[ 13 * kd + 0 ] = 1; /*O */

    /*H2O2 */
    ncf[ 14 * kd + 1 ] = 2; /*H */
    ncf[ 14 * kd + 0 ] = 2; /*O */

    /*O2 */
    ncf[ 15 * kd + 0 ] = 2; /*O */

    /*CH3 */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 1 ] = 3; /*H */

    /*CH2O */
    ncf[ 17 * kd + 1 ] = 2; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CO */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 0 ] = 2; /*O */

    /*C2H2 */
    ncf[ 20 * kd + 2 ] = 2; /*C */
    ncf[ 20 * kd + 1 ] = 2; /*H */

    /*C2H6 */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 6; /*H */

    /*CH2CO */
    ncf[ 22 * kd + 2 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */
    ncf[ 22 * kd + 0 ] = 1; /*O */

    /*aC3H5 */
    ncf[ 23 * kd + 2 ] = 3; /*C */
    ncf[ 23 * kd + 1 ] = 5; /*H */

    /*C6H5 */
    ncf[ 24 * kd + 2 ] = 6; /*C */
    ncf[ 24 * kd + 1 ] = 5; /*H */

    /*C6H5CH2 */
    ncf[ 25 * kd + 2 ] = 7; /*C */
    ncf[ 25 * kd + 1 ] = 7; /*H */

    /*C6H4O2 */
    ncf[ 26 * kd + 2 ] = 6; /*C */
    ncf[ 26 * kd + 1 ] = 4; /*H */
    ncf[ 26 * kd + 0 ] = 2; /*O */

    /*C6H5CHO */
    ncf[ 27 * kd + 2 ] = 7; /*C */
    ncf[ 27 * kd + 1 ] = 6; /*H */
    ncf[ 27 * kd + 0 ] = 1; /*O */

    /*N2 */
    ncf[ 28 * kd + 3 ] = 2; /*N */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE( double *  a, double *  b, double *  e)
{

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double *  T, double *  C, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[29]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);
}

#ifdef AMREX_USE_CUDA
/*GPU version of productionRate: no more use of thermo namespace vectors */
/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE inline void  productionRate(double * wdot, double * sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];


    for (int i = 0; i < 29; ++i) {
        wdot[i] = 0.0;
    }

    return;
}

AMREX_GPU_HOST_DEVICE inline void comp_qfqr(double *  qf, double * qr, double * sc, double * tc, double invT)
{

    return;
}
#endif


#ifndef AMREX_USE_CUDA
static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[0];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[0];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif




/*compute the production rate for each species pointwise on CPU */
void productionRate(double *  wdot, double *  sc, double T)
{
  productionrate_(wdot, sc, &T);
}

void comp_k_f(double *  tc, double invT, double *  k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<0; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double *  tc, double invT, double *  Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[29];
    gibbs(g_RT, tc);


#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<0; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;


    return;
}

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
{

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 29; ++i) {
        mixture += sc[i];
    }

    double Corr[0];
    for (int i = 0; i < 0; ++i) {
        Corr[i] = 1.0;
    }

    for (int i=0; i<0; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}
#endif


#ifndef AMREX_USE_CUDA
/*compute the production rate for each species */
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[0*npt], Kc_s[0*npt], mixture[npt], g_RT[29*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
}

void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
    }
}

void vcomp_gibbs(int npt, double *  g_RT, double *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[29];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
        g_RT[7*npt+i] = g[7];
        g_RT[8*npt+i] = g[8];
        g_RT[9*npt+i] = g[9];
        g_RT[10*npt+i] = g[10];
        g_RT[11*npt+i] = g[11];
        g_RT[12*npt+i] = g[12];
        g_RT[13*npt+i] = g[13];
        g_RT[14*npt+i] = g[14];
        g_RT[15*npt+i] = g[15];
        g_RT[16*npt+i] = g[16];
        g_RT[17*npt+i] = g[17];
        g_RT[18*npt+i] = g[18];
        g_RT[19*npt+i] = g[19];
        g_RT[20*npt+i] = g[20];
        g_RT[21*npt+i] = g[21];
        g_RT[22*npt+i] = g[22];
        g_RT[23*npt+i] = g[23];
        g_RT[24*npt+i] = g[24];
        g_RT[25*npt+i] = g[25];
        g_RT[26*npt+i] = g[26];
        g_RT[27*npt+i] = g[27];
        g_RT[28*npt+i] = g[28];
    }
}

void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

    }
}

void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
		double *  k_f_s, double *  Kc_s,
		double *  tc, double *  invT, double *  T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;
    }
}
#endif

/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void DWDOT_PRECOND(double *  J, double *  sc, double *  Tp, int * HP)
{
    double c[29];

    for (int k=0; k<29; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<29; k++) {
        J[870+k] *= 1.e-6;
        J[k*30+29] *= 1.e6;
    }

    return;
}

/*compute an approx to the SPS Jacobian */
AMREX_GPU_HOST_DEVICE void SLJ_PRECOND_CSC(double *  Jsps, int * indx, int * len, double * sc, double * Tp, int * HP, double * gamma)
{
    double c[29];
    double J[900];
    double mwt[29];

    get_mw(mwt);

    for (int k=0; k<29; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* Change of coord */
    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<29; k++) {
        J[870+k] = 1.e-6 * J[870+k] * mwt[k];
        J[k*30+29] = 1.e6 * J[k*30+29] / mwt[k];
    }
    /* dTdot/dT */
    /* dwdot[l]/[k] */
    for (int k=0; k<29; k++) {
        for (int l=0; l<29; l++) {
            /* DIAG elem */
            if (k == l){
                J[ 30 * k + l] =  J[ 30 * k + l] * mwt[l] / mwt[k];
            /* NOT DIAG and not last column nor last row */
            } else {
                J[ 30 * k + l] =  J[ 30 * k + l] * mwt[l] / mwt[k];
            }
        }
    }

    for (int k=0; k<(*len); k++) {
        Jsps[k] = J[indx[k]];
    }

    return;
}

/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)
{
    double c[29];

    for (int k=0; k<29; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<29; k++) {
        J[870+k] *= 1.e-6;
        J[k*30+29] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    double c[29];
    double J[900];

    for (int k=0; k<29; k++) {
        c[k] = 1.0/ 29.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(J[ 30 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of simplified Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_PRECOND( int * nJdata, int * consP)
{
    double c[29];
    double J[900];

    for (int k=0; k<29; k++) {
        c[k] = 1.0/ 29.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 30 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}


#ifndef AMREX_USE_CUDA
/*compute the sparsity pattern of the simplified precond Jacobian on CPU */
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    double c[29];
    double J[900];

    for (int k=0; k<29; k++) {
        c[k] = 1.0/ 29.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 30*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[30*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 30*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}
#else

/*compute the sparsity pattern of the simplified precond Jacobian on GPU */
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_PRECOND(int * rowPtr, int * colIndx, int * consP)
{
    double c[29];
    double J[900];

    for (int k=0; k<29; k++) {
        c[k] = 1.0/ 29.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l=0; l<30; l++) {
        for (int k=0; k<30; k++) {
            if (k == l) {
                colIndx[nJdata_tmp-1] = l+1; 
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[30*k + l] != 0.0) {
                    colIndx[nJdata_tmp-1] = k+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        rowPtr[l+1] = nJdata_tmp;
    }

    return;
}
#endif

/*compute the sparsity pattern of the Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    double c[29];
    double J[900];
    int offset_row;
    int offset_col;

    for (int k=0; k<29; k++) {
        c[k] = 1.0/ 29.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 30;
        offset_col = nc * 30;
        for (int k=0; k<30; k++) {
            for (int l=0; l<30; l++) {
                if(J[30*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}


#ifdef AMREX_USE_CUDA
/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE
void aJacobian(double * J, double * sc, double T, int consP)
{


    for (int i=0; i<900; i++) {
        J[i] = 0.0;
    }

    double wdot[29];
    for (int k=0; k<29; k++) {
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
    for (int k = 0; k < 29; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[29];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[29];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[29];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[29], dcRdT[29], e_RT[29];
    double * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 29; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[870+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 29; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 29; ++m) {
            dehmixdc += eh_RT[m]*J[k*30+m];
        }
        J[k*30+29] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[899] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


#ifndef AMREX_USE_CUDA
/*compute the reaction Jacobian on CPU */
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<900; i++) {
        J[i] = 0.0;
    }

    double wdot[29];
    for (int k=0; k<29; k++) {
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
    for (int k = 0; k < 29; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[29];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[29];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[29];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[29], dcRdT[29], e_RT[29];
    double * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 29; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[870+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 29; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 29; ++m) {
            dehmixdc += eh_RT[m]*J[k*30+m];
        }
        J[k*30+29] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[899] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}
#endif


/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP)
{
    for (int i=0; i<900; i++) {
        J[i] = 0.0;
    }

    double wdot[29];
    for (int k=0; k<29; k++) {
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
    for (int k = 0; k < 29; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[29];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[29];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[29];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    double c_R[29], dcRdT[29], e_RT[29];
    double * eh_RT;
    if (HP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 29; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[870+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 29; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 29; ++m) {
            dehmixdc += eh_RT[m]*J[k*30+m];
        }
        J[k*30+29] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[899] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void dcvpRdT(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            +7.12034480e-02
            +1.73991910e-04 * tc[1]
            -5.05836180e-07 * tc[2]
            +3.03945480e-10 * tc[3];
        /*species 1: C2H4 */
        species[1] =
            -7.57052247e-03
            +1.14198058e-04 * tc[1]
            -2.07476626e-07 * tc[2]
            +1.07953749e-10 * tc[3];
        /*species 2: CH4 */
        species[2] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 3: H2 */
        species[3] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
        /*species 4: C3H6 */
        species[4] =
            +2.09251800e-02
            +8.97358800e-06 * tc[1]
            -5.00673600e-08 * tc[2]
            +2.86325840e-11 * tc[3];
        /*species 5: C4H81 */
        species[5] =
            +3.08533800e-02
            +1.01730494e-05 * tc[1]
            -7.39646640e-08 * tc[2]
            +4.44407720e-11 * tc[3];
        /*species 6: iC4H8 */
        species[6] =
            +2.59029570e-02
            +1.63970708e-05 * tc[1]
            -6.65797770e-08 * tc[2]
            +3.55834320e-11 * tc[3];
        /*species 7: C6H6 */
        species[7] =
            +5.84276130e-02
            -5.89717100e-05 * tc[1]
            -2.08171320e-08 * tc[2]
            +3.28501012e-11 * tc[3];
        /*species 8: C6H5CH3 */
        species[8] =
            +2.10994380e-02
            +1.70732036e-04 * tc[1]
            -3.97831980e-07 * tc[2]
            +2.23826416e-10 * tc[3];
        /*species 9: H */
        species[9] =
            +7.05332819e-13
            -3.99183928e-15 * tc[1]
            +6.90244896e-18 * tc[2]
            -3.71092933e-21 * tc[3];
        /*species 10: O */
        species[10] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 11: OH */
        species[11] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 12: HO2 */
        species[12] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 13: H2O */
        species[13] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 14: H2O2 */
        species[14] =
            -5.42822417e-04
            +3.34671402e-05 * tc[1]
            -6.47312439e-08 * tc[2]
            +3.44981745e-11 * tc[3];
        /*species 15: O2 */
        species[15] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 16: CH3 */
        species[16] =
            +2.01095175e-03
            +1.14604371e-05 * tc[1]
            -2.06135228e-08 * tc[2]
            +1.01754294e-11 * tc[3];
        /*species 17: CH2O */
        species[17] =
            -9.90833369e-03
            +7.46440016e-05 * tc[1]
            -1.13785578e-07 * tc[2]
            +5.27090608e-11 * tc[3];
        /*species 18: CO */
        species[18] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 19: CO2 */
        species[19] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 20: C2H2 */
        species[20] =
            +2.33615629e-02
            -7.10343630e-05 * tc[1]
            +8.40457311e-08 * tc[2]
            -3.40029190e-11 * tc[3];
        /*species 21: C2H6 */
        species[21] =
            -5.50154270e-03
            +1.19887658e-04 * tc[1]
            -2.12539886e-07 * tc[2]
            +1.07474308e-10 * tc[3];
        /*species 22: CH2CO */
        species[22] =
            +1.81188721e-02
            -3.47894948e-05 * tc[1]
            +2.80319270e-08 * tc[2]
            -8.05830460e-12 * tc[3];
        /*species 23: aC3H5 */
        species[23] =
            +1.98138210e-02
            +2.49941200e-05 * tc[1]
            -1.00066665e-07 * tc[2]
            +6.33862840e-11 * tc[3];
        /*species 24: C6H5 */
        species[24] =
            +5.21789680e-02
            -5.11168540e-05 * tc[1]
            -2.11983363e-08 * tc[2]
            +3.03335900e-11 * tc[3];
        /*species 25: C6H5CH2 */
        species[25] =
            +3.85128320e-02
            +6.57229840e-05 * tc[1]
            -2.30918163e-07 * tc[2]
            +1.41692272e-10 * tc[3];
        /*species 26: C6H4O2 */
        species[26] =
            +5.78424450e-02
            -7.64288780e-05 * tc[1]
            +1.38937968e-08 * tc[2]
            +1.45186604e-11 * tc[3];
        /*species 27: C6H5CHO */
        species[27] =
            +6.63692450e-02
            -6.96327060e-05 * tc[1]
            -1.88998131e-08 * tc[2]
            +3.43228404e-11 * tc[3];
        /*species 28: N2 */
        species[28] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
    } else {
        /*species 0: POSF10325 */
        species[0] =
            +6.61883880e-02
            -4.83170000e-05 * tc[1]
            +1.07596491e-08 * tc[2]
            -4.09075920e-13 * tc[3];
        /*species 1: C2H4 */
        species[1] =
            +1.46454151e-02
            -1.34215583e-05 * tc[1]
            +4.41668769e-09 * tc[2]
            -5.02824244e-13 * tc[3];
        /*species 2: CH4 */
        species[2] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 3: H2 */
        species[3] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
        /*species 4: C3H6 */
        species[4] =
            +1.49083400e-02
            -9.89979800e-06 * tc[1]
            +2.16360660e-09 * tc[2]
            -1.50648160e-13 * tc[3];
        /*species 5: C4H81 */
        species[5] =
            +3.43505070e-02
            -3.17663940e-05 * tc[1]
            +9.92689860e-09 * tc[2]
            -1.01444180e-12 * tc[3];
        /*species 6: iC4H8 */
        species[6] =
            +2.96114870e-02
            -2.61542580e-05 * tc[1]
            +7.97158020e-09 * tc[2]
            -8.05388520e-13 * tc[3];
        /*species 7: C6H6 */
        species[7] =
            +2.38544330e-02
            -1.76255452e-05 * tc[1]
            +3.62970630e-09 * tc[2]
            -7.28860120e-14 * tc[3];
        /*species 8: C6H5CH3 */
        species[8] =
            +2.66912870e-02
            -1.93677010e-05 * tc[1]
            +4.72158870e-09 * tc[2]
            -3.78654404e-13 * tc[3];
        /*species 9: H */
        species[9] =
            -2.30842973e-11
            +3.23123896e-14 * tc[1]
            -1.42054571e-17 * tc[2]
            +1.99278943e-21 * tc[3];
        /*species 10: O */
        species[10] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 11: OH */
        species[11] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 12: HO2 */
        species[12] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 13: H2O */
        species[13] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 14: H2O2 */
        species[14] =
            +4.90831694e-03
            -3.80278450e-06 * tc[1]
            +1.11355796e-09 * tc[2]
            -1.15163322e-13 * tc[3];
        /*species 15: O2 */
        species[15] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 16: CH3 */
        species[16] =
            +7.23990037e-03
            -5.97428696e-06 * tc[1]
            +1.78705393e-09 * tc[2]
            -1.86861758e-13 * tc[3];
        /*species 17: CH2O */
        species[17] =
            +9.20000082e-03
            -8.84517626e-06 * tc[1]
            +3.01923636e-09 * tc[2]
            -3.53542256e-13 * tc[3];
        /*species 18: CO */
        species[18] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 19: CO2 */
        species[19] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 20: C2H2 */
        species[20] =
            +5.96166664e-03
            -4.74589704e-06 * tc[1]
            +1.40223651e-09 * tc[2]
            -1.44494085e-13 * tc[3];
        /*species 21: C2H6 */
        species[21] =
            +2.16852677e-02
            -2.00512134e-05 * tc[1]
            +6.64236003e-09 * tc[2]
            -7.60011560e-13 * tc[3];
        /*species 22: CH2CO */
        species[22] =
            +9.00359745e-03
            -8.33879270e-06 * tc[1]
            +2.77003765e-09 * tc[2]
            -3.17935280e-13 * tc[3];
        /*species 23: aC3H5 */
        species[23] =
            +1.43247310e-02
            -1.13563264e-05 * tc[1]
            +3.32424030e-09 * tc[2]
            -3.61455548e-13 * tc[3];
        /*species 24: C6H5 */
        species[24] =
            +2.22416300e-02
            -1.74399956e-05 * tc[1]
            +4.13663550e-09 * tc[2]
            -2.12584224e-13 * tc[3];
        /*species 25: C6H5CH2 */
        species[25] =
            +2.34938730e-02
            -1.70750734e-05 * tc[1]
            +4.16725230e-09 * tc[2]
            -3.34457680e-13 * tc[3];
        /*species 26: C6H4O2 */
        species[26] =
            +2.36149950e-02
            -2.04691520e-05 * tc[1]
            +5.85965220e-09 * tc[2]
            -5.09840880e-13 * tc[3];
        /*species 27: C6H5CHO */
        species[27] =
            +2.56804190e-02
            -2.09334580e-05 * tc[1]
            +5.82402900e-09 * tc[2]
            -5.39351680e-13 * tc[3];
        /*species 28: N2 */
        species[28] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
AMREX_GPU_HOST_DEVICE void progressRate(double *  qdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

#ifndef AMREX_USE_CUDA
    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }
#endif


    return;
}


/*compute the progress rate for each reaction */
AMREX_GPU_HOST_DEVICE void progressRateFR(double *  q_f, double *  q_r, double *  sc, double T)
{
    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *  kc, double *  g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void gibbs(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            -4.044634400000000e+04 * invT
            -1.736510830000000e+01
            -3.168577700000000e+00 * tc[0]
            -3.560172400000000e-02 * tc[1]
            -1.449932583333333e-05 * tc[2]
            +1.405100500000000e-08 * tc[3]
            -3.799318500000000e-12 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +5.089775930000000e+03 * invT
            -1.381294799999999e-01
            -3.959201480000000e+00 * tc[0]
            +3.785261235000000e-03 * tc[1]
            -9.516504866666667e-06 * tc[2]
            +5.763239608333333e-09 * tc[3]
            -1.349421865000000e-12 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            +1.074826000000000e+03 * invT
            -1.465203300000000e+01
            -1.493307000000000e+00 * tc[0]
            -1.046259000000000e-02 * tc[1]
            -7.477990000000000e-07 * tc[2]
            +1.390760000000000e-09 * tc[3]
            -3.579073000000000e-13 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            -1.790400400000000e+03 * invT
            -1.988133100000000e+01
            -1.181138000000000e+00 * tc[0]
            -1.542669000000000e-02 * tc[1]
            -8.477541166666667e-07 * tc[2]
            +2.054574000000000e-09 * tc[3]
            -5.555096499999999e-13 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            -4.037306900000000e+03 * invT
            -1.002924750000000e+01
            -2.647140500000000e+00 * tc[0]
            -1.295147850000000e-02 * tc[1]
            -1.366422566666667e-06 * tc[2]
            +1.849438250000000e-09 * tc[3]
            -4.447929000000000e-13 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            +9.181777300000000e+03 * invT
            -4.873360540000000e+01
            +4.843773400000000e+00 * tc[0]
            -2.921380650000000e-02 * tc[1]
            +4.914309166666667e-06 * tc[2]
            +5.782536666666667e-10 * tc[3]
            -4.106262650000000e-13 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            +4.075630000000000e+03 * invT
            -1.866694370000000e+01
            -1.615266300000000e+00 * tc[0]
            -1.054971900000000e-02 * tc[1]
            -1.422766966666667e-05 * tc[2]
            +1.105088833333333e-08 * tc[3]
            -2.797830200000000e-12 * tc[4];
        /*species 9: H */
        species[9] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
        /*species 10: O */
        species[10] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 11: OH */
        species[11] =
            +3.381538120000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 13: H2O */
        species[13] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 15: O2 */
        species[15] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +1.644499880000000e+04 * invT
            +2.069026070000000e+00
            -3.673590400000000e+00 * tc[0]
            -1.005475875000000e-03 * tc[1]
            -9.550364266666668e-07 * tc[2]
            +5.725978541666666e-10 * tc[3]
            -1.271928670000000e-13 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.430895670000000e+04 * invT
            +4.190910250000000e+00
            -4.793723150000000e+00 * tc[0]
            +4.954166845000000e-03 * tc[1]
            -6.220333466666666e-06 * tc[2]
            +3.160710508333333e-09 * tc[3]
            -6.588632600000000e-13 * tc[4];
        /*species 18: CO */
        species[18] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 19: CO2 */
        species[19] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            +2.642898070000000e+04 * invT
            -1.313102400600000e+01
            -8.086810940000000e-01 * tc[0]
            -1.168078145000000e-02 * tc[1]
            +5.919530250000000e-06 * tc[2]
            -2.334603641666667e-09 * tc[3]
            +4.250364870000000e-13 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.152220550000000e+04 * invT
            +1.624601760000000e+00
            -4.291424920000000e+00 * tc[0]
            +2.750771350000000e-03 * tc[1]
            -9.990638133333334e-06 * tc[2]
            +5.903885708333334e-09 * tc[3]
            -1.343428855000000e-12 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            -7.270000000000000e+03 * invT
            -1.007981170000000e+01
            -2.135836300000000e+00 * tc[0]
            -9.059436050000000e-03 * tc[1]
            +2.899124566666666e-06 * tc[2]
            -7.786646400000000e-10 * tc[3]
            +1.007288075000000e-13 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +1.924562900000000e+04 * invT
            -1.581003050000000e+01
            -1.363183500000000e+00 * tc[0]
            -9.906910499999999e-03 * tc[1]
            -2.082843333333333e-06 * tc[2]
            +2.779629583333333e-09 * tc[3]
            -7.923285500000000e-13 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            +3.977959000000000e+04 * invT
            -4.502568030000000e+01
            +3.693145300000000e+00 * tc[0]
            -2.608948400000000e-02 * tc[1]
            +4.259737833333333e-06 * tc[2]
            +5.888426750000000e-10 * tc[3]
            -3.791698750000000e-13 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            +2.330702700000000e+04 * invT
            -2.306770460000000e+01
            -4.811154000000000e-01 * tc[0]
            -1.925641600000000e-02 * tc[1]
            -5.476915333333333e-06 * tc[2]
            +6.414393416666667e-09 * tc[3]
            -1.771153400000000e-12 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            -1.761104700000000e+04 * invT
            -3.019144305000000e+01
            +9.519300500000000e-01 * tc[0]
            -2.892122250000000e-02 * tc[1]
            +6.369073166666667e-06 * tc[2]
            -3.859388000000000e-10 * tc[3]
            -1.814832550000000e-13 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            -6.116934900000000e+03 * invT
            -4.339446840000000e+01
            +3.162733400000000e+00 * tc[0]
            -3.318462250000000e-02 * tc[1]
            +5.802725500000000e-06 * tc[2]
            +5.249948083333333e-10 * tc[3]
            -4.290355050000000e-13 * tc[4];
        /*species 28: N2 */
        species[28] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
    } else {
        /*species 0: POSF10325 */
        species[0] =
            -4.877881600000000e+04 * invT
            +1.190462710000000e+02
            -2.322816300000000e+01 * tc[0]
            -3.309419400000000e-02 * tc[1]
            +4.026416666666667e-06 * tc[2]
            -2.988791416666667e-10 * tc[3]
            +5.113449000000000e-15 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +4.939886140000000e+03 * invT
            -8.269258140000002e+00
            -2.036111160000000e+00 * tc[0]
            -7.322707550000000e-03 * tc[1]
            +1.118463191666667e-06 * tc[2]
            -1.226857691666667e-10 * tc[3]
            +6.285303050000000e-15 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            -9.235703000000000e+02 * invT
            +2.004560700000000e+01
            -6.732257000000000e+00 * tc[0]
            -7.454170000000000e-03 * tc[1]
            +8.249831666666666e-07 * tc[2]
            -6.010018333333334e-11 * tc[3]
            +1.883102000000000e-15 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            -2.139723100000000e+03 * invT
            -1.348961690000000e+01
            -2.053584100000000e+00 * tc[0]
            -1.717525350000000e-02 * tc[1]
            +2.647199500000000e-06 * tc[2]
            -2.757471833333334e-10 * tc[3]
            +1.268052250000000e-14 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            -5.006675800000000e+03 * invT
            +3.393792100000000e+00
            -4.460947000000000e+00 * tc[0]
            -1.480574350000000e-02 * tc[1]
            +2.179521500000000e-06 * tc[2]
            -2.214327833333333e-10 * tc[3]
            +1.006735650000000e-14 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            +5.204346200000000e+03 * invT
            +3.825378950000000e+01
            -9.138124500000000e+00 * tc[0]
            -1.192721650000000e-02 * tc[1]
            +1.468795433333333e-06 * tc[2]
            -1.008251750000000e-10 * tc[3]
            +9.110751500000001e-16 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            -6.976490800000000e+02 * invT
            +5.966881900000000e+01
            -1.294003400000000e+01 * tc[0]
            -1.334564350000000e-02 * tc[1]
            +1.613975083333333e-06 * tc[2]
            -1.311552416666667e-10 * tc[3]
            +4.733180050000000e-15 * tc[4];
        /*species 9: H */
        species[9] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 10: O */
        species[10] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 11: OH */
        species[11] =
            +3.718857740000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 13: H2O */
        species[13] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 15: O2 */
        species[15] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +1.677558430000000e+04 * invT
            -6.194354070000000e+00
            -2.285717720000000e+00 * tc[0]
            -3.619950185000000e-03 * tc[1]
            +4.978572466666667e-07 * tc[2]
            -4.964038700000000e-11 * tc[3]
            +2.335771970000000e-15 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.399583230000000e+04 * invT
            -1.189563292000000e+01
            -1.760690080000000e+00 * tc[0]
            -4.600000410000000e-03 * tc[1]
            +7.370980216666666e-07 * tc[2]
            -8.386767666666666e-11 * tc[3]
            +4.419278200000001e-15 * tc[4];
        /*species 18: CO */
        species[18] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 19: CO2 */
        species[19] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            +2.593599920000000e+04 * invT
            +5.377850850000001e+00
            -4.147569640000000e+00 * tc[0]
            -2.980833320000000e-03 * tc[1]
            +3.954914200000000e-07 * tc[2]
            -3.895101425000000e-11 * tc[3]
            +1.806176065000000e-15 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.142639320000000e+04 * invT
            -1.404372920000000e+01
            -1.071881500000000e+00 * tc[0]
            -1.084263385000000e-02 * tc[1]
            +1.670934450000000e-06 * tc[2]
            -1.845100008333333e-10 * tc[3]
            +9.500144500000000e-15 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            -7.778500000000000e+03 * invT
            +3.879050115000000e+00
            -4.511297320000000e+00 * tc[0]
            -4.501798725000000e-03 * tc[1]
            +6.948993916666666e-07 * tc[2]
            -7.694549016666667e-11 * tc[3]
            +3.974191005000000e-15 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +1.748244900000000e+04 * invT
            +1.774383770000000e+01
            -6.500787700000000e+00 * tc[0]
            -7.162365500000000e-03 * tc[1]
            +9.463605333333332e-07 * tc[2]
            -9.234000833333333e-11 * tc[3]
            +4.518194349999999e-15 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            +3.626104700000000e+04 * invT
            +3.155195400000000e+01
            -8.597310999999999e+00 * tc[0]
            -1.112081500000000e-02 * tc[1]
            +1.453332966666666e-06 * tc[2]
            -1.149065416666667e-10 * tc[3]
            +2.657302800000000e-15 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            +1.856420300000000e+04 * invT
            +6.570956900000000e+01
            -1.404398000000000e+01 * tc[0]
            -1.174693650000000e-02 * tc[1]
            +1.422922783333333e-06 * tc[2]
            -1.157570083333333e-10 * tc[3]
            +4.180721000000000e-15 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            -2.108577000000000e+04 * invT
            +4.803129300000000e+01
            -1.173084000000000e+01 * tc[0]
            -1.180749750000000e-02 * tc[1]
            +1.705762666666667e-06 * tc[2]
            -1.627681166666667e-10 * tc[3]
            +6.373011000000000e-15 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            -1.101974400000000e+04 * invT
            +6.161653300000000e+01
            -1.365073700000000e+01 * tc[0]
            -1.284020950000000e-02 * tc[1]
            +1.744454833333333e-06 * tc[2]
            -1.617785833333333e-10 * tc[3]
            +6.741896000000000e-15 * tc[4];
        /*species 28: N2 */
        species[28] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void helmholtz(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            -4.04463440e+04 * invT
            -1.83651083e+01
            -3.16857770e+00 * tc[0]
            -3.56017240e-02 * tc[1]
            -1.44993258e-05 * tc[2]
            +1.40510050e-08 * tc[3]
            -3.79931850e-12 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +5.08977593e+03 * invT
            -1.13812948e+00
            -3.95920148e+00 * tc[0]
            +3.78526124e-03 * tc[1]
            -9.51650487e-06 * tc[2]
            +5.76323961e-09 * tc[3]
            -1.34942187e-12 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.17935173e+02 * invT
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            +1.07482600e+03 * invT
            -1.56520330e+01
            -1.49330700e+00 * tc[0]
            -1.04625900e-02 * tc[1]
            -7.47799000e-07 * tc[2]
            +1.39076000e-09 * tc[3]
            -3.57907300e-13 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            -1.79040040e+03 * invT
            -2.08813310e+01
            -1.18113800e+00 * tc[0]
            -1.54266900e-02 * tc[1]
            -8.47754117e-07 * tc[2]
            +2.05457400e-09 * tc[3]
            -5.55509650e-13 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            -4.03730690e+03 * invT
            -1.10292475e+01
            -2.64714050e+00 * tc[0]
            -1.29514785e-02 * tc[1]
            -1.36642257e-06 * tc[2]
            +1.84943825e-09 * tc[3]
            -4.44792900e-13 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            +9.18177730e+03 * invT
            -4.97336054e+01
            +4.84377340e+00 * tc[0]
            -2.92138065e-02 * tc[1]
            +4.91430917e-06 * tc[2]
            +5.78253667e-10 * tc[3]
            -4.10626265e-13 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            +4.07563000e+03 * invT
            -1.96669437e+01
            -1.61526630e+00 * tc[0]
            -1.05497190e-02 * tc[1]
            -1.42276697e-05 * tc[2]
            +1.10508883e-08 * tc[3]
            -2.79783020e-12 * tc[4];
        /*species 9: H */
        species[9] =
            +2.54736599e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 10: O */
        species[10] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 11: OH */
        species[11] =
            +3.38153812e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 13: H2O */
        species[13] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            -1.77025821e+04 * invT
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 15: O2 */
        species[15] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +1.64449988e+04 * invT
            +1.06902607e+00
            -3.67359040e+00 * tc[0]
            -1.00547588e-03 * tc[1]
            -9.55036427e-07 * tc[2]
            +5.72597854e-10 * tc[3]
            -1.27192867e-13 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.43089567e+04 * invT
            +3.19091025e+00
            -4.79372315e+00 * tc[0]
            +4.95416684e-03 * tc[1]
            -6.22033347e-06 * tc[2]
            +3.16071051e-09 * tc[3]
            -6.58863260e-13 * tc[4];
        /*species 18: CO */
        species[18] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 19: CO2 */
        species[19] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            +2.64289807e+04 * invT
            -1.41310240e+01
            -8.08681094e-01 * tc[0]
            -1.16807815e-02 * tc[1]
            +5.91953025e-06 * tc[2]
            -2.33460364e-09 * tc[3]
            +4.25036487e-13 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.15222055e+04 * invT
            +6.24601760e-01
            -4.29142492e+00 * tc[0]
            +2.75077135e-03 * tc[1]
            -9.99063813e-06 * tc[2]
            +5.90388571e-09 * tc[3]
            -1.34342886e-12 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            -7.27000000e+03 * invT
            -1.10798117e+01
            -2.13583630e+00 * tc[0]
            -9.05943605e-03 * tc[1]
            +2.89912457e-06 * tc[2]
            -7.78664640e-10 * tc[3]
            +1.00728807e-13 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +1.92456290e+04 * invT
            -1.68100305e+01
            -1.36318350e+00 * tc[0]
            -9.90691050e-03 * tc[1]
            -2.08284333e-06 * tc[2]
            +2.77962958e-09 * tc[3]
            -7.92328550e-13 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            +3.97795900e+04 * invT
            -4.60256803e+01
            +3.69314530e+00 * tc[0]
            -2.60894840e-02 * tc[1]
            +4.25973783e-06 * tc[2]
            +5.88842675e-10 * tc[3]
            -3.79169875e-13 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            +2.33070270e+04 * invT
            -2.40677046e+01
            -4.81115400e-01 * tc[0]
            -1.92564160e-02 * tc[1]
            -5.47691533e-06 * tc[2]
            +6.41439342e-09 * tc[3]
            -1.77115340e-12 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            -1.76110470e+04 * invT
            -3.11914431e+01
            +9.51930050e-01 * tc[0]
            -2.89212225e-02 * tc[1]
            +6.36907317e-06 * tc[2]
            -3.85938800e-10 * tc[3]
            -1.81483255e-13 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            -6.11693490e+03 * invT
            -4.43944684e+01
            +3.16273340e+00 * tc[0]
            -3.31846225e-02 * tc[1]
            +5.80272550e-06 * tc[2]
            +5.24994808e-10 * tc[3]
            -4.29035505e-13 * tc[4];
        /*species 28: N2 */
        species[28] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
    } else {
        /*species 0: POSF10325 */
        species[0] =
            -4.87788160e+04 * invT
            +1.18046271e+02
            -2.32281630e+01 * tc[0]
            -3.30941940e-02 * tc[1]
            +4.02641667e-06 * tc[2]
            -2.98879142e-10 * tc[3]
            +5.11344900e-15 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +4.93988614e+03 * invT
            -9.26925814e+00
            -2.03611116e+00 * tc[0]
            -7.32270755e-03 * tc[1]
            +1.11846319e-06 * tc[2]
            -1.22685769e-10 * tc[3]
            +6.28530305e-15 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 3: H2 */
        species[3] =
            -9.50158922e+02 * invT
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            -9.23570300e+02 * invT
            +1.90456070e+01
            -6.73225700e+00 * tc[0]
            -7.45417000e-03 * tc[1]
            +8.24983167e-07 * tc[2]
            -6.01001833e-11 * tc[3]
            +1.88310200e-15 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            -2.13972310e+03 * invT
            -1.44896169e+01
            -2.05358410e+00 * tc[0]
            -1.71752535e-02 * tc[1]
            +2.64719950e-06 * tc[2]
            -2.75747183e-10 * tc[3]
            +1.26805225e-14 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            -5.00667580e+03 * invT
            +2.39379210e+00
            -4.46094700e+00 * tc[0]
            -1.48057435e-02 * tc[1]
            +2.17952150e-06 * tc[2]
            -2.21432783e-10 * tc[3]
            +1.00673565e-14 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            +5.20434620e+03 * invT
            +3.72537895e+01
            -9.13812450e+00 * tc[0]
            -1.19272165e-02 * tc[1]
            +1.46879543e-06 * tc[2]
            -1.00825175e-10 * tc[3]
            +9.11075150e-16 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            -6.97649080e+02 * invT
            +5.86688190e+01
            -1.29400340e+01 * tc[0]
            -1.33456435e-02 * tc[1]
            +1.61397508e-06 * tc[2]
            -1.31155242e-10 * tc[3]
            +4.73318005e-15 * tc[4];
        /*species 9: H */
        species[9] =
            +2.54736599e+04 * invT
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 10: O */
        species[10] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 11: OH */
        species[11] =
            +3.71885774e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 13: H2O */
        species[13] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            -1.78617877e+04 * invT
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 15: O2 */
        species[15] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +1.67755843e+04 * invT
            -7.19435407e+00
            -2.28571772e+00 * tc[0]
            -3.61995018e-03 * tc[1]
            +4.97857247e-07 * tc[2]
            -4.96403870e-11 * tc[3]
            +2.33577197e-15 * tc[4];
        /*species 17: CH2O */
        species[17] =
            -1.39958323e+04 * invT
            -1.28956329e+01
            -1.76069008e+00 * tc[0]
            -4.60000041e-03 * tc[1]
            +7.37098022e-07 * tc[2]
            -8.38676767e-11 * tc[3]
            +4.41927820e-15 * tc[4];
        /*species 18: CO */
        species[18] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 19: CO2 */
        species[19] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            +2.59359992e+04 * invT
            +4.37785085e+00
            -4.14756964e+00 * tc[0]
            -2.98083332e-03 * tc[1]
            +3.95491420e-07 * tc[2]
            -3.89510143e-11 * tc[3]
            +1.80617607e-15 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            -1.14263932e+04 * invT
            -1.50437292e+01
            -1.07188150e+00 * tc[0]
            -1.08426339e-02 * tc[1]
            +1.67093445e-06 * tc[2]
            -1.84510001e-10 * tc[3]
            +9.50014450e-15 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            -7.77850000e+03 * invT
            +2.87905011e+00
            -4.51129732e+00 * tc[0]
            -4.50179872e-03 * tc[1]
            +6.94899392e-07 * tc[2]
            -7.69454902e-11 * tc[3]
            +3.97419100e-15 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +1.74824490e+04 * invT
            +1.67438377e+01
            -6.50078770e+00 * tc[0]
            -7.16236550e-03 * tc[1]
            +9.46360533e-07 * tc[2]
            -9.23400083e-11 * tc[3]
            +4.51819435e-15 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            +3.62610470e+04 * invT
            +3.05519540e+01
            -8.59731100e+00 * tc[0]
            -1.11208150e-02 * tc[1]
            +1.45333297e-06 * tc[2]
            -1.14906542e-10 * tc[3]
            +2.65730280e-15 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            +1.85642030e+04 * invT
            +6.47095690e+01
            -1.40439800e+01 * tc[0]
            -1.17469365e-02 * tc[1]
            +1.42292278e-06 * tc[2]
            -1.15757008e-10 * tc[3]
            +4.18072100e-15 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            -2.10857700e+04 * invT
            +4.70312930e+01
            -1.17308400e+01 * tc[0]
            -1.18074975e-02 * tc[1]
            +1.70576267e-06 * tc[2]
            -1.62768117e-10 * tc[3]
            +6.37301100e-15 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            -1.10197440e+04 * invT
            +6.06165330e+01
            -1.36507370e+01 * tc[0]
            -1.28402095e-02 * tc[1]
            +1.74445483e-06 * tc[2]
            -1.61778583e-10 * tc[3]
            +6.74189600e-15 * tc[4];
        /*species 28: N2 */
        species[28] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void cv_R(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            +2.16857770e+00
            +7.12034480e-02 * tc[1]
            +8.69959550e-05 * tc[2]
            -1.68612060e-07 * tc[3]
            +7.59863700e-11 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +2.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 3: H2 */
        species[3] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            +4.93307000e-01
            +2.09251800e-02 * tc[1]
            +4.48679400e-06 * tc[2]
            -1.66891200e-08 * tc[3]
            +7.15814600e-12 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            +1.81138000e-01
            +3.08533800e-02 * tc[1]
            +5.08652470e-06 * tc[2]
            -2.46548880e-08 * tc[3]
            +1.11101930e-11 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            +1.64714050e+00
            +2.59029570e-02 * tc[1]
            +8.19853540e-06 * tc[2]
            -2.21932590e-08 * tc[3]
            +8.89585800e-12 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            -5.84377340e+00
            +5.84276130e-02 * tc[1]
            -2.94858550e-05 * tc[2]
            -6.93904400e-09 * tc[3]
            +8.21252530e-12 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            +6.15266300e-01
            +2.10994380e-02 * tc[1]
            +8.53660180e-05 * tc[2]
            -1.32610660e-07 * tc[3]
            +5.59566040e-11 * tc[4];
        /*species 9: H */
        species[9] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 10: O */
        species[10] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 11: OH */
        species[11] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 13: H2O */
        species[13] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 15: O2 */
        species[15] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +2.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +3.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 18: CO */
        species[18] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            -1.91318906e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +3.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            +1.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +3.63183500e-01
            +1.98138210e-02 * tc[1]
            +1.24970600e-05 * tc[2]
            -3.33555550e-08 * tc[3]
            +1.58465710e-11 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            -4.69314530e+00
            +5.21789680e-02 * tc[1]
            -2.55584270e-05 * tc[2]
            -7.06611210e-09 * tc[3]
            +7.58339750e-12 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            -5.18884600e-01
            +3.85128320e-02 * tc[1]
            +3.28614920e-05 * tc[2]
            -7.69727210e-08 * tc[3]
            +3.54230680e-11 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            -1.95193005e+00
            +5.78424450e-02 * tc[1]
            -3.82144390e-05 * tc[2]
            +4.63126560e-09 * tc[3]
            +3.62966510e-12 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            -4.16273340e+00
            +6.63692450e-02 * tc[1]
            -3.48163530e-05 * tc[2]
            -6.29993770e-09 * tc[3]
            +8.58071010e-12 * tc[4];
        /*species 28: N2 */
        species[28] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: POSF10325 */
        species[0] =
            +2.22281630e+01
            +6.61883880e-02 * tc[1]
            -2.41585000e-05 * tc[2]
            +3.58654970e-09 * tc[3]
            -1.02268980e-13 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +1.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 3: H2 */
        species[3] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            +5.73225700e+00
            +1.49083400e-02 * tc[1]
            -4.94989900e-06 * tc[2]
            +7.21202200e-10 * tc[3]
            -3.76620400e-14 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            +1.05358410e+00
            +3.43505070e-02 * tc[1]
            -1.58831970e-05 * tc[2]
            +3.30896620e-09 * tc[3]
            -2.53610450e-13 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            +3.46094700e+00
            +2.96114870e-02 * tc[1]
            -1.30771290e-05 * tc[2]
            +2.65719340e-09 * tc[3]
            -2.01347130e-13 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            +8.13812450e+00
            +2.38544330e-02 * tc[1]
            -8.81277260e-06 * tc[2]
            +1.20990210e-09 * tc[3]
            -1.82215030e-14 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            +1.19400340e+01
            +2.66912870e-02 * tc[1]
            -9.68385050e-06 * tc[2]
            +1.57386290e-09 * tc[3]
            -9.46636010e-14 * tc[4];
        /*species 9: H */
        species[9] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 10: O */
        species[10] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 11: OH */
        species[11] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 13: H2O */
        species[13] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 15: O2 */
        species[15] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +1.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +7.60690080e-01
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 18: CO */
        species[18] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            +3.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +7.18815000e-02
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            +3.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +5.50078770e+00
            +1.43247310e-02 * tc[1]
            -5.67816320e-06 * tc[2]
            +1.10808010e-09 * tc[3]
            -9.03638870e-14 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            +7.59731100e+00
            +2.22416300e-02 * tc[1]
            -8.71999780e-06 * tc[2]
            +1.37887850e-09 * tc[3]
            -5.31460560e-14 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            +1.30439800e+01
            +2.34938730e-02 * tc[1]
            -8.53753670e-06 * tc[2]
            +1.38908410e-09 * tc[3]
            -8.36144200e-14 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            +1.07308400e+01
            +2.36149950e-02 * tc[1]
            -1.02345760e-05 * tc[2]
            +1.95321740e-09 * tc[3]
            -1.27460220e-13 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            +1.26507370e+01
            +2.56804190e-02 * tc[1]
            -1.04667290e-05 * tc[2]
            +1.94134300e-09 * tc[3]
            -1.34837920e-13 * tc[4];
        /*species 28: N2 */
        species[28] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void cp_R(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            +3.16857770e+00
            +7.12034480e-02 * tc[1]
            +8.69959550e-05 * tc[2]
            -1.68612060e-07 * tc[3]
            +7.59863700e-11 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 3: H2 */
        species[3] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            +1.49330700e+00
            +2.09251800e-02 * tc[1]
            +4.48679400e-06 * tc[2]
            -1.66891200e-08 * tc[3]
            +7.15814600e-12 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            +1.18113800e+00
            +3.08533800e-02 * tc[1]
            +5.08652470e-06 * tc[2]
            -2.46548880e-08 * tc[3]
            +1.11101930e-11 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            +2.64714050e+00
            +2.59029570e-02 * tc[1]
            +8.19853540e-06 * tc[2]
            -2.21932590e-08 * tc[3]
            +8.89585800e-12 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            -4.84377340e+00
            +5.84276130e-02 * tc[1]
            -2.94858550e-05 * tc[2]
            -6.93904400e-09 * tc[3]
            +8.21252530e-12 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            +1.61526630e+00
            +2.10994380e-02 * tc[1]
            +8.53660180e-05 * tc[2]
            -1.32610660e-07 * tc[3]
            +5.59566040e-11 * tc[4];
        /*species 9: H */
        species[9] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 10: O */
        species[10] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 11: OH */
        species[11] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 13: H2O */
        species[13] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 15: O2 */
        species[15] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 18: CO */
        species[18] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            +2.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +1.36318350e+00
            +1.98138210e-02 * tc[1]
            +1.24970600e-05 * tc[2]
            -3.33555550e-08 * tc[3]
            +1.58465710e-11 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            -3.69314530e+00
            +5.21789680e-02 * tc[1]
            -2.55584270e-05 * tc[2]
            -7.06611210e-09 * tc[3]
            +7.58339750e-12 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            +4.81115400e-01
            +3.85128320e-02 * tc[1]
            +3.28614920e-05 * tc[2]
            -7.69727210e-08 * tc[3]
            +3.54230680e-11 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            -9.51930050e-01
            +5.78424450e-02 * tc[1]
            -3.82144390e-05 * tc[2]
            +4.63126560e-09 * tc[3]
            +3.62966510e-12 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            -3.16273340e+00
            +6.63692450e-02 * tc[1]
            -3.48163530e-05 * tc[2]
            -6.29993770e-09 * tc[3]
            +8.58071010e-12 * tc[4];
        /*species 28: N2 */
        species[28] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: POSF10325 */
        species[0] =
            +2.32281630e+01
            +6.61883880e-02 * tc[1]
            -2.41585000e-05 * tc[2]
            +3.58654970e-09 * tc[3]
            -1.02268980e-13 * tc[4];
        /*species 1: C2H4 */
        species[1] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 3: H2 */
        species[3] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 4: C3H6 */
        species[4] =
            +6.73225700e+00
            +1.49083400e-02 * tc[1]
            -4.94989900e-06 * tc[2]
            +7.21202200e-10 * tc[3]
            -3.76620400e-14 * tc[4];
        /*species 5: C4H81 */
        species[5] =
            +2.05358410e+00
            +3.43505070e-02 * tc[1]
            -1.58831970e-05 * tc[2]
            +3.30896620e-09 * tc[3]
            -2.53610450e-13 * tc[4];
        /*species 6: iC4H8 */
        species[6] =
            +4.46094700e+00
            +2.96114870e-02 * tc[1]
            -1.30771290e-05 * tc[2]
            +2.65719340e-09 * tc[3]
            -2.01347130e-13 * tc[4];
        /*species 7: C6H6 */
        species[7] =
            +9.13812450e+00
            +2.38544330e-02 * tc[1]
            -8.81277260e-06 * tc[2]
            +1.20990210e-09 * tc[3]
            -1.82215030e-14 * tc[4];
        /*species 8: C6H5CH3 */
        species[8] =
            +1.29400340e+01
            +2.66912870e-02 * tc[1]
            -9.68385050e-06 * tc[2]
            +1.57386290e-09 * tc[3]
            -9.46636010e-14 * tc[4];
        /*species 9: H */
        species[9] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 10: O */
        species[10] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 11: OH */
        species[11] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 12: HO2 */
        species[12] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 13: H2O */
        species[13] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 14: H2O2 */
        species[14] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 15: O2 */
        species[15] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 16: CH3 */
        species[16] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 17: CH2O */
        species[17] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 18: CO */
        species[18] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 19: CO2 */
        species[19] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 20: C2H2 */
        species[20] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 22: CH2CO */
        species[22] =
            +4.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 23: aC3H5 */
        species[23] =
            +6.50078770e+00
            +1.43247310e-02 * tc[1]
            -5.67816320e-06 * tc[2]
            +1.10808010e-09 * tc[3]
            -9.03638870e-14 * tc[4];
        /*species 24: C6H5 */
        species[24] =
            +8.59731100e+00
            +2.22416300e-02 * tc[1]
            -8.71999780e-06 * tc[2]
            +1.37887850e-09 * tc[3]
            -5.31460560e-14 * tc[4];
        /*species 25: C6H5CH2 */
        species[25] =
            +1.40439800e+01
            +2.34938730e-02 * tc[1]
            -8.53753670e-06 * tc[2]
            +1.38908410e-09 * tc[3]
            -8.36144200e-14 * tc[4];
        /*species 26: C6H4O2 */
        species[26] =
            +1.17308400e+01
            +2.36149950e-02 * tc[1]
            -1.02345760e-05 * tc[2]
            +1.95321740e-09 * tc[3]
            -1.27460220e-13 * tc[4];
        /*species 27: C6H5CHO */
        species[27] =
            +1.36507370e+01
            +2.56804190e-02 * tc[1]
            -1.04667290e-05 * tc[2]
            +1.94134300e-09 * tc[3]
            -1.34837920e-13 * tc[4];
        /*species 28: N2 */
        species[28] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesInternalEnergy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            +2.16857770e+00
            +3.56017240e-02 * tc[1]
            +2.89986517e-05 * tc[2]
            -4.21530150e-08 * tc[3]
            +1.51972740e-11 * tc[4]
            -4.04463440e+04 * invT;
        /*species 1: C2H4 */
        species[1] =
            +2.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 2: CH4 */
        species[2] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 3: H2 */
        species[3] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 4: C3H6 */
        species[4] =
            +4.93307000e-01
            +1.04625900e-02 * tc[1]
            +1.49559800e-06 * tc[2]
            -4.17228000e-09 * tc[3]
            +1.43162920e-12 * tc[4]
            +1.07482600e+03 * invT;
        /*species 5: C4H81 */
        species[5] =
            +1.81138000e-01
            +1.54266900e-02 * tc[1]
            +1.69550823e-06 * tc[2]
            -6.16372200e-09 * tc[3]
            +2.22203860e-12 * tc[4]
            -1.79040040e+03 * invT;
        /*species 6: iC4H8 */
        species[6] =
            +1.64714050e+00
            +1.29514785e-02 * tc[1]
            +2.73284513e-06 * tc[2]
            -5.54831475e-09 * tc[3]
            +1.77917160e-12 * tc[4]
            -4.03730690e+03 * invT;
        /*species 7: C6H6 */
        species[7] =
            -5.84377340e+00
            +2.92138065e-02 * tc[1]
            -9.82861833e-06 * tc[2]
            -1.73476100e-09 * tc[3]
            +1.64250506e-12 * tc[4]
            +9.18177730e+03 * invT;
        /*species 8: C6H5CH3 */
        species[8] =
            +6.15266300e-01
            +1.05497190e-02 * tc[1]
            +2.84553393e-05 * tc[2]
            -3.31526650e-08 * tc[3]
            +1.11913208e-11 * tc[4]
            +4.07563000e+03 * invT;
        /*species 9: H */
        species[9] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 10: O */
        species[10] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 11: OH */
        species[11] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 12: HO2 */
        species[12] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 13: H2O */
        species[13] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 14: H2O2 */
        species[14] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 15: O2 */
        species[15] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 16: CH3 */
        species[16] =
            +2.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 17: CH2O */
        species[17] =
            +3.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 18: CO */
        species[18] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 20: C2H2 */
        species[20] =
            -1.91318906e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +3.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 22: CH2CO */
        species[22] =
            +1.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.27000000e+03 * invT;
        /*species 23: aC3H5 */
        species[23] =
            +3.63183500e-01
            +9.90691050e-03 * tc[1]
            +4.16568667e-06 * tc[2]
            -8.33888875e-09 * tc[3]
            +3.16931420e-12 * tc[4]
            +1.92456290e+04 * invT;
        /*species 24: C6H5 */
        species[24] =
            -4.69314530e+00
            +2.60894840e-02 * tc[1]
            -8.51947567e-06 * tc[2]
            -1.76652803e-09 * tc[3]
            +1.51667950e-12 * tc[4]
            +3.97795900e+04 * invT;
        /*species 25: C6H5CH2 */
        species[25] =
            -5.18884600e-01
            +1.92564160e-02 * tc[1]
            +1.09538307e-05 * tc[2]
            -1.92431803e-08 * tc[3]
            +7.08461360e-12 * tc[4]
            +2.33070270e+04 * invT;
        /*species 26: C6H4O2 */
        species[26] =
            -1.95193005e+00
            +2.89212225e-02 * tc[1]
            -1.27381463e-05 * tc[2]
            +1.15781640e-09 * tc[3]
            +7.25933020e-13 * tc[4]
            -1.76110470e+04 * invT;
        /*species 27: C6H5CHO */
        species[27] =
            -4.16273340e+00
            +3.31846225e-02 * tc[1]
            -1.16054510e-05 * tc[2]
            -1.57498443e-09 * tc[3]
            +1.71614202e-12 * tc[4]
            -6.11693490e+03 * invT;
        /*species 28: N2 */
        species[28] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: POSF10325 */
        species[0] =
            +2.22281630e+01
            +3.30941940e-02 * tc[1]
            -8.05283333e-06 * tc[2]
            +8.96637425e-10 * tc[3]
            -2.04537960e-14 * tc[4]
            -4.87788160e+04 * invT;
        /*species 1: C2H4 */
        species[1] =
            +1.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 2: CH4 */
        species[2] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 3: H2 */
        species[3] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 4: C3H6 */
        species[4] =
            +5.73225700e+00
            +7.45417000e-03 * tc[1]
            -1.64996633e-06 * tc[2]
            +1.80300550e-10 * tc[3]
            -7.53240800e-15 * tc[4]
            -9.23570300e+02 * invT;
        /*species 5: C4H81 */
        species[5] =
            +1.05358410e+00
            +1.71752535e-02 * tc[1]
            -5.29439900e-06 * tc[2]
            +8.27241550e-10 * tc[3]
            -5.07220900e-14 * tc[4]
            -2.13972310e+03 * invT;
        /*species 6: iC4H8 */
        species[6] =
            +3.46094700e+00
            +1.48057435e-02 * tc[1]
            -4.35904300e-06 * tc[2]
            +6.64298350e-10 * tc[3]
            -4.02694260e-14 * tc[4]
            -5.00667580e+03 * invT;
        /*species 7: C6H6 */
        species[7] =
            +8.13812450e+00
            +1.19272165e-02 * tc[1]
            -2.93759087e-06 * tc[2]
            +3.02475525e-10 * tc[3]
            -3.64430060e-15 * tc[4]
            +5.20434620e+03 * invT;
        /*species 8: C6H5CH3 */
        species[8] =
            +1.19400340e+01
            +1.33456435e-02 * tc[1]
            -3.22795017e-06 * tc[2]
            +3.93465725e-10 * tc[3]
            -1.89327202e-14 * tc[4]
            -6.97649080e+02 * invT;
        /*species 9: H */
        species[9] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 10: O */
        species[10] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 11: OH */
        species[11] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 12: HO2 */
        species[12] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 13: H2O */
        species[13] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 14: H2O2 */
        species[14] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 15: O2 */
        species[15] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 16: CH3 */
        species[16] =
            +1.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 17: CH2O */
        species[17] =
            +7.60690080e-01
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 18: CO */
        species[18] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 20: C2H2 */
        species[20] =
            +3.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +7.18815000e-02
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 22: CH2CO */
        species[22] =
            +3.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.77850000e+03 * invT;
        /*species 23: aC3H5 */
        species[23] =
            +5.50078770e+00
            +7.16236550e-03 * tc[1]
            -1.89272107e-06 * tc[2]
            +2.77020025e-10 * tc[3]
            -1.80727774e-14 * tc[4]
            +1.74824490e+04 * invT;
        /*species 24: C6H5 */
        species[24] =
            +7.59731100e+00
            +1.11208150e-02 * tc[1]
            -2.90666593e-06 * tc[2]
            +3.44719625e-10 * tc[3]
            -1.06292112e-14 * tc[4]
            +3.62610470e+04 * invT;
        /*species 25: C6H5CH2 */
        species[25] =
            +1.30439800e+01
            +1.17469365e-02 * tc[1]
            -2.84584557e-06 * tc[2]
            +3.47271025e-10 * tc[3]
            -1.67228840e-14 * tc[4]
            +1.85642030e+04 * invT;
        /*species 26: C6H4O2 */
        species[26] =
            +1.07308400e+01
            +1.18074975e-02 * tc[1]
            -3.41152533e-06 * tc[2]
            +4.88304350e-10 * tc[3]
            -2.54920440e-14 * tc[4]
            -2.10857700e+04 * invT;
        /*species 27: C6H5CHO */
        species[27] =
            +1.26507370e+01
            +1.28402095e-02 * tc[1]
            -3.48890967e-06 * tc[2]
            +4.85335750e-10 * tc[3]
            -2.69675840e-14 * tc[4]
            -1.10197440e+04 * invT;
        /*species 28: N2 */
        species[28] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesEnthalpy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            +3.16857770e+00
            +3.56017240e-02 * tc[1]
            +2.89986517e-05 * tc[2]
            -4.21530150e-08 * tc[3]
            +1.51972740e-11 * tc[4]
            -4.04463440e+04 * invT;
        /*species 1: C2H4 */
        species[1] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 * invT;
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 3: H2 */
        species[3] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 4: C3H6 */
        species[4] =
            +1.49330700e+00
            +1.04625900e-02 * tc[1]
            +1.49559800e-06 * tc[2]
            -4.17228000e-09 * tc[3]
            +1.43162920e-12 * tc[4]
            +1.07482600e+03 * invT;
        /*species 5: C4H81 */
        species[5] =
            +1.18113800e+00
            +1.54266900e-02 * tc[1]
            +1.69550823e-06 * tc[2]
            -6.16372200e-09 * tc[3]
            +2.22203860e-12 * tc[4]
            -1.79040040e+03 * invT;
        /*species 6: iC4H8 */
        species[6] =
            +2.64714050e+00
            +1.29514785e-02 * tc[1]
            +2.73284513e-06 * tc[2]
            -5.54831475e-09 * tc[3]
            +1.77917160e-12 * tc[4]
            -4.03730690e+03 * invT;
        /*species 7: C6H6 */
        species[7] =
            -4.84377340e+00
            +2.92138065e-02 * tc[1]
            -9.82861833e-06 * tc[2]
            -1.73476100e-09 * tc[3]
            +1.64250506e-12 * tc[4]
            +9.18177730e+03 * invT;
        /*species 8: C6H5CH3 */
        species[8] =
            +1.61526630e+00
            +1.05497190e-02 * tc[1]
            +2.84553393e-05 * tc[2]
            -3.31526650e-08 * tc[3]
            +1.11913208e-11 * tc[4]
            +4.07563000e+03 * invT;
        /*species 9: H */
        species[9] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 10: O */
        species[10] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 11: OH */
        species[11] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 12: HO2 */
        species[12] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 13: H2O */
        species[13] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 14: H2O2 */
        species[14] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 15: O2 */
        species[15] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 16: CH3 */
        species[16] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 * invT;
        /*species 17: CH2O */
        species[17] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 * invT;
        /*species 18: CO */
        species[18] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 20: C2H2 */
        species[20] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 * invT;
        /*species 22: CH2CO */
        species[22] =
            +2.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.27000000e+03 * invT;
        /*species 23: aC3H5 */
        species[23] =
            +1.36318350e+00
            +9.90691050e-03 * tc[1]
            +4.16568667e-06 * tc[2]
            -8.33888875e-09 * tc[3]
            +3.16931420e-12 * tc[4]
            +1.92456290e+04 * invT;
        /*species 24: C6H5 */
        species[24] =
            -3.69314530e+00
            +2.60894840e-02 * tc[1]
            -8.51947567e-06 * tc[2]
            -1.76652803e-09 * tc[3]
            +1.51667950e-12 * tc[4]
            +3.97795900e+04 * invT;
        /*species 25: C6H5CH2 */
        species[25] =
            +4.81115400e-01
            +1.92564160e-02 * tc[1]
            +1.09538307e-05 * tc[2]
            -1.92431803e-08 * tc[3]
            +7.08461360e-12 * tc[4]
            +2.33070270e+04 * invT;
        /*species 26: C6H4O2 */
        species[26] =
            -9.51930050e-01
            +2.89212225e-02 * tc[1]
            -1.27381463e-05 * tc[2]
            +1.15781640e-09 * tc[3]
            +7.25933020e-13 * tc[4]
            -1.76110470e+04 * invT;
        /*species 27: C6H5CHO */
        species[27] =
            -3.16273340e+00
            +3.31846225e-02 * tc[1]
            -1.16054510e-05 * tc[2]
            -1.57498443e-09 * tc[3]
            +1.71614202e-12 * tc[4]
            -6.11693490e+03 * invT;
        /*species 28: N2 */
        species[28] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: POSF10325 */
        species[0] =
            +2.32281630e+01
            +3.30941940e-02 * tc[1]
            -8.05283333e-06 * tc[2]
            +8.96637425e-10 * tc[3]
            -2.04537960e-14 * tc[4]
            -4.87788160e+04 * invT;
        /*species 1: C2H4 */
        species[1] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 * invT;
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 3: H2 */
        species[3] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 4: C3H6 */
        species[4] =
            +6.73225700e+00
            +7.45417000e-03 * tc[1]
            -1.64996633e-06 * tc[2]
            +1.80300550e-10 * tc[3]
            -7.53240800e-15 * tc[4]
            -9.23570300e+02 * invT;
        /*species 5: C4H81 */
        species[5] =
            +2.05358410e+00
            +1.71752535e-02 * tc[1]
            -5.29439900e-06 * tc[2]
            +8.27241550e-10 * tc[3]
            -5.07220900e-14 * tc[4]
            -2.13972310e+03 * invT;
        /*species 6: iC4H8 */
        species[6] =
            +4.46094700e+00
            +1.48057435e-02 * tc[1]
            -4.35904300e-06 * tc[2]
            +6.64298350e-10 * tc[3]
            -4.02694260e-14 * tc[4]
            -5.00667580e+03 * invT;
        /*species 7: C6H6 */
        species[7] =
            +9.13812450e+00
            +1.19272165e-02 * tc[1]
            -2.93759087e-06 * tc[2]
            +3.02475525e-10 * tc[3]
            -3.64430060e-15 * tc[4]
            +5.20434620e+03 * invT;
        /*species 8: C6H5CH3 */
        species[8] =
            +1.29400340e+01
            +1.33456435e-02 * tc[1]
            -3.22795017e-06 * tc[2]
            +3.93465725e-10 * tc[3]
            -1.89327202e-14 * tc[4]
            -6.97649080e+02 * invT;
        /*species 9: H */
        species[9] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 10: O */
        species[10] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 11: OH */
        species[11] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 12: HO2 */
        species[12] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 13: H2O */
        species[13] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 14: H2O2 */
        species[14] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 15: O2 */
        species[15] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 16: CH3 */
        species[16] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 * invT;
        /*species 17: CH2O */
        species[17] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 * invT;
        /*species 18: CO */
        species[18] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 19: CO2 */
        species[19] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 20: C2H2 */
        species[20] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 * invT;
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 * invT;
        /*species 22: CH2CO */
        species[22] =
            +4.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.77850000e+03 * invT;
        /*species 23: aC3H5 */
        species[23] =
            +6.50078770e+00
            +7.16236550e-03 * tc[1]
            -1.89272107e-06 * tc[2]
            +2.77020025e-10 * tc[3]
            -1.80727774e-14 * tc[4]
            +1.74824490e+04 * invT;
        /*species 24: C6H5 */
        species[24] =
            +8.59731100e+00
            +1.11208150e-02 * tc[1]
            -2.90666593e-06 * tc[2]
            +3.44719625e-10 * tc[3]
            -1.06292112e-14 * tc[4]
            +3.62610470e+04 * invT;
        /*species 25: C6H5CH2 */
        species[25] =
            +1.40439800e+01
            +1.17469365e-02 * tc[1]
            -2.84584557e-06 * tc[2]
            +3.47271025e-10 * tc[3]
            -1.67228840e-14 * tc[4]
            +1.85642030e+04 * invT;
        /*species 26: C6H4O2 */
        species[26] =
            +1.17308400e+01
            +1.18074975e-02 * tc[1]
            -3.41152533e-06 * tc[2]
            +4.88304350e-10 * tc[3]
            -2.54920440e-14 * tc[4]
            -2.10857700e+04 * invT;
        /*species 27: C6H5CHO */
        species[27] =
            +1.36507370e+01
            +1.28402095e-02 * tc[1]
            -3.48890967e-06 * tc[2]
            +4.85335750e-10 * tc[3]
            -2.69675840e-14 * tc[4]
            -1.10197440e+04 * invT;
        /*species 28: N2 */
        species[28] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesEntropy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: POSF10325 */
        species[0] =
            +3.16857770e+00 * tc[0]
            +7.12034480e-02 * tc[1]
            +4.34979775e-05 * tc[2]
            -5.62040200e-08 * tc[3]
            +1.89965925e-11 * tc[4]
            +2.05336860e+01 ;
        /*species 1: C2H4 */
        species[1] =
            +3.95920148e+00 * tc[0]
            -7.57052247e-03 * tc[1]
            +2.85495146e-05 * tc[2]
            -2.30529584e-08 * tc[3]
            +6.74710933e-12 * tc[4]
            +4.09733096e+00 ;
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 3: H2 */
        species[3] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
        /*species 4: C3H6 */
        species[4] =
            +1.49330700e+00 * tc[0]
            +2.09251800e-02 * tc[1]
            +2.24339700e-06 * tc[2]
            -5.56304000e-09 * tc[3]
            +1.78953650e-12 * tc[4]
            +1.61453400e+01 ;
        /*species 5: C4H81 */
        species[5] =
            +1.18113800e+00 * tc[0]
            +3.08533800e-02 * tc[1]
            +2.54326235e-06 * tc[2]
            -8.21829600e-09 * tc[3]
            +2.77754825e-12 * tc[4]
            +2.10624690e+01 ;
        /*species 6: iC4H8 */
        species[6] =
            +2.64714050e+00 * tc[0]
            +2.59029570e-02 * tc[1]
            +4.09926770e-06 * tc[2]
            -7.39775300e-09 * tc[3]
            +2.22396450e-12 * tc[4]
            +1.26763880e+01 ;
        /*species 7: C6H6 */
        species[7] =
            -4.84377340e+00 * tc[0]
            +5.84276130e-02 * tc[1]
            -1.47429275e-05 * tc[2]
            -2.31301467e-09 * tc[3]
            +2.05313133e-12 * tc[4]
            +4.38898320e+01 ;
        /*species 8: C6H5CH3 */
        species[8] =
            +1.61526630e+00 * tc[0]
            +2.10994380e-02 * tc[1]
            +4.26830090e-05 * tc[2]
            -4.42035533e-08 * tc[3]
            +1.39891510e-11 * tc[4]
            +2.02822100e+01 ;
        /*species 9: H */
        species[9] =
            +2.50000000e+00 * tc[0]
            +7.05332819e-13 * tc[1]
            -9.97959820e-16 * tc[2]
            +7.66938773e-19 * tc[3]
            -2.31933083e-22 * tc[4]
            -4.46682853e-01 ;
        /*species 10: O */
        species[10] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 11: OH */
        species[11] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 12: HO2 */
        species[12] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 13: H2O */
        species[13] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 14: H2O2 */
        species[14] =
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 15: O2 */
        species[15] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 16: CH3 */
        species[16] =
            +3.67359040e+00 * tc[0]
            +2.01095175e-03 * tc[1]
            +2.86510928e-06 * tc[2]
            -2.29039142e-09 * tc[3]
            +6.35964335e-13 * tc[4]
            +1.60456433e+00 ;
        /*species 17: CH2O */
        species[17] =
            +4.79372315e+00 * tc[0]
            -9.90833369e-03 * tc[1]
            +1.86610004e-05 * tc[2]
            -1.26428420e-08 * tc[3]
            +3.29431630e-12 * tc[4]
            +6.02812900e-01 ;
        /*species 18: CO */
        species[18] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 19: CO2 */
        species[19] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 20: C2H2 */
        species[20] =
            +8.08681094e-01 * tc[0]
            +2.33615629e-02 * tc[1]
            -1.77585907e-05 * tc[2]
            +9.33841457e-09 * tc[3]
            -2.12518243e-12 * tc[4]
            +1.39397051e+01 ;
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00 * tc[0]
            -5.50154270e-03 * tc[1]
            +2.99719144e-05 * tc[2]
            -2.36155428e-08 * tc[3]
            +6.71714427e-12 * tc[4]
            +2.66682316e+00 ;
        /*species 22: CH2CO */
        species[22] =
            +2.13583630e+00 * tc[0]
            +1.81188721e-02 * tc[1]
            -8.69737370e-06 * tc[2]
            +3.11465856e-09 * tc[3]
            -5.03644037e-13 * tc[4]
            +1.22156480e+01 ;
        /*species 23: aC3H5 */
        species[23] =
            +1.36318350e+00 * tc[0]
            +1.98138210e-02 * tc[1]
            +6.24853000e-06 * tc[2]
            -1.11185183e-08 * tc[3]
            +3.96164275e-12 * tc[4]
            +1.71732140e+01 ;
        /*species 24: C6H5 */
        species[24] =
            -3.69314530e+00 * tc[0]
            +5.21789680e-02 * tc[1]
            -1.27792135e-05 * tc[2]
            -2.35537070e-09 * tc[3]
            +1.89584937e-12 * tc[4]
            +4.13325350e+01 ;
        /*species 25: C6H5CH2 */
        species[25] =
            +4.81115400e-01 * tc[0]
            +3.85128320e-02 * tc[1]
            +1.64307460e-05 * tc[2]
            -2.56575737e-08 * tc[3]
            +8.85576700e-12 * tc[4]
            +2.35488200e+01 ;
        /*species 26: C6H4O2 */
        species[26] =
            -9.51930050e-01 * tc[0]
            +5.78424450e-02 * tc[1]
            -1.91072195e-05 * tc[2]
            +1.54375520e-09 * tc[3]
            +9.07416275e-13 * tc[4]
            +2.92395130e+01 ;
        /*species 27: C6H5CHO */
        species[27] =
            -3.16273340e+00 * tc[0]
            +6.63692450e-02 * tc[1]
            -1.74081765e-05 * tc[2]
            -2.09997923e-09 * tc[3]
            +2.14517753e-12 * tc[4]
            +4.02317350e+01 ;
        /*species 28: N2 */
        species[28] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: POSF10325 */
        species[0] =
            +2.32281630e+01 * tc[0]
            +6.61883880e-02 * tc[1]
            -1.20792500e-05 * tc[2]
            +1.19551657e-09 * tc[3]
            -2.55672450e-14 * tc[4]
            -9.58181080e+01 ;
        /*species 1: C2H4 */
        species[1] =
            +2.03611116e+00 * tc[0]
            +1.46454151e-02 * tc[1]
            -3.35538958e-06 * tc[2]
            +4.90743077e-10 * tc[3]
            -3.14265152e-14 * tc[4]
            +1.03053693e+01 ;
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 3: H2 */
        species[3] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
        /*species 4: C3H6 */
        species[4] =
            +6.73225700e+00 * tc[0]
            +1.49083400e-02 * tc[1]
            -2.47494950e-06 * tc[2]
            +2.40400733e-10 * tc[3]
            -9.41551000e-15 * tc[4]
            -1.33133500e+01 ;
        /*species 5: C4H81 */
        species[5] =
            +2.05358410e+00 * tc[0]
            +3.43505070e-02 * tc[1]
            -7.94159850e-06 * tc[2]
            +1.10298873e-09 * tc[3]
            -6.34026125e-14 * tc[4]
            +1.55432010e+01 ;
        /*species 6: iC4H8 */
        species[6] =
            +4.46094700e+00 * tc[0]
            +2.96114870e-02 * tc[1]
            -6.53856450e-06 * tc[2]
            +8.85731133e-10 * tc[3]
            -5.03367825e-14 * tc[4]
            +1.06715490e+00 ;
        /*species 7: C6H6 */
        species[7] =
            +9.13812450e+00 * tc[0]
            +2.38544330e-02 * tc[1]
            -4.40638630e-06 * tc[2]
            +4.03300700e-10 * tc[3]
            -4.55537575e-15 * tc[4]
            -2.91156650e+01 ;
        /*species 8: C6H5CH3 */
        species[8] =
            +1.29400340e+01 * tc[0]
            +2.66912870e-02 * tc[1]
            -4.84192525e-06 * tc[2]
            +5.24620967e-10 * tc[3]
            -2.36659002e-14 * tc[4]
            -4.67287850e+01 ;
        /*species 9: H */
        species[9] =
            +2.50000001e+00 * tc[0]
            -2.30842973e-11 * tc[1]
            +8.07809740e-15 * tc[2]
            -1.57838412e-18 * tc[3]
            +1.24549339e-22 * tc[4]
            -4.46682914e-01 ;
        /*species 10: O */
        species[10] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 11: OH */
        species[11] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 12: HO2 */
        species[12] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 13: H2O */
        species[13] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 14: H2O2 */
        species[14] =
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 15: O2 */
        species[15] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 16: CH3 */
        species[16] =
            +2.28571772e+00 * tc[0]
            +7.23990037e-03 * tc[1]
            -1.49357174e-06 * tc[2]
            +1.98561548e-10 * tc[3]
            -1.16788599e-14 * tc[4]
            +8.48007179e+00 ;
        /*species 17: CH2O */
        species[17] =
            +1.76069008e+00 * tc[0]
            +9.20000082e-03 * tc[1]
            -2.21129406e-06 * tc[2]
            +3.35470707e-10 * tc[3]
            -2.20963910e-14 * tc[4]
            +1.36563230e+01 ;
        /*species 18: CO */
        species[18] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 19: CO2 */
        species[19] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 20: C2H2 */
        species[20] =
            +4.14756964e+00 * tc[0]
            +5.96166664e-03 * tc[1]
            -1.18647426e-06 * tc[2]
            +1.55804057e-10 * tc[3]
            -9.03088033e-15 * tc[4]
            -1.23028121e+00 ;
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00 * tc[0]
            +2.16852677e-02 * tc[1]
            -5.01280335e-06 * tc[2]
            +7.38040003e-10 * tc[3]
            -4.75007225e-14 * tc[4]
            +1.51156107e+01 ;
        /*species 22: CH2CO */
        species[22] =
            +4.51129732e+00 * tc[0]
            +9.00359745e-03 * tc[1]
            -2.08469817e-06 * tc[2]
            +3.07781961e-10 * tc[3]
            -1.98709550e-14 * tc[4]
            +6.32247205e-01 ;
        /*species 23: aC3H5 */
        species[23] =
            +6.50078770e+00 * tc[0]
            +1.43247310e-02 * tc[1]
            -2.83908160e-06 * tc[2]
            +3.69360033e-10 * tc[3]
            -2.25909717e-14 * tc[4]
            -1.12430500e+01 ;
        /*species 24: C6H5 */
        species[24] =
            +8.59731100e+00 * tc[0]
            +2.22416300e-02 * tc[1]
            -4.35999890e-06 * tc[2]
            +4.59626167e-10 * tc[3]
            -1.32865140e-14 * tc[4]
            -2.29546430e+01 ;
        /*species 25: C6H5CH2 */
        species[25] =
            +1.40439800e+01 * tc[0]
            +2.34938730e-02 * tc[1]
            -4.26876835e-06 * tc[2]
            +4.63028033e-10 * tc[3]
            -2.09036050e-14 * tc[4]
            -5.16655890e+01 ;
        /*species 26: C6H4O2 */
        species[26] =
            +1.17308400e+01 * tc[0]
            +2.36149950e-02 * tc[1]
            -5.11728800e-06 * tc[2]
            +6.51072467e-10 * tc[3]
            -3.18650550e-14 * tc[4]
            -3.63004530e+01 ;
        /*species 27: C6H5CHO */
        species[27] =
            +1.36507370e+01 * tc[0]
            +2.56804190e-02 * tc[1]
            -5.23336450e-06 * tc[2]
            +6.47114333e-10 * tc[3]
            -3.37094800e-14 * tc[4]
            -4.79657960e+01 ;
        /*species 28: N2 */
        species[28] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }
    return;
}


/*save atomic weights into array */
void atomicWeight(double *  awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */
    awt[4] = 39.948000; /*AR */

    return;
}


/* get temperature given internal energy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, &emin);
    CKUBMS(&tmax, y, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,&e1);
        CKCVBS(&t1,y,&cv);
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
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double hin  = *h;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, &hmin);
    CKHBMS(&tmax, y, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,&h1);
        CKCPBS(&t1,y,&cp);
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


/*compute the critical parameters for each species */
void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i)
{

    double   EPS[29];
    double   SIG[29];
    double    wt[29];
    double avogadro = 6.02214199e23;
    double boltzmann = 1.3806503e-16; //we work in CGS
    double Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    get_mw(wt);

    /*species 0: POSF10325 */
    Tci[0] = 1.316 * EPS[0] ; 
    ai[0] = (5.55 * pow(avogadro,2.0) * EPS[0]*boltzmann * pow(1e-8*SIG[0],3.0) ) / (pow(wt[0],2.0)); 
    bi[0] = 0.855 * avogadro * pow(1e-8*SIG[0],3.0) / (wt[0]); 
    acentric_i[0] = 0.0 ;

    /*species 1: C2H4 */
    /*Imported from NIST */
    Tci[1] = 282.340000 ; 
    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(28.054000,2.0) * 50.410000); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (28.054000 * 50.410000); 
    acentric_i[1] = 0.087000 ;

    /*species 2: CH4 */
    /*Imported from NIST */
    Tci[2] = 190.560000 ; 
    ai[2] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[2],2.0) / (pow(16.043030,2.0) * 45.990000); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (16.043030 * 45.990000); 
    acentric_i[2] = 0.011000 ;

    /*species 3: H2 */
    /*Imported from NIST */
    Tci[3] = 33.145000 ; 
    ai[3] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[3],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (2.015880 * 12.964000); 
    acentric_i[3] = -0.219000 ;

    /*species 4: C3H6 */
    Tci[4] = 1.316 * EPS[4] ; 
    ai[4] = (5.55 * pow(avogadro,2.0) * EPS[4]*boltzmann * pow(1e-8*SIG[4],3.0) ) / (pow(wt[4],2.0)); 
    bi[4] = 0.855 * avogadro * pow(1e-8*SIG[4],3.0) / (wt[4]); 
    acentric_i[4] = 0.0 ;

    /*species 5: C4H81 */
    Tci[5] = 1.316 * EPS[5] ; 
    ai[5] = (5.55 * pow(avogadro,2.0) * EPS[5]*boltzmann * pow(1e-8*SIG[5],3.0) ) / (pow(wt[5],2.0)); 
    bi[5] = 0.855 * avogadro * pow(1e-8*SIG[5],3.0) / (wt[5]); 
    acentric_i[5] = 0.0 ;

    /*species 6: iC4H8 */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: C6H6 */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: C6H5CH3 */
    Tci[8] = 1.316 * EPS[8] ; 
    ai[8] = (5.55 * pow(avogadro,2.0) * EPS[8]*boltzmann * pow(1e-8*SIG[8],3.0) ) / (pow(wt[8],2.0)); 
    bi[8] = 0.855 * avogadro * pow(1e-8*SIG[8],3.0) / (wt[8]); 
    acentric_i[8] = 0.0 ;

    /*species 9: H */
    Tci[9] = 1.316 * EPS[9] ; 
    ai[9] = (5.55 * pow(avogadro,2.0) * EPS[9]*boltzmann * pow(1e-8*SIG[9],3.0) ) / (pow(wt[9],2.0)); 
    bi[9] = 0.855 * avogadro * pow(1e-8*SIG[9],3.0) / (wt[9]); 
    acentric_i[9] = 0.0 ;

    /*species 10: O */
    Tci[10] = 1.316 * EPS[10] ; 
    ai[10] = (5.55 * pow(avogadro,2.0) * EPS[10]*boltzmann * pow(1e-8*SIG[10],3.0) ) / (pow(wt[10],2.0)); 
    bi[10] = 0.855 * avogadro * pow(1e-8*SIG[10],3.0) / (wt[10]); 
    acentric_i[10] = 0.0 ;

    /*species 11: OH */
    Tci[11] = 1.316 * EPS[11] ; 
    ai[11] = (5.55 * pow(avogadro,2.0) * EPS[11]*boltzmann * pow(1e-8*SIG[11],3.0) ) / (pow(wt[11],2.0)); 
    bi[11] = 0.855 * avogadro * pow(1e-8*SIG[11],3.0) / (wt[11]); 
    acentric_i[11] = 0.0 ;

    /*species 12: HO2 */
    Tci[12] = 1.316 * EPS[12] ; 
    ai[12] = (5.55 * pow(avogadro,2.0) * EPS[12]*boltzmann * pow(1e-8*SIG[12],3.0) ) / (pow(wt[12],2.0)); 
    bi[12] = 0.855 * avogadro * pow(1e-8*SIG[12],3.0) / (wt[12]); 
    acentric_i[12] = 0.0 ;

    /*species 13: H2O */
    /*Imported from NIST */
    Tci[13] = 647.096000 ; 
    ai[13] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[13],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[13] = 0.08664 * Rcst * Tci[13] / (18.015340 * 220.640000); 
    acentric_i[13] = 0.344300 ;

    /*species 14: H2O2 */
    Tci[14] = 1.316 * EPS[14] ; 
    ai[14] = (5.55 * pow(avogadro,2.0) * EPS[14]*boltzmann * pow(1e-8*SIG[14],3.0) ) / (pow(wt[14],2.0)); 
    bi[14] = 0.855 * avogadro * pow(1e-8*SIG[14],3.0) / (wt[14]); 
    acentric_i[14] = 0.0 ;

    /*species 15: O2 */
    /*Imported from NIST */
    Tci[15] = 154.581000 ; 
    ai[15] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[15],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[15] = 0.08664 * Rcst * Tci[15] / (31.998800 * 50.430466); 
    acentric_i[15] = 0.022200 ;

    /*species 16: CH3 */
    Tci[16] = 1.316 * EPS[16] ; 
    ai[16] = (5.55 * pow(avogadro,2.0) * EPS[16]*boltzmann * pow(1e-8*SIG[16],3.0) ) / (pow(wt[16],2.0)); 
    bi[16] = 0.855 * avogadro * pow(1e-8*SIG[16],3.0) / (wt[16]); 
    acentric_i[16] = 0.0 ;

    /*species 17: CH2O */
    Tci[17] = 1.316 * EPS[17] ; 
    ai[17] = (5.55 * pow(avogadro,2.0) * EPS[17]*boltzmann * pow(1e-8*SIG[17],3.0) ) / (pow(wt[17],2.0)); 
    bi[17] = 0.855 * avogadro * pow(1e-8*SIG[17],3.0) / (wt[17]); 
    acentric_i[17] = 0.0 ;

    /*species 18: CO */
    /*Imported from NIST */
    Tci[18] = 132.850000 ; 
    ai[18] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[18],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[18] = 0.08664 * Rcst * Tci[18] / (28.010000 * 34.940000); 
    acentric_i[18] = 0.045000 ;

    /*species 19: CO2 */
    /*Imported from NIST */
    Tci[19] = 304.120000 ; 
    ai[19] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[19],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[19] = 0.08664 * Rcst * Tci[19] / (44.009950 * 73.740000); 
    acentric_i[19] = 0.225000 ;

    /*species 20: C2H2 */
    /*Imported from NIST */
    Tci[20] = 308.300000 ; 
    ai[20] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[20],2.0) / (pow(26.038000,2.0) * 61.140000); 
    bi[20] = 0.08664 * Rcst * Tci[20] / (26.038000 * 61.140000); 
    acentric_i[20] = 0.189000 ;

    /*species 21: C2H6 */
    /*Imported from NIST */
    Tci[21] = 305.320000 ; 
    ai[21] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[21],2.0) / (pow(30.070120,2.0) * 48.720000); 
    bi[21] = 0.08664 * Rcst * Tci[21] / (30.070120 * 48.720000); 
    acentric_i[21] = 0.099000 ;

    /*species 22: CH2CO */
    Tci[22] = 1.316 * EPS[22] ; 
    ai[22] = (5.55 * pow(avogadro,2.0) * EPS[22]*boltzmann * pow(1e-8*SIG[22],3.0) ) / (pow(wt[22],2.0)); 
    bi[22] = 0.855 * avogadro * pow(1e-8*SIG[22],3.0) / (wt[22]); 
    acentric_i[22] = 0.0 ;

    /*species 23: aC3H5 */
    Tci[23] = 1.316 * EPS[23] ; 
    ai[23] = (5.55 * pow(avogadro,2.0) * EPS[23]*boltzmann * pow(1e-8*SIG[23],3.0) ) / (pow(wt[23],2.0)); 
    bi[23] = 0.855 * avogadro * pow(1e-8*SIG[23],3.0) / (wt[23]); 
    acentric_i[23] = 0.0 ;

    /*species 24: C6H5 */
    Tci[24] = 1.316 * EPS[24] ; 
    ai[24] = (5.55 * pow(avogadro,2.0) * EPS[24]*boltzmann * pow(1e-8*SIG[24],3.0) ) / (pow(wt[24],2.0)); 
    bi[24] = 0.855 * avogadro * pow(1e-8*SIG[24],3.0) / (wt[24]); 
    acentric_i[24] = 0.0 ;

    /*species 25: C6H5CH2 */
    Tci[25] = 1.316 * EPS[25] ; 
    ai[25] = (5.55 * pow(avogadro,2.0) * EPS[25]*boltzmann * pow(1e-8*SIG[25],3.0) ) / (pow(wt[25],2.0)); 
    bi[25] = 0.855 * avogadro * pow(1e-8*SIG[25],3.0) / (wt[25]); 
    acentric_i[25] = 0.0 ;

    /*species 26: C6H4O2 */
    Tci[26] = 1.316 * EPS[26] ; 
    ai[26] = (5.55 * pow(avogadro,2.0) * EPS[26]*boltzmann * pow(1e-8*SIG[26],3.0) ) / (pow(wt[26],2.0)); 
    bi[26] = 0.855 * avogadro * pow(1e-8*SIG[26],3.0) / (wt[26]); 
    acentric_i[26] = 0.0 ;

    /*species 27: C6H5CHO */
    Tci[27] = 1.316 * EPS[27] ; 
    ai[27] = (5.55 * pow(avogadro,2.0) * EPS[27]*boltzmann * pow(1e-8*SIG[27],3.0) ) / (pow(wt[27],2.0)); 
    bi[27] = 0.855 * avogadro * pow(1e-8*SIG[27],3.0) / (wt[27]); 
    acentric_i[27] = 0.0 ;

    /*species 28: N2 */
    /*Imported from NIST */
    Tci[28] = 126.192000 ; 
    ai[28] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[28],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[28] = 0.08664 * Rcst * Tci[28] / (28.013400 * 33.958000); 
    acentric_i[28] = 0.037200 ;

    return;
}


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 118;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 16994;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 29;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(double* WT ) {
    WT[0] = 1.54297990E+02;
    WT[1] = 2.80541800E+01;
    WT[2] = 1.60430300E+01;
    WT[3] = 2.01594000E+00;
    WT[4] = 4.20812700E+01;
    WT[5] = 5.61083600E+01;
    WT[6] = 5.61083600E+01;
    WT[7] = 7.81147200E+01;
    WT[8] = 9.21418100E+01;
    WT[9] = 1.00797000E+00;
    WT[10] = 1.59994000E+01;
    WT[11] = 1.70073700E+01;
    WT[12] = 3.30067700E+01;
    WT[13] = 1.80153400E+01;
    WT[14] = 3.40147400E+01;
    WT[15] = 3.19988000E+01;
    WT[16] = 1.50350600E+01;
    WT[17] = 3.00264900E+01;
    WT[18] = 2.80105500E+01;
    WT[19] = 4.40099500E+01;
    WT[20] = 2.60382400E+01;
    WT[21] = 3.00701200E+01;
    WT[22] = 4.20376400E+01;
    WT[23] = 4.10733000E+01;
    WT[24] = 7.71067500E+01;
    WT[25] = 9.11338400E+01;
    WT[26] = 1.08097580E+02;
    WT[27] = 1.06125270E+02;
    WT[28] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(double* EPS ) {
    EPS[0] = 7.50460000E+02;
    EPS[1] = 2.80800000E+02;
    EPS[2] = 1.41400000E+02;
    EPS[3] = 3.80000000E+01;
    EPS[4] = 2.66800000E+02;
    EPS[5] = 3.57000000E+02;
    EPS[6] = 3.57000000E+02;
    EPS[7] = 4.64800000E+02;
    EPS[8] = 4.95300000E+02;
    EPS[9] = 1.45000000E+02;
    EPS[10] = 8.00000000E+01;
    EPS[11] = 8.00000000E+01;
    EPS[12] = 1.07400000E+02;
    EPS[13] = 5.72400000E+02;
    EPS[14] = 1.07400000E+02;
    EPS[15] = 1.07400000E+02;
    EPS[16] = 1.44000000E+02;
    EPS[17] = 4.98000000E+02;
    EPS[18] = 9.81000000E+01;
    EPS[19] = 2.44000000E+02;
    EPS[20] = 2.09000000E+02;
    EPS[21] = 2.52300000E+02;
    EPS[22] = 4.36000000E+02;
    EPS[23] = 2.66800000E+02;
    EPS[24] = 4.64800000E+02;
    EPS[25] = 4.95300000E+02;
    EPS[26] = 4.85000000E+02;
    EPS[27] = 5.93000000E+02;
    EPS[28] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(double* SIG ) {
    SIG[0] = 6.83400000E+00;
    SIG[1] = 3.97100000E+00;
    SIG[2] = 3.74600000E+00;
    SIG[3] = 2.92000000E+00;
    SIG[4] = 4.98200000E+00;
    SIG[5] = 5.17600000E+00;
    SIG[6] = 5.17600000E+00;
    SIG[7] = 5.29000000E+00;
    SIG[8] = 5.68000000E+00;
    SIG[9] = 2.05000000E+00;
    SIG[10] = 2.75000000E+00;
    SIG[11] = 2.75000000E+00;
    SIG[12] = 3.45800000E+00;
    SIG[13] = 2.60500000E+00;
    SIG[14] = 3.45800000E+00;
    SIG[15] = 3.45800000E+00;
    SIG[16] = 3.80000000E+00;
    SIG[17] = 3.59000000E+00;
    SIG[18] = 3.65000000E+00;
    SIG[19] = 3.76300000E+00;
    SIG[20] = 4.10000000E+00;
    SIG[21] = 4.30200000E+00;
    SIG[22] = 3.97000000E+00;
    SIG[23] = 4.98200000E+00;
    SIG[24] = 5.29000000E+00;
    SIG[25] = 5.68000000E+00;
    SIG[26] = 5.42500000E+00;
    SIG[27] = 5.47000000E+00;
    SIG[28] = 3.62100000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(double* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 4.30000000E-01;
    DIP[9] = 0.00000000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[11] = 0.00000000E+00;
    DIP[12] = 0.00000000E+00;
    DIP[13] = 1.84400000E+00;
    DIP[14] = 0.00000000E+00;
    DIP[15] = 0.00000000E+00;
    DIP[16] = 0.00000000E+00;
    DIP[17] = 0.00000000E+00;
    DIP[18] = 0.00000000E+00;
    DIP[19] = 0.00000000E+00;
    DIP[20] = 0.00000000E+00;
    DIP[21] = 0.00000000E+00;
    DIP[22] = 0.00000000E+00;
    DIP[23] = 0.00000000E+00;
    DIP[24] = 0.00000000E+00;
    DIP[25] = 4.30000000E-01;
    DIP[26] = 4.00000000E-01;
    DIP[27] = 2.80000000E+00;
    DIP[28] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(double* POL ) {
    POL[0] = 0.00000000E+00;
    POL[1] = 0.00000000E+00;
    POL[2] = 2.60000000E+00;
    POL[3] = 7.90000000E-01;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 1.03200000E+01;
    POL[8] = 1.23000000E+01;
    POL[9] = 0.00000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 0.00000000E+00;
    POL[12] = 0.00000000E+00;
    POL[13] = 0.00000000E+00;
    POL[14] = 0.00000000E+00;
    POL[15] = 1.60000000E+00;
    POL[16] = 0.00000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 1.95000000E+00;
    POL[19] = 2.65000000E+00;
    POL[20] = 0.00000000E+00;
    POL[21] = 0.00000000E+00;
    POL[22] = 0.00000000E+00;
    POL[23] = 0.00000000E+00;
    POL[24] = 1.03200000E+01;
    POL[25] = 1.23000000E+01;
    POL[26] = 0.00000000E+00;
    POL[27] = 0.00000000E+00;
    POL[28] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(double* ZROT ) {
    ZROT[0] = 1.00000000E+00;
    ZROT[1] = 1.50000000E+00;
    ZROT[2] = 1.30000000E+01;
    ZROT[3] = 2.80000000E+02;
    ZROT[4] = 1.00000000E+00;
    ZROT[5] = 1.00000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[7] = 0.00000000E+00;
    ZROT[8] = 1.00000000E+00;
    ZROT[9] = 0.00000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[11] = 0.00000000E+00;
    ZROT[12] = 1.00000000E+00;
    ZROT[13] = 4.00000000E+00;
    ZROT[14] = 3.80000000E+00;
    ZROT[15] = 3.80000000E+00;
    ZROT[16] = 0.00000000E+00;
    ZROT[17] = 2.00000000E+00;
    ZROT[18] = 1.80000000E+00;
    ZROT[19] = 2.10000000E+00;
    ZROT[20] = 2.50000000E+00;
    ZROT[21] = 1.50000000E+00;
    ZROT[22] = 2.00000000E+00;
    ZROT[23] = 1.00000000E+00;
    ZROT[24] = 0.00000000E+00;
    ZROT[25] = 1.00000000E+00;
    ZROT[26] = 1.00000000E+00;
    ZROT[27] = 1.00000000E+00;
    ZROT[28] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 2;
    NLIN[1] = 2;
    NLIN[2] = 2;
    NLIN[3] = 1;
    NLIN[4] = 2;
    NLIN[5] = 2;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 2;
    NLIN[9] = 0;
    NLIN[10] = 0;
    NLIN[11] = 1;
    NLIN[12] = 2;
    NLIN[13] = 2;
    NLIN[14] = 2;
    NLIN[15] = 1;
    NLIN[16] = 1;
    NLIN[17] = 2;
    NLIN[18] = 1;
    NLIN[19] = 1;
    NLIN[20] = 1;
    NLIN[21] = 2;
    NLIN[22] = 2;
    NLIN[23] = 2;
    NLIN[24] = 2;
    NLIN[25] = 2;
    NLIN[26] = 2;
    NLIN[27] = 2;
    NLIN[28] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(double* COFETA) {
    COFETA[0] = -8.80525237E+00;
    COFETA[1] = -2.43954160E+00;
    COFETA[2] = 5.73087700E-01;
    COFETA[3] = -3.16172064E-02;
    COFETA[4] = -2.50655444E+01;
    COFETA[5] = 5.33982977E+00;
    COFETA[6] = -5.89982992E-01;
    COFETA[7] = 2.47780650E-02;
    COFETA[8] = -2.00094664E+01;
    COFETA[9] = 3.57220167E+00;
    COFETA[10] = -3.87936446E-01;
    COFETA[11] = 1.71483254E-02;
    COFETA[12] = -1.38347699E+01;
    COFETA[13] = 1.00106621E+00;
    COFETA[14] = -4.98105694E-02;
    COFETA[15] = 2.31450475E-03;
    COFETA[16] = -2.51406631E+01;
    COFETA[17] = 5.30723075E+00;
    COFETA[18] = -5.89742369E-01;
    COFETA[19] = 2.49294033E-02;
    COFETA[20] = -2.48316231E+01;
    COFETA[21] = 4.94595777E+00;
    COFETA[22] = -5.12278955E-01;
    COFETA[23] = 2.03286378E-02;
    COFETA[24] = -2.48316231E+01;
    COFETA[25] = 4.94595777E+00;
    COFETA[26] = -5.12278955E-01;
    COFETA[27] = 2.03286378E-02;
    COFETA[28] = -2.15748728E+01;
    COFETA[29] = 3.36866449E+00;
    COFETA[30] = -2.66238062E-01;
    COFETA[31] = 8.00671416E-03;
    COFETA[32] = -2.02907587E+01;
    COFETA[33] = 2.74306418E+00;
    COFETA[34] = -1.73194603E-01;
    COFETA[35] = 3.50387394E-03;
    COFETA[36] = -2.04078397E+01;
    COFETA[37] = 3.65436395E+00;
    COFETA[38] = -3.98339635E-01;
    COFETA[39] = 1.75883009E-02;
    COFETA[40] = -1.50926240E+01;
    COFETA[41] = 1.92606504E+00;
    COFETA[42] = -1.73487476E-01;
    COFETA[43] = 7.82572931E-03;
    COFETA[44] = -1.50620763E+01;
    COFETA[45] = 1.92606504E+00;
    COFETA[46] = -1.73487476E-01;
    COFETA[47] = 7.82572931E-03;
    COFETA[48] = -1.71463238E+01;
    COFETA[49] = 2.68036374E+00;
    COFETA[50] = -2.72570227E-01;
    COFETA[51] = 1.21650964E-02;
    COFETA[52] = -1.05420863E+01;
    COFETA[53] = -1.37777096E+00;
    COFETA[54] = 4.20502308E-01;
    COFETA[55] = -2.40627230E-02;
    COFETA[56] = -1.71312832E+01;
    COFETA[57] = 2.68036374E+00;
    COFETA[58] = -2.72570227E-01;
    COFETA[59] = 1.21650964E-02;
    COFETA[60] = -1.71618309E+01;
    COFETA[61] = 2.68036374E+00;
    COFETA[62] = -2.72570227E-01;
    COFETA[63] = 1.21650964E-02;
    COFETA[64] = -2.02316497E+01;
    COFETA[65] = 3.63241793E+00;
    COFETA[66] = -3.95581049E-01;
    COFETA[67] = 1.74725495E-02;
    COFETA[68] = -1.98330577E+01;
    COFETA[69] = 2.69480162E+00;
    COFETA[70] = -1.65880845E-01;
    COFETA[71] = 3.14504769E-03;
    COFETA[72] = -1.66188336E+01;
    COFETA[73] = 2.40307799E+00;
    COFETA[74] = -2.36167638E-01;
    COFETA[75] = 1.05714061E-02;
    COFETA[76] = -2.40014975E+01;
    COFETA[77] = 5.14359547E+00;
    COFETA[78] = -5.74269731E-01;
    COFETA[79] = 2.44937679E-02;
    COFETA[80] = -2.33666446E+01;
    COFETA[81] = 4.80350223E+00;
    COFETA[82] = -5.38341336E-01;
    COFETA[83] = 2.32747213E-02;
    COFETA[84] = -2.46410937E+01;
    COFETA[85] = 5.19497183E+00;
    COFETA[86] = -5.78827339E-01;
    COFETA[87] = 2.46050921E-02;
    COFETA[88] = -2.23395647E+01;
    COFETA[89] = 3.86433912E+00;
    COFETA[90] = -3.41553983E-01;
    COFETA[91] = 1.17083447E-02;
    COFETA[92] = -2.51527853E+01;
    COFETA[93] = 5.30723075E+00;
    COFETA[94] = -5.89742369E-01;
    COFETA[95] = 2.49294033E-02;
    COFETA[96] = -2.15813667E+01;
    COFETA[97] = 3.36866449E+00;
    COFETA[98] = -2.66238062E-01;
    COFETA[99] = 8.00671416E-03;
    COFETA[100] = -2.02962585E+01;
    COFETA[101] = 2.74306418E+00;
    COFETA[102] = -1.73194603E-01;
    COFETA[103] = 3.50387394E-03;
    COFETA[104] = -2.05908562E+01;
    COFETA[105] = 2.96167485E+00;
    COFETA[106] = -2.05616587E-01;
    COFETA[107] = 5.06998104E-03;
    COFETA[108] = -1.45584504E+01;
    COFETA[109] = 2.59853727E-01;
    COFETA[110] = 1.82996694E-01;
    COFETA[111] = -1.32040000E-02;
    COFETA[112] = -1.65695594E+01;
    COFETA[113] = 2.39056562E+00;
    COFETA[114] = -2.34558144E-01;
    COFETA[115] = 1.05024037E-02;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(double* COFLAM) {
    COFLAM[0] = -1.47353053E+01;
    COFLAM[1] = 5.55326618E+00;
    COFLAM[2] = -3.20904193E-01;
    COFLAM[3] = 1.98082505E-03;
    COFLAM[4] = -1.46152839E+01;
    COFLAM[5] = 6.36251406E+00;
    COFLAM[6] = -5.03832130E-01;
    COFLAM[7] = 1.26121050E-02;
    COFLAM[8] = 1.33091602E+01;
    COFLAM[9] = -4.96140261E+00;
    COFLAM[10] = 1.03295505E+00;
    COFLAM[11] = -5.63420090E-02;
    COFLAM[12] = 9.24084392E+00;
    COFLAM[13] = -4.69568931E-01;
    COFLAM[14] = 1.15980279E-01;
    COFLAM[15] = -2.61033830E-03;
    COFLAM[16] = -1.70514683E+01;
    COFLAM[17] = 7.35906811E+00;
    COFLAM[18] = -6.52820880E-01;
    COFLAM[19] = 1.99982154E-02;
    COFLAM[20] = -1.45831906E+01;
    COFLAM[21] = 5.90742506E+00;
    COFLAM[22] = -3.94364916E-01;
    COFLAM[23] = 5.56233506E-03;
    COFLAM[24] = -1.27545639E+01;
    COFLAM[25] = 5.23633389E+00;
    COFLAM[26] = -3.15210955E-01;
    COFLAM[27] = 2.58612951E-03;
    COFLAM[28] = -3.57113176E+01;
    COFLAM[29] = 1.46486345E+01;
    COFLAM[30] = -1.61123262E+00;
    COFLAM[31] = 6.19430346E-02;
    COFLAM[32] = -3.13392011E+01;
    COFLAM[33] = 1.27590033E+01;
    COFLAM[34] = -1.34469180E+00;
    COFLAM[35] = 4.96068611E-02;
    COFLAM[36] = -8.57929284E-01;
    COFLAM[37] = 3.65436395E+00;
    COFLAM[38] = -3.98339635E-01;
    COFLAM[39] = 1.75883009E-02;
    COFLAM[40] = 1.69267361E+00;
    COFLAM[41] = 1.92606504E+00;
    COFLAM[42] = -1.73487476E-01;
    COFLAM[43] = 7.82572931E-03;
    COFLAM[44] = 1.50119731E+01;
    COFLAM[45] = -3.63267854E+00;
    COFLAM[46] = 5.92839101E-01;
    COFLAM[47] = -2.62920439E-02;
    COFLAM[48] = -1.12960913E+00;
    COFLAM[49] = 2.34014305E+00;
    COFLAM[50] = -1.63245030E-01;
    COFLAM[51] = 5.80319600E-03;
    COFLAM[52] = 2.33729817E+01;
    COFLAM[53] = -8.96536433E+00;
    COFLAM[54] = 1.52828665E+00;
    COFLAM[55] = -7.58551778E-02;
    COFLAM[56] = 8.83996545E-01;
    COFLAM[57] = 1.31525428E+00;
    COFLAM[58] = 1.91774576E-02;
    COFLAM[59] = -4.41642722E-03;
    COFLAM[60] = -1.93718739E+00;
    COFLAM[61] = 2.89110219E+00;
    COFLAM[62] = -2.71096923E-01;
    COFLAM[63] = 1.15344907E-02;
    COFLAM[64] = 1.39937901E+01;
    COFLAM[65] = -4.64256494E+00;
    COFLAM[66] = 9.07728674E-01;
    COFLAM[67] = -4.77274469E-02;
    COFLAM[68] = 5.39305086E+00;
    COFLAM[69] = -2.39312375E+00;
    COFLAM[70] = 7.39585221E-01;
    COFLAM[71] = -4.58435589E-02;
    COFLAM[72] = 1.18777264E+01;
    COFLAM[73] = -3.15463949E+00;
    COFLAM[74] = 6.01973268E-01;
    COFLAM[75] = -3.03211261E-02;
    COFLAM[76] = -1.13649314E+01;
    COFLAM[77] = 5.88177395E+00;
    COFLAM[78] = -5.68651819E-01;
    COFLAM[79] = 2.03561485E-02;
    COFLAM[80] = -7.70164502E+00;
    COFLAM[81] = 4.56884453E+00;
    COFLAM[82] = -4.04747583E-01;
    COFLAM[83] = 1.40841060E-02;
    COFLAM[84] = -1.09902209E+01;
    COFLAM[85] = 4.70647707E+00;
    COFLAM[86] = -2.52272495E-01;
    COFLAM[87] = 1.75193258E-04;
    COFLAM[88] = -8.32871231E+00;
    COFLAM[89] = 3.97067262E+00;
    COFLAM[90] = -2.21252287E-01;
    COFLAM[91] = 1.47870329E-03;
    COFLAM[92] = -2.14189975E+01;
    COFLAM[93] = 9.40841118E+00;
    COFLAM[94] = -9.66247514E-01;
    COFLAM[95] = 3.55085385E-02;
    COFLAM[96] = -3.38121733E+01;
    COFLAM[97] = 1.38837806E+01;
    COFLAM[98] = -1.51167380E+00;
    COFLAM[99] = 5.75748164E-02;
    COFLAM[100] = -2.79630575E+01;
    COFLAM[101] = 1.15247530E+01;
    COFLAM[102] = -1.19585643E+00;
    COFLAM[103] = 4.35928746E-02;
    COFLAM[104] = -2.41168238E+01;
    COFLAM[105] = 1.00711211E+01;
    COFLAM[106] = -1.01480297E+00;
    COFLAM[107] = 3.58690084E-02;
    COFLAM[108] = -2.56251957E+01;
    COFLAM[109] = 1.04275397E+01;
    COFLAM[110] = -1.03208759E+00;
    COFLAM[111] = 3.55971767E-02;
    COFLAM[112] = 1.29306158E+01;
    COFLAM[113] = -3.52817362E+00;
    COFLAM[114] = 6.45499013E-01;
    COFLAM[115] = -3.19375299E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(double* COFD) {
    COFD[0] = -1.38113115E+01;
    COFD[1] = 6.10910971E-01;
    COFD[2] = 2.74613049E-01;
    COFD[3] = -1.72629105E-02;
    COFD[4] = -2.11258476E+01;
    COFD[5] = 4.49446616E+00;
    COFD[6] = -2.96370856E-01;
    COFD[7] = 9.98580587E-03;
    COFD[8] = -2.26451895E+01;
    COFD[9] = 5.54454211E+00;
    COFD[10] = -4.71329438E-01;
    COFD[11] = 1.90999657E-02;
    COFD[12] = -1.82795157E+01;
    COFD[13] = 4.61015109E+00;
    COFD[14] = -3.82545459E-01;
    COFD[15] = 1.65444749E-02;
    COFD[16] = -2.17017559E+01;
    COFD[17] = 4.61366606E+00;
    COFD[18] = -3.15066862E-01;
    COFD[19] = 1.09216645E-02;
    COFD[20] = -2.01398625E+01;
    COFD[21] = 3.77377367E+00;
    COFD[22] = -1.86478209E-01;
    COFD[23] = 4.58817582E-03;
    COFD[24] = -2.01398625E+01;
    COFD[25] = 3.77377367E+00;
    COFD[26] = -1.86478209E-01;
    COFD[27] = 4.58817582E-03;
    COFD[28] = -1.81032125E+01;
    COFD[29] = 2.75595461E+00;
    COFD[30] = -3.54899742E-02;
    COFD[31] = -2.67520318E-03;
    COFD[32] = -1.76727730E+01;
    COFD[33] = 2.50076717E+00;
    COFD[34] = 2.22033801E-03;
    COFD[35] = -4.48335559E-03;
    COFD[36] = -2.09357329E+01;
    COFD[37] = 5.52503719E+00;
    COFD[38] = -4.67605898E-01;
    COFD[39] = 1.88908964E-02;
    COFD[40] = -2.16332288E+01;
    COFD[41] = 5.40860960E+00;
    COFD[42] = -4.73010822E-01;
    COFD[43] = 1.99398149E-02;
    COFD[44] = -2.16608259E+01;
    COFD[45] = 5.40860960E+00;
    COFD[46] = -4.73010822E-01;
    COFD[47] = 1.99398149E-02;
    COFD[48] = -2.28139639E+01;
    COFD[49] = 5.61234608E+00;
    COFD[50] = -4.91375180E-01;
    COFD[51] = 2.04179832E-02;
    COFD[52] = -1.51900988E+01;
    COFD[53] = 1.91174168E+00;
    COFD[54] = 8.81438912E-02;
    COFD[55] = -8.55341166E-03;
    COFD[56] = -2.28263210E+01;
    COFD[57] = 5.61234608E+00;
    COFD[58] = -4.91375180E-01;
    COFD[59] = 2.04179832E-02;
    COFD[60] = -2.28011547E+01;
    COFD[61] = 5.61234608E+00;
    COFD[62] = -4.91375180E-01;
    COFD[63] = 2.04179832E-02;
    COFD[64] = -2.26092158E+01;
    COFD[65] = 5.53027841E+00;
    COFD[66] = -4.68614364E-01;
    COFD[67] = 1.89477473E-02;
    COFD[68] = -1.68489747E+01;
    COFD[69] = 2.48068244E+00;
    COFD[70] = 5.18653536E-03;
    COFD[71] = -4.62524223E-03;
    COFD[72] = -2.26357059E+01;
    COFD[73] = 5.58177744E+00;
    COFD[74] = -4.90187864E-01;
    COFD[75] = 2.04753987E-02;
    COFD[76] = -2.18641291E+01;
    COFD[77] = 4.80496160E+00;
    COFD[78] = -3.45476597E-01;
    COFD[79] = 1.24586877E-02;
    COFD[80] = -2.23333626E+01;
    COFD[81] = 5.13339632E+00;
    COFD[82] = -3.97833815E-01;
    COFD[83] = 1.51137805E-02;
    COFD[84] = -2.16832171E+01;
    COFD[85] = 4.73756385E+00;
    COFD[86] = -3.34700069E-01;
    COFD[87] = 1.19121851E-02;
    COFD[88] = -1.83142602E+01;
    COFD[89] = 3.06214549E+00;
    COFD[90] = -8.03954556E-02;
    COFD[91] = -5.30273859E-04;
    COFD[92] = -2.16922067E+01;
    COFD[93] = 4.61366606E+00;
    COFD[94] = -3.15066862E-01;
    COFD[95] = 1.09216645E-02;
    COFD[96] = -1.80988918E+01;
    COFD[97] = 2.75595461E+00;
    COFD[98] = -3.54899742E-02;
    COFD[99] = -2.67520318E-03;
    COFD[100] = -1.76693225E+01;
    COFD[101] = 2.50076717E+00;
    COFD[102] = 2.22033801E-03;
    COFD[103] = -4.48335559E-03;
    COFD[104] = -1.78388304E+01;
    COFD[105] = 2.57553304E+00;
    COFD[106] = -8.86366178E-03;
    COFD[107] = -3.95164836E-03;
    COFD[108] = -1.60462384E+01;
    COFD[109] = 1.75146402E+00;
    COFD[110] = 1.11279769E-01;
    COFD[111] = -9.64032668E-03;
    COFD[112] = -2.26156323E+01;
    COFD[113] = 5.57760228E+00;
    COFD[114] = -4.89797685E-01;
    COFD[115] = 2.04643263E-02;
    COFD[116] = -2.11258476E+01;
    COFD[117] = 4.49446616E+00;
    COFD[118] = -2.96370856E-01;
    COFD[119] = 9.98580587E-03;
    COFD[120] = -2.19327397E+01;
    COFD[121] = 5.60638188E+00;
    COFD[122] = -4.91272522E-01;
    COFD[123] = 2.04396264E-02;
    COFD[124] = -1.98075055E+01;
    COFD[125] = 5.02169524E+00;
    COFD[126] = -4.31582804E-01;
    COFD[127] = 1.84953568E-02;
    COFD[128] = -1.42229194E+01;
    COFD[129] = 3.38669384E+00;
    COFD[130] = -2.28784122E-01;
    COFD[131] = 1.00790953E-02;
    COFD[132] = -2.21933982E+01;
    COFD[133] = 5.59472344E+00;
    COFD[134] = -4.91421518E-01;
    COFD[135] = 2.05117088E-02;
    COFD[136] = -2.25171771E+01;
    COFD[137] = 5.58249828E+00;
    COFD[138] = -4.78873376E-01;
    COFD[139] = 1.95316774E-02;
    COFD[140] = -2.25171771E+01;
    COFD[141] = 5.58249828E+00;
    COFD[142] = -4.78873376E-01;
    COFD[143] = 1.95316774E-02;
    COFD[144] = -2.23733906E+01;
    COFD[145] = 5.38363480E+00;
    COFD[146] = -4.40406530E-01;
    COFD[147] = 1.73594979E-02;
    COFD[148] = -2.23389643E+01;
    COFD[149] = 5.29923580E+00;
    COFD[150] = -4.25997743E-01;
    COFD[151] = 1.65974476E-02;
    COFD[152] = -1.82251914E+01;
    COFD[153] = 5.05237312E+00;
    COFD[154] = -4.35182396E-01;
    COFD[155] = 1.86363074E-02;
    COFD[156] = -1.74792112E+01;
    COFD[157] = 4.29676909E+00;
    COFD[158] = -3.44085306E-01;
    COFD[159] = 1.49671135E-02;
    COFD[160] = -1.74984476E+01;
    COFD[161] = 4.29676909E+00;
    COFD[162] = -3.44085306E-01;
    COFD[163] = 1.49671135E-02;
    COFD[164] = -1.89616623E+01;
    COFD[165] = 4.68595732E+00;
    COFD[166] = -3.91842840E-01;
    COFD[167] = 1.69262542E-02;
    COFD[168] = -2.08812333E+01;
    COFD[169] = 5.08859217E+00;
    COFD[170] = -3.90525428E-01;
    COFD[171] = 1.47376395E-02;
    COFD[172] = -1.89685165E+01;
    COFD[173] = 4.68595732E+00;
    COFD[174] = -3.91842840E-01;
    COFD[175] = 1.69262542E-02;
    COFD[176] = -1.89544778E+01;
    COFD[177] = 4.68595732E+00;
    COFD[178] = -3.91842840E-01;
    COFD[179] = 1.69262542E-02;
    COFD[180] = -1.98646734E+01;
    COFD[181] = 5.04367502E+00;
    COFD[182] = -4.34153325E-01;
    COFD[183] = 1.85956055E-02;
    COFD[184] = -2.16379567E+01;
    COFD[185] = 5.29019717E+00;
    COFD[186] = -4.24502606E-01;
    COFD[187] = 1.65197343E-02;
    COFD[188] = -1.86157761E+01;
    COFD[189] = 4.55689508E+00;
    COFD[190] = -3.75937921E-01;
    COFD[191] = 1.62703488E-02;
    COFD[192] = -2.16802612E+01;
    COFD[193] = 5.52918296E+00;
    COFD[194] = -4.85360709E-01;
    COFD[195] = 2.03448006E-02;
    COFD[196] = -2.12121370E+01;
    COFD[197] = 5.39823225E+00;
    COFD[198] = -4.72294645E-01;
    COFD[199] = 1.99340225E-02;
    COFD[200] = -2.18273547E+01;
    COFD[201] = 5.55753905E+00;
    COFD[202] = -4.88136714E-01;
    COFD[203] = 2.04294957E-02;
    COFD[204] = -2.20453723E+01;
    COFD[205] = 5.44448440E+00;
    COFD[206] = -4.51529024E-01;
    COFD[207] = 1.79698119E-02;
    COFD[208] = -2.21885140E+01;
    COFD[209] = 5.59472344E+00;
    COFD[210] = -4.91421518E-01;
    COFD[211] = 2.05117088E-02;
    COFD[212] = -2.23716664E+01;
    COFD[213] = 5.38363480E+00;
    COFD[214] = -4.40406530E-01;
    COFD[215] = 1.73594979E-02;
    COFD[216] = -2.23376752E+01;
    COFD[217] = 5.29923580E+00;
    COFD[218] = -4.25997743E-01;
    COFD[219] = 1.65974476E-02;
    COFD[220] = -2.23523338E+01;
    COFD[221] = 5.32879676E+00;
    COFD[222] = -4.30989895E-01;
    COFD[223] = 1.68597740E-02;
    COFD[224] = -2.18341320E+01;
    COFD[225] = 5.02186857E+00;
    COFD[226] = -3.79735818E-01;
    COFD[227] = 1.41857239E-02;
    COFD[228] = -1.85864144E+01;
    COFD[229] = 4.54915847E+00;
    COFD[230] = -3.75000738E-01;
    COFD[231] = 1.62324821E-02;
    COFD[232] = -2.26451895E+01;
    COFD[233] = 5.54454211E+00;
    COFD[234] = -4.71329438E-01;
    COFD[235] = 1.90999657E-02;
    COFD[236] = -1.98075055E+01;
    COFD[237] = 5.02169524E+00;
    COFD[238] = -4.31582804E-01;
    COFD[239] = 1.84953568E-02;
    COFD[240] = -1.72167708E+01;
    COFD[241] = 4.16886779E+00;
    COFD[242] = -3.28518156E-01;
    COFD[243] = 1.43341626E-02;
    COFD[244] = -1.24693568E+01;
    COFD[245] = 2.76686648E+00;
    COFD[246] = -1.49120141E-01;
    COFD[247] = 6.66220432E-03;
    COFD[248] = -1.99269592E+01;
    COFD[249] = 4.95514826E+00;
    COFD[250] = -4.23691395E-01;
    COFD[251] = 1.81828318E-02;
    COFD[252] = -2.09912990E+01;
    COFD[253] = 5.28557747E+00;
    COFD[254] = -4.61402384E-01;
    COFD[255] = 1.96111546E-02;
    COFD[256] = -2.09912990E+01;
    COFD[257] = 5.28557747E+00;
    COFD[258] = -4.61402384E-01;
    COFD[259] = 1.96111546E-02;
    COFD[260] = -2.17262575E+01;
    COFD[261] = 5.48539572E+00;
    COFD[262] = -4.80731929E-01;
    COFD[263] = 2.01857298E-02;
    COFD[264] = -2.20314708E+01;
    COFD[265] = 5.54982225E+00;
    COFD[266] = -4.87416516E-01;
    COFD[267] = 2.04093655E-02;
    COFD[268] = -1.57199037E+01;
    COFD[269] = 4.19936335E+00;
    COFD[270] = -3.32311009E-01;
    COFD[271] = 1.44921003E-02;
    COFD[272] = -1.50270339E+01;
    COFD[273] = 3.46140064E+00;
    COFD[274] = -2.38440092E-01;
    COFD[275] = 1.04960087E-02;
    COFD[276] = -1.50420953E+01;
    COFD[277] = 3.46140064E+00;
    COFD[278] = -2.38440092E-01;
    COFD[279] = 1.04960087E-02;
    COFD[280] = -1.62775714E+01;
    COFD[281] = 3.79163564E+00;
    COFD[282] = -2.80257365E-01;
    COFD[283] = 1.22656902E-02;
    COFD[284] = -2.14087397E+01;
    COFD[285] = 5.57282008E+00;
    COFD[286] = -4.76690890E-01;
    COFD[287] = 1.94000719E-02;
    COFD[288] = -1.62824412E+01;
    COFD[289] = 3.79163564E+00;
    COFD[290] = -2.80257365E-01;
    COFD[291] = 1.22656902E-02;
    COFD[292] = -1.62724462E+01;
    COFD[293] = 3.79163564E+00;
    COFD[294] = -2.80257365E-01;
    COFD[295] = 1.22656902E-02;
    COFD[296] = -1.72738845E+01;
    COFD[297] = 4.19029808E+00;
    COFD[298] = -3.31177076E-01;
    COFD[299] = 1.44446234E-02;
    COFD[300] = -2.14082453E+01;
    COFD[301] = 5.55346617E+00;
    COFD[302] = -4.87783156E-01;
    COFD[303] = 2.04210886E-02;
    COFD[304] = -1.59525102E+01;
    COFD[305] = 3.66023858E+00;
    COFD[306] = -2.63401043E-01;
    COFD[307] = 1.15432000E-02;
    COFD[308] = -1.92867554E+01;
    COFD[309] = 4.83375900E+00;
    COFD[310] = -4.09146560E-01;
    COFD[311] = 1.76006599E-02;
    COFD[312] = -1.87897298E+01;
    COFD[313] = 4.66162351E+00;
    COFD[314] = -3.88920477E-01;
    COFD[315] = 1.68089648E-02;
    COFD[316] = -1.94823660E+01;
    COFD[317] = 4.87333294E+00;
    COFD[318] = -4.13769241E-01;
    COFD[319] = 1.77802244E-02;
    COFD[320] = -2.11309207E+01;
    COFD[321] = 5.41773516E+00;
    COFD[322] = -4.73414338E-01;
    COFD[323] = 1.99258685E-02;
    COFD[324] = -1.99235839E+01;
    COFD[325] = 4.95514826E+00;
    COFD[326] = -4.23691395E-01;
    COFD[327] = 1.81828318E-02;
    COFD[328] = -2.17251451E+01;
    COFD[329] = 5.48539572E+00;
    COFD[330] = -4.80731929E-01;
    COFD[331] = 2.01857298E-02;
    COFD[332] = -2.20306514E+01;
    COFD[333] = 5.54982225E+00;
    COFD[334] = -4.87416516E-01;
    COFD[335] = 2.04093655E-02;
    COFD[336] = -2.19259008E+01;
    COFD[337] = 5.53120252E+00;
    COFD[338] = -4.85560309E-01;
    COFD[339] = 2.03509814E-02;
    COFD[340] = -2.22862449E+01;
    COFD[341] = 5.58519781E+00;
    COFD[342] = -4.83910410E-01;
    COFD[343] = 1.99365500E-02;
    COFD[344] = -1.59327297E+01;
    COFD[345] = 3.65620899E+00;
    COFD[346] = -2.62933804E-01;
    COFD[347] = 1.15253223E-02;
    COFD[348] = -1.82795157E+01;
    COFD[349] = 4.61015109E+00;
    COFD[350] = -3.82545459E-01;
    COFD[351] = 1.65444749E-02;
    COFD[352] = -1.42229194E+01;
    COFD[353] = 3.38669384E+00;
    COFD[354] = -2.28784122E-01;
    COFD[355] = 1.00790953E-02;
    COFD[356] = -1.24693568E+01;
    COFD[357] = 2.76686648E+00;
    COFD[358] = -1.49120141E-01;
    COFD[359] = 6.66220432E-03;
    COFD[360] = -1.03270606E+01;
    COFD[361] = 2.19285409E+00;
    COFD[362] = -7.54492786E-02;
    COFD[363] = 3.51398213E-03;
    COFD[364] = -1.43135319E+01;
    COFD[365] = 3.31177824E+00;
    COFD[366] = -2.18945280E-01;
    COFD[367] = 9.64764419E-03;
    COFD[368] = -1.52614188E+01;
    COFD[369] = 3.64565939E+00;
    COFD[370] = -2.61726871E-01;
    COFD[371] = 1.14799244E-02;
    COFD[372] = -1.52614188E+01;
    COFD[373] = 3.64565939E+00;
    COFD[374] = -2.61726871E-01;
    COFD[375] = 1.14799244E-02;
    COFD[376] = -1.63065147E+01;
    COFD[377] = 4.02718684E+00;
    COFD[378] = -3.10842073E-01;
    COFD[379] = 1.35952684E-02;
    COFD[380] = -1.66210564E+01;
    COFD[381] = 4.10820924E+00;
    COFD[382] = -3.21095223E-01;
    COFD[383] = 1.40301968E-02;
    COFD[384] = -1.14366381E+01;
    COFD[385] = 2.78323501E+00;
    COFD[386] = -1.51214064E-01;
    COFD[387] = 6.75150012E-03;
    COFD[388] = -1.09595712E+01;
    COFD[389] = 2.30836460E+00;
    COFD[390] = -8.76339315E-02;
    COFD[391] = 3.90878445E-03;
    COFD[392] = -1.09628982E+01;
    COFD[393] = 2.30836460E+00;
    COFD[394] = -8.76339315E-02;
    COFD[395] = 3.90878445E-03;
    COFD[396] = -1.18998012E+01;
    COFD[397] = 2.57507000E+00;
    COFD[398] = -1.24033737E-01;
    COFD[399] = 5.56694959E-03;
    COFD[400] = -1.71982995E+01;
    COFD[401] = 4.63881404E+00;
    COFD[402] = -3.86139633E-01;
    COFD[403] = 1.66955081E-02;
    COFD[404] = -1.19006548E+01;
    COFD[405] = 2.57507000E+00;
    COFD[406] = -1.24033737E-01;
    COFD[407] = 5.56694959E-03;
    COFD[408] = -1.18988955E+01;
    COFD[409] = 2.57507000E+00;
    COFD[410] = -1.24033737E-01;
    COFD[411] = 5.56694959E-03;
    COFD[412] = -1.25141260E+01;
    COFD[413] = 2.77873601E+00;
    COFD[414] = -1.50637360E-01;
    COFD[415] = 6.72684281E-03;
    COFD[416] = -1.60528285E+01;
    COFD[417] = 4.11188603E+00;
    COFD[418] = -3.21540884E-01;
    COFD[419] = 1.40482564E-02;
    COFD[420] = -1.17159737E+01;
    COFD[421] = 2.48123210E+00;
    COFD[422] = -1.11322604E-01;
    COFD[423] = 4.99282389E-03;
    COFD[424] = -1.37794315E+01;
    COFD[425] = 3.23973858E+00;
    COFD[426] = -2.09989036E-01;
    COFD[427] = 9.27667906E-03;
    COFD[428] = -1.34709807E+01;
    COFD[429] = 3.09379603E+00;
    COFD[430] = -1.91268635E-01;
    COFD[431] = 8.47480224E-03;
    COFD[432] = -1.39924781E+01;
    COFD[433] = 3.26384506E+00;
    COFD[434] = -2.12947087E-01;
    COFD[435] = 9.39743888E-03;
    COFD[436] = -1.57034851E+01;
    COFD[437] = 3.93614244E+00;
    COFD[438] = -2.99111497E-01;
    COFD[439] = 1.30888229E-02;
    COFD[440] = -1.43129712E+01;
    COFD[441] = 3.31177824E+00;
    COFD[442] = -2.18945280E-01;
    COFD[443] = 9.64764419E-03;
    COFD[444] = -1.63063503E+01;
    COFD[445] = 4.02718684E+00;
    COFD[446] = -3.10842073E-01;
    COFD[447] = 1.35952684E-02;
    COFD[448] = -1.66209380E+01;
    COFD[449] = 4.10820924E+00;
    COFD[450] = -3.21095223E-01;
    COFD[451] = 1.40301968E-02;
    COFD[452] = -1.64898075E+01;
    COFD[453] = 4.08133132E+00;
    COFD[454] = -3.17684493E-01;
    COFD[455] = 1.38851127E-02;
    COFD[456] = -1.73937806E+01;
    COFD[457] = 4.40038200E+00;
    COFD[458] = -3.56980932E-01;
    COFD[459] = 1.55040067E-02;
    COFD[460] = -1.16906297E+01;
    COFD[461] = 2.47469981E+00;
    COFD[462] = -1.10436257E-01;
    COFD[463] = 4.95273813E-03;
    COFD[464] = -2.17017559E+01;
    COFD[465] = 4.61366606E+00;
    COFD[466] = -3.15066862E-01;
    COFD[467] = 1.09216645E-02;
    COFD[468] = -2.21933982E+01;
    COFD[469] = 5.59472344E+00;
    COFD[470] = -4.91421518E-01;
    COFD[471] = 2.05117088E-02;
    COFD[472] = -1.99269592E+01;
    COFD[473] = 4.95514826E+00;
    COFD[474] = -4.23691395E-01;
    COFD[475] = 1.81828318E-02;
    COFD[476] = -1.43135319E+01;
    COFD[477] = 3.31177824E+00;
    COFD[478] = -2.18945280E-01;
    COFD[479] = 9.64764419E-03;
    COFD[480] = -2.23964038E+01;
    COFD[481] = 5.56066804E+00;
    COFD[482] = -4.88405706E-01;
    COFD[483] = 2.04357330E-02;
    COFD[484] = -2.28126594E+01;
    COFD[485] = 5.58523510E+00;
    COFD[486] = -4.81201481E-01;
    COFD[487] = 1.97107111E-02;
    COFD[488] = -2.28126594E+01;
    COFD[489] = 5.58523510E+00;
    COFD[490] = -4.81201481E-01;
    COFD[491] = 1.97107111E-02;
    COFD[492] = -2.27843627E+01;
    COFD[493] = 5.43128350E+00;
    COFD[494] = -4.49151750E-01;
    COFD[495] = 1.78402439E-02;
    COFD[496] = -2.27981313E+01;
    COFD[497] = 5.36784270E+00;
    COFD[498] = -4.37680260E-01;
    COFD[499] = 1.72142944E-02;
    COFD[500] = -1.83542556E+01;
    COFD[501] = 4.98756925E+00;
    COFD[502] = -4.27526072E-01;
    COFD[503] = 1.83341755E-02;
    COFD[504] = -1.76808721E+01;
    COFD[505] = 4.24719726E+00;
    COFD[506] = -3.38206061E-01;
    COFD[507] = 1.47350654E-02;
    COFD[508] = -1.77028170E+01;
    COFD[509] = 4.24719726E+00;
    COFD[510] = -3.38206061E-01;
    COFD[511] = 1.47350654E-02;
    COFD[512] = -1.91261963E+01;
    COFD[513] = 4.61801405E+00;
    COFD[514] = -3.83535652E-01;
    COFD[515] = 1.65862513E-02;
    COFD[516] = -2.13884087E+01;
    COFD[517] = 5.17440955E+00;
    COFD[518] = -4.04678430E-01;
    COFD[519] = 1.54706350E-02;
    COFD[520] = -1.91345696E+01;
    COFD[521] = 4.61801405E+00;
    COFD[522] = -3.83535652E-01;
    COFD[523] = 1.65862513E-02;
    COFD[524] = -1.91174465E+01;
    COFD[525] = 4.61801405E+00;
    COFD[526] = -3.83535652E-01;
    COFD[527] = 1.65862513E-02;
    COFD[528] = -1.99835686E+01;
    COFD[529] = 4.97875278E+00;
    COFD[530] = -4.26485475E-01;
    COFD[531] = 1.82931933E-02;
    COFD[532] = -2.20998738E+01;
    COFD[533] = 5.36053938E+00;
    COFD[534] = -4.36434519E-01;
    COFD[535] = 1.71484255E-02;
    COFD[536] = -1.87733838E+01;
    COFD[537] = 4.49191492E+00;
    COFD[538] = -3.68041771E-01;
    COFD[539] = 1.59498676E-02;
    COFD[540] = -2.18653077E+01;
    COFD[541] = 5.47368915E+00;
    COFD[542] = -4.79424291E-01;
    COFD[543] = 2.01372920E-02;
    COFD[544] = -2.14369874E+01;
    COFD[545] = 5.37331605E+00;
    COFD[546] = -4.70491203E-01;
    COFD[547] = 1.99134666E-02;
    COFD[548] = -2.20033797E+01;
    COFD[549] = 5.51276597E+00;
    COFD[550] = -4.83701824E-01;
    COFD[551] = 2.02915297E-02;
    COFD[552] = -2.24615468E+01;
    COFD[553] = 5.49330641E+00;
    COFD[554] = -4.60498247E-01;
    COFD[555] = 1.84639199E-02;
    COFD[556] = -2.23903059E+01;
    COFD[557] = 5.56066804E+00;
    COFD[558] = -4.88405706E-01;
    COFD[559] = 2.04357330E-02;
    COFD[560] = -2.27820795E+01;
    COFD[561] = 5.43128350E+00;
    COFD[562] = -4.49151750E-01;
    COFD[563] = 1.78402439E-02;
    COFD[564] = -2.27964005E+01;
    COFD[565] = 5.36784270E+00;
    COFD[566] = -4.37680260E-01;
    COFD[567] = 1.72142944E-02;
    COFD[568] = -2.28129363E+01;
    COFD[569] = 5.39331915E+00;
    COFD[570] = -4.42116537E-01;
    COFD[571] = 1.74516613E-02;
    COFD[572] = -2.23688893E+01;
    COFD[573] = 5.11898303E+00;
    COFD[574] = -3.95452058E-01;
    COFD[575] = 1.49902966E-02;
    COFD[576] = -1.87483158E+01;
    COFD[577] = 4.48550694E+00;
    COFD[578] = -3.67277454E-01;
    COFD[579] = 1.59194755E-02;
    COFD[580] = -2.01398625E+01;
    COFD[581] = 3.77377367E+00;
    COFD[582] = -1.86478209E-01;
    COFD[583] = 4.58817582E-03;
    COFD[584] = -2.25171771E+01;
    COFD[585] = 5.58249828E+00;
    COFD[586] = -4.78873376E-01;
    COFD[587] = 1.95316774E-02;
    COFD[588] = -2.09912990E+01;
    COFD[589] = 5.28557747E+00;
    COFD[590] = -4.61402384E-01;
    COFD[591] = 1.96111546E-02;
    COFD[592] = -1.52614188E+01;
    COFD[593] = 3.64565939E+00;
    COFD[594] = -2.61726871E-01;
    COFD[595] = 1.14799244E-02;
    COFD[596] = -2.28126594E+01;
    COFD[597] = 5.58523510E+00;
    COFD[598] = -4.81201481E-01;
    COFD[599] = 1.97107111E-02;
    COFD[600] = -2.27786045E+01;
    COFD[601] = 5.40563818E+00;
    COFD[602] = -4.44444322E-01;
    COFD[603] = 1.75813146E-02;
    COFD[604] = -2.27786045E+01;
    COFD[605] = 5.40563818E+00;
    COFD[606] = -4.44444322E-01;
    COFD[607] = 1.75813146E-02;
    COFD[608] = -2.22466531E+01;
    COFD[609] = 5.02889918E+00;
    COFD[610] = -3.80860199E-01;
    COFD[611] = 1.42428137E-02;
    COFD[612] = -2.20760694E+01;
    COFD[613] = 4.88495888E+00;
    COFD[614] = -3.58090775E-01;
    COFD[615] = 1.30937796E-02;
    COFD[616] = -1.94091156E+01;
    COFD[617] = 5.32291505E+00;
    COFD[618] = -4.65883522E-01;
    COFD[619] = 1.97916109E-02;
    COFD[620] = -1.87878849E+01;
    COFD[621] = 4.61260432E+00;
    COFD[622] = -3.82854484E-01;
    COFD[623] = 1.65575163E-02;
    COFD[624] = -1.88114917E+01;
    COFD[625] = 4.61260432E+00;
    COFD[626] = -3.82854484E-01;
    COFD[627] = 1.65575163E-02;
    COFD[628] = -2.02566224E+01;
    COFD[629] = 4.97613338E+00;
    COFD[630] = -4.26175206E-01;
    COFD[631] = 1.82809270E-02;
    COFD[632] = -2.03437836E+01;
    COFD[633] = 4.57152878E+00;
    COFD[634] = -3.08371263E-01;
    COFD[635] = 1.05838559E-02;
    COFD[636] = -2.02660394E+01;
    COFD[637] = 4.97613338E+00;
    COFD[638] = -4.26175206E-01;
    COFD[639] = 1.82809270E-02;
    COFD[640] = -2.02468029E+01;
    COFD[641] = 4.97613338E+00;
    COFD[642] = -4.26175206E-01;
    COFD[643] = 1.82809270E-02;
    COFD[644] = -2.10572534E+01;
    COFD[645] = 5.31360223E+00;
    COFD[646] = -4.64787000E-01;
    COFD[647] = 1.97483720E-02;
    COFD[648] = -2.13352637E+01;
    COFD[649] = 4.87252053E+00;
    COFD[650] = -3.56127804E-01;
    COFD[651] = 1.29948788E-02;
    COFD[652] = -1.98832231E+01;
    COFD[653] = 4.84731557E+00;
    COFD[654] = -4.10638352E-01;
    COFD[655] = 1.76543886E-02;
    COFD[656] = -2.25018756E+01;
    COFD[657] = 5.59178974E+00;
    COFD[658] = -4.85668031E-01;
    COFD[659] = 2.00491907E-02;
    COFD[660] = -2.22817344E+01;
    COFD[661] = 5.59185582E+00;
    COFD[662] = -4.91155812E-01;
    COFD[663] = 2.05043018E-02;
    COFD[664] = -2.25119517E+01;
    COFD[665] = 5.58206320E+00;
    COFD[666] = -4.82956809E-01;
    COFD[667] = 1.98731634E-02;
    COFD[668] = -2.20355148E+01;
    COFD[669] = 5.14570932E+00;
    COFD[670] = -3.99877142E-01;
    COFD[671] = 1.52199557E-02;
    COFD[672] = -2.28056965E+01;
    COFD[673] = 5.58523510E+00;
    COFD[674] = -4.81201481E-01;
    COFD[675] = 1.97107111E-02;
    COFD[676] = -2.22439282E+01;
    COFD[677] = 5.02889918E+00;
    COFD[678] = -3.80860199E-01;
    COFD[679] = 1.42428137E-02;
    COFD[680] = -2.20739808E+01;
    COFD[681] = 4.88495888E+00;
    COFD[682] = -3.58090775E-01;
    COFD[683] = 1.30937796E-02;
    COFD[684] = -2.21500241E+01;
    COFD[685] = 4.93303891E+00;
    COFD[686] = -3.65678823E-01;
    COFD[687] = 1.34760462E-02;
    COFD[688] = -2.12908209E+01;
    COFD[689] = 4.48269001E+00;
    COFD[690] = -2.94539894E-01;
    COFD[691] = 9.89462789E-03;
    COFD[692] = -1.98610390E+01;
    COFD[693] = 4.84231878E+00;
    COFD[694] = -4.10101001E-01;
    COFD[695] = 1.76356687E-02;
    COFD[696] = -2.01398625E+01;
    COFD[697] = 3.77377367E+00;
    COFD[698] = -1.86478209E-01;
    COFD[699] = 4.58817582E-03;
    COFD[700] = -2.25171771E+01;
    COFD[701] = 5.58249828E+00;
    COFD[702] = -4.78873376E-01;
    COFD[703] = 1.95316774E-02;
    COFD[704] = -2.09912990E+01;
    COFD[705] = 5.28557747E+00;
    COFD[706] = -4.61402384E-01;
    COFD[707] = 1.96111546E-02;
    COFD[708] = -1.52614188E+01;
    COFD[709] = 3.64565939E+00;
    COFD[710] = -2.61726871E-01;
    COFD[711] = 1.14799244E-02;
    COFD[712] = -2.28126594E+01;
    COFD[713] = 5.58523510E+00;
    COFD[714] = -4.81201481E-01;
    COFD[715] = 1.97107111E-02;
    COFD[716] = -2.27786045E+01;
    COFD[717] = 5.40563818E+00;
    COFD[718] = -4.44444322E-01;
    COFD[719] = 1.75813146E-02;
    COFD[720] = -2.27786045E+01;
    COFD[721] = 5.40563818E+00;
    COFD[722] = -4.44444322E-01;
    COFD[723] = 1.75813146E-02;
    COFD[724] = -2.22466531E+01;
    COFD[725] = 5.02889918E+00;
    COFD[726] = -3.80860199E-01;
    COFD[727] = 1.42428137E-02;
    COFD[728] = -2.20760694E+01;
    COFD[729] = 4.88495888E+00;
    COFD[730] = -3.58090775E-01;
    COFD[731] = 1.30937796E-02;
    COFD[732] = -1.94091156E+01;
    COFD[733] = 5.32291505E+00;
    COFD[734] = -4.65883522E-01;
    COFD[735] = 1.97916109E-02;
    COFD[736] = -1.87878849E+01;
    COFD[737] = 4.61260432E+00;
    COFD[738] = -3.82854484E-01;
    COFD[739] = 1.65575163E-02;
    COFD[740] = -1.88114917E+01;
    COFD[741] = 4.61260432E+00;
    COFD[742] = -3.82854484E-01;
    COFD[743] = 1.65575163E-02;
    COFD[744] = -2.02566224E+01;
    COFD[745] = 4.97613338E+00;
    COFD[746] = -4.26175206E-01;
    COFD[747] = 1.82809270E-02;
    COFD[748] = -2.03437836E+01;
    COFD[749] = 4.57152878E+00;
    COFD[750] = -3.08371263E-01;
    COFD[751] = 1.05838559E-02;
    COFD[752] = -2.02660394E+01;
    COFD[753] = 4.97613338E+00;
    COFD[754] = -4.26175206E-01;
    COFD[755] = 1.82809270E-02;
    COFD[756] = -2.02468029E+01;
    COFD[757] = 4.97613338E+00;
    COFD[758] = -4.26175206E-01;
    COFD[759] = 1.82809270E-02;
    COFD[760] = -2.10572534E+01;
    COFD[761] = 5.31360223E+00;
    COFD[762] = -4.64787000E-01;
    COFD[763] = 1.97483720E-02;
    COFD[764] = -2.13352637E+01;
    COFD[765] = 4.87252053E+00;
    COFD[766] = -3.56127804E-01;
    COFD[767] = 1.29948788E-02;
    COFD[768] = -1.98832231E+01;
    COFD[769] = 4.84731557E+00;
    COFD[770] = -4.10638352E-01;
    COFD[771] = 1.76543886E-02;
    COFD[772] = -2.25018756E+01;
    COFD[773] = 5.59178974E+00;
    COFD[774] = -4.85668031E-01;
    COFD[775] = 2.00491907E-02;
    COFD[776] = -2.22817344E+01;
    COFD[777] = 5.59185582E+00;
    COFD[778] = -4.91155812E-01;
    COFD[779] = 2.05043018E-02;
    COFD[780] = -2.25119517E+01;
    COFD[781] = 5.58206320E+00;
    COFD[782] = -4.82956809E-01;
    COFD[783] = 1.98731634E-02;
    COFD[784] = -2.20355148E+01;
    COFD[785] = 5.14570932E+00;
    COFD[786] = -3.99877142E-01;
    COFD[787] = 1.52199557E-02;
    COFD[788] = -2.28056965E+01;
    COFD[789] = 5.58523510E+00;
    COFD[790] = -4.81201481E-01;
    COFD[791] = 1.97107111E-02;
    COFD[792] = -2.22439282E+01;
    COFD[793] = 5.02889918E+00;
    COFD[794] = -3.80860199E-01;
    COFD[795] = 1.42428137E-02;
    COFD[796] = -2.20739808E+01;
    COFD[797] = 4.88495888E+00;
    COFD[798] = -3.58090775E-01;
    COFD[799] = 1.30937796E-02;
    COFD[800] = -2.21500241E+01;
    COFD[801] = 4.93303891E+00;
    COFD[802] = -3.65678823E-01;
    COFD[803] = 1.34760462E-02;
    COFD[804] = -2.12908209E+01;
    COFD[805] = 4.48269001E+00;
    COFD[806] = -2.94539894E-01;
    COFD[807] = 9.89462789E-03;
    COFD[808] = -1.98610390E+01;
    COFD[809] = 4.84231878E+00;
    COFD[810] = -4.10101001E-01;
    COFD[811] = 1.76356687E-02;
    COFD[812] = -1.81032125E+01;
    COFD[813] = 2.75595461E+00;
    COFD[814] = -3.54899742E-02;
    COFD[815] = -2.67520318E-03;
    COFD[816] = -2.23733906E+01;
    COFD[817] = 5.38363480E+00;
    COFD[818] = -4.40406530E-01;
    COFD[819] = 1.73594979E-02;
    COFD[820] = -2.17262575E+01;
    COFD[821] = 5.48539572E+00;
    COFD[822] = -4.80731929E-01;
    COFD[823] = 2.01857298E-02;
    COFD[824] = -1.63065147E+01;
    COFD[825] = 4.02718684E+00;
    COFD[826] = -3.10842073E-01;
    COFD[827] = 1.35952684E-02;
    COFD[828] = -2.27843627E+01;
    COFD[829] = 5.43128350E+00;
    COFD[830] = -4.49151750E-01;
    COFD[831] = 1.78402439E-02;
    COFD[832] = -2.22466531E+01;
    COFD[833] = 5.02889918E+00;
    COFD[834] = -3.80860199E-01;
    COFD[835] = 1.42428137E-02;
    COFD[836] = -2.22466531E+01;
    COFD[837] = 5.02889918E+00;
    COFD[838] = -3.80860199E-01;
    COFD[839] = 1.42428137E-02;
    COFD[840] = -2.12039763E+01;
    COFD[841] = 4.42962851E+00;
    COFD[842] = -2.86319550E-01;
    COFD[843] = 9.48633214E-03;
    COFD[844] = -2.09413003E+01;
    COFD[845] = 4.24466287E+00;
    COFD[846] = -2.57969164E-01;
    COFD[847] = 8.08847972E-03;
    COFD[848] = -2.01030164E+01;
    COFD[849] = 5.51409582E+00;
    COFD[850] = -4.83844807E-01;
    COFD[851] = 2.02965712E-02;
    COFD[852] = -1.97761100E+01;
    COFD[853] = 4.93652913E+00;
    COFD[854] = -4.21485300E-01;
    COFD[855] = 1.80955487E-02;
    COFD[856] = -1.98013312E+01;
    COFD[857] = 4.93652913E+00;
    COFD[858] = -4.21485300E-01;
    COFD[859] = 1.80955487E-02;
    COFD[860] = -2.12131123E+01;
    COFD[861] = 5.26829491E+00;
    COFD[862] = -4.59325230E-01;
    COFD[863] = 1.95275380E-02;
    COFD[864] = -1.74926607E+01;
    COFD[865] = 3.15361887E+00;
    COFD[866] = -9.38468185E-02;
    COFD[867] = 1.12007553E-04;
    COFD[868] = -2.12236380E+01;
    COFD[869] = 5.26829491E+00;
    COFD[870] = -4.59325230E-01;
    COFD[871] = 1.95275380E-02;
    COFD[872] = -2.12021613E+01;
    COFD[873] = 5.26829491E+00;
    COFD[874] = -4.59325230E-01;
    COFD[875] = 1.95275380E-02;
    COFD[876] = -2.17763744E+01;
    COFD[877] = 5.50642944E+00;
    COFD[878] = -4.83016239E-01;
    COFD[879] = 2.02671529E-02;
    COFD[880] = -2.01549977E+01;
    COFD[881] = 4.23025405E+00;
    COFD[882] = -2.55784321E-01;
    COFD[883] = 7.98147639E-03;
    COFD[884] = -2.08713448E+01;
    COFD[885] = 5.15826153E+00;
    COFD[886] = -4.46699782E-01;
    COFD[887] = 1.90459433E-02;
    COFD[888] = -2.26266184E+01;
    COFD[889] = 5.50195433E+00;
    COFD[890] = -4.62690815E-01;
    COFD[891] = 1.86010137E-02;
    COFD[892] = -2.25690703E+01;
    COFD[893] = 5.58516153E+00;
    COFD[894] = -4.80440252E-01;
    COFD[895] = 1.96479586E-02;
    COFD[896] = -2.26092626E+01;
    COFD[897] = 5.48733210E+00;
    COFD[898] = -4.59335784E-01;
    COFD[899] = 1.83981530E-02;
    COFD[900] = -2.10753161E+01;
    COFD[901] = 4.59015462E+00;
    COFD[902] = -3.11301450E-01;
    COFD[903] = 1.07307395E-02;
    COFD[904] = -2.27764511E+01;
    COFD[905] = 5.43128350E+00;
    COFD[906] = -4.49151750E-01;
    COFD[907] = 1.78402439E-02;
    COFD[908] = -2.12007188E+01;
    COFD[909] = 4.42962851E+00;
    COFD[910] = -2.86319550E-01;
    COFD[911] = 9.48633214E-03;
    COFD[912] = -2.09387694E+01;
    COFD[913] = 4.24466287E+00;
    COFD[914] = -2.57969164E-01;
    COFD[915] = 8.08847972E-03;
    COFD[916] = -2.10514354E+01;
    COFD[917] = 4.30471901E+00;
    COFD[918] = -2.67147055E-01;
    COFD[919] = 8.53981315E-03;
    COFD[920] = -1.94550832E+01;
    COFD[921] = 3.52672843E+00;
    COFD[922] = -1.49232080E-01;
    COFD[923] = 2.77504261E-03;
    COFD[924] = -2.08547637E+01;
    COFD[925] = 5.15581320E+00;
    COFD[926] = -4.46543388E-01;
    COFD[927] = 1.90458118E-02;
    COFD[928] = -1.76727730E+01;
    COFD[929] = 2.50076717E+00;
    COFD[930] = 2.22033801E-03;
    COFD[931] = -4.48335559E-03;
    COFD[932] = -2.23389643E+01;
    COFD[933] = 5.29923580E+00;
    COFD[934] = -4.25997743E-01;
    COFD[935] = 1.65974476E-02;
    COFD[936] = -2.20314708E+01;
    COFD[937] = 5.54982225E+00;
    COFD[938] = -4.87416516E-01;
    COFD[939] = 2.04093655E-02;
    COFD[940] = -1.66210564E+01;
    COFD[941] = 4.10820924E+00;
    COFD[942] = -3.21095223E-01;
    COFD[943] = 1.40301968E-02;
    COFD[944] = -2.27981313E+01;
    COFD[945] = 5.36784270E+00;
    COFD[946] = -4.37680260E-01;
    COFD[947] = 1.72142944E-02;
    COFD[948] = -2.20760694E+01;
    COFD[949] = 4.88495888E+00;
    COFD[950] = -3.58090775E-01;
    COFD[951] = 1.30937796E-02;
    COFD[952] = -2.20760694E+01;
    COFD[953] = 4.88495888E+00;
    COFD[954] = -3.58090775E-01;
    COFD[955] = 1.30937796E-02;
    COFD[956] = -2.09413003E+01;
    COFD[957] = 4.24466287E+00;
    COFD[958] = -2.57969164E-01;
    COFD[959] = 8.08847972E-03;
    COFD[960] = -2.05967421E+01;
    COFD[961] = 4.02612453E+00;
    COFD[962] = -2.25065111E-01;
    COFD[963] = 6.48488943E-03;
    COFD[964] = -2.03845633E+01;
    COFD[965] = 5.56619223E+00;
    COFD[966] = -4.88860312E-01;
    COFD[967] = 2.04450332E-02;
    COFD[968] = -2.01219251E+01;
    COFD[969] = 5.01930029E+00;
    COFD[970] = -4.31305027E-01;
    COFD[971] = 1.84846322E-02;
    COFD[972] = -2.01478340E+01;
    COFD[973] = 5.01930029E+00;
    COFD[974] = -4.31305027E-01;
    COFD[975] = 1.84846322E-02;
    COFD[976] = -2.15623852E+01;
    COFD[977] = 5.34961360E+00;
    COFD[978] = -4.68762622E-01;
    COFD[979] = 1.98932200E-02;
    COFD[980] = -1.83871985E+01;
    COFD[981] = 3.54020956E+00;
    COFD[982] = -1.52304845E-01;
    COFD[983] = 2.96268912E-03;
    COFD[984] = -2.15734148E+01;
    COFD[985] = 5.34961360E+00;
    COFD[986] = -4.68762622E-01;
    COFD[987] = 1.98932200E-02;
    COFD[988] = -2.15536156E+01;
    COFD[989] = 5.35030481E+00;
    COFD[990] = -4.68818560E-01;
    COFD[991] = 1.98942796E-02;
    COFD[992] = -2.20578790E+01;
    COFD[993] = 5.56193702E+00;
    COFD[994] = -4.88512557E-01;
    COFD[995] = 2.04380728E-02;
    COFD[996] = -1.98066211E+01;
    COFD[997] = 4.01137806E+00;
    COFD[998] = -2.22734310E-01;
    COFD[999] = 6.36741567E-03;
    COFD[1000] = -2.12018521E+01;
    COFD[1001] = 5.23328793E+00;
    COFD[1002] = -4.55218686E-01;
    COFD[1003] = 1.93667472E-02;
    COFD[1004] = -2.26941847E+01;
    COFD[1005] = 5.45651759E+00;
    COFD[1006] = -4.53709578E-01;
    COFD[1007] = 1.80890894E-02;
    COFD[1008] = -2.26781083E+01;
    COFD[1009] = 5.56399361E+00;
    COFD[1010] = -4.75020089E-01;
    COFD[1011] = 1.93065281E-02;
    COFD[1012] = -2.26225827E+01;
    COFD[1013] = 5.42359285E+00;
    COFD[1014] = -4.47762010E-01;
    COFD[1015] = 1.77644050E-02;
    COFD[1016] = -2.08711113E+01;
    COFD[1017] = 4.43077544E+00;
    COFD[1018] = -2.86496187E-01;
    COFD[1019] = 9.49507497E-03;
    COFD[1020] = -2.27897781E+01;
    COFD[1021] = 5.36784270E+00;
    COFD[1022] = -4.37680260E-01;
    COFD[1023] = 1.72142944E-02;
    COFD[1024] = -2.09377754E+01;
    COFD[1025] = 4.24466287E+00;
    COFD[1026] = -2.57969164E-01;
    COFD[1027] = 8.08847972E-03;
    COFD[1028] = -2.05939847E+01;
    COFD[1029] = 4.02612453E+00;
    COFD[1030] = -2.25065111E-01;
    COFD[1031] = 6.48488943E-03;
    COFD[1032] = -2.07475551E+01;
    COFD[1033] = 4.10156709E+00;
    COFD[1034] = -2.36431182E-01;
    COFD[1035] = 7.03931247E-03;
    COFD[1036] = -1.93438403E+01;
    COFD[1037] = 3.42798339E+00;
    COFD[1038] = -1.35210942E-01;
    COFD[1039] = 2.12316045E-03;
    COFD[1040] = -2.11737258E+01;
    COFD[1041] = 5.22581725E+00;
    COFD[1042] = -4.54351699E-01;
    COFD[1043] = 1.93332012E-02;
    COFD[1044] = -2.09357329E+01;
    COFD[1045] = 5.52503719E+00;
    COFD[1046] = -4.67605898E-01;
    COFD[1047] = 1.88908964E-02;
    COFD[1048] = -1.82251914E+01;
    COFD[1049] = 5.05237312E+00;
    COFD[1050] = -4.35182396E-01;
    COFD[1051] = 1.86363074E-02;
    COFD[1052] = -1.57199037E+01;
    COFD[1053] = 4.19936335E+00;
    COFD[1054] = -3.32311009E-01;
    COFD[1055] = 1.44921003E-02;
    COFD[1056] = -1.14366381E+01;
    COFD[1057] = 2.78323501E+00;
    COFD[1058] = -1.51214064E-01;
    COFD[1059] = 6.75150012E-03;
    COFD[1060] = -1.83542556E+01;
    COFD[1061] = 4.98756925E+00;
    COFD[1062] = -4.27526072E-01;
    COFD[1063] = 1.83341755E-02;
    COFD[1064] = -1.94091156E+01;
    COFD[1065] = 5.32291505E+00;
    COFD[1066] = -4.65883522E-01;
    COFD[1067] = 1.97916109E-02;
    COFD[1068] = -1.94091156E+01;
    COFD[1069] = 5.32291505E+00;
    COFD[1070] = -4.65883522E-01;
    COFD[1071] = 1.97916109E-02;
    COFD[1072] = -2.01030164E+01;
    COFD[1073] = 5.51409582E+00;
    COFD[1074] = -4.83844807E-01;
    COFD[1075] = 2.02965712E-02;
    COFD[1076] = -2.03845633E+01;
    COFD[1077] = 5.56619223E+00;
    COFD[1078] = -4.88860312E-01;
    COFD[1079] = 2.04450332E-02;
    COFD[1080] = -1.47968712E+01;
    COFD[1081] = 4.23027636E+00;
    COFD[1082] = -3.36139991E-01;
    COFD[1083] = 1.46507621E-02;
    COFD[1084] = -1.34230272E+01;
    COFD[1085] = 3.48624238E+00;
    COFD[1086] = -2.41554467E-01;
    COFD[1087] = 1.06263545E-02;
    COFD[1088] = -1.34247866E+01;
    COFD[1089] = 3.48624238E+00;
    COFD[1090] = -2.41554467E-01;
    COFD[1091] = 1.06263545E-02;
    COFD[1092] = -1.46554748E+01;
    COFD[1093] = 3.83606243E+00;
    COFD[1094] = -2.86076532E-01;
    COFD[1095] = 1.25205829E-02;
    COFD[1096] = -1.95739570E+01;
    COFD[1097] = 5.61113230E+00;
    COFD[1098] = -4.90190187E-01;
    COFD[1099] = 2.03260675E-02;
    COFD[1100] = -1.46559141E+01;
    COFD[1101] = 3.83606243E+00;
    COFD[1102] = -2.86076532E-01;
    COFD[1103] = 1.25205829E-02;
    COFD[1104] = -1.46550083E+01;
    COFD[1105] = 3.83606243E+00;
    COFD[1106] = -2.86076532E-01;
    COFD[1107] = 1.25205829E-02;
    COFD[1108] = -1.57994893E+01;
    COFD[1109] = 4.22225052E+00;
    COFD[1110] = -3.35156428E-01;
    COFD[1111] = 1.46104855E-02;
    COFD[1112] = -1.97550088E+01;
    COFD[1113] = 5.56931926E+00;
    COFD[1114] = -4.89105511E-01;
    COFD[1115] = 2.04493129E-02;
    COFD[1116] = -1.43151174E+01;
    COFD[1117] = 3.68038508E+00;
    COFD[1118] = -2.65779346E-01;
    COFD[1119] = 1.16360771E-02;
    COFD[1120] = -1.76147026E+01;
    COFD[1121] = 4.86049500E+00;
    COFD[1122] = -4.12200578E-01;
    COFD[1123] = 1.77160971E-02;
    COFD[1124] = -1.72232223E+01;
    COFD[1125] = 4.69060745E+00;
    COFD[1126] = -3.92369888E-01;
    COFD[1127] = 1.69459661E-02;
    COFD[1128] = -1.79344949E+01;
    COFD[1129] = 4.91373893E+00;
    COFD[1130] = -4.18747629E-01;
    COFD[1131] = 1.79856610E-02;
    COFD[1132] = -1.94688688E+01;
    COFD[1133] = 5.43830787E+00;
    COFD[1134] = -4.75472880E-01;
    COFD[1135] = 1.99909996E-02;
    COFD[1136] = -1.83539686E+01;
    COFD[1137] = 4.98756925E+00;
    COFD[1138] = -4.27526072E-01;
    COFD[1139] = 1.83341755E-02;
    COFD[1140] = -2.01029332E+01;
    COFD[1141] = 5.51409582E+00;
    COFD[1142] = -4.83844807E-01;
    COFD[1143] = 2.02965712E-02;
    COFD[1144] = -2.03845035E+01;
    COFD[1145] = 5.56619223E+00;
    COFD[1146] = -4.88860312E-01;
    COFD[1147] = 2.04450332E-02;
    COFD[1148] = -2.02678158E+01;
    COFD[1149] = 5.55236751E+00;
    COFD[1150] = -4.87675676E-01;
    COFD[1151] = 2.04178098E-02;
    COFD[1152] = -2.05542181E+01;
    COFD[1153] = 5.59742840E+00;
    COFD[1154] = -4.86964048E-01;
    COFD[1155] = 2.01281667E-02;
    COFD[1156] = -1.42894441E+01;
    COFD[1157] = 3.67490723E+00;
    COFD[1158] = -2.65114792E-01;
    COFD[1159] = 1.16092671E-02;
    COFD[1160] = -2.16332288E+01;
    COFD[1161] = 5.40860960E+00;
    COFD[1162] = -4.73010822E-01;
    COFD[1163] = 1.99398149E-02;
    COFD[1164] = -1.74792112E+01;
    COFD[1165] = 4.29676909E+00;
    COFD[1166] = -3.44085306E-01;
    COFD[1167] = 1.49671135E-02;
    COFD[1168] = -1.50270339E+01;
    COFD[1169] = 3.46140064E+00;
    COFD[1170] = -2.38440092E-01;
    COFD[1171] = 1.04960087E-02;
    COFD[1172] = -1.09595712E+01;
    COFD[1173] = 2.30836460E+00;
    COFD[1174] = -8.76339315E-02;
    COFD[1175] = 3.90878445E-03;
    COFD[1176] = -1.76808721E+01;
    COFD[1177] = 4.24719726E+00;
    COFD[1178] = -3.38206061E-01;
    COFD[1179] = 1.47350654E-02;
    COFD[1180] = -1.87878849E+01;
    COFD[1181] = 4.61260432E+00;
    COFD[1182] = -3.82854484E-01;
    COFD[1183] = 1.65575163E-02;
    COFD[1184] = -1.87878849E+01;
    COFD[1185] = 4.61260432E+00;
    COFD[1186] = -3.82854484E-01;
    COFD[1187] = 1.65575163E-02;
    COFD[1188] = -1.97761100E+01;
    COFD[1189] = 4.93652913E+00;
    COFD[1190] = -4.21485300E-01;
    COFD[1191] = 1.80955487E-02;
    COFD[1192] = -2.01219251E+01;
    COFD[1193] = 5.01930029E+00;
    COFD[1194] = -4.31305027E-01;
    COFD[1195] = 1.84846322E-02;
    COFD[1196] = -1.34230272E+01;
    COFD[1197] = 3.48624238E+00;
    COFD[1198] = -2.41554467E-01;
    COFD[1199] = 1.06263545E-02;
    COFD[1200] = -1.32093628E+01;
    COFD[1201] = 2.90778936E+00;
    COFD[1202] = -1.67388544E-01;
    COFD[1203] = 7.45220609E-03;
    COFD[1204] = -1.32244035E+01;
    COFD[1205] = 2.90778936E+00;
    COFD[1206] = -1.67388544E-01;
    COFD[1207] = 7.45220609E-03;
    COFD[1208] = -1.43190389E+01;
    COFD[1209] = 3.17651319E+00;
    COFD[1210] = -2.02028974E-01;
    COFD[1211] = 8.94232502E-03;
    COFD[1212] = -1.94093572E+01;
    COFD[1213] = 5.16013126E+00;
    COFD[1214] = -4.46824543E-01;
    COFD[1215] = 1.90464887E-02;
    COFD[1216] = -1.43238998E+01;
    COFD[1217] = 3.17651319E+00;
    COFD[1218] = -2.02028974E-01;
    COFD[1219] = 8.94232502E-03;
    COFD[1220] = -1.43139231E+01;
    COFD[1221] = 3.17651319E+00;
    COFD[1222] = -2.02028974E-01;
    COFD[1223] = 8.94232502E-03;
    COFD[1224] = -1.50766130E+01;
    COFD[1225] = 3.47945612E+00;
    COFD[1226] = -2.40703722E-01;
    COFD[1227] = 1.05907441E-02;
    COFD[1228] = -1.94373127E+01;
    COFD[1229] = 5.02567894E+00;
    COFD[1230] = -4.32045169E-01;
    COFD[1231] = 1.85132214E-02;
    COFD[1232] = -1.40999008E+01;
    COFD[1233] = 3.08120012E+00;
    COFD[1234] = -1.89629903E-01;
    COFD[1235] = 8.40361952E-03;
    COFD[1236] = -1.70534856E+01;
    COFD[1237] = 4.14240922E+00;
    COFD[1238] = -3.25239774E-01;
    COFD[1239] = 1.41980687E-02;
    COFD[1240] = -1.65488358E+01;
    COFD[1241] = 3.95035840E+00;
    COFD[1242] = -3.00959418E-01;
    COFD[1243] = 1.31692593E-02;
    COFD[1244] = -1.72556499E+01;
    COFD[1245] = 4.17889917E+00;
    COFD[1246] = -3.29752510E-01;
    COFD[1247] = 1.43850275E-02;
    COFD[1248] = -1.90883268E+01;
    COFD[1249] = 4.84384483E+00;
    COFD[1250] = -4.10265575E-01;
    COFD[1251] = 1.76414287E-02;
    COFD[1252] = -1.76775033E+01;
    COFD[1253] = 4.24719726E+00;
    COFD[1254] = -3.38206061E-01;
    COFD[1255] = 1.47350654E-02;
    COFD[1256] = -1.97750001E+01;
    COFD[1257] = 4.93652913E+00;
    COFD[1258] = -4.21485300E-01;
    COFD[1259] = 1.80955487E-02;
    COFD[1260] = -2.01211076E+01;
    COFD[1261] = 5.01930029E+00;
    COFD[1262] = -4.31305027E-01;
    COFD[1263] = 1.84846322E-02;
    COFD[1264] = -1.99915387E+01;
    COFD[1265] = 4.99136452E+00;
    COFD[1266] = -4.27975210E-01;
    COFD[1267] = 1.83519204E-02;
    COFD[1268] = -2.06365262E+01;
    COFD[1269] = 5.20105041E+00;
    COFD[1270] = -4.51465908E-01;
    COFD[1271] = 1.92210125E-02;
    COFD[1272] = -1.40756935E+01;
    COFD[1273] = 3.07549274E+00;
    COFD[1274] = -1.88889344E-01;
    COFD[1275] = 8.37152866E-03;
    COFD[1276] = -2.16608259E+01;
    COFD[1277] = 5.40860960E+00;
    COFD[1278] = -4.73010822E-01;
    COFD[1279] = 1.99398149E-02;
    COFD[1280] = -1.74984476E+01;
    COFD[1281] = 4.29676909E+00;
    COFD[1282] = -3.44085306E-01;
    COFD[1283] = 1.49671135E-02;
    COFD[1284] = -1.50420953E+01;
    COFD[1285] = 3.46140064E+00;
    COFD[1286] = -2.38440092E-01;
    COFD[1287] = 1.04960087E-02;
    COFD[1288] = -1.09628982E+01;
    COFD[1289] = 2.30836460E+00;
    COFD[1290] = -8.76339315E-02;
    COFD[1291] = 3.90878445E-03;
    COFD[1292] = -1.77028170E+01;
    COFD[1293] = 4.24719726E+00;
    COFD[1294] = -3.38206061E-01;
    COFD[1295] = 1.47350654E-02;
    COFD[1296] = -1.88114917E+01;
    COFD[1297] = 4.61260432E+00;
    COFD[1298] = -3.82854484E-01;
    COFD[1299] = 1.65575163E-02;
    COFD[1300] = -1.88114917E+01;
    COFD[1301] = 4.61260432E+00;
    COFD[1302] = -3.82854484E-01;
    COFD[1303] = 1.65575163E-02;
    COFD[1304] = -1.98013312E+01;
    COFD[1305] = 4.93652913E+00;
    COFD[1306] = -4.21485300E-01;
    COFD[1307] = 1.80955487E-02;
    COFD[1308] = -2.01478340E+01;
    COFD[1309] = 5.01930029E+00;
    COFD[1310] = -4.31305027E-01;
    COFD[1311] = 1.84846322E-02;
    COFD[1312] = -1.34247866E+01;
    COFD[1313] = 3.48624238E+00;
    COFD[1314] = -2.41554467E-01;
    COFD[1315] = 1.06263545E-02;
    COFD[1316] = -1.32244035E+01;
    COFD[1317] = 2.90778936E+00;
    COFD[1318] = -1.67388544E-01;
    COFD[1319] = 7.45220609E-03;
    COFD[1320] = -1.32399106E+01;
    COFD[1321] = 2.90778936E+00;
    COFD[1322] = -1.67388544E-01;
    COFD[1323] = 7.45220609E-03;
    COFD[1324] = -1.43394069E+01;
    COFD[1325] = 3.17651319E+00;
    COFD[1326] = -2.02028974E-01;
    COFD[1327] = 8.94232502E-03;
    COFD[1328] = -1.94253036E+01;
    COFD[1329] = 5.16013126E+00;
    COFD[1330] = -4.46824543E-01;
    COFD[1331] = 1.90464887E-02;
    COFD[1332] = -1.43444709E+01;
    COFD[1333] = 3.17651319E+00;
    COFD[1334] = -2.02028974E-01;
    COFD[1335] = 8.94232502E-03;
    COFD[1336] = -1.43340796E+01;
    COFD[1337] = 3.17651319E+00;
    COFD[1338] = -2.02028974E-01;
    COFD[1339] = 8.94232502E-03;
    COFD[1340] = -1.50911794E+01;
    COFD[1341] = 3.47945612E+00;
    COFD[1342] = -2.40703722E-01;
    COFD[1343] = 1.05907441E-02;
    COFD[1344] = -1.94570287E+01;
    COFD[1345] = 5.02567894E+00;
    COFD[1346] = -4.32045169E-01;
    COFD[1347] = 1.85132214E-02;
    COFD[1348] = -1.41191261E+01;
    COFD[1349] = 3.08120012E+00;
    COFD[1350] = -1.89629903E-01;
    COFD[1351] = 8.40361952E-03;
    COFD[1352] = -1.70757047E+01;
    COFD[1353] = 4.14240922E+00;
    COFD[1354] = -3.25239774E-01;
    COFD[1355] = 1.41980687E-02;
    COFD[1356] = -1.65675362E+01;
    COFD[1357] = 3.95035840E+00;
    COFD[1358] = -3.00959418E-01;
    COFD[1359] = 1.31692593E-02;
    COFD[1360] = -1.72753760E+01;
    COFD[1361] = 4.17889917E+00;
    COFD[1362] = -3.29752510E-01;
    COFD[1363] = 1.43850275E-02;
    COFD[1364] = -1.91102652E+01;
    COFD[1365] = 4.84384483E+00;
    COFD[1366] = -4.10265575E-01;
    COFD[1367] = 1.76414287E-02;
    COFD[1368] = -1.76992976E+01;
    COFD[1369] = 4.24719726E+00;
    COFD[1370] = -3.38206061E-01;
    COFD[1371] = 1.47350654E-02;
    COFD[1372] = -1.98001639E+01;
    COFD[1373] = 4.93652913E+00;
    COFD[1374] = -4.21485300E-01;
    COFD[1375] = 1.80955487E-02;
    COFD[1376] = -2.01469731E+01;
    COFD[1377] = 5.01930029E+00;
    COFD[1378] = -4.31305027E-01;
    COFD[1379] = 1.84846322E-02;
    COFD[1380] = -2.00180417E+01;
    COFD[1381] = 4.99136452E+00;
    COFD[1382] = -4.27975210E-01;
    COFD[1383] = 1.83519204E-02;
    COFD[1384] = -2.06629641E+01;
    COFD[1385] = 5.20105041E+00;
    COFD[1386] = -4.51465908E-01;
    COFD[1387] = 1.92210125E-02;
    COFD[1388] = -1.40949196E+01;
    COFD[1389] = 3.07549274E+00;
    COFD[1390] = -1.88889344E-01;
    COFD[1391] = 8.37152866E-03;
    COFD[1392] = -2.28139639E+01;
    COFD[1393] = 5.61234608E+00;
    COFD[1394] = -4.91375180E-01;
    COFD[1395] = 2.04179832E-02;
    COFD[1396] = -1.89616623E+01;
    COFD[1397] = 4.68595732E+00;
    COFD[1398] = -3.91842840E-01;
    COFD[1399] = 1.69262542E-02;
    COFD[1400] = -1.62775714E+01;
    COFD[1401] = 3.79163564E+00;
    COFD[1402] = -2.80257365E-01;
    COFD[1403] = 1.22656902E-02;
    COFD[1404] = -1.18998012E+01;
    COFD[1405] = 2.57507000E+00;
    COFD[1406] = -1.24033737E-01;
    COFD[1407] = 5.56694959E-03;
    COFD[1408] = -1.91261963E+01;
    COFD[1409] = 4.61801405E+00;
    COFD[1410] = -3.83535652E-01;
    COFD[1411] = 1.65862513E-02;
    COFD[1412] = -2.02566224E+01;
    COFD[1413] = 4.97613338E+00;
    COFD[1414] = -4.26175206E-01;
    COFD[1415] = 1.82809270E-02;
    COFD[1416] = -2.02566224E+01;
    COFD[1417] = 4.97613338E+00;
    COFD[1418] = -4.26175206E-01;
    COFD[1419] = 1.82809270E-02;
    COFD[1420] = -2.12131123E+01;
    COFD[1421] = 5.26829491E+00;
    COFD[1422] = -4.59325230E-01;
    COFD[1423] = 1.95275380E-02;
    COFD[1424] = -2.15623852E+01;
    COFD[1425] = 5.34961360E+00;
    COFD[1426] = -4.68762622E-01;
    COFD[1427] = 1.98932200E-02;
    COFD[1428] = -1.46554748E+01;
    COFD[1429] = 3.83606243E+00;
    COFD[1430] = -2.86076532E-01;
    COFD[1431] = 1.25205829E-02;
    COFD[1432] = -1.43190389E+01;
    COFD[1433] = 3.17651319E+00;
    COFD[1434] = -2.02028974E-01;
    COFD[1435] = 8.94232502E-03;
    COFD[1436] = -1.43394069E+01;
    COFD[1437] = 3.17651319E+00;
    COFD[1438] = -2.02028974E-01;
    COFD[1439] = 8.94232502E-03;
    COFD[1440] = -1.55666415E+01;
    COFD[1441] = 3.48070094E+00;
    COFD[1442] = -2.40859499E-01;
    COFD[1443] = 1.05972514E-02;
    COFD[1444] = -2.06463744E+01;
    COFD[1445] = 5.41688482E+00;
    COFD[1446] = -4.73387188E-01;
    COFD[1447] = 1.99280175E-02;
    COFD[1448] = -1.55741053E+01;
    COFD[1449] = 3.48070094E+00;
    COFD[1450] = -2.40859499E-01;
    COFD[1451] = 1.05972514E-02;
    COFD[1452] = -1.55588279E+01;
    COFD[1453] = 3.48070094E+00;
    COFD[1454] = -2.40859499E-01;
    COFD[1455] = 1.05972514E-02;
    COFD[1456] = -1.63542394E+01;
    COFD[1457] = 3.82388595E+00;
    COFD[1458] = -2.84480724E-01;
    COFD[1459] = 1.24506311E-02;
    COFD[1460] = -2.08367725E+01;
    COFD[1461] = 5.35267674E+00;
    COFD[1462] = -4.69010505E-01;
    COFD[1463] = 1.98979152E-02;
    COFD[1464] = -1.52792891E+01;
    COFD[1465] = 3.36790500E+00;
    COFD[1466] = -2.26321740E-01;
    COFD[1467] = 9.97135055E-03;
    COFD[1468] = -1.84777607E+01;
    COFD[1469] = 4.49330851E+00;
    COFD[1470] = -3.68208715E-01;
    COFD[1471] = 1.59565402E-02;
    COFD[1472] = -1.78903913E+01;
    COFD[1473] = 4.29613154E+00;
    COFD[1474] = -3.44012526E-01;
    COFD[1475] = 1.49643715E-02;
    COFD[1476] = -1.86499071E+01;
    COFD[1477] = 4.53572533E+00;
    COFD[1478] = -3.73386925E-01;
    COFD[1479] = 1.61678881E-02;
    COFD[1480] = -2.05272328E+01;
    COFD[1481] = 5.18417470E+00;
    COFD[1482] = -4.49491573E-01;
    COFD[1483] = 1.91438508E-02;
    COFD[1484] = -1.91208314E+01;
    COFD[1485] = 4.61801405E+00;
    COFD[1486] = -3.83535652E-01;
    COFD[1487] = 1.65862513E-02;
    COFD[1488] = -2.12111746E+01;
    COFD[1489] = 5.26829491E+00;
    COFD[1490] = -4.59325230E-01;
    COFD[1491] = 1.95275380E-02;
    COFD[1492] = -2.15609288E+01;
    COFD[1493] = 5.34961360E+00;
    COFD[1494] = -4.68762622E-01;
    COFD[1495] = 1.98932200E-02;
    COFD[1496] = -2.14647019E+01;
    COFD[1497] = 5.33054865E+00;
    COFD[1498] = -4.66761560E-01;
    COFD[1499] = 1.98252926E-02;
    COFD[1500] = -2.19080065E+01;
    COFD[1501] = 5.44858060E+00;
    COFD[1502] = -4.76668834E-01;
    COFD[1503] = 2.00375810E-02;
    COFD[1504] = -1.52486273E+01;
    COFD[1505] = 3.35922578E+00;
    COFD[1506] = -2.25181399E-01;
    COFD[1507] = 9.92132878E-03;
    COFD[1508] = -1.51900988E+01;
    COFD[1509] = 1.91174168E+00;
    COFD[1510] = 8.81438912E-02;
    COFD[1511] = -8.55341166E-03;
    COFD[1512] = -2.08812333E+01;
    COFD[1513] = 5.08859217E+00;
    COFD[1514] = -3.90525428E-01;
    COFD[1515] = 1.47376395E-02;
    COFD[1516] = -2.14087397E+01;
    COFD[1517] = 5.57282008E+00;
    COFD[1518] = -4.76690890E-01;
    COFD[1519] = 1.94000719E-02;
    COFD[1520] = -1.71982995E+01;
    COFD[1521] = 4.63881404E+00;
    COFD[1522] = -3.86139633E-01;
    COFD[1523] = 1.66955081E-02;
    COFD[1524] = -2.13884087E+01;
    COFD[1525] = 5.17440955E+00;
    COFD[1526] = -4.04678430E-01;
    COFD[1527] = 1.54706350E-02;
    COFD[1528] = -2.03437836E+01;
    COFD[1529] = 4.57152878E+00;
    COFD[1530] = -3.08371263E-01;
    COFD[1531] = 1.05838559E-02;
    COFD[1532] = -2.03437836E+01;
    COFD[1533] = 4.57152878E+00;
    COFD[1534] = -3.08371263E-01;
    COFD[1535] = 1.05838559E-02;
    COFD[1536] = -1.74926607E+01;
    COFD[1537] = 3.15361887E+00;
    COFD[1538] = -9.38468185E-02;
    COFD[1539] = 1.12007553E-04;
    COFD[1540] = -1.83871985E+01;
    COFD[1541] = 3.54020956E+00;
    COFD[1542] = -1.52304845E-01;
    COFD[1543] = 2.96268912E-03;
    COFD[1544] = -1.95739570E+01;
    COFD[1545] = 5.61113230E+00;
    COFD[1546] = -4.90190187E-01;
    COFD[1547] = 2.03260675E-02;
    COFD[1548] = -1.94093572E+01;
    COFD[1549] = 5.16013126E+00;
    COFD[1550] = -4.46824543E-01;
    COFD[1551] = 1.90464887E-02;
    COFD[1552] = -1.94253036E+01;
    COFD[1553] = 5.16013126E+00;
    COFD[1554] = -4.46824543E-01;
    COFD[1555] = 1.90464887E-02;
    COFD[1556] = -2.06463744E+01;
    COFD[1557] = 5.41688482E+00;
    COFD[1558] = -4.73387188E-01;
    COFD[1559] = 1.99280175E-02;
    COFD[1560] = -1.19157919E+01;
    COFD[1561] = 9.28955130E-01;
    COFD[1562] = 2.42107090E-01;
    COFD[1563] = -1.59823963E-02;
    COFD[1564] = -2.06516336E+01;
    COFD[1565] = 5.41688482E+00;
    COFD[1566] = -4.73387188E-01;
    COFD[1567] = 1.99280175E-02;
    COFD[1568] = -2.12652533E+01;
    COFD[1569] = 5.59961818E+00;
    COFD[1570] = -4.91624858E-01;
    COFD[1571] = 2.05035550E-02;
    COFD[1572] = -2.12831323E+01;
    COFD[1573] = 5.61184117E+00;
    COFD[1574] = -4.90532156E-01;
    COFD[1575] = 2.03507922E-02;
    COFD[1576] = -1.77563250E+01;
    COFD[1577] = 3.57475686E+00;
    COFD[1578] = -1.56396297E-01;
    COFD[1579] = 3.12157721E-03;
    COFD[1580] = -2.11388331E+01;
    COFD[1581] = 5.55529675E+00;
    COFD[1582] = -4.87942518E-01;
    COFD[1583] = 2.04249054E-02;
    COFD[1584] = -2.07653719E+01;
    COFD[1585] = 5.01092022E+00;
    COFD[1586] = -3.77985635E-01;
    COFD[1587] = 1.40968645E-02;
    COFD[1588] = -2.15095980E+01;
    COFD[1589] = 5.46737673E+00;
    COFD[1590] = -4.55696085E-01;
    COFD[1591] = 1.81982625E-02;
    COFD[1592] = -2.12661865E+01;
    COFD[1593] = 5.24930667E+00;
    COFD[1594] = -4.17435088E-01;
    COFD[1595] = 1.61434424E-02;
    COFD[1596] = -1.87383952E+01;
    COFD[1597] = 3.96926341E+00;
    COFD[1598] = -2.16412264E-01;
    COFD[1599] = 6.06012078E-03;
    COFD[1600] = -2.13847439E+01;
    COFD[1601] = 5.17440955E+00;
    COFD[1602] = -4.04678430E-01;
    COFD[1603] = 1.54706350E-02;
    COFD[1604] = -1.74914373E+01;
    COFD[1605] = 3.15361887E+00;
    COFD[1606] = -9.38468185E-02;
    COFD[1607] = 1.12007553E-04;
    COFD[1608] = -1.83862949E+01;
    COFD[1609] = 3.54020956E+00;
    COFD[1610] = -1.52304845E-01;
    COFD[1611] = 2.96268912E-03;
    COFD[1612] = -1.84877055E+01;
    COFD[1613] = 3.61246265E+00;
    COFD[1614] = -1.63179700E-01;
    COFD[1615] = 3.49175937E-03;
    COFD[1616] = -1.54850731E+01;
    COFD[1617] = 2.22760630E+00;
    COFD[1618] = 3.88962143E-02;
    COFD[1619] = -6.06812513E-03;
    COFD[1620] = -2.10643259E+01;
    COFD[1621] = 5.53614847E+00;
    COFD[1622] = -4.86046736E-01;
    COFD[1623] = 2.03659188E-02;
    COFD[1624] = -2.28263210E+01;
    COFD[1625] = 5.61234608E+00;
    COFD[1626] = -4.91375180E-01;
    COFD[1627] = 2.04179832E-02;
    COFD[1628] = -1.89685165E+01;
    COFD[1629] = 4.68595732E+00;
    COFD[1630] = -3.91842840E-01;
    COFD[1631] = 1.69262542E-02;
    COFD[1632] = -1.62824412E+01;
    COFD[1633] = 3.79163564E+00;
    COFD[1634] = -2.80257365E-01;
    COFD[1635] = 1.22656902E-02;
    COFD[1636] = -1.19006548E+01;
    COFD[1637] = 2.57507000E+00;
    COFD[1638] = -1.24033737E-01;
    COFD[1639] = 5.56694959E-03;
    COFD[1640] = -1.91345696E+01;
    COFD[1641] = 4.61801405E+00;
    COFD[1642] = -3.83535652E-01;
    COFD[1643] = 1.65862513E-02;
    COFD[1644] = -2.02660394E+01;
    COFD[1645] = 4.97613338E+00;
    COFD[1646] = -4.26175206E-01;
    COFD[1647] = 1.82809270E-02;
    COFD[1648] = -2.02660394E+01;
    COFD[1649] = 4.97613338E+00;
    COFD[1650] = -4.26175206E-01;
    COFD[1651] = 1.82809270E-02;
    COFD[1652] = -2.12236380E+01;
    COFD[1653] = 5.26829491E+00;
    COFD[1654] = -4.59325230E-01;
    COFD[1655] = 1.95275380E-02;
    COFD[1656] = -2.15734148E+01;
    COFD[1657] = 5.34961360E+00;
    COFD[1658] = -4.68762622E-01;
    COFD[1659] = 1.98932200E-02;
    COFD[1660] = -1.46559141E+01;
    COFD[1661] = 3.83606243E+00;
    COFD[1662] = -2.86076532E-01;
    COFD[1663] = 1.25205829E-02;
    COFD[1664] = -1.43238998E+01;
    COFD[1665] = 3.17651319E+00;
    COFD[1666] = -2.02028974E-01;
    COFD[1667] = 8.94232502E-03;
    COFD[1668] = -1.43444709E+01;
    COFD[1669] = 3.17651319E+00;
    COFD[1670] = -2.02028974E-01;
    COFD[1671] = 8.94232502E-03;
    COFD[1672] = -1.55741053E+01;
    COFD[1673] = 3.48070094E+00;
    COFD[1674] = -2.40859499E-01;
    COFD[1675] = 1.05972514E-02;
    COFD[1676] = -2.06516336E+01;
    COFD[1677] = 5.41688482E+00;
    COFD[1678] = -4.73387188E-01;
    COFD[1679] = 1.99280175E-02;
    COFD[1680] = -1.55816822E+01;
    COFD[1681] = 3.48070094E+00;
    COFD[1682] = -2.40859499E-01;
    COFD[1683] = 1.05972514E-02;
    COFD[1684] = -1.55661750E+01;
    COFD[1685] = 3.48070094E+00;
    COFD[1686] = -2.40859499E-01;
    COFD[1687] = 1.05972514E-02;
    COFD[1688] = -1.63588981E+01;
    COFD[1689] = 3.82388595E+00;
    COFD[1690] = -2.84480724E-01;
    COFD[1691] = 1.24506311E-02;
    COFD[1692] = -2.08438809E+01;
    COFD[1693] = 5.35267674E+00;
    COFD[1694] = -4.69010505E-01;
    COFD[1695] = 1.98979152E-02;
    COFD[1696] = -1.52861376E+01;
    COFD[1697] = 3.36790500E+00;
    COFD[1698] = -2.26321740E-01;
    COFD[1699] = 9.97135055E-03;
    COFD[1700] = -1.84863000E+01;
    COFD[1701] = 4.49330851E+00;
    COFD[1702] = -3.68208715E-01;
    COFD[1703] = 1.59565402E-02;
    COFD[1704] = -1.78969684E+01;
    COFD[1705] = 4.29613154E+00;
    COFD[1706] = -3.44012526E-01;
    COFD[1707] = 1.49643715E-02;
    COFD[1708] = -1.86570209E+01;
    COFD[1709] = 4.53572533E+00;
    COFD[1710] = -3.73386925E-01;
    COFD[1711] = 1.61678881E-02;
    COFD[1712] = -2.05356023E+01;
    COFD[1713] = 5.18417470E+00;
    COFD[1714] = -4.49491573E-01;
    COFD[1715] = 1.91438508E-02;
    COFD[1716] = -1.91291147E+01;
    COFD[1717] = 4.61801405E+00;
    COFD[1718] = -3.83535652E-01;
    COFD[1719] = 1.65862513E-02;
    COFD[1720] = -2.12216591E+01;
    COFD[1721] = 5.26829491E+00;
    COFD[1722] = -4.59325230E-01;
    COFD[1723] = 1.95275380E-02;
    COFD[1724] = -2.15719260E+01;
    COFD[1725] = 5.34961360E+00;
    COFD[1726] = -4.68762622E-01;
    COFD[1727] = 1.98932200E-02;
    COFD[1728] = -2.14761835E+01;
    COFD[1729] = 5.33054865E+00;
    COFD[1730] = -4.66761560E-01;
    COFD[1731] = 1.98252926E-02;
    COFD[1732] = -2.19194378E+01;
    COFD[1733] = 5.44858060E+00;
    COFD[1734] = -4.76668834E-01;
    COFD[1735] = 2.00375810E-02;
    COFD[1736] = -1.52554761E+01;
    COFD[1737] = 3.35922578E+00;
    COFD[1738] = -2.25181399E-01;
    COFD[1739] = 9.92132878E-03;
    COFD[1740] = -2.28011547E+01;
    COFD[1741] = 5.61234608E+00;
    COFD[1742] = -4.91375180E-01;
    COFD[1743] = 2.04179832E-02;
    COFD[1744] = -1.89544778E+01;
    COFD[1745] = 4.68595732E+00;
    COFD[1746] = -3.91842840E-01;
    COFD[1747] = 1.69262542E-02;
    COFD[1748] = -1.62724462E+01;
    COFD[1749] = 3.79163564E+00;
    COFD[1750] = -2.80257365E-01;
    COFD[1751] = 1.22656902E-02;
    COFD[1752] = -1.18988955E+01;
    COFD[1753] = 2.57507000E+00;
    COFD[1754] = -1.24033737E-01;
    COFD[1755] = 5.56694959E-03;
    COFD[1756] = -1.91174465E+01;
    COFD[1757] = 4.61801405E+00;
    COFD[1758] = -3.83535652E-01;
    COFD[1759] = 1.65862513E-02;
    COFD[1760] = -2.02468029E+01;
    COFD[1761] = 4.97613338E+00;
    COFD[1762] = -4.26175206E-01;
    COFD[1763] = 1.82809270E-02;
    COFD[1764] = -2.02468029E+01;
    COFD[1765] = 4.97613338E+00;
    COFD[1766] = -4.26175206E-01;
    COFD[1767] = 1.82809270E-02;
    COFD[1768] = -2.12021613E+01;
    COFD[1769] = 5.26829491E+00;
    COFD[1770] = -4.59325230E-01;
    COFD[1771] = 1.95275380E-02;
    COFD[1772] = -2.15536156E+01;
    COFD[1773] = 5.35030481E+00;
    COFD[1774] = -4.68818560E-01;
    COFD[1775] = 1.98942796E-02;
    COFD[1776] = -1.46550083E+01;
    COFD[1777] = 3.83606243E+00;
    COFD[1778] = -2.86076532E-01;
    COFD[1779] = 1.25205829E-02;
    COFD[1780] = -1.43139231E+01;
    COFD[1781] = 3.17651319E+00;
    COFD[1782] = -2.02028974E-01;
    COFD[1783] = 8.94232502E-03;
    COFD[1784] = -1.43340796E+01;
    COFD[1785] = 3.17651319E+00;
    COFD[1786] = -2.02028974E-01;
    COFD[1787] = 8.94232502E-03;
    COFD[1788] = -1.55588279E+01;
    COFD[1789] = 3.48070094E+00;
    COFD[1790] = -2.40859499E-01;
    COFD[1791] = 1.05972514E-02;
    COFD[1792] = -2.12652533E+01;
    COFD[1793] = 5.59961818E+00;
    COFD[1794] = -4.91624858E-01;
    COFD[1795] = 2.05035550E-02;
    COFD[1796] = -1.55661750E+01;
    COFD[1797] = 3.48070094E+00;
    COFD[1798] = -2.40859499E-01;
    COFD[1799] = 1.05972514E-02;
    COFD[1800] = -1.55511344E+01;
    COFD[1801] = 3.48070094E+00;
    COFD[1802] = -2.40859499E-01;
    COFD[1803] = 1.05972514E-02;
    COFD[1804] = -1.63493345E+01;
    COFD[1805] = 3.82388595E+00;
    COFD[1806] = -2.84480724E-01;
    COFD[1807] = 1.24506311E-02;
    COFD[1808] = -2.08293255E+01;
    COFD[1809] = 5.35267674E+00;
    COFD[1810] = -4.69010505E-01;
    COFD[1811] = 1.98979152E-02;
    COFD[1812] = -1.52721107E+01;
    COFD[1813] = 3.36790500E+00;
    COFD[1814] = -2.26321740E-01;
    COFD[1815] = 9.97135055E-03;
    COFD[1816] = -1.84688406E+01;
    COFD[1817] = 4.49330851E+00;
    COFD[1818] = -3.68208715E-01;
    COFD[1819] = 1.59565402E-02;
    COFD[1820] = -1.78834935E+01;
    COFD[1821] = 4.29613154E+00;
    COFD[1822] = -3.44012526E-01;
    COFD[1823] = 1.49643715E-02;
    COFD[1824] = -1.86424545E+01;
    COFD[1825] = 4.53572533E+00;
    COFD[1826] = -3.73386925E-01;
    COFD[1827] = 1.61678881E-02;
    COFD[1828] = -2.05184870E+01;
    COFD[1829] = 5.18417470E+00;
    COFD[1830] = -4.49491573E-01;
    COFD[1831] = 1.91438508E-02;
    COFD[1832] = -1.91121742E+01;
    COFD[1833] = 4.61801405E+00;
    COFD[1834] = -3.83535652E-01;
    COFD[1835] = 1.65862513E-02;
    COFD[1836] = -2.12002655E+01;
    COFD[1837] = 5.26829491E+00;
    COFD[1838] = -4.59325230E-01;
    COFD[1839] = 1.95275380E-02;
    COFD[1840] = -2.15521921E+01;
    COFD[1841] = 5.35030481E+00;
    COFD[1842] = -4.68818560E-01;
    COFD[1843] = 1.98942796E-02;
    COFD[1844] = -2.14573159E+01;
    COFD[1845] = 5.33204143E+00;
    COFD[1846] = -4.66932795E-01;
    COFD[1847] = 1.98318389E-02;
    COFD[1848] = -2.20898696E+01;
    COFD[1849] = 5.51001635E+00;
    COFD[1850] = -4.83405862E-01;
    COFD[1851] = 2.02810787E-02;
    COFD[1852] = -1.52414485E+01;
    COFD[1853] = 3.35922578E+00;
    COFD[1854] = -2.25181399E-01;
    COFD[1855] = 9.92132878E-03;
    COFD[1856] = -2.26092158E+01;
    COFD[1857] = 5.53027841E+00;
    COFD[1858] = -4.68614364E-01;
    COFD[1859] = 1.89477473E-02;
    COFD[1860] = -1.98646734E+01;
    COFD[1861] = 5.04367502E+00;
    COFD[1862] = -4.34153325E-01;
    COFD[1863] = 1.85956055E-02;
    COFD[1864] = -1.72738845E+01;
    COFD[1865] = 4.19029808E+00;
    COFD[1866] = -3.31177076E-01;
    COFD[1867] = 1.44446234E-02;
    COFD[1868] = -1.25141260E+01;
    COFD[1869] = 2.77873601E+00;
    COFD[1870] = -1.50637360E-01;
    COFD[1871] = 6.72684281E-03;
    COFD[1872] = -1.99835686E+01;
    COFD[1873] = 4.97875278E+00;
    COFD[1874] = -4.26485475E-01;
    COFD[1875] = 1.82931933E-02;
    COFD[1876] = -2.10572534E+01;
    COFD[1877] = 5.31360223E+00;
    COFD[1878] = -4.64787000E-01;
    COFD[1879] = 1.97483720E-02;
    COFD[1880] = -2.10572534E+01;
    COFD[1881] = 5.31360223E+00;
    COFD[1882] = -4.64787000E-01;
    COFD[1883] = 1.97483720E-02;
    COFD[1884] = -2.17763744E+01;
    COFD[1885] = 5.50642944E+00;
    COFD[1886] = -4.83016239E-01;
    COFD[1887] = 2.02671529E-02;
    COFD[1888] = -2.20578790E+01;
    COFD[1889] = 5.56193702E+00;
    COFD[1890] = -4.88512557E-01;
    COFD[1891] = 2.04380728E-02;
    COFD[1892] = -1.57994893E+01;
    COFD[1893] = 4.22225052E+00;
    COFD[1894] = -3.35156428E-01;
    COFD[1895] = 1.46104855E-02;
    COFD[1896] = -1.50766130E+01;
    COFD[1897] = 3.47945612E+00;
    COFD[1898] = -2.40703722E-01;
    COFD[1899] = 1.05907441E-02;
    COFD[1900] = -1.50911794E+01;
    COFD[1901] = 3.47945612E+00;
    COFD[1902] = -2.40703722E-01;
    COFD[1903] = 1.05907441E-02;
    COFD[1904] = -1.63542394E+01;
    COFD[1905] = 3.82388595E+00;
    COFD[1906] = -2.84480724E-01;
    COFD[1907] = 1.24506311E-02;
    COFD[1908] = -2.12831323E+01;
    COFD[1909] = 5.61184117E+00;
    COFD[1910] = -4.90532156E-01;
    COFD[1911] = 2.03507922E-02;
    COFD[1912] = -1.63588981E+01;
    COFD[1913] = 3.82388595E+00;
    COFD[1914] = -2.84480724E-01;
    COFD[1915] = 1.24506311E-02;
    COFD[1916] = -1.63493345E+01;
    COFD[1917] = 3.82388595E+00;
    COFD[1918] = -2.84480724E-01;
    COFD[1919] = 1.24506311E-02;
    COFD[1920] = -1.73374529E+01;
    COFD[1921] = 4.21416723E+00;
    COFD[1922] = -3.34163932E-01;
    COFD[1923] = 1.45697432E-02;
    COFD[1924] = -2.14449559E+01;
    COFD[1925] = 5.56531152E+00;
    COFD[1926] = -4.88789821E-01;
    COFD[1927] = 2.04437116E-02;
    COFD[1928] = -1.59863030E+01;
    COFD[1929] = 3.67388294E+00;
    COFD[1930] = -2.64990709E-01;
    COFD[1931] = 1.16042706E-02;
    COFD[1932] = -1.93276434E+01;
    COFD[1933] = 4.85015581E+00;
    COFD[1934] = -4.10945109E-01;
    COFD[1935] = 1.76651398E-02;
    COFD[1936] = -1.88463816E+01;
    COFD[1937] = 4.68393046E+00;
    COFD[1938] = -3.91610863E-01;
    COFD[1939] = 1.69174645E-02;
    COFD[1940] = -1.95552142E+01;
    COFD[1941] = 4.90255048E+00;
    COFD[1942] = -4.17368501E-01;
    COFD[1943] = 1.79287358E-02;
    COFD[1944] = -2.11606963E+01;
    COFD[1945] = 5.42846112E+00;
    COFD[1946] = -4.74321870E-01;
    COFD[1947] = 1.99459749E-02;
    COFD[1948] = -1.99803490E+01;
    COFD[1949] = 4.97875278E+00;
    COFD[1950] = -4.26485475E-01;
    COFD[1951] = 1.82931933E-02;
    COFD[1952] = -2.17753205E+01;
    COFD[1953] = 5.50642944E+00;
    COFD[1954] = -4.83016239E-01;
    COFD[1955] = 2.02671529E-02;
    COFD[1956] = -2.20571038E+01;
    COFD[1957] = 5.56193702E+00;
    COFD[1958] = -4.88512557E-01;
    COFD[1959] = 2.04380728E-02;
    COFD[1960] = -2.19586553E+01;
    COFD[1961] = 5.54610272E+00;
    COFD[1962] = -4.87043162E-01;
    COFD[1963] = 2.03974796E-02;
    COFD[1964] = -2.22630459E+01;
    COFD[1965] = 5.60031657E+00;
    COFD[1966] = -4.87628580E-01;
    COFD[1967] = 2.01686065E-02;
    COFD[1968] = -1.59633387E+01;
    COFD[1969] = 3.66853818E+00;
    COFD[1970] = -2.64346221E-01;
    COFD[1971] = 1.15784613E-02;
    COFD[1972] = -1.68489747E+01;
    COFD[1973] = 2.48068244E+00;
    COFD[1974] = 5.18653536E-03;
    COFD[1975] = -4.62524223E-03;
    COFD[1976] = -2.16379567E+01;
    COFD[1977] = 5.29019717E+00;
    COFD[1978] = -4.24502606E-01;
    COFD[1979] = 1.65197343E-02;
    COFD[1980] = -2.14082453E+01;
    COFD[1981] = 5.55346617E+00;
    COFD[1982] = -4.87783156E-01;
    COFD[1983] = 2.04210886E-02;
    COFD[1984] = -1.60528285E+01;
    COFD[1985] = 4.11188603E+00;
    COFD[1986] = -3.21540884E-01;
    COFD[1987] = 1.40482564E-02;
    COFD[1988] = -2.20998738E+01;
    COFD[1989] = 5.36053938E+00;
    COFD[1990] = -4.36434519E-01;
    COFD[1991] = 1.71484255E-02;
    COFD[1992] = -2.13352637E+01;
    COFD[1993] = 4.87252053E+00;
    COFD[1994] = -3.56127804E-01;
    COFD[1995] = 1.29948788E-02;
    COFD[1996] = -2.13352637E+01;
    COFD[1997] = 4.87252053E+00;
    COFD[1998] = -3.56127804E-01;
    COFD[1999] = 1.29948788E-02;
    COFD[2000] = -2.01549977E+01;
    COFD[2001] = 4.23025405E+00;
    COFD[2002] = -2.55784321E-01;
    COFD[2003] = 7.98147639E-03;
    COFD[2004] = -1.98066211E+01;
    COFD[2005] = 4.01137806E+00;
    COFD[2006] = -2.22734310E-01;
    COFD[2007] = 6.36741567E-03;
    COFD[2008] = -1.97550088E+01;
    COFD[2009] = 5.56931926E+00;
    COFD[2010] = -4.89105511E-01;
    COFD[2011] = 2.04493129E-02;
    COFD[2012] = -1.94373127E+01;
    COFD[2013] = 5.02567894E+00;
    COFD[2014] = -4.32045169E-01;
    COFD[2015] = 1.85132214E-02;
    COFD[2016] = -1.94570287E+01;
    COFD[2017] = 5.02567894E+00;
    COFD[2018] = -4.32045169E-01;
    COFD[2019] = 1.85132214E-02;
    COFD[2020] = -2.08367725E+01;
    COFD[2021] = 5.35267674E+00;
    COFD[2022] = -4.69010505E-01;
    COFD[2023] = 1.98979152E-02;
    COFD[2024] = -1.77563250E+01;
    COFD[2025] = 3.57475686E+00;
    COFD[2026] = -1.56396297E-01;
    COFD[2027] = 3.12157721E-03;
    COFD[2028] = -2.08438809E+01;
    COFD[2029] = 5.35267674E+00;
    COFD[2030] = -4.69010505E-01;
    COFD[2031] = 1.98979152E-02;
    COFD[2032] = -2.08293255E+01;
    COFD[2033] = 5.35267674E+00;
    COFD[2034] = -4.69010505E-01;
    COFD[2035] = 1.98979152E-02;
    COFD[2036] = -2.14449559E+01;
    COFD[2037] = 5.56531152E+00;
    COFD[2038] = -4.88789821E-01;
    COFD[2039] = 2.04437116E-02;
    COFD[2040] = -1.90499441E+01;
    COFD[2041] = 3.99221757E+00;
    COFD[2042] = -2.19854880E-01;
    COFD[2043] = 6.22736279E-03;
    COFD[2044] = -2.05128705E+01;
    COFD[2045] = 5.23843909E+00;
    COFD[2046] = -4.55815614E-01;
    COFD[2047] = 1.93898040E-02;
    COFD[2048] = -2.19317743E+01;
    COFD[2049] = 5.45216133E+00;
    COFD[2050] = -4.52916925E-01;
    COFD[2051] = 1.80456400E-02;
    COFD[2052] = -2.20036369E+01;
    COFD[2053] = 5.55935694E+00;
    COFD[2054] = -4.74154740E-01;
    COFD[2055] = 1.92584304E-02;
    COFD[2056] = -2.19399793E+01;
    COFD[2057] = 5.41841631E+00;
    COFD[2058] = -4.46818971E-01;
    COFD[2059] = 1.77127652E-02;
    COFD[2060] = -2.01015340E+01;
    COFD[2061] = 4.41511629E+00;
    COFD[2062] = -2.84086963E-01;
    COFD[2063] = 9.37586971E-03;
    COFD[2064] = -2.20947902E+01;
    COFD[2065] = 5.36053938E+00;
    COFD[2066] = -4.36434519E-01;
    COFD[2067] = 1.71484255E-02;
    COFD[2068] = -2.01531862E+01;
    COFD[2069] = 4.23025405E+00;
    COFD[2070] = -2.55784321E-01;
    COFD[2071] = 7.98147639E-03;
    COFD[2072] = -1.98052638E+01;
    COFD[2073] = 4.01137806E+00;
    COFD[2074] = -2.22734310E-01;
    COFD[2075] = 6.36741567E-03;
    COFD[2076] = -1.99254109E+01;
    COFD[2077] = 4.08598993E+00;
    COFD[2078] = -2.33972521E-01;
    COFD[2079] = 6.91534316E-03;
    COFD[2080] = -1.85959819E+01;
    COFD[2081] = 3.43982423E+00;
    COFD[2082] = -1.36303972E-01;
    COFD[2083] = 2.15194160E-03;
    COFD[2084] = -2.04833713E+01;
    COFD[2085] = 5.23112374E+00;
    COFD[2086] = -4.54967682E-01;
    COFD[2087] = 1.93570423E-02;
    COFD[2088] = -2.26357059E+01;
    COFD[2089] = 5.58177744E+00;
    COFD[2090] = -4.90187864E-01;
    COFD[2091] = 2.04753987E-02;
    COFD[2092] = -1.86157761E+01;
    COFD[2093] = 4.55689508E+00;
    COFD[2094] = -3.75937921E-01;
    COFD[2095] = 1.62703488E-02;
    COFD[2096] = -1.59525102E+01;
    COFD[2097] = 3.66023858E+00;
    COFD[2098] = -2.63401043E-01;
    COFD[2099] = 1.15432000E-02;
    COFD[2100] = -1.17159737E+01;
    COFD[2101] = 2.48123210E+00;
    COFD[2102] = -1.11322604E-01;
    COFD[2103] = 4.99282389E-03;
    COFD[2104] = -1.87733838E+01;
    COFD[2105] = 4.49191492E+00;
    COFD[2106] = -3.68041771E-01;
    COFD[2107] = 1.59498676E-02;
    COFD[2108] = -1.98832231E+01;
    COFD[2109] = 4.84731557E+00;
    COFD[2110] = -4.10638352E-01;
    COFD[2111] = 1.76543886E-02;
    COFD[2112] = -1.98832231E+01;
    COFD[2113] = 4.84731557E+00;
    COFD[2114] = -4.10638352E-01;
    COFD[2115] = 1.76543886E-02;
    COFD[2116] = -2.08713448E+01;
    COFD[2117] = 5.15826153E+00;
    COFD[2118] = -4.46699782E-01;
    COFD[2119] = 1.90459433E-02;
    COFD[2120] = -2.12018521E+01;
    COFD[2121] = 5.23328793E+00;
    COFD[2122] = -4.55218686E-01;
    COFD[2123] = 1.93667472E-02;
    COFD[2124] = -1.43151174E+01;
    COFD[2125] = 3.68038508E+00;
    COFD[2126] = -2.65779346E-01;
    COFD[2127] = 1.16360771E-02;
    COFD[2128] = -1.40999008E+01;
    COFD[2129] = 3.08120012E+00;
    COFD[2130] = -1.89629903E-01;
    COFD[2131] = 8.40361952E-03;
    COFD[2132] = -1.41191261E+01;
    COFD[2133] = 3.08120012E+00;
    COFD[2134] = -1.89629903E-01;
    COFD[2135] = 8.40361952E-03;
    COFD[2136] = -1.52792891E+01;
    COFD[2137] = 3.36790500E+00;
    COFD[2138] = -2.26321740E-01;
    COFD[2139] = 9.97135055E-03;
    COFD[2140] = -2.11388331E+01;
    COFD[2141] = 5.55529675E+00;
    COFD[2142] = -4.87942518E-01;
    COFD[2143] = 2.04249054E-02;
    COFD[2144] = -1.52861376E+01;
    COFD[2145] = 3.36790500E+00;
    COFD[2146] = -2.26321740E-01;
    COFD[2147] = 9.97135055E-03;
    COFD[2148] = -1.52721107E+01;
    COFD[2149] = 3.36790500E+00;
    COFD[2150] = -2.26321740E-01;
    COFD[2151] = 9.97135055E-03;
    COFD[2152] = -1.59863030E+01;
    COFD[2153] = 3.67388294E+00;
    COFD[2154] = -2.64990709E-01;
    COFD[2155] = 1.16042706E-02;
    COFD[2156] = -2.05128705E+01;
    COFD[2157] = 5.23843909E+00;
    COFD[2158] = -4.55815614E-01;
    COFD[2159] = 1.93898040E-02;
    COFD[2160] = -1.50233475E+01;
    COFD[2161] = 3.26660767E+00;
    COFD[2162] = -2.13287177E-01;
    COFD[2163] = 9.41137857E-03;
    COFD[2164] = -1.81735763E+01;
    COFD[2165] = 4.38391495E+00;
    COFD[2166] = -3.54941287E-01;
    COFD[2167] = 1.54195107E-02;
    COFD[2168] = -1.76285640E+01;
    COFD[2169] = 4.19935698E+00;
    COFD[2170] = -3.32310212E-01;
    COFD[2171] = 1.44920670E-02;
    COFD[2172] = -1.83538377E+01;
    COFD[2173] = 4.42828044E+00;
    COFD[2174] = -3.60417833E-01;
    COFD[2175] = 1.56455103E-02;
    COFD[2176] = -2.02922701E+01;
    COFD[2177] = 5.11106992E+00;
    COFD[2178] = -4.42047129E-01;
    COFD[2179] = 1.89042990E-02;
    COFD[2180] = -1.87685041E+01;
    COFD[2181] = 4.49191492E+00;
    COFD[2182] = -3.68041771E-01;
    COFD[2183] = 1.59498676E-02;
    COFD[2184] = -2.08696226E+01;
    COFD[2185] = 5.15826153E+00;
    COFD[2186] = -4.46699782E-01;
    COFD[2187] = 1.90459433E-02;
    COFD[2188] = -2.12005646E+01;
    COFD[2189] = 5.23328793E+00;
    COFD[2190] = -4.55218686E-01;
    COFD[2191] = 1.93667472E-02;
    COFD[2192] = -2.10856216E+01;
    COFD[2193] = 5.20655472E+00;
    COFD[2194] = -4.52109091E-01;
    COFD[2195] = 1.92461018E-02;
    COFD[2196] = -2.17920541E+01;
    COFD[2197] = 5.41750012E+00;
    COFD[2198] = -4.73406713E-01;
    COFD[2199] = 1.99264525E-02;
    COFD[2200] = -1.50031687E+01;
    COFD[2201] = 3.26223357E+00;
    COFD[2202] = -2.12746642E-01;
    COFD[2203] = 9.38912883E-03;
    COFD[2204] = -2.18641291E+01;
    COFD[2205] = 4.80496160E+00;
    COFD[2206] = -3.45476597E-01;
    COFD[2207] = 1.24586877E-02;
    COFD[2208] = -2.16802612E+01;
    COFD[2209] = 5.52918296E+00;
    COFD[2210] = -4.85360709E-01;
    COFD[2211] = 2.03448006E-02;
    COFD[2212] = -1.92867554E+01;
    COFD[2213] = 4.83375900E+00;
    COFD[2214] = -4.09146560E-01;
    COFD[2215] = 1.76006599E-02;
    COFD[2216] = -1.37794315E+01;
    COFD[2217] = 3.23973858E+00;
    COFD[2218] = -2.09989036E-01;
    COFD[2219] = 9.27667906E-03;
    COFD[2220] = -2.18653077E+01;
    COFD[2221] = 5.47368915E+00;
    COFD[2222] = -4.79424291E-01;
    COFD[2223] = 2.01372920E-02;
    COFD[2224] = -2.25018756E+01;
    COFD[2225] = 5.59178974E+00;
    COFD[2226] = -4.85668031E-01;
    COFD[2227] = 2.00491907E-02;
    COFD[2228] = -2.25018756E+01;
    COFD[2229] = 5.59178974E+00;
    COFD[2230] = -4.85668031E-01;
    COFD[2231] = 2.00491907E-02;
    COFD[2232] = -2.26266184E+01;
    COFD[2233] = 5.50195433E+00;
    COFD[2234] = -4.62690815E-01;
    COFD[2235] = 1.86010137E-02;
    COFD[2236] = -2.26941847E+01;
    COFD[2237] = 5.45651759E+00;
    COFD[2238] = -4.53709578E-01;
    COFD[2239] = 1.80890894E-02;
    COFD[2240] = -1.76147026E+01;
    COFD[2241] = 4.86049500E+00;
    COFD[2242] = -4.12200578E-01;
    COFD[2243] = 1.77160971E-02;
    COFD[2244] = -1.70534856E+01;
    COFD[2245] = 4.14240922E+00;
    COFD[2246] = -3.25239774E-01;
    COFD[2247] = 1.41980687E-02;
    COFD[2248] = -1.70757047E+01;
    COFD[2249] = 4.14240922E+00;
    COFD[2250] = -3.25239774E-01;
    COFD[2251] = 1.41980687E-02;
    COFD[2252] = -1.84777607E+01;
    COFD[2253] = 4.49330851E+00;
    COFD[2254] = -3.68208715E-01;
    COFD[2255] = 1.59565402E-02;
    COFD[2256] = -2.07653719E+01;
    COFD[2257] = 5.01092022E+00;
    COFD[2258] = -3.77985635E-01;
    COFD[2259] = 1.40968645E-02;
    COFD[2260] = -1.84863000E+01;
    COFD[2261] = 4.49330851E+00;
    COFD[2262] = -3.68208715E-01;
    COFD[2263] = 1.59565402E-02;
    COFD[2264] = -1.84688406E+01;
    COFD[2265] = 4.49330851E+00;
    COFD[2266] = -3.68208715E-01;
    COFD[2267] = 1.59565402E-02;
    COFD[2268] = -1.93276434E+01;
    COFD[2269] = 4.85015581E+00;
    COFD[2270] = -4.10945109E-01;
    COFD[2271] = 1.76651398E-02;
    COFD[2272] = -2.19317743E+01;
    COFD[2273] = 5.45216133E+00;
    COFD[2274] = -4.52916925E-01;
    COFD[2275] = 1.80456400E-02;
    COFD[2276] = -1.81735763E+01;
    COFD[2277] = 4.38391495E+00;
    COFD[2278] = -3.54941287E-01;
    COFD[2279] = 1.54195107E-02;
    COFD[2280] = -2.13425698E+01;
    COFD[2281] = 5.40460130E+00;
    COFD[2282] = -4.72718910E-01;
    COFD[2283] = 1.99362717E-02;
    COFD[2284] = -2.09191285E+01;
    COFD[2285] = 5.30153901E+00;
    COFD[2286] = -4.63335119E-01;
    COFD[2287] = 1.96897053E-02;
    COFD[2288] = -2.14326461E+01;
    COFD[2289] = 5.41729961E+00;
    COFD[2290] = -4.73400281E-01;
    COFD[2291] = 1.99269567E-02;
    COFD[2292] = -2.22116706E+01;
    COFD[2293] = 5.54251230E+00;
    COFD[2294] = -4.70946314E-01;
    COFD[2295] = 1.90785869E-02;
    COFD[2296] = -2.18590741E+01;
    COFD[2297] = 5.47368915E+00;
    COFD[2298] = -4.79424291E-01;
    COFD[2299] = 2.01372920E-02;
    COFD[2300] = -2.26242685E+01;
    COFD[2301] = 5.50195433E+00;
    COFD[2302] = -4.62690815E-01;
    COFD[2303] = 1.86010137E-02;
    COFD[2304] = -2.26924003E+01;
    COFD[2305] = 5.45651759E+00;
    COFD[2306] = -4.53709578E-01;
    COFD[2307] = 1.80890894E-02;
    COFD[2308] = -2.26924062E+01;
    COFD[2309] = 5.47738722E+00;
    COFD[2310] = -4.57512679E-01;
    COFD[2311] = 1.82977612E-02;
    COFD[2312] = -2.22623522E+01;
    COFD[2313] = 5.18930931E+00;
    COFD[2314] = -4.07143570E-01;
    COFD[2315] = 1.55986909E-02;
    COFD[2316] = -1.81432461E+01;
    COFD[2317] = 4.37565431E+00;
    COFD[2318] = -3.53906025E-01;
    COFD[2319] = 1.53760786E-02;
    COFD[2320] = -2.23333626E+01;
    COFD[2321] = 5.13339632E+00;
    COFD[2322] = -3.97833815E-01;
    COFD[2323] = 1.51137805E-02;
    COFD[2324] = -2.12121370E+01;
    COFD[2325] = 5.39823225E+00;
    COFD[2326] = -4.72294645E-01;
    COFD[2327] = 1.99340225E-02;
    COFD[2328] = -1.87897298E+01;
    COFD[2329] = 4.66162351E+00;
    COFD[2330] = -3.88920477E-01;
    COFD[2331] = 1.68089648E-02;
    COFD[2332] = -1.34709807E+01;
    COFD[2333] = 3.09379603E+00;
    COFD[2334] = -1.91268635E-01;
    COFD[2335] = 8.47480224E-03;
    COFD[2336] = -2.14369874E+01;
    COFD[2337] = 5.37331605E+00;
    COFD[2338] = -4.70491203E-01;
    COFD[2339] = 1.99134666E-02;
    COFD[2340] = -2.22817344E+01;
    COFD[2341] = 5.59185582E+00;
    COFD[2342] = -4.91155812E-01;
    COFD[2343] = 2.05043018E-02;
    COFD[2344] = -2.22817344E+01;
    COFD[2345] = 5.59185582E+00;
    COFD[2346] = -4.91155812E-01;
    COFD[2347] = 2.05043018E-02;
    COFD[2348] = -2.25690703E+01;
    COFD[2349] = 5.58516153E+00;
    COFD[2350] = -4.80440252E-01;
    COFD[2351] = 1.96479586E-02;
    COFD[2352] = -2.26781083E+01;
    COFD[2353] = 5.56399361E+00;
    COFD[2354] = -4.75020089E-01;
    COFD[2355] = 1.93065281E-02;
    COFD[2356] = -1.72232223E+01;
    COFD[2357] = 4.69060745E+00;
    COFD[2358] = -3.92369888E-01;
    COFD[2359] = 1.69459661E-02;
    COFD[2360] = -1.65488358E+01;
    COFD[2361] = 3.95035840E+00;
    COFD[2362] = -3.00959418E-01;
    COFD[2363] = 1.31692593E-02;
    COFD[2364] = -1.65675362E+01;
    COFD[2365] = 3.95035840E+00;
    COFD[2366] = -3.00959418E-01;
    COFD[2367] = 1.31692593E-02;
    COFD[2368] = -1.78903913E+01;
    COFD[2369] = 4.29613154E+00;
    COFD[2370] = -3.44012526E-01;
    COFD[2371] = 1.49643715E-02;
    COFD[2372] = -2.15095980E+01;
    COFD[2373] = 5.46737673E+00;
    COFD[2374] = -4.55696085E-01;
    COFD[2375] = 1.81982625E-02;
    COFD[2376] = -1.78969684E+01;
    COFD[2377] = 4.29613154E+00;
    COFD[2378] = -3.44012526E-01;
    COFD[2379] = 1.49643715E-02;
    COFD[2380] = -1.78834935E+01;
    COFD[2381] = 4.29613154E+00;
    COFD[2382] = -3.44012526E-01;
    COFD[2383] = 1.49643715E-02;
    COFD[2384] = -1.88463816E+01;
    COFD[2385] = 4.68393046E+00;
    COFD[2386] = -3.91610863E-01;
    COFD[2387] = 1.69174645E-02;
    COFD[2388] = -2.20036369E+01;
    COFD[2389] = 5.55935694E+00;
    COFD[2390] = -4.74154740E-01;
    COFD[2391] = 1.92584304E-02;
    COFD[2392] = -1.76285640E+01;
    COFD[2393] = 4.19935698E+00;
    COFD[2394] = -3.32310212E-01;
    COFD[2395] = 1.44920670E-02;
    COFD[2396] = -2.09191285E+01;
    COFD[2397] = 5.30153901E+00;
    COFD[2398] = -4.63335119E-01;
    COFD[2399] = 1.96897053E-02;
    COFD[2400] = -2.03766950E+01;
    COFD[2401] = 5.13263469E+00;
    COFD[2402] = -4.44457285E-01;
    COFD[2403] = 1.89932102E-02;
    COFD[2404] = -2.10944088E+01;
    COFD[2405] = 5.34286099E+00;
    COFD[2406] = -4.68100992E-01;
    COFD[2407] = 1.98731399E-02;
    COFD[2408] = -2.21065306E+01;
    COFD[2409] = 5.58360799E+00;
    COFD[2410] = -4.82701436E-01;
    COFD[2411] = 1.98437922E-02;
    COFD[2412] = -2.14323189E+01;
    COFD[2413] = 5.37331605E+00;
    COFD[2414] = -4.70491203E-01;
    COFD[2415] = 1.99134666E-02;
    COFD[2416] = -2.25674389E+01;
    COFD[2417] = 5.58516153E+00;
    COFD[2418] = -4.80440252E-01;
    COFD[2419] = 1.96479586E-02;
    COFD[2420] = -2.26768913E+01;
    COFD[2421] = 5.56399361E+00;
    COFD[2422] = -4.75020089E-01;
    COFD[2423] = 1.93065281E-02;
    COFD[2424] = -2.26583214E+01;
    COFD[2425] = 5.57930518E+00;
    COFD[2426] = -4.77998103E-01;
    COFD[2427] = 1.94753934E-02;
    COFD[2428] = -2.25093649E+01;
    COFD[2429] = 5.43187491E+00;
    COFD[2430] = -4.49258240E-01;
    COFD[2431] = 1.78460456E-02;
    COFD[2432] = -1.76002031E+01;
    COFD[2433] = 4.19171952E+00;
    COFD[2434] = -3.31354810E-01;
    COFD[2435] = 1.44520623E-02;
    COFD[2436] = -2.16832171E+01;
    COFD[2437] = 4.73756385E+00;
    COFD[2438] = -3.34700069E-01;
    COFD[2439] = 1.19121851E-02;
    COFD[2440] = -2.18273547E+01;
    COFD[2441] = 5.55753905E+00;
    COFD[2442] = -4.88136714E-01;
    COFD[2443] = 2.04294957E-02;
    COFD[2444] = -1.94823660E+01;
    COFD[2445] = 4.87333294E+00;
    COFD[2446] = -4.13769241E-01;
    COFD[2447] = 1.77802244E-02;
    COFD[2448] = -1.39924781E+01;
    COFD[2449] = 3.26384506E+00;
    COFD[2450] = -2.12947087E-01;
    COFD[2451] = 9.39743888E-03;
    COFD[2452] = -2.20033797E+01;
    COFD[2453] = 5.51276597E+00;
    COFD[2454] = -4.83701824E-01;
    COFD[2455] = 2.02915297E-02;
    COFD[2456] = -2.25119517E+01;
    COFD[2457] = 5.58206320E+00;
    COFD[2458] = -4.82956809E-01;
    COFD[2459] = 1.98731634E-02;
    COFD[2460] = -2.25119517E+01;
    COFD[2461] = 5.58206320E+00;
    COFD[2462] = -4.82956809E-01;
    COFD[2463] = 1.98731634E-02;
    COFD[2464] = -2.26092626E+01;
    COFD[2465] = 5.48733210E+00;
    COFD[2466] = -4.59335784E-01;
    COFD[2467] = 1.83981530E-02;
    COFD[2468] = -2.26225827E+01;
    COFD[2469] = 5.42359285E+00;
    COFD[2470] = -4.47762010E-01;
    COFD[2471] = 1.77644050E-02;
    COFD[2472] = -1.79344949E+01;
    COFD[2473] = 4.91373893E+00;
    COFD[2474] = -4.18747629E-01;
    COFD[2475] = 1.79856610E-02;
    COFD[2476] = -1.72556499E+01;
    COFD[2477] = 4.17889917E+00;
    COFD[2478] = -3.29752510E-01;
    COFD[2479] = 1.43850275E-02;
    COFD[2480] = -1.72753760E+01;
    COFD[2481] = 4.17889917E+00;
    COFD[2482] = -3.29752510E-01;
    COFD[2483] = 1.43850275E-02;
    COFD[2484] = -1.86499071E+01;
    COFD[2485] = 4.53572533E+00;
    COFD[2486] = -3.73386925E-01;
    COFD[2487] = 1.61678881E-02;
    COFD[2488] = -2.12661865E+01;
    COFD[2489] = 5.24930667E+00;
    COFD[2490] = -4.17435088E-01;
    COFD[2491] = 1.61434424E-02;
    COFD[2492] = -1.86570209E+01;
    COFD[2493] = 4.53572533E+00;
    COFD[2494] = -3.73386925E-01;
    COFD[2495] = 1.61678881E-02;
    COFD[2496] = -1.86424545E+01;
    COFD[2497] = 4.53572533E+00;
    COFD[2498] = -3.73386925E-01;
    COFD[2499] = 1.61678881E-02;
    COFD[2500] = -1.95552142E+01;
    COFD[2501] = 4.90255048E+00;
    COFD[2502] = -4.17368501E-01;
    COFD[2503] = 1.79287358E-02;
    COFD[2504] = -2.19399793E+01;
    COFD[2505] = 5.41841631E+00;
    COFD[2506] = -4.46818971E-01;
    COFD[2507] = 1.77127652E-02;
    COFD[2508] = -1.83538377E+01;
    COFD[2509] = 4.42828044E+00;
    COFD[2510] = -3.60417833E-01;
    COFD[2511] = 1.56455103E-02;
    COFD[2512] = -2.14326461E+01;
    COFD[2513] = 5.41729961E+00;
    COFD[2514] = -4.73400281E-01;
    COFD[2515] = 1.99269567E-02;
    COFD[2516] = -2.10944088E+01;
    COFD[2517] = 5.34286099E+00;
    COFD[2518] = -4.68100992E-01;
    COFD[2519] = 1.98731399E-02;
    COFD[2520] = -2.15746136E+01;
    COFD[2521] = 5.44803850E+00;
    COFD[2522] = -4.76610560E-01;
    COFD[2523] = 2.00355294E-02;
    COFD[2524] = -2.22159630E+01;
    COFD[2525] = 5.51722375E+00;
    COFD[2526] = -4.66081431E-01;
    COFD[2527] = 1.88044011E-02;
    COFD[2528] = -2.19982918E+01;
    COFD[2529] = 5.51276597E+00;
    COFD[2530] = -4.83701824E-01;
    COFD[2531] = 2.02915297E-02;
    COFD[2532] = -2.26074491E+01;
    COFD[2533] = 5.48733210E+00;
    COFD[2534] = -4.59335784E-01;
    COFD[2535] = 1.83981530E-02;
    COFD[2536] = -2.26212238E+01;
    COFD[2537] = 5.42359285E+00;
    COFD[2538] = -4.47762010E-01;
    COFD[2539] = 1.77644050E-02;
    COFD[2540] = -2.26194619E+01;
    COFD[2541] = 5.44502013E+00;
    COFD[2542] = -4.51625610E-01;
    COFD[2543] = 1.79750820E-02;
    COFD[2544] = -2.22607146E+01;
    COFD[2545] = 5.20705035E+00;
    COFD[2546] = -4.10106772E-01;
    COFD[2547] = 1.57532724E-02;
    COFD[2548] = -1.83249299E+01;
    COFD[2549] = 4.42045763E+00;
    COFD[2550] = -3.59451578E-01;
    COFD[2551] = 1.56056164E-02;
    COFD[2552] = -1.83142602E+01;
    COFD[2553] = 3.06214549E+00;
    COFD[2554] = -8.03954556E-02;
    COFD[2555] = -5.30273859E-04;
    COFD[2556] = -2.20453723E+01;
    COFD[2557] = 5.44448440E+00;
    COFD[2558] = -4.51529024E-01;
    COFD[2559] = 1.79698119E-02;
    COFD[2560] = -2.11309207E+01;
    COFD[2561] = 5.41773516E+00;
    COFD[2562] = -4.73414338E-01;
    COFD[2563] = 1.99258685E-02;
    COFD[2564] = -1.57034851E+01;
    COFD[2565] = 3.93614244E+00;
    COFD[2566] = -2.99111497E-01;
    COFD[2567] = 1.30888229E-02;
    COFD[2568] = -2.24615468E+01;
    COFD[2569] = 5.49330641E+00;
    COFD[2570] = -4.60498247E-01;
    COFD[2571] = 1.84639199E-02;
    COFD[2572] = -2.20355148E+01;
    COFD[2573] = 5.14570932E+00;
    COFD[2574] = -3.99877142E-01;
    COFD[2575] = 1.52199557E-02;
    COFD[2576] = -2.20355148E+01;
    COFD[2577] = 5.14570932E+00;
    COFD[2578] = -3.99877142E-01;
    COFD[2579] = 1.52199557E-02;
    COFD[2580] = -2.10753161E+01;
    COFD[2581] = 4.59015462E+00;
    COFD[2582] = -3.11301450E-01;
    COFD[2583] = 1.07307395E-02;
    COFD[2584] = -2.08711113E+01;
    COFD[2585] = 4.43077544E+00;
    COFD[2586] = -2.86496187E-01;
    COFD[2587] = 9.49507497E-03;
    COFD[2588] = -1.94688688E+01;
    COFD[2589] = 5.43830787E+00;
    COFD[2590] = -4.75472880E-01;
    COFD[2591] = 1.99909996E-02;
    COFD[2592] = -1.90883268E+01;
    COFD[2593] = 4.84384483E+00;
    COFD[2594] = -4.10265575E-01;
    COFD[2595] = 1.76414287E-02;
    COFD[2596] = -1.91102652E+01;
    COFD[2597] = 4.84384483E+00;
    COFD[2598] = -4.10265575E-01;
    COFD[2599] = 1.76414287E-02;
    COFD[2600] = -2.05272328E+01;
    COFD[2601] = 5.18417470E+00;
    COFD[2602] = -4.49491573E-01;
    COFD[2603] = 1.91438508E-02;
    COFD[2604] = -1.87383952E+01;
    COFD[2605] = 3.96926341E+00;
    COFD[2606] = -2.16412264E-01;
    COFD[2607] = 6.06012078E-03;
    COFD[2608] = -2.05356023E+01;
    COFD[2609] = 5.18417470E+00;
    COFD[2610] = -4.49491573E-01;
    COFD[2611] = 1.91438508E-02;
    COFD[2612] = -2.05184870E+01;
    COFD[2613] = 5.18417470E+00;
    COFD[2614] = -4.49491573E-01;
    COFD[2615] = 1.91438508E-02;
    COFD[2616] = -2.11606963E+01;
    COFD[2617] = 5.42846112E+00;
    COFD[2618] = -4.74321870E-01;
    COFD[2619] = 1.99459749E-02;
    COFD[2620] = -2.01015340E+01;
    COFD[2621] = 4.41511629E+00;
    COFD[2622] = -2.84086963E-01;
    COFD[2623] = 9.37586971E-03;
    COFD[2624] = -2.02922701E+01;
    COFD[2625] = 5.11106992E+00;
    COFD[2626] = -4.42047129E-01;
    COFD[2627] = 1.89042990E-02;
    COFD[2628] = -2.22116706E+01;
    COFD[2629] = 5.54251230E+00;
    COFD[2630] = -4.70946314E-01;
    COFD[2631] = 1.90785869E-02;
    COFD[2632] = -2.21065306E+01;
    COFD[2633] = 5.58360799E+00;
    COFD[2634] = -4.82701436E-01;
    COFD[2635] = 1.98437922E-02;
    COFD[2636] = -2.22159630E+01;
    COFD[2637] = 5.51722375E+00;
    COFD[2638] = -4.66081431E-01;
    COFD[2639] = 1.88044011E-02;
    COFD[2640] = -2.09002742E+01;
    COFD[2641] = 4.72895031E+00;
    COFD[2642] = -3.33332771E-01;
    COFD[2643] = 1.18431478E-02;
    COFD[2644] = -2.24554521E+01;
    COFD[2645] = 5.49330641E+00;
    COFD[2646] = -4.60498247E-01;
    COFD[2647] = 1.84639199E-02;
    COFD[2648] = -2.10730345E+01;
    COFD[2649] = 4.59015462E+00;
    COFD[2650] = -3.11301450E-01;
    COFD[2651] = 1.07307395E-02;
    COFD[2652] = -2.08693817E+01;
    COFD[2653] = 4.43077544E+00;
    COFD[2654] = -2.86496187E-01;
    COFD[2655] = 9.49507497E-03;
    COFD[2656] = -2.09501825E+01;
    COFD[2657] = 4.48561185E+00;
    COFD[2658] = -2.94994032E-01;
    COFD[2659] = 9.91723779E-03;
    COFD[2660] = -1.96907838E+01;
    COFD[2661] = 3.86787656E+00;
    COFD[2662] = -2.00935789E-01;
    COFD[2663] = 5.30029143E-03;
    COFD[2664] = -2.02646611E+01;
    COFD[2665] = 5.10426133E+00;
    COFD[2666] = -4.41256919E-01;
    COFD[2667] = 1.88737290E-02;
    COFD[2668] = -2.16922067E+01;
    COFD[2669] = 4.61366606E+00;
    COFD[2670] = -3.15066862E-01;
    COFD[2671] = 1.09216645E-02;
    COFD[2672] = -2.21885140E+01;
    COFD[2673] = 5.59472344E+00;
    COFD[2674] = -4.91421518E-01;
    COFD[2675] = 2.05117088E-02;
    COFD[2676] = -1.99235839E+01;
    COFD[2677] = 4.95514826E+00;
    COFD[2678] = -4.23691395E-01;
    COFD[2679] = 1.81828318E-02;
    COFD[2680] = -1.43129712E+01;
    COFD[2681] = 3.31177824E+00;
    COFD[2682] = -2.18945280E-01;
    COFD[2683] = 9.64764419E-03;
    COFD[2684] = -2.23903059E+01;
    COFD[2685] = 5.56066804E+00;
    COFD[2686] = -4.88405706E-01;
    COFD[2687] = 2.04357330E-02;
    COFD[2688] = -2.28056965E+01;
    COFD[2689] = 5.58523510E+00;
    COFD[2690] = -4.81201481E-01;
    COFD[2691] = 1.97107111E-02;
    COFD[2692] = -2.28056965E+01;
    COFD[2693] = 5.58523510E+00;
    COFD[2694] = -4.81201481E-01;
    COFD[2695] = 1.97107111E-02;
    COFD[2696] = -2.27764511E+01;
    COFD[2697] = 5.43128350E+00;
    COFD[2698] = -4.49151750E-01;
    COFD[2699] = 1.78402439E-02;
    COFD[2700] = -2.27897781E+01;
    COFD[2701] = 5.36784270E+00;
    COFD[2702] = -4.37680260E-01;
    COFD[2703] = 1.72142944E-02;
    COFD[2704] = -1.83539686E+01;
    COFD[2705] = 4.98756925E+00;
    COFD[2706] = -4.27526072E-01;
    COFD[2707] = 1.83341755E-02;
    COFD[2708] = -1.76775033E+01;
    COFD[2709] = 4.24719726E+00;
    COFD[2710] = -3.38206061E-01;
    COFD[2711] = 1.47350654E-02;
    COFD[2712] = -1.76992976E+01;
    COFD[2713] = 4.24719726E+00;
    COFD[2714] = -3.38206061E-01;
    COFD[2715] = 1.47350654E-02;
    COFD[2716] = -1.91208314E+01;
    COFD[2717] = 4.61801405E+00;
    COFD[2718] = -3.83535652E-01;
    COFD[2719] = 1.65862513E-02;
    COFD[2720] = -2.13847439E+01;
    COFD[2721] = 5.17440955E+00;
    COFD[2722] = -4.04678430E-01;
    COFD[2723] = 1.54706350E-02;
    COFD[2724] = -1.91291147E+01;
    COFD[2725] = 4.61801405E+00;
    COFD[2726] = -3.83535652E-01;
    COFD[2727] = 1.65862513E-02;
    COFD[2728] = -1.91121742E+01;
    COFD[2729] = 4.61801405E+00;
    COFD[2730] = -3.83535652E-01;
    COFD[2731] = 1.65862513E-02;
    COFD[2732] = -1.99803490E+01;
    COFD[2733] = 4.97875278E+00;
    COFD[2734] = -4.26485475E-01;
    COFD[2735] = 1.82931933E-02;
    COFD[2736] = -2.20947902E+01;
    COFD[2737] = 5.36053938E+00;
    COFD[2738] = -4.36434519E-01;
    COFD[2739] = 1.71484255E-02;
    COFD[2740] = -1.87685041E+01;
    COFD[2741] = 4.49191492E+00;
    COFD[2742] = -3.68041771E-01;
    COFD[2743] = 1.59498676E-02;
    COFD[2744] = -2.18590741E+01;
    COFD[2745] = 5.47368915E+00;
    COFD[2746] = -4.79424291E-01;
    COFD[2747] = 2.01372920E-02;
    COFD[2748] = -2.14323189E+01;
    COFD[2749] = 5.37331605E+00;
    COFD[2750] = -4.70491203E-01;
    COFD[2751] = 1.99134666E-02;
    COFD[2752] = -2.19982918E+01;
    COFD[2753] = 5.51276597E+00;
    COFD[2754] = -4.83701824E-01;
    COFD[2755] = 2.02915297E-02;
    COFD[2756] = -2.24554521E+01;
    COFD[2757] = 5.49330641E+00;
    COFD[2758] = -4.60498247E-01;
    COFD[2759] = 1.84639199E-02;
    COFD[2760] = -2.23842815E+01;
    COFD[2761] = 5.56066804E+00;
    COFD[2762] = -4.88405706E-01;
    COFD[2763] = 2.04357330E-02;
    COFD[2764] = -2.27742038E+01;
    COFD[2765] = 5.43128350E+00;
    COFD[2766] = -4.49151750E-01;
    COFD[2767] = 1.78402439E-02;
    COFD[2768] = -2.27880759E+01;
    COFD[2769] = 5.36784270E+00;
    COFD[2770] = -4.37680260E-01;
    COFD[2771] = 1.72142944E-02;
    COFD[2772] = -2.28041813E+01;
    COFD[2773] = 5.39331915E+00;
    COFD[2774] = -4.42116537E-01;
    COFD[2775] = 1.74516613E-02;
    COFD[2776] = -2.23601792E+01;
    COFD[2777] = 5.11898303E+00;
    COFD[2778] = -3.95452058E-01;
    COFD[2779] = 1.49902966E-02;
    COFD[2780] = -1.87434358E+01;
    COFD[2781] = 4.48550694E+00;
    COFD[2782] = -3.67277454E-01;
    COFD[2783] = 1.59194755E-02;
    COFD[2784] = -1.80988918E+01;
    COFD[2785] = 2.75595461E+00;
    COFD[2786] = -3.54899742E-02;
    COFD[2787] = -2.67520318E-03;
    COFD[2788] = -2.23716664E+01;
    COFD[2789] = 5.38363480E+00;
    COFD[2790] = -4.40406530E-01;
    COFD[2791] = 1.73594979E-02;
    COFD[2792] = -2.17251451E+01;
    COFD[2793] = 5.48539572E+00;
    COFD[2794] = -4.80731929E-01;
    COFD[2795] = 2.01857298E-02;
    COFD[2796] = -1.63063503E+01;
    COFD[2797] = 4.02718684E+00;
    COFD[2798] = -3.10842073E-01;
    COFD[2799] = 1.35952684E-02;
    COFD[2800] = -2.27820795E+01;
    COFD[2801] = 5.43128350E+00;
    COFD[2802] = -4.49151750E-01;
    COFD[2803] = 1.78402439E-02;
    COFD[2804] = -2.22439282E+01;
    COFD[2805] = 5.02889918E+00;
    COFD[2806] = -3.80860199E-01;
    COFD[2807] = 1.42428137E-02;
    COFD[2808] = -2.22439282E+01;
    COFD[2809] = 5.02889918E+00;
    COFD[2810] = -3.80860199E-01;
    COFD[2811] = 1.42428137E-02;
    COFD[2812] = -2.12007188E+01;
    COFD[2813] = 4.42962851E+00;
    COFD[2814] = -2.86319550E-01;
    COFD[2815] = 9.48633214E-03;
    COFD[2816] = -2.09377754E+01;
    COFD[2817] = 4.24466287E+00;
    COFD[2818] = -2.57969164E-01;
    COFD[2819] = 8.08847972E-03;
    COFD[2820] = -2.01029332E+01;
    COFD[2821] = 5.51409582E+00;
    COFD[2822] = -4.83844807E-01;
    COFD[2823] = 2.02965712E-02;
    COFD[2824] = -1.97750001E+01;
    COFD[2825] = 4.93652913E+00;
    COFD[2826] = -4.21485300E-01;
    COFD[2827] = 1.80955487E-02;
    COFD[2828] = -1.98001639E+01;
    COFD[2829] = 4.93652913E+00;
    COFD[2830] = -4.21485300E-01;
    COFD[2831] = 1.80955487E-02;
    COFD[2832] = -2.12111746E+01;
    COFD[2833] = 5.26829491E+00;
    COFD[2834] = -4.59325230E-01;
    COFD[2835] = 1.95275380E-02;
    COFD[2836] = -1.74914373E+01;
    COFD[2837] = 3.15361887E+00;
    COFD[2838] = -9.38468185E-02;
    COFD[2839] = 1.12007553E-04;
    COFD[2840] = -2.12216591E+01;
    COFD[2841] = 5.26829491E+00;
    COFD[2842] = -4.59325230E-01;
    COFD[2843] = 1.95275380E-02;
    COFD[2844] = -2.12002655E+01;
    COFD[2845] = 5.26829491E+00;
    COFD[2846] = -4.59325230E-01;
    COFD[2847] = 1.95275380E-02;
    COFD[2848] = -2.17753205E+01;
    COFD[2849] = 5.50642944E+00;
    COFD[2850] = -4.83016239E-01;
    COFD[2851] = 2.02671529E-02;
    COFD[2852] = -2.01531862E+01;
    COFD[2853] = 4.23025405E+00;
    COFD[2854] = -2.55784321E-01;
    COFD[2855] = 7.98147639E-03;
    COFD[2856] = -2.08696226E+01;
    COFD[2857] = 5.15826153E+00;
    COFD[2858] = -4.46699782E-01;
    COFD[2859] = 1.90459433E-02;
    COFD[2860] = -2.26242685E+01;
    COFD[2861] = 5.50195433E+00;
    COFD[2862] = -4.62690815E-01;
    COFD[2863] = 1.86010137E-02;
    COFD[2864] = -2.25674389E+01;
    COFD[2865] = 5.58516153E+00;
    COFD[2866] = -4.80440252E-01;
    COFD[2867] = 1.96479586E-02;
    COFD[2868] = -2.26074491E+01;
    COFD[2869] = 5.48733210E+00;
    COFD[2870] = -4.59335784E-01;
    COFD[2871] = 1.83981530E-02;
    COFD[2872] = -2.10730345E+01;
    COFD[2873] = 4.59015462E+00;
    COFD[2874] = -3.11301450E-01;
    COFD[2875] = 1.07307395E-02;
    COFD[2876] = -2.27742038E+01;
    COFD[2877] = 5.43128350E+00;
    COFD[2878] = -4.49151750E-01;
    COFD[2879] = 1.78402439E-02;
    COFD[2880] = -2.11974824E+01;
    COFD[2881] = 4.42962851E+00;
    COFD[2882] = -2.86319550E-01;
    COFD[2883] = 9.48633214E-03;
    COFD[2884] = -2.09352623E+01;
    COFD[2885] = 4.24466287E+00;
    COFD[2886] = -2.57969164E-01;
    COFD[2887] = 8.08847972E-03;
    COFD[2888] = -2.10476554E+01;
    COFD[2889] = 4.30471901E+00;
    COFD[2890] = -2.67147055E-01;
    COFD[2891] = 8.53981315E-03;
    COFD[2892] = -1.94513324E+01;
    COFD[2893] = 3.52672843E+00;
    COFD[2894] = -1.49232080E-01;
    COFD[2895] = 2.77504261E-03;
    COFD[2896] = -2.08530414E+01;
    COFD[2897] = 5.15581320E+00;
    COFD[2898] = -4.46543388E-01;
    COFD[2899] = 1.90458118E-02;
    COFD[2900] = -1.76693225E+01;
    COFD[2901] = 2.50076717E+00;
    COFD[2902] = 2.22033801E-03;
    COFD[2903] = -4.48335559E-03;
    COFD[2904] = -2.23376752E+01;
    COFD[2905] = 5.29923580E+00;
    COFD[2906] = -4.25997743E-01;
    COFD[2907] = 1.65974476E-02;
    COFD[2908] = -2.20306514E+01;
    COFD[2909] = 5.54982225E+00;
    COFD[2910] = -4.87416516E-01;
    COFD[2911] = 2.04093655E-02;
    COFD[2912] = -1.66209380E+01;
    COFD[2913] = 4.10820924E+00;
    COFD[2914] = -3.21095223E-01;
    COFD[2915] = 1.40301968E-02;
    COFD[2916] = -2.27964005E+01;
    COFD[2917] = 5.36784270E+00;
    COFD[2918] = -4.37680260E-01;
    COFD[2919] = 1.72142944E-02;
    COFD[2920] = -2.20739808E+01;
    COFD[2921] = 4.88495888E+00;
    COFD[2922] = -3.58090775E-01;
    COFD[2923] = 1.30937796E-02;
    COFD[2924] = -2.20739808E+01;
    COFD[2925] = 4.88495888E+00;
    COFD[2926] = -3.58090775E-01;
    COFD[2927] = 1.30937796E-02;
    COFD[2928] = -2.09387694E+01;
    COFD[2929] = 4.24466287E+00;
    COFD[2930] = -2.57969164E-01;
    COFD[2931] = 8.08847972E-03;
    COFD[2932] = -2.05939847E+01;
    COFD[2933] = 4.02612453E+00;
    COFD[2934] = -2.25065111E-01;
    COFD[2935] = 6.48488943E-03;
    COFD[2936] = -2.03845035E+01;
    COFD[2937] = 5.56619223E+00;
    COFD[2938] = -4.88860312E-01;
    COFD[2939] = 2.04450332E-02;
    COFD[2940] = -2.01211076E+01;
    COFD[2941] = 5.01930029E+00;
    COFD[2942] = -4.31305027E-01;
    COFD[2943] = 1.84846322E-02;
    COFD[2944] = -2.01469731E+01;
    COFD[2945] = 5.01930029E+00;
    COFD[2946] = -4.31305027E-01;
    COFD[2947] = 1.84846322E-02;
    COFD[2948] = -2.15609288E+01;
    COFD[2949] = 5.34961360E+00;
    COFD[2950] = -4.68762622E-01;
    COFD[2951] = 1.98932200E-02;
    COFD[2952] = -1.83862949E+01;
    COFD[2953] = 3.54020956E+00;
    COFD[2954] = -1.52304845E-01;
    COFD[2955] = 2.96268912E-03;
    COFD[2956] = -2.15719260E+01;
    COFD[2957] = 5.34961360E+00;
    COFD[2958] = -4.68762622E-01;
    COFD[2959] = 1.98932200E-02;
    COFD[2960] = -2.15521921E+01;
    COFD[2961] = 5.35030481E+00;
    COFD[2962] = -4.68818560E-01;
    COFD[2963] = 1.98942796E-02;
    COFD[2964] = -2.20571038E+01;
    COFD[2965] = 5.56193702E+00;
    COFD[2966] = -4.88512557E-01;
    COFD[2967] = 2.04380728E-02;
    COFD[2968] = -1.98052638E+01;
    COFD[2969] = 4.01137806E+00;
    COFD[2970] = -2.22734310E-01;
    COFD[2971] = 6.36741567E-03;
    COFD[2972] = -2.12005646E+01;
    COFD[2973] = 5.23328793E+00;
    COFD[2974] = -4.55218686E-01;
    COFD[2975] = 1.93667472E-02;
    COFD[2976] = -2.26924003E+01;
    COFD[2977] = 5.45651759E+00;
    COFD[2978] = -4.53709578E-01;
    COFD[2979] = 1.80890894E-02;
    COFD[2980] = -2.26768913E+01;
    COFD[2981] = 5.56399361E+00;
    COFD[2982] = -4.75020089E-01;
    COFD[2983] = 1.93065281E-02;
    COFD[2984] = -2.26212238E+01;
    COFD[2985] = 5.42359285E+00;
    COFD[2986] = -4.47762010E-01;
    COFD[2987] = 1.77644050E-02;
    COFD[2988] = -2.08693817E+01;
    COFD[2989] = 4.43077544E+00;
    COFD[2990] = -2.86496187E-01;
    COFD[2991] = 9.49507497E-03;
    COFD[2992] = -2.27880759E+01;
    COFD[2993] = 5.36784270E+00;
    COFD[2994] = -4.37680260E-01;
    COFD[2995] = 1.72142944E-02;
    COFD[2996] = -2.09352623E+01;
    COFD[2997] = 4.24466287E+00;
    COFD[2998] = -2.57969164E-01;
    COFD[2999] = 8.08847972E-03;
    COFD[3000] = -2.05912423E+01;
    COFD[3001] = 4.02612453E+00;
    COFD[3002] = -2.25065111E-01;
    COFD[3003] = 6.48488943E-03;
    COFD[3004] = -2.07445786E+01;
    COFD[3005] = 4.10156709E+00;
    COFD[3006] = -2.36431182E-01;
    COFD[3007] = 7.03931247E-03;
    COFD[3008] = -1.93408889E+01;
    COFD[3009] = 3.42798339E+00;
    COFD[3010] = -1.35210942E-01;
    COFD[3011] = 2.12316045E-03;
    COFD[3012] = -2.11724381E+01;
    COFD[3013] = 5.22581725E+00;
    COFD[3014] = -4.54351699E-01;
    COFD[3015] = 1.93332012E-02;
    COFD[3016] = -1.78388304E+01;
    COFD[3017] = 2.57553304E+00;
    COFD[3018] = -8.86366178E-03;
    COFD[3019] = -3.95164836E-03;
    COFD[3020] = -2.23523338E+01;
    COFD[3021] = 5.32879676E+00;
    COFD[3022] = -4.30989895E-01;
    COFD[3023] = 1.68597740E-02;
    COFD[3024] = -2.19259008E+01;
    COFD[3025] = 5.53120252E+00;
    COFD[3026] = -4.85560309E-01;
    COFD[3027] = 2.03509814E-02;
    COFD[3028] = -1.64898075E+01;
    COFD[3029] = 4.08133132E+00;
    COFD[3030] = -3.17684493E-01;
    COFD[3031] = 1.38851127E-02;
    COFD[3032] = -2.28129363E+01;
    COFD[3033] = 5.39331915E+00;
    COFD[3034] = -4.42116537E-01;
    COFD[3035] = 1.74516613E-02;
    COFD[3036] = -2.21500241E+01;
    COFD[3037] = 4.93303891E+00;
    COFD[3038] = -3.65678823E-01;
    COFD[3039] = 1.34760462E-02;
    COFD[3040] = -2.21500241E+01;
    COFD[3041] = 4.93303891E+00;
    COFD[3042] = -3.65678823E-01;
    COFD[3043] = 1.34760462E-02;
    COFD[3044] = -2.10514354E+01;
    COFD[3045] = 4.30471901E+00;
    COFD[3046] = -2.67147055E-01;
    COFD[3047] = 8.53981315E-03;
    COFD[3048] = -2.07475551E+01;
    COFD[3049] = 4.10156709E+00;
    COFD[3050] = -2.36431182E-01;
    COFD[3051] = 7.03931247E-03;
    COFD[3052] = -2.02678158E+01;
    COFD[3053] = 5.55236751E+00;
    COFD[3054] = -4.87675676E-01;
    COFD[3055] = 2.04178098E-02;
    COFD[3056] = -1.99915387E+01;
    COFD[3057] = 4.99136452E+00;
    COFD[3058] = -4.27975210E-01;
    COFD[3059] = 1.83519204E-02;
    COFD[3060] = -2.00180417E+01;
    COFD[3061] = 4.99136452E+00;
    COFD[3062] = -4.27975210E-01;
    COFD[3063] = 1.83519204E-02;
    COFD[3064] = -2.14647019E+01;
    COFD[3065] = 5.33054865E+00;
    COFD[3066] = -4.66761560E-01;
    COFD[3067] = 1.98252926E-02;
    COFD[3068] = -1.84877055E+01;
    COFD[3069] = 3.61246265E+00;
    COFD[3070] = -1.63179700E-01;
    COFD[3071] = 3.49175937E-03;
    COFD[3072] = -2.14761835E+01;
    COFD[3073] = 5.33054865E+00;
    COFD[3074] = -4.66761560E-01;
    COFD[3075] = 1.98252926E-02;
    COFD[3076] = -2.14573159E+01;
    COFD[3077] = 5.33204143E+00;
    COFD[3078] = -4.66932795E-01;
    COFD[3079] = 1.98318389E-02;
    COFD[3080] = -2.19586553E+01;
    COFD[3081] = 5.54610272E+00;
    COFD[3082] = -4.87043162E-01;
    COFD[3083] = 2.03974796E-02;
    COFD[3084] = -1.99254109E+01;
    COFD[3085] = 4.08598993E+00;
    COFD[3086] = -2.33972521E-01;
    COFD[3087] = 6.91534316E-03;
    COFD[3088] = -2.10856216E+01;
    COFD[3089] = 5.20655472E+00;
    COFD[3090] = -4.52109091E-01;
    COFD[3091] = 1.92461018E-02;
    COFD[3092] = -2.26924062E+01;
    COFD[3093] = 5.47738722E+00;
    COFD[3094] = -4.57512679E-01;
    COFD[3095] = 1.82977612E-02;
    COFD[3096] = -2.26583214E+01;
    COFD[3097] = 5.57930518E+00;
    COFD[3098] = -4.77998103E-01;
    COFD[3099] = 1.94753934E-02;
    COFD[3100] = -2.26194619E+01;
    COFD[3101] = 5.44502013E+00;
    COFD[3102] = -4.51625610E-01;
    COFD[3103] = 1.79750820E-02;
    COFD[3104] = -2.09501825E+01;
    COFD[3105] = 4.48561185E+00;
    COFD[3106] = -2.94994032E-01;
    COFD[3107] = 9.91723779E-03;
    COFD[3108] = -2.28041813E+01;
    COFD[3109] = 5.39331915E+00;
    COFD[3110] = -4.42116537E-01;
    COFD[3111] = 1.74516613E-02;
    COFD[3112] = -2.10476554E+01;
    COFD[3113] = 4.30471901E+00;
    COFD[3114] = -2.67147055E-01;
    COFD[3115] = 8.53981315E-03;
    COFD[3116] = -2.07445786E+01;
    COFD[3117] = 4.10156709E+00;
    COFD[3118] = -2.36431182E-01;
    COFD[3119] = 7.03931247E-03;
    COFD[3120] = -2.08861598E+01;
    COFD[3121] = 4.17077470E+00;
    COFD[3122] = -2.46899919E-01;
    COFD[3123] = 7.55092148E-03;
    COFD[3124] = -1.95094733E+01;
    COFD[3125] = 3.50834090E+00;
    COFD[3126] = -1.47169195E-01;
    COFD[3127] = 2.69999047E-03;
    COFD[3128] = -2.10567656E+01;
    COFD[3129] = 5.19883613E+00;
    COFD[3130] = -4.51205592E-01;
    COFD[3131] = 1.92107800E-02;
    COFD[3132] = -1.60462384E+01;
    COFD[3133] = 1.75146402E+00;
    COFD[3134] = 1.11279769E-01;
    COFD[3135] = -9.64032668E-03;
    COFD[3136] = -2.18341320E+01;
    COFD[3137] = 5.02186857E+00;
    COFD[3138] = -3.79735818E-01;
    COFD[3139] = 1.41857239E-02;
    COFD[3140] = -2.22862449E+01;
    COFD[3141] = 5.58519781E+00;
    COFD[3142] = -4.83910410E-01;
    COFD[3143] = 1.99365500E-02;
    COFD[3144] = -1.73937806E+01;
    COFD[3145] = 4.40038200E+00;
    COFD[3146] = -3.56980932E-01;
    COFD[3147] = 1.55040067E-02;
    COFD[3148] = -2.23688893E+01;
    COFD[3149] = 5.11898303E+00;
    COFD[3150] = -3.95452058E-01;
    COFD[3151] = 1.49902966E-02;
    COFD[3152] = -2.12908209E+01;
    COFD[3153] = 4.48269001E+00;
    COFD[3154] = -2.94539894E-01;
    COFD[3155] = 9.89462789E-03;
    COFD[3156] = -2.12908209E+01;
    COFD[3157] = 4.48269001E+00;
    COFD[3158] = -2.94539894E-01;
    COFD[3159] = 9.89462789E-03;
    COFD[3160] = -1.94550832E+01;
    COFD[3161] = 3.52672843E+00;
    COFD[3162] = -1.49232080E-01;
    COFD[3163] = 2.77504261E-03;
    COFD[3164] = -1.93438403E+01;
    COFD[3165] = 3.42798339E+00;
    COFD[3166] = -1.35210942E-01;
    COFD[3167] = 2.12316045E-03;
    COFD[3168] = -2.05542181E+01;
    COFD[3169] = 5.59742840E+00;
    COFD[3170] = -4.86964048E-01;
    COFD[3171] = 2.01281667E-02;
    COFD[3172] = -2.06365262E+01;
    COFD[3173] = 5.20105041E+00;
    COFD[3174] = -4.51465908E-01;
    COFD[3175] = 1.92210125E-02;
    COFD[3176] = -2.06629641E+01;
    COFD[3177] = 5.20105041E+00;
    COFD[3178] = -4.51465908E-01;
    COFD[3179] = 1.92210125E-02;
    COFD[3180] = -2.19080065E+01;
    COFD[3181] = 5.44858060E+00;
    COFD[3182] = -4.76668834E-01;
    COFD[3183] = 2.00375810E-02;
    COFD[3184] = -1.54850731E+01;
    COFD[3185] = 2.22760630E+00;
    COFD[3186] = 3.88962143E-02;
    COFD[3187] = -6.06812513E-03;
    COFD[3188] = -2.19194378E+01;
    COFD[3189] = 5.44858060E+00;
    COFD[3190] = -4.76668834E-01;
    COFD[3191] = 2.00375810E-02;
    COFD[3192] = -2.20898696E+01;
    COFD[3193] = 5.51001635E+00;
    COFD[3194] = -4.83405862E-01;
    COFD[3195] = 2.02810787E-02;
    COFD[3196] = -2.22630459E+01;
    COFD[3197] = 5.60031657E+00;
    COFD[3198] = -4.87628580E-01;
    COFD[3199] = 2.01686065E-02;
    COFD[3200] = -1.85959819E+01;
    COFD[3201] = 3.43982423E+00;
    COFD[3202] = -1.36303972E-01;
    COFD[3203] = 2.15194160E-03;
    COFD[3204] = -2.17920541E+01;
    COFD[3205] = 5.41750012E+00;
    COFD[3206] = -4.73406713E-01;
    COFD[3207] = 1.99264525E-02;
    COFD[3208] = -2.22623522E+01;
    COFD[3209] = 5.18930931E+00;
    COFD[3210] = -4.07143570E-01;
    COFD[3211] = 1.55986909E-02;
    COFD[3212] = -2.25093649E+01;
    COFD[3213] = 5.43187491E+00;
    COFD[3214] = -4.49258240E-01;
    COFD[3215] = 1.78460456E-02;
    COFD[3216] = -2.22607146E+01;
    COFD[3217] = 5.20705035E+00;
    COFD[3218] = -4.10106772E-01;
    COFD[3219] = 1.57532724E-02;
    COFD[3220] = -1.96907838E+01;
    COFD[3221] = 3.86787656E+00;
    COFD[3222] = -2.00935789E-01;
    COFD[3223] = 5.30029143E-03;
    COFD[3224] = -2.23601792E+01;
    COFD[3225] = 5.11898303E+00;
    COFD[3226] = -3.95452058E-01;
    COFD[3227] = 1.49902966E-02;
    COFD[3228] = -1.94513324E+01;
    COFD[3229] = 3.52672843E+00;
    COFD[3230] = -1.49232080E-01;
    COFD[3231] = 2.77504261E-03;
    COFD[3232] = -1.93408889E+01;
    COFD[3233] = 3.42798339E+00;
    COFD[3234] = -1.35210942E-01;
    COFD[3235] = 2.12316045E-03;
    COFD[3236] = -1.95094733E+01;
    COFD[3237] = 3.50834090E+00;
    COFD[3238] = -1.47169195E-01;
    COFD[3239] = 2.69999047E-03;
    COFD[3240] = -1.70798825E+01;
    COFD[3241] = 2.40238329E+00;
    COFD[3242] = 1.33288432E-02;
    COFD[3243] = -4.88004423E-03;
    COFD[3244] = -2.17720676E+01;
    COFD[3245] = 5.41516167E+00;
    COFD[3246] = -4.73408139E-01;
    COFD[3247] = 1.99387805E-02;
    COFD[3248] = -2.26156323E+01;
    COFD[3249] = 5.57760228E+00;
    COFD[3250] = -4.89797685E-01;
    COFD[3251] = 2.04643263E-02;
    COFD[3252] = -1.85864144E+01;
    COFD[3253] = 4.54915847E+00;
    COFD[3254] = -3.75000738E-01;
    COFD[3255] = 1.62324821E-02;
    COFD[3256] = -1.59327297E+01;
    COFD[3257] = 3.65620899E+00;
    COFD[3258] = -2.62933804E-01;
    COFD[3259] = 1.15253223E-02;
    COFD[3260] = -1.16906297E+01;
    COFD[3261] = 2.47469981E+00;
    COFD[3262] = -1.10436257E-01;
    COFD[3263] = 4.95273813E-03;
    COFD[3264] = -1.87483158E+01;
    COFD[3265] = 4.48550694E+00;
    COFD[3266] = -3.67277454E-01;
    COFD[3267] = 1.59194755E-02;
    COFD[3268] = -1.98610390E+01;
    COFD[3269] = 4.84231878E+00;
    COFD[3270] = -4.10101001E-01;
    COFD[3271] = 1.76356687E-02;
    COFD[3272] = -1.98610390E+01;
    COFD[3273] = 4.84231878E+00;
    COFD[3274] = -4.10101001E-01;
    COFD[3275] = 1.76356687E-02;
    COFD[3276] = -2.08547637E+01;
    COFD[3277] = 5.15581320E+00;
    COFD[3278] = -4.46543388E-01;
    COFD[3279] = 1.90458118E-02;
    COFD[3280] = -2.11737258E+01;
    COFD[3281] = 5.22581725E+00;
    COFD[3282] = -4.54351699E-01;
    COFD[3283] = 1.93332012E-02;
    COFD[3284] = -1.42894441E+01;
    COFD[3285] = 3.67490723E+00;
    COFD[3286] = -2.65114792E-01;
    COFD[3287] = 1.16092671E-02;
    COFD[3288] = -1.40756935E+01;
    COFD[3289] = 3.07549274E+00;
    COFD[3290] = -1.88889344E-01;
    COFD[3291] = 8.37152866E-03;
    COFD[3292] = -1.40949196E+01;
    COFD[3293] = 3.07549274E+00;
    COFD[3294] = -1.88889344E-01;
    COFD[3295] = 8.37152866E-03;
    COFD[3296] = -1.52486273E+01;
    COFD[3297] = 3.35922578E+00;
    COFD[3298] = -2.25181399E-01;
    COFD[3299] = 9.92132878E-03;
    COFD[3300] = -2.10643259E+01;
    COFD[3301] = 5.53614847E+00;
    COFD[3302] = -4.86046736E-01;
    COFD[3303] = 2.03659188E-02;
    COFD[3304] = -1.52554761E+01;
    COFD[3305] = 3.35922578E+00;
    COFD[3306] = -2.25181399E-01;
    COFD[3307] = 9.92132878E-03;
    COFD[3308] = -1.52414485E+01;
    COFD[3309] = 3.35922578E+00;
    COFD[3310] = -2.25181399E-01;
    COFD[3311] = 9.92132878E-03;
    COFD[3312] = -1.59633387E+01;
    COFD[3313] = 3.66853818E+00;
    COFD[3314] = -2.64346221E-01;
    COFD[3315] = 1.15784613E-02;
    COFD[3316] = -2.04833713E+01;
    COFD[3317] = 5.23112374E+00;
    COFD[3318] = -4.54967682E-01;
    COFD[3319] = 1.93570423E-02;
    COFD[3320] = -1.50031687E+01;
    COFD[3321] = 3.26223357E+00;
    COFD[3322] = -2.12746642E-01;
    COFD[3323] = 9.38912883E-03;
    COFD[3324] = -1.81432461E+01;
    COFD[3325] = 4.37565431E+00;
    COFD[3326] = -3.53906025E-01;
    COFD[3327] = 1.53760786E-02;
    COFD[3328] = -1.76002031E+01;
    COFD[3329] = 4.19171952E+00;
    COFD[3330] = -3.31354810E-01;
    COFD[3331] = 1.44520623E-02;
    COFD[3332] = -1.83249299E+01;
    COFD[3333] = 4.42045763E+00;
    COFD[3334] = -3.59451578E-01;
    COFD[3335] = 1.56056164E-02;
    COFD[3336] = -2.02646611E+01;
    COFD[3337] = 5.10426133E+00;
    COFD[3338] = -4.41256919E-01;
    COFD[3339] = 1.88737290E-02;
    COFD[3340] = -1.87434358E+01;
    COFD[3341] = 4.48550694E+00;
    COFD[3342] = -3.67277454E-01;
    COFD[3343] = 1.59194755E-02;
    COFD[3344] = -2.08530414E+01;
    COFD[3345] = 5.15581320E+00;
    COFD[3346] = -4.46543388E-01;
    COFD[3347] = 1.90458118E-02;
    COFD[3348] = -2.11724381E+01;
    COFD[3349] = 5.22581725E+00;
    COFD[3350] = -4.54351699E-01;
    COFD[3351] = 1.93332012E-02;
    COFD[3352] = -2.10567656E+01;
    COFD[3353] = 5.19883613E+00;
    COFD[3354] = -4.51205592E-01;
    COFD[3355] = 1.92107800E-02;
    COFD[3356] = -2.17720676E+01;
    COFD[3357] = 5.41516167E+00;
    COFD[3358] = -4.73408139E-01;
    COFD[3359] = 1.99387805E-02;
    COFD[3360] = -1.49828430E+01;
    COFD[3361] = 3.25781069E+00;
    COFD[3362] = -2.12199367E-01;
    COFD[3363] = 9.36657283E-03;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 4;
    KTDIF[1] = 10;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(double* COFTD) {
    COFTD[0] = 3.22793099E-02;
    COFTD[1] = 8.16227034E-04;
    COFTD[2] = -4.14297547E-07;
    COFTD[3] = 6.67378576E-11;
    COFTD[4] = 2.49017478E-01;
    COFTD[5] = 4.29036573E-04;
    COFTD[6] = -2.42668617E-07;
    COFTD[7] = 4.20801371E-11;
    COFTD[8] = 3.39557243E-01;
    COFTD[9] = 1.79335036E-04;
    COFTD[10] = -1.10135705E-07;
    COFTD[11] = 2.06427239E-11;
    COFTD[12] = 0.00000000E+00;
    COFTD[13] = 0.00000000E+00;
    COFTD[14] = 0.00000000E+00;
    COFTD[15] = 0.00000000E+00;
    COFTD[16] = 2.72587742E-01;
    COFTD[17] = 4.31830859E-04;
    COFTD[18] = -2.45662623E-07;
    COFTD[19] = 4.27920337E-11;
    COFTD[20] = 2.11501250E-01;
    COFTD[21] = 5.48075431E-04;
    COFTD[22] = -3.01849605E-07;
    COFTD[23] = 5.13164012E-11;
    COFTD[24] = 2.11501250E-01;
    COFTD[25] = 5.48075431E-04;
    COFTD[26] = -3.01849605E-07;
    COFTD[27] = 5.13164012E-11;
    COFTD[28] = 1.50554111E-01;
    COFTD[29] = 6.51990284E-04;
    COFTD[30] = -3.48888001E-07;
    COFTD[31] = 5.81440411E-11;
    COFTD[32] = 1.35289561E-01;
    COFTD[33] = 6.79033126E-04;
    COFTD[34] = -3.60786540E-07;
    COFTD[35] = 5.98432251E-11;
    COFTD[36] = -1.44152190E-01;
    COFTD[37] = -7.99993584E-05;
    COFTD[38] = 4.89707442E-08;
    COFTD[39] = -9.14277269E-12;
    COFTD[40] = 4.06682492E-01;
    COFTD[41] = 3.84705248E-05;
    COFTD[42] = -2.54846868E-08;
    COFTD[43] = 5.86302354E-12;
    COFTD[44] = 4.12895615E-01;
    COFTD[45] = 3.90582612E-05;
    COFTD[46] = -2.58740310E-08;
    COFTD[47] = 5.95259633E-12;
    COFTD[48] = 4.28230888E-01;
    COFTD[49] = 1.20873273E-04;
    COFTD[50] = -7.70268349E-08;
    COFTD[51] = 1.52678954E-11;
    COFTD[52] = 2.27469146E-02;
    COFTD[53] = 6.73078907E-04;
    COFTD[54] = -3.40935843E-07;
    COFTD[55] = 5.48499211E-11;
    COFTD[56] = 4.29789463E-01;
    COFTD[57] = 1.21313199E-04;
    COFTD[58] = -7.73071792E-08;
    COFTD[59] = 1.53234639E-11;
    COFTD[60] = 4.26579943E-01;
    COFTD[61] = 1.20407274E-04;
    COFTD[62] = -7.67298757E-08;
    COFTD[63] = 1.52090336E-11;
    COFTD[64] = 3.31191185E-01;
    COFTD[65] = 1.81326714E-04;
    COFTD[66] = -1.11096391E-07;
    COFTD[67] = 2.07635959E-11;
    COFTD[68] = 1.22693382E-01;
    COFTD[69] = 6.21278143E-04;
    COFTD[70] = -3.29965208E-07;
    COFTD[71] = 5.47161548E-11;
    COFTD[72] = 4.30605547E-01;
    COFTD[73] = 9.35961902E-05;
    COFTD[74] = -6.03983623E-08;
    COFTD[75] = 1.23115170E-11;
    COFTD[76] = 2.93191523E-01;
    COFTD[77] = 4.01430006E-04;
    COFTD[78] = -2.30705763E-07;
    COFTD[79] = 4.05176586E-11;
    COFTD[80] = 3.05613225E-01;
    COFTD[81] = 3.24505886E-04;
    COFTD[82] = -1.89889572E-07;
    COFTD[83] = 3.38663465E-11;
    COFTD[84] = 2.74036956E-01;
    COFTD[85] = 3.96249742E-04;
    COFTD[86] = -2.26857964E-07;
    COFTD[87] = 3.97176979E-11;
    COFTD[88] = 1.59288984E-01;
    COFTD[89] = 6.02833801E-04;
    COFTD[90] = -3.24837576E-07;
    COFTD[91] = 5.43909010E-11;
    COFTD[92] = 2.71946054E-01;
    COFTD[93] = 4.30814303E-04;
    COFTD[94] = -2.45084319E-07;
    COFTD[95] = 4.26912987E-11;
    COFTD[96] = 1.50452493E-01;
    COFTD[97] = 6.51550218E-04;
    COFTD[98] = -3.48652516E-07;
    COFTD[99] = 5.81047963E-11;
    COFTD[100] = 1.35224069E-01;
    COFTD[101] = 6.78704414E-04;
    COFTD[102] = -3.60611888E-07;
    COFTD[103] = 5.98142558E-11;
    COFTD[104] = 1.41501819E-01;
    COFTD[105] = 6.76400843E-04;
    COFTD[106] = -3.60211338E-07;
    COFTD[107] = 5.98385258E-11;
    COFTD[108] = 7.24978506E-02;
    COFTD[109] = 7.61676489E-04;
    COFTD[110] = -3.93710165E-07;
    COFTD[111] = 6.41510700E-11;
    COFTD[112] = 4.31331269E-01;
    COFTD[113] = 9.20536800E-05;
    COFTD[114] = -5.94509616E-08;
    COFTD[115] = 1.21437993E-11;
    COFTD[116] = -2.27026268E-01;
    COFTD[117] = 9.23446255E-04;
    COFTD[118] = -3.86383917E-07;
    COFTD[119] = 5.52022104E-11;
    COFTD[120] = -5.08744745E-02;
    COFTD[121] = 8.54342586E-04;
    COFTD[122] = -4.15926453E-07;
    COFTD[123] = 6.53063261E-11;
    COFTD[124] = 1.05124122E-01;
    COFTD[125] = 6.50665957E-04;
    COFTD[126] = -3.42564538E-07;
    COFTD[127] = 5.64804120E-11;
    COFTD[128] = 1.44152190E-01;
    COFTD[129] = 7.99993584E-05;
    COFTD[130] = -4.89707442E-08;
    COFTD[131] = 9.14277269E-12;
    COFTD[132] = -4.05742102E-02;
    COFTD[133] = 8.66054242E-04;
    COFTD[134] = -4.24254167E-07;
    COFTD[135] = 6.68660357E-11;
    COFTD[136] = -1.04603726E-01;
    COFTD[137] = 9.18547452E-04;
    COFTD[138] = -4.33740076E-07;
    COFTD[139] = 6.68534946E-11;
    COFTD[140] = -1.04603726E-01;
    COFTD[141] = 9.18547452E-04;
    COFTD[142] = -4.33740076E-07;
    COFTD[143] = 6.68534946E-11;
    COFTD[144] = -1.56336482E-01;
    COFTD[145] = 9.43814119E-04;
    COFTD[146] = -4.29437786E-07;
    COFTD[147] = 6.47088533E-11;
    COFTD[148] = -1.67919551E-01;
    COFTD[149] = 9.47918492E-04;
    COFTD[150] = -4.27148814E-07;
    COFTD[151] = 6.39808439E-11;
    COFTD[152] = 0.00000000E+00;
    COFTD[153] = 0.00000000E+00;
    COFTD[154] = 0.00000000E+00;
    COFTD[155] = 0.00000000E+00;
    COFTD[156] = 2.35283119E-01;
    COFTD[157] = 4.65670599E-04;
    COFTD[158] = -2.60939824E-07;
    COFTD[159] = 4.49271822E-11;
    COFTD[160] = 2.37053352E-01;
    COFTD[161] = 4.69174231E-04;
    COFTD[162] = -2.62903094E-07;
    COFTD[163] = 4.52652072E-11;
    COFTD[164] = 1.80186965E-01;
    COFTD[165] = 6.02882805E-04;
    COFTD[166] = -3.27063140E-07;
    COFTD[167] = 5.50170790E-11;
    COFTD[168] = -1.74352698E-01;
    COFTD[169] = 8.62246873E-04;
    COFTD[170] = -3.79545489E-07;
    COFTD[171] = 5.60262093E-11;
    COFTD[172] = 1.80513677E-01;
    COFTD[173] = 6.03975942E-04;
    COFTD[174] = -3.27656165E-07;
    COFTD[175] = 5.51168351E-11;
    COFTD[176] = 1.79840299E-01;
    COFTD[177] = 6.01722902E-04;
    COFTD[178] = -3.26433894E-07;
    COFTD[179] = 5.49112302E-11;
    COFTD[180] = 1.00039110E-01;
    COFTD[181] = 6.50468660E-04;
    COFTD[182] = -3.41778999E-07;
    COFTD[183] = 5.62779132E-11;
    COFTD[184] = -1.61357564E-01;
    COFTD[185] = 9.05920260E-04;
    COFTD[186] = -4.07879153E-07;
    COFTD[187] = 6.10626290E-11;
    COFTD[188] = 2.00119897E-01;
    COFTD[189] = 5.64793704E-04;
    COFTD[190] = -3.09445484E-07;
    COFTD[191] = 5.24139335E-11;
    COFTD[192] = -2.00309448E-02;
    COFTD[193] = 8.50440115E-04;
    COFTD[194] = -4.21064468E-07;
    COFTD[195] = 6.67959710E-11;
    COFTD[196] = 1.63245097E-02;
    COFTD[197] = 7.90133388E-04;
    COFTD[198] = -3.98292458E-07;
    COFTD[199] = 6.38851432E-11;
    COFTD[200] = -2.72323768E-02;
    COFTD[201] = 8.39184413E-04;
    COFTD[202] = -4.13849924E-07;
    COFTD[203] = 6.54928043E-11;
    COFTD[204] = -1.41640506E-01;
    COFTD[205] = 9.21404324E-04;
    COFTD[206] = -4.23210110E-07;
    COFTD[207] = 6.41400322E-11;
    COFTD[208] = -4.05265093E-02;
    COFTD[209] = 8.65036069E-04;
    COFTD[210] = -4.23755394E-07;
    COFTD[211] = 6.67874248E-11;
    COFTD[212] = -1.56283739E-01;
    COFTD[213] = 9.43495709E-04;
    COFTD[214] = -4.29292909E-07;
    COFTD[215] = 6.46870228E-11;
    COFTD[216] = -1.67878917E-01;
    COFTD[217] = 9.47689110E-04;
    COFTD[218] = -4.27045451E-07;
    COFTD[219] = 6.39653615E-11;
    COFTD[220] = -1.64884446E-01;
    COFTD[221] = 9.51009610E-04;
    COFTD[222] = -4.29932210E-07;
    COFTD[223] = 6.45267379E-11;
    COFTD[224] = -1.96491231E-01;
    COFTD[225] = 9.44114020E-04;
    COFTD[226] = -4.13070492E-07;
    COFTD[227] = 6.07423579E-11;
    COFTD[228] = 2.01521643E-01;
    COFTD[229] = 5.62744089E-04;
    COFTD[230] = -3.08519239E-07;
    COFTD[231] = 5.22805986E-11;
}

