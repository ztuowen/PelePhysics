#include <math.h>
#include <iostream>
#include <cstring>
#include <chrono>


#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver     */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <cvode/cvode_spils.h>         /* access to CVSpils interface */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>
#include <cvode/cvode_impl.h>

#include <AMReX_Print.H>

/**********************************/

/* Functions Called by the Solver */
static int cF_RHS_unit(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

static int cJac_unit(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int reactor_unit_init(const int* cvode_iE); 

int react_unit(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in, 
		realtype *P_in, realtype *dt_react, realtype *time);

void reactor_unit_close();

static int check_flag_unit(void *flagvalue, const char *funcname, int opt);

static void PrintFinalStats_unit(void *cvodeMem);

extern "C" {
    void get_t_given_ey_(double * e, double * y, double * t, int *ierr); 
    void get_t_given_hy_(double * h, double * y, double * t, int *ierr); 
    void ckcvms_(double * T, double * cvms);
    void ckcpms_(double * T, double * cvms);
    void ckums_(double * T,  double * ums);
    void ckhms_(double * T,  double * ums);
    void ckwc_(double * T, double * C, double * wdot);
    void ckwt_(double * wt);
    void dwdot_(double * J, double * sc, double * Tp, int * consP);
    void ckindx_(int * mm, int * kk, int * ii, int * nfit); 

}
/**********************************/
/* Main Kernel fct called in solver RHS */
void fKernelSpec_unit(realtype *dt, realtype *yvec_d, realtype *ydot_d,
		double *rhoX_init, double *rhoXsrc_ext, double *rYs);


