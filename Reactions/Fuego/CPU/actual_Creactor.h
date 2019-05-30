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

#ifdef USE_KLU 
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#include <sundials/sundials_sparse.h>
#include <sunmatrix/sunmatrix_sparse.h>
#endif


#include <AMReX_Print.H>

/**********************************/

    struct CVUserData {
#ifdef USE_KLU 
      int NNZ; 
      SUNMatrix *PS;
      realtype **Jdata = NULL;
      int **rowVals = NULL;
      int **colPtrs = NULL;
      int *indx = NULL;
      realtype *JSPSmat = NULL;
      klu_common *Common;
      klu_symbolic **Symbolic;
      klu_numeric **Numeric;
#else
      realtype **(**Jbd);
      realtype **(**P);
      sunindextype *(**pivot);
#endif
      int NEQ;
      int NCELLS;
      /* NRG */
      double *rhoX_init = NULL;
      double *rhoXsrc_ext = NULL;
      double *rYsrc       = NULL;
      /* hacks */
      bool FirstTimePrecond;
      bool reactor_cvode_initialized;

      int iDense_Creact;
      int iJac_Creact;
      int iE_Creact;
    };
    
    typedef CVUserData *UserData;

    UserData AllocUserData(int NEQ, int NCELLS, int iE, int iR, int iJ);

//    /**********************************/
//    /* Attr */
//    void setJac(int iJac) {
//	    iJac_Creact = iJac;
//    }
//
//    void setReacType(int iR) {
//	    iE_Creact = iR;
//    }
//
//    void setIntegType(int iD) {
//	    iDense_Creact = iD;
//    }
//
//    int iJac_return() const {
//	    return iJac_Creact;
//    }
//
//    int iE_return() const {
//	    return iE_Creact;
//    }
//
//    int iDense_return() const {
//	    return iDense_Creact;
//    }
    
//protected:
//    int iDense_Creact;
//    int iJac_Creact;
//    int iE_Creact;
//};
    

/**********************************/
/* Functions Called by the Solver */
 int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

 int cJac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

#ifdef USE_KLU 
 int cJac_KLU(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

 int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

#ifdef USE_KLU 
 int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);
 int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);
#endif

 int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);


/**********************************/
/* Functions Called by the Program */
int react(realtype *rY_in, realtype *rY_src_in, 
    		realtype *rX_in, realtype *rX_src_in, 
    		realtype *P_in, realtype *dt_react, realtype *time,
		const int* cvode_iE, const int* Ncells);

void react_info(const int* cvode_iE, const int* Ncells);
    
void reactor_close(void *cvodeMem, SUNLinearSolver LS, SUNMatrix A, void *data, N_Vector y, N_Vector tol);
    
    
/**********************************/
/* Helper functions */
int check_flag(void *flagvalue, const char *funcname, int opt);

void PrintFinalStats(void *cvodeMem, realtype Temp, UserData data);

void check_state(N_Vector yvec, int NEQ, int NCELLS,bool ok_to_int);


/**********************************/
/* Main Kernel fct called in solver RHS */
void fKernelSpec(int neq, int ncells, realtype *dt, realtype *yvec_d, realtype *ydot_d,
    		double *rhoX_init, double *rhoXsrc_ext, double *rYs, int iE);


