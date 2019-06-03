#include <CPU/actual_Creactor.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>

/**********************************/
/* Dummy MAIN CVODE call routine to print params once */
void react_info(const int* cvode_iE, const int* Ncells){

	int mm, ii, nfit;
	int NEQ;

	/* FIRST STEP get NEQ */
	CKINDX(&mm,&NEQ,&ii,&nfit);
        printf("Nb of spec is %d \n", NEQ);

	/* Args */
	printf("Ncells in one solve is %d\n",*Ncells);
	printf("Reactor type is %d\n",*cvode_iE);

	/* ParmParse from the inputs file */ 
	/* Reuse dummy variables mm and ii */
	amrex::ParmParse pp("ns");
	pp.query("cvode_iJac",mm);
	pp.query("cvode_iDense", ii);

	if (ii == 1) {
            printf("\n--> Using a Direct Dense Solver \n");
#ifdef USE_KLU 
	} else if (ii == 5) {
            printf("\n--> Using a Sparse Dense Solver \n");
#endif
	} else if (ii == 99) {
            printf("\n--> Using an Iterative Solver \n");
	} else {
	    amrex::Abort("\n--> ... what are you doing ? \n");
	}

	if (mm == 0) {
            printf("\n--> Without Analytical J\n");
#ifdef USE_KLU 
	    if (ii == 5) {
	        amrex::Abort("\n--> SPARSE SOLVER SHOULD HAVE AN AJ...\n");
	    }
#endif
	} else {
            printf("\n--> With Analytical J ... ");
	    if (ii == 99) {
#ifdef USE_KLU 
                printf(" KLU enabled so GMRES is SPARSE \n");
	        int nJdata;
                int HP;
                if (*cvode_iE == 1) {
                    HP = 0;
                } else {
                    HP = 1;
                }
		SPARSITY_INFO_PRECOND(&nJdata,&HP);
                printf("--> SPARSE Preconditioner -- non zero entries %d represents %f %% fill pattern.\n", nJdata, nJdata/float((NEQ+1) * (NEQ+1)) *100.0);
#else
                printf(" \n");
#endif
	    } else {
                printf(" \n");
	    }
	}

	/* Ok we're done ...*/
        printf(" --> DONE WITH INITIALIZATION (CPU) \n\n");
}

/* Main CVODE call routine */
int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
		realtype *P_in, 
                realtype *dt_react, realtype *time, 
		const int* cvode_iE, const int* Ncells){

        /* declarations */
	realtype time_init, time_out, dummy_time, temperature_save ;
	realtype reltol;
	N_Vector atol = NULL;
	realtype *ratol;
	int mm, ii, nfit;
	int NEQ, NCELLS, neq_tot;
	int flag;
        int iverbose       = 2;
        /* FOR CVODE */
        N_Vector y         = NULL;
        SUNLinearSolver LS = NULL;
        SUNMatrix A        = NULL;
        void *cvode_mem    = NULL;
        UserData data      = NULL;

	/* Print initial time and expected output time */
        time_init = *time;
	time_out  = *time + (*dt_react);
        if (iverbose > 3) {
	    printf("BEG : time curr is %14.6e and dt_react is %14.6e and final time should be %14.6e \n", time_init, *dt_react, time_out);
	}

	/* FIRST STEP get NEQ */
	CKINDX(&mm,&NEQ,&ii,&nfit);

	/* ParmParse from the inputs file */ 
	amrex::ParmParse pp("ns");
	pp.query("cvode_iJac",mm);
	pp.query("cvode_iDense", ii);

	/* Args */

	NCELLS         = *Ncells;
        neq_tot        = (NEQ + 1) * NCELLS;

	/* Definition of main vector */
	y = N_VNew_Serial(neq_tot);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	cvode_mem = CVodeCreate(CV_BDF);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

        /* Does not work for more than 1 cell right now */
	data = AllocUserData(NEQ, NCELLS, *cvode_iE, ii, mm);
	if(check_flag((void *)data, "AllocUserData", 2)) return(1);

	/* Set the pointer to user-defined data */
	flag = CVodeSetUserData(cvode_mem, data);
	if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);   
	/* Define vectors to be used later */
	data->rhoX_init = (double *) malloc(NCELLS*sizeof(double));
	data->rhoXsrc_ext = (double *) malloc( NCELLS*sizeof(double));
	data->rYsrc = (double *)  malloc((NCELLS*NEQ)*sizeof(double));

	/* Fill */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetArrayPointer(y);
	// rhoY,T
	std::memcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS));
	temperature_save = rY_in[NEQ];
	// rhoY_src_ext
	std::memcpy(data->rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double));
	// rhoE/rhoH
	std::memcpy(data->rhoX_init, rX_in, sizeof(realtype) * NCELLS);
	std::memcpy(data->rhoXsrc_ext, rX_src_in, sizeof(realtype) * NCELLS);

	/* check if state is ready to integrate */
        bool actual_ok_to_react = true;
	check_state(y, NEQ, NCELLS, actual_ok_to_react);
	if (!actual_ok_to_react)  { 
	    amrex::Abort("\n Check_state failed: state is out of react bounds \n");
	}

        realtype time_tmp = time_init; //0.0e+0;
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function, the inital time, and 
	 * initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, cF_RHS, time_tmp, y);
	if (check_flag(&flag, "CVodeInit", 1)){
		amrex::Abort("CVode init error");
	}
	
	/* Definition of tolerances: one for each species */
	reltol = 1.0e-10;
        atol  = N_VNew_Serial(neq_tot);
	ratol = N_VGetArrayPointer(atol);
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = 1.0e-10;
        }
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, atol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	flag = CVodeSetInitStep(cvode_mem, 1.0e-09);
	if (check_flag(&flag, "CVodeSetInitStep", 1)) return(1);

	//flag = CVodeSetNonlinConvCoef(cvode_mem, 1.0e-03);
	//if (check_flag(&flag, "CVodeSetNonlinConvCoef", 1)) return(1);

	flag = CVodeSetMaxNonlinIters(cvode_mem, 50);
	if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

	flag = CVodeSetMaxErrTestFails(cvode_mem, 100);
	if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return(1);

	if (data->iDense_Creact == 1) {
            /* Create dense SUNMatrix for use in linear solves */
	    A = SUNDenseMatrix(neq_tot, neq_tot);
	    if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

	    /* Create dense SUNLinearSolver object for use by CVode */
	    LS = SUNDenseLinearSolver(y, A);
	    if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	    if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

#ifdef USE_KLU 
	} else if (data->iDense_Creact == 5) {

	    /* Create sparse SUNMatrix for use in linear solves */
	    A = SUNSparseMatrix(neq_tot, neq_tot, (data->NNZ)*NCELLS, CSC_MAT);
            if(check_flag((void *)A, "SUNSparseMatrix", 0)) return(1);

	    /* Create KLU solver object for use by CVode */
	    LS = SUNLinSol_KLU(y, A);
	    if(check_flag((void *)LS, "SUNLinSol_KLU", 0)) return(1);

	    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
	    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
	    if(check_flag(&flag, "CVodeSetLinearSolver", 1)) return(1);
#endif

	} else if (data->iDense_Creact == 99) {

            /* Create the linear solver object */
	    if (data->iJac_Creact == 0) {
	        LS = SUNSPGMR(y, PREC_NONE, 0);
	        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
	    } else {
	        LS = SUNSPGMR(y, PREC_LEFT, 0);
	        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
	    }

	    /* Set CVSpils linear solver to LS */
	    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
	    if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);
	} else {
	    int nJdata;
	    int HP = 0;
	    SPARSITY_INFO(&nJdata, &HP, 1);
            printf("--> SPARSE solver -- non zero entries %d represents %f %% sparsity pattern.", nJdata, nJdata/float((NEQ+1) * (NEQ+1)) *100.0);
	    amrex::Abort("\n \n");
	}

	if (data->iJac_Creact == 1) {
	    if (data->iDense_Creact == 99) {
	        /* Set the JAcobian-times-vector function */
	        flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
	        if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);
	        /* Set the preconditioner solve and setup functions */
#ifdef USE_KLU 
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond_sparse, PSolve_sparse);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#else
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#endif
#ifdef USE_KLU 
	    } else if (data->iDense_Creact == 5){
		/* Set the user-supplied Jacobian routine Jac */
		flag = CVodeSetJacFn(cvode_mem, cJac_KLU);
		if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1); 
#endif
	    } else {
	        /* Set the user-supplied Jacobian routine Jac */
                flag = CVodeSetJacFn(cvode_mem, cJac);
		if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1);
	    }
	}

        /* Set the max number of time steps */ 
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

        /* Set the max order */ 
        //flag = CVodeSetMaxOrd(cvode_mem, 2);
	//if(check_flag(&flag, "CVodeSetMaxOrd", 1)) return(1);

	data->reactor_cvode_initialized = true;

        /* Actual Prgm */

	flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
	/* ONE STEP MODE FOR DEBUGGING */
	//flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_ONE_STEP);
	if (check_flag(&flag, "CVode", 1)){
		amrex::Abort("CVode integration error");
	}

	*dt_react = dummy_time - time_init;
	*time = time_init + *dt_react;
	if (iverbose > 3) {
	    printf("END : time curr is %14.6e and dt_react is %14.6e \n", dummy_time, *dt_react);
	}

	/* Pack data to return in main routine external */
	std::memcpy(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype));
	for  (int i = 0; i < NCELLS; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}

	/* VERBOSE MODE */
	if (iverbose > 10) {
            for (int tid = 0; tid < NCELLS; tid ++) {
	        double rhov, energy, temp, energy2;
	        double MF[NEQ];
                //realtype activity[NEQ], cdot[NEQ], molecular_weight[NEQ];
	        int  lierr;
	        rhov = 0.0;
                int offset = tid * (NEQ + 1); 
                for (int k = 0; k < NEQ; k ++) {
	    	rhov =  rhov + rY_in[offset + k];
	        }
                //CKWT(molecular_weight);
                for (int k = 0; k < NEQ; k ++) {
	    	    MF[k] = rY_in[offset + k]/rhov;
	            //activity[k] = rY_in[offset + k]/(molecular_weight[k]);
	        }
	        energy = rX_in[tid]/rhov ;
	        if (*cvode_iE == 1) { 
	            GET_T_GIVEN_EY(&energy, MF, &temp, &lierr);
	            CKHBMS(&temp, MF, &energy2);
	            CKUBMS(&temp, MF, &energy);
	    	    CKPY(&rhov, &temp, MF, P_in);
	            //printf("T,e,h,p,rho 'I'? 
		    printf("%4.16e %4.16e %4.16e %4.16e %4.16e %4.16e %4.16e %4.16e %4.16e %4.16e %4.16e \n",*time, rY_in[offset + NEQ],energy, energy2, *P_in, rhov, MF[0], MF[1], MF[5], MF[3], MF[2]);
	        } else {
	            GET_T_GIVEN_HY(&energy, MF, &temp, &lierr);
	            CKHBMS(&temp, MF, &energy);
	            CKUBMS(&temp, MF, &energy2);
	    	    CKPY(&rhov, &temp, MF, P_in);
	            printf("e,h,p,rho 'II'? %4.16e %4.16e %4.16e %4.16e\n",energy2, energy, *P_in, rhov);
	        }
	    }
	    printf("\nAdditional verbose info --\n");
	    PrintFinalStats(cvode_mem, temperature_save, data);
	}

	// Clean up
        reactor_close(cvode_mem, LS, A, data, y, atol);

        long int nfe;
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	return nfe;
}


/* RHS routine */
 int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
		void *user_data){

        /* Make local copies of pointers in user_data (cell M)*/
        UserData data_wk;
        data_wk = (UserData) user_data;   

	realtype *y_d      = N_VGetArrayPointer(y_in);
	realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

	fKernelSpec(data_wk->NEQ, data_wk->NCELLS, &t, y_d, ydot_d, 
			    data_wk->rhoX_init, data_wk->rhoXsrc_ext, data_wk->rYsrc,
			    data_wk->iE_Creact);

	return(0);
}


/*
 * kernels
 */

/* RHS source terms evaluation */
void fKernelSpec(int NEQ, int NCELLS, realtype *dt, realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs, int iE)
{
  int tid;
  for (tid = 0; tid < NCELLS; tid ++) {
      realtype massfrac[NEQ],activity[NEQ];
      realtype Xi[NEQ], cXi[NEQ];
      realtype cdot[NEQ], molecular_weight[NEQ];
      realtype temp, energy;
      int lierr;

      int offset = tid * (NEQ + 1); 

      /* MW CGS */
      CKWT(molecular_weight);
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NEQ; i++){
          rho = rho + yvec_d[offset + i];
      }
      /* temp */
      temp = yvec_d[offset + NEQ];
      /* Yks, C CGS*/
      for (int i = 0; i < NEQ; i++){
          massfrac[i] = yvec_d[offset + i] / rho;
	  activity[i] = yvec_d[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      energy = (rhoX_init[tid] + rhoXsrc_ext[tid]*(*dt)) /rho;

      /* Fuego calls on device */
      if (iE == 1){
          GET_T_GIVEN_EY(&energy, massfrac, &temp, &lierr);
          CKUMS(&temp, Xi);
          CKCVMS(&temp, cXi);
      } else {
          GET_T_GIVEN_HY(&energy, massfrac, &temp, &lierr);
          CKHMS(&temp, Xi);
          CKCPMS(&temp, cXi);
      }
      CKWC(&temp, activity, cdot);
      int cX = 0.0;
      for (int i = 0; i < NEQ; i++){
          cX = cX + massfrac[i] * cXi[i];
      }

      /* Fill ydot vect */
      ydot_d[offset + NEQ] = rhoXsrc_ext[tid];
      for (int i = 0; i < NEQ; i++){
          ydot_d[offset + i] = cdot[i] * molecular_weight[i] + rYs[tid * (NEQ) + i];
          ydot_d[offset + NEQ] = ydot_d[offset + NEQ]  - ydot_d[offset + i] * Xi[i];
      }
      ydot_d[offset + NEQ] = ydot_d[offset + NEQ] /(rho * cX);
  }
}


/* Analytical Jacobian evaluation */
 int cJac(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk;
  data_wk = (UserData) user_data;   
  realtype *ydata  = N_VGetArrayPointer(u);

  int tid;
  for (tid = 0; tid < data_wk->NCELLS; tid ++) {
      realtype *J_col_k;
      realtype  temp; 
      realtype activity[data_wk->NEQ], molecular_weight[data_wk->NEQ];
      realtype Jmat_tmp[(data_wk->NEQ+1)*(data_wk->NEQ+1)];

      int offset = tid * (data_wk->NEQ + 1); 

      /* MW CGS */
      CKWT(molecular_weight);
      /* temp */
      temp = ydata[offset + data_wk->NEQ];
      for (int i = 0; i < data_wk->NEQ; i++){
          activity[i] = ydata[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      int consP;
      if (data_wk->iE_Creact == 1) {
	  consP = 0;
          DWDOT(Jmat_tmp, activity, &temp, &consP);
      } else {
          consP = 1;
          DWDOT(Jmat_tmp, activity, &temp, &consP);
      }
      /* fill the sunMat */
      for (int k = 0; k < data_wk->NEQ; k++){
	  J_col_k = SM_COLUMN_D(J,offset + k);
	  for (int i = 0; i < data_wk->NEQ; i++){
	        J_col_k[offset + i] = Jmat_tmp[k*(data_wk->NEQ+1)+i] * molecular_weight[i] / molecular_weight[k]; 
          }
	  J_col_k[offset + data_wk->NEQ] = Jmat_tmp[k*(data_wk->NEQ+1)+data_wk->NEQ] / molecular_weight[k]; 
      }
      J_col_k = SM_COLUMN_D(J,offset + data_wk->NEQ);
      for (int i = 0; i < data_wk->NEQ; i++){
          J_col_k[offset + i] = Jmat_tmp[data_wk->NEQ*(data_wk->NEQ+1)+i] * molecular_weight[i]; 
      }
  }

  return(0);

}


#ifdef USE_KLU 
/* Analytical SPARSE Jacobian evaluation */
 int cJac_KLU(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  /* Make local copies of pointers to input data (big M) */
  realtype *ydata  = N_VGetArrayPointer(u);
  sunindextype *colptrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype *rowvals_tmp = SUNSparseMatrix_IndexValues(J);
  realtype *data = SUNSparseMatrix_Data(J);

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk;
  data_wk = (UserData) user_data;   

  /* MW CGS */
  realtype molecular_weight[data_wk->NEQ];
  CKWT(molecular_weight);

  /* Fixed RowVals */
  for (int i=0;i<data_wk->NNZ;i++) {
      rowvals_tmp[i] = data_wk->rowVals[0][i];
  }
  /* Fixed colPtrs */
  colptrs_tmp[0] = data_wk->colPtrs[0][0];
  for (int i=0;i<data_wk->NCELLS*(data_wk->NEQ + 1);i++) {
      colptrs_tmp[i+1] = data_wk->colPtrs[0][i+1];
  }

  /* Temp vectors */
  realtype temp_save_lcl, temp;
  realtype activity[data_wk->NEQ];
  realtype Jmat_tmp[(data_wk->NEQ+1)*(data_wk->NEQ+1)];
  int tid, offset, nbVals, idx;
  temp_save_lcl = 0.0;
  for (tid = 0; tid < data_wk->NCELLS; tid ++) {
      offset = tid * (data_wk->NEQ + 1); 
      /* temp */
      temp = ydata[offset + data_wk->NEQ];
      /* Do we recompute Jac ? */
      if (fabs(temp - temp_save_lcl) > 1.0) {
          for (int i = 0; i < data_wk->NEQ; i++){
              activity[i] = ydata[offset + i]/(molecular_weight[i]);
          }
          /* NRG CGS */
          int consP;
          if (data_wk->iE_Creact == 1) {
              consP = 0;
              DWDOT(Jmat_tmp, activity, &temp, &consP);
          } else {
              consP = 1;
              DWDOT(Jmat_tmp, activity, &temp, &consP);
          }
	  temp_save_lcl = temp;
	  /* rescale */
          for (int i = 0; i < data_wk->NEQ; i++) {
              for (int k = 0; k < data_wk->NEQ; k++) {
                  Jmat_tmp[k*(data_wk->NEQ+1) + i] = Jmat_tmp[k*(data_wk->NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
              }
              Jmat_tmp[i*(data_wk->NEQ+1) + data_wk->NEQ] = Jmat_tmp[i*(data_wk->NEQ+1) + data_wk->NEQ] / molecular_weight[i];
          }
          for (int i = 0; i < data_wk->NEQ; i++) {
              Jmat_tmp[data_wk->NEQ*(data_wk->NEQ+1) + i] = Jmat_tmp[data_wk->NEQ*(data_wk->NEQ+1) + i] * molecular_weight[i];
          }
      }
      /* Go from Dense to Sparse */
      for (int i = 1; i < data_wk->NEQ+2; i++) {
	  nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i - 1];
	  for (int j = 0; j < nbVals; j++) {
	          idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i - 1] + j ];
	              data[ data_wk->colPtrs[0][offset + i - 1] + j ] = Jmat_tmp[(i - 1) * (data_wk->NEQ + 1) + idx];
	  }
      }
  }

  return(0);

}
#endif


#ifdef USE_KLU 
/* Preconditioner setup routine for GMRES solver when KLU sparse mode is activated 
 * Generate and preprocess P
*/
 int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{

  int ok,tid;
  /* Make local copies of pointers in user_data, and of pointer to u's data */
  UserData data_wk;
  data_wk = (UserData) user_data;   
  //realtype **(**Jbd);
  //Jbd = (data_wk->Jbd);
  realtype *udata;
  udata = N_VGetArrayPointer(u);


  /* MW CGS */
  realtype molecular_weight[data_wk->NEQ];
  CKWT(molecular_weight);

  /* Formalism */
  int consP;
  if (data_wk->iE_Creact == 1) { 
       consP = 0;
  } else {
       consP = 1;
  }

  //start_Jcomp = std::chrono::system_clock::now();
  if (jok) {
        /* jok = SUNTRUE: Copy Jbd to P */
        *jcurPtr = SUNFALSE;
  } else {
        /* Temp vectors */
        realtype temp, temp_save_lcl;
        realtype activity[data_wk->NEQ];
	int offset;
        temp_save_lcl = 0.0;
        for (tid = 0; tid < data_wk->NCELLS; tid ++) {
            offset = tid * (data_wk->NEQ + 1); 
            /* temp */
            temp = udata[offset + data_wk->NEQ];
            /* Do we recompute Jac ? */
            //if (fabs(temp - temp_save_lcl) > 1.0) {
                for (int i = 0; i < data_wk->NEQ; i++){
		    activity[i] = udata[offset + i]/(molecular_weight[i]);
                }
		SLJ_PRECOND_CSC(data_wk->JSPSmat,data_wk->indx,&(data_wk->NNZ), activity, &temp, &consP, &gamma);
	        //temp_save_lcl = temp;
	    //} else {
	    //    /* if not: copy the one from prev cell */
	    //    for (int k = 1; k < data_wk->NNZ; k++) { 
	    //        data_wk->Jdata[tid][k] = data_wk->Jdata[tid-1][k];
	    //    }
	    //}
	}

        *jcurPtr = SUNTRUE;
  }

  int nbVals;
  for (int i = 1; i < data_wk->NEQ+2; i++) {
      /* nb non zeros elem should be the same for all cells */
      nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i-1];
      for (int j = 0; j < nbVals; j++) {
    	  /* row of non zero elem should be the same for all cells */
    	  int idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i-1] + j ];
          /* Scale by -gamma */
          /* Add identity matrix */
          for (tid = 0; tid < data_wk->NCELLS; tid ++) {
    	      if (idx == (i-1)) {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = 1.0 - gamma * data_wk->JSPSmat[ data_wk->colPtrs[tid][i-1] + j ]; 
    	      } else {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = - gamma * data_wk->JSPSmat[ data_wk->colPtrs[tid][i-1] + j ]; 
    	      }
          }
      }
  }
  
  if (!(data_wk->FirstTimePrecond)) {
      for (tid = 0; tid < data_wk->NCELLS; tid ++) {
          ok = klu_refactor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], data_wk->Numeric[tid], &(data_wk->Common[tid]));
      }
  } else {
      for (tid = 0; tid < data_wk->NCELLS; tid ++) {
          data_wk->Numeric[tid] = klu_factor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], &(data_wk->Common[tid])) ; 
      }
      data_wk->FirstTimePrecond = false;
  }

  return(0);
}
#else


/* Preconditioner setup routine for GMRES solver when no sparse mode is activated 
 * Generate and preprocess P
*/
 int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{

  UserData data_wk;
  data_wk = (UserData) user_data;   
  realtype **(**P), **(**Jbd);
  sunindextype *(**pivot), ierr;
  realtype *udata; //, **a, **j;
  realtype temp;
  realtype Jmat[(data_wk->NEQ+1)*(data_wk->NEQ+1)];
  realtype activity[data_wk->NEQ];

  /* MW CGS */
  realtype molecular_weight[data_wk->NEQ];
  CKWT(molecular_weight);

  /* Make local copies of pointers in user_data, and of pointer to u's data */
  P = (data_wk->P);
  Jbd = (data_wk->Jbd);
  pivot = (data_wk->pivot);
  udata = N_VGetArrayPointer(u);

  if (jok) {
      /* jok = SUNTRUE: Copy Jbd to P */
      denseCopy(Jbd[0][0], P[0][0], data_wk->NEQ+1, data_wk->NEQ+1);
      *jcurPtr = SUNFALSE;
  } else {
      /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */
      /* Make local copies of problem variables, for efficiency. */
      for (int i = 0; i < data_wk->NEQ; i++){
          activity[i] = udata[i]/(molecular_weight[i]);
      }
      temp = udata[data_wk->NEQ];

      // C in mol/cm3
      int consP;
      if (data_wk->iE_Creact == 1) { 
          consP = 0;
      } else {
          consP = 1;
      }
      DWDOT_PRECOND(Jmat, activity, &temp, &consP);
      //DWDOT(Jmat, activity, &temp, &consP);
      /* Compute Jacobian.  Load into P. */
      denseScale(0.0, Jbd[0][0], data_wk->NEQ+1, data_wk->NEQ+1);
      for (int i = 0; i < data_wk->NEQ; i++) {
          for (int k = 0; k < data_wk->NEQ; k++) {
              (Jbd[0][0])[k][i] = Jmat[k*(data_wk->NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
          }
          (Jbd[0][0])[i][data_wk->NEQ] = Jmat[i*(data_wk->NEQ+1) + data_wk->NEQ] / molecular_weight[i];
      }
      for (int i = 0; i < data_wk->NEQ; i++) {
          (Jbd[0][0])[data_wk->NEQ][i] = Jmat[data_wk->NEQ*(data_wk->NEQ+1) + i] * molecular_weight[i];
      }
      (Jbd[0][0])[data_wk->NEQ][data_wk->NEQ] = Jmat[(data_wk->NEQ+1)*(data_wk->NEQ+1)-1];

      denseCopy(Jbd[0][0], P[0][0], data_wk->NEQ+1, data_wk->NEQ+1);

      *jcurPtr = SUNTRUE;
  }
  
  /* Scale by -gamma */
  denseScale(-gamma, P[0][0], data_wk->NEQ+1, data_wk->NEQ+1);
  //denseScale(0.0, P[0][0], data_wk->NEQ+1, data_wk->NEQ+1);
  
  /* Add identity matrix and do LU decompositions on blocks in place. */
  denseAddIdentity(P[0][0], data_wk->NEQ+1);
  ierr = denseGETRF(P[0][0], data_wk->NEQ+1, data_wk->NEQ+1, pivot[0][0]);
  if (ierr != 0) return(1);
  
  return(0);
}
#endif


#ifdef USE_KLU 
/* PSolve for GMRES solver when KLU sparse mode is activated */
 int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{

  UserData data_wk;
  data_wk = (UserData) user_data;

  realtype *zdata;
  zdata = N_VGetArrayPointer(z);

  N_VScale(1.0, r, z);

  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  int tid, offset_beg, offset_end;
  realtype zdata_cell[data_wk->NEQ+1];
  for (tid = 0; tid < data_wk->NCELLS; tid ++) {
      offset_beg = tid * (data_wk->NEQ + 1); 
      offset_end = (tid + 1) * (data_wk->NEQ + 1);
      std::memcpy(zdata_cell, zdata+offset_beg, (data_wk->NEQ+1)*sizeof(realtype));
      klu_solve(data_wk->Symbolic[tid], data_wk->Numeric[tid], data_wk->NEQ+1, 1, zdata_cell, &(data_wk->Common[tid])) ; 
      std::memcpy(zdata+offset_beg, zdata_cell, (data_wk->NEQ+1)*sizeof(realtype));
  }

  return(0);
}
#else


/* PSolve for GMRES solver when no sparse mode is activated */
 int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  realtype **(**P);
  sunindextype *(**pivot);
  realtype *zdata, *v;
  UserData data_wk;

  /* Extract the P and pivot arrays from user_data. */

  data_wk = (UserData) user_data;
  P = data_wk->P;
  pivot = data_wk->pivot;
  zdata = N_VGetArrayPointer(z);
  
  N_VScale(1.0, r, z);
  
  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  v = zdata;
  denseGETRS(P[0][0], data_wk->NEQ+1, pivot[0][0], v);

  return(0);
}
#endif


/* 
 * OTHERS
*/


void check_state(N_Vector yvec, int NEQ, int NCELLS, bool actual_ok_to_react) 
{
  realtype *ydata;
  ydata = N_VGetArrayPointer(yvec);

  double rho, Temp;
  int offset;
  for (int tid = 0; tid < NCELLS; tid ++) {
      rho = 0.0;
      offset = tid * (NEQ + 1); 
      for (int k = 0; k < NEQ; k ++) {
          rho =  rho + ydata[offset + k];
          //printf("yspec ? %3.16e \n", ydata[offset + k]);
      }
      //printf("rho ? %3.16e \n", rho);
      Temp = ydata[offset + NEQ];
      if ((rho < 1.0e-10) || (rho > 1.e10)) {
          actual_ok_to_react = false;
      }
      if ((Temp < 200.0) || (Temp > 5000.0)) {
          actual_ok_to_react = false; 
      }
  }

}

/* Get and print some final statistics */
void PrintFinalStats(void *cvodeMem, realtype Temp, UserData data_wk)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn;
  long int nli, npe, nps, ncfl, netfails;
  int flag;
  realtype hlast, hinused, hcur;

  //data_wk = (UserData) data;   

  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netfails);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetLastStep(cvodeMem, &hlast);
  check_flag(&flag, "CVodeGetLastStep", 1);
  flag = CVodeGetActualInitStep(cvodeMem, &hinused);
  check_flag(&flag, "CVodeGetActualInitStep", 1);
  flag = CVodeGetCurrentTime(cvodeMem, &hcur);
  check_flag(&flag, "CVodeGetCurrentTime", 1);

  if (data_wk->iDense_Creact == 1){
      flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
      check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
      flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
      check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  } else if (data_wk->iDense_Creact == 99){
      flag = CVSpilsGetNumRhsEvals(cvodeMem, &nfeLS);
      check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);
      flag = CVSpilsGetNumJtimesEvals(cvodeMem, &nje);
      check_flag(&flag, "CVSpilsGetNumJTSetupEvals", 1);
      flag = CVSpilsGetNumPrecEvals(cvodeMem, &npe);
      check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
      flag = CVSpilsGetNumPrecSolves(cvodeMem, &nps);
      check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
      flag = CVSpilsGetNumLinIters(cvodeMem, &nli);
      check_flag(&flag, "CVSpilsGetNumLinIters", 1);
      flag = CVSpilsGetNumConvFails(cvodeMem, &ncfl); 
      check_flag(&flag, "CVSpilsGetNumConvFails", 1);
  }

  printf("-- Final Statistics --\n");
  printf("NonLinear (Newton) related --\n");
  printf("    DT(dt, dtcur), RHS, Iterations, ErrTestFails, LinSolvSetups = %f %-6ld(%14.6e %14.6e) %-6ld %-6ld %-6ld %-6ld \n",
  Temp, nst, hlast, hcur, nfe, nni, netfails, nsetups);

  if (data_wk->iDense_Creact == 1){
      printf("Linear (Dense Direct Solve) related --\n");
      printf("    FD RHS, NumJacEvals                           = %f %-6ld %-6ld \n", Temp, nfeLS, nje);
  } else if (data_wk->iDense_Creact == 99){
      // LinSolvSetups actually reflects the number of time the LinSolver has been called. 
      // NonLinIterations can be taken without the need for LinItes
      printf("Linear (Krylov GMRES Solve) related --\n");
      printf("    RHSeval, jtvEval, NumPrecEvals, NumPrecSolves = %f %-6ld %-6ld %-6ld %-6ld \n", 
            	      Temp, nfeLS, nje, npe, nps);
      printf("    Iterations, ConvFails = %f %-6ld %-6ld \n", 
         	      Temp, nli, ncfl );
  }
}


/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}


/* Alloc Data for CVODE */
UserData AllocUserData(int neq, int ncells, int iE, int iR, int iJ)
{
  UserData data_wk;
  data_wk = (UserData) malloc(sizeof *data_wk);

 (data_wk->NEQ)    = neq;
 (data_wk->NCELLS) = ncells;  
 (data_wk->FirstTimePrecond) = true;
 (data_wk->reactor_cvode_initialized) = false;

 (data_wk->iE_Creact)      = iE;
 (data_wk->iDense_Creact)  = iR;
 (data_wk->iJac_Creact)    = iJ;

#ifndef USE_KLU 
if (data_wk->iDense_Creact == 99) {
      /* Precond data */
      (data_wk->P) = new realtype***[ncells];
      (data_wk->Jbd) = new realtype***[ncells];
      (data_wk->pivot) = new sunindextype**[ncells];
      for(int i = 0; i < ncells; ++i) {
              (data_wk->P)[i] = new realtype**[ncells];
              (data_wk->Jbd)[i] = new realtype**[ncells];
              (data_wk->pivot)[i] = new sunindextype*[ncells];
      }

      for(int i = 0; i < ncells; ++i) {
          (data_wk->P)[i][i] = newDenseMat(neq+1, neq+1);
          (data_wk->Jbd)[i][i] = newDenseMat(neq+1, neq+1);
          (data_wk->pivot)[i][i] = newIndexArray(neq+1);
      }
 } 
#else 
  /* Sparse Direct and Sparse (It) Precond data */
  data_wk->colPtrs = new int*[ncells];
  data_wk->rowVals = new int*[ncells];
  data_wk->Jdata = new realtype*[ncells];

  int HP;
  if (data_wk->iE_Creact == 1) {
      HP = 0;
  } else {
      HP = 1;
  }
  if (data_wk->iDense_Creact == 5) {
      /* Sparse Matrix for Direct Sparse KLU solver */
      (data_wk->PS) = new SUNMatrix[1];
      SPARSITY_INFO(&(data_wk->NNZ),&HP,ncells);
      //printf("--> SPARSE solver -- non zero entries %d represents %f %% fill pattern.\n", data_wk->NNZ, data_wk->NNZ/float((neq+1) * (neq+1) * ncells * ncells) *100.0);
          (data_wk->PS)[0] = SUNSparseMatrix((neq+1)*ncells, (neq+1)*ncells, data_wk->NNZ, CSC_MAT);
          data_wk->colPtrs[0] = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[0]); 
          data_wk->rowVals[0] = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[0]);
          data_wk->Jdata[0] = SUNSparseMatrix_Data((data_wk->PS)[0]);
          SPARSITY_PREPROC(data_wk->rowVals[0],data_wk->colPtrs[0],&HP,ncells);
      //}

  } else if (data_wk->iDense_Creact == 99) {
      /* KLU internal storage */
      data_wk->Common = new klu_common[ncells];
      data_wk->Symbolic = new klu_symbolic*[ncells];
      data_wk->Numeric = new klu_numeric*[ncells];
      /* Sparse Matrices for It Sparse KLU block-solve */
      data_wk->PS = new SUNMatrix[ncells];
      /* Nb of non zero elements*/
      SPARSITY_INFO_PRECOND(&(data_wk->NNZ),&HP);
      //printf("--> SPARSE Preconditioner -- non zero entries %d represents %f %% fill pattern.\n", data_wk->NNZ, data_wk->NNZ/float((neq+1) * (neq+1)) *100.0);
      data_wk->indx = new int[data_wk->NNZ];
      data_wk->JSPSmat = new realtype[data_wk->NNZ];
      for(int i = 0; i < ncells; ++i) {
          (data_wk->PS)[i] = SUNSparseMatrix(neq+1, neq+1, data_wk->NNZ, CSC_MAT);
          data_wk->colPtrs[i] = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[i]); 
          data_wk->rowVals[i] = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[i]);
          data_wk->Jdata[i] = SUNSparseMatrix_Data((data_wk->PS)[i]);
          SPARSITY_PREPROC_PRECOND(data_wk->rowVals[i],data_wk->colPtrs[i],data_wk->indx,&HP);
          klu_defaults (&(data_wk->Common[i]));
          data_wk->Symbolic[i] = klu_analyze (neq+1, data_wk->colPtrs[i], data_wk->rowVals[i], &(data_wk->Common[i])) ; 
      }
  }
#endif

  return(data_wk);
}



/* Free memory */
void reactor_close(void *cvodeMem, SUNLinearSolver LS, SUNMatrix A, void *data, N_Vector y, N_Vector v_tol){

  UserData data_wk;
  data_wk = (UserData) data;

  CVodeFree(&cvodeMem);
  SUNLinSolFree(LS);
  if (data_wk->iDense_Creact == 1) {
    SUNMatDestroy(A);
  }
  N_VDestroy(y); 
  /* hacks */
  data_wk->reactor_cvode_initialized = false;
  data_wk->FirstTimePrecond = true;
  /* NRG */  
  free(data_wk->rhoX_init);
  free(data_wk->rhoXsrc_ext);
  free(data_wk->rYsrc);

#ifndef USE_KLU 
  if ( data_wk->iDense_Creact == 99) {
      for(int i = 0; i < data_wk->NCELLS; ++i) {
          destroyMat((data_wk->P)[i][i]);
          destroyMat((data_wk->Jbd)[i][i]);
          destroyArray((data_wk->pivot)[i][i]);
      }
  }
#endif

  free(data_wk);

  N_VDestroy(v_tol);          /* Free the atol vector */
}


