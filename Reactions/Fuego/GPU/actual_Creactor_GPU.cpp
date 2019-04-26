#include <actual_Creactor_GPU.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>

/**********************************/
/* Global Variables */
  N_Vector y         = NULL;
  SUNLinearSolver LS = NULL;
  int NEQ       = 0;
  int NCELLS    = 0;
  int iDense_Creact = 1;
  int iJac_Creact   = 0;
  int iE_Creact     = 1;
  int iverbose      = 1;
  void *cvode_mem   = NULL;
  double *rhoe_init = NULL;
  double *rhoh_init = NULL;
  double *rhoesrc_ext = NULL;
  double *rhohsrc_ext = NULL;
  double *rYsrc = NULL;
  double *dt_save;
  double *gamma_d;
  double *temp_old = NULL;
  cusparseMatDescr_t descrA;
  cusolverSpHandle_t cusolverHandle;
  cusparseHandle_t cusparseHandle;
  csrqrInfo_t info;
  void *buffer_qr = NULL;
  size_t workspaceInBytes, internalDataInBytes;
  cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
  cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
  cudaError_t cudaStat1 = cudaSuccess;
  //UserData user_data;

/**********************************/
/* Definitions */
int reactor_init(const int* cvode_iE,const int* Ncells){ 

	int flag, neq_tot ;
	realtype reltol, time;
	N_Vector atol;
	realtype *ratol;
	int mm, ii, nfit;

	CKINDX(&mm,&NEQ,&ii,&nfit);
        if (iverbose > 0) {
	    printf("Nb of spec is %d \n", NEQ);
	}

	/* ParmParse from the inputs file */ 
	amrex::ParmParse pp("ns");
	pp.query("cvode_iJac",iJac_Creact);
	pp.query("cvode_iDense", iDense_Creact);

	/* Args */
	iE_Creact      = *cvode_iE;
	NCELLS         = *Ncells;
        neq_tot        = (NEQ + 1) * NCELLS;

        if (iverbose > 0) {
	    printf("Ncells in one solve ? %d\n",NCELLS);
	}

        /* User data */
        UserData user_data;
        cudaMallocManaged(&user_data, sizeof(struct CVodeUserData));
        user_data->ncells_d[0] = NCELLS;
        user_data->neqs_per_cell[0] = NEQ;
        user_data->flagP = iE_Creact; 

        if ((iDense_Creact == 99) && (iJac_Creact == 1)) { 
            int HP;
            if (iE_Creact == 1) {
                HP = 0;
            } else {
                HP = 1;
            }
            /* Precond data */ 
            if (iverbose > 0) {
                printf("Alloc stuff for Precond \n");
                // Find sparsity pattern to fill structure of sparse matrix
                SPARSITY_INFO_PRECOND(&(user_data->NNZ),&HP);
                printf("--> SPARSE Preconditioner -- non zero entries %d represents %f %% fill pattern.\n", user_data->NNZ, user_data->NNZ/float((NEQ+1) * (NEQ+1)) *100.0);
            }
            cudaMallocManaged(&(user_data->csr_row_count_d), (NEQ+2) * sizeof(int));
            cudaMallocManaged(&(user_data->csr_col_index_d), user_data->NNZ * sizeof(int));
            cudaMallocManaged(&(user_data->csr_jac_d), user_data->NNZ * NCELLS * sizeof(double));
            cudaMallocManaged(&(user_data->csr_val_d), user_data->NNZ * NCELLS * sizeof(double));

            SPARSITY_PREPROC_PRECOND_GPU(user_data->csr_row_count_d, user_data->csr_col_index_d, &HP);
            if (iverbose > 1) {
                for (int i=0; i<NEQ+1; i++) {
                    printf("\n row %d csr_row_count %d \n", i, user_data->csr_row_count_d[i+1]);
                }
            }

            // Create Sparse batch QR solver
            // qr info and matrix descriptor
            cusolver_status = cusolverSpCreate(&cusolverHandle);
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

            cusparse_status = cusparseCreateMatDescr(&descrA); 
            assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

            cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ONE);
 
            cusparse_status = cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ONE);
            assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

            cusolver_status = cusolverSpCreateCsrqrInfo(&info);
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

            // symbolic analysis
            cusolver_status = cusolverSpXcsrqrAnalysisBatched(cusolverHandle,
                                                      NEQ+1, // size per subsystem
                                                      NEQ+1, // size per subsystem
                                                      user_data->NNZ,
                                                      descrA,
                                                      user_data->csr_row_count_d,
                                                      user_data->csr_col_index_d,
                                                      info);
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);


            // allocate working space 
            cusolver_status = cusolverSpDcsrqrBufferInfoBatched(cusolverHandle,
                                                      NEQ+1, // size per subsystem
                                                      NEQ+1, // size per subsystem
                                                      user_data->NNZ,
                                                      descrA,
                                                      user_data->csr_val_d,
                                                      user_data->csr_row_count_d,
                                                      user_data->csr_col_index_d,
                                                      NCELLS,
                                                      info,
                                                      &internalDataInBytes,
                                                      &workspaceInBytes);
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
            
            cudaStat1 = cudaMalloc((void**)&buffer_qr, workspaceInBytes);
            assert(cudaStat1 == cudaSuccess);

        }

	/* Initialize chemistry onto the device */
        initialize_chemistry_device(user_data);

	/* Definition of main vector */
	y = N_VNew_Cuda(neq_tot);
	if(check_flag((void*)y, "N_VNew_Cuda", 0)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	//if ((*cvode_meth == 2) && (*cvode_itmeth == 2))
	//{
	cvode_mem = CVodeCreate(CV_BDF);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
	//} else {
	//    amrex::Abort("\n--> Weird inputs to CVodeCreate. Viable options are CV_BDF (=2), CV_NEWTON (=2)\n");
	//}

        flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));

        time = 0.0e+0;
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function, the inital time, and 
	 * initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, cF_RHS, time, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	
	/* Definition of tolerances: one for each species */
	reltol = 1.0e-10;
        atol  = N_VNew_Cuda(neq_tot);
	ratol = N_VGetHostArrayPointer_Cuda(atol);
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = 1.0e-10;
        }
	N_VCopyToDevice_Cuda(atol);
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, atol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	//if (iDense_Creact == 1) {
            //cuSolver_method LinearSolverMethod = QR;
            //flag = cv_cuSolver_SetLinearSolver(cvode_mem, LinearSolverMethod, false, 0);
            //flag = cv_cuSolver_CSR_SetSizes(cvode_mem, NEQ+1, (NEQ+1)*(NEQ+1),NCELLS);
            //flag = cv_cuSolver_SetJacFun(cvode_mem, &fun_csr_jac);
            //int csr_row_count[NEQ+2];
            //int csr_col_index[(NEQ+1)*(NEQ+1)];
            //csr_row_count[0] = 1;
            //for (int i=0; i<NEQ+1; i++) {
            //    csr_row_count[i+1] = csr_row_count[i] + (NEQ+1);
            //    printf("\n csr_row_count ? %d \n", csr_row_count[i+1]);
            //    for (int j=0; j<NEQ+1; j++) {
            //        csr_col_index[i*(NEQ+1) + j] = j+1 ;
            //        printf("csr_col_index ? %d", csr_col_index[i*(NEQ+1) + j]); 
            //    }
            //}
            //flag = cv_cuSolver_SystemInitialize(cvode_mem, &csr_row_count[0], &csr_col_index[0]);
	//} else 
	if (iDense_Creact == 99) {
            printf("\n--> Using an Iterative Solver \n");

            /* Create the linear solver object */
            if (iJac_Creact == 0) { 
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
	    amrex::Abort("\n--> Only solver implemented is iterative GMRES ...\n");
	}

	if (iJac_Creact == 0) {
            printf("\n--> Without Analytical J\n");
	} else {
            printf("\n--> With Analytical J\n");
	    if (iDense_Creact == 99) {
                if (iverbose > 0) {
                    printf("\n    (99)\n");
		}
	        /* Set the JAcobian-times-vector function */
	        flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
	        if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);

	        /* Set the preconditioner solve and setup functions */
	        //flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
	        flag = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
            }
	}

        /* Set the max number of time steps */ 
	flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

        /* Set the max order */
        flag = CVodeSetMaxOrd(cvode_mem, 5);
        if(check_flag(&flag, "CVodeSetMaxOrd", 1)) return(1);

	/* Define vectors to be used later in creact */
	// GPU stuff: might want to rethink this and put everything in userdata
	cudaMallocManaged(&dt_save, 1*sizeof(double));
	cudaMallocManaged(&gamma_d, 1*sizeof(double));
	if (iE_Creact == 1) { 
	    cudaMallocManaged(&rhoe_init, NCELLS*sizeof(double));
	    cudaMallocManaged(&rhoesrc_ext, NCELLS*sizeof(double));
	} else {
	    cudaMallocManaged(&rhoh_init, NCELLS*sizeof(double));
	    cudaMallocManaged(&rhohsrc_ext, NCELLS*sizeof(double));
	}
	cudaMallocManaged(&rYsrc, (NCELLS*NEQ)*sizeof(double));
	// ReInit stuff
	cudaMallocManaged(&temp_old, NCELLS*sizeof(double));
	for  (int i = 0; i < NCELLS; i++) {
		temp_old[i] = 0.0;
	}

	N_VDestroy(atol);          /* Free the atol vector */

	/* Ok we're done ...*/
        if (iverbose > 0) {
	    printf(" --> DONE WITH INITIALIZATION (GPU) %d \n", iE_Creact);
	}

	return(0);
}

/* Main routine for external looping */
int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
		realtype *P_in, 
                realtype *dt_react, realtype *time, int *Init) {

	realtype time_init, time_out ;
	int flag;

	//cudaError_t cuda_status = cudaSuccess;

        time_init = *time;
	time_out  = *time + *dt_react;

        if (iverbose > 3) {
	    printf("BEG : time curr is %14.6e and dt_react is %14.6e and final time should be %14.6e \n", time_init, *dt_react, time_out);
	}

	/* Get Device MemCpy of in arrays */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
	// rhoY,T
	cudaMemcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS), cudaMemcpyHostToDevice);
	// rhoY_src_ext
	cudaMemcpy(rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double), cudaMemcpyHostToDevice);
	// rhoE/rhoH
	if (iE_Creact == 1) { 
	    cudaMemcpy(rhoe_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	    cudaMemcpy(rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	} else {
	    cudaMemcpy(rhoh_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	    cudaMemcpy(rhohsrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
	}

	/* Call CVODE: ReInit for convergence */
        if (iverbose > 1) {
            printf("\n -------------------------------------\n");
	}
	//if (*Init == 1) {
            //if (iverbose > 1) {
            //    printf("ReInit always \n");
	    //}
	    CVodeReInit(cvode_mem, time_init, y);
	//} else {
	//    // Really bad here when pseveral cells are packed together
	//    double delta_T_max = 0.0;
	//    int offset;
	//    for  (int tid = 0; tid < NCELLS; tid++) {
	//	offset = tid * (NEQ + 1);
	//        temp_old[tid] = fabs(rY_in[offset + NEQ] - temp_old[tid]);
	//	if (delta_T_max < temp_old[tid]) {
	//	    delta_T_max = temp_old[tid];
	//	}
	//    }
        //    if (delta_T_max > 50.0) {
        //        if (iverbose > 1) {
	//            printf("ReInit delta_T = %f \n", delta_T_max);
	//	}
	//        CVodeReInit(cvode_mem, time_init, y);
	//    } else {
        //        if (iverbose > 1) {
	//            printf("ReInit Partial delta_T = %f \n", delta_T_max);
	//	}
	//        CVodeReInitPartial(cvode_mem, time_init, y);
	//    }
        //}
	flag = CVode(cvode_mem, time_out, y, &time_init, CV_NORMAL);
	if (check_flag(&flag, "CVode", 1)) return(1);

	/* Pack data to return in main routine external */
	cudaMemcpy(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype), cudaMemcpyDeviceToHost);
	for  (int i = 0; i < NCELLS; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}

	//if (*Init != 1) {
	//    int offset;
	//    for  (int tid = 0; tid < NCELLS; tid++) {
	//	offset = tid * (NEQ + 1);
	//        temp_old[tid] = rY_in[offset + NEQ];
	//    }
	//}

	/* If in debug mode: print stats */
        if (iverbose > 1) {
           int ierr;
	   double htmp;
           printf("\n......cvode done:\n");
           ierr = CVodeGetLastStep(cvode_mem, &htmp);
           printf(" -- last successful step size = %4.8e \n",htmp);
           ierr = CVodeGetCurrentStep(cvode_mem, &htmp);
           printf(" -- next step to try = %4.8e \n",htmp);
           ierr = CVodeGetCurrentTime(cvode_mem, &htmp);
           printf(" -- integrated time reached = %4.8e \n",htmp);
	   long int itmp, itmp2; 
           ierr = CVodeGetNumSteps(cvode_mem, &itmp);
           printf(" -- number of time steps (nst) = %-6ld \n",itmp);
           ierr = CVodeGetNumRhsEvals(cvode_mem, &itmp);
           ierr = CVSpilsGetNumRhsEvals(cvode_mem, &itmp2);
           itmp = itmp + itmp2;
           printf(" -- number of fRHS EVAL (nfe+nfels) = %-6ld \n", itmp);
           ierr = CVDlsGetNumJacEvals(cvode_mem, &itmp);
           printf(" -- number of Jac EVAL = %-6ld \n", itmp);
	   int itmp3; 
           ierr = CVodeGetLastOrder(cvode_mem, &itmp3);
           printf(" -- method order last used = %d \n", itmp3);
           ierr = CVodeGetCurrentOrder(cvode_mem, &itmp3);
           printf(" -- method order to be used = %d \n", itmp3);
           if (iDense_Creact == 99){
               ierr = CVSpilsGetNumPrecSolves(cvode_mem, &itmp);
	       printf(" -- number of Precond Solves %-6ld \n", itmp);
           }
           ierr = CVodeGetNumNonlinSolvIters(cvode_mem, &itmp);
           printf(" -- number of Newton iterations (nni) %-6ld \n", itmp);
           ierr = CVodeGetNumNonlinSolvConvFails(cvode_mem, &itmp);
	   printf(" -- number of Newton failures %-6ld \n", itmp);
           printf(" -------------------------------------\n");
	}

        long int nfe;
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	return nfe;
}

/* RHS routine used in CVODE */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
		void *user_data){

        //std::chrono::time_point<std::chrono::system_clock> start, end;		
	//start = std::chrono::system_clock::now();
	cudaError_t cuda_status = cudaSuccess;

	/* Get Device pointers for Kernel call */
	realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y_in);
	realtype *ydot_d      = N_VGetDeviceArrayPointer_Cuda(ydot_in);
	// UV/HP
	cudaMemcpy(dt_save, &t, sizeof(double), cudaMemcpyHostToDevice);

        if (iE_Creact == 1) {
	    /* GPU tests */
	    unsigned block = 32;
	    unsigned grid = NCELLS/32 + 1;
	    fKernelSpec<<<grid,block>>>(dt_save, user_data, 
			    yvec_d, ydot_d, 
			    rhoe_init, rhoesrc_ext, rYsrc);
	    cuda_status = cudaDeviceSynchronize();
	    //std::cout << "In fun_rhs, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
	    assert(cuda_status == cudaSuccess);
	} else {
	    unsigned block = 32;
	    unsigned grid = NCELLS/32 + 1;
	    fKernelSpec<<<grid,block>>>(dt_save, user_data,
			    yvec_d, ydot_d, 
			    rhoh_init, rhohsrc_ext, rYsrc);
	    cuda_status = cudaDeviceSynchronize();
	    //std::cout << "In fun_rhs, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
	    assert(cuda_status == cudaSuccess);
	}
	//end = std::chrono::system_clock::now();
	//std::chrono::duration<double> elapsed_seconds = end - start;
	//std::cout << " RHS duration " << elapsed_seconds.count() << std::endl;
	return(0);
}

static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
               booleantype *jcurPtr, realtype gamma, void *user_data) {

        // allocate working space 
        UserData udata = static_cast<CVodeUserData*>(user_data);

        cudaError_t cuda_status = cudaSuccess;

        /* Get Device pointers for Kernel call */
        realtype *u_d      = N_VGetDeviceArrayPointer_Cuda(u);
        realtype *udot_d   = N_VGetDeviceArrayPointer_Cuda(fu);

	cudaMemcpy(dt_save, &tn, sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(gamma_d, &gamma, sizeof(double), cudaMemcpyHostToDevice);

        if (jok) {
            //printf(" jok is OK \n");
	    unsigned block = 32;
	    unsigned grid = NCELLS/32 + 1;
	    fKernelComputeAJ<<<grid,block>>>(user_data, u_d, udot_d,gamma_d, udata->csr_val_d);
            cuda_status = cudaDeviceSynchronize();  
            assert(cuda_status == cudaSuccess);

            *jcurPtr = SUNFALSE;
        } else {
            //printf(" jok is NOT OK \n");
	    unsigned block = 32;
	    unsigned grid = NCELLS/32 + 1;
	    fKernelComputeAJ<<<grid,block>>>(user_data, u_d, udot_d,gamma_d, udata->csr_val_d);
            cuda_status = cudaDeviceSynchronize();  
            assert(cuda_status == cudaSuccess);

            *jcurPtr = SUNTRUE;
        }

        //printf(" ... after kernels, NNZ ? %d \n", udata->NNZ);
        //for (int i = 0; i < udata->NNZ; i++) {
        //     printf(" FIRST CELL csr_val_d[ %d ] = %14.6e \n", i, udata->csr_val_d[i]); //, csr_val[i]);
        //}

        cusolver_status = cusolverSpDcsrqrBufferInfoBatched(cusolverHandle,NEQ+1,NEQ+1, 
                                (udata->NNZ),
                                descrA,
                                udata->csr_val_d,
                                udata->csr_row_count_d,
                                udata->csr_col_index_d,
                                NCELLS,
                                info,
                                &internalDataInBytes,
                                &workspaceInBytes);

        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

	return(0);
}

static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{

        cudaError_t cuda_status = cudaSuccess;

        UserData udata = static_cast<CVodeUserData*>(user_data);

        /* Get Device pointers for Kernel call */
        realtype *u_d      = N_VGetDeviceArrayPointer_Cuda(u);
        realtype *udot_d   = N_VGetDeviceArrayPointer_Cuda(fu);

        realtype *z_d      = N_VGetDeviceArrayPointer_Cuda(z);
        realtype *r_d      = N_VGetDeviceArrayPointer_Cuda(r);

        cusolver_status = cusolverSpDcsrqrsvBatched(cusolverHandle,NEQ+1,NEQ+1,
                               (udata->NNZ),
                               descrA,
                               udata->csr_val_d,
                               udata->csr_row_count_d,
                               udata->csr_col_index_d,
                               r_d, 
                               z_d,
                               NCELLS,
                               info,
                               buffer_qr);


        /* Checks */
        N_VCopyFromDevice_Cuda(z);
        N_VCopyFromDevice_Cuda(r);

        realtype *z_h      = N_VGetHostArrayPointer_Cuda(z);
        realtype *r_h      = N_VGetHostArrayPointer_Cuda(r);

        if (iverbose > 4) {
            for(int batchId = 0 ; batchId < NCELLS; batchId++){
                // measure |bj - Aj*xj|
                double *csrValAj = (udata->csr_val_d) + batchId * (udata->NNZ);
                double *xj       = z_h + batchId * (NEQ+1);
                double *bj       = r_h + batchId * (NEQ+1);
                // sup| bj - Aj*xj|
                double sup_res = 0;
                for(int row = 0 ; row < (NEQ+1) ; row++){
                    const int start = udata->csr_row_count_d[row] - 1;
                    const int end = udata->csr_row_count_d[row +1] - 1;
                    //printf("     row %d =  %d values \n", row, end - start);
                    double Ax = 0.0; // Aj(row,:)*xj
                    for(int colidx = start ; colidx < end ; colidx++){
                        const int col = udata->csr_col_index_d[colidx] - 1;
                        const double Areg = csrValAj[colidx];
                        const double xreg = xj[col];
                        //printf("        Value %d = col is %d and A is %E \n", colidx, col, Areg);
                        Ax = Ax + Areg * xreg;
                    }
                    double r = bj[row] - Ax;
                    sup_res = (sup_res > fabs(r))? sup_res : fabs(r);
                }
                printf("batchId %d: sup|bj - Aj*xj| = %E \n", batchId, sup_res);
            }
        }

        return(0);
}


/*
 * CUDA kernels
 */
//__global__ void fKernelCompute(void *user_data, realtype *u_d, realtype *udot_d, double *gamma)
//{
//  UserData udata = static_cast<CVodeUserData*>(user_data);
//
//  int tid = blockDim.x * blockIdx.x + threadIdx.x;
//
//  if (tid < udata->ncells_d[0]) {
//      /* offsets */
//      int u_offset = tid * (udata->neqs_per_cell[0] + 1); 
//      int jac_offset = tid * udata->NNZ;
//      realtype* u_curr = u_d + u_offset;
//      realtype* csr_jac_cell = udata->csr_jac_d + jac_offset;
//
//      /* Fill the Sps Mat */
//      int nbVals;
//      for (int i = 1; i < udata->neqs_per_cell[0]+2; i++) {
//          nbVals = udata->csr_row_count_d[i]-udata->csr_row_count_d[i-1];
//          printf("Row %d : we have %d nonzero values \n", i-1, nbVals);
//          for (int j = 0; j < nbVals; j++) {
//    	      int idx = udata->csr_col_index_d[ udata->csr_row_count_d[i-1] + j ];
//              /* Scale by -gamma */
//              /* Add identity matrix */
//    	      if (idx == (i-1)) {
//                  udata->csr_jac_cell[ udata->csr_row_count_d[i-1] + j ] = 1.0 - gamma * Jmat[ idx ][ idx ]; 
//    	      } else {
//                  udata->csr_jac_cell[ udata->csr_row_count_d[i-1] + j ] = - gamma * Jmat[ idx ][ i-1 ]; 
//    	      }
//          }
//      }
//
//  }
//
//}

__global__ void fKernelComputeAJ(void *user_data, realtype *u_d, realtype *udot_d, double * gamma, double * csr_val_arg)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  int tid = blockDim.x * blockIdx.x + threadIdx.x;
          
  if (tid < udata->ncells_d[0]) {
      /* local tmp vars */
      realtype activity[21];
      realtype molecular_weight[21];
      realtype temp;
      realtype Jmat[484];

      /* offsets */
      int u_offset = tid * (udata->neqs_per_cell[0] + 1); 
      int jac_offset = tid * udata->NNZ;
      realtype* u_curr = u_d + u_offset;
      //realtype* csr_jac_cell = udata->csr_jac_d + jac_offset;
      //realtype* csr_val_cell = udata->csr_val_d + jac_offset;
      realtype* csr_jac_cell = udata->csr_jac_d + jac_offset;
      realtype* csr_val_cell = csr_val_arg + jac_offset;

      /* MW CGS */
      molecularWeight_d(molecular_weight);
      /* temp */
      temp = u_curr[udata->neqs_per_cell[0]];
      /* Yks, C CGS*/
      for (int i = 0; i < udata->neqs_per_cell[0]; i++){
          activity[i] = u_d[i]/(molecular_weight[i]);
      }
      /* Fuego calls on device 
       * NB to be more accurate should use energy to
       * recompute temp ...      */
      if (udata->flagP == 1){
          int consP = 0 ;
          dwdot_d(Jmat, activity, &temp, &consP, user_data);
      } else {
          int consP = 1 ;
          dwdot_d(Jmat, activity, &temp, &consP, user_data);
      }
      /* renorm the DenseMat */
      for (int i = 0; i < udata->neqs_per_cell[0]; i++){
	  for (int k = 0; k < udata->neqs_per_cell[0]; k++){
              Jmat[k*(udata->neqs_per_cell[0]+1)+i] = Jmat[k*(udata->neqs_per_cell[0]+1)+i] * molecular_weight[i] / molecular_weight[k];
	  }
	  Jmat[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] = Jmat[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] / molecular_weight[i]; 
          Jmat[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] = Jmat[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] * molecular_weight[i]; 
      }
      /* Fill the Sps Mat */
      int nbVals;
      for (int i = 1; i < udata->neqs_per_cell[0]+2; i++) {
          nbVals = udata->csr_row_count_d[i]-udata->csr_row_count_d[i-1];
          for (int j = 0; j < nbVals; j++) {
    	      int idx = udata->csr_col_index_d[ udata->csr_row_count_d[i-1] + j - 1 ] - 1;
              /* Scale by -gamma */
              /* Add identity matrix */
    	      if (idx == (i-1)) {
                  csr_val_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = 1.0 - (*gamma) * Jmat[ idx * (udata->neqs_per_cell[0]+1) + idx ]; 
                  csr_jac_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = Jmat[ idx * (udata->neqs_per_cell[0]+1) + idx ]; 
    	      } else {
                  csr_val_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = - (*gamma) * Jmat[ idx * (udata->neqs_per_cell[0]+1) + i-1 ]; 
                  csr_jac_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = Jmat[ idx * (udata->neqs_per_cell[0]+1) + i-1 ]; 
    	      }
          }
      }

  }

}

__global__ void fKernelSpec(realtype *dt, void *user_data, 
		            realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs)
{
  int tid;
  UserData udata = static_cast<CVodeUserData*>(user_data);

  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < udata->ncells_d[0]) {
      realtype massfrac[21],activity[21];
      realtype Xi[21], cXi[21];
      realtype cdot[21], molecular_weight[21];
      realtype temp, energy;
      int lierr;

      int offset = tid * (udata->neqs_per_cell[0] + 1); 

      /* MW CGS */
      molecularWeight_d(molecular_weight);

      /* rho */ 
      realtype rho = 0.0;
      for (int i = 0; i < udata->neqs_per_cell[0]; i++){
          rho = rho + yvec_d[offset + i];
      }

      /* temp */
      temp = yvec_d[offset + udata->neqs_per_cell[0]];

      /* Yks, C CGS*/
      for (int i = 0; i < udata->neqs_per_cell[0]; i++){
          massfrac[i] = yvec_d[offset + i] / rho;
	  activity[i] = yvec_d[offset + i]/(molecular_weight[i]);
      }

      /* NRG CGS */
      energy = (rhoX_init[tid] + rhoXsrc_ext[tid]*(*dt)) /rho;

      /* Fuego calls on device */
      if (udata->flagP == 1){
          get_t_given_ey_d_(&energy, massfrac, &temp, &lierr);
          ckums_d(&temp, Xi);
          ckcvms_d(&temp, cXi);
      } else {
          get_t_given_hy_d_(&energy, massfrac, &temp, &lierr);
          ckhms_d(&temp, Xi);
          ckcpms_d(&temp, cXi);
      }
      ckwc_d(&temp, activity, cdot, user_data);
      int cX = 0.0;
      for (int i = 0; i < udata->neqs_per_cell[0]; i++){
          cX = cX + massfrac[i] * cXi[i];
      }

      /* Fill ydot vect */
      ydot_d[offset + udata->neqs_per_cell[0]] = rhoXsrc_ext[tid];
      for (int i = 0; i < udata->neqs_per_cell[0]; i++){
          ydot_d[offset + i] = cdot[i] * molecular_weight[i] + rYs[tid * (udata->neqs_per_cell[0]) + i];
          ydot_d[offset + udata->neqs_per_cell[0]] = ydot_d[offset + udata->neqs_per_cell[0]]  - ydot_d[offset + i] * Xi[i];
      }
      ydot_d[offset + udata->neqs_per_cell[0]] = ydot_d[offset + udata->neqs_per_cell[0]] /(rho * cX);
  }
}


__global__ void fKernelJacCSR(realtype t, void *user_data,
                                          realtype *yvec_d, realtype *ydot_d,
                                          realtype* csr_jac,
                                          const int size, const int nnz, 
                                          const int nbatched)
{

    UserData udata = static_cast<CVodeUserData*>(user_data);
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < nbatched) {
        realtype activity[21];
        realtype molecular_weight[21];
        realtype temp;
        realtype Jmat[484];

        int jac_offset = tid * nnz;
        int y_offset = tid * size;

        realtype* csr_jac_cell = csr_jac + jac_offset;
        realtype* actual_y = yvec_d + y_offset;

        /* MW CGS */
        molecularWeight_d(molecular_weight);
        /* rho */ 
        //realtype rho = 0.0;
        //for (int i = 0; i < udata->neqs_per_cell[0]; i++){
        //    rho = rho + actual_y[i];
        //}
        /* temp */
        temp = actual_y[udata->neqs_per_cell[0]];
        /* Yks, C CGS*/
        for (int i = 0; i < udata->neqs_per_cell[0]; i++){
	    activity[i] = actual_y[i]/(molecular_weight[i]);
        }
        /* Fuego calls on device 
         * NB to be more accurate should use energy to
         * recompute temp ...      */
        if (udata->flagP == 1){
            int consP = 0 ;
            dwdot_d(Jmat, activity, &temp, &consP, user_data);
        } else {
            int consP = 1 ;
            dwdot_d(Jmat, activity, &temp, &consP, user_data);
        }
        /* fill the sunMat */
        for (int k = 0; k < udata->neqs_per_cell[0]; k++){
	    for (int i = 0; i < udata->neqs_per_cell[0]; i++){
                csr_jac_cell[k*(udata->neqs_per_cell[0]+1)+i] = Jmat[i*(udata->neqs_per_cell[0]+1)+k] * molecular_weight[k] / molecular_weight[i];
	    }
	    csr_jac_cell[k*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] = Jmat[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+k] * molecular_weight[k]; 
        }
        for (int i = 0; i < udata->neqs_per_cell[0]; i++){
            csr_jac_cell[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] = Jmat[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] / molecular_weight[i]; 
        }
    }
}


 /* Free and destroy memory */
void reactor_close(){

	SUNLinSolFree(LS);
	N_VDestroy(y);          /* Free the y vector */
	CVodeFree(&cvode_mem);
	if (iE_Creact == 1) { 
	    cudaFree(rhoe_init);
	    cudaFree(rhoesrc_ext);
	} else {
	    cudaFree(rhoh_init);
	    cudaFree(rhohsrc_ext);
	}
	cudaFree(rYsrc);
}

/* Get and print some final statistics */
static void PrintFinalStats(void *cvodeMem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvodeMem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_flag(void *flagvalue, const char *funcname, int opt)
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

/* 
 * Non device functions
 */


/*
 * CUDA device functions
 * DRM19
 */
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
__device__ void ckwc_d(double * T, double * C, double * wdot, void *user_data)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate_d(wdot, C, *T, user_data);

    /*convert to chemkin units */
    for (id = 0; id < 21; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}

/*compute the production rate for each species */
__device__ void productionRate_d(double * wdot, double * sc, double T, void *user_data)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    double qdot, q_f[84], q_r[84];
    comp_qfqr_d(q_f, q_r, sc, tc, invT, user_data);

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

__device__ void comp_k_f_d(double * tc, double invT, double * k_f, double * Corr, double * sc, void *user_data)
{

    UserData udata = static_cast<CVodeUserData*>(user_data);

    for (int i=0; i<84; ++i) {
        k_f[i] = udata->prefactor_units[i] * udata->fwd_A[i]
                    * exp(udata->fwd_beta[i] * tc[0] - udata->activation_units[i] * udata->fwd_Ea[i] * invT);
    }
    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 21; ++i) {
        mixture += sc[i];
    }

    /* troe */
    {
        double alpha[8];
        alpha[0] = mixture + (udata->TB[0][0] - 1)*sc[0] + (udata->TB[0][1] - 1)*sc[5] + (udata->TB[0][2] - 1)*sc[10] + (udata->TB[0][3] - 1)*sc[11] + (udata->TB[0][4] - 1)*sc[12] + (udata->TB[0][5] - 1)*sc[18] + (udata->TB[0][6] - 1)*sc[20];
        alpha[1] = mixture + (udata->TB[1][0] - 1)*sc[0] + (udata->TB[1][1] - 1)*sc[5] + (udata->TB[1][2] - 1)*sc[10] + (udata->TB[1][3] - 1)*sc[11] + (udata->TB[1][4] - 1)*sc[12] + (udata->TB[1][5] - 1)*sc[18] + (udata->TB[1][6] - 1)*sc[20];
        alpha[2] = mixture + (udata->TB[2][0] - 1)*sc[0] + (udata->TB[2][1] - 1)*sc[5] + (udata->TB[2][2] - 1)*sc[10] + (udata->TB[2][3] - 1)*sc[11] + (udata->TB[2][4] - 1)*sc[12] + (udata->TB[2][5] - 1)*sc[18] + (udata->TB[2][6] - 1)*sc[20];
        alpha[3] = mixture + (udata->TB[3][0] - 1)*sc[0] + (udata->TB[3][1] - 1)*sc[5] + (udata->TB[3][2] - 1)*sc[10] + (udata->TB[3][3] - 1)*sc[11] + (udata->TB[3][4] - 1)*sc[12] + (udata->TB[3][5] - 1)*sc[18];
        alpha[4] = mixture + (udata->TB[4][0] - 1)*sc[0] + (udata->TB[4][1] - 1)*sc[5] + (udata->TB[4][2] - 1)*sc[10] + (udata->TB[4][3] - 1)*sc[11] + (udata->TB[4][4] - 1)*sc[12] + (udata->TB[4][5] - 1)*sc[18] + (udata->TB[4][6] - 1)*sc[20];
        alpha[5] = mixture + (udata->TB[5][0] - 1)*sc[0] + (udata->TB[5][1] - 1)*sc[5] + (udata->TB[5][2] - 1)*sc[10] + (udata->TB[5][3] - 1)*sc[11] + (udata->TB[5][4] - 1)*sc[12] + (udata->TB[5][5] - 1)*sc[18] + (udata->TB[5][6] - 1)*sc[20];
        alpha[6] = mixture + (udata->TB[6][0] - 1)*sc[0] + (udata->TB[6][1] - 1)*sc[5] + (udata->TB[6][2] - 1)*sc[10] + (udata->TB[6][3] - 1)*sc[11] + (udata->TB[6][4] - 1)*sc[12] + (udata->TB[6][5] - 1)*sc[18] + (udata->TB[6][6] - 1)*sc[20];
        alpha[7] = mixture + (udata->TB[7][0] - 1)*sc[0] + (udata->TB[7][1] - 1)*sc[5] + (udata->TB[7][2] - 1)*sc[10] + (udata->TB[7][3] - 1)*sc[11] + (udata->TB[7][4] - 1)*sc[12] + (udata->TB[7][5] - 1)*sc[18] + (udata->TB[7][6] - 1)*sc[20];
        for (int i=0; i<8; i++)
        {
            double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
            redP = alpha[i-0] / k_f[i] * udata->phase_units[i] * udata->low_A[i] * exp(udata->low_beta[i] * tc[0] - udata->activation_units[i] * udata->low_Ea[i] *invT);
            F = redP / (1.0 + redP);
            logPred = log10(redP);
            logFcent = log10(
                (fabs(udata->troe_Tsss[i]) > 1.e-100 ? (1.-udata->troe_a[i])*exp(-tc[1]/udata->troe_Tsss[i]) : 0.) 
                + (fabs(udata->troe_Ts[i]) > 1.e-100 ? udata->troe_a[i] * exp(-tc[1]/udata->troe_Ts[i]) : 0.) 
                + (udata->troe_len[i] == 4 ? exp(-udata->troe_Tss[i] * invT) : 0.) );
            troe_c = -.4 - .67 * logFcent;
            troe_n = .75 - 1.27 * logFcent;
            troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
            F_troe = pow(10., logFcent / (1.0 + troe*troe));
            Corr[i] = F * F_troe;
        }
    }

    /* simple three-body correction */
    {
        double alpha;
        alpha = mixture + (udata->TB[8][0] - 1)*sc[0] + (udata->TB[8][1] - 1)*sc[5] + (udata->TB[8][2] - 1)*sc[10] + (udata->TB[8][3] - 1)*sc[11] + (udata->TB[8][4] - 1)*sc[12] + (udata->TB[8][5] - 1)*sc[18] + (udata->TB[8][6] - 1)*sc[20];
        Corr[8] = alpha;
        alpha = mixture + (udata->TB[9][0] - 1)*sc[0] + (udata->TB[9][1] - 1)*sc[3] + (udata->TB[9][2] - 1)*sc[5] + (udata->TB[9][3] - 1)*sc[10] + (udata->TB[9][4] - 1)*sc[11] + (udata->TB[9][5] - 1)*sc[12] + (udata->TB[9][6] - 1)*sc[18] + (udata->TB[9][7] - 1)*sc[20];
        Corr[9] = alpha;
        alpha = mixture + (udata->TB[10][0] - 1)*sc[3] + (udata->TB[10][1] - 1)*sc[5] + (udata->TB[10][2] - 1)*sc[11] + (udata->TB[10][3] - 1)*sc[12] + (udata->TB[10][4] - 1)*sc[18] + (udata->TB[10][5] - 1)*sc[19] + (udata->TB[10][6] - 1)*sc[20];
        Corr[10] = alpha;
        alpha = mixture + (udata->TB[11][0] - 1)*sc[0] + (udata->TB[11][1] - 1)*sc[5] + (udata->TB[11][2] - 1)*sc[10] + (udata->TB[11][3] - 1)*sc[12] + (udata->TB[11][4] - 1)*sc[18] + (udata->TB[11][5] - 1)*sc[20];
        Corr[11] = alpha;
        alpha = mixture + (udata->TB[12][0] - 1)*sc[0] + (udata->TB[12][1] - 1)*sc[5] + (udata->TB[12][2] - 1)*sc[10] + (udata->TB[12][3] - 1)*sc[18] + (udata->TB[12][4] - 1)*sc[20];
        Corr[12] = alpha;
        alpha = mixture + (udata->TB[13][0] - 1)*sc[0] + (udata->TB[13][1] - 1)*sc[5] + (udata->TB[13][2] - 1)*sc[10] + (udata->TB[13][3] - 1)*sc[11] + (udata->TB[13][4] - 1)*sc[12] + (udata->TB[13][5] - 1)*sc[18];
        Corr[13] = alpha;
    }
    return;
}

__device__ void comp_Kc_d(double * tc, double invT, double * Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[21];
    gibbs_d(g_RT, tc);

    Kc[0] = g_RT[1] + g_RT[7] - g_RT[9];
    Kc[1] = g_RT[1] + g_RT[9] - g_RT[10];
    Kc[2] = g_RT[1] + g_RT[13] - g_RT[14];
    Kc[3] = g_RT[1] + g_RT[14] - g_RT[15];
    Kc[4] = g_RT[1] + g_RT[16] - g_RT[17];
    Kc[5] = g_RT[1] + g_RT[17] - g_RT[18];
    Kc[6] = g_RT[0] + g_RT[11] - g_RT[14];
    Kc[7] = 2*g_RT[9] - g_RT[18];
    Kc[8] = g_RT[1] + g_RT[2] - g_RT[4];
    Kc[9] = g_RT[2] + g_RT[11] - g_RT[12];
    Kc[10] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[11] = -g_RT[0] + 2*g_RT[1];
    Kc[12] = g_RT[1] + g_RT[4] - g_RT[5];
    Kc[13] = -g_RT[1] - g_RT[11] + g_RT[13];
    Kc[14] = g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4];
    Kc[15] = g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[16] = -g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13];
    Kc[17] = -g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13];
    Kc[18] = -g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14];
    Kc[19] = g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10];
    Kc[20] = g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13];
    Kc[21] = -g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13];
    Kc[22] = g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14];
    Kc[23] = g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16];
    Kc[24] = g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17];
    Kc[25] = g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18];
    Kc[26] = -g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12];
    Kc[27] = g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14];
    Kc[28] = g_RT[1] + 2*g_RT[3] - g_RT[3] - g_RT[6];
    Kc[29] = g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6];
    Kc[30] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19];
    Kc[31] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20];
    Kc[32] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4];
    Kc[33] = g_RT[0] - 2*g_RT[0] + 2*g_RT[1];
    Kc[34] = -g_RT[0] + 2*g_RT[1] + g_RT[5] - g_RT[5];
    Kc[35] = -g_RT[0] + 2*g_RT[1] + g_RT[12] - g_RT[12];
    Kc[36] = -g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6];
    Kc[37] = g_RT[1] - 2*g_RT[4] + g_RT[6];
    Kc[38] = -g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10];
    Kc[39] = -g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13];
    Kc[40] = -g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14];
    Kc[41] = g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15];
    Kc[42] = -g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18];
    Kc[43] = g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5];
    Kc[44] = -g_RT[2] + 2*g_RT[4] - g_RT[5];
    Kc[45] = -g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[46] = -g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14];
    Kc[47] = -g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14];
    Kc[48] = g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9];
    Kc[49] = g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9];
    Kc[50] = g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10];
    Kc[51] = -g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12];
    Kc[52] = g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13];
    Kc[53] = g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14];
    Kc[54] = g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18];
    Kc[55] = -g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14];
    Kc[56] = -g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10];
    Kc[57] = -g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15];
    Kc[58] = -g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12];
    Kc[59] = g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13];
    Kc[60] = g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9];
    Kc[61] = -g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16];
    Kc[62] = g_RT[7] - 2*g_RT[9] + g_RT[10];
    Kc[63] = -g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19];
    Kc[64] = -g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20];
    Kc[65] = -g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11];
    Kc[66] = g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11];
    Kc[67] = g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9];
    Kc[68] = g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8];
    Kc[69] = -g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16];
    Kc[70] = g_RT[8] - 2*g_RT[9] + g_RT[10];
    Kc[71] = -g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11];
    Kc[72] = -g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12];
    Kc[73] = g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14];
    Kc[74] = -g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15];
    Kc[75] = g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14];
    Kc[76] = -g_RT[1] + 2*g_RT[9] - g_RT[17];
    Kc[77] = g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13];
    Kc[78] = g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14];
    Kc[79] = g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18];
    Kc[80] = -g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13];
    Kc[81] = g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13];
    Kc[82] = g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15];
    Kc[83] = g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17];

    for (int i=0; i<84; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refCinv;
    Kc[2] *= refCinv;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refCinv;
    Kc[6] *= refCinv;
    Kc[7] *= refCinv;
    Kc[8] *= refCinv;
    Kc[9] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refCinv;
    Kc[12] *= refCinv;
    Kc[13] *= refC;
    Kc[28] *= refCinv;
    Kc[29] *= refCinv;
    Kc[30] *= refCinv;
    Kc[31] *= refCinv;
    Kc[33] *= refCinv;
    Kc[34] *= refCinv;
    Kc[35] *= refCinv;
    Kc[65] *= refC;
    Kc[80] *= refC;

    return;
}

__device__ void comp_qfqr_d(double *  qf, double * qr, double * sc, double * tc, double invT, void *user_data)
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

    double T = tc[1];
    double Corr[84];
    for (int i = 0; i < 84; ++i) {
        Corr[i] = 1.0;
    }
    double k_f_save[84];
    double Kc_save[84];
    comp_k_f_d(tc,invT,k_f_save,Corr,sc,user_data);
    comp_Kc_d(tc,invT,Kc_save);


    for (int i=0; i<84; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}

/*compute the reaction Jacobian */
__device__ void dwdot_d(double * J, double * sc, double * Tp, int * consP, void *user_data)
{
    double c[21];

    for (int k=0; k<21; k++) {
        c[k] = 1.e6 * sc[k];
    }

    ajacobian_d(J, c, *Tp, *consP, user_data);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<21; k++) {
        J[462+k] *= 1.e-6;
        J[k*22+21] *= 1.e6;
    }

    return;
}

void initialize_chemistry_device(UserData user_data)
{
    // (0):  H + CH2 (+M) <=> CH3 (+M)
    user_data->fwd_A[0]     = 25000000000000000;
    user_data->fwd_beta[0]  = -0.80000000000000004;
    user_data->fwd_Ea[0]    = 0;
    user_data->low_A[0]     = 3.2000000000000002e+27;
    user_data->low_beta[0]  = -3.1400000000000001;
    user_data->low_Ea[0]    = 1230;
    user_data->troe_a[0]    = 0.68000000000000005;
    user_data->troe_Tsss[0] = 78;
    user_data->troe_Ts[0]   = 1995;
    user_data->troe_Tss[0]  = 5590;
    user_data->troe_len[0]  = 4;
    user_data->prefactor_units[0]  = 1.0000000000000002e-06;
    user_data->activation_units[0] = 0.50321666580471969;
    user_data->phase_units[0]      = 1e-12;
    user_data->is_PD[0] = 1;
    user_data->nTB[0] = 7;
    cudaMallocManaged(&user_data->TB[0], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[0], 7 * sizeof(int));
    user_data->TBid[0][0] = 0; user_data->TB[0][0] = 2; // H2
    user_data->TBid[0][1] = 5; user_data->TB[0][1] = 6; // H2O
    user_data->TBid[0][2] = 10; user_data->TB[0][2] = 2; // CH4
    user_data->TBid[0][3] = 11; user_data->TB[0][3] = 1.5; // CO
    user_data->TBid[0][4] = 12; user_data->TB[0][4] = 2; // CO2
    user_data->TBid[0][5] = 18; user_data->TB[0][5] = 3; // C2H6
    user_data->TBid[0][6] = 20; user_data->TB[0][6] = 0.69999999999999996; // AR

    // (1):  H + CH3 (+M) <=> CH4 (+M)
    user_data->fwd_A[1]     = 12700000000000000;
    user_data->fwd_beta[1]  = -0.63;
    user_data->fwd_Ea[1]    = 383;
    user_data->low_A[1]     = 2.4769999999999999e+33;
    user_data->low_beta[1]  = -4.7599999999999998;
    user_data->low_Ea[1]    = 2440;
    user_data->troe_a[1]    = 0.78300000000000003;
    user_data->troe_Tsss[1] = 74;
    user_data->troe_Ts[1]   = 2941;
    user_data->troe_Tss[1]  = 6964;
    user_data->troe_len[1]  = 4;
    user_data->prefactor_units[1]  = 1.0000000000000002e-06;
    user_data->activation_units[1] = 0.50321666580471969;
    user_data->phase_units[1]      = 1e-12;
    user_data->is_PD[1] = 1;
    user_data->nTB[1] = 7;
    cudaMallocManaged(&user_data->TB[1], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[1], 7 * sizeof(int));
    user_data->TBid[1][0] = 0; user_data->TB[1][0] = 2; // H2
    user_data->TBid[1][1] = 5; user_data->TB[1][1] = 6; // H2O
    user_data->TBid[1][2] = 10; user_data->TB[1][2] = 2; // CH4
    user_data->TBid[1][3] = 11; user_data->TB[1][3] = 1.5; // CO
    user_data->TBid[1][4] = 12; user_data->TB[1][4] = 2; // CO2
    user_data->TBid[1][5] = 18; user_data->TB[1][5] = 3; // C2H6
    user_data->TBid[1][6] = 20; user_data->TB[1][6] = 0.69999999999999996; // AR

    // (2):  H + HCO (+M) <=> CH2O (+M)
    user_data->fwd_A[2]     = 1090000000000;
    user_data->fwd_beta[2]  = 0.47999999999999998;
    user_data->fwd_Ea[2]    = -260;
    user_data->low_A[2]     = 1.35e+24;
    user_data->low_beta[2]  = -2.5699999999999998;
    user_data->low_Ea[2]    = 1425;
    user_data->troe_a[2]    = 0.78239999999999998;
    user_data->troe_Tsss[2] = 271;
    user_data->troe_Ts[2]   = 2755;
    user_data->troe_Tss[2]  = 6570;
    user_data->troe_len[2]  = 4;
    user_data->prefactor_units[2]  = 1.0000000000000002e-06;
    user_data->activation_units[2] = 0.50321666580471969;
    user_data->phase_units[2]      = 1e-12;
    user_data->is_PD[2] = 1;
    user_data->nTB[2] = 7;
    cudaMallocManaged(&user_data->TB[2], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[2], 7 * sizeof(int));
    user_data->TBid[2][0] = 0; user_data->TB[2][0] = 2; // H2
    user_data->TBid[2][1] = 5; user_data->TB[2][1] = 6; // H2O
    user_data->TBid[2][2] = 10; user_data->TB[2][2] = 2; // CH4
    user_data->TBid[2][3] = 11; user_data->TB[2][3] = 1.5; // CO
    user_data->TBid[2][4] = 12; user_data->TB[2][4] = 2; // CO2
    user_data->TBid[2][5] = 18; user_data->TB[2][5] = 3; // C2H6
    user_data->TBid[2][6] = 20; user_data->TB[2][6] = 0.69999999999999996; // AR

    // (3):  H + CH2O (+M) <=> CH3O (+M)
    user_data->fwd_A[3]     = 540000000000;
    user_data->fwd_beta[3]  = 0.45400000000000001;
    user_data->fwd_Ea[3]    = 2600;
    user_data->low_A[3]     = 2.2e+30;
    user_data->low_beta[3]  = -4.7999999999999998;
    user_data->low_Ea[3]    = 5560;
    user_data->troe_a[3]    = 0.75800000000000001;
    user_data->troe_Tsss[3] = 94;
    user_data->troe_Ts[3]   = 1555;
    user_data->troe_Tss[3]  = 4200;
    user_data->troe_len[3]  = 4;
    user_data->prefactor_units[3]  = 1.0000000000000002e-06;
    user_data->activation_units[3] = 0.50321666580471969;
    user_data->phase_units[3]      = 1e-12;
    user_data->is_PD[3] = 1;
    user_data->nTB[3] = 6;
    cudaMallocManaged(&user_data->TB[3], 6 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[3], 6 * sizeof(int));
    user_data->TBid[3][0] = 0; user_data->TB[3][0] = 2; // H2
    user_data->TBid[3][1] = 5; user_data->TB[3][1] = 6; // H2O
    user_data->TBid[3][2] = 10; user_data->TB[3][2] = 2; // CH4
    user_data->TBid[3][3] = 11; user_data->TB[3][3] = 1.5; // CO
    user_data->TBid[3][4] = 12; user_data->TB[3][4] = 2; // CO2
    user_data->TBid[3][5] = 18; user_data->TB[3][5] = 3; // C2H6

    // (4):  H + C2H4 (+M) <=> C2H5 (+M)
    user_data->fwd_A[4]     = 1080000000000;
    user_data->fwd_beta[4]  = 0.45400000000000001;
    user_data->fwd_Ea[4]    = 1820;
    user_data->low_A[4]     = 1.1999999999999999e+42;
    user_data->low_beta[4]  = -7.6200000000000001;
    user_data->low_Ea[4]    = 6970;
    user_data->troe_a[4]    = 0.97529999999999994;
    user_data->troe_Tsss[4] = 210;
    user_data->troe_Ts[4]   = 984;
    user_data->troe_Tss[4]  = 4374;
    user_data->troe_len[4]  = 4;
    user_data->prefactor_units[4]  = 1.0000000000000002e-06;
    user_data->activation_units[4] = 0.50321666580471969;
    user_data->phase_units[4]      = 1e-12;
    user_data->is_PD[4] = 1;
    user_data->nTB[4] = 7;
    cudaMallocManaged(&user_data->TB[4], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[4], 7 * sizeof(int));
    user_data->TBid[4][0] = 0; user_data->TB[4][0] = 2; // H2
    user_data->TBid[4][1] = 5; user_data->TB[4][1] = 6; // H2O
    user_data->TBid[4][2] = 10; user_data->TB[4][2] = 2; // CH4
    user_data->TBid[4][3] = 11; user_data->TB[4][3] = 1.5; // CO
    user_data->TBid[4][4] = 12; user_data->TB[4][4] = 2; // CO2
    user_data->TBid[4][5] = 18; user_data->TB[4][5] = 3; // C2H6
    user_data->TBid[4][6] = 20; user_data->TB[4][6] = 0.69999999999999996; // AR

    // (5):  H + C2H5 (+M) <=> C2H6 (+M)
    user_data->fwd_A[5]     = 5.21e+17;
    user_data->fwd_beta[5]  = -0.98999999999999999;
    user_data->fwd_Ea[5]    = 1580;
    user_data->low_A[5]     = 1.9900000000000001e+41;
    user_data->low_beta[5]  = -7.0800000000000001;
    user_data->low_Ea[5]    = 6685;
    user_data->troe_a[5]    = 0.84219999999999995;
    user_data->troe_Tsss[5] = 125;
    user_data->troe_Ts[5]   = 2219;
    user_data->troe_Tss[5]  = 6882;
    user_data->troe_len[5]  = 4;
    user_data->prefactor_units[5]  = 1.0000000000000002e-06;
    user_data->activation_units[5] = 0.50321666580471969;
    user_data->phase_units[5]      = 1e-12;
    user_data->is_PD[5] = 1;
    user_data->nTB[5] = 7;
    cudaMallocManaged(&user_data->TB[5], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[5], 7 * sizeof(int));
    user_data->TBid[5][0] = 0; user_data->TB[5][0] = 2; // H2
    user_data->TBid[5][1] = 5; user_data->TB[5][1] = 6; // H2O
    user_data->TBid[5][2] = 10; user_data->TB[5][2] = 2; // CH4
    user_data->TBid[5][3] = 11; user_data->TB[5][3] = 1.5; // CO
    user_data->TBid[5][4] = 12; user_data->TB[5][4] = 2; // CO2
    user_data->TBid[5][5] = 18; user_data->TB[5][5] = 3; // C2H6
    user_data->TBid[5][6] = 20; user_data->TB[5][6] = 0.69999999999999996; // AR

    // (6):  H2 + CO (+M) <=> CH2O (+M)
    user_data->fwd_A[6]     = 43000000;
    user_data->fwd_beta[6]  = 1.5;
    user_data->fwd_Ea[6]    = 79600;
    user_data->low_A[6]     = 5.0699999999999998e+27;
    user_data->low_beta[6]  = -3.4199999999999999;
    user_data->low_Ea[6]    = 84350;
    user_data->troe_a[6]    = 0.93200000000000005;
    user_data->troe_Tsss[6] = 197;
    user_data->troe_Ts[6]   = 1540;
    user_data->troe_Tss[6]  = 10300;
    user_data->troe_len[6]  = 4;
    user_data->prefactor_units[6]  = 1.0000000000000002e-06;
    user_data->activation_units[6] = 0.50321666580471969;
    user_data->phase_units[6]      = 1e-12;
    user_data->is_PD[6] = 1;
    user_data->nTB[6] = 7;
    cudaMallocManaged(&user_data->TB[6], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[6], 7 * sizeof(int));
    user_data->TBid[6][0] = 0; user_data->TB[6][0] = 2; // H2
    user_data->TBid[6][1] = 5; user_data->TB[6][1] = 6; // H2O
    user_data->TBid[6][2] = 10; user_data->TB[6][2] = 2; // CH4
    user_data->TBid[6][3] = 11; user_data->TB[6][3] = 1.5; // CO
    user_data->TBid[6][4] = 12; user_data->TB[6][4] = 2; // CO2
    user_data->TBid[6][5] = 18; user_data->TB[6][5] = 3; // C2H6
    user_data->TBid[6][6] = 20; user_data->TB[6][6] = 0.69999999999999996; // AR

    // (7):  2 CH3 (+M) <=> C2H6 (+M)
    user_data->fwd_A[7]     = 21200000000000000;
    user_data->fwd_beta[7]  = -0.96999999999999997;
    user_data->fwd_Ea[7]    = 620;
    user_data->low_A[7]     = 1.7700000000000001e+50;
    user_data->low_beta[7]  = -9.6699999999999999;
    user_data->low_Ea[7]    = 6220;
    user_data->troe_a[7]    = 0.53249999999999997;
    user_data->troe_Tsss[7] = 151;
    user_data->troe_Ts[7]   = 1038;
    user_data->troe_Tss[7]  = 4970;
    user_data->troe_len[7]  = 4;
    user_data->prefactor_units[7]  = 1.0000000000000002e-06;
    user_data->activation_units[7] = 0.50321666580471969;
    user_data->phase_units[7]      = 1e-12;
    user_data->is_PD[7] = 1;
    user_data->nTB[7] = 7;
    cudaMallocManaged(&user_data->TB[7], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[7], 7 * sizeof(int));
    user_data->TBid[7][0] = 0; user_data->TB[7][0] = 2; // H2
    user_data->TBid[7][1] = 5; user_data->TB[7][1] = 6; // H2O
    user_data->TBid[7][2] = 10; user_data->TB[7][2] = 2; // CH4
    user_data->TBid[7][3] = 11; user_data->TB[7][3] = 1.5; // CO
    user_data->TBid[7][4] = 12; user_data->TB[7][4] = 2; // CO2
    user_data->TBid[7][5] = 18; user_data->TB[7][5] = 3; // C2H6
    user_data->TBid[7][6] = 20; user_data->TB[7][6] = 0.69999999999999996; // AR

    // (8):  O + H + M <=> OH + M
    user_data->fwd_A[8]     = 5e+17;
    user_data->fwd_beta[8]  = -1;
    user_data->fwd_Ea[8]    = 0;
    user_data->prefactor_units[8]  = 1.0000000000000002e-12;
    user_data->activation_units[8] = 0.50321666580471969;
    user_data->phase_units[8]      = 1e-12;
    user_data->is_PD[8] = 0;
    user_data->nTB[8] = 7;
    cudaMallocManaged(&user_data->TB[8], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[8], 7 * sizeof(int));
    user_data->TBid[8][0] = 0; user_data->TB[8][0] = 2; // H2
    user_data->TBid[8][1] = 5; user_data->TB[8][1] = 6; // H2O
    user_data->TBid[8][2] = 10; user_data->TB[8][2] = 2; // CH4
    user_data->TBid[8][3] = 11; user_data->TB[8][3] = 1.5; // CO
    user_data->TBid[8][4] = 12; user_data->TB[8][4] = 2; // CO2
    user_data->TBid[8][5] = 18; user_data->TB[8][5] = 3; // C2H6
    user_data->TBid[8][6] = 20; user_data->TB[8][6] = 0.69999999999999996; // AR

    // (9):  O + CO + M <=> CO2 + M
    user_data->fwd_A[9]     = 602000000000000;
    user_data->fwd_beta[9]  = 0;
    user_data->fwd_Ea[9]    = 3000;
    user_data->prefactor_units[9]  = 1.0000000000000002e-12;
    user_data->activation_units[9] = 0.50321666580471969;
    user_data->phase_units[9]      = 1e-12;
    user_data->is_PD[9] = 0;
    user_data->nTB[9] = 8;
    cudaMallocManaged(&user_data->TB[9], 8 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[9], 8 * sizeof(int));
    user_data->TBid[9][0] = 0; user_data->TB[9][0] = 2; // H2
    user_data->TBid[9][1] = 3; user_data->TB[9][1] = 6; // O2
    user_data->TBid[9][2] = 5; user_data->TB[9][2] = 6; // H2O
    user_data->TBid[9][3] = 10; user_data->TB[9][3] = 2; // CH4
    user_data->TBid[9][4] = 11; user_data->TB[9][4] = 1.5; // CO
    user_data->TBid[9][5] = 12; user_data->TB[9][5] = 3.5; // CO2
    user_data->TBid[9][6] = 18; user_data->TB[9][6] = 3; // C2H6
    user_data->TBid[9][7] = 20; user_data->TB[9][7] = 0.5; // AR

    // (10):  H + O2 + M <=> HO2 + M
    user_data->fwd_A[10]     = 2.8e+18;
    user_data->fwd_beta[10]  = -0.85999999999999999;
    user_data->fwd_Ea[10]    = 0;
    user_data->prefactor_units[10]  = 1.0000000000000002e-12;
    user_data->activation_units[10] = 0.50321666580471969;
    user_data->phase_units[10]      = 1e-12;
    user_data->is_PD[10] = 0;
    user_data->nTB[10] = 7;
    cudaMallocManaged(&user_data->TB[10], 7 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[10], 7 * sizeof(int));
    user_data->TBid[10][0] = 3; user_data->TB[10][0] = 0; // O2
    user_data->TBid[10][1] = 5; user_data->TB[10][1] = 0; // H2O
    user_data->TBid[10][2] = 11; user_data->TB[10][2] = 0.75; // CO
    user_data->TBid[10][3] = 12; user_data->TB[10][3] = 1.5; // CO2
    user_data->TBid[10][4] = 18; user_data->TB[10][4] = 1.5; // C2H6
    user_data->TBid[10][5] = 19; user_data->TB[10][5] = 0; // N2
    user_data->TBid[10][6] = 20; user_data->TB[10][6] = 0; // AR

    // (11):  2 H + M <=> H2 + M
    user_data->fwd_A[11]     = 1e+18;
    user_data->fwd_beta[11]  = -1;
    user_data->fwd_Ea[11]    = 0;
    user_data->prefactor_units[11]  = 1.0000000000000002e-12;
    user_data->activation_units[11] = 0.50321666580471969;
    user_data->phase_units[11]      = 1e-12;
    user_data->is_PD[11] = 0;
    user_data->nTB[11] = 6;
    cudaMallocManaged(&user_data->TB[11], 6 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[11], 6 * sizeof(int));
    user_data->TBid[11][0] = 0; user_data->TB[11][0] = 0; // H2
    user_data->TBid[11][1] = 5; user_data->TB[11][1] = 0; // H2O
    user_data->TBid[11][2] = 10; user_data->TB[11][2] = 2; // CH4
    user_data->TBid[11][3] = 12; user_data->TB[11][3] = 0; // CO2
    user_data->TBid[11][4] = 18; user_data->TB[11][4] = 3; // C2H6
    user_data->TBid[11][5] = 20; user_data->TB[11][5] = 0.63; // AR

    // (12):  H + OH + M <=> H2O + M
    user_data->fwd_A[12]     = 2.2e+22;
    user_data->fwd_beta[12]  = -2;
    user_data->fwd_Ea[12]    = 0;
    user_data->prefactor_units[12]  = 1.0000000000000002e-12;
    user_data->activation_units[12] = 0.50321666580471969;
    user_data->phase_units[12]      = 1e-12;
    user_data->is_PD[12] = 0;
    user_data->nTB[12] = 5;
    cudaMallocManaged(&user_data->TB[12], 5 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[12], 5 * sizeof(int));
    user_data->TBid[12][0] = 0; user_data->TB[12][0] = 0.72999999999999998; // H2
    user_data->TBid[12][1] = 5; user_data->TB[12][1] = 3.6499999999999999; // H2O
    user_data->TBid[12][2] = 10; user_data->TB[12][2] = 2; // CH4
    user_data->TBid[12][3] = 18; user_data->TB[12][3] = 3; // C2H6
    user_data->TBid[12][4] = 20; user_data->TB[12][4] = 0.38; // AR

    // (13):  HCO + M <=> H + CO + M
    user_data->fwd_A[13]     = 1.87e+17;
    user_data->fwd_beta[13]  = -1;
    user_data->fwd_Ea[13]    = 17000;
    user_data->prefactor_units[13]  = 1.0000000000000002e-06;
    user_data->activation_units[13] = 0.50321666580471969;
    user_data->phase_units[13]      = 1e-6;
    user_data->is_PD[13] = 0;
    user_data->nTB[13] = 6;
    cudaMallocManaged(&user_data->TB[13], 6 * sizeof(double));
    cudaMallocManaged(&user_data->TBid[13], 6 * sizeof(int));
    user_data->TBid[13][0] = 0; user_data->TB[13][0] = 2; // H2
    user_data->TBid[13][1] = 5; user_data->TB[13][1] = 0; // H2O
    user_data->TBid[13][2] = 10; user_data->TB[13][2] = 2; // CH4
    user_data->TBid[13][3] = 11; user_data->TB[13][3] = 1.5; // CO
    user_data->TBid[13][4] = 12; user_data->TB[13][4] = 2; // CO2
    user_data->TBid[13][5] = 18; user_data->TB[13][5] = 3; // C2H6

    // (14):  O + H2 <=> H + OH
    user_data->fwd_A[14]     = 50000;
    user_data->fwd_beta[14]  = 2.6699999999999999;
    user_data->fwd_Ea[14]    = 6290;
    user_data->prefactor_units[14]  = 1.0000000000000002e-06;
    user_data->activation_units[14] = 0.50321666580471969;
    user_data->phase_units[14]      = 1e-12;
    user_data->is_PD[14] = 0;
    user_data->nTB[14] = 0;

    // (15):  O + HO2 <=> OH + O2
    user_data->fwd_A[15]     = 20000000000000;
    user_data->fwd_beta[15]  = 0;
    user_data->fwd_Ea[15]    = 0;
    user_data->prefactor_units[15]  = 1.0000000000000002e-06;
    user_data->activation_units[15] = 0.50321666580471969;
    user_data->phase_units[15]      = 1e-12;
    user_data->is_PD[15] = 0;
    user_data->nTB[15] = 0;

    // (16):  O + CH2 <=> H + HCO
    user_data->fwd_A[16]     = 80000000000000;
    user_data->fwd_beta[16]  = 0;
    user_data->fwd_Ea[16]    = 0;
    user_data->prefactor_units[16]  = 1.0000000000000002e-06;
    user_data->activation_units[16] = 0.50321666580471969;
    user_data->phase_units[16]      = 1e-12;
    user_data->is_PD[16] = 0;
    user_data->nTB[16] = 0;

    // (17):  O + CH2(S) <=> H + HCO
    user_data->fwd_A[17]     = 15000000000000;
    user_data->fwd_beta[17]  = 0;
    user_data->fwd_Ea[17]    = 0;
    user_data->prefactor_units[17]  = 1.0000000000000002e-06;
    user_data->activation_units[17] = 0.50321666580471969;
    user_data->phase_units[17]      = 1e-12;
    user_data->is_PD[17] = 0;
    user_data->nTB[17] = 0;

    // (18):  O + CH3 <=> H + CH2O
    user_data->fwd_A[18]     = 84300000000000;
    user_data->fwd_beta[18]  = 0;
    user_data->fwd_Ea[18]    = 0;
    user_data->prefactor_units[18]  = 1.0000000000000002e-06;
    user_data->activation_units[18] = 0.50321666580471969;
    user_data->phase_units[18]      = 1e-12;
    user_data->is_PD[18] = 0;
    user_data->nTB[18] = 0;

    // (19):  O + CH4 <=> OH + CH3
    user_data->fwd_A[19]     = 1020000000;
    user_data->fwd_beta[19]  = 1.5;
    user_data->fwd_Ea[19]    = 8600;
    user_data->prefactor_units[19]  = 1.0000000000000002e-06;
    user_data->activation_units[19] = 0.50321666580471969;
    user_data->phase_units[19]      = 1e-12;
    user_data->is_PD[19] = 0;
    user_data->nTB[19] = 0;

    // (20):  O + HCO <=> OH + CO
    user_data->fwd_A[20]     = 30000000000000;
    user_data->fwd_beta[20]  = 0;
    user_data->fwd_Ea[20]    = 0;
    user_data->prefactor_units[20]  = 1.0000000000000002e-06;
    user_data->activation_units[20] = 0.50321666580471969;
    user_data->phase_units[20]      = 1e-12;
    user_data->is_PD[20] = 0;
    user_data->nTB[20] = 0;

    // (21):  O + HCO <=> H + CO2
    user_data->fwd_A[21]     = 30000000000000;
    user_data->fwd_beta[21]  = 0;
    user_data->fwd_Ea[21]    = 0;
    user_data->prefactor_units[21]  = 1.0000000000000002e-06;
    user_data->activation_units[21] = 0.50321666580471969;
    user_data->phase_units[21]      = 1e-12;
    user_data->is_PD[21] = 0;
    user_data->nTB[21] = 0;

    // (22):  O + CH2O <=> OH + HCO
    user_data->fwd_A[22]     = 39000000000000;
    user_data->fwd_beta[22]  = 0;
    user_data->fwd_Ea[22]    = 3540;
    user_data->prefactor_units[22]  = 1.0000000000000002e-06;
    user_data->activation_units[22] = 0.50321666580471969;
    user_data->phase_units[22]      = 1e-12;
    user_data->is_PD[22] = 0;
    user_data->nTB[22] = 0;

    // (23):  O + C2H4 <=> CH3 + HCO
    user_data->fwd_A[23]     = 19200000;
    user_data->fwd_beta[23]  = 1.8300000000000001;
    user_data->fwd_Ea[23]    = 220;
    user_data->prefactor_units[23]  = 1.0000000000000002e-06;
    user_data->activation_units[23] = 0.50321666580471969;
    user_data->phase_units[23]      = 1e-12;
    user_data->is_PD[23] = 0;
    user_data->nTB[23] = 0;

    // (24):  O + C2H5 <=> CH3 + CH2O
    user_data->fwd_A[24]     = 132000000000000;
    user_data->fwd_beta[24]  = 0;
    user_data->fwd_Ea[24]    = 0;
    user_data->prefactor_units[24]  = 1.0000000000000002e-06;
    user_data->activation_units[24] = 0.50321666580471969;
    user_data->phase_units[24]      = 1e-12;
    user_data->is_PD[24] = 0;
    user_data->nTB[24] = 0;

    // (25):  O + C2H6 <=> OH + C2H5
    user_data->fwd_A[25]     = 89800000;
    user_data->fwd_beta[25]  = 1.9199999999999999;
    user_data->fwd_Ea[25]    = 5690;
    user_data->prefactor_units[25]  = 1.0000000000000002e-06;
    user_data->activation_units[25] = 0.50321666580471969;
    user_data->phase_units[25]      = 1e-12;
    user_data->is_PD[25] = 0;
    user_data->nTB[25] = 0;

    // (26):  O2 + CO <=> O + CO2
    user_data->fwd_A[26]     = 2500000000000;
    user_data->fwd_beta[26]  = 0;
    user_data->fwd_Ea[26]    = 47800;
    user_data->prefactor_units[26]  = 1.0000000000000002e-06;
    user_data->activation_units[26] = 0.50321666580471969;
    user_data->phase_units[26]      = 1e-12;
    user_data->is_PD[26] = 0;
    user_data->nTB[26] = 0;

    // (27):  O2 + CH2O <=> HO2 + HCO
    user_data->fwd_A[27]     = 100000000000000;
    user_data->fwd_beta[27]  = 0;
    user_data->fwd_Ea[27]    = 40000;
    user_data->prefactor_units[27]  = 1.0000000000000002e-06;
    user_data->activation_units[27] = 0.50321666580471969;
    user_data->phase_units[27]      = 1e-12;
    user_data->is_PD[27] = 0;
    user_data->nTB[27] = 0;

    // (28):  H + 2 O2 <=> HO2 + O2
    user_data->fwd_A[28]     = 3e+20;
    user_data->fwd_beta[28]  = -1.72;
    user_data->fwd_Ea[28]    = 0;
    user_data->prefactor_units[28]  = 1.0000000000000002e-12;
    user_data->activation_units[28] = 0.50321666580471969;
    user_data->phase_units[28]      = 1e-18;
    user_data->is_PD[28] = 0;
    user_data->nTB[28] = 0;

    // (29):  H + O2 + H2O <=> HO2 + H2O
    user_data->fwd_A[29]     = 9.38e+18;
    user_data->fwd_beta[29]  = -0.76000000000000001;
    user_data->fwd_Ea[29]    = 0;
    user_data->prefactor_units[29]  = 1.0000000000000002e-12;
    user_data->activation_units[29] = 0.50321666580471969;
    user_data->phase_units[29]      = 1e-18;
    user_data->is_PD[29] = 0;
    user_data->nTB[29] = 0;

    // (30):  H + O2 + N2 <=> HO2 + N2
    user_data->fwd_A[30]     = 3.75e+20;
    user_data->fwd_beta[30]  = -1.72;
    user_data->fwd_Ea[30]    = 0;
    user_data->prefactor_units[30]  = 1.0000000000000002e-12;
    user_data->activation_units[30] = 0.50321666580471969;
    user_data->phase_units[30]      = 1e-18;
    user_data->is_PD[30] = 0;
    user_data->nTB[30] = 0;

    // (31):  H + O2 + AR <=> HO2 + AR
    user_data->fwd_A[31]     = 7e+17;
    user_data->fwd_beta[31]  = -0.80000000000000004;
    user_data->fwd_Ea[31]    = 0;
    user_data->prefactor_units[31]  = 1.0000000000000002e-12;
    user_data->activation_units[31] = 0.50321666580471969;
    user_data->phase_units[31]      = 1e-18;
    user_data->is_PD[31] = 0;
    user_data->nTB[31] = 0;

    // (32):  H + O2 <=> O + OH
    user_data->fwd_A[32]     = 83000000000000;
    user_data->fwd_beta[32]  = 0;
    user_data->fwd_Ea[32]    = 14413;
    user_data->prefactor_units[32]  = 1.0000000000000002e-06;
    user_data->activation_units[32] = 0.50321666580471969;
    user_data->phase_units[32]      = 1e-12;
    user_data->is_PD[32] = 0;
    user_data->nTB[32] = 0;

    // (33):  2 H + H2 <=> 2 H2
    user_data->fwd_A[33]     = 90000000000000000;
    user_data->fwd_beta[33]  = -0.59999999999999998;
    user_data->fwd_Ea[33]    = 0;
    user_data->prefactor_units[33]  = 1.0000000000000002e-12;
    user_data->activation_units[33] = 0.50321666580471969;
    user_data->phase_units[33]      = 1e-18;
    user_data->is_PD[33] = 0;
    user_data->nTB[33] = 0;

    // (34):  2 H + H2O <=> H2 + H2O
    user_data->fwd_A[34]     = 6e+19;
    user_data->fwd_beta[34]  = -1.25;
    user_data->fwd_Ea[34]    = 0;
    user_data->prefactor_units[34]  = 1.0000000000000002e-12;
    user_data->activation_units[34] = 0.50321666580471969;
    user_data->phase_units[34]      = 1e-18;
    user_data->is_PD[34] = 0;
    user_data->nTB[34] = 0;

    // (35):  2 H + CO2 <=> H2 + CO2
    user_data->fwd_A[35]     = 5.5e+20;
    user_data->fwd_beta[35]  = -2;
    user_data->fwd_Ea[35]    = 0;
    user_data->prefactor_units[35]  = 1.0000000000000002e-12;
    user_data->activation_units[35] = 0.50321666580471969;
    user_data->phase_units[35]      = 1e-18;
    user_data->is_PD[35] = 0;
    user_data->nTB[35] = 0;

    // (36):  H + HO2 <=> O2 + H2
    user_data->fwd_A[36]     = 28000000000000;
    user_data->fwd_beta[36]  = 0;
    user_data->fwd_Ea[36]    = 1068;
    user_data->prefactor_units[36]  = 1.0000000000000002e-06;
    user_data->activation_units[36] = 0.50321666580471969;
    user_data->phase_units[36]      = 1e-12;
    user_data->is_PD[36] = 0;
    user_data->nTB[36] = 0;

    // (37):  H + HO2 <=> 2 OH
    user_data->fwd_A[37]     = 134000000000000;
    user_data->fwd_beta[37]  = 0;
    user_data->fwd_Ea[37]    = 635;
    user_data->prefactor_units[37]  = 1.0000000000000002e-06;
    user_data->activation_units[37] = 0.50321666580471969;
    user_data->phase_units[37]      = 1e-12;
    user_data->is_PD[37] = 0;
    user_data->nTB[37] = 0;

    // (38):  H + CH4 <=> CH3 + H2
    user_data->fwd_A[38]     = 660000000;
    user_data->fwd_beta[38]  = 1.6200000000000001;
    user_data->fwd_Ea[38]    = 10840;
    user_data->prefactor_units[38]  = 1.0000000000000002e-06;
    user_data->activation_units[38] = 0.50321666580471969;
    user_data->phase_units[38]      = 1e-12;
    user_data->is_PD[38] = 0;
    user_data->nTB[38] = 0;

    // (39):  H + HCO <=> H2 + CO
    user_data->fwd_A[39]     = 73400000000000;
    user_data->fwd_beta[39]  = 0;
    user_data->fwd_Ea[39]    = 0;
    user_data->prefactor_units[39]  = 1.0000000000000002e-06;
    user_data->activation_units[39] = 0.50321666580471969;
    user_data->phase_units[39]      = 1e-12;
    user_data->is_PD[39] = 0;
    user_data->nTB[39] = 0;

    // (40):  H + CH2O <=> HCO + H2
    user_data->fwd_A[40]     = 23000000000;
    user_data->fwd_beta[40]  = 1.05;
    user_data->fwd_Ea[40]    = 3275;
    user_data->prefactor_units[40]  = 1.0000000000000002e-06;
    user_data->activation_units[40] = 0.50321666580471969;
    user_data->phase_units[40]      = 1e-12;
    user_data->is_PD[40] = 0;
    user_data->nTB[40] = 0;

    // (41):  H + CH3O <=> OH + CH3
    user_data->fwd_A[41]     = 32000000000000;
    user_data->fwd_beta[41]  = 0;
    user_data->fwd_Ea[41]    = 0;
    user_data->prefactor_units[41]  = 1.0000000000000002e-06;
    user_data->activation_units[41] = 0.50321666580471969;
    user_data->phase_units[41]      = 1e-12;
    user_data->is_PD[41] = 0;
    user_data->nTB[41] = 0;

    // (42):  H + C2H6 <=> C2H5 + H2
    user_data->fwd_A[42]     = 115000000;
    user_data->fwd_beta[42]  = 1.8999999999999999;
    user_data->fwd_Ea[42]    = 7530;
    user_data->prefactor_units[42]  = 1.0000000000000002e-06;
    user_data->activation_units[42] = 0.50321666580471969;
    user_data->phase_units[42]      = 1e-12;
    user_data->is_PD[42] = 0;
    user_data->nTB[42] = 0;

    // (43):  OH + H2 <=> H + H2O
    user_data->fwd_A[43]     = 216000000;
    user_data->fwd_beta[43]  = 1.51;
    user_data->fwd_Ea[43]    = 3430;
    user_data->prefactor_units[43]  = 1.0000000000000002e-06;
    user_data->activation_units[43] = 0.50321666580471969;
    user_data->phase_units[43]      = 1e-12;
    user_data->is_PD[43] = 0;
    user_data->nTB[43] = 0;

    // (44):  2 OH <=> O + H2O
    user_data->fwd_A[44]     = 35700;
    user_data->fwd_beta[44]  = 2.3999999999999999;
    user_data->fwd_Ea[44]    = -2110;
    user_data->prefactor_units[44]  = 1.0000000000000002e-06;
    user_data->activation_units[44] = 0.50321666580471969;
    user_data->phase_units[44]      = 1e-12;
    user_data->is_PD[44] = 0;
    user_data->nTB[44] = 0;

    // (45):  OH + HO2 <=> O2 + H2O
    user_data->fwd_A[45]     = 29000000000000;
    user_data->fwd_beta[45]  = 0;
    user_data->fwd_Ea[45]    = -500;
    user_data->prefactor_units[45]  = 1.0000000000000002e-06;
    user_data->activation_units[45] = 0.50321666580471969;
    user_data->phase_units[45]      = 1e-12;
    user_data->is_PD[45] = 0;
    user_data->nTB[45] = 0;

    // (46):  OH + CH2 <=> H + CH2O
    user_data->fwd_A[46]     = 20000000000000;
    user_data->fwd_beta[46]  = 0;
    user_data->fwd_Ea[46]    = 0;
    user_data->prefactor_units[46]  = 1.0000000000000002e-06;
    user_data->activation_units[46] = 0.50321666580471969;
    user_data->phase_units[46]      = 1e-12;
    user_data->is_PD[46] = 0;
    user_data->nTB[46] = 0;

    // (47):  OH + CH2(S) <=> H + CH2O
    user_data->fwd_A[47]     = 30000000000000;
    user_data->fwd_beta[47]  = 0;
    user_data->fwd_Ea[47]    = 0;
    user_data->prefactor_units[47]  = 1.0000000000000002e-06;
    user_data->activation_units[47] = 0.50321666580471969;
    user_data->phase_units[47]      = 1e-12;
    user_data->is_PD[47] = 0;
    user_data->nTB[47] = 0;

    // (48):  OH + CH3 <=> CH2 + H2O
    user_data->fwd_A[48]     = 56000000;
    user_data->fwd_beta[48]  = 1.6000000000000001;
    user_data->fwd_Ea[48]    = 5420;
    user_data->prefactor_units[48]  = 1.0000000000000002e-06;
    user_data->activation_units[48] = 0.50321666580471969;
    user_data->phase_units[48]      = 1e-12;
    user_data->is_PD[48] = 0;
    user_data->nTB[48] = 0;

    // (49):  OH + CH3 <=> CH2(S) + H2O
    user_data->fwd_A[49]     = 25010000000000;
    user_data->fwd_beta[49]  = 0;
    user_data->fwd_Ea[49]    = 0;
    user_data->prefactor_units[49]  = 1.0000000000000002e-06;
    user_data->activation_units[49] = 0.50321666580471969;
    user_data->phase_units[49]      = 1e-12;
    user_data->is_PD[49] = 0;
    user_data->nTB[49] = 0;

    // (50):  OH + CH4 <=> CH3 + H2O
    user_data->fwd_A[50]     = 100000000;
    user_data->fwd_beta[50]  = 1.6000000000000001;
    user_data->fwd_Ea[50]    = 3120;
    user_data->prefactor_units[50]  = 1.0000000000000002e-06;
    user_data->activation_units[50] = 0.50321666580471969;
    user_data->phase_units[50]      = 1e-12;
    user_data->is_PD[50] = 0;
    user_data->nTB[50] = 0;

    // (51):  OH + CO <=> H + CO2
    user_data->fwd_A[51]     = 47600000;
    user_data->fwd_beta[51]  = 1.228;
    user_data->fwd_Ea[51]    = 70;
    user_data->prefactor_units[51]  = 1.0000000000000002e-06;
    user_data->activation_units[51] = 0.50321666580471969;
    user_data->phase_units[51]      = 1e-12;
    user_data->is_PD[51] = 0;
    user_data->nTB[51] = 0;

    // (52):  OH + HCO <=> H2O + CO
    user_data->fwd_A[52]     = 50000000000000;
    user_data->fwd_beta[52]  = 0;
    user_data->fwd_Ea[52]    = 0;
    user_data->prefactor_units[52]  = 1.0000000000000002e-06;
    user_data->activation_units[52] = 0.50321666580471969;
    user_data->phase_units[52]      = 1e-12;
    user_data->is_PD[52] = 0;
    user_data->nTB[52] = 0;

    // (53):  OH + CH2O <=> HCO + H2O
    user_data->fwd_A[53]     = 3430000000;
    user_data->fwd_beta[53]  = 1.1799999999999999;
    user_data->fwd_Ea[53]    = -447;
    user_data->prefactor_units[53]  = 1.0000000000000002e-06;
    user_data->activation_units[53] = 0.50321666580471969;
    user_data->phase_units[53]      = 1e-12;
    user_data->is_PD[53] = 0;
    user_data->nTB[53] = 0;

    // (54):  OH + C2H6 <=> C2H5 + H2O
    user_data->fwd_A[54]     = 3540000;
    user_data->fwd_beta[54]  = 2.1200000000000001;
    user_data->fwd_Ea[54]    = 870;
    user_data->prefactor_units[54]  = 1.0000000000000002e-06;
    user_data->activation_units[54] = 0.50321666580471969;
    user_data->phase_units[54]      = 1e-12;
    user_data->is_PD[54] = 0;
    user_data->nTB[54] = 0;

    // (55):  HO2 + CH2 <=> OH + CH2O
    user_data->fwd_A[55]     = 20000000000000;
    user_data->fwd_beta[55]  = 0;
    user_data->fwd_Ea[55]    = 0;
    user_data->prefactor_units[55]  = 1.0000000000000002e-06;
    user_data->activation_units[55] = 0.50321666580471969;
    user_data->phase_units[55]      = 1e-12;
    user_data->is_PD[55] = 0;
    user_data->nTB[55] = 0;

    // (56):  HO2 + CH3 <=> O2 + CH4
    user_data->fwd_A[56]     = 1000000000000;
    user_data->fwd_beta[56]  = 0;
    user_data->fwd_Ea[56]    = 0;
    user_data->prefactor_units[56]  = 1.0000000000000002e-06;
    user_data->activation_units[56] = 0.50321666580471969;
    user_data->phase_units[56]      = 1e-12;
    user_data->is_PD[56] = 0;
    user_data->nTB[56] = 0;

    // (57):  HO2 + CH3 <=> OH + CH3O
    user_data->fwd_A[57]     = 20000000000000;
    user_data->fwd_beta[57]  = 0;
    user_data->fwd_Ea[57]    = 0;
    user_data->prefactor_units[57]  = 1.0000000000000002e-06;
    user_data->activation_units[57] = 0.50321666580471969;
    user_data->phase_units[57]      = 1e-12;
    user_data->is_PD[57] = 0;
    user_data->nTB[57] = 0;

    // (58):  HO2 + CO <=> OH + CO2
    user_data->fwd_A[58]     = 150000000000000;
    user_data->fwd_beta[58]  = 0;
    user_data->fwd_Ea[58]    = 23600;
    user_data->prefactor_units[58]  = 1.0000000000000002e-06;
    user_data->activation_units[58] = 0.50321666580471969;
    user_data->phase_units[58]      = 1e-12;
    user_data->is_PD[58] = 0;
    user_data->nTB[58] = 0;

    // (59):  CH2 + O2 <=> OH + HCO
    user_data->fwd_A[59]     = 13200000000000;
    user_data->fwd_beta[59]  = 0;
    user_data->fwd_Ea[59]    = 1500;
    user_data->prefactor_units[59]  = 1.0000000000000002e-06;
    user_data->activation_units[59] = 0.50321666580471969;
    user_data->phase_units[59]      = 1e-12;
    user_data->is_PD[59] = 0;
    user_data->nTB[59] = 0;

    // (60):  CH2 + H2 <=> H + CH3
    user_data->fwd_A[60]     = 500000;
    user_data->fwd_beta[60]  = 2;
    user_data->fwd_Ea[60]    = 7230;
    user_data->prefactor_units[60]  = 1.0000000000000002e-06;
    user_data->activation_units[60] = 0.50321666580471969;
    user_data->phase_units[60]      = 1e-12;
    user_data->is_PD[60] = 0;
    user_data->nTB[60] = 0;

    // (61):  CH2 + CH3 <=> H + C2H4
    user_data->fwd_A[61]     = 40000000000000;
    user_data->fwd_beta[61]  = 0;
    user_data->fwd_Ea[61]    = 0;
    user_data->prefactor_units[61]  = 1.0000000000000002e-06;
    user_data->activation_units[61] = 0.50321666580471969;
    user_data->phase_units[61]      = 1e-12;
    user_data->is_PD[61] = 0;
    user_data->nTB[61] = 0;

    // (62):  CH2 + CH4 <=> 2 CH3
    user_data->fwd_A[62]     = 2460000;
    user_data->fwd_beta[62]  = 2;
    user_data->fwd_Ea[62]    = 8270;
    user_data->prefactor_units[62]  = 1.0000000000000002e-06;
    user_data->activation_units[62] = 0.50321666580471969;
    user_data->phase_units[62]      = 1e-12;
    user_data->is_PD[62] = 0;
    user_data->nTB[62] = 0;

    // (63):  CH2(S) + N2 <=> CH2 + N2
    user_data->fwd_A[63]     = 15000000000000;
    user_data->fwd_beta[63]  = 0;
    user_data->fwd_Ea[63]    = 600;
    user_data->prefactor_units[63]  = 1.0000000000000002e-06;
    user_data->activation_units[63] = 0.50321666580471969;
    user_data->phase_units[63]      = 1e-12;
    user_data->is_PD[63] = 0;
    user_data->nTB[63] = 0;

    // (64):  CH2(S) + AR <=> CH2 + AR
    user_data->fwd_A[64]     = 9000000000000;
    user_data->fwd_beta[64]  = 0;
    user_data->fwd_Ea[64]    = 600;
    user_data->prefactor_units[64]  = 1.0000000000000002e-06;
    user_data->activation_units[64] = 0.50321666580471969;
    user_data->phase_units[64]      = 1e-12;
    user_data->is_PD[64] = 0;
    user_data->nTB[64] = 0;

    // (65):  CH2(S) + O2 <=> H + OH + CO
    user_data->fwd_A[65]     = 28000000000000;
    user_data->fwd_beta[65]  = 0;
    user_data->fwd_Ea[65]    = 0;
    user_data->prefactor_units[65]  = 1.0000000000000002e-06;
    user_data->activation_units[65] = 0.50321666580471969;
    user_data->phase_units[65]      = 1e-12;
    user_data->is_PD[65] = 0;
    user_data->nTB[65] = 0;

    // (66):  CH2(S) + O2 <=> CO + H2O
    user_data->fwd_A[66]     = 12000000000000;
    user_data->fwd_beta[66]  = 0;
    user_data->fwd_Ea[66]    = 0;
    user_data->prefactor_units[66]  = 1.0000000000000002e-06;
    user_data->activation_units[66] = 0.50321666580471969;
    user_data->phase_units[66]      = 1e-12;
    user_data->is_PD[66] = 0;
    user_data->nTB[66] = 0;

    // (67):  CH2(S) + H2 <=> CH3 + H
    user_data->fwd_A[67]     = 70000000000000;
    user_data->fwd_beta[67]  = 0;
    user_data->fwd_Ea[67]    = 0;
    user_data->prefactor_units[67]  = 1.0000000000000002e-06;
    user_data->activation_units[67] = 0.50321666580471969;
    user_data->phase_units[67]      = 1e-12;
    user_data->is_PD[67] = 0;
    user_data->nTB[67] = 0;

    // (68):  CH2(S) + H2O <=> CH2 + H2O
    user_data->fwd_A[68]     = 30000000000000;
    user_data->fwd_beta[68]  = 0;
    user_data->fwd_Ea[68]    = 0;
    user_data->prefactor_units[68]  = 1.0000000000000002e-06;
    user_data->activation_units[68] = 0.50321666580471969;
    user_data->phase_units[68]      = 1e-12;
    user_data->is_PD[68] = 0;
    user_data->nTB[68] = 0;

    // (69):  CH2(S) + CH3 <=> H + C2H4
    user_data->fwd_A[69]     = 12000000000000;
    user_data->fwd_beta[69]  = 0;
    user_data->fwd_Ea[69]    = -570;
    user_data->prefactor_units[69]  = 1.0000000000000002e-06;
    user_data->activation_units[69] = 0.50321666580471969;
    user_data->phase_units[69]      = 1e-12;
    user_data->is_PD[69] = 0;
    user_data->nTB[69] = 0;

    // (70):  CH2(S) + CH4 <=> 2 CH3
    user_data->fwd_A[70]     = 16000000000000;
    user_data->fwd_beta[70]  = 0;
    user_data->fwd_Ea[70]    = -570;
    user_data->prefactor_units[70]  = 1.0000000000000002e-06;
    user_data->activation_units[70] = 0.50321666580471969;
    user_data->phase_units[70]      = 1e-12;
    user_data->is_PD[70] = 0;
    user_data->nTB[70] = 0;

    // (71):  CH2(S) + CO <=> CH2 + CO
    user_data->fwd_A[71]     = 9000000000000;
    user_data->fwd_beta[71]  = 0;
    user_data->fwd_Ea[71]    = 0;
    user_data->prefactor_units[71]  = 1.0000000000000002e-06;
    user_data->activation_units[71] = 0.50321666580471969;
    user_data->phase_units[71]      = 1e-12;
    user_data->is_PD[71] = 0;
    user_data->nTB[71] = 0;

    // (72):  CH2(S) + CO2 <=> CH2 + CO2
    user_data->fwd_A[72]     = 7000000000000;
    user_data->fwd_beta[72]  = 0;
    user_data->fwd_Ea[72]    = 0;
    user_data->prefactor_units[72]  = 1.0000000000000002e-06;
    user_data->activation_units[72] = 0.50321666580471969;
    user_data->phase_units[72]      = 1e-12;
    user_data->is_PD[72] = 0;
    user_data->nTB[72] = 0;

    // (73):  CH2(S) + CO2 <=> CO + CH2O
    user_data->fwd_A[73]     = 14000000000000;
    user_data->fwd_beta[73]  = 0;
    user_data->fwd_Ea[73]    = 0;
    user_data->prefactor_units[73]  = 1.0000000000000002e-06;
    user_data->activation_units[73] = 0.50321666580471969;
    user_data->phase_units[73]      = 1e-12;
    user_data->is_PD[73] = 0;
    user_data->nTB[73] = 0;

    // (74):  CH3 + O2 <=> O + CH3O
    user_data->fwd_A[74]     = 26750000000000;
    user_data->fwd_beta[74]  = 0;
    user_data->fwd_Ea[74]    = 28800;
    user_data->prefactor_units[74]  = 1.0000000000000002e-06;
    user_data->activation_units[74] = 0.50321666580471969;
    user_data->phase_units[74]      = 1e-12;
    user_data->is_PD[74] = 0;
    user_data->nTB[74] = 0;

    // (75):  CH3 + O2 <=> OH + CH2O
    user_data->fwd_A[75]     = 36000000000;
    user_data->fwd_beta[75]  = 0;
    user_data->fwd_Ea[75]    = 8940;
    user_data->prefactor_units[75]  = 1.0000000000000002e-06;
    user_data->activation_units[75] = 0.50321666580471969;
    user_data->phase_units[75]      = 1e-12;
    user_data->is_PD[75] = 0;
    user_data->nTB[75] = 0;

    // (76):  2 CH3 <=> H + C2H5
    user_data->fwd_A[76]     = 4990000000000;
    user_data->fwd_beta[76]  = 0.10000000000000001;
    user_data->fwd_Ea[76]    = 10600;
    user_data->prefactor_units[76]  = 1.0000000000000002e-06;
    user_data->activation_units[76] = 0.50321666580471969;
    user_data->phase_units[76]      = 1e-12;
    user_data->is_PD[76] = 0;
    user_data->nTB[76] = 0;

    // (77):  CH3 + HCO <=> CH4 + CO
    user_data->fwd_A[77]     = 26480000000000;
    user_data->fwd_beta[77]  = 0;
    user_data->fwd_Ea[77]    = 0;
    user_data->prefactor_units[77]  = 1.0000000000000002e-06;
    user_data->activation_units[77] = 0.50321666580471969;
    user_data->phase_units[77]      = 1e-12;
    user_data->is_PD[77] = 0;
    user_data->nTB[77] = 0;

    // (78):  CH3 + CH2O <=> HCO + CH4
    user_data->fwd_A[78]     = 3320;
    user_data->fwd_beta[78]  = 2.8100000000000001;
    user_data->fwd_Ea[78]    = 5860;
    user_data->prefactor_units[78]  = 1.0000000000000002e-06;
    user_data->activation_units[78] = 0.50321666580471969;
    user_data->phase_units[78]      = 1e-12;
    user_data->is_PD[78] = 0;
    user_data->nTB[78] = 0;

    // (79):  CH3 + C2H6 <=> C2H5 + CH4
    user_data->fwd_A[79]     = 6140000;
    user_data->fwd_beta[79]  = 1.74;
    user_data->fwd_Ea[79]    = 10450;
    user_data->prefactor_units[79]  = 1.0000000000000002e-06;
    user_data->activation_units[79] = 0.50321666580471969;
    user_data->phase_units[79]      = 1e-12;
    user_data->is_PD[79] = 0;
    user_data->nTB[79] = 0;

    // (80):  HCO + H2O <=> H + CO + H2O
    user_data->fwd_A[80]     = 2.244e+18;
    user_data->fwd_beta[80]  = -1;
    user_data->fwd_Ea[80]    = 17000;
    user_data->prefactor_units[80]  = 1.0000000000000002e-06;
    user_data->activation_units[80] = 0.50321666580471969;
    user_data->phase_units[80]      = 1e-12;
    user_data->is_PD[80] = 0;
    user_data->nTB[80] = 0;

    // (81):  HCO + O2 <=> HO2 + CO
    user_data->fwd_A[81]     = 7600000000000;
    user_data->fwd_beta[81]  = 0;
    user_data->fwd_Ea[81]    = 400;
    user_data->prefactor_units[81]  = 1.0000000000000002e-06;
    user_data->activation_units[81] = 0.50321666580471969;
    user_data->phase_units[81]      = 1e-12;
    user_data->is_PD[81] = 0;
    user_data->nTB[81] = 0;

    // (82):  CH3O + O2 <=> HO2 + CH2O
    user_data->fwd_A[82]     = 4.2799999999999999e-13;
    user_data->fwd_beta[82]  = 7.5999999999999996;
    user_data->fwd_Ea[82]    = -3530;
    user_data->prefactor_units[82]  = 1.0000000000000002e-06;
    user_data->activation_units[82] = 0.50321666580471969;
    user_data->phase_units[82]      = 1e-12;
    user_data->is_PD[82] = 0;
    user_data->nTB[82] = 0;

    // (83):  C2H5 + O2 <=> HO2 + C2H4
    user_data->fwd_A[83]     = 840000000000;
    user_data->fwd_beta[83]  = 0;
    user_data->fwd_Ea[83]    = 3875;
    user_data->prefactor_units[83]  = 1.0000000000000002e-06;
    user_data->activation_units[83] = 0.50321666580471969;
    user_data->phase_units[83]      = 1e-12;
    user_data->is_PD[83] = 0;
    user_data->nTB[83] = 0;

    return;
}
//void finalize_chemistry_device(void *user_data)
//{

    //UserData udata = static_cast<CVodeUserData*>(user_data);

    //for (int i=0; i<84; i++) {
    //    if (nTB[i] != 0) {
    //        nTB[i] = 0;
    //        free(TB[i]);
    //        free(TBid[i]);
    //    }
    //}
    //return;
//}


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
__device__ void ajacobian_d(double * J, double * sc, double T, int consP, void *user_data)
{

    UserData udata = static_cast<CVodeUserData*>(user_data);

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
    alpha = mixture + (udata->TB[0][0] - 1)*sc[0] + (udata->TB[0][1] - 1)*sc[5] + (udata->TB[0][2] - 1)*sc[10] + (udata->TB[0][3] - 1)*sc[11] + (udata->TB[0][4] - 1)*sc[12] + (udata->TB[0][5] - 1)*sc[18] + (udata->TB[0][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = udata->prefactor_units[0] * udata->fwd_A[0]
                * exp(udata->fwd_beta[0] * tc[0] - udata->activation_units[0] * udata->fwd_Ea[0] * invT);
    dlnkfdT = udata->fwd_beta[0] * invT + udata->activation_units[0] * udata->fwd_Ea[0] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[0] * exp(udata->low_beta[0] * tc[0] - udata->activation_units[0] * udata->low_Ea[0] * invT);
    Pr = udata->phase_units[0] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[0] * invT + udata->activation_units[0] * udata->low_Ea[0] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[0]) > 1.e-100 ? (1.-udata->troe_a[0])*exp(-T/udata->troe_Tsss[0]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[0]) > 1.e-100 ? udata->troe_a[0] * exp(-T/udata->troe_Ts[0]) : 0.);
    Fcent3 = (udata->troe_len[0] == 4 ? exp(-udata->troe_Tss[0] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[0]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[0] : 0.)
      + (fabs(udata->troe_Ts[0]) > 1.e-100 ? -Fcent2/udata->troe_Ts[0] : 0.)
      + (udata->troe_len[0] == 4 ? Fcent3*udata->troe_Tss[0]*invT2 : 0.) );
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
        dqdci = (udata->TB[0][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
        J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[7];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
        J[31] += dqdci;               /* dwdot[CH3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[0][1] - 1)*dcdc_fac;
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
        dqdci = (udata->TB[0][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[0][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[249] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[251] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[0][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[271] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[273] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[0][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[403] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[405] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[0][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[447] -= dqdci;              /* dwdot[CH2]/d[AR] */
        J[449] += dqdci;              /* dwdot[CH3]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[0][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[7];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[0][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f*sc[1];
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac - k_r;
        dqdc[10] = udata->TB[0][2]*dcdc_fac;
        dqdc[11] = udata->TB[0][3]*dcdc_fac;
        dqdc[12] = udata->TB[0][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = udata->TB[0][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = udata->TB[0][6]*dcdc_fac;
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
    alpha = mixture + (udata->TB[1][0] - 1)*sc[0] + (udata->TB[1][1] - 1)*sc[5] + (udata->TB[1][2] - 1)*sc[10] + (udata->TB[1][3] - 1)*sc[11] + (udata->TB[1][4] - 1)*sc[12] + (udata->TB[1][5] - 1)*sc[18] + (udata->TB[1][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = udata->prefactor_units[1] * udata->fwd_A[1]
                * exp(udata->fwd_beta[1] * tc[0] - udata->activation_units[1] * udata->fwd_Ea[1] * invT);
    dlnkfdT = udata->fwd_beta[1] * invT + udata->activation_units[1] * udata->fwd_Ea[1] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[1] * exp(udata->low_beta[1] * tc[0] - udata->activation_units[1] * udata->low_Ea[1] * invT);
    Pr = udata->phase_units[1] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[1] * invT + udata->activation_units[1] * udata->low_Ea[1] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[1]) > 1.e-100 ? (1.-udata->troe_a[1])*exp(-T/udata->troe_Tsss[1]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[1]) > 1.e-100 ? udata->troe_a[1] * exp(-T/udata->troe_Ts[1]) : 0.);
    Fcent3 = (udata->troe_len[1] == 4 ? exp(-udata->troe_Tss[1] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[1]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[1] : 0.)
      + (fabs(udata->troe_Ts[1]) > 1.e-100 ? -Fcent2/udata->troe_Ts[1] : 0.)
      + (udata->troe_len[1] == 4 ? Fcent3*udata->troe_Tss[1]*invT2 : 0.) );
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
        dqdci = (udata->TB[1][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[9] -= dqdci;                /* dwdot[CH3]/d[H2] */
        J[10] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[9];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[32] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[1][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[120] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[1];
        J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[1][2] - 1)*dcdc_fac - k_r;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[1][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[1][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[273] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[274] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[1][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[1][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[449] -= dqdci;              /* dwdot[CH3]/d[AR] */
        J[450] += dqdci;              /* dwdot[CH4]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[9];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[1][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*sc[1];
        dqdc[10] = udata->TB[1][2]*dcdc_fac - k_r;
        dqdc[11] = udata->TB[1][3]*dcdc_fac;
        dqdc[12] = udata->TB[1][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = udata->TB[1][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = udata->TB[1][6]*dcdc_fac;
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
    alpha = mixture + (udata->TB[2][0] - 1)*sc[0] + (udata->TB[2][1] - 1)*sc[5] + (udata->TB[2][2] - 1)*sc[10] + (udata->TB[2][3] - 1)*sc[11] + (udata->TB[2][4] - 1)*sc[12] + (udata->TB[2][5] - 1)*sc[18] + (udata->TB[2][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = udata->prefactor_units[2] * udata->fwd_A[2]
                * exp(udata->fwd_beta[2] * tc[0] - udata->activation_units[2] * udata->fwd_Ea[2] * invT);
    dlnkfdT = udata->fwd_beta[2] * invT + udata->activation_units[2] * udata->fwd_Ea[2] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[2] * exp(udata->low_beta[2] * tc[0] - udata->activation_units[2] * udata->low_Ea[2] * invT);
    Pr = udata->phase_units[2] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[2] * invT + udata->activation_units[2] * udata->low_Ea[2] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[2]) > 1.e-100 ? (1.-udata->troe_a[2])*exp(-T/udata->troe_Tsss[2]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[2]) > 1.e-100 ? udata->troe_a[2] * exp(-T/udata->troe_Ts[2]) : 0.);
    Fcent3 = (udata->troe_len[2] == 4 ? exp(-udata->troe_Tss[2] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[2]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[2] : 0.)
      + (fabs(udata->troe_Ts[2]) > 1.e-100 ? -Fcent2/udata->troe_Ts[2] : 0.)
      + (udata->troe_len[2] == 4 ? Fcent3*udata->troe_Tss[2]*invT2 : 0.) );
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
        dqdci = (udata->TB[2][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[13];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[2][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[2][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[2][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[2][4] - 1)*dcdc_fac;
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
        dqdci = (udata->TB[2][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[2][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[453] -= dqdci;              /* dwdot[HCO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[2][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[13];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[2][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = udata->TB[2][2]*dcdc_fac;
        dqdc[11] = udata->TB[2][3]*dcdc_fac;
        dqdc[12] = udata->TB[2][4]*dcdc_fac;
        dqdc[13] = dcdc_fac + k_f*sc[1];
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = udata->TB[2][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = udata->TB[2][6]*dcdc_fac;
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
    alpha = mixture + (udata->TB[3][0] - 1)*sc[0] + (udata->TB[3][1] - 1)*sc[5] + (udata->TB[3][2] - 1)*sc[10] + (udata->TB[3][3] - 1)*sc[11] + (udata->TB[3][4] - 1)*sc[12] + (udata->TB[3][5] - 1)*sc[18];
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = udata->prefactor_units[3] * udata->fwd_A[3]
                * exp(udata->fwd_beta[3] * tc[0] - udata->activation_units[3] * udata->fwd_Ea[3] * invT);
    dlnkfdT = udata->fwd_beta[3] * invT + udata->activation_units[3] * udata->fwd_Ea[3] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[3] * exp(udata->low_beta[3] * tc[0] - udata->activation_units[3] * udata->low_Ea[3] * invT);
    Pr = udata->phase_units[3] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[3] * invT + udata->activation_units[3] * udata->low_Ea[3] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[3]) > 1.e-100 ? (1.-udata->troe_a[3])*exp(-T/udata->troe_Tsss[3]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[3]) > 1.e-100 ? udata->troe_a[3] * exp(-T/udata->troe_Ts[3]) : 0.);
    Fcent3 = (udata->troe_len[3] == 4 ? exp(-udata->troe_Tss[3] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[3]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[3] : 0.)
      + (fabs(udata->troe_Ts[3]) > 1.e-100 ? -Fcent2/udata->troe_Ts[3] : 0.)
      + (udata->troe_len[3] == 4 ? Fcent3*udata->troe_Tss[3]*invT2 : 0.) );
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
        dqdci = (udata->TB[3][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        J[15] += dqdci;               /* dwdot[CH3O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[14];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
        J[37] += dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[3][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[125] += dqdci;              /* dwdot[CH3O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[3][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[235] += dqdci;              /* dwdot[CH3O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[3][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[256] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        J[257] += dqdci;              /* dwdot[CH3O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[3][4] - 1)*dcdc_fac;
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
        dqdci = (udata->TB[3][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[410] -= dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[411] += dqdci;              /* dwdot[CH3O]/d[C2H6] */
    }
    else {
        dqdc[0] = udata->TB[3][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[14];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[3][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = udata->TB[3][2]*dcdc_fac;
        dqdc[11] = udata->TB[3][3]*dcdc_fac;
        dqdc[12] = udata->TB[3][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac + k_f*sc[1];
        dqdc[15] = dcdc_fac - k_r;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = udata->TB[3][5]*dcdc_fac;
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
    alpha = mixture + (udata->TB[4][0] - 1)*sc[0] + (udata->TB[4][1] - 1)*sc[5] + (udata->TB[4][2] - 1)*sc[10] + (udata->TB[4][3] - 1)*sc[11] + (udata->TB[4][4] - 1)*sc[12] + (udata->TB[4][5] - 1)*sc[18] + (udata->TB[4][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = udata->prefactor_units[4] * udata->fwd_A[4]
                * exp(udata->fwd_beta[4] * tc[0] - udata->activation_units[4] * udata->fwd_Ea[4] * invT);
    dlnkfdT = udata->fwd_beta[4] * invT + udata->activation_units[4] * udata->fwd_Ea[4] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[4] * exp(udata->low_beta[4] * tc[0] - udata->activation_units[4] * udata->low_Ea[4] * invT);
    Pr = udata->phase_units[4] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[4] * invT + udata->activation_units[4] * udata->low_Ea[4] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[4]) > 1.e-100 ? (1.-udata->troe_a[4])*exp(-T/udata->troe_Tsss[4]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[4]) > 1.e-100 ? udata->troe_a[4] * exp(-T/udata->troe_Ts[4]) : 0.);
    Fcent3 = (udata->troe_len[4] == 4 ? exp(-udata->troe_Tss[4] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[4]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[4] : 0.)
      + (fabs(udata->troe_Ts[4]) > 1.e-100 ? -Fcent2/udata->troe_Ts[4] : 0.)
      + (udata->troe_len[4] == 4 ? Fcent3*udata->troe_Tss[4]*invT2 : 0.) );
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
        dqdci = (udata->TB[4][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[16];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[38] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[4][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[126] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[4][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[236] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[4][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[258] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[259] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[4][4] - 1)*dcdc_fac;
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
        dqdci = (udata->TB[4][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[4][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[456] -= dqdci;              /* dwdot[C2H4]/d[AR] */
        J[457] += dqdci;              /* dwdot[C2H5]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[4][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[16];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[4][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = udata->TB[4][2]*dcdc_fac;
        dqdc[11] = udata->TB[4][3]*dcdc_fac;
        dqdc[12] = udata->TB[4][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f*sc[1];
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = udata->TB[4][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = udata->TB[4][6]*dcdc_fac;
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
    alpha = mixture + (udata->TB[5][0] - 1)*sc[0] + (udata->TB[5][1] - 1)*sc[5] + (udata->TB[5][2] - 1)*sc[10] + (udata->TB[5][3] - 1)*sc[11] + (udata->TB[5][4] - 1)*sc[12] + (udata->TB[5][5] - 1)*sc[18] + (udata->TB[5][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = udata->prefactor_units[5] * udata->fwd_A[5]
                * exp(udata->fwd_beta[5] * tc[0] - udata->activation_units[5] * udata->fwd_Ea[5] * invT);
    dlnkfdT = udata->fwd_beta[5] * invT + udata->activation_units[5] * udata->fwd_Ea[5] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[5] * exp(udata->low_beta[5] * tc[0] - udata->activation_units[5] * udata->low_Ea[5] * invT);
    Pr = udata->phase_units[5] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[5] * invT + udata->activation_units[5] * udata->low_Ea[5] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[5]) > 1.e-100 ? (1.-udata->troe_a[5])*exp(-T/udata->troe_Tsss[5]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[5]) > 1.e-100 ? udata->troe_a[5] * exp(-T/udata->troe_Ts[5]) : 0.);
    Fcent3 = (udata->troe_len[5] == 4 ? exp(-udata->troe_Tss[5] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[5]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[5] : 0.)
      + (fabs(udata->troe_Ts[5]) > 1.e-100 ? -Fcent2/udata->troe_Ts[5] : 0.)
      + (udata->troe_len[5] == 4 ? Fcent3*udata->troe_Tss[5]*invT2 : 0.) );
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
        dqdci = (udata->TB[5][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[39] -= dqdci;               /* dwdot[C2H5]/d[H] */
        J[40] += dqdci;               /* dwdot[C2H6]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[5][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[127] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[5][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[237] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[5][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[259] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[5][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[281] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H5] */
        dqdci =  + k_f*sc[1];
        J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[392] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[5][5] - 1)*dcdc_fac - k_r;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[413] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[5][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[457] -= dqdci;              /* dwdot[C2H5]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[5][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[17];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[5][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = udata->TB[5][2]*dcdc_fac;
        dqdc[11] = udata->TB[5][3]*dcdc_fac;
        dqdc[12] = udata->TB[5][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[1];
        dqdc[18] = udata->TB[5][5]*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = udata->TB[5][6]*dcdc_fac;
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
    alpha = mixture + (udata->TB[6][0] - 1)*sc[0] + (udata->TB[6][1] - 1)*sc[5] + (udata->TB[6][2] - 1)*sc[10] + (udata->TB[6][3] - 1)*sc[11] + (udata->TB[6][4] - 1)*sc[12] + (udata->TB[6][5] - 1)*sc[18] + (udata->TB[6][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = udata->prefactor_units[6] * udata->fwd_A[6]
                * exp(udata->fwd_beta[6] * tc[0] - udata->activation_units[6] * udata->fwd_Ea[6] * invT);
    dlnkfdT = udata->fwd_beta[6] * invT + udata->activation_units[6] * udata->fwd_Ea[6] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[6] * exp(udata->low_beta[6] * tc[0] - udata->activation_units[6] * udata->low_Ea[6] * invT);
    Pr = udata->phase_units[6] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[6] * invT + udata->activation_units[6] * udata->low_Ea[6] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[6]) > 1.e-100 ? (1.-udata->troe_a[6])*exp(-T/udata->troe_Tsss[6]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[6]) > 1.e-100 ? udata->troe_a[6] * exp(-T/udata->troe_Ts[6]) : 0.);
    Fcent3 = (udata->troe_len[6] == 4 ? exp(-udata->troe_Tss[6] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[6]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[6] : 0.)
      + (fabs(udata->troe_Ts[6]) > 1.e-100 ? -Fcent2/udata->troe_Ts[6] : 0.)
      + (udata->troe_len[6] == 4 ? Fcent3*udata->troe_Tss[6]*invT2 : 0.) );
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
        dqdci = (udata->TB[6][0] - 1)*dcdc_fac + k_f*sc[11];
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[6][1] - 1)*dcdc_fac;
        J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[6][2] - 1)*dcdc_fac;
        J[220] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[6][3] - 1)*dcdc_fac + k_f*sc[0];
        J[242] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[6][4] - 1)*dcdc_fac;
        J[264] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[308] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[319] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[6][5] - 1)*dcdc_fac;
        J[396] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[6][6] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H2]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[6][0]*dcdc_fac + k_f*sc[11];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[6][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = udata->TB[6][2]*dcdc_fac;
        dqdc[11] = udata->TB[6][3]*dcdc_fac + k_f*sc[0];
        dqdc[12] = udata->TB[6][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = udata->TB[6][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = udata->TB[6][6]*dcdc_fac;
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
    alpha = mixture + (udata->TB[7][0] - 1)*sc[0] + (udata->TB[7][1] - 1)*sc[5] + (udata->TB[7][2] - 1)*sc[10] + (udata->TB[7][3] - 1)*sc[11] + (udata->TB[7][4] - 1)*sc[12] + (udata->TB[7][5] - 1)*sc[18] + (udata->TB[7][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[9]*sc[9];
    k_f = udata->prefactor_units[7] * udata->fwd_A[7]
                * exp(udata->fwd_beta[7] * tc[0] - udata->activation_units[7] * udata->fwd_Ea[7] * invT);
    dlnkfdT = udata->fwd_beta[7] * invT + udata->activation_units[7] * udata->fwd_Ea[7] * invT2;
    /* pressure-fall-off */
    k_0 = udata->low_A[7] * exp(udata->low_beta[7] * tc[0] - udata->activation_units[7] * udata->low_Ea[7] * invT);
    Pr = udata->phase_units[7] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = udata->low_beta[7] * invT + udata->activation_units[7] * udata->low_Ea[7] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(udata->troe_Tsss[7]) > 1.e-100 ? (1.-udata->troe_a[7])*exp(-T/udata->troe_Tsss[7]) : 0.);
    Fcent2 = (fabs(udata->troe_Ts[7]) > 1.e-100 ? udata->troe_a[7] * exp(-T/udata->troe_Ts[7]) : 0.);
    Fcent3 = (udata->troe_len[7] == 4 ? exp(-udata->troe_Tss[7] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(udata->troe_Tsss[7]) > 1.e-100 ? -Fcent1/udata->troe_Tsss[7] : 0.)
      + (fabs(udata->troe_Ts[7]) > 1.e-100 ? -Fcent2/udata->troe_Ts[7] : 0.)
      + (udata->troe_len[7] == 4 ? Fcent3*udata->troe_Tss[7]*invT2 : 0.) );
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
        dqdci = (udata->TB[7][0] - 1)*dcdc_fac;
        J[9] += -2 * dqdci;           /* dwdot[CH3]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[7][1] - 1)*dcdc_fac;
        J[119] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2*sc[9];
        J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[216] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[7][2] - 1)*dcdc_fac;
        J[229] += -2 * dqdci;         /* dwdot[CH3]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[7][3] - 1)*dcdc_fac;
        J[251] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[7][4] - 1)*dcdc_fac;
        J[273] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[7][5] - 1)*dcdc_fac - k_r;
        J[405] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[7][6] - 1)*dcdc_fac;
        J[449] += -2 * dqdci;         /* dwdot[CH3]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[7][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = udata->TB[7][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*2*sc[9];
        dqdc[10] = udata->TB[7][2]*dcdc_fac;
        dqdc[11] = udata->TB[7][3]*dcdc_fac;
        dqdc[12] = udata->TB[7][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = udata->TB[7][5]*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = udata->TB[7][6]*dcdc_fac;
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
    alpha = mixture + (udata->TB[8][0] - 1)*sc[0] + (udata->TB[8][1] - 1)*sc[5] + (udata->TB[8][2] - 1)*sc[10] + (udata->TB[8][3] - 1)*sc[11] + (udata->TB[8][4] - 1)*sc[12] + (udata->TB[8][5] - 1)*sc[18] + (udata->TB[8][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = udata->prefactor_units[8] * udata->fwd_A[8]
                * exp(udata->fwd_beta[8] * tc[0] - udata->activation_units[8] * udata->fwd_Ea[8] * invT);
    dlnkfdT = udata->fwd_beta[8] * invT + udata->activation_units[8] * udata->fwd_Ea[8] * invT2;
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
        dqdci = (udata->TB[8][0] - 1)*q_nocor;
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
        dqdci = (udata->TB[8][1] - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[114] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[8][2] - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[8][3] - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[246] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[8][4] - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[8][5] - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[8][6] - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[444] += dqdci;              /* dwdot[OH]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[8][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = udata->TB[8][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = udata->TB[8][2]*q_nocor;
        dqdc[11] = udata->TB[8][3]*q_nocor;
        dqdc[12] = udata->TB[8][4]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = udata->TB[8][5]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = udata->TB[8][6]*q_nocor;
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
    alpha = mixture + (udata->TB[9][0] - 1)*sc[0] + (udata->TB[9][1] - 1)*sc[3] + (udata->TB[9][2] - 1)*sc[5] + (udata->TB[9][3] - 1)*sc[10] + (udata->TB[9][4] - 1)*sc[11] + (udata->TB[9][5] - 1)*sc[12] + (udata->TB[9][6] - 1)*sc[18] + (udata->TB[9][7] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = udata->prefactor_units[9] * udata->fwd_A[9]
                * exp(udata->fwd_beta[9] * tc[0] - udata->activation_units[9] * udata->fwd_Ea[9] * invT);
    dlnkfdT = udata->fwd_beta[9] * invT + udata->activation_units[9] * udata->fwd_Ea[9] * invT2;
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
        dqdci = (udata->TB[9][0] - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[12] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[11];
        J[46] -= dqdci;               /* dwdot[O]/d[O] */
        J[55] -= dqdci;               /* dwdot[CO]/d[O] */
        J[56] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[O2] */
        dqdci = (udata->TB[9][1] - 1)*q_nocor;
        J[68] -= dqdci;               /* dwdot[O]/d[O2] */
        J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
        J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[9][2] - 1)*q_nocor;
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[122] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[9][3] - 1)*q_nocor;
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[232] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[9][4] - 1)*q_nocor + k_f*sc[2];
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[9][5] - 1)*q_nocor - k_r;
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[9][6] - 1)*q_nocor;
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[408] += dqdci;              /* dwdot[CO2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[9][7] - 1)*q_nocor;
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[452] += dqdci;              /* dwdot[CO2]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[9][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[11];
        dqdc[3] = udata->TB[9][1]*q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = udata->TB[9][2]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = udata->TB[9][3]*q_nocor;
        dqdc[11] = udata->TB[9][4]*q_nocor + k_f*sc[2];
        dqdc[12] = udata->TB[9][5]*q_nocor - k_r;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = udata->TB[9][6]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = udata->TB[9][7]*q_nocor;
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
    alpha = mixture + (udata->TB[10][0] - 1)*sc[3] + (udata->TB[10][1] - 1)*sc[5] + (udata->TB[10][2] - 1)*sc[11] + (udata->TB[10][3] - 1)*sc[12] + (udata->TB[10][4] - 1)*sc[18] + (udata->TB[10][5] - 1)*sc[19] + (udata->TB[10][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = udata->prefactor_units[10] * udata->fwd_A[10]
                * exp(udata->fwd_beta[10] * tc[0] - udata->activation_units[10] * udata->fwd_Ea[10] * invT);
    dlnkfdT = udata->fwd_beta[10] * invT + udata->activation_units[10] * udata->fwd_Ea[10] * invT2;
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
        dqdci = (udata->TB[10][0] - 1)*q_nocor + k_f*sc[1];
        J[67] -= dqdci;               /* dwdot[H]/d[O2] */
        J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[10][1] - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (udata->TB[10][2] - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[10][3] - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[270] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[10][4] - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[399] -= dqdci;              /* dwdot[O2]/d[C2H6] */
        J[402] += dqdci;              /* dwdot[HO2]/d[C2H6] */
        /* d()/d[N2] */
        dqdci = (udata->TB[10][5] - 1)*q_nocor;
        J[419] -= dqdci;              /* dwdot[H]/d[N2] */
        J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
        /* d()/d[AR] */
        dqdci = (udata->TB[10][6] - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] = udata->TB[10][0]*q_nocor + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[5] = udata->TB[10][1]*q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = udata->TB[10][2]*q_nocor;
        dqdc[12] = udata->TB[10][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = udata->TB[10][4]*q_nocor;
        dqdc[19] = udata->TB[10][5]*q_nocor;
        dqdc[20] = udata->TB[10][6]*q_nocor;
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
    alpha = mixture + (udata->TB[11][0] - 1)*sc[0] + (udata->TB[11][1] - 1)*sc[5] + (udata->TB[11][2] - 1)*sc[10] + (udata->TB[11][3] - 1)*sc[12] + (udata->TB[11][4] - 1)*sc[18] + (udata->TB[11][5] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[1];
    k_f = udata->prefactor_units[11] * udata->fwd_A[11]
                * exp(udata->fwd_beta[11] * tc[0] - udata->activation_units[11] * udata->fwd_Ea[11] * invT);
    dlnkfdT = udata->fwd_beta[11] * invT + udata->activation_units[11] * udata->fwd_Ea[11] * invT2;
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
        dqdci = (udata->TB[11][0] - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2*sc[1];
        J[22] += dqdci;               /* dwdot[H2]/d[H] */
        J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[11][1] - 1)*q_nocor;
        J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[11][2] - 1)*q_nocor;
        J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[221] += -2 * dqdci;         /* dwdot[H]/d[CH4] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[11][3] - 1)*q_nocor;
        J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[11][4] - 1)*q_nocor;
        J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[397] += -2 * dqdci;         /* dwdot[H]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[11][5] - 1)*q_nocor;
        J[440] += dqdci;              /* dwdot[H2]/d[AR] */
        J[441] += -2 * dqdci;         /* dwdot[H]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[11][0]*q_nocor - k_r;
        dqdc[1] = q_nocor + k_f*2*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = udata->TB[11][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = udata->TB[11][2]*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = udata->TB[11][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = udata->TB[11][4]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = udata->TB[11][5]*q_nocor;
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
    alpha = mixture + (udata->TB[12][0] - 1)*sc[0] + (udata->TB[12][1] - 1)*sc[5] + (udata->TB[12][2] - 1)*sc[10] + (udata->TB[12][3] - 1)*sc[18] + (udata->TB[12][4] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = udata->prefactor_units[12] * udata->fwd_A[12]
                * exp(udata->fwd_beta[12] * tc[0] - udata->activation_units[12] * udata->fwd_Ea[12] * invT);
    dlnkfdT = udata->fwd_beta[12] * invT + udata->activation_units[12] * udata->fwd_Ea[12] * invT2;
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
        dqdci = (udata->TB[12][0] - 1)*q_nocor;
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
        dqdci = (udata->TB[12][1] - 1)*q_nocor - k_r;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[12][2] - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[12][3] - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (udata->TB[12][4] - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[444] -= dqdci;              /* dwdot[OH]/d[AR] */
        J[445] += dqdci;              /* dwdot[H2O]/d[AR] */
    }
    else {
        dqdc[0] = udata->TB[12][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = udata->TB[12][1]*q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = udata->TB[12][2]*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = udata->TB[12][3]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = udata->TB[12][4]*q_nocor;
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
    alpha = mixture + (udata->TB[13][0] - 1)*sc[0] + (udata->TB[13][1] - 1)*sc[5] + (udata->TB[13][2] - 1)*sc[10] + (udata->TB[13][3] - 1)*sc[11] + (udata->TB[13][4] - 1)*sc[12] + (udata->TB[13][5] - 1)*sc[18];
    /* forward */
    phi_f = sc[13];
    k_f = udata->prefactor_units[13] * udata->fwd_A[13]
                * exp(udata->fwd_beta[13] * tc[0] - udata->activation_units[13] * udata->fwd_Ea[13] * invT);
    dlnkfdT = udata->fwd_beta[13] * invT + udata->activation_units[13] * udata->fwd_Ea[13] * invT2;
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
        dqdci = (udata->TB[13][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[11] += dqdci;               /* dwdot[CO]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[11];
        J[23] += dqdci;               /* dwdot[H]/d[H] */
        J[33] += dqdci;               /* dwdot[CO]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2O] */
        dqdci = (udata->TB[13][1] - 1)*q_nocor;
        J[111] += dqdci;              /* dwdot[H]/d[H2O] */
        J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (udata->TB[13][2] - 1)*q_nocor;
        J[221] += dqdci;              /* dwdot[H]/d[CH4] */
        J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (udata->TB[13][3] - 1)*q_nocor - k_r*sc[1];
        J[243] += dqdci;              /* dwdot[H]/d[CO] */
        J[253] += dqdci;              /* dwdot[CO]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (udata->TB[13][4] - 1)*q_nocor;
        J[265] += dqdci;              /* dwdot[H]/d[CO2] */
        J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[287] += dqdci;              /* dwdot[H]/d[HCO] */
        J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[C2H6] */
        dqdci = (udata->TB[13][5] - 1)*q_nocor;
        J[397] += dqdci;              /* dwdot[H]/d[C2H6] */
        J[407] += dqdci;              /* dwdot[CO]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
    }
    else {
        dqdc[0] = udata->TB[13][0]*q_nocor;
        dqdc[1] = q_nocor - k_r*sc[11];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = udata->TB[13][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = udata->TB[13][2]*q_nocor;
        dqdc[11] = udata->TB[13][3]*q_nocor - k_r*sc[1];
        dqdc[12] = udata->TB[13][4]*q_nocor;
        dqdc[13] = q_nocor + k_f;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = udata->TB[13][5]*q_nocor;
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
    k_f = udata->prefactor_units[14] * udata->fwd_A[14]
                * exp(udata->fwd_beta[14] * tc[0] - udata->activation_units[14] * udata->fwd_Ea[14] * invT);
    dlnkfdT = udata->fwd_beta[14] * invT + udata->activation_units[14] * udata->fwd_Ea[14] * invT2;
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
    k_f = udata->prefactor_units[15] * udata->fwd_A[15]
                * exp(udata->fwd_beta[15] * tc[0] - udata->activation_units[15] * udata->fwd_Ea[15] * invT);
    dlnkfdT = udata->fwd_beta[15] * invT + udata->activation_units[15] * udata->fwd_Ea[15] * invT2;
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
    k_f = udata->prefactor_units[16] * udata->fwd_A[16]
                * exp(udata->fwd_beta[16] * tc[0] - udata->activation_units[16] * udata->fwd_Ea[16] * invT);
    dlnkfdT = udata->fwd_beta[16] * invT + udata->activation_units[16] * udata->fwd_Ea[16] * invT2;
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
    k_f = udata->prefactor_units[17] * udata->fwd_A[17]
                * exp(udata->fwd_beta[17] * tc[0] - udata->activation_units[17] * udata->fwd_Ea[17] * invT);
    dlnkfdT = udata->fwd_beta[17] * invT + udata->activation_units[17] * udata->fwd_Ea[17] * invT2;
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
    k_f = udata->prefactor_units[18] * udata->fwd_A[18]
                * exp(udata->fwd_beta[18] * tc[0] - udata->activation_units[18] * udata->fwd_Ea[18] * invT);
    dlnkfdT = udata->fwd_beta[18] * invT + udata->activation_units[18] * udata->fwd_Ea[18] * invT2;
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
    k_f = udata->prefactor_units[19] * udata->fwd_A[19]
                * exp(udata->fwd_beta[19] * tc[0] - udata->activation_units[19] * udata->fwd_Ea[19] * invT);
    dlnkfdT = udata->fwd_beta[19] * invT + udata->activation_units[19] * udata->fwd_Ea[19] * invT2;
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
    k_f = udata->prefactor_units[20] * udata->fwd_A[20]
                * exp(udata->fwd_beta[20] * tc[0] - udata->activation_units[20] * udata->fwd_Ea[20] * invT);
    dlnkfdT = udata->fwd_beta[20] * invT + udata->activation_units[20] * udata->fwd_Ea[20] * invT2;
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
    k_f = udata->prefactor_units[21] * udata->fwd_A[21]
                * exp(udata->fwd_beta[21] * tc[0] - udata->activation_units[21] * udata->fwd_Ea[21] * invT);
    dlnkfdT = udata->fwd_beta[21] * invT + udata->activation_units[21] * udata->fwd_Ea[21] * invT2;
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
    k_f = udata->prefactor_units[22] * udata->fwd_A[22]
                * exp(udata->fwd_beta[22] * tc[0] - udata->activation_units[22] * udata->fwd_Ea[22] * invT);
    dlnkfdT = udata->fwd_beta[22] * invT + udata->activation_units[22] * udata->fwd_Ea[22] * invT2;
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
    k_f = udata->prefactor_units[23] * udata->fwd_A[23]
                * exp(udata->fwd_beta[23] * tc[0] - udata->activation_units[23] * udata->fwd_Ea[23] * invT);
    dlnkfdT = udata->fwd_beta[23] * invT + udata->activation_units[23] * udata->fwd_Ea[23] * invT2;
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
    k_f = udata->prefactor_units[24] * udata->fwd_A[24]
                * exp(udata->fwd_beta[24] * tc[0] - udata->activation_units[24] * udata->fwd_Ea[24] * invT);
    dlnkfdT = udata->fwd_beta[24] * invT + udata->activation_units[24] * udata->fwd_Ea[24] * invT2;
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
    k_f = udata->prefactor_units[25] * udata->fwd_A[25]
                * exp(udata->fwd_beta[25] * tc[0] - udata->activation_units[25] * udata->fwd_Ea[25] * invT);
    dlnkfdT = udata->fwd_beta[25] * invT + udata->activation_units[25] * udata->fwd_Ea[25] * invT2;
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
    k_f = udata->prefactor_units[26] * udata->fwd_A[26]
                * exp(udata->fwd_beta[26] * tc[0] - udata->activation_units[26] * udata->fwd_Ea[26] * invT);
    dlnkfdT = udata->fwd_beta[26] * invT + udata->activation_units[26] * udata->fwd_Ea[26] * invT2;
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
    k_f = udata->prefactor_units[27] * udata->fwd_A[27]
                * exp(udata->fwd_beta[27] * tc[0] - udata->activation_units[27] * udata->fwd_Ea[27] * invT);
    dlnkfdT = udata->fwd_beta[27] * invT + udata->activation_units[27] * udata->fwd_Ea[27] * invT2;
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
    k_f = udata->prefactor_units[28] * udata->fwd_A[28]
                * exp(udata->fwd_beta[28] * tc[0] - udata->activation_units[28] * udata->fwd_Ea[28] * invT);
    dlnkfdT = udata->fwd_beta[28] * invT + udata->activation_units[28] * udata->fwd_Ea[28] * invT2;
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
    k_f = udata->prefactor_units[29] * udata->fwd_A[29]
                * exp(udata->fwd_beta[29] * tc[0] - udata->activation_units[29] * udata->fwd_Ea[29] * invT);
    dlnkfdT = udata->fwd_beta[29] * invT + udata->activation_units[29] * udata->fwd_Ea[29] * invT2;
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
    k_f = udata->prefactor_units[30] * udata->fwd_A[30]
                * exp(udata->fwd_beta[30] * tc[0] - udata->activation_units[30] * udata->fwd_Ea[30] * invT);
    dlnkfdT = udata->fwd_beta[30] * invT + udata->activation_units[30] * udata->fwd_Ea[30] * invT2;
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
    k_f = udata->prefactor_units[31] * udata->fwd_A[31]
                * exp(udata->fwd_beta[31] * tc[0] - udata->activation_units[31] * udata->fwd_Ea[31] * invT);
    dlnkfdT = udata->fwd_beta[31] * invT + udata->activation_units[31] * udata->fwd_Ea[31] * invT2;
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
    k_f = udata->prefactor_units[32] * udata->fwd_A[32]
                * exp(udata->fwd_beta[32] * tc[0] - udata->activation_units[32] * udata->fwd_Ea[32] * invT);
    dlnkfdT = udata->fwd_beta[32] * invT + udata->activation_units[32] * udata->fwd_Ea[32] * invT2;
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
    k_f = udata->prefactor_units[33] * udata->fwd_A[33]
                * exp(udata->fwd_beta[33] * tc[0] - udata->activation_units[33] * udata->fwd_Ea[33] * invT);
    dlnkfdT = udata->fwd_beta[33] * invT + udata->activation_units[33] * udata->fwd_Ea[33] * invT2;
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
    k_f = udata->prefactor_units[34] * udata->fwd_A[34]
                * exp(udata->fwd_beta[34] * tc[0] - udata->activation_units[34] * udata->fwd_Ea[34] * invT);
    dlnkfdT = udata->fwd_beta[34] * invT + udata->activation_units[34] * udata->fwd_Ea[34] * invT2;
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
    k_f = udata->prefactor_units[35] * udata->fwd_A[35]
                * exp(udata->fwd_beta[35] * tc[0] - udata->activation_units[35] * udata->fwd_Ea[35] * invT);
    dlnkfdT = udata->fwd_beta[35] * invT + udata->activation_units[35] * udata->fwd_Ea[35] * invT2;
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
    k_f = udata->prefactor_units[36] * udata->fwd_A[36]
                * exp(udata->fwd_beta[36] * tc[0] - udata->activation_units[36] * udata->fwd_Ea[36] * invT);
    dlnkfdT = udata->fwd_beta[36] * invT + udata->activation_units[36] * udata->fwd_Ea[36] * invT2;
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
    k_f = udata->prefactor_units[37] * udata->fwd_A[37]
                * exp(udata->fwd_beta[37] * tc[0] - udata->activation_units[37] * udata->fwd_Ea[37] * invT);
    dlnkfdT = udata->fwd_beta[37] * invT + udata->activation_units[37] * udata->fwd_Ea[37] * invT2;
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
    k_f = udata->prefactor_units[38] * udata->fwd_A[38]
                * exp(udata->fwd_beta[38] * tc[0] - udata->activation_units[38] * udata->fwd_Ea[38] * invT);
    dlnkfdT = udata->fwd_beta[38] * invT + udata->activation_units[38] * udata->fwd_Ea[38] * invT2;
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
    k_f = udata->prefactor_units[39] * udata->fwd_A[39]
                * exp(udata->fwd_beta[39] * tc[0] - udata->activation_units[39] * udata->fwd_Ea[39] * invT);
    dlnkfdT = udata->fwd_beta[39] * invT + udata->activation_units[39] * udata->fwd_Ea[39] * invT2;
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
    k_f = udata->prefactor_units[40] * udata->fwd_A[40]
                * exp(udata->fwd_beta[40] * tc[0] - udata->activation_units[40] * udata->fwd_Ea[40] * invT);
    dlnkfdT = udata->fwd_beta[40] * invT + udata->activation_units[40] * udata->fwd_Ea[40] * invT2;
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
    k_f = udata->prefactor_units[41] * udata->fwd_A[41]
                * exp(udata->fwd_beta[41] * tc[0] - udata->activation_units[41] * udata->fwd_Ea[41] * invT);
    dlnkfdT = udata->fwd_beta[41] * invT + udata->activation_units[41] * udata->fwd_Ea[41] * invT2;
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
    k_f = udata->prefactor_units[42] * udata->fwd_A[42]
                * exp(udata->fwd_beta[42] * tc[0] - udata->activation_units[42] * udata->fwd_Ea[42] * invT);
    dlnkfdT = udata->fwd_beta[42] * invT + udata->activation_units[42] * udata->fwd_Ea[42] * invT2;
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
    k_f = udata->prefactor_units[43] * udata->fwd_A[43]
                * exp(udata->fwd_beta[43] * tc[0] - udata->activation_units[43] * udata->fwd_Ea[43] * invT);
    dlnkfdT = udata->fwd_beta[43] * invT + udata->activation_units[43] * udata->fwd_Ea[43] * invT2;
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
    k_f = udata->prefactor_units[44] * udata->fwd_A[44]
                * exp(udata->fwd_beta[44] * tc[0] - udata->activation_units[44] * udata->fwd_Ea[44] * invT);
    dlnkfdT = udata->fwd_beta[44] * invT + udata->activation_units[44] * udata->fwd_Ea[44] * invT2;
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
    k_f = udata->prefactor_units[45] * udata->fwd_A[45]
                * exp(udata->fwd_beta[45] * tc[0] - udata->activation_units[45] * udata->fwd_Ea[45] * invT);
    dlnkfdT = udata->fwd_beta[45] * invT + udata->activation_units[45] * udata->fwd_Ea[45] * invT2;
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
    k_f = udata->prefactor_units[46] * udata->fwd_A[46]
                * exp(udata->fwd_beta[46] * tc[0] - udata->activation_units[46] * udata->fwd_Ea[46] * invT);
    dlnkfdT = udata->fwd_beta[46] * invT + udata->activation_units[46] * udata->fwd_Ea[46] * invT2;
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
    k_f = udata->prefactor_units[47] * udata->fwd_A[47]
                * exp(udata->fwd_beta[47] * tc[0] - udata->activation_units[47] * udata->fwd_Ea[47] * invT);
    dlnkfdT = udata->fwd_beta[47] * invT + udata->activation_units[47] * udata->fwd_Ea[47] * invT2;
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
    k_f = udata->prefactor_units[48] * udata->fwd_A[48]
                * exp(udata->fwd_beta[48] * tc[0] - udata->activation_units[48] * udata->fwd_Ea[48] * invT);
    dlnkfdT = udata->fwd_beta[48] * invT + udata->activation_units[48] * udata->fwd_Ea[48] * invT2;
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
    k_f = udata->prefactor_units[49] * udata->fwd_A[49]
                * exp(udata->fwd_beta[49] * tc[0] - udata->activation_units[49] * udata->fwd_Ea[49] * invT);
    dlnkfdT = udata->fwd_beta[49] * invT + udata->activation_units[49] * udata->fwd_Ea[49] * invT2;
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
    k_f = udata->prefactor_units[50] * udata->fwd_A[50]
                * exp(udata->fwd_beta[50] * tc[0] - udata->activation_units[50] * udata->fwd_Ea[50] * invT);
    dlnkfdT = udata->fwd_beta[50] * invT + udata->activation_units[50] * udata->fwd_Ea[50] * invT2;
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
    k_f = udata->prefactor_units[51] * udata->fwd_A[51]
                * exp(udata->fwd_beta[51] * tc[0] - udata->activation_units[51] * udata->fwd_Ea[51] * invT);
    dlnkfdT = udata->fwd_beta[51] * invT + udata->activation_units[51] * udata->fwd_Ea[51] * invT2;
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
    k_f = udata->prefactor_units[52] * udata->fwd_A[52]
                * exp(udata->fwd_beta[52] * tc[0] - udata->activation_units[52] * udata->fwd_Ea[52] * invT);
    dlnkfdT = udata->fwd_beta[52] * invT + udata->activation_units[52] * udata->fwd_Ea[52] * invT2;
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
    k_f = udata->prefactor_units[53] * udata->fwd_A[53]
                * exp(udata->fwd_beta[53] * tc[0] - udata->activation_units[53] * udata->fwd_Ea[53] * invT);
    dlnkfdT = udata->fwd_beta[53] * invT + udata->activation_units[53] * udata->fwd_Ea[53] * invT2;
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
    k_f = udata->prefactor_units[54] * udata->fwd_A[54]
                * exp(udata->fwd_beta[54] * tc[0] - udata->activation_units[54] * udata->fwd_Ea[54] * invT);
    dlnkfdT = udata->fwd_beta[54] * invT + udata->activation_units[54] * udata->fwd_Ea[54] * invT2;
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
    k_f = udata->prefactor_units[55] * udata->fwd_A[55]
                * exp(udata->fwd_beta[55] * tc[0] - udata->activation_units[55] * udata->fwd_Ea[55] * invT);
    dlnkfdT = udata->fwd_beta[55] * invT + udata->activation_units[55] * udata->fwd_Ea[55] * invT2;
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
    k_f = udata->prefactor_units[56] * udata->fwd_A[56]
                * exp(udata->fwd_beta[56] * tc[0] - udata->activation_units[56] * udata->fwd_Ea[56] * invT);
    dlnkfdT = udata->fwd_beta[56] * invT + udata->activation_units[56] * udata->fwd_Ea[56] * invT2;
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
    k_f = udata->prefactor_units[57] * udata->fwd_A[57]
                * exp(udata->fwd_beta[57] * tc[0] - udata->activation_units[57] * udata->fwd_Ea[57] * invT);
    dlnkfdT = udata->fwd_beta[57] * invT + udata->activation_units[57] * udata->fwd_Ea[57] * invT2;
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
    k_f = udata->prefactor_units[58] * udata->fwd_A[58]
                * exp(udata->fwd_beta[58] * tc[0] - udata->activation_units[58] * udata->fwd_Ea[58] * invT);
    dlnkfdT = udata->fwd_beta[58] * invT + udata->activation_units[58] * udata->fwd_Ea[58] * invT2;
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
    k_f = udata->prefactor_units[59] * udata->fwd_A[59]
                * exp(udata->fwd_beta[59] * tc[0] - udata->activation_units[59] * udata->fwd_Ea[59] * invT);
    dlnkfdT = udata->fwd_beta[59] * invT + udata->activation_units[59] * udata->fwd_Ea[59] * invT2;
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
    k_f = udata->prefactor_units[60] * udata->fwd_A[60]
                * exp(udata->fwd_beta[60] * tc[0] - udata->activation_units[60] * udata->fwd_Ea[60] * invT);
    dlnkfdT = udata->fwd_beta[60] * invT + udata->activation_units[60] * udata->fwd_Ea[60] * invT2;
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
    k_f = udata->prefactor_units[61] * udata->fwd_A[61]
                * exp(udata->fwd_beta[61] * tc[0] - udata->activation_units[61] * udata->fwd_Ea[61] * invT);
    dlnkfdT = udata->fwd_beta[61] * invT + udata->activation_units[61] * udata->fwd_Ea[61] * invT2;
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
    k_f = udata->prefactor_units[62] * udata->fwd_A[62]
                * exp(udata->fwd_beta[62] * tc[0] - udata->activation_units[62] * udata->fwd_Ea[62] * invT);
    dlnkfdT = udata->fwd_beta[62] * invT + udata->activation_units[62] * udata->fwd_Ea[62] * invT2;
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
    k_f = udata->prefactor_units[63] * udata->fwd_A[63]
                * exp(udata->fwd_beta[63] * tc[0] - udata->activation_units[63] * udata->fwd_Ea[63] * invT);
    dlnkfdT = udata->fwd_beta[63] * invT + udata->activation_units[63] * udata->fwd_Ea[63] * invT2;
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
    k_f = udata->prefactor_units[64] * udata->fwd_A[64]
                * exp(udata->fwd_beta[64] * tc[0] - udata->activation_units[64] * udata->fwd_Ea[64] * invT);
    dlnkfdT = udata->fwd_beta[64] * invT + udata->activation_units[64] * udata->fwd_Ea[64] * invT2;
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
    k_f = udata->prefactor_units[65] * udata->fwd_A[65]
                * exp(udata->fwd_beta[65] * tc[0] - udata->activation_units[65] * udata->fwd_Ea[65] * invT);
    dlnkfdT = udata->fwd_beta[65] * invT + udata->activation_units[65] * udata->fwd_Ea[65] * invT2;
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
    k_f = udata->prefactor_units[66] * udata->fwd_A[66]
                * exp(udata->fwd_beta[66] * tc[0] - udata->activation_units[66] * udata->fwd_Ea[66] * invT);
    dlnkfdT = udata->fwd_beta[66] * invT + udata->activation_units[66] * udata->fwd_Ea[66] * invT2;
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
    k_f = udata->prefactor_units[67] * udata->fwd_A[67]
                * exp(udata->fwd_beta[67] * tc[0] - udata->activation_units[67] * udata->fwd_Ea[67] * invT);
    dlnkfdT = udata->fwd_beta[67] * invT + udata->activation_units[67] * udata->fwd_Ea[67] * invT2;
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
    k_f = udata->prefactor_units[68] * udata->fwd_A[68]
                * exp(udata->fwd_beta[68] * tc[0] - udata->activation_units[68] * udata->fwd_Ea[68] * invT);
    dlnkfdT = udata->fwd_beta[68] * invT + udata->activation_units[68] * udata->fwd_Ea[68] * invT2;
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
    k_f = udata->prefactor_units[69] * udata->fwd_A[69]
                * exp(udata->fwd_beta[69] * tc[0] - udata->activation_units[69] * udata->fwd_Ea[69] * invT);
    dlnkfdT = udata->fwd_beta[69] * invT + udata->activation_units[69] * udata->fwd_Ea[69] * invT2;
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
    k_f = udata->prefactor_units[70] * udata->fwd_A[70]
                * exp(udata->fwd_beta[70] * tc[0] - udata->activation_units[70] * udata->fwd_Ea[70] * invT);
    dlnkfdT = udata->fwd_beta[70] * invT + udata->activation_units[70] * udata->fwd_Ea[70] * invT2;
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
    k_f = udata->prefactor_units[71] * udata->fwd_A[71]
                * exp(udata->fwd_beta[71] * tc[0] - udata->activation_units[71] * udata->fwd_Ea[71] * invT);
    dlnkfdT = udata->fwd_beta[71] * invT + udata->activation_units[71] * udata->fwd_Ea[71] * invT2;
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
    k_f = udata->prefactor_units[72] * udata->fwd_A[72]
                * exp(udata->fwd_beta[72] * tc[0] - udata->activation_units[72] * udata->fwd_Ea[72] * invT);
    dlnkfdT = udata->fwd_beta[72] * invT + udata->activation_units[72] * udata->fwd_Ea[72] * invT2;
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
    k_f = udata->prefactor_units[73] * udata->fwd_A[73]
                * exp(udata->fwd_beta[73] * tc[0] - udata->activation_units[73] * udata->fwd_Ea[73] * invT);
    dlnkfdT = udata->fwd_beta[73] * invT + udata->activation_units[73] * udata->fwd_Ea[73] * invT2;
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
    k_f = udata->prefactor_units[74] * udata->fwd_A[74]
                * exp(udata->fwd_beta[74] * tc[0] - udata->activation_units[74] * udata->fwd_Ea[74] * invT);
    dlnkfdT = udata->fwd_beta[74] * invT + udata->activation_units[74] * udata->fwd_Ea[74] * invT2;
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
    k_f = udata->prefactor_units[75] * udata->fwd_A[75]
                * exp(udata->fwd_beta[75] * tc[0] - udata->activation_units[75] * udata->fwd_Ea[75] * invT);
    dlnkfdT = udata->fwd_beta[75] * invT + udata->activation_units[75] * udata->fwd_Ea[75] * invT2;
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
    k_f = udata->prefactor_units[76] * udata->fwd_A[76]
                * exp(udata->fwd_beta[76] * tc[0] - udata->activation_units[76] * udata->fwd_Ea[76] * invT);
    dlnkfdT = udata->fwd_beta[76] * invT + udata->activation_units[76] * udata->fwd_Ea[76] * invT2;
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
    k_f = udata->prefactor_units[77] * udata->fwd_A[77]
                * exp(udata->fwd_beta[77] * tc[0] - udata->activation_units[77] * udata->fwd_Ea[77] * invT);
    dlnkfdT = udata->fwd_beta[77] * invT + udata->activation_units[77] * udata->fwd_Ea[77] * invT2;
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
    k_f = udata->prefactor_units[78] * udata->fwd_A[78]
                * exp(udata->fwd_beta[78] * tc[0] - udata->activation_units[78] * udata->fwd_Ea[78] * invT);
    dlnkfdT = udata->fwd_beta[78] * invT + udata->activation_units[78] * udata->fwd_Ea[78] * invT2;
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
    k_f = udata->prefactor_units[79] * udata->fwd_A[79]
                * exp(udata->fwd_beta[79] * tc[0] - udata->activation_units[79] * udata->fwd_Ea[79] * invT);
    dlnkfdT = udata->fwd_beta[79] * invT + udata->activation_units[79] * udata->fwd_Ea[79] * invT2;
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
    k_f = udata->prefactor_units[80] * udata->fwd_A[80]
                * exp(udata->fwd_beta[80] * tc[0] - udata->activation_units[80] * udata->fwd_Ea[80] * invT);
    dlnkfdT = udata->fwd_beta[80] * invT + udata->activation_units[80] * udata->fwd_Ea[80] * invT2;
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
    k_f = udata->prefactor_units[81] * udata->fwd_A[81]
                * exp(udata->fwd_beta[81] * tc[0] - udata->activation_units[81] * udata->fwd_Ea[81] * invT);
    dlnkfdT = udata->fwd_beta[81] * invT + udata->activation_units[81] * udata->fwd_Ea[81] * invT2;
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
    k_f = udata->prefactor_units[82] * udata->fwd_A[82]
                * exp(udata->fwd_beta[82] * tc[0] - udata->activation_units[82] * udata->fwd_Ea[82] * invT);
    dlnkfdT = udata->fwd_beta[82] * invT + udata->activation_units[82] * udata->fwd_Ea[82] * invT2;
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
    k_f = udata->prefactor_units[83] * udata->fwd_A[83]
                * exp(udata->fwd_beta[83] * tc[0] - udata->activation_units[83] * udata->fwd_Ea[83] * invT);
    dlnkfdT = udata->fwd_beta[83] * invT + udata->activation_units[83] * udata->fwd_Ea[83] * invT2;
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
