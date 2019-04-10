--- Options are :

- cvode_iJac
# 1 Using Analytical J
# 0 Using FD J

- cvode_iDense
# 1  = Using a Direct Dense solver
# 99 = Using a GMRES Iterative solver
        !! No Jac so cvode_iJac here differentiate between
	   No preconditioning (0), Left preconditioned (1)
        !! When Left preconditioned (1) + KLU ENABLED then
           The Preconditioning is sparse and optimized
# 5  = ( IF KLU ENABLED ) Using a Sparse Direct solver
        !! THis one can only be used with cvode_iJac = 1


--- RECAP Posssible choices of inputs:
        cvode_iDense     cvode_iJac 
        ---------------------------
             Direct Dense solver
        ---------------------------
             1               0
             1               1
        ---------------------------
             Direct Sparse solver
              ( KLU ENABLED )
        ---------------------------
             5               1
        ---------------------------
             GMRES Iterative solver
        ---------------------------
             99              0
             99              1
        ---------------------------
             GMRES Iterative solver
              ( KLU ENABLED )
        ---------------------------
             99              1
        
