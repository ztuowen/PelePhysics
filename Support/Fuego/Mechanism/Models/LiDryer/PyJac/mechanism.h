#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 8
/* Species Indexes
0  H2
1  O2
2  H2O
3  H
4  O
5  OH
6  HO2
7  H2O2
8  N2
*/

//Number of species
#define NSP 9
//Number of variables. NN = NSP + 1 (temperature)
#define NN 10
//Number of forward reactions
#define FWD_RATES 21
//Number of reversible reactions
#define REV_RATES 21
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 6

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

