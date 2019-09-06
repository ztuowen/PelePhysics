#ifndef _kernel_H_
#define _kernel_H_
#include <AMReX_GpuDevice.H>
#include "chemistry_file.H"

#define NUM_SPECIES 21
#define NUM_REACTIONS 84

using namespace amrex;

AMREX_GPU_DEVICE 
void 
W_spec_d(int i, int j, int k,
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& Y,
		amrex::Array4<amrex::Real> const& cdots)
{

	amrex::GpuArray<amrex::Real,NUM_SPECIES> sc;
	amrex::GpuArray<amrex::Real,NUM_REACTIONS> qf;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> wdot;
	amrex::Real imw[NUM_SPECIES];
	amrex::Real invT, logT, k_f, qdot, Corr;

	get_imw(imw);

	invT = 1/temp(i,j,k);
	logT = log(temp(i,j,k));

	for (int n = 0; n < NUM_SPECIES; n++) {
	        sc[n] = Y(i,j,k,n) * rho(i,j,k) * imw[n];
	}

        /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
        qf[0] = sc[1]*sc[7];
        /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
        qf[1] = sc[1]*sc[9];
        /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
        qf[2] = sc[1]*sc[13];
        /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
        qf[3] = sc[1]*sc[14];
        /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
        qf[4] = sc[1]*sc[16];
        /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
        qf[5] = sc[1]*sc[17];
        /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
        qf[6] = sc[0]*sc[11];
        /*reaction 8: 2 CH3 (+M) <=> C2H6 (+M) */
        qf[7] = sc[9]*sc[9];
        /*reaction 9: O + H + M <=> OH + M */
        qf[8] = sc[1]*sc[2];
        /*reaction 10: O + CO + M <=> CO2 + M */
        qf[9] = sc[2]*sc[11];
        /*reaction 11: H + O2 + M <=> HO2 + M */
        qf[10] = sc[1]*sc[3];
        /*reaction 12: 2 H + M <=> H2 + M */
        qf[11] = sc[1]*sc[1];
        /*reaction 13: H + OH + M <=> H2O + M */
        qf[12] = sc[1]*sc[4];
        /*reaction 14: HCO + M <=> H + CO + M */
        qf[13] = sc[13];
        /*reaction 15: O + H2 <=> H + OH */
        qf[14] = sc[0]*sc[2];
        /*reaction 16: O + HO2 <=> OH + O2 */
        qf[15] = sc[2]*sc[6];
        /*reaction 17: O + CH2 <=> H + HCO */
        qf[16] = sc[2]*sc[7];
        /*reaction 18: O + CH2(S) <=> H + HCO */
        qf[17] = sc[2]*sc[8];
        /*reaction 19: O + CH3 <=> H + CH2O */
        qf[18] = sc[2]*sc[9];
        /*reaction 20: O + CH4 <=> OH + CH3 */
        qf[19] = sc[2]*sc[10];
        /*reaction 21: O + HCO <=> OH + CO */
        qf[20] = sc[2]*sc[13];
        /*reaction 22: O + HCO <=> H + CO2 */
        qf[21] = sc[2]*sc[13];
        /*reaction 23: O + CH2O <=> OH + HCO */
        qf[22] = sc[2]*sc[14];
        /*reaction 24: O + C2H4 <=> CH3 + HCO */
        qf[23] = sc[2]*sc[16];
        /*reaction 25: O + C2H5 <=> CH3 + CH2O */
        qf[24] = sc[2]*sc[17];
        /*reaction 26: O + C2H6 <=> OH + C2H5 */
        qf[25] = sc[2]*sc[18];
        /*reaction 27: O2 + CO <=> O + CO2 */
        qf[26] = sc[3]*sc[11];
        /*reaction 28: O2 + CH2O <=> HO2 + HCO */
        qf[27] = sc[3]*sc[14];
        /*reaction 29: H + 2 O2 <=> HO2 + O2 */
        qf[28] = sc[1]*sc[3]*sc[3];
        /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
        qf[29] = sc[1]*sc[3]*sc[5];
        /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
        qf[30] = sc[1]*sc[3]*sc[19];
        /*reaction 32: H + O2 + AR <=> HO2 + AR */
        qf[31] = sc[1]*sc[3]*sc[20];
        /*reaction 33: H + O2 <=> O + OH */
        qf[32] = sc[1]*sc[3];
        /*reaction 34: 2 H + H2 <=> 2 H2 */
        qf[33] = sc[0]*sc[1]*sc[1];
        /*reaction 35: 2 H + H2O <=> H2 + H2O */
        qf[34] = sc[1]*sc[1]*sc[5];
        /*reaction 36: 2 H + CO2 <=> H2 + CO2 */
        qf[35] = sc[1]*sc[1]*sc[12];
        /*reaction 37: H + HO2 <=> O2 + H2 */
        qf[36] = sc[1]*sc[6];
        /*reaction 38: H + HO2 <=> 2 OH */
        qf[37] = sc[1]*sc[6];
        /*reaction 39: H + CH4 <=> CH3 + H2 */
        qf[38] = sc[1]*sc[10];
        /*reaction 40: H + HCO <=> H2 + CO */
        qf[39] = sc[1]*sc[13];
        /*reaction 41: H + CH2O <=> HCO + H2 */
        qf[40] = sc[1]*sc[14];
        /*reaction 42: H + CH3O <=> OH + CH3 */
        qf[41] = sc[1]*sc[15];
        /*reaction 43: H + C2H6 <=> C2H5 + H2 */
        qf[42] = sc[1]*sc[18];
        /*reaction 44: OH + H2 <=> H + H2O */
        qf[43] = sc[0]*sc[4];
        /*reaction 45: 2 OH <=> O + H2O */
        qf[44] = sc[4]*sc[4];
        /*reaction 46: OH + HO2 <=> O2 + H2O */
        qf[45] = sc[4]*sc[6];
        /*reaction 47: OH + CH2 <=> H + CH2O */
        qf[46] = sc[4]*sc[7];
        /*reaction 48: OH + CH2(S) <=> H + CH2O */
        qf[47] = sc[4]*sc[8];
        /*reaction 49: OH + CH3 <=> CH2 + H2O */
        qf[48] = sc[4]*sc[9];
        /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
        qf[49] = sc[4]*sc[9];
        /*reaction 51: OH + CH4 <=> CH3 + H2O */
        qf[50] = sc[4]*sc[10];
        /*reaction 52: OH + CO <=> H + CO2 */
        qf[51] = sc[4]*sc[11];
        /*reaction 53: OH + HCO <=> H2O + CO */
        qf[52] = sc[4]*sc[13];
        /*reaction 54: OH + CH2O <=> HCO + H2O */
        qf[53] = sc[4]*sc[14];
        /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
        qf[54] = sc[4]*sc[18];
        /*reaction 56: HO2 + CH2 <=> OH + CH2O */
        qf[55] = sc[6]*sc[7];
        /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
        qf[56] = sc[6]*sc[9];
        /*reaction 58: HO2 + CH3 <=> OH + CH3O */
        qf[57] = sc[6]*sc[9];
        /*reaction 59: HO2 + CO <=> OH + CO2 */
        qf[58] = sc[6]*sc[11];
        /*reaction 60: CH2 + O2 <=> OH + HCO */
        qf[59] = sc[3]*sc[7];
        /*reaction 61: CH2 + H2 <=> H + CH3 */
        qf[60] = sc[0]*sc[7];
        /*reaction 62: CH2 + CH3 <=> H + C2H4 */
        qf[61] = sc[7]*sc[9];
        /*reaction 63: CH2 + CH4 <=> 2 CH3 */
        qf[62] = sc[7]*sc[10];
        /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
        qf[63] = sc[8]*sc[19];
        /*reaction 65: CH2(S) + AR <=> CH2 + AR */
        qf[64] = sc[8]*sc[20];
        /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
        qf[65] = sc[3]*sc[8];
        /*reaction 67: CH2(S) + O2 <=> CO + H2O */
        qf[66] = sc[3]*sc[8];
        /*reaction 68: CH2(S) + H2 <=> CH3 + H */
        qf[67] = sc[0]*sc[8];
        /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
        qf[68] = sc[5]*sc[8];
        /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
        qf[69] = sc[8]*sc[9];
        /*reaction 71: CH2(S) + CH4 <=> 2 CH3 */
        qf[70] = sc[8]*sc[10];
        /*reaction 72: CH2(S) + CO <=> CH2 + CO */
        qf[71] = sc[8]*sc[11];
        /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
        qf[72] = sc[8]*sc[12];
        /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
        qf[73] = sc[8]*sc[12];
        /*reaction 75: CH3 + O2 <=> O + CH3O */
        qf[74] = sc[3]*sc[9];
        /*reaction 76: CH3 + O2 <=> OH + CH2O */
        qf[75] = sc[3]*sc[9];
        /*reaction 77: 2 CH3 <=> H + C2H5 */
        qf[76] = sc[9]*sc[9];
        /*reaction 78: CH3 + HCO <=> CH4 + CO */
        qf[77] = sc[9]*sc[13];
        /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
        qf[78] = sc[9]*sc[14];
        /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
        qf[79] = sc[9]*sc[18];
        /*reaction 81: HCO + H2O <=> H + CO + H2O */
        qf[80] = sc[5]*sc[13];
        /*reaction 82: HCO + O2 <=> HO2 + CO */
        qf[81] = sc[3]*sc[13];
        /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
        qf[82] = sc[3]*sc[15];
        /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
        qf[83] = sc[3]*sc[17];

        // (0):  O + H + M <=> OH + M
        qf[8] *= 1.0000000000000002e-12 * 5e+17 
                   * exp(-1 * logT - 0.50321666580471969 * 0 * invT);
        // (1):  O + H2 <=> H + OH
        qf[14] *= 1.0000000000000002e-06 * 50000 
               * exp(2.6699999999999999 * logT - 0.50321666580471969 * 6290 * invT);
        // (2):  O + HO2 <=> OH + O2
        qf[15] *= 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (3):  O + CH2 <=> H + HCO
        qf[16] *= 1.0000000000000002e-06 * 80000000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (4):  O + CH2(S) <=> H + HCO
        qf[17] *= 1.0000000000000002e-06 * 15000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (5):  O + CH3 <=> H + CH2O
        qf[18] *= 1.0000000000000002e-06 * 84300000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (6):  O + CH4 <=> OH + CH3
        qf[19] *= 1.0000000000000002e-06 * 1020000000 
                   * exp(1.5 * logT - 0.50321666580471969 * 8600 * invT);
        // (7):  O + CO + M <=> CO2 + M
        qf[9] *= 1.0000000000000002e-12 * 602000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 3000 * invT);
        // (8):  O + HCO <=> OH + CO
        qf[20] *= 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (9):  O + HCO <=> H + CO2
        qf[21] *= 1.0000000000000002e-06 * 30000000000000
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (10):  O + CH2O <=> OH + HCO
        qf[22] *= 1.0000000000000002e-06 * 39000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 3540 * invT);
        // (11):  O + C2H4 <=> CH3 + HCO
        qf[23] *= 1.0000000000000002e-06 * 19200000 
                   * exp(1.8300000000000001 * logT - 0.50321666580471969 * 220 * invT);
        // (12):  O + C2H5 <=> CH3 + CH2O
        qf[24] *= 1.0000000000000002e-06 * 132000000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (13):  O + C2H6 <=> OH + C2H5
        qf[25] *= 1.0000000000000002e-06 * 89800000 
                   * exp(1.9199999999999999 * logT - 0.50321666580471969 * 5690 * invT);
        // (14):  O2 + CO <=> O + CO2
        qf[26] *= 1.0000000000000002e-06 * 2500000000000 
                   * exp(0 * logT - 0.50321666580471969 * 47800 * invT);
        // (15):  O2 + CH2O <=> HO2 + HCO
        qf[27] *= 1.0000000000000002e-06 * 100000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 40000 * invT);
        // (16):  H + O2 + M <=> HO2 + M
        qf[10] *= 1.0000000000000002e-12 * 2.8e+18 
                   * exp(-0.85999999999999999 * logT - 0.50321666580471969 * 0 * invT);
        // (17):  H + 2 O2 <=> HO2 + O2
        qf[28] *= 1.0000000000000002e-12 * 3e+20 
                   * exp(-1.72 * logT - 0.50321666580471969 * 0 * invT);
        // (18):  H + O2 + H2O <=> HO2 + H2O
        qf[29] *= 1.0000000000000002e-12 * 9.38e+18 
                   * exp(-0.76000000000000001 * logT - 0.50321666580471969 * 0 * invT);
        // (19):  H + O2 + N2 <=> HO2 + N2
	qf[30] *= 1.0000000000000002e-12 * 3.75e+20
                   * exp(-1.72 * logT - 0.50321666580471969 * 0 * invT);
        // (20):  H + O2 + AR <=> HO2 + AR
        qf[31] *= 1.0000000000000002e-12 * 7e+17 
                   * exp(-0.80000000000000004 * logT - 0.50321666580471969 * 0 * invT);
        // (21):  H + O2 <=> O + OH
        qf[32] *= 1.0000000000000002e-06 * 83000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 14413 * invT);
        // (22):  2 H + M <=> H2 + M
        qf[11] *= 1.0000000000000002e-12 * 1e+18 
                   * exp(-1 * logT - 0.50321666580471969 * 0 * invT);
        // (23):  2 H + H2 <=> 2 H2
        qf[33] *= 1.0000000000000002e-12 * 90000000000000000 
                   * exp(-0.59999999999999998 * logT - 0.50321666580471969 * 0 * invT);
        // (24):  2 H + H2O <=> H2 + H2O
        qf[34] *= 1.0000000000000002e-12 * 6e+19 
                   * exp(-1.25 * logT - 0.50321666580471969 * 0 * invT);
        // (25):  2 H + CO2 <=> H2 + CO2
        qf[35] *= 1.0000000000000002e-12 * 5.5e+20 
                   * exp(-2 * logT - 0.50321666580471969 * 0 * invT);
        // (26):  H + OH + M <=> H2O + M
        qf[12] *= 1.0000000000000002e-12 * 2.2e+22 
                   * exp(-2 * logT - 0.50321666580471969 * 0 * invT);
        // (27):  H + HO2 <=> O2 + H2
        qf[36] *= 1.0000000000000002e-06 * 28000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 1068 * invT);
        // (28):  H + HO2 <=> 2 OH
        qf[37] *= 1.0000000000000002e-06 * 134000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 635 * invT);
        // (29):  H + CH2 (+M) <=> CH3 (+M)
        qf[0] *= 1.0000000000000002e-06 * 25000000000000000 
                   * exp(-0.80000000000000004 * logT - 0.50321666580471969 * 0 * invT);
        // (30):  H + CH3 (+M) <=> CH4 (+M)
        qf[1] *= 1.0000000000000002e-06 * 12700000000000000 
                   * exp(-0.63 * logT - 0.50321666580471969 * 383 * invT);
        // (31):  H + CH4 <=> CH3 + H2
        qf[38] *= 1.0000000000000002e-06 * 660000000 
                   * exp(1.6200000000000001 * logT - 0.50321666580471969 * 10840 * invT);
        // (32):  H + HCO (+M) <=> CH2O (+M)
        qf[2] *= 1.0000000000000002e-06 * 1090000000000 
                   * exp(0.47999999999999998 * logT - 0.50321666580471969 * -260 * invT);
        // (33):  H + HCO <=> H2 + CO
        qf[39] *= 1.0000000000000002e-06 * 73400000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (34):  H + CH2O (+M) <=> CH3O (+M)
        qf[3] *= 1.0000000000000002e-06 * 540000000000 
                   * exp(0.45400000000000001 * logT - 0.50321666580471969 * 2600 * invT);
        // (35):  H + CH2O <=> HCO + H2
        qf[40] *= 1.0000000000000002e-06 * 23000000000 
                   * exp(1.05 * logT - 0.50321666580471969 * 3275 * invT);
        // (36):  H + CH3O <=> OH + CH3
        qf[41] *= 1.0000000000000002e-06 * 32000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (37):  H + C2H4 (+M) <=> C2H5 (+M)
        qf[4] *= 1.0000000000000002e-06 * 1080000000000 
                   * exp(0.45400000000000001 * logT - 0.50321666580471969 * 1820 * invT);
        // (38):  H + C2H5 (+M) <=> C2H6 (+M)
        qf[5] *= 1.0000000000000002e-06 * 5.21e+17 
                   * exp(-0.98999999999999999 * logT - 0.50321666580471969 * 1580 * invT);
        // (39):  H + C2H6 <=> C2H5 + H2
        qf[42] *= 1.0000000000000002e-06 * 115000000 
                   * exp(1.8999999999999999 * logT - 0.50321666580471969 * 7530 * invT);
        // (40):  H2 + CO (+M) <=> CH2O (+M)
        qf[6] *= 1.0000000000000002e-06 * 43000000 
                   * exp(1.5 * logT - 0.50321666580471969 * 79600 * invT);
        // (41):  OH + H2 <=> H + H2O
        qf[43] *= 1.0000000000000002e-06 * 216000000 
                   * exp(1.51 * logT - 0.50321666580471969 * 3430 * invT);
        // (42):  2 OH <=> O + H2O
        qf[44] *= 1.0000000000000002e-06 * 35700 
                   * exp(2.3999999999999999 * logT - 0.50321666580471969 * -2110 * invT);
        // (43):  OH + HO2 <=> O2 + H2O
        qf[45] *= 1.0000000000000002e-06 * 29000000000000 
                   * exp(0 * logT - 0.50321666580471969 * -500 * invT);
        // (44):  OH + CH2 <=> H + CH2O
        k_f = 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[46] *= Corr * k_f;
        // (45):  OH + CH2(S) <=> H + CH2O
        k_f = 1.0000000000000002e-06 * 30000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[47] *= Corr * k_f;
        // (46):  OH + CH3 <=> CH2 + H2O
        k_f = 1.0000000000000002e-06 * 56000000 
                   * exp(1.6000000000000001 * logT - 0.50321666580471969 * 5420 * invT);
        Corr  = 1.0;
        qf[48] *= Corr * k_f;
        // (47):  OH + CH3 <=> CH2(S) + H2O
        k_f = 1.0000000000000002e-06 * 25010000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[49] *= Corr * k_f;
        // (48):  OH + CH4 <=> CH3 + H2O
        k_f = 1.0000000000000002e-06 * 100000000 
                   * exp(1.6000000000000001 * logT - 0.50321666580471969 * 3120 * invT);
        Corr  = 1.0;
        qf[50] *= Corr * k_f;
        // (49):  OH + CO <=> H + CO2
        k_f = 1.0000000000000002e-06 * 47600000 
                   * exp(1.228 * logT - 0.50321666580471969 * 70 * invT);
        Corr  = 1.0;
        qf[51] *= Corr * k_f;
        // (50):  OH + HCO <=> H2O + CO
        k_f = 1.0000000000000002e-06 * 50000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[52] *= Corr * k_f;
        // (51):  OH + CH2O <=> HCO + H2O
        k_f = 1.0000000000000002e-06 * 3430000000 
                   * exp(1.1799999999999999 * logT - 0.50321666580471969 * -447 * invT);
        Corr  = 1.0;
        qf[53] *= Corr * k_f;
        // (52):  OH + C2H6 <=> C2H5 + H2O
        k_f = 1.0000000000000002e-06 * 3540000 
                   * exp(2.1200000000000001 * logT - 0.50321666580471969 * 870 * invT);
        Corr  = 1.0;
        qf[54] *= Corr * k_f;
        // (53):  HO2 + CH2 <=> OH + CH2O
        k_f = 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[55] *= Corr * k_f;
        // (54):  HO2 + CH3 <=> O2 + CH4
        k_f = 1.0000000000000002e-06 * 1000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[56] *= Corr * k_f;
        // (55):  HO2 + CH3 <=> OH + CH3O
        k_f = 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[57] *= Corr * k_f;
        // (56):  HO2 + CO <=> OH + CO2
        k_f = 1.0000000000000002e-06 * 150000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 23600 * invT);
        Corr  = 1.0;
        qf[58] *= Corr * k_f;
        // (57):  CH2 + O2 <=> OH + HCO
        k_f = 1.0000000000000002e-06 * 13200000000000 
                   * exp(0 * logT - 0.50321666580471969 * 1500 * invT);
        Corr  = 1.0;
        qf[59] *= Corr * k_f;
        // (58):  CH2 + H2 <=> H + CH3
        k_f = 1.0000000000000002e-06 * 500000 
                   * exp(2 * logT - 0.50321666580471969 * 7230 * invT);
        Corr  = 1.0;
        qf[60] *= Corr * k_f;
        // (59):  CH2 + CH3 <=> H + C2H4
        k_f = 1.0000000000000002e-06 * 40000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[61] *= Corr * k_f;
        // (60):  CH2 + CH4 <=> 2 CH3
        k_f = 1.0000000000000002e-06 * 2460000 
                   * exp(2 * logT - 0.50321666580471969 * 8270 * invT);
        Corr  = 1.0;
        qf[62] *= Corr * k_f;
        // (61):  CH2(S) + N2 <=> CH2 + N2
        k_f = 1.0000000000000002e-06 * 15000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 600 * invT);
        Corr  = 1.0;
        qf[63] *= Corr * k_f;
        // (62):  CH2(S) + AR <=> CH2 + AR
        k_f = 1.0000000000000002e-06 * 9000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 600 * invT);
        Corr  = 1.0;
        qf[64] *= Corr * k_f;
        // (63):  CH2(S) + O2 <=> H + OH + CO
        k_f = 1.0000000000000002e-06 * 28000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[65] *= Corr * k_f;
        // (64):  CH2(S) + O2 <=> CO + H2O
        k_f = 1.0000000000000002e-06 * 12000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[66] *= Corr * k_f;
        // (65):  CH2(S) + H2 <=> CH3 + H
        k_f = 1.0000000000000002e-06 * 70000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[67] *= Corr * k_f;
        // (66):  CH2(S) + H2O <=> CH2 + H2O
        k_f = 1.0000000000000002e-06 * 30000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[68] *= Corr * k_f;
        // (67):  CH2(S) + CH3 <=> H + C2H4
        k_f = 1.0000000000000002e-06 * 12000000000000 
                   * exp(0 * logT - 0.50321666580471969 * -570 * invT);
        Corr  = 1.0;
        qf[69] *= Corr * k_f;
        // (68):  CH2(S) + CH4 <=> 2 CH3
        k_f = 1.0000000000000002e-06 * 16000000000000 
                   * exp(0 * logT - 0.50321666580471969 * -570 * invT);
        Corr  = 1.0;
        qf[70] *= Corr * k_f;
        // (69):  CH2(S) + CO <=> CH2 + CO
        k_f = 1.0000000000000002e-06 * 9000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[71] *= Corr * k_f;
        // (70):  CH2(S) + CO2 <=> CH2 + CO2
        k_f = 1.0000000000000002e-06 * 7000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[72] *= Corr * k_f;
        // (71):  CH2(S) + CO2 <=> CO + CH2O
        k_f = 1.0000000000000002e-06 * 14000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[73] *= Corr * k_f;
        // (72):  CH3 + O2 <=> O + CH3O
        k_f = 1.0000000000000002e-06 * 26750000000000 
                   * exp(0 * logT - 0.50321666580471969 * 28800 * invT);
        Corr  = 1.0;
        qf[74] *= Corr * k_f;
        // (73):  CH3 + O2 <=> OH + CH2O
        qf[75] *= 1.0000000000000002e-06 * 36000000000 
                   * exp(0 * logT - 0.50321666580471969 * 8940 * invT);
        // (74):  2 CH3 (+M) <=> C2H6 (+M)
        qf[7] *= 1.0000000000000002e-06 * 21200000000000000 
                   * exp(-0.96999999999999997 * logT - 0.50321666580471969 * 620 * invT);
        // (75):  2 CH3 <=> H + C2H5
        qf[76] *= 1.0000000000000002e-06 * 4990000000000 
                   * exp(0.10000000000000001 * logT - 0.50321666580471969 * 10600 * invT);
        // (76):  CH3 + HCO <=> CH4 + CO
        qf[77] *= 1.0000000000000002e-06 * 26480000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (77):  CH3 + CH2O <=> HCO + CH4
        qf[78] *= 1.0000000000000002e-06 * 3320 
                   * exp(2.8100000000000001 * logT - 0.50321666580471969 * 5860 * invT);
        // (78):  CH3 + C2H6 <=> C2H5 + CH4
        qf[79] *= 1.0000000000000002e-06 * 6140000 
                   * exp(1.74 * logT - 0.50321666580471969 * 10450 * invT);
        // (79):  HCO + H2O <=> H + CO + H2O
        qf[80] *= 1.0000000000000002e-06 * 2.244e+18 
                   * exp(-1 * logT - 0.50321666580471969 * 17000 * invT);
        // (80):  HCO + M <=> H + CO + M
        qf[13] *= 1.0000000000000002e-06 * 1.87e+17 
                   * exp(-1 * logT - 0.50321666580471969 * 17000 * invT);
        // (81):  HCO + O2 <=> HO2 + CO
        qf[81] *= 1.0000000000000002e-06 * 7600000000000 
                   * exp(0 * logT - 0.50321666580471969 * 400 * invT);
        // (82):  CH3O + O2 <=> HO2 + CH2O
        qf[82] *= 1.0000000000000002e-06 * 4.2799999999999999e-13 
                   * exp(7.5999999999999996 * logT - 0.50321666580471969 * -3530 * invT);
        // (83):  C2H5 + O2 <=> HO2 + C2H4
        qf[83] *= 1.0000000000000002e-06 * 840000000000 
                   * exp(0 * logT - 0.50321666580471969 * 3875 * invT);

        for (int i = 0; i < 21; ++i) {
            wdot[i] = 0.0;
        }

        qdot = qf[0];
        wdot[1] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[1];
        wdot[1] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[2];
        wdot[1]  -= qdot;
        wdot[13] -= qdot;

        qdot = qf[3];
        wdot[1]  -= qdot;
        wdot[14] -= qdot;

        qdot = qf[4];
        wdot[1]  -= qdot;
        wdot[16] -= qdot;

        qdot = qf[5];
        wdot[1] -= qdot;
        wdot[17] -= qdot;

        qdot = qf[6];
        wdot[0] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[7];
        wdot[9] -= 2 * qdot;

        qdot = qf[8];
        wdot[1] -= qdot;
        wdot[2] -= qdot;

        qdot = qf[9];
        wdot[2] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[10];
        wdot[1] -= qdot;
        wdot[3] -= qdot;

        qdot = qf[11];
        wdot[1] -= 2 * qdot;

        qdot = qf[12];
        wdot[1] -= qdot;
        wdot[4] -= qdot;

        qdot = qf[13];
        wdot[13] -= qdot;

        qdot = qf[14];
        wdot[0] -= qdot;
        wdot[2] -= qdot;

        qdot = qf[15];
        wdot[2] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[16];
        wdot[2] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[17];
        wdot[2] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[18];
        wdot[2] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[19];
        wdot[2] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[20];
        wdot[2] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[21];
        wdot[2] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[22];
        wdot[2] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[23];
        wdot[2] -= qdot;
        wdot[16] -= qdot;

        qdot = qf[24];
        wdot[2] -= qdot;
        wdot[17] -= qdot;

        qdot = qf[25];
        wdot[2] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[26];
        wdot[3] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[27];
        wdot[3] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[28];
        wdot[1] -= qdot;
        wdot[3] -= 2 * qdot;

        qdot = qf[29];
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[5] -= qdot;

        qdot = qf[30];
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[19] -= qdot;

        qdot = qf[31];
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[20] -= qdot;

        qdot = qf[32];
        wdot[1] -= qdot;
        wdot[3] -= qdot;

        qdot = qf[33];
        wdot[0] -= qdot;
        wdot[1] -= 2 * qdot;

        qdot = qf[34];
        wdot[1] -= 2 * qdot;
        wdot[5] -= qdot;

        qdot = qf[35];
        wdot[1] -= 2 * qdot;
        wdot[12] -= qdot;

        qdot = qf[36];
        wdot[1] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[37];
        wdot[1] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[38];
        wdot[1] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[39];
        wdot[1] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[40];
        wdot[1] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[41];
        wdot[1] -= qdot;
        wdot[15] -= qdot;

        qdot = qf[42];
        wdot[1] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[43];
        wdot[0] -= qdot;
        wdot[4] -= qdot;

        qdot = qf[44];
        wdot[4] -= 2 * qdot;

        qdot = qf[45];
        wdot[4] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[46];
        wdot[4] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[47];
        wdot[4] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[48];
        wdot[4] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[49];
        wdot[4] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[50];
        wdot[4] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[51];
        wdot[4] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[52];
        wdot[4] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[53];
        wdot[4] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[54];
        wdot[4] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[55];
        wdot[6] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[56];
        wdot[6] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[57];
        wdot[6] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[58];
        wdot[6] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[59];
        wdot[3] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[60];
        wdot[0] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[61];
        wdot[7] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[62];
        wdot[7] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[63];
        wdot[8] -= qdot;
        wdot[19] -= qdot;

        qdot = qf[64];
        wdot[8] -= qdot;
        wdot[20] -= qdot;

        qdot = qf[65];
        wdot[3] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[66];
        wdot[3] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[67];
        wdot[0] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[68];
        wdot[5] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[69];
        wdot[8] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[70];
        wdot[8] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[71];
        wdot[8] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[72];
        wdot[8] -= qdot;
        wdot[12] -= qdot;

        qdot = qf[73];
        wdot[8] -= qdot;
        wdot[12] -= qdot;

        qdot = qf[74];
        wdot[3] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[75];
        wdot[3] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[76];
        wdot[9] -= 2 * qdot;

        qdot = qf[77];
        wdot[9] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[78];
        wdot[9] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[79];
        wdot[9] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[80];
        wdot[5] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[81];
        wdot[3] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[82];
        wdot[3] -= qdot;
        wdot[15] -= qdot;

        qdot = qf[83];
        wdot[3] -= qdot;
        wdot[17] -= qdot;

	for (int n = 0; n < NUM_SPECIES; n++) {
		cdots(i,j,k,n) = -wdot[n];
	}
}

#endif
