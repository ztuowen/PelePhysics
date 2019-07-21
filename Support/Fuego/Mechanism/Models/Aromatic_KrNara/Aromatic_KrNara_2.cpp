#include "chemistry_file.H"

/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(double* COFETA) {
    for (int i = 0 ; i < 631 ; i++) {
	    COFETA[i] =0.0;
    }
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(double* COFLAM) {
    for (int i = 0 ; i < 631 ; i++) {
	    COFLAM[i] =0.0;
    }
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 1;
    KTDIF[1] = 5;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(double* COFTD) {
    for (int i = 0 ; i < 1263 ; i++) {
	    COFTD[i] =0.0;
    }
}

