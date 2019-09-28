#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#include <Transport_F.H>
#include <onePt_F.H>
#include <PlotFileFromMF.H>
#include <mechanism.h>
#include <chemistry_file.H>
#include <actual_reactor.H>

std::string inputs_name = "";

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);
    {

      ParmParse pp;
    
      std::string probin_file = "probin";
      pp.query("probin_file",probin_file);
      int probin_file_length = probin_file.length();
      std::vector<int> probin_file_name(probin_file_length);

      for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

      extern_init(&(probin_file_name[0]),&probin_file_length);

#ifdef _OPENMP
#pragma omp parallel
#endif  
      extern_init_reactor();

    
      std::vector<int> npts(3,1);
      for (int i = 0; i < BL_SPACEDIM; ++i) {
	//npts[i] = 128;
      }
    
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

      std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
      for (int i=0; i<BL_SPACEDIM; ++i) {
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
      }
    
      int max_size = 32;
      pp.query("max_size",max_size);
      BoxArray ba(domain);
      ba.maxSize(max_size);

      int num_grow = 0;
      DistributionMapping dmap(ba);
      MultiFab mass_frac(ba,dmap,NUM_SPECIES,num_grow);
      MultiFab temperature(ba,dmap,1,num_grow);
      MultiFab internal_energy(ba,dmap,1,num_grow);

      IntVect tilesize(D_DECL(10240,8,32));

      Real T = 2360.71152280225;
      Array<Real,NUM_SPECIES> rhoY={ 1.385750517830831E-003,  3.439324932465339E-007,  2.939621541435051E-008,
                                     1.529589462521231E-006,  2.504950310424716E-004,  1.614391793151254E-005,
                                     5.618952763518058E-004,  7.904437268168153E-008,  1.964038191110560E-005,
                                     6.715516612026511E-004,  0.000000000000000E+000,  4.996348205941698E-012,
                                     4.490662043323236E-013,  0.000000000000000E+000,  0.000000000000000E+000,
                                     0.000000000000000E+000,  0.000000000000000E+000,  1.378418241889268E-012,
                                     0.000000000000000E+000,  0.000000000000000E+000,  4.617045057908065E-004};

      Array<Real,NUM_SPECIES> FrhoY = { -45.9715016276571      , 3.173207943712875E-002,
                                        -1.211735435128485E-003, -5.548817821158521E-003,  -43.7301905977032     ,
                                        -9.096045747696510E-002,  0.337605209731932     ,  8.135029186239005E-003,
                                        0.943763235459477      ,-0.646263505985702      ,-1.916363537409386E-006,
                                        1.698395916474832E-004 , 1.305803790052794E-005 , 1.444739232509387E-002,
                                        -9.07367952512635      ,-6.690623492013716E-005 , 0.222474510456600     ,
                                        2.828132632384628E-005 , 5.802097075854978E-004 , 2.915751802574431E-002,
                                        -15.3182145674259};
      Real Frhoh_cgs = -271871146560.329 * 10;

      //Frhoh_cgs = 0;
      //FrhoY = {0};
      
      Array<Real,NUM_SPECIES> FrhoY_cgs;
      for (int n=0; n<NUM_SPECIES; ++n) FrhoY_cgs[n] = FrhoY[n] * 1.e-3;

      Real dt = 3.304704155e-06;

      Real R=0;
      for (int n=0; n<NUM_SPECIES; ++n) R+=rhoY[n];

      Array<Real,NUM_SPECIES+1> vect;
      for (int n=0; n<NUM_SPECIES; ++n) vect[n] = rhoY[n] * 1.e-3;
      vect[NUM_SPECIES] = T;

      Array<Real,NUM_SPECIES> Y;
      for (int n=0; n<NUM_SPECIES; ++n) Y[n] = rhoY[n] / R;
      Real h_cgs;
      CKHBMS(&T,&(Y[0]),&h_cgs);
      
      Real R_cgs = R * 1.e-3;
      Real rhoh_cgs = h_cgs * R_cgs;

      double pressure = 1.0; // dummy FIXME
      double time_init = 0;
      int i=0, j=0, k=0;

      Real fc = react(&(vect[0]), &(FrhoY_cgs[0]), &rhoh_cgs, &Frhoh_cgs, &pressure, &dt, &time_init, &i, &j, &k);

      Print() << "fc = " << fc << " final T: " << vect[NUM_SPECIES] << std::endl;
      
      extern_close();

    }

    Finalize();

    return 0;
}
