#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"
#include <GPU_misc.H>

#include <main_F.H>
#include <PlotFileFromMF.H>

std::string inputs_name = "";

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    Real timer_tot = amrex::second();

    Real timer_init = 0.;
    Real timer_advance = 0.;
    Real timer_print = 0.;
    Real timer_print_tmp = 0.;

    {

      timer_init = amrex::second();

      ParmParse pp;
    
      std::string probin_file = "probin";
      pp.query("probin_file",probin_file);
      int probin_file_length = probin_file.length();
      std::vector<int> probin_file_name(probin_file_length);

      for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

      int fuel_idx = NC12H26_ID;
      int oxy_idx  = O2_ID;
      int bath_idx = N2_ID;

      extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx,&bath_idx);
    
      std::vector<int> npts(3,1);
      for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 8;
      }
      npts[1] = 32;
    
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

      std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
      for (int i=0; i<BL_SPACEDIM; ++i) {
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
      }
    
      int max_size = 8;
      pp.query("max_size",max_size);
      BoxArray ba(domain);
      ba.maxSize(max_size);

      int num_spec;
      num_spec = NUM_SPECIES;

      DistributionMapping dm{ba};

      int num_grow = 0;
      MultiFab mass_frac(ba,dm,num_spec,num_grow);
      MultiFab temperature(ba,dm,1,num_grow);
      MultiFab density(ba,dm,1,num_grow);

      IntVect tilesize(D_DECL(10240,8,32));

      BL_PROFILE_VAR("initialize_data()", InitData);
    
      int box_count =0;
      for (MFIter mfi(mass_frac,tilesize); mfi.isValid(); ++mfi) {
	const Box& box = mfi.tilebox();
	initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			BL_TO_FORTRAN_N_3D(temperature[mfi],0),
			BL_TO_FORTRAN_N_3D(density[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
	box_count +=1;
      }

      BL_PROFILE_VAR_STOP(InitData);

      printf("That many boxes: %d \n", box_count);

      timer_init = amrex::second() - timer_init;

      timer_print = amrex::second();
      ParmParse ppa("amr");
      std::string pltfile("plt");  
      ppa.query("plot_file",pltfile);
      std::string outfile = amrex::Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
      PlotFileFromMF(temperature,outfile);
      timer_print = amrex::second() - timer_print;

      BL_PROFILE_VAR("advance()", Advance);

      timer_advance = amrex::second();

      MultiFab wdots(ba,dm,num_spec,num_grow);
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

	const Box& box = mfi.tilebox();

	auto  mf      = mass_frac.array(mfi);
	auto  temp    = temperature.array(mfi);
	auto  rho     = density.array(mfi); 
	auto  cdots   = wdots.array(mfi);

	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_RTY2W(i, j, k, rho, temp, mf, cdots);
		// put the gunction here
	    });

      }

      BL_PROFILE_VAR_STOP(Advance);

      timer_advance = amrex::second() - timer_advance;

      timer_print_tmp = amrex::second();
      outfile = amrex::Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(wdots,outfile);
      timer_print = amrex::second() - timer_print_tmp + timer_print;

      extern_close();

    }

    timer_tot = amrex::second() - timer_tot;

    ParallelDescriptor::ReduceRealMax({timer_tot, timer_init, timer_advance, timer_print},
                                     ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run Time total        = " << timer_tot     << "\n"
                   << "Run Time init         = " << timer_init    << "\n"
                   << "Run Time advance      = " << timer_advance << "\n"
                   << "Run Time print plt    = " << timer_print << "\n";

    BL_PROFILE_VAR_STOP(pmain);

    Finalize();

    return 0;
}
