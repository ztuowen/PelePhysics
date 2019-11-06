#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <AMReX_GpuDevice.H>
#include <kernel.H>

template <typename L>
void ForMarc (Box const& box, int nc, L f) noexcept
{
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
#ifdef AMREX_USE_GPU
    Gpu::ExecutionConfig ec;
    ec.numBlocks.x = 2560;
    ec.numBlocks.x = 5120;
    ec.numBlocks.y = 1;
    ec.numBlocks.z = 1;
    ec.numThreads.x = nc;
    ec.numThreads.y = 1;
    ec.numThreads.z = 1;
    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept {
      int block_idx_x = blockIdx.x;
      int grid_dim_x = gridDim.x;
      int thread_idx_x = threadIdx.x;
#else
      int block_idx_x = 0;
      int grid_dim_x = 1;
      int thread_idx_x = 0;
#endif
      for (int icell = block_idx_x, stride = grid_dim_x; icell < ncells; icell += stride)
      {
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;
        int n = thread_idx_x;
        f(i,j,k,n);
      }
#ifdef AMREX_USE_GPU
    });
    AMREX_GPU_ERROR_CHECK();
#endif
}

void checker(int idx,
             int block_dim_x,
             int block_dim_y,
             int block_dim_z)
{  

  int num_spec = 21;
  int num_reac = 84;
  int num_threads = block_dim_x * block_dim_y * block_dim_z;
  
  {
    int npass = (num_spec + num_threads - 1) / num_threads;
    for (int pass = 0; pass<npass; ++pass) {
      int sid = pass*num_threads + idx;
      if (sid < num_spec) {
        Print() << "C  idx, sid: " << idx << " ," << sid << std::endl;
      }
    }
  }

  {
    int npass = (num_reac + num_threads - 1) / num_threads;

    for (int pass = 0; pass<npass; ++pass) {
      int rid = pass*num_threads + idx;
      if (rid < num_reac) {
        Print() << "Q  idx, rid: " << idx << " ," << rid << std::endl;
      }
    }
  }

  {
    int npass = (num_spec + num_threads - 1) / num_threads;
    for (int pass = 0; pass<npass; ++pass) {
      int sid = pass*num_threads + idx;

      if (sid<num_spec) {
        Print() << "W  idx, rid: " << idx << " ," << sid << std::endl;
      }
    }
  }
}



int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);
    {
      ParmParse pp;

      int nc=256;
      pp.query("nc",nc);
      Vector<int> n_cells(BL_SPACEDIM,nc);
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(n_cells[0]-1,n_cells[1]-1,n_cells[2]-1)));

      int max_size = 64;
      pp.query("max_size",max_size);
      BoxArray ba(domain);
      ba.maxSize(max_size);

      int num_spec = NUM_SPECIES;
      int num_reac = NUM_REACTIONS;

      DistributionMapping dm(ba);

      int num_grow = 0;
      MultiFab mass_frac(ba,dm,num_spec,num_grow);
      MultiFab temperature(ba,dm,1,num_grow);
      MultiFab density(ba,dm,1,num_grow);

      //checker(0,32,1,1);
#if 1

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("INIT");
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          const auto& temp = temperature.array(mfi);
          const auto& Y    = mass_frac.array(mfi);
          const auto& rho  = density.array(mfi);
          For(bx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
            for (int n=0; n<num_spec; ++n) {
              Y(i,j,k,n) = 1./num_spec;
            }
            temp(i,j,k) = 450;
            rho(i,j,k) = 0.75;
          });
        }
      }

      MultiFab wdot(ba,dm,num_spec,num_grow);
      wdot.setVal(0);

      kinit();

      #ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("PREFETCH");
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          temperature.prefetchToDevice(mfi);
          density.prefetchToDevice(mfi);
          mass_frac.prefetchToDevice(mfi);
        }
      }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("COMPUTE_W");
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& box = mfi.tilebox();
          const auto& temp = temperature.array(mfi);
          const auto& Y    = mass_frac.array(mfi);
          const auto& rho  = density.array(mfi);
          const auto& w    = wdot.array(mfi);
          int numPts = box.numPts();

          ForMarc(box, 84,
          [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
          {
            Real wtmp[21];
            W_spec_d(rho(i,j,k),temp(i,j,k),&(Y(i,j,k,0)),numPts,wtmp);

#ifdef AMREX_USE_GPU
            int idx = threadIdx.x;
            int block_dim_x = blockDim.x;
            int block_dim_y = blockDim.y;
            int block_dim_z = blockDim.z;
#else
            int idx = 0;
            int block_dim_x = 1;
            int block_dim_y = 1;
            int block_dim_z = 1;
#endif
            int num_threads = block_dim_x * block_dim_y * block_dim_z;
            for (int sid = idx; sid < num_spec; sid+=num_threads) {
              w(i,j,k,sid) = wtmp[sid];
            }
          });
        }
      }

      VisMF::Write(wdot,"WDOT");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("COMPUTE_W1");
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& box = mfi.tilebox();
          const auto& temp = temperature.array(mfi);
          const auto& Y    = mass_frac.array(mfi);
          const auto& rho  = density.array(mfi);
          const auto& w    = wdot.array(mfi);
          int numPts = box.numPts();

          For(box,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
            Real wtmp[21];
            W_spec(rho(i,j,k),temp(i,j,k),&(Y(i,j,k,0)),numPts,wtmp);
            for (int n=0; n<num_spec; ++n) {
              w(i,j,k,n) = wtmp[n];
            }
          });
        }
      }

      VisMF::Write(wdot,"WDOT1");

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        BL_PROFILE("COMPUTE_W_MAX");
        for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& box = mfi.tilebox();
          const auto& temp = temperature.array(mfi);
          const auto& Y    = mass_frac.array(mfi);
          const auto& rho  = density.array(mfi);
          const auto& w    = wdot.array(mfi);
          int numPts = box.numPts();

          For(box,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
            Real wtmp[21];
            W_specMAX(rho(i,j,k),temp(i,j,k),&(Y(i,j,k,0)),numPts,wtmp);
            for (int n=0; n<num_spec; ++n) {
              w(i,j,k,n) = wtmp[n];
            }
          });
        }
      }

      VisMF::Write(wdot,"WDOTM");

#endif
    }

    Finalize();

    return 0;
}

