
/*------------------------------------------------------------------------
 * Copyright (C) 2011 For the list of authors, see file AUTHORS.
 *
 * This file is part of SOFI2D.
 *
 * SOFI2D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 *
 * SOFI2D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOFI2D. See file COPYING and/or
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 --------------------------------------------------------------------------*/

/*  ----------------------------------------------------------------------
 *  This is program SOFI2D.
 *  Parallel 2-D Viscoelastic Finite Difference Seismic Modeling *  using the Standard Staggered Grid (SSG)
 *
 *  PLEASE DO NOT DISTRIBUTE. PLEASE REFER OTHER PEOPLE TO :
 *
 *  Prof. Dr. Thomas Bohlen, Karlsruhe Institute of Technology,
 *  Geophysical Institute,
 *  Hertzstr. 16, 76187 Karlsruhe, Germany
 *  Phone/Fax: +49 (0)721 608 44416
 *  mailto:thomas.bohlen@kit.edu,
 *  http://www.gpi.kit.edu/
 *  http://www.gpi.kit.edu/SOFI2D.php
 *
 *  If you want to publish synthetic data calculated with this program please
 *  give a reference to the following paper:
 *  Bohlen, T., 2002, Parallel 3-D viscoelastic finite-difference seismic modelling,
 *  Computers @ Geopsciences, Vol. 28, No. 8, 887-889.
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar_struct.h"
#include "globvar_struct_fwi.h"
#include "logging.h"
#include "macros.h"
#include "enums.h"
#ifdef EBUG
#include "debug_buffers.h"
#endif

#include <unistd.h>

int main(int argc, char **argv)
{
    /* variables in main */
    int ishot, nshots, snapcheck;   /* Added ishot and nshots for multiple shots */
    int iter;
    clock_t cpu_time1 = 0, cpu_time = 0 ;
    FILE *log_fp = NULL;
    char ext[10];
    double time1 = 0.0, time2 = 0.0, time9 = 0.0;
    float *hc = NULL;
    st_acquisition acq = { };
    st_boundary nb;
    st_velocity vel;
    st_stress stress;
    st_model mod, testmod;
    st_model_av mod_av;
    st_visc_mem visco_mem;
    st_seismogram section, section_obs;
    st_signals signals = { };
    int ntr_loc = 0, nsrc_loc = 0;
    st_pml_coeff pml_coeff;
    st_pml_wfd pml_wfd;
    st_buffer velbuff, stressbuff;
    st_freq_velocity fourier_vel_fw, fourier_vel_back;
    /*finv is an instance of vector function created in util;*/
    double L2 = 0.0;
    int ntast=1;
    int lsnap = 0, nsnap = 0;
    int iteration = 0, cdf = 0, groupnum = 0, nf = 0;
    int it_group = 0;
    int ncplx= 0;


    /* FWI variables */
    /*float alpha_SL_old;
     * int fwi_run = 1;
     * int gradient_optimization = 1;
     * int steplength_search = 0; */
    /* Variables for step length calculation */
    //int step3 = 0, countstep;signals
    //int step1, step2, step3 = 0; // itests, iteste, stepmax, countstep;
    //float scalefac;

    /* declare struct for global variables */
    GlobVar gv = {.MPID = -1,.OUTNTIMESTEPINFO = 100,.NDT = 1,.IDX = 1,.IDY = 1 };

    /* declare struct for inversion variables */
    GlobVarInv vinv = {.ITERMAX = 1,.DTINV = 1,.WORKFLOW_STAGE = 1,.LBFGS_ITER_START = 1,
        .FWI_RUN = 1,.GRADIENT_OPTIMIZATION = 1,.FREQ_NR = 1};

    /* declare struct for acquisition variables */
    AcqVar acq = { };

    /* declare struct for wavefield variables */
    MemWavefield mpw = { };

    /* declare struct for model variables */
    MemModel mpm = { };

    /* declare struct for FWI wavefield and model variables */
//    if (gv.MODE == FWI) {
    MemInv minv = { };
//    }

    /* declare struct for performance measures */
    Perform perf = { };

    /* initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &gv.MPID);signals

    /* initialize logging */
    log_init(NULL);
    log_banner(LOG_SOFI);

    time1 = MPI_Wtime();
    cpu_time1 = clock();

    if (gv.MPID == 0) {
        if (argc != 2) {
            log_fatal
                ("Unexpected number of commandline arguments; ssignalsingle argument required: name of json parameter file.\n");
        }

        const char *fileinp = argv[1];

        /* check if parameter file can be opened */
        if (access(fileinp, R_OK) != 0) {
            log_fatal("Cannot open/read json parameter file %s.\n", fileinp);
        }

        /* check suffix of parameter file */
        if (!STRSTRCOMP(fileinp, ".json")) {
            log_fatal("Parameter file %s has no .json suffix.\n", fileinp);
        }

        /* read json parameter file */
        read_par_json(fileinp, &gv, &vinv);

    }

    /* exchange parameters between MPI processes */
    exchange_par(&gv, &vinv);

    /* read FWI workflow */
    if (gv.MODE == FWI) {
        if (vinv.USE_WORKFLOW) {
            read_workflow(&gv, &vinv);
        }

        switch (vinv.TIME_FILT) {
            case 1:
                vinv.F_LOW_PASS = vinv.F_LOW_PASS_START;
                break;
            /*read frequencies from file */
            case 2:
                filter_frequencies(&vinv);
                vinv.F_LOW_PASS = vinv.F_LOW_PASS_EXT[vinv.FREQ_NR];
                break;
        }
        log_info("F_LOW_PASS, F_LOW_PASS_START, ORDER: %f, %f, %d\n", vinv.F_LOW_PASS, vinv.F_LOW_PASS_START, vinv.ORDER);
        //if (vinv.TIME_FILT == 2) {
            /*read frequencies from file */
        //    filter_frequencies(&vinv);
            /* start with first low pass frequency */
        //    vinv.F_LOW_PASS_START = vinv.F_LOW_PASS[1];
        //} else {
        //    vinv.F_LOW_PASS = vector(1, 1);
        //    vinv.F_LOW_PASS[1] = vinv.F_LOW_PASS_START;
        //}
    }

    /* set logging verbosity */

    /* check file system/output directories */
    check_fs(&gv, &vinv);

    sprintf(ext, ".%i", gv.MPID);
    strcat(gv.LOG_FILE, ext);

    /* set up logging output */
    switch (gv.LOG) {
      case 0:                  /* logging to stdout/stderr for all ranks */
          log_fp = NULL;
          log_set_output(NULL);
          log_infoc(0, "Log messages sent to stdout/stderr on all MPI ranks.\n");
          break;
      case 1:                  /* logging to file for all ranks */
          log_infoc(0, "Now redirecting log messages to log file on all MPI ranks.\n");
          if ((log_fp = fopen(gv.LOG_FILE, "w")) == NULL) {
              log_fatal("Opening log file %s for writing failed.\n", gv.LOG_FILE);
          }
          log_set_output(log_fp);
          log_info("This is the log file %s generated by PE %d.\n", gv.LOG_FILE, gv.MPID);
          break;
      case 2:                  /* logging to stdout/stderr on rank 0, logging to file for all other ranks */
          if (0 == gv.MPID) {
              log_fp = NULL;
              log_set_output(NULL);
              log_info("Now redirecting log messages to log file on all MPI ranks except rank 0.\n");
          } else {
              if ((log_fp = fopen(gv.LOG_FILE, "w")) == NULL) {
                  log_fatal("Opening log file %s for writing failed.\n", gv.LOG_FILE);
              }
              log_set_output(log_fp);
              log_info("This is the log file %s generated by PE %d.\n", gv.LOG_FILE, gv.MPID);
          }
          break;
      default:
          log_warn("Unknown value %d for parameter LOG encountered; using LOG=0.\n", gv.LOG);
          gv.LOG = 0;
          log_fp = NULL;
          log_set_output(NULL);
          log_infoc(0, "Log messages sent to stdout/stderr on all MPI ranks.\n");
          break;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* domain decomposition */
    initproc(&gv);

    /* reading acquisition parameters */
    nshots = acq_read(&acq, &gv);

    /* output of parameters */
    if (gv.MPID == 0) {
        write_par(&gv, &vinv);
    }

    /* memory allocation of buffers */
    initmem(&mpm, &mpw, &minv, &gv, &vinv);

    /* initialize FD operators */
    initfd(&gv);

    /* Holberg coefficients for FD operators */
    hc = holbergcoeff(&gv);

    MPI_Barrier(MPI_COMM_WORLD);

    /* create model grids */
    readmod(&mpm, &minv, &gv, &vinv);


    MPI_Barrier(MPI_COMM_WORLD);

    /* check if the FD run will be stable and free of numerical dispersion */
   checkfd(hc, acq.srcpos, acq.nsrc, acq.recpos, &gv);

    /* calculate damping coefficients for CPMLs */
    if (gv.ABS_TYPE == 1) {
        PML_pro(&mpm, &gv);
    }

    /* calculate 2-D array for exponential damping of reflections at the edges of the numerical mesh */
    if (gv.ABS_TYPE == 2) {
        absorb(mpm.absorb_coeff, &gv);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*prepmod(&mpm, &gv);
     * 
     * MPI_Barrier(MPI_COMM_WORLD); */

    time2 = MPI_Wtime();
    log_infoc(0, "Starting time stepping around real time %4.2fs.\n", time2 - time1);

    /*---------------------------  Start inversion  ----------------------*/

    /*------------------------  loop over iterations  --------------------*/

    for (iter = 1; iter <= vinv.ITERMAX; iter++) {  /* fullwaveform iteration loop */
        if(gv.MPID==0) log_info("OUT1: FWI_RUN: %d, STEPLENGTH_SEARCH: %d, GRADIENT_OPTIMIZATION: %d\n", vinv.FWI_RUN, vinv.STEPLENGTH_SEARCH, vinv.GRADIENT_OPTIMIZATION);
        if (gv.MPID == 0) {
            time2 = MPI_Wtime();
            log_info("------------------------------------------------------------------\n");
            if (gv.MODE == FWI) {
                log_info("                   TD-FWI ITERATION %d \t of %d \n", iter, vinv.ITERMAX);
            } else {
                log_info("                        FD-SIMULATION \n");
            }
            log_info("------------------------------------------------------------------\n");
        }
        // At each iteration the workflow is applied
        if (gv.MODE == FWI) {
            if (vinv.USE_WORKFLOW) {
                apply_workflow(iter, &gv, &vinv);
            }

            init_grad(iter, &minv, &gv, &vinv);

        }

        if(gv.MPID==0) log_info("OUT2: FWI_RUN: %d, STEPLENGTH_SEARCH: %d, GRADIENT_OPTIMIZATION: %d\n", vinv.FWI_RUN, vinv.STEPLENGTH_SEARCH, vinv.GRADIENT_OPTIMIZATION);
        /*-----------------------------------------------------*/
        /*  While loop for Wolfe step length search            */
        /*-----------------------------------------------------*/
        while (vinv.FWI_RUN || vinv.STEPLENGTH_SEARCH || vinv.GRADIENT_OPTIMIZATION) {

            if(gv.MPID==0) log_info("FWI_RUN: %d\n", vinv.FWI_RUN);
            /*-----------------------------------------------------*/
            /*              Calculate Misfit and gradient          */
            /*-----------------------------------------------------*/
            if (vinv.FWI_RUN) {
                /* For the calculation of the material parameters between gridpoints
                 * they have to be averaged. For this, values lying at 0 and NX+1,
                 * for example, are required on the local grid. These are now copied from the
                 * neighbouring grids */
                prepmod(&mpm, &gv);

                MPI_Barrier(MPI_COMM_WORLD);

                if(gv.MODE == FWI) init_inv(iter, &minv, &gv, &vinv);

                /*----------------------  loop over multiple shots  ------------------*/
                for (ishot = 1; ishot <= nshots; ishot++) {
                    
                    gv.SOURCE_SHAPE = gv.SOURCE_SHAPE_OLD;
                    
                    
                    /*------------------------------------------------------------------------------*/
                    /*----------- Start of inversion of source time function -----------------------*/                       /* Do not Excute STF if this is a step length search run for Wolfe condition
                     * Therefore (gradient_optimization==1) is added. */
                    /*if (((INV_STF == 1) && ((iter == 1) || (do_stf == 1))) && (gradient_optimization == 1)) {
                        stf_inv();
                    }*/                        /*------------------------------------------------------------------------------*/ 

                    snapcheck = initsrc(ishot, nshots, &acq, &gv);
                    
                    /*------------------------------------------------------------------------------*/
                    /*---------- Start of Time Domain Filtering ------------------------------------*/
                    /*time domain filtering of the source signal */
                    /*if (((TIME_FILT == 1) || (TIME_FILT == 2)) && (SOURCE_SHAPE != 6) && (INV_STF == 0))
                    {
                        log_info("Time Domain Filter applied: Lowpass with corner frequency of %.2f Hz, order %d\n",F_LOW_PASS, ORDER);
                        timedomain_filt(signals, F_LOW_PASS, ORDER, nsrc_loc, ns, 1);
                    }

                    MPI_Barrier(MPI_COMM_WORLD);*/
                    /*------------------------------------------------------------------------------*/

                    /* initialize wavefield with zero */
                    zero_wavefield(iter, &mpw, &minv, &gv, &vinv);

                    /* determine block index boundaries for inner area and frame */
                    subgrid_bounds(1, gv.NX, 1, gv.NY, &gv);
                    /* forward propagation */
                    /*pshot = 0;
                    if (gv.METHOD) {
                        dummy = (1 / (finv[nf - 1] * gv.TAST * gv.DT));
                        ntast = (dummy);
                        if (!ntast)
                            ntast = 1;

                        if (gv.STFI) {
                            stfi(&acq, &nb, &nb_fix, &vel, &stress, &mod, &mod_av, &visco_mem, &section, &section_obs,
                                &signals, nsrc_loc, ntr_loc, &pml_coeff, &pml_wfd, &stressbuff, &velbuff,
                                &fourier_vel_fw, finv, &L2, nf, ntast, lsnap, nsnap, ishot, cdf, iteration, groupnum,
                                it_group, ncplx, &gv);
                        }
                    }*/
                /*****************************************************
             * Timestep loop (simulation)
             *****************************************************/
            /*timeloop(&nb, &nb_fix, &vel, &stress, &mod, &mod_av, &section,
                     acq.srcpos_loc, acq.recpos_loc, &signals, nsrc_loc,
                     mod.absorb_coeff, &pml_coeff, &pml_wfd, &stressbuff, &velbuff, &fourier_vel_fw,
                     finv, nf, ntast, ntr_loc, lsnap, nsnap, 0, ishot, 1, &gv); */

                    /*------------------------------------------------------------------------------*/
                    /*---------- Start loop over timesteps (forward model) -------------------------*/
                    time_loop(iter, ishot, snapcheck, hc, acq.srcpos_loc, acq.signals, acq.signals,
                              acq.nsrc_loc, 0, &acq, &mpm, &mpw, &minv, &gv, &vinv, &perf);

                    /* gather and output seismograms if applicable */
                    saveseis(ishot, &acq, &minv, &gv, &vinv);

                    /*------------------------------------------------------------------------------*/
                    /*---------- Inversion: Start inversion process --------------------------------*/
                    if (gv.MODE == FWI) {
                        inversion(iter, ishot, snapcheck, hc, &acq, &mpm, &mpw, &minv, &gv, &vinv, &perf);
                    }
                }
                /*----------------  end of loop over multiple shots  -----------------*/
                
                vinv.FWI_RUN = 0;
                vinv.STEPLENGTH_SEARCH = 0;
                vinv.GRADIENT_OPTIMIZATION = 0;
            }
        }
        vinv.FWI_RUN = 1;
    }

    /* deallocate memory */
    /*freemem(&mpm, &mpw, &gv);*/

    if (gv.SEISMO)
        free_imatrix(acq.recpos, 1, 3, 1, gv.NTRG);

    MPI_Barrier(MPI_COMM_WORLD);

    if (gv.MPID == 0) {
        perf.time_av_v_update = perf.time_av_v_update / (double)perf.infocounter;
        perf.time_av_s_update = perf.time_av_s_update / (double)perf.infocounter;
        perf.time_av_v_exchange = perf.time_av_v_exchange / (double)perf.infocounter;
        perf.time_av_s_exchange = perf.time_av_s_exchange / (double)perf.infocounter;
        perf.time_av_timestep = perf.time_av_timestep / (double)perf.infocounter;
        log_info("Approximate average times for\n");
        log_info("  velocity update: .. %.6lfs.\n", perf.time_av_v_update);
        log_info("  stress update: .... %.6lfs.\n", perf.time_av_s_update);
        log_info("  velocity exchange:  %.6lfs.\n", perf.time_av_v_exchange);
        log_info("  stress exchange: .. %.6lfs.\n", perf.time_av_s_exchange);
        log_info("  time step: ........ %.6lfs.\n", perf.time_av_timestep);
        cpu_time = clock() - cpu_time1;
        log_info("CPU time of program per PE: %.3lfs.\n", (double)cpu_time / (double)CLOCKS_PER_SEC);
        time9 = MPI_Wtime();
        log_info("Total real time of program: %.3lfs.\n", time9 - time1);
    }

    /* finalize logging */
    log_finalize();

    if (log_fp)
        fclose(log_fp);

    /* finalize MPI */
    MPI_Finalize();

    return EXIT_SUCCESS;
}
