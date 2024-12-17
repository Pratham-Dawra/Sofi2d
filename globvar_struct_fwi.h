
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

/*------------------------------------------------------------------------
 * globvar_struct.h - global variables of viscoelastic 2D FD program
 * generally, for the names of the global variables uppercase letters 
 * are used.
 * ----------------------------------------------------------------------*/

#ifndef GLOBVAR_STRUCT_FWI_H_INCLUDED
#define GLOBVAR_STRUCT_FWI_H_INCLUDED

//#include "macros.h"
//#include "enums.h"
//#include "memw_struct.h"

typedef struct {

    /* Inversion parameters */
    int ITERMAX;                // maximum number of iterations
    char DATA_DIR[STRING_SIZE]; // folder containing real data
    int ADJOINT_TYPE;           // choose seismic component for back-propagation (x & y, y only or x only)
    char MISFIT_LOG_FILE[STRING_SIZE];  // file name to save misfits
    int VELOCITY;               // switch to minimize particle velocities (instead of displacements)

    int INV_VP_ITER;            // iteration at which inversion for vp starts
    int INV_VS_ITER;            // iteration at which inversion for vs starts
    int INV_RHO_ITER;           // iteration at which inversion for density starts

    /* Output of inversion results */
    char INV_MODELFILE[STRING_SIZE];    // file name to output inverted model
    int NFSTART;                // iteration step for first model output
    int NF;                     // iteration increment for model output
    char JACOBIAN[STRING_SIZE]; // file name to output gradients
    int NFSTART_JAC;            // first iteration step to save gradients
    int NF_JAC;                 // gradients are saved every nf_jac iteration step

    /* Workflow setting */
    int USE_WORKFLOW;           // switch to use workflow file
    int WORKFLOW_STAGE;         // 
    char FILE_WORKFLOW[STRING_SIZE];    // file containing workflow parameters
    float **WORKFLOW;           //
    int WORKFLOW_LINES;         //
    char WORKFLOW_HEADER[STRING_SIZE];  //
//    int WORKFLOW_MAX_VAR        // maximum number of parameters

    /* Ápproximate Hessian */
    int EPRECOND;               // switch to activate approximated Hessian for preconditioning
    float EPSILON_WE;           // defines a water level to stabilize the approximated Hessian
//    float EPSILON_WE_SH;                          !!! NO SH !!!
    int EPRECOND_PER_SHOT;      // switch to calculate approximated Hessian for each shot
//    int EPRECOND_PER_SHOT_SH;                     !!! NO SH !!!
    int EPRECOND_ITER;
    int EPRECOND_MAX;           // maximum EPRECOND in workflow.

    /* PCG and L-BFGS */
    int GRAD_METHOD;            // choose method gradient calculation (PCG or L-BFGS)
    int WOLFE_CONDITION;        // switch to use step length search based on Wolfe condition
    int WOLFE_NUM_TEST;         // maximum number of test calculations
    int WOLFE_TRY_OLD_STEPLENGTH;   // switch to try old steplength first
    float WOLFE_C1_SL;          // c1 parameter of Wolfe condition
    float WOLFE_C2_SL;          // c2 parameter of Wolfe condition
    /* Variables for L-BFGS */
    int N_LBFGS;                // 
    int LBFGS_STEP_LENGTH;      // switch to use second step length search
    float LBFGS_SCALE_GRADIENTS;    // NOT EXPLAINED IN MANUAL !!!
    int LBFGS_ITER_START;       //
    float ALPHA_SL_OLD;         // switch
    int FWI_RUN;                // switch
    int GRADIENT_OPTIMIZATION;  //switch
    int STEPLENGTH_SEARCH;      // switch
    /* Variables for step length calculation */
    int STEP3;                  //
    int COUNTSTEP;              //

    /* for Wolfe condition */
    /*int steplength_search = 0;
     * int FWI_run = 1;
     * int gradient_optimization = 1;
     * float alpha_SL_min = 0, alpha_SL_max = 0, alpha_SL = 1.0;
     * float alpha_SL_old;
     * float L2_SL_old = 0, L2_SL_new = 0;
     * float c1_SL = 1e-4, c2_SL = 0.9;
     * int wolfe_status;
     * int wolfe_sum_FWI = 0;
     * int wolfe_found_lower_L2 = 0;
     * float alpha_SL_FS;
     * float L2_SL_FS;
     * int use_wolfe_failsafe = 0;
     * int wolfe_SLS_failed = 0; */

    /* step length estimation */
    float EPS_SCALE;            // scalar for step length calculation
    int STEPMAX;                // number of attemps to find step length
    float SCALEFAC;             // factor to divide steplength by
    int TESTSHOT_START;         // lowest shot number for calculation of misfit norm of data residuals
    int TESTSHOT_END;           // highest shot number for calculation of misfit norm of data residuals
    int TESTSHOT_INCR;          // shot increment for calculation of misfit norm of data residuals
    int NO_OF_TESTSHOTS;        // number of test shots
    int ITESTSHOT;              //
    float EPSILON;              //
    float EPSILON_u;            //
    float EPSILON_rho;          //

    /* Misfit definition */
    int LNORM;                  // norm used for misfit definition
    int NORMALIZE;              // switch to normalize synthetic and measured data before calculation of residuals
    int DTINV;                  // increment of time samples used for misfit calculation
    int NTDTINV;                // number of timesteps used for misfit calculation (=ceil(NT/DTINV))
    float WATERLEVEL_LNORM8;    // finite water-level to keep division regular
    float L2;                   //
    float L2_ALL_SHOTS;         //
    float ENERGY;               //
    float ENERGY_ALL_SHOTS;     //
    int WRITE_DIFF;             // write difference between measured and synthetic seismogramms to disk

    /* Abort criterion */
    float PRO;                  // abort criterium per FWI stage [%]
    int MIN_ITER;               // minimum number of itereations to run (per frequency)

    /* Source time function inversion */
    int INV_STF;                // invert for source time function
    int N_STF;                  // incerment between STF iteration steps
    int N_STF_START;            // iteration step where STF starts 
    char PARA[STRING_SIZE];     // STF optimization approach
    int TAPER_STF;              // switch to taper source signal
    int TRKILL_STF;             // switch to apply trace killing for STF
    char TRKILL_FILE_STF[STRING_SIZE];  // file listing traces to be killed
    int TRKILL_STF_OFFSET;      // activate offset based trace killing for STF
    int TRKILL_STF_OFFSET_INVERT;   // invert sense of trace killing (traces between lower and upper limit are kept)
    float TRKILL_STF_OFFSET_LOWER;  // lower offset limit for trace killing [m]
    float TRKILL_STF_OFFSET_UPPER;  // upper offset limit for trace killing [m]
    int STF_FULL;               //         !!! NOT IN MANUAL - needed ???
    int KILLED_TRACES;          //
    int KILLED_TRACES_TESTSHOTS;    //

    /* Gradient smoothing */
    int GRAD_FILTER;            // switch to apply gradient smoothing
    int FILT_SIZE_GRAD;         // filter size [grid cells]
    int GRAD_FILT_WAVELENGTH;   // wavelength dependent filter size []
    float A;                    // weighting factor

    /* Gradient taper */
    int SWS_TAPER_GRAD_VERT;    // switch to apply vertical taper
    int SWS_TAPER_GRAD_HOR;     // switch to apply horizontal taper
    int GRADT1;                 // first corner point of taper [grid cells]
    int GRADT2;                 // second corner point of taper [grid cells]
    int GRADT3;                 // third corner point of taper [grid cells]
    int GRADT4;                 // fourth corner point of taper [grid cells]
    int SWS_TAPER_GRAD_SOURCES; // switch to apply cylindrical taper around sources
    int SWS_TAPER_CIRCULAR_PER_SHOT;    // switch to apply cylindrical taper around each single source
    int SRTSHAPE;               // switch to define shape of taper
    float SRTRADIUS;            // radius of taper [m]
    int FILTSIZE;               // size where taper is zero around source [grid cells]
    int SWS_TAPER_FILE;         // switch to use externally defined taper from file
    int SWS_TAPER_FILE_PER_SHOT;    // switch to use externally defined taper for each shot from file
    char TAPER_FILE_NAME[STRING_SIZE];  // externally defined taper file
    int SWS_TESTSHOT;           //

    /* Spatial filter of gradients */
    int SPATFILTER;             // switch to apply spatial Gaussian filter
    int SPAT_FILT_SIZE;         //
    int SPAT_FILT_1;            //
    int SPAT_FILT_ITER;         //

    /* Frequency filtering */
    int TIME_FILT;              // switch to apply frequency filter
    char FREQ_FILE[STRING_SIZE];    // file providing individual frequency
    int NFREQ;                  // number of frequencies in FREQ_FILE
    int FREQ_NR;                // actual frequency to read
    float F_LOW_PASS;           // actual frequency of low pass filter
    float *F_LOW_PASS_EXT;      // vector containing frequencies of low pass filter
    //int ZERO_PHASE   -> mentioned in manual but doesn't exist in code
    float F_LOW_PASS_START;     // lower limit of butterworth low pass filter [Hz]
    float F_LOW_PASS_END;       // upper limit of butterworth low pass filter [Hz]
    float F_LOW_PASS_INCR;      // bandwidth increase per inversion step [Hz]
    int ORDER;                  // half order of butterworth low pass filter
    float F_HIGH_PASS;          // upper limit of an additional butterworth high pass filter [Hz]
    int WRITE_FILTERED_DATA;    // write the time filtered measured data to disk

    /* Time windowing */
    int TIMEWIN;                // switch to apply time windowing in the time series
    int TW_IND;                 // switch to read entire window from file
    char PICKS_FILE[STRING_SIZE];   // file containg picks for time windowing (one file per shot)
    float TWLENGTH_PLUS;        // length of the time window after picked time []
    float TWLENGTH_MINUS;       // length of the time window before picked time []
    float GAMMA;                // damping factor

    /* Trace killing */
    int TRKILL;                 // switch to apply trace killing
    char TRKILL_FILE[STRING_SIZE];  // file containing trace kill table
    int TRKILL_OFFSET;          // activate offset based trace killing
    float TRKILL_OFFSET_LOWER;  // lower offset limit for trace killing [m]
    float TRKILL_OFFSET_UPPER;  // upper offset limit for trace killing [m]

    /* Limits for model parameters */
    float VPUPPERLIM;           // maximum value for vp
    float VPLOWERLIM;           // minimum value for vp
    float VSUPPERLIM;           // maximum value for vs
    float VSLOWERLIM;           // minimum value for vs
    float RHOUPPERLIM;          // maximum value for density
    float RHOLOWERLIM;          // minimum value for density

    /* Limit update of model parameters */
    int S;                      // individual limits for model parameter updates
    float S_VP;                 // maximum deviation from starting model [%]
    float S_VS;                 // maximum deviation from starting model [%]
    float S_RHO;                // maximum deviation from starting model [%]
    float VP_AVG;               // average velocity (P-wave)
    float VS_AVG;               // average velocity (S-wave)
    float RHO_AVG;              // average density
    float C_VP;                 // squared average velocity (P-wave)
    float C_VS;                 // squared average velocity (S-wave)
    float C_RHO;                // squared average density

    float VP_VS_RATIO;          // minimum vp/vs-ratio

    int MODEL_FILTER;           // switch to apply smoothing of vp and vs models
    int FILT_SIZE;              // filter length [grid cells]

    //int SOURCE_SHAPE_SH;                          !!! NO SH !!!
    // int ACOUSTIC;                                !!! move to globvar_struct_fw
    // int WAVETYPE; // 1=PSV, 2=SH, 3=PSV & SH     !!! NO SH !!!
    // int JOINT_INVERSION_PSV_SH_TYPE;             !!! NO SH !!!
    // int JOINT_EQUAL_WEIGHTING;                   !!! NO SH !!!
    // float PHI;   -> PLANE_WAVE_ANGLE in GLOB_VAR
    // float JOINT_INVERSION_PSV_SH_ALPHA_VS;       !!! NO SH !!!
    // float JOINT_INVERSION_PSV_SH_ALPHA_RHO;      !!! NO SH !!!
    // int TAPER;                                   !!! NOT USED IN CODE !!!
    // int TAPERLENGTH;                             !!! NOT USED IN CODE !!!
    // int PARAMETERIZATION;                       // choose parameterization (velocities, impedance or Lamé)  ??? drop ???
    // int FORWARD_ONLY;                            !!! make enum (in GlobVar) !!!
    // int INVTYPE;                                // not mentioned in manual - makes no sense in code !!!

} GlobVarInv;

#endif
