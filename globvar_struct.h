
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

#ifndef GLOBVAR_STRUCT_H_INCLUDED
#define GLOBVAR_STRUCT_H_INCLUDED

#include "macros.h"
#include "enums.h"
#include "memw_struct.h"

typedef void (*FDop_s_fct)(int i, int j, float *vxx, float *vyx, float *vxy, float *vyy, MemWavefield *mpw);
typedef void (*FDop_v_fct)(int i, int j, float *sxx_x, float *sxy_x, float *sxy_y, float *syy_y, MemWavefield *mpw);
typedef void (*FDop_ac_s_fct)(int i, int j, float *vxx, float *vyy, MemWavefield *mpw);
typedef void (*FDop_ac_v_fct)(int i, int j, float *sxx_x, float *syy_y, MemWavefield *mpw);

typedef struct {

    // FD Params
    RUNMODE MODE;               // run mode (Forward only, FWI) 
    WEQTYPE WEQ;                // wave equation
    float DH;                   // spacial increment [m]
    int FDORDER;                // spatial FD order
    int ND;                     // gv->FDORDER / 2
    float TIME;                 // time (of modelling) [s]
    float DT;                   // time increment (of modelling) [s]
    float DTDH;                 // DT/DH
    int FDORDER_TIME;           // temporal FD order
    int NT;                     // number of timesteps (=iround(TIME/DT))
    FDop_s_fct FDOP_S;          // function pointer for FD operator
    FDop_v_fct FDOP_V;          // function pointer for FD operator
    FDop_ac_s_fct FDOP_AC_S;    // function pointer for FD operator (acoustic cases)
    FDop_ac_v_fct FDOP_AC_V;    // function pointer for FD operator (acoustic cases)
    int MAXRELERROR;            // switch of maximum relative group velocity error

    // MPI-variables
    int MPID;                   // ID of processor
    int MPID_SHOT;                  // ID of processor in COMM_SHOT domain
    int NPROC;                  // number of processors (=NPROCX*NPROCY; also number of MPI processes)
    int NPROCX;                 // number of processors in x-direction
    int NPROCY;                 // number of processors in y-direction
    int POS[3];                 // processor location in the logical proc array (1=x, 2=y, slot 0 unused)
    int INDEX[5];               // ID of neighboring processes (1=left, 2=right, 3=upper, 4=lower, slot 0 unused)
    int GGRID[5];               // global grid points that this process ID actually handles (1=left, 2=right, ...)
    int *GX;                    // subgrid array in x-direction
    int *GY;                    // subgrid array in y-diection
    int NX;                     // number of grid points in x-direction for this processor (local)
    int NY;                     // number of grid points in y-direction for this processor (local)
    int SNAPIDX[5];             // start/end indices for snapshots for this processor (1=left, 2=right, ...)
    char *BUFF_ADDR;            // buffer for buffering messages
    int BUFFSIZE;               // size of buffer for buffering messages
    const int TAG1;             // MPI tag
    const int TAG2;             // MPI tag
    const int TAG3;             // MPI tag
    const int TAG4;             // MPI tag
    const int TAG5;             // MPI tag
    const int TAG6;             // MPI TAG1

    // Boundary
    int FREE_SURF;              // switch to apply free surface at the top of the model
    int BOUNDARY;               // switch to apply periodic boundary condition at edges
    int ABS_TYPE;               // type of the absorbing boundary
    int FW;                     // width of absorbing frame [gridpoints]
    float DAMPING;              // attenuation at the edges of the grid [%]
    // PML-Parameters
    float NPOWER;               // exponent for calculation of damping profile
    float K_MAX_CPML;           // 
    float FPML;                 // dominant signal frequency (usually FC) [Hz]
    float VPPML;                // attenuation velocity within the PML boundary [m/s]

    // Models
    int READMOD;                // switch to read model parameters from MFILE
    char MFILE[STRING_SIZE];    // model file name
    int NXG;                    // number of grid points in x-direction (global)
    int NYG;                    // number of grid points in y-direction (global)
    float VPMIN;                // min P-wave velocity
    float VPMAX;                // max P-wave velocity
    float VSMIN;                // min S-wave velocity
    float VSMAX;                // max S-wave velocity
    // Attenuation
    int L;                      // number of relaxation parameters
    float F_REF;                // reference frequency for dispersion [Hz]
    float *FL;                  // frequency of each relaxation parameters [Hz]
    float TAU;                  // ratio of retardation and relaxation time

    // Source
    int SOURCE_TYPE;            // type of source
    int SOURCE_SHAPE;           // shape of source-signal
    int SOURCE_SHAPE_OLD;       // old SOURCE_SHAPE needed for STF (FWI)
    int SOURCE_TOPO;            // switch to place sources along topography
    char SIGNAL_FILE[STRING_SIZE];  // name of external signal file
    int SRCREC;                 // switch to read source parameters from external source file
    char SOURCE_FILE[STRING_SIZE];  // name of source parameter file
    int RUN_MULTIPLE_SHOTS;     // multiple shots modeled simultaneously (0) or individually; added for multiple shots
    float PLANE_WAVE_DEPTH;     // depth of plane wave excitation [meter]
    float PLANE_WAVE_ANGLE;     // dip of plane wave from vertical [degree]
    float TS;                   // duration of source signal [in second]
    float XS;                   // Source location
    float YS;                   // Source location

    // Receiver
    int READREC;                // switch to read receiver positions from file
    char REC_FILE[STRING_SIZE]; // name of external receiver file
    int REC_TOPO;               // switch to place receivers along topography
    float REFREC[4];            // reference point for receiver coordinate system
    float XREC1;                // x-position of first receiver [m]
    float XREC2;                // x-position of last receiver [m]
    float YREC1;                // y-position of first receiver [m]
    float YREC2;                // y-position of last receiver [m]
    float NGEOPH;               // distance between two adjacent receivers [gridpoints]; in auto mode NGEOPH will be 
                                // calculated from model dimensions, type integer is incorrect
    int NTRG;                   // global number of receiver
    int NTR;                    // number of receivers per PE
    int REC_ARRAY;              // number of receivers in 1D receiver array
    float REC_ARRAY_DEPTH;      // depth of first plane [m] 
    float REC_ARRAY_DIST;       // increment between receiver planes [m]
    int DRX;                    // increment between receivers in each plane [gridpoints]

    // Seismograms
    int SEISMO;                 // switch to output components of seismograms
    int NDT;                    // sampling rate of seismograms [timesteps DT]
    int NS;                     // number of samples of seismogram
    int SEIS_FORMAT;            // data output format for seismograms
    char SEIS_FILE[STRING_SIZE];    // name of output file of seismograms
    float **SECTIONVX;          // buffer for seismogram output (vx-component)
    float **SECTIONVY;          // buffer for seismogram output (vy-component)
    float **SECTIONP;           // buffer for seismogram output (p-component)
    float **SECTIONCURL;        // buffer for seismogram output (curl-component)
    float **SECTIONDIV;         // buffer for seismogram output (div-component)
    float **SEISMO_FULLDATA;    // buffer for merge of seismograms

    // Snapshots
    int SNAP;                   // switch to output of snapshots
    int SNAP_FORMAT;            // data output format for snapshots
    char SNAP_FILE[STRING_SIZE];    // name of output file of snapshots
    float TSNAP1;               // first snapshot [s] (IFOS: SNAPSHOT_START)
    float TSNAP2;               // last snapshot [s] (IFOS: SNAPSHOT_END)
    float TSNAPINC;             // increment between snapshots [s] (IFOS: SNAPSHOT_INCR)
    int SNAPSHOT_START;         // first shot to write snapshots
    int SNAPSHOT_END;           // last shot to write snapshots
    int SNAPSHOT_INCR;          // shot increment to write snapshots
    int IDX;                    // increment in x-direction [gridpoints]
    int IDY;                    // increment in y-direction [gridpoints]
    int SNAPIDCHECK;            // reset IDX and IDY if MODE=FWI

    // Others
    int WRITE_MODELFILES;       // switch to output model files
    int SIGOUT;                 // switch to output source wavelet
    int SIGOUT_FORMAT;          // data output format for source wavelet
    char SIGOUT_FILE[STRING_SIZE];  // name of output file of source wavelet
    int LOG;                    // switch to output logging information
    char LOG_VERBOSITY[STRING_SIZE];    // log output level (verbosity)       
    char LOG_FILE[STRING_SIZE]; // name of output file of logging information
    int OUTNTIMESTEPINFO;       // every OUTNTIMESTEPINFO th timestep, information on the time step will be given to screen/file

    // Inversion parameters
    int STFI;                   // switch to apply source time function inversion (STFI)
    int STFI_CALC;              // switch to calculate STFI
    int ITMIN;                  // minimum number of iteration
    int ITMAX;                  // maximum number of iteration

} GlobVar;

#endif
