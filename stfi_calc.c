
/*------------------------------------------------------------------------
 * Copyright (C) 2015 For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS3D.
 * 
 * IFOS3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS3D. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
--------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
*   Calculating source time function inversion (STFI) filter 
*  ----------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "fd.h"
#include "util.h"
#include "logging.h"
#include "kiss_fftr.h"

void stfi_calc(int ishot, st_seismogram *section, st_seismogram *section_obs, st_signals *signals, int ntr_loc,
               int nfft_min, GlobVar *gv)
{
    /* ishot,                   // shot number
     * float **section,         // modelled trace
     * float **section_obs,     // target trace
     * float **signals,         // filter
     * int ntr_loc,             // number of traces per MPI process
     * int nfft_min,            // min. length for Fourier transform
     * st_acquisition *acq      // struct containing acquisition parameters
     * GlobVar *gv);            // global variable struct */

    float **trace = NULL, **target = NULL;

    if (0 == gv->MYID_SHOT)
        log_info("STFI: Calculating filter for shot %d ...\n", ishot);

    switch (gv->ADJOINT_TYPE) {
      case 1:
          trace = section->vx;
          target = section_obs->vx;
          break;
      case 3:
          trace = section->vz;
          target = section_obs->vz;
          break;
      case 2:
          /* FALLTHRU */
      case 4:
          trace = section->vy;
          target = section_obs->vy;
          break;
      case 5:
          trace = section->p;
          target = section_obs->p;
          break;
    }

    // find next power of two for FFT; kissfft does not need a power of 2 
    // but definitely an even number of samples for a real-to-complex FFT
    int nfft = 2;
    while (nfft < gv->NS || nfft < nfft_min) {
        nfft *= 2;
    }

    int ncplx = nfft / 2 + 1;

    // allocate scalar buffers for frequency domain
    kiss_fft_scalar *trace_pad = (kiss_fft_scalar *) calloc(nfft, sizeof(kiss_fft_scalar));
    kiss_fft_scalar *target_pad = (kiss_fft_scalar *) calloc(nfft, sizeof(kiss_fft_scalar));
    // allocate complex buffers for frequency domain
    kiss_fft_cpx *trace_frq = (kiss_fft_cpx *) calloc(ncplx, sizeof(kiss_fft_cpx));
    kiss_fft_cpx *target_frq = (kiss_fft_cpx *) calloc(ncplx, sizeof(kiss_fft_cpx));

    // allocate scalar and complex buffers for filter calculation
    kiss_fft_scalar *AC = (kiss_fft_scalar *) calloc(ncplx, sizeof(kiss_fft_scalar));
    kiss_fft_scalar *tmp_AC = (kiss_fft_scalar *) calloc(ncplx, sizeof(kiss_fft_scalar));
    kiss_fft_cpx *CC = (kiss_fft_cpx *) calloc(ncplx, sizeof(kiss_fft_cpx));
    kiss_fft_cpx *tmp_CC = (kiss_fft_cpx *) calloc(ncplx, sizeof(kiss_fft_cpx));
    kiss_fft_scalar tmp_Err = 0;
    kiss_fft_scalar tmp_Err_max = 0;
    kiss_fft_scalar Err_max = 0;

    // forward transform state
    kiss_fftr_cfg kiss_fftr_state = kiss_fftr_alloc(nfft, 0, 0, 0);

    for (int j = 1; j <= ntr_loc; j++) {

        // pad input trace
        memcpy(trace_pad, &(trace[j][0]), gv->NS * sizeof(kiss_fft_scalar));

        // pad target trace
        memcpy(target_pad, &(target[j][0]), gv->NS * sizeof(kiss_fft_scalar));

        // perform forward Fourier transforms
        kiss_fftr(kiss_fftr_state, trace_pad, trace_frq);
        kiss_fftr(kiss_fftr_state, target_pad, target_frq);

        // calculate cross- and auto-correlations and determine maximum
        for (int i = 0; i < ncplx; i++) {
            tmp_CC[i].r += target_frq[i].r * trace_frq[i].r + target_frq[i].i * trace_frq[i].i;
            tmp_CC[i].i += target_frq[i].i * trace_frq[i].r - target_frq[i].r * trace_frq[i].i;
            tmp_Err = pow(trace_frq[i].r, 2) + pow(trace_frq[i].i, 2);
            tmp_Err_max += tmp_Err;
            tmp_AC[i] += tmp_Err;
        }
    }

    MPI_Barrier(gv->COMM_SHOT);
    MPI_Allreduce(tmp_CC, CC, 2 * ncplx, MPI_FLOAT, MPI_SUM, gv->COMM_SHOT);
    MPI_Allreduce(tmp_AC, AC, ncplx, MPI_FLOAT, MPI_SUM, gv->COMM_SHOT);
    MPI_Allreduce(&tmp_Err_max, &Err_max, 1, MPI_FLOAT, MPI_SUM, gv->COMM_SHOT);
    kiss_fft_scalar fac = Err_max * gv->STFI_EPS;

    // design matching filter
    for (int i = 0; i < ncplx; i++) {
        signals->stfi_filt_new[(ishot - 1) * ncplx + i].r = CC[i].r / (AC[i] + fac);
        signals->stfi_filt_new[(ishot - 1) * ncplx + i].i = CC[i].i / (AC[i] + fac);
    }

    free(trace_pad);
    trace_pad = NULL;
    free(target_pad);
    target_pad = NULL;

    free(trace_frq);
    trace_frq = NULL;
    free(target_frq);
    target_frq = NULL;

    free(AC);
    AC = NULL;
    free(tmp_AC);
    tmp_AC = NULL;
    free(CC);
    CC = NULL;
    free(tmp_CC);
    tmp_CC = NULL;

    free(kiss_fftr_state);
    kiss_fftr_state = NULL;

    return;
}
