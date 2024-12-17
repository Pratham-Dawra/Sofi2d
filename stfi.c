
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
 *   Source time function inversion (STFI)
 *  ----------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "fd.h"
#include "util.h"
#include "logging.h"
#include "seismo_shift.h"
#include "kiss_fftr.h"

void stfi(AcqVar*acq, st_seismogram *section,
          st_seismogram *section_obs, st_signals *signals, int nsrc_loc, int ntr_loc,
          st_buffer *stressbuff, const float *finv, double *L2, int nf, int ishot, int cdf,
          int iteration, int groupnum, int it_group, int ncplx, GlobVar *gv, int iter, int snapcheck,
          float *hc,float **srcpos_loc, float **signalx, float **signaly, int nsrc, int sw,
          MemModel *mpm, MemWavefield *mpw, MemInv *minv, GlobVarInv *vinv, Perform *perf)
{

    if (2 == it_group || iteration == gv->ITMIN || gv->RAND_SHOT) {
        gv->STFI_CALC = 1;
    }

    /* Source Time Function Inversion (STFI) */
    if (gv->STFI_CALC) {

        if (0 == gv->MYID_SHOT) log_info("Starting STFI for shot %d ...\n", ishot);

        /*****************************************************
        * Timestep loop (simulation for STFI)
        *****************************************************/
        time_loop(iter, ishot, snapcheck, hc, srcpos_loc, signalx, signaly, nsrc, sw, acq, mpm, mpw, minv ,gv ,vinv, perf);

        /* STFI: Calculate filter */
        seismo_shift(section, ntr_loc, acq->srcpos[5][ishot], gv);
        calc_res(ishot, section, section_obs, signals, ntr_loc, finv, nf, L2, 0, groupnum, acq, gv);
    }

    if (2 == gv->STFI && 1 != cdf && nsrc_loc > 0) {
        stfi_merge(ishot, signals, ncplx);
    }
    /*Copy new derived filter in filter application buffer */
    memcpy(signals->stfi_filt, signals->stfi_filt_new, acq->nsrc * ncplx * sizeof(kiss_fft_cpx));

    if (nsrc_loc > 0) {
        /* STFI: Apply filter */
        stfi_apply(ishot, signals, nsrc_loc, 1024, gv);
    }
    if (gv->SIGOUT && (1 == cdf || gv->RAND_SHOT || 2 == gv->STFI)) {
        savesig(signals->fw, acq, nsrc_loc, ishot, iteration, 1, gv);
    }
    if (2 == gv->STFI)
        gv->STFI_CALC = 1;

    /* initialize wavefield with zero */

    zero_wavefield( iter, mpw, minv, gv, vinv);
    return;
}
