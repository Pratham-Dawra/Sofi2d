
#include <stdlib.h>
#include <string.h>

#include "fd.h"
#include "util.h"
#include "logging.h"
#include "kiss_fftr.h"

void stfi_merge(int ishot, st_signals *signals, int ncplx)
{
    /* ishot,                   // shot number
     * float **signals,         // input/output source signature (overwrite)
     * int ncplx,               // size of complex buffer */

    double tmp_filt_re, tmp_filt_im;

    log_info("STFI: Merging filter for shot %d ...\n", ishot);

    // merging old and new filter to input traces
    for (int i = 0; i < ncplx; i++) {
        tmp_filt_re =
            signals->stfi_filt[(ishot - 1) * ncplx + i].r * signals->stfi_filt_new[(ishot - 1) * ncplx + i].r -
            signals->stfi_filt[(ishot - 1) * ncplx + i].i * signals->stfi_filt_new[(ishot - 1) * ncplx + i].i;
        tmp_filt_im =
            signals->stfi_filt[(ishot - 1) * ncplx + i].r * signals->stfi_filt_new[(ishot - 1) * ncplx + i].i +
            signals->stfi_filt[(ishot - 1) * ncplx + i].i * signals->stfi_filt_new[(ishot - 1) * ncplx + i].r;
        signals->stfi_filt_new[(ishot - 1) * ncplx + i].r = tmp_filt_re;
        signals->stfi_filt_new[(ishot - 1) * ncplx + i].i = tmp_filt_im;
    }

    return;
}
