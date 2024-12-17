
#include <stdlib.h>
#include <string.h>

#include "fd.h"
#include "util.h"
#include "logging.h"
#include "kiss_fftr.h"
//#include "stfi.h"

void stfi_apply(int ishot, st_signals *signals, int nsrc_loc, int nfft_min, GlobVar *gv)
{
    /* ishot,                   // shot number
     * float **signals,         // input/output source signature (overwrite)
     * int nsrc_loc,            // number of sources per MPI process
     * int nfft_min,            // min. length for Fourier transform
     * GlobVar *gv);            // global variable struct */

    log_info("STFI: Applying filter for shot %d ...\n", ishot);

    // find next power of two for FFT; kissfft does not need a power of 2 
    // but definitely an even number of samples for a real-to-complex FFT
    int nfft = 2;
    while (nfft < gv->NS || nfft < nfft_min) {
        nfft *= 2;
    }

    int ncplx = nfft / 2 + 1;

    // allocate scalar buffers for frequency domain
    kiss_fft_scalar *trace_pad = (kiss_fft_scalar *) calloc(nfft, sizeof(kiss_fft_scalar));
    kiss_fft_scalar *trace_filt_pad = (kiss_fft_scalar *) calloc(nfft, sizeof(kiss_fft_scalar));
    // allocate complex buffers for frequency domain
    kiss_fft_cpx *trace_frq = (kiss_fft_cpx *) calloc(ncplx, sizeof(kiss_fft_cpx));
    kiss_fft_cpx *trace_filt_frq = (kiss_fft_cpx *) calloc(ncplx, sizeof(kiss_fft_cpx));

    // forward transform state
    kiss_fftr_cfg kiss_fftr_state = kiss_fftr_alloc(nfft, 0, 0, 0);
    kiss_fftr_cfg kiss_fftri_state = kiss_fftr_alloc(nfft, 1, 0, 0);

    for (int j = 1; j <= nsrc_loc; j++) {

        // pad input trace
        memcpy(trace_pad, &(signals->fw[j][0]), gv->NS * sizeof(kiss_fft_scalar));

        // perform forward Fourier transforms
        kiss_fftr(kiss_fftr_state, trace_pad, trace_frq);

        // apply matching filter to input traces
        for (int i = 0; i < ncplx; i++) {
            trace_filt_frq[i].r = trace_frq[(j - 1) * ncplx + i].r * signals->stfi_filt[(ishot - 1) * ncplx + i].r -
                trace_frq[(j - 1) * ncplx + i].i * signals->stfi_filt[(ishot - 1) * ncplx + i].i;
            trace_filt_frq[i].i = trace_frq[(j - 1) * ncplx + i].r * signals->stfi_filt[(ishot - 1) * ncplx + i].i +
                trace_frq[(j - 1) * ncplx + i].i * signals->stfi_filt[(ishot - 1) * ncplx + i].r;
        }

        // perform inverse Fourier transforms
        kiss_fftri(kiss_fftri_state, trace_filt_frq, trace_filt_pad);

        // reduce trace to original number of samples
        for (int k = 0; k <= gv->NS-1; k++) {
            signals->fw[j][k] = trace_filt_pad[k - 1] / nfft;
        }
    }

    free(trace_pad);
    trace_pad = NULL;
    free(trace_filt_pad);
    trace_filt_pad = NULL;

    free(trace_frq);
    trace_frq = NULL;
    free(trace_filt_frq);
    trace_filt_frq = NULL;

    free(kiss_fftr_state);
    kiss_fftr_state = NULL;
    free(kiss_fftri_state);
    kiss_fftri_state = NULL;

    return;
}
