/* vim: set et ts=4 sw=4 sts=4 ft=c:
 *
 * Copyright (C) 2016 Christopher Snowhill.  All rights reserved.
 * https://github.com/kode54/fft-resampler
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

/* Welcome to FFT convolution, with your host, kode54! Today, we'll be following
 * along with some code I mostly picked up from examples I discovered scattered
 * around various parts of the Internet. I chose to use kissfft because it's
 * fairly fast, and also under a less restrictive license than FFTW. */

#ifndef _SIMPLE_CONVOLVER_H_
#define _SIMPLE_CONVOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Fully opaque convolver state created and returned here, otherwise NULL on
 * failure. Users are welcome to change this to pass in a const pointer to an
 * impulse and its size, which will be copied and no longer needed upon return.
 * It is assumed that there will be one impulse per input channel, and that
 * each impulse will have one channel per output.
 *
 * Mode may be one of the following:
 * - 0: single 1 channel impulse, same number of input and output channels
 * - 1: single multi channel impulse, individual channels mapped 1:1 from
 *      inputs to outputs, same number of input and output channels
 * - 2: multiple multi-channel impulses, one impulse per input channel,
 *      containing one channel per output channel, and the results are
 *      summed together. */
void * convolver_create(const float * const* impulses, int impulse_size, int input_channels, int output_channels, int mode);

/* This function is for re-importing a modified impulse set into an existing
 * instance, with the same number of channels per input and output, so the
 * same number of impulses and channels per impulse. Useful if you are
 * creating filter impulses, and wish to restage with a newly generated filter
 * set. */
void convolver_restage(void *, const float * const* impulses);

/* Pass an instance of the convolver here to clean up when you're done with it */
void convolver_delete(void *);

/* This will clear the intermediate buffers of the convolver, useful for
 * restarting a stream with the same filter parameters. */
void convolver_clear(void *);

/* This returns the number of samples necessary for the convolver to generate
 * a block of output. When there is output buffered, it is a bad idea to fill
 * more input samples in. */
int convolver_get_free_count(void *);

/* Write the samples here, one n-channel set at a time. */
void convolver_write(void *, const float *);

/* This will return how many samples are buffered for output. */
int convolver_ready(void *);

/* This will return one n-channel set of samples at a time, or else silence
 * if the buffer is empty. */
void convolver_read(void *, float *);

#ifdef __cplusplus
}
#endif

#endif
