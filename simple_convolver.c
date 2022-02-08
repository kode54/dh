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

#include "simple_convolver.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Here comes the magic import header! */

#ifdef USE_FFTW
#include <fftw3.h>
#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#include <mm_malloc.h>
#else
#include "kissfft/kiss_fftr.h"
#endif

#ifdef __APPLE__
#define _mm_malloc(a, b) _memalign_malloc(a, b)
static void *_memalign_malloc(size_t size, size_t align) {
	void *ret = NULL;
	if(posix_memalign(&ret, align, size) != 0) {
		return NULL;
	}
	return ret;
}
#define _mm_free(a) free(a)
#endif

#if !defined(USE_FFTW) && defined(__APPLE__)
static inline int _malloc_dspsplitcomplex(DSPSplitComplex *out, size_t fftlen) {
	fftlen = (fftlen / 2) + 1;
	out->realp = _mm_malloc(sizeof(float) * fftlen, 16);
	out->imagp = _mm_malloc(sizeof(float) * fftlen, 16);
	if(out->realp == NULL || out->imagp == NULL) return -1;
	return 0;
}

static inline void _free_dspsplitcomplex(DSPSplitComplex *cpx) {
	if(cpx->realp) _mm_free(cpx->realp);
	if(cpx->imagp) _mm_free(cpx->imagp);
}
#endif

/* Only the simplest of state information is necessary for this, and it's
 * designed around a simple usage case. As many samples as you can generate,
 * input one at a time, and output pulled one at a time as well. */

typedef struct convolver_state {
	int fftlen; /* size of FFT */
	int impulselen; /* size of impulse */
	int fftlenover2; /* half size of FFT, rounded up */
#if !defined(USE_FFTW) && defined(__APPLE__)
	int fftlenlog2; /* log2 of FFT size */
#endif
	int stepsize; /* size of overlapping steps */
	int buffered_in; /* how many input samples buffered */
	int buffered_out; /* output samples buffered */
	int inputs; /* Input channels */
	int outputs; /* Output channels */
	int mode; /* Mode */
#ifdef USE_FFTW
	fftwf_plan *p_fw, p_bw; /* forward and backwards plans */
	fftwf_complex *f_in, *f_out, **f_ir; /* input, output, and impulse in frequency domain */
#elif defined(__APPLE__)
	FFTSetup setup; /* setup */
	DSPSplitComplex f_in, f_out, *f_ir; /* input, output, and impulse in frequency domain */
#else
	kiss_fftr_cfg cfg_fw, cfg_bw; /* forward and backwards instances */
	kiss_fft_cpx *f_in, *f_out, **f_ir; /* input, output, and impulse in frequency domain */
#endif
	float *revspace, **outspace, **inspace; /* reverse, output, and input work space */
} convolver_state;

/* Fully opaque convolver state created and returned here, otherwise NULL on
 * failure. Users are welcome to change this to pass in a const pointer to an
 * impulse and its size, which will be copied and no longer needed upon return.
 * It is assumed that there will be one impulse per input channel, and that
 * each impulse will have one channel per output. */

void *convolver_create(const float *const *impulse, int impulse_size, int input_channels, int output_channels, int mode) {
	convolver_state *state;
	int fftlen, total_channels, i, j, k;

	if(mode < 0 || mode > 2)
		return 0;

	if((mode == 0 || mode == 1) && input_channels != output_channels)
		return 0;

	state = (convolver_state *)calloc(1, sizeof(convolver_state));

	if(!state)
		return NULL;

	state->mode = mode;
	state->inputs = input_channels;
	state->outputs = output_channels;
	if(mode == 0)
		total_channels = 1;
	else if(mode == 1)
		total_channels = input_channels;
	else if(mode == 2)
		total_channels = input_channels * output_channels;

	state->stepsize = 512;
	state->impulselen = impulse_size;

	fftlen = state->impulselen + state->stepsize + 1;

	{
		// round up to a power of two
		int pow = 1;
		while(fftlen > 2) {
			pow++;
			fftlen /= 2;
		}
		fftlen = 2 << pow;
		state->fftlenover2 = 1 << pow;
#ifndef USE_FFTW
		state->fftlenlog2 = pow + 1;
#endif
	}

	state->fftlen = fftlen;
	state->buffered_in = 0;
	state->buffered_out = 0;

	/* Prepare arrays for multiple inputs */
	/* And we use kissfft's aligned malloc functions/macros to allocate these things. */

#ifdef USE_FFTW
	if((state->f_in = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (fftlen / 2 + 1))) == NULL)
#elif defined(__APPLE__)
	if(_malloc_dspsplitcomplex(&state->f_in, fftlen) < 0)
#else
	if((state->f_in = (kiss_fft_cpx *)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx) * (fftlen / 2 + 1))) == NULL)
#endif
		goto error;

#ifdef USE_FFTW
	if((state->f_out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (fftlen / 2 + 1))) == NULL)
#elif defined(__APPLE__)
	if(_malloc_dspsplitcomplex(&state->f_out, fftlen) < 0)
#else
	if((state->f_out = (kiss_fft_cpx *)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx) * (fftlen / 2 + 1))) == NULL)
#endif
		goto error;

#ifdef USE_FFTW
	if((state->f_ir = (fftwf_complex **)calloc(sizeof(fftwf_complex *), total_channels)) == NULL)
#elif defined(__APPLE__)
	if((state->f_ir = (DSPSplitComplex *)calloc(sizeof(DSPSplitComplex), total_channels)) == NULL)
#else
	if((state->f_ir = (kiss_fft_cpx **)calloc(sizeof(kiss_fft_cpx *), total_channels)) == NULL)
#endif
		goto error;
	for(i = 0; i < total_channels; ++i) {
#ifdef USE_FFTW
		if((state->f_ir[i] = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (fftlen / 2 + 1))) == NULL)
#elif defined(__APPLE__)
		if(_malloc_dspsplitcomplex(&state->f_ir[i], fftlen) < 0)
#else
		if((state->f_ir[i] = (kiss_fft_cpx *)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx) * (fftlen / 2 + 1))) == NULL)
#endif
			goto error;
	}

#ifdef USE_FFTW
	if((state->revspace = (float *)fftwf_malloc(sizeof(float) * fftlen)) == NULL)
#elif defined(__APPLE__)
	if((state->revspace = (float *)_mm_malloc(sizeof(float) * fftlen, 16)) == NULL)
#else
	if((state->revspace = (float *)KISS_FFT_MALLOC(sizeof(float) * fftlen)) == NULL)
#endif
		goto error;

	if((state->outspace = (float **)calloc(sizeof(float *), output_channels)) == NULL)
		goto error;
	for(i = 0; i < output_channels; ++i) {
#ifdef USE_FFTW
		if((state->outspace[i] = (float *)fftwf_malloc(sizeof(float) * fftlen)) == NULL)
#elif defined(__APPLE__)
		if((state->outspace[i] = (float *)_mm_malloc(sizeof(float) * fftlen, 16)) == NULL)
#else
		if((state->outspace[i] = (float *)calloc(sizeof(float), fftlen)) == NULL)
#endif
			goto error;
#if defined(USE_FFTW) || defined(__APPLE__)
		memset(state->outspace[i], 0, sizeof(float) * fftlen);
#endif
	}

	if((state->inspace = (float **)calloc(sizeof(float *), input_channels)) == NULL)
		goto error;
	for(i = 0; i < input_channels; ++i) {
#ifdef USE_FFTW
		if((state->inspace[i] = (float *)fftwf_malloc(sizeof(float) * fftlen)) == NULL)
#elif defined(__APPLE__)
		if((state->inspace[i] = (float *)_mm_malloc(sizeof(float) * fftlen, 16)) == NULL)
#else
		if((state->inspace[i] = (float *)calloc(sizeof(float), fftlen)) == NULL)
#endif
			goto error;
#if defined(USE_FFTW) || defined(__APPLE__)
		memset(state->inspace[i], 0, sizeof(float) * fftlen);
#endif
	}

#ifdef USE_FFTW
	if((state->p_fw = (fftwf_plan *)calloc(sizeof(fftwf_plan), input_channels)) == NULL)
		goto error;

	for(i = 0; i < input_channels; ++i) {
		if((state->p_fw[i] = fftwf_plan_dft_r2c_1d(fftlen, state->inspace[i], state->f_in, FFTW_ESTIMATE)) == NULL)
			goto error;
	}
	if((state->p_bw = fftwf_plan_dft_c2r_1d(fftlen, state->f_out, state->revspace, FFTW_ESTIMATE)) == NULL)
		goto error;
#elif defined(__APPLE__)
	if((state->setup = vDSP_create_fftsetup(state->fftlenlog2, FFT_RADIX2)) == NULL)
		goto error;
#else
	if((state->cfg_fw = kiss_fftr_alloc(fftlen, 0, NULL, NULL)) == NULL)
		goto error;
	if((state->cfg_bw = kiss_fftr_alloc(fftlen, 1, NULL, NULL)) == NULL)
		goto error;
#endif

	convolver_restage(state, impulse);

	return state;

error:
	convolver_delete(state);
	return NULL;
}

/* Restage the convolver with a new impulse set, same size/parameters */
void convolver_restage(void *state_, const float *const *impulse) {
	convolver_state *state = (convolver_state *)state_;

	float *impulse_temp;

	int impulse_count;
	int channels_per_impulse;
	int fftlen = state->fftlen;
	int impulse_size = state->impulselen;
	int i, j, k;
	int lenover2 = state->fftlenover2;

#ifdef USE_FFTW
	fftwf_plan p;
	fftwf_complex **f_ir = state->f_ir;
#elif defined(__APPLE__)
	int log2n = state->fftlenlog2;
	FFTSetup setup = state->setup;
	DSPSplitComplex *f_ir = state->f_ir;
#else
	kiss_fftr_cfg cfg_fw = state->cfg_fw;
	kiss_fft_cpx **f_ir = state->f_ir;
#endif

	if(state->mode == 0 || state->mode == 1)
		impulse_count = 1;
	else if(state->mode == 2)
		impulse_count = state->inputs;

	if(state->mode == 0)
		channels_per_impulse = 1;
	else if(state->mode == 1 || state->mode == 2)
		channels_per_impulse = state->outputs;

	/* Since the FFT requires a full input for every transformaton, we allocate
	 * a temporary buffer, which we fill with the impulse, then pad with silence. */

	if((impulse_temp = (float *)malloc(sizeof(float) * fftlen)) == NULL)
		return;

	memset(impulse_temp + impulse_size, 0, sizeof(float) * (fftlen - impulse_size));

	for(i = 0; i < impulse_count; ++i) {
		for(j = 0; j < channels_per_impulse; ++j) {
			for(k = 0; k < impulse_size; ++k) {
				impulse_temp[k] = impulse[i][j + k * channels_per_impulse];
			}

			/* Our first actual transformation, which is cached for the life of this convolver. */
#ifdef USE_FFTW
			p = fftwf_plan_dft_r2c_1d(fftlen, impulse_temp, f_ir[i * channels_per_impulse + j], FFTW_ESTIMATE);
			if(p) {
				fftwf_execute(p);
				fftwf_destroy_plan(p);
			}
#elif defined(__APPLE__)
			vDSP_ctoz((DSPComplex *)impulse_temp, 2, &f_ir[i * channels_per_impulse + j], 1, lenover2);
			vDSP_fft_zrip(setup, &f_ir[i * channels_per_impulse + j], 1, log2n, FFT_FORWARD);
#else
			kiss_fftr(cfg_fw, impulse_temp, f_ir[i * channels_per_impulse + j]);
#endif
		}
	}

	free(impulse_temp);
}

/* Delete our opaque state, by freeing all of its member structures, then the
 * top level structure itself. */

void convolver_delete(void *state_) {
	if(state_) {
		int i, input_channels, output_channels, total_channels;
		convolver_state *state = (convolver_state *)state_;
		input_channels = state->inputs;
		output_channels = state->outputs;
		if(state->mode == 0)
			total_channels = 1;
		else if(state->mode == 1)
			total_channels = input_channels;
		else if(state->mode == 2)
			total_channels = input_channels * output_channels;
#ifdef USE_FFTW
		if(state->p_fw) {
			for(i = 0; i < input_channels; ++i) {
				if(state->p_fw[i])
					fftwf_destroy_plan(state->p_fw[i]);
			}
			free(state->p_fw);
		}
		if(state->p_bw)
			fftwf_destroy_plan(state->p_bw);
#elif defined(__APPLE__)
		if(state->setup)
			vDSP_destroy_fftsetup(state->setup);
#else
		if(state->cfg_fw)
			kiss_fftr_free(state->cfg_fw);
		if(state->cfg_bw)
			kiss_fftr_free(state->cfg_bw);
#endif
		if(state->f_ir) {
			for(i = 0; i < total_channels; ++i) {
#if !defined(USE_FFTW) && defined(__APPLE__)
				_free_dspsplitcomplex(&state->f_ir[i]);
#else
				if(state->f_ir[i])
#ifdef USE_FFTW
					fftwf_free(state->f_ir[i]);
#else
					KISS_FFT_FREE(state->f_ir[i]);
#endif
#endif
			}
			free(state->f_ir);
		}
#if !defined(USE_FFTW) && defined(__APPLE__)
		_free_dspsplitcomplex(&state->f_out);
#else
		if(state->f_out)
#ifdef USE_FFTW
			fftwf_free(state->f_out);
#else
			KISS_FFT_FREE(state->f_out);
#endif
#endif
#if !defined(USE_FFTW) && defined(__APPLE__)
		_free_dspsplitcomplex(&state->f_in);
#else
		if(state->f_in)
#ifdef USE_FFTW
			fftwf_free(state->f_in);
#else
			KISS_FFT_FREE(state->f_in);
#endif
#endif
		if(state->revspace)
#ifdef USE_FFTW
			fftwf_free(state->revspace);
#elif defined(__APPLE__)
			_mm_free(state->revspace);
#else
			KISS_FFT_FREE(state->revspace);
#endif
		if(state->outspace) {
			for(i = 0; i < output_channels; ++i) {
				if(state->outspace[i])
#ifdef USE_FFTW
					fftwf_free(state->outspace[i]);
#elif defined(__APPLE__)
					_mm_free(state->outspace[i]);
#else
					free(state->outspace[i]);
#endif
			}
			free(state->outspace);
		}
		if(state->inspace) {
			for(i = 0; i < input_channels; ++i) {
				if(state->inspace[i])
#ifdef USE_FFTW
					fftwf_free(state->inspace[i]);
#elif defined(__APPLE__)
					_mm_free(state->inspace[i]);
#else
					free(state->inspace[i]);
#endif
			}
			free(state->inspace);
		}
		free(state);
	}
}

/* This resets the state between uses, if you need to restart output on startup. */

void convolver_clear(void *state_) {
	if(state_) {
		/* Clearing for a new use setup only requires resetting the input and
		 * output buffers, not actually changing any of the FFT state. */

		int i, input_channels, output_channels, fftlen;
		convolver_state *state = (convolver_state *)state_;
		input_channels = state->inputs;
		output_channels = state->outputs;
		fftlen = state->fftlen;
		state->buffered_in = 0;
		state->buffered_out = 0;
		for(i = 0; i < input_channels; ++i)
			memset(state->inspace[i], 0, sizeof(float) * fftlen);
		for(i = 0; i < output_channels; ++i)
			memset(state->outspace[i], 0, sizeof(float) * fftlen);
	}
}

/* Input sample data is fed in here, one sample at a time. */

static void convolver_write(void *state_, const float *input_samples, int count) {
	if(state_) {
		convolver_state *state = (convolver_state *)state_;

		int i, j, k, input_channels;
		input_channels = state->inputs;

		for(j = 0; j < count; ++j) {
			for(i = 0; i < input_channels; ++i)
				state->inspace[i][state->buffered_in] = input_samples[i];
			input_samples += input_channels;
			++state->buffered_in;
		}

		for(j = 0; j < input_channels; ++j) {
			memset(&state->inspace[j][state->buffered_in], 0, (state->fftlen - state->buffered_in) * sizeof(float));
		}

		/* And every stepsize samples buffered, it convolves a new block of samples. */

		{
			int output_channels = state->outputs;
			int fftlen;
			int index;
			float fftlen_if;
			int lenover2 = state->fftlenover2;
#ifdef USE_FFTW
			fftwf_complex *f_in = state->f_in;
			fftwf_complex *f_out = state->f_out;
#elif defined(__APPLE__)
			DSPSplitComplex *f_in = &state->f_in;
			DSPSplitComplex *f_out = &state->f_out;
			FFTSetup setup = state->setup;
			int log2n = state->fftlenlog2;
			float scale;
#else
			kiss_fft_cpx *f_in = state->f_in;
			kiss_fft_cpx *f_out = state->f_out;
#endif
			float *revspace = state->revspace;

#if !defined(USE_FFTW) && defined(__APPLE__)
			fftlen = state->fftlen;
#endif

			if(state->mode == 0 || state->mode == 1) {
				for(i = 0; i < input_channels; ++i) {
					int index = i * state->mode;
#ifdef USE_FFTW
					fftwf_complex *f_ir = state->f_ir[index];
#elif defined(__APPLE__)
					DSPSplitComplex *f_ir = &state->f_ir[index];
					float preserveIRNyq;
					float preserveSigNyq;
#else
					kiss_fft_cpx *f_ir = state->f_ir[index];
#endif
					float *outspace;
#if !defined(USE_FFTW) && defined(__APPLE__)
					vDSP_ctoz((DSPComplex *)(state->inspace[i]), 2, f_in, 1, lenover2);

					vDSP_fft_zrip(setup, f_in, 1, log2n, FFT_FORWARD);

					preserveIRNyq = f_ir->imagp[0];
					f_ir->imagp[0] = 0;
					preserveSigNyq = f_in->imagp[0];
					f_in->imagp[0] = 0;

					vDSP_zvmul(f_in, 1, f_ir, 1, f_out, 1, lenover2, 1);

					f_out->imagp[0] = preserveIRNyq * preserveSigNyq;
					f_ir->imagp[0] = preserveIRNyq;

					vDSP_fft_zrip(setup, f_out, 1, log2n, FFT_INVERSE);

					vDSP_ztoc(f_out, 1, (DSPComplex *)revspace, 2, lenover2);

					scale = 1.0 / (4.0 * (float)fftlen);
					vDSP_vsmul(revspace, 1, &scale, revspace, 1, fftlen);

					outspace = state->outspace[i];
					vDSP_vadd(revspace, 1, outspace, 1, outspace, 1, fftlen);
#else
#ifdef USE_FFTW
					fftwf_execute(state->p_fw[i]);
#else
					kiss_fftr(state->cfg_fw, state->inspace[i], f_in);
#endif

					for(k = 0, fftlen = lenover2; k < fftlen; ++k) {
#ifdef USE_FFTW
						float re = f_ir[k][0] * f_in[k][0] - f_ir[k][1] * f_in[k][1];
						float im = f_ir[k][1] * f_in[k][0] + f_ir[k][0] * f_in[k][1];
						f_out[k][0] = re;
						f_out[k][1] = im;
#else
						float re = f_ir[k].r * f_in[k].r - f_ir[k].i * f_in[k].i;
						float im = f_ir[k].i * f_in[k].r + f_ir[k].r * f_in[k].i;
						f_out[k].r = re;
						f_out[k].i = im;
#endif
					}

#ifdef USE_FFTW
					fftwf_execute(state->p_bw);
#else
					kiss_fftri(state->cfg_bw, f_out, revspace);
#endif

					outspace = state->outspace[i];
					for(j = 0, fftlen = state->fftlen, fftlen_if = 1.0f / (float)fftlen; j < fftlen; ++j)
						outspace[j] += revspace[j] * fftlen_if;
#endif
				}
			} else if(state->mode == 2) {
				for(i = 0; i < input_channels; ++i) {
					/* First the input samples are transformed to frequency domain, like
					 * the cached impulse was in the setup function. */

#ifdef USE_FFTW
					fftwf_execute(state->p_fw[i]);
#elif defined(__APPLE__)
					float preserveIRNyq;
					float preserveSigNyq;

					vDSP_ctoz((DSPComplex *)(state->inspace[i]), 2, f_in, 1, lenover2);

					vDSP_fft_zrip(setup, f_in, 1, log2n, FFT_FORWARD);

					preserveSigNyq = f_in->imagp[0];
					f_in->imagp[0] = 0;
#else
					kiss_fftr(state->cfg_fw, state->inspace[i], f_in);
#endif

					for(j = 0; j < output_channels; ++j) {
						/* Then we cross multiply the products of the frequency domain, the
						 * real and imaginary values, into output real and imaginary pairs. */

						int index = i * output_channels + j;
#ifdef USE_FFTW
						fftwf_complex *f_ir = state->f_ir[index];
#elif defined(__APPLE__)
						DSPSplitComplex *f_ir = &state->f_ir[index];
						float scale;
#else
						kiss_fft_cpx *f_ir = state->f_ir[index];
#endif
						float *outspace;
#if !defined(USE_FFTW) && defined(__APPLE__)
						preserveIRNyq = f_ir->imagp[0];
						f_ir->imagp[0] = 0;

						vDSP_zvmul(f_in, 1, f_ir, 1, f_out, 1, lenover2, 1);

						f_ir->imagp[0] = preserveIRNyq;
						f_out->imagp[0] = preserveIRNyq * preserveSigNyq;
#else
						for(k = 0, fftlen = lenover2; k < fftlen; ++k) {
#ifdef USE_FFTW
							float re = f_ir[k][0] * f_in[k][0] - f_ir[k][1] * f_in[k][1];
							float im = f_ir[k][1] * f_in[k][0] + f_ir[k][0] * f_in[k][1];
							f_out[k][0] = re;
							f_out[k][1] = im;
#else
							float re = f_ir[k].r * f_in[k].r - f_ir[k].i * f_in[k].i;
							float im = f_ir[k].i * f_in[k].r + f_ir[k].r * f_in[k].i;
							f_out[k].r = re;
							f_out[k].i = im;
#endif
						}
#endif

						/* Then we transform back from frequency to time domain. */

#ifdef USE_FFTW
						fftwf_execute(state->p_bw);
#elif defined(__APPLE__)
						vDSP_fft_zrip(setup, f_out, 1, log2n, FFT_INVERSE);
						vDSP_ztoc(f_out, 1, (DSPComplex *)revspace, 2, lenover2);
#else
						kiss_fftri(state->cfg_bw, f_out, revspace);
#endif

						/* Then we add the entire revspace block onto our output, dividing
						 * each value by the total number of samples in the buffer. Remember,
						 * since there is some overlap, this addition step is important. */

#if !defined(USE_FFTW) && defined(__APPLE__)
						scale = 1.0 / (4.0 * (float)fftlen);
						vDSP_vsmul(revspace, 1, &scale, revspace, 1, fftlen);

						outspace = state->outspace[j];
						vDSP_vadd(revspace, 1, outspace, 1, outspace, 1, fftlen);
#else
						outspace = state->outspace[j];
						for(k = 0, fftlen = state->fftlen, fftlen_if = 1.0f / (float)fftlen; k < fftlen; ++k)
							outspace[k] += revspace[k] * fftlen_if;
#endif
					}
				}
			}

			/* Output samples are now buffered and ready for retrieval. */

			state->buffered_out = state->buffered_in;
			state->buffered_in = 0;
		}
	}
}

/* Call this to process samples */

void convolver_run(void *state_, const float *input_samples, float *output_samples, int count) {
	if(state_) {
		convolver_state *state = (convolver_state *)state_;

		while(count > 0) {
			int count_to_do = count;
			if(count_to_do > state->stepsize)
				count_to_do = state->stepsize;

			convolver_write(state_, input_samples, count_to_do);

			input_samples += count_to_do * state->inputs;

			int i, j, output_channels = state->outputs;

			for(j = 0; j < count_to_do; ++j) {
				for(i = 0; i < output_channels; ++i) {
					float sample = state->outspace[i][j];

					output_samples[i] = sample;
				}

				output_samples += output_channels;
			}

			for(i = 0; i < output_channels; ++i) {
				float *outspace = state->outspace[i];
				memmove(outspace, outspace + state->buffered_out, (state->fftlen - state->buffered_out) * sizeof(float));
				memset(outspace + state->fftlen - state->buffered_out, 0, state->buffered_out * sizeof(float));
			}

			count -= state->buffered_out;
			state->buffered_out = 0;
		}
	}
}
