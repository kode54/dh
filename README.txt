Dirac to Headphones convolver
-----------------------------

This code and supplied data can translate surround input
to stereo headphone output, using a supplied set of user
captured impulse responses. Currently, it bundles the
secret sauce, but you are welcome to supply your own set
of impulse responses.

Currently supports three FFT libraries:

1) KissFFT, bundled.
2) FFTW 3, if FFTW=1 is passed to Makefile
3) Apple vDSP, the fastest on supported hardware
