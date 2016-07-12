CFLAGS = -O2

DH2_OBJS = dh2.o

ST_OBJS = sample_trim.o

CONV_OBJS = simple_convolver.o kissfft/kiss_fft.o kissfft/kiss_fftr.o 

all: dh2

dh2 : $(DH2_OBJS) $(CONV_OBJS)
	$(CC) -o $@ $^

samples.h : sample_trim
	./sample_trim > samples.h

sample_trim : $(ST_OBJS)
	$(CC) -o $@ $^

dh2.o : dh2.c samples.h
	$(CC) -c $(CFLAGS) -o $@ dh2.c

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $*.c

clean:
	rm -f $(DH2_OBJS) $(ST_OBJS) $(CONV_OBJS) dh2 sample_trim samples.h > /dev/null

