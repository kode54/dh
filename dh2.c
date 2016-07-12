#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "simple_convolver.h"

#include "samples.h"

unsigned int get_le16(const unsigned char * ptr)
{
	return ptr[0] + (ptr[1] << 8);
}

unsigned int get_le32(const unsigned char * ptr)
{
	return ptr[0] + (ptr[1] << 8) + (ptr[2] << 16) + (ptr[3] << 24);
}

unsigned int get_be32(const unsigned char * ptr)
{
	return ptr[3] + (ptr[2] << 8) + (ptr[1] << 16) + (ptr[0] << 24);
}

int main(int argc, char ** argv)
{
	FILE *f, *g;
	unsigned char buffer[1024];
	unsigned int riff_size, fmt_size, data_offset, data_size;
	unsigned int id;
	unsigned int sample_rate, sample_count, samples_out;
	unsigned int sample, preset, samples_written;
	unsigned int i;

	const speaker_preset * set;

	void *conv;

	float samples[6];

	if (argc != 3)
	{
		fprintf(stderr, "Usage:\tdh2 <input.wav> <output.raw>\n");
		return 1;
	}

	f = fopen(argv[1], "rb");

	if ( !f )
	{
		fprintf(stderr, "Unable to open %s.\n", argv[1]);
		return 1;
	}

	if ( fread(buffer, 1, 8, f) != 8 )
	{
		fclose(f);
		fprintf(stderr, "Unable to read WAV header.\n");
		return 1;
	}

	if ( get_be32( buffer ) != 'RIFF' )
	{
		fclose(f);
		fprintf(stderr, "Not a RIFF file.\n");
		return 1;
	}

	riff_size = get_le32( buffer + 4 );

	if ( riff_size < 4 )
	{
		fclose(f);
		fprintf(stderr, "RIFF too small.\n");
		return 1;
	}

	fread(buffer, 1, 4, f);
	riff_size -= 4;

	if ( get_be32( buffer ) != 'WAVE' )
	{
		fclose(f);
		fprintf(stderr, "Not WAVE format.\n");
		return 1;
	}

	fmt_size = 0;
	data_size = 0;

	while ( riff_size >= 8 )
	{
		fread(buffer, 1, 8, f);
		riff_size -= 8;

		id = get_be32(buffer);
		if (id == 'fmt ')
		{
			if (fmt_size)
			{
				fclose(f);
				fprintf(stderr, "Multiple fmt chunks found.\n");
				return 1;
			}

			fmt_size = get_le32(buffer + 4);
			if (fmt_size & 1) ++fmt_size;
			fread(buffer, 1, fmt_size, f);
			riff_size -= fmt_size;

			sample_rate = get_le32(buffer + 4);
		}
		else if (id == 'data')
		{
			if (data_size)
			{
				fclose(f);
				fprintf(stderr, "Multiple data chunks found.\n");
				return 1;
			}

			data_size = get_le32(buffer + 4);
			if (data_size & 1) data_size++;
			riff_size -= data_size;
			data_offset = ftell(f);
			fseek(f, data_size, SEEK_CUR);
		}
		else
		{
			int size = get_le32(buffer + 4);
			if (size & 1) size++;
			riff_size -= size;
			fseek(f, size, SEEK_CUR);
		}
	}

	if (!fmt_size || !data_size)
	{
		fclose(f);
		if (!fmt_size)
		{
			fprintf(stderr, "Missing fmt chunk.\n");
		}
		if (!data_size)
		{
			fprintf(stderr, "Missing data chunk.\n");
		}
		return 1;
	}

	g = fopen(argv[2], "wb");
	if (!g)
	{
		fclose(f);
		fprintf(stderr, "Unable to open %s for writing.\n", argv[2]);
		return 1;
	}

	fseek(f, data_offset, SEEK_SET);

	sample_count = data_size / (4 * 6);

	set = speaker_presets[1];

	for (preset = 0; preset < speaker_preset_count; ++preset)
	{
		if (set[preset].frequency == sample_rate) break;
	}

	conv = convolver_create(set[preset].impulses->impulse, set[preset].impulses->count, 6, 2, 2);

	sample = 0;
	samples_written = 0;

	while (samples_written < sample_count)
	{
		unsigned int samples_free;
		unsigned int samples_ready = convolver_ready(conv);
		while (samples_ready)
		{
			convolver_read(conv, samples);
			fwrite(samples, 4, 2, g);
			++samples_written;
			--samples_ready;
			if (samples_written == sample_count) break;
		}
		if (samples_written == sample_count) break;
		samples_free = convolver_get_free_count(conv);
		while (samples_free)
		{
			if (sample < sample_count)
			{
				fread(samples, 4, 6, f);
				++sample;
			}
			else
			{
				memset(samples, 0, 6*4);
			}
			convolver_write(conv, samples);
			--samples_free;
		}
	}

	fclose(g);
	fclose(f);

	return 0;
}
