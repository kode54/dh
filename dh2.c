#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simple_convolver.h"

#include "samples.h"

unsigned int get_le16(const unsigned char *ptr) {
	return ptr[0] + (ptr[1] << 8);
}

unsigned int get_le32(const unsigned char *ptr) {
	return ptr[0] + (ptr[1] << 8) + (ptr[2] << 16) + (ptr[3] << 24);
}

unsigned int get_be32(const unsigned char *ptr) {
	return ptr[3] + (ptr[2] << 8) + (ptr[1] << 16) + (ptr[0] << 24);
}

int main(int argc, char **argv) {
	FILE *f, *g;
	unsigned char buffer[1024];
	unsigned int riff_size, fmt_size, data_offset, data_size;
	unsigned int id;
	unsigned int sample_rate, sample_count, samples_out;
	unsigned int sample, preset, samples_written;
	unsigned int i;

	const speaker_preset *set;

	void *conv;

	float inbuffer[1024 * 6];
	float outbuffer[1024 * 2];

	if(argc != 3) {
		fprintf(stderr, "Usage:\tdh2 <input.wav> <output.raw>\n");
		return 1;
	}

	f = fopen(argv[1], "rb");

	if(!f) {
		fprintf(stderr, "Unable to open %s.\n", argv[1]);
		return 1;
	}

	if(fread(buffer, 1, 8, f) != 8) {
		fclose(f);
		fprintf(stderr, "Unable to read WAV header.\n");
		return 1;
	}

	if(get_be32(buffer) != 'RIFF') {
		fclose(f);
		fprintf(stderr, "Not a RIFF file.\n");
		return 1;
	}

	riff_size = get_le32(buffer + 4);

	if(riff_size < 4) {
		fclose(f);
		fprintf(stderr, "RIFF too small.\n");
		return 1;
	}

	fread(buffer, 1, 4, f);
	riff_size -= 4;

	if(get_be32(buffer) != 'WAVE') {
		fclose(f);
		fprintf(stderr, "Not WAVE format.\n");
		return 1;
	}

	fmt_size = 0;
	data_size = 0;

	while(riff_size >= 8) {
		fread(buffer, 1, 8, f);
		riff_size -= 8;

		id = get_be32(buffer);
		if(id == 'fmt ') {
			if(fmt_size) {
				fclose(f);
				fprintf(stderr, "Multiple fmt chunks found.\n");
				return 1;
			}

			fmt_size = get_le32(buffer + 4);
			if(fmt_size & 1) ++fmt_size;
			fread(buffer, 1, fmt_size, f);
			riff_size -= fmt_size;

			sample_rate = get_le32(buffer + 4);
		} else if(id == 'data') {
			if(data_size) {
				fclose(f);
				fprintf(stderr, "Multiple data chunks found.\n");
				return 1;
			}

			data_size = get_le32(buffer + 4);
			if(data_size & 1) data_size++;
			riff_size -= data_size;
			data_offset = ftell(f);
			fseek(f, data_size, SEEK_CUR);
		} else {
			int size = get_le32(buffer + 4);
			if(size & 1) size++;
			riff_size -= size;
			fseek(f, size, SEEK_CUR);
		}
	}

	if(!fmt_size || !data_size) {
		fclose(f);
		if(!fmt_size) {
			fprintf(stderr, "Missing fmt chunk.\n");
		}
		if(!data_size) {
			fprintf(stderr, "Missing data chunk.\n");
		}
		return 1;
	}

	g = fopen(argv[2], "wb");
	if(!g) {
		fclose(f);
		fprintf(stderr, "Unable to open %s for writing.\n", argv[2]);
		return 1;
	}

	fseek(f, data_offset, SEEK_SET);

	sample_count = data_size / (4 * 6);

	set = speaker_presets[1];

	for(preset = 0; preset < speaker_preset_count; ++preset) {
		if(set[preset].frequency == sample_rate) break;
	}

	conv = convolver_create(set[preset].impulses->impulse, set[preset].impulses->count, 6, 2, 2);

	sample = 0;
	samples_written = 0;

	while(samples_written < sample_count) {
		size_t samples_in = 0;
		while(samples_in < 1024 && sample < sample_count) {
			size_t samples_to_read = sample_count - sample;
			if(samples_to_read > 1024 - samples_in)
				samples_to_read = 1024 - samples_in;
			fread(inbuffer + samples_in * 6, 4 * 6, samples_to_read, f);
			sample += samples_to_read;
			samples_in += samples_to_read;
		}

		convolver_run(conv, inbuffer, outbuffer, samples_in);

		fwrite(outbuffer, 2 * 4, samples_in, g);

		samples_written += samples_in;
	}

	fclose(g);
	fclose(f);

	return 0;
}
