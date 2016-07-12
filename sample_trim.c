#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _countof(d) (sizeof((d))/sizeof(((d)[0])))

static const int actual_frequencies[] = { 8000, 11025, 16000, 22050, 32000, 44100, 48000 };
static const int frequencies[] = { 8, 11, 16, 22, 32, 441, 48 };
static const int frequency_count = _countof(frequencies);
static const char* speakers[] = { "FL", "FR", "FC", "LFE", "BL", "BR" };
static const int speaker_count = _countof(speakers);

static int sample_counts[frequency_count];

static const int impulse_wav_size = 524332;

unsigned int get_le32(const unsigned char * ptr)
{
	return ptr[0] + (ptr[1] << 8) + (ptr[2] << 16) + (ptr[3] << 24);
}

unsigned int get_be32(const unsigned char * ptr)
{
	return ptr[3] + (ptr[2] << 8) + (ptr[1] << 16) + (ptr[0] << 24);
}

unsigned int filter_sample(const unsigned char * ptr)
{
	unsigned int sample = get_le32(ptr);
	if (sample == 0x80000000) sample = 0;
	return sample;
}

int find_data(const unsigned char * buffer, unsigned int size, int * data_offset, int * data_size)
{
	const unsigned char * ptr = buffer;
	unsigned int tag;
	unsigned int size_;
	if (size < 8) return -1;
	tag = get_be32(ptr);
	if (tag != 'RIFF') return -1;
	size_ = get_le32(ptr + 4);
	if (size_ + 8 > size) return -1;
	size = size_;
	ptr += 8;
	if (size < 4) return -1;
	tag = get_be32(ptr);
	if (tag != 'WAVE') return -1;
	size -= 4;
	ptr += 4;
	do {
		tag = get_be32(ptr);
		size_ = get_le32(ptr + 4);
		if (size_ + 8 > size) return -1;
		if (tag == 'data') {
			*data_offset = (ptr - buffer) + 8;
			*data_size = size_;
			return 0;
		}
		if (size_ & 1) size_ += 1;
		ptr += size_ + 8;
		size -= size_ + 8;
	} while (size > 0);
	return -1;
}

int signum(float value)
{
	return (value > 0) - (value < 0);
}

void fprint_float(FILE *f, float value)
{
	int sign = signum(value);
	float fraction;
	int whole;
	value = fabs(value);
	fraction = fmod(value, 1.0);
	whole = (int)(value -= fraction);
	if (sign < 0) fputc('-', f);
	fprintf(f, "%u", whole);
	if (fraction)
	{
		fputc('.', f);
		do
		{
			value = fraction * 10.0;
			fraction = fmod(value, 1.0);
			whole = (int)(value -= fraction);
			fputc('0' + whole, f);
		}
		while (fraction);
	}
}

int main(void)
{
	int level, frequency, speaker, sample;
	int min_sample, max_sample, sample_count;
	int data_offset, data_size;
	FILE *f;
	char name[128];
	unsigned char * buffer[speaker_count];

	fprintf(stdout, "typedef struct speaker_impulses\n{\n\tunsigned int count;\n\tconst float * impulse[%u];\n} speaker_impulses;\n\n", speaker_count);

	for (speaker = 0; speaker < speaker_count; ++speaker)
	{
		buffer[speaker] = (unsigned char *) malloc(impulse_wav_size);
	}

	for (level = 1; level <= 3; ++level)
	{
		for (frequency = 0; frequency < frequency_count; ++frequency)
		{
			min_sample = 0x7fffffff;
			max_sample = 0;
			for (speaker = 0; speaker < speaker_count; ++speaker)
			{
				sprintf(name, "samples/processed/sample_%u_%s_dh%u.wav", frequencies[frequency], speakers[speaker], level);
				f = fopen(name, "rb");
				fread(buffer[speaker], 1, impulse_wav_size, f);
				fclose(f);

				if (find_data(buffer[speaker], impulse_wav_size, &data_offset, &data_size) < 0)
				{
					fprintf(stderr, "Invalid sample: %s\n", name);
					return 1;
				}

				for (sample = 0, sample_count = data_size / 8; sample < sample_count; ++sample)
				{
					unsigned int sample1 = filter_sample(buffer[speaker] + data_offset + sample * 8);
					unsigned int sample2 = filter_sample(buffer[speaker] + data_offset + sample * 8 + 4);
					if (sample1 || sample2) break;
				}

				if (sample < min_sample) min_sample = sample;

				for (sample = sample_count - 1; sample >= 0; --sample)
				{
					unsigned int sample1 = filter_sample(buffer[speaker] + data_offset + sample * 8);
					unsigned int sample2 = filter_sample(buffer[speaker] + data_offset + sample * 8 + 4);
					if (sample1 || sample2) break;
				}

				if (sample > max_sample) max_sample = sample;
			}

			sample_counts[frequency] = max_sample - min_sample + 1;

			for (speaker = 0; speaker < speaker_count; ++speaker)
			{
				find_data(buffer[speaker], impulse_wav_size, &data_offset, &data_size);
				fprintf(stdout, "static const unsigned int impulse_l%u_%u_%s[%u * 2] = {\n", level, frequencies[frequency], speakers[speaker], max_sample - min_sample + 1);
				for (sample = min_sample; sample <= max_sample; ++sample)
				{
					unsigned int i;
					if (((sample - min_sample) & 3) == 0) fprintf(stdout, "\t");
					i = filter_sample(buffer[speaker] + data_offset + sample * 8);
					fprintf(stdout, "0x%08x, ", i);
					i = filter_sample(buffer[speaker] + data_offset + sample * 8 + 4);
					fprintf(stdout, "0x%08x, ", i);
					if (((sample - min_sample) & 3) == 3) fprintf(stdout, "\n");
				}
				if ((max_sample - min_sample + 1) & 3) fprintf(stdout, "\n");
				fprintf(stdout, "};\n\n");
			}
		}
	}

	for (level = 1; level <= 3; ++level)
	{
		for (frequency = 0; frequency < frequency_count; ++frequency)
		{
			fprintf(stdout, "static const speaker_impulses impulses_l%u_%u = {\n\t%u,\n\t{\n", level, frequencies[frequency], sample_counts[frequency]);

			for (speaker = 0; speaker < speaker_count; ++speaker)
			{
				fprintf(stdout, "\t\t(const float *)&impulse_l%u_%u_%s%s\n", level, frequencies[frequency], speakers[speaker], speaker < speaker_count - 1 ? "," : "");
			}

			fprintf(stdout, "\t}\n};\n\n");
		}
	}

	fprintf(stdout, "typedef struct speaker_preset\n{\n\tunsigned int frequency;\n\tconst speaker_impulses *impulses;\n} speaker_preset;\n\n");

	fprintf(stdout, "static const int speaker_preset_count = %u;\n\n", frequency_count);

	fprintf(stdout, "static const speaker_preset speaker_presets[3][%u] =\n{\n", frequency_count);

	for (level = 1; level <= 3; ++level)
	{
		fprintf(stdout, "\t{\n");
		for (frequency = frequency_count - 1; frequency >= 0; --frequency)
		{
			fprintf(stdout, "\t\t{ %u, &impulses_l%u_%u }%s\n", actual_frequencies[frequency], level, frequencies[frequency], frequency > 0 ? "," : "");
		}
		fprintf(stdout, "\t}%s\n", level < 3 ? "," : "");
	}

	fprintf(stdout, "};\n");

	return 0;
}
