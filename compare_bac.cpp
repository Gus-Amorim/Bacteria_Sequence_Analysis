#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <vector> 

int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010
#define NUM_THREADS		4

void Init()
{
	M2 = 1;
	for (int i = 0; i < LEN - 2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
		M2 *= AA_NUMBER;
	M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
	M = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
}

class Bacteria
{
private:
	std::vector<long> second;
	long one_l[AA_NUMBER];
	long indexs;
	long total;
	long total_l;
	long complement;

	void InitVectors()
	{
		vector.resize(M);
		second.resize(M1);
		std::fill(vector.begin(), vector.end(), 0);
		std::fill(second.begin(), second.end(), 0);
		memset(one_l, 0, AA_NUMBER * sizeof(long));
		total = 0;
		total_l = 0;
		complement = 0;
	}

	void init_buffer(char* buffer)
	{
		complement++;
		indexs = 0;
		for (int i=0; i<LEN-1; i++)
		{
			short enc = encode(buffer[i]);
			one_l[enc]++;
			total_l++;
			indexs = indexs * AA_NUMBER + enc;
		}
		second[indexs]++;
	}

	void cont_buffer(char ch)
	{
		short enc = encode(ch);
		one_l[enc]++;
		total_l++;
		long index = indexs * AA_NUMBER + enc;
		vector[index]++;
		total++;
		indexs = (indexs % M2) * AA_NUMBER + enc;
		second[indexs]++;
	}

public:
	std::vector<long> vector;

	Bacteria(const char* filename)
	{
		FILE* bacteria_file;
		errno_t OK = fopen_s(&bacteria_file, filename, "r");

		if (OK != 0)
		{
			fprintf(stderr, "Error: failed to open file %s\n", filename);
			exit(1);
		}

		InitVectors();
		const int BUFFER_SIZE = 200;
		char buffer[BUFFER_SIZE];

		while (fgets(buffer, sizeof(buffer), bacteria_file) != NULL)
		{
			if (buffer[0] == '>')
			{
				continue;
			}
			for (size_t i = 0; i < strlen(buffer); i++)
			{
				if (buffer[i] != '\n')
				{
					this->cont_buffer(buffer[i]);
				}
			}
		}
		fclose(bacteria_file);
	}

	~Bacteria()
	{
		second.clear();
		vector.clear();
	}	

	double CompareBacteria(Bacteria* b2)
	{
		double correlation = 0;
		double vector_len1 = 0;
		double vector_len2 = 0;

		double inv_total_b1 = 1.0 / this->total;
		double inv_total_plus_complement_b1 = 1.0 / (this->total + this->complement);
		double inv_total_l_b1 = 1.0 / this->total_l;

		double inv_total_b2 = 1.0 / b2->total;
		double inv_total_plus_complement_b2 = 1.0 / (b2->total + b2->complement);
		double inv_total_l_b2 = 1.0 / b2->total_l;

		long* second_b1 = &this->second[0];
		long* second_b2 = &b2->second[0];
		long* vector_b1 = &this->vector[0];
		long* vector_b2 = &b2->vector[0];

		for (long i = 0; i < M; i++)
		{
			long i_div_AA_NUMBER = i / AA_NUMBER;
			long i_mod_AA_NUMBER = i % AA_NUMBER;
			long i_div_M1 = i / M1;
			long i_mod_M1 = i % M1;

			double p1_b1 = second_b1[i_div_AA_NUMBER] * inv_total_plus_complement_b1;
			double p2_b1 = this->one_l[i_mod_AA_NUMBER] * inv_total_l_b1;
			double p3_b1 = second_b1[i_mod_M1] * inv_total_plus_complement_b1;
			double p4_b1 = this->one_l[i_div_M1] * inv_total_l_b1;
			double stochastic1 = this->total * (p1_b1 * p2_b1 + p3_b1 * p4_b1) * 0.5;

			double vec_i = vector_b1[i];
			double t1 = (stochastic1 > EPSILON) ? (vec_i - stochastic1) / stochastic1 : 0.0;
			vector_len1 += t1 * t1;

			double p1_b2 = second_b2[i_div_AA_NUMBER] * inv_total_plus_complement_b2;
			double p2_b2 = b2->one_l[i_mod_AA_NUMBER] * inv_total_l_b2;
			double p3_b2 = second_b2[i_mod_M1] * inv_total_plus_complement_b2;
			double p4_b2 = b2->one_l[i_div_M1] * inv_total_l_b2;
			double stochastic2 = b2->total * (p1_b2 * p2_b2 + p3_b2 * p4_b2) * 0.5;

			double t2 = (stochastic2 > EPSILON) ? (vector_b2[i] - stochastic2) / stochastic2 : 0.0;
			vector_len2 += t2 * t2;

			correlation += t1 * t2;
		}
		return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
	}
};

void ReadInputFile(const char* input_name)
{
	FILE* input_file;
	errno_t OK = fopen_s(&input_file, input_name, "r");

	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
		exit(1);
	}

	fscanf_s(input_file, "%d", &number_bacteria);
	bacteria_name = new char* [number_bacteria];

	for (long i = 0; i < number_bacteria; i++)
	{
		char name[10];
		fscanf_s(input_file, "%s", name, 10);
		bacteria_name[i] = new char[20];
		sprintf_s(bacteria_name[i], 20, "data/%s.faa", name);
	}
	fclose(input_file);
}

void CompareAllBacteria()
{
	std::vector<Bacteria*> bacteriaList(number_bacteria);

	for (int i = 0; i < number_bacteria; i++)
	{
		bacteriaList[i] = new Bacteria(bacteria_name[i]);
	}

#pragma omp parallel for schedule(dynamic) num_threads(NUM_THREADS)
	for (int i = 0; i < number_bacteria - 1; i++)
	{
		Bacteria* b1 = bacteriaList[i];

		for (int j = i + 1; j < number_bacteria; j++)
		{
			Bacteria* b2 = bacteriaList[j];
			double correlation = b1->CompareBacteria(b2);
#pragma omp critical
			{
				printf("%03d %03d -> %.10lf\n", i, j, correlation);
			}
		}
	}

	for (int i = 0; i < number_bacteria; i++)
	{
		delete bacteriaList[i];
	}
}

int main(int argc, char* argv[])
{
	time_t t1 = time(NULL);

	Init();
	ReadInputFile("list.txt");
	CompareAllBacteria();

	time_t t2 = time(NULL);
	printf("time elapsed: %lld seconds\n", t2 - t1);
	return 0;
}