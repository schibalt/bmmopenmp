// BlockedMatrixMultiply.cpp : Defines the entry point for the console application.

#define BUFSIZE 26
//#define ARRAY_SIZE 128
#define __STDC_FORMAT_MACROS

#include <inttypes.h>
#include <pthread.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdint-gcc.h>
#include <vector>
#include <fstream>
#include <sys/time.h>
#include <ctgmath>

#ifdef WIN32
#include <Windows.h>
#endif

using std::min;
using std::stringstream;
using std::ofstream;
using std::ifstream;
using std::string;
using std::vector;
using std::endl;
using std::pow;

int ARRAY_SIZE;

uint64_t GetTimeMs64() {

#ifdef WIN32
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed
	 * since January 1, 1601 (UTC) and copy it
	 * to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	uint64_t ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10; // From 100 nano seconds (10^-7)

	return ret;
#else
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	uint64_t ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
//	ret /= 1000;
	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);

	return ret;
#endif
}

vector<string> &split(const string &s, char delim, vector<string> &elems) {

	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {

		elems.push_back(item);
	}
	return elems;
}

vector<string> split(const string &s, char delim) {

	vector<string> elems;
	split(s, delim, elems);
	return elems;
}

vector<string> instvector(ifstream* matrixfile, string matrixfilename) {

	string firstline, line;
	bool print = false;
	//process matrix a
	if ((*matrixfile).is_open()) {

		getline(*matrixfile, firstline);
		while (getline(*matrixfile, line))
			printf("an extra line was found in %s\n", matrixfilename.c_str());
		(*matrixfile).close();
	} else
		printf("Unable to open %s\n", matrixfilename.c_str());

	vector<string> matrixvector = split(firstline, ' ');

	int elementcounter = 0;

	for (vector<string>::iterator it = matrixvector.begin();
			it != matrixvector.end(); ++it) {
		elementcounter++;

		if (print)
			printf("%i\t%s\n", elementcounter, (*it).c_str());
	}
	if (print)
		printf("\n");
	return matrixvector;
	//process matrix a
}

void printmatrix(int arraysize, float** matrix, string name) {

	printf("\nprinting %d x %d %s matrix\n\n", arraysize, arraysize,
			name.c_str());
	for (int x = 0; x < arraysize; x++) {

		for (int y = 0; y < arraysize; y++)
			printf("%f ", matrix[x][y]);
		printf("\n");
	}
}

int main(int argc, char* argv[]) {

	string matrixafilename = argv[2];
	string matrixbfilename = argv[3];
	string answermatrixfilename = argv[4];

	ifstream matrixafile(matrixafilename.c_str());
	ifstream matrixbfile(matrixbfilename.c_str());
	ifstream answermatrixfile(answermatrixfilename.c_str());

	string firstline, line;
	vector<string> matrixavector = instvector(&matrixafile, matrixafilename);
	vector<string> matrixbvector = instvector(&matrixbfile, matrixbfilename);
	vector<string> answermatrixvector = instvector(&answermatrixfile,
			answermatrixfilename);

	int nelements = matrixavector.size();
	ARRAY_SIZE = sqrt(nelements);

	char * pEnd;
	int block_size = strtol(argv[1], &pEnd, 10);

	//inst the arrays
	float** A = new float*[ARRAY_SIZE];
	float** B = new float*[ARRAY_SIZE];
	float** C = new float*[ARRAY_SIZE];
	float** answer = new float*[ARRAY_SIZE];

	for (int i = 0; i < ARRAY_SIZE; ++i) {

		A[i] = new float[ARRAY_SIZE];
		B[i] = new float[ARRAY_SIZE];
		C[i] = new float[ARRAY_SIZE];
		answer[i] = new float[ARRAY_SIZE];

		for (int j = 0; j < ARRAY_SIZE; j++) {

			A[i][j] = atof(matrixavector[i * ARRAY_SIZE + j].c_str());
			B[i][j] = atof(matrixbvector[i * ARRAY_SIZE + j].c_str());
			C[i][j] = 0;
			answer[i][j] = atof(answermatrixvector[i * ARRAY_SIZE + j].c_str());
		}
	}
	// arrays inst

//	printmatrix(ARRAY_SIZE, A, "A");
//	printmatrix(ARRAY_SIZE, B, "B");
//	printmatrix(ARRAY_SIZE, C, "C");
//	printmatrix(ARRAY_SIZE, answer, "answer");

	printf("max %d threads\n", omp_get_max_threads());

	int numthreads;

	if (block_size <= 0)
		block_size = sqrt(omp_get_max_threads());

	numthreads = pow(block_size, 2);
	omp_set_num_threads(numthreads);
	uint64_t start = GetTimeMs64();

	for (int i = 0; i < ARRAY_SIZE; i += block_size) {

		for (int j = 0; j < ARRAY_SIZE; j += block_size) {

#pragma omp parallel
			for (int k = 0; k < ARRAY_SIZE; k++) {
//#pragma omp critical
//								printf("%d threads\n", omp_get_num_threads());

				int myidx = omp_get_thread_num() % block_size;
				int myidy = omp_get_thread_num() / block_size;

//#pragma omp critical
//				printf(
//						"Hello from thread %02i i = %02d j = %02d k = %02d myidx = %02d myidy = %d\n",
//						omp_get_thread_num(), i, j, k, myidx, myidy);
				float result = C[i + myidx][j + myidy]
						+ A[i + myidx][k] * B[k][j + myidy];
#pragma omp critical
				C[i + myidx][j + myidy] += result;
			} //k

		} //j
	} //i

	uint64_t end = GetTimeMs64();
	uint64_t runtime = end - start;

	printmatrix(ARRAY_SIZE, C, "C");
	printmatrix(ARRAY_SIZE, answer, "answer");
	char buffer[50];
	sprintf(buffer,
			"run time %" PRIu64 " μs with block size %i and %i threads\n",
			runtime, block_size, numthreads);
	printf("%s", buffer);
	ofstream outfile;
	outfile.open("bmmompexperiments.txt", ofstream::app);

	outfile << buffer;
	outfile.close();

// De-Allocate memory to prevent memory leak
	for (int i = 0; i < ARRAY_SIZE; ++i) {

		delete[] A[i];
		delete[] B[i];
		delete[] C[i];
		delete[] answer[i];
	}

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] answer;

	return 0;
}
