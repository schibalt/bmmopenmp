// BlockedMatrixMultiply.cpp : Defines the entry point for the console application.

#define BUFSIZE 26
//#define ARRAY_SIZE 128

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

vector<string> split(const string &s, char delim);

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
	ret /= 1000;

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

void printmatrix(int ARRAY_SIZE, float** matrix) {

	for (int x = 0; x < ARRAY_SIZE; x++) {

		for (int y = 0; y < ARRAY_SIZE; y++)
			printf("%f ", matrix[x][y]);
		printf("\n\n");
	}
	printf("\n");
}

int main(int argc, char* argv[]) {

//	ofstream outfile;
//	outfile.open("bmmexperiments.txt");

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
	int ARRAY_SIZE = sqrt(nelements);

	uint64_t starttime = GetTimeMs64();
	char * pEnd;
	long block_size = strtol(argv[1], &pEnd, 10);

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

	printmatrix(ARRAY_SIZE, A);
	printmatrix(ARRAY_SIZE, B);
	printmatrix(ARRAY_SIZE, C);
	printmatrix(ARRAY_SIZE, answer);

	/*
	 printf("\n\n");
	 //this is the normal, working algorithm
	 for (int i = 0; i < ARRAY_SIZE; i++) {
	 for (int j = 0; j < ARRAY_SIZE; j++) {
	 *C[i, j] = 0;
	 for (int k = 0; k < ARRAY_SIZE; k++) {
	 int a = *A[i, k];
	 int b = *B[k, j];
	 int c = *C[i, j];
	 int product = c + (a * b);
	 *C[i, j] = product;
	 }
	 }
	 }*/

	//print
	//				printf("\n\n");
	//				for (int x = 0; x < ARRAY_SIZE; x++) {
	//
	//					for (int y = 0; y < ARRAY_SIZE; y++) {
	//						 if (*C[x] != ARRAY_SIZE)
	//						 exit(EXIT_FAILURE);
	//						 exit(1);
	//						printf("%f", C[x][y]);
	//					}
	//					printf("\n\n");
	//				}
	for (int i = 0; i < ARRAY_SIZE; i += block_size) {

		for (int j = 0; j < ARRAY_SIZE; j += block_size) {

			   omp_set_num_threads(16);
#pragma omp parallel for
			for (int k = 0; k < ARRAY_SIZE; k++) {

				int myidx = omp_get_thread_num() % ARRAY_SIZE;
				int myidy = omp_get_thread_num() / ARRAY_SIZE;
				C[i + myidx][j + myidy] = C[i + myidx][j + myidy]
						+ A[i + myidx][k] * B[k][j + myidy];
			}

//			uint64_t blockstarttime = GetTimeMs64();

			/* old stuff
			 for (int i = 0; i < ARRAY_SIZE; i++) {

			 for (int jj = j; jj < min((int) (j + block_size), ARRAY_SIZE);
			 jj++) {

			 for (int kk = i;
			 kk < min((int) (i + block_size), ARRAY_SIZE); kk++)

			 C[i][jj] = C[i][jj] + (A[i][kk] * B[kk][jj]);
			 }
			 }
			 * 			 */
			/*uint64_t blockendtime = GetTimeMs64();
			 uint64_t blockruntime = blockendtime - blockstarttime;

			 int blocknumber = j / block_size
			 + i / (pow(block_size, 2) / ARRAY_SIZE);
			 stringstream ss;
			 ss << "\tblock " << blocknumber << " runtime " << blockruntime
			 << " microseconds (" << blockruntime / 1000000
			 << " seconds/" << blockruntime / 1000000 / 60
			 << " minutes)\n";

			 string logElement = ss.str();
			 outfile << logElement;}
			 *
			 */
		}
	}

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

	uint64_t endtime = GetTimeMs64();
	uint64_t runtime = endtime - starttime;

//			struct tm now;
	char timebuf[BUFSIZE];
//			size_t err;

// Convert to an ASCII representation.
//			err = strftime(timebuf, BUFSIZE, "%c", &now);

	stringstream ss;
	ss << timebuf << ": \t runtime " << runtime << " microseconds ("
			<< runtime / 1000000 << " seconds/" << runtime / 1000000 / 60
			<< " minutes) for block size " << block_size << "\n";

//	outfile.close();

	return 0;
}
