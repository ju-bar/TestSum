// TestSum.cpp : This is where you find "main"
// J. Barthel, Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
// Copyright (c) 2018

#include "pch.h"
#include "summation.h"
#include <iostream>
#include <chrono>
#include <string>
#include <vector>

#ifndef __TESTSUM_PRM
#define __TESTSUM_PRM
#define _TESTSUM_NUM			4 // number of test functions
#define _TESTSUM_SCALE_MIN		10000 // min size of the array
#define _TESTSUM_SCALE_STEP		10 // step size of the array size
#define _TESTSUM_SCALE_NUM		5 // number of scales by factors of _TESTSUM_SCALE_STEP
#define _TESTSUM_REPEAT_NUM		50 // number of test repeats for statistics
#endif

using namespace std;

int main()
{
	int nerr = 0;
    cout << "\n";
	cout << "Program: TestSum\n";
	cout << "  by J. Barthel, FZ-Juelich, 2018\n";
	cout << "\n";
	cout << "- initializations \n";
	size_t nitems = (size_t)_TESTSUM_SCALE_MIN;
	size_t idx = 0;
	unsigned int urnd = 0;
	vector<string> stestname; // list of test names
	stestname.push_back("float straight");
	stestname.push_back("double straight");
	stestname.push_back("float kahan");
	stestname.push_back("float dncs2");
	float* fsum = NULL;
	float* ftim = NULL;
	float* farr = NULL;
	float* fsumcur = NULL;
	float* ftimcur = NULL;
	chrono::high_resolution_clock _clock; // clock with ns step
	auto cl_start = _clock.now(); // clock status
	auto cl_stop = _clock.now();
	auto cl_dif = cl_stop - cl_start;
	fsum = (float*)calloc(_TESTSUM_SCALE_NUM*_TESTSUM_NUM, sizeof(float)); // list of summation results
	ftim = (float*)calloc(_TESTSUM_SCALE_NUM*_TESTSUM_NUM, sizeof(float)); // list of timing results
	if (NULL == fsum || NULL == ftim) {
		cerr << "Error: failed to allocate memory for results.\n";
		nerr = 1; goto _Exit;
	}
	//
	// loop over number of scales
	cout << "- starting tests ... \n";
	for (int isca = 0; isca < _TESTSUM_SCALE_NUM; isca++) {
		cout << "  - current scale (" << isca+1 << "/" << _TESTSUM_SCALE_NUM << "): " << nitems << "\n";
		if (NULL != farr) { free(farr); } // clear data array
		farr = (float*)calloc(nitems, sizeof(float)); // allocate new data array and preset with zeroes
		if (NULL == farr) {
			cerr << "Error: failed to allocate memory for current scale (" << nitems << ").\n";
			nerr = 2; goto _Exit;
		}
		// prepare data with random numbers between 0 and 1
		cout << "  - setting array data: random numbers between 0 and 1.\n";
		for (idx = 0; idx < nitems; idx++) {
			rand_s(&urnd);
			farr[idx] = (float)( (double)urnd / (double)(UINT_MAX - 1) );
		}
		// loop over tests
		for (int itest = 0; itest < _TESTSUM_NUM; itest++) {
			fsumcur = &fsum[isca + itest * _TESTSUM_SCALE_NUM]; // slightly strange array setup, but I want lists depending on array size.
			ftimcur = &ftim[isca + itest * _TESTSUM_SCALE_NUM];
			cout << "  - running " << stestname[itest].c_str() << " ... \n";
			cl_start = _clock.now();
			for (int irep = 0; irep < _TESTSUM_REPEAT_NUM; irep++) {
				switch (itest) {
				case 0:
					fstrsum(farr, nitems, fsumcur);
					break;
				case 1:
					fdstrsum(farr, nitems, fsumcur);
					break;
				case 2:
					fkahan(farr, nitems, fsumcur);
					break;
				case 3:
					fdncs2(farr, nitems, fsumcur);
					break;
				default:
					cerr << "Error: invalid test ID.\n";
					nerr = 3; goto _Exit;
					break;
				}
			} // repeats
			cl_stop = _clock.now();
			cl_dif = cl_stop - cl_start;
			*ftimcur = (float)(chrono::duration <double, nano>(cl_dif).count()/((double)_TESTSUM_REPEAT_NUM*nitems));
			cout << "      result: " << *fsumcur << "  (" << *ftimcur << " ns/item)\n";
		} // tests
		// increase the scale
		nitems *= (size_t)_TESTSUM_SCALE_STEP;
	} // scale


_Exit:
	if (NULL != fsum) { free(fsum); }
	if (NULL != ftim) { free(ftim); }
	if (NULL != farr) { free(farr); }
	cout << "Done.\n";
	cout << "\n";
	exit(nerr);
}

