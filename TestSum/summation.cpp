// summation.cpp : This is where you find summation routines
// J. Barthel, Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
// Copyright (c) 2018

#include "pch.h"
#include "summation.h"
#include <math.h>
#include <malloc.h>
#include <memory.h>


using namespace std;


void fstrsum(float *a, size_t n, float *s)
{
	*s = 0.f;
	if (n > 0) {
		for (size_t i = 0; i < n; i++) {
			*s += a[i];
		}
	}
	return;
}


void fdstrsum(float *a, size_t n, float *s)
{
	*s = 0.f;
	double dtmp = 0.0;
	if (n > 0) {
		for (size_t i = 0; i < n; i++) {
			dtmp += (double)a[i];
		}
		*s = (float)dtmp;
	}
	return;
}


void fkahan(float *a, size_t n, float *s)
{
	*s = 0.f;
	if (n > 0) {
		float c = 0.f, t = 0.f, y = 0.f;
		for (size_t i = 0; i < n; i++) {
			y = a[i] - c;
			t = *s + y; // temp result
			c = (t - *s) - y; // new error
			*s = t; // updated sum
		}
	}
}


void fdncs2m(float *a, size_t n, float *s)
{
	// roll out to 5
	if (n <= 0) { *s = 0.; return; }
	if (n == 1) { *s = a[0]; return; }
	if (n == 2) { *s = a[0] + a[1]; return; }
	if (n == 3) { *s = a[0] + a[1] + a[2]; return; }
	if (n == 4) { *s = a[0] + a[1] + a[2] + a[3]; return; }
	if (n == 5) { *s = a[0] + a[1] + a[2] + a[3] + a[4]; return; }
	// recurse shift 2
	*s = 0.f;
	size_t n0 = 0, n1 = 1, n2 = 2, nc = 0, idx = 0, itmp = 0;
	// calculate number of strides through the buffer
	size_t ntmp = (size_t)ceil((double)n / (double)_SUMMATION_BUFFER);
	float r = 0.0f, t = 0.0f;
	float* dst = s; // preset destination with single result number
	if (ntmp > 1) { // there will be more than one stride -> more than one number
		dst = (float*)calloc(ntmp, sizeof(float)); // allocate new destination buffer
	}
	else { // butterfly on "a" directly, using s as output target
		/*for (idx = 0; idx < n; idx++) {
			*s += a[idx];
		}
		t = *s;
		*s = 0.0f;*/
		r = 0.0f;
		n2 = 2; n1 = 1;
		while (n2 <= n) {
			for (idx = n2 - 1; idx < n; idx += n2) {
				a[idx] += a[idx - n1];
			}
			if (n1 <= n % n2) { // handle left-over component
				r += a[idx - n1];
			}
			n1 = n2;
			n2 = n2 << 1;
		}
		*s += (a[n1 - 1] + r);
		return;
	}
	float* tmp = a; // pre-link to beginning of input temp buffer
	while (n0 < n) { // loop over strides, using dst as output target
		nc = __min(_SUMMATION_BUFFER, n - n0); // number of copied items
		if (nc == _SUMMATION_BUFFER) { // butterfly on full 2^M buffer length. This is repeated ntmp-1 times.
			tmp = &a[n0]; // link to offset in a
			n2 = 2; n1 = 1;
			while (n2 < _SUMMATION_BUFFER) {
				for (idx = n2 - 1; idx < _SUMMATION_BUFFER; idx += n2) {
					tmp[idx] += tmp[idx - n1];
				}
				n1 = n2;
				n2 = n2 << 1;
			}
			dst[itmp] += (tmp[n1 - 1] + tmp[n2 - 1]); // store intermediate result in stride slot of destination
		}
		else { // butterfly on remaining buffer length (not power of two), this happens only once!
			// roll-out to 5
			if (1 == nc) { dst[itmp] = a[n0]; }
			else if (2 == nc) { dst[itmp] = a[n0] + a[1 + n0]; }
			else if (3 == nc) { dst[itmp] = a[n0] + a[1 + n0] + a[2 + n0]; }
			else if (4 == nc) { dst[itmp] = a[n0] + a[1 + n0] + a[2 + n0] + a[3 + n0]; }
			else if (5 == nc) { dst[itmp] = a[n0] + a[1 + n0] + a[2 + n0] + a[3 + n0] + a[4 + n0]; }
			else if (5 < nc) {
				tmp = &a[n0]; // link to offset in a
				/*for (idx = 0; idx < nc; idx++) {
					dst[itmp] += tmp[idx];
				}
				t = dst[itmp];
				dst[itmp] = 0.0f;*/
				n2 = 2; n1 = 1;
				while (n2 <= nc) {
					for (idx = n2 - 1; idx < nc; idx += n2) {
						tmp[idx] += tmp[idx - n1];
					}
					if (n1 <= nc % n2) { // handle left-over component
						r += tmp[idx - n1];
					}
					n1 = n2;
					n2 = n2 << 1;
				}
				dst[itmp] += (tmp[n1 - 1] + r);
			}
		}
		n0 += _SUMMATION_BUFFER;
		itmp++;
	}
	if (ntmp > 1) { // recurse dst if more than one buffer stride happened
		fdncs2m(dst, ntmp, s); // sum on dst buffer to s
		free(dst); // release dst buffer memory
	}
	return; // s received the result
}


void fdncs2(float *a, size_t n, float *s)
{
	// roll out to 5
	if (n <= 0) { *s = 0.; return; }
	if (n == 1) { *s = a[0]; return; }
	if (n == 2) { *s = a[0] + a[1]; return; }
	if (n == 3) { *s = a[0] + a[1] + a[2]; return; }
	if (n == 4) { *s = a[0] + a[1] + a[2] + a[3]; return; }
	if (n == 5) { *s = a[0] + a[1] + a[2] + a[3] + a[4]; return; }
	// recurse shift 2
	*s = 0.f;
	size_t n0 = 0, n1 = 1, n2 = 2, nc = 0, idx = 0, itmp = 0;
	// calculate number of strides through the buffer
	size_t ntmp = (size_t)ceil((double)n / (double)_SUMMATION_BUFFER);
	size_t nbb = sizeof(float)*_SUMMATION_BUFFER;
	size_t nbc = 0;
	float r = 0.0f, t = 0.0f;
	float* tmp = (float*)malloc(nbb); // alloc working buffer
	float* dst = s; // preset destination with single result number
	if (ntmp > 1) { // there will be more than one stride -> more than one number
		dst = (float*)calloc(ntmp, sizeof(float)); // allocate new destination buffer
	}
	else { // small n -> butterfly on tmp (copy of a), using s as output target
		/*for (idx = 0; idx < n; idx++) {
			*s += a[idx];
		}
		t = *s;
		*s = 0.0f;*/
		memcpy(tmp, a, sizeof(float)*n); // prepare summation buffer
		r = 0.0f;
		n2 = 2; n1 = 1;
		while (n2 <= n) {
			for (idx = n2 - 1; idx < n; idx += n2) {
				tmp[idx] += tmp[idx - n1];
			}
			if (n1 <= n % n2) { // handle left-over component
				r += tmp[idx - n1];
			}
			n1 = n2;
			n2 = n2 << 1;
		}
		*s += (tmp[n1 - 1] + r);
		free(tmp);
		return;
	}
	while (n0 < n) { // loop over strides, using dst as output target
		nc = __min(_SUMMATION_BUFFER, n - n0); // number of copied items
		if (nc == _SUMMATION_BUFFER) { // butterfly on full 2^M buffer length. This is repeated ntmp-1 times.
			memcpy(tmp, &a[n0], nbb); // prepare summation buffer
			n2 = 2; n1 = 1;
			while (n2 < _SUMMATION_BUFFER) {
				for (idx = n2 - 1; idx < _SUMMATION_BUFFER; idx += n2) {
					tmp[idx] += tmp[idx - n1];
				}
				n1 = n2;
				n2 = n2 << 1;
			}
			dst[itmp] += (tmp[n1 - 1] + tmp[n2 - 1]); // store intermediate result in stride slot of destination
		}
		else { // butterfly on remaining buffer length (not power of two), this happens only once!
			// roll-out to 5
			if (1 == nc) { dst[itmp] = a[n0]; }
			else if (2 == nc) { dst[itmp] = a[n0] + a[1 + n0]; }
			else if (3 == nc) { dst[itmp] = a[n0] + a[1 + n0] + a[2 + n0]; }
			else if (4 == nc) { dst[itmp] = a[n0] + a[1 + n0] + a[2 + n0] + a[3 + n0]; }
			else if (5 == nc) { dst[itmp] = a[n0] + a[1 + n0] + a[2 + n0] + a[3 + n0] + a[4 + n0]; }
			else if (5 < nc) {
				memcpy(tmp, &a[n0], sizeof(float)*nc); // prepare summation buffer
				/*for (idx = 0; idx < nc; idx++) {
					dst[itmp] += tmp[idx];
				}
				t = dst[itmp];
				dst[itmp] = 0.0f;*/
				n2 = 2; n1 = 1;
				while (n2 <= nc) {
					for (idx = n2 - 1; idx < nc; idx += n2) {
						tmp[idx] += tmp[idx - n1];
					}
					if (n1 <= nc % n2) { // handle left-over component
						r += tmp[idx - n1];
					}
					n1 = n2;
					n2 = n2 << 1;
				}
				dst[itmp] += (tmp[n1 - 1] + r);
			}
		}
		n0 += _SUMMATION_BUFFER;
		itmp++;
	}
	free(tmp); // release tmp buffer memory
	if (ntmp > 1) { // recurse dst if more than one buffer stride happened
		fdncs2m(dst, ntmp, s); // sum on dst buffer to s
		free(dst); // release dst buffer memory
	}
	return; // s received the result
}
