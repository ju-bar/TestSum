// "summation.h"
//
// Declaration of summation routines over a memory block
//
// J. Barthel, Forschungszentrum Jülich GmbH, 52425 Jülich, Germany
// Copyright (c) 2018
#pragma once

#define _SUMMATION_BUFFER		0x1000	// number of temporary buffer items for summation


// straight float summation without attempting error correction
void fstrsum(float* a, size_t n, float *s);

// straight float summation with double accumulator
void fdstrsum(float* a, size_t n, float *s);

// kahan summation with float error correction
// see <https://en.wikipedia.org/wiki/Kahan_summation_algorithm>
void fkahan(float* a, size_t n, float *s);

// strided 2-fold butterfly summation
// ! Warning: This routine may modify the input a.
//   Make a backup of the data or provide a copy if you still need it !
void fdncs2m(float* a, size_t n, float *s);

// strided 2-fold butterfly summation without input modification
void fdncs2(float* a, size_t n, float *s);
