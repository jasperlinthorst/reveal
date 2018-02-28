/////////////////////////////////////////////////////////////////
// Defaults.h
//
// Default constants for use in PROBCONS.  The emission
// probabilities were computed using the program used to build
// the BLOSUM62 matrix from the BLOCKS 5.0 dataset.  Transition
// parameters were obtained via unsupervised EM training on the
// BALIBASE 2.0 benchmark alignment database.
/////////////////////////////////////////////////////////////////

#ifndef DEFAULTS_H
#define DEFAULTS_H

#include <string>

using namespace std;

float initDistrib1Default[] = { 0.33333333f, 0.33333333f, 0.33333333f };

float gapOpen1Default[] = { 0.013652682f, 0.013652682f };
float gapExtend1Default[] = { 0.9744453f, 0.9744453f };

float initDistrib2Default[] = { 0.2f, 0.2f, 0.2f, 0.2f, 0.2f };

float gapOpen2Default[] = { 0.012986835f, 0.012986835f, 0.0018214881f, 0.0018214881f};
float gapExtend2Default[] = { 0.7126062401851738f, 0.7126062401851738f, 0.99656342579062f, 0.99656342579062f};

string alphabetDefault = "ACGTN";

float emitSingleDefault[5] = {
	0.2f, 0.2f, 0.2f, 0.2f, 0.2f 
};

float emitPairsDefault[5][5] = {
	{0.12064298095701059f, 0.0f, 0.0f, 0.0f, 0.0f}, 
	{0.010367271172731285f*2, 0.12064298095701059f, 0.0f, 0.0f, 0.0f},
	{0.01862247669752685f*2, 0.010367271172731285f*2, 0.12064298095701059f, 0.0f, 0.0f},
	{0.010367271172731285f*2, 0.01862247669752685f*2, 0.010367271172731285f*2, 0.12064298095701059f, 0.0f},
	{0.04f*2, 0.04f*2, 0.04f*2, 0.04f*2, 0.04f}
};

#endif