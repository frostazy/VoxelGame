/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#include "HyVoxelPrivatePCH.h"

#include "Common/WhiteNoise.h"
#include "Math.h"
#include <stdlib.h>

using namespace HyVoxel;

int CWhiteNoise::samples[WHITE_NOISE_DIM*WHITE_NOISE_DIM*WHITE_NOISE_DIM];

void CWhiteNoise::initialize()
{
	for (int i = 0; i < WHITE_NOISE_DIM*WHITE_NOISE_DIM*WHITE_NOISE_DIM; i++)
	{
		samples[i] = rand();
	}
}

