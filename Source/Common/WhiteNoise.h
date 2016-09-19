/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

namespace HyVoxel
{
	/// Implements deterministic 3D white noise for ints
	class CWhiteNoise
	{
	private:
		static const int WHITE_NOISE_DIM = 128;
	public:
		/// Initialize the Noise sample buffer
		static void initialize();
		/// Returns an integer noise value for the supplied 3D point. The value is in the range from 0 to RAND_MAX - 1.
		static int getValue(unsigned int x, unsigned int y, unsigned int z)
		{
			x %= WHITE_NOISE_DIM;
			y %= WHITE_NOISE_DIM;
			z %= WHITE_NOISE_DIM;
			return samples[y*WHITE_NOISE_DIM*WHITE_NOISE_DIM + z*WHITE_NOISE_DIM + x];
		}
	private:
		static int samples[WHITE_NOISE_DIM*WHITE_NOISE_DIM*WHITE_NOISE_DIM];
	};
}