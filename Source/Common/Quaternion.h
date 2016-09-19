/*
  Copyright (C) 2008 Alex Diener

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Alex Diener adiener@sacredsoftware.net
*/

#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include "Common/MatrixAlg.h"
#include "Common/Vector.h"

namespace HyVoxel
{
	namespace Algebra
	{
		typedef struct Quaternion Quaternion;

		/// A Quaternion object
		struct Quaternion
		{
			float x;
			float y;
			float z;
			float w;
		};

		/// Load identity quaternion
		void Quaternion_loadIdentity(Quaternion * quaternion);
		/// Returns identity quaternion
		Quaternion Quaternion_identity();
		/// Initializes quaternion
		Quaternion Quaternion_withValues(float x, float y, float z, float w);

		/// Initializes quaternion from a vector
		Quaternion Quaternion_fromVector(struct Vector vector);
		/// Returns vector from quaternion
		struct Vector Quaternion_toVector(Quaternion quaternion);
		/// Creates quaternion from axis angle
		Quaternion Quaternion_fromAxisAngle(struct Vector axis, float angle);
		/// Converts quaternion to axis angle
		void Quaternion_toAxisAngle(Quaternion quaternion, struct Vector * axis, float * angle);
		/// Converts quaternion to Matrix
		struct Matrix Quaternion_toMatrix(Quaternion quaternion);

		/// Normalizes quaternion
		void Quaternion_normalize(Quaternion * quaternion);
		/// Returns normalized quaternion
		Quaternion Quaternion_normalized(Quaternion quaternion);

		/// Multiplies quaternion
		void Quaternion_multiply(Quaternion * quaternion1, Quaternion quaternion2);
		/// Returns mutiplication of quaternion
		Quaternion Quaternion_multiplied(Quaternion quaternion1, Quaternion quaternion2);
		/// Interpolates between two quaternions
		Quaternion Quaternion_slerp(Quaternion start, Quaternion end, float alpha);

		/// Rotates quaternion
		void Quaternion_rotate(Quaternion * quaternion, struct Vector axis, float angle);
		/// Returns rotated quaternion
		Quaternion Quaternion_rotated(Quaternion quaternion, struct Vector axis, float angle);

		/// Inverts quaternion
		void Quaternion_invert(Quaternion * quaternion);
		/// Returns inverted quaternion
		Quaternion Quaternion_inverted(Quaternion quaternion);

		/// Returns quaternion to vector multiplication
		struct Vector Quaternion_multiplyVector(Quaternion quaternion, struct Vector vector);
	}
}

#endif
