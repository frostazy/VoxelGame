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

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "HyMath.h"

namespace HyVoxel
{
	/// 3D algebra functions and data types
	namespace Algebra
	{
		/// Initializes a vector
		inline Vector Vector_withValues(float x, float y, float z)
		{
			Vector vector;
			vector.x = x;
			vector.y = y;
			vector.z = z;
			return vector;
		}

		/// Normalizes vector
		bool Vector_normalize(Vector * vector);
		/// Returns normalized vector
		Vector Vector_normalized(Vector vector);

		/// Returns the magnitude of the vector
		float Vector_magnitude(Vector vector);
		/// Returns the squared magnitude
		float Vector_magnitudeSquared(Vector vector);
		/// Returns sum of two vectors
		Vector Vector_add(Vector vector1, Vector vector2);
		/// Returns substraction of two vectors
		Vector Vector_subtract(Vector vector1, Vector vector2);
		/// Returns dot product between two vectors
		float Vector_dot(Vector vector1, Vector vector2);
		/// Returns cross product between two vectors
		Vector Vector_cross(Vector vector1, Vector vector2);
		Vector Vector_mutiply(float value, Vector vector);

		const double piover180 = 0.0174532925;
		const float piover180f = 0.0174532925f;

		Vector closestPtPointTriangle(Vector p, Vector a, Vector b, Vector c);
		bool testSphereTriangle(Vector center, float r, Vector a, Vector b, Vector c, Vector &p);
		double sign2D(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y);

		bool pointInsideTriangle3D(Vector p, Vector v0, Vector v1, Vector v2);
		int segmentIntersectTriangle3D(Vector p0, Vector p1, Vector v0, Vector v1, Vector v2, Vector* p);
	}
}
#endif
