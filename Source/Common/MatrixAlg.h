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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "Common/Quaternion.h"
#include "Common/Vector.h"

namespace HyVoxel
{
	namespace Algebra
	{
		//typedef struct Matrix Matrix;

		/// Loads identity matrix
		void Matrix_loadIdentity(Matrix * matrix);
		/// Returns identity matrix
		Matrix Matrix_identity();

		/// Initializes matrix with values
		Matrix Matrix_withValues(float m0,  float m4,  float m8,  float m12,
								 float m1,  float m5,  float m9,  float m13,
								 float m2,  float m6,  float m10, float m14,
								 float m3,  float m7,  float m11, float m15);
		/// Creates matrix from direction vectors
		Matrix Matrix_fromDirectionVectors(struct Vector right, struct Vector up, struct Vector front);

		/// Mutiplies matrices
		void Matrix_multiply(Matrix * matrix1, Matrix matrix2);
		/// Returns result of matrix multiplication
		Matrix Matrix_multiplied(Matrix matrix1, Matrix matrix2);

		/// Translates matrix
		void Matrix_translate(Matrix * matrix1, float x, float y, float z);
		/// Returns translated matrix
		Matrix Matrix_translated(Matrix matrix1, float x, float y, float z);

		/// Scales matrix
		void Matrix_scale(Matrix * matrix, float x, float y, float z);
		/// Returns scales matrix
		Matrix Matrix_scaled(Matrix matrix, float x, float y, float z);

		/// Rotates matrix
		void Matrix_rotate(Matrix * matrix, struct Vector axis, float angle);
		/// Returns rotated matrix
		Matrix Matrix_rotated(Matrix matrix, struct Vector axis, float angle);

		/// Applies a shear transofmration to matrix on X
		void Matrix_shearX(Matrix * matrix, float y, float z);
		/// Returns a shear transofmration to matrix on X
		Matrix Matrix_shearedX(Matrix matrix, float y, float z);

		/// Applies a shear transofmration to matrix on Y
		void Matrix_shearY(Matrix * matrix, float x, float z);
		/// Returns a shear transofmration to matrix on Y
		Matrix Matrix_shearedY(Matrix matrix, float x, float z);

		/// Applies a shear transofmration to matrix on Z
		void Matrix_shearZ(Matrix * matrix, float x, float y);
		/// Returns a shear transofmration to matrix on Z
		Matrix Matrix_shearedZ(Matrix matrix, float x, float y);

		/// Applies perspective transformation to matrix
		void Matrix_applyPerspective(Matrix * matrix, float fovY, float aspect, float zNear, float zFar);
		/// Returns the result of applying perspective transformation to a matrix
		Matrix Matrix_perspective(Matrix matrix, float fovY, float aspect, float zNear, float zFar);

		/// Transposes a matrix
		void Matrix_transpose(Matrix * matrix);
		/// Returns matrix transpose
		Matrix Matrix_transposed(Matrix matrix);

		/// Returns matrix determinant
		float Matrix_determinant(Matrix matrix);

		/// Inverts a matrix
		void Matrix_invert(Matrix * matrix);
		/// Reurns inverted matrix
		Matrix Matrix_inverted(Matrix matrix);

		/// Returns vector from multiplication of matrix to another vector
		inline struct Vector Matrix_multiplyVector(Matrix matrix, Vector vector)
		{
			Vector result;
			result.x = ((matrix.m[0] * vector.x) + (matrix.m[4] * vector.y) + (matrix.m[8]  * vector.z) + matrix.m[12]);
			result.y = ((matrix.m[1] * vector.x) + (matrix.m[5] * vector.y) + (matrix.m[9]  * vector.z) + matrix.m[13]);
			result.z = ((matrix.m[2] * vector.x) + (matrix.m[6] * vector.y) + (matrix.m[10] * vector.z) + matrix.m[14]);
			return result;
		}
	}
}

#endif
