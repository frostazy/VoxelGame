#pragma once

// Please don't include this file directly.
// Please include "HyVoxelLib.h"

#include <cmath>
#include <algorithm>

namespace HyVoxel {

	#pragma pack(push, 8)

	// Magic numbers for numerical precision.
	#define DELTA			(0.00001f)

	#define HY_PI       3.14159265358979323846   // pi


	// TODO: liang
	// rename to "Vector3f"
	struct Vector
	{
		float x;
		float y;
		float z;

		Vector() {}
		Vector(float inX, float inY, float inZ) : x(inX), y(inY), z(inZ) {}
		Vector(float v) : x(v), y(v), z(v) {}

		inline bool operator==(const Vector& v) const { return x==v.x && y==v.y && z==v.z; }

		inline Vector operator+(const Vector& v) const { return Vector(x + v.x, y + v.y, z + v.z); }
		inline Vector operator+(float bias) const { return Vector(x + bias, y + bias, z + bias); }
		inline Vector operator+=(const Vector& v) {
			x += v.x; y += v.y; z += v.z;
			return *this;
		}

		inline Vector operator-(const Vector& v) const { return Vector(x - v.x, y - v.y, z - v.z); }
		inline Vector operator-(float bias) const { return Vector(x - bias, y - bias, z - bias); }
		inline Vector operator-=(const Vector& v) {
			x -= v.x; y -= v.y; z -= v.z;
			return *this;
		}

		inline Vector operator-() const { return Vector(-x, -y, -z); }

		inline Vector operator*(const Vector& v) const { return Vector(x * v.x, y * v.y, z * v.z); }
		inline Vector operator*(float scale) const { return Vector(x * scale, y * scale, z * scale); }

		inline Vector operator/(const Vector& v) const { return Vector(x / v.x, y / v.y, z / v.z); }
		inline Vector operator/(float scale) const {
			const float rscale = 1.f / scale;
			return Vector(x * rscale, y * rscale, z * rscale);
		}
		inline Vector operator/=(float v) {
			const float rv = 1.0f / v;
			x *= rv; y *= rv; z *= rv;
			return *this;
		}

		inline bool Normalize(float tolerance = 1.e-8f) {
			const float squareSum = x*x + y*y + z*z;
			if (squareSum > tolerance)
			{
				const float scale = 1.0f / std::sqrt(squareSum);
				x *= scale; y *= scale; z *= scale;
				return true;
			}
			return false;
		}

		inline Vector GetSafeNormal(float tolerance = 1.e-8f) const {
			const float squareSum = x*x + y*y + z*z;

			// Not sure if it's safe to add tolerance in there. Might introduce too many errors
			if (squareSum == 1.f)
			{
				return *this;
			}
			else if (squareSum < tolerance)
			{
				return Vector(0.0f);
			}
			const float scale = 1.0f / std::sqrt(squareSum);
			return Vector(x*scale, y*scale, z*scale);
		}

		inline float GetAbsMax() { return std::max(std::max(std::abs(x), std::abs(y)), std::abs(z)); }

		Vector RotateAngleAxis(const float angleDeg, const Vector& axis) const;

		inline static float DotProduct(const Vector& vector1, const Vector& vector2) {
			return ((vector1.x * vector2.x) + (vector1.y * vector2.y) + (vector1.z * vector2.z));
		}

		inline static Vector CrossProduct(const Vector& vector1, const Vector& vector2) {
			return Vector (
				(vector1.y * vector2.z) - (vector1.z * vector2.y),
				(vector1.z * vector2.x) - (vector1.x * vector2.z),
				(vector1.x * vector2.y) - (vector1.y * vector2.x)
			);
		}

		inline static float Dist(const Vector& v1, const Vector& v2) {
			return std::sqrt(std::pow(v2.x - v1.x, 2) + std::pow(v2.y - v1.y, 2) + std::pow(v2.z - v1.z, 2));
		}

		inline static void CreateOrthonormalBasis(Vector& xAxis, Vector& yAxis, Vector& zAxis)
		{
			xAxis -= zAxis * DotProduct(xAxis, zAxis) / DotProduct(zAxis, zAxis);
			yAxis -= zAxis * DotProduct(yAxis, zAxis) / DotProduct(zAxis, zAxis);

			if (DotProduct(xAxis, xAxis) < DELTA*DELTA)
			{
				xAxis = CrossProduct(yAxis, zAxis);
			}

			if (DotProduct(yAxis, yAxis) < DELTA*DELTA)
			{
				yAxis = CrossProduct(xAxis, zAxis);
			}

			xAxis.Normalize(); yAxis.Normalize(); zAxis.Normalize();
		}



		#if 0
		Vector& operator=(const HyVoxel::Vector& v)
		{
			x = v.X;
			y = v.Y;
			z = v.Z;
			return *this;
		}
		#endif
	};

	inline Vector operator*(float scale, const Vector& v) { return v.operator*(scale); }

	inline Vector Vector::RotateAngleAxis(const float angleDeg, const Vector& axis) const
	{
		float rad = angleDeg / 180.0f * (float)HY_PI;
		float s = std::sinf(rad);
		float c = std::cosf(rad);

		const float xx = axis.x * axis.x;
		const float yy = axis.y * axis.y;
		const float zz = axis.z * axis.z;

		const float xy = axis.x * axis.y;
		const float yz = axis.y * axis.z;
		const float zx = axis.z * axis.x;

		const float xs = axis.x * s;
		const float ys = axis.y * s;
		const float zs = axis.z * s;

		const float OMC = 1.f - c;

		return Vector(
			(OMC * xx + c) * x + (OMC * xy - zs) * y + (OMC * zx + ys) * z,
			(OMC * xy + zs) * x + (OMC * yy + c) * y + (OMC * yz - xs) * z,
			(OMC * zx - ys) * x + (OMC * yz + xs) * y + (OMC * zz + c) * z
		);
	}


	struct Vector2f
	{
		float x;
		float y;

		Vector2f(float inX, float inY) : x(inX), y(inY) { }

		inline Vector2f operator+(const Vector2f& v) const { return Vector2f(x + v.x, y + v.y); }
		inline Vector2f operator+(float bias) const { return Vector2f(x + bias, y + bias); }

		inline Vector2f operator-(const Vector2f& v) const { return Vector2f(x - v.x, y - v.y); }
		inline Vector2f operator-(float bias) const { return Vector2f(x - bias, y - bias); }

		inline Vector2f operator*(const Vector2f& v) const { return Vector2f(x * v.x, y * v.y); }
		inline Vector2f operator*(float scale) const { return Vector2f(x * scale, y * scale); }

		inline Vector2f operator/(const Vector2f& v) const { return Vector2f(x / v.x, y / v.y); }
		inline Vector2f operator/(float scale) const {
			const float rscale = 1.f / scale;
			return Vector2f(x * rscale, y * rscale);
		}

	};

	// TODO: liang
	// rename to "Matrix4x4f"
	struct Matrix
	{
		float m[16];

		Matrix() {}
		#if 0
		Matrix(const FMatrix& mtx)
		{
			typedef float Float4x4[4][4];
			const Float4x4& M = *((const Float4x4*)&mtx);
			m[0] = M[0][0]; m[1] = M[0][1]; m[2] = M[0][2]; m[3] = M[0][3];
			m[4] = M[1][0]; m[5] = M[1][1]; m[6] = M[1][2]; m[7] = M[1][3];
			m[8] = M[2][0]; m[9] = M[2][1]; m[10] = M[2][2]; m[11] = M[2][3];
			m[12] = M[3][0]; m[13] = M[3][1]; m[14] = M[3][2]; m[15] = M[3][3];
		}
		#endif
	};

	#pragma pack(pop)

}