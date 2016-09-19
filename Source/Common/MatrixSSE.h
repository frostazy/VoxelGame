/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#if (defined(_WINDOWS) || defined(__SSE2__))

#pragma once

#include "HyVoxelConfig.h"

namespace HyVoxel
{
	namespace Algebra
	{
		/// Simplified 10-Element QEF Matrix for Simplification Error Calculations (in SSE)
		class QEFMatrixSSE
		{
		public:
			QEFMatrixSSE();
			QEFMatrixSSE(const double c);
			QEFMatrixSSE(const QEFMatrixSSE& q);
			QEFMatrixSSE(const double plane[]);

			const double* getRawPtr() const;
			void reset();

			double operator[](const int c) const;
			const QEFMatrixSSE operator+(const QEFMatrixSSE& q) const;
			QEFMatrixSSE& operator+=(const QEFMatrixSSE& q);

		private:
			__declspec(align(16)) double m[10];
		};

		/// Positional Error for simplified QEF Matrix
		double vertexErrorSSE(const QEFMatrixSSE& q, const double x, const double y, const double z);

		/// SSE Vector Operations
		class VectorSSE
		{
		public:
			static float dot(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1);
			static float dot3(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2);
			static void cross(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout);
			static void cross3(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2, float& xout, float& yout, float& zout);
			static float crossDot(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1);
			static float crossDot3(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2);
			static float length(const float x, const float y, const float z);
			static float lengthSqr(const float x, const float y, const float z);
			static void normalize(const float x, const float y, const float z, float& xout, float& yout, float& zout);
			static void normalize2(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout);
			static void add(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout);
			static void subtract(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout);
			static void multiply(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout);
			static void divide(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout);

			static double dot(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1);
			static double dot3(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);
			static void cross(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, double& xout, double& yout, double& zout);
			static void cross3(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, double& xout, double& yout, double& zout);
			static double crossDot(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1);
			static double crossDot3(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);
			static double length(const double x, const double y, const double z);
			static double lengthSqr(const double x, const double y, const double z);
			static void normalize(const double x, const double y, const double z, double& xout, double& yout, double& zout);
		};

		/// SSE Inverse-Square-Root
		float invSQR(const float x);
		double invSQR(const double x);

		/// SSE Triangle Centroid
		void triCentroid(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2, float& xout, float& yout, float& zout);

		/// SSE (Actual) Triangle Area
		double triArea(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);
	};
};

#endif