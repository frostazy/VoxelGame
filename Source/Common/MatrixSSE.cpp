#if (defined(_WINDOWS) || defined(__SSE2__))

#include "HyVoxelPrivatePCH.h"

#include "MatrixSSE.h"
#include <emmintrin.h>
#include <xmmintrin.h>

HyVoxel::Algebra::QEFMatrixSSE::QEFMatrixSSE()
{
}

HyVoxel::Algebra::QEFMatrixSSE::QEFMatrixSSE(const double c)
{
	if (c == 0.0)
	{
		memset(this->m, 0, 10 * sizeof(double));
	}
	else
	{
		for (int i = 0; i < 10; i++)
		{
			this->m[i] = c;
		}
	}
}

HyVoxel::Algebra::QEFMatrixSSE::QEFMatrixSSE(const QEFMatrixSSE& q)
{
	memcpy(this->m, q.m, 10 * sizeof(double));
}

HyVoxel::Algebra::QEFMatrixSSE::QEFMatrixSSE(const double plane[])
{
	__m128d* p = (__m128d*)plane;
	__m128d* s = (__m128d*)this->m;

	// aa, ab
	__m128d aa = _mm_unpacklo_pd(*p, *p);									// ab, ab -> aa
	*s++ = _mm_mul_pd(aa, *p);

	// ac, ad
	*s++ = _mm_mul_pd(aa, *(p + 1));

	// bb, bc
	__m128d bb = _mm_unpackhi_pd(*p, *p);									// ab, ab -> bb
	__m128d bc = _mm_shuffle_pd(*p, *(p + 1), _MM_SHUFFLE2(0, 1));			// ab, cd -> bc
	*s++ = _mm_mul_pd(bb, bc);

	// bd, cc
	__m128d dc = _mm_shuffle_pd(*(p + 1), *(p + 1), 1);						// cd, cd -> dc
	*s++ = _mm_mul_pd(bc, dc);

	// cd, dd
	__m128d dd = _mm_unpackhi_pd(*(p + 1), *(p + 1));						// cd, cd -> dd
	*s = _mm_mul_pd(dd, *(p + 1));
}

const double* HyVoxel::Algebra::QEFMatrixSSE::getRawPtr() const
{
	return this->m;
}

void HyVoxel::Algebra::QEFMatrixSSE::reset()
{
	memset(this->m, 0, 10 * sizeof(double));
}

double HyVoxel::Algebra::QEFMatrixSSE::operator[](const int c) const
{
	return this->m[c];
}

const HyVoxel::Algebra::QEFMatrixSSE HyVoxel::Algebra::QEFMatrixSSE::operator+(const QEFMatrixSSE& q) const
{
	QEFMatrixSSE r;

	__m128d* s = (__m128d*)r.getRawPtr();
	__m128d* p1 = (__m128d*)this->m;
	__m128d* p2 = (__m128d*)q.getRawPtr();

	// manually unrolling the loop increases perf by a factor of ~2
	s[0] = _mm_add_pd(p1[0], p2[0]);
	s[1] = _mm_add_pd(p1[1], p2[1]);
	s[2] = _mm_add_pd(p1[2], p2[2]);
	s[3] = _mm_add_pd(p1[3], p2[3]);
	s[4] = _mm_add_pd(p1[4], p2[4]);

	return r;
}

HyVoxel::Algebra::QEFMatrixSSE& HyVoxel::Algebra::QEFMatrixSSE::operator+=(const QEFMatrixSSE& q)
{
	__m128d* p1 = (__m128d*)this->m;
	__m128d* p2 = (__m128d*)q.getRawPtr();

	p1[0] = _mm_add_pd(p1[0], p2[0]);
	p1[1] = _mm_add_pd(p1[1], p2[1]);
	p1[2] = _mm_add_pd(p1[2], p2[2]);
	p1[3] = _mm_add_pd(p1[3], p2[3]);
	p1[4] = _mm_add_pd(p1[4], p2[4]);

	return *this;
}

double HyVoxel::Algebra::vertexErrorSSE(const QEFMatrixSSE& q, const double x, const double y, const double z)
{
	__m128d* s = (__m128d*)q.getRawPtr();
	__m128d p[5] = 
	{
		_mm_set_pd(2.0*y,	x),
		_mm_set_pd(2.0,		2.0*z),
		_mm_set_pd(2.0*y,	y),
		_mm_set_pd(z,		2.0),
		_mm_set_pd(1.0,		2.0*z)
	};
	__m128d xx = _mm_set_pd(x, x);
	p[0] = _mm_mul_pd(p[0], xx);
	p[1] = _mm_mul_pd(p[1], xx);
	__m128d zy = _mm_set_pd(z, y);
	p[2] = _mm_mul_pd(p[2], zy);
	p[3] = _mm_mul_pd(p[3], zy);

	__m128d et, error = _mm_mul_pd(s[0], p[0]);

	et = _mm_mul_pd(s[1], p[1]);
	error = _mm_add_pd(error, et);

	et = _mm_mul_pd(s[2], p[2]);
	error = _mm_add_pd(error, et);

	et = _mm_mul_pd(s[3], p[3]);
	error = _mm_add_pd(error, et);

	et = _mm_mul_pd(s[4], p[4]);
	error = _mm_add_pd(error, et);
	error = _mm_add_sd(_mm_unpackhi_pd(error, error), error);

	return *reinterpret_cast<double*>(&error);
}

namespace HyVoxel
{
	namespace Algebra
	{
		inline __m128 ftov(const float x, const float y, const float z)
		{
			return _mm_set_ps(0.f, z, y, x);
		}

		inline void vtof(const __m128& v, float& x, float& y, float& z)
		{
			union
			{
				__m128 v;
				float f[4];
			} u;
			u.v = v;
			x = u.f[0];
			y = u.f[1];
			z = u.f[2];
		}

		inline __m128 dotv(const __m128 v0, const __m128 v1)
		{
			__m128 t0, t1;
			t0 = _mm_mul_ps(v0, v1);
			t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 3, 0, 1));
			t0 = _mm_add_ps(t0, t1);
			t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(0, 1, 2, 3));
			return _mm_add_ps(t0, t1);
		}
	};
};

float HyVoxel::Algebra::VectorSSE::dot(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1)
{
	const __m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);
	return _mm_cvtss_f32(dotv(v0, v1));
}

float HyVoxel::Algebra::VectorSSE::dot3(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2)
{
	const __m128 v = ftov(x0, y0, z0);
	__m128 v0 = ftov(x1, y1, z1);
	__m128 v1 = ftov(x2, y2, z2);
	v0 = _mm_sub_ps(v0, v);
	v1 = _mm_sub_ps(v1, v);
	return _mm_cvtss_f32(dotv(v0, v1));
}

void HyVoxel::Algebra::VectorSSE::cross(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout)
{
	__m128 t0, t1;
	const __m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);
	t0 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 1, 0, 2)));
	t1 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1)));
	t0 =_mm_sub_ps(t0, t1);
	vtof(t0, xout, yout, zout);
}

void HyVoxel::Algebra::VectorSSE::cross3(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2, float& xout, float& yout, float& zout)
{
	__m128 t0, t1;
	const __m128 v = ftov(x0, y0, z0);
	__m128 v0 = ftov(x1, y1, z1);
	__m128 v1 = ftov(x2, y2, z2);
	v0 = _mm_sub_ps(v0, v);
	v1 = _mm_sub_ps(v1, v);
	t0 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 1, 0, 2)));
	t1 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1)));
	t0 =_mm_sub_ps(t0, t1);
	vtof(t0, xout, yout, zout);
}

float HyVoxel::Algebra::VectorSSE::crossDot(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1)
{
	__m128 t0, t1;
	const __m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);

	// cross
	t0 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 1, 0, 2)));
	t1 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1)));
	t0 =_mm_sub_ps(t0, t1);

	// dot
	return _mm_cvtss_f32(dotv(t0, t0));
}

float HyVoxel::Algebra::VectorSSE::crossDot3(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2)
{
	__m128 t0, t1;
	const __m128 v = ftov(x0, y0, z0);
	__m128 v0 = ftov(x1, y1, z1);
	__m128 v1 = ftov(x2, y2, z2);
	v0 = _mm_sub_ps(v0, v);
	v1 = _mm_sub_ps(v1, v);

	// cross
	t0 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 1, 0, 2)));
	t1 = _mm_mul_ps(_mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 0, 2, 1)));
	t0 =_mm_sub_ps(t0, t1);

	// dot
	return _mm_cvtss_f32(dotv(t0, t0));
}

float HyVoxel::Algebra::VectorSSE::length(const float x, const float y, const float z)
{
	const __m128 v = ftov(x, y, z);
	return _mm_cvtss_f32(_mm_sqrt_ss(dotv(v, v)));
}

float HyVoxel::Algebra::VectorSSE::lengthSqr(const float x, const float y, const float z)
{
	const __m128 v = ftov(x, y, z);
	return _mm_cvtss_f32(dotv(v, v));
}

void HyVoxel::Algebra::VectorSSE::normalize(const float x, const float y, const float z, float& xout, float& yout, float& zout)
{
	__m128 v = ftov(x, y, z);
	v = _mm_mul_ps(v, _mm_rsqrt_ps(dotv(v, v)));
	vtof(v, xout, yout, zout);
}

void HyVoxel::Algebra::VectorSSE::normalize2(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout)
{
	__m128 v = ftov(x1, y1, z1);
	const __m128 v0 = ftov(x0, y0, z0);
	v = _mm_sub_ps(v, v0);
	v = _mm_mul_ps(v, _mm_rsqrt_ps(dotv(v, v)));
	vtof(v, xout, yout, zout);
}

void HyVoxel::Algebra::VectorSSE::add(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout)
{
	const __m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);
	__m128 t = _mm_add_ps(v0, v1);
	vtof(t, xout, yout, zout);
}

void HyVoxel::Algebra::VectorSSE::subtract(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout)
{
	const __m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);
	__m128 t = _mm_sub_ps(v0, v1);
	vtof(t, xout, yout, zout);
}

void HyVoxel::Algebra::VectorSSE::multiply(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout)
{
	const __m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);
	__m128 t = _mm_mul_ps(v0, v1);
	vtof(t, xout, yout, zout);
}

void HyVoxel::Algebra::VectorSSE::divide(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, float& xout, float& yout, float& zout)
{
	const __m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);
	__m128 t = _mm_div_ps(v0, v1);
	vtof(t, xout, yout, zout);
}

double HyVoxel::Algebra::VectorSSE::dot(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1)
{
	__m128d t0, t1;
	const __m128d v0xy = _mm_set_pd(y0, x0);
	const __m128d v0zw = _mm_set_pd(0.0, z0);
	const __m128d v1xy = _mm_set_pd(y1, x1);
	const __m128d v1zw = _mm_set_pd(0.0, z1);
	t0 = _mm_mul_pd(v0xy, v1xy);
	t1 = _mm_mul_sd(v0zw, v1zw);
	t0 = _mm_add_sd(t0, t1);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	return *reinterpret_cast<double*>(&t1);
}

double HyVoxel::Algebra::VectorSSE::dot3(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	__m128d v0xy = _mm_set_pd(y1, x1);
	__m128d v0zw = _mm_set_pd(0.0, z1);
	__m128d v1xy = _mm_set_pd(y2, x2);
	__m128d v1zw = _mm_set_pd(0.0, z2);
	__m128d t0 = _mm_set_pd(y0, x0);
	__m128d t1 = _mm_set_pd(0.0, z0);
	v0xy = _mm_sub_pd(v0xy, t0);
	v0zw = _mm_sub_sd(v0zw, t1);
	v1xy = _mm_sub_pd(v1xy, t0);
	v1zw = _mm_sub_sd(v1zw, t1);
	t0 = _mm_mul_pd(v0xy, v1xy);
	t1 = _mm_mul_sd(v0zw, v1zw);
	t0 = _mm_add_sd(t0, t1);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	return *reinterpret_cast<double*>(&t1);
}

void HyVoxel::Algebra::VectorSSE::cross(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, double& xout, double& yout, double& zout)
{
	static const __m128d sign = _mm_set_pd(0.0, -0.0);
	__m128d t0, t1;
	const __m128d v0xy = _mm_set_pd(y0, x0);
	const __m128d v0zw = _mm_set_pd(0.0, z0);
	const __m128d v1xy = _mm_set_pd(y1, x1);
	const __m128d v1zw = _mm_set_pd(0.0, z1);
	
	t0 = _mm_mul_pd(_mm_unpacklo_pd(v0zw, v0zw), v1xy);
	t1 = _mm_mul_pd(_mm_unpacklo_pd(v1zw, v1zw), v0xy);
	t1 = _mm_sub_pd(t0, t1);
	t1 = _mm_shuffle_pd(t1, t1, 1);
	t1 = _mm_xor_pd(t1, sign);
	t0 = _mm_mul_pd(v0xy, _mm_shuffle_pd(v1xy, v1xy, 1));
	t0 = _mm_sub_sd(t0, _mm_unpackhi_pd(t0, t0));

	xout = *reinterpret_cast<double*>(&t1);
    
// #ifdef _WIN32
#if 0
	yout = *reinterpret_cast<double*>(&(_mm_unpackhi_pd(t1, t1)));
#else
    __m128d rvalue = _mm_unpackhi_pd(t1, t1);
    yout = *reinterpret_cast<double*>(&rvalue);
#endif
    
	zout = *reinterpret_cast<double*>(&t0);
}

void HyVoxel::Algebra::VectorSSE::cross3(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, double& xout, double& yout, double& zout)
{
	static const __m128d sign = _mm_set_pd(0.0, -0.0);
	__m128d v0xy = _mm_set_pd(y1, x1);
	__m128d v0zw = _mm_set_pd(0.0, z1);
	__m128d v1xy = _mm_set_pd(y2, x2);
	__m128d v1zw = _mm_set_pd(0.0, z2);
	__m128d t0 = _mm_set_pd(y0, x0);
	__m128d t1 = _mm_set_pd(0.0, z0);
	v0xy = _mm_sub_pd(v0xy, t0);
	v0zw = _mm_sub_sd(v0zw, t1);
	v1xy = _mm_sub_pd(v1xy, t0);
	v1zw = _mm_sub_sd(v1zw, t1);

	t0 = _mm_mul_pd(_mm_unpacklo_pd(v0zw, v0zw), v1xy);
	t1 = _mm_mul_pd(_mm_unpacklo_pd(v1zw, v1zw), v0xy);
	t1 = _mm_sub_pd(t0, t1);
	t1 = _mm_shuffle_pd(t1, t1, 1);
	t1 = _mm_xor_pd(t1, sign);
	t0 = _mm_mul_pd(v0xy, _mm_shuffle_pd(v1xy, v1xy, 1));
	t0 = _mm_sub_sd(t0, _mm_unpackhi_pd(t0, t0));

	xout = *reinterpret_cast<double*>(&t1);

// #ifdef _WIN32
#if 0
	yout = *reinterpret_cast<double*>(&(_mm_unpackhi_pd(t1, t1)));
#else
    __m128d rvalue = _mm_unpackhi_pd(t1, t1);
    yout = *reinterpret_cast<double*>(&rvalue);
#endif
    
	zout = *reinterpret_cast<double*>(&t0);
}

double HyVoxel::Algebra::VectorSSE::crossDot(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1)
{
	__m128d t0, t1;
	static const __m128d sign = _mm_set_pd(0.0, -0.0);
	__m128d v0xy = _mm_set_pd(y0, x0);
	__m128d v0zw = _mm_set_pd(0.0, z0);
	__m128d v1xy = _mm_set_pd(y1, x1);
	__m128d v1zw = _mm_set_pd(0.0, z1);

	// cross
	t0 = _mm_mul_pd(_mm_unpacklo_pd(v0zw, v0zw), v1xy);
	t1 = _mm_mul_pd(_mm_unpacklo_pd(v1zw, v1zw), v0xy);
	t1 = _mm_sub_pd(t0, t1);
	t1 = _mm_shuffle_pd(t1, t1, 1);
	t1 = _mm_xor_pd(t1, sign);
	t0 = _mm_mul_pd(v0xy, _mm_shuffle_pd(v1xy, v1xy, 1));
	t0 = _mm_sub_sd(t0, _mm_unpackhi_pd(t0, t0));

	// dot
	t1 = _mm_mul_pd(t1, t1);
	t0 = _mm_mul_sd(t0, t0);
	t0 = _mm_add_sd(t1, t0);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	return *reinterpret_cast<double*>(&t1);
}

double HyVoxel::Algebra::VectorSSE::crossDot3(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	static const __m128d sign = _mm_set_pd(0.0, -0.0);
	__m128d v0xy = _mm_set_pd(y1, x1);
	__m128d v0zw = _mm_set_pd(0.0, z1);
	__m128d v1xy = _mm_set_pd(y2, x2);
	__m128d v1zw = _mm_set_pd(0.0, z2);
	__m128d t0 = _mm_set_pd(y0, x0);
	__m128d t1 = _mm_set_pd(0.0, z0);
	v0xy = _mm_sub_pd(v0xy, t0);
	v0zw = _mm_sub_sd(v0zw, t1);
	v1xy = _mm_sub_pd(v1xy, t0);
	v1zw = _mm_sub_sd(v1zw, t1);

	// cross
	t0 = _mm_mul_pd(_mm_unpacklo_pd(v0zw, v0zw), v1xy);
	t1 = _mm_mul_pd(_mm_unpacklo_pd(v1zw, v1zw), v0xy);
	t1 = _mm_sub_pd(t0, t1);
	t1 = _mm_shuffle_pd(t1, t1, 1);
	t1 = _mm_xor_pd(t1, sign);
	t0 = _mm_mul_pd(v0xy, _mm_shuffle_pd(v1xy, v1xy, 1));
	t0 = _mm_sub_sd(t0, _mm_unpackhi_pd(t0, t0));

	// dot
	t1 = _mm_mul_pd(t1, t1);
	t0 = _mm_mul_sd(t0, t0);
	t0 = _mm_add_sd(t1, t0);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	return *reinterpret_cast<double*>(&t1);
}

double HyVoxel::Algebra::VectorSSE::length(const double x, const double y, const double z)
{
	__m128d t0, t1;
	const __m128d v0xy = _mm_set_pd(y, x);
	const __m128d v0zw = _mm_set_pd(0.0, z);
	t0 = _mm_mul_pd(v0xy, v0xy);
	t1 = _mm_mul_sd(v0zw, v0zw);
	t0 = _mm_add_sd(t0, t1);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	t1 = _mm_sqrt_sd(t1, t1);
	return *reinterpret_cast<double*>(&t1);
}

double HyVoxel::Algebra::VectorSSE::lengthSqr(const double x, const double y, const double z)
{
	__m128d t0, t1;
	const __m128d v0xy = _mm_set_pd(y, x);
	const __m128d v0zw = _mm_set_pd(0.0, z);
	t0 = _mm_mul_pd(v0xy, v0xy);
	t1 = _mm_mul_sd(v0zw, v0zw);
	t0 = _mm_add_sd(t0, t1);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	return *reinterpret_cast<double*>(&t1);
}

void HyVoxel::Algebra::VectorSSE::normalize(const double x, const double y, const double z, double& xout, double& yout, double& zout)
{
	__m128d t0, t1;
	__m128d v0xy = _mm_set_pd(y, x);
	__m128d v0zw = _mm_set_pd(0.0, z);
	t0 = _mm_mul_pd(v0xy, v0xy);
	t1 = _mm_mul_sd(v0zw, v0zw);
	t0 = _mm_add_sd(t0, t1);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	t1 = _mm_sqrt_sd(t1, t1);
	v0xy = _mm_div_pd(v0xy, _mm_unpacklo_pd(t1, t1));
	v0zw = _mm_div_sd(v0zw, t1);

	xout = *reinterpret_cast<double*>(&v0xy);
	
// #ifdef _WIN32
#if 0
	yout = *reinterpret_cast<double*>(&(_mm_unpackhi_pd(v0xy, v0xy)));
#else
    __m128d rvalue = _mm_unpackhi_pd(v0xy, v0xy);
    yout = *reinterpret_cast<double*>(&rvalue);
#endif
    
	zout = *reinterpret_cast<double*>(&v0zw);
}

float HyVoxel::Algebra::invSQR(const float x)
{
	__m128 s = _mm_set_ps(0.f, 0.f, 0.f, x);
	s = _mm_rsqrt_ss(s);
	return *reinterpret_cast<float*>(&s);
}

double HyVoxel::Algebra::invSQR(const double x)
{
	__m128d s = _mm_set_pd(1.0, x);
	s = _mm_sqrt_sd(s, s);
	s = _mm_div_sd(_mm_shuffle_pd(s, s, 1), s);
	return *reinterpret_cast<double*>(&s);
}

void HyVoxel::Algebra::triCentroid(const float x0, const float y0, const float z0, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2, float& xout, float& yout, float& zout)
{
	static const __m128 d = _mm_set_ps1(3.f);
	__m128 v0 = ftov(x0, y0, z0);
	const __m128 v1 = ftov(x1, y1, z1);
	const __m128 v2 = ftov(x2, y2, z2);
	v0 = _mm_add_ps(v0, v1);
	v0 = _mm_add_ps(v0, v2);
	v0 = _mm_div_ps(v0, d);
	vtof(v0, xout, yout, zout);
}

double HyVoxel::Algebra::triArea(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
	static const __m128d sign = _mm_set_pd(0.0, -0.0);
	__m128d v0xy = _mm_set_pd(y1, x1);
	__m128d v0zw = _mm_set_pd(0.0, z1);
	__m128d v1xy = _mm_set_pd(y2, x2);
	__m128d v1zw = _mm_set_pd(0.0, z2);
	__m128d t0 = _mm_set_pd(y0, x0);
	__m128d t1 = _mm_set_pd(0.0, z0);
	v0xy = _mm_sub_pd(v0xy, t0);
	v0zw = _mm_sub_sd(v0zw, t1);
	v1xy = _mm_sub_pd(v1xy, t0);
	v1zw = _mm_sub_sd(v1zw, t1);

	// cross
	t0 = _mm_mul_pd(_mm_unpacklo_pd(v0zw, v0zw), v1xy);
	t1 = _mm_mul_pd(_mm_unpacklo_pd(v1zw, v1zw), v0xy);
	t1 = _mm_sub_pd(t0, t1);
	t1 = _mm_shuffle_pd(t1, t1, 1);
	t1 = _mm_xor_pd(t1, sign);
	t0 = _mm_mul_pd(v0xy, _mm_shuffle_pd(v1xy, v1xy, 1));
	t0 = _mm_sub_sd(t0, _mm_unpackhi_pd(t0, t0));

	// len
	t1 = _mm_mul_pd(t1, t1);
	t0 = _mm_mul_sd(t0, t0);
	t0 = _mm_add_sd(t1, t0);
	t1 = _mm_add_sd(_mm_unpackhi_pd(t0, t0), t0);
	t1 = _mm_sqrt_sd(t1, t1);

	// half
	t0 = _mm_set1_pd(0.5);
	t1 = _mm_mul_sd(t0, t1);
	return *reinterpret_cast<double*>(&t1);
}

#endif