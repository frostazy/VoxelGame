/*** ********************************************************************
Matrix class implementation.
Mar 13 2008, HE Zhao
http://hezhao.net
he@hezhao.net
************************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cstring>

namespace HyVoxel
{
	namespace Algebra
	{
		/// A matrix for minimizing Quadratic Error functions
		class QEFMatrix
		{
		public:
			void reset()
			{
				memset(m, 0, 16 * sizeof(double));
			}
			const double* getRawPtr() const
			{
				return this->m;
			}


			QEFMatrix()
			{
				memset(m, 0, 16 * sizeof(double));
			}
			QEFMatrix(double c)
			{
				for (int i = 0; i < 16; i++)
				{
					m[i] = c;
				}
			}
			QEFMatrix(double m11, double m12, double m13, double m14,
					  double m21, double m22, double m23, double m24,
					  double m31, double m32, double m33, double m34,
					  double m41, double m42, double m43, double m44)
			{
				m[0] = m11;
				m[1] = m12;
				m[2] = m13;
				m[3] = m14;
				m[4] = m21;
				m[5] = m22;
				m[6] = m23;
				m[7] = m24;
				m[8] = m31;
				m[9] = m32;
				m[10] = m33;
				m[11] = m34;
				m[12] = m41;
				m[13] = m42;
				m[14] = m43;
				m[15] = m44;
			}
			QEFMatrix(const QEFMatrix& q)
			{
				memcpy(m, q.m, 16 * sizeof(double));
			}

			/*
			* matrix constructor for a specific plane
			*
			* plane - a 4-dimension array
			*/
			QEFMatrix(const double plane[])
			{
				const double& a = plane[0];
				const double& b = plane[1];
				const double& c = plane[2];
				const double& d = plane[3];
				/*m[0]	= a*a;
				m[1]	= a*b;
				m[2]	= a*c;
				m[3]	= a*d;
				m[4]	= b*a;
				m[5]	= b*b;
				m[6]	= b*c;
				m[7]	= b*d;
				m[8]	= c*a;
				m[9]	= c*b;
				m[10]	= c*c;
				m[11]	= c*d;
				m[12]	= d*a;
				m[13]	= d*b;
				m[14]	= d*c;
				m[15]	= d*d;*/
				m[0]	= a*a;
				m[5]	= b*b;
				m[10]	= c*c;
				m[15]	= d*d;
				m[1]	= m[4]	= a*b;
				m[2]	= m[8]	= a*c;
				m[3]	= m[12]	= a*d;
				m[6]	= m[9]	= b*c;
				m[7]	= m[13]	= b*d;
				m[11]	= m[14]	= c*d;
			}

			QEFMatrix(double plane1[], double plane2[])
			{
				double a1 = plane1[0];
				double b1 = plane1[1];
				double c1 = plane1[2];
				double d1 = plane1[3];
				double a2 = plane2[0];
				double b2 = plane2[1];
				double c2 = plane2[2];
				double d2 = plane2[3];
				m[0] = a1*a1 + a2*a2;
				m[1] = a1*b1 + a2*b2;
				m[2] = a1*c1 + a2*c2;
				m[3] = a1*d1 + a2*d2;
				m[4] = a1*b1 + a2*b2;
				m[5] = b1*b1 + b2*b2;
				m[6] = b1*c1 + b2*c2;
				m[7] = b1*d1 + b2*d2;
				m[8] = a1*c1 + a2*c2;
				m[9] = b1*c1 + b2*c2;
				m[10] = c1*c1 + c2*c2;
				m[11] = c1*d1 + c2*d2;
				m[12] = a1*d1 + a2*d2;
				m[13] = b1*d1 + b2*d2;
				m[14] = c1*d1 + c2*d2;
				m[15] = d1*d1 + d2*d2;
			}

			double operator[](int c) const
			{
				return m[c];
			}

			/*
			 * calculate determinant of a 3*3 matrix
			 *
			 * a(mn)  - index of matrix element
			 */
			double det(int a11, int a12, int a13,
					   int a21, int a22, int a23,
					   int a31, int a32, int a33)
			{
				double det = m[a11]*m[a22]*m[a33] + m[a13]*m[a21]*m[a32] + m[a12]*m[a23]*m[a31]
							 - m[a13]*m[a22]*m[a31] - m[a11]*m[a23]*m[a32] - m[a12]*m[a21]*m[a33];
				return det;
			}

			/* For test only */
			void print()
			{
				/*
				for (int i = 0; i < 16; )
				{
					printf("%.8lf ", m[i++]);
					if (i % 4 == 0)
					{
						printf("\n");
					}
				}
				*/
			}

			inline const QEFMatrix operator+(const QEFMatrix& n) const
			{
				return QEFMatrix(
						   m[0]+n[0],   m[1]+n[1],   m[2]+n[2],   m[3]+n[3],
						   m[4]+n[4],   m[5]+n[5],   m[6]+n[6],   m[7]+n[7],
						   m[8]+n[8],   m[9]+n[9],  m[10]+n[10], m[11]+n[11],
						   m[12]+n[12], m[13]+n[13], m[14]+n[14], m[15]+n[15]);
			}

			inline QEFMatrix& operator+=(const QEFMatrix& n)
			{
				m[0]+=n[0];
				m[1]+=n[1];
				m[2]+=n[2];
				m[3]+=n[3];
				m[4]+=n[4];
				m[5]+=n[5];
				m[6]+=n[6];
				m[7]+=n[7];
				m[8]+=n[8];
				m[9]+=n[9];
				m[10]+=n[10];
				m[11]+=n[11];
				m[12]+=n[12];
				m[13]+=n[13];
				m[14]+=n[14];
				m[15]+=n[15];
				return *this;
			}

			inline const QEFMatrix operator+(const double& c) const
			{
				return QEFMatrix(c+m[0],		c+m[1],		c+m[2],		c+m[3],
								 c+m[4],		c+m[5],		c+m[6],		c+m[7],
								 c+m[8],		c+m[9],		c+m[10],	c+m[11],
								 c+m[12],	c+m[13],	c+m[14],	c+m[15]);
			}

			inline QEFMatrix& operator+=(const double& c)
			{
				for (int i = 0; i < 16; i++)
				{
					m[i] += c;
				}
				return *this;
			}

			inline const QEFMatrix operator*(const double& c) const
			{
				return QEFMatrix(c*m[0],		c*m[1],		c*m[2],		c*m[3],
								 c*m[4],		c*m[5],		c*m[6],		c*m[7],
								 c*m[8],		c*m[9],		c*m[10],	c*m[11],
								 c*m[12],	c*m[13],	c*m[14],	c*m[15]);
			}

			inline QEFMatrix& operator*=(const double& c)
			{
				for (int i = 0; i < 16; i++)
				{
					m[i] *= c;
				}
				return *this;
			}

			inline QEFMatrix& operator=(const double& c)
			{
				for (int i = 0; i < 16; i++)
				{
					m[i] = c;
				}
				return *this;
			}

		public:
			double m[16];
		};

		inline double vertex_error(const QEFMatrix& q, const double& x, const double& y, const double& z)
		{
			return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[5]*y*y
				   + 2*q[6]*y*z + 2*q[7]*y + q[10]*z*z + 2*q[11]*z + q[15];
		};

		inline bool solveQEF(QEFMatrix& q_bar, double* vx, double* vy, double* vz, double *error)
		{
			if (q_bar[1] != q_bar[4] || q_bar[2] != q_bar[8] || q_bar[6] != q_bar[9] ||
					q_bar[3] != q_bar[12] || q_bar[7] != q_bar[13] || q_bar[11] != q_bar[14])
			{
				return false;
			}

			QEFMatrix q_delta = QEFMatrix(q_bar[0], q_bar[1],  q_bar[2],  q_bar[3],
										  q_bar[4], q_bar[5],  q_bar[6],  q_bar[7],
										  q_bar[8], q_bar[9], q_bar[10], q_bar[11],
										  0,        0,	      0,        1);

			double det = q_delta.det(0, 1, 2, 4, 5, 6, 8, 9, 10);
			if ((det != 0))
			{
				*vx = -1/det*(q_delta.det(1, 2, 3, 5, 6, 7, 9, 10, 11));
				*vy =  1/det*(q_delta.det(0, 2, 3, 4, 6, 7, 8, 10, 11));
				*vz = -1/det*(q_delta.det(0, 1, 3, 4, 5, 7, 8, 9, 11));
			}
			*error = vertex_error(q_bar, *vx, *vy, *vz);
			return true;
		}
	}

}

#endif
