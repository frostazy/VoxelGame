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

#include "HyVoxelPrivatePCH.h"

#include "Common/Vector.h"

#include <math.h>

using namespace HyVoxel;

/*
Vector HyVoxel::Algebra::Vector_withValues(float x, float y, float z) {
	Vector vector;

	vector.x = x;
	vector.y = y;
	vector.z = z;
	return vector;
}
*/

bool HyVoxel::Algebra::Vector_normalize(Vector * vector)
{
	float magnitude;

	magnitude = sqrtf((vector->x * vector->x) + (vector->y * vector->y) + (vector->z * vector->z));
	if (magnitude > 0.0000001f)
	{
		vector->x /= magnitude;
		vector->y /= magnitude;
		vector->z /= magnitude;
		return true;
	}
	else
	{
		return false;
	}
}

Vector HyVoxel::Algebra::Vector_normalized(Vector vector)
{
	Vector_normalize(&vector);
	return vector;
}

float HyVoxel::Algebra::Vector_magnitude(Vector vector)
{
	return sqrtf((vector.x * vector.x) + (vector.y * vector.y) + (vector.z * vector.z));
}

float HyVoxel::Algebra::Vector_magnitudeSquared(Vector vector)
{
	return ((vector.x * vector.x) + (vector.y * vector.y) + (vector.z * vector.z));
}

Vector HyVoxel::Algebra::Vector_add(Vector vector1, Vector vector2)
{
	return Vector_withValues((vector1.x + vector2.x), (vector1.y + vector2.y), (vector1.z + vector2.z));
}

Vector HyVoxel::Algebra::Vector_subtract(Vector vector1, Vector vector2)
{
	return Vector_withValues((vector1.x - vector2.x), (vector1.y - vector2.y), (vector1.z - vector2.z));
}

float HyVoxel::Algebra::Vector_dot(Vector vector1, Vector vector2)
{
	return ((vector1.x * vector2.x) + (vector1.y * vector2.y) + (vector1.z * vector2.z));
}

Vector HyVoxel::Algebra::Vector_cross(Vector vector1, Vector vector2)
{
	Vector result;

	result.x = ((vector1.y * vector2.z) - (vector1.z * vector2.y));
	result.y = ((vector1.z * vector2.x) - (vector1.x * vector2.z));
	result.z = ((vector1.x * vector2.y) - (vector1.y * vector2.x));
	return result;
}

Vector HyVoxel::Algebra::Vector_mutiply(float value, Vector vector)
{
	return Vector_withValues(value*vector.x, value*vector.y, value*vector.z);
}

Vector HyVoxel::Algebra::closestPtPointTriangle(Vector p, Vector a, Vector b, Vector c)
{
	Vector ab = Vector_withValues(b.x - a.x, b.y - a.y, b.z - a.z);
	Vector ac = Vector_withValues(c.x - a.x, c.y - a.y, c.z - a.z);
	Vector bc = Vector_withValues(c.x - b.x, c.y - b.y, c.z - b.z);

	// Compute parametric position s for projection P' of P on AB,
	// P' = A + s*AB, s = snom/(snom+sdenom)
	Vector pa = Vector_subtract(p, a);
	Vector pc = Vector_subtract(p, c);
	float snom = Vector_dot(pa, ab);
	float sdenom = Vector_dot(Vector_subtract(p, b), Vector_subtract(a, b));

	// Compute parametric position t for projection P' of P on AC,
	// P' = A + t*AC, s = tnom/(tnom+tdenom)
	float tnom = Vector_dot(pa, ac), tdenom = Vector_dot(pc, Vector_subtract(a, c));

	if (snom <= 0.0f && tnom <= 0.0f)
	{
		return a;    // Vertex region early out
	}

	// Compute parametric position u for projection P' of P on BC,
	// P' = B + u*BC, u = unom/(unom+udenom)
	float unom = Vector_dot(Vector_subtract(p, b), bc);
	float udenom = Vector_dot(Vector_subtract(p, c), Vector_subtract(b, c));

	if (sdenom <= 0.0f && unom <= 0.0f)
	{
		return b;    // Vertex region early out
	}
	if (tdenom <= 0.0f && udenom <= 0.0f)
	{
		return c;    // Vertex region early out
	}

	// P is outside (or on) AB if the triple scalar product [N PA PB] <= 0
	Vector n = Vector_cross(Vector_subtract(b, a), Vector_subtract(c, a));
	float vc = Vector_dot(n, Vector_cross(Vector_subtract(a, p), Vector_subtract(b, p)));
	// If P outside AB and within feature region of AB,
	// return projection of P onto AB
	if (vc <= 0.0f && snom >= 0.0f && sdenom >= 0.0f)
	{
		return Vector_add(a, Vector_mutiply(snom / (snom + sdenom), ab));
	}

	// P is outside (or on) BC if the triple scalar product [N PB PC] <= 0
	float va = Vector_dot(n, Vector_cross(Vector_subtract(b, p), Vector_subtract(c, p)));
	// If P outside BC and within feature region of BC,
	// return projection of P onto BC
	if (va <= 0.0f && unom >= 0.0f && udenom >= 0.0f)
	{
		return Vector_add(b, Vector_mutiply(unom / (unom + udenom), bc));
	}

	// P is outside (or on) CA if the triple scalar product [N PC PA] <= 0
	float vb = Vector_dot(n, Vector_cross(Vector_subtract(c, p), Vector_subtract(a, p)));
	// If P outside CA and within feature region of CA,
	// return projection of P onto CA
	if (vb <= 0.0f && tnom >= 0.0f && tdenom >= 0.0f)
	{
		return Vector_add(a, Vector_mutiply(tnom / (tnom + tdenom), ac));
	}

	// P must project inside face region. Compute Q using barycentric coordinates
	float u = va / (va + vb + vc);
	float v = vb / (va + vb + vc);
	float w = 1.0f - u - v; // = vc / (va + vb + vc)
	return Vector_add(
			   Vector_mutiply(u, a),
			   Vector_add(
				   Vector_mutiply(v, b),
				   Vector_mutiply(w, c)));
}

// Returns true if sphere s intersects triangle ABC, false otherwise.
// The point p on abc closest to the sphere center is also returned
bool HyVoxel::Algebra::testSphereTriangle(Vector center, float r, Vector a, Vector b, Vector c, Vector &p)
{
	// Find point P on triangle ABC closest to sphere center
	p = closestPtPointTriangle(center, a, b, c);

	// Sphere and triangle intersect if the (squared) distance from sphere
	// center to point p is less than the (squared) sphere radius
	Vector v = Vector_subtract(p, center);
	return Vector_dot(v, v) <= r * r;
}

double HyVoxel::Algebra::sign2D(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y)
{
	return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
}

bool HyVoxel::Algebra::pointInsideTriangle3D(Vector p, Vector v0, Vector v1, Vector v2)
{
	Vector u, v, w;

	// get triangle edge vectors and plane normal
	u = Vector_subtract(v1, v0);
	v = Vector_subtract(v2, v0);

	float  uu, uv, vv, wu, wv, D;
	uu = Vector_dot(u, u);
	uv = Vector_dot(u, v);
	vv = Vector_dot(v, v);
	w = Vector_subtract(p, v0);
	wu = Vector_dot(w, u);
	wv = Vector_dot(w, v);
	D = uv * uv - uu * vv;

	float s, t;
	s = (uv * wv - vv * wu) / D;
	if (s < 0.0 || s > 1.0)        // P is outside T
	{
		return false;
	}
	t = (uv * wu - uu * wv) / D;
	if (t < 0.0 || (s + t) > 1.0)  // P is outside T
	{
		return false;
	}

	return true;                   // P is in T
}

//results
//0 the segment does not intersect the triangle
//1 the segment intersects the triangle
//the segment is in the same plane that the triangle
//2 the segment is inside the triangle
//3 p0 is inside the triangle
//4 p1 is inside the triangle

int HyVoxel::Algebra::segmentIntersectTriangle3D(Vector p0, Vector p1, Vector v0, Vector v1, Vector v2, Vector* p)
{
	Vector u, v, n;     // triangle vectors
	Vector dir, w0, w;  // ray vectors
	float     r, a, b;	// params to calc ray-plane intersect

	// get triangle edge vectors and plane normal
	u = Vector_subtract(v1, v0);
	v = Vector_subtract(v2, v0);
	n = Vector_cross(u, v);

	dir = Vector_subtract(p1, p0);  // ray direction vector
	w0 = Vector_subtract(p0, v0);
	a = -Vector_dot(n,w0);
	b = Vector_dot(n, dir);
	if (fabs(b) < 0.00000001)       // ray is  parallel to triangle plane
	{
		if (a == 0)
		{
			if (pointInsideTriangle3D(p0, v0, v1, v2))
				if (pointInsideTriangle3D(p1, v0, v1, v2))
				{
					return 2;
				}
				else
				{
					return 3;
				}
			else if (pointInsideTriangle3D(p1, v0, v1, v2))
			{
				return 4;
			}
			else
			{
				return 0;
			}
		}
		else
		{
			return 0;    // ray disjoint from plane
		}
	}

	// get intersect point of ray with triangle plane
	r = a / b;
	if (r < 0.0)                    // ray goes away from triangle
	{
		return 0;    // => no intersect
	}
	// for a segment, also test if (r > 1.0) => no intersect

	*p = Vector_add(p0, Vector_mutiply(r, dir));   // intersect point of ray and plane

	// is P inside T?
	float  uu, uv, vv, wu, wv, D;
	uu = Vector_dot(u, u);
	uv = Vector_dot(u, v);
	vv = Vector_dot(v, v);
	w = Vector_subtract(*p, v0);
	wu = Vector_dot(w, u);
	wv = Vector_dot(w, v);
	D = uv * uv - uu * vv;

	float s, t;
	s = (uv * wv - vv * wu) / D;
	if (s < 0.0 || s > 1.0)         // P is outside T
	{
		return 0;
	}
	t = (uv * wu - uu * wv) / D;
	if (t < 0.0 || (s + t) > 1.0)  // P is outside T
	{
		return 0;
	}

	return 1;                       // P is in T
}