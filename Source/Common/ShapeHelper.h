// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "HyVoxelLib.h"
#include "HyContainerImpl.h"

#include <vector>
#include <algorithm>

namespace HyVoxel {

	class Edge
	{
	public:
		Edge(int start, int end);

		int start;
		int end;
		Edge Reverse();
		bool operator== (const Edge& edge) const;
	};

	class Face
	{
	public:
		Face(hyvector<int> vi);

		hyvector<int> vi;
		Vector Normal(const hyvector<Vector>& verts) const;
		Vector Center(const hyvector<Vector>& verts) const;
		float MinEdge(const hyvector<Vector>& verts) const;
		hyvector<Edge> GetEdge() const;
		bool valid;
		void Remove();
		bool operator== (const Face &face) const;
	};

	class ShapeHelper
	{
	public:
		ShapeHelper();

		static void AssignUV(const hyvector<Vector>& verts, const hyvector<Face>& faces, hyvector<Vector2f>& UV0, float size);
		static void AssignGlobalUV(const hyvector<Vector>& verts, const hyvector<Face>& faces, hyvector<Vector2f>& uv0, float size);

		static void Scale(hyvector<Vector>& verts, const hyvector<Face>& faces,
			const Vector& center, const Vector& vec);

		static void ScaleFace(hyvector<Vector>& verts, const hyvector<Face>& faces,
			const Face& face, const Vector& vec);

		static hyvector<Face> SubDivideFace(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Face& face, int cuts);

		static void Translate(hyvector<Vector>& verts, const hyvector<Face>& faces,
			const Vector& vec, const hyvector<int>& vi);

		static Face ExtrudeFace(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Face& face, float translate_forwards);

		static Vector ExtrudeConeFace(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Face& face, float translate_forwards);

		static void RemoveFace(hyvector<Face>& faces, const Face& face);
		static Face AddReversedFace(hyvector<Face>& faces, const Face & face);
		static Face CopyFace(hyvector<Vector>& verts, hyvector<Face>& faces, hyvector<Vector2f>& uvs,
			const Face& face, const hyvector<Vector>& old_verts, const hyvector<Vector2f>& old_uvs, float forward);
		//正方体，需要指定 中心，X，Y，Z方向单位向量，半径
		static void AddCube(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Vector& center, const Vector& x, const Vector& y, const Vector& z, float r);
		//圆柱体，需要指定 底面圆中心，高线方向，底面半径，高，圆边数，总边数 num+2
		static void AddCylinder(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Vector& center, const Vector& normal, float radius, float height, int num);
		//圆台
		static void AddAltar(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Vector& center, const Vector& normal, float r_bottom, float r_top, float height, int num);
		//圆锥
		static void AddCone(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Vector& center, const Vector& normal, float radius, float height, int num);
		//球形
		static void AddSphere(hyvector<Vector>& verts, hyvector<Face>& faces,
			const Vector& center, float radius, int num);
		static void AddTriangle(hyvector<Vector>& verts, hyvector<int>& tris,
			const Vector& v0, const Vector& v1, const Vector& v2);
		//返回平面圆形各顶点
		static hyvector<Vector> Circle(const Vector& center, const Vector& normal, float radius, int num);
		//多边形至三角形
		static void ToTris(const hyvector<Vector>& verts, const hyvector<Face>& faces, hyvector<int>& otris);
		//分离所有共顶点面的顶点
		static void SeperateVert(hyvector<Vector>& verts, hyvector<Face>& faces, hyvector<Vector>& normals);
		//磨边
		static void Bevel(hyvector<Vector>& verts, hyvector<Face>& faces, float bevel_size);
	};
	
}
