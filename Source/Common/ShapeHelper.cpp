// Fill out your copyright notice in the Description page of Project Settings.

#include "ShapeHelper.h"

#include <cmath>
#include <algorithm>

using namespace std;

namespace HyVoxel {

	ShapeHelper::ShapeHelper()
	{
	}

	void ShapeHelper::AssignUV(const hyvector<Vector>& verts, const hyvector<Face>& faces, hyvector<Vector2f>& uv0, float size)
	{
		uv0 = hyvector<Vector2f>();
		for (auto vert : verts)
			uv0.push_back(Vector2f(0.0f, 0.0f));
		for (auto face : faces) {
			Vector normal = face.Normal(verts);
			Vector U = (verts[face.vi[1]] - verts[face.vi[0]]).GetSafeNormal(0.01);
			Vector V = Vector::CrossProduct(normal, U).GetSafeNormal(0.01);
			for (int i = 1; i < face.vi.size(); i++) {
				Vector vec = verts[face.vi[i]] - verts[face.vi[0]];
				float u = Vector::DotProduct(vec, U) / size;
				float v = Vector::DotProduct(vec, V) / size;
				uv0[face.vi[i]] = Vector2f(u, v);
			}
		}
	}

	void ShapeHelper::AssignGlobalUV(const hyvector<Vector>& verts, const hyvector<Face>& faces, hyvector<Vector2f>& uv0, float size)
	{
		uv0 = hyvector<Vector2f>();
		for (auto vert : verts)
			uv0.push_back(Vector2f(0.0, 0.0));
		for (auto face : faces) {
			Vector normal = face.Normal(verts);
			for (auto vi : face.vi)
				if (abs(normal.z) >= 0.7)
					uv0[vi] = Vector2f(verts[vi].x, verts[vi].y) / size;
				else if (abs(normal.x) >= abs(normal.y))
					uv0[vi] = Vector2f(verts[vi].y, verts[vi].z) / size;
				else
					uv0[vi] = Vector2f(verts[vi].x, verts[vi].z) / size;
		}
	}



	void ShapeHelper::Scale(hyvector<Vector>& verts, const hyvector<Face>& faces,
		const Vector& center, const Vector& vec)
	{
		for (auto v : verts)
			v = center + (v - center) * vec;
	}

	void ShapeHelper::ScaleFace(hyvector<Vector>& verts, const hyvector<Face>& faces,
		const Face& face, const Vector& vec)
	{
		Vector center = face.Center(verts);
		for (auto vi : face.vi)
			verts[vi] = center + (verts[vi] - center) * vec;
	}

	hyvector<Face> ShapeHelper::SubDivideFace(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Face& face, int cuts)
	{

		hyvector<Face> new_faces;
		if (cuts == 1) {
			new_faces.push_back(face);
			return new_faces;
		}
		int start_vert = verts.size();
		//hyvector<hyvector<Vector>> new_verts;
		for (int i = 0; i < cuts + 1; i++)
			for (int j = 0; j < cuts + 1; j++) {
				Vector v1 = verts[face.vi[0]] + (verts[face.vi[3]] - verts[face.vi[0]]) * i / cuts;
				Vector v2 = verts[face.vi[1]] + (verts[face.vi[2]] - verts[face.vi[1]]) * i / cuts;
				Vector vert = v1 + (v2 - v1) * j / cuts;
				verts.push_back(vert);
			}
		// Remove old one
		int fi;
		if ((fi = faces.Find(face)) != -1) {
			faces[fi].Remove();
		}

		// Divide neighboring edges
		//hyvector<int> neighbor_vi0;
		//hyvector<int> neighbor_vi1;
		//hyvector<int> neighbor_vi2;
		//hyvector<int> neighbor_vi3;
		for (int nfi = 0; nfi < faces.size(); nfi++) {
			if (!faces[nfi].valid)
				continue;
			int nfvi;
			if ((nfvi = faces[nfi].vi.Find(face.vi[0])) != -1) {
				//UE_LOG(LogTemp, Log, TEXT("nvfi:0,%d,%d"), nfi, nfvi)
				//neighbor_vi0.push_back(nfi);
				faces[nfi].vi[nfvi] = start_vert;
			}
			if ((nfvi = faces[nfi].vi.Find(face.vi[1])) != -1) {
				//neighbor_vi1.push_back(nfi);
				faces[nfi].vi[nfvi] = start_vert + cuts;
			}
			if ((nfvi = faces[nfi].vi.Find(face.vi[2])) != -1) {
				//neighbor_vi2.push_back(nfi);
				faces[nfi].vi[nfvi] = start_vert + cuts*(cuts + 2);
			}
			if ((nfvi = faces[nfi].vi.Find(face.vi[3])) != -1) {
				//neighbor_vi3.push_back(nfi);
				faces[nfi].vi[nfvi] = start_vert + cuts*(cuts + 1);
			}
		}
		hyvector<int> f;
		for (int i = 0; i < cuts + 1; i++)
			f.push_back(start_vert + i*(cuts + 1));
		faces.push_back(Face(f));
		f.clear();
		for (int i = 0; i < cuts + 1; i++)
			f.push_back(start_vert + cuts*(cuts + 1) + i);
		faces.push_back(Face(f));
		f.clear();
		for (int i = 0; i < cuts + 1; i++)
			f.push_back(start_vert + cuts*(cuts + 2) - i*(cuts + 1));
		faces.push_back(Face(f));
		f.clear();
		for (int i = 0; i < cuts + 1; i++)
			f.push_back(start_vert + cuts - i);
		faces.push_back(Face(f));
		// Unused
		//for (auto nfi : neighbor_vi0)
		//	if (neighbor_vi3.Find(nfi) != -1) {
		//		int pos = faces[nfi].vi.Find(start_vert);
		//		for (int i = 1; i < cuts; i++)
		//			faces[nfi].vi.Insert(start_vert + i*(cuts + 1), ++pos);
		//	}
		//for (auto nfi : neighbor_vi3)
		//	if (neighbor_vi2.Find(nfi) != -1) {
		//		int pos = faces[nfi].vi.Find(start_vert + cuts*(cuts + 1));
		//		for (int i = 1; i < cuts; i++)
		//			faces[nfi].vi.Insert(start_vert + cuts*(cuts + 1) + i, ++pos);
		//	}
		//for (auto nfi : neighbor_vi2)
		//	if (neighbor_vi1.Find(nfi) != -1) {
		//		int pos = faces[nfi].vi.Find(start_vert + cuts*(cuts + 2));
		//		for (int i = 1; i < cuts; i++)
		//			faces[nfi].vi.Insert(start_vert + cuts*(cuts + 2) - i*(cuts + 1), ++pos);
		//	}
		//for (auto nfi : neighbor_vi1)
		//	if (neighbor_vi0.Find(nfi) != -1) {
		//		int pos = faces[nfi].vi.Find(start_vert + cuts);
		//		for (int i = 1; i < cuts; i++)
		//			faces[nfi].vi.Insert(start_vert + cuts - i, ++pos);
		//	}
		for (int i = 0; i < cuts; i++)
			for (int j = 0; j < cuts; j++) {
				int ff[4] = {
					start_vert + i*(cuts + 1) + j,
					start_vert + i*(cuts + 1) + j + 1,
					start_vert + (i + 1)*(cuts + 1) + j + 1,
					start_vert + (i + 1)*(cuts + 1) + j
				};
				hyvector<int> f;
				f.Append(ff, 4);
				Face new_face = Face(f);
				new_faces.push_back(new_face);
				faces.push_back(new_face);
			}
		return new_faces;
	}

	void ShapeHelper::Translate(hyvector<Vector>& verts, const hyvector<Face>& faces,
		const Vector& vec, const hyvector<int>& vis)
	{
		for (auto vi : vis)
			verts[vi] += vec;
	}

	Face ShapeHelper::ExtrudeFace(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Face& face, float translate_forwards)
	{
		hyvector<int> vi = face.vi;
		// there is an issure that const Face& face may be released in the process ??
		Face f = face;
		int n = vi.size();
		int m = verts.size();
		Vector forward_vec = face.Normal(verts) * translate_forwards;
		hyvector<int> fn;
		for (int i = 0; i < n; i++) {
			verts.push_back(verts[face.vi[i]] + forward_vec);
			fn.push_back(m + i);
		}
		for (int i = 0; i < n - 1; i++) {
			hyvector<int> ff;
			ff.push_back(m + i + 1);
			ff.push_back(m + i);
			ff.push_back(vi[i]);
			ff.push_back(vi[i + 1]);
			faces.push_back(Face(ff));
		}
		hyvector<int> ff;
		ff.push_back(m);
		ff.push_back(m + n - 1);
		ff.push_back(vi[n - 1]);
		ff.push_back(vi[0]);
		faces.push_back(Face(ff));

		int fi;
		if ((fi = faces.Find(f)) != -1) {
			faces[fi].Remove();
		}
		Face new_face = Face(fn);
		faces.push_back(new_face);
		return new_face;
	}

	Vector ShapeHelper::ExtrudeConeFace(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Face& face, float translate_forwards)
	{
		int n = face.vi.size();
		Vector vec = face.Center(verts) + face.Normal(verts) * translate_forwards;
		hyvector<int> f;
		verts.push_back(vec);
		for (int i = 0; i < n - 1; i++) {
			hyvector<int> ff;
			ff.push_back(verts.size() - 1);
			ff.push_back(face.vi[i]);
			ff.push_back(face.vi[i + 1]);
			faces.push_back(Face(ff));
		}
		hyvector<int> ff;
		ff.push_back(verts.size() - 1);
		ff.push_back(face.vi[n - 1]);
		ff.push_back(face.vi[0]);
		faces.push_back(Face(ff));

		int fi;
		if ((fi = faces.Find(face)) != -1) {
			faces[fi].Remove();
		}
		Face new_face = Face(f);
		faces.push_back(new_face);
		return vec;
	}

	void ShapeHelper::RemoveFace(hyvector<Face>& faces, const Face& face)
	{
		int fi;
		if ((fi = faces.Find(face)) != -1) {
			faces[fi].Remove();
		}
	}

	Face ShapeHelper::AddReversedFace(hyvector<Face>& faces, const Face & face)
	{
		hyvector<int> fvi;
		for (int i = face.vi.size() - 1; i >= 0; i--)
			fvi.push_back(face.vi[i]);
		Face newFace = Face(fvi);
		faces.push_back(newFace);
		return newFace;
	}

	Face ShapeHelper::CopyFace(hyvector<Vector>& verts, hyvector<Face>& faces, hyvector<Vector2f>& uvs,
		const Face& face, const hyvector<Vector>& old_verts, const hyvector<Vector2f>& old_uvs, float forward)
	{
		bool add_uvs = (old_uvs.size() != 0);
		hyvector<int> fvi;
		for (auto vi : face.vi) {
			int i;
			if ((i = verts.Find(old_verts[vi])) != -1 && forward == 0.0) {
				fvi.push_back(i);
				if (add_uvs)
					uvs[i] = old_uvs[vi];
			}
			else {
				fvi.push_back(verts.size());
				if (forward == 0.0)
					verts.push_back(old_verts[vi]);
				else
					verts.push_back(old_verts[vi] + face.Normal(old_verts)*forward);
				if (add_uvs)
					uvs.push_back(old_uvs[vi]);
			}
		}
		Face new_face = Face(fvi);
		faces.push_back(new_face);
		return new_face;
	}

	void ShapeHelper::AddCube(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Vector& center, const Vector& x, const Vector& y, const Vector& z, float r)
	{
		int start_vert = verts.size();
		verts.push_back(center + (-x - y - z)*r);
		verts.push_back(center + (x - y - z)*r);
		verts.push_back(center + (x - y + z)*r);
		verts.push_back(center + (-x - y + z)*r);
		verts.push_back(center + (-x + y - z)*r);
		verts.push_back(center + (x + y - z)*r);
		verts.push_back(center + (x + y + z)*r);
		verts.push_back(center + (-x + y + z)*r);
		const int FACES[6][4] = {
			{ start_vert + 2, start_vert + 6, start_vert + 5, start_vert + 1 },
			{ start_vert + 7, start_vert + 3, start_vert + 0, start_vert + 4 },
			{ start_vert + 3, start_vert + 2, start_vert + 1, start_vert + 0 },
			{ start_vert + 6, start_vert + 7, start_vert + 4, start_vert + 5 },
			{ start_vert + 7, start_vert + 6, start_vert + 2, start_vert + 3 },
			{ start_vert + 0, start_vert + 1, start_vert + 5, start_vert + 4 }
		};
		for (int i = 0; i < 6; i++) {
			hyvector<int> f;
			f.Append(FACES[i], 4);
			faces.push_back(Face(f));
		}
	}


	void ShapeHelper::AddCylinder(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Vector& center, const Vector& normal, float radius, float height, int num)
	{
		bool with_tops = true;
		int start_vert = verts.size();
		hyvector<Vector> verts_c = Circle(center, normal, radius, num);
		verts += verts_c;
		verts.push_back(verts_c[0] + normal * height);
		//push_back side faces
		for (int i = 0; i < num - 1; i++) {
			verts.push_back(verts_c[i + 1] + normal * height);
			int FACES[4] = {
				start_vert + num + i,
				start_vert + num + i + 1,
				start_vert + i + 1,
				start_vert + i
			};
			hyvector<int> f;
			f.Append(FACES, 4);
			faces.push_back(Face(f));
		}
		int FACES[4] = {
			start_vert + 2 * num - 1,
			start_vert + num,
			start_vert + 0,
			start_vert + num - 1
		};
		hyvector<int> f;
		f.Append(FACES, 4);
		faces.push_back(Face(f));
		if (with_tops) {
			//push_back bottom face
			f.clear();
			for (int i = 0; i < num; i++)
				f.push_back(start_vert + i);
			faces.push_back(Face(f));
			//push_back top face
			f.clear();
			for (int i = num - 1; i >= 0; i--)
				f.push_back(start_vert + num + i);
			faces.push_back(Face(f));
		}
	}

	void ShapeHelper::AddAltar(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Vector& center, const Vector& normal, float r_bottom, float r_top, float height, int num)
	{
		bool with_tops = true;
		int start_vert = verts.size();
		verts += Circle(center, normal, r_bottom, num);
		verts += Circle(center + normal * height, normal, r_top, num);
		//push_back side faces
		for (int i = 0; i < num - 1; i++) {
			int FACES[4] = {
				start_vert + num + i,
				start_vert + num + i + 1,
				start_vert + i + 1,
				start_vert + i
			};
			hyvector<int> f;
			f.Append(FACES, 4);
			faces.push_back(Face(f));
		}
		int FACES[4] = {
			start_vert + 2 * num - 1,
			start_vert + num,
			start_vert + 0,
			start_vert + num - 1
		};
		hyvector<int> f;
		f.Append(FACES, 4);
		faces.push_back(Face(f));
		if (with_tops) {
			//push_back bottom face
			f.clear();
			for (int i = 0; i < num; i++)
				f.push_back(start_vert + i);
			faces.push_back(Face(f));
			//push_back top face
			f.clear();
			for (int i = num - 1; i >= 0; i--)
				f.push_back(start_vert + num + i);
			faces.push_back(Face(f));
		}
	}

	void ShapeHelper::AddCone(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Vector& center, const Vector& normal, float radius, float height, int num)
	{
		bool with_tops = true;
		int start_vert = verts.size();
		verts += Circle(center, normal, radius, num);
		verts.push_back(center + normal * height);
		//push_back side faces
		for (int i = 0; i < num - 1; i++) {
			int FACES[3] = {
				start_vert + num,
				start_vert + i + 1,
				start_vert + i
			};
			hyvector<int> f;
			f.Append(FACES, 3);
			faces.push_back(Face(f));
		}
		int FACES[3] = {
			start_vert + num,
			start_vert + 0,
			start_vert + num - 1
		};
		hyvector<int> f;
		f.Append(FACES, 3);
		faces.push_back(Face(f));
		if (with_tops) {
			//push_back bottom face
			f.clear();
			for (int i = 0; i < num; i++)
				f.push_back(start_vert + i);
			faces.push_back(Face(f));
		}
	}

	void ShapeHelper::AddSphere(hyvector<Vector>& verts, hyvector<Face>& faces,
		const Vector& center, float radius, int num)
	{
		int start_vert = verts.size();
		Vector normal = Vector(0.0, 0.0, 1.0);
		verts += Circle(center, normal, radius, 3 * num);
		hyvector<int> f1, f2;
		for (int i = 0; i < 3 * num; i++) {
			f1.push_back(start_vert + i);
			f2.push_back(start_vert + 3 * num - 1 - i);
		}
		Face face1 = Face(f1);
		Face face2 = Face(f2);
		faces.push_back(face1);
		faces.push_back(face2);
		float arc = 3.141592654 / 2 / (float)num;
		for (int i = 0; i < num - 1; i++) {
			float h = (sin((float)(i + 1)*arc) - sin((float)i*arc))*radius;
			float r = cos((float)(i + 1)*arc) / cos((float)i*arc);
			face1 = ExtrudeFace(verts, faces, face1, h);
			ScaleFace(verts, faces, face1, Vector(r, r, r));
			face2 = ExtrudeFace(verts, faces, face2, h);
			ScaleFace(verts, faces, face2, Vector(r, r, r));
		}
		float h = (1 - cos(arc))*radius;
		ExtrudeConeFace(verts, faces, face1, h);
		ExtrudeConeFace(verts, faces, face2, h);
	}

	void ShapeHelper::AddTriangle(hyvector<Vector>& verts, hyvector<int>& tris,
		const Vector& v0, const Vector& v1, const Vector& v2)
	{
		verts.clear();
		verts.push_back(v0);
		verts.push_back(v1);
		verts.push_back(v2);

		tris.clear();
		tris.push_back(0);
		tris.push_back(1);
		tris.push_back(2);
	}

	hyvector<Vector> ShapeHelper::Circle(const Vector& center, const Vector& normal, float radius, int num)
	{
		Vector cnormal = normal.GetSafeNormal(0.01);
		Vector X = Vector(1.0, 0.0, 0.0);
		Vector Y = Vector(0.0, 1.0, 0.0);
		Vector::CreateOrthonormalBasis(X, Y, cnormal);
		//make sure normal == X cross Y
		if (Vector::DotProduct(Vector::CrossProduct(X, Y), cnormal) < 0)
			std::swap(X, Y);
		float arc = 2 * 3.141592654 / num;
		hyvector<Vector> points;
		for (int i = 0; i < num; i++) {
			Vector p = center + cos((float)i*arc)*radius*X + sin((float)i*arc)*radius*Y;
			points.push_back(p);
			//UE_LOG(LogTemp, Log, TEXT("p:%f,%f,%f"), p.X, p.Y, p.Z)
		}
		return points;
	}

	void ShapeHelper::ToTris(const hyvector<Vector>& verts, const hyvector<Face>& faces, hyvector<int>& otris)
	{
		hyvector<Face> new_faces;
		for (auto face : faces) {
			if (face.valid && face.vi.size() >= 3) {
				for (int i = 0; i < face.vi.size() - 2; i++) {
					int TRIS[3] = {
						face.vi[0], face.vi[i + 1], face.vi[i + 2],
					};
					hyvector<int> f;
					f.Append(TRIS, 3);
					new_faces.push_back(Face(f));
				}
			}
		}
		otris.clear();
		for (auto face : new_faces) {
			otris.push_back(face.vi[0]);
			otris.push_back(face.vi[1]);
			otris.push_back(face.vi[2]);
		}
	}

	void ShapeHelper::SeperateVert(hyvector<Vector>& verts, hyvector<Face>& faces, hyvector<Vector>& normals)
	{
		hyvector<Face> new_faces;
		hyvector<Vector> new_verts;
		normals.clear();
		for (auto face : faces) {
			if (face.valid) {
				hyvector<int> f;
				for (auto vi : face.vi) {
					f.push_back(new_verts.size());
					new_verts.push_back(verts[vi]);
					normals.push_back(face.Normal(verts));
				}
				new_faces.push_back(Face(f));
			}
		}
		verts = new_verts;
		faces = new_faces;
	}

	void ShapeHelper::Bevel(hyvector<Vector>& verts, hyvector<Face>& faces, float bevel_ratio)
	{
		hyvector<hyvector<Edge>> neighbor_e_of_v;
		hyvector<Edge> edges;
		hyvector<hyvector<Edge>> new_edges_of_e;
		int start_vert = verts.size();
		int start_face = faces.size();
		for (auto v : verts) {
			neighbor_e_of_v.push_back(hyvector<Edge>());
		}

		//UE_LOG(LogTemp, Log, TEXT("bevel:%d,%d"), start_vert, start_face)

		for (int fi = 0; fi < start_face; fi++) {
			Face face = faces[fi];
			if (!face.valid)
				continue;

			// Add a smaller faces along original one
			hyvector<int> f;
			for (int i = 0; i < face.vi.size(); i++)
				f.push_back(start_vert + i);
			Face new_face = Face(f);
			faces.push_back(new_face);

			// Record the transformed edge of original one
			// Original edge is (edge[i].start, edge[i].end), new edge is (start_vert+i, start_vert+i+1)
			hyvector<Edge> fe = face.GetEdge();
			for (int i = 0; i < face.vi.size(); i++) {
				int ei;
				Edge e = fe[i];
				if ((ei = edges.Find(e)) != -1) {
					continue;
				}
				else if ((ei = edges.Find(e.Reverse())) != -1) {
					new_edges_of_e[ei].push_back(Edge(start_vert + i, start_vert + (i + 1) % face.vi.size()));
				}
				else {
					edges.push_back(e);
					hyvector<Edge> es;
					new_edges_of_e.push_back(es);
					new_edges_of_e.back().push_back(Edge(start_vert + i, start_vert + (i + 1) % face.vi.size()));
					neighbor_e_of_v[e.start].push_back(e);
					neighbor_e_of_v[e.end].push_back(e);
				}
			}

			// Add vertexes for new faces
			Vector center = face.Center(verts);
			for (auto vi : face.vi) {
				Vector new_vert = verts[vi] + (center - verts[vi])*bevel_ratio;
				verts.push_back(new_vert);
			}
			start_vert = verts.size();
		}

		// Add faces along original edges
		for (int ei = 0; ei < edges.size(); ei++) {
			if (new_edges_of_e[ei].size() < 2)
				continue;
			hyvector<int> f;
			f.push_back(new_edges_of_e[ei][1].end);
			f.push_back(new_edges_of_e[ei][1].start);
			f.push_back(new_edges_of_e[ei][0].end);
			f.push_back(new_edges_of_e[ei][0].start);
			Face new_face = Face(f);
			faces.push_back(new_face);
		}

		// Add faces along original vertexes
		for (int vi = 0; vi < neighbor_e_of_v.size(); vi++) {
			Vector center_vert = Vector(0.0, 0.0, 0.0);
			if (neighbor_e_of_v[vi].size() < 3)
				//if (neighbor_e_of_v[vi].size() < 2)		
				continue;
			for (auto e : neighbor_e_of_v[vi]) {
				int ei = edges.Find(e);
				//UE_LOG(LogTemp, Log, TEXT("bevel:%d,%d,%d"), vi, neighbor_e_of_v.size(), ei)
				if (ei == -1)
					continue;
				if (vi == e.start)
					center_vert += verts[new_edges_of_e[ei][0].start];
				else
					center_vert += verts[new_edges_of_e[ei][0].end];
			}
			center_vert /= neighbor_e_of_v[vi].size();
			verts.push_back(center_vert);
			int new_vi = verts.size() - 1;
			for (auto e : neighbor_e_of_v[vi]) {
				hyvector<int> f;
				f.push_back(new_vi);
				int ei = edges.Find(e);
				if (ei == -1)
					continue;
				if (new_edges_of_e[ei].size() < 2)
					continue;
				if (vi == e.start) {
					f.push_back(new_edges_of_e[ei][1].end);
					f.push_back(new_edges_of_e[ei][0].start);
				}
				else {
					f.push_back(new_edges_of_e[ei][0].end);
					f.push_back(new_edges_of_e[ei][1].start);
				}
				Face new_face = Face(f);
				faces.push_back(new_face);
			}
		}

		for (int i = 0; i < start_face; i++)
			faces[i].Remove();
	}




	Face::Face(hyvector<int> vi)
	{
		this->vi = vi;
		this->valid = true;
	}

	Vector Face::Normal(const hyvector<Vector>& verts) const
	{
		if (vi.size() > 2)
			return Vector::CrossProduct(verts[vi[2]] - verts[vi[0]],
				verts[vi[1]] - verts[vi[0]]).GetSafeNormal(0.01);
		else return Vector(0.0, 0.0, 1.0);
	}

	Vector Face::Center(const hyvector<Vector>& verts) const
	{
		Vector center = Vector(0.0, 0.0, 0.0);
		for (int i = 0; i < vi.size(); i++)
			center += verts[vi[i]];
		center /= vi.size();
		return center;
	}

	float Face::MinEdge(const hyvector<Vector>& verts) const
	{
		float length = Vector::Dist(verts[vi[0]], verts[vi[vi.size() - 1]]);
		for (int i = 0; i < vi.size() - 1; i++)
			if (Vector::Dist(verts[vi[i]], verts[vi[i + 1]]) < length)
				length = Vector::Dist(verts[vi[i]], verts[vi[i + 1]]) < length;
		return length / 2;
	}

	hyvector<Edge> Face::GetEdge() const
	{
		hyvector<Edge> edges;
		for (int i = 0; i < vi.size(); i++)
			edges.push_back(Edge(vi[i], vi[(i + 1) % vi.size()]));
		return edges;
	}

	void Face::Remove()
	{
		valid = false;
	}

	bool Face::operator==(const Face& face) const
	{
		return this->vi == face.vi;
	}

	Edge::Edge(int start, int end)
	{
		this->start = start;
		this->end = end;
	}

	Edge Edge::Reverse()
	{
		return Edge(end, start);
	}

	bool Edge::operator==(const Edge& edge) const
	{
		return (this->start == edge.start) && (this->end == edge.end);
	}

}


// TODO: liang
//  to be moved into somewhere like help function interface.
#include "ProcModeling/Spaceship/SpaceshipGenerator.h"
namespace HyVoxel
{
	void HySpaceshipImpl::GenerateSphere(const Vector& center, float radius, int num, OutputMesh& omesh)
	{
		OutVectorImpl<Vector>* verts = new OutVectorImpl<Vector>();
		OutVectorImpl<int>* tris = new OutVectorImpl<int>();

		hyvector<Face> faces;
		ShapeHelper::AddSphere(*verts, faces, center, radius, num);

		ShapeHelper::ToTris(*verts, faces, *tris);

		verts->SetOutVectorValues();
		tris->SetOutVectorValues();
		omesh.verts = verts;
		omesh.normals = NULL;
		omesh.tris = tris;
	}

	void HySpaceshipImpl::GenerateTriangle(const Vector& v0, const Vector& v1, const Vector& v2, OutputMesh& omesh)
	{
		OutVectorImpl<Vector>* verts = new OutVectorImpl<Vector>();
		OutVectorImpl<int>* tris = new OutVectorImpl<int>();

		ShapeHelper::AddTriangle(*verts, *tris, v0, v1, v2);

		verts->SetOutVectorValues();
		tris->SetOutVectorValues();
		omesh.verts = verts;
		omesh.normals = NULL;
		omesh.tris = tris;
	}
}