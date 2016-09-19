// Fill out your copyright notice in the Description page of Project Settings.

#include "SpaceshipGenerator.h"

#include <time.h>  
#include <cmath> 

namespace HyVoxel {

	float Random(float start, float end)
	{
		return start + (end - start) * rand() / float(RAND_MAX);
	}

	float GetAspectRatio(const Face& face, hyvector<Vector>verts)
	{
		if (face.vi.size() < 4)
			return 1000.0;
		float distMax = Vector::Dist(verts[face.vi[0]], verts[face.vi[face.vi.size() - 1]]);
		float distMin = Vector::Dist(verts[face.vi[0]], verts[face.vi[face.vi.size() - 1]]);
		for (int i = 0; i < face.vi.size() - 1; i++) {
			float dist = Vector::Dist(verts[face.vi[i]], verts[face.vi[i + 1]]);
			if (dist > distMax)
				distMax = dist;
			if (dist < distMin)
				distMin = dist;
		}
		return distMax / distMin;
	}

	// Extrudes a face along its normal by translate_forwards units.
	// Returns the new face, and optionally fills out extruded_face_list
	// with all the additional side faces created from the extrusion.


	// Similar to extrude_face, except corrigates the geometry to create "ribs".
	// Returns the new face.
	Face RibbedExtrudeFace(hyvector<Vector>& verts, hyvector<Face>& faces, const Face& face,
		float forwards, int ribs = 3, float scale = 0.9)
	{
		if (faces.Find(face) == -1)
			return face;
		float forwardsPerRib = forwards / float(ribs);
		Face newFace = face;
		for (int i = 0; i < ribs; i++) {
			newFace = ShapeHelper::ExtrudeFace(verts, faces, newFace, forwardsPerRib * 0.25);
			newFace = ShapeHelper::ExtrudeFace(verts, faces, newFace, 0.0);
			ShapeHelper::ScaleFace(verts, faces, newFace, Vector(scale, scale, scale));
			newFace = ShapeHelper::ExtrudeFace(verts, faces, newFace, forwardsPerRib * 0.5);
			newFace = ShapeHelper::ExtrudeFace(verts, faces, newFace, 0.0);
			ShapeHelper::ScaleFace(verts, faces, newFace, Vector(1 / scale, 1 / scale, 1 / scale));
			newFace = ShapeHelper::ExtrudeFace(verts, faces, newFace, forwardsPerRib * 0.25);
		}
		return newFace;
	}

	void ShipCabinImpl::InitCore(ShipComponentImpl& component, Vector shipCenter, float shipCoreSize)
	{
		ShapeHelper::AddCube(*component.verts, component.faces, shipCenter,
			Vector(1.0, 0.0, 0.0), Vector(0.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0), shipCoreSize / 2);
		type = MAIN;
		length = shipCoreSize;
		fi[FRONT] = 0;
		fi[BACK] = 1;
		fi[LEFT] = 2;
		fi[RIGHT] = 3;
		fi[TOP] = 4;
		fi[BOTTOM] = 5;
		vi[0] = 6;
		vi[1] = 2;
		vi[2] = 1;
		vi[3] = 5;
		vi[4] = 7;
		vi[5] = 3;
		vi[6] = 0;
		vi[7] = 4;
		for (auto ft : fType)
			ft = UNSPECIFIED;
	}

	void ShipCabinImpl::Extrude(ShipComponentImpl& component, ShipCabinImpl& newCabin, Direction dir,
		float cabinLength, float sizeControl, bool useRibbing)
	{
		newCabin.type = ShipCabinImpl::MAIN;
		if (this->successors[0] == NULL)
			this->successors[0] = &newCabin;
		else
			this->successors[1] = &newCabin;
		newCabin.predecessor = this;

		if (useRibbing && type != RIBBED && (dir == FRONT || dir == BACK) && sizeControl > 1.0 &&
			Random(0.0, 1.0) < 0.1) {
			this->fType[dir] = INVALID;
			float scale = Random(0.7, 0.9);
			int ribs = Random(4, 6);
			RibbedExtrudeFace(*component.verts, component.faces, component.faces[(this->fi[dir])], 2 * cabinLength, ribs, scale);

			int fnum = component.faces.size();
			newCabin.fi[dir] = fnum - 1;
			newCabin.fType[dir] = UNSPECIFIED;
			Face& f1 = component.faces[this->fi[dir]];
			Face& f2 = component.faces[newCabin.fi[dir]];
			newCabin.vi[0] = f1.vi[0];
			newCabin.vi[1] = f1.vi[1];
			newCabin.vi[2] = f1.vi[2];
			newCabin.vi[3] = f1.vi[3];
			newCabin.vi[4] = f2.vi[0];
			newCabin.vi[5] = f2.vi[1];
			newCabin.vi[6] = f2.vi[2];
			newCabin.vi[7] = f2.vi[3];
			newCabin.length = 2 * cabinLength;
			return;
		}

		this->fType[dir] = INVALID;
		ShapeHelper::ExtrudeFace(*component.verts, component.faces, component.faces[(this->fi[dir])], cabinLength);

		int fnum = component.faces.size();
		newCabin.fi[dir] = fnum - 1;
		newCabin.fType[dir] = UNSPECIFIED;
		Face& f1 = component.faces[this->fi[dir]];
		Face& f2 = component.faces[newCabin.fi[dir]];
		newCabin.vi[0] = f1.vi[0];
		newCabin.vi[1] = f1.vi[1];
		newCabin.vi[2] = f1.vi[2];
		newCabin.vi[3] = f1.vi[3];
		newCabin.vi[4] = f2.vi[0];
		newCabin.vi[5] = f2.vi[1];
		newCabin.vi[6] = f2.vi[2];
		newCabin.vi[7] = f2.vi[3];
		newCabin.length = cabinLength;
		int fis[4];
		fis[0] = fnum - 2;
		fis[1] = fnum - 3;
		fis[2] = fnum - 4;
		fis[3] = fnum - 5;
		if (dir == FRONT || dir == BACK) {
			for (auto& fi : fis) {
				if (component.faces[fi].Normal(*component.verts).y < -0.1) {
					newCabin.fType[LEFT] = UNSPECIFIED;
					newCabin.fi[LEFT] = fi;
				}
				else if (component.faces[fi].Normal(*component.verts).y > 0.1) {
					newCabin.fType[RIGHT] = UNSPECIFIED;
					newCabin.fi[RIGHT] = fi;
				}
				else if (component.faces[fi].Normal(*component.verts).z > 0) {
					newCabin.fType[TOP] = UNSPECIFIED;
					newCabin.fi[TOP] = fi;
				}
				else if (component.faces[fi].Normal(*component.verts).z < 0) {
					newCabin.fType[BOTTOM] = UNSPECIFIED;
					newCabin.fi[BOTTOM] = fi;
				}
			}
		}
		else {
			for (auto& fi : fis) {
				if (component.faces[fi].Normal(*component.verts).x < -0.2) {
					newCabin.fType[BACK] = UNSPECIFIED;
					newCabin.fi[BACK] = fi;
				}
				else if (component.faces[fi].Normal(*component.verts).x > 0.2) {
					newCabin.fType[FRONT] = UNSPECIFIED;
					newCabin.fi[FRONT] = fi;
				}
			}
		}
		return;
	}

	float ShipCabinImpl::RandomAdjust(ShipComponentImpl& component, Direction dir, float sizeControl, float sizeExpect,
		bool useScaling, bool useTranslation, bool useRotation)
	{
		if (type == RIBBED)
			return sizeControl;

		hyvector<int>& vi = component.faces[fi[dir]].vi;

		// Maybe apply some scaling
		if (Random(0.0, 1.0) > 0.3) {
			float s = Random(1.3, 1.5);
			if (sizeControl > sizeExpect * 1.1)
				s = 1 / s;
			else if (sizeControl > sizeExpect / 1.1 && Random(0.0, 1.0) > 0.5)
				s = 1 / s;
			sizeControl *= s;
			float sy = s * Random(0.9, 1.1);
			float sz = s * Random(0.9, 1.1);
			ShapeHelper::ScaleFace(*component.verts, component.faces, component.faces[fi[dir]], Vector(1, sy, sz));
		}

		// Maybe apply some sideways translation
		if (Random(0.0, 1.0) > 0.3) {
			Vector sideways = Vector(0.0, 0.0, Random(0.1, 0.5) * length);
			if (Random(0.0, 1.0) > 0.5)
				sideways = -sideways;
			hyvector<int> vis;
			vis.push_back(vi[0]);
			vis.push_back(vi[1]);
			vis.push_back(vi[2]);
			vis.push_back(vi[3]);
			ShapeHelper::Translate(*component.verts, component.faces, sideways * sizeControl, vis);
		}

		// Maybe add some rotation around Y axis
		float angle = 0.0;
		if (Random(0.0, 1.0) > 0.3)
			angle = 5.0;
		if (Random(0.0, 1.0) > 0.7)
			angle = 10.0;
		if (Random(0.0, 1.0) > 0.5)
			angle = -angle;
		(*component.verts)[vi[0]] = (*component.verts)[vi[0]].RotateAngleAxis(angle * sizeControl, Vector(0.0, 1.0, 0.0));
		(*component.verts)[vi[1]] = (*component.verts)[vi[1]].RotateAngleAxis(angle * sizeControl, Vector(0.0, 1.0, 0.0));
		(*component.verts)[vi[2]] = (*component.verts)[vi[2]].RotateAngleAxis(angle * sizeControl, Vector(0.0, 1.0, 0.0));
		(*component.verts)[vi[3]] = (*component.verts)[vi[3]].RotateAngleAxis(angle * sizeControl, Vector(0.0, 1.0, 0.0));

		return sizeControl;
	}

	void ShipCabinImpl::Draw(ShipComponentImpl& component, const hyvector<Vector>& refVerts, float scale, float margin)
	{
		int vnum = (*component.verts).size();
		int fnum = component.faces.size();
		hyvector<Vector> nv;
		for (auto& i : this->vi)
			nv.push_back(refVerts[i]);
		// Scaling
		Vector center(0.0, 0.0, 0.0);
		for (auto& v : nv)
			center += v;
		center /= 8;
		for (auto& v : nv) {
			v = center*(1.0 - scale) + v*scale;
			(*component.verts).push_back(v);
		}
		const int FACES[6][4] = {
			{ vnum + 3, vnum + 2, vnum + 1, vnum + 0 },
			{ vnum + 4, vnum + 5, vnum + 6, vnum + 7 },
			{ vnum + 0, vnum + 1, vnum + 5, vnum + 4 },
			{ vnum + 1, vnum + 2, vnum + 6, vnum + 5 },
			{ vnum + 2, vnum + 3, vnum + 7, vnum + 6 },
			{ vnum + 3, vnum + 0, vnum + 4, vnum + 7 }
		};
		for (int i = 0; i < 6; i++) {
			hyvector<int> f;
			f.Append(FACES[i], 4);
			component.faces.push_back(Face(f));
		}
	}

	SpaceshipGenerator::SpaceshipGenerator(int num_hull_segments_min,
		int num_hull_segments_max,
		bool create_asymmetry_segments,
		int num_asymmetry_segments_min,
		int num_asymmetry_segments_max,
		bool create_face_detail,
		bool apply_bevel_modifier)
		: num_hull_segments_min(num_hull_segments_min)
		, num_hull_segments_max(num_hull_segments_max)
		, num_asymmetry_segments_min(num_asymmetry_segments_min)
		, num_asymmetry_segments_max(num_asymmetry_segments_max)
		, create_face_detail(create_face_detail)
		, apply_bevel_modifier(apply_bevel_modifier)
	{}

	bool SpaceshipGenerator::Generate(Vector center, Vector rotation, float size, unsigned int seed)
	{
		for (int i = 0; i < Spaceship::COMPONENT_COUNT; ++i)
		{
			ShipComponentImpl& comp = components[i];
			comp.verts = new OutVectorImpl<Vector>();
			comp.normals = new OutVectorImpl<Vector>();
			comp.tris = new OutVectorImpl<int>();
			comp.uvs = new OutVectorImpl<Vector2f>();
		}

		shipCoreSize = size;
		shipCenter = center;

		if (seed == 0)
			seed = (unsigned int)time(NULL);
		srand(seed);
		//UE_LOG(LogTemp, Log, TEXT("start"));

		ShipComponentImpl& hulls = components[HULL];
		ShipComponentImpl& accessories = components[ACCESSORY];
		ShipComponentImpl& spheres = components[HULL_SPHERE];
		ShipComponentImpl& engines = components[ENGINE_EMIT];
		ShipComponentImpl& lights = components[LIGHT];
		ShipComponentImpl& windows = components[WINDOW_EMIT];
		ShipComponentImpl& internals = components[CABIN];

		// Let's start with a unit BMesh cube scaled randomly
		cabins[cabinNum++].InitCore(hulls, shipCenter, shipCoreSize);
		cabinNum++;


		// Extrude out the hull along the X axis, adding some semi - random perturbations
		int frontSegments = rand() % (num_hull_segments_max - num_hull_segments_min + 1) + num_hull_segments_min;
		int backSegments = rand() % (num_hull_segments_max - num_hull_segments_min + 1) + num_hull_segments_min;

		// Extrude rear face first
		int refc = 0;
		float cabinLength = Random(0.4, 0.7) * shipCoreSize;
		float sizeControl = 1.0;
		for (int j = 0; j < backSegments - 1; j++) {

			// Most of the time, extrude out the face with some random deviations
			cabins[refc].Extrude(hulls, cabins[cabinNum], BACK, cabinLength * sizeControl, sizeControl);
			refc = cabinNum;
			cabinNum++;
			sizeControl = cabins[refc].RandomAdjust(hulls, BACK, sizeControl, 1.2);
		}

		// The rear cabin, make it large enough
		float scaleBack = Random(0.7, 1.3);
		float scaleBackEnd = Random(0.7, 1.0);
		ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, hulls.faces[cabins[refc].fi[BACK]], Vector(1, scaleBack, scaleBack));
		cabinLength = Random(0.3, 0.8) * cabinLength * sizeControl;
		cabins[refc].Extrude(hulls, cabins[cabinNum], BACK, cabinLength * sizeControl, sizeControl, false);
		refc = cabinNum;
		cabinNum++;
		cabins[refc].type = ShipCabinImpl::BACK_END;

		cabins[refc].Extrude(hulls, cabins[cabinNum], BACK, cabinLength * sizeControl, sizeControl, false);
		refc = cabinNum;
		cabinNum++;
		ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, hulls.faces[cabins[refc].fi[BACK]],
			Vector(1, scaleBackEnd, scaleBackEnd));
		cabins[refc].type = ShipCabinImpl::BACK_END;

		refc = 0;
		cabinLength = Random(0.4, 0.7) * shipCoreSize;
		sizeControl = 1.0;
		for (int j = 0; j < frontSegments - 1; j++) {

			// Most of the time, extrude out the face with some random deviations
			cabins[refc].Extrude(hulls, cabins[cabinNum], FRONT, cabinLength * sizeControl, sizeControl);
			refc = cabinNum;
			cabinNum++;
			sizeControl = cabins[refc].RandomAdjust(hulls, FRONT, sizeControl, 0.8);
		}

		cabinLength = Random(1.0, 1.5) * cabinLength * sizeControl;

		cabins[refc].Extrude(hulls, cabins[cabinNum], FRONT, cabinLength * sizeControl, sizeControl, false);
		refc = cabinNum;
		cabinNum++;
		cabins[refc].type = ShipCabinImpl::FRONT_END;

		cabins[refc].Extrude(hulls, cabins[cabinNum], FRONT, cabinLength * sizeControl, sizeControl, false);
		refc = cabinNum;
		cabinNum++;
		float scaleFront = Random(0.3, 0.5);
		ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, hulls.faces[cabins[refc].fi[FRONT]],
			Vector(1, scaleFront, scaleFront));
		cabins[refc].type = ShipCabinImpl::FRONT_END;

		// Add some large asynmmetrical sections of the hull that stick out
		int mainCabinNum = cabinNum;
		if (create_asymmetry_segments) {
			for (int i = 0; i < mainCabinNum; i++) {
				if (cabins[i].type == ShipCabinImpl::MAIN) {
					Direction dirs[4] = {
						LEFT,
						RIGHT,
						TOP,
						BOTTOM
					};
					for (auto& dir : dirs) {
						if (cabins[i].fType[dir] == ShipCabinImpl::UNSPECIFIED) {
							Face sideFace = hulls.faces[cabins[i].fi[dir]];
							if (GetAspectRatio(sideFace, *hulls.verts) < 3.0 && Random(0.0, 1.0) < 0.35) {
								if (cabins[i].predecessor != NULL)
									if (cabins[i].predecessor->fType[dir] == ShipCabinImpl::INVALID)
										continue;

								// Extrude!
								cabins[i].fType[dir] = ShipCabinImpl::INVALID;
								int segments = (int)Random(num_asymmetry_segments_min, num_asymmetry_segments_max + 1);

								//if (apply_longer_asymmetry_segments > 0)
								//	if (Random(0.0, 1.0) > 0.75) {
								//		segments += 2;
								//		apply_longer_asymmetry_segments -= 1;
								//	}

								for (int j = 0; j < segments; j++) {
									int cabinLength = Random(0.1 * shipCoreSize, 0.4 * shipCoreSize);
									cabins[i].Extrude(hulls, cabins[cabinNum], dir, cabinLength, 1.0, false);
									if (j == 0)
										cabins[cabinNum].type = ShipCabinImpl::SIDE_BASE;
									else if (j == segments - 1)
										cabins[cabinNum].type = ShipCabinImpl::SIDE_END;
									else
										cabins[cabinNum].type = ShipCabinImpl::SIDE;

									// maybe apply some scaling
									if (j == 0 || j == segments - 1 || Random(0.0, 1.0) > 0.25) {
										sideFace = hulls.faces[cabins[cabinNum].fi[dir]];
										float s = 1.0 / Random(1.5, 1.8);
										float sx = 1.0 / Random(1.1, 1.3);
										ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, sideFace, Vector(sx, s, s));
									}
									cabinNum++;
								}
							}
						}
					}
				}
			}
		}

		hyvector<Face> main_engine_faces;
		hyvector<Face> engine_faces;
		hyvector<Face> grid_faces;
		hyvector<Face> antenna_faces;
		hyvector<Face> weapon_faces;
		hyvector<Face> sphere_faces;
		hyvector<Face> disc_faces;
		hyvector<Face> cylinder_faces;



		for (int i = 0; i < cabinNum; i++) {
			if (cabins[i].type == ShipCabinImpl::BACK_END) {
				main_engine_faces.push_back(hulls.faces[cabins[i].fi[BACK]]);
				cabins[i].fType[BACK] = ShipCabinImpl::ENGINE;
			}
			else if (cabins[i].type == ShipCabinImpl::MAIN) {
				if (cabins[i].fType[LEFT] == ShipCabinImpl::UNSPECIFIED) {
					if (Random(0.0, 1.0) > 0.5) {
						if (cabins[i].predecessor != NULL)
							if (cabins[i].predecessor->fType[LEFT] == ShipCabinImpl::SPHERE)
								continue;
						sphere_faces.push_back(hulls.faces[cabins[i].fi[LEFT]]);
						cabins[i].fType[LEFT] = ShipCabinImpl::SPHERE;
					}
				}
				if (cabins[i].fType[RIGHT] == ShipCabinImpl::UNSPECIFIED) {
					if (Random(0.0, 1.0) > 0.5) {
						if (cabins[i].predecessor != NULL)
							if (cabins[i].predecessor->fType[RIGHT] == ShipCabinImpl::SPHERE)
								continue;
						sphere_faces.push_back(hulls.faces[cabins[i].fi[RIGHT]]);
						cabins[i].fType[RIGHT] = ShipCabinImpl::SPHERE;
					}
				}
				if (cabins[i].fType[BOTTOM] == ShipCabinImpl::UNSPECIFIED)
					if (Random(0.0, 1.0) < 0.5 &&
						hulls.faces[cabins[i].fi[BOTTOM]].Normal(*hulls.verts).z < -0.8) {
						disc_faces.push_back(hulls.faces[cabins[i].fi[BOTTOM]]);
						cabins[i].fType[BOTTOM] = ShipCabinImpl::DISC;
					}
			}
			else if (cabins[i].type == ShipCabinImpl::SIDE_BASE) {
				if (cabins[i].fType[BACK] == ShipCabinImpl::UNSPECIFIED)
					if (Random(0.0, 1.0) < 0.5) {
						engine_faces.push_back(hulls.faces[cabins[i].fi[BACK]]);
						cabins[i].fType[BACK] = ShipCabinImpl::ENGINE;
					}
			}
		}

		// Now we've categorized, let's actually add the detail

		for (auto& face : main_engine_faces)
			AddMainEngineToFace(face);

		for (auto& face : engine_faces)
			AddEngineToFace(face);

		//for (auto& face : grid_faces)
		//	AddGridToFace(face);

		for (auto& face : antenna_faces)
			AddAntennaToFace(face);

		for (auto& face : weapon_faces)
			AddWeaponToFace(face);

		for (auto& face : sphere_faces)
			AddSphereToFace(face);

		for (auto& face : disc_faces)
			AddDiscToFace(face);

		for (auto& face : cylinder_faces)
			AddCylinderToFace(face);


		for (int i = 0; i < cabinNum; i++)
			cabins[i].Draw(internals, *hulls.verts, 0.9);
		ShapeHelper::Bevel(*internals.verts, internals.faces, 0.5);
		ShapeHelper::SeperateVert(*internals.verts, internals.faces, *internals.normals);
		ShapeHelper::ToTris(*internals.verts, internals.faces, *internals.tris);

		ShapeHelper::Bevel(*hulls.verts, hulls.faces, 0.5);
		ShapeHelper::SeperateVert(*hulls.verts, hulls.faces, *hulls.normals);
		ShapeHelper::AssignGlobalUV(*hulls.verts, hulls.faces, *hulls.uvs, shipCoreSize);
		ShapeHelper::ToTris(*hulls.verts, hulls.faces, *hulls.tris);

		ShapeHelper::SeperateVert(*accessories.verts, accessories.faces, *accessories.normals);
		ShapeHelper::AssignGlobalUV(*accessories.verts, accessories.faces, *accessories.uvs, shipCoreSize);
		ShapeHelper::ToTris(*accessories.verts, accessories.faces, *accessories.tris);

		ShapeHelper::SeperateVert(*spheres.verts, spheres.faces, *spheres.normals);
		ShapeHelper::AssignGlobalUV(*spheres.verts, spheres.faces, *spheres.uvs, shipCoreSize);
		ShapeHelper::ToTris(*spheres.verts, spheres.faces, *spheres.tris);

		ShapeHelper::Bevel(*engines.verts, engines.faces, 0.5);
		ShapeHelper::SeperateVert(*engines.verts, engines.faces, *engines.normals);
		ShapeHelper::ToTris(*engines.verts, engines.faces, *engines.tris);

		ShapeHelper::ToTris(*lights.verts, lights.faces, *lights.tris);

		//Add window faces along the hull
		ShapeHelper::Bevel((*components[WINDOW_EMIT].verts), components[WINDOW_EMIT].faces, 0.5);
		for (auto face : components[HULL_SPHERE].faces) {
			Vector normal = face.Normal((*components[HULL_SPHERE].verts));
			if (abs(normal.z) < 0.7) {
				ShapeHelper::CopyFace((*components[WINDOW_EMIT].verts), components[WINDOW_EMIT].faces, (*components[WINDOW_EMIT].uvs), face, (*components[HULL_SPHERE].verts), hyvector<Vector2f>(), 0.0);
			}
		}
		ShapeHelper::SeperateVert((*components[WINDOW_EMIT].verts), components[WINDOW_EMIT].faces, (*components[WINDOW_EMIT].normals));
		for (auto face : components[WINDOW_EMIT].faces) {
			Vector normal = face.Normal((*components[WINDOW_EMIT].verts));
			for (auto vi : face.vi)
				(*components[WINDOW_EMIT].verts)[vi] += normal*0.01;
		}
		ShapeHelper::AssignGlobalUV((*components[WINDOW_EMIT].verts), components[WINDOW_EMIT].faces, (*components[WINDOW_EMIT].uvs), shipCoreSize * 2.5);
		ShapeHelper::ToTris((*components[WINDOW_EMIT].verts), components[WINDOW_EMIT].faces, *components[WINDOW_EMIT].tris);

		return true;

	}

	SpaceshipGenerator::~SpaceshipGenerator()
	{
	}

	// Extrudes a face along its normal by translate_forwards units.
	// Returns the new face, and optionally fills out extruded_face_list
	// with all the additional side faces created from the extrusion.


	void SpaceshipGenerator::AddMainEngineToFace(const Face& face) {
		// Given a face, splits it into a uniform grid and extrudes each grid face
		// out and back in again, making an exhaust shape.
		// The more square the face is, the more grid divisions it might have
		if (!face.valid)
			return;
		ShipComponentImpl& hulls = components[HULL];
		ShipComponentImpl& engines = components[ENGINE_EMIT];
		int num_cuts_max = int(6 - GetAspectRatio(face, *hulls.verts));
		int num_cuts = (int)Random(1.5, num_cuts_max);
		float exhaust_main_length = Random(0.5, 1.0) * shipCoreSize;
		float exhaust_length = Random(0.2, 0.3) * shipCoreSize;
		float scale_outer = 1 / Random(1.3, 1.6);
		float scale_inner = Random(1.1, 1.2);
		hyvector<Face> faces_divided = ShapeHelper::SubDivideFace(*hulls.verts, hulls.faces, face, num_cuts);
		for (auto face_divided : faces_divided) {
			Face new_face = face_divided;
			new_face = ShapeHelper::ExtrudeFace(*hulls.verts, hulls.faces, new_face, exhaust_main_length / num_cuts);
			ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, new_face, Vector(1.2, 1.2, 1.2));
			new_face = ShapeHelper::ExtrudeFace(*hulls.verts, hulls.faces, new_face, exhaust_length / num_cuts);
			ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, new_face, Vector(scale_outer, scale_outer, scale_outer));
			new_face = ShapeHelper::ExtrudeFace(*hulls.verts, hulls.faces, new_face, -exhaust_main_length / num_cuts);
			ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, new_face, Vector(scale_inner, scale_inner, scale_inner));

			new_face = ShapeHelper::CopyFace(*engines.verts, engines.faces, *engines.uvs, new_face, *hulls.verts, *hulls.uvs, 0.0);
			ShapeHelper::ScaleFace(*engines.verts, engines.faces, new_face, Vector(0.5, 0.5, 0.5));
			ShapeHelper::AddReversedFace(engines.faces, new_face);
			new_face = ShapeHelper::ExtrudeFace(*engines.verts, engines.faces, new_face, exhaust_length / num_cuts);
		}
	}

	void SpaceshipGenerator::AddEngineToFace(const Face& face) {
		// Given a face, splits it into a uniform grid and extrudes each grid face
		// out and back in again, making an exhaust shape.
		// The more square the face is, the more grid divisions it might have
		if (!face.valid)
			return;
		ShipComponentImpl& hulls = components[HULL];
		ShipComponentImpl& engines = components[ENGINE_EMIT];
		float exhaust_length = Random(0.1, 0.15) * shipCoreSize;
		float scale_outer = 1 / Random(1.5, 2.0);
		float scale_inner = Random(1.1, 1.2);
		Face new_face = face;
		new_face = ShapeHelper::ExtrudeFace(*hulls.verts, hulls.faces, new_face, exhaust_length);
		ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, new_face, Vector(scale_outer, scale_outer, scale_outer));
		new_face = ShapeHelper::ExtrudeFace(*hulls.verts, hulls.faces, new_face, -exhaust_length);
		ShapeHelper::ScaleFace(*hulls.verts, hulls.faces, new_face, Vector(scale_inner, scale_inner, scale_inner));

		//ShapeHelper::CopyFace(engines.verts, engines.faces, engines.uvs, new_face, hulls.verts, hulls.uvs, 0.0);
		new_face = ShapeHelper::CopyFace(*engines.verts, engines.faces, *engines.uvs, new_face, *hulls.verts, *hulls.uvs, -0.25 * exhaust_length);
		ShapeHelper::ScaleFace(*engines.verts, engines.faces, new_face, Vector(0.5, 0.5, 0.5));
		ShapeHelper::AddReversedFace(engines.faces, new_face);
		new_face = ShapeHelper::ExtrudeFace(*engines.verts, engines.faces, new_face, 0.5 * exhaust_length);
	}

	void SpaceshipGenerator::AddGridToFace(const Face& face) {
		if (!face.valid)
			return;
		int num_cuts = (int)Random(2, 5);
		float grid_length = Random(0.025 * shipCoreSize, 0.15 * shipCoreSize);
		float scale = 0.8;
		hyvector<Face> faces_divided = ShapeHelper::SubDivideFace((*components[HULL].verts), components[HULL].faces, face, num_cuts);
		for (auto face_divided : faces_divided) {
			//material_index = Material.hull_lights if random() > 0.5 else Material.hull
			Face new_face = ShapeHelper::ExtrudeFace((*components[HULL].verts), components[HULL].faces, face_divided, grid_length);
			//for extruded_face in extruded_face_list :
			//if abs(face.normal.z) < 0.707 : # side face
			//extruded_face.material_index = material_index
			ShapeHelper::ScaleFace((*components[HULL].verts), components[HULL].faces, new_face, Vector(scale, scale, scale));
		}
	}

	void SpaceshipGenerator::AddAntennaToFace(const Face& face) {
		if (!face.valid || face.vi.size() < 4)
			return;
		int horizontal_step = (int)Random(6, 12);
		int vertical_step = (int)Random(6, 12);
		for (int h = 2; h < horizontal_step - 1; h++) {
			for (int v = 2; v < vertical_step - 1; v++) {
				if (Random(0.0, 1.0) > 0.9) {
					Vector v1 = (*components[HULL].verts)[face.vi[0]] + ((*components[HULL].verts)[face.vi[3]] - (*components[HULL].verts)[face.vi[0]]) * h / horizontal_step;
					Vector v2 = (*components[HULL].verts)[face.vi[1]] + ((*components[HULL].verts)[face.vi[2]] - (*components[HULL].verts)[face.vi[1]]) * h / horizontal_step;
					Vector pos = v1 + (v2 - v1) * v / vertical_step;
					float height = Random(0.8, 3.0)*face.MinEdge((*components[HULL].verts));
					bool add_light = false;
					if (height > 2.0*face.MinEdge((*components[HULL].verts)))
						add_light = true;
					if (height < 0.2 * shipCoreSize)
						continue;
					float base_height = Random(0.05, 0.15)*height;
					float base_radius = Random(0.4, 4);
					int num_segments = (int)Random(4.0, 7.0);
					// Spire
					ShapeHelper::AddAltar((*components[ACCESSORY].verts), components[ACCESSORY].faces, pos, face.Normal((*components[HULL].verts)), base_radius, 0.1*base_radius, height, num_segments);
					// Base
					ShapeHelper::AddAltar((*components[ACCESSORY].verts), components[ACCESSORY].faces, pos, face.Normal((*components[HULL].verts)), base_radius*Random(1.5, 2.0), base_radius*Random(1.0, 1.5), base_height, num_segments);
					// Light
					if (add_light)
						ShapeHelper::AddCylinder((*components[LIGHT].verts), components[LIGHT].faces, pos + face.Normal((*components[HULL].verts))*height*0.99, face.Normal((*components[HULL].verts)), 0.15*base_radius, 0.5*base_radius, num_segments);
				}
			}
		}
	}

	void SpaceshipGenerator::AddWeaponToFace(const Face& face) {
		if (!face.valid)
			//TO DO: ZY
			return;
	}

	void SpaceshipGenerator::AddSphereToFace(const Face& face) {
		if (!face.valid)
			return;
		ShipComponentImpl& hulls = components[HULL];
		ShipComponentImpl& accessories = components[ACCESSORY];
		float radius = Random(0.5, 1.2) * face.MinEdge(*hulls.verts);
		int num_segments = 6;
		Vector spCenter = face.Center(*hulls.verts) - 0.15 * radius * face.Normal(*hulls.verts);
		ShapeHelper::AddSphere(*accessories.verts, accessories.faces, spCenter, radius, num_segments);
	}

	void SpaceshipGenerator::AddDiscToFace(const Face& face) {
		if (!face.valid)
			return;
		ShipComponentImpl& hulls = components[HULL];
		ShipComponentImpl& accessories = components[ACCESSORY];
		ShipComponentImpl& engines = components[ENGINE_EMIT];
		float height = 0.25 * face.MinEdge(*hulls.verts);
		ShapeHelper::AddAltar(*accessories.verts, accessories.faces, face.Center(*hulls.verts) - height * face.Normal(*hulls.verts),
			face.Normal(*hulls.verts), height * 4, height * 3, height * 2, 8);
		Face new_face = ShapeHelper::CopyFace(*engines.verts, engines.faces, *engines.uvs, accessories.faces[accessories.faces.size() - 1], *accessories.verts, hyvector<Vector2f>(), -0.5 * height);
		ShapeHelper::AddReversedFace(engines.faces, new_face);
		new_face = ShapeHelper::ExtrudeFace(*engines.verts, engines.faces, new_face, 0.6 * height);
		new_face = ShapeHelper::ExtrudeFace(*engines.verts, engines.faces, new_face, 0.0);
		ShapeHelper::ScaleFace(*engines.verts, engines.faces, new_face, Vector(0.5, 0.5, 0.5));
		new_face = ShapeHelper::ExtrudeFace(*engines.verts, engines.faces, new_face, -0.6 * height);
	}

	void SpaceshipGenerator::AddCylinderToFace(const Face& face) {
		if (!face.valid || face.vi.size() < 4 || GetAspectRatio(face, (*components[HULL].verts)) > 2)
			return;
		int step = (int)Random(5, 9);
		int num_segments = (int)Random(6, 13);
		float cylinder_depth = 2 * face.MinEdge((*components[HULL].verts)) / (step + 2);
		float cylinder_size = cylinder_depth * 0.6;
		// UE_LOG(LogTemp, Log, TEXT("cylinder:%f,%f,%d"), cylinder_size, cylinder_depth, num_segments);
		for (int h = 2; h < step - 1; h++) {
			for (int v = 2; v < step - 1; v++) {
				Vector v1 = (*components[HULL].verts)[face.vi[0]] + ((*components[HULL].verts)[face.vi[3]] - (*components[HULL].verts)[face.vi[0]]) * h / step;
				Vector v2 = (*components[HULL].verts)[face.vi[1]] + ((*components[HULL].verts)[face.vi[2]] - (*components[HULL].verts)[face.vi[1]]) * h / step;
				Vector pos = v1 + (v2 - v1) * v / step;
				Vector cylinder_dir = face.Normal((*components[HULL].verts)).RotateAngleAxis(90, Vector(0.0, 1.0, 0.0));
				ShapeHelper::AddCylinder((*components[ACCESSORY].verts), components[ACCESSORY].faces, pos, cylinder_dir, cylinder_size, cylinder_depth, num_segments);
				if (Random(0.0, 1.0) > 0.75)
					ShapeHelper::AddCylinder((*components[LIGHT].verts), components[LIGHT].faces, pos + face.Normal((*components[HULL].verts))*cylinder_size, cylinder_dir, cylinder_size*0.1, cylinder_depth*0.2, num_segments);
			}
		}
	}

	float SpaceshipGenerator::Random(float start, float end)
	{
		return start + (end - start) * rand() / float(RAND_MAX);
	}

	void HySpaceshipImpl::GenerateSpaceship(Vector center, Vector rotation, float spaceship_size, unsigned int seed, Spaceship& outship)
	{
		SpaceshipGenerator generator;
		generator.Generate(center, rotation, spaceship_size, seed);

		for (int i = 0; i < Spaceship::COMPONENT_COUNT; ++i)
		{
			ShipComponent& comp = outship.components[i];
			ShipComponentImpl& compi = generator.components[i];
			compi.verts->SetOutVectorValues(); comp.verts = compi.verts;
			compi.normals->SetOutVectorValues(); comp.normals = compi.normals;
			compi.tris->SetOutVectorValues(); comp.tris = compi.tris;
			compi.uvs->SetOutVectorValues(); comp.uvs = compi.uvs;
		}
	}

}
