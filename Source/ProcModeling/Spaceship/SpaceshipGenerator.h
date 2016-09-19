// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "HyVoxelLib.h"
#include "Common/ShapeHelper.h"

namespace HyVoxel {

	//struct ShipComponent {
	//	hyvector<Face> faces;
	//	hyvector<Vector> verts;
	//	hyvector<int> tris;
	//	hyvector<Vector> normals;
	//	hyvector<Vector2f> uvs;
	//};

	struct ShipComponentImpl
	{
		hyvector<Face>				faces;
		OutVectorImpl<Vector>*		verts = NULL;
		OutVectorImpl<Vector>*		normals = NULL;
		OutVectorImpl<int>*			tris = NULL;
		OutVectorImpl<Vector2f>*	uvs = NULL;
	};

	struct ShipCabinImpl {

		enum CType {
			MAIN = 0,			// Possible antenna, sphere, disc, cylinder added
			FRONT_END,		// Possible antenna added
			BACK_END,		// Main engine added
			RIBBED,			// Nothing added
			SIDE_BASE,		// Possible engine added; no sphere
			SIDE,			// Possible engine added; no sphere
			SIDE_END,		// Possible engine added; no sphere, antenna
							//...
		}
		type = MAIN;

		enum FType {
			INVALID = 0,
			UNSPECIFIED,
			SIDED,
			ENGINE,
			GRID,
			ANTENNA,
			WEAPON,
			SPHERE,
			DISC,
			CYLINDER
		}
		fType[6] = { INVALID };;

		ShipCabinImpl*				successors[6] = { NULL };
		ShipCabinImpl*				predecessor = NULL;
		int							fi[6] = { -1 };
		int							vi[8] = { -1 };
		float						length = -1;
		float						sizeControl = 1.0;

		void InitCore(ShipComponentImpl& component, Vector shipCenter, float shipCoreSize);

		void Extrude(ShipComponentImpl& component, ShipCabinImpl& newCabin, Direction dir,
			float cabinLength, float sizeControl = 1.0, bool useRibbing = true);

		float RandomAdjust(ShipComponentImpl& component, Direction dir, float sizeControl = 1.0, float sizeExpect = 1.0,
			bool useScaling = true, bool useTranslation = true, bool useRotation = true);

		void Draw(ShipComponentImpl& component, const hyvector<Vector>& refVerts,
			float scale = 1.0, float margin = 0.0);

		Direction Opposite(Direction dir)
		{
			return Direction((dir + 3) % 6);
		}
	};

	class SpaceshipGenerator
	{
	public:

		static const int HULL = 0;
		static const int CABIN = 1;
		static const int ACCESSORY = 2;
		static const int ENGINE_EMIT = 3;
		static const int WINDOW_EMIT = 4;
		static const int HULL_SPHERE = 5;
		static const int LIGHT = 6;

		int num_hull_segments_min;
		int num_hull_segments_max;
		bool create_asymmetry_segments;
		int num_asymmetry_segments_min;
		int num_asymmetry_segments_max;
		bool create_face_detail;
		bool apply_bevel_modifier;
		//bool allow_horizontal_symmetry;
		//bool allow_vertical_symmetry;

		SpaceshipGenerator(
			int num_hull_segments_min = 4,
			int num_hull_segments_max = 6,
			bool create_asymmetry_segments = true,
			int num_asymmetry_segments_min = 1,
			int num_asymmetry_segments_max = 3,
			bool create_face_detail = true,
			bool apply_bevel_modifier = true);
		~SpaceshipGenerator();

		ShipComponentImpl components[Spaceship::COMPONENT_COUNT];
		ShipCabinImpl cabins[Spaceship::CABIN_COUNT];
		int cabinNum = 0;

		float shipCoreSize;
		Vector shipCenter;

		float Random(float start, float end);
		//float GetAspectRatio(const Face& face, std::vector<Vector> verts);
		//Face RibbedExtrudeFace(const Face& face, float translate_forwards, int num_ribs, float rib_scale);

		void AddMainEngineToFace(const Face& face);
		void AddEngineToFace(const Face& face);
		void AddGridToFace(const Face& face);
		void AddAntennaToFace(const Face& face);
		void AddWeaponToFace(const Face& face);
		void AddSphereToFace(const Face& face);
		void AddDiscToFace(const Face& face);
		void AddCylinderToFace(const Face& face);

		bool Generate(Vector center, Vector rotation, float size, unsigned int seed = 0);
	};

	class HySpaceshipImpl : public IHySpaceshipInterface
	{
	public:
		virtual void Release() override { delete this; }

		virtual void GenerateSpaceship(Vector center, Vector rotation, float spaceship_size, unsigned int seed, Spaceship& outship) override;

		// TODO: liang
		//  to be moved into somewhere like help function interface.
		virtual void GenerateSphere(const Vector& center, float radius, int num, OutputMesh& omesh) override;
        virtual void GenerateTriangle(const Vector& v0, const Vector& v1, const Vector& v2, OutputMesh& omesh) override;
	};
}
