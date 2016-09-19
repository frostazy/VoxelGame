/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

#include "HyVoxelConfig.h"

// we support per-material physics properties.  currently, the following 10 basic materials are defined.
namespace HyVoxel
{
	namespace Physics
	{
		// base material classes used for realistic physics simulation
		enum materialClasses
		{
			MAT_CORK = 0,
			MAT_WOOD,
			MAT_ICE,
			MAT_WATER,
			MAT_PLASTIC,
			MAT_CONCRETE,
			MAT_ROCK,
			MAT_METAL,
			MAT_HELIUM,
			MAT_DEFAULT
		};
		static const unsigned char cMATCLASSES = 10;

		// common physics properties
		const float CONST_GRAVITY = -100.f;		// ~1x gravity [dm^2/s]
		const float AIR_DENSITY = 0.000001205f;	// air density [Mg/dm^3]
		const float TERM_VEL = 550.f;			// ~terminal velocity [dm/s]

		// physics material library (one entry per materialClass)
		class CPhysicsMaterialLibrary
		{
		public:
			CPhysicsMaterialLibrary() { initMaterials(); };
		private:
			struct matData
			{
				double density;				// approx. material density [Mg/dm^3]
				float* friction = NULL;		// coeff of friction [dimensionless]
				float* rfriction = NULL;	// coeff of rolling friction [dimensionless]
				float restitution;			// coeff of restitution (i.e. elasticity) [dimensionless]
				double cstrength;			// compressive strength (breaking pressure -- fracture) [Mg/(dm*s^2)]
				matData()
				{
					friction = (float*)VF_RAWALLOC(cMATCLASSES*sizeof(float));
					rfriction = (float*)VF_RAWALLOC(cMATCLASSES*sizeof(float));
				}
				~matData()
				{
					VF_FREE(friction);
					VF_FREE(rfriction);
				}
				matData& operator =(const matData& a)
				{
					density = a.density;
					restitution = a.restitution;
					cstrength = a.cstrength;
					memcpy(friction, a.friction, cMATCLASSES*sizeof(float));
					memcpy(rfriction, a.rfriction, cMATCLASSES*sizeof(float));
					return *this;
				}
			};
			void initMaterials();
		public:
			TVFMap<int, matData> m_materials;
		};
	};
};
