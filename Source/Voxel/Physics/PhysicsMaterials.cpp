#include "HyVoxelPrivatePCH.h"

#include "PhysicsMaterials.h"

void HyVoxel::Physics::CPhysicsMaterialLibrary::initMaterials()
{
	// Constructs the physics material library.
	// TODO :: this will eventually be read from a (plain text?) file

	// The Coefficient of Friction is defined between contact pairs.
	// We must define [MATCOUNT] coefficients per material.

	// CORK ----------------------------------------------------------
	m_materials[MAT_CORK].density = 0.00025f;
	m_materials[MAT_CORK].restitution = 0.6f;
	m_materials[MAT_CORK].restitution = 0.6f;
	m_materials[MAT_CORK].cstrength = 2.f;

	m_materials[MAT_CORK].friction[MAT_CORK] = 0.7f;
	m_materials[MAT_CORK].friction[MAT_WOOD] = m_materials[MAT_WOOD].friction[MAT_CORK] = 0.7f;
	m_materials[MAT_CORK].friction[MAT_ICE] = m_materials[MAT_ICE].friction[MAT_CORK] = 0.1f;
	m_materials[MAT_CORK].friction[MAT_WATER] = m_materials[MAT_WATER].friction[MAT_CORK] = 1.0f;
	m_materials[MAT_CORK].friction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].friction[MAT_CORK] = 0.7f;
	m_materials[MAT_CORK].friction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].friction[MAT_CORK] = 0.7f;
	m_materials[MAT_CORK].friction[MAT_ROCK] = m_materials[MAT_ROCK].friction[MAT_CORK] = 0.7f;
	m_materials[MAT_CORK].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_CORK] = 0.7f;
	m_materials[MAT_CORK].friction[MAT_DEFAULT] = m_materials[MAT_CORK].friction[MAT_ROCK];

	m_materials[MAT_CORK].rfriction[MAT_CORK] = 0.1f;
	m_materials[MAT_CORK].rfriction[MAT_WOOD] = m_materials[MAT_WOOD].rfriction[MAT_CORK] = 0.1f;
	m_materials[MAT_CORK].rfriction[MAT_ICE] = m_materials[MAT_ICE].rfriction[MAT_CORK] = 0.0001f;
	m_materials[MAT_CORK].rfriction[MAT_WATER] = m_materials[MAT_WATER].rfriction[MAT_CORK] = 1.0f;
	m_materials[MAT_CORK].rfriction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].rfriction[MAT_CORK] = 0.1f;
	m_materials[MAT_CORK].rfriction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].rfriction[MAT_CORK] = 0.1f;
	m_materials[MAT_CORK].rfriction[MAT_ROCK] = m_materials[MAT_ROCK].rfriction[MAT_CORK] = 0.1f;
	m_materials[MAT_CORK].rfriction[MAT_METAL] = m_materials[MAT_METAL].rfriction[MAT_CORK] = 0.1f;
	m_materials[MAT_CORK].rfriction[MAT_DEFAULT] = m_materials[MAT_CORK].rfriction[MAT_ROCK];

	// WOOD ----------------------------------------------------------
	m_materials[MAT_WOOD].density = 0.0007f;
	m_materials[MAT_WOOD].restitution = 0.6f;
	m_materials[MAT_WOOD].cstrength = 1035.f;

	m_materials[MAT_WOOD].friction[MAT_WOOD] = 0.3f;
	m_materials[MAT_WOOD].friction[MAT_ICE] = m_materials[MAT_ICE].friction[MAT_WOOD] = 0.1f;
	m_materials[MAT_WOOD].friction[MAT_WATER] = m_materials[MAT_WATER].friction[MAT_WOOD] = 1.0f;
	m_materials[MAT_WOOD].friction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].friction[MAT_WOOD] = 0.3f;
	m_materials[MAT_WOOD].friction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].friction[MAT_WOOD] = 0.6f;
	m_materials[MAT_WOOD].friction[MAT_ROCK] = m_materials[MAT_ROCK].friction[MAT_WOOD] = 0.4f;
	m_materials[MAT_WOOD].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_WOOD] = 0.4f;
	m_materials[MAT_WOOD].friction[MAT_DEFAULT] = m_materials[MAT_WOOD].friction[MAT_ROCK];

	m_materials[MAT_WOOD].rfriction[MAT_WOOD] = 0.04f;
	m_materials[MAT_WOOD].rfriction[MAT_ICE] = m_materials[MAT_ICE].rfriction[MAT_WOOD] = 0.0001f;
	m_materials[MAT_WOOD].rfriction[MAT_WATER] = m_materials[MAT_WATER].rfriction[MAT_WOOD] = 1.0f;
	m_materials[MAT_WOOD].rfriction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].rfriction[MAT_WOOD] = 0.04f;
	m_materials[MAT_WOOD].rfriction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].rfriction[MAT_WOOD] = 0.05f;
	m_materials[MAT_WOOD].rfriction[MAT_ROCK] = m_materials[MAT_ROCK].rfriction[MAT_WOOD] = 0.04f;
	m_materials[MAT_WOOD].rfriction[MAT_METAL] = m_materials[MAT_METAL].rfriction[MAT_WOOD] = 0.04f;
	m_materials[MAT_WOOD].rfriction[MAT_DEFAULT] = m_materials[MAT_WOOD].rfriction[MAT_ROCK];

	// ICE ----------------------------------------------------------
	m_materials[MAT_ICE].density = 0.0009f;
	m_materials[MAT_ICE].restitution = 0.65f;
	m_materials[MAT_ICE].cstrength = 750.f;

	m_materials[MAT_ICE].friction[MAT_ICE] = 0.01f;
	m_materials[MAT_ICE].friction[MAT_WATER] = m_materials[MAT_WATER].friction[MAT_ICE] = 1.0f;
	m_materials[MAT_ICE].friction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].friction[MAT_ICE] = 0.1f;
	m_materials[MAT_ICE].friction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].friction[MAT_ICE] = 0.1f;
	m_materials[MAT_ICE].friction[MAT_ROCK] = m_materials[MAT_ROCK].friction[MAT_ICE] = 0.1f;
	m_materials[MAT_ICE].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_ICE] = 0.1f;
	m_materials[MAT_ICE].friction[MAT_DEFAULT] = m_materials[MAT_ICE].friction[MAT_ROCK];

	m_materials[MAT_ICE].rfriction[MAT_ICE] = 0.00001f;
	m_materials[MAT_ICE].rfriction[MAT_WATER] = m_materials[MAT_WATER].rfriction[MAT_ICE] = 1.0f;
	m_materials[MAT_ICE].rfriction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].rfriction[MAT_ICE] = 0.0001f;
	m_materials[MAT_ICE].rfriction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].rfriction[MAT_ICE] = 0.0001f;
	m_materials[MAT_ICE].rfriction[MAT_ROCK] = m_materials[MAT_ROCK].rfriction[MAT_ICE] = 0.0001f;
	m_materials[MAT_ICE].rfriction[MAT_METAL] = m_materials[MAT_METAL].rfriction[MAT_ICE] = 0.0001f;
	m_materials[MAT_ICE].rfriction[MAT_DEFAULT] = m_materials[MAT_ICE].rfriction[MAT_ROCK];

	// WATER ----------------------------------------------------------
	m_materials[MAT_WATER].density = 0.001f;
	m_materials[MAT_WATER].restitution = 0.01f;
	m_materials[MAT_WATER].cstrength = 0.f;

	m_materials[MAT_WATER].friction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].friction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].friction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].friction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].friction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].friction[MAT_ROCK] = m_materials[MAT_ROCK].friction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].friction[MAT_DEFAULT] = m_materials[MAT_WATER].friction[MAT_ROCK];

	m_materials[MAT_WATER].rfriction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].rfriction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].rfriction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].rfriction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].rfriction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].rfriction[MAT_ROCK] = m_materials[MAT_ROCK].rfriction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].rfriction[MAT_METAL] = m_materials[MAT_METAL].rfriction[MAT_WATER] = 1.0f;
	m_materials[MAT_WATER].rfriction[MAT_DEFAULT] = m_materials[MAT_WATER].rfriction[MAT_ROCK];

	// PLASTIC ----------------------------------------------------------
	m_materials[MAT_PLASTIC].density = 0.0012f;
	m_materials[MAT_PLASTIC].restitution = 0.7f;
	m_materials[MAT_PLASTIC].cstrength = 10000.f;

	m_materials[MAT_PLASTIC].friction[MAT_PLASTIC] = 0.6f;
	m_materials[MAT_PLASTIC].friction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].friction[MAT_PLASTIC] = 0.5f;
	m_materials[MAT_PLASTIC].friction[MAT_ROCK] = m_materials[MAT_ROCK].friction[MAT_PLASTIC] = 0.4f;
	m_materials[MAT_PLASTIC].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_PLASTIC] = 0.3f;
	m_materials[MAT_PLASTIC].friction[MAT_DEFAULT] = m_materials[MAT_PLASTIC].friction[MAT_ROCK];

	m_materials[MAT_PLASTIC].rfriction[MAT_PLASTIC] = 0.04f;
	m_materials[MAT_PLASTIC].rfriction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].rfriction[MAT_PLASTIC] = 0.07f;
	m_materials[MAT_PLASTIC].rfriction[MAT_ROCK] = m_materials[MAT_ROCK].rfriction[MAT_PLASTIC] = 0.06f;
	m_materials[MAT_PLASTIC].rfriction[MAT_METAL] = m_materials[MAT_METAL].rfriction[MAT_PLASTIC] = 0.06f;
	m_materials[MAT_PLASTIC].rfriction[MAT_DEFAULT] = m_materials[MAT_PLASTIC].rfriction[MAT_ROCK];

	// CONCRETE ----------------------------------------------------------
	m_materials[MAT_CONCRETE].density = 0.002f;
	m_materials[MAT_CONCRETE].restitution = 0.65f;
	m_materials[MAT_CONCRETE].cstrength = 3000.f;

	m_materials[MAT_CONCRETE].friction[MAT_CONCRETE] = 0.7f;
	m_materials[MAT_CONCRETE].friction[MAT_ROCK] = m_materials[MAT_ROCK].friction[MAT_CONCRETE] = 0.7f;
	m_materials[MAT_CONCRETE].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_CONCRETE] = 0.6f;
	m_materials[MAT_CONCRETE].friction[MAT_DEFAULT] = m_materials[MAT_CONCRETE].friction[MAT_ROCK];

	m_materials[MAT_CONCRETE].rfriction[MAT_CONCRETE] = 0.07f;
	m_materials[MAT_CONCRETE].rfriction[MAT_ROCK] = m_materials[MAT_ROCK].rfriction[MAT_CONCRETE] = 0.07f;
	m_materials[MAT_CONCRETE].rfriction[MAT_METAL] = m_materials[MAT_METAL].rfriction[MAT_CONCRETE] = 0.06f;
	m_materials[MAT_CONCRETE].rfriction[MAT_DEFAULT] = m_materials[MAT_CONCRETE].rfriction[MAT_ROCK];

	// ROCK ----------------------------------------------------------
	m_materials[MAT_ROCK].density = 0.004f;
	m_materials[MAT_ROCK].restitution = 0.725f;
	m_materials[MAT_ROCK].cstrength = 13000.f;

	m_materials[MAT_ROCK].friction[MAT_ROCK] = 0.6f;
	m_materials[MAT_ROCK].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_ROCK] = 0.6f;
	m_materials[MAT_ROCK].friction[MAT_DEFAULT] = m_materials[MAT_ROCK].friction[MAT_ROCK];

	m_materials[MAT_ROCK].rfriction[MAT_ROCK] = 0.06f;
	m_materials[MAT_ROCK].rfriction[MAT_METAL] = m_materials[MAT_METAL].rfriction[MAT_ROCK] = 0.06f;
	m_materials[MAT_ROCK].rfriction[MAT_DEFAULT] = m_materials[MAT_ROCK].rfriction[MAT_ROCK];

	// METAL ----------------------------------------------------------
	m_materials[MAT_METAL].density = 0.01f;
	m_materials[MAT_METAL].restitution = 0.6f;
	m_materials[MAT_METAL].cstrength = 22500.f;

	m_materials[MAT_METAL].friction[MAT_METAL] = 0.5f;
	m_materials[MAT_METAL].friction[MAT_DEFAULT] = m_materials[MAT_METAL].friction[MAT_ROCK];

	m_materials[MAT_METAL].rfriction[MAT_METAL] = 0.005f;
	m_materials[MAT_METAL].rfriction[MAT_DEFAULT] = m_materials[MAT_METAL].rfriction[MAT_ROCK];

	// HELIUM ----------------------------------------------------------
	m_materials[MAT_HELIUM].density = 0.0000001785;
	m_materials[MAT_HELIUM].restitution = 0.01f;
	m_materials[MAT_HELIUM].cstrength = 0.f;

	m_materials[MAT_HELIUM].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_CORK] = m_materials[MAT_CORK].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_WOOD] = m_materials[MAT_WOOD].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_ICE] = m_materials[MAT_ICE].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_WATER] = m_materials[MAT_WATER].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_PLASTIC] = m_materials[MAT_PLASTIC].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_CONCRETE] = m_materials[MAT_CONCRETE].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_ROCK] = m_materials[MAT_ROCK].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_METAL] = m_materials[MAT_METAL].friction[MAT_HELIUM] = 0.01f;
	m_materials[MAT_HELIUM].friction[MAT_DEFAULT] = m_materials[MAT_HELIUM].friction[MAT_ROCK];

	// DEFAULT/UNDEFINED ----------------------------------------------------------
	m_materials[MAT_DEFAULT] = m_materials[MAT_ROCK];
}
