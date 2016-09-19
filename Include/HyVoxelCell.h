#pragma once

// Please don't include this file directly.
// Please include "HyVoxelLib.h"

#include <cstdlib>
#include <cstdint>

namespace HyVoxel {

	#pragma pack(push, 8)
	typedef uint8_t				VoxelType;
	const VoxelType VOXEL_HAS_MATERIAL = 0x01;
	const VoxelType VOXEL_HAS_VECTOR = 0x02;
	const VoxelType VOXEL_HAS_ROTATION = 0x04;

	typedef uint16_t			MaterialId;
	const int VOXEL_BITS_MATERIAL = 16;
	/// A special material ID that specified no information is defined for the voxel
	const int VOXEL_NIL = ~(-1 << VOXEL_BITS_MATERIAL);
	/// A special material ID that specified material information is undefined, only voxel coordinates should be considered
	const int VOXEL_ONLY_COORDS = VOXEL_NIL - 1;
	const int VOXEL_EMPTY = 0;

	// Voxel data
	struct Voxel
	{
		VoxelType type;
		MaterialId material;
		unsigned char vx;
		unsigned char vy;
		unsigned char vz;
		unsigned char rx;
		unsigned char ry;
		unsigned char rz;
	};

	/// A 64bit integer that identifies a single world octree Cell
	typedef uint64_t		CellId;

	/// Bits to encode octree X coordinate in Cell Id
	const int CELL_BITS_X = 24;
	/// Bits to encode octree Y coordinate in Cell Id
	const int CELL_BITS_Y = 12;
	/// Bits to encode octree Z coordinate in Cell Id
	const int CELL_BITS_Z = 24;

	const unsigned int MASK_4BIT = 0x0000000F;
	const unsigned int MASK_23BIT = 0x007FFFFF;
	const unsigned int MASK_11BIT = 0x000007FF;
	const CellId MASK_LOD = 0xF000000000000000;
	const CellId MASK_X = 0x0FFFFFF000000000;
	const CellId MASK_Y = 0x0000000FFF000000;
	const CellId MASK_Z = 0x0000000000FFFFFF;

	/// Compose a CellId from the location of a Cell in the octree
	inline CellId packCellId(int level, int x, int y, int z)
	{
		uint64_t _lvl = level & MASK_4BIT;
		uint64_t _x = std::abs(x) & MASK_23BIT;
		if (x < 0)
			_x |= 0x1 << (CELL_BITS_X - 1);
		uint64_t _y = std::abs(y) & MASK_11BIT;
		if (y < 0)
			_y |= 0x1 << (CELL_BITS_Y - 1);
		uint64_t _z = std::abs(z) & MASK_23BIT;
		if (z < 0)
			_z |= 0x1 << (CELL_BITS_Z - 1);

		return _lvl << (CELL_BITS_X + CELL_BITS_Y + CELL_BITS_Z) |
			_x << (CELL_BITS_Y + CELL_BITS_Z) |
			_y << (CELL_BITS_Z) |
			_z;
	}

	/// Obtain location of a Cell in the octree from its CellId identifier
	inline void unpackCellId(CellId id, int &level, int &x, int &y, int &z)
	{
		level = (int)((id & MASK_LOD) >> (CELL_BITS_X + CELL_BITS_Y + CELL_BITS_Z));

		uint32_t _x = (uint32_t)((id & MASK_X) >> (CELL_BITS_Y + CELL_BITS_Z));
		uint32_t _y = (uint32_t)((id & MASK_Y) >> (CELL_BITS_Z));
		uint32_t _z = (uint32_t)(id & MASK_Z);

		x = (_x & MASK_23BIT) * ((_x & (0x1 << (CELL_BITS_X - 1))) ? -1 : 1);
		y = (_y & MASK_11BIT) * ((_y & (0x1 << (CELL_BITS_Y - 1))) ? -1 : 1);
		z = (_z & MASK_23BIT) * ((_z & (0x1 << (CELL_BITS_Z - 1))) ? -1 : 1);
	}

	/// Return the index into the cell blocks for a given x, y, z offset.
	template <int xDim, int yDim, int zDim, typename IdxType> IdxType getCellBlockIndex(int x, int y, int z)
	{
		return IdxType((xDim * yDim * z) + (yDim * x) + y);
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	class VoxelData
	{
	public:
		enum { cSize = xDim * yDim * zDim };
		enum { cAxes = 3 };
		enum { cDimX = xDim };
		enum { cDimY = yDim };
		enum { cDimZ = zDim };
		struct Index
		{
			Index() : index(0) {}
			Index(int idx) : index(static_cast<IdxType>(idx)) {}
			Index(int x, int y, int z) : index(IdxType(getCellBlockIndex<xDim, yDim, zDim, IdxType>(x, y, z))) {}
			operator IdxType() const
			{
				return index;
			}
			IdxType index;
		};

		void setVoxel(const Index& index, const Voxel& voxel);
		void getVoxel(const Index& index, Voxel& voxel) const;

		void setType(const Index& index, VoxelType type);
		void addType(const Index& index, VoxelType type);
		void removeType(const Index& index, VoxelType type);
		VoxelType getType(const Index& index) const;
		bool hasType(const Index& index, VoxelType type) const;

		void setMaterial(const Index& index, MaterialId material);
		void setVector(const Index& index, double dx, double dy, double dz);
		void setVector(const Index& index, unsigned char dx, unsigned char dy, unsigned char dz);
		void setRotation(const Index& index, double dx, double dy, double dz);
		void setRotation(const Index& index, unsigned char dx, unsigned char dy, unsigned char dz);
		MaterialId getMaterial(const Index& index) const;
		void getVector(const Index& index, double& dx, double& dy, double& dz) const;
		void getVector(const Index& index, float& dx, float& dy, float& dz) const;
		void getVector(const Index& index, unsigned char& dx, unsigned char& dy, unsigned char& dz) const;
		void getRotation(const Index& index, double& dx, double& dy, double& dz) const;
		void getRotation(const Index& index, float& dx, float& dy, float& dz) const;
		void getRotation(const Index& index, unsigned char& dx, unsigned char& dy, unsigned char& dz) const;

		/// Initializes the blocks.
		void clear();
		void copy(const VoxelData<xDim, yDim, zDim, IdxType>& source);

		/// Raw data interfaces.
		VoxelType* getTypes() { return types; }
		MaterialId* getMaterials() { return materials; }
		const MaterialId* getMaterials() const { return materials; }
		static int getMaterialsSize() { return sizeof(MaterialId) * cSize; }
		unsigned char* getVectors() { return (unsigned char*)vectors; }
		const unsigned char* getVectors() const { return (unsigned char*)vectors; }
		static int getVectorsSize() { return sizeof(unsigned char) * cAxes * cSize; }
		unsigned char* getRotations() { return (unsigned char*)vectors; }
		const unsigned char* getRotations() const { return (unsigned char*)rotations; }
		static int getRotationsSize() { return sizeof(unsigned char) * cAxes * cSize; }
	private:
		VoxelType types[cSize];
		MaterialId materials[cSize];
		unsigned char vectors[cSize][cAxes];
		unsigned char rotations[cSize][cAxes];
	public:
		enum { cSerialSize = cSize*(sizeof(VoxelType) + sizeof(MaterialId) + cAxes + cAxes) };

		// TODO: liang
		// void serialize(unsigned char* buffer) const;
		// void unserialize(unsigned char* buffer);
	};

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::setVoxel(const Index& index, const Voxel& voxel)
	{
		types[index] = voxel.type;
		materials[index] = voxel.material;
		vectors[index][0] = voxel.vx;
		vectors[index][1] = voxel.vy;
		vectors[index][2] = voxel.vz;
		rotations[index][0] = voxel.rx;
		rotations[index][1] = voxel.ry;
		rotations[index][2] = voxel.rz;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::getVoxel(const Index& index, Voxel& voxel) const
	{
		voxel.type = types[index];
		voxel.material = materials[index];
		voxel.vx = vectors[index][0];
		voxel.vy = vectors[index][1];
		voxel.vz = vectors[index][2];
		voxel.rx = rotations[index][0];
		voxel.ry = rotations[index][1];
		voxel.rz = rotations[index][2];
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::setType(const Index& index, VoxelType type)
	{
		types[index] = type;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::addType(const Index& index, VoxelType type)
	{
		types[index] |= type;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::removeType(const Index& index, VoxelType type)
	{
		types[index] &= ~type;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	VoxelType VoxelData<xDim, yDim, zDim, IdxType>::getType(const Index& index) const
	{
		return types[index];
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	bool VoxelData<xDim, yDim, zDim, IdxType>::hasType(const Index& index, VoxelType type) const
	{
		return (types[index] & type) != 0;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::setMaterial(const Index& index, MaterialId material)
	{
		materials[index] = material;
		types[index] |= VOXEL_HAS_MATERIAL;
		/*
		if (material != VOXEL_NIL && material != VOXEL_ONLY_COORDS && material != VOXEL_EMPTY)
		types[index] |= VOXEL_HAS_MATERIAL;
		else
		types[index] &= ~VOXEL_HAS_MATERIAL;
		*/
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::setVector(const Index& index, double dx, double dy, double dz)
	{
		// dx, dy and dx range from 0 to 1. Must be converted to a 8 bit integer so they can be stored in a voxel
		vectors[index][0] = encodeVectorCoord(dx);
		vectors[index][1] = encodeVectorCoord(dy);
		vectors[index][2] = encodeVectorCoord(dz);

		types[index] |= VOXEL_HAS_VECTOR;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::setVector(const Index& index, unsigned char dx, unsigned char dy, unsigned char dz)
	{
		vectors[index][0] = dx;
		vectors[index][1] = dy;
		vectors[index][2] = dz;

		types[index] |= VOXEL_HAS_VECTOR;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::setRotation(const Index& index, double dx, double dy, double dz)
	{
		// dx, dy and dx range from 0 to 1. Must be converted to a 8 bit integer so they can be stored in a voxel
		rotations[index][0] = (unsigned char)(255.0*dx);
		rotations[index][1] = (unsigned char)(255.0*dy);
		rotations[index][2] = (unsigned char)(255.0*dz);
		types[index] |= VOXEL_HAS_VECTOR;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::setRotation(const Index& index, unsigned char dx, unsigned char dy, unsigned char dz)
	{
		rotations[index][0] = dx;
		rotations[index][1] = dy;
		rotations[index][2] = dz;
		types[index] |= VOXEL_HAS_VECTOR;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	MaterialId VoxelData<xDim, yDim, zDim, IdxType>::getMaterial(const Index& index) const
	{
		return materials[index];
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::getVector(const Index& index, double& dx, double& dy, double& dz) const
	{
		dx = decodeVectorCoord(vectors[index][0]);
		dy = decodeVectorCoord(vectors[index][1]);
		dz = decodeVectorCoord(vectors[index][2]);
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::getVector(const Index& index, float& dx, float& dy, float& dz) const
	{
		dx = (float)decodeVectorCoord(vectors[index][0]);
		dy = (float)decodeVectorCoord(vectors[index][1]);
		dz = (float)decodeVectorCoord(vectors[index][2]);
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::getVector(const Index& index, unsigned char& dx, unsigned char& dy, unsigned char& dz) const
	{
		dx = vectors[index][0];
		dy = vectors[index][1];
		dz = vectors[index][2];
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::getRotation(const Index& index, double& dx, double& dy, double& dz) const
	{
		dx = rotations[index][0] / 255.0;
		dy = rotations[index][1] / 255.0;
		dz = rotations[index][2] / 255.0;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::getRotation(const Index& index, float& dx, float& dy, float& dz) const
	{
		dx = rotations[index][0] / 255.0f;
		dy = rotations[index][1] / 255.0f;
		dz = rotations[index][2] / 255.0f;
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::getRotation(const Index& index, unsigned char& dx, unsigned char& dy, unsigned char& dz) const
	{
		dx = rotations[index][0];
		dy = rotations[index][1];
		dz = rotations[index][2];
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::clear()
	{
		memset(types, 0, cSize * sizeof(VoxelType));
		//memset(materials, 0, cSize*sizeof(MaterialId));
		//memset(vectors, 0, 3*cSize);
		//memset(rotations, 0, 3*cSize);
	}

	template <int xDim, int yDim, int zDim, typename IdxType>
	void VoxelData<xDim, yDim, zDim, IdxType>::copy(const VoxelData<xDim, yDim, zDim, IdxType>& source)
	{
		memcpy(types, source.types, cSize * sizeof(VoxelType));
		memcpy(materials, source.materials, cSize * sizeof(MaterialId));
		memcpy(vectors, source.vectors, 3 * cSize);
		memcpy(rotations, source.rotations, 3 * cSize);
	}

	// TODO: liang
	// rename
	/// Number of voxels along one Cell Dim
	const int BLOCK_DIMENSION = 40;
	/// Number of voxels a Cell will overlap to its neighbors
	const int BLOCK_MARGIN = 2;
	/// Actual number of voxels in one Cell Dim once the margin is considered
	const int BLOCK_SIZE = (BLOCK_DIMENSION + 2 * BLOCK_MARGIN);

	// TODO: liang
	// rename
	typedef VoxelData<BLOCK_DIMENSION, BLOCK_DIMENSION, BLOCK_DIMENSION, unsigned short> BlockVoxelData;
	typedef VoxelData<BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE, unsigned int>	ContourVoxelData;

	const int CELL_WIDTH_WITH_MARGIN1 = BLOCK_DIMENSION + 2;
	typedef VoxelData<CELL_WIDTH_WITH_MARGIN1, CELL_WIDTH_WITH_MARGIN1, CELL_WIDTH_WITH_MARGIN1, unsigned int>	CellVoxelMargin1;

	#pragma pack(pop)
}
