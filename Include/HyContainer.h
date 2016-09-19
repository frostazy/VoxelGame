#pragma once

// Please don't include this file directly.
// Please include "HyVoxelLib.h"

namespace HyVoxel {

	#pragma pack(push, 8)

	template<typename T>
	class OutVector
	{
	public:
		int length = 0;
		T*	pdata = nullptr;

		virtual void Release() = 0;
	};

	struct Vector;

	struct VoxelizationInputMesh
	{
		HyVoxel::Vector*	verts = nullptr;
		int					vertNum = 0;
		int*				indices = nullptr;
		int					indexNum = 0;
	};

	struct OutputMesh
	{
		OutVector<HyVoxel::Vector>*	verts = nullptr;
		OutVector<HyVoxel::Vector>*	normals = nullptr;
		OutVector<int>*				tris = nullptr;

		~OutputMesh()
		{
			if (verts)
			{
				verts->Release();
				verts = nullptr;
			}
			if (normals)
			{
				normals->Release();
				normals = nullptr;
			}
			if (tris)
			{
				tris->Release();
				tris = nullptr;
			}
		}
	};

	#pragma pack(pop)
}
