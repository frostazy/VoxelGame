/************************************************************
* (C) Voxel Farm Inc. 2015
*/


#pragma once

#include <stdint.h>

namespace HyVoxel
{

	/****************************
	*							*
	*	Union Find				*
	*	(Fast Classification)	*
	*							*
	****************************/
	template <typename idxType>
	class CUnionFind
	{
	public:
		CUnionFind() {};
		CUnionFind(idxType N);
		~CUnionFind();
		idxType Find(idxType id);
		void Union(idxType id1, idxType id2);
		void Reset();
	protected:
		struct node
		{
			idxType parent;
			idxType rank;
		};
		node* uf;
		idxType size;
	};

	template <typename idxType>
	CUnionFind<idxType>::CUnionFind(idxType N)
	{
		this->size = N;
		this->uf = VF_ALLOC(node, N);
		for (idxType i = 0; i < N; i++)
		{
			uf[i].parent = i;
			uf[i].rank = 0;
		}
	}

	template <typename idxType>
	CUnionFind<idxType>::~CUnionFind()
	{
		VF_FREE(uf);
	}

	template <typename idxType>
	idxType CUnionFind<idxType>::Find(idxType id)
	{
		if (uf[id].parent != id)
		{
			uf[id].parent = Find(uf[id].parent);
		}
		return uf[id].parent;
	}

	template <typename idxType>
	void CUnionFind<idxType>::Union(idxType id1, idxType id2)
	{
		if (id1 == id2)
		{
			return;
		}

		idxType id1Root = Find(id1);
		idxType id2Root = Find(id2);

		if (id1Root == id2Root)
		{
			return;
		}

		if (uf[id1Root].rank < uf[id2Root].rank)
		{
			uf[id1Root].parent = id2Root;
		}

		else if (uf[id1Root].rank > uf[id2Root].rank)
		{
			uf[id2Root].parent = id1Root;
		}

		else
		{
			uf[id2Root].parent = id1Root;
			uf[id1Root].rank++;
		}
	}

	template <typename idxType>
	void CUnionFind<idxType>::Reset()
	{
		for (idxType i = 0; i < this->size; i++)
		{
			uf[i].parent = i;
			uf[i].rank = 0;
		}
	}

	typedef CUnionFind<uint16_t> CUnionFind16;
	typedef CUnionFind<uint32_t> CUnionFind32;
}
