/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

#include "matrix.h"

namespace HyVoxel
{
	/// Specifies a vertex in a polygonal mesh
	struct FQ_Vertex
	{
		/// 3D coordinates of the vertex
		double x, y, z;
		/// Vertex type. Only vertices of the same type can be collapsed during optimization.
		int type;
		/// Pointer to the mesh faces adjacent to the vertex
		int faces;
		/// Whether the vertex has been deleted
		bool deleted;
		/// An external piece of data that can be stored in the vertex
		int external;
		/// Contains contouring flags to differentiate border vertices and roaming vector vertices
		int flags;
		/// Maximum error the vertex will tolerate during optimization
		double maxerror;
	};

	/// A face in the mesh. The first three componens are the vertex indices.
	/// The fourth component is the material type.
	typedef int FQ_Face[4];

	/// An entry in a linked list of face adjacencies for a vertex
	struct CFaceLink
	{
		/// Identifier to the adjancent face
		int faceId;
		/// Pointer to the next entry in the list, -1 if it is the tail.
		int next;
	};

	/// An efficient mesh representation that features Quadratic Error function simplification using a Multiple-Choice Algorithm
	class CFastQuadrics
	{
	public:
		/// Array of vertices in the mesh
		FQ_Vertex*	vertices;
		/// Array of faces in the mesh
		FQ_Face*	faces;
		/// Array of all nodes in adjancency lists
		CFaceLink*	faceIndex;
		/// Number of vertices in the mesh
		int vertexCount;
		/// Number of faces in the mesh
		int faceCount;
		/// Number of adjancency links in the mesh
		int faceLinkCount;
		/// Set to true to skip building adjancency lists
		bool fastBuild;
	public:
		CFastQuadrics(void);
		~CFastQuadrics(void);
	public:
		/// Allocates memory for the vertex and face arrays
		void allocate(
			/// Number of vertices to be allocated
			int vertCount,
			/// Number of faces to be allocated
			int faceCount);
		
		/// Releases memory held by this mesh
		void release();

		/// Allows one to resize the container to hold vertexCount verts and faceCount faces
		void resize();

		/// Swaps one CFastQuadrics with another
		void swap(CFastQuadrics& mesh);

		/// Copies the mesh parameter to the mesh
		void copy(CFastQuadrics& mesh);

		/// Appends the mesh parameter to the mesh
		void append(CFastQuadrics& mesh);

		/// Adds a vertex to the mesh
		int addVertex(
			/// X coordinate for the new vertex
			double x,
			/// Y coordinate for the new vertex
			double y,
			/// Z coordinate for the new vertex
			double z,
			/// A type holder for the vertex. Can be chosen by the application.
			int type,
			/// Maximum simplification error the vertex will allow
			double maxerror = 1000000.0,
			/// An external application-specific buffer
			int external = 0,
			/// The vertex flags determined during contouring
			int flags = 0);

		/// Adds a triangle to the mesh
		void addFace(
			/// Index for the first vertex in the triangle
			int vid0,
			/// Index for the second vertex in the triangle
			int vid1,
			/// Index for the third vertex in the triangle
			int vid2,
			/// A m aterial indentifier. Can be chosen by the application.
			int material);
	private:
		void addFaceToIndex(FQ_Vertex& v, int fid);
	};
}