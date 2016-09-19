/************************************************************
* (C) Voxel Farm Inc. 2015
*/
#include "HyVoxelPrivatePCH.h"

#include "HyVoxelConfig.h"
#include "FastQuadrics.h"

using namespace HyVoxel;

CFastQuadrics::CFastQuadrics(void)
{
	vertexCount = 0;
	faceCount = 0;
	faceLinkCount = 0;
	vertices = NULL;
	faces = NULL;
	faceIndex = NULL;
	fastBuild = false;
}

CFastQuadrics::~CFastQuadrics(void)
{
	release();
}

void CFastQuadrics::allocate(int vertCount, int faceCount)
{
	release();
	vertices = VF_ALLOC(FQ_Vertex, vertCount);
	faces = VF_ALLOC(FQ_Face, faceCount);
	if (fastBuild == false)	faceIndex = VF_ALLOC(CFaceLink, 3 * faceCount);
}

void CFastQuadrics::release()
{
	faceCount = 0;
	vertexCount = 0;
	faceLinkCount = 0;
	if(vertices != NULL)
	{
		VF_FREE(vertices);
	}
	if(faces != NULL)
	{
		VF_FREE(faces);
	}
	if(faceIndex != NULL)
	{
		VF_FREE(faceIndex);
	}
	vertices = NULL;
	faces = NULL;
	faceIndex = NULL;
}

void CFastQuadrics::copy(CFastQuadrics& mesh)
{
	if (mesh.vertexCount > 0 && mesh.faceCount > 0)
	{
		allocate(mesh.vertexCount, mesh.faceCount);

		vertexCount = mesh.vertexCount;
		memcpy(vertices, mesh.vertices, sizeof(FQ_Vertex)*vertexCount);
		faceCount = mesh.faceCount;
		memcpy(faces, mesh.faces, sizeof(FQ_Face)*faceCount);
		if (fastBuild == false)
		{
			faceLinkCount = mesh.faceLinkCount;
			memcpy(faceIndex, mesh.faceIndex, sizeof(CFaceLink)*faceLinkCount);
		}
		else
			faceLinkCount = 0;
	}
}

void CFastQuadrics::append(CFastQuadrics& mesh)
{
	if (mesh.vertexCount > 0 && mesh.faceCount > 0)
	{
		CFastQuadrics tempMesh;
		tempMesh.allocate(vertexCount + mesh.vertexCount, faceCount + mesh.faceCount);

		tempMesh.vertexCount = vertexCount + mesh.vertexCount;
		tempMesh.faceCount = faceCount + mesh.faceCount;
		tempMesh.faceLinkCount = faceLinkCount + mesh.faceLinkCount;

		//copy first mesh
		FQ_Vertex* currentVertex = tempMesh.vertices;
		FQ_Face* currentFace = tempMesh.faces;
		CFaceLink* currentLink = tempMesh.faceIndex;

		memcpy(currentVertex, vertices, sizeof(FQ_Vertex)*vertexCount);
		memcpy(currentFace, faces, sizeof(FQ_Face)*faceCount);
		memcpy(currentLink, faceIndex, sizeof(CFaceLink)*faceLinkCount);

		//copy second mesh
		currentVertex += vertexCount;
		currentFace += faceCount;
		currentLink += faceLinkCount;

		memcpy(currentVertex, mesh.vertices, sizeof(FQ_Vertex)*mesh.vertexCount);
		memcpy(currentFace, mesh.faces, sizeof(FQ_Face)*mesh.faceCount);
		memcpy(currentLink, mesh.faceIndex, sizeof(CFaceLink)*mesh.faceLinkCount);

		//fix the indexes
		for (int v=0; v<mesh.vertexCount; v++)
		{
			currentVertex->faces += faceLinkCount;
			currentVertex++;
		}

		for (int f=0; f<mesh.faceCount; f++)
		{
			(*currentFace)[0] += vertexCount;
			(*currentFace)[1] += vertexCount;
			(*currentFace)[2] += vertexCount;

			currentFace++;
		}

		for (int l=0; l<mesh.faceLinkCount; l++)
		{
			currentLink->faceId += faceCount;
			if (currentLink->next >= 0)
			{
				currentLink->next += faceLinkCount;
			}

			currentLink++;
		}

		//copy back
		allocate(tempMesh.vertexCount, tempMesh.faceCount);

		vertexCount = tempMesh.vertexCount;
		memcpy(vertices, tempMesh.vertices, sizeof(FQ_Vertex)*vertexCount);
		faceCount = tempMesh.faceCount;
		memcpy(faces, tempMesh.faces, sizeof(FQ_Face)*faceCount);
		faceLinkCount = tempMesh.faceLinkCount;
		memcpy(faceIndex, tempMesh.faceIndex, sizeof(CFaceLink)*faceLinkCount);
	}
}

int CFastQuadrics::addVertex(double x, double y, double z, int type, double maxerror, int external, int flags)
{
	// Get reference to new vertex in buffer
	FQ_Vertex& v = vertices[vertexCount];

	// Set vertex properties
	v.x = x;
	v.y = y;
	v.z = z;
	v.type = type;
	v.external = external;
	v.flags = flags;
	v.deleted = false;
	v.maxerror = maxerror;
	v.faces = -1;
	vertexCount++;

	return (vertexCount - 1);
}

void CFastQuadrics::addFaceToIndex(FQ_Vertex& v, int fid)
{
	CFaceLink* link = NULL;

	// Get the face list for the vertex
	int* last = &v.faces;

	// Move to the end of the list
	while (*last != -1)
	{
		link = &faceIndex[*last];
		last = &(link->next);
	}

	// Add new link to the chain
	*last = faceLinkCount;

	// Get reference to link element
	link = &faceIndex[faceLinkCount];

	// Set link data
	link->faceId = fid;
	link->next = -1;

	faceLinkCount++;
}

void CFastQuadrics::addFace(int vid0, int vid1, int vid2, int material)
{
	// Get references to vertices
	FQ_Vertex& v0 = vertices[vid0];
	FQ_Vertex& v1 = vertices[vid1];
	FQ_Vertex& v2 = vertices[vid2];

	// Get reference to next face in buffer
	FQ_Face& f = faces[faceCount];

	// Set face data
	f[0] = vid0;
	f[1] = vid1;
	f[2] = vid2;
	f[3] = material;

	// Update face indices for vertices
	if (!fastBuild)
	{
		addFaceToIndex(v0, faceCount);
		addFaceToIndex(v1, faceCount);
		addFaceToIndex(v2, faceCount);
	}

	faceCount++;
}

void CFastQuadrics::resize()
{
	vertices = (FQ_Vertex*)VF_RAWREALLOC(vertices, vertexCount*sizeof(FQ_Vertex));
	faces = (FQ_Face*)VF_RAWREALLOC(faces, faceCount*sizeof(FQ_Face));
	if (fastBuild == false)	faceIndex = (CFaceLink*)VF_RAWREALLOC(faceIndex, faceLinkCount*sizeof(CFaceLink));
}

void CFastQuadrics::swap(CFastQuadrics& mesh)
{
	FQ_Vertex* tv = this->vertices;
	this->vertices = mesh.vertices;
	mesh.vertices = tv;

	FQ_Face* tf = this->faces;
	this->faces = mesh.faces;
	mesh.faces = tf;

	int tvc = this->vertexCount;
	this->vertexCount = mesh.vertexCount;
	mesh.vertexCount = tvc;

	int tfc = this->faceCount;
	this->faceCount = mesh.faceCount;
	mesh.faceCount = tfc;

	if (this->fastBuild == false && mesh.fastBuild == false)
	{
		CFaceLink* tl = this->faceIndex;
		this->faceIndex = mesh.faceIndex;
		mesh.faceIndex = tl;

		int tlc = this->faceLinkCount;
		this->faceLinkCount = mesh.faceLinkCount;
		mesh.faceLinkCount = tlc;
	}
	else
	{
		this->fastBuild = true;
		mesh.fastBuild = true;
	}
}