#include "HyVoxelPrivatePCH.h"

#include "Contour.h"

#include "Voxel/HyVoxelImpl.h"

#include "Common/FastQuadrics.h"
#include "Common/MatrixSSE.h"
#include "Common/UnionFind.h"
#include "Common/HyContainerImpl.h"

using namespace HyVoxel::Algebra;

namespace HyVoxel {

CMaterialLibrary gMaterialLibrary;

static int getMaterialMedium(CMaterial& materialDefinition, int level)
{
	int medium = materialDefinition.medium;
	/*
	if (medium == CellMesh::MEDIUM_FOLIAGE && level >= BILLBOARD_MAX_LEVEL_TYPE1 - 2)
	medium = CellMesh::MEDIUM_SOLID;
	*/
	return medium;
}


ContourThreadContext::ContourThreadContext()
{
	// Create voxel work buffer
	voxelData = VF_NEW ContourVoxelData();

	// Create work buffer to store vertices
	points = VF_ALLOC(Vector, BLOCK_CUBE_SIZE);

	// Create work buffer for contouring point index
	pointIndex = VF_ALLOC(int, BLOCK_CUBE_SIZE);

	// Allocate work buffer for octree nodes
	nodes = VF_ALLOC(OctreeNode, CONTOUR_OCTREE_SIZE*CONTOUR_OCTREE_SIZE*CONTOUR_OCTREE_SIZE);

	// Allocate work buffer for quads
	quadIndex = VF_ALLOC(int, 3 * 4 * BLOCK_CUBE_SIZE);

	// Allocate work buffer for quadtypes
	quadTypes = VF_ALLOC(int, 3 * 4 * BLOCK_CUBE_SIZE);

	// Allocate work buffer for vertex types
	vertexTypes = VF_ALLOC(int, BLOCK_CUBE_SIZE);

	// Allocate work buffer for edgesigns
	edgeSign = VF_ALLOC(bool, 3 * BLOCK_CUBE_SIZE);
	materials = VF_ALLOC(MaterialId, 3 * BLOCK_CUBE_SIZE);

	// static simplification buffers
	comp_fVert = VF_ALLOC(bool, BLOCK_CUBE_SIZE);
	comp_sVert = VF_ALLOC(bool, BLOCK_CUBE_SIZE);
	comp_lVert = VF_ALLOC(bool, BLOCK_CUBE_SIZE);
	comp_vPlanes = VF_ALLOC(bool, 3 * 4 * BLOCK_CUBE_SIZE);
	comp_mMap = VF_ALLOC(unsigned short, BLOCK_CUBE_SIZE);

	// these need to be aligned for SIMD ops
	comp_quadrics = VF_ALLOC_ALIGNED(QEFMatrixSSE, BLOCK_CUBE_SIZE);
	comp_qPlanes = VF_ALLOC_ALIGNED(double, 3 * 4 * 4 * BLOCK_CUBE_SIZE);
}

ContourThreadContext::~ContourThreadContext()
{
	VF_DELETE voxelData;
	VF_FREE(points);
	VF_FREE(pointIndex);
	VF_FREE(nodes);
	VF_FREE(quadIndex);
	VF_FREE(quadTypes);
	VF_FREE(vertexTypes);
	VF_FREE(edgeSign);
	VF_FREE(materials);

	VF_FREE(comp_fVert);
	VF_FREE(comp_sVert);
	VF_FREE(comp_lVert);
	VF_FREE(comp_vPlanes);
	VF_FREE(comp_mMap);

	// must use aligned-free
	VF_FREE_ALIGNED(comp_quadrics);
	VF_FREE_ALIGNED(comp_qPlanes);
}

void HyVoxelInterfaceImpl::Contour(CellId cellId, IContourThreadContext* contourThreadContext,
	HyVoxel::OutputMesh& outMesh)
{
	ContourThreadContext& tc = static_cast<ContourThreadContextImpl*>(contourThreadContext)->ctx;

	CellMesh newCellMesh;
	GenerateCellMesh(&tc, &gMaterialLibrary, &newCellMesh, true);

	const float DECIMETER_TO_CENTIMETER = 10.0f;
	const float scale = CELL_WIDTH * (1 << CELL_LOD_MIN) * DECIMETER_TO_CENTIMETER;

	int level, cellX, cellY, cellZ;
	unpackCellId(cellId, level, cellX, cellY, cellZ);
	HyVoxel::Vector offset(cellX*scale, cellZ*scale, cellY*scale);

	OutVectorImpl<HyVoxel::Vector>* overts = new OutVectorImpl<HyVoxel::Vector>();
	OutVectorImpl<HyVoxel::Vector>* onormals = new OutVectorImpl<HyVoxel::Vector>();
	OutVectorImpl<int>* otris = new OutVectorImpl<int>();

	std::vector<HyVoxel::Vector>& verts = *overts;
	std::vector<HyVoxel::Vector>& normals = *onormals;
	std::vector<int>& tris = *otris;

	CFastQuadrics& fq = newCellMesh.fqs[CellMesh::MEDIUM_SOLID];
	for (int itFace = 0; itFace < fq.faceCount; ++itFace)
	{
		FQ_Face& face = fq.faces[itFace];
		FQ_Vertex& vfV0 = fq.vertices[face[0]];
		FQ_Vertex& vfV1 = fq.vertices[face[1]];
		FQ_Vertex& vfV2 = fq.vertices[face[2]];

		HyVoxel::Vector v0(vfV0.x, vfV0.z, vfV0.y);
		HyVoxel::Vector v1(vfV1.x, vfV1.z, vfV1.y);
		HyVoxel::Vector v2(vfV2.x, vfV2.z, vfV2.y);

		v0 = v0*scale + offset;
		v1 = v1*scale + offset;
		v2 = v2*scale + offset;

		HyVoxel::Vector n = HyVoxel::Vector::CrossProduct(v2 - v0, v1 - v0);
		n.Normalize();

		verts.push_back(v0); normals.push_back(n);
		tris.push_back(verts.size() - 1);
		verts.push_back(v1); normals.push_back(n);
		tris.push_back(verts.size() - 1);
		verts.push_back(v2); normals.push_back(n);
		tris.push_back(verts.size() - 1);
	}

	overts->SetOutVectorValues(); outMesh.verts = overts;
	onormals->SetOutVectorValues(); outMesh.normals = onormals;
	otris->SetOutVectorValues(); outMesh.tris = otris;
}

int contourEx(
	int level,
	ContourVoxelData* data,
	CMaterialLibrary* materialLibrary,
	Vector* points,
	int* pointIndex,
	int* quadIndex,
	int* vertexType,
	bool* edgeSign,
	MaterialId* materials,
	int& vertexCount)
{
	int quadCount = 0;
	vertexCount = 0;
	int lastVoxelIndex = BLOCK_SIZE - BLOCK_MARGIN - 1;

	// Initialize index for points
	memset(pointIndex, 0xFF, BLOCK_CUBE_SIZE*sizeof(int));

	// Contour
	for (int zi = 0; zi < BLOCK_SIZE - 1; zi++)
	{
		for (int xi = 0; xi < BLOCK_SIZE - 1; xi++)
		{
			for (int yi = 0; yi < BLOCK_SIZE - 1; yi++)
			{
				// Compute index for voxel in buffer
				ContourVoxelData::Index idx(xi, yi, zi);

				// Read voxel
				int mv = 0;
				if (data->hasType(idx, VOXEL_HAS_MATERIAL))
				{
					mv = data->getMaterial(idx);
				}

				// Perform Dual Contouring.
				// Inspect three edges for a voxel.
				for (int edge = 0; edge < 3; edge++)
				{
					int mxi = xi + VoxelEdgeEndpoints[edge][0];
					int myi = yi + VoxelEdgeEndpoints[edge][1];
					int mzi = zi + VoxelEdgeEndpoints[edge][2];

					if ((edge == 0) && ((yi == 0) || (zi == 0)) ||
						(edge == 1) && ((xi == 0) || (zi == 0)) ||
						(edge == 2) && ((xi == 0) || (yi == 0)))
					{
						continue;
					}

					// Get edge endpoint value for the voxel
					ContourVoxelData::Index idx2(mxi, myi, mzi);

					// Get material for edge endpoint voxel
					int mv2 = 0;
					if (data->hasType(idx2, VOXEL_HAS_MATERIAL))
					{
						mv2 = data->getMaterial(idx2);
					}

					if (mv != mv2)
					{
						// frontier is with air
						if (mv == 0 || mv2 == 0)
						{
							// Pick crossing's material
							int materialId = (mv != 0) ? mv : mv2;

							// Make sure the material is within bounds
							if (materialId >= materialLibrary->materialCount)
							{
								continue;
							}

							// Quad will be produced by taking the four neighboring voxels to the edge
							for (int i = 0; i < 4; i++)
							{
								int nxi = xi + VoxelEdgeLinks[edge][i][0];
								int nyi = yi + VoxelEdgeLinks[edge][i][1];
								int nzi = zi + VoxelEdgeLinks[edge][i][2];

								// Get neighbor voxel and its floating point
								ContourVoxelData::Index bidx(nxi, nyi, nzi);

								// Set coordinates to vertex
								int vid = pointIndex[bidx];
								if (vid == -1)
								{
									double dx, dy, dz;
									if (data->hasType(bidx, VOXEL_HAS_VECTOR))
									{
										data->getVector(bidx, dx, dy, dz);
									}
									else
									{
										dx = 0.5;
										dy = 0.5;
										dz = 0.5;
									}

									// Compute point coordinates within cell
									float px = (float)(((nxi + dx) - BLOCK_MARGIN) / BLOCK_DIMENSION);
									float py = (float)(((nyi + dy) - BLOCK_MARGIN) / BLOCK_DIMENSION);
									float pz = (float)(((nzi + dz) - BLOCK_MARGIN) / BLOCK_DIMENSION);

									vid = vertexCount;
									pointIndex[bidx] = vid;
									vertexType[vid] = 0;
									vertexCount++;

									Vector& v = points[vid];
									v.x = px;
									v.y = py;
									v.z = pz;
								}
								quadIndex[4 * quadCount + i] = vid;
								int& vtype = vertexType[vid];

								if (nxi < BLOCK_MARGIN)
								{
									vtype |= 0x1000;
									if ((nxi == BLOCK_MARGIN - 1) && (xi == BLOCK_MARGIN))
									{
										vtype |= 0x001;
									}
								}
								else if ((nxi == lastVoxelIndex) && (xi == lastVoxelIndex + 1))
								{
									vtype |= 0x010;
								}
								else if (nxi > lastVoxelIndex)
								{
									vtype |= 0x1000;
									if ((nxi == lastVoxelIndex + 1) && (xi == lastVoxelIndex + 1) && (edge != 0))
									{
										vtype |= 0x100;
									}
								}

								if (nyi < BLOCK_MARGIN)
								{
									vtype |= 0x1000;
									if ((nyi == BLOCK_MARGIN - 1) && (yi == BLOCK_MARGIN))
									{
										vtype |= 0x002;
									}
								}
								else if ((nyi == lastVoxelIndex) && (yi == lastVoxelIndex + 1))
								{
									vtype |= 0x020;
								}
								else if (nyi > lastVoxelIndex)
								{
									vtype |= 0x1000;
									if ((nyi == lastVoxelIndex + 1) && (yi == lastVoxelIndex + 1) && (edge != 1))
									{
										vtype |= 0x200;
									}
								}

								if (nzi < BLOCK_MARGIN)
								{
									vtype |= 0x1000;
									if ((nzi == BLOCK_MARGIN - 1) && (zi == BLOCK_MARGIN))
									{
										vtype |= 0x004;
									}
								}
								else if ((nzi == lastVoxelIndex) && (zi == lastVoxelIndex + 1))
								{
									vtype |= 0x040;
								}
								else if (nzi > lastVoxelIndex)
								{
									vtype |= 0x1000;
									if ((nzi == lastVoxelIndex + 1) && (zi == lastVoxelIndex + 1) && (edge != 2))
									{
										vtype |= 0x400;
									}
								}
							}

							edgeSign[quadCount] = mv == 0;

							// Set material info
							materials[quadCount] = materialId;
							quadCount++;
						}
						// See if there is a material boundary in the voxel data or a seam
						else
						{
							int medium = materialLibrary->materialIndex[mv].medium;
							int medium2 = materialLibrary->materialIndex[mv2].medium;

							if (medium != medium2)
							{
								// Quad will be produced by taking the four neighboring voxels to the edge
								int quads[4];
								for (int i = 0; i < 4; i++)
								{
									int nxi = xi + VoxelEdgeLinks[edge][i][0];
									int nyi = yi + VoxelEdgeLinks[edge][i][1];
									int nzi = zi + VoxelEdgeLinks[edge][i][2];

									// Get neighbor voxel and its floating point
									ContourVoxelData::Index bidx(nxi, nyi, nzi);

									// Set coordinates to vertex
									int vid = pointIndex[bidx];
									if (vid == -1)
									{
										double dx, dy, dz;
										if (data->hasType(bidx, VOXEL_HAS_VECTOR))
										{
											data->getVector(bidx, dx, dy, dz);
										}
										else
										{
											dx = 0.5;
											dy = 0.5;
											dz = 0.5;
										}

										// Compute point coordinates within cell
										float px = ((float)(nxi + dx) - BLOCK_MARGIN) / (BLOCK_SIZE - 2 * BLOCK_MARGIN);
										float py = ((float)(nyi + dy) - BLOCK_MARGIN) / (BLOCK_SIZE - 2 * BLOCK_MARGIN);
										float pz = ((float)(nzi + dz) - BLOCK_MARGIN) / (BLOCK_SIZE - 2 * BLOCK_MARGIN);

										vid = vertexCount;
										pointIndex[bidx] = vid;
										vertexType[vid] = 0;
										vertexCount++;
										Vector& v = points[vid];
										v.x = px;
										v.y = py;
										v.z = pz;
									}
									quads[i] = vid;
									int& vtype = vertexType[vid];

									if (nxi < BLOCK_MARGIN)
									{
										vtype |= 0x1000;
										if ((nxi == BLOCK_MARGIN - 1) && (xi == BLOCK_MARGIN))
										{
											vtype |= 0x001;
										}
									}
									else if ((nxi == lastVoxelIndex) && (xi == lastVoxelIndex + 1))
									{
										vtype |= 0x010;
									}
									else if (nxi > lastVoxelIndex)
									{
										vtype |= 0x1000;
										if ((nxi == lastVoxelIndex + 1) && (xi == lastVoxelIndex + 1) && (edge != 0))
										{
											vtype |= 0x100;
										}
									}

									if (nyi < BLOCK_MARGIN)
									{
										vtype |= 0x1000;
										if ((nyi == BLOCK_MARGIN - 1) && (yi == BLOCK_MARGIN))
										{
											vtype |= 0x002;
										}
									}
									else if ((nyi == lastVoxelIndex) && (yi == lastVoxelIndex + 1))
									{
										vtype |= 0x020;
									}
									else if (nyi > lastVoxelIndex)
									{
										vtype |= 0x1000;
										if ((nyi == lastVoxelIndex + 1) && (yi == lastVoxelIndex + 1) && (edge != 1))
										{
											vtype |= 0x200;
										}
									}

									if (nzi < BLOCK_MARGIN)
									{
										vtype |= 0x1000;
										if ((nzi == BLOCK_MARGIN - 1) && (zi == BLOCK_MARGIN))
										{
											vtype |= 0x004;
										}
									}
									else if ((nzi == lastVoxelIndex) && (zi == lastVoxelIndex + 1))
									{
										vtype |= 0x040;
									}
									else if (nzi > lastVoxelIndex)
									{
										vtype |= 0x1000;
										if ((nzi == lastVoxelIndex + 1) && (zi == lastVoxelIndex + 1) && (edge != 2))
										{
											vtype |= 0x400;
										}
									}
								}

								// frontier is between two mediums
								int matIds[2] = { mv, mv2 };
								int mediums[2] = { medium, medium2 };

								for (int side = 0; side < 2; side++)
								{
									// Pick crossing's material
									int materialId = matIds[side];

									// Make sure the material is within bounds
									if (materialId >= materialLibrary->materialCount)
									{
										continue;
									}

									quadIndex[4 * quadCount] = quads[0];
									quadIndex[4 * quadCount + 1] = quads[1];
									quadIndex[4 * quadCount + 2] = quads[2];
									quadIndex[4 * quadCount + 3] = quads[3];

									// Set quad orientation.
									edgeSign[quadCount] = side != 0;

									// Set material info
									materials[quadCount] = materialId;
									quadCount++;
								}
							}
						}
					}
				}
			}
		}
	}

	return quadCount;
}

inline double calculate_vertex_error(const CFastQuadrics* mesh, const QEFMatrixSSE* quadrics, const int id_v1, const int id_v2, bool& collapseOnV1)
{
	QEFMatrixSSE q_bar = (quadrics[id_v1] + quadrics[id_v2]);
	double error1 = vertexErrorSSE(q_bar, mesh->vertices[id_v1].x, mesh->vertices[id_v1].y, mesh->vertices[id_v1].z);
	double error2 = vertexErrorSSE(q_bar, mesh->vertices[id_v2].x, mesh->vertices[id_v2].y, mesh->vertices[id_v2].z);
	if (error1 < error2)
	{
		collapseOnV1 = true;
		return error1;
	}
	else
	{
		collapseOnV1 = false;
		return error2;
	}
}

inline double calculate_vertex_error_fv(const CFastQuadrics* mesh, const QEFMatrixSSE* quadrics, const int id_v1, const int id_v2)
{
	QEFMatrixSSE q_bar = (quadrics[id_v1] + quadrics[id_v2]);
	double error = vertexErrorSSE(q_bar, mesh->vertices[id_v1].x, mesh->vertices[id_v1].y, mesh->vertices[id_v1].z);
	return error;
}

inline double fisqrt(const double& x)
{
	// fast inverse sqrt (~5% faster than 1.0/sqrt(x)) on this machine
	static const float th = 1.5f;
	long i;
	float q, r;

	r = (float)x;
	q = r*0.5f;

	i = *(long*)&r;
	i = 0x5F3759DF - (i >> 1);
	r = *(float*)&i;
	r = r*(th - r*r*q);

	if (std::isnan(r))
		return 1.0;

	return (double)r;
}

void update_face_quadric(const CFastQuadrics* mesh, QEFMatrixSSE* quadrics, const int& vid, double* planes, bool* validPlanes)
{
	QEFMatrixSSE& m = quadrics[vid];
	m.reset();
	double *plane, M;
	int* index = &mesh->vertices[vid].faces;
	while (*index != -1)
	{
		CFaceLink& fi = mesh->faceIndex[*index];
		FQ_Face& adjF = mesh->faces[fi.faceId];
		if (adjF[0] != -1)
		{
			plane = &planes[4 * fi.faceId];
			if (!validPlanes[fi.faceId])
			{
				const FQ_Vertex& v0 = mesh->vertices[adjF[0]];
				const FQ_Vertex& v1 = mesh->vertices[adjF[1]];
				const FQ_Vertex& v2 = mesh->vertices[adjF[2]];
				const double& x0 = v0.x;	const double& y0 = v0.y;	const double& z0 = v0.z;
				const double& x1 = v1.x;	const double& y1 = v1.y;	const double& z1 = v1.z;
				const double& x2 = v2.x;	const double& y2 = v2.y;	const double& z2 = v2.z;
				double& a = plane[0] = (y1 - y0)*(z2 - z0) - (z1 - z0)*(y2 - y0);
				double& b = plane[1] = (z1 - z0)*(x2 - x0) - (x1 - x0)*(z2 - z0);
				double& c = plane[2] = (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0);
				M = fisqrt(a*a + b*b + c*c);
				a *= M;	b *= M;	c *= M;
				plane[3] = -1 * (a*x0 + b*y0 + c*z0);
				validPlanes[fi.faceId] = true;
			}
			m += QEFMatrixSSE(plane);
		}
		index = &fi.next;
	}
}

// vec1 normalized;  vec2 not
inline bool borderCP(float* vec1, float* vec2)
{
	float cpLen = VectorSSE::crossDot(vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
	if (cpLen > 0.001f)	return false;
	else				return true;
}

inline unsigned __int64 edgeID(const int id1, const int id2)
{
	union
	{
		unsigned __int64 edgeID;
		unsigned __int32 vIDs[2];
	} IDtoE;
	IDtoE.vIDs[0] = std::min(id1, id2);
	IDtoE.vIDs[1] = std::max(id1, id2);
	return IDtoE.edgeID;
}


static bool validateCollapse(CFastQuadrics* mesh, FQ_Vertex* vSource, FQ_Vertex* vDest, const unsigned short* matMap = NULL, const int edgeType = 0)
{
	// we can't collapse an edge if...
	// ... it actually has no common edge
	bool commonEdge = false;
	float collapseDir[3];
	if (matMap != NULL)
	{
		VectorSSE::normalize2((float)vSource->x, (float)vSource->y, (float)vSource->z, (float)vDest->x, (float)vDest->y, (float)vDest->z, collapseDir[0], collapseDir[1], collapseDir[2]);
	}

	TVFMap<unsigned __int64, int> edgeMap;
	double oldNorm[3], newNorm[3], dotp, area;
	int* fidx = &vSource->faces;
	while (*fidx != -1)
	{
		CFaceLink& fl = mesh->faceIndex[*fidx];
		FQ_Face& f = mesh->faces[fl.faceId];
		if (f[0] != -1)
		{
			FQ_Vertex* v1 = &mesh->vertices[f[0]];
			FQ_Vertex* v2 = &mesh->vertices[f[1]];
			FQ_Vertex* v3 = &mesh->vertices[f[2]];
			if (!((v1 == vSource || v2 == vSource || v3 == vSource) && (v1 == vDest || v2 == vDest || v3 == vDest)))
			{
				VectorSSE::cross3(v1->x, v1->y, v1->z, v2->x, v2->y, v2->z, v3->x, v3->y, v3->z, oldNorm[0], oldNorm[1], oldNorm[2]);

				// ... it will affect a boundary other than the boundary being tested
				if (v1 == vSource)
				{
					v1 = vDest;
					if (edgeType != 0)
					{
						int vf = v2->flags & 0xFF;
						if (vf > 0 && vf != edgeType)
						{
							return false;
						}
						vf = v3->flags & 0xFF;
						if (vf > 0 && vf != edgeType)
						{
							return false;
						}
					}
				}
				else if (v2 == vSource)
				{
					v2 = vDest;
					if (edgeType != 0)
					{
						int vf = v1->flags & 0xFF;
						if (vf > 0 && vf != edgeType)
						{
							return false;
						}
						vf = v3->flags & 0xFF;
						if (vf > 0 && vf != edgeType)
						{
							return false;
						}
					}
				}
				else
				{
					v3 = vDest;
					if (edgeType != 0)
					{
						int vf = v1->flags & 0xFF;
						if (vf > 0 && vf != edgeType)
						{
							return false;
						}
						vf = v2->flags & 0xFF;
						if (vf > 0 && vf != edgeType)
						{
							return false;
						}
					}
				}

				VectorSSE::cross3(v1->x, v1->y, v1->z, v2->x, v2->y, v2->z, v3->x, v3->y, v3->z, newNorm[0], newNorm[1], newNorm[2]);

				// ... if we flip the normal
				dotp = VectorSSE::dot(oldNorm[0], oldNorm[1], oldNorm[2], newNorm[0], newNorm[1], newNorm[2]);
				if (dotp < 0.0)
				{
					return false;
				}

				// ... if we create a null
				area = VectorSSE::dot(newNorm[0], newNorm[1], newNorm[2], newNorm[0], newNorm[1], newNorm[2]);
				if (area < 1.0e-8)
				{
					return false;
				}

				// ... if we disturb a material boundary
				if (matMap != NULL)
				{
					if (v1 == vDest)
					{
						if (matMap[f[1]] >= 0x8000)
						{
							unsigned __int64 eid = edgeID(f[0], f[1]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								float edgeDir[3];
								VectorSSE::normalize2((float)v2->x, (float)v2->y, (float)v2->z, (float)v1->x, (float)v1->y, (float)v1->z, edgeDir[0], edgeDir[1], edgeDir[2]);
								if (!borderCP(collapseDir, edgeDir))	return false;
							}
						}
						if (matMap[f[2]] >= 0x8000)
						{
							unsigned __int64 eid = edgeID(f[0], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								float edgeDir[3];
								VectorSSE::normalize2((float)v3->x, (float)v3->y, (float)v3->z, (float)v1->x, (float)v1->y, (float)v1->z, edgeDir[0], edgeDir[1], edgeDir[2]);
								if (!borderCP(collapseDir, edgeDir))	return false;
							}
						}
					}
					else if (v2 == vDest)
					{
						if (matMap[f[0]] >= 0x8000)
						{
							unsigned __int64 eid = edgeID(f[0], f[1]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								float edgeDir[3];
								VectorSSE::normalize2((float)v1->x, (float)v1->y, (float)v1->z, (float)v2->x, (float)v2->y, (float)v2->z, edgeDir[0], edgeDir[1], edgeDir[2]);
								if (!borderCP(collapseDir, edgeDir))	return false;
							}
						}
						if (matMap[f[2]] >= 0x8000)
						{
							unsigned __int64 eid = edgeID(f[1], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								float edgeDir[3];
								VectorSSE::normalize2((float)v3->x, (float)v3->y, (float)v3->z, (float)v2->x, (float)v2->y, (float)v2->z, edgeDir[0], edgeDir[1], edgeDir[2]);
								if (!borderCP(collapseDir, edgeDir))	return false;
							}
						}
					}
					else
					{
						if (matMap[f[0]] >= 0x8000)
						{
							unsigned __int64 eid = edgeID(f[0], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								float edgeDir[3];
								VectorSSE::normalize2((float)v1->x, (float)v1->y, (float)v1->z, (float)v3->x, (float)v3->y, (float)v3->z, edgeDir[0], edgeDir[1], edgeDir[2]);
								if (!borderCP(collapseDir, edgeDir))	return false;
							}
						}
						if (matMap[f[1]] >= 0x8000)
						{
							unsigned __int64 eid = edgeID(f[1], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								float edgeDir[3];
								VectorSSE::normalize2((float)v2->x, (float)v2->y, (float)v2->z, (float)v3->x, (float)v3->y, (float)v3->z, edgeDir[0], edgeDir[1], edgeDir[2]);
								if (!borderCP(collapseDir, edgeDir))	return false;
							}
						}
					}
				}
			}
			else
			{
				commonEdge = true;
				if (matMap != NULL)
				{
					if (v1 != vSource && v1 != vDest && matMap[f[0]] >= 0x8000)
					{
						if (v2 == vSource)
						{
							unsigned __int64 eid = edgeID(f[0], f[1]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								return false;
							}
						}
						else
						{
							unsigned __int64 eid = edgeID(f[0], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								return false;
							}
						}
					}
					else if (v2 != vSource && v2 != vDest && matMap[f[1]] >= 0x8000)
					{
						if (v1 == vSource)
						{
							unsigned __int64 eid = edgeID(f[0], f[1]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								return false;
							}
						}
						else
						{
							unsigned __int64 eid = edgeID(f[1], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								return false;
							}
						}
					}
					else if (v3 != vSource && v3 != vDest && matMap[f[2]] >= 0x8000)
					{
						if (v1 == vSource)
						{
							unsigned __int64 eid = edgeID(f[0], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								return false;
							}
						}
						else
						{
							unsigned __int64 eid = edgeID(f[1], f[2]);
							if (edgeMap.count(eid) == 0)
							{
								edgeMap[eid] = f[3] & 0xFF;
							}
							else if ((f[3] & 0xFF) != edgeMap[eid])
							{
								return false;
							}
						}
					}
				}
				if (edgeType != 0)
				{
					if (v1 == vSource || v1 == vDest)
					{
						if (v2 == vSource || v2 == vDest)
						{
							if ((v3->flags & 0xFF) != 0)
							{
								return false;
							}
						}
						else
						{
							if ((v2->flags & 0xFF) != 0)
							{
								return false;
							}
						}
					}
					else if (v2 == vSource || v2 == vDest)
					{
						if (v1 == vSource || v1 == vDest)
						{
							if ((v3->flags & 0xFF) != 0)
							{
								return false;
							}
						}
						else
						{
							if ((v1->flags & 0xFF) != 0)
							{
								return false;
							}
						}
					}
					else
					{
						if (v1 == vSource || v1 == vDest)
						{
							if ((v2->flags & 0xFF) != 0)
							{
								return false;
							}
						}
						else
						{
							if ((v1->flags & 0xFF) != 0)
							{
								return false;
							}
						}
					}
				}
			}
		}
		fidx = &fl.next;
	}
	return commonEdge;
}


void contractMeshEdgeSeams(CFastQuadrics* mesh, QEFMatrixSSE* quadrics, bool* validPlanes, double* planes, const unsigned short* matMap, double* matError, const CMaterialLibrary* materialLibrary)
{
	const unsigned short* mm[2] = { NULL, matMap };
	static const unsigned short allVF[18] =
	{
		0x01,			// xlo
		0x10,			// xhi
		0x02,			// ylo
		0x20,			// yhi
		0x04,			// zlo
		0x40,			// zhi
		0x01 | 0x02,		// xlo ylo
		0x01 | 0x20,		// xlo yhi
		0x10 | 0x02,		// xhi ylo
		0x10 | 0x20,		// xhi yhi
		0x01 | 0x04,		// xlo zlo
		0x01 | 0x40,		// xlo zhi
		0x10 | 0x04,		// xhi zlo
		0x10 | 0x40,		// xhi zhi
		0x02 | 0x04,		// ylo zlo
		0x02 | 0x40,		// ylo zhi
		0x20 | 0x04,		// yhi zlo
		0x20 | 0x40		// yhi zhi
	};
	// sort verts by seam type
	// note 1:
	//		judging by recent dumps, the ordering is consistent (i.e. face links, etc. match) on each side of the boundary
	//		no need to sort when this is the case.
	// note 2:
	//		we CANNOT merge a X border to a XY border (or similar)
	//		as not all X and XY verts are common in border cells
	//		(this will cause holes)
	std::deque<int> seamVerts[18];
	for (int i = 0; i < mesh->vertexCount; i++)
	{
		const FQ_Vertex& v = mesh->vertices[i];
		if (v.deleted || v.faces == -1)
		{
			continue;
		}
		const int vflags = v.flags & 0xFF;
		for (int j = 0; j < 18; j++)
		{
			if (vflags == allVF[j])
			{
				seamVerts[j].push_back(i);
			}
		}
	}
	for (int i = 17; i >= 0; i--)
	{
		// we cannot contract away a single edge (this will cause holes)
		if (seamVerts[i].size() < 3)
		{
			continue;
		}
		bool updateVQ = true;
		// for this seam type...
		while (!seamVerts[i].empty())
		{
			// find the next pairing of vert ids
			while (mesh->vertices[seamVerts[i].front()].deleted)
			{
				seamVerts[i].pop_front();
			}
			if (seamVerts[i].empty())
			{
				continue;
			}
			int id1 = seamVerts[i].front();
			seamVerts[i].pop_front();
			if (seamVerts[i].empty())
			{
				continue;
			}
			while (mesh->vertices[seamVerts[i].front()].deleted)
			{
				seamVerts[i].pop_front();
			}
			if (seamVerts[i].empty())
			{
				continue;
			}
			int* index = &mesh->vertices[id1].faces;
			bool nfound = false;
			std::deque<int>::iterator nptr;
			while (*index != -1 && !nfound)
			{
				CFaceLink& link = mesh->faceIndex[*index];
				FQ_Face& f = mesh->faces[link.faceId];
				if (f[0] != -1)
				{
					for (nptr = seamVerts[i].begin(); nptr != seamVerts[i].end(); ++nptr)
					{
						if (!mesh->vertices[*nptr].deleted)
						{
							if (f[0] == *nptr || f[1] == *nptr || f[2] == *nptr)
							{
								nfound = true;
								break;
							}
						}
					}
				}
				index = &link.next;
			}
			if (!nfound)
			{
				updateVQ = true;
				continue;
			}
			int id2 = *nptr;
			if (nptr != seamVerts[i].begin())
			{
				seamVerts[i].erase(nptr);
				seamVerts[i].push_front(id2);
			}

			FQ_Vertex* v1 = &mesh->vertices[id1];
			FQ_Vertex* v2 = &mesh->vertices[id2];

			// check max edge len (if necessary)
			if (materialLibrary != NULL)
			{
				const double edgeLen = VectorSSE::lengthSqr(v1->x - v2->x, v1->y - v2->y, v1->z - v2->z);
				if (edgeLen > materialLibrary->materialIndex[std::min(v1->type, v2->type)].maxEdgeLength)
				{
					updateVQ = true;
					continue;
				}
			}

			// update quadric info
			if (updateVQ)
			{
				update_face_quadric(mesh, quadrics, id1, planes, validPlanes);
			}
			update_face_quadric(mesh, quadrics, id2, planes, validPlanes);
			QEFMatrixSSE vnQuadric = quadrics[id1] + quadrics[id2];

			// check mat borders (direction v1->v2)
			bool collapseOK = false;
			if (matMap[id1] == matMap[id2] || matMap[id1] < 0x8000)
			{
				const double vError = vertexErrorSSE(vnQuadric, v2->x, v2->y, v2->z);
				if (vError < matError[id1])
				{
					collapseOK = validateCollapse(mesh, v1, v2, mm[matMap[id2] >= 0x8000], allVF[i]);
				}
			}
			// try (v2->v1 on failure)
			if (!collapseOK)
			{
				if (matMap[id1] == matMap[id2] || matMap[id2] < 0x8000)
				{
					const double vError = vertexErrorSSE(vnQuadric, v1->x, v1->y, v1->z);
					if (vError < matError[id2])
					{
						collapseOK = validateCollapse(mesh, v2, v1, mm[matMap[id1] >= 0x8000], allVF[i]);
					}
				}
				if (collapseOK)
				{
					// next id is id2, replace it with id1
					seamVerts[i].front() = id1;
					FQ_Vertex* t = v1;
					v1 = v2;
					v2 = t;
					int ti = id1;
					id1 = id2;
					id2 = ti;
				}
			}
			// collapse verts
			if (collapseOK)
			{
				updateVQ = true;
				// 'v1' is deleted
				int* index = &v1->faces;
				while (*index != -1)
				{
					CFaceLink& link = mesh->faceIndex[*index];
					FQ_Face& f = mesh->faces[link.faceId];
					if (f[0] != -1)
					{
						validPlanes[link.faceId] = false;
						for (int j = 0; j < 3; j++)
						{
							if (f[j] == id1)
								f[j] = id2;
						}
						if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
							f[0] = -1;
					}
					index = &link.next;
				}
				// append adjacency lists
				index = &v2->faces;
				CFaceLink* tail = NULL;
				while (*index != -1)
				{
					tail = &mesh->faceIndex[*index];
					index = &(tail->next);
				}
				if (tail != NULL)
				{
					tail->next = v1->faces;
				}
				else
				{
					v2->faces = v1->faces;
				}
				// update 'absorbed' vertex
				v1->deleted = true;
				v1->external = id2;
				v2->flags |= v1->flags;
				v1->type = std::min(v1->type, v2->type);
				matError[id2] = std::min(matError[id1], matError[id2]);
			}
			// no collapse, quadrics already computed in prior step
			else
			{
				updateVQ = false;
			}
		}
	}
}

inline void updateFaceLinks(CFastQuadrics* mesh, const int idSrc, const int idDest, bool* lockVerts, bool* skipVerts, bool* validPlanes)
{
	int *index, *childIndex;
	index = &mesh->vertices[idDest].faces;
	CFaceLink* prevLink = NULL;
	while (*index != -1)
	{
		CFaceLink& fl = mesh->faceIndex[*index];
		FQ_Face& f = mesh->faces[fl.faceId];

		// redundant face
		if (f[0] == -1)
		{
			if (prevLink != NULL)
			{
				prevLink->next = fl.next;
			}
			else
			{
				mesh->vertices[idDest].faces = fl.next;
			}
		}
		// deleted/repeated face
		else if (((f[0] == idDest && f[1] == idDest) ||
			(f[0] == idDest && f[2] == idDest) ||
			(f[1] == idDest && f[2] == idDest)) ||
			(f[0] == idSrc || f[1] == idSrc || f[2] == idSrc) &&
			(f[0] == idDest || f[1] == idDest || f[2] == idDest))
		{
			if (prevLink != NULL)
			{
				prevLink->next = fl.next;
			}
			else
			{
				mesh->vertices[idDest].faces = fl.next;
			}
			for (int vc = 0; vc < 3; vc++)
			{
				CFaceLink* prevCLink = NULL;
				childIndex = &mesh->vertices[f[vc]].faces;
				while (*childIndex != -1)
				{
					// update neighbour indices
					CFaceLink& cfl = mesh->faceIndex[*childIndex];
					FQ_Face& fc = mesh->faces[cfl.faceId];
					if (cfl.faceId != fl.faceId)
					{
						prevCLink = &cfl;
						for (int fvi = 0; fvi < 3; fvi++)
						{
							if (fc[fvi] == idSrc)
							{
								fc[fvi] = idDest;
							}
						}
					}
					childIndex = &cfl.next;
				}
			}
			f[0] = -1;
		}
		// ok face
		else
		{
			prevLink = &fl;
			for (int fvi = 0; fvi < 3; fvi++)
			{
				lockVerts[f[fvi]] = true;
				skipVerts[f[fvi]] = false;
			}
		}
		index = &fl.next;
	}

	// signal vertices and quadrics in need of update
	mesh->vertices[idDest].external = -2;
	index = &mesh->vertices[idSrc].faces;
	while (*index != -1)
	{
		CFaceLink& fl = mesh->faceIndex[*index];
		FQ_Face& f = mesh->faces[fl.faceId];

		if (f[0] != -1)
		{
			validPlanes[fl.faceId] = false;
			for (int fvi = 0; fvi < 3; fvi++)
			{
				mesh->vertices[f[fvi]].external = -2;
			}
		}
		index = &fl.next;
	}
}

void contractMeshEdgeTest(CFastQuadrics* mesh, const CMaterialLibrary* materialLibrary, const int level, const int medium, const double maxError, const bool onlyBoundaries, ContourThreadContext* tc)
{
	// sanity check
	if (mesh->faceCount == 0)
	{
		return;
	}
	static const int maxHomoLOD = materialLibrary->maxHomogeneousLOD;

	QEFMatrixSSE* quadrics = tc->comp_quadrics;
	bool* fixVerts = tc->comp_fVert;
	bool* skipVerts = tc->comp_sVert;
	bool* lockVerts = tc->comp_lVert;
	bool* validPlanes = tc->comp_vPlanes;
	double* planes = tc->comp_qPlanes;
	unsigned short* matMap = tc->comp_mMap;

	static const int QPSize = 3 * 4 * 4 * BLOCK_CUBE_SIZE;
	bool ownsExtraBuffers = false;
	bool* edgeVert = NULL;
	double* matError = NULL;
	if (2 * mesh->vertexCount < BLOCK_CUBE_SIZE &&
		4 * mesh->faceCount + mesh->vertexCount < QPSize)
	{
		edgeVert = tc->comp_fVert + mesh->vertexCount;
		matError = tc->comp_qPlanes + 4 * mesh->faceCount;
	}
	else
	{
		ownsExtraBuffers = true;
		edgeVert = VF_ALLOC(bool, mesh->vertexCount);
		matError = VF_ALLOC(double, mesh->vertexCount);
	}
	memset(edgeVert, 0, mesh->vertexCount*sizeof(bool));
	memset(validPlanes, 0, mesh->faceCount*sizeof(bool));
	memset(fixVerts, 0, mesh->vertexCount*sizeof(bool));
	memset(matMap, 0, mesh->vertexCount*sizeof(unsigned short));

	// invalidate vertex quadrics; set perma-fixed verts based on contouring type and/or material properties
	int *index;
	for (int i = 0; i < mesh->vertexCount; i++)
	{
		FQ_Vertex& v = mesh->vertices[i];
		v.external = -2;
		if (medium == CellMesh::MEDIUM_FOLIAGE)
		{
			if (level <= maxHomoLOD)	matError[i] = 0.0005;
			else						matError[i] = 1.0e10;
			v.type = 65536;
			index = &(v.faces);
			while (*index != -1)
			{
				CFaceLink& link = mesh->faceIndex[*index];
				FQ_Face& f = mesh->faces[link.faceId];
				if (f[0] != -1)
				{
					unsigned char m = (f[3] & 0xFF);
					v.type = std::min(v.type, (int)m);
				}
				index = &link.next;
			}
		}
		else if (level <= maxHomoLOD)
		{
			TVFSet<unsigned short> vertMats;
			unsigned short& material = matMap[i];
			double& merror = matError[i];
			merror = 1.0e20;
			index = &(v.faces);
			bool blend = true;
			while (*index != -1)
			{
				CFaceLink& link = mesh->faceIndex[*index];
				FQ_Face& f = mesh->faces[link.faceId];
				if (f[0] != -1)
				{
					unsigned short m = (f[3] & 0xFF);
					if (vertMats.count(m) == 0)
					{
						vertMats.insert(m);
						material += m;
						merror = std::min(merror, materialLibrary->materialIndex[m].simplificationError);
						blend &= (materialLibrary->materialIndex[m].blend && !materialLibrary->materialIndex[m].homogeneous);
					}
				}
				index = &link.next;
			}
			merror *= 0.01;
			merror += maxError;
			v.type = vertMats.empty() ? 0 : (*vertMats.begin());	// min id
			if (!blend && vertMats.size() > 1) material |= 0x8000;
		}
		else
		{
			double& merror = matError[i];
			merror = 1.0e20;
			v.type = 65536;
			index = &(v.faces);
			while (*index != -1)
			{
				CFaceLink& link = mesh->faceIndex[*index];
				FQ_Face& f = mesh->faces[link.faceId];
				if (f[0] != -1)
				{
					unsigned char m = (f[3] & 0xFF);
					v.type = std::min(v.type, (int)m);
					merror = std::min(merror, materialLibrary->materialIndex[m].simplificationError);
				}
				index = &link.next;
			}
			merror *= 0.05;
			merror += maxError;
		}
	}

	// seam simplification
	if (medium != CellMesh::MEDIUM_FOLIAGE)
	{
		if (medium == CellMesh::MEDIUM_WATER)
		{
			// check maxEdgeLength
			contractMeshEdgeSeams(mesh, quadrics, validPlanes, planes, matMap, matError, materialLibrary);
		}
		else
		{
			// ignore maxEdgeLength
			contractMeshEdgeSeams(mesh, quadrics, validPlanes, planes, matMap, matError, NULL);
		}
	}

	// no cell simplification
	if (onlyBoundaries)
	{
		if (ownsExtraBuffers)
		{
			VF_FREE(edgeVert);
			VF_FREE(matError);
		}
		return;
	}

	// cell simplification
	bool secondPass = false;
	memset(fixVerts, 0, mesh->vertexCount*sizeof(bool));
	const unsigned short* mm[2] = { NULL, matMap };
	if (level <= maxHomoLOD && medium == CellMesh::MEDIUM_SOLID)
	{
		secondPass = true;

		// we initially configure all verts to be skipped;  we then selectively unlock verts neighbouring fixed vertices
		// this grows the simplification around fixed points, rather than simplifying pseudo-randomly [ (0,0,0) -> (1,1,1) ] and becoming stuck on them
		memset(skipVerts, 1, mesh->vertexCount*sizeof(bool));
		for (int i = 0; i < mesh->vertexCount; i++)
		{
			FQ_Vertex& v = mesh->vertices[i];
			if (v.deleted)
			{
				fixVerts[i] = true;
			}
			else if (v.faces == -1)
			{
				fixVerts[i] = true;
				v.deleted = true;
				v.external = -1;
			}
			else if (v.flags != 0)
			{
				fixVerts[i] = true;
				index = &v.faces;
				while (*index != -1)
				{
					CFaceLink& link = mesh->faceIndex[*index];
					FQ_Face& f = mesh->faces[link.faceId];
					if (f[0] != -1)
					{
						for (int j = 0; j < 3; j++)
						{
							skipVerts[f[j]] = false;
							edgeVert[f[j]] = true;
						}
					}
					index = &link.next;
				}
			}
		}
	}
	else
	{
		memset(skipVerts, 0, mesh->vertexCount*sizeof(bool));
		for (int i = 0; i < mesh->vertexCount; i++)
		{
			FQ_Vertex& v = mesh->vertices[i];
			if (v.faces == -1)
			{
				fixVerts[i] = true;
				v.deleted = true;
				v.external = -1;
			}
			else if (v.deleted || v.flags != 0)
			{
				fixVerts[i] = true;
			}
		}
	}

	// main loop
	static const int maxIterations = 150; // for safety, convergence takes around ~20-40 iterations per cell
	bool locked, collapseOnV1, collapseOK, mcollapseOnV1, done = false;
	double minError, error;
	int vidx, midx, iter = 0;

	// first pass ----------------------------------------------------------------------------------------------------------------------
	//	for each vertex, find the minimum error edge collapse
	//	if said collapse does not invalidate geometry, perform it
	//		pros: fast, only checks one set of geom per active vert
	//		cons: will choose the same move repeatedly if the neighbourhood doesn't change, and thus, miss some pretty obvious collapses
	while (!done && ++iter < maxIterations)
	{
		done = true;
		memset(lockVerts, 0, mesh->vertexCount*sizeof(bool));
		for (int i = 0; i < mesh->vertexCount; i++)
		{
			if (lockVerts[i] || fixVerts[i] || skipVerts[i])
			{
				continue;
			}

			// find neighbour vertex with lowest error;  if found, attempt to collapse the common edge
			FQ_Vertex& v = mesh->vertices[i];
			locked = false;
			minError = matError[i];
			midx = -1;
			index = &v.faces;
			while (*index != -1 && !locked)
			{
				CFaceLink& link = mesh->faceIndex[*index];
				FQ_Face& f = mesh->faces[link.faceId];
				if (f[0] != -1)
				{
					for (int vi = 0; vi < 3; vi++)
					{
						vidx = f[vi];
						if (vidx == i)
						{
							continue;
						}

						if (lockVerts[vidx])
						{
							locked = true;
							break;
						}

						// it may be better to distinguish verts which have no moves in the prior step
						// from vert pairs we don't wish to retest (right now, this acts as a global skip,
						// but it would perhaps be better to:
						//		a) skip verts with no moves should their neighbourhood remain static
						//		b) skip test *pairs* which failed geom tests in the previous step
						if (skipVerts[vidx])
						{
							continue;
						}

						if (matMap[i] != matMap[vidx] && (matMap[i] >= 0x8000 || matMap[vidx] >= 0x8000))	// contract any no-blend edges
																											//if (matMap[i] != matMap[vidx])													// contract only identical edges
						{
							continue;
						}

						FQ_Vertex& vn = mesh->vertices[vidx];

						// old maxEdgeLength test for water medium
						if (medium == CellMesh::MEDIUM_WATER)
						{
							int matIndex = f[3] & 0xFF;
							double edgeLength = VectorSSE::dot(v.x - vn.x, v.y - vn.y, v.z - vn.z, v.x - vn.x, v.y - vn.y, v.z - vn.z);
							if (edgeLength > materialLibrary->materialIndex[matIndex].maxEdgeLength)
							{
								continue;
							}
						}

						if (v.external == -2)
						{
							update_face_quadric(mesh, quadrics, i, planes, validPlanes);
							v.external = -1;
						}
						if (vn.external == -2)
						{
							update_face_quadric(mesh, quadrics, vidx, planes, validPlanes);
							vn.external = -1;
						}
						if (fixVerts[vidx])
						{
							collapseOnV1 = false;
							error = calculate_vertex_error_fv(mesh, quadrics, vidx, i);
							if (error < matError[i])
							{
								if (error > 0.0)	error *= 0.5;
								else				error *= 2.0;
							}
						}
						else
						{
							error = calculate_vertex_error(mesh, quadrics, i, vidx, collapseOnV1);
							if (error < matError[i])
							{
								if (collapseOnV1 && edgeVert[i])
								{
									if (!edgeVert[vidx])
									{
										double dfc = abs(v.x - 0.5) + abs(v.y - 0.5) + abs(v.z - 0.5);
										if (error > 0.0)	error += error*dfc;
										else				error -= error*dfc;
									}
									else
									{
										double dfc1 = abs(v.x - 0.5) + abs(v.y - 0.5) + abs(v.z - 0.5);
										double dfc2 = abs(vn.x - 0.5) + abs(vn.y - 0.5) + abs(vn.z - 0.5);
										if (dfc2 < dfc1)
										{
											if (error > 0.0)	error += error;
											else				error -= error;
										}
									}
								}
								else if (!collapseOnV1 && edgeVert[vidx])
								{
									if (!edgeVert[i])
									{
										double dfc = abs(vn.x - 0.5) + abs(vn.y - 0.5) + abs(vn.z - 0.5);
										if (error > 0.0)	error += error*dfc;
										else				error -= error*dfc;
									}
									else
									{
										double dfc1 = abs(vn.x - 0.5) + abs(vn.y - 0.5) + abs(vn.z - 0.5);
										double dfc2 = abs(v.x - 0.5) + abs(v.y - 0.5) + abs(v.z - 0.5);
										if (dfc2 < dfc1)
										{
											if (error > 0.0)	error += error;
											else				error -= error;
										}
									}
								}
							}
						}
						if (error < minError)
						{
							midx = vidx;
							minError = error;
							mcollapseOnV1 = collapseOnV1;
						}
					}
				}
				index = &link.next;
			}

			if (locked)
			{
				continue;
			}

			// check for collapse
			if (midx != -1)
			{
				if (mcollapseOnV1)			collapseOK = validateCollapse(mesh, &mesh->vertices[midx], &mesh->vertices[i], mm[matMap[i] >= 0x8000]);	// collapsing to v1 coordinates, check the faces of v2
				else						collapseOK = validateCollapse(mesh, &mesh->vertices[i], &mesh->vertices[midx], mm[matMap[i] >= 0x8000]);	// collapsing to v2 coordinates, check the faces of v1

				if (collapseOK)
				{
					// signal that another iteration is required
					done = false;

					int ids[2]; // source, dest
					if (mcollapseOnV1)
					{
						ids[0] = midx;
						ids[1] = i;
					}
					else
					{
						ids[0] = i;
						ids[1] = midx;
					}

					FQ_Vertex& vSrc = mesh->vertices[ids[0]];
					FQ_Vertex& vDest = mesh->vertices[ids[1]];

					fixVerts[ids[0]] = true;
					vSrc.deleted = true;
					vSrc.external = ids[1];
					vDest.flags |= vSrc.flags;
					matError[ids[1]] = std::min(matError[i], matError[midx]);

					// append adjacency lists
					index = &(vDest.faces);
					CFaceLink* tail = NULL;
					while (*index != -1)
					{
						tail = &mesh->faceIndex[*index];
						index = &(tail->next);
					}
					if (tail != NULL)
					{
						tail->next = vSrc.faces;
					}
					else
					{
						vDest.faces = vSrc.faces;
					}

					// update face links.  prune deleted faces.  swap references from vSrc to vDest.
					updateFaceLinks(mesh, ids[0], ids[1], lockVerts, skipVerts, validPlanes);
				}
				else
				{
					skipVerts[i] = true;
				}
			}
			else
			{
				skipVerts[i] = true;
			}
		}
	}

	// second pass ----------------------------------------------------------------------------------------------------------------------
	//	only *required* for simplification which begins with the border edges
	//	(otherwise, disconnected pieces will always be skipped)
	//	for each vertex, find the minimum error ALLOWED edge collapse and perform it
	//		pros: the collapse will always be allowed, thus we shouldn't 'miss' any edges
	//		cons: we need to check all geom with error < minerror, which can be expensive
	if (secondPass)
	{
		done = false;
		iter = 0;
		memset(skipVerts, 0, mesh->vertexCount*sizeof(bool));
		while (!done && ++iter < maxIterations)
		{
			done = true;
			memset(lockVerts, 0, mesh->vertexCount*sizeof(bool));
			for (int i = 0; i < mesh->vertexCount; i++)
			{
				if (lockVerts[i] || fixVerts[i] || skipVerts[i])
				{
					continue;
				}

				// find neighbour vertex with lowest error;  if found, attempt to collapse the common edge
				FQ_Vertex& v = mesh->vertices[i];
				locked = false;
				minError = matError[i];
				midx = -1;
				index = &v.faces;
				const unsigned short* MM = mm[matMap[i] >= 0x8000];
				while (*index != -1 && !locked)
				{
					CFaceLink& link = mesh->faceIndex[*index];
					FQ_Face& f = mesh->faces[link.faceId];
					if (f[0] != -1)
					{
						int matIndex = f[3] & 0xFF;
						for (int vi = 0; vi < 3; vi++)
						{
							vidx = f[vi];
							if (vidx == i)
							{
								continue;
							}

							if (lockVerts[vidx])
							{
								locked = true;
								break;
							}

							if (skipVerts[vidx])
							{
								continue;
							}

							if (matMap[i] != matMap[vidx] && (matMap[i] >= 0x8000 || matMap[vidx] >= 0x8000))
								//if (matMap[i] != matMap[vidx])
							{
								continue;
							}

							FQ_Vertex& vn = mesh->vertices[vidx];

							// old maxEdgeLength test for water medium
							if (medium == CellMesh::MEDIUM_WATER)
							{
								int matIndex = f[3] & 0xFF;
								double edgeLength = VectorSSE::dot(v.x - vn.x, v.y - vn.y, v.z - vn.z, v.x - vn.x, v.y - vn.y, v.z - vn.z);
								if (edgeLength > materialLibrary->materialIndex[matIndex].maxEdgeLength)
								{
									continue;
								}
							}

							if (v.external == -2)
							{
								update_face_quadric(mesh, quadrics, i, planes, validPlanes);
								v.external = -1;
							}
							if (vn.external == -2)
							{
								update_face_quadric(mesh, quadrics, vidx, planes, validPlanes);
								vn.external = -1;
							}
							if (fixVerts[vidx])
							{
								collapseOnV1 = false;
								error = calculate_vertex_error_fv(mesh, quadrics, vidx, i);
								if (error < matError[i])
								{
									collapseOK = validateCollapse(mesh, &mesh->vertices[i], &mesh->vertices[vidx], MM);
									if (!collapseOK)	continue;
								}
							}
							else
							{
								error = calculate_vertex_error(mesh, quadrics, i, vidx, collapseOnV1);
								if (error < matError[i])
								{
									if (collapseOnV1)			collapseOK = validateCollapse(mesh, &mesh->vertices[vidx], &mesh->vertices[i], MM);
									else						collapseOK = validateCollapse(mesh, &mesh->vertices[i], &mesh->vertices[vidx], MM);
									if (!collapseOK)			continue;
								}
							}
							if (error < minError)
							{
								midx = vidx;
								minError = error;
								mcollapseOnV1 = collapseOnV1;
							}
						}
					}
					index = &link.next;
				}

				if (locked)
				{
					continue;
				}

				// check for collapse
				if (midx != -1)
				{
					// signal that another iteration is required
					done = false;

					int ids[2]; // source, dest
					if (mcollapseOnV1)
					{
						ids[0] = midx;
						ids[1] = i;
					}
					else
					{
						ids[0] = i;
						ids[1] = midx;
					}

					FQ_Vertex& vSrc = mesh->vertices[ids[0]];
					FQ_Vertex& vDest = mesh->vertices[ids[1]];

					fixVerts[ids[0]] = true;
					vSrc.deleted = true;
					vSrc.external = ids[1];
					vDest.flags |= vSrc.flags;
					matError[ids[1]] = std::min(matError[i], matError[midx]);

					// append adjacency lists
					index = &(vDest.faces);
					CFaceLink* tail = NULL;
					while (*index != -1)
					{
						tail = &mesh->faceIndex[*index];
						index = &(tail->next);
					}
					if (tail != NULL)
					{
						tail->next = vSrc.faces;
					}
					else
					{
						vDest.faces = vSrc.faces;
					}

					// update face links.  prune deleted faces.  swap references from vSrc to vDest.
					updateFaceLinks(mesh, ids[0], ids[1], lockVerts, skipVerts, validPlanes);
				}
				else
				{
					skipVerts[i] = true;
				}
			}
		}
	}

	if (ownsExtraBuffers)
	{
		VF_FREE(edgeVert);
		VF_FREE(matError);
	}
}

struct isoMeshData
{
	bool	isolated;				// records whether or not this mesh is isolated from the cell boundaries
	int		vertexCount;			// records the number of unique vertices in this mesh
	int		idxC, idxS, idxE;		// tracks the current, initial and final indices in the vertex map
	float	xMin, yMin, zMin;		// records the lower bounds of the mesh

	isoMeshData()
	{
		this->isolated = true;
		this->vertexCount = 0;
		this->xMin = FLT_MAX;
		this->yMin = FLT_MAX;
		this->zMin = FLT_MAX;
	};
};

static bool trintersect(const float& v1x, const float& v1y, const float& v1z,	// vert1
	const float& v2x, const float& v2y, const float& v2z,	// vert2
	const float& v3x, const float& v3y, const float& v3z,	// vert3
	const float& rx, const float& ry, const float& rz,		// ray origin
	const float& rdx, const float& rdy, const float& rdz,	// ray direction
	float &int_dist)										// intersection distance
{
	// c.f. Moller–Trumbore intersection
	static const float epsilon = 0.000001f;
	const float e1x = v2x - v1x;
	const float e1y = v2y - v1y;
	const float e1z = v2z - v1z;
	const float e2x = v3x - v1x;
	const float e2y = v3y - v1y;
	const float e2z = v3z - v1z;

	float px = rdy*e2z - rdz*e2y;
	float py = rdz*e2x - rdx*e2z;
	float pz = rdx*e2y - rdy*e2x;

	float det = px*e1x + py*e1y + pz*e1z;

	if (det > -epsilon && det < epsilon)
	{
		return false;
	}

	float inv_det = 1.f / det;

	float tx = rx - v1x;
	float ty = ry - v1y;
	float tz = rz - v1z;

	float u = inv_det * (tx*px + ty*py + tz*pz);

	if (u < 0.f || u > 1.f)
	{
		return false;
	}

	float qx = ty*e1z - tz*e1y;
	float qy = tz*e1x - tx*e1z;
	float qz = tx*e1y - ty*e1x;

	float v = inv_det * (rdx*qx + rdy*qy + rdz*qz);

	if (v < 0.f || v > 1.f)
	{
		return false;
	}

	int_dist = inv_det * (qx*e2x + qy*e2y + qz*e2z);

	if (int_dist > epsilon)
	{
		return true;
	}
	else
	{
		return false;
	}
}

static void removePockets(CFastQuadrics* mesh)
{
	// allocate classification buffers :: we require one buffer to index isolated mesh vertices, and another to track face visitation to optimize the raytests
	unsigned char* faceVisitation = (unsigned char*)calloc(mesh->faceCount, sizeof(unsigned char));
	FQ_Vertex** vertexIndices = VF_ALLOC(FQ_Vertex*, mesh->vertexCount);

	// step 1 :: classify connected meshes using UnionFind
	CUnionFind32 classify(mesh->vertexCount);
	for (int i = 0; i < mesh->vertexCount; i++)
	{
		if (mesh->vertices[i].deleted)
		{
			continue;
		}

		int* index = &mesh->vertices[i].faces;

		// this vertex has no faces, mark it for deletion regardless of whether or not we delete its associated mesh
		if (*index == -1)
		{
			mesh->vertices[i].deleted = true;
		}
		while (*index != -1)
		{
			CFaceLink& fi = mesh->faceIndex[*index];
			FQ_Face& adjFace = mesh->faces[fi.faceId];
			if (adjFace[0] != -1 && adjFace[3] != -1)
			{
				classify.Union(i, adjFace[0]);
				classify.Union(i, adjFace[1]);
				classify.Union(i, adjFace[2]);
			}
			index = &fi.next;
		}
	}

	// step 2 :: identify border types, and mesh AABBs and vertex counts
	// important note: this requires the vertex type property to be unchanged from that of contouring
	int classifiedType;
	TVFMap<int, isoMeshData*> meshInfo;
	for (int i = 0; i < mesh->vertexCount; i++)
	{
		if (mesh->vertices[i].deleted)
		{
			continue;
		}

		classifiedType = classify.Find(i);
		if (meshInfo.count(classifiedType) == 0)
		{
			meshInfo[classifiedType] = VF_NEW isoMeshData();
		}
		if (mesh->vertices[i].type != 0)
		{
			meshInfo[classifiedType]->isolated = false;
		}

		// update type lower bounds and vertex counts (only necessary if type is still labelled as isolated)
		if (meshInfo[classifiedType]->isolated)
		{
			meshInfo[classifiedType]->vertexCount++;
			meshInfo[classifiedType]->xMin = std::min(meshInfo[classifiedType]->xMin, (float)mesh->vertices[i].x);
			meshInfo[classifiedType]->yMin = std::min(meshInfo[classifiedType]->yMin, (float)mesh->vertices[i].y);
			meshInfo[classifiedType]->zMin = std::min(meshInfo[classifiedType]->zMin, (float)mesh->vertices[i].z);
		}
	}

	// step 3 :: determine array index positions
	int currentPos = 0;
	for (TVFMap<int, isoMeshData*>::iterator i = meshInfo.begin(); i != meshInfo.end(); ++i)
	{
		if (i->second->isolated)
		{
			i->second->idxS = currentPos;
			i->second->idxC = currentPos;
			currentPos += i->second->vertexCount;
			i->second->idxE = currentPos;
		}
	}

	// step 4 :: fill vertex array
	for (int i = 0; i < mesh->vertexCount; i++)
	{
		if (mesh->vertices[i].deleted)
		{
			continue;
		}

		classifiedType = classify.Find(i);
		if (meshInfo[classifiedType]->isolated)
		{
			vertexIndices[meshInfo[classifiedType]->idxC++] = &mesh->vertices[i];
		}
	}

	// step 5 :: compute intersections
	float tx, ty, tz, tNx, tNy, tNz, atNx, atNy, atNz;
	float rx, ry, rz, rDx, rDy, rDz, curDist, minDist;
	int minFaceID;
	for (TVFMap<int, isoMeshData*>::iterator i = meshInfo.begin(); i != meshInfo.end(); ++i)
	{
		if (i->second->isolated)
		{
			minDist = FLT_MAX;
			tx = -FLT_MAX;
			i->second->idxC = i->second->idxS;
			while (true)
			{
				// no suitable faces found
				if (i->second->idxC == i->second->idxE)
				{
					break;
				}

				// determine ray information
				FQ_Vertex* v = vertexIndices[i->second->idxC++];
				int* index = &v->faces;
				while (*index != -1)
				{
					CFaceLink& fi = mesh->faceIndex[*index];
					FQ_Face& adjFace = mesh->faces[fi.faceId];
					if (adjFace[0] != -1 && adjFace[3] != -1)
					{
						FQ_Vertex& v1 = mesh->vertices[adjFace[0]];
						FQ_Vertex& v2 = mesh->vertices[adjFace[1]];
						FQ_Vertex& v3 = mesh->vertices[adjFace[2]];

						// face centroid
						triCentroid(
							(float)v1.x, (float)v1.y, (float)v1.z,
							(float)v2.x, (float)v2.y, (float)v2.z,
							(float)v3.x, (float)v3.y, (float)v3.z,
							tx, ty, tz);

						// face normal
						VectorSSE::cross3(
							(float)v1.x, (float)v1.y, (float)v1.z,
							(float)v2.x, (float)v2.y, (float)v2.z,
							(float)v3.x, (float)v3.y, (float)v3.z,
							tNx, tNy, tNz);

						break;
					}
					index = &fi.next;
				}

				// a valid start face was found, raycast
				if (tx > -FLT_MAX)
				{
					atNx = abs(tNx);
					atNy = abs(tNy);
					atNz = abs(tNz);

					// compute ray origin + direction
					if (atNx > atNy)
					{
						if (atNx > atNz)
						{
							rx = i->second->xMin - 0.00001f;
							ry = ty;
							rz = tz;
							rDx = 1.f;
							rDy = 0.f;
							rDz = 0.f;
						}
						else
						{
							rx = tx;
							ry = ty;
							rz = i->second->zMin - 0.00001f;
							rDx = 0.f;
							rDy = 0.f;
							rDz = 1.f;
						}
					}
					else
					{
						if (atNy > atNz)
						{
							rx = tx;
							ry = i->second->yMin - 0.00001f;
							rz = tz;
							rDx = 0.f;
							rDy = 1.f;
							rDz = 0.f;
						}
						else
						{
							rx = tx;
							ry = ty;
							rz = i->second->zMin - 0.00001f;
							rDx = 0.f;
							rDy = 0.f;
							rDz = 1.f;
						}
					}

					// intersect each face of the mesh
					for (int j = i->second->idxS; j < i->second->idxE; j++)
					{
						FQ_Vertex* v = vertexIndices[j];
						int* index = &v->faces;
						while (*index != -1)
						{
							CFaceLink& fi = mesh->faceIndex[*index];
							FQ_Face& adjFace = mesh->faces[fi.faceId];
							if (faceVisitation[fi.faceId] == 0 && adjFace[0] != -1 && adjFace[3] != -1)
							{
								faceVisitation[fi.faceId] = 1;
								FQ_Vertex& v1 = mesh->vertices[adjFace[0]];
								FQ_Vertex& v2 = mesh->vertices[adjFace[1]];
								FQ_Vertex& v3 = mesh->vertices[adjFace[2]];
								if (trintersect(/*v1*/(float)v1.x, (float)v1.y, (float)v1.z,
									/*v2*/(float)v2.x, (float)v2.y, (float)v2.z,
									/*v3*/(float)v3.x, (float)v3.y, (float)v3.z,
									/*ray*/rx, ry, rz, rDx, rDy, rDz, curDist)
									&& curDist < minDist)
								{
									minDist = curDist;
									minFaceID = fi.faceId;
								}
							}
							index = &fi.next;
						}
					}
					break;
				}
			}

			// we have found at least one valid ray intersection, check its normal wrt the casted ray
			if (minDist < FLT_MAX)
			{
				FQ_Vertex& v1 = mesh->vertices[mesh->faces[minFaceID][0]];
				FQ_Vertex& v2 = mesh->vertices[mesh->faces[minFaceID][1]];
				FQ_Vertex& v3 = mesh->vertices[mesh->faces[minFaceID][2]];

				// face normal
				VectorSSE::cross3(
					(float)v1.x, (float)v1.y, (float)v1.z,
					(float)v2.x, (float)v2.y, (float)v2.z,
					(float)v3.x, (float)v3.y, (float)v3.z,
					tNx, tNy, tNz);

				// comp with casted ray dir;  mark mesh to be retained if the conditional is true
				if (rDx > 0.f && tNx <= 0.f)
				{
					i->second->isolated = false;
				}
				else if (rDy > 0.f && tNy <= 0.f)
				{
					i->second->isolated = false;
				}
				else if (rDz > 0.f && tNz <= 0.f)
				{
					i->second->isolated = false;
				}
			}
		}
	}

	// step 6 :: mark faces and vertices for deletion
	for (TVFMap<int, isoMeshData*>::const_iterator i = meshInfo.begin(); i != meshInfo.end(); ++i)
	{
		if (i->second->isolated)
		{
			for (int j = i->second->idxS; j < i->second->idxE; j++)
			{
				FQ_Vertex* v = vertexIndices[j];
				v->deleted = true;
				int* index = &v->faces;
				while (*index != -1)
				{
					CFaceLink& fi = mesh->faceIndex[*index];
					FQ_Face& adjFace = mesh->faces[fi.faceId];
					adjFace[0] = adjFace[3] = -1;
					index = &fi.next;
				}
			}
		}
	}

	// cleanup
	for (TVFMap<int, isoMeshData*>::iterator i = meshInfo.begin(); i != meshInfo.end(); ++i)
	{
		VF_DELETE i->second;
	}
	VF_FREE(faceVisitation);
	VF_FREE(vertexIndices);
}

void contractMesh(CFastQuadrics& mesh, CMaterialLibrary& materialLibrary, int level, int medium, double maxError, bool onlyBoundaries, ContourThreadContext* tc)
{
	contractMeshEdgeTest(&mesh, &materialLibrary, level, medium, maxError, onlyBoundaries, tc);

	//pocket removal
	removePockets(&mesh);
}

int GenerateCellMesh(
	/// A ContourThreadContext instance containing the work buffers for the contouring operation
	ContourThreadContext* tc,
	/// Material library to be used in contouring
	CMaterialLibrary* materialLibrary,
	/// The CellMesh object that will receive the polygonal output
	CellMesh* cellData,
	/// Signals whether simplify the meshes
	bool compress,
	/// Signals whether to simplify only the cell boundaries
	bool compressOnlyBoundaries
	)
{
	CellId cell = cellData->cellId;
	ContourVoxelData* data = tc->voxelData;
	Vector* points = tc->points;
	int* pointIndex = tc->pointIndex;
	int* quadIndex = tc->quadIndex;
	int* vertexTypes = tc->vertexTypes;
	bool* edgeSign = tc->edgeSign;
	MaterialId* materials = tc->materials;

	// Unpack cell coordinates
	int level, xc, yc, zc;
	unpackCellId(cell, level, xc, yc, zc);

	int vertexCount = 0;
	int quadCount = contourEx(
		level, data, materialLibrary, points,
		pointIndex, quadIndex, vertexTypes, edgeSign, materials, vertexCount);

	// See if contouring returned any quads
	if (quadCount > 0)
	{
		// Compute cell size in world units
		double scale = (1 << level)*CELL_WIDTH;

		// Classify quads
		int quadMediumCount[CellMesh::MEDIUM_MAX];
		for (int i = 0; i < CellMesh::MEDIUM_MAX; ++i)
		{
			quadMediumCount[i] = 0;
		}
		for (int i = 0; i < quadCount; i++)
		{
			int mat = materials[i];
			CMaterial& materialDefinition = materialLibrary->materialIndex[mat];
			int medium = getMaterialMedium(materialDefinition, level);
			quadMediumCount[medium]++;
		}

		// Initialize mesh
		int totalFaces = 0;
		for (int medium = 0; medium < CellMesh::MEDIUM_MAX; medium++)
		{
			if (quadMediumCount[medium] == 0)
			{
				continue;
			}

			CFastQuadrics hiResMesh;
			hiResMesh.allocate(vertexCount, 2 * quadMediumCount[medium]);
			if (hiResMesh.faces == NULL)
			{
				printf("bad alloc");
				exit(EXIT_FAILURE);
			}

			// Add vertices returned by contouring
			for (int i = 0; i < vertexCount; i++)
			{
				Vector& p = points[i];
				hiResMesh.addVertex(p.x, p.y, p.z, 0, 1.0e6, 0, vertexTypes[i]);
			}

			// Add triangles returned by contouring.
			// Each quad contributes two triangles
			for (int i = 0; i < quadCount; i++)
			{
				// Read material from contouring info
				int mat = materials[i];

				CMaterial& materialDefinition = materialLibrary->materialIndex[mat];
				int quadMedium = getMaterialMedium(materialDefinition, level);
				if (quadMedium != medium)
				{
					continue;
				}

				// Read four vertex IDs for the current quad from contouring octree
				int id0 = quadIndex[4 * i];
				int id1 = quadIndex[4 * i + 1];
				int id2 = quadIndex[4 * i + 2];
				int id3 = quadIndex[4 * i + 3];

				//if the quad is fully outside of the cell mark it for future deletion
				int materialMask = 0;
				if ((vertexTypes[id0] & vertexTypes[id1] & vertexTypes[id2] & vertexTypes[id3] & 0x1000) != 0)
				{
					materialMask = 0x100;
				}

				// Add triangles depending on edge sign
				int triangles[2][3];
				if (edgeSign[i])
				{
					triangles[0][0] = id0;
					triangles[0][1] = id1;
					triangles[0][2] = id2;

					triangles[1][0] = id0;
					triangles[1][1] = id2;
					triangles[1][2] = id3;
				}
				else
				{
					triangles[0][0] = id0;
					triangles[0][1] = id2;
					triangles[0][2] = id1;

					triangles[1][0] = id0;
					triangles[1][1] = id3;
					triangles[1][2] = id2;
				}
				for (int t = 0; t < 2; t++)
				{
					if (triangles[t][0] != triangles[t][1] &&
						triangles[t][0] != triangles[t][2] &&
						triangles[t][1] != triangles[t][2])
					{
						hiResMesh.addFace(triangles[t][0], triangles[t][1], triangles[t][2], (mat | materialMask));
					}
				}
			}

			//seams
			//set data for the seams
			for (int front = 0; front < 2; front++)
			{
				for (int axis = 0; axis < 3; axis++)
				{
					//mask for the seam
					int mask = (1 << axis) << (8 * front);

					int faceCount = 0;
					int vertexCount = 0;
					bool resetPoinIndex = true;

					//count the faces and vertices in the seam
					for (int i = 0; i < hiResMesh.faceCount; i++)
					{
						FQ_Face& f = hiResMesh.faces[i];

						if ((f[3] & 0x100) != 0)
						{
							continue;
						}

						bool inSeamFace = false;
						for (int vi = 0; vi < 3; vi++)
						{
							FQ_Vertex& v = hiResMesh.vertices[f[vi]];
							if ((v.flags & mask) != 0)
							{
								inSeamFace = true;
								break;
							}
						}

						if (inSeamFace)
						{
							if (resetPoinIndex)
							{
								memset(pointIndex, 0xFF, BLOCK_CUBE_SIZE*sizeof(int));
								resetPoinIndex = false;
							}

							faceCount++;
							for (int vi = 0; vi < 3; vi++)
							{
								int& idx = pointIndex[f[vi]];
								if (idx == -1)
								{
									vertexCount++;
									idx = 1;
								}
							}
						}
					}

					CFastQuadrics& seamMesh = cellData->seamMesh[medium][front][axis];
					if (faceCount > 0)
					{
						//allocate memory for the seam
						seamMesh.allocate(vertexCount, faceCount);

						//reset the buffer of indexes
						memset(pointIndex, 0xFF, BLOCK_CUBE_SIZE*sizeof(int));

						//add faces and vertices to the seam
						for (int i = 0; i < hiResMesh.faceCount; i++)
						{
							FQ_Face& f = hiResMesh.faces[i];
							if ((f[3] & 0x100) != 0)
							{
								continue;
							}

							bool inSeamFace = false;
							for (int vi = 0; vi < 3; vi++)
							{
								FQ_Vertex& v = hiResMesh.vertices[f[vi]];
								if ((v.flags & mask) != 0)
								{
									inSeamFace = true;
									break;
								}
							}

							if (inSeamFace)
							{
								int triangle[3];
								for (int vi = 0; vi < 3; vi++)
								{
									int& idx = pointIndex[f[vi]];
									if (idx == -1)
									{
										FQ_Vertex& v = hiResMesh.vertices[f[vi]];
										idx = seamMesh.addVertex(v.x, v.y, v.z, v.type, v.maxerror, f[vi], v.flags);
									}
									triangle[vi] = idx;
								}
								seamMesh.addFace(triangle[0], triangle[1], triangle[2], f[3]);
							}
						}
					}
					else
					{
						seamMesh.release();
					}
				}
			}

			int errorExp = std::max(0, level - CELL_LOD_MIN - 3);
			if (compress)
			{
				double error = 0.00000001*errorExp;
				contractMesh(hiResMesh, *materialLibrary, level, medium, error, compressOnlyBoundaries, tc);
			}
			else
			{
				// assign vert (mat) types here if there's no simplification (otherwise, the simplification will do this step inline)
				// Remember last material for a case where no material is known
				int lastValidMaterial = 0;

				// Iterate through all vertices in mesh
				for (int i = 0; i < hiResMesh.vertexCount; i++)
				{
					FQ_Vertex& v = hiResMesh.vertices[i];
					int* index = &v.faces;
					int matId = INT_MAX;

					// Find material of a neighboring face and use it for the vertex
					bool materialFound = false;
					while (*index != -1)
					{
						CFaceLink& fi = hiResMesh.faceIndex[*index];
						FQ_Face& adjFace = hiResMesh.faces[fi.faceId];

						if (adjFace[0] != -1 && adjFace[3] != -1)
						{
							matId = std::min(adjFace[3], matId);
							materialFound = true;
						}
						index = &fi.next;
					}
					if (materialFound)
					{
						// Store found material in the type property for the vertex
						v.type = matId;

						// Rememeber last material
						lastValidMaterial = matId;
					}
					else
					{
						// Store last found in the type property for the vertex
						v.type = lastValidMaterial;
					}
				}
			}

			// copy to cellData mesh
			if (hiResMesh.faces == NULL)
			{
				printf("Invalid Mesh");
				exit(EXIT_FAILURE);
			}

			memset(pointIndex, 0xFF, BLOCK_CUBE_SIZE*sizeof(int));

			int faceCount = 0;
			int vertexCount = 0;

			for (int i = 0; i < hiResMesh.faceCount; i++)
			{
				FQ_Face& f = hiResMesh.faces[i];
				if (f[0] != -1)
				{
					if ((f[3] & 0x100) == 0)
					{
						for (int vi = 0; vi < 3; vi++)
						{
							int& idx = pointIndex[f[vi]];
							if (idx == -1)
							{
								vertexCount++;
								idx = 1;
							}
						}

						faceCount++;
					}
					else
					{
						f[0] = -1;
						f[1] = -1;
						f[2] = -1;
						f[3] = -1;
					}
				}
			}

			CFastQuadrics& fq = cellData->fqs[medium];
			fq.allocate(vertexCount, faceCount);
			memset(pointIndex, 0xFF, BLOCK_CUBE_SIZE*sizeof(int));

			for (int i = 0; i < hiResMesh.faceCount; i++)
			{
				FQ_Face& f = hiResMesh.faces[i];
				if (f[0] != -1)
				{
					int triangle[3];
					for (int vi = 0; vi < 3; vi++)
					{
						int& idx = pointIndex[f[vi]];
						if (idx == -1)
						{
							FQ_Vertex& v = hiResMesh.vertices[f[vi]];
							idx = fq.addVertex(v.x, v.y, v.z, v.type, v.maxerror, 0, v.flags);
							v.external = idx;
						}

						triangle[vi] = idx;
					}

					fq.addFace(triangle[0], triangle[1], triangle[2], f[3]);
				}
			}

			// update the operation list
			int maxSteps = hiResMesh.vertexCount;	// setting this to a fixed constant is... unwise.
			for (int front = 0; front < 2; front++)
			{
				for (int axis = 0; axis < 3; axis++)
				{
					CFastQuadrics& seamMesh = cellData->seamMesh[medium][front][axis];
					if (seamMesh.faceCount > 0)
					{
						for (int i = 0; i < seamMesh.vertexCount; i++)
						{
							FQ_Vertex& v = seamMesh.vertices[i];

							int index = v.external;
							int steps = 0;

							while ((index >= 0) && (steps < maxSteps))
							{
								steps++;
								FQ_Vertex& nv = hiResMesh.vertices[index];
								index = nv.external;
								if (!nv.deleted)
								{
									break;
								}
							}

							if (steps < maxSteps)
							{
								v.external = index;
							}
							else
							{
								v.external = -1;
							}
						}
					}
				}
			}
			totalFaces += cellData->fqs[medium].faceCount;
		}

		return totalFaces;
	}
	else
	{
		return 0;
	}
}



}

