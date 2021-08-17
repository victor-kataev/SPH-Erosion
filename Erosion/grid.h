#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <list>
#include <omp.h>
#include <cstring>
#include <iomanip>
#include <float.h>

#include "voxel.h"

struct Triangle
{
	Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c)
		: A(a), B(b), C(c)
	{
		norm = glm::normalize(glm::cross(B - A, C - A));
	}

	glm::vec3 A;
	glm::vec3 B;
	glm::vec3 C;
	glm::vec3 norm;
};

class Grid
{
private:
	int m_DimX;
	int m_DimY;
	int m_DimZ;
	glm::vec3 m_Dim;
	int m_empty_voxels;
	int m_filled_voxels;

	Voxel* m_Grid;
	//FluidSystem m_Fluid;

	std::vector<float> surfaceParts;
	std::vector<float> fluidParts;
	std::vector<unsigned int> indices;

	struct Heightfield
	{
		size_t dimX = 0;
		size_t dimY = 0;

		unsigned char* map = nullptr;
	};

	Heightfield m_heightfield;

	void destroyGrid()
	{
		delete[] m_Grid;
		surfaceParts.clear();
		fluidParts.clear();
	}

	void createGrid(int dimx, int dimy, int dimz)
	{
		m_Grid = new Voxel[dimx * dimy * dimz];
		m_Dim = glm::vec3(dimx, dimy, dimz);
		m_DimX = m_Dim.x;
		m_DimY = m_Dim.y;
		m_DimZ = m_Dim.z;
		m_empty_voxels = 0;
		m_filled_voxels = 0;
	}

public:
	Grid(int dimX = 512, int dimY = 512, int dimZ = 512)
		: m_DimX(dimX), m_DimY(dimY), m_DimZ(dimZ)
	{
		m_Dim = glm::vec3(dimX, dimY, dimZ);
		m_Grid = new Voxel[dimX * dimY * dimZ];

		m_heightfield.dimX = 512;
		m_heightfield.dimY = 512;
		m_heightfield.map = new unsigned char[512 * 512]; //hard coded
		memset(m_heightfield.map, 0, 512 * 512);
	}

	Voxel GetVoxel(int x, int y, int z) const
	{
		return m_Grid[x + m_DimX * (y + m_DimY * z)];
	}

	void SetVoxel(int x, int y, int z, VoxelType type)
	{
		//m_Grid[x + m_DimX * (y + m_DimY * z)] = val;
		m_Grid[x + m_DimX * (y + m_DimY * z)].type = type;
		m_Grid[x + m_DimX * (y + m_DimY * z)].position = glm::vec3(x, y, z);

		//todo
	}

	void LoadHeightfield(unsigned char* img)
	{
		memcpy(m_heightfield.map, img, (size_t)m_heightfield.dimX * m_heightfield.dimY);
		HeightFieldMax();
	}

	unsigned char GetHeightfieldAt(int x, int y)
	{
		return m_heightfield.map[m_heightfield.dimY * x + y];
	}

	void HeightFieldMax() const
	{
		int max = 0;
		for (int i = 0; i < m_heightfield.dimX * m_heightfield.dimY; i++)
			if (m_heightfield.map[i] > max)
				max = m_heightfield.map[i];
		//std::cout << "Heightfield max = " << max << " " << __LINE__ << std::endl;
	}

	void genIndices()
	{
		for(int z = 0, j = m_Dim.z-1; z < m_Dim.z && j >= 0; z++, j--)
			for (int x = 0, i = m_Dim.x-1; x < m_Dim.x && i >= 0; x++, i--)
			{
				if (x + 1 < m_Dim.x && z + 1 < m_Dim.z)
				{
					indices.push_back(z * m_Dim.x + x); //position in vertex buffer
					indices.push_back(z * m_Dim.x + (x+1));
					indices.push_back((z+1) * m_Dim.x + x);
				}
				if (j - 1 >= 0 && i - 1 >= 0)
				{
					indices.push_back(j * m_Dim.x + i);
					indices.push_back(j * m_Dim.x + (i - 1));
					indices.push_back((j - 1) * m_Dim.x + i);
				}
			}
	}

	void UpdateGrid(int dimx, int dimy, int dimz)
	{
		destroyGrid();
		createGrid(dimx, dimy, dimz);
		for (int z = 0; z < m_Dim.z; z++)
			for (int x = 0; x < m_Dim.x; x++)
			{
				int y = GetHeightfieldAt(x, z);
				if (dimy <= y)
				{
					y = dimy-1;
				}
				SetVoxel(x, y, z, VoxelType::VOXEL_MAT);
				surfaceParts.push_back(x);
				surfaceParts.push_back(y);
				surfaceParts.push_back(z);

				//normals
				glm::vec3 u(0.0);
				glm::vec3 d(0.0);
				glm::vec3 r(0.0);
				glm::vec3 l(0.0);

				if (x - 1 >= 0)
					l = glm::vec3(x - (x - 1), y - GetHeightfieldAt(x - 1, z), z - z);
				if (x + 1 < m_Dim.x)
					r = glm::vec3((x + 1) - x, GetHeightfieldAt(x + 1, z) - y, z - z);
				if (z - 1 >= 0)
					u = glm::vec3(x - x, y - GetHeightfieldAt(x, z - 1), z - (z - 1));
				if (z + 1 < m_Dim.y)
					d = glm::vec3(x - x, GetHeightfieldAt(x, z + 1) - y, (z + 1) - z);

				glm::vec3 normal = glm::normalize(glm::cross(u, l) + glm::cross(u, r) + glm::cross(d, l) + glm::cross(d, r));
				surfaceParts.push_back(normal.x);
				surfaceParts.push_back(normal.y);
				surfaceParts.push_back(normal.z);
			}
		genIndices();
	}

	bool rayIntersectsTriangle(const glm::vec3& pos, const glm::vec3& dir, const Triangle tri, float* t)
	{
		float u, v;
		glm::vec3 E1 = tri.B - tri.A;
		glm::vec3 E2 = tri.C - tri.A;
		glm::vec3 N = glm::cross(E1, E2);
		//glm::vec3 N = tri.norm;
		float det = -glm::dot(dir, N);
		float invdet = 1.0 / det;
		glm::vec3 AO = pos - tri.A;
		glm::vec3 DAO = glm::cross(AO, dir);
		u = glm::dot(E2, DAO) * invdet;
		v = -glm::dot(E1, DAO) * invdet;
		*t = glm::dot(AO, N) * invdet;
		return (abs(det) >= 1e-6 && *t >= 0.0 && u >= 0.0 && v >= 0.0 && (u + v) <= 1.0);
	}

	std::vector<Triangle> getCellTriangles(glm::vec2 cellIndex)
	{
		glm::vec3 A(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1]), cellIndex[1]);
		glm::vec3 AA(cellIndex[0] + 1, GetHeightfieldAt(cellIndex[0] + 1, cellIndex[1] + 1), cellIndex[1] + 1);
		glm::vec3 B(cellIndex[0] + 1, GetHeightfieldAt(cellIndex[0] + 1, cellIndex[1]), cellIndex[1]);
		glm::vec3 C(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1] + 1), cellIndex[1] + 1);
		
		std::vector<Triangle> ret;
		ret.emplace_back(C, B, A); //ABC = normal is opposite
		ret.emplace_back(AA, B, C);

		return ret;
	}

	//3D-DDA (modified)
	glm::vec2 findAdjacentCell(const glm::vec3& pos, const glm::vec3& dir, glm::vec2 cellIndex)
	{
		float t_x, t_z;
		glm::vec2 deltaT;
		glm::vec2 cellDim(1, 1);
		glm::vec2 rayOrigin(pos.x, pos.z);
		glm::vec2 dir2D(dir.x, dir.z);

		if (dir2D[0] < 0)
		{
			deltaT[0] = -cellDim[0] / dir2D[0];
			t_x = (floor(rayOrigin[0] / cellDim[0]) * cellDim[0] - rayOrigin[0]) / dir2D[0];
		}
		else if (dir2D[0] > 0)
		{
			deltaT[0] = cellDim[0] / dir2D[0];
			t_x = ((floor(rayOrigin[0] / cellDim[0]) + 1) * cellDim[0] - rayOrigin[0]) / dir2D[0];
		}
		else
		{
			deltaT[0] = 0;
			t_x = INFINITY;
		}

		if (dir2D[1] < 0)
		{
			deltaT[1] = -cellDim[1] / dir2D[1];
			t_z = (floor(rayOrigin[1] / cellDim[1]) * cellDim[1] - rayOrigin[1]) / dir2D[1];
		}
		else if (dir2D[1] > 0)
		{
			deltaT[1] = cellDim[1] / dir2D[1];
			t_z = ((floor(rayOrigin[1] / cellDim[1]) + 1) * cellDim[1] - rayOrigin[1]) / dir2D[1];
		}
		else
		{
			deltaT[1] = 0;
			t_z = INFINITY;
		}

		if (t_x < t_z)
		{
			if (dir2D[0] < 0)
				cellIndex[0]--;
			else if (dir2D[0] > 0)
				cellIndex[0]++;
		}
		else
		{
			if (dir2D[1] < 0)
				cellIndex[1]--;
			else if (dir2D[1] > 0)
				cellIndex[1]++;
		}

		if (cellIndex[0] < 0 || cellIndex[0] >= m_DimX || cellIndex[1] < 0 || cellIndex[1] >= m_DimZ)
			cellIndex = glm::vec2(-1, -1);

		return cellIndex;
	}

	bool mappedOnTriangle(const Triangle tri, const glm::vec3& pos, const glm::vec3& dir, glm::vec3& cp, glm::vec3& norm)
	{
		float t;

		if (rayIntersectsTriangle(pos, dir, tri, &t))
		{
			//t += 0.000002;
			norm = dir;
			cp = pos + t * norm;
			return true;
		}
		return false;
	}

	bool mappedBetweenTriangles(const Triangle & tri1, const Triangle & tri2, const glm::vec3 & pos, glm::vec3& cp, glm::vec3& norm)
	{
		float t;

		if (rayIntersectsTriangle(pos, glm::normalize(tri1.norm + tri2.norm), tri1, &t))
		{
			//t += 0.000002;
			norm = glm::normalize(tri1.norm + tri2.norm);
			cp = pos + t * norm;
			return true;
		}
		else if (rayIntersectsTriangle(pos, glm::normalize(tri1.norm + tri2.norm), tri2, &t))
		{
			//t += 0.000002;
			norm = glm::normalize(tri1.norm + tri2.norm);
			cp = pos + t * norm;
			return true;
		}

		return false;
	}

	float min3(float a, float b, float c)
	{
		if (a < b)
			if (a < c)
				return a;
			else
				return c;
		else if (b < c)
			return b;
		else
			return c;
	}

	bool cornerCaseABC(const Triangle ABC, const Triangle AABC, const glm::vec3& posCurr, const glm::vec3& posNext, glm::vec2 cellIndex, glm::vec3& cp, glm::vec3& norm)
	{
		glm::vec2 posCurr2D = { posCurr.x, posCurr.z };
		glm::vec2 origin = { floor(posCurr.x), floor(posCurr.z) };
		glm::vec2 right = { origin.x + 1, origin.y };
		glm::vec2 down = { origin.x, origin.y + 1 };
		float dorigin = glm::length(origin - posCurr2D);
		float dright = glm::length(right - posCurr2D);
		float ddown = glm::length(down - posCurr2D);
		float dmin = min3(dorigin, dright, ddown);
		std::vector<Triangle> triangles;
		std::vector<Triangle> tmp;

		if (dmin == dorigin)
		{
			triangles.push_back(ABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y - 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y - 1));
			triangles.push_back(tmp[1]); //AABC
			//std::cout << "origin\n";
		}
		else if (dmin == dright)
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y - 1));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y - 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "right\n";
		}
		else //dmin == ddown
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y + 1));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y + 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "down\n";
		}

		glm::vec3 n(0.0);
		for (const auto& tr : triangles)
			n += tr.norm;
		n = glm::normalize(n);

		for (const auto& tr : triangles)
			if (mappedOnTriangle(tr, posNext, n, cp, norm))
				return true;
		return false;
	}

	bool cornerCaseAABC(const Triangle ABC, const Triangle AABC, const glm::vec3& posCurr, const glm::vec3& posNext, glm::vec2 cellIndex, glm::vec3& cp, glm::vec3& norm)
	{
		glm::vec2 posCurr2D = { posCurr.x, posCurr.z };
		glm::vec2 origin = { floor(posCurr.x) + 1, floor(posCurr.z) + 1 };
		glm::vec2 left = { origin.x - 1, origin.y };
		glm::vec2 up = { origin.x, origin.y - 1 };
		float dorigin = glm::length(origin - posCurr2D);
		float dleft = glm::length(left - posCurr2D);
		float dup = glm::length(up - posCurr2D);
		float dmin = min3(dorigin, dup, dleft);
		std::vector<Triangle> triangles;
		std::vector<Triangle> tmp;

		if (dmin == dorigin)
		{
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y + 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y + 1));
			triangles.push_back(tmp[0]); //ABC
			//std::cout << "origin AABC\n";
			//std::cout << "posCurr: " << posCurr.x << " " << posCurr.z << std::endl;
			//std::cout << "cellIndex: " << cellIndex.x << " " << cellIndex.y << std::endl;
		}
		else if (dmin == dup)
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y - 1));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y - 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "up\n";
		}
		else //dmin == dleft
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y + 1));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y + 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "left\n";
		}

		glm::vec3 n(0.0);
		for (const auto& tr : triangles)
			n += tr.norm;
		n = glm::normalize(n);

		for (const auto& tr : triangles)
			if (mappedOnTriangle(tr, posNext, n, cp, norm))
				return true;
		return false;
	}

	bool pointLaysOnTriangle(const glm::vec3 & p, Triangle & tri)
	{
		float S = fabs(glm::dot(tri.B - tri.A, tri.C - tri.A)) / 2.0f;
		float a = fabs(glm::dot(tri.B - p, tri.C - p)) / 2.0f / S;
		float b = fabs(glm::dot(tri.C - p, tri.A - p)) / 2.0f / S;
		float c = 1 - a - b;
		return a >= 0 && a <= 1 && b >= 0 && b <= 1 && c >= 0 && c <= 1;
	}

	bool pointLaysOnPlane(const glm::vec3& p, const Triangle& tri)
	{
		const float myEpsilon = 1.192092896e-05F;
		float d = glm::dot(tri.norm, tri.B);
		float res = tri.norm.x * tri.C.x + tri.norm.y * tri.C.y + tri.norm.z * tri.C.z - d;
		//std::cout << "plane res: " << res << std::endl;
		if (fabs(res) < myEpsilon)
			return true;
		return false;
	}

	bool collision(const glm::vec3& posCurr, const glm::vec3& posNext, const glm::vec3& velNext, glm::vec3 & contactP, glm::vec3 & norm)
	{
		float t; //where t is R.Origin + t*R.Dir
		glm::vec3 dir = glm::normalize(velNext);
		glm::vec3 dirOpposite = -dir;

		//std::cout << posCurr.x << " " << posCurr.y << " " << posCurr.z << std::endl;
		glm::vec2 cellIndex(floor(posCurr.x), floor(posCurr.z));
		glm::vec2 cellIndexNext(floor(posNext.x), floor(posNext.z));
		if (cellIndex[0] < 0 || cellIndex[0] >= m_DimX-1 || cellIndex[1] < 0 || cellIndex[1] >= m_DimZ-1
			|| cellIndexNext[0] < 0 || cellIndexNext[0] >= m_DimX-1 || cellIndexNext[1] < 0 || cellIndexNext[1] >= m_DimZ-1)
			return false;
		
		//if in the same cell
		if (cellIndex == cellIndexNext)
		{
			//std::cout << "in the same cell\n";
			//std::cout << "dir: " << dir.x << " " << dir.y << " " << dir.z << std::endl;
			//std::cout << "diroposite: " << dirOposite.x << " " << dirOposite.y << " " << dirOposite.z << std::endl;
			std::vector<Triangle> cellTriangles = getCellTriangles(cellIndex);
			Triangle ABC = cellTriangles[0];
			Triangle AABC = cellTriangles[1];

			glm::vec3 CPtmp;
			float dABC = INFINITY, dAABC = INFINITY;
			if (rayIntersectsTriangle(posNext, dirOpposite, ABC, &t))
			{
				CPtmp = posNext + t * dirOpposite;
				dABC = glm::length(CPtmp - posNext);
				//std::cout << "true1\n";
			}

			if (rayIntersectsTriangle(posNext, dirOpposite, AABC, &t))
			{
				CPtmp = posNext + t * dirOpposite;
				dAABC = glm::length(CPtmp - posNext);
				//std::cout << "true2\n";
			}

			
			//the cell is built out of two triangles
			//check intersection with both and try to map on appropriate one
			if (dABC < dAABC || (fabs(dABC - dAABC) < FLT_EPSILON && dABC != INFINITY)) //or dABC == dAABC
			{
				//std::cout << "ABC case\n";
				//try to map on the same triangle
				if (rayIntersectsTriangle(posNext, ABC.norm, ABC, &t))
				{
					//t += 0.000002;
					norm = ABC.norm;
					contactP = posNext + t * norm;
					//std::cout << "mapped on ABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
					return true;
				}
				else
				{
					glm::vec2 adjCellIndex = findAdjacentCell(posNext, glm::normalize(dir + ABC.norm), cellIndex);
					if (adjCellIndex == glm::vec2(-1, -1))
						return false;
					std::vector<Triangle> adjCellTriangles = getCellTriangles(adjCellIndex);
					Triangle ABC_adj = adjCellTriangles[0];
					Triangle AABC_adj = adjCellTriangles[1];
					glm::vec2 cellDir = { adjCellIndex.x - cellIndex.x, adjCellIndex.y - cellIndex.y };

					//intersection with hypotenuse
					if (cellDir.x > 0 || cellDir.y > 0)
					{
						//std::cout << "hypotenuse\n";
						if (mappedBetweenTriangles(ABC, AABC, posNext, contactP, norm))
						{
							//std::cout << "mapped between ABC AABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{

						//intersection with cathenus
						if (mappedBetweenTriangles(ABC, AABC_adj, posNext, contactP, norm))
						{
							//std::cout << "mapped between ABC AABC_adj, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC " << __LINE__ << std::endl;
							return true;
						}
					}
				}
			}
			else if (dAABC < dABC)
			{
				//std::cout << "AABC case\n";
				//try to map on the same triangle
				if (rayIntersectsTriangle(posNext, AABC.norm, AABC, &t))
				{
					//t += 0.000002;
					norm = AABC.norm;
					contactP = posNext + t * norm;
					//std::cout << "mapped on AABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
					return true;
				}
				else
				{
					glm::vec2 adjCellIndex = findAdjacentCell(posNext, glm::normalize(dir + AABC.norm), cellIndex);
					if (adjCellIndex == glm::vec2(-1, -1))
						return false;
					std::vector<Triangle> adjCellTriangles = getCellTriangles(adjCellIndex);
					Triangle ABC_adj = adjCellTriangles[0];
					Triangle AABC_adj = adjCellTriangles[1];
					glm::vec2 cellDir = { adjCellIndex.x - cellIndex.x, adjCellIndex.y - cellIndex.y };

					//intersection with hypotenuse
					if (cellDir.x < 0 || cellDir.y < 0)
					{
						//std::cout << "hypotenuse\n";
						if (mappedBetweenTriangles(ABC, AABC, posNext, contactP, norm))
						{
							//std::cout << "mapped between AABC ABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseAABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{

						//intersection with cathenus
						if (mappedBetweenTriangles(AABC, ABC_adj, posNext, contactP, norm))
						{
							//std::cout << "mapped between AABC ABC_adj, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseAABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC " << __LINE__ << std::endl;
							return true;
						}
					}
				}
			}
			//particle did not go through collision did not occur
			return false;
		}
		else
		{
			//std::cout << "in different cells\n";
			//std::cout << "cellIndex: " << cellIndex.x << " " << cellIndex.y << " | cellIndexNext: " << cellIndexNext.x << " " << cellIndexNext.y << std::endl;
			std::vector<Triangle> cellTriangles = getCellTriangles(cellIndex);
			Triangle ABC_poscurr = cellTriangles[0];
			Triangle AABC_poscurr = cellTriangles[1];
			std::vector<Triangle> cellNextTriangles = getCellTriangles(cellIndexNext);
			Triangle ABC_posnext = cellNextTriangles[0];
			Triangle AABC_posnext = cellNextTriangles[1];
			glm::vec2 cellDir = { cellIndexNext - cellIndex };

			//diagonal corner case
			if (cellDir.x != 0 && cellDir.y != 0)
			{
				if (cellDir.x + cellDir.y == -2)
				{
					if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
					{
						//std::cout << "corner case ABC diagonal " << __LINE__ << std::endl;
						return true;
					}
				}
				else if(cellDir.x + cellDir.y == 2)
				{
					if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
					{
						//std::cout << "corner case AABC diagonal " << __LINE__ << std::endl;
						return true;
					}
				}
				else
				{
					if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
					{
						//std::cout << "corner case ABC diagonal " << __LINE__ << std::endl;
						return true;
					}
				}
			}
			else
			{
				if (rayIntersectsTriangle(posNext, dirOpposite, ABC_poscurr, &t))
				{
					if (cellDir.x == -1 || cellDir.y == -1)
					{
						if (mappedBetweenTriangles(ABC_poscurr, AABC_posnext, posNext, contactP, norm))
						{
							//std::cout << "mapped between ABC_poscurr AABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{
						if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}

				}
				else if (rayIntersectsTriangle(posNext, dirOpposite, AABC_poscurr, &t))
				{
					if (cellDir.x == 1 || cellDir.y == 1)
					{
						if (mappedBetweenTriangles(AABC_poscurr, ABC_posnext, posNext, contactP, norm))
						{
							//std::cout << "mapped between AABC_poscurr ABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{
						if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}
				}
				else
				{
					glm::vec3 CPtmp;
					float dABC_posnext = INFINITY, dAABC_posnext = INFINITY;
					if (rayIntersectsTriangle(posCurr, dir, ABC_posnext, &t))
					{
						CPtmp = posCurr + t * dir;
						dABC_posnext = glm::length(CPtmp - posCurr);
					}

					if (rayIntersectsTriangle(posCurr, dir, AABC_posnext, &t))
					{
						CPtmp = posCurr + t * dir;
						dAABC_posnext = glm::length(CPtmp - posCurr);
					}


					if (dABC_posnext < dAABC_posnext || (fabs(dABC_posnext - dAABC_posnext) < FLT_EPSILON && dABC_posnext != INFINITY))
					{
						if (rayIntersectsTriangle(posNext, ABC_posnext.norm, ABC_posnext, &t))
						{
							//t += 0.000002;
							norm = ABC_posnext.norm;
							contactP = posNext + t * norm;
							//std::cout << "mapped on ABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cellDir.x == 1 || cellDir.y == 1)
						{
							if (mappedBetweenTriangles(ABC_posnext, AABC_poscurr, posNext, contactP, norm))
							{
								//std::cout << "mapped between ABC_posnext AABC_poscurr, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
								return true;
							}

							if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case AABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
						else
						{
							if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
					}
					else if (dAABC_posnext < dABC_posnext)
					{
						if (rayIntersectsTriangle(posNext, AABC_posnext.norm, AABC_posnext, &t))
						{
							//t += 0.000002;
							norm = AABC_posnext.norm;
							contactP = posNext + t * norm;
							//std::cout << "mapped on AABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cellDir.x == -1 || cellDir.y == -1)
						{
							if (mappedBetweenTriangles(AABC_posnext, ABC_poscurr, posNext, contactP, norm))
							{
								//std::cout << "mapped between AABC_posnext ABC_poscurr, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
								return true;
							}

							if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
						else
						{
							if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
					}
				}
			}
			return false;
		}
	}

	std::vector<unsigned int> GetIndices()
	{
		return indices;
	}

	size_t GetIndicesSize() const
	{
		return indices.size();
	}

	glm::vec3 GetDim() const
	{
		return m_Dim;
	}

	std::vector<float> GetSurfaceParts()
	{
		return surfaceParts;
	}

	size_t GetSurfacePartsSize() const
	{
		return surfaceParts.size();
	}

	std::vector<float> GetFluidParts()
	{
		return fluidParts;
	}

	size_t GetFluidPartsSize() const
	{
		return fluidParts.size();
	}
};

