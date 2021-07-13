#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <list>
#include <omp.h>
#include <cstring>

#include "voxel.h"

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
		std::cout << "Heightfield max = " << max << std::endl;
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
		/*for (int y = 0; y < m_Fluid.m_Dim.y; y++)
			for(int z = 0; z < m_Fluid.m_Dim.z; z++)
				for(int x = 0; x < m_Fluid.m_Dim.x; x++)
				{
					Voxel& fluid_voxel = m_Fluid.m_Volume[x + (int)m_Fluid.m_Dim.x * (y + (int)m_Fluid.m_Dim.y * z)];
					if (fluid_voxel.position.x >= 0 && fluid_voxel.position.x < m_Dim.x &&
						fluid_voxel.position.y >= 0 && fluid_voxel.position.y < m_Dim.y &&
						fluid_voxel.position.z >= 0 && fluid_voxel.position.z < m_Dim.z)
					{
						fluidParts.push_back((float)fluid_voxel.position.x / 100.0f);
						fluidParts.push_back((float)fluid_voxel.position.y / 100.0f);
						fluidParts.push_back((float)fluid_voxel.position.z / 100.0f);
					}
				}*/
	}

	bool rayIntersectsTriangle(const glm::vec3& pos, const glm::vec3& dir, const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, float* t, float* u, float* v, glm::vec3& N)
	{
		glm::vec3 E1 = B - A;
		glm::vec3 E2 = C - A;
		N = glm::cross(E1, E2);
		//glm::vec3 dir = glm::normalize(currParticle.velocity);
		float det = -glm::dot(dir, N);
		float invdet = 1.0 / det;
		glm::vec3 AO = pos - A;
		glm::vec3 DAO = glm::cross(AO, dir);
		*u = glm::dot(E2, DAO) * invdet;
		*v = -glm::dot(E1, DAO) * invdet;
		*t = glm::dot(AO, N) * invdet;
		return (abs(det) >= 1e-6 && *t >= 0.0 && *u >= 0.0 && *v >= 0.0 && (*u + *v) <= 1.0);
	}

	//bool testCellForCollision(const glm::ivec2& cellIndex, Voxel& currParticle, const glm::vec3& nextPos)
	//{
	//	float t, u, v; //where t is R.Origin + t*R.Dir;  u,v - barycentric coords
	//	glm::vec3 N; //normal

	//	//triangle points in 3D
	//	glm::vec3 A(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1]), cellIndex[1]);
	//	glm::vec3 B(cellIndex[0] + 1, GetHeightfieldAt(cellIndex[0] + 1, cellIndex[1]), cellIndex[1]);
	//	glm::vec3 C(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1] + 1), cellIndex[1] + 1);
	//	
	//	//first triangle of the cell
	//	if (rayIntersectsTriangle(currParticle, A, B, C, &t, &u, &v, N))
	//	{
	//		glm::vec3 contactP = currParticle.position + t * glm::normalize(currParticle.velocity); //R.Origin + t*R.Dir
	//		float currNextDistance = sqrt(pow(currParticle.position.x - nextPos.x, 2) + pow(currParticle.position.y - nextPos.y, 2) + pow(currParticle.position.z - nextPos.z, 2));
	//		float currIntersecDistance = sqrt(pow(currParticle.position.x - contactP.x, 2) + pow(currParticle.position.y - intersectionPoint.y, 2) + pow(currParticle.position.z - intersectionPoint.z, 2));
	//		//if distance between current particle and next step intersects the triangle
	//		if (currIntersecDistance <= currNextDistance)
	//		{
	//			//depth of penetration
	//			float d = sqrt(pow(intersectionPoint.x - nextPos.x, 2) + pow(intersectionPoint.y - nextPos.y, 2) + pow(intersectionPoint.z - nextPos.z, 2));
	//			glm::vec3 n = N;
	//			n = glm::normalize(n);

	//			currParticle.position = intersectionPoint;
	//			currParticle.velocity = currParticle.velocity - (1 + FluidSystem::Cr * (d / glm::length(currParticle.velocity))) * (currParticle.velocity * N) * N;
	//			return true;
	//		}
	//		else
	//		{
	//			//in order to not check the intersection with another triangle
	//			return false;
	//		}
	//	}


	//	//second triangle of the cell
	//	A = { A.x + 1, A.y, A.z + 1 };
	//	if (rayIntersectsTriangle(currParticle, A, B, C, &t, &u, &v, N))
	//	{
	//		glm::vec3 intersectionPoint = currParticle.position + t * glm::normalize(currParticle.velocity); //R.Origin + t*R.Dir
	//		float currNextDistance = sqrt(pow(currParticle.position.x - nextPos.x, 2) + pow(currParticle.position.y - nextPos.y, 2) + pow(currParticle.position.z - nextPos.z, 2));
	//		float currIntersecDistance = sqrt(pow(currParticle.position.x - intersectionPoint.x, 2) + pow(currParticle.position.y - intersectionPoint.y, 2) + pow(currParticle.position.z - intersectionPoint.z, 2));
	//		//if distance between current particle and next step intersects the triangle
	//		if (currIntersecDistance <= currNextDistance)
	//		{
	//			float d = sqrt(pow(intersectionPoint.x - nextPos.x, 2) + pow(intersectionPoint.y - nextPos.y, 2) + pow(intersectionPoint.z - nextPos.z, 2));

	//			currParticle.position = intersectionPoint;
	//			currParticle.velocity = currParticle.velocity - (1 + FluidSystem::Cr * (d / glm::length(currParticle.velocity))) * (currParticle.velocity * N) * N;
	//			return true;
	//		}
	//	}

	//	return false;
	//}

	//3D-DDA
	bool collision(const glm::vec3& posCurr, glm::vec3& posNext, glm::vec3& velNext, glm::vec3 & contactP, glm::vec3 & norm)
	{
		glm::ivec2 cellIndex(floor(posCurr.x), floor(posCurr.z));
		glm::ivec2 cellIndexNext(floor(posNext.x), floor(posNext.z));
		if (cellIndex[0] < 0 || cellIndex[0] >= m_DimX-1 || cellIndex[1] < 0 || cellIndex[1] >= m_DimZ-1
			|| cellIndexNext[0] < 0 || cellIndexNext[0] >= m_DimX-1 || cellIndexNext[1] < 0 || cellIndexNext[1] >= m_DimZ-1)
			return false;
		
		//if in the same cell
		if (cellIndex == cellIndexNext)
		{
			float t, u, v; //where t is R.Origin + t*R.Dir;  u,v - barycentric coords
			glm::vec3 N; //normal
			glm::vec3 dir = glm::normalize(velNext);

			//triangle points in 3D
			glm::vec3 A(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1]), cellIndex[1]);
			glm::vec3 AA(cellIndex[0]+1, GetHeightfieldAt(cellIndex[0]+1, cellIndex[1]+1), cellIndex[1]+1);
			glm::vec3 B(cellIndex[0] + 1, GetHeightfieldAt(cellIndex[0] + 1, cellIndex[1]), cellIndex[1]);
			glm::vec3 C(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1] + 1), cellIndex[1] + 1);

			//the cell is built of two triangles
			//check intersection with both and try to map on appropriate one
			if (rayIntersectsTriangle(posCurr, dir, C, B, A, &t, &u, &v, N))
			{
				norm = glm::normalize(N);
				//try to map on the same triangle
				if (rayIntersectsTriangle(posNext, norm, C, B, A, &t, &u, &v, N))
				{
					t += 0.0003;
					contactP = posNext + t * norm;
					return true;
				}
				else
				{	//try to map on second triangle
					if (rayIntersectsTriangle(posNext, norm, AA, B, C, &t, &u, &v, N))
					{
						t += 0.0003;
						contactP = posNext + t * norm;
						return true;
					}
				}
				//particle didn't go through collision didn not occur
				return false;
			}
			else if (rayIntersectsTriangle(posCurr, dir, AA, B, C, &t, &u, &v, N))
			{
				norm = glm::normalize(N);
				if (rayIntersectsTriangle(posNext, norm, AA, B, C, &t, &u, &v, N))
				{
					t += 0.0003;
					contactP = posNext + t * norm;
					return true;
				}
				else
				{
					if (rayIntersectsTriangle(posNext, norm, C, B, A, &t, &u, &v, N))
					{
						t += 0.0003;
						contactP = posNext + t * norm;
						return true;
					}
				}
				return false;
			}
			return false;
		}

		glm::vec2 vel2D(velNext.x, velNext.z);
		glm::vec2 rayDirection = glm::normalize(vel2D);
		glm::vec2 rayOrigin(posCurr.x, posCurr.z);
		glm::vec2 cellDim(1, 1);
		glm::vec2 deltaT;
		float t_x, t_z;

		if (rayDirection[0] < 0)
		{
			deltaT[0] = -cellDim[0] / rayDirection[0];
			t_x = (floor(rayOrigin[0] / cellDim[0]) * cellDim[0] - rayOrigin[0]) / rayDirection[0];
		}
		else if (rayDirection[0] > 0)
		{
			deltaT[0] = cellDim[0] / rayDirection[0];
			t_x = ((floor(rayOrigin[0] / cellDim[0]) + 1) * cellDim[0] - rayOrigin[0]) / rayDirection[0];
		}
		else
		{
			deltaT[0] = 0;
			t_x = INFINITY;
		}

		if (rayDirection[1] < 0)
		{
			deltaT[1] = -cellDim[1] / rayDirection[1];
			t_z = (floor(rayOrigin[1] / cellDim[1]) * cellDim[1] - rayOrigin[1]) / rayDirection[1];
		}
		else if (rayDirection[1] > 0)
		{
			deltaT[1] = cellDim[1] / rayDirection[1];
			t_z = ((floor(rayOrigin[1] / cellDim[1]) + 1) * cellDim[1] - rayOrigin[1]) / rayDirection[1];
		}
		else
		{
			deltaT[1] = 0;
			t_z = INFINITY;
		}

		float t = 0;
		std::vector<glm::ivec2> cellsToTest;
		cellsToTest.push_back(cellIndex);
		while (1)
		{
			if (t_x < t_z)
			{
				t = t_x;
				t_x += deltaT[0];
				if (rayDirection[0] < 0)
					cellIndex[0]--;
				else if(rayDirection[0] > 0)
					cellIndex[0]++;
			}
			else
			{
				t = t_z;
				t_z += deltaT[1];
				if (rayDirection[1] < 0)
					cellIndex[1]--;
				else if (rayDirection[1] > 0)
					cellIndex[1]++;
			}

			if (cellIndex == cellIndexNext)
			{
				cellsToTest.push_back(cellIndex);
				break;
			}

			if (cellIndex[0] < 0 || cellIndex[0] >= m_DimX || cellIndex[1] < 0 || cellIndex[1] >= m_DimZ)
				break;

			cellsToTest.push_back(cellIndex);
		}

		for (const auto& cell : cellsToTest)
		{
			float t, u, v; //where t is R.Origin + t*R.Dir;  u,v - barycentric coords
			glm::vec3 N; //normal
			glm::vec3 dir = glm::normalize(velNext);

			glm::vec3 A(cell[0], GetHeightfieldAt(cell[0], cell[1]), cell[1]);
			glm::vec3 AA(cell[0]+1, GetHeightfieldAt(cell[0]+1, cell[1]+1), cell[1]+1);
			glm::vec3 B(cell[0] + 1, GetHeightfieldAt(cell[0] + 1, cell[1]), cell[1]);
			glm::vec3 C(cell[0], GetHeightfieldAt(cell[0], cell[1] + 1), cell[1] + 1);

			//1st triangle
			if (rayIntersectsTriangle(posCurr, dir, C, B, A, &t, &u, &v, N))
			{
				norm = glm::normalize(N);
				if (rayIntersectsTriangle(posNext, norm, C, B, A, &t, &u, &v, N))
				{
					t += 0.0003;
					contactP = posNext + t * norm;
					return true;
				}
				else
				{
					glm::vec3 A1(cellIndexNext[0], GetHeightfieldAt(cellIndexNext[0], cellIndexNext[1]), cellIndexNext[1]);
					glm::vec3 AA1(cellIndexNext[0] + 1, GetHeightfieldAt(cellIndexNext[0] + 1, cellIndexNext[1] + 1), cellIndexNext[1] + 1);
					glm::vec3 B1(cellIndexNext[0] + 1, GetHeightfieldAt(cellIndexNext[0] + 1, cellIndexNext[1]), cellIndexNext[1]);
					glm::vec3 C1(cellIndexNext[0], GetHeightfieldAt(cellIndexNext[0], cellIndexNext[1] + 1), cellIndexNext[1] + 1);

					//1st triangle
					if (rayIntersectsTriangle(posNext, norm, C1, B1, A1, &t, &u, &v, N))
					{
						t += 0.0003;
						contactP = posNext + t * norm;
						return true;
					}
					//2nd triangle
					else if (rayIntersectsTriangle(posNext, norm, AA1, B1, C1, &t, &u, &v, N))
					{
						t += 0.0003;
						contactP = posNext + t * norm;
						return true;
					}
				}
				return false;
			}
			else if (rayIntersectsTriangle(posCurr, dir, AA, B, C, &t, &u, &v, N))
			{
				norm = glm::normalize(N);
				if (rayIntersectsTriangle(posNext, norm, AA, B, C, &t, &u, &v, N))
				{
					t += 0.0003;
					contactP = posNext + t * norm;
					return true;
				}
				else
				{
					glm::vec3 A1(cellIndexNext[0], GetHeightfieldAt(cellIndexNext[0], cellIndexNext[1]), cellIndexNext[1]);
					glm::vec3 AA1(cellIndexNext[0] + 1, GetHeightfieldAt(cellIndexNext[0] + 1, cellIndexNext[1] + 1), cellIndexNext[1] + 1);
					glm::vec3 B1(cellIndexNext[0] + 1, GetHeightfieldAt(cellIndexNext[0] + 1, cellIndexNext[1]), cellIndexNext[1]);
					glm::vec3 C1(cellIndexNext[0], GetHeightfieldAt(cellIndexNext[0], cellIndexNext[1] + 1), cellIndexNext[1] + 1);

					//1st triangle
					if (rayIntersectsTriangle(posNext, norm, C1, B1, A1, &t, &u, &v, N))
					{
						t += 0.0003;
						contactP = posNext + t * norm;
						return true;
					}
					//2nd triangle
					else if (rayIntersectsTriangle(posNext, norm, AA1, B1, C1, &t, &u, &v, N))
					{
						t += 0.0003;
						contactP = posNext + t * norm;
						return true;
					}
				}
				return false;
			}
		}

		return false;
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

