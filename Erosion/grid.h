#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <list>
#include <omp.h>
#include <cstring>

#include "fluid_system.h"

class Graph
{
public:
	Graph(int dimx, int dimz)
		: m_DimX(dimx), m_DimZ(dimz)
	{
		m_adjacents = new std::vector<std::pair<int, int>>*[dimz];
		for(int i = 0; i < dimx; i++)
			m_adjacents[i] = new std::vector<std::pair<int, int>>[dimx];
	}

	~Graph()
	{
		for (int i = 0; i < m_DimZ; i++)
			delete[] m_adjacents[i];
		delete[] m_adjacents;
	}

	void AddEdge(int x1, int z1, int x2, int z2)
	{
		if (x2 >= m_DimX || x2 < 0 || z2 >= m_DimZ || z2 < 0)
			return;

		std::pair<int, int> coords2 = std::make_pair<int, int>((int)x2, (int)z2);
		auto found = std::find(m_adjacents[x1][z1].begin(), m_adjacents[x1][z1].end(), coords2);
		if (found != m_adjacents[x1][z1].end())
			return;
		std::pair<int, int> coords1 = std::make_pair<int, int>((int)x1, (int)z1);
		m_adjacents[x1][z1].push_back(coords2);
		m_adjacents[x2][z2].push_back(coords1);
	}

	std::vector<std::pair<int, int>> GetAdjacents(int x, int z) const
	{
		return m_adjacents[x][z];
	}

private:
	std::vector<std::pair<int, int>> **m_adjacents;
	int m_DimX;
	int m_DimZ;

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
	FluidSystem m_Fluid;

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

	void LoadFluid(const FluidSystem& fs)
	{
		m_Fluid = fs;
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

	//void buildGraph()
	//{
	//	
	//	for(int z = 0; z < m_Dim.z; z++)
	//		for (int x = 0; x < m_Dim.x; x++)
	//		{
	//			m_Graph.AddEdge(x, z, x + 1, z);
	//			m_Graph.AddEdge(x, z, x, z + 1);
	//			m_Graph.AddEdge(x, z, x - 1, z);
	//			m_Graph.AddEdge(x, z, x, z - 1);
	//		}
	//}

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

	//void generateIndices()
	//{
	//	bool **visited = new bool*[m_DimZ];
	//	for (int i = 0; i < m_DimZ; i++)
	//	{
	//		visited[i] = new bool[m_DimX];
	//		for (int j = 0; j < m_DimX; j++)
	//			visited[i][j] = false;
	//	}

	//	std::list<std::pair<int, int>> queue;

	//	visited[0][0] = true;
	//	queue.push_back(std::pair<int, int>(0, 0));
	//	std::list<unsigned int> triangle;

	//	while (!queue.empty())
	//	{
	//		std::pair<int, int> current = queue.front();
	//		queue.pop_front();
	//		triangle.push_back(current.second * m_DimX + current.first);
	//		
	//		if (triangle.size() == 3)
	//		{
	//			for (unsigned int index : triangle)
	//				indices.push_back(index);

	//			triangle.pop_front();
	//		}

	//		for (const auto& adj : m_Graph.GetAdjacents(current.first, current.second))
	//		{
	//			if (!visited[adj.first][adj.second])
	//			{
	//				queue.push_back(adj);
	//				visited[adj.first][adj.second] = true;
	//			}
	//		}
	//	}
	//}

	void UpdateGridDepricated(int dimx, int dimy, int dimz)
	{
		destroyGrid();
		createGrid(dimx, dimy, dimz);
		for (int z = 0; z < m_Dim.z; z++)
			for (int x = 0; x < m_Dim.x; x++)
				for (int y = 0; y < m_Dim.y; y++)
					if (y <= GetHeightfieldAt(x, z))
					{
						SetVoxel(x, y, z, VoxelType::VOXEL_MAT);
						//SetVoxel(x, y, z, 1);
						m_filled_voxels++;

						if (y == GetHeightfieldAt(x, z) || y == m_Dim.y - 1)
						{
							surfaceParts.push_back((float)x / 100.0f);
							surfaceParts.push_back((float)y / 100.0f);
							surfaceParts.push_back((float)z / 100.0f);
						}
					}
					else
					{
						SetVoxel(x, y, z, VoxelType::VOXEL_AIR);
						//SetVoxel(x, y, z, 0);
						m_empty_voxels++;
					}
		genIndices();
		std::cout << "empty voxels = " << m_empty_voxels << " --- filled voxels = " << m_filled_voxels << std::endl;
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
				surfaceParts.push_back((float)x / 100.0f);
				surfaceParts.push_back((float)y / 100.0f);
				surfaceParts.push_back((float)z / 100.0f);
			}
		genIndices();
		for (int y = 0; y < m_Fluid.m_Dim.y; y++)
			for(int z = 0; z < m_Fluid.m_Dim.z; z++)
				for(int x = 0; x < m_Fluid.m_Dim.x; x++)
				{
					Voxel fluid_voxel = m_Fluid.m_Volume[x + (int)m_Fluid.m_Dim.x * (y + (int)m_Fluid.m_Dim.y * z)];
					if (fluid_voxel.position.x >= 0 && fluid_voxel.position.x < m_Dim.x &&
						fluid_voxel.position.y >= 0 && fluid_voxel.position.y < m_Dim.y &&
						fluid_voxel.position.z >= 0 && fluid_voxel.position.z < m_Dim.z)
					{
						fluidParts.push_back((float)fluid_voxel.position.x / 100.0f);
						fluidParts.push_back((float)fluid_voxel.position.y / 100.0f);
						fluidParts.push_back((float)fluid_voxel.position.z / 100.0f);
					}
				}
	}

	unsigned int* GetIndices()
	{
		return indices.data();
	}

	size_t GetIndicesSize() const
	{
		return indices.size();
	}

	glm::vec3 GetDim() const
	{
		return m_Dim;
	}

	float* GetSurfaceParts()
	{
		return surfaceParts.data();
	}

	size_t GetSurfacePartsSize() const
	{
		return surfaceParts.size();
	}

	float* GetFluidParts()
	{
		return fluidParts.data();
	}

	size_t GetFluidPartsSize() const
	{
		return fluidParts.size();
	}
};

