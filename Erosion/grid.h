#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <omp.h>

struct Voxel
{
	float density;

	Voxel()
	{
		density = 0;
	}

	Voxel(float d)
	{
		density = d;
	}
};


class Grid
{
private:
	int m_DimX;
	int m_DimY;
	int m_DimZ;
	glm::vec3 m_Dim;

	bool* m_Grid;
	std::vector<float> coords;

	struct Heightfield
	{
		int dimX = 0;
		int dimY = 0;

		unsigned char* map = nullptr;
	};

	Heightfield m_heightfield;

	void destroyGrid()
	{
		delete[] m_Grid;
		coords.clear();
	}

	void createGrid(int dimx, int dimy, int dimz)
	{
		m_Grid = new bool[dimx * dimy * dimz];
		m_Dim = glm::vec3(dimx, dimy, dimz);
		m_DimX = m_Dim.x;
		m_DimY = m_Dim.y;
		m_DimZ = m_Dim.z;
	}

public:
	Grid(int dimX = 512, int dimY = 512, int dimZ = 512)
		: m_DimX(dimX), m_DimY(dimY), m_DimZ(dimZ)
	{
		m_Dim = glm::vec3(dimX, dimY, dimZ);
		m_Grid = new bool[dimX * dimY * dimZ];

		m_heightfield.dimX = dimX;
		m_heightfield.dimY = dimY;
		m_heightfield.map = new unsigned char[dimX * dimZ];
		memset(m_heightfield.map, 0, dimX * dimZ);
		//for (int i = 0; i < dimX * dimZ; i++)
		//	m_heightfield.map[i] = 0;
	}

	Voxel GetVoxel(int x, int y, int z) const
	{
		return m_Grid[x + m_DimX * (y + m_DimY * z)];
	}

	void SetVoxel(int x, int y, int z, bool density)
	{
		m_Grid[x + m_DimX * (y + m_DimY * z)] = density;
	}

	void LoadHeightfield(unsigned char* img)
	{
		memcpy(m_heightfield.map, img, m_heightfield.dimX * m_heightfield.dimY);
	}

	unsigned char GetHeightfieldAt(int x, int y)
	{
		return m_heightfield.map[m_DimY * x + y];
	}

	void UpdateGrid(int dimx, int dimy, int dimz)
	{
		destroyGrid();
		createGrid(dimx, dimy, dimz);
		for (int y = 0; y < m_Dim.y; y++)
			for (int z = 0; z < m_Dim.z; z++)
				for (int x = 0; x < m_Dim.x; x++)
					if (y <= GetHeightfieldAt(x, z))
					{
						SetVoxel(x, y, z, 1);
						coords.push_back((float)x / 100.0f + 100);
						coords.push_back((float)y / 100.0f + 100);
						coords.push_back((float)z / 100.0f + 100);
					}
					else
						SetVoxel(x, y, z, 0);
	}

	glm::vec3 GetDim() const
	{
		return m_Dim;
	}

	float*  GetCoords()
	{
		return coords.data();
	}

	size_t GetCoordsSize() const
	{
		return coords.size();
	}
};

