#pragma once
#include "voxel.h"
#include <glm/glm.hpp>


class FluidSystem
{
public:
	FluidSystem() = default;

	FluidSystem(glm::vec3 o, glm::vec3 d)
		: m_Origin(o), m_Dim(d)
	{
		m_Volume = new Voxel[d.x * d.y * d.z];
		for(int y = 0; y < m_Dim.y; y++)
			for(int z = 0; z < m_Dim.z; z++)
				for (int x = 0; x < m_Dim.x; x++)
				{
					m_Volume->position = glm::vec3(x, y, z);
					m_Volume->type = VoxelType::VOXEL_WAT;
					m_Volume->velocity = glm::vec3(1.0);
				}
	}

	~FluidSystem()
	{
		delete[] m_Volume;
	}
	
	void operator=(const FluidSystem& other)
	{
		m_Dim = other.m_Dim;
		m_Origin = other.m_Origin;
		m_Volume = new Voxel[m_Dim.x * m_Dim.y * m_Dim.z];
		for (int i = 0; i < m_Dim.x * m_Dim.y * m_Dim.z; i++)
			m_Volume[i] = other.m_Volume[i];

	}

	Voxel* GetVolume() const
	{
		return m_Volume;
	}

	glm::vec3 m_Origin;

private:
	glm::vec3 m_Dim;
	Voxel *m_Volume;
};