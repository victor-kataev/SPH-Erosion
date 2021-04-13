#pragma once
#include "voxel.h"
#include <glm/glm.hpp>


#define GRAVITY 9.8

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
					m_Volume[x + (int)m_Dim.x * (y + (int)m_Dim.y * z)].position = glm::vec3(x, y, z) + m_Origin;
					m_Volume[x + (int)m_Dim.x * (y + (int)m_Dim.y * z)].type = VoxelType::VOXEL_WAT;
					m_Volume[x + (int)m_Dim.x * (y + (int)m_Dim.y * z)].velocity = glm::vec3(0.0, -GRAVITY, 0.0);
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

	void Run()
	{
		glm::vec3 gravity(0.0, -1 * GRAVITY, 0.0);
		glm::vec3 delta;
		for (int y = 0; y < m_Dim.y; y++)
			for (int z = 0; z < m_Dim.z; z++)
				for (int x = 0; x < m_Dim.x; x++)
				{
					size_t idx = x + (int)m_Dim.x * (y + (int)m_Dim.y * z);
					delta = 0.01f * (gravity + m_Volume[idx].velocity);
					m_Volume[idx].position += delta;
				}
	}


	glm::vec3 m_Origin;
	glm::vec3 m_Dim;
	Voxel *m_Volume;

	static float Cr;

};


float FluidSystem::Cr = 0.3;

typedef unsigned int uint;


#define MAX_PARAM  50
#define MAX_BUF	   25

struct FBufs
{
	inline ~FBufs() {
		for (int i = 0; i < MAX_BUF; i++) {
			if (mcpu[i] != nullptr) {
				free(mcpu[i]);
			}
		}
	}


	//inline Vector3DF* bufVector3DF(int n) { return (Vector3DF*)mcpu[n]; }
	//inline Float3* bufFloat3(int n) { return (Float3*)mcpu[n]; }
	inline float* bufF(int n) { return (float*)mcpu[n]; }
	inline uint* bufI(int n) { return (uint*)mcpu[n]; }
	inline char* bufC(int n) { return (char*)mcpu[n]; }

	char* mcpu[MAX_BUF] = { nullptr };
};


class FluidSystemSPH
{
private:
	//Simulation parameters
	float m_Param[MAX_PARAM];
	FBufs m_Fluid;
};