#pragma once
#include "voxel.h"
#include <glm/glm.hpp>

#include "mesh.h"
#include "sphere.h"


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
					m_Volume[x + (int)m_Dim.x * (y + (int)m_Dim.y * z)].velocity = glm::vec3(0.2, -GRAVITY, 0.0);
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
#define PI 3.141592f

#define FPOS 0
#define FVEL 1
#define FDENSITY 2
#define FMASS 3
#define FPRESS 4

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


struct FluidParticle
{
	glm::vec3 Position;
	glm::vec3 Velocity;
	glm::vec3 Acceleration;
	float Mass;
	float Density;
	float Pressure;
	glm::vec3 PressureForce;
	glm::vec3 ViscosityForce;
	glm::vec3 GravityForce;
	glm::vec3 SurfaceForce;
};

class FluidSystemSPH
{
public:
	FluidSystemSPH()
	{

	}

	void Initialize(int nParts)
	{
		num = nParts;

		for (int i = 0; i < cbrt(num); i++)
			for (int j = 0; j < cbrt(num); j++)
				for (int k = 0; k < cbrt(num); k++)
				{
					float x = -0.2 + i * 0.025;
					float y = -0.05 + j * 0.025;
					float z = -0.15 + k * 0.025;

					FluidParticle particle;
					particle.Position = glm::vec3(x + 1, y + 256, z + 1);
					particle.Velocity = glm::vec3(0.0);
					particle.Acceleration = glm::vec3(0.0);
					particle.Mass = m;
					m_Particles.push_back(particle);
				}

		vol = num * m / p0;
		smoothRadius = cbrt(3 * vol * x / (4 * PI * num));
		h = smoothRadius;
		l = sqrt(p0 / x);
	}

	void Run(Mesh terrain)
	{

		//compute density and pressure
#pragma omp parallel for collapse(2)
		for (int i = 0; i < num; i++)
		{
			FluidParticle& currPart = m_Particles[i];
			
			float density = 0;
			for (int j = 0; j < num; j++)
			{
				FluidParticle& neighbPart = m_Particles[j];

				float distance = glm::length(currPart.Position - neighbPart.Position);
				if (distance <= smoothRadius && distance != 0)
					density += neighbPart.Mass * kernDefault(currPart.Position - neighbPart.Position);
			}
			currPart.Density = density;
			currPart.Pressure = k * (currPart.Density - p0);
		}


		//compute forces
#pragma omp parallel for collapse(3)
		for (int i = 0; i < num; i++)
		{
			FluidParticle& currPart = m_Particles[i];

			//internal forces
			glm::vec3 fPress(0.0);
			glm::vec3 fVisc(0.0);
			for (int j = 0; j < num; j++)
			{
				FluidParticle& neighbPart = m_Particles[j];

				float distance = glm::length(currPart.Position - neighbPart.Position);
				if (distance <= smoothRadius && distance != 0)
				{
					fPress += (currPart.Pressure / (float)pow(currPart.Density, 2) + neighbPart.Pressure / (float)pow(neighbPart.Density, 2)) * neighbPart.Mass * gradPressure(currPart.Position - neighbPart.Position);
					fVisc += (neighbPart.Velocity - currPart.Velocity) * neighbPart.Mass / neighbPart.Density * laplVisc(currPart.Position - neighbPart.Position);
				}

			}
			currPart.PressureForce = -currPart.Density * fPress;
			currPart.ViscosityForce = visc * fVisc;

			//external forces
			currPart.GravityForce = currPart.Density * g;
			glm::vec3 n(0.0); //inward surface normal
			float colorFieldLapl = 0.0;
			for (int j = 0; j < num; j++)
			{
				FluidParticle& neighbPart = m_Particles[j];

				n += neighbPart.Mass / neighbPart.Density * gradDefault(currPart.Position - neighbPart.Position);
				colorFieldLapl += neighbPart.Mass / neighbPart.Density * laplDefault(currPart.Position - neighbPart.Position);
			}

			currPart.SurfaceForce = glm::vec3(0.0);
			if (glm::length(n) >= l)
				currPart.SurfaceForce = -surf_tens * colorFieldLapl * (n / glm::length(n));
		}

		advance(terrain);
		m_Time += dt;
	}

	void Draw()
	{
		for (int i = 0; i < num; i++)
		{
			float r = s * cbrt(3 * m_Particles[i].Mass / (4 * PI * m_Particles[i].Density));
			m_Spheres.emplace_back(30, 30, r, m_Particles[i].Position);
			m_Spheres[i].Draw();
		}
		m_Spheres.clear();
	}

private:
	//Simulation parameters
	//float m_Param[MAX_PARAM];
	//FBufs m_Fluid;
	glm::vec3 g = glm::vec3(0.0, -9.82, 0.0);
	const float dt = 0.01;
	float m_Time = 0.0;
	const float p0 = 998.29;
	const float m = 0.02;
	const float visc = 3.5;
	const float surf_tens = 0.0728;
	//const float l = 7.065;
	float l;
	const float k = 3.0;
	const float cR = 0;
	const int x = 20;
	//const float h = 0.0457;
	float h;
	int num = 10;
	float vol;
	float fluidDensity = 1000.0;
	float smoothRadius;
	float s = 1;

private:
	std::vector<FluidParticle> m_Particles;
	std::vector<Sphere> m_Spheres;

private:

	std::vector<float> getVertices()
	{
		std::vector<float> output;
		for (const auto& particle : m_Particles)
		{
			output.push_back(particle.Position.x);
			output.push_back(particle.Position.y);
			output.push_back(particle.Position.z);
		}
		return output;
	}

	//leap-frog + collision handling
	void advance(Mesh terrain)
	{
#pragma omp parallel for
		for (int i = 0; i < num; i++)
		{
			FluidParticle& currPart = m_Particles[i];
			glm::vec3 F;
			glm::vec3 fInternal;
			glm::vec3 fExternal;
			glm::vec3 velNext;
			glm::vec3 posNext;

			fInternal = currPart.PressureForce + currPart.ViscosityForce;
			fExternal = currPart.GravityForce + currPart.SurfaceForce;
			F = fInternal + fExternal;
			currPart.Acceleration = F / currPart.Density;

			//if initial velocity offset not initialized
			if (m_Time == 0.0)
			{
				currPart.Velocity = currPart.Velocity - 0.5f * dt * currPart.Acceleration;
			}

			velNext = currPart.Velocity + dt * currPart.Acceleration;
			posNext = currPart.Position + dt * velNext;
			handleCollision(currPart.Position, posNext, velNext, terrain);

			currPart.Velocity = velNext;
			currPart.Position = posNext;
		}
	}

	
	void handleCollision(glm::vec3 & posCurr,  const glm::vec3& posNext, glm::vec3 & velNext, const Mesh terrain)
	{
		float t, u, v; //where t is R.Origin + t*R.Dir;  u,v - barycentric coords
		glm::vec3 N; //normal
		glm::vec3 dir = glm::normalize(posNext - posCurr);

		//check collision with each triangle
		for (int i = 0; i < terrain.Indices.size(); i += 3)
		{
			float a1 = terrain.Vertices[terrain.Indices[i] * 3]; //times 3, bc 1 vertex is composed out of 3 floats
			float a2 = terrain.Vertices[terrain.Indices[i] * 3 + 1];
			float a3 = terrain.Vertices[terrain.Indices[i] * 3 + 2];
			float b1 = terrain.Vertices[terrain.Indices[i+1] * 3];
			float b2 = terrain.Vertices[terrain.Indices[i+1] * 3 + 1];
			float b3 = terrain.Vertices[terrain.Indices[i+1] * 3 + 2];
			float c1 = terrain.Vertices[terrain.Indices[i+2] * 3];
			float c2 = terrain.Vertices[terrain.Indices[i+2] * 3 + 1];
			float c3 = terrain.Vertices[terrain.Indices[i+2] * 3 + 2];

			//triangle
			glm::vec3 A(a1, a2, a3);
			glm::vec3 B(b1, b2, b3);
			glm::vec3 C(c1, c2, c3);

			if (rayIntersectsTriangle(posCurr, dir, A, B, C, &t, &u, &v, N))
			{
				glm::vec3 contactP = posCurr + t * dir;
				float distanceToTravel = glm::length(posNext - posCurr);
				float distanceToIntersect = glm::length(contactP - posCurr);
				if (distanceToTravel > distanceToIntersect)
				{
					float d = glm::length(posNext - contactP); //penetration depth 
					posCurr = contactP;
					velNext = velNext - (1 + cR * (d / (dt * glm::length(velNext)))) * (velNext * N) * N;
					return;
				}
			}
		}
	}

	bool rayIntersectsTriangle(const glm::vec3 & pos, const glm::vec3 & dir, const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, float* t, float* u, float* v, glm::vec3& N)
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

	float kernDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		if (len >= 0 && len <= h)
			return (315.0 / (64.0 * PI * pow(h, 9))) * pow((h*h - len*len), 3);
		else
			return 0;
	}

	glm::vec3 gradDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		return -945.0f / (32.0f * PI * (float)pow(h, 9)) * r * (float)pow((h * h - len * len), 2);
	}

	float laplDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		return (-945.0f / 32.0f * PI * (float)pow(h, 9)) * (h*h - len*len) * (3 * h*h - 7 * len*len);
	}

	glm::vec3 gradPressure(const glm::vec3& r)
	{

		float len = glm::length(r);
		return (float)(-45.0f / PI * (float)pow(h, 6)) * (r / len) * (float)pow((h - len), 2);

	}

	float laplVisc(const glm::vec3& r)
	{
		float len = glm::length(r);
		return (45.0f / PI * (float)pow(h, 6)) * (h - len);
	}
};