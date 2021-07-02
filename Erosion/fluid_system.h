#pragma once
#include <memory>
#include <glm/gtc/matrix_transform.hpp>

#include "voxel.h"
#include "mesh.h"
#include "sphere.h"
#include "shader.h"


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
	glm::vec3 SurfaceNormal;
};

class FluidSystemSPH
{
public:
	FluidSystemSPH()
	{
		m_Origin = glm::vec3(0.0, 0.0, 0.0);
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
					particle.Position = glm::vec3(x + m_Origin.x, y + m_Origin.y, z + m_Origin.z);
					particle.Velocity = glm::vec3(0.0, 0.0, 0.0);
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
				if (distance <= smoothRadius)
					density += neighbPart.Mass * kernDefault(currPart.Position - neighbPart.Position);
			}
			currPart.Density = density;
			currPart.Pressure = k * (currPart.Density - p0);
		}


		//internal forces
#pragma omp parallel for collapse(2)
		for (int i = 0; i < num; i++)
		{
			FluidParticle& cur = m_Particles[i];

			glm::vec3 fPress(0.0);
			glm::vec3 fVisc(0.0);
			glm::vec3 n(0.0); //inward surface normal
			for (int j = 0; j < num; j++)
			{
				FluidParticle& neighb = m_Particles[j];
				float distance = glm::length(cur.Position - neighb.Position);
				if (distance <= smoothRadius && i != j)
				{
					fPress += (cur.Pressure / (cur.Density * cur.Density) + neighb.Pressure / (neighb.Density * neighb.Density)) * neighb.Mass * gradPressure(cur.Position - neighb.Position);
					fVisc += (neighb.Velocity - cur.Velocity) * (neighb.Mass / neighb.Density) * laplVisc(cur.Position - neighb.Position);
					n += (neighb.Mass / neighb.Density) * gradDefault(cur.Position - neighb.Position);
				}

			}
			cur.PressureForce = (-1) * cur.Density * fPress;
			cur.ViscosityForce = visc * fVisc;
			cur.SurfaceNormal = n;
		}

		//external forces
#pragma omp parallel for collapse(2)
		for (int i = 0; i < num; i++)
		{
			FluidParticle& currPart = m_Particles[i];
			currPart.GravityForce = currPart.Density * g;
			float colorFieldLapl = 0.0;
			for (int j = 0; j < num; j++)
			{
				FluidParticle& neighbPart = m_Particles[j];

				float distance = glm::length(currPart.Position - neighbPart.Position);
				if (distance <= smoothRadius)
				{
					colorFieldLapl += (neighbPart.Mass / neighbPart.Density) * laplDefault(currPart.Position - neighbPart.Position);
				}
			}

			//currPart.SurfaceForce = glm::vec3(0.0);
			//if (glm::length(n) >= l)
			//float ln = glm::length(n);
			currPart.SurfaceForce = -surf_tens * colorFieldLapl * currPart.SurfaceNormal;
		}

		advance(terrain);
		m_Time += m_Dt;
	}

	void Draw(const Shader& shader)
	{
		if (!m_Sphere)
			m_Sphere = std::make_unique<Sphere>(10, 10, 1, glm::vec3(0.0, 0.0, 0.0));

		for (int i = 0; i < num; i++)
		{
			float r = s * cbrt(3 * m_Particles[i].Mass / (4 * PI * m_Particles[i].Density));
			// r = 0.01f;
			glm::mat4 model = glm::mat4(1.0);
			model = glm::translate(model, m_Particles[i].Position);
			model = glm::scale(model, glm::vec3(r));
			shader.setMat4("model", model);
			m_Sphere->Draw();
		}
	}

	void SetOrigin(const glm::vec3& o )
	{
		m_Origin = o;
	}

	glm::vec3 GetOrigin() const
	{
		return m_Origin;
	}

	void SetDeltaTime(float dt)
	{
		m_Dt = dt;
	}

	float GetDeltaTime() const
	{
		return m_Dt;
	}

	void PrintCoords() const
	{
		for (int i = 0; i < num; i++)
			std::cout << "[" << i << "] " << m_Particles[i].Position.x << ' ' << m_Particles[i].Position.y << ' ' << m_Particles[i].Position.z << std::endl;
	}

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
//#pragma omp parallel for
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
			glm::vec3 acc = F / currPart.Density;

			//if initial velocity offset not initialized
			//if (m_Time == 0.0)
			//{
			//	currPart.Velocity = currPart.Velocity - 0.5f * m_Dt * currPart.Acceleration;
			//}

			velNext = currPart.Velocity + m_Dt * acc;
			posNext = currPart.Position + m_Dt * velNext;
			handleCollision(currPart.Position, posNext, velNext, terrain);

			currPart.Velocity = velNext;
			currPart.Position = posNext;
			
		}
	}

	
	void handleCollision(const glm::vec3 & posCurr,  glm::vec3& posNext, glm::vec3 & velNext, const Mesh terrain)
	{
		float t, u, v; //where t is R.Origin + t*R.Dir;  u,v - barycentric coords
		glm::vec3 N; //normal
		//glm::vec3 dir = glm::normalize(posNext - posCurr);
		glm::vec3 dir = glm::normalize(velNext); 
		std::cout << dir.x << " " << dir.y << " " << dir.z << std::endl;


		//check collision with each triangle
		for (int i = 0; i < terrain.Indices.size(); i += 3)
		{
			float a1 = terrain.Vertices[terrain.Indices[i] * 6]; // 6 - offset
			float a2 = terrain.Vertices[terrain.Indices[i] * 6 + 1];
			float a3 = terrain.Vertices[terrain.Indices[i] * 6 + 2];
			float b1 = terrain.Vertices[terrain.Indices[i+1] * 6];
			float b2 = terrain.Vertices[terrain.Indices[i+1] * 6 + 1];
			float b3 = terrain.Vertices[terrain.Indices[i+1] * 6 + 2];
			float c1 = terrain.Vertices[terrain.Indices[i+2] * 6];
			float c2 = terrain.Vertices[terrain.Indices[i+2] * 6 + 1];
			float c3 = terrain.Vertices[terrain.Indices[i+2] * 6 + 2];

			//triangle
			glm::vec3 A(a1, a2, a3);
			glm::vec3 B(b1, b2, b3);
			glm::vec3 C(c1, c2, c3);

			if (rayIntersectsTriangle(posCurr, dir, C, B, A, &t, &u, &v, N))
			{
				glm::vec3 dirN = N;
				//std::cout << "before int: " << N.x << " " << N.y << " " << N.z << std::endl;
				if (rayIntersectsTriangle(posNext, N, C, B, A, &t, &u, &v, N))
				{
					//std::cout << "Intersec with itself: " << N.x << " " << N.y << " " << N.z << std::endl;
					t += 0.00009; 
					glm::vec3 contactP = posNext + t * N;
					float d = glm::length(contactP - posNext);
					velNext = velNext - (float)(1 + cR * (d / (m_Dt * glm::length(velNext)))) * glm::dot(velNext, N) * N;
					posNext = contactP;
					return;
				}
				else
				{
					//iterate over all triangles again 
					for (int i = 0; i < terrain.Indices.size(); i += 3)
					{
						float a1 = terrain.Vertices[terrain.Indices[i] * 6]; // 6 - offset
						float a2 = terrain.Vertices[terrain.Indices[i] * 6 + 1];
						float a3 = terrain.Vertices[terrain.Indices[i] * 6 + 2];
						float b1 = terrain.Vertices[terrain.Indices[i + 1] * 6];
						float b2 = terrain.Vertices[terrain.Indices[i + 1] * 6 + 1];
						float b3 = terrain.Vertices[terrain.Indices[i + 1] * 6 + 2];
						float c1 = terrain.Vertices[terrain.Indices[i + 2] * 6];
						float c2 = terrain.Vertices[terrain.Indices[i + 2] * 6 + 1];
						float c3 = terrain.Vertices[terrain.Indices[i + 2] * 6 + 2];

						//triangle
						glm::vec3 A(a1, a2, a3);
						glm::vec3 B(b1, b2, b3);
						glm::vec3 C(c1, c2, c3);
						if (rayIntersectsTriangle(posNext, dirN, C, B, A, &t, &u, &v, N))
						{
							//std::cout << "Intersec with adj: " << N.x << " " << N.y << " " << N.z << std::endl;
							t += 10.2;
							glm::vec3 contactP = posNext + t * dirN;
							float d = glm::length(contactP - posNext);
							velNext = velNext - (float)(1 + cR * d / (m_Dt * glm::length(velNext))) * glm::dot(velNext, dirN) * dirN;
							posNext = contactP;
							return;
						}
					}
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
		float invdet = 1.0 / det ;
		glm::vec3 AO = pos - A;
		glm::vec3 DAO = glm::cross(AO, dir);
		*u = glm::dot(E2, DAO) * invdet;
		*v = -glm::dot(E1, DAO) * invdet;
		*t = glm::dot(AO, N) * invdet;
		return (abs(det) >= 1e-6 && *t >= 0.0 && *u >= 0.0 && *v >= 0.0 && (*u + *v) <= 1.0);
	}

	bool RayIntersectsTriangle(glm::vec3 rayOrigin, glm::vec3 rayVector, glm::vec3 & A, glm::vec3& B, glm::vec3& C, glm::vec3& outIntersectionPoint, glm::vec3 & N)
	{
		const float EPSILON = 0.0000001;
		glm::vec3 vertex0 = A;
		glm::vec3 vertex1 = B;
		glm::vec3 vertex2 = C;
		glm::vec3 edge1, edge2, h, s, q;
		float a, f, u, v;
		edge1 = vertex1 - vertex0;
		edge2 = vertex2 - vertex0;
		N = glm::cross(edge1, edge2);
		h = glm::cross(rayVector,edge2);
		a = glm::dot(edge1, h);
		if (a > -EPSILON && a < EPSILON)
			return false;    // This ray is parallel to this triangle.
		f = 1.0 / a;
		s = rayOrigin - vertex0;
		u = f * glm::dot(s, h);
		if (u < 0.0 || u > 1.0)
			return false;
		q = glm::cross(s, edge1);
		v = f * glm::dot(rayVector, q);
		if (v < 0.0 || u + v > 1.0)
			return false;
		// At this stage we can compute t to find out where the intersection point is on the line.
		float t = f * glm::dot(edge2, q);
		if (t > EPSILON) // ray intersection
		{
			outIntersectionPoint = rayOrigin + rayVector * t;
			return true;
		}
		else // This means that there is a line intersection but not a ray intersection.
			return false;
	}

	float kernDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		if (len >= 0 && len <= h)
			return (float)(315.0f / (64.0f * PI * powf(h, 9.0f))) * powf((h*h - len*len), 3.0f);
		else
			return 0;
	}

	glm::vec3 gradDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		return -945.0f / (32.0f * PI * (float)powf(h, 9.0f)) * r * (float)powf((h * h - len * len), 2.0f);
	}

	float laplDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		return (-945.0f / (32.0f * PI * powf(h, 9.0f))) * (h*h - len*len) * (3 * h*h - 7 * len*len);
	}

	glm::vec3 gradPressure(const glm::vec3& r)
	{

		//float len = glm::length(r);
		//return (float)(-45.0f / (PI * (float)pow(h, 6))) * (r / len) * (float)pow((h - len), 2);


		float dist = glm::length(r);
		if (dist < 10e-5) {
			return -(glm::normalize(glm::vec3(1.0f, 1.0f, 1.0f)) * (float)(45.f / (PI * powf(h, 6.0f))) * powf(h - dist, 2.0f));
		}
		else {
			glm::vec3 returnValue = -(glm::normalize(r) * (float)(45.f / (PI * powf(h, 6.0f))) * powf(h - dist, 2.0f));

			return returnValue;
		}

	}

	float laplVisc(const glm::vec3& r)
	{
		float len = glm::length(r);
		return (float)(45.0f / (PI * powf(h, 6.0f))) * (h - len);
	}

	

private:
		//Simulation parameters
		//float m_Param[MAX_PARAM];
		//FBufs m_Fluid;
		glm::vec3 g = glm::vec3(0.0, -9.82f, 0.0);
		float m_Dt = 0.0f;
		float m_Time = 0.0f;
		const float p0 = 998.29f;
		const float m = 0.02f;
		const float visc = 3.5f;
		const float surf_tens = 0.0728f;
		//const float l = 7.065;
		float l;
		const float k = 3.0f;
		const float cR = 0.1f;
		const int x = 20;
		//const float h = 0.0457;
		float h;
		int num = 10;
		float vol;
		float smoothRadius;
		float s = 1;

private:
	std::vector<FluidParticle> m_Particles;
	std::unique_ptr<Sphere> m_Sphere;
	glm::vec3 m_Origin;
};