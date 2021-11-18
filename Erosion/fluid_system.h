#pragma once
#include <memory>
#include <glm/gtc/matrix_transform.hpp>

#include "mesh.h"
#include "sphere.h"
#include "grid.h"


#define GRAVITY 9.8



typedef unsigned int uint;


#define MAX_PARAM  50
#define MAX_BUF	   25
#define M_PI 3.14159265358979323846

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
		init_num = nParts;
		id = 0;

		for (int i = 0; i < cbrt(num); i++)
			for (int j = 0; j < cbrt(num); j++)
				for (int k = 0; k < cbrt(num); k++)
				{
					float x = -0.2 + i * 0.025;
					float y = -0.05 + j * 0.025;
					float z = -0.15 + k * 0.025;

					FluidParticle particle;
					particle.Id = id++;
					particle.Position = glm::vec3(x + m_Origin.x, y + m_Origin.y, z + m_Origin.z);
					particle.Velocity = glm::vec3(0.0, 0.0, 0.0);
					particle.Acceleration = glm::vec3(0.0);
					//particle.Mass = MASS;
					m_Particles.push_back(particle);
				}

		vol = num * MASS / p0;
		//smoothRadius = cbrt(3 * vol * x / (4 * PI * num));
		smoothRadius = h;
		deltaS = h / 2.0f;
		//h = smoothRadius;
		//l = sqrt(p0 / x);
	}

	void Run(Grid& grid)
	{

		//compute density and pressure
#pragma omp parallel for collapse(2)
		for (int i = 0; i < m_Particles.size(); i++)
		{
			FluidParticle& currPart = m_Particles[i];

			float density = 0;
			int cnt = 0;
			for (int j = 0; j < m_Particles.size(); j++)
			{
				FluidParticle& neighbPart = m_Particles[j];

				float distance = glm::length(currPart.Position - neighbPart.Position);
				if (distance <= smoothRadius)
				{
					cnt++;
					density += MASS * kernDefault(currPart.Position - neighbPart.Position);
				}
			}
			if (density < 1000 && cnt == 1)
				density = 328.0f;
			else if (density < 1000 && cnt == 2)
				density = 1005.0f;
			else if (density < 1000 && cnt == 3)
				density = 1020.0f;
			else if(density < 1000)
				density = 1002.0f;
			
			currPart.Density = density;
			currPart.Pressure = k * (currPart.Density - p0);
		}


		//internal forces
#pragma omp parallel for collapse(2)
		for (int i = 0; i < m_Particles.size(); i++)
		{
			FluidParticle& cur = m_Particles[i];

			glm::vec3 fPress(0.0);
			glm::vec3 fVisc(0.0);
			glm::vec3 n(0.0); //inward surface normal
			for (int j = 0; j < m_Particles.size(); j++)
			{
				FluidParticle& neighb = m_Particles[j];
				float distance = glm::length(cur.Position - neighb.Position);
				if (distance <= smoothRadius && i != j)
				{
					cur.NeighbId = neighb.Id;
					glm::vec3 grad = gradPressure(cur.Position - neighb.Position);
					fPress += (cur.Pressure / (cur.Density * cur.Density) + neighb.Pressure / (neighb.Density * neighb.Density)) * MASS * grad;
					//std::cout << "pressure: " << fPress.x << " " << fPress.y << " " << fPress.z;
					//std::cout << "	kern: " << grad.x << " " << grad.y << " " << grad.z << std::endl;
					fVisc += (neighb.Velocity - cur.Velocity) * (MASS / neighb.Density) * laplVisc(cur.Position - neighb.Position);
					//std::cout << "visc: " << fVisc.x << " " << fVisc.y << " " << fVisc.z << std::endl;
					n += (MASS / neighb.Density) * gradDefault(cur.Position - neighb.Position);
				}

			}
			cur.PressureForce = -(fPress*cur.Density);
			cur.ViscosityForce = fVisc * visc;
			cur.SurfaceNormal = n;
		}


		//external forces
#pragma omp parallel for collapse(2)
		for (int i = 0; i < m_Particles.size(); i++)
		{
			FluidParticle& currPart = m_Particles[i];
			currPart.GravityForce = currPart.Density * g;
			float colorFieldLapl = 0.0;
			for (int j = 0; j < m_Particles.size(); j++)
			{
				FluidParticle& neighbPart = m_Particles[j];

				float distance = glm::length(currPart.Position - neighbPart.Position);
				if (distance <= smoothRadius)
					colorFieldLapl += (MASS / neighbPart.Density) * laplDefault(currPart.Position - neighbPart.Position);
			}

			//currPart.SurfaceForce = glm::vec3(0.0);
			//if (glm::length(n) >= l)
			//float ln = glm::length(n);
			currPart.SurfaceForce = -surf_tens * colorFieldLapl * currPart.SurfaceNormal;
		}


		advance(grid);
		//m_Time += m_Dt;
	}

	void Draw(const Shader& shader, int selected_part)
	{
		shader.setVec3("dirLight.dir", glm::vec3(-0.1, -0.7, 0.2));
		shader.setVec3("dirLight.color", glm::vec3(1.0));
		shader.setVec3("material.ka", glm::vec3(0.2f));
		shader.setVec3("material.kd", glm::vec3(0.7f));
		shader.setVec3("material.ks", glm::vec3(1.0));
		shader.setFloat("material.ksh", 4.0f);

		if (!m_Sphere)
			m_Sphere = std::make_unique<Sphere>(10, 10, 1, glm::vec3(0.0, 0.0, 0.0));

		for (int i = 0; i < m_Particles.size(); i++)
		{
			float r = s * cbrt(3 * MASS / (4 * PI * m_Particles[i].Density));
			r = 0.01f;
			glm::mat4 model = glm::mat4(1.0);
			model = glm::translate(model, m_Particles[i].Position);
			model = glm::scale(model, glm::vec3(r));
			shader.setMat4("model", model);
			if (i == selected_part)
				shader.setVec3("myColor", glm::vec3(1.0, 1.0, 0.0));
			else
				shader.setVec3("myColor", glm::vec3(0.0, 0.2, 0.5));
				//shader.setVec3("myColor", glm::vec3(glm::length(m_Particles[i].PressureForce), 0.0, 1.0));
			m_Sphere->Draw();
		}

		for (int i = 0; i < m_BParticles.size(); i++)
		{
			float r = 0.01f;
			glm::mat4 model = glm::mat4(1.0);
			model = glm::translate(model, m_BParticles[i].Position);
			model = glm::scale(model, glm::vec3(r));
			shader.setMat4("model", model);
			shader.setVec3("myColor", glm::vec3(0.8, 0.8, 0.8));
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
		deltaT = dt;
	}

	float GetDeltaTime() const
	{
		return deltaT;
	}


	void PrintCoords() const
	{
		for (int i = 0; i < num; i++)
			std::cout << "[" << i << "] " << m_Particles[i].Position.x << ' ' << m_Particles[i].Position.y << ' ' << m_Particles[i].Position.z << std::endl;
	}

	void AddParticles(int n) {
		float offset = -(cbrt(n) / 2 * 0.15f);
		for (int i = 0; i < cbrt(n); i++) {
			for (int j = 0; j < cbrt(n); j++) {
				for (int k = 0; k < cbrt(n); k++) {
					float x = -0.2 + i * 0.025;
					float y = -0.05 + j * 0.025;
					float z = -0.15 + k * 0.025;
					FluidParticle tmp;
					tmp.Id = id++;
					tmp.Position = glm::vec3(x+ m_Origin.x, y+ m_Origin.y, z+ m_Origin.z);
					tmp.Velocity = glm::vec3(0.0, 0.0, 0.0);
					tmp.Acceleration = glm::vec3(0.0);
					//tmp.Mass = MASS;
					m_Particles.push_back(tmp);
				}
			}
		}
		num += n;
	}

	void Reset()
	{
		m_Particles.clear();
		m_BParticles.clear();
		num = 0;
		id = 0;
		AddParticles(init_num);
	}

	float* GetMass()
	{
		return &MASS;
	}

	float* GetVisc()
	{
		return &visc;
	}

	float* GetSurfTen()
	{
		return &surf_tens;
	}

	float* Getp0()
	{
		return &p0;
	}

	glm::vec3* GetGrav()
	{
		return &g;
	}

	float* GetDamping()
	{
		return &damping;
	}

	FluidParticle GetParticle(int id)
	{
		return m_Particles[id];
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
	void advance(Grid & grid)
	{
		omp_lock_t writelock;
		omp_init_lock(&writelock);
#pragma omp parallel for
		for (int i = 0; i < m_Particles.size(); i++)
		{
			FluidParticle& currPart = m_Particles[i];
			glm::vec3 F;
			glm::vec3 fInternal;
			glm::vec3 fExternal;
			glm::vec3 fBoundary;
			glm::vec3 velNext;
			glm::vec3 posNext;

			fInternal = currPart.PressureForce + currPart.ViscosityForce;
			fExternal = currPart.GravityForce + currPart.SurfaceForce;

			F = fInternal + fExternal;
			glm::vec3 acc = F / currPart.Density;
			currPart.Acceleration = acc; //for debug purposes

			//if initial velocity offset is not initialized
			//if (m_Time == 0.0)
			//{
			//	currPart.Velocity = currPart.Velocity - 0.5f * m_Dt * currPart.Acceleration;
			//}
			
			velNext = currPart.Velocity + acc * deltaT;
			posNext = currPart.Position + velNext * deltaT;

			grid.SeedCell(currPart.Position, m_BParticles, deltaS, &writelock);
			//fBoundary = grid.CalculateBoundaryForce(currPart);

			glm::vec3 contactP;
			glm::vec3 norm;
			if (grid.collision(currPart.Position, posNext, velNext, contactP, norm) && deltaT != 0)
			{
				float d = glm::length(posNext - contactP);
				if(glm::length(velNext))
					velNext = velNext - (float)(1 + cR * (d / (deltaT * glm::length(velNext)))) * glm::dot(velNext, norm) * norm;
				glm::vec3 velNextRev = -velNext;
				float k = glm::dot(glm::normalize(velNextRev), norm);
				float v_n_len = k * glm::length(velNext);
				glm::vec3 v_n = v_n_len * (norm);
				glm::vec3 v_t = velNext + v_n;
				glm::vec3 reflected = v_n + (-v_t);
				velNext = velNext + damping* reflected;
				glm::vec3 v = -velNext;
				//velNext = velNext + damping * v;
				posNext = contactP;
			}

			/*if (collisionS(posNext, contactP, norm) && deltaT != 0)
			{
				float d = glm::length(posNext - contactP);
				velNext = velNext - norm * (float)(1 + 0.5f * d / (deltaT * glm::length(velNext))) * glm::dot(velNext, norm);
				posNext = contactP;
			}*/

			currPart.Velocity = velNext;
			currPart.Position = posNext;
			
		}
		omp_destroy_lock(&writelock);
	}

	bool collisionS(glm::vec3 pos, glm::vec3& contactP, glm::vec3& normal) {
		if (abs(pos.x) < len && abs(pos.y) < len && abs(pos.z) < len) {
			return false;
		}
		float x = pos.x;
		float y = pos.y;
		float z = pos.z;
		char max = 'x';
		float maxP = abs(x);
		if (maxP < abs(y)) {
			max = 'y';
			maxP = abs(y);
		}
		if (maxP < abs(z)) {
			max = 'z';
			maxP = abs(z);
		}
		contactP = pos;
		switch (max) {
		case 'x':
			if (x < -len) {
				contactP.x = -len;
				normal = glm::vec3(1, 0, 0);
			}
			else {
				contactP.x = len;
				normal = glm::vec3(-1, 0, 0);
			}
			break;
		case'y':
			if (y < -len) {
				contactP.y = -len;
				//print(contactP);
				normal = glm::vec3(0, 1, 0);
			}
			else {
				contactP.y = len;
				normal = glm::vec3(0, -1, 0);
			}
			break;
		case'z':
			if (z < -len) {
				contactP.z = -len;
				normal = glm::vec3(0, 0, 1);
			}
			else {
				contactP.z = len;
				normal = glm::vec3(0, 0, -1);
			}
			break;
		}
		return true;
	}


	float kernDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		if (len > h) 
			return 0.0f;
		return (float)(315.0f / (64.0f * PI * powf(h, 9.0f))) * powf((h*h - len*len), 3.0f);
	}

	glm::vec3 gradDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		return -(r * (float)(945.0f / (32.0f * PI * powf(h, 9.0f))) * powf((h*h - len*len), 2.0f));
	}

	float laplDefault(const glm::vec3& r)
	{
		float len = glm::length(r);
		return -(945.0f / (32.0f * PI * powf(h, 9.0f))) * (h*h - len*len) * (3 * h*h - 7 * len*len);
	}

	glm::vec3 gradPressure(const glm::vec3& r)
	{

		float dist = glm::length(r);
		if (dist < 10e-5) {
			return glm::vec3(1.0f) * (float)(-45.f / (PI * powf(h, 6.0f))) * powf(h - dist, 2.0f);
		}
		else {
			glm::vec3 returnValue = (float)(-45.f / (PI * powf(h, 6.0f))) * glm::normalize(r) * powf(h - dist, 2.0f);

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
		float deltaT = 0.0f;
		float m_Time = 0.0f;
		float p0 = 998.29f;
		float MASS = 0.02f;
		float visc = 3.5f;
		float surf_tens = 0.0728f;
		float k = 3.0f;
		//const float l = 7.065;
		float l;
		float cR = 0.5f;
		//const int x = 20;
		float h = 0.0457f;
		int num = 10;
		int init_num;
		float vol;
		float smoothRadius;
		float s = 1;
		float len = 0.2f;
		int id;
		float damping = 0.1f;
		float deltaS;

private:
	std::vector<FluidParticle> m_Particles;
	std::vector<FluidParticle> m_BParticles;
	std::unique_ptr<Sphere> m_Sphere;
	glm::vec3 m_Origin;
};