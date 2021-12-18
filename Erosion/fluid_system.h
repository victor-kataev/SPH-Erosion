#pragma once
#include <memory>
#include <glm/gtc/matrix_transform.hpp>

#include "mesh.h"
#include "sphere.h"
#include "grid.h"


#define GRAVITY 9.8



typedef unsigned int uint;



#define SEDIMENT_MAX 0.7f
#define SOLID_DENSITY 1842.0f
#define DAMPING 110.0f
#define C_2_MASS(C) (C*MASS_C_COEFFICIENT)
#define MASS_2_C (m*INV_MASS_C_COEFFICIENT)


#define FLUID_MASS 0.02f
#define FLUID_BASE_DENSITY 998.29f

float MASS_C_COEFFICIENT = SOLID_DENSITY * FLUID_MASS / FLUID_BASE_DENSITY;
float INV_MASS_C_COEFFICIENT = 1.0f / (SOLID_DENSITY * FLUID_MASS / FLUID_BASE_DENSITY);
float SEDIMENT_MAX_MASS = C_2_MASS(SEDIMENT_MAX);

#define PI 3.141592f
float kd = 6.0f * PI * 1.78e-5;
float ks = 1.119e6;

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
		table_size = 2 * nParts;
		//id = 0;
		m_ParticlesSpatialHash.resize(table_size);
		m_ParticleNeighbors.resize(num);

		for (int i = 0; i < cbrt(num); i++)
			for (int j = 0; j < cbrt(num); j++)
				for (int k = 0; k < cbrt(num); k++)
				{
					float x = -0.2 + i * 0.025;
					float y = -0.05 + j * 0.025;
					float z = -0.15 + k * 0.025;

					FluidParticle particle;
					particle.Id = FluidParticle::IdCount++;
					particle.Position = glm::vec3(x + m_Origin.x, y + m_Origin.y, z + m_Origin.z);
					particle.Velocity = glm::vec3(0.0, 0.0, 0.0);
					particle.Acceleration = glm::vec3(0.0);
					particle.sedim = 0.0f;
					particle.sedim_delta = 0.0f;
					//particle.Mass = MASS;
					m_Particles.push_back(particle);
					
				}

		//vol = num * MASS / p0;
		//smoothRadius = cbrt(3 * vol * x / (4 * PI * num));
		smoothRadius = h;
		deltaS = h / 2.0f;
		//h = smoothRadius;
		//l = sqrt(p0 / x);
	}

	void Run(Grid& grid)
	{
		rebuildTable();
		//compute density and pressure
#pragma omp parallel for collapse(2)
		for (int i = 0; i < m_Particles.size(); i++)
		{
			FluidParticle& currPart = m_Particles[i];

			float density = 0;
			int cnt = 0;
			//for (int j = 0; j < m_Particles.size(); j++)
			for (int j = 0; j < m_ParticleNeighbors[i].size(); j++)
			{
				int idx = m_ParticleNeighbors[i][j];
				FluidParticle& neighbPart = m_Particles[idx];

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
			//for (int j = 0; j < m_Particles.size(); j++)
			for (int j = 0; j < m_ParticleNeighbors[i].size(); j++)
			{
				int idx = m_ParticleNeighbors[i][j];
				FluidParticle& neighb = m_Particles[idx];
				float distance = glm::length(cur.Position - neighb.Position);
				if (distance <= smoothRadius && i != idx)
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
			currPart.GravityForce = currPart.Density * gravityVector;
			float colorFieldLapl = 0.0;
			//for (int j = 0; j < m_Particles.size(); j++)
			for (int j = 0; j < m_ParticleNeighbors[i].size(); j++)
			{
				int idx = m_ParticleNeighbors[i][j];
				FluidParticle& neighbPart = m_Particles[idx];

				float distance = glm::length(currPart.Position - neighbPart.Position);
				if (distance <= smoothRadius)
					colorFieldLapl += (MASS / neighbPart.Density) * laplDefault(currPart.Position - neighbPart.Position);
			}

			//currPart.SurfaceForce = glm::vec3(0.0);
			//if (glm::length(n) >= l)
			//float ln = glm::length(n);
			currPart.SurfaceForce = -surf_tens * colorFieldLapl * currPart.SurfaceNormal;
		}

		computeBoundaryForces(grid);
		computeSedimentTransfer();
		//computeErosion();
		//computeDeposition();

		advance(grid);
		//m_Time += m_Dt;
	}

	void Draw(const Shader& shader, int selected_part)
	{
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
				shader.setVec3("myColor", glm::vec3(0.0, 0.0, 0.8));
				//shader.setVec3("myColor", glm::vec3(glm::length(m_Particles[i].PressureForce), 0.0, 1.0));
			m_Sphere->Draw();
		}

		//draw boundary particles
		/*for (int i = 0; i < m_BParticles.size(); i++)
		{
			float r = 0.01f;
			glm::mat4 model = glm::mat4(1.0);
			model = glm::translate(model, m_BParticles[i].Position);
			model = glm::scale(model, glm::vec3(r));
			shader.setMat4("model", model);
			shader.setVec3("myColor", glm::vec3(0.8, 0.8, 0.8));
			m_Sphere->Draw();
		}*/

		//hightlight nearest boundary particles
		if(render_boundary)
			for (const auto& bp : m_NearestBParticles)
			{
				float r = 0.01f;
				glm::mat4 model = glm::mat4(1.0);
				model = glm::translate(model, glm::vec3(bp.Position.x, bp.Position.y, bp.Position.z));
				model = glm::scale(model, glm::vec3(r));
				shader.setMat4("model", model);
				shader.setVec3("myColor", glm::vec3(0.6, 0.0, 0.0));
				m_Sphere->Draw();
			}

		m_NearestBParticles.clear();

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
					tmp.Id = FluidParticle::IdCount++;
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
		AddParticles(init_num);
	}

	float* GetMu()
	{
		return &mu;
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
		return &gravityVector;
	}

	float* GetDamping()
	{
		return &damping;
	}

	bool* GetRenderBoundary()
	{
		return &render_boundary;
	}

	float* GetKs()
	{
		return &Ks;
	}

	float* GetKd()
	{
		return &Kd;
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

	void computeErosion()
	{

	}

	void computeDeposition()
	{
	
	}

	//advecton + diffustion
	void computeSedimentTransfer()
	{
		for (int i = 0; i < m_Particles.size(); i++)
		{
			FluidParticle& currPart = m_Particles[i];

			glm::vec3 vSettling = 2.0f / 9.0f * 0.001f * 0.001f * (float)((SOLID_DENSITY - currPart.Density) / visc) * gravityVector;
			float dC = 0.0f;
			float q = 0.0f;

			for (const auto& j : m_ParticleNeighbors[i])
			{
				FluidParticle& neigh = m_Particles[j];
				glm::vec3 rij = glm::normalize(currPart.Position - neigh.Position);
				float rij_len = glm::length(currPart.Position - neigh.Position);
				float v_r = glm::dot(vSettling, rij);


				//donor-acceptor
				if (v_r >= 0.0f)
				{
					//neighbor is donor
					if (neigh.sedim > 0.0f && currPart.sedim < SEDIMENT_MAX)
					{
						v_r *= richardson_zaki(currPart.sedim);
						q = MASS * neigh.sedim / neigh.Density;
						q *= -v_r* gradCubicSpline(rij_len);
						dC += q;
					}
				}
				else
				{
					//current is donor
					if (currPart.sedim > 0.0f && neigh.sedim < SEDIMENT_MAX)
					{
						v_r *= richardson_zaki(neigh.sedim);
						q = MASS * currPart.sedim / currPart.Density;
						q *= -v_r * gradCubicSpline(rij_len);
						dC += q;
					}
				}
			
				//diffusion
				dC += MASS / (currPart.Density * neigh.Density) * 0.1f * (currPart.sedim - neigh.sedim) * gradCubicSpline(rij_len);

				if (dC <= 0.0f) 
					currPart.sedim_delta += dC;
				else
					neigh.sedim_delta -= dC;
				dC = 0.0f;

			}
		}
	}

	void computeBoundaryForces(Grid& grid)
	{
		omp_lock_t writelock;
		omp_init_lock(&writelock);

#pragma omp parallel for
		for (int i = 0; i < m_Particles.size(); i++)
		{
			FluidParticle& currPart = m_Particles[i];
			usetfp nearest_boundary;
			glm::vec3 fBoundary(0.0f);


			omp_set_lock(&writelock);
			grid.SeedCell(currPart.Position, m_BParticles, deltaS);
			nearest_boundary = grid.FindNearestBoundary(currPart.Position, smoothRadius);
			omp_unset_lock(&writelock);

			if (!nearest_boundary.empty())
			{
				float shortest_dist = UINT_MAX;
				FluidParticle closestbp;
				for (const auto& bp : nearest_boundary)
				{
					//no-slip
					fBoundary += deltaS * deltaS * (-mu) * currPart.Velocity * laplVisc(currPart.Position - bp.Position);
					float dist = glm::length(bp.Position - currPart.Position);
					if (dist < shortest_dist)
					{
						shortest_dist = dist;
						closestbp = bp;
					}
				}

				//no penetration
				/*if (glm::dot(glm::normalize(closestbp.Position - currPart.Position), closestbp.SurfaceNormal) >= 0)
				{
					fBoundary += (Ks * shortest_dist - glm::dot(currPart.Velocity, closestbp.SurfaceNormal) * Kd) * closestbp.SurfaceNormal;
					currPart.shortest = shortest_dist;//debug only
					currPart.underSurf = true;
					currPart.lastDetectedBoundaryPos = closestbp.Position;
					currPart.lastDetectedBoundaryNorm = closestbp.SurfaceNormal;
				}
				else
				{
					currPart.underSurf = false;
				}*/

			}
			/*else if (nearest_boundary.empty() && currPart.underSurf)
			{
				nearest_boundary = grid.FindNearestBoundary(currPart.Position, smoothRadius*10);
				FluidParticle closestbp;
				float shortest_dist = UINT_MAX;
				for (const auto& bp : nearest_boundary)
				{
					float dist = glm::length(bp.Position - currPart.Position);
					if (dist < shortest_dist)
					{
						shortest_dist = dist;
						closestbp = bp;
					}
				}
				if(!nearest_boundary.empty())
					fBoundary += (Ks * shortest_dist - glm::dot(currPart.Velocity, closestbp.SurfaceNormal) * Kd) * closestbp.SurfaceNormal;
			}*/
			currPart.fBoundary = fBoundary;//debug only

			omp_set_lock(&writelock);
			m_NearestBParticles.merge(nearest_boundary);
			omp_unset_lock(&writelock);
		}
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
			glm::vec3 velNext;
			glm::vec3 posNext;




			fInternal = currPart.PressureForce + currPart.ViscosityForce;
			fExternal = currPart.GravityForce + currPart.SurfaceForce + currPart.fBoundary;

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

			
			glm::vec3 contactP;
			glm::vec3 norm;
			if (grid.collision(currPart.Position, posNext, velNext, contactP, norm) && deltaT != 0)
			{
				float d = glm::length(posNext - contactP);
				if(glm::length(velNext))
					velNext = velNext - (float)(1 + cR * (d / (deltaT * glm::length(velNext)))) * glm::dot(velNext, norm) * norm;
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

	float gradCubicSpline(float rij)
	{
		float h1 = 1.0f / h;
		float q = rij * h1;

		//kernel normalizing factor
		float factor = 1.0f / PI * h1 * h1 * h1;
		float tmp2 = 2.0f - q;
		float val;
		if (rij > 1e-12)
			if (q > 2.0f)
				val = 0.0f;
			else if (q > 1.0f)
				val = -0.75f * tmp2 * tmp2;
			else
				val = -3.0f * q * (1 - 0.75 * q);
		else
			val = 0.0f;

		return val * factor;
	}

	size_t hash(const glm::vec3& pos)
	{
		static float cc = 1.0f / h;
		int rx = floor(pos.x * cc);
		int ry = floor(pos.y * cc);
		int rz = floor(pos.z * cc);
		unsigned long res = (rx * 73856093) ^ (ry * 19349663) ^ (rz * 83492791);

		return res % table_size;
	}

	void rebuildTable()
	{
		m_ParticlesSpatialHash.clear();
		m_ParticleNeighbors.clear();
		m_ParticlesSpatialHash.resize(table_size);
		m_ParticleNeighbors.resize(num);

		for (int i = 0; i < m_Particles.size(); i++)
			m_ParticlesSpatialHash[hash(m_Particles[i].Position)].push_back(i);

#pragma omp parallel for collapse(4)
		for (int i = 0; i < m_Particles.size(); i++)
		{
			glm::vec3& r_q = m_Particles[i].Position;
			glm::vec3 bbmin = r_q - glm::vec3(h);
			glm::vec3 bbmax = r_q + glm::vec3(h);

			std::unordered_set<int> tmp;

			float a, b, c;
			for(int y = 0; y < 3; y++)
				for(int z = 0; z < 3; z++)
					for (int x = 0; x < 3; x++)
					{
						if (x == 0)
							a = bbmin.x;
						else if (x == 1)
							a = bbmin.x + h;
						else
							a = bbmax.x;
						
						if (y == 0)
							b = bbmin.y;
						else if (y == 1)
							b = bbmin.y + h;
						else
							b = bbmax.y;
						
						if (z == 0)
							c = bbmin.z;
						else if (z == 1)
							c = bbmin.z + h;
						else
							c = bbmax.z;
							
						std::vector<int>& cellParts = m_ParticlesSpatialHash[hash(glm::vec3(a, b, c))];
						for (const auto part_idx : cellParts)
						{
							glm::vec3& r_j = m_Particles[part_idx].Position;
							if (glm::length(r_q - r_j) <= h)
								tmp.insert(part_idx);
						}
					}
			m_ParticleNeighbors[i].insert(m_ParticleNeighbors[i].begin(), tmp.begin(), tmp.end());
		}
	}

	float richardson_zaki(float C)
	{
		if (C > SEDIMENT_MAX)
			return 0.0f;
		float e = 4.5;
		return (1 - powf(C / SEDIMENT_MAX, e));
	}

private:
		//Simulation parameters
		//float m_Param[MAX_PARAM];
		//FBufs m_Fluid;
	bool render_boundary = false;

		glm::vec3 gravityVector = glm::vec3(0.0, -9.82f, 0.0);
		float deltaT = 0.01f;
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
		int num;
		int init_num;
		float vol;
		float smoothRadius;
		float s = 1;
		float len = 0.2f;
		float damping = 0.1f;
		float deltaS;
		float mu = 0.27f;
		float Ks = 1.119e6;
		//float Ks = 45000;
		float Kd = 300.0f;
		int table_size;

private:
	std::vector<FluidParticle> m_Particles;
	std::vector<std::vector<int>> m_ParticlesSpatialHash;
	std::vector<std::vector<int>> m_ParticleNeighbors;
	std::vector<FluidParticle> m_BParticles;
	usetfp m_NearestBParticles;
	std::unique_ptr<Sphere> m_Sphere;
	glm::vec3 m_Origin;
};