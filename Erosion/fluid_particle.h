#pragma once

#include <glm/glm.hpp>

struct FluidParticle
{
	unsigned long Id;
	static unsigned long IdCount;
	glm::vec3 Position;
	glm::vec3 Velocity;
	glm::vec3 Acceleration;
	//float Mass;
	float Density;
	float Pressure;
	glm::vec3 PressureForce;
	glm::vec3 ViscosityForce;
	glm::vec3 GravityForce;
	glm::vec3 SurfaceForce;
	glm::vec3 SurfaceNormal;
	glm::vec3 fBoundary;
	float shortest;
	int NeighbId;
	int cnt;
	bool underSurf = false;
	glm::vec3 lastDetectedBoundaryPos;
	glm::vec3 lastDetectedBoundaryNorm;
	float sedim;
	float sedim_delta;
	mutable float dM;
	float sedim_ratio;
	char triangle; // 'A' - ABC, 'B' - AABC
	double lifetime;
};

unsigned long FluidParticle::IdCount = 0;
