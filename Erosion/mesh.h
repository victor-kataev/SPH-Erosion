#pragma once

#include <vector>

struct Mesh
{
	Mesh(const std::vector<float>& verts, const std::vector<unsigned int>& indices)
		: Vertices(verts), Indices(indices)
	{}

	std::vector<float> Vertices;
	std::vector<unsigned int> Indices;
};