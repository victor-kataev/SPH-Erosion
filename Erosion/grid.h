#pragma once
#include <glm/glm.hpp>
#include <vector>
#include <list>
#include <omp.h>
#include <cstring>
#include <iomanip>
#include <float.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdlib.h>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "voxel.h"
#include "shader.h"


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
	float dM;
	float sedim_ratio;
};

unsigned long FluidParticle::IdCount = 0;

struct FluidParticleHasher
{
	size_t operator() (const FluidParticle& obj) const
	{
		return std::hash<unsigned long>() (obj.Id);
	}
};

struct FluidParticleComparator
{
	bool operator() (const FluidParticle& obj1, const FluidParticle& obj2) const
	{
		return obj1.Id == obj2.Id;
	}
};

typedef std::unordered_set<FluidParticle, FluidParticleHasher, FluidParticleComparator> usetfp;


struct Triangle
{
	Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c)
		: A(a), B(b), C(c)
	{
		norm = glm::normalize(glm::cross(B - A, C - A));
	}

	glm::vec3 A;
	glm::vec3 B;
	glm::vec3 C;
	glm::vec3 norm;
};

bool cmpX(const FluidParticle& a, const FluidParticle& b)
{
	return a.Position.x < b.Position.x;
}

bool cmpY(const FluidParticle& a, const FluidParticle& b)
{
	return a.Position.y < b.Position.y;
}

bool cmpZ(const FluidParticle& a, const FluidParticle& b)
{
	return a.Position.z < b.Position.z;
}


struct KdNode
{
	FluidParticle* particle;
	KdNode* left_child;
	KdNode* right_child;
};

class Kdtree
{
private:
	std::vector<FluidParticle> m_Particles;
	KdNode* root;

	KdNode* buildTree(const std::vector<FluidParticle>::iterator& begin, const std::vector<FluidParticle>::iterator& end, int depth)
	{
		if (begin == end)
			return nullptr;

		int axis = depth % 3;

		switch (axis)
		{
		case 0: { std::sort(begin, end, cmpX); break; }
		case 1: { std::sort(begin, end, cmpY); break; }
		case 2: { std::sort(begin, end, cmpZ); break; }
		}

		int median = (end - begin) / 2;
		KdNode* node = new KdNode();
		node->particle = &(*(begin + median));
		node->left_child = buildTree(begin, begin + median, depth + 1);
		node->right_child = buildTree(begin + median+1, end, depth + 1);
		return node;
	}

	void nn(const glm::vec3& q, KdNode* node, int cd, usetfp& nearest, float sr)
	{
		if (node == nullptr)
			return;

		if (glm::length(q - node->particle->Position) <= sr)
			nearest.insert(*node->particle);

		if (q[cd] < node->particle->Position[cd])
		{
			nn(q, node->left_child, (cd + 1) % 3, nearest, sr);
			if(abs(q[cd] - node->particle->Position[cd]) <= sr)
				nn(q, node->right_child, (cd + 1) % 3, nearest, sr);
		}
		else
		{
			nn(q, node->right_child, (cd + 1) % 3, nearest, sr);
			if (abs(q[cd] - node->particle->Position[cd]) <= sr)
				nn(q, node->left_child, (cd + 1) % 3, nearest, sr);
		}
	}

public:
	Kdtree() = default;

	Kdtree(std::vector<FluidParticle> particles)
	{
		m_Particles = particles;
		root = buildTree(m_Particles.begin(), m_Particles.end(), 0);
	}

	usetfp NearestNeighbors(glm::vec3& point, float sr)
	{
		usetfp nearest;
		nn(point, root, 0, nearest, sr);
		return nearest;
	}
};

class Grid
{
private:
	glm::vec3 m_Dim;
	glm::vec2 m_CellSize;
	ImVec4 m_Color;

	unsigned int VAO, VBO, EBO;

	Voxel* m_Grid;
	//float* m_Grid;
	//std::unordered_map<int, bool> m_SeededCells;
	std::unordered_map<int, Kdtree*> m_SeededCells;

	std::vector<float> vertexData;
	std::vector<unsigned int> indexData;

	struct Heightfield
	{
		size_t dimX = 0;
		size_t dimY = 0;

		unsigned char* map = nullptr;
	};

	Heightfield m_Heightfield;

	void destroyGrid()
	{
		if (m_Grid)
			delete[] m_Grid;
		vertexData.clear();
		indexData.clear();
	}

	void setupGraphics()
	{
		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &EBO);
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(float), &vertexData[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexData.size() * sizeof(unsigned int), &indexData[0], GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
	}

	void createGrid()
	{
		m_Grid = new Voxel[m_Dim.x * m_Dim.y * m_Dim.z];

		for (int z = 0; z < m_Dim.z; z++)
			for (int x = 0; x < m_Dim.x; x++)
			{
				int y = GetHeightfieldAt(x, z);
				if (y >= m_Dim.y)
					y = m_Dim.y - 1;

				SetVoxel(x, y, z);
				//m_Grid[x + m_DimX * (y + m_DimY * z) + 0] = x * m_CellSize[0];
				//m_Grid[x + m_DimX * (y + m_DimY * z) + 1] = y;
				//m_Grid[x + m_DimX * (y + m_DimY * z) + 2] = z * m_CellSize[1];
				vertexData.push_back(x * m_CellSize[0]);
				vertexData.push_back(y);
				vertexData.push_back(z * m_CellSize[1]);

				//normals
				glm::vec3 u(0.0);
				glm::vec3 d(0.0);
				glm::vec3 r(0.0);
				glm::vec3 l(0.0);

				glm::vec3 n1(0.0);
				glm::vec3 n2(0.0);
				glm::vec3 n3(0.0);
				glm::vec3 n4(0.0);
				glm::vec3 n5(0.0);
				glm::vec3 n6(0.0);

				if (x - 1 >= 0 && z - 1 >= 0)
				{
					l = glm::vec3(x - 1 - x, GetHeightfieldAt(x - 1, z) - y, z - z);
					u = glm::vec3(x - x, GetHeightfieldAt(x, z - 1) - y, z - 1 - z);
					//n1 = glm::cross(l, u);
					n1 = glm::cross(u, l);
				}

				if (x + 1 < m_Dim.x && z - 1 >= 0)
				{
					u = glm::vec3(x - x, GetHeightfieldAt(x, z - 1)-y, z - 1 - z);
					glm::vec3 ur = glm::vec3(x + 1 - x, GetHeightfieldAt(x + 1, z - 1) - y, z - 1 - z);
					r = glm::vec3(x + 1 - x, GetHeightfieldAt(x + 1, z) - y, z - z);

					//n2 = glm::cross(u, ur);
					n2 = glm::cross(ur, u);
					//n3 = glm::cross(ur, r);
					n3 = glm::cross(r, ur);
				}

				if (x + 1 < m_Dim.x && z + 1 < m_Dim.z)
				{
					r = glm::vec3(x + 1 - x, GetHeightfieldAt(x + 1, z) - y, z - z);
					d = glm::vec3(x - x, GetHeightfieldAt(x, z + 1) - y, z + 1 - z);
					//n4 = glm::cross(r, d);
					n4 = glm::cross(d, r);
				}

				if (x - 1 >= 0 && z + 1 < m_Dim.z)
				{
					d = glm::vec3(x - x, GetHeightfieldAt(x, z + 1) - y, z + 1 - z);
					l = glm::vec3(x - 1 - x, GetHeightfieldAt(x - 1, z) - y, z - z);
					glm::vec3 dl = glm::vec3(x - 1 - x, GetHeightfieldAt(x - 1, z + 1) - y, z + 1 - z);
					//n5 = glm::cross(d, dl);
					n5 = glm::cross(dl, d);
					//n6 = glm::cross(dl, l);
					n6 = glm::cross(l, dl);
				}

				//if (x - 1 >= 0)
				//	l = glm::vec3(x - (x - 1), y - GetHeightfieldAt(x - 1, z), z - z);
				//if (x + 1 < m_Dim.x)
				//	r = glm::vec3((x + 1) - x, GetHeightfieldAt(x + 1, z) - y, z - z);
				//if (z - 1 >= 0)
				//	u = glm::vec3(x - x, y - GetHeightfieldAt(x, z - 1), z - (z - 1));
				//if (z + 1 < m_Dim.y)
				//	d = glm::vec3(x - x, GetHeightfieldAt(x, z + 1) - y, (z + 1) - z);

				//glm::vec3 normal = glm::normalize(glm::cross(u, l) + glm::cross(u, r) + glm::cross(d, l) + glm::cross(d, r));
				glm::vec3 normal = glm::normalize(n1 + n2 + n3 + n4 + n5 + n6);
				vertexData.push_back(normal.x);
				vertexData.push_back(normal.y);
				vertexData.push_back(normal.z);
			}
		generateIndices();
	}

	void generateIndices()
	{
		for (int z = 0, j = m_Dim.z - 1; z < m_Dim.z && j >= 0; z++, j--)
			for (int x = 0, i = m_Dim.x - 1; x < m_Dim.x && i >= 0; x++, i--)
			{
				if (x + 1 < m_Dim.x && z + 1 < m_Dim.z)
				{
					indexData.push_back(z * m_Dim.x + x); //position in vertex buffer
					indexData.push_back(z * m_Dim.x + (x + 1));
					indexData.push_back((z + 1) * m_Dim.x + x);
				}
				if (j - 1 >= 0 && i - 1 >= 0)
				{
					indexData.push_back(j * m_Dim.x + i);
					indexData.push_back(j * m_Dim.x + (i - 1));
					indexData.push_back((j - 1) * m_Dim.x + i);
				}
			}
	}

	bool rayIntersectsTriangle(const glm::vec3& pos, const glm::vec3& dir, const Triangle tri, float* t)
	{
		float u, v;
		glm::vec3 E1 = tri.B - tri.A;
		glm::vec3 E2 = tri.C - tri.A;
		glm::vec3 N = glm::cross(E1, E2);
		//glm::vec3 N = tri.norm;
		float det = -glm::dot(dir, N);
		float invdet = 1.0 / det;
		glm::vec3 AO = pos - tri.A;
		glm::vec3 DAO = glm::cross(AO, dir);
		u = glm::dot(E2, DAO) * invdet;
		v = -glm::dot(E1, DAO) * invdet;
		*t = glm::dot(AO, N) * invdet;
		return (abs(det) >= 1e-6 && *t >= 0.0 && u >= 0.0 && v >= 0.0 && (u + v) <= 1.0);
	}

	void loadHeightMapFromPicture(const char* pic_path)
	{
		int width, height, channels;
		unsigned char* img = stbi_load(pic_path, &width, &height, &channels, 1);
		if (img == NULL)
		{
			printf("Error in loading the image\n");
			exit(1);
		}
		printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);

		m_Heightfield.dimX = width;
		m_Heightfield.dimY = height;
		m_Heightfield.map = new unsigned char[width * height];
		memset(m_Heightfield.map, 0, width * height);
		memcpy(m_Heightfield.map, img, (size_t)m_Heightfield.dimX * m_Heightfield.dimY);

		stbi_image_free(img);
	}

	std::vector<glm::vec2> findCellsWithinRadius(const glm::vec3& pos, const float r)
	{
		float x_min, x_max, z_min, z_max;
		std::vector<glm::vec2> output;
		output.emplace_back(floor(pos.x), floor(pos.z));
		x_min = floor(pos.x - r);
		x_max = floor(pos.x + r);
		z_min = floor(pos.z - r);
		z_max = floor(pos.z + r);

		for (int x = x_min; x <= x_max; x++)
			for (int z = z_min; z <= z_max; z++)
				if (
					(x < m_Dim.x && z < m_Dim.z && x >= 0 && z >= 0)
					&&
					(glm::length(glm::vec2(x, z) - glm::vec2(pos.x, pos.z)) <= r
					|| glm::length(glm::vec2(x + 1, z) - glm::vec2(pos.x, pos.z)) <= r
					|| glm::length(glm::vec2(x, z + 1) - glm::vec2(pos.x, pos.z)) <= r
					|| glm::length(glm::vec2(x + 1, z + 1) - glm::vec2(pos.x, pos.z)) <= r

					|| glm::length(glm::vec2(x, pos.z) - glm::vec2(pos.x, pos.z)) <= r//???
					|| glm::length(glm::vec2(pos.x, z) - glm::vec2(pos.x, pos.z)) <= r//???
					|| glm::length(glm::vec2(x + 1, pos.z) - glm::vec2(pos.x, pos.z)) <= r//???
					|| glm::length(glm::vec2(pos.x, z + 1) - glm::vec2(pos.x, pos.z)) <= r)//???
					)
				{
					if(glm::vec2(x, z) != glm::vec2(floor(pos.x), floor(pos.z)))
						output.emplace_back(x, z);
				}
		return output;
	}

public:
	Grid(const char* pic_path, glm::vec3& dim)
		: m_Dim(dim),
		m_CellSize(1.0f, 1.0f),
		m_Color(0.31f, 0.23f, 0.16f, 1.00f)
	{
		loadHeightMapFromPicture(pic_path);
		createGrid();
		setupGraphics();
	}

	//Voxel GetVoxel(int x, int y, int z) const
	//{
	//	return m_Grid[x + m_DimX * (y + m_DimY * z)];
	//}

	void SetVoxel(int x, int y, int z)
	{
		//m_Grid[x + m_DimX * (y + m_DimY * z)] = val;
		//m_Grid[x + m_DimX * (y + m_DimY * z)].type = type;
		m_Grid[x + (int)m_Dim.x * (y + (int)m_Dim.y * z)].position = glm::vec3(x * m_CellSize[0], y, z * m_CellSize[1]);

		//todo
	}

	unsigned char GetHeightfieldAt(int x, int y)
	{
		return m_Heightfield.map[m_Heightfield.dimY * x + y];
	}

	//for debug
	void HeightFieldMax() const
	{
		int max = 0;
		for (int i = 0; i < m_Heightfield.dimX * m_Heightfield.dimY; i++)
			if (m_Heightfield.map[i] > max)
				max = m_Heightfield.map[i];
		//std::cout << "Heightfield max = " << max << " " << __LINE__ << std::endl;
	}


	void UpdateGrid(glm::vec3& dim)
	{
		m_Dim = dim;
		destroyGrid();
		createGrid();

		glDeleteBuffers(1, &VBO);
		glDeleteBuffers(1, &EBO);
		glDeleteVertexArrays(1, &VAO);
		setupGraphics();
	}

	void UpdateHeightMap(const char* pic_path)
	{
		delete[] m_Heightfield.map;
		loadHeightMapFromPicture(pic_path);
	}

	std::vector<Triangle> getCellTriangles(glm::vec2 cellIndex)
	{
		glm::vec3 A(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1]), cellIndex[1]);
		glm::vec3 AA(cellIndex[0] + 1, GetHeightfieldAt(cellIndex[0] + 1, cellIndex[1] + 1), cellIndex[1] + 1);
		glm::vec3 B(cellIndex[0] + 1, GetHeightfieldAt(cellIndex[0] + 1, cellIndex[1]), cellIndex[1]);
		glm::vec3 C(cellIndex[0], GetHeightfieldAt(cellIndex[0], cellIndex[1] + 1), cellIndex[1] + 1);

		std::vector<Triangle> ret;
		ret.emplace_back(C, B, A); //ABC = normal is opposite
		ret.emplace_back(AA, B, C);

		return ret;
	}

	//3D-DDA (modified)
	glm::vec2 findAdjacentCell(const glm::vec3& pos, const glm::vec3& dir, glm::vec2 cellIndex)
	{
		float t_x, t_z;
		glm::vec2 deltaT;
		glm::vec2 cellDim(1, 1);
		glm::vec2 rayOrigin(pos.x, pos.z);
		glm::vec2 dir2D(dir.x, dir.z);

		if (dir2D[0] < 0)
		{
			deltaT[0] = -cellDim[0] / dir2D[0];
			t_x = (floor(rayOrigin[0] / cellDim[0]) * cellDim[0] - rayOrigin[0]) / dir2D[0];
		}
		else if (dir2D[0] > 0)
		{
			deltaT[0] = cellDim[0] / dir2D[0];
			t_x = ((floor(rayOrigin[0] / cellDim[0]) + 1) * cellDim[0] - rayOrigin[0]) / dir2D[0];
		}
		else
		{
			deltaT[0] = 0;
			t_x = INFINITY;
		}

		if (dir2D[1] < 0)
		{
			deltaT[1] = -cellDim[1] / dir2D[1];
			t_z = (floor(rayOrigin[1] / cellDim[1]) * cellDim[1] - rayOrigin[1]) / dir2D[1];
		}
		else if (dir2D[1] > 0)
		{
			deltaT[1] = cellDim[1] / dir2D[1];
			t_z = ((floor(rayOrigin[1] / cellDim[1]) + 1) * cellDim[1] - rayOrigin[1]) / dir2D[1];
		}
		else
		{
			deltaT[1] = 0;
			t_z = INFINITY;
		}

		if (t_x < t_z)
		{
			if (dir2D[0] < 0)
				cellIndex[0]--;
			else if (dir2D[0] > 0)
				cellIndex[0]++;
		}
		else
		{
			if (dir2D[1] < 0)
				cellIndex[1]--;
			else if (dir2D[1] > 0)
				cellIndex[1]++;
		}

		if (cellIndex[0] < 0 || cellIndex[0] >= m_Dim.x || cellIndex[1] < 0 || cellIndex[1] >= m_Dim.z)
			cellIndex = glm::vec2(-1, -1);

		return cellIndex;
	}

	bool mappedOnTriangle(const Triangle tri, const glm::vec3& pos, const glm::vec3& dir, glm::vec3& cp, glm::vec3& norm)
	{
		float t;

		if (rayIntersectsTriangle(pos, dir, tri, &t))
		{
			//t += 0.000002;
			norm = dir;
			cp = pos + t * norm;
			return true;
		}
		return false;
	}

	bool mappedBetweenTriangles(const Triangle& tri1, const Triangle& tri2, const glm::vec3& pos, glm::vec3& cp, glm::vec3& norm)
	{
		float t;

		if (rayIntersectsTriangle(pos, glm::normalize(tri1.norm + tri2.norm), tri1, &t))
		{
			//t += 0.000002;
			norm = glm::normalize(tri1.norm + tri2.norm);
			cp = pos + t * norm;
			return true;
		}
		else if (rayIntersectsTriangle(pos, glm::normalize(tri1.norm + tri2.norm), tri2, &t))
		{
			//t += 0.000002;
			norm = glm::normalize(tri1.norm + tri2.norm);
			cp = pos + t * norm;
			return true;
		}

		return false;
	}

	float min3(float a, float b, float c)
	{
		if (a < b)
			if (a < c)
				return a;
			else
				return c;
		else if (b < c)
			return b;
		else
			return c;
	}

	bool cornerCaseABC(const Triangle ABC, const Triangle AABC, const glm::vec3& posCurr, const glm::vec3& posNext, glm::vec2 cellIndex, glm::vec3& cp, glm::vec3& norm)
	{
		glm::vec2 posCurr2D = { posCurr.x, posCurr.z };
		glm::vec2 origin = { floor(posCurr.x), floor(posCurr.z) };
		glm::vec2 right = { origin.x + 1, origin.y };
		glm::vec2 down = { origin.x, origin.y + 1 };
		float dorigin = glm::length(origin - posCurr2D);
		float dright = glm::length(right - posCurr2D);
		float ddown = glm::length(down - posCurr2D);
		float dmin = min3(dorigin, dright, ddown);
		std::vector<Triangle> triangles;
		std::vector<Triangle> tmp;

		if (dmin == dorigin)
		{
			triangles.push_back(ABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y - 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y - 1));
			triangles.push_back(tmp[1]); //AABC
			//std::cout << "origin\n";
		}
		else if (dmin == dright)
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y - 1));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y - 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "right\n";
		}
		else //dmin == ddown
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y + 1));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y + 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "down\n";
		}

		glm::vec3 n(0.0);
		for (const auto& tr : triangles)
			n += tr.norm;
		n = glm::normalize(n);

		for (const auto& tr : triangles)
			if (mappedOnTriangle(tr, posNext, n, cp, norm))
				return true;
		return false;
	}

	bool cornerCaseAABC(const Triangle ABC, const Triangle AABC, const glm::vec3& posCurr, const glm::vec3& posNext, glm::vec2 cellIndex, glm::vec3& cp, glm::vec3& norm)
	{
		glm::vec2 posCurr2D = { posCurr.x, posCurr.z };
		glm::vec2 origin = { floor(posCurr.x) + 1, floor(posCurr.z) + 1 };
		glm::vec2 left = { origin.x - 1, origin.y };
		glm::vec2 up = { origin.x, origin.y - 1 };
		float dorigin = glm::length(origin - posCurr2D);
		float dleft = glm::length(left - posCurr2D);
		float dup = glm::length(up - posCurr2D);
		float dmin = min3(dorigin, dup, dleft);
		std::vector<Triangle> triangles;
		std::vector<Triangle> tmp;

		if (dmin == dorigin)
		{
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y + 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y + 1));
			triangles.push_back(tmp[0]); //ABC
			//std::cout << "origin AABC\n";
			//std::cout << "posCurr: " << posCurr.x << " " << posCurr.z << std::endl;
			//std::cout << "cellIndex: " << cellIndex.x << " " << cellIndex.y << std::endl;
		}
		else if (dmin == dup)
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y - 1));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x + 1, cellIndex.y - 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "up\n";
		}
		else //dmin == dleft
		{
			triangles.push_back(ABC);
			triangles.push_back(AABC);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y));
			triangles.push_back(tmp[1]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x, cellIndex.y + 1));
			triangles.push_back(tmp[0]);
			tmp = getCellTriangles(glm::vec2(cellIndex.x - 1, cellIndex.y + 1));
			triangles.insert(triangles.end(), tmp.begin(), tmp.end());
			//std::cout << "left\n";
		}

		glm::vec3 n(0.0);
		for (const auto& tr : triangles)
			n += tr.norm;
		n = glm::normalize(n);

		for (const auto& tr : triangles)
			if (mappedOnTriangle(tr, posNext, n, cp, norm))
				return true;
		return false;
	}

	bool pointLaysOnTriangle(const glm::vec3& p, Triangle& tri)
	{
		float S = fabs(glm::dot(tri.B - tri.A, tri.C - tri.A)) / 2.0f;
		float a = fabs(glm::dot(tri.B - p, tri.C - p)) / 2.0f / S;
		float b = fabs(glm::dot(tri.C - p, tri.A - p)) / 2.0f / S;
		float c = 1 - a - b;
		return a >= 0 && a <= 1 && b >= 0 && b <= 1 && c >= 0 && c <= 1;
	}

	bool pointLaysOnPlane(const glm::vec3& p, const Triangle& tri)
	{
		const float myEpsilon = 1.192092896e-05F;
		float d = glm::dot(tri.norm, tri.B);
		float res = tri.norm.x * tri.C.x + tri.norm.y * tri.C.y + tri.norm.z * tri.C.z - d;
		//std::cout << "plane res: " << res << std::endl;
		if (fabs(res) < myEpsilon)
			return true;
		return false;
	}

	bool collision(const glm::vec3& posCurr, const glm::vec3& posNext, const glm::vec3& velNext, glm::vec3& contactP, glm::vec3& norm)
	{
		float t; //where t is R.Origin + t*R.Dir
		glm::vec3 dir = glm::normalize(velNext);
		glm::vec3 dirOpposite = -dir;

		////std::cout << posCurr.x << " " << posCurr.y << " " << posCurr.z << std::endl;
		glm::vec2 cellIndex(floor(posCurr.x), floor(posCurr.z));
		glm::vec2 cellIndexNext(floor(posNext.x), floor(posNext.z));
		if (cellIndex[0] < 0 || cellIndex[0] >= m_Dim.x - 1 || cellIndex[1] < 0 || cellIndex[1] >= m_Dim.z - 1
			|| cellIndexNext[0] < 0 || cellIndexNext[0] >= m_Dim.x - 1 || cellIndexNext[1] < 0 || cellIndexNext[1] >= m_Dim.z - 1)
			return false;

		//if in the same cell
		if (cellIndex == cellIndexNext)
		{
			////std::cout << "in the same cell\n";
			////std::cout << "dir: " << dir.x << " " << dir.y << " " << dir.z << std::endl;
			////std::cout << "diroposite: " << dirOpposite.x << " " << dirOpposite.y << " " << dirOpposite.z << std::endl;
			std::vector<Triangle> cellTriangles = getCellTriangles(cellIndex);
			Triangle ABC = cellTriangles[0];
			Triangle AABC = cellTriangles[1];

			glm::vec3 CPtmp;
			float dABC = INFINITY, dAABC = INFINITY;
			if (rayIntersectsTriangle(posNext, dirOpposite, ABC, &t))
			{
				CPtmp = posNext + t * dirOpposite;
				dABC = glm::length(CPtmp - posNext);
				//std::cout << "true1\n";
			}

			if (rayIntersectsTriangle(posNext, dirOpposite, AABC, &t))
			{
				CPtmp = posNext + t * dirOpposite;
				dAABC = glm::length(CPtmp - posNext);
				//std::cout << "true2\n";
			}


			//the cell is built out of two triangles
			//check intersection with both and try to map on appropriate one
			if (dABC < dAABC || (fabs(dABC - dAABC) < FLT_EPSILON && dABC != INFINITY)) //or dABC == dAABC
			{
				//std::cout << "ABC case\n";
				//try to map on the same triangle
				if (rayIntersectsTriangle(posNext, ABC.norm, ABC, &t))
				{
					//t += 0.000002;
					norm = ABC.norm;
					contactP = posNext + t * norm;
					//std::cout << "mapped on ABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
					return true;
				}
				else
				{
					glm::vec2 adjCellIndex = findAdjacentCell(posNext, glm::normalize(dir + ABC.norm), cellIndex);
					if (adjCellIndex == glm::vec2(-1, -1))
						return false;
					std::vector<Triangle> adjCellTriangles = getCellTriangles(adjCellIndex);
					Triangle ABC_adj = adjCellTriangles[0];
					Triangle AABC_adj = adjCellTriangles[1];
					glm::vec2 cellDir = { adjCellIndex.x - cellIndex.x, adjCellIndex.y - cellIndex.y };

					//intersection with hypotenuse
					if (cellDir.x > 0 || cellDir.y > 0)
					{
						//std::cout << "hypotenuse\n";
						if (mappedBetweenTriangles(ABC, AABC, posNext, contactP, norm))
						{
							//std::cout << "mapped between ABC AABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{

						//intersection with cathenus
						if (mappedBetweenTriangles(ABC, AABC_adj, posNext, contactP, norm))
						{
							//std::cout << "mapped between ABC AABC_adj, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC " << __LINE__ << std::endl;
							return true;
						}
					}
				}
			}
			else if (dAABC < dABC)
			{
				//std::cout << "AABC case\n";
				//try to map on the same triangle
				if (rayIntersectsTriangle(posNext, AABC.norm, AABC, &t))
				{
					//t += 0.000002;
					norm = AABC.norm;
					contactP = posNext + t * norm;
					//std::cout << "mapped on AABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
					return true;
				}
				else
				{
					glm::vec2 adjCellIndex = findAdjacentCell(posNext, glm::normalize(dir + AABC.norm), cellIndex);
					if (adjCellIndex == glm::vec2(-1, -1))
						return false;
					std::vector<Triangle> adjCellTriangles = getCellTriangles(adjCellIndex);
					Triangle ABC_adj = adjCellTriangles[0];
					Triangle AABC_adj = adjCellTriangles[1];
					glm::vec2 cellDir = { adjCellIndex.x - cellIndex.x, adjCellIndex.y - cellIndex.y };

					//intersection with hypotenuse
					if (cellDir.x < 0 || cellDir.y < 0)
					{
						//std::cout << "hypotenuse\n";
						if (mappedBetweenTriangles(ABC, AABC, posNext, contactP, norm))
						{
							//std::cout << "mapped between AABC ABC, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseAABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{

						//intersection with cathenus
						if (mappedBetweenTriangles(AABC, ABC_adj, posNext, contactP, norm))
						{
							//std::cout << "mapped between AABC ABC_adj, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						//corner case
						if (cornerCaseAABC(ABC, AABC, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC " << __LINE__ << std::endl;
							return true;
						}
					}
				}
			}
			//particle did not go through collision did not occur
			return false;
		}
		else
		{
			////std::cout << "in different cells\n";
			////std::cout << "cellIndex: " << cellIndex.x << " " << cellIndex.y << " | cellIndexNext: " << cellIndexNext.x << " " << cellIndexNext.y << std::endl;
			std::vector<Triangle> cellTriangles = getCellTriangles(cellIndex);
			Triangle ABC_poscurr = cellTriangles[0];
			Triangle AABC_poscurr = cellTriangles[1];
			std::vector<Triangle> cellNextTriangles = getCellTriangles(cellIndexNext);
			Triangle ABC_posnext = cellNextTriangles[0];
			Triangle AABC_posnext = cellNextTriangles[1];
			glm::vec2 cellDir = { cellIndexNext - cellIndex };

			//diagonal corner case
			if (cellDir.x != 0 && cellDir.y != 0)
			{
				if (cellDir.x + cellDir.y == -2)
				{
					if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
					{
						//std::cout << "corner case ABC diagonal " << __LINE__ << std::endl;
						return true;
					}
				}
				else if (cellDir.x + cellDir.y == 2)
				{
					if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
					{
						//std::cout << "corner case AABC diagonal " << __LINE__ << std::endl;
						return true;
					}
				}
				else
				{
					if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
					{
						//std::cout << "corner case ABC diagonal " << __LINE__ << std::endl;
						return true;
					}
				}
			}
			else
			{
				if (rayIntersectsTriangle(posNext, dirOpposite, ABC_poscurr, &t))
				{
					if (cellDir.x == -1 || cellDir.y == -1)
					{
						if (mappedBetweenTriangles(ABC_poscurr, AABC_posnext, posNext, contactP, norm))
						{
							//std::cout << "mapped between ABC_poscurr AABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{
						if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}

				}
				else if (rayIntersectsTriangle(posNext, dirOpposite, AABC_poscurr, &t))
				{
					if (cellDir.x == 1 || cellDir.y == 1)
					{
						if (mappedBetweenTriangles(AABC_poscurr, ABC_posnext, posNext, contactP, norm))
						{
							//std::cout << "mapped between AABC_poscurr ABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}
					else
					{
						if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
						{
							//std::cout << "corner case AABC celldiff " << __LINE__ << std::endl;
							return true;
						}
					}
				}
				else
				{
					glm::vec3 CPtmp;
					float dABC_posnext = INFINITY, dAABC_posnext = INFINITY;
					if (rayIntersectsTriangle(posCurr, dir, ABC_posnext, &t))
					{
						CPtmp = posCurr + t * dir;
						dABC_posnext = glm::length(CPtmp - posCurr);
					}

					if (rayIntersectsTriangle(posCurr, dir, AABC_posnext, &t))
					{
						CPtmp = posCurr + t * dir;
						dAABC_posnext = glm::length(CPtmp - posCurr);
					}


					if (dABC_posnext < dAABC_posnext || (fabs(dABC_posnext - dAABC_posnext) < FLT_EPSILON && dABC_posnext != INFINITY))
					{
						if (rayIntersectsTriangle(posNext, ABC_posnext.norm, ABC_posnext, &t))
						{
							//t += 0.000002;
							norm = ABC_posnext.norm;
							contactP = posNext + t * norm;
							//std::cout << "mapped on ABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cellDir.x == 1 || cellDir.y == 1)
						{
							if (mappedBetweenTriangles(ABC_posnext, AABC_poscurr, posNext, contactP, norm))
							{
								//std::cout << "mapped between ABC_posnext AABC_poscurr, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
								return true;
							}

							if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case AABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
						else
						{
							if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
					}
					else if (dAABC_posnext < dABC_posnext)
					{
						if (rayIntersectsTriangle(posNext, AABC_posnext.norm, AABC_posnext, &t))
						{
							//t += 0.000002;
							norm = AABC_posnext.norm;
							contactP = posNext + t * norm;
							//std::cout << "mapped on AABC_posnext, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
							return true;
						}

						if (cellDir.x == -1 || cellDir.y == -1)
						{
							if (mappedBetweenTriangles(AABC_posnext, ABC_poscurr, posNext, contactP, norm))
							{
								//std::cout << "mapped between AABC_posnext ABC_poscurr, cellindex: " << cellIndex.x << " " << cellIndex.y << " " << __LINE__ << std::endl;
								return true;
							}

							if (cornerCaseABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
						else
						{
							if (cornerCaseAABC(ABC_poscurr, AABC_poscurr, posCurr, posNext, cellIndex, contactP, norm))
							{
								//std::cout << "corner case ABC celldiff " << __LINE__ << std::endl;
								return true;
							}
						}
					}
				}
			}
			return false;
		}
	}

	void SeedCell(const glm::vec3 partpos, std::vector<FluidParticle>& boundary, float deltaS)
	{
		glm::vec2 ccell = {floor(partpos.x), floor(partpos.z)};
		if (ccell[0] < 0 || ccell[0] >= m_Dim.x || ccell[1] < 0 || ccell[1] >= m_Dim.z)
			return;

		std::vector<glm::vec2> cellsWithinRadius = findCellsWithinRadius(partpos, deltaS*2);

		int cellIdx;
		std::vector<FluidParticle> boundaryParts;
		for (const auto& cell : cellsWithinRadius)
		{
			cellIdx = cell[0] * 1000 + cell[1];

			if (m_SeededCells.find(cellIdx) == m_SeededCells.end())
			{
				fillCell(cell, boundaryParts, deltaS);
				m_SeededCells[cellIdx] = new Kdtree(boundaryParts);
				boundary.insert(boundary.end(), boundaryParts.begin(), boundaryParts.end());
			}
		}
	}

	void fillCell(const glm::vec2& cell, std::vector<FluidParticle>& boundary, float deltaS)
	{
		std::vector<Triangle> cellTriangles = getCellTriangles(cell);
		Triangle ABC = cellTriangles[0];
		Triangle AABC = cellTriangles[1];

		float num = glm::length(ABC.C - ABC.A) / deltaS;
		float step = 1.0f / num;

		float t;
		for (t = 0; t < 1; t += step)
		{
			glm::vec3 a = ABC.A + t * (ABC.B - ABC.A);
			glm::vec3 b = ABC.A + t * (ABC.C - ABC.A);
			//glm::vec3 c = ABC.A + t * (AABC.A - AABC.C);
			float n = glm::length(b - a) / deltaS;
			float s = 1.0f / n;

			for (float u = 0; u <= 1.0f; u += s)
			{
				glm::vec3 d = b + u * (a - b);
				//glm::vec3 d = a + u * (b - a);

				FluidParticle bp;
				bp.Id = FluidParticle::IdCount++;
				bp.Position = d;
				bp.Velocity = glm::vec3(0.0);
				bp.SurfaceNormal = ABC.norm;
				boundary.push_back(bp);
			}
		}

		num = glm::length(AABC.A - AABC.C) / deltaS;
		step = 1.0f / num;
		for (t = step; t < 1; t += step)
		{
			glm::vec3 a = ABC.A + t * (ABC.B - ABC.A);
			glm::vec3 c = ABC.A + t * (AABC.A - AABC.C);
			float n = glm::length(c - a) / deltaS;
			float s = 1.0f / n;
			for (float u = 0; u < 1; u += s)
			{
				glm::vec3 e = c + u * (a - c);
				//glm::vec3 e = a + u * (c - a);
				FluidParticle bp;
				bp.Id = FluidParticle::IdCount++;
				bp.Position = e;
				bp.Velocity = glm::vec3(0.0);
				bp.SurfaceNormal = AABC.norm;
				boundary.push_back(bp);
			}
		}
	}

	usetfp FindNearestBoundary(glm::vec3& pos, float sr)
	{
		/*glm::vec2 cell = {floor(pos.x), floor(pos.z)};
		if (cell[0] < 0 || cell[0] >= m_Dim.x || cell[1] < 0 || cell[1] >= m_Dim.z)
			return {};
			*/ //could be usefull later

		std::vector<glm::vec2> cellsWithinRadius = findCellsWithinRadius(pos, sr);

		int cellIdx;
		usetfp nearestParticles;
		for (const auto& cell : cellsWithinRadius)
		{
			cellIdx = cell[0] * 1000 + cell[1];

			if (m_SeededCells.find(cellIdx) != m_SeededCells.end())
				nearestParticles.merge(m_SeededCells[cellIdx]->NearestNeighbors(pos, sr));
		}
		return nearestParticles;
	}

	void Draw(Shader& shader)
	{
		glm::mat4 model = glm::mat4(1.0f);
		shader.setMat4("model", model);
		shader.setVec3("material.ka", glm::vec3(0.2f));
		shader.setVec3("material.kd", glm::vec3(0.7f));
		shader.setVec3("material.ks", glm::vec3(1.0));
		shader.setFloat("material.ksh", 4.0f);

		shader.setVec3("myColor", glm::vec3(m_Color.x, m_Color.y, m_Color.z));
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, indexData.size(), GL_UNSIGNED_INT, 0);

	}

	void SetColor(ImVec4& color)
	{
		m_Color = color;
	}

	ImVec4 GetColor() const
	{
		return m_Color;
	}
};

