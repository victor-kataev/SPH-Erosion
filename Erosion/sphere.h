#pragma once

#include <GL/glew.h>


#include <vector>
#define PI 3.141592f


struct Sphere
{
	Sphere(int sectorCount, int stackCount, float r, glm::vec3 center)
		: m_Radius(r)
	{
		float x, y, z, xy;
		float nx, ny, nz, lengthInv = 1.0f / r;
		
		float sectorStep = 2 * PI / sectorCount;
		float stackStep = PI / stackCount;
		float sectorAngle, stackAngle;

		for(int i = 0; i <= stackCount; i++)
		{
			stackAngle = PI / 2 - i * stackStep;
			xy = r * cosf(stackAngle);
			z = r * sinf(stackAngle) + center.z;

			for (int j = 0; j <= sectorCount; j++)
			{
				sectorAngle = j * sectorStep;

				x = xy * cosf(sectorAngle) + center.x;
				y = xy * sinf(sectorAngle) + center.y;
				m_Vertices.push_back(x);
				m_Vertices.push_back(y);
				m_Vertices.push_back(z);

				nx = (x - center.x) * lengthInv;
				ny = (y - center.y) * lengthInv;
				nz = (z - center.z) * lengthInv;
				m_Vertices.push_back(nx);
				m_Vertices.push_back(ny);
				m_Vertices.push_back(nz);
			}
		}

		int k1, k2;
		for (int i = 0; i < stackCount; ++i)
		{
			k1 = i * (sectorCount + 1);     
			k2 = k1 + sectorCount + 1; 

			for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
			{
				if (i != 0)
				{
					m_Indices.push_back(k1);
					m_Indices.push_back(k2);
					m_Indices.push_back(k1 + 1);
				}

				if (i != (stackCount - 1))
				{
					m_Indices.push_back(k1 + 1);
					m_Indices.push_back(k2);
					m_Indices.push_back(k2 + 1);
				}

				m_LineIndices.push_back(k1);
				m_LineIndices.push_back(k2);
				if (i != 0)
				{
					m_LineIndices.push_back(k1);
					m_LineIndices.push_back(k1 + 1);
				}
			}
		}

		SetUpBuffers();
	}

	~Sphere()
	{
		glDeleteVertexArrays(1, &VAO);
		glDeleteBuffers(1, &VBO);

		m_Vertices.clear();
		m_Indices.clear();
		m_LineIndices.clear();
	}

	void Draw() const
	{
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, m_Indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	}

	void SetUpBuffers()
	{
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);

		glGenBuffers(1, &VBO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, m_Vertices.size() * sizeof(float), &m_Vertices[0], GL_STATIC_DRAW);
		glGenBuffers(1, &EBO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_Indices.size() * sizeof(unsigned int), &m_Indices[0], GL_STATIC_DRAW);

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
	}

	std::vector<float> GetVertices() const
	{
		return m_Vertices;
	}

	std::vector<unsigned int> GetIndices() const
	{
		return m_Indices;
	}

	float m_Radius;
	std::vector<float> m_Vertices;
	std::vector<unsigned int> m_Indices;
	std::vector<unsigned int> m_LineIndices;
	unsigned int VAO;
	unsigned int VBO;
	unsigned int EBO;
};