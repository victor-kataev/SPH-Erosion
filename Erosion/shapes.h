#pragma once

#include <GL/glew.h>

#include "mesh.h"


class Shape
{
public:
	Shape() = default;

	void CreateCube()
	{
		std::vector<float> verts{
            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f, //back
             0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
             0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
             0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
            -0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,

            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f, 1.0f, //front
             0.5f, -0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
             0.5f,  0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
             0.5f,  0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
            -0.5f,  0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f, 1.0f,

            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f, //left
            -0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
            -0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,

             0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f, //right
             0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
             0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
             0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
             0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
             0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,

            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f, //bottom
             0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
             0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
             0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
            -0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f
	    };

        std::vector<unsigned int> inds;
        for (int i = 0; i < 30; i++)
            inds.push_back(i);

        cubeMesh.Vertices = verts;
        cubeMesh.Indices = inds;
            
        glGenVertexArrays(1, &cubeVAO);
        glBindVertexArray(cubeVAO);

        glGenBuffers(1, &cubeVBO);
        glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
        glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), &cubeMesh.Vertices[0], GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
	}
	
    void DrawCube()
    {
        glBindVertexArray(cubeVAO);
        glDrawArrays(GL_TRIANGLES, 0, cubeMesh.Vertices.size() / 6);
        glBindVertexArray(0);
    }

    void CreateBowl()
    {
        std::vector<float> verts;
        std::vector<unsigned int> inds;

        int dimx = 21;
        int dimz = 8;
        float sectorL = 90.0f / (dimx / 2);
        for(int z = 0; z < dimz; z++)
            for (int x = 0; x < dimx; x++)
            {
                float y = sin(glm::radians(-sectorL * x));
                //verts.push_back((x-dimx/2+1.5)/10.0f);
                verts.push_back(x - dimx / 2);
                verts.push_back(y * 5.0f + 200);
                //verts.push_back(1.0);
                //verts.push_back((z-dimz/2)/5.0f);
                verts.push_back(z - dimz/2);

                //normals
                glm::vec3 u(0.0);
                glm::vec3 d(0.0);
                glm::vec3 r(0.0);
                glm::vec3 l(0.0);

                if (x - 1 >= 0)
                    l = glm::vec3(x - (x - 1), y - sin(glm::radians(-sectorL * (x-1))), z - z);
                if (x + 1 < 11)
                    r = glm::vec3((x + 1) - x, sin(glm::radians(-sectorL * (x + 1))) - y, z - z);
                if (z - 1 >= 0)
                    u = glm::vec3(x - x, y - y, z - (z - 1));
                if (z + 1 < 8)
                    d = glm::vec3(x - x, y - y, (z + 1) - z);

                glm::vec3 normal = glm::normalize(glm::cross(u, l) + glm::cross(u, r) + glm::cross(d, l) + glm::cross(d, r));
                verts.push_back(normal.x);
                verts.push_back(normal.y);
                verts.push_back(normal.z);
            }

        for (int z = 0, j = dimz-1; z < dimz && j >= 0; z++, j--)
            for (int x = 0, i = dimx-1; x < dimx && i >= 0; x++, i--)
            {
                if (x + 1 < dimx && z + 1 < dimz)
                {
                    inds.push_back(z * dimx + x); //position in vertex buffer
                    inds.push_back(z * dimx + (x + 1));
                    inds.push_back((z + 1) * dimx + x);
                }
                if (j - 1 >= 0 && i - 1 >= 0)
                {
                    inds.push_back(j * dimx + i);
                    inds.push_back(j * dimx + (i - 1));
                    inds.push_back((j - 1) * dimx + i);
                }
            }
       
        bowlMesh.Vertices = verts;
        bowlMesh.Indices = inds;

        glGenVertexArrays(1, &bowlVAO);
        glBindVertexArray(bowlVAO);

        glGenBuffers(1, &bowlVBO);
        glBindBuffer(GL_ARRAY_BUFFER, bowlVBO);
        glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), &verts[0], GL_STATIC_DRAW);

        glGenBuffers(1, &bowlEBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bowlEBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, inds.size() * sizeof(unsigned int), &inds[0], GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);
    }

    void DrawBowl()
    {
        glBindVertexArray(bowlVAO);
        glDrawElements(GL_TRIANGLES, bowlMesh.Indices.size(), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
    }
	
	Mesh cubeMesh;
    Mesh bowlMesh;


private:
    unsigned int cubeVAO;
    unsigned int cubeVBO;
    unsigned int cubeEBO;
    unsigned int bowlVAO;
    unsigned int bowlVBO;
    unsigned int bowlEBO;

};