#pragma once
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <vector>


#define DEBUG_SEDIMENTATION_FLAG	0x01
#define DEBUG_DEPOSITION_FLAG		0x01 << 1
#define DEBUG_EROSION_FLAG			0x01 << 2



struct FluidParticleLight
{
	unsigned long id;
	float sedim;
	float sedim_delta;
	float sedim_ratio;
	float dM;
	char triangle;
};


//particle - particle interaction
struct PPInteraction
{
	int ij;
	float dC;
	FluidParticleLight particle1;
	FluidParticleLight particle2;
};


//an attempt to produce a UI debugger
//singleton
class Debugger
{
public:
	static Debugger* Get()
	{
		if (!m_Instance)
			m_Instance = new Debugger();
		return m_Instance;
	}

	Debugger(const Debugger&) = delete; //copy constructor

	~Debugger()
	{
		if (m_Instance)
			delete m_Instance;
	}

	//in sedimentation only dC and sedim_delta are changed
	void PushBackPostSedimentation(int ij, float dC, const FluidParticle& fp1, const FluidParticle& fp2)
	{
		FluidParticleLight p1, p2;
		p1.id = fp1.Id;
		p1.sedim = fp1.sedim;
		p1.sedim_delta = fp1.sedim_delta;
		p1.sedim_ratio = fp1.sedim_ratio;
		
		p2.id = fp2.Id;
		p2.sedim = fp2.sedim;
		p2.sedim_delta = fp2.sedim_delta;
		p2.sedim_ratio = fp2.sedim_ratio;

		m_PostSedimentation[ij] = { ij, dC, p1, p2 };
	}

	void DisplayDebugWindow(uint8_t debug_flags) const
	{
		ImGui::Begin("Erosion Debugger");
		if (debug_flags & DEBUG_SEDIMENTATION_FLAG)
			displaySedimentation();
		ImGui::End();
	}

	void ClearBuffers()
	{
		m_PostSedimentation.clear();
	}

	void PostSedimentationBufferInit(int size)
	{
		m_PostSedimentation.resize(size);
	}

	std::vector<PPInteraction>& GetPostSedimentationBuffer()
	{
		return m_PostSedimentation;
	}


private:
	static Debugger* m_Instance;
	

	std::vector<PPInteraction> m_PostSedimentation;

	void displaySedimentation() const
	{
		if (!ImGui::CollapsingHeader("Post sedimentation"))
			return;

		for (int i = 0; i < m_PostSedimentation.size(); i++)
		{
			const auto& data = m_PostSedimentation[i];
			if (ImGui::TreeNode((void*)(intptr_t)i, "dC[%d]", i))
			{
				ImGui::Text("ij: %d", data.ij);
				ImGui::Text("dC: %d", data.dC);
				ImGui::Text("part1_id: %d", data.particle1.id);
				ImGui::TreePop();
			}
		}
	}

	Debugger()
	{

	}
};

Debugger* Debugger::m_Instance = nullptr;