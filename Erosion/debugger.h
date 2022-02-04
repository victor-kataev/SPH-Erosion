#pragma once
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <vector>


#define DEBUG_SEDIMENTATION_DISPLAY_FLAG	0x01
#define DEBUG_DEPOSITION_DISPLAY_FLAG		0x01 << 1
#define DEBUG_EROSION_DISPLAY_FLAG			0x01 << 2



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

	void PushBackPostDeposition(int ij, float dC, const FluidParticle& fp, const FluidParticle& bp)
	{
		FluidParticleLight p1, p2;
		p1.id = fp.Id;
		p1.sedim = fp.sedim;
		p1.sedim_delta = fp.sedim_delta;
		p1.sedim_ratio = fp.sedim_ratio;

		p2.id = bp.Id;
		p2.dM = bp.dM;

		m_PostDeposition[ij] = { ij, dC, p1, p2 };
	}

	void DisplayDebugWindow(uint8_t display_flags) const
	{
		ImGui::Begin("Erosion Debugger");
		if (DEBUG_SEDIMENTATION_DISPLAY_FLAG & display_flags)
			displaySedimentation();
		if (DEBUG_DEPOSITION_DISPLAY_FLAG & display_flags)
			displayDeposition();
		ImGui::End();
	}

	//void ClearBuffers()
	//{
	//	m_PostSedimentation.clear();
	//}

	void PostSedimentationBufferInit(int size)
	{
		m_PostSedimentation.clear();
		m_PostSedimentation.resize(size);
	}

	void PostDepositionBufferInit(int size)
	{
		m_PostDeposition.clear();
		m_PostDeposition.resize(size);
	}

	//std::vector<PPInteraction>& GetPostSedimentationBuffer()
	//{
	//	return m_PostSedimentation;
	//}


private:
	static Debugger* m_Instance;
	

	std::vector<PPInteraction> m_PostSedimentation;
	std::vector<PPInteraction> m_PostDeposition;

	void displaySedimentation() const
	{
		if (!ImGui::CollapsingHeader("Post SEDIMENTATION"))
			return;

		for (int i = 0; i < m_PostSedimentation.size(); i++)
		{
			const auto& data = m_PostSedimentation[i];
			if (ImGui::TreeNode((void*)(intptr_t)i, "dC[%d]", i))
			{
				//ImGui::Text("ij: %d", data.ij);
				ImGui::Text("dC: %.15f", data.dC);
				ImGui::Text("particle_1");
				ImGui::Text("\tid: %d", data.particle1.id);
				ImGui::Text("\tsedim: %.10f", data.particle1.sedim);
				ImGui::Text("\tsedim_delta: %.15f", data.particle1.sedim_delta);
				ImGui::Text("\tsedim_ratio: %.10f", data.particle1.sedim_ratio);
				ImGui::Text("particle_2");
				ImGui::Text("\tid: %d", data.particle2.id);
				ImGui::Text("\tsedim: %.10f", data.particle2.sedim);
				ImGui::Text("\tsedim_delta: %.15f", data.particle2.sedim_delta);
				ImGui::Text("\tsedim_ratio: %.10f", data.particle2.sedim_ratio);
				ImGui::TreePop();
			}
		}
	}

	void displayDeposition() const
	{
		if (!ImGui::CollapsingHeader("Post DEPOSITION"))
			return;

		for (int i = 0; i < m_PostDeposition.size(); i++)
		{
			const auto& data = m_PostDeposition[i];
			const char* format;
			if (data.dC == -777.0f)
				format = "dC_BP[%d] (skipped) (sedim)";
			if (data.dC == -666.0f)
				format = "dC_BP[%d] (skipped) (settling)";
			else
				format = "dC_BP[%d]";
			if (ImGui::TreeNode((void*)(intptr_t)i, format, i))
			{
				//ImGui::Text("ij: %d", data.ij);
				ImGui::Text("dC_BP: %.15f (%e)", data.dC, data.dC);
				ImGui::Text("fluid_part");
				ImGui::Text("\tid: %d", data.particle1.id);
				ImGui::Text("\tsedim: %.10f", data.particle1.sedim);
				ImGui::Text("\tsedim_delta: %.15f", data.particle1.sedim_delta);
				ImGui::Text("\tsedim_ratio: %.10f", data.particle1.sedim_ratio);
				ImGui::Text("boundary_part");
				ImGui::Text("\tid: %d", data.particle2.id);
				ImGui::Text("\tdM: %.15f", data.particle2.dM);
				ImGui::TreePop();
			}
		}
	}

	Debugger()
	{

	}
};

Debugger* Debugger::m_Instance = nullptr;