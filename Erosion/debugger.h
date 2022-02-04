#pragma once
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <vector>


#define SEDIMENTATION_DISPLAY_FLAG	0x01
#define DEPOSITION_DISPLAY_FLAG		0x01 << 1
#define EROSION_DISPLAY_FLAG		0x01 << 2
#define SEDIMENT_FLOW_DISPLAY_FLAG	0x01 << 3



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

	void PushBackPostSedimentFlowSphBoundary(int ij, float dC, const FluidParticle& fp, const FluidParticle& bp)
	{
		FluidParticleLight p1, p2;
		p1.id = fp.Id;
		p1.sedim = fp.sedim;
		p1.sedim_delta = fp.sedim_delta;
		p1.sedim_ratio = fp.sedim_ratio;

		p2.id = bp.Id;
		p2.dM = bp.dM;

		m_PostSedimentFlowSphBoundary[ij] = { ij, dC, p1, p2 };
	}
	
	void PushBackPostSedimentFlowSphSph(int ij, float dC, const FluidParticle& fp1, const FluidParticle& fp2)
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

		m_PostSedimentFlowSphSph[ij] = { ij, dC, p1, p2 };
	}

	void PushBackPostErosion(int ij, float dC, const FluidParticle& fp, const FluidParticle& bp, float m)
	{
		FluidParticleLight p1, p2;
		p1.id = fp.Id;
		p1.sedim = fp.sedim;
		p1.sedim_delta = fp.sedim_delta;
		p1.sedim_ratio = fp.sedim_ratio;
		p1.dM = m;

		p2.id = bp.Id;
		p2.dM = bp.dM;

		m_PostErosion[ij] = { ij, dC, p1, p2 };
	}

	void DisplayDebugWindow(uint8_t display_flags) const
	{
		ImGui::Begin("Erosion Debugger");
		if (SEDIMENTATION_DISPLAY_FLAG & display_flags)
			displaySedimentation();
		if (DEPOSITION_DISPLAY_FLAG & display_flags)
			displayDeposition();
		if (SEDIMENT_FLOW_DISPLAY_FLAG & display_flags)
			displaySedimentFlow();
		if (EROSION_DISPLAY_FLAG & display_flags)
			displayErosion();
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

	void PostSedimentFlowSphSphInit(int size)
	{
		m_PostSedimentFlowSphSph.clear();
		m_PostSedimentFlowSphSph.resize(size);
	}
	
	void PostSedimentFlowSphBoundaryInit(int size)
	{
		m_PostSedimentFlowSphBoundary.clear();
		m_PostSedimentFlowSphBoundary.resize(size);
	}

	void PostErosionInit(int size)
	{
		m_PostErosion.clear();
		m_PostErosion.resize(size);
	}

	//std::vector<PPInteraction>& GetPostSedimentationBuffer()
	//{
	//	return m_PostSedimentation;
	//}


private:
	static Debugger* m_Instance;
	

	std::vector<PPInteraction> m_PostSedimentation;
	std::vector<PPInteraction> m_PostDeposition;
	std::vector<PPInteraction> m_PostErosion;
	std::vector<PPInteraction> m_PostSedimentFlowSphBoundary;
	std::vector<PPInteraction> m_PostSedimentFlowSphSph;

	void displaySedimentation() const
	{
		if (!ImGui::CollapsingHeader("Post SEDIMENTATION"))
			return;

		ImGui::Text("--Total %d--", m_PostSedimentation.size());
		for (int i = 0; i < m_PostSedimentation.size(); i++)
		{
			const auto& data = m_PostSedimentation[i];
			const char* format;
			if (data.dC == -777.0f)
				format = "dC[%d] (skipped)";
			else
				format = "dC[%d]";
			if (ImGui::TreeNode((void*)(intptr_t)i, format, i))
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

		ImGui::Text("--Total %d--", m_PostDeposition.size());
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
				ImGui::Text("boundary_part");
				ImGui::Text("\tid: %d", data.particle2.id);
				ImGui::Text("\tdM: %.15f", data.particle2.dM);
				ImGui::Text("fluid_part");
				ImGui::Text("\tid: %d", data.particle1.id);
				ImGui::Text("\tsedim: %.10f", data.particle1.sedim);
				ImGui::Text("\tsedim_delta: %.15f", data.particle1.sedim_delta);
				ImGui::Text("\tsedim_ratio: %.10f", data.particle1.sedim_ratio);
				ImGui::TreePop();
			}
		}
	}

	void displaySedimentFlow() const
	{
		if (!ImGui::CollapsingHeader("Post SEDIMENT FLOW"))
			return;

		if (ImGui::TreeNode("sph-sph"))
		{
			ImGui::Text("--Total %d--", m_PostSedimentFlowSphSph.size());
			for (int i = 0; i < m_PostSedimentFlowSphSph.size(); i++)
			{
				const auto& data = m_PostSedimentFlowSphSph[i];
				const char* format;
				if (data.dC == -777.0f)
					format = "dC[%d] (skipped)";
				else
					format = "dC[%d]";

				if (ImGui::TreeNode((void*)(intptr_t)i, format, i))
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

			ImGui::TreePop();
		}
		
		if (ImGui::TreeNode("sph-boundary"))
		{
			ImGui::Text("--Total %d--", m_PostSedimentFlowSphBoundary.size());
			for (int i = 0; i < m_PostSedimentFlowSphBoundary.size(); i++)
			{
				const auto& data = m_PostSedimentFlowSphBoundary[i];
				const char* format;
				if (data.dC == -777.0f)
					format = "dC[%d] (skipped)";
				else
					format = "dC[%d]";

				if (ImGui::TreeNode((void*)(intptr_t)i, format, i))
				{
					//ImGui::Text("ij: %d", data.ij);
					ImGui::Text("dC_BP: %.15f (%e)", data.dC, data.dC);
					ImGui::Text("boundary_part");
					ImGui::Text("\tid: %d", data.particle2.id);
					ImGui::Text("\tdM: %.15f", data.particle2.dM);
					ImGui::Text("fluid_part");
					ImGui::Text("\tid: %d", data.particle1.id);
					ImGui::Text("\tsedim: %.10f", data.particle1.sedim);
					ImGui::Text("\tsedim_delta: %.15f", data.particle1.sedim_delta);
					ImGui::Text("\tsedim_ratio: %.10f", data.particle1.sedim_ratio);
					ImGui::TreePop();
				}
			}

			ImGui::TreePop();
		}
	}

	void displayErosion() const
	{
		if (!ImGui::CollapsingHeader("Post EROSION"))
			return;

		ImGui::Text("--Total %d--", m_PostErosion.size());
		for (int i = 0; i < m_PostErosion.size(); i++)
		{
			const auto& data = m_PostErosion[i];
			const char* format;
			if (data.dC == -666.0f)
				format = "[%d] (skipped)";
			else
				format = "[%d]";
			if (ImGui::TreeNode((void*)(intptr_t)i, format, i))
			{
				//ImGui::Text("ij: %d", data.ij);
				ImGui::Text("dC_BP: %.15f (%e)", data.dC, data.dC);
				ImGui::Text("boundary_part");
				ImGui::Text("\tid: %d", data.particle2.id);
				ImGui::Text("\tdM: %.15f", data.particle2.dM);
				ImGui::Text("\tm: %.15f", data.particle1.dM);
				ImGui::Text("fluid_part");
				ImGui::Text("\tid: %d", data.particle1.id);
				ImGui::Text("\tsedim: %.10f", data.particle1.sedim);
				ImGui::Text("\tsedim_delta: %.15f", data.particle1.sedim_delta);
				ImGui::Text("\tsedim_ratio: %.10f", data.particle1.sedim_ratio);
				ImGui::TreePop();
			}
		}
	}

	Debugger()
	{

	}
};

Debugger* Debugger::m_Instance = nullptr;