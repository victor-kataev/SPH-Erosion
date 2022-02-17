#pragma once
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <vector>
#include <unordered_map>
#include <map>

#include "fluid_particle.h"

//#define UI_DEBUG



#define SEDIMENTATION_DISPLAY_FLAG				0x01
#define DEPOSITION_DISPLAY_FLAG					0x01 << 1
#define EROSION_DISPLAY_FLAG					0x01 << 2
#define SEDIMENT_FLOW_DISPLAY_FLAG				0x01 << 3
#define CELLS_DISPLAY_FLAG						0x01 << 4
#define PARTICLES_SEDIMENTATION_DISPLAY_FLAG	0x01 << 5
#define PARTICLES_DEPOSITION_DISPLAY_FLAG		0x01 << 6
#define PARTICLES_SEDIMENT_FLOW_DISPLAY_FLAG	0x01 << 7
#define PARTICLES_EROSION_DISPLAY_FLAG			0x01 << 8



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

struct CellData
{
	CellData() {};

	glm::vec2 cell;
	std::pair<float, float> triangleMass; //first - ABC, second - AABC
	
	Triangle ABCbeforeUpdate;
	Triangle AABCbeforeUpdate;

	Triangle ABCafterUpdate;
	Triangle AABCafterUpdate;
};

struct ParticleData
{
	FluidParticleLight particle;
	std::vector<int> dCs;
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
	void SedimentationParticleInteraction(int ij, float dC, const FluidParticle& fp1, const FluidParticle& fp2)
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
		m_FluidParticlesAfterSedimentation[fp1.Id].dCs.push_back(ij);
	}

	void DepositionParticleInteraction(int ij, float dC, const FluidParticle& fp, const FluidParticle& bp)
	{
		FluidParticleLight p1, p2;
		p1.id = fp.Id;
		p1.sedim = fp.sedim;
		p1.sedim_delta = fp.sedim_delta;
		p1.sedim_ratio = fp.sedim_ratio;

		p2.id = bp.Id;
		p2.dM = bp.dM;

		m_PostDeposition[ij] = { ij, dC, p1, p2 };
		m_FluidParticlesAfterDeposition[fp.Id].dCs.push_back(ij);
	}

	void SedimentFlowSphBoundaryInteraction(int ij, float dC, const FluidParticle& fp, const FluidParticle& bp)
	{
		FluidParticleLight p1, p2;
		p1.id = fp.Id;
		p1.sedim = fp.sedim;
		p1.sedim_delta = fp.sedim_delta;
		p1.sedim_ratio = fp.sedim_ratio;

		p2.id = bp.Id;
		p2.dM = bp.dM;

		m_PostSedimentFlowSphBoundary[ij] = { ij, dC, p1, p2 };
		m_FluidParticlesAfterSedimentFlowSphBoundary[fp.Id].dCs.push_back(ij);
		m_BoundaryParticlesAfterSedimentFlow[bp.Id].dCs.push_back(ij);
	}
	
	void SedimentFlowSphSphInteraction(int ij, float dC, const FluidParticle& fp1, const FluidParticle& fp2)
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
		m_FluidParticlesAfterSedimentFlowSphSph[fp1.Id].dCs.push_back(ij);
	}

	void InsertPostErosionInteraction(int ij, float dC, const FluidParticle& fp, const FluidParticle& bp, float m)
	{
		FluidParticleLight p1, p2;
		p1.id = fp.Id;
		p1.sedim = fp.sedim;
		p1.sedim_delta = fp.sedim_delta;
		p1.sedim_ratio = fp.sedim_ratio;
		p1.dM = m;

		p2.id = bp.Id;
		p2.dM = bp.dM;
		p2.triangle = bp.triangle;

		m_PostErosion_dCs[ij] = { ij, dC, p1, p2 };
		m_FluidParticlesAfterErosion[fp.Id].dCs.push_back(ij);
		m_BoundaryParticlesAfterErosion[bp.Id].dCs.push_back(ij);
	}

	void CellsClear()
	{
		m_Cells.clear();
	}

	void PushBackMass(const glm::vec2& cell, const std::pair<float, float>& triangleMass)
	{
		int cellIdx = cell[0] * 1000 + cell[1];
		m_Cells[cellIdx].cell = cell;
		m_Cells[cellIdx].triangleMass = triangleMass;
	}

	void ABCGeometryBeforeUpdate(const glm::vec2& cell, const Triangle& abc)
	{
		int cellIdx = cell[0] * 1000 + cell[1];
		m_Cells[cellIdx].ABCbeforeUpdate = abc;
	}

	void AABCGeometryBeforeUpdate(const glm::vec2& cell, const Triangle& aabc)
	{
		int cellIdx = cell[0] * 1000 + cell[1];
		m_Cells[cellIdx].AABCbeforeUpdate = aabc;
	}

	void ABCGeometryAfterUpdate(const glm::vec2& cell, const Triangle& abc)
	{
		int cellIdx = cell[0] * 1000 + cell[1];
		m_Cells[cellIdx].ABCafterUpdate = abc;
	}

	void AABCGeometryAfterUpdate(const glm::vec2& cell, const Triangle& aabc)
	{
		int cellIdx = cell[0] * 1000 + cell[1];
		m_Cells[cellIdx].AABCafterUpdate = aabc;
	}


	void InsertFluidParticleAfterSedimentation(const FluidParticle& fp)
	{
		FluidParticleLight p;
		p.id = fp.Id;
		p.sedim = fp.sedim;
		p.sedim_delta = fp.sedim_delta;
		p.sedim_ratio = fp.sedim_ratio;

		m_FluidParticlesAfterSedimentation[fp.Id].particle = p;
	}

	void InsertFluidParticleAfterDeposition(const FluidParticle& fp)
	{
		FluidParticleLight p;
		p.id = fp.Id;
		p.sedim = fp.sedim;
		p.sedim_delta = fp.sedim_delta;
		p.sedim_ratio = fp.sedim_ratio;

		m_FluidParticlesAfterDeposition[fp.Id].particle = p;
	}
	
	void InsertFluidParticleAfterSedimentFlowSphSph(const FluidParticle& fp)
	{
		FluidParticleLight p;
		p.id = fp.Id;
		p.sedim = fp.sedim;
		p.sedim_delta = fp.sedim_delta;
		p.sedim_ratio = fp.sedim_ratio;

		m_FluidParticlesAfterSedimentFlowSphSph[fp.Id].particle = p;
	}

	void InsertFluidParticleAfterSedimentFlowSphBoundary(const FluidParticle& fp)
	{
		FluidParticleLight p;
		p.id = fp.Id;
		p.sedim = fp.sedim;
		p.sedim_delta = fp.sedim_delta;
		p.sedim_ratio = fp.sedim_ratio;

		m_FluidParticlesAfterSedimentFlowSphBoundary[fp.Id].particle = p;
	}
	
	void InsertBoundaryParticleAfterSedimentFlow(const FluidParticle& bp)
	{
		FluidParticleLight p;
		p.id = bp.Id;
		p.dM = bp.dM;

		m_BoundaryParticlesAfterSedimentFlow[bp.Id].particle = p;
	}

	void InsertFluidParticleAfterErosion(const FluidParticle& fp)
	{
		FluidParticleLight p;
		p.id = fp.Id;
		p.sedim = fp.sedim;
		p.sedim_delta = fp.sedim_delta;
		p.sedim_ratio = fp.sedim_ratio;

		m_FluidParticlesAfterErosion[fp.Id].particle = p;
	}

	void InsertBoundaryParticleAfterErosion(const FluidParticle& bp)
	{
		FluidParticleLight p;
		p.id = bp.Id;
		p.dM = bp.dM;
		p.triangle = bp.triangle;

		m_BoundaryParticlesAfterErosion[bp.Id].particle = p;
	}

	void DisplayDebugWindow(uint16_t display_flags, bool& pause) const
	{
		ImGui::Begin("Erosion Debugger");
		if (ImGui::Button("Next frame"))
			pause = false;

		ImGui::Text("--Pairs Interaction--");

		if (SEDIMENTATION_DISPLAY_FLAG & display_flags)
			displaySedimentation();
		if (DEPOSITION_DISPLAY_FLAG & display_flags)
			displayDeposition();
		if (SEDIMENT_FLOW_DISPLAY_FLAG & display_flags)
			displaySedimentFlow();
		if (EROSION_DISPLAY_FLAG & display_flags)
			displayErosion();
		if (CELLS_DISPLAY_FLAG & display_flags)
			displayCells();

		ImGui::Text("--Particles--");
		if (PARTICLES_SEDIMENTATION_DISPLAY_FLAG & display_flags)
			displayParticlesSedimentation();
		if (PARTICLES_DEPOSITION_DISPLAY_FLAG & display_flags)
			displayParticlesDeposition();
		if (PARTICLES_SEDIMENT_FLOW_DISPLAY_FLAG & display_flags)
			displayParticlesSedimentFlow();
		if (PARTICLES_EROSION_DISPLAY_FLAG & display_flags)
			displayParticlesErosion();

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

		//todo
		m_FluidParticlesAfterSedimentation.clear();
	}


	void PostDepositionBufferInit(int size)
	{
		m_PostDeposition.clear();
		m_PostDeposition.resize(size);

		//todo
		m_FluidParticlesAfterDeposition.clear();
	}

	void PostSedimentFlowSphSphInit(int size)
	{
		m_PostSedimentFlowSphSph.clear();
		m_PostSedimentFlowSphSph.resize(size);

		//todo
		m_FluidParticlesAfterSedimentFlowSphSph.clear();
	}
	
	void PostSedimentFlowSphBoundaryInit(int size)
	{
		m_PostSedimentFlowSphBoundary.clear();
		m_PostSedimentFlowSphBoundary.resize(size);

		//todo
		m_FluidParticlesAfterSedimentFlowSphBoundary.clear();
		m_BoundaryParticlesAfterSedimentFlow.clear();
	}

	void PostErosionInit(int size)
	{
		m_PostErosion_dCs.clear();
		m_PostErosion_dCs.resize(size);

		//todo
		m_FluidParticlesAfterErosion.clear();
		m_BoundaryParticlesAfterErosion.clear();
	}

	//std::vector<PPInteraction>& GetPostSedimentationBuffer()
	//{
	//	return m_PostSedimentation;
	//}


private:
	static Debugger* m_Instance;
	

	std::vector<PPInteraction> m_PostSedimentation;
	std::vector<PPInteraction> m_PostDeposition;
	std::vector<PPInteraction> m_PostErosion_dCs;
	std::vector<PPInteraction> m_PostSedimentFlowSphBoundary;
	std::vector<PPInteraction> m_PostSedimentFlowSphSph;

	std::map<int, CellData> m_Cells;

	std::map<unsigned long, ParticleData> m_FluidParticlesAfterSedimentation;
	std::map<unsigned long, ParticleData> m_FluidParticlesAfterDeposition;

	std::map<unsigned long, ParticleData> m_FluidParticlesAfterSedimentFlowSphSph;
	std::map<unsigned long, ParticleData> m_FluidParticlesAfterSedimentFlowSphBoundary;
	std::map<unsigned long, ParticleData> m_BoundaryParticlesAfterSedimentFlow;
	
	std::map<unsigned long, ParticleData> m_FluidParticlesAfterErosion;
	std::map<unsigned long, ParticleData> m_BoundaryParticlesAfterErosion;

	void displayParticlesSedimentation() const
	{
		if (!ImGui::CollapsingHeader("After SEDIMENTATION"))
			return;

		for (const auto& part : m_FluidParticlesAfterSedimentation)
		{
			if (ImGui::TreeNode((void*)(intptr_t)part.first, "(%d)", part.first))
			{
				ImGui::Text("id: %d", part.second.particle.id);
				ImGui::Text("sedim: %.15f", part.second.particle.sedim);
				ImGui::Text("sedim_delta: %.15f", part.second.particle.sedim_delta);
				ImGui::Text("sedim_ratio: %f", part.second.particle.sedim_ratio);
				
				for (int i = 0; i < part.second.dCs.size(); i++)
				{
					const auto& dC_idx = part.second.dCs[i];
					const auto& data = m_PostSedimentation[dC_idx];
					const char* format;
					if (data.dC == -777.0f)
						format = "dC[%d] (skipped)";
					else
						format = "dC[%d]";
					if (ImGui::TreeNode((void*)(intptr_t)i, format, dC_idx))
					{
						
						ImGui::Text("ij: %d", data.ij);
						ImGui::Text("[dC]: %.15f", data.dC);
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
		}
	}

	void displayParticlesDeposition() const
	{
		if (!ImGui::CollapsingHeader("After DEPOSITION"))
			return;

		for (const auto& part : m_FluidParticlesAfterDeposition)
		{
			if (ImGui::TreeNode((void*)(intptr_t)part.first, "(%d)", part.first))
			{
				ImGui::Text("id: %d", part.second.particle.id);
				ImGui::Text("sedim: %.15f", part.second.particle.sedim);
				ImGui::Text("sedim_delta: %.15f", part.second.particle.sedim_delta);
				ImGui::Text("sedim_ratio: %f", part.second.particle.sedim_ratio);
				
				for (int i = 0; i < part.second.dCs.size(); i++)
				{
					const auto& dC_idx = part.second.dCs[i];
					const auto& data = m_PostDeposition[dC_idx];
					const char* format;
					if (data.dC == -777.0f)
						format = "dC_BP[%d] (skipped) (sedim)";
					else if (data.dC == -666.0f)
						format = "dC_BP[%d] (skipped) (settling)";
					else
						format = "dC_BP[%d]";
					if (ImGui::TreeNode((void*)(intptr_t)i, format, dC_idx))
					{
						ImGui::Text("ij: %d", data.ij);
						ImGui::Text("[dC_BP]: %.15f (%e)", data.dC, data.dC);
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
	}

	void displayParticlesSedimentFlow() const
	{
		if (!ImGui::CollapsingHeader("After SEDIMENT FLOW"))
			return;

		if (ImGui::TreeNode("sph-sph"))
		{
			for (const auto& pair: m_FluidParticlesAfterSedimentFlowSphSph)
			{
				if (ImGui::TreeNode((void*)(intptr_t)pair.first, "(%d)", pair.first))
				{
					ImGui::Text("id: %d", pair.second.particle.id);
					ImGui::Text("sedim: %.15f", pair.second.particle.sedim);
					ImGui::Text("sedim_delta: %.15f", pair.second.particle.sedim_delta);
					ImGui::Text("sedim_ratio: %f", pair.second.particle.sedim_ratio);

					for (int i = 0; i < pair.second.dCs.size(); i++)
					{
						int dC_idx = pair.second.dCs[i];
						const auto& data = m_PostSedimentFlowSphSph[dC_idx];
						const char* format;
						if (data.dC == -777.0f)
							format = "dC[%d] (skipped)";
						else
							format = "dC[%d]";

						if (ImGui::TreeNode((void*)(intptr_t)i, format, i))
						{
							ImGui::Text("ij: %d", data.ij);
							ImGui::Text("[[dC]]: %.15f", data.dC);
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
			}
			ImGui::TreePop();
		}
		
		if (ImGui::TreeNode("sph-boundary"))
		{
			if (ImGui::TreeNode("sph"))
			{
				for (const auto& pair : m_FluidParticlesAfterSedimentFlowSphBoundary)
				{
					if (ImGui::TreeNode((void*)(intptr_t)pair.first, "(%d)", pair.first))
					{
						ImGui::Text("id: %d", pair.second.particle.id);
						ImGui::Text("sedim: %.15f", pair.second.particle.sedim);
						ImGui::Text("sedim_delta: %.15f", pair.second.particle.sedim_delta);
						ImGui::Text("sedim_ratio: %f", pair.second.particle.sedim_ratio);

						for (int i = 0; i < pair.second.dCs.size(); i++)
						{
							int dC_idx = pair.second.dCs[i];
							const auto& data = m_PostSedimentFlowSphBoundary[dC_idx];
							const char* format;
							if (data.dC == -777.0f)
								format = "dC_BP[%d] (skipped)";
							else
								format = "dC_BP[%d]";
							if (ImGui::TreeNode((void*)(intptr_t)i, format, dC_idx))
							{
								ImGui::Text("ij: %d", data.ij);
								ImGui::Text("[[dC_BP]]: %.15f (%e)", data.dC, data.dC);
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
				ImGui::TreePop();
			}

			if (ImGui::TreeNode("boundary"))
			{
				ImGui::Text("--Total-- %d", m_BoundaryParticlesAfterSedimentFlow.size());
				for (const auto& pair : m_BoundaryParticlesAfterSedimentFlow)
				{
					if (ImGui::TreeNode((void*)(intptr_t)pair.first, "(%d)", pair.first))
					{
						ImGui::Text("id: %d", pair.second.particle.id);
						ImGui::Text("dM: %.15f", pair.second.particle.dM);

						for (int i = 0; i < pair.second.dCs.size(); i++)
						{
							int dC_idx = pair.second.dCs[i];
							const auto& data = m_PostSedimentFlowSphBoundary[dC_idx];
							const char* format;
							if (data.dC == -777.0f)
								format = "dC_BP[%d] (skipped)";
							else
								format = "dC_BP[%d]";
							if (ImGui::TreeNode((void*)(intptr_t)i, format, dC_idx))
							{
								ImGui::Text("ij: %d", data.ij);
								ImGui::Text("[[dC_BP]]: %.15f (%e)", data.dC, data.dC);
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
				ImGui::TreePop();
			}
			ImGui::TreePop();
		}
	}

	void displayParticlesErosion() const
	{
		if (!ImGui::CollapsingHeader("After EROSION"))
			return;

		if (ImGui::TreeNode("sph"))
		{
			for (const auto& pair : m_FluidParticlesAfterErosion)
			{
				if (ImGui::TreeNode((void*)(intptr_t)pair.first, "(%d)", pair.first))
				{
					ImGui::Text("id: %d", pair.second.particle.id);
					ImGui::Text("sedim: %.15f", pair.second.particle.sedim);
					ImGui::Text("sedim_delta: %.15f", pair.second.particle.sedim_delta);
					ImGui::Text("sedim_ratio: %f", pair.second.particle.sedim_ratio);
					for (int i = 0; i < pair.second.dCs.size(); i++)
					{
						int dC_idx = pair.second.dCs[i];
						const auto& data = m_PostErosion_dCs[dC_idx];
						const char* format;
						if (data.dC == -666.0f)
							format = "dC_BP[%d] (skipped: minVrel > vRel)";
						else
							format = "dC_BP[%d]";
						if (ImGui::TreeNode((void*)(intptr_t)i, format, dC_idx))
						{
							ImGui::Text("ij: %d", data.ij);
							ImGui::Text("[[dC_BP]]: %.15f (%e)", data.dC, data.dC);
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
					ImGui::TreePop();
				}
			}
			ImGui::TreePop();
		}
		
		if (ImGui::TreeNode("boundary"))
		{
			ImGui::Text("--Total-- %d", m_BoundaryParticlesAfterErosion.size());
			for (const auto& pair : m_BoundaryParticlesAfterErosion)
			{
				if (ImGui::TreeNode((void*)(intptr_t)pair.first, "(%d)", pair.first))
				{
					ImGui::Text("id: %d", pair.second.particle.id);
					ImGui::Text("dM: %.15f", pair.second.particle.dM);
					ImGui::Text("triangle: %c", pair.second.particle.triangle);

					for (int i = 0; i < pair.second.dCs.size(); i++)
					{
						int dC_idx = pair.second.dCs[i];
						const auto& data = m_PostErosion_dCs[dC_idx];
						const char* format;
						if (data.dC == -666.0f)
							format = "dC_BP[%d] (skipped: minVrel > vRel)";
						else
							format = "dC_BP[%d]";
						if (ImGui::TreeNode((void*)(intptr_t)i, format, dC_idx))
						{
							ImGui::Text("ij: %d", data.ij);
							ImGui::Text("[[dC_BP]]: %.15f (%e)", data.dC, data.dC);
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
					ImGui::TreePop();
				}
			}
			ImGui::TreePop();
		}
	}

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
			else if (data.dC == -666.0f)
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
					format = "dC_BP[%d] (skipped)";
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

			ImGui::TreePop();
		}
	}

	void displayErosion() const
	{
		if (!ImGui::CollapsingHeader("Post EROSION"))
			return;

		ImGui::Text("--Total %d--", m_PostErosion_dCs.size());
		for (int i = 0; i < m_PostErosion_dCs.size(); i++)
		{
			const auto& data = m_PostErosion_dCs[i];
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
				ImGui::Text("\ttriangle: %f", data.particle2.triangle);
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

	void displayCells() const
	{
		if (!ImGui::CollapsingHeader("CELLS"))
			return;


		for (const auto& cell : m_Cells)
		{
			int cellIdx = cell.first;
			const CellData& cd = cell.second;
			if (ImGui::TreeNode((void*)(intptr_t)cellIdx, "(%d, %d)", (int)cd.cell[0], (int)cd.cell[1]))
			{
				ImGui::Text("ABC mass: %.15f", cd.triangleMass.first);
				ImGui::Text("AABC mass: %.15f", cd.triangleMass.second);
				if (ImGui::TreeNode("ABC before-after"))
				{
					ImGui::Text("A(%.2f, %.2f): %.15f", cd.ABCbeforeUpdate.C[0], cd.ABCbeforeUpdate.C[2], cd.ABCbeforeUpdate.C[1]);
					ImGui::Text("A(%.2f, %.2f): %.15f", cd.ABCafterUpdate.C[0],  cd.ABCafterUpdate.C[2],  cd.ABCafterUpdate.C[1]);
					ImGui::Text("B(%.2f, %.2f): %.15f", cd.ABCbeforeUpdate.B[0], cd.ABCbeforeUpdate.B[2], cd.ABCbeforeUpdate.B[1]);
					ImGui::Text("B(%.2f, %.2f): %.15f", cd.ABCafterUpdate.B[0],  cd.ABCafterUpdate.B[2],  cd.ABCafterUpdate.B[1]);
					ImGui::Text("C(%.2f, %.2f): %.15f", cd.ABCbeforeUpdate.A[0], cd.ABCbeforeUpdate.A[2], cd.ABCbeforeUpdate.A[1]);
					ImGui::Text("C(%.2f, %.2f): %.15f", cd.ABCafterUpdate.A[0],  cd.ABCafterUpdate.A[2],  cd.ABCafterUpdate.A[1]);

					ImGui::TreePop();
				}

				if (ImGui::TreeNode("AABC before-after"))
				{
					ImGui::Text("AA(%.2f, %.2f): %.15f", cd.AABCbeforeUpdate.A[0], cd.AABCbeforeUpdate.A[2], cd.AABCbeforeUpdate.A[1]);
					ImGui::Text("AA(%.2f, %.2f): %.15f", cd.AABCafterUpdate.A[0],  cd.AABCafterUpdate.A[2],  cd.AABCafterUpdate.A[1]);
					ImGui::Text("B (%.2f, %.2f): %.15f", cd.AABCbeforeUpdate.B[0], cd.AABCbeforeUpdate.B[2], cd.AABCbeforeUpdate.B[1]);
					ImGui::Text("B (%.2f, %.2f): %.15f", cd.AABCafterUpdate.B[0],  cd.AABCafterUpdate.B[2],  cd.AABCafterUpdate.B[1]);
					ImGui::Text("C (%.2f, %.2f): %.15f", cd.AABCbeforeUpdate.C[0], cd.AABCbeforeUpdate.C[2], cd.AABCbeforeUpdate.C[1]);
					ImGui::Text("C (%.2f, %.2f): %.15f", cd.AABCafterUpdate.C[0],  cd.AABCafterUpdate.C[2],  cd.AABCafterUpdate.C[1]);

					ImGui::TreePop();
				}

				ImGui::TreePop();
			}
		}
	}

	Debugger()
	{

	}
};

Debugger* Debugger::m_Instance = nullptr;