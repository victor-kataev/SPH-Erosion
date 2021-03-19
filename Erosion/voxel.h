enum VoxelType
{
	VOXEL_AIR = 0,
	VOXEL_WAT,
	VOXEL_MAT
};

struct Voxel
{
	float density;
	glm::vec3 position;
	glm::vec3 velocity;
	VoxelType type;

	Voxel()
	{
		density = 0;
		position = glm::vec3(0.0);
		velocity = glm::vec3(0.0);
		type = VoxelType::VOXEL_AIR;
	}

	Voxel(VoxelType t)
	{
		type = t;
	}
};