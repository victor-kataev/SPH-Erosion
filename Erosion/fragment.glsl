#version 330 core

in vec3 Normal;
in vec3 FragPos;

out vec4 color;

struct Material
{
	vec3 ka;
	vec3 kd;
	vec3 ks;
	float ksh;
};

struct Light
{
	vec3 dir;
	vec3 color;
};

uniform vec3 myColor;

uniform vec3 viewerPos;
uniform Material material;
uniform Light dirLight;

void main()
{
	vec3 lightDir = normalize(-dirLight.dir);
	
	float diffuseImpact = max(dot(lightDir, Normal), 0.0);

	vec3 reflected = reflect(-lightDir, Normal);
	vec3 viewerDir = normalize(viewerPos - FragPos);
	float specularImpact = pow(max(dot(reflected, viewerDir), 0.0), material.ksh);

	vec3 ambient = dirLight.color * material.ka * myColor;
	vec3 diffuse = dirLight.color * material.kd * diffuseImpact * myColor;
	vec3 specular = dirLight.color * material.ks * specularImpact * myColor;

	vec3 result = ambient + diffuse + specular;

	color = vec4(result, 1.0);
}