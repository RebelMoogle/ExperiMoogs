#version 330

uniform vec4 shaderColor;
in vec3 fragNormal;
in vec2 fragTexcoord;
in vec3 lightDir;
in vec3 viewDir;
in vec3 reflectDir;
out vec4 outputColor;

void main()
{
// illuminated streamlines.

	//outputColor = vec4(shaderColor.rgb, 1.0) * vec4(normalize(fragNormal.xyz), 1);
	//vec4(fragTexcoord.x, 0.5, 0.5, 1.0- fragTexcoord.x);
	// more color stages??#
	//outputColor = vec4(fragTexcoord.y, 1.0-fragTexcoord.y, 0.0, 1.0);

	float lambertian = max(dot(normalize(lightDir), normalize(fragNormal.xyz)), 0.0);
	float specular = 0.0;

	float specAngle = max(dot(reflectDir, viewDir), 0.0);
	specular = pow(specAngle, 4.0);
  
  	outputColor = vec4(0.1,0.1,0.1,0.0) + vec4(lambertian*shaderColor.rgb + specular, 1.0);
}