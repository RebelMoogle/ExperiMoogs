#version 330

uniform mat4 modelViewMatrix;
uniform mat4 modelViewProjectionMatrix;
uniform vec4 shaderColor;
uniform vec3 camDirection;
uniform vec3 camPosition;
in vec4 position;
in vec4 normal;
in vec2 texcoord;

out vec3 fragNormal;
out vec2 fragTexcoord;
out vec4 fragColor;

out vec3 lightDir;
out vec3 viewDir;
out vec3 reflectDir;


void main(){

	gl_Position = modelViewProjectionMatrix * position;
	fragTexcoord = texcoord;
	//fragNormal = normal;

	lightDir = normalize(position.xyz - camPosition); //vec3(1., 1., 0.); 
	fragNormal = cross(cross(normal.xyz, lightDir), normal.xyz);
	reflectDir = reflect(lightDir, normalize(fragNormal));
	viewDir = camDirection;

}