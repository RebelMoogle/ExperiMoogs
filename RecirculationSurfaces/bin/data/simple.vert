#version 330

uniform mat4 modelViewMatrix;
uniform mat4 modelViewProjectionMatrix;
uniform vec4 shaderColor;
uniform vec3 camDirection;
uniform vec3 camPosition;
in vec4 position;
in vec4 normal;
in vec2 texcoord;

out vec4 fragNormal;
out vec2 fragTexcoord;
out vec4 fragColor;

out vec3 lightDir;
out vec3 viewDir;
out vec3 reflectDir;


void main(){

	gl_Position = modelViewProjectionMatrix * position;
	fragTexcoord = texcoord;
	fragNormal = normal;

	lightDir = gl_Position.xyz - (modelViewProjectionMatrix*vec4(camPosition, 1.0)).xyz; //vec3(1., 1., 0.); 
	reflectDir = reflect(lightDir, normalize(normal.xyz));
	viewDir = camDirection;

}