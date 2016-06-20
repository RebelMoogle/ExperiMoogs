#version 330

uniform mat4 modelViewProjectionMatrix;
uniform vec4 shaderColor;
in vec4 position;
in vec4 color;

out vec4 fragColor;



void main(){

	gl_Position = modelViewProjectionMatrix * position;
	fragColor = color;

}