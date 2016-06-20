#version 330 //point cloud shader

uniform vec4 shaderColor;
in vec4 fragColor;

out vec4 outputColor;

void main()
{
  
  	outputColor = vec4(length(fragColor.rgb), 0.0, 0.0, fragColor.a);
}