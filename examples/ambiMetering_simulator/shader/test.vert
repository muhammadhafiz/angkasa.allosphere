#version 120

varying vec4 color;
varying vec3 normal, lightDir, eyeVec;
uniform sampler2D texture0;


void main() {
    vec4 textureColor = texture2D(texture0, gl_TexCoord[0].st);

    color = gl_Color;
    //color = vec4(1,1,1,1);
    vec4 vertex = gl_ModelViewMatrix * gl_Vertex;
    normal = gl_NormalMatrix * gl_Normal;
    vec3 V = vertex.xyz;
    eyeVec = normalize(-V);
    lightDir = normalize(vec3(gl_LightSource[0].position.xyz - V));
    gl_TexCoord[0] = gl_MultiTexCoord0;

    //gl_PointSize = abs(gl_Vertex.x* 20.0);
    //gl_PointSize = length(textureColor) * 5.0;
    gl_PointSize = textureColor.r * 10.0;
    //gl_PointSize = 5.0;

    #if defined(OMNI)
      gl_Position = omni_render(vertex);
    #else
      gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    #endif

}
