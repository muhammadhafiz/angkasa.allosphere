//#version 120
uniform sampler2D texture0;
uniform sampler2D texture1;

void main() {
    gl_TexCoord[0] = gl_MultiTexCoord0;
    vec4 textureColor = texture2D(texture0, gl_TexCoord[0].st);
    vec4 textureColor1 = texture2D(texture1, gl_TexCoord[0].st);
    vec4 vertex = gl_ModelViewMatrix * gl_Vertex;

// ALLOSPHERE:
    gl_PointSize = (textureColor.r * 5000.0);

    #if defined(OMNI)
      gl_Position = omni_render(vertex);
    #else
      gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    #endif
}
