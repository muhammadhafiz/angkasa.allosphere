//#version 120
uniform sampler2D texture0;
uniform sampler2D texture1;

void main() {
    gl_TexCoord[0] = gl_MultiTexCoord0;
    vec4 textureColor = texture2D(texture0, gl_TexCoord[0].st);
    vec4 textureColor1 = texture2D(texture1, gl_TexCoord[0].st);

    vec4 vertex = gl_ModelViewMatrix * gl_Vertex;
//    gl_PointSize = 0.01 + (length(textureColor) * 5.0);
//    gl_PointSize = 0.01 + (textureColor.r * 5.0);
    gl_PointSize = 0.0001+ (textureColor.r * 20.0);

    #if defined(OMNI)
      gl_Position = omni_render(vertex);
    #else
      gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    #endif
}
