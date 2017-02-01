varying vec4 color;
varying vec3 normal, lightDir, eyeVec;

//uniform sampler2D texture0;
uniform sampler2D texture0;
uniform float lighting;
uniform float texture;
uniform vec4 ambient, diffuse;
uniform vec3 meterVal;

void main() {
    vec4 colorMixed = color;
    //if (texture > 0.0) {
      vec4 textureColor = texture2D(texture0, gl_TexCoord[0].st);
      colorMixed = mix(color, textureColor, texture);
    //}

    // colors
    vec4 c_ambient = colorMixed * 0.2;
    vec4 c_diffuse = colorMixed * 0.6;
    vec4 c_specular = vec4(1.0);

    // calc
    vec3 N = normalize(normal);
    vec3 L = lightDir;
    float lambert = max(dot(N, L), 0.0);
    vec3 E = eyeVec;
    vec3 R = reflect(L, N);
    // float shininess = 0.9 + 1e-20; // [?] why this value?
    float shininess = 3.0; // [?] why this value?
    float spec = pow(max(dot(R, E), 0.0), shininess);

    // result
    vec4 ambient = colorMixed;
    vec4 diffuse = c_diffuse * lambert;
    vec4 specular = c_specular * spec * lambert;
    vec4 final_color = ambient + diffuse + specular;
    // added this clamp otherwise colors get wack
    final_color = clamp(final_color, 0.0, 1.0);
    final_color = mix(colorMixed, final_color, lighting);
    // final_color.a = color.a;

    gl_FragColor = textureColor;// vec4(1,1,1,0.1); //final_color;textureColor; 

}
