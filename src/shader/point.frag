uniform sampler2D texture0;
uniform sampler2D texture1;

void main() {
    vec4 textureColor = texture2D(texture0, gl_TexCoord[0].st);
    vec4 textureColor1 = texture2D(texture1, gl_TexCoord[0].st);

    //gl_FragColor = textureColor1; //vec4(1,1,1,0.1); //final_color;
    //gl_FragColor = (textureColor1*textureColor1)*5.0;
    gl_FragColor = (textureColor1*textureColor1*textureColor1) * 5.0;
}
