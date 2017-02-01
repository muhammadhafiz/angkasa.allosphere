#include "allocore/graphics/al_Shader.hpp"
#include "allocore/al_Allocore.hpp"
#include "alloutil/al_ShaderManager.hpp"

using namespace al;

class MyApp : public App{
public:
  ShaderManager sm;
	Mesh sphere;

	MyApp(){
    addSphere(sphere,1, 30, 30);
    sphere.primitive(Graphics::POINTS);
    for(int i=0; i<sphere.vertices().size(); i++){
      float f = float(i)/sphere.vertices().size();
      sphere.color(RGB(f,0,0));
    }

    nav().pos(0,0,10);
		initWindow();
	}

  void loadShaders(){
    sm.addShaderFile("point", "point.vert", "point.frag");
  }

  void onCreate(const ViewpointWindow& w){
    loadShaders();
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SMOOTH);
  }

	void onAnimate(double dt){
    if(sm.poll()){
      loadShaders();
    }
	}

	void onDraw(Graphics& g){
    ShaderProgram* s = sm.get("point");
    s->begin();
    g.draw(sphere);
    // g.pointSize(6.0);
    s->end();
	}
};

int main(){
	MyApp().start();
}
