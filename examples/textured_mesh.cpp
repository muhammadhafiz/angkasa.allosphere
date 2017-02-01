// MAT201B
// Fall 2015
// Author(s): Tim Wood
//
// Shows how to:
// - Load an image into a texture
// - Bind a texture and draw on a mesh
//

#include "allocore/io/al_App.hpp"
using namespace al;
using namespace std;

struct AlloApp : App {
  Light light;
  Mesh m;
  Texture texture;

  AlloApp() {
    nav().pos(0, 0, 30);
    // light.pos(0, 0, 10);

    addSphereWithTexcoords(m);
    m.primitive(Graphics::POINTS);
    // m.generateNormals(); // don't need normals if we don't use lights

    initWindow();

    // load image into texture print out error and exit if failure
    Image image;
    SearchPaths searchPaths;
    searchPaths.addSearchPath(".");
    string filename = searchPaths.find("surround.jpg").filepath();
    if (image.load(filename)) {
      cout << "Read image from " << filename << endl;
    } else {
      cout << "Failed to read image from " << filename << "!!!" << endl;
      exit(-1);
    }
    texture.allocate(image.array());
  }

  virtual void onDraw(Graphics& g, const Viewpoint& v) {
    // light();  // disable light so texture shows up as is
    g.pushMatrix();
    g.scale(10);
    // g.rotate(180, 0, 0, 1);  // rotate sphere to show texture as we'd expect
    g.pointSize(20.0);
    texture.bind();
    g.draw(m);
    texture.unbind();
    g.popMatrix();
  }

  virtual void onAnimate(double dt) {}

  // onMouseMove is when the mouse moves
  virtual void onMouseMove(const ViewpointWindow& w, const Mouse& m) {
    // normalize mouse position from -1.0 to 1.0
    float x = float(m.x()) / w.width() * 2.f - 1.f;
    float y = (float(m.y()) / w.height() * 2.f - 1.f) * -1.f;

    // move light with mouse
    light.pos(Vec3f(x, y, 1.f) * 10.f);
  }

  // other mouse callbacks
  //
  virtual void onMouseDown(const ViewpointWindow& w, const Mouse& m) {}
  virtual void onMouseUp(const ViewpointWindow& w, const Mouse& m) {}
  virtual void onMouseDrag(const ViewpointWindow& w, const Mouse& m) {}
};
int main(int argc, char* argv[]) {
  AlloApp app;
  app.start();
  return 0;
}
