

#include <allocore/al_Allocore.hpp>
#include <allocore/types/al_SingleRWRingBuffer.hpp>
// #include <alloutil/al_OmniStereoGraphicsRenderer.hpp>
#include "../OmniRender.hpp"
#include "Cuttlebone/Cuttlebone.hpp"
#include "Gamma/scl.h"

#include "../soundfilebuffered.hpp"
#include "allospherespeakerlayouts.h"

#include "../state.hpp"

#include "allocore/graphics/al_Shader.hpp"
#include "alloutil/al_ShaderManager.hpp"

// #define AUDIO_BLOCK_SIZE 512

using namespace al;
using namespace std;

//#define SIN gam::scl::sinT9
//#define COS gam::scl::cosT8

#define SIN sin
#define COS cos

class MeterParams {
public:
	MeterParams() :
	    mSphere(true),
	    mPeak(false),
	    mDecibel(false),
	    mGrayscale(false),
	    mPeakGain(3.0),
	    mRmsGain(20.0),
	    mDbRange(80.0)
	{}

	bool mSphere;
	bool mPeak;
	bool mDecibel;
	bool mGrayscale;


	float mPeakGain;
	float mRmsGain;
	double mDbRange;
};

class AmbiMeterRenderer : public OmniStereoGraphicsRenderer
{
public:

	AmbiMeterRenderer() :
	    mMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
	    mTextureBuffer(((SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float) ),
	    mSpriteTex(16,16, Graphics::LUMINANCE, Graphics::FLOAT)
	{
		addSphereWithTexcoords(mMesh, 1, SPATIAL_SAMPLING);
		mMesh.primitive(Graphics::POINTS);

		addSphereWithTexcoords(mWireframeMesh, 1, SPATIAL_SAMPLING);
		mWireframeMesh.primitive(Graphics::LINES);

		gaussianSprite(mSpriteTex);

//		mMesh.vertex(-10,-10,-10); mMesh.color(RGB(1,0,0)); mMesh.texCoord(0, 0);
//	    mMesh.vertex( 10,-10,-10); mMesh.color(RGB(1,1,0)); mMesh.texCoord(1, 0);
//	    mMesh.vertex( 0, 10,-10); mMesh.color(RGB(0,1,1)); mMesh.texCoord(0.5, 1);

//		lens().near(0.1).far(25).fovy(45);
//		nav().pos(0,0,4);
//		nav().quat().fromAxisAngle(0. * M_2PI, 0, 1, 0);

		// mTexture.create();
		// mTexture.resize(32, 32);

		initWindow(Window::Dim(0,0, 600, 400), "Angkasa AlloSphere Renderer", VIS_FPS);
		mStateTaker.start();

//		light.ambient(Color(0.4, 0.4, 0.4, 1.0));
//	    light.pos(5, 5, 5);
		params.mSphere = true;
	}

	// Create a Gaussian "bump" function to use for the sprite
	void gaussianSprite(Texture &t){
		int Nx = t.width();
		int Ny = t.height();
		float * pixels = t.data<float>();

		for(int j=0; j<Ny; ++j)
		{
			float y = float(j)/(Ny-1)*2-1;
			for(int i=0; i<Nx; ++i)
			{
				float x = float(i)/(Nx-1)*2-1;
				float m = exp(-3*(x*x + y*y));
				pixels[j*Nx + i] = m;
			}
		}
	}

	// Load shaders whenever shaderManager polls
	virtual void loadShaders(){
		shaderManager.vertLibCode = "#version 120\n";
		shaderManager.vertLibCode.append("#define OMNI 1\n" + mOmni.glsl());
		shaderManager.addShaderFile("al_point", "al_point.vert", "al_point.frag");
	}

	virtual void onAnimate(al_sec dt) override {
		static bool initd = false;
		if(shaderManager.poll() || !initd){
			loadShaders();
			initd = true;
		}
		// light.pos(nav().pos());
		pose = nav();
		// std::cout << dt << std::endl;
	}

	virtual void onDraw(Graphics &g) override {

		State state;
		bool newState = false;
		while (mStateTaker.get(state) > 0) { newState = true;} // Pop all states in queue

		if(newState) {
			memcpy(mMeterValues.data.ptr, state.rmsTexture,
			       SPATIAL_SAMPLING *SPATIAL_SAMPLING * 3 * sizeof(float));
			mTexture.submit(mMeterValues, true);

			mTexture.filter(Texture::LINEAR);
			// std::cout << "New State" << std::endl;
		}

//		g.polygonMode(Graphics::LINE);
//		light();
		shader().begin();
		omni().uniforms(shader());
	    // omni().clearColor() = Color(0.4);

		shader().uniform("lighting", 0.0);
	    shader().uniform("texture", 0.0);

		g.blendAdd();
		// Draw wireframe mesh
		g.lineWidth(0.25);
		g.color(1, 1, 1, 0.5);
		g.draw(mWireframeMesh);
		shader().end();

		// // Read the texture buffer
		// int bytesToRead = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
		// while(mTextureBuffer.readSpace() >= bytesToRead) {
		// 	mTextureBuffer.read(mMeterValues.data.ptr, bytesToRead);
		// 	mTexture.submit(mMeterValues, true);
		// 	mTexture.filter(Texture::LINEAR);
		// }

		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glEnable(GL_POINT_SMOOTH);
		glEnable(GL_POINT_SPRITE);
		glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);

		ShaderProgram* s = shaderManager.get("al_point");
		s->begin();
		omni().uniforms(*s);
		s->uniform("texture0", 0);
		s->uniform("texture1", 1);
		mTexture.bind(0); //binds to textureID
		mSpriteTex.bind(1);
		g.draw(mMesh);
		mSpriteTex.unbind(1);
		mTexture.unbind(0);
		s->end();

		glDisable(GL_POINT_SPRITE);


	}

	virtual void start() {
      // mStateTaker.start();                        // non-blocking
      OmniStereoGraphicsRenderer::start();  // blocks
    }

//	virtual void onKeyDown(const Keyboard& k) override {
//		if (k.key() == 'g') {
//			params.mGrayscale = !params.mGrayscale;
//		} else if (k.key() == 'h') {
//			params.mDecibel = !params.mDecibel;
//		} else if (k.key() == 'j') {
//			params.mSphere = !params.mSphere;
//		} else if (k.key() == 'k') {
//			params.mPeak = !params.mPeak;
//		}
//	}

private:
	Mesh mMesh;
	Mesh mWireframeMesh;
	Texture mTexture;
	Texture mSpriteTex;
	Array mMeterValues;
	Array mNewRmsMeterValues;
	SingleRWRingBuffer mTextureBuffer;
	ShaderManager shaderManager;

	// Light light;

	cuttlebone::Taker<State> mStateTaker;

	MeterParams params;
};


int main(int argc, char *argv[])
{
	AmbiMeterRenderer app;
	app.start();
}
