

#include <allocore/al_Allocore.hpp>
#include <allocore/types/al_SingleRWRingBuffer.hpp>
#include <alloutil/al_OmniStereoGraphicsRenderer.hpp>
#include "Cuttlebone/Cuttlebone.hpp"
#include "Gamma/scl.h"

#include "soundfilebuffered.hpp"
#include "allospherespeakerlayouts.h"

#include "../ambiMetering_renderer/state.hpp"

#define AUDIO_BLOCK_SIZE 512

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
	    mTextureBuffer(((SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float) )
	{
		addSphereWithTexcoords(mMesh, 1, 128);
		addCube(mFrontalMesh, false, 0.1);

//		mMesh.vertex(-10,-10,-10); mMesh.color(RGB(1,0,0)); mMesh.texCoord(0, 0);
//	    mMesh.vertex( 10,-10,-10); mMesh.color(RGB(1,1,0)); mMesh.texCoord(1, 0);
//	    mMesh.vertex( 0, 10,-10); mMesh.color(RGB(0,1,1)); mMesh.texCoord(0.5, 1);

//		lens().near(0.1).far(25).fovy(45);
//		nav().pos(0,0,4);
//		nav().quat().fromAxisAngle(0. * M_2PI, 0, 1, 0);

		// mTexture.create();
		// mTexture.resize(32, 32);

		initWindow(Window::Dim(0,0, 600, 400), "Ambisonics Metering", 30);
		mStateTaker.start();

//		light.ambient(Color(0.4, 0.4, 0.4, 1.0));
//	    light.pos(5, 5, 5);
		params.mSphere = true;
	}

	virtual void onAnimate(al_sec dt) override {
		// light.pos(nav().pos());
		pose = nav();
		// std::cout << dt << std::endl;
	}

	virtual void onDraw(Graphics &g) override {

//		g.polygonMode(Graphics::LINE);
		State s;
//		light();
		omni().uniforms(shader());
	    omni().clearColor() = Color(0.4);

		shader().uniform("lighting", 0.1);
	    shader().uniform("texture", 0.9);

		bool newState = false;
		while (mStateTaker.get(s) > 0) { newState = true;} // Pop all states in queue

		if(newState) {
			memcpy(mMeterValues.data.ptr, s.rmsTexture,
			       SPATIAL_SAMPLING *SPATIAL_SAMPLING * 3 * sizeof(float));
			mTexture.submit(mMeterValues, true);

			mTexture.filter(Texture::LINEAR);
			// std::cout << "New State" << std::endl;
		}

		g.polygonMode(Graphics::FILL);
		g.antialiasing(Graphics::NICEST);
		g.depthTesting(true);
		g.blending(true);

//			mTexture.filter(Texture::LINEAR);
		if (params.mSphere) {
			mTexture.bind();
			g.draw(mMesh);
			mTexture.unbind();
			g.pushMatrix();
			g.translate(1.5, 0, 0);
			g.draw(mFrontalMesh);
			g.popMatrix();
		} else {
			mTexture.quad(g, 2, 2, -2, -1);
			g.pushMatrix();
			g.rotate(180, 0, 1, 0);
			mTexture.quad(g, 2, 2, -2, -1);
			g.popMatrix();
		}
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
	Mesh mFrontalMesh;
	Texture mTexture;
	Array mMeterValues;
	SingleRWRingBuffer mTextureBuffer;

	Light light;

	cuttlebone::Taker<State> mStateTaker;

	MeterParams params;
};


int main(int argc, char *argv[])
{
	AmbiMeterRenderer app;
	app.start();
}
