#include <allocore/al_Allocore.hpp>
#include <allocore/types/al_SingleRWRingBuffer.hpp>
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

class AmbiMeterRenderer : public App
{
public:

	AmbiMeterRenderer() :
	    mMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
	    mNewRmsMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
	    mNewPeakMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
	    mTextureBuffer(((SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float) ),
	    mStateMaker("127.0.0.1"),

			mSpriteTex(16,16, Graphics::LUMINANCE, Graphics::FLOAT)
			// mStateMaker("192.168.10.255")
	{
		addSphereWithTexcoords(mMesh, 1, 64);
		// std::cout << mMesh.vertices().size() << '\n';
		// mMesh.primitive(Graphics::POINTS);
		// mMesh.vertex(0,0,0);
		// mMesh.color(0.5, 0.5, 0.5, 1);
		// addCube(mFrontalMesh, false, 0.1);
		// mMesh.color(0.5, 0., 0., 1);

		// Create a Gaussian "bump" function to use for the sprite
		int Nx = mSpriteTex.width();
		int Ny = mSpriteTex.height();
		float * pixels = mSpriteTex.data<float>();

		for(int j=0; j<Ny; ++j){ float y = float(j)/(Ny-1)*2-1;
		for(int i=0; i<Nx; ++i){ float x = float(i)/(Nx-1)*2-1;
			float m = exp(-3*(x*x + y*y));
			pixels[j*Nx + i] = m;
		}}

		lens().near(0.1).far(25).fovy(45);
		nav().pos(0,0,4);
		nav().quat().fromAxisAngle(0. * M_2PI, 0, 1, 0);

		initWindow(Window::Dim(0,0, 600, 400), "Spatiotemporal Granulation: Stretch", 30);

		SearchPaths sp;
		sp.addSearchPath("/Users/hafiz/Documents/thesis/media/samples/4chan/");
		// sp.addSearchPath("../Ambisonics/");
		// string filename = "Fireworks_4chan_48.wav";
		string filename = "Fireworks_4chan_stretch.wav";
		string fullPath = sp.find(filename).filepath();

		if (fullPath.size() == 0) {
			cout << "ERROR: File not found!" << endl;
			return;
		} else {
			mSoundFile = new SoundFileBuffered(fullPath);
			if (!mSoundFile->opened() || mSoundFile->channels() != 4) {
				cout << "ERROR: Soundfile invalid." << endl;
			}
			mSoundFile->setReadCallback(sampleSpace, this);
		}


//		SpeakerLayout speakerLayout = AllosphereSpeakerLayouts::threeRings54();
//		SpeakerLayout speakerLayout = AllosphereSpeakerLayouts::threeRings27();
//		SpeakerLayout speakerLayout = AllosphereSpeakerLayouts::centerRing();
//		SpeakerLayout speakerLayout = AllosphereSpeakerLayouts::topAndBottom();
//		SpeakerLayout speakerLayout = AllosphereSpeakerLayouts::sparse();
		SpeakerLayout speakerLayout = HeadsetSpeakerLayout();

		// Create spatializer
		spatializer = new AmbisonicsSpatializer(speakerLayout, 3, 1);
		spatializer->numSpeakers(speakerLayout.numSpeakers());
		spatializer->numFrames(AUDIO_BLOCK_SIZE);

		mRmsSamples = mSoundFile->frameRate() / 25;
		mRmsCounter = 0;
		mRmsAccum.resize(SPATIAL_SAMPLING * SPATIAL_SAMPLING);
		mPeakHold.resize(SPATIAL_SAMPLING * SPATIAL_SAMPLING);
		for (int i = 0; i < SPATIAL_SAMPLING * SPATIAL_SAMPLING; i++) {
			mRmsAccum[i] = 0.0;
			mPeakHold[i] = 0.0;
		}

		initAudio(mSoundFile->frameRate(), AUDIO_BLOCK_SIZE, 2, 0);
		mStateMaker.start();

	}

	virtual void onCreate(const ViewpointWindow &win) override {

	}

	virtual void onDraw(Graphics &g) override {
		glEnable(GL_POINT_SPRITE);
		glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);///
		g.blendAdd();

		glEnable(GL_POINT_SMOOTH);
		// glEnable( GL_BLEND );
		// glBlendFunc( Graphics::SRC_ALPHA, Graphics::DST_ALPHA);
		g.pointSize(6.0);
		g.polygonMode(Graphics::POINT);
		g.antialiasing(Graphics::NICEST);
		g.depthTesting(true);
		g.blending(true);
		// g.pointSize(3);
		// g.shadeModel(Graphics::SMOOTH);
		// g.nicest();

		// g.polygonMode(Graphics::LINE);
		// g.pointAtten(0, 0, 1);

		// // Enable blending to hide texture edges
		// g.blendAdd();

		int bytesToRead = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
		while(mTextureBuffer.readSpace() >= bytesToRead) {
			mTextureBuffer.read(mMeterValues.data.ptr, bytesToRead);
			mTexture.submit(mMeterValues, true);
			// mTexture = mSpriteTex;
			// mTexture.filter(Texture::NEAREST);
			mTexture.filter(Texture::LINEAR);
			// float *meterval = (float *) obj->mMeterValues.data.ptr;
			// meterval = (float *) obj->mMeterValues.data.ptr;

			// float meterval = mMeterValues[0];
			// mMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),

			// std::cout << (float *) mMeterValues.data.ptr << '\n';
		}
		// mMeterValues.read(mMeterValues.data.ptr, bytesToRead);
		// std::cout << mMeterValues.elem<float>(0,0,0) << '\n';


		float buf[3];
		mMeterValues.read(buf, 0, 0);
		std::cout << buf[0] << '\n';

		if (params.mSphere) {
			// mTexture.bind();
			mSpriteTex.bind();
			g.draw(mMesh);
			mSpriteTex.unbind();
			// mTexture.unbind();
		}

		// for (int i = 0; i< 64; i++) {
		// 	g.pointsize(map_size(meterval[i]))
		// 	g.color(some_function(i, meshed[i].vertices()[0].x))
		// 	g.draw(meshes[i])
		// }

		// mSpriteTex.bind()
		// for (int i = 0; i< 64; i++) {
		// 	g.pointsize(map_size(meterval[i]))
		// 	g.color(some_function(i, meshed[i].vertices()[0].x))
		// 	g.draw(meshes[i])
		// }
		// mSpriteTex.unbind()
		//
		//
		//
		// // setup
		// addSphereWithTexcoords(mMesh, 1, 64);
		// Mesh meshes[64];
		//
		// for(int i =0; i< 64; i++){
		// 	Mesh& m = meshes[i]
		// 	m.reset();
		// 	m.vertex(mMesh.vertices()[i]);
		// }
		// //
		// // for (i ~ 64) {
		// // 	Mesh& m = meshes[i]
		// // 	m.reset();
		// // 	m.vertex(mMesh.vertices()[i]);
		// // }
		//
		//
		// // draw
		// //
		// // sprite.bind()
		// // for () {
		// // 	g.pointsize(map_size(meterval[i]))
		// // 	g.color(some_function(i, meshed[i].vertices()[0].x))
		// // 	g.draw(meshes[i])
		// // }
		// // sprite.unbind()
		//
		// // glDisable(GL_POINT_SMOOTH);
		// // glDisable(GL_POINT_SPRITE);
	}

	virtual void onSound(AudioIOData &io) {

		int numFrames = io.framesPerBuffer();
		// std::cout << numFrames << '\n';
		assert(AUDIO_BLOCK_SIZE == numFrames);

		spatializer->prepare();

		float * ambiChans = spatializer->ambiChans();

		int framesRead = mSoundFile->read(readBuffer, AUDIO_BLOCK_SIZE);
		// std::cout << framesRead << '\n';
		if (framesRead != AUDIO_BLOCK_SIZE) {
			// TODO Handle overrun
			// std::cout << "framesRead != AUDIO_BLOCK_SIZE" << '\n';
		}

		float gain = 1.0;
		for (int chan = 0; chan < 4; chan++) { //for every ambi channel
			for (int j = 0; j < framesRead; j++) { //for every audio block
				ambiChans[j + chan*AUDIO_BLOCK_SIZE] = (readBuffer[j*4 + chan] * gain); //interleaved
				// std::cout << j + chan*AUDIO_BLOCK_SIZE << '\n';
			}
		}
		spatializer->finalize(io);
	}

	static void sampleSpace(float *buf, int numChannels, int numFrames, void *userData) { ///
		AmbiMeterRenderer *obj = static_cast<AmbiMeterRenderer *>(userData);

		float sqrt2 = 1.0/sqrt(2.0);
		float *w = buf;
		float *x = buf + 1;
		float *y = buf + 2;
		float *z = buf + 3;
		float *rmsBuffer = (float *) obj->mNewRmsMeterValues.data.ptr;
		float *peakBuffer = (float *) obj->mNewPeakMeterValues.data.ptr;

		for (int frame = 0; frame < numFrames; frame++) {

			obj->mRmsCounter++;
			for(int elevIndex = 0; elevIndex < SPATIAL_SAMPLING; elevIndex++) {
				float elev = M_PI/2.0 - (M_PI * (elevIndex + 0.5)/(float) SPATIAL_SAMPLING);
				// std::cout << elev << '\n';
				float cosElev = COS(elev);
				for(int azimuthIndex = 0; azimuthIndex < SPATIAL_SAMPLING; azimuthIndex++) {
					float azimuth = M_2PI * (azimuthIndex + 0.5)/(float) SPATIAL_SAMPLING;
					float sampledSpace = (*w * sqrt2) + (*x * COS(azimuth) * cosElev)
					        + (*y * SIN(azimuth) * cosElev)
					        + (*z * cosElev);

					obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] += sampledSpace * sampledSpace;

					if (obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] < fabs(sampledSpace)) {
						obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = fabs(sampledSpace);//Hold the peak
					}

					if (obj->mRmsCounter == obj->mRmsSamples) { // Time for a new texture
						obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = sqrt(obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] / obj->mRmsSamples); //RMS
						float rmsIntensity;
						float peakIntensity;
						if (obj->params.mDecibel) {
							peakIntensity = (20 * log10(obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mPeakGain) + obj->params.mDbRange)/ obj->params.mDbRange;
							rmsIntensity = (20 * log10(obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mRmsGain) + obj->params.mDbRange)/ obj->params.mDbRange;
						} else {
							peakIntensity = obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mPeakGain;
							rmsIntensity = obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mRmsGain;
					}

						if (obj->params.mGrayscale) {
							*rmsBuffer++ = rmsIntensity; ///
							*rmsBuffer++ = rmsIntensity;
							*rmsBuffer++ = rmsIntensity;
						} else {
							float r,g,b;
							if (rmsIntensity <= 0.5) { ///
								g = 2 * rmsIntensity;
								b = 2 * (0.5 - rmsIntensity);
								r = 0.0;
							} else {
								r = 2 * (rmsIntensity - 0.5);
								g = 2 * ( 1- rmsIntensity);
								b = 0.0;
							}
							*rmsBuffer++ = sqrt(r); ///
							*rmsBuffer++ = sqrt(g);
							*rmsBuffer++ = sqrt(b);
						}
//						obj->mState.rmsTexture[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex];

//						obj->mState.peakTexture[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex];

						obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = 0.0;
						obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = 0.0;
						int bytesToWrite = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
						if (obj->mTextureBuffer.writeSpace() >= bytesToWrite) {
							obj->mTextureBuffer.write(obj->mNewRmsMeterValues.data.ptr, bytesToWrite);
						} else {
							//							cout << "Texture Buffer overrun!" << endl;
						}
					}
				}
			}

			w += 4;
			x += 4;
			y += 4;
			z += 4;

			if (obj->mRmsCounter == obj->mRmsSamples) {
				obj->mRmsCounter = 0;
				rmsBuffer = (float *) obj->mNewRmsMeterValues.data.ptr;
				peakBuffer = (float *) obj->mNewPeakMeterValues.data.ptr;
				memcpy(obj->mState.peakTexture, peakBuffer, SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float));
				memcpy(obj->mState.rmsTexture, rmsBuffer, SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float));
				obj->mStateMaker.set(obj->mState);
			}
		}
	}

	virtual void onKeyDown(const Keyboard& k) override {
		if (k.key() == 'g') {
			params.mGrayscale = !params.mGrayscale;
		} else if (k.key() == 'h') {
			params.mDecibel = !params.mDecibel;
		} else if (k.key() == 'j') {
			params.mSphere = !params.mSphere;
		} else if (k.key() == 'k') {
			params.mPeak = !params.mPeak;
		}
	}

private:
	Mesh mMesh;
	Mesh mFrontalMesh;
	Texture mTexture;
	Array mMeterValues;
	Array mNewRmsMeterValues;
	Array mNewPeakMeterValues;
	SingleRWRingBuffer mTextureBuffer;
	Texture mSpriteTex;


	AmbisonicsSpatializer *spatializer;

	int mRmsCounter;
	int mRmsSamples;
	vector<float> mRmsAccum;
	vector<float> mPeakHold;
	SoundFileBuffered *mSoundFile;
	float readBuffer[AUDIO_BLOCK_SIZE * 4];

	MeterParams params;

	State mState;
	cuttlebone::Maker<State> mStateMaker;
};


int main(int argc, char *argv[])
{
	AmbiMeterRenderer app;
	app.start();
}
