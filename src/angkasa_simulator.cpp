#include "allocore/al_Allocore.hpp"
#include <allocore/types/al_SingleRWRingBuffer.hpp>
#include "Cuttlebone/Cuttlebone.hpp"
#include "Gamma/scl.h"

#include "soundfilebuffered.hpp"
#include "allosphere/allospherespeakerlayouts.h"
#include "state.hpp"

#include "allocore/graphics/al_Shader.hpp"
#include "alloutil/al_ShaderManager.hpp"

#define AUDIO_BLOCK_SIZE 64 // default: 512

using namespace al;
using namespace std;

#define SIN sin
#define COS cos

//#define SIN gam::scl::sinT9
//#define COS gam::scl::cosT8

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

class MyApp : public App
{
public:
  ShaderManager sm;
	Mesh sphere;

	MyApp():
    mMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
    mNewRmsMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
    mNewPeakMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
    mTextureBuffer(((SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float) ),
    mStateMaker("127.0.0.1"),
    mSpriteTex(16,16, Graphics::LUMINANCE, Graphics::FLOAT)
  {
    addSphereWithTexcoords(sphere, 1, SPATIAL_SAMPLING);
    sphere.primitive(Graphics::POINTS);

		// Create a Gaussian "bump" function to use for the sprite
		int Nx = mSpriteTex.width();
		int Ny = mSpriteTex.height();
		float * pixels = mSpriteTex.data<float>();

		for(int j=0; j<Ny; ++j){ float y = float(j)/(Ny-1)*2-1;
		for(int i=0; i<Nx; ++i){ float x = float(i)/(Nx-1)*2-1;
			float m = exp(-3*(x*x + y*y));
			pixels[j*Nx + i] = m;
		}}

		nav().pos(0,0,5);
    initWindow(Window::Dim(0,0, 600, 400), "Spatiotemporal Granulation", 30);

    // Find Ambisonics files
    SearchPaths sp;
		// sp.addSearchPath("/Users/hafiz/code/allo_SG/src/data/samples/");
		sp.addSearchPath(".", true);
		// sp.addSearchPath("/Users/hafiz/Documents/thesis/media/samples/4chan/");
    // string filename = "Fireworks_4chan_stretch.wav";
    // string filename = "Fireworks_4chan_48.wav";
		string filename = "Fireworks_0.1_constantStretch.wav";
		// string filename = "Fireworks_SGStretch_10secs.wav";
		// string filename = "crossSynth_organFireworks.wav";
		// string filename = "crossSynth_chinookFireworks.wav";
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
    //
    // Create spatializer
    SpeakerLayout speakerLayout = HeadsetSpeakerLayout();
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

  void loadShaders(){
    sm.addShaderFile("point", "point.vert", "point.frag");
  }

  void onCreate(const ViewpointWindow& w){

		loadShaders();
  }

	void onAnimate(double dt){
    if(sm.poll()){
      loadShaders();
    }
	}

	void onDraw(Graphics& g){
		// Tell GPU to render a screen-aligned textured quad at each vertex
		g.depthTesting(false);
		g.blending(true);
    g.blendAdd();

		int bytesToRead = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
    while(mTextureBuffer.readSpace() >= bytesToRead) {
      mTextureBuffer.read(mMeterValues.data.ptr, bytesToRead);
      mTexture.submit(mMeterValues, true);
      mTexture.filter(Texture::LINEAR);
    }

    ShaderProgram* s = sm.get("point");
    s->begin();
    s->uniform("texture0", 0);
		s->uniform("texture1", 1);
    // std::cout << buf[0] << '\n';
		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glEnable(GL_POINT_SMOOTH);
    mTexture.bind(0); //binds to textureID
		glEnable(GL_POINT_SPRITE);
		glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
		mSpriteTex.bind(1);
		// g.pointSize(10);
    g.draw(sphere);
		// g.rotate(90, 0, 0, 1);  // rotate sphere to show texture as we'd expect
		mSpriteTex.unbind(1);
		glDisable(GL_POINT_SPRITE);
    mTexture.unbind(0);
		glDisable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glDisable(GL_POINT_SMOOTH);
    // g.pointSize(6.0);
    s->end();
    mTexture.quad(g, 1, 1, 3, 0);
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
		}

		float gain = 1.0;
		for (int chan = 0; chan < 4; chan++) { //for every ambi channel
			for (int j = 0; j < framesRead; j++) { //for every audio block
				ambiChans[j + chan*AUDIO_BLOCK_SIZE] = (readBuffer[j*4 + chan] * gain); //interleaved
			}
		}
		spatializer->finalize(io);
	}

  	static void sampleSpace(float *buf, int numChannels, int numFrames, void *userData) { ///
  		MyApp *obj = static_cast<MyApp *>(userData);

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
  						}
              else {
  							float r,g,b;
  							if (rmsIntensity <= 0.5) { ///
									r = 0.0;
  								g = 2 * rmsIntensity;
  								b = 2 * (0.5 - rmsIntensity);
  							} else {
  								r = 2 * (rmsIntensity - 0.5);
  								g = 2 * ( 1- rmsIntensity);
  								b = 0.0;
  							}
  							*rmsBuffer++ = sqrt(r); ///
  							*rmsBuffer++ = sqrt(g);
  							*rmsBuffer++ = sqrt(b);
  						}


							// //map color to hue, so a continuous spectrum based on amplitude of grain
							//
  						// if (obj->params.mGrayscale) {
  						// 	*rmsBuffer++ = rmsIntensity; ///
  						// 	*rmsBuffer++ = rmsIntensity;
  						// 	*rmsBuffer++ = rmsIntensity;
  						// }
              // else {
  						// 	float r,g,b;
  						// 	if (rmsIntensity <= 0.5) { ///
							// 		r = 0.0;
  						// 		g = 2 * rmsIntensity;
  						// 		b = 2 * (0.5 - rmsIntensity);
  						// 	} else {
  						// 		r = 2 * (rmsIntensity - 0.5);
  						// 		g = 2 * ( 1- rmsIntensity);
  						// 		b = 0.0;
  						// 	}
  						// 	*rmsBuffer++ = sqrt(r); ///
  						// 	*rmsBuffer++ = sqrt(g);
  						// 	*rmsBuffer++ = sqrt(b);
  						// }
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



int main(){
	MyApp().start();
}


//TODO:
//- extract 1 grain, and encode into higher order. play only that grain
//-- extract another grain from a different position of same recording, and stretch them (or phase shift) at a different rate
//-- real time TPV
//-- extract 1 grain from another recording, and mix together
//- extract one section of the space, and only synthesize the grains in that sampleSpace
//-- slowly reduce the section of the space, until only 1 grain is heard (torch light)
//- RMS cross synthesis
//- Spectral cross synthesis
//-- port from ofx. trigger multiple grains in one direction
//comp: start with 1 grain from one position, triggered over and over again
// spatially morph the sound to include other grains from other positions
// trigger a granulation engine, with corresponding azimuth elevation values
// start increasing grain width, and overlapping
//constant stretch, from 1x, to reverse slow 0.01x.
//split the space up into 2 halves, and twist around axis based on envelope of signal
//each grain is then phase shifted in time, to create a "spatial delay"/ video screen tearinge-esque(?)
//red/blue separation is based on amplitude of signal, or spatial phase shift of grains
//spatialization: mix ambi and vbap? so that only a minimum number of speakers are used to play sounds from one direction?
