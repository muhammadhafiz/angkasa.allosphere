

#include <allocore/al_Allocore.hpp>
#include <allocore/types/al_SingleRWRingBuffer.hpp>
#include "Cuttlebone/Cuttlebone.hpp"
#include "Gamma/scl.h"

#include "soundfilebuffered.hpp"
#include "allospherespeakerlayouts.h"

#define AUDIO_BLOCK_SIZE 512

using namespace al;
using namespace std;

#define SPATIAL_SAMPLING 7

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

class AmbiMeter : public App
{
public:

	AmbiMeter() :
	    mMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
	    mNewMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
	    mTextureBuffer(((SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float) )
	{
		addSphereWithTexcoords(mMesh, 1, 128);
		addCube(mFrontalMesh, false, 0.1);

		lens().near(0.1).far(25).fovy(45);
		nav().pos(0,0,4);
		nav().quat().fromAxisAngle(0. * M_2PI, 0, 1, 0);

		initWindow(Window::Dim(0,0, 600, 400), "Ambisonics Metering", 30);

		SearchPaths sp;
		sp.addSearchPath("/Users/hafiz/Documents/thesis/media/samples/4chan/");
//		sp.addSearchPath("/home/andres/Documents/src/Allostuff/AlloSystem-github/Ambisonics/Ambisonia/");
		// sp.addSearchPath("../Ambisonics/");
		// sp.addSearchPath("../Ambisonics/NYC Ambisonics Sampler/");


		string filename = "buzzard_short.wav";
//		string filename = "NYC AMB Long Island City 01 B-Format.wav";
//		string filename = "NYC AMB Subway 02 B-Format.wav";

		string fullPath = sp.find(filename).filepath();

		if (fullPath.size() == 0) {
			cout << "ERROR: File not found!" << endl;
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

	}

	virtual void onCreate(const ViewpointWindow &win) override {

	}

	virtual void onDraw(Graphics &g) override {

//		g.polygonMode(Graphics::LINE);

		int bytesToRead = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
		while(mTextureBuffer.readSpace() >= bytesToRead) {
			mTextureBuffer.read(mMeterValues.data.ptr, bytesToRead);
			mTexture.submit(mMeterValues, true);

			mTexture.filter(Texture::NEAREST);
//			mTexture.filter(Texture::LINEAR);
		}
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

	virtual void onSound(AudioIOData &io) {

		int numFrames = io.framesPerBuffer();

		assert(AUDIO_BLOCK_SIZE == numFrames);

		//		AmbisonicsSpatializer* spatializer = ud->spatializer;
		spatializer->prepare(io);

		float * ambiChans = spatializer->ambiChans();

		int framesRead = mSoundFile->read(readBuffer, AUDIO_BLOCK_SIZE);

		if (framesRead != AUDIO_BLOCK_SIZE) {
			// TODO Handle overrun
		}

		float gain = 1.0;
		for (int chan = 0; chan < 4; chan++) {
			for (int j = 0; j < framesRead; j++) {
//				*ambiChans++ = (ud->read_buffer[j*4 + i] * ud->gain);
				ambiChans[j + chan*AUDIO_BLOCK_SIZE] = (readBuffer[j*4 + chan] * gain);
				//io.out(i, j) = (ud->read_buffer[j*4 + i] * ud->gain);
			}
		}
		spatializer->finalize(io);
	}

	static void sampleSpace(float *buf, int numChannels, int numFrames, void *userData) {
		AmbiMeter *obj = static_cast<AmbiMeter *>(userData);

		float sqrt2 = 1.0/sqrt(2.0);
		float *w = buf;
		float *x = buf + 1;
		float *y = buf + 2;
		float *z = buf + 3;
		float *d = (float *) obj->mNewMeterValues.data.ptr;

		for (int frame = 0; frame < numFrames; frame++) {

			obj->mRmsCounter++;
			for(int elevIndex = 0; elevIndex < SPATIAL_SAMPLING; elevIndex++) {
				float elev = M_PI/2.0 - (M_PI * (elevIndex + 0.5)/(float) SPATIAL_SAMPLING);
				float cosElev = COS(elev);
				for(int azimuthIndex = 0; azimuthIndex < SPATIAL_SAMPLING; azimuthIndex++) {
					float azimuth = M_2PI * (azimuthIndex + 0.5)/(float) SPATIAL_SAMPLING;
					float sampledSpace = (*w * sqrt2) + *x*COS(azimuth) * cosElev
					        + *y * SIN(azimuth) * cosElev
					        + *z * cosElev;

					obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] += sampledSpace * sampledSpace;
					if (obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] < fabs(sampledSpace)) {
						obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = fabs(sampledSpace);
					}
					if (obj->mRmsCounter == obj->mRmsSamples) { // Time for a new texture
						obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = sqrt(obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] / obj->mRmsSamples);
						float intensity;
						if (obj->params.mDecibel) {
							if (obj->params.mPeak) {
								intensity = (20 * log10(obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mPeakGain) + obj->params.mDbRange)/ obj->params.mDbRange;
							} else {
								intensity = (20 * log10(obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mRmsGain) + obj->params.mDbRange)/ obj->params.mDbRange;
							}
						} else {
							if (obj->params.mPeak) {
								intensity = obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mPeakGain;
							} else {
								intensity = obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] * obj->params.mRmsGain;
							}
						}

						if (obj->params.mGrayscale) {
							*d++ = intensity;
							*d++ = intensity;
							*d++ = intensity;
						} else {
							float r,g,b;
							if (intensity <= 0.5) {
								g = 2 * intensity;
								b = 2 * (0.5 - intensity);
								r = 0.0;
							} else {
								r = 2 * (intensity - 0.5);
								g = 2 * ( 1- intensity);
								b = 0.0;
							}
							*d++ = sqrt(r);
							*d++ = sqrt(g);
							*d++ = sqrt(b);
						}
						obj->mRmsAccum[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = 0.0;
						obj->mPeakHold[azimuthIndex * SPATIAL_SAMPLING + elevIndex] = 0.0;
						int bytesToWrite = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
						if (obj->mTextureBuffer.writeSpace() >= bytesToWrite) {
							obj->mTextureBuffer.write(obj->mNewMeterValues.data.ptr, bytesToWrite);
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
				d = (float *) obj->mNewMeterValues.data.ptr;
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
	Array mNewMeterValues;
	SingleRWRingBuffer mTextureBuffer;


	AmbisonicsSpatializer *spatializer;

	int mRmsCounter;
	int mRmsSamples;
	vector<float> mRmsAccum;
	vector<float> mPeakHold;
	SoundFileBuffered *mSoundFile;
	float readBuffer[AUDIO_BLOCK_SIZE * 4];

	MeterParams params;
};


int main(int argc, char *argv[])
{
	AmbiMeter app;
	app.start();
}
