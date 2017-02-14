// Angkasa: Spatiotemporal Granulation [Allosphere]
// Author: Muhammad Hafiz Wan Rosli
// January 1st 2017

#include "allocore/al_Allocore.hpp"
#include <allocore/types/al_SingleRWRingBuffer.hpp>
#include "Cuttlebone/Cuttlebone.hpp"
#include "Gamma/scl.h"

#include "soundfilebuffered.hpp"
#include "allosphere/allospherespeakerlayouts.h"
#include "state.hpp"

#include "allocore/graphics/al_Shader.hpp"
#include "allocore/ui/al_ParameterMIDI.hpp"
#include "alloutil/al_ShaderManager.hpp"

using namespace al;
using namespace std;

#define AUDIO_BLOCK_SIZE 1024 // default: 512
#define SIN sin
#define COS cos
// #define SIN gam::scl::sinT9
// #define COS gam::scl::cosT8

ParameterMIDI parameterMIDI; // KORG nanoKONTROL2
Parameter Slider0("RMS Gain", "", 0.85, "", 0.0, 1.0);
Parameter Slider1("RMS Threshold", "", 0.0, "", 0.0, 0.5);
Parameter Slider2("Master Mix", "", 0.0, "", 0.0, 1.0);
Parameter Slider7("Master Gain", "", 0.0, "", 0.0, 1.0);
Parameter recButton("PrintOut Vals", "", 0.0, "", 0.0, 1.0);

class MeterParams {
public:
	MeterParams() :
	mSphere(true),
	mDecibel(false),
	mGrayscale(false),
	mRmsGain(1.0), //20.0
	mCrossFade(0.0),
	mDbRange(80.0)
	{}

		bool mSphere;
		bool mDecibel;
		bool mGrayscale;

		float mRmsGain;
		float mCrossFade;
		double mDbRange;
	};

	class Angkasa : public App
	{
	public:
		Angkasa():
		mMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
		mNewRmsMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
		mTextureBuffer(((SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float)),
		mSpriteTex(16,16, Graphics::LUMINANCE, Graphics::FLOAT),
		mStateMaker("127.0.0.1")
		// mStateMaker("192.168.10.255")
		{
			addSphereWithTexcoords(mMesh, 1, SPATIAL_SAMPLING);
			mMesh.primitive(Graphics::POINTS);

			addSphereWithTexcoords(mWireframeMesh, 1, SPATIAL_SAMPLING);
			mWireframeMesh.primitive(Graphics::LINES);

			gaussianSprite(mSpriteTex);

			nav().pos(0,0,5);
			initWindow(Window::Dim(0,0, 600, 400), "Angkasa" /*| Spatiotemporal Granulation*/, visFps);

			// Find Ambisonics files
			SearchPaths sp;
			sp.addSearchPath(".", true);
			// 0 = carrier, i.e. source, 1 = modulator, i.e. temporal modifier
			// filename[0] = "norm_gulls.wav";
			// filename[0] = "norm_gamelan_bapangSelisir.wav";
			// filename[0] = "norm_insects.wav";
			// filename[0] = "norm_steamtrain.wav";
			filename[0] = "norm_organ.wav";
			filename[1] = "norm_fireworks.wav";
			for (int i = 0; i < 2; i++){
				fullPath[i] = sp.find(filename[i]).filepath();
			}
			if (fullPath[0].size() == 0 || fullPath[1].size() == 0) {
				cout << "ERROR: File not found!" << endl;
				return;
			} else {
				mSoundFile0 = new SoundFileBuffered(fullPath[0],true, AUDIO_BLOCK_SIZE * 4);
				mSoundFile1 = new SoundFileBuffered(fullPath[1],true, AUDIO_BLOCK_SIZE * 4);

				if (!mSoundFile0->opened() || mSoundFile0->channels() != 4) {
					cout << "ERROR: Soundfile invalid." << endl;
				}
			}

			// Choice of Allosphere Speker Layouts
			SpeakerLayout speakerLayout = HeadsetSpeakerLayout();

			// Create spatializer
			spatializer = new AmbisonicsSpatializer(speakerLayout, 3, 1);
			spatializer->numSpeakers(speakerLayout.numSpeakers());
			spatializer->numFrames(AUDIO_BLOCK_SIZE);

			mRmsSamples = mSoundFile0->frameRate() / visFps;
			mRmsCounter = 0;
			mRmsAccum.resize(SPATIAL_SAMPLING * SPATIAL_SAMPLING);
			for (int i = 0; i < SPATIAL_SAMPLING * SPATIAL_SAMPLING; i++) {
				mRmsAccum[i] = 0.0;
			}

			initAudio(mSoundFile0->frameRate(), AUDIO_BLOCK_SIZE, 2, 0);
			mStateMaker.start();
		}

		// Draw waveform
		static void addWaveformDisplay(Mesh& m, float *buffer, float R, float G, float B) {
			float bufferSize = AUDIO_BLOCK_SIZE;

			int zoomOut = 2;  // # audio samples per OpenGL Mesh vertex
			m.primitive(Graphics::LINE_STRIP);
			int n = 0;  // current sample #
			for (int i = 0; i < bufferSize / zoomOut; ++i) {
				float max = -1, min = 1;
				for (int j = 0; j < zoomOut; ++j) {
					if (buffer[n] > max) max = buffer[n];
					if (buffer[n] < min) min = buffer[n];
					n++;
				}
				float verticalScale = 0.5;
				m.color(RGB(R,G,B));
				m.vertex((float)i / (bufferSize / zoomOut), max * verticalScale, 0);
				m.vertex((float)i / (bufferSize / zoomOut), min * verticalScale, 0);
			}
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
			shaderManager.addShaderFile("point", "point.vert", "point.frag");
		}

		virtual void onCreate(const ViewpointWindow &win) override {
			loadShaders();
		}

		virtual void onAnimate(double dt) override{
			if(shaderManager.poll()){
				loadShaders();
			}
			params.mRmsGain = Slider0.get();
			params.mCrossFade = Slider2.get();
			// Print out values for saving parameters
			if(recButton.get() == true){
				cout << "Sample 0: " << filename[0] << endl;
				cout << "Sample 1: " << filename[1] << endl;
				cout << "Slider0 val: " << Slider0.get() << endl;
				cout << "Slider1 val: " << Slider1.get() << endl;
				cout << "Slider2 val: " << Slider2.get() << endl;
				cout << "Slider7 val: " << Slider7.get() << endl;
			}
			// Master mix for playing different Ambi files
			masterMix = Slider2.get();
			if(masterMix < 0.5){
				signal_mix[0] =  1 - (masterMix*2);
				signal_mix[1] = masterMix*2;
				signal_mix[2] = 0;
			} else if(masterMix > 0.5){
				signal_mix[0] = 0;
				signal_mix[1] = (1- masterMix)*2;
				signal_mix[2] = (masterMix*2) - 1;
			}
		}

		virtual void onDraw(Graphics &g) override {
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SMOOTH);
			glEnable(GL_POINT_SPRITE);
			glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
			g.blendAdd();
			// Draw wireframe mesh
			g.lineWidth(0.1);
			g.color(1,1,1,0.25);
			g.draw(mWireframeMesh);

			// Read the texture buffer
			int bytesToRead = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
			while(mTextureBuffer.readSpace() >= bytesToRead) {
				mTextureBuffer.read(mMeterValues.data.ptr, bytesToRead);
				mTexture.submit(mMeterValues, true);
				mTexture.filter(Texture::LINEAR);
			}
			ShaderProgram* s = shaderManager.get("point");
			s->begin();
			s->uniform("texture0", 0);
			s->uniform("texture1", 1);
			mTexture.bind(0); //binds to textureID
			mSpriteTex.bind(1);
			g.draw(mMesh);
			mSpriteTex.unbind(1);
			glDisable(GL_POINT_SPRITE);
			mTexture.unbind(0);
			s->end();

			mTexture.quad(g, 1, 1, 3, 0);

			// Draw waveform(s)
			Mesh waveform0;
			Mesh waveform1;
			addWaveformDisplay(waveform0, originalWonly, 1, 0, 0);
			addWaveformDisplay(waveform1, reconstructWonly, 0, 1, 0);
			g.translate(1.0, -1.0, 0);
			g.draw(waveform0);
			g.draw(waveform1);
		}

		virtual void onSound(AudioIOData &io) override{
			spatializer->prepare();
			float * ambiChans = spatializer->ambiChans();

			int numFrames = io.framesPerBuffer();
			assert(AUDIO_BLOCK_SIZE == numFrames);

			int framesRead[2];
			framesRead[0] = mSoundFile0->read(readBuffer[0], AUDIO_BLOCK_SIZE);
			framesRead[1] = mSoundFile1->read(readBuffer[1], AUDIO_BLOCK_SIZE);
			if (framesRead[0] != AUDIO_BLOCK_SIZE || framesRead[1] != AUDIO_BLOCK_SIZE) {
				cout << "buffer overrun! framesRead[0]: " << framesRead[0] << " frames" << endl;
				cout << "buffer overrun! framesRead[1]: " << framesRead[1] << " frames" << endl;
			}

			float sqrt2 = 1.0/sqrt(2.0);
			float reconstructAmbiNormFactor = SPATIAL_SAMPLING * SPATIAL_SAMPLING;

			// Pointer for buffer 0
			float *w_0 = readBuffer[0];
			float *x_0 = readBuffer[0] + 1;
			float *y_0 = readBuffer[0] + 2;
			float *z_0 = readBuffer[0] + 3;

			// Pointer for buffer 1
			float *w_1 = readBuffer[1];
			float *x_1 = readBuffer[1] + 1;
			float *y_1 = readBuffer[1] + 2;
			float *z_1 = readBuffer[1] + 3;

			// Pointer for reconstructed buffer
			float *w_r = reconstructAmbiBuffer;
			float *x_r = reconstructAmbiBuffer + 1;
			float *y_r = reconstructAmbiBuffer + 2;
			float *z_r = reconstructAmbiBuffer + 3;

			float *rmsBuffer = (float *) mNewRmsMeterValues.data.ptr;

			for (int frame = 0; frame < numFrames; frame++) { // For each frame in the buffer
				// mRmsCounter++;
				float maxVal= 0;

				for(int elevIndex = 0; elevIndex < SPATIAL_SAMPLING; elevIndex++) { // For sampled elevation
					float elev = M_PI/2.0 - (M_PI * (elevIndex + 0.5)/(float) SPATIAL_SAMPLING); // -1.51844 to 1.51844
					float cosElev = COS(elev);
					for(int azimuthIndex = 0; azimuthIndex < SPATIAL_SAMPLING; azimuthIndex++) { // For each sampled azimuth
						float azimuth = M_2PI * (azimuthIndex + 0.5)/(float) SPATIAL_SAMPLING; // 0.10472 to 6.17847
						int which_STgrain = azimuthIndex * SPATIAL_SAMPLING + elevIndex; // Current grain
						float STgrain_mix;
						float rmsIntensity;
						float r,g,b;

						// Calculate the STgrain for sample_0
						float STgrain_0= (*w_0 * sqrt2)
													 + (*x_0 * COS(azimuth) * cosElev)
													 + (*y_0 * SIN(azimuth) * cosElev)
													 + (*z_0 * SIN(elev));

						// Calculate the STgrain for sample_1
						float STgrain_1= (*w_1 * sqrt2)
													 + (*x_1 * COS(azimuth) * cosElev)
													 + (*y_1 * SIN(azimuth) * cosElev)
													 + (*z_1 * SIN(elev));

						// Accumulate STgrain^2 of each STgrain for sample_1 (for RMS)
						crossSynth_RMSAccum[which_STgrain] += STgrain_1 * STgrain_1;

						// Calculate the interpolated RMS of each STgrain of sample_1
						interp_crossSynth_RMS = ipl::linear((float) (frame)/ (numFrames-1), prev_crossSynth_RMS[which_STgrain], crossSynth_RMS[which_STgrain]);

						if(frame == numFrames-1){ // At the end of audio block
							prev_crossSynth_RMS[which_STgrain] = crossSynth_RMS[which_STgrain];
							crossSynth_RMS[which_STgrain] = sqrt(crossSynth_RMSAccum[which_STgrain]/ (float)numFrames);// * 1.2;
							crossSynth_RMSAccum[which_STgrain]= 0.0;
						}
						STgrain_mix = (STgrain_0 * signal_mix[0])
												+ (STgrain_1 * signal_mix[1])
												+ (STgrain_0 * (interp_crossSynth_RMS* 4.0) * signal_mix[2]);
						// STgrain_mix = STgrain_0 * (interp_crossSynth_RMS* 4.0);
						// STgrain_mix = STgrain_0;


						// Accumulate STgrain_mix^2
						mRmsAccum[which_STgrain] += STgrain_mix * STgrain_mix;
						mInterp_RMS = ipl::linear((float) (frame)/ (numFrames-1), mPrev_RMS[which_STgrain], mCurrent_RMS[which_STgrain]);

						rmsIntensity = mInterp_RMS * params.mRmsGain;
						//control the exponent for source width control (audio & gfx)

						rmsThreshold = Slider1.get();
						if (rmsIntensity <= rmsThreshold) {
							r = 0.0;
							g = 2 * rmsIntensity;
							b = 2 * (rmsThreshold - rmsIntensity);
							rmsIntensity = rmsThreshold;
							// STslice[frame][which_STgrain] = 0.0;
						} else {
							active_STgrains++;
							r = 2 * (rmsIntensity - rmsThreshold);
							g = 2 * ( 1- rmsIntensity);
							b = 0.0;
							// STslice[frame][which_STgrain] = rmsIntensity * STgrain_mix* 2.0;
							// STslice[frame][which_STgrain] = STgrain_mix;
						}

						STslice[frame][which_STgrain] = (rmsIntensity- rmsThreshold) * STgrain_mix * 2.0;
						// STslice[frame][whichGrain] = ((rmsIntensity- rmsThreshold)* (1.0/ (1.0-rmsThreshold))) * STgrain_mix;

						if(frame == numFrames-1){ // At the end of audio block (255/ 511 maybe?)
							mPrev_RMS[which_STgrain] = mCurrent_RMS[which_STgrain];
							mCurrent_RMS[which_STgrain] = sqrt(mRmsAccum[which_STgrain] / (float) numFrames);
							mRmsAccum[which_STgrain] = 0.0;

							*rmsBuffer++ = sqrt(r); //rmsIntensity- rmsThreshold;
							*rmsBuffer++ = sqrt(g); //rmsIntensity- rmsThreshold;
							*rmsBuffer++ = sqrt(b); //rmsIntensity- rmsThreshold;

							int bytesToWrite = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
							if (mTextureBuffer.writeSpace() >= bytesToWrite) {
								mTextureBuffer.write(mNewRmsMeterValues.data.ptr, bytesToWrite);
							} else {
								// cout << "Texture Buffer overrun!" << endl;
							}
						}

						// Project each grain on to corresponding direction
						float W_reconstruct = STslice[frame][which_STgrain] * sqrt2;
						float X_reconstruct = STslice[frame][which_STgrain] * COS(azimuth) * cosElev;
						float Y_reconstruct = STslice[frame][which_STgrain] * SIN(azimuth) * cosElev;
						float Z_reconstruct = STslice[frame][which_STgrain] * SIN(elev);

						// Add each grain
						*w_r += W_reconstruct;
						*x_r += X_reconstruct;
						*y_r += Y_reconstruct;
						*z_r += Z_reconstruct;

						// if(*w_r > maxVal){
						// 	maxVal = *w_r ;
						// 	cout << maxVal << endl;
						// }
					} // End of Azi
				} // End of elev

				// Normalize for each reconstructed ambi channel
				// *w_r /= (reconstructAmbiNormFactor/2.0);
				// float amplitudeControl = 1 - (active_STgrains+1/ 101);
				// float amplitudeControl = (active_STgrains+1)/ 2.0;
				// float amplitudeControl = (reconstructAmbiNormFactor/2)/ (101-active_STgrains);
				// float amplitudeControl = (active_STgrains+1) * 2/ (active_STgrains+1);
				// float amplitudeControl = 50.0/ ((101-active_STgrains)* 2.0));
				// cout << amplitudeControl << endl;

				// 100 grains -> w_r/50 	( 100/2)
				// 100 grains -> w_r/ (reconstructAmbiNormFactor/2)/ (101-active_STgrains) = 1/ (50/ 1) = 0.002
				// 50 grain -> w_r/ (reconstructAmbiNormFactor/2)/ (101-active_STgrains) = 1/ (50/ 51) = 1.0
				// 1 grain -> w_r/ (reconstructAmbiNormFactor/2)/ (101-active_STgrains) = 1/ (50/ 100) = 2

				// 100 GRAINS -> 50  50/ (101-active_STgrains)^2
				// 1 GRAIN -> 0.05  50/ (101-active_STgrains)^2
				// ((active_STgrains+1) * 5) / (active_STgrains/10)

				// *w_r /= active_STgrains;
				// cout << active_STgrains << endl;
				// float amplitudeControl = pow((active_STgrains+1), 2) / (reconstructAmbiNormFactor);

				*w_r /= reconstructAmbiNormFactor/ 2.0;
				*x_r /= reconstructAmbiNormFactor/ 4.0;
				*y_r /= reconstructAmbiNormFactor/ 4.0;
				*z_r /= reconstructAmbiNormFactor/ 2.0;
				//
				// *w_r /= (float)(active_STgrains+1)/1.0;
				// *x_r /= (float)(active_STgrains+1)/2.0;
				// *y_r /= (float)(active_STgrains+1)/2.0;
				// *z_r /= (float)(active_STgrains+1)/1.0;
				//
				// *w_r /= (reconstructAmbiNormFactor/ amplitudeControl);
				// *x_r /= (reconstructAmbiNormFactor/ (amplitudeControl * 2));
				// *y_r /= (reconstructAmbiNormFactor/ (amplitudeControl * 2));
				// *z_r /= (reconstructAmbiNormFactor/ amplitudeControl);


				// Temporary: Draw the waveforms
				originalWonly[frame] = *w_1; //*w_0
				reconstructWonly[frame] = *w_r;//reconstructAmbiBuffer[frame*4] * interp_crossSynth_RMS[frame]* 2.0; //*r_z/ (reconstructAmbiNorm/2); //2,4,4,2

				// OUTPUT
				float masterGain = Slider7.get(); // Master Gain
				for (int chan = 0; chan < 4; chan++) { // For every Ambi Channel
					ambiChans[frame+ chan*AUDIO_BLOCK_SIZE] = reconstructAmbiBuffer[frame*4 + chan]* masterGain;
				}
				*w_r = 0.0; *x_r = 0.0; *y_r = 0.0; *z_r = 0.0;

				w_0+=4; x_0+=4; y_0+=4; z_0+=4;
				w_1+=4; x_1+=4; y_1+=4; z_1+=4;
				w_r+=4; x_r+=4; y_r+=4; z_r+=4;

				// if (mRmsCounter == mRmsSamples) {
				// 	mRmsCounter = 0;
				// 	rmsBuffer = (float *) mNewRmsMeterValues.data.ptr;
				// 	memcpy(mState.rmsTexture, rmsBuffer, SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float));
				// 	mStateMaker.set(mState);
				// }
				active_STgrains = 0;
			} // End of audio frames
			spatializer->finalize(io);
		} // End of onSound


		virtual void onKeyDown(const Keyboard& k) override {
			if (k.key() == 'g') {
				params.mGrayscale = !params.mGrayscale;
			}
		}

	private:
		Mesh mMesh;
		Mesh mWireframeMesh;
		Texture mTexture;
		Texture mSpriteTex;
		Array mMeterValues;
		Array mNewRmsMeterValues;
		SingleRWRingBuffer mTextureBuffer;

		AmbisonicsSpatializer *spatializer;
		ShaderManager shaderManager;

		int mRmsCounter;
		int mRmsSamples;
		vector<float> mRmsAccum;

		float mInterp_RMS;
		float mPrev_RMS[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float mCurrent_RMS[SPATIAL_SAMPLING*SPATIAL_SAMPLING];

		SoundFileBuffered *mSoundFile0;
		SoundFileBuffered *mSoundFile1;

		float readBuffer[2][AUDIO_BLOCK_SIZE * 4];
		float reconstructAmbiBuffer[AUDIO_BLOCK_SIZE * 4];

		float reconstructWonly[AUDIO_BLOCK_SIZE];
		float originalWonly[AUDIO_BLOCK_SIZE];

		float STslice[AUDIO_BLOCK_SIZE][SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		// float STslice[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float crossSynth_RMSAccum[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float crossSynth_RMS[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float prev_crossSynth_RMS[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float interp_crossSynth_RMS;

		float masterMix;
		float signal_mix[3];
		float rmsThreshold;
		int active_STgrains;
		int visFps = 30;

		string filename[2];
		string fullPath[2];

		MeterParams params;

		State mState;
		cuttlebone::Maker<State> mStateMaker;
	};


	int main(int argc, char *argv[])
	{			float maxVal= 0;

		parameterMIDI.init();
		parameterMIDI.connectControl(Slider0, 0, 1);
		parameterMIDI.connectControl(Slider1, 1, 1);
		parameterMIDI.connectControl(Slider2, 2, 1);
		parameterMIDI.connectControl(Slider7, 7, 1);
		parameterMIDI.connectControl(recButton, 45, 1);

		Angkasa app;
		app.start();
	}
