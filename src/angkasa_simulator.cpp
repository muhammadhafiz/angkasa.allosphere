// Angkasa: Spatiotemporal Granulation [Allosphere]
// Author: Muhammad Hafiz Wan Rosli
// January 1st 2017

#include "allocore/al_Allocore.hpp"
#include <allocore/types/al_SingleRWRingBuffer.hpp>
#include "Cuttlebone/Cuttlebone.hpp"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"

#include "alloaudio/al_SoundfileBuffered.hpp"
#include "allosphere/allospherespeakerlayouts.h"
#include "allosphere/allospherespeakerlayouts.cpp"
#include "state.hpp"

#include "allocore/graphics/al_Shader.hpp"
#include "allocore/ui/al_ParameterMIDI.hpp"
#include "alloutil/al_ShaderManager.hpp"

using namespace al;
using namespace std;

#define AUDIO_BLOCK_SIZE 512 // default: 512
#define NUM_OF_FILES 8 // num of loaded files
#define SIN sin
#define COS cos
// #define SIN gam::scl::sinT9
// #define COS gam::scl::cosT8

ParameterMIDI parameterMIDI; // KORG nanoKONTROL2
Parameter Slider0("RMS Gain", "", 0.05, "", 0.0, 1.0);
Parameter Slider1("RMS Threshold", "", 0.0, "", 0.0, 0.2);
Parameter Slider2("Master Mix", "", 0.0, "", 0.0, 1.0);

Parameter Slider6("Grain Size", "", AUDIO_BLOCK_SIZE, "", 128, AUDIO_BLOCK_SIZE);
Parameter Slider7("Master Gain", "", 0.0, "", 0.0, 1.0);

Parameter recButton("PrintOut Vals", "", 0.0, "", 0.0, 1.0);
Parameter track_backward("Change Carrier (backward)", "", 0, "", 0, 1);
Parameter track_forward("Change Carrier (forward)", "", 0, "", 0, 1);
Parameter marker_backward("Change Modulator (backward)", "", 0, "", 0, 1);
Parameter marker_forward("Change Modulator (forward)", "", 0, "", 0, 1);

class MeterParams {
public:
	MeterParams() :
	mSphere(true),
	mGrayscale(false),
	mMasterGain(0),
	mRmsGain(1.0),
	mCrossFade(0.0),
	mCarrierFile(0),
	mModulatorFile(1)
	{}
		bool mSphere;
		bool mGrayscale;

		float mMasterGain;
		float mRmsGain;
		float mCrossFade;

		int mCarrierFile;
		int mModulatorFile;

	};

	class Angkasa : public App
	{
	public:
		Angkasa():
		mMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
		mNewRmsMeterValues(3, AlloFloat32Ty, SPATIAL_SAMPLING, SPATIAL_SAMPLING),
		mTextureBuffer(((SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float)),
		mSpriteTex(16,16, Graphics::LUMINANCE, Graphics::FLOAT),
		mStateMaker(defaultBroadcastIP())
		{
			addSphereWithTexcoords(mMesh, 1, SPATIAL_SAMPLING);
			mMesh.primitive(Graphics::POINTS);

			addSphereWithTexcoords(mWireframeMesh, 1, SPATIAL_SAMPLING);
			mWireframeMesh.primitive(Graphics::LINES);

			gaussianSprite(mSpriteTex);

			nav().pos(0,0,5);
			initWindow(Window::Dim(0,0, 600, 400), "Angkasa", visFps);

			// Find Ambisonics files
			SearchPaths sp;
			sp.addSearchPath(".", true);
			// 0 = carrier, i.e. source, 1 = modulator
			filename[0] = "norm_organ.wav";
			filename[1] = "norm_fireworks.wav";
			filename[2] = "norm_insects.wav";
			filename[3] = "norm_gamelan_bapangSelisir.wav";
			filename[4] = "norm_gulls.wav";
			filename[5] = "norm_steamtrain.wav";
			filename[6] = "rain_barrett.wav";
			filename[7] = "fireworks_0.1_constantStretch.wav";
			for (int i = 0; i < NUM_OF_FILES; i++){
				fullPath[i] = sp.find(filename[i]).filepath();

				if (fullPath[i].size() == 0 || fullPath[i].size() == 0) {
					cout << "ERROR: File: "<< i << " not found!" << endl;
					return;
				} else {
					mSoundFile[i] = new SoundFileBuffered(fullPath[i],true, AUDIO_BLOCK_SIZE * 4);
					cout << "[Soundfile "<< i << "]"
							 << " Framerate: "<< mSoundFile[i]->frameRate()
							 << ", Frames: " << mSoundFile[i]->frames()
							 << ", Channels: " << mSoundFile[i]->channels()
							 << ", Samples: " << mSoundFile[i]->samples()
							 << endl;
					if (!mSoundFile[i]->opened() || mSoundFile[i]->channels() != 4) {
						cout << "ERROR: Soundfile " << i << " is invalid." << endl;
					}
				}
			} cout << endl;

			// Choice of Allosphere Speaker Layouts
			SpeakerLayout speakerLayout = HeadsetSpeakerLayout();
			if(sim()) speakerLayout = AllosphereSpeakerLayouts::threeRings54();

			// Create spatializer
			spatializer = new AmbisonicsSpatializer(speakerLayout, 3, 1);
			spatializer->numSpeakers(speakerLayout.numSpeakers());
			spatializer->numFrames(AUDIO_BLOCK_SIZE);

			mRmsSamples = mSoundFile[0]->frameRate() / visFps;
			mRmsCounter = 0;
			mRmsAccum.resize(SPATIAL_SAMPLING * SPATIAL_SAMPLING);
			for (int i = 0; i < SPATIAL_SAMPLING * SPATIAL_SAMPLING; i++) {
				mRmsAccum[i] = 0.0;
			}

			if(sim()){
				initAudio(mSoundFile[0]->frameRate(), AUDIO_BLOCK_SIZE, 60, 60);
				AudioDevice indev("ECHO X5", AudioDevice::INPUT);
				AudioDevice outdev("ECHO X5", AudioDevice::OUTPUT);
				indev.print();
				outdev.print();
				audioIO().deviceIn(indev);
				audioIO().deviceOut(outdev);
			} else initAudio(mSoundFile[0]->frameRate(), AUDIO_BLOCK_SIZE, 2, 0);
			mStateMaker.start();
		}

		// Check if using Allosphere
		bool sim(){
			std::string hostname = Socket::hostName();
			return (hostname == "gr01" || hostname == "audio.10g");
		}

		// Broadcast to Allosphere or local IP
		const char* defaultBroadcastIP(){
			if(sim()) return "192.168.10.255";
			else return "127.0.0.1";
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
			params.mMasterGain = Slider7.get();

			if (track_backward.get() || track_forward.get()){
				params.mCarrierFile-= track_backward.get();
				params.mCarrierFile+= track_forward.get();
				if (params.mCarrierFile > NUM_OF_FILES-1){
					params.mCarrierFile = 0;
				}else if (params.mCarrierFile < 0){
					params.mCarrierFile = NUM_OF_FILES-1;
				} cout << "Carrier: " << params.mCarrierFile << endl;
			}

			if (marker_backward.get() || marker_forward.get()){
				params.mModulatorFile-= marker_backward.get();
				params.mModulatorFile+= marker_forward.get();
				if (params.mModulatorFile > NUM_OF_FILES-1){
					params.mModulatorFile = 0;
				}else if (params.mModulatorFile < 0){
					params.mModulatorFile = NUM_OF_FILES-1;
				} cout << "Modulator: " << params.mModulatorFile << endl;
			}

			// Print out values for saving parameters
			if(recButton.get() == true){
				cout << "Sample 0: " << filename[params.mCarrierFile] << endl;
				cout << "Sample 1: " << filename[params.mModulatorFile] << endl;
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
			g.color(1,1,1,0.1);
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
			mTexture.bind(0); // Binds to textureID
			mSpriteTex.bind(1);
			g.draw(mMesh);
			mSpriteTex.unbind(1);
			glDisable(GL_POINT_SPRITE);
			mTexture.unbind(0);
			s->end();

			// Temporarily draw quad for texture
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

			int numFrames = AUDIO_BLOCK_SIZE;
			// int numFrames = Slider6.get();
			// assert(AUDIO_BLOCK_SIZE == numFrames);

			int framesRead[2];
			float sqrt2 = 1.0/ sqrt(2.0);
			float windowFunction[numFrames];
			gam::tbl::hann(windowFunction, numFrames);

			// int grainSize = AUDIO_BLOCK_SIZE/2;
			// framesRead[0] = mSoundFile0->read(readBuffer[0], grainSize); io.buffersize());
			framesRead[0] = mSoundFile[params.mCarrierFile]->read(readBuffer[0], numFrames);
			framesRead[1] = mSoundFile[params.mModulatorFile]->read(readBuffer[1], numFrames);
			if (framesRead[0] != numFrames || framesRead[1] != numFrames) {
				cout << "buffer overrun! framesRead[0]: " << framesRead[0] << " frames" << endl;
				cout << "buffer overrun! framesRead[1]: " << framesRead[1] << " frames" << endl;
			}

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

			// Pointer for output RMS
			float *rmsBuffer = (float *) mNewRmsMeterValues.data.ptr;

			for (int frame = 0; frame < numFrames; frame++) { // For each frame in the buffer
				mRmsCounter++;
				for(int elevIndex = 0; elevIndex < SPATIAL_SAMPLING; elevIndex++) { // For sampled elevation
					float elev = M_PI/2.0 - (M_PI * (elevIndex + 0.5)/(float) SPATIAL_SAMPLING); // -1.51844 to 1.51844
					float cosElev = COS(elev);
					for(int azimuthIndex = 0; azimuthIndex < SPATIAL_SAMPLING; azimuthIndex++) { // For each sampled azimuth
						float azimuth = M_2PI * (azimuthIndex + 0.5)/(float) SPATIAL_SAMPLING; // 0.10472 to 6.17847
						int STgrain_index = azimuthIndex * SPATIAL_SAMPLING + elevIndex;
						float STgrain_mix[AUDIO_BLOCK_SIZE];
						float rmsScaled;
						float r,g,b;
						bool STgrain_ON;

						// Calculate the STgrain for sample_0
						float STgrain_0 = (*w_0 * sqrt2)
													  + (*x_0 * COS(azimuth) * cosElev)
													  + (*y_0 * SIN(azimuth) * cosElev)
													  + (*z_0 * SIN(elev));

						// Calculate the STgrain for sample_1
						float STgrain_1 = (*w_1 * sqrt2)
													  + (*x_1 * COS(azimuth) * cosElev)
													  + (*y_1 * SIN(azimuth) * cosElev)
													  + (*z_1 * SIN(elev));

						// Accumulate STgrain^2 of each STgrain for sample_1 (for RMS)
						crossSynth_RMSAccum[STgrain_index] += STgrain_1 * STgrain_1;

						// Calculate the interpolated RMS of each STgrain of sample_1
						interp_crossSynth_RMS = ipl::linear((float) (frame)/ (numFrames-1), prev_crossSynth_RMS[STgrain_index], crossSynth_RMS[STgrain_index]);

						if(frame == numFrames-1){ // At the end of audio block
							prev_crossSynth_RMS[STgrain_index] = crossSynth_RMS[STgrain_index];
							crossSynth_RMS[STgrain_index] = sqrt(crossSynth_RMSAccum[STgrain_index]/ (float)numFrames);// * 1.2;
							crossSynth_RMSAccum[STgrain_index]= 0.0;
						}
						// STgrain_mix[frame] = STgrain_0;
						STgrain_mix[frame] = (STgrain_0 * signal_mix[0])
															 + (STgrain_1 * signal_mix[1])
															 + (STgrain_0 * (interp_crossSynth_RMS* 4.0) * signal_mix[2]);


						// Accumulate STgrain_mix^2
						mRmsAccum[STgrain_index] += STgrain_mix[frame] * STgrain_mix[frame];
						rmsScaled = mRms[STgrain_index] * params.mRmsGain;

						rmsThreshold = Slider1.get();
						if (rmsScaled <= rmsThreshold) {
							STgrain_ON = false;
							r = 0.0;
							g = 2 * rmsScaled;
							b = 2 * (rmsThreshold - rmsScaled);
						} else {
							STgrain_ON = true;
							active_STgrains++;
							r = 2 * (rmsScaled - rmsThreshold);
							g = 2 * ( 1- rmsScaled);
							b = 0.0;
						}

						// Trigger STgrains based on conditionals, in this case, RMS threshold
						if(STgrain_ON == true){
							STslice[frame][STgrain_index] = STgrain_mix[frame] * windowFunction[frame];
						}else{
							STslice[frame][STgrain_index] = 0.0;
						}

						if(frame == numFrames-1){ // At the end of audio block (255/ 511 maybe?)
							mRms[STgrain_index] = sqrt(mRmsAccum[STgrain_index] / (float) numFrames);
							mRmsAccum[STgrain_index] = 0.0;

							*rmsBuffer++ = sqrt(r); // rmsScaled- rmsThreshold;
							*rmsBuffer++ = sqrt(g); // rmsScaled- rmsThreshold;
							*rmsBuffer++ = sqrt(b); // rmsScaled- rmsThreshold;

							int bytesToWrite = SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float);
							if (mTextureBuffer.writeSpace() >= bytesToWrite) {
								mTextureBuffer.write(mNewRmsMeterValues.data.ptr, bytesToWrite);
							} else {
								// cout << "Texture Buffer overrun!" << endl;
							}
						}

						// Project each grain on to corresponding direction
						float W_reconstruct = STslice[frame][STgrain_index] * sqrt2;
						float X_reconstruct = STslice[frame][STgrain_index] * COS(azimuth) * cosElev;
						float Y_reconstruct = STslice[frame][STgrain_index] * SIN(azimuth) * cosElev;
						float Z_reconstruct = STslice[frame][STgrain_index] * SIN(elev);

						// Add each ST_grain
						*w_r += W_reconstruct;
						*x_r += X_reconstruct;
						*y_r += Y_reconstruct;
						*z_r += Z_reconstruct;

					} // End of Azi
				} // End of Elev

				// Normalize (based on acive grains) for each reconstructed ambi channel
				*w_r /= (float)(active_STgrains+1)/ 1.0;
				*x_r /= (float)(active_STgrains+1)/ 2.0;
				*y_r /= (float)(active_STgrains+1)/ 2.0;
				*z_r /= (float)(active_STgrains+1)/ 1.0;

				// Temporary: Draw the waveforms
				originalWonly[frame] = *w_0; //*w_0
				reconstructWonly[frame] = *w_r;//reconstructAmbiBuffer[frame*4] * interp_crossSynth_RMS[frame]* 2.0; //*r_z/ (reconstructAmbiNorm/2); //2,4,4,2

				// OUTPUT
				for (int chan = 0; chan < 4; chan++) { // For every Ambi Channel
					ambiChans[frame+ chan*AUDIO_BLOCK_SIZE] = reconstructAmbiBuffer[frame*4 + chan]* params.mMasterGain;
					// if (ambiChans[frame+ chan*AUDIO_BLOCK_SIZE] > maxVal){
					// 	maxVal = ambiChans[frame+ chan*AUDIO_BLOCK_SIZE];
					// 	cout << maxVal << endl;
					// }
				}

				*w_r = 0.0; *x_r = 0.0; *y_r = 0.0; *z_r = 0.0;

				// Increment pointer by 4 (4 channel interleave Ambi files)
				w_0+=4; x_0+=4; y_0+=4; z_0+=4;
				w_1+=4; x_1+=4; y_1+=4; z_1+=4;
				w_r+=4; x_r+=4; y_r+=4; z_r+=4;

				if (mRmsCounter == mRmsSamples) {
					mRmsCounter = 0;
					rmsBuffer = (float *) mNewRmsMeterValues.data.ptr;
					memcpy(mState.rmsTexture, rmsBuffer, SPATIAL_SAMPLING * SPATIAL_SAMPLING * 3 * sizeof(float));
					mStateMaker.set(mState);
				}
				active_STgrains = 0;
			} // End of audio frames
			spatializer->finalize(io);
		} // End of onSound

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

		vector<float> mRmsAccum;
		int mRmsCounter;
		int mRmsSamples;
		float mRms[SPATIAL_SAMPLING*SPATIAL_SAMPLING];

		SoundFileBuffered *mSoundFile[8];

		float readBuffer[2][AUDIO_BLOCK_SIZE * 4];
		float reconstructAmbiBuffer[AUDIO_BLOCK_SIZE * 4];

		float reconstructWonly[AUDIO_BLOCK_SIZE];
		float originalWonly[AUDIO_BLOCK_SIZE];

		float reconstructAmbiNormFactor = SPATIAL_SAMPLING * SPATIAL_SAMPLING;
		float STslice[AUDIO_BLOCK_SIZE][SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		// float STslice[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float crossSynth_RMSAccum[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float crossSynth_RMS[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float prev_crossSynth_RMS[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		float interp_crossSynth_RMS;

		float maxVal = 0.0;
		float masterMix;
		float signal_mix[3];
		float rmsThreshold;
		int active_STgrains;
		int visFps = 30;

		string filename[NUM_OF_FILES];
		string fullPath[NUM_OF_FILES];

		MeterParams params;

		State mState;
		cuttlebone::Maker<State> mStateMaker;
	};


	int main(int argc, char *argv[])
	{

		Angkasa app;

		int midiDevice = 0;
		if(app.sim()) midiDevice = 3;
		parameterMIDI.init(midiDevice, false);
		parameterMIDI.connectControl(Slider0, 0, 1);
		parameterMIDI.connectControl(Slider1, 1, 1);
		parameterMIDI.connectControl(Slider2, 2, 1);

		parameterMIDI.connectControl(Slider6, 6, 1);
		parameterMIDI.connectControl(Slider7, 7, 1);

		parameterMIDI.connectControl(recButton, 45, 1);
		parameterMIDI.connectControl(track_backward, 58, 1);
		parameterMIDI.connectControl(track_forward, 59, 1);
		parameterMIDI.connectControl(marker_backward, 61, 1);
		parameterMIDI.connectControl(marker_forward, 62, 1);

		app.start();
	}
