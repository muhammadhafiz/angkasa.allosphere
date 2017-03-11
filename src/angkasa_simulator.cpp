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

#include "granulate.h"

using namespace al;
using namespace std;

// #define AUDIO_BLOCK_SIZE 1024 // default: 512
#define NUM_OF_FILES 3 // num of loaded files
#define TOTAL_SPATIAL_SAMPLING SPATIAL_SAMPLING * SPATIAL_SAMPLING
#define SIN sin
#define COS cos
// #define SIN gam::scl::sinT9
// #define COS gam::scl::cosT8

ParameterMIDI parameterMIDI; // KORG nanoKONTROL2
Parameter Slider0("RMS Gain", "", 1.0, "", 0.0, 1.0);
Parameter Slider1("RMS Threshold", "", 0.0, "", 0.0, 0.3);
Parameter Slider2("Master Mix", "", 0.0, "", 0.0, 1.0);

Parameter Slider7("Master Gain", "", 0.5, "", 0.0, 1.0);

Parameter recButton("PrintOut Vals", "", 0.0, "", 0.0, 1.0);
Parameter track_backward("Change Carrier (backward)", "", 0, "", 0, 1);
Parameter track_forward("Change Carrier (forward)", "", 0, "", 0, 1);
Parameter marker_backward("Change Modulator (backward)", "", 0, "", 0, 1);
Parameter marker_forward("Change Modulator (forward)", "", 0, "", 0, 1);

Parameter reset_grains("Reset Grains", "", 0.0, "", 0.0, 1.0);
Parameter Slider16("Grain Voices", "", 1, "", 1, 20);
Parameter Slider17("Grain Duration", "", 1, "", 1, 300);
Parameter Slider18("Grain Ramp", "", 0, "", 0, 100);
Parameter Slider19("Grain Offset", "", 0, "", -149, 150);
Parameter Slider20("Grain Delay", "", 0, "", 0, 300);
Parameter Slider21("Grain Stretch", "", 1, "", 1, 5);
Parameter Slider22("Grain Randomness", "", 0., "", 0., 1.);

class MeterParams {
public:
	MeterParams() :
	mMasterGain(0),
	mRmsGain(1.0),
	mCrossFade(0.0),
	mCarrierFile(0),
	mModulatorFile(1)
	{}
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
		mTextureBuffer(((TOTAL_SPATIAL_SAMPLING * 3 * 4) + 1 ) * sizeof(float)),
		mSpriteTex(16,16, Graphics::LUMINANCE, Graphics::FLOAT),
		mStateMaker(defaultBroadcastIP())
		{
			addSphereWithTexcoords(mMesh, 1, SPATIAL_SAMPLING);
			mMesh.primitive(Graphics::POINTS);

			addSphereWithTexcoords(mWireframeMesh, 1, SPATIAL_SAMPLING);
			mWireframeMesh.primitive(Graphics::LINES);

			gaussianSprite(mSpriteTex);

			nav().pos(0,0,5);
			initWindow(Window::Dim(0,0, 600, 400), "Angkasa", VIS_FPS);

			// Find Ambisonics files
			SearchPaths sp;
			sp.addSearchPath(".", true);
			// 0 = carrier, i.e. source, 1 = modulator
			// Problem when source and modulator is the same file. Try to avoid!
			filename[0] = "norm_organ.wav";
			filename[1] = "norm_insects_long.wav";
			// filename[2] = "norm_insects.wav";
			filename[2] = "norm_gamelan_bapangSelisir.wav";
			// filename[4] = "norm_gulls.wav";
			// filename[5] = "norm_steamtrain.wav";
			// filename[6] = "norm_440_45-45.wav";
			// filename[7] = "norm_440_0-0_660_90-0.wav";
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

			mRmsSamples = mSoundFile[0]->frameRate() / VIS_FPS;
			// mRmsCounter = 0;
			mRmsAccum.resize(SPATIAL_SAMPLING * SPATIAL_SAMPLING);
			// TODO: use memset to write zeros instead
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

			currentFile = params.mCarrierFile;
			changeSample = false;
			framesRead[0]= 0;
			framesRead[1]= 0;
			channelWeights[0] = 2.0;
			channelWeights[1] = 4.0;
			channelWeights[2] = 4.0;
			channelWeights[3] = 2.0;

			// GRANULATOR
			preProjectedFile = new SoundFile(fullPath[0]);
			preProjectedFile->openRead();
			float preProjectedBuffer[preProjectedFile->samples()];
			preProjectedFile->read( preProjectedBuffer, preProjectedFile->samples() );

			g_N = 1;
			g_duration = 1;
			g_ramp = 0;
			g_offset = 0;
			g_delay = 0;
			g_stretch = 1;
			g_random = 0.;
			numSamplesInFile = mSoundFile[0]->frames();

			for (int i = 0; i < TOTAL_SPATIAL_SAMPLING; i++) {
				numSamplesReadFromFile[i] = 0;
				grani[i].setSampleRate(preProjectedFile->frameRate());
				grani[i].setRandomFactor(g_random);
				grani[i].setStretch(g_stretch);
				grani[i].setGrainParameters(g_duration, g_ramp, g_offset, g_delay);
				grani[i].openFile( preProjectedFile->frames() );
				grani[i].setVoices(g_N);
			}

			float sqrt2 = 1.0/ sqrt(2.0);
			for(int elevIndex = 0; elevIndex < SPATIAL_SAMPLING; elevIndex++) { // For sampled elevation
				float elev = M_PI/2.0 - (M_PI * (elevIndex + 0.5)/(float) SPATIAL_SAMPLING); // -1.51844 to 1.51844
				float sinElev = SIN(elev);
				float cosElev = COS(elev);
				for(int azimuthIndex = 0; azimuthIndex < SPATIAL_SAMPLING; azimuthIndex++) { // For each sampled azimuth
					float azimuth = M_2PI * (azimuthIndex + 0.5)/(float) SPATIAL_SAMPLING; // 0.10472 to 6.17847
					float sinAzi = SIN(azimuth);
					float cosAzi = COS(azimuth);
					int STgrain_index = azimuthIndex * SPATIAL_SAMPLING + elevIndex;

					float *w_pp = preProjectedBuffer;
					float *x_pp = preProjectedBuffer + 1;
					float *y_pp = preProjectedBuffer + 2;
					float *z_pp = preProjectedBuffer + 3;

					for (int frame =0; frame < preProjectedFile->frames(); frame++) {
						float projectedSignal = (*w_pp * sqrt2)
																	+ (*x_pp * cosAzi * cosElev)
																	+ (*y_pp * sinAzi * cosElev)
																	+ (*z_pp * sinElev);
						grani[STgrain_index].writeData(projectedSignal);
						w_pp+=4; x_pp+=4; y_pp+=4; z_pp+=4;
					}
				}
			}
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

			// Master mix for playing different Ambi files
			masterMix = params.mCrossFade;
			if(masterMix < 0.5){
				signal_mix[0] =  1 - (masterMix*2);
				signal_mix[1] = masterMix*2;
				signal_mix[2] = 0;
			} else if(masterMix > 0.5){
				signal_mix[0] = 0;
				signal_mix[1] = (1- masterMix)*2;
				signal_mix[2] = (masterMix*2) - 1;
			}

			// std::cout << "RMS GAIN " << params.mRmsGain
			// 					<< "RMS Threshold" << Slider1.get()
			// 					<< std::endl;


			std::cout << "Carrier: " << params.mCarrierFile
								<< "| Modulator: " << params.mModulatorFile
								<< "| Voices: " << g_N
								<< "| Duration: " << g_duration
								<< "| Ramp: " << g_ramp
								<< "| Offset: " << g_offset
								<< "| Delay: " << g_delay
								<< "| Stretch: " << g_stretch
								<< "| Random: " << g_random
								<< std::endl;
		}

		virtual void onDraw(Graphics &g) override {
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SMOOTH);
			glEnable(GL_POINT_SPRITE);
			glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
			g.blendAdd();
			// Draw wireframe mesh
			g.lineWidth(0.25);
			g.color(1, 1, 1, 0.5);
			g.draw(mWireframeMesh);

			// Read the texture buffer
			int bytesToRead = TOTAL_SPATIAL_SAMPLING * 3 * sizeof(float);
			while(mTextureBuffer.readSpace() >= bytesToRead) {
				mTextureBuffer.read(mMeterValues.data.ptr, bytesToRead);
				mTexture.submit(mMeterValues, true);
				// mTexture.filter(Texture::LINEAR);
			}

			ShaderProgram* s = shaderManager.get("point");
			s->begin();
			s->uniform("texture0", 0);
			s->uniform("texture1", 1);
			mTexture.bind(0);
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

			// Granulation parameters
			g_N 				= Slider16.get();
			g_duration 	= Slider17.get();
			g_ramp 			= Slider18.get();
			g_offset 		= Slider19.get();
			g_delay 		= Slider20.get();
			g_stretch 	= Slider21.get();
			g_random 		= Slider22.get();
			// std::cout << g_duration << std::endl;
			for (int i = 0; i< TOTAL_SPATIAL_SAMPLING; i++){
				grani[i].setVoices(g_N);
				grani[i].setRandomFactor(g_random);
				grani[i].setStretch(g_stretch);
				grani[i].setGrainParameters(g_duration, g_ramp, g_offset, g_delay);
			}

			// Set up spatializer
			spatializer->prepare();
			float * ambiChans = spatializer->ambiChans();
			memset(reconstructAmbiBuffer, 0, sizeof(float) * AUDIO_BLOCK_SIZE * 4);

			int numOfGatedGrains[AUDIO_BLOCK_SIZE];
			memset(numOfGatedGrains, 0, sizeof(float) * AUDIO_BLOCK_SIZE);

			// Read into buffer per CB
			framesRead[0] = mSoundFile[params.mCarrierFile]->read(readBuffer[0], AUDIO_BLOCK_SIZE);
			framesRead[1] = mSoundFile[params.mModulatorFile]->read(readBuffer[1], AUDIO_BLOCK_SIZE);
			if (framesRead[0] != AUDIO_BLOCK_SIZE || framesRead[1] != AUDIO_BLOCK_SIZE) {
				cout << "buffer overrun! framesRead[0]: " << framesRead[0] << " frames" << endl;
				cout << "buffer overrun! framesRead[1]: " << framesRead[1] << " frames" << endl;
			}

			// Pointer for output RMS, NxN
			float *rmsBuffer = (float *) mNewRmsMeterValues.data.ptr;

			float sqrt2 = 1.0/ sqrt(2.0);
			for(int elevIndex = 0; elevIndex < SPATIAL_SAMPLING; elevIndex++) { // For sampled elevation
				float elev = M_PI/2.0 - (M_PI * (elevIndex + 0.5)/(float) SPATIAL_SAMPLING); // -1.51844 to 1.51844
				float sinElev = SIN(elev);
				float cosElev = COS(elev);
				for(int azimuthIndex = 0; azimuthIndex < SPATIAL_SAMPLING; azimuthIndex++) { // For each sampled azimuth
					float azimuth = M_2PI * (azimuthIndex + 0.5)/(float) SPATIAL_SAMPLING; // 0.10472 to 6.17847
					float sinAzi = SIN(azimuth);
					float cosAzi = COS(azimuth);

					int STgrain_index = azimuthIndex * SPATIAL_SAMPLING + elevIndex;
					numSamplesReadFromFile[STgrain_index] = 0;
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

					grani[STgrain_index].setThreshold(params.mRmsGain, Slider1.get());

					for (int frame = 0; frame < AUDIO_BLOCK_SIZE; frame++) {

						if (params.mCarrierFile != currentFile){
							changeSample = true;
						}
						if(changeSample == true){
							if (numSamplesReadFromFile[STgrain_index] < numSamplesInFile){
								// cout << "changesample!" << endl;
								float STgrain_mix = (*w_0 * sqrt2)
															  	+ (*x_0 * COS(azimuth) * cosElev)
																  + (*y_0 * SIN(azimuth) * cosElev)
															  	+ (*z_0 * SIN(elev));

								grani[STgrain_index].writeData(STgrain_mix);
								numSamplesReadFromFile[STgrain_index]++;
							}
							else if (numSamplesReadFromFile[STgrain_index] >= numSamplesInFile){
								numSamplesReadFromFile[STgrain_index] = 0;
								changeSample = false;
								break;
								cout << "reset!" << endl;
							}
						}


						// Tick granulation engine
						STslice[STgrain_index] = grani[STgrain_index].tick();

						// Calculate number of grains that is above RMS threshold
						numOfGatedGrains[frame] += grani[STgrain_index].getNumOfGatedGrains();

						// Accumulate granulated output
						mRmsAccum[STgrain_index] += STslice[STgrain_index] * STslice[STgrain_index];

						// Project each grain on to corresponding direction
						float W_reconstruct = STslice[STgrain_index] * sqrt2;
						float X_reconstruct = STslice[STgrain_index] * cosAzi * cosElev;
						float Y_reconstruct = STslice[STgrain_index] * sinAzi * cosElev;
						float Z_reconstruct = STslice[STgrain_index] * sinElev;

						// Add each ST_grain
						*w_r += W_reconstruct;
						*x_r += X_reconstruct;
						*y_r += Y_reconstruct;
						*z_r += Z_reconstruct;

						// Increment pointer by 4 (4 channel interleave Ambi files)
						w_0+=4; x_0+=4; y_0+=4; z_0+=4;
						w_r+=4; x_r+=4; y_r+=4; z_r+=4;
						currentFile = params.mCarrierFile;
					} // End of audio frames

					mRms[STgrain_index] = sqrt(mRmsAccum[STgrain_index] / (float) AUDIO_BLOCK_SIZE);
					mRmsAccum[STgrain_index] = 0.0;

					*rmsBuffer++ = mRms[STgrain_index];
					*rmsBuffer++ = mRms[STgrain_index];
					*rmsBuffer++ = mRms[STgrain_index];
				} // End of Azi
			} // End of Elev

			// OUTPUT
			float * w_0 = readBuffer[0];
			float * w_r = reconstructAmbiBuffer;
			for (int frame = 0; frame < AUDIO_BLOCK_SIZE; frame++){
				// cout << "frame "<< frame << " has a total grain of " << numOfGatedGrains[frame] << endl;
				// AUDIO OUTPUT
				for (int chan = 0; chan < 4; chan++) { // For every Ambi Channel
						// reconstructAmbiBuffer[frame*4 + chan] /= (float)(TOTAL_SPATIAL_SAMPLING) / channelWeights[chan];
						reconstructAmbiBuffer[frame*4 + chan] /= (float)(numOfGatedGrains[frame]+1) / channelWeights[chan];
						ambiChans[frame+ chan*AUDIO_BLOCK_SIZE] = reconstructAmbiBuffer[frame*4 + chan]* params.mMasterGain;
					}

				// VISUAL OUTPUT
				mRmsCounter++;
				if (mRmsCounter == mRmsSamples) {
					mRmsCounter = 0;
					rmsBuffer = (float *) mNewRmsMeterValues.data.ptr;
					memcpy(mState.rmsTexture, rmsBuffer, TOTAL_SPATIAL_SAMPLING * 3 * sizeof(float));
					mStateMaker.set(mState);

					int bytesToWrite = TOTAL_SPATIAL_SAMPLING * 3 * sizeof(float);
					if (mTextureBuffer.writeSpace() >= bytesToWrite) {
						mTextureBuffer.write(mNewRmsMeterValues.data.ptr, bytesToWrite);
					} else {
						// cout << "Texture Buffer overrun!" << endl;
					}
				}

				// TEMP: Draw the waveforms
				originalWonly[frame] = *w_0;
				w_0 += 4;
			 	reconstructWonly[frame] = *w_r;
				w_r += 4;
			}

			spatializer->finalize(io);
		} // End of onSound

		virtual void onKeyDown(const Keyboard& k) override {
			if (k.key() == '1') {
				params.mCarrierFile = 0;
			} else if (k.key() == '2'){
				params.mCarrierFile = 1;
			} else if (k.key() == '3'){
				params.mCarrierFile = 2;
			} else if (k.key() == '4'){
				params.mCarrierFile = 3;
			} else if (k.key() == '5'){
				params.mCarrierFile = 4;
			} else if (k.key() == '6'){
				params.mCarrierFile = 5;
			} else if (k.key() == '7'){
				params.mCarrierFile = 6;
			} else if (k.key() == '8'){
				params.mCarrierFile = 7;

			} else if (k.key() == 'r'){
				params.mModulatorFile = 0;
			} else if (k.key() == 't'){
				params.mModulatorFile = 1;
			} else if (k.key() == 'y'){
				params.mModulatorFile = 2;
			} else if (k.key() == 'u'){
				params.mModulatorFile = 3;
			} else if (k.key() == 'i'){
				params.mModulatorFile = 4;
			} else if (k.key() == 'o'){
				params.mModulatorFile = 5;
			} else if (k.key() == 'p'){
				params.mModulatorFile = 6;
			} else if (k.key() == '['){
				params.mModulatorFile = 7;
			}
		}

	private:
		Mesh mMesh;
		Mesh mWireframeMesh;
		Texture mTexture;
		Texture mSpriteTex;
		al::Array mMeterValues;
		al::Array mNewRmsMeterValues;
		SingleRWRingBuffer mTextureBuffer;

		AmbisonicsSpatializer *spatializer;
		ShaderManager shaderManager;
		MeterParams params;

		int mRmsCounter;
		int mRmsSamples;
		vector<float> mRmsAccum;
		float mRms[SPATIAL_SAMPLING*SPATIAL_SAMPLING];

		SoundFileBuffered *mSoundFile[8];
		SoundFile *preProjectedFile;
		string filename[NUM_OF_FILES];
		string fullPath[NUM_OF_FILES];
		int framesRead[2];
		float channelWeights[4];

		// int numOfGatedGrains[AUDIO_BLOCK_SIZE];
		float STslice[TOTAL_SPATIAL_SAMPLING];
		float readBuffer[2][AUDIO_BLOCK_SIZE * 4];
		float reconstructAmbiBuffer[AUDIO_BLOCK_SIZE * 4];
		float reconstructWonly[AUDIO_BLOCK_SIZE];
		float originalWonly[AUDIO_BLOCK_SIZE];

		float masterMix;
		float signal_mix[3];
		float rmsThreshold;
		bool changeSample;
		int currentFile;
		int numSamplesInFile;
		int numSamplesReadFromFile[SPATIAL_SAMPLING*SPATIAL_SAMPLING];

		Granulate grani[TOTAL_SPATIAL_SAMPLING];
		int g_N;
		int g_duration;
		int g_ramp;
		int g_offset;
		int g_delay;
		int g_stretch;
		float g_random;

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

		// parameterMIDI.connectControl(Slider6, 6, 1);
		parameterMIDI.connectControl(Slider7, 7, 1);

		parameterMIDI.connectControl(recButton, 45, 1);
		parameterMIDI.connectControl(track_backward, 58, 1);
		parameterMIDI.connectControl(track_forward, 59, 1);
		parameterMIDI.connectControl(marker_backward, 61, 1);
		parameterMIDI.connectControl(marker_forward, 62, 1);

		parameterMIDI.connectControl(reset_grains, 46, 1);
		parameterMIDI.connectControl(Slider16, 16, 1);
		parameterMIDI.connectControl(Slider17, 17, 1);
		parameterMIDI.connectControl(Slider18, 18, 1);
		parameterMIDI.connectControl(Slider19, 19, 1);
		parameterMIDI.connectControl(Slider20, 20, 1);
		parameterMIDI.connectControl(Slider21, 21, 1);
		parameterMIDI.connectControl(Slider22, 22, 1);

		app.start();
	}
