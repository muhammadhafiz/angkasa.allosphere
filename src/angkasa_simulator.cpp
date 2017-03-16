// Angkasa: Spatiotemporal Granulation [Allosphere]
// Author: Muhammad Hafiz Wan Rosli
// January 1st 2017


#include "allocore/al_Allocore.hpp"
#include <allocore/types/al_SingleRWRingBuffer.hpp>
#include "Cuttlebone/Cuttlebone.hpp"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "alloGLV/al_ControlGLV.hpp"
#include "GLV/glv.h"

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
#define NUM_DREHBANK_KNOBS TOTAL_SPATIAL_SAMPLING
#define SIN sin
#define COS cos
// #define SIN gam::scl::sinT9
// #define COS gam::scl::cosT8


class MeterParams {
public:
	MeterParams() :
	mResetAll(0),
	mOverwriteSample(0),
	mMasterGain(0),
	mRmsGain(1.0),
	mRmsThreshold(0.0),
	mRmsScaler(1.0),
	mCrossFade(0.0),
	mCarrierFile(0),
	mModulatorFile(1)
	{}
		bool mResetAll;
		bool mOverwriteSample;
		float mMasterGain;
		float mRmsGain;
		float mRmsThreshold;
		float mRmsScaler;
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
		mSpriteTex(16, 16, Graphics::LUMINANCE, Graphics::FLOAT),
		mStateMaker(defaultBroadcastIP())
		{
			addSphereWithTexcoords(mMesh, 1, SPATIAL_SAMPLING);
			mMesh.primitive(Graphics::POINTS);

			addSphereWithTexcoords(mWireframeMesh, 1, SPATIAL_SAMPLING);
			mWireframeMesh.primitive(Graphics::LINES);

			gaussianSprite(mSpriteTex);

			nav().pos(0,0,6);
			initWindow(Window::Dim(0,0, 600, 400), "Angkasa", VIS_FPS);

			// Find Ambisonics files
			SearchPaths sp;
			sp.addSearchPath(".", true);
			// 0 = carrier, i.e. source, 1 = modulator
			// Problem when source and modulator is the same file. Try to avoid!
			// filename[0] = "norm_organ_fade.wav";
			filename[2] = "norm_organ.wav";
			filename[0] = "norm_insects_1min.wav";
			// filename[1] = "norm_insects_long.wav";
			// filename[1] = "norm_insects.wav";
			filename[1] = "norm_gamelan_bapangSelisir_long.wav";
			// filename[0] = "gambangan_30b.wav";
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
			mRmsAccum.resize(TOTAL_SPATIAL_SAMPLING);
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
			// changeSample = false;
			params.mOverwriteSample = false;
			framesRead[0]= 0;
			framesRead[1]= 0;
			framesRead[2]= 0;
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
			// g_random = 0.;
			g_randomDur = 0.;
			g_randomOffset = 0.;
			g_randomDelay = 0.;
			g_randomPointer = 0.;


			numSamplesInFile = mSoundFile[0]->frames();

			for (int i = 0; i < TOTAL_SPATIAL_SAMPLING; i++) {
				prevFile[i] = params.mCarrierFile;
				numSamplesReadFromFile[i] = 0;
				grani[i].setSampleRate(preProjectedFile->frameRate());
				// grani[i].setRandomFactor(g_random);
				grani[i].setRandomFactor(g_randomDur, g_randomDelay, g_randomOffset, g_randomPointer);
				grani[i].setStretch(g_stretch);
				grani[i].setGrainParameters(g_duration, g_ramp, g_offset, g_delay);
				grani[i].openFile( preProjectedFile->frames() );
				// cout << preProjectedFile->frames() << endl;
				grani[i].openFile( preProjectedFile->frameRate()*3 );
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

			// MIDI controller (Doepfer Drehbank)
			drehbank = new ParameterMIDI;
			int drehbank_midiPort = 0;
			if(sim()) drehbank_midiPort = 4;
			drehbank->init(drehbank_midiPort, false);

			// g_stretch = 17 -> 32 (16 knobs, 9 STG each)
			// g_offset = 33 -> 48 (16 knobs, 9 STG each)
			// g_delay = 49 -> 64 (16 knobs, 9 STG each)
			// ceil((i+1/9)+16)
			// knobsPerSTG = TOTAL_SPATIAL_SAMPLING/ 16; //9
			// knobsOffset = 16;

			// for (int i = 0; i< TOTAL_SPATIAL_SAMPLING; i++){
			// 	whichKnobSets = ceil(((float)i+1)/knobsPerSTG);
			// 	g_offset_knob = whichKnobSets + knobsOffset;
			// 	g_delay_knob = whichKnobSets + (knobsOffset * 2);
			// 	g_stretch_knob = whichKnobSets + (knobsOffset * 3);
			//
			// 	// cout << g_offset_knob << endl;
			// 	g_offset = drehbankKnob[g_offset_knob]->get();
			// 	g_delay = drehbankKnob[g_delay_knob]->get();
			// 	g_stretch = drehbankKnob[g_stretch_knob]->get();
			// 	//
			// 	grani[i].setStretch(g_stretch);
			// 	grani[i].setGrainParameters(g_duration, g_ramp, g_offset, g_delay);
			// }


			drehbankKnob[0] = new Parameter("Grain Voices", "", 1, "", 1, 10);
			drehbankKnob[1] = new Parameter("Grain Duration", "", 1, "", 1, 300);
			drehbankKnob[2] = new Parameter("Grain Ramp", "", 0, "", 0, 100);
			// drehbankKnob[3] = new Parameter("Grain Randomness", "", 0., "", 0., 0.1);
			drehbankKnob[4] = new Parameter("Grain Offset", "", 0, "", -149, 150);
			drehbankKnob[5] = new Parameter("Grain Delay", "", 0, "", 0, 300);
			drehbankKnob[6] = new Parameter("Grain Stretch", "", 1, "", 1, 10);
			drehbankKnob[8] = new Parameter("RMS Gain", "", 1.0, "", 0.0, 1.0);
			drehbankKnob[9] = new Parameter("RMS Threshold", "", 0.0, "", 0.0, 0.3);
			drehbankKnob[25] = new Parameter("RMS Scaler", "", 1.0, "", 0.0, 1.0);
			drehbankKnob[14] = new Parameter("Subwoofer Mix", "", 0.5, "", 0., 1.0);
			drehbankKnob[15] = new Parameter("Master Gain", "", 0.0, "", 0.0, 1.0);
			drehbankKnob[17] = new Parameter("Grain Randomness_dur", "", 0., "", 0., 1.0);
			drehbankKnob[20] = new Parameter("Grain Randomness_offset", "", 0., "", 0., 1.0);
			drehbankKnob[21] = new Parameter("Grain Randomness_delay", "", 0., "", 0., 1.0);
			drehbankKnob[22] = new Parameter("Grain Randomness_pointer", "", 0., "", 0., 1.0);

			drehbankKnob[46] = new Parameter("Overwrite Sample", "", 0., "", 0., 1.0); // If not 0.0
			drehbankKnob[47] = new Parameter("RESET POINTER", "", 0., "", 0., 1.0); // If not 0.0

			drehbankKnob[60] = new Parameter("Signal Mix 0", "", 0., "", 0., 1.0);
			drehbankKnob[61] = new Parameter("Signal Mix 1", "", 0., "", 0., 1.0);
			// drehbankKnob[62] = new Parameter("Signal Mix 2", "", 0., "", 0., 1.0);


			drehbank->connectControl(*drehbankKnob[0], 0, 1);
			drehbank->connectControl(*drehbankKnob[1], 1, 1);
			drehbank->connectControl(*drehbankKnob[2], 2, 1);
			// drehbank->connectControl(*drehbankKnob[3], 3, 1);
			drehbank->connectControl(*drehbankKnob[4], 4, 1);
			drehbank->connectControl(*drehbankKnob[5], 5, 1);
			drehbank->connectControl(*drehbankKnob[6], 6, 1);
			drehbank->connectControl(*drehbankKnob[8], 8, 1);
			drehbank->connectControl(*drehbankKnob[9], 9, 1);
			drehbank->connectControl(*drehbankKnob[25], 25, 1);
			drehbank->connectControl(*drehbankKnob[14], 14, 1);
			drehbank->connectControl(*drehbankKnob[15], 15, 1);
			drehbank->connectControl(*drehbankKnob[17], 17, 1);
			drehbank->connectControl(*drehbankKnob[20], 20, 1);
			drehbank->connectControl(*drehbankKnob[21], 21, 1);
			drehbank->connectControl(*drehbankKnob[22], 22, 1);
			drehbank->connectControl(*drehbankKnob[46], 46, 1);
			drehbank->connectControl(*drehbankKnob[47], 47, 1);
			drehbank->connectControl(*drehbankKnob[60], 60, 1);
			drehbank->connectControl(*drehbankKnob[61], 61, 1);
			// drehbank->connectControl(*drehbankKnob[62], 62, 1);
			// drehbank->connectControl(*drehbankKnob[63], 63, 1);


			// string pName_offset = "STG_offset";
			// string pName_delay = "STG_delay";
			// string pName_stretch = "STG_stretch";
			//
			// for(int i = 0; i< TOTAL_SPATIAL_SAMPLING; i++){
			// 	whichKnobSets = floor(((float)i)/knobsPerSTG);
			// 	g_offset_knob = whichKnobSets + knobsOffset;
			// 	g_delay_knob = whichKnobSets + (knobsOffset * 2);
			// 	g_stretch_knob = whichKnobSets + (knobsOffset * 3);
			//
			// 	pName_offset += std::to_string(i);
			// 	pName_delay += std::to_string(i);
			// 	pName_stretch += std::to_string(i);
			//
			// 	drehbankKnob[g_offset_knob] = new Parameter(pName_offset, "", 0, "", -149, 150);
			// 	drehbankKnob[g_delay_knob] = new Parameter(pName_delay, "", 0, "", 0, 300);
			// 	drehbankKnob[g_stretch_knob] = new Parameter(pName_stretch, "", 0, "", 1, 10);
			//
			// 	drehbank->connectControl(*drehbankKnob[g_offset_knob], g_offset_knob, 1);
			// 	drehbank->connectControl(*drehbankKnob[g_delay_knob], g_delay_knob, 1);
			// 	drehbank->connectControl(*drehbankKnob[g_stretch_knob], g_stretch_knob, 1);
			// }

			// MIDI controller (nanoKontrol2)
			// nanoKontrol2 = new ParameterMIDI;
			// int nanoKontrol2_midiPort = 0;
			// if(sim()) nanoKontrol2_midiPort = 3;
			// nanoKontrol2->init(nanoKontrol2_midiPort, false);
			//
			// Slider0 = new Parameter("RMS Gain", "", 1.0, "", 0.0, 1.0);
			// Slider1 = new Parameter("RMS Threshold", "", 0.0, "", 0.0, 0.3);
			// Slider2 = new Parameter("Master Mix", "", 0.0, "", 0.0, 1.0);
			// Slider7 = new Parameter("Master Gain", "", 0.5, "", 0.0, 1.0);
			// Slider16 = new Parameter("Grain Voices", "", 1, "", 1, 20);
			// Slider17 = new Parameter("Grain Duration", "", 1, "", 1, 300);
			// Slider18 = new Parameter("Grain Ramp", "", 0, "", 0, 100);
			// Slider19 = new Parameter("Grain Offset", "", 0, "", -149, 150);
			// Slider20 = new Parameter("Grain Delay", "", 0, "", 0, 300);
			// Slider21 = new Parameter("Grain Stretch", "", 1, "", 1, 5);
			// Slider22 = new Parameter("Grain Randomness", "", 0., "", 0., 0.1);
			//
			// recButton = new Parameter("PrintOut Vals", "", 0.0, "", 0.0, 1.0);
			// track_backward = new Parameter("Change Carrier (backward)", "", 0, "", 0, 1);
			// track_forward = new Parameter("Change Carrier (forward)", "", 0, "", 0, 1);
			// marker_backward = new Parameter("Change Modulator (backward)", "", 0, "", 0, 1);
			// marker_forward = new Parameter("Change Modulator (forward)", "", 0, "", 0, 1);
			//
			// nanoKontrol2->connectControl(*Slider0, 0, 1);
			// nanoKontrol2->connectControl(*Slider1, 1, 1);
			// nanoKontrol2->connectControl(*Slider2, 2, 1);
			// nanoKontrol2->connectControl(*Slider7, 7, 1);
			// nanoKontrol2->connectControl(*Slider16, 16, 1);
			// nanoKontrol2->connectControl(*Slider17, 17, 1);
			// nanoKontrol2->connectControl(*Slider18, 18, 1);
			// nanoKontrol2->connectControl(*Slider19, 19, 1);
			// nanoKontrol2->connectControl(*Slider20, 20, 1);
			// nanoKontrol2->connectControl(*Slider21, 21, 1);
			// nanoKontrol2->connectControl(*Slider22, 22, 1);
			// nanoKontrol2->connectControl(*recButton, 45, 1);
			// nanoKontrol2->connectControl(*track_backward, 58, 1);
			// nanoKontrol2->connectControl(*track_forward, 59, 1);
			// nanoKontrol2->connectControl(*marker_backward, 61, 1);
			// nanoKontrol2->connectControl(*marker_forward, 62, 1);

			// GUI
			// GLVBinding gui;
			// glv::Slider slider;
			// glv::Table layout;

			// cout << drehbankKnob[1]->getName() << endl;
			gui = new GLVBinding;
			layout[0] = new glv::Table;
			layout[1] = new glv::Table;

			// Connect GUI to window
			gui->bindTo(window());

			// Configure GUI
			// gui->stretch(1,100);
			// gui->style().color.set(glv::Color(0.7), 0.5);
			gui->style().color.set(glv::Color(0.75), 0.75);

			layout[0]->arrangement(">p");
			layout[1]->arrangement("d");

			// layout[0]->arrange();

			// for(int i = 0; i < 9; i++){
			// 	int index = i;
			// 	if(i == 3) continue;
			// 	else if(i > 6 && i < 9) index+=5;
			// 	// else if(i == 9) index+= 6;
			// 	// cout << index << endl;
			//
			// 	slider[index] = new glv::Slider;
			// 	slider[index]->interval(drehbankKnob[index]->min(), drehbankKnob[index]->max());
			// 	*layout[0] << *slider[index];
			// 	*layout[0] << new glv::Label(drehbankKnob[index]->getName());
			// 	// layout[0]->arrange();
			// }


			slider[19] = new glv::Slider;
			// slider[19]->interval(0, numSamplesInFile-1);
			slider[19]->interval(0, ( preProjectedFile->frameRate()*3 )-1);
			*layout[0] << *slider[19];
			*layout[0] << new glv::Label("Pointer Position");

			slider[0]= new glv::Slider;
			slider[0]->interval(drehbankKnob[0]->min(), drehbankKnob[0]->max());
			*layout[0] << *slider[0];
			*layout[0] << new glv::Label("Grain Voices");

			slider[1]= new glv::Slider;
			slider[1]->interval(drehbankKnob[1]->min(), drehbankKnob[1]->max());
			*layout[0] << *slider[1];
			*layout[0] << new glv::Label("Grain Duration");

			slider[2]= new glv::Slider;
			slider[2]->interval(drehbankKnob[2]->min(), drehbankKnob[2]->max());
			*layout[0] << *slider[2];
			*layout[0] << new glv::Label("Grain Ramp");

			slider[4]= new glv::Slider;
			slider[4]->interval(drehbankKnob[4]->min(), drehbankKnob[4]->max());
			*layout[0] << *slider[4];
			*layout[0] << new glv::Label("Grain Offset");

			slider[5]= new glv::Slider;
			slider[5]->interval(drehbankKnob[5]->min(), drehbankKnob[5]->max());
			*layout[0] << *slider[5];
			*layout[0] << new glv::Label("Grain Delay");

			slider[6]= new glv::Slider;
			slider[6]->interval(drehbankKnob[6]->min(), drehbankKnob[6]->max());
			*layout[0] << *slider[6];
			*layout[0] << new glv::Label("Grain Stretch");

			slider[17]= new glv::Slider;
			slider[17]->interval(drehbankKnob[17]->min(), drehbankKnob[17]->max());
			*layout[0] << *slider[17];
			*layout[0] << new glv::Label("Randomize Grain Dur");

			slider[20]= new glv::Slider;
			slider[20]->interval(drehbankKnob[20]->min(), drehbankKnob[20]->max());
			*layout[0] << *slider[20];
			*layout[0] << new glv::Label("Randomize Grain Offset");

			slider[21]= new glv::Slider;
			slider[21]->interval(drehbankKnob[21]->min(), drehbankKnob[21]->max());
			*layout[0] << *slider[21];
			*layout[0] << new glv::Label("Randomize Grain Delay");

			slider[22]= new glv::Slider;
			slider[22]->interval(drehbankKnob[22]->min(), drehbankKnob[22]->max());
			*layout[0] << *slider[22];
			*layout[0] << new glv::Label("Randomize Grain Pointer");

			slider[8]= new glv::Slider;
			slider[8]->interval(drehbankKnob[8]->min(), drehbankKnob[8]->max());
			*layout[0] << *slider[8];
			*layout[0] << new glv::Label("RMS Gain");

			slider[9]= new glv::Slider;
			slider[9]->interval(drehbankKnob[9]->min(), drehbankKnob[9]->max());
			*layout[0] << *slider[9];
			*layout[0] << new glv::Label("RMS Threshold");

			slider[25]= new glv::Slider;
			slider[25]->interval(drehbankKnob[25]->min(), drehbankKnob[25]->max());
			*layout[0] << *slider[25];
			*layout[0] << new glv::Label("RMS Scaler");

			slider[60]= new glv::Slider;
			slider[60]->interval(drehbankKnob[60]->min(), drehbankKnob[60]->max());
			*layout[0] << *slider[60];
			*layout[0] << new glv::Label("Signal 0 Mix");

			slider[61]= new glv::Slider;
			slider[61]->interval(drehbankKnob[61]->min(), drehbankKnob[61]->max());
			*layout[0] << *slider[61];
			*layout[0] << new glv::Label("Signal 1 Mix");

			slider[14]= new glv::Slider;
			slider[14]->interval(drehbankKnob[14]->min(), drehbankKnob[14]->max());
			*layout[0] << *slider[14];
			*layout[0] << new glv::Label("Subwoofer Mix");

			slider[15]= new glv::Slider;
			slider[15]->interval(drehbankKnob[15]->min(), drehbankKnob[15]->max());
			*layout[0] << *slider[15];
			*layout[0] << new glv::Label("Master Gain");


			// slider[29]= new glv::Slider;
			// slider[29]->interval(drehbankKnob[29]->min(), drehbankKnob[29]->max());
			// *layout[0] << *slider[29];
			// *layout[0] << new glv::Label("RMS Scaler");
			//
			// slider[15]= new glv::Slider;
			// slider[15]->interval(drehbankKnob[15]->min(), drehbankKnob[15]->max());
			// *layout[0] << *slider[15];
			// *layout[0] << new glv::Label(drehbankKnob[15]->getName());
			//
			//
			// slider[60] = new glv::Slider;
			// slider[61] = new glv::Slider;
			// slider[62] = new glv::Slider;
			// slider[63] = new glv::Slider;
			// slider[60]->interval(drehbankKnob[60]->min(), drehbankKnob[60]->max());
			// slider[61]->interval(drehbankKnob[61]->min(), drehbankKnob[61]->max());
			// slider[62]->interval(drehbankKnob[62]->min(), drehbankKnob[62]->max());
			// slider[63]->interval(drehbankKnob[63]->min(), drehbankKnob[63]->max());
			// *layout[0] << *slider[60];
			// *layout[0] << new glv::Label("Signal 0");
			// *layout[0] << *slider[61];
			// *layout[0] << new glv::Label("Signal 1");
			// *layout[0] << *slider[62];
			// *layout[0] << new glv::Label("Signal 2");
			// *layout[0] << *slider[63];
			// *layout[0] << new glv::Label("Subwoofer");
			// layout[0]->arrange();
			// *gui << *layout[0];

			// nd = new glv::NumberDialer;
			// *layout[0] << *nd;

			double plotWidth = 400;
			// data.resize(glv::Data::DOUBLE, 4, AUDIO_BLOCK_SIZE);
			data.resize(glv::Data::DOUBLE, 4, AUDIO_BLOCK_SIZE);
			// for (int i = 0; i< AUDIO_BLOCK_SIZE; i++){
			// 	// data.assign(1, 0, i, i+1);
			// 	// data.assign(0, 1, i, i+1);
			// 	data.assign(0.9, 0, i);
			// 	data.assign(0.25, 1, i);
			// 	data.assign(-0.25, 2, i);
			// 	data.assign(-0.9, 3, i);
			// }

			plotFunction1D[0] = new glv::PlotFunction1D(glv::Color(0.75, 0.75, 0.75, 1.0));
			plotFunction1D[1] = new glv::PlotFunction1D(glv::Color(0.75, 0.00, 0.00, 1.0));
			plotFunction1D[2] = new glv::PlotFunction1D(glv::Color(0.00, 0.75, 0.00, 1.0));
			plotFunction1D[3] = new glv::PlotFunction1D(glv::Color(0.25, 0.00, 1.00, 1.0));

			plotTest[0] = new glv::Plot(glv::Rect(    0, 0* plotWidth/4, plotWidth,  plotWidth/4), *plotFunction1D[0]);
			plotTest[1] = new glv::Plot(glv::Rect(    0, 0* plotWidth/4, plotWidth,  plotWidth/4), *plotFunction1D[1]);
			plotTest[2] = new glv::Plot(glv::Rect(    0, 0* plotWidth/4, plotWidth,  plotWidth/4), *plotFunction1D[2]);
			plotTest[3] = new glv::Plot(glv::Rect(    0, 0* plotWidth/4, plotWidth,  plotWidth/4), *plotFunction1D[3]);
			// glv::Plot v1__(glv::Rect(    0,0*d/8, d,  d/8), *plotFunction1D[0]);
			// glv::Plot v1__(glv::Rect(100,20), *new glv::PlotFunction1D(glv::Color(0.5,0,0)));
			plotTest[0]->data() = data.slice(0).stride(data.size(0)).shape(1, data.size(1,2));
			plotTest[1]->data() = data.slice(1).stride(data.size(0)).shape(1, data.size(1,2));
			plotTest[2]->data() = data.slice(2).stride(data.size(0)).shape(1, data.size(1,2));
			plotTest[3]->data() = data.slice(3).stride(data.size(0)).shape(1, data.size(1,2));

			// plotTest[0]->major(data.size()/64);
			// plotTest[1]->major(data.size()/64);
			// plotTest[2]->major(data.size()/64);
			// plotTest[3]->major(data.size()/64);
			// cout << data.size() << endl;
			plotTest[0]->minor(AUDIO_BLOCK_SIZE/ 4, 0);
			plotTest[1]->minor(AUDIO_BLOCK_SIZE/ 4, 0);
			plotTest[2]->minor(AUDIO_BLOCK_SIZE/ 4, 0);
			plotTest[3]->minor(AUDIO_BLOCK_SIZE/ 4, 0);

			plotTest[0]->major(AUDIO_BLOCK_SIZE, 0);
			plotTest[1]->major(AUDIO_BLOCK_SIZE, 0);
			plotTest[2]->major(AUDIO_BLOCK_SIZE, 0);
			plotTest[3]->major(AUDIO_BLOCK_SIZE, 0);

			// plotTest[0]->showGrid(false, 1);
			// plotTest[1]->showGrid(false, 1);
			// plotTest[2]->showGrid(false, 1);
			// plotTest[3]->showGrid(false, 1);

			plotTest[0]->showAxis(false, -1);
			plotTest[1]->showAxis(false, -1);
			plotTest[2]->showAxis(false, -1);
			plotTest[3]->showAxis(false, -1);

			plotTest[0]->range(0, data.size(1,2), 0).range(-1.2, 1.2, 1);
			plotTest[1]->range(0, data.size(1,2), 0).range(-1.2, 1.2, 1);
			plotTest[2]->range(0, data.size(1,2), 0).range(-1.2, 1.2, 1);
			plotTest[3]->range(0, data.size(1,2), 0).range(-1.2, 1.2, 1);

			plotTest[0]->lockScroll(true, -1).lockZoom(true, -1);
			plotTest[1]->lockScroll(true, -1).lockZoom(true, -1);
			plotTest[2]->lockScroll(true, -1).lockZoom(true, -1);
			plotTest[3]->lockScroll(true, -1).lockZoom(true, -1);

			// topTest = new glv::GLV;

			// *layout[0] << *plotFunction1D[0];
			// *layout[0] << *plotTest[0];
			// *layout[0] << *plotTest[1];
			// *layout[0] << *plotTest[2];
			// *layout[0] << *plotTest[3];
			// layout[1]->arrangement(">p");
			*layout[1] << *plotTest[0] << *plotTest[1] << *plotTest[2] << *plotTest[3];

			// button[0] = new glv::Button(glv::Rect(100));
			button[0] = new glv::Button(glv::Rect(20, 20), false);
			button[1] = new glv::Button(glv::Rect(20, 20), false);
			*layout[1] << *button[0];
			*layout[1] << *button[1];

			layout[1]->pos(glv::Place::TR).anchor(0.787, 0);
			layout[0]->arrange();
			layout[1]->arrange();

			*gui << *layout[0];
			*gui << *layout[1];
			// *gui << *button[0] << *button[1];
			// button[0].pos();

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

			// if (track_backward->get() || track_forward->get()){
			// 	params.mCarrierFile-= track_backward->get();
			// 	params.mCarrierFile+= track_forward->get();
			// 	if (params.mCarrierFile > NUM_OF_FILES-1){
			// 		params.mCarrierFile = 0;
			// 	}else if (params.mCarrierFile < 0){
			// 		params.mCarrierFile = NUM_OF_FILES-1;
			// 	} cout << "Carrier: " << params.mCarrierFile << endl;
			// }
			//
			// if (marker_backward->get() || marker_forward->get()){
			// 	params.mModulatorFile-= marker_backward->get();
			// 	params.mModulatorFile+= marker_forward->get();
			// 	if (params.mModulatorFile > NUM_OF_FILES-1){
			// 		params.mModulatorFile = 0;
			// 	}else if (params.mModulatorFile < 0){
			// 		params.mModulatorFile = NUM_OF_FILES-1;
			// 	} cout << "Modulator: " << params.mModulatorFile << endl;
			// }

			// Master mix for playing different Ambi files
			// masterMix = params.mCrossFade;
			// if(masterMix < 0.5){
			// 	signal_mix[0] =  1 - (masterMix*2);
			// 	signal_mix[1] = masterMix*2;
			// 	signal_mix[2] = 0;
			// } else if(masterMix > 0.5){
			// 	signal_mix[0] = 0;
			// 	signal_mix[1] = (1- masterMix)*2;
			// 	signal_mix[2] = (masterMix*2) - 1;
			// }

			// std::cout << "RMS GAIN " << params.mRmsGain
			// 					<< "RMS Threshold" << params.mRmsThreshold
			// 					<< std::endl;


			// std::cout << "Car: " << params.mCarrierFile
			// 					<< "| V: " << g_N
			// 					<< "| Du: " << g_duration
			// 					<< "| R_Du: " << g_randomDur
			// 					<< "| R: " << g_ramp
			// 					<< "| O: " << g_offset
			// 					<< "| R_O: " << g_randomOffset
			// 					<< "| De: " << g_delay
			// 					<< "| R_De: " << g_randomDelay
			// 					<< "| S: " << g_stretch
			// 					<< "| R_P: " << g_randomPointer
			// 					// << "| Random: " << g_random
			// 					<< std::endl;

			// std::cout << "Carrier: " << params.mCarrierFile
			// 					<< "| Voices: " << g_N
			// 					<< "| Duration: " << g_duration
			// 					<< "| Ramp: " << g_ramp
			// 					// << "| Offset: " << g_offset
			// 					// << "| Delay: " << g_delay
			// 					// << "| Stretch: " << g_stretch
			// 					// << "| Random: " << g_random
			// 					<< std::endl;
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
			// Mesh waveform0;
			// Mesh waveform1;
			// addWaveformDisplay(waveform0, originalWonly, 1, 0, 0);
			// addWaveformDisplay(waveform1, reconstructWonly, 0, 1, 0);
			// g.translate(1.0, -1.0, 0);
			// g.draw(waveform0);
			// g.draw(waveform1);

			// GUI
			// What the slider shows
			slider[0]->setValue(g_N);
			slider[1]->setValue(g_duration);
			slider[2]->setValue(g_ramp);
			slider[4]->setValue(g_offset);
			slider[5]->setValue(g_delay);
			slider[6]->setValue(g_stretch);
			slider[8]->setValue(params.mRmsGain);
			slider[9]->setValue(params.mRmsThreshold);
			slider[14]->setValue(sub_mix);
			slider[15]->setValue(params.mMasterGain);
			slider[19]->setValue(grani[0].getGrainPointer());
			slider[25]->setValue(params.mRmsScaler);
			slider[60]->setValue(signal_mix[0]);
			slider[61]->setValue(signal_mix[1]);

			slider[17]->setValue(g_randomDur);
			slider[20]->setValue(g_randomOffset);
			slider[21]->setValue(g_randomDelay);
			slider[22]->setValue(g_randomPointer);

			button[0]->setValue(params.mOverwriteSample);
			button[1]->setValue(params.mResetAll);

			for (int i = 0; i< AUDIO_BLOCK_SIZE; i++){
				for (int j= 0; j< 4; j++){
				data.assign(gui_reconstruct[(i*4)+j], j, i);
				}
			}
			// nd->setValue(params.mCarrierFile);
			// cout << grani[0].getGrainPointer() << endl;

		}

		virtual void onSound(AudioIOData &io) override{

			// Set up spatializer
			spatializer->prepare();
			float * ambiChans = spatializer->ambiChans();
			memset(reconstructAmbiBuffer, 0, sizeof(float) * AUDIO_BLOCK_SIZE * 4);

			int numOfGatedGrains[AUDIO_BLOCK_SIZE];
			memset(numOfGatedGrains, 0, sizeof(float) * AUDIO_BLOCK_SIZE);

			memset(subChan, 0, sizeof(float) * AUDIO_BLOCK_SIZE);
			// Read into buffer per CB
			// framesRead[0] = mSoundFile[params.mCarrierFile]->read(readBuffer[0], AUDIO_BLOCK_SIZE);
			// framesRead[1] = mSoundFile[params.mModulatorFile]->read(readBuffer[1], AUDIO_BLOCK_SIZE);
			// if (framesRead[0] != AUDIO_BLOCK_SIZE || framesRead[1] != AUDIO_BLOCK_SIZE) {
			// 	cout << "buffer overrun! framesRead[0]: " << framesRead[0] << " frames" << endl;
			// 	cout << "buffer overrun! framesRead[1]: " << framesRead[1] << " frames" << endl;
			// }

			// Read into buffer per CB
			// framesRead[0] = mSoundFile[params.mCarrierFile]->read(readBuffer[0], AUDIO_BLOCK_SIZE);
			framesRead[0] = mSoundFile[0]->read(readBuffer[0], AUDIO_BLOCK_SIZE);
			framesRead[1] = mSoundFile[1]->read(readBuffer[1], AUDIO_BLOCK_SIZE);
			framesRead[2] = mSoundFile[2]->read(readBuffer[2], AUDIO_BLOCK_SIZE);

			// params.mRmsGain = Slider0->get();
			// params.mRmsThreshold = Slider1->get();
			// params.mMasterGain = Slider7->get();
			//


			signal_mix[0] = drehbankKnob[60]->get();
			signal_mix[1] = drehbankKnob[61]->get();
			// signal_mix[2] = drehbankKnob[62]->get();
			sub_mix = drehbankKnob[14]->get();
			// // Granulation parameters
			// g_N 				= Slider16->get();
			// g_duration 	= Slider17->get();
			// g_ramp 			= Slider18->get();
			// g_offset 		= Slider19->get();
			// g_delay 		= Slider20->get();
			// g_stretch 	= Slider21->get();
			// g_random 		= Slider22->get();

			params.mRmsGain = drehbankKnob[8]->get();
			params.mRmsThreshold = drehbankKnob[9]->get();
			params.mRmsScaler = drehbankKnob[25]->get();
			params.mMasterGain = drehbankKnob[15]->get();


			if (drehbankKnob[46]->get() != 0.0) {
				params.mOverwriteSample = true;
			} else if (drehbankKnob[46]->get() == 0.0) {
				params.mOverwriteSample = false;
			}

			if (drehbankKnob[47]->get() != 0.0) {
				params.mResetAll = true;
			} else if (drehbankKnob[47]->get() == 0.0) {
				params.mResetAll = false;
			}


			// Granulation parameters
			// g_N 				= Slider16->get();
			g_N 				= drehbankKnob[0]->get();
			g_duration 	= drehbankKnob[1]->get();
			g_randomDur = drehbankKnob[17]->get();
			g_ramp 			= drehbankKnob[2]->get();
			// g_random 		= drehbankKnob[3]->get();
			g_offset 		= drehbankKnob[4]->get();
			g_randomOffset = drehbankKnob[20]->get();
			g_delay 		= drehbankKnob[5]->get();
			g_randomDelay = drehbankKnob[21]->get();
			g_stretch 	= drehbankKnob[6]->get();
			g_randomPointer = drehbankKnob[22]->get();

			for (int i = 0; i< TOTAL_SPATIAL_SAMPLING; i++){
				grani[i].setRandomFactor(g_randomDur, g_randomDelay, g_randomOffset, g_randomPointer);
				grani[i].setStretch(g_stretch);
				grani[i].setGrainParameters(g_duration, g_ramp, g_offset, g_delay);
				grani[i].setVoices(g_N);
			}


			// for (int i = 0; i< TOTAL_SPATIAL_SAMPLING; i++){
			// 	whichKnobSets = floor(((float)i)/knobsPerSTG);
			// 	g_offset_knob = whichKnobSets + knobsOffset;
			// 	g_delay_knob = whichKnobSets + (knobsOffset * 2);
			// 	g_stretch_knob = whichKnobSets + (knobsOffset * 3);
			//
			// 	g_offset = drehbankKnob[g_offset_knob]->get();
			// 	g_delay = drehbankKnob[g_delay_knob]->get();
			// 	g_stretch = drehbankKnob[g_stretch_knob]->get();
			// 	// cout << g_stretch_knob << endl;
			// 	grani[i].setVoices(g_N);
			// 	grani[i].setRandomFactor(g_random);
			// 	grani[i].setStretch(g_stretch);
			// 	grani[i].setGrainParameters(g_duration, g_ramp, g_offset, g_delay);
			// }

			for (int frame = 0; frame < AUDIO_BLOCK_SIZE*4; frame++) {
				mixedReadBuffer[frame] = (readBuffer[0][frame]* signal_mix[0]) + (readBuffer[1][frame]* signal_mix[1]) + (readBuffer[2][frame]* signal_mix[2]);
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
					// numSamplesReadFromFile[STgrain_index] = 0;
					float *w_0 = mixedReadBuffer;
					float *x_0 = mixedReadBuffer + 1;
					float *y_0 = mixedReadBuffer + 2;
					float *z_0 = mixedReadBuffer + 3;

					// // // Pointer for buffer 0
					// float *w_0 = readBuffer[0];
					// float *x_0 = readBuffer[0] + 1;
					// float *y_0 = readBuffer[0] + 2;
					// float *z_0 = readBuffer[0] + 3;

					// // Pointer for buffer 1
					// float *w_1 = readBuffer[1];
					// float *x_1 = readBuffer[1] + 1;
					// float *y_1 = readBuffer[1] + 2;
					// float *z_1 = readBuffer[1] + 3;

					// Pointer for reconstructed buffer
					float *w_r = reconstructAmbiBuffer;
					float *x_r = reconstructAmbiBuffer + 1;
					float *y_r = reconstructAmbiBuffer + 2;
					float *z_r = reconstructAmbiBuffer + 3;

					grani[STgrain_index].setThreshold(params.mRmsGain, params.mRmsThreshold * params.mRmsScaler);

					for (int frame = 0; frame < AUDIO_BLOCK_SIZE; frame++) {
						if (params.mOverwriteSample == true){
							float STgrain_mix = (*w_0 * sqrt2)
																+ (*x_0 * cosAzi * cosElev)
																+ (*y_0 * sinAzi * cosElev)
																+ (*z_0 * sinElev);

							grani[STgrain_index].writeData(STgrain_mix);
							numSamplesReadFromFile[STgrain_index]++;

							if (numSamplesReadFromFile[STgrain_index] >= numSamplesInFile){
									numSamplesReadFromFile[STgrain_index] = 0;
							}
						}

						if (params.mResetAll == true){
							grani[STgrain_index].resetGlobalPointer();
							// grani[STgrain_index].reset();
							// cout << "RESET!" << endl;
						}

						// Tick granulation engine
						STslice[STgrain_index] = grani[STgrain_index].tick();

						// Calculate number of grains that is above RMS threshold for normalization
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
						// currentFile = params.mCarrierFile;
					} // End of audio frames

					mRms[STgrain_index] = sqrt(mRmsAccum[STgrain_index] / (float) AUDIO_BLOCK_SIZE);
					mRmsAccum[STgrain_index] = 0.0;

					*rmsBuffer++ = mRms[STgrain_index];
					*rmsBuffer++ = mRms[STgrain_index];
					*rmsBuffer++ = mRms[STgrain_index];
				} // End of Azi
			} // End of Elev
			params.mResetAll = false;

			// OUTPUT
			// float * w_0 = readBuffer[0];
			// float * w_r = reconstructAmbiBuffer;
			// Pointer for reconstructed buffer
			// float *gui_w_r = reconstructAmbiBuffer;
			// float *gui_x_r = reconstructAmbiBuffer + 1;
			// float *gui_y_r = reconstructAmbiBuffer + 2;
			// float *gui_z_r = reconstructAmbiBuffer + 3;
			for (int frame = 0; frame < AUDIO_BLOCK_SIZE; frame++){
				// cout << "frame "<< frame << " has a total grain of " << numOfGatedGrains[frame] << endl;
				// AUDIO OUTPUT
				for (int chan = 0; chan < 4; chan++) { // For every Ambi Channel
						// reconstructAmbiBuffer[frame*4 + chan] /= (float)(TOTAL_SPATIAL_SAMPLING) / channelWeights[chan];
						reconstructAmbiBuffer[frame*4 + chan] /= (float)(numOfGatedGrains[frame]+1) / channelWeights[chan];
						ambiChans[frame+ chan*AUDIO_BLOCK_SIZE] = reconstructAmbiBuffer[frame*4 + chan]* params.mMasterGain;
						subChan[frame] += ambiChans[frame+ chan*AUDIO_BLOCK_SIZE];
					}

				if(sim()){
					// subChan[frame] /= 4;
					float subValue = subChan[frame] * sub_mix;
					io.out(47, frame) = subValue;
					// subChan[frame] = 0;
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
				// cout << frame << endl;
				// for (int i = 0; i< AUDIO_BLOCK_SIZE; i++){
				// 	// data.assign(1, 0, i, i+1);
				// 	// data.assign(0, 1, i, i+1);

				// glv::Data& d = data;
				// d.assign(*gui_w_r, 0, frame);
				// d.assign(*gui_w_r, 1, frame);
				// d.assign(*gui_w_r, 2, frame);
				// d.assign(*gui_w_r, 3, frame);
				// }

				// TEMP: Draw the waveforms
				// originalWonly[frame] = w_0;
				// w_0 += 4;
				// gui_reconstruct[frame] = *gui_w_r;
				// gui_reconstruct[frame+1] = *gui_x_r;
				// gui_reconstruct[frame+2] = *gui_y_r;
				// gui_reconstruct[frame+3] = *gui_z_r;
				// gui_w_r += 4; gui_x_r +=4; gui_y_r += 4; gui_z_r +=4;
			}
			memcpy(gui_reconstruct, reconstructAmbiBuffer, sizeof(float) * AUDIO_BLOCK_SIZE * 4);
			spatializer->finalize(io);
		} // End of onSound

		virtual void onKeyDown(const Keyboard& k) override {
			// k.print();
			// if (k.key() == 8) {
			// 	params.mOverwriteSample = !params.mOverwriteSample;
			// } else if (k.key() == 92){
			// 	params.mResetAll = true;
			// 	// params.mResetAll = !params.mResetAll;
			// 	// params.mResetAll = false;
			// }

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

		// bool mOverwriteSample;
		int mRmsCounter;
		int mRmsSamples;
		vector<float> mRmsAccum;
		float mRms[SPATIAL_SAMPLING*SPATIAL_SAMPLING];

		SoundFileBuffered *mSoundFile[8];
		SoundFile *preProjectedFile;
		string filename[NUM_OF_FILES];
		string fullPath[NUM_OF_FILES];
		int framesRead[3];
		float channelWeights[4];

		// int numOfGatedGrains[AUDIO_BLOCK_SIZE];
		float STslice[TOTAL_SPATIAL_SAMPLING];
		float mixedReadBuffer[AUDIO_BLOCK_SIZE * 4];
		float readBuffer[3][AUDIO_BLOCK_SIZE * 4];
		float reconstructAmbiBuffer[AUDIO_BLOCK_SIZE * 4];
		float gui_reconstruct[AUDIO_BLOCK_SIZE * 4];
		float subChan[AUDIO_BLOCK_SIZE];
		// float reconstructWonly[AUDIO_BLOCK_SIZE];
		// float originalWonly[AUDIO_BLOCK_SIZE];

		float masterMix;
		float sub_mix;
		float signal_mix[3];
		float rmsThreshold;
		// bool changeSample;
		int numSamplesInFile;
		int currentFile;
		int prevFile[SPATIAL_SAMPLING*SPATIAL_SAMPLING];
		int numSamplesReadFromFile[SPATIAL_SAMPLING*SPATIAL_SAMPLING];

		Granulate grani[TOTAL_SPATIAL_SAMPLING];
		int g_N;
		int g_duration;
		int g_ramp;
		int g_offset;
		int g_delay;
		int g_stretch;
		// float g_random;
		float g_randomDur;
		float g_randomDelay;
		float g_randomOffset;
		float g_randomPointer;

		State mState;
		cuttlebone::Maker<State> mStateMaker;

		// MIDI controller (Doepfer Drehbank)
		int knobsPerSTG;
		int knobsOffset;
		int whichKnobSets;
		int g_offset_knob;
		int g_delay_knob;
		int g_stretch_knob;

		ParameterMIDI *drehbank;
		Parameter *drehbankKnob[NUM_DREHBANK_KNOBS];

		// // MIDI controller (nanoKontrol2)
		// ParameterMIDI *nanoKontrol2; // KORG nanoKONTROL2
		// Parameter *Slider0;
		// Parameter *Slider1;
		// Parameter *Slider2;
		// Parameter *Slider7;
		// Parameter *Slider16;
		// Parameter *Slider17;
		// Parameter *Slider18;
		// Parameter *Slider19;
		// Parameter *Slider20;
		// Parameter *Slider21;
		// Parameter *Slider22;
		// Parameter *recButton;

		// Parameter *track_backward;
		// Parameter *track_forward;
		// Parameter *marker_backward;
		// Parameter *marker_forward;

		// GUI
		GLVBinding *gui;
		glv::Slider *slider[NUM_DREHBANK_KNOBS];
		glv::Button *button[2];
		glv::NumberDialer *nd;
		glv::Table *layout[2];
		glv::PlotFunction1D *plotFunction1D[4];
		glv::Plot *plotTest[4];
		glv::Data data;
	};


	int main(int argc, char *argv[])
	{

		Angkasa app;

		cout 	<< "\nControls (Doepfer Drehbank): " << endl
					<< "Knob[0]: Grain Voices" << endl
					<< "Knob[1]: Grain Duration" << endl
					<< "Knob[2]: Grain Ramp" << endl
					<< "Knob[4]: Grain Offset" << endl
					<< "Knob[5]: Grain Delay" << endl
					<< "Knob[6]: Grain Stretch" << endl
					<< "Knob[8]: RMS Gain" << endl
					<< "Knob[9]: RMS Threshold" << endl
					<< "Knob[25]: RMS Scaler" << endl
					<< "Knob[14]: Subwoofer Mix" << endl
					<< "Knob[15]: Master Gain" << endl
					<< "Knob[17]: Grain Randomness_dur" << endl
					<< "Knob[20]: Grain Randomness_offset" << endl
					<< "Knob[21]: Grain Randomness_delay" << endl
					<< "Knob[22]: Grain Randomness_pointer" << endl
					<< "Knob[46]: Overwrite Sample" << endl
					<< "Knob[47]: RESET POINTER" << endl
					<< "Knob[60]: Signal Mix 0" << endl
					<< "Knob[61]: Signal Mix 1" << endl;

		app.start();

	}
