#ifndef STATE_HPP
#define STATE_HPP

#define AUDIO_BLOCK_SIZE 1024 // default: 512
#define SPATIAL_SAMPLING 15//15//30//12//18
#define VIS_FPS 60


struct State {
	float rmsTexture[SPATIAL_SAMPLING*SPATIAL_SAMPLING*3];
};



#endif // STATE_HPP
