#ifndef STATE_HPP
#define STATE_HPP


#define SPATIAL_SAMPLING 10//15//30//12//18

struct State {
	float rmsTexture[SPATIAL_SAMPLING*SPATIAL_SAMPLING*3];
	// float peakTexture[SPATIAL_SAMPLING*SPATIAL_SAMPLING*3];
};



#endif // STATE_HPP
