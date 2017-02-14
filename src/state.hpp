#ifndef STATE_HPP
#define STATE_HPP


#define SPATIAL_SAMPLING 14//18

struct State {
	float rmsTexture[SPATIAL_SAMPLING*SPATIAL_SAMPLING*3];
	float peakTexture[SPATIAL_SAMPLING*SPATIAL_SAMPLING*3];
};



#endif // STATE_HPP
