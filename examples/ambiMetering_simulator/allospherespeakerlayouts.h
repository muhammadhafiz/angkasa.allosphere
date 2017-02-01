#ifndef ALLOSPHERESPEAKERLAYOUTS_H
#define ALLOSPHERESPEAKERLAYOUTS_H

#include "allocore/sound/al_AudioScene.hpp"
#include "alloutil/al_AlloSphereSpeakerLayout.hpp"


namespace  al {

class AllosphereSpeakerLayouts
{
public:
  AllosphereSpeakerLayouts();
  
  static SpeakerLayout threeRings54();
  static SpeakerLayout threeRings27();
  static SpeakerLayout centerRing();
  static SpeakerLayout topAndBottom();
  static SpeakerLayout sparse();
};
}

#endif // ALLOSPHERESPEAKERLAYOUTS_H
