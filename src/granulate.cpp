#include "granulate.h"
// #include "FileRead.h"
#include <cmath>


Granulate :: Granulate( void )
{
  this->setGrainParameters(); // use default values
  this->setRandomFactor();
  gStretch_ = 0;
  stretchCounter_ = 0;
  gain_ = 1.0;
}

Granulate :: Granulate( unsigned int nVoices, int projectedBufferSize )
{
  this->setGrainParameters(); // use default values
  this->setRandomFactor();
  gStretch_ = 0;
  stretchCounter_ = 0;
  this->openFile( projectedBufferSize );
  this->setVoices( nVoices );
}

Granulate :: ~Granulate( void )
{
}

void Granulate :: setSampleRate( const float sr )
{
  sampleRate_ = sr;

  if(sr != sampleRate_ ) {
    std::cout << "Different Sample Rate!!!" << std::endl;
    // recompute stuff based on the new sample rate
  }
}

void Granulate :: writeData( double ringBufferValue ) {
  ring_(ringBufferValue);
}

int Granulate :: getNumOfGatedGrains() {
  return gNumOfGatedActiveGrains_;
}

void Granulate :: setThreshold ( float rmsGain, float rmsThreshold ) {
  gRmsGain_ = rmsGain;
  gRmsThreshold_ = rmsThreshold;
}

void Granulate :: setStretch( unsigned int stretchFactor )
{
  if ( stretchFactor <= 1 )
    gStretch_ = 0;
  else if ( gStretch_ >= 1000 )
    gStretch_ = 1000;
  else
    gStretch_ = stretchFactor - 1;
}

void Granulate :: setGrainParameters( unsigned int duration, unsigned int rampPercent,
                                      int offset, unsigned int delay )
{
  gDuration_ = duration;
  if ( gDuration_ == 0 ) {
    gDuration_ = 1;
    // oStream_ << "Granulate::setGrainParameters: duration argument cannot be zero ... setting to 1 millisecond.";
    // handleError( StkError::WARNING );
  }

  gRampPercent_ = rampPercent;
  if ( gRampPercent_ > 100 ) {
    gRampPercent_ = 100;
    // oStream_ << "Granulate::setGrainParameters: rampPercent argument cannot be greater than 100 ... setting to 100.";
    // handleError( StkError::WARNING );
  }

  gOffset_ = offset;
  gDelay_ = delay;
}

void Granulate :: setRandomFactor( double randomness )
{
  if ( randomness < 0.0 ) gRandomFactor_ = 0.0;
  else if ( randomness > 1.0 ) gRandomFactor_ = 0.97;

  gRandomFactor_ = 0.97 * randomness;
};

void Granulate :: openFile( int projectedBufferSize )
{

  data_ = &ring_;
  data_->resize( projectedBufferSize );
  ring_.resize( projectedBufferSize );
  std::memset( data_->elems(), 0.0, projectedBufferSize * sizeof(float));
  std::memset( ring_.elems(), 0.0, projectedBufferSize * sizeof(float));

  lastFrame_ = 0.0;

  this->reset();

#if defined(_STK_DEBUG_)
  std::ostringstream message;
  message << "Granulate::openFile: file = " << fileName << ", file frames = " << file.fileSize() << '.';
  handleError( message.str(), StkError::DEBUG_PRINT );
#endif

}

void Granulate :: getProjectedSignal( Array<double> projectedSignal )
{
  // std::memcpy(ring_.elems(), projectedSignal, file.frames() * numChannels_ * sizeof(float));
}

void Granulate :: reset( void )
{
  gPointer_ = 0;

  // Reset grain parameters.
  size_t count;
  size_t nVoices = (unsigned int)grains_.size();
  for ( unsigned int i=0; i<grains_.size(); i++ ) {
    grains_[i].repeats = 0;
    count = ( i * gDuration_ * 0.001 * sampleRate_ / nVoices );
    grains_[i].counter = count;
    grains_[i].state = GRAIN_STOPPED;
  }

    lastFrame_ = 0.0;

  // std::cout << nVoices << std::endl;
}

void Granulate :: setVoices( unsigned int nVoices )
{
  // std::cout << nVoices << std::endl;
#if defined(_STK_DEBUG_)
  std::ostringstream message;
  message << "Granulate::setVoices: nVoices = " << nVoices << ", existing voices = " << grains_.size() << '.';
  handleError( message.str(), StkError::DEBUG_PRINT );
#endif

  size_t oldSize = grains_.size();
  grains_.resize( nVoices );

  // Initialize new grain voices.
  size_t count;
  // for ( size_t i=0; i<nVoices; i++ ) { //weird singing effect
  for ( size_t i=oldSize; i<nVoices; i++ ) {
    grains_[i].repeats = 0;
    count = ( i * gDuration_ * 0.001 * sampleRate_ / nVoices );
    grains_[i].counter = count;
    grains_[i].pointer = gPointer_;
    grains_[i].state = GRAIN_STOPPED;
  }

  gain_ = 1.0 / grains_.size();
}

void Granulate :: calculateGrain( Granulate::Grain& grain )
{
  // std::cout << "ATTACK!!!" << std::endl;
  size_t grainSize = (  gDuration_ * 0.001 * sampleRate_  );
  grain.rms = 0;
  grain.rmsAccum = 0;
  for (int k = 0; k< grainSize; k++){
    grain.rmsAccum += data_->operator[]( grain.startPointer + k ) *
                      data_->operator[]( grain.startPointer + k ) ;
  }
  grain.rms = sqrt(grain.rmsAccum / (float) gDuration_) * gRmsGain_;
  // std::cout << " rmsAccum: "<< grain.rmsAccum << std::endl;

  if (grain.rms <= gRmsThreshold_){
    grain.gatedActive = false;
  }else{
    grain.gatedActive = true;
  }

  if ( grain.repeats > 0 ) {
    grain.repeats--;
    grain.pointer = grain.startPointer;
    if ( grain.attackCount > 0 ) {
      grain.eScaler = 0.0;
      grain.eRate = -grain.eRate;
      grain.counter = grain.attackCount;
      grain.state = GRAIN_FADEIN;
    }
    else {
      grain.counter = grain.sustainCount;
      grain.state = GRAIN_SUSTAIN;
    }
    return;
  }

  // Calculate duration and envelope parameters.
  double seconds = gDuration_ * 0.001;
  seconds += ( seconds * gRandomFactor_ * noise() );
  // seconds += ( seconds * gRandomFactor_ * ( 2.0 * rand() / (RAND_MAX + 1.0) - 1.0 ) );
  unsigned long count = (unsigned long) ( seconds * sampleRate_ );
  grain.attackCount = (unsigned int) ( gRampPercent_ * 0.005 * count );
  grain.decayCount = grain.attackCount;
  grain.sustainCount = count - 2 * grain.attackCount;
  grain.eScaler = 0.0;
  if ( grain.attackCount > 0 ) {
    grain.eRate = 1.0 / grain.attackCount;
    grain.counter = grain.attackCount;
    grain.state = GRAIN_FADEIN;
  }
  else {
    grain.counter = grain.sustainCount;
    grain.state = GRAIN_SUSTAIN;
  }

  // Calculate delay parameter.
  seconds = gDelay_ * 0.001;
  seconds += ( seconds * gRandomFactor_ * noise() );
  // seconds += ( seconds * gRandomFactor_ * ( 2.0 * rand() / (RAND_MAX + 1.0) - 1.0 ) );
  count = (unsigned long) ( seconds * sampleRate_ );
  grain.delayCount = count;

  // Save stretch parameter.
  grain.repeats = gStretch_;

  // Calculate offset parameter.
  seconds = gOffset_ * 0.001;
  seconds += ( seconds * gRandomFactor_ * std::abs( noise() ) );
  // seconds += ( seconds * gRandomFactor_ * std::abs( ( 2.0 * rand() / (RAND_MAX + 1.0) - 1.0 ) ) );
  int offset = (int) ( seconds * sampleRate_ );

  // Add some randomization to the pointer start position.
  seconds = gDuration_ * 0.001 * gRandomFactor_ * noise();
  // seconds = gDuration_ * 0.001 * gRandomFactor_ * ( 2.0 * rand() / (RAND_MAX + 1.0) - 1.0 );
  offset += (int) ( seconds * sampleRate_ );
  grain.pointer += offset;
  while ( grain.pointer >= data_->size() ) grain.pointer -= data_->size();
  if ( grain.pointer <  0 ) grain.pointer = 0;
  grain.startPointer = grain.pointer;
}


double Granulate :: tick()
{
  // std::cout << grains_.size() << std::endl;
  unsigned int i, j;
  lastFrame_ = 0.0;

  if ( data_->size() == 0 ) return 0.0;

  double sample;
  gNumOfGatedActiveGrains_ = 0;
  for ( i=0; i<grains_.size(); i++ ) {
    if ( grains_[i].counter == 0 ) { // Update the grain state.
      switch ( grains_[i].state ) {

      case GRAIN_STOPPED:
      // We're done waiting between grains ... setup for new grain
        this->calculateGrain( grains_[i] );
        break;

      case GRAIN_FADEIN:
        // We're done ramping up the envelope
        if ( grains_[i].sustainCount > 0 ) {
          grains_[i].counter = grains_[i].sustainCount;
          grains_[i].state = GRAIN_SUSTAIN;
          break;
        }
        // else no sustain state (i.e. perfect triangle window)

      case GRAIN_SUSTAIN:
        // We're done with flat part of envelope ... setup to ramp down
        if ( grains_[i].decayCount > 0 ) {
          grains_[i].counter = grains_[i].decayCount;
          grains_[i].eRate = -grains_[i].eRate;
          grains_[i].state = GRAIN_FADEOUT;
          break;
        }
        // else no fade out state (gRampPercent = 0)

      case GRAIN_FADEOUT:
        // We're done ramping down ... setup for wait between grains
        if ( grains_[i].delayCount > 0 ) {
          grains_[i].counter = grains_[i].delayCount;
          grains_[i].state = GRAIN_STOPPED;
          break;
        }
        // else no delay (gDelay = 0)

        this->calculateGrain( grains_[i] );
      }
  }

    if ( grains_[i].state > 0 ) {
        // sample = ring_[ nChannels * grains_[i].pointer + j ];
        sample = data_->operator[]( grains_[i].pointer );
        if ( grains_[i].state == GRAIN_FADEIN || grains_[i].state == GRAIN_FADEOUT ) {
          sample *= grains_[i].eScaler;
          grains_[i].eScaler += grains_[i].eRate;
        }

        if (grains_[i].gatedActive == true){
          // gNumOfGatedActiveGrains_++;
          gNumOfGatedActiveGrains_ = 1;
          // lastFrame_ += (sample * gain_);
          lastFrame_ += sample;
        }else{
          lastFrame_ += 0;
        }

      // Increment and check pointer limits.
      grains_[i].pointer++;
      if ( grains_[i].pointer >= data_->size() )
        grains_[i].pointer = 0;
    }
    grains_[i].counter--;
  }
  // std::cout << gNumOfGatedActiveGrains_ << std::endl;

  // Increment our global file pointer at the stretch rate.
  if ( stretchCounter_++ == gStretch_ ) {
    gPointer_++;
    if ( (unsigned long) gPointer_ >= data_->size() ) gPointer_ = 0;
    stretchCounter_ = 0;
  }

  return lastFrame_;
}
