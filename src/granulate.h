// #ifndef STK_GRANULATE_H
// #define STK_GRANULATE_H

#include "allocore/al_Allocore.hpp"
#include "alloaudio/al_SoundfileBuffered.hpp"
#include "allocore/ui/al_ParameterMIDI.hpp"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/Containers.h"
#include "Gamma/Noise.h"

#include <vector>

using namespace gam;

class Granulate
{
 public:
  //! Default constructor.
  Granulate( void );

  //! Constructor taking input audio file and number of voices arguments.
  Granulate( unsigned int nVoices, int projectedBufferSize );

  //! Class destructor.
  ~Granulate( void );

  void setSampleRate( const float sr );

  //! Write current frame to ring buffer
  void writeData( double ringBufferValue );

  //! Write current frame to ring buffer
  int getNumOfGatedGrains();

  //! Set the treshold for triggering grains
  void setThreshold (double rmsGain, double rmsThreshold);

  //! Load a monophonic soundfile to be "granulated".
  /*!
    An StkError will be thrown if the file is not found, its format
    is unknown or unsupported, or the file has more than one channel.
  */
  void openFile( int projectedBufferSize );

  // Setup for granulation
  // void setup( int projectedBufferSize );

  //! Get grain pointer position
  int getGrainPointer( );

  //! Project Signal for each granular stream
  void getProjectedSignal( Array<double> projectedSignal );

  //! Reset all STG's global pointers to 0
  void resetGlobalPointer( void );

  //! Reset the file pointer and all existing grains to the file start.
  /*!
    Multiple grains are offset from one another in time by grain
    duration / nVoices.
  */
  void reset( void );

  //! Set the number of simultaneous grain "voices" to use.
  /*!
    Multiple grains are offset from one another in time by grain
    duration / nVoices.  For this reason, it is best to set the grain
    parameters before calling this function (during initialization).
  */
  void setVoices( unsigned int nVoices = 1 );

  //! Set the stretch factor used for grain playback (1 - 1000).
  /*!
    Granular synthesis allows for time-stetching without affecting
    the original pitch of a sound.  A stretch factor of 4 will produce
    a resulting sound of length 4 times the orignal sound.  The
    default parameter of 1 produces no stretching.
  */
  void setStretch( unsigned int stretchFactor = 1 );

  //! Set global grain parameters used to determine individual grain settings.
  /*!
    Each grain is defined as having a length of \e duration
    milliseconds which must be greater than zero.  For values of \e
    rampPercent (0 - 100) greater than zero, a linear envelope will be
    applied to each grain.  If \e rampPercent = 100, the resultant
    grain "window" is triangular while \e rampPercent = 50 produces a
    trapezoidal window.  In addition, each grain can have a time delay
    of length \e delay and a grain pointer increment of length \e
    offset, which can be negative, before the next ramp onset (in
    milliseconds).  The \e offset parameter controls grain pointer
    jumps between enveloped grain segments, while the \e delay
    parameter causes grain calculations to pause between grains.  The
    actual values calculated for each grain will be randomized by a
    factor set using the setRandomFactor() function.
  */
  void setGrainParameters( unsigned int duration = 30, unsigned int rampPercent = 50,
                           int offset = 0, unsigned int delay = 0 );

  //! This factor is used when setting individual grain parameters (0.0 - 1.0).
  /*!
    This random factor is applied when all grain state durations
    are calculated.  If set to 0.0, no randomness occurs.  When
    randomness = 1.0, a grain segment of length \e duration will be
    randomly augmented by up to +- \e duration seconds (i.e., a 30
    millisecond length will be augmented by an extra length of up to
    +30 or -30 milliseconds).
   */
  // void setRandomFactor( double randomness = 0.1 );
  void setRandomFactor( double randomness_dur = 0.1, double randomness_delay = 0.1, double randomness_offset = 0.1, double randomness_pointer = 0.1);

  //! Return the specified channel value of the last computed frame.
  /*!
    The \c channel argument must be less than the number of output
    channels, which can be determined with the channelsOut() function
    (the first channel is specified by 0).  However, range checking is
    only performed if _STK_DEBUG_ is defined during compilation, in
    which case an out-of-range value will trigger an StkError
    exception. \sa lastFrame()
  */
  double lastOut();

  //! Compute one sample frame and return the specified \c channel value.
  double tick();

  //! Fill the StkFrames object with computed sample frames, starting at the specified channel.
  /*!
    The \c channel argument plus the number of output channels must
    be less than the number of channels in the StkFrames argument (the
    first channel is specified by 0).  However, range checking is only
    performed if _STK_DEBUG_ is defined during compilation, in which
    case an out-of-range value will trigger an StkError exception.
  */
  Array<double>& tick( Array<double>& frames, unsigned int channel = 0 );
  //! StkFrames: The default constructor initializes the frame data structure to size zero.

  enum GrainState {
    GRAIN_STOPPED,
    GRAIN_FADEIN,
    GRAIN_SUSTAIN,
    GRAIN_FADEOUT
  };

 protected:

  struct Grain {
    double eScaler;
    double eRate;
    unsigned long attackCount;
    unsigned long sustainCount;
    unsigned long decayCount;
    unsigned long delayCount;
    unsigned long counter;
    //unsigned long pointer;
    double pointer;
    unsigned long startPointer;
    unsigned int repeats;
    float rms;
    float rmsAccum;
    bool gatedActive;
    GrainState state;

    // Default constructor.
    Grain()
      :eScaler(0.0), eRate(0.0), attackCount(0), sustainCount(0), decayCount(0),
       delayCount(0), counter(0), pointer(0), startPointer(0), repeats(0), rms(0),
       rmsAccum(0), gatedActive(false), state(GRAIN_STOPPED) {}
  };

  void calculateGrain( Granulate::Grain& grain );

  // Array<double> * data_;
  Array<double> * data_;
  Ring<double> ring_;

  std::vector<Grain> grains_;
  // Noise noise;
  NoiseWhite<> noise;//( 2.0 * rand() / (RAND_MAX + 1.0) - 1.0 )
  // float noise;
  //long gPointer_;
  double gPointer_;

  // Global grain parameters.
  unsigned int gNumOfGatedActiveGrains_;
  unsigned int gDuration_;
  unsigned int gRampPercent_;
  unsigned int gDelay_;
  unsigned int gStretch_;
  unsigned int stretchCounter_;
  float gRmsThreshold_;
  float gRmsGain_;
  int gOffset_;
  // double gRandomFactor_;
  double gRandomFactor_dur_;
  double gRandomFactor_delay_;
  double gRandomFactor_offset_;
  double gRandomFactor_pointer_;
  double gain_;

  float sampleRate_;
  // unsigned int numChannels_;
  //!Generator- STK lastFrame_: Return an StkFrames reference to the last output sample frame.
  double lastFrame_;

};

inline double Granulate :: lastOut( )
{
  return lastFrame_;
}

// inline Array<double>& Granulate :: tick( Array<double>& frames, unsigned int channel )
// {
//   unsigned int nChannels = lastFrame_.channels();
// #if defined(_STK_DEBUG_)
//   if ( channel > frames.channels() - nChannels ) {
//     oStream_ << "Granulate::tick(): channel and StkFrames arguments are incompatible!";
//     handleError( StkError::FUNCTION_ARGUMENT );
//   }
// #endif
//
//   double *samples = &frames[channel];
//   unsigned int j, hop = frames.channels() - nChannels;
//   for ( unsigned int i=0; i<frames.frames(); i++, samples += hop ) {
//     *samples++ = tick();
//     for ( j=1; j<nChannels; j++ )
//       *samples++ = lastFrame_[j];
//   }
//
//   return frames;
// }


// #endif
