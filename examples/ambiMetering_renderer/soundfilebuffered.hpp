#ifndef SOUNDFILEBUFFERED_H
#define SOUNDFILEBUFFERED_H

#include <mutex>
#include <thread>
#include <condition_variable>
#include <memory>

#include "Gamma/SoundFile.h"
#include "allocore/types/al_SingleRWRingBuffer.hpp"


namespace al
{

class SoundFileBuffered
{
public:
  SoundFileBuffered(std::string fullPath, bool loop = false, int bufferFrames = 1024);
  ~SoundFileBuffered();
  
  int read(float *buffer, int numFrames);
  
  bool opened() const;								///< Returns whether the sound file is open
  gam::SoundFile::EncodingType encoding() const;	///< Get encoding type
  gam::SoundFile::Format format() const;			///< Get format
  double frameRate() const;							///< Get frames/second
  int frames() const;								///< Get number of frames
  int channels() const;								///< Get number of channels
  int samples() const;								///< Get number of samples ( = frames x channels)
  
  typedef void (*CallbackFunc)(float *buffer, int numChannels,
                               int numFrames, void * userData);
  
  void setReadCallback(CallbackFunc func, void *userData);
  
private:
  bool mRunning;
  bool mLoop;
  std::mutex mLock;
  std::condition_variable mCondVar;
  std::thread *mReaderThread;
  SingleRWRingBuffer *mRingBuffer;
  int mBufferFrames;
  
  gam::SoundFile mSf;
  CallbackFunc mReadCallback;
  void *mCallbackData;
  
  static void readFunction(SoundFileBuffered *obj);
};

} // namespace al

#endif // SOUNDFILEBUFFERED_H
