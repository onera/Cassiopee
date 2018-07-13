# Sound support

try: import pyaudio
except: raise ImportError, "Sound module requires pyaudio."

import wave

# Globals
audioHandle = None
musicFileHandle = None
poolSize = 3
soundPool = [None]*poolSize

# Init sound module
def initSound():
    global audioHandle
    audioHandle = pyaudio.PyAudio()
    return audioHandle

def closeSound():
    audioHandle.terminate()

# def music callback (internal)
def musicCallback__(in_data, frame_count, time_info, status):
    data = musicFileHandle.readframes(frame_count)
    if data == '': musicFileHandle.rewind(); musicFileHandle.readframes(frame_count)
    return (data, pyaudio.paContinue)

# Play music from a wav file
# Retourne le musicHandle
def playMusic(fileName):
    global musicFileHandle
    musicFileHandle = wave.open(fileName, 'rb')
    if audioHandle is None: initSound()
    stream = audioHandle.open(format=audioHandle.get_format_from_width(musicFileHandle.getsampwidth()),
                    channels=musicFileHandle.getnchannels(),
                    rate=musicFileHandle.getframerate(),
                    output=True,
                    stream_callback=musicCallback__)
    stream.start_stream()
    return stream

# Stop music
def stopMusic(musicHandle):
    musicHandle.stop_stream()

# End music and close everything
def endMusic(musicHandle):
    global musicFileHandle
    musicHandle.stop_stream()
    musicHandle.close()
    musicFileHandle.close()
    musicFileHandle = None

# Register sound (store in memory)
# Retourne soundHandle
def registerSound(fileName):
    wf = wave.open(fileName, 'rb')
    nchannels = wf.getnchannels()
    sampleWidth = wf.getsampwidth()
    frameRate = wf.getframerate()
    nframes = wf.getnframes()
    gdata = wf.readframes(nframes)
    sizeFrame = len(gdata)/nframes
    wf.close()
    # stream, pt, nchannels, sampleWidth, frameRate, sizeFrame, gdata
    return [None, 0, nchannels, sampleWidth, frameRate, sizeFrame, gdata]

def soundCallback0__(in_data, frame_count, time_info, status):
    h = soundPool[0]
    pt = h[1]; s = h[5]; gdata = h[6]
    data = gdata[s*pt:s*pt+s*frame_count]
    pt += frame_count; h[1] = pt
    if s*pt > len(gdata): h[1] = 0
    return (data, pyaudio.paContinue)
def soundCallback1__(in_data, frame_count, time_info, status):
    h = soundPool[1]
    pt = h[1]; s = h[5]; gdata = h[6]
    data = gdata[s*pt:s*pt+s*frame_count]
    pt += frame_count; h[1] = pt
    if s*pt > len(gdata): h[1] = 0
    return (data, pyaudio.paContinue)
def soundCallback2__(in_data, frame_count, time_info, status):
    h = soundPool[2]
    pt = h[1]; s = h[5]; gdata = h[6]
    data = gdata[s*pt:s*pt+s*frame_count]
    pt += frame_count; h[1] = pt
    if s*pt > len(gdata): h[1] = 0
    return (data, pyaudio.paContinue)

def playSound(soundHandle):
    if audioHandle is None: initSound()
    # Cherche un pool de libre
    i = -1
    for j in xrange(poolSize): 
        if soundPool[j] is None: i = j; break
    if i == -1: return None # full
    #print 'found free pool=%d'%i
    if i == 0: callback = soundCallback0__
    elif i == 1: callback = soundCallback1__
    elif i == 2: callback = soundCallback2__
    else: callback = soundCallback0__
    soundPool[i] = [None, 0, soundHandle[2], soundHandle[3], soundHandle[4], 
    soundHandle[5], soundHandle[6]]
    stream = audioHandle.open(format=audioHandle.get_format_from_width(soundHandle[3]),
                              channels=soundHandle[2],
                              rate=soundHandle[4],
                              output=True,
                              stream_callback=callback)
    soundPool[i][0] = stream
    stream.start_stream()
    return stream

# Close all sounds (if finished)
def closeAllSounds():
    for i in xrange(poolSize):
        h = soundPool[i]
        #print 'Checking pool=',i
        if h is not None:
            s = h[0]
            if not s.is_active(): 
                s.close()
                soundPool[i] = None
                #print "closing pool=%d"%i
