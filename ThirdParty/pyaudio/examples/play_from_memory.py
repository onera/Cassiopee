"""PyAudio Example: Play a wave file (callback version)."""

import pyaudio
import wave
import time
import sys

slot = 0
pt = 0

# Small file, read all!
fileName = '380_gunshot_single-mike-koenig.wav'
wf = wave.open(fileName, 'rb')
nchannels = wf.getnchannels()
sampleWidth = wf.getsampwidth()
frameRate = wf.getframerate()
nframes = wf.getnframes()
gdata = wf.readframes(nframes)
wf.close()

wf = wave.open(fileName, 'rb')

# instantiate PyAudio (1)
p = pyaudio.PyAudio()

# define callback (2)
def callback(in_data, frame_count, time_info, status):
    global pt
    data = gdata[4*pt:4*pt+4*frame_count]
    pt += frame_count
    if 4*pt > len(gdata): pt = 0
    return (data, pyaudio.paContinue)


def callback2(in_data, frame_count, time_info, status):
    data = wf.readframes(frame_count)
    return (data, pyaudio.paContinue)

def gunshot():
    stream = p.open(format=p.get_format_from_width(sampleWidth),
                    channels=nchannels,
                    rate=frameRate,
                    output=True,
                    stream_callback=callback)
    stream.start_stream()
    return stream

def gunshot2():
    stream = p.open(format=p.get_format_from_width(sampleWidth),
                    channels=nchannels,
                    rate=frameRate,
                    output=True,
                    stream_callback=callback2)
    stream.start_stream()
    return stream

import sys
streams = []
while 1 != 2:
    a = sys.stdin.read(1)
    #print '>>',a, len(streams)
    l = len(streams)
    if l < 4: s = gunshot(); streams.append(s); l += 1
    for i in xrange(l):
        j = l-i-1
        if streams[j].is_active() == False: 
            streams[j].close()
            del streams[j]
    time.sleep(0.03)

wf.close()
# close PyAudio (7)
p.terminate()
