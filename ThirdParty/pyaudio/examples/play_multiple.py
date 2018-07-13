"""PyAudio Example: Play a wave file (callback version)."""

import pyaudio
import wave
import time
import sys

slot = 0

fileName = '380_gunshot_single-mike-koenig.wav'
#fileName = 'shotgun-mossberg590-RA_The_Sun_God-451502290.wav'
wf = wave.open(fileName, 'rb')
wf2 = wave.open(fileName, 'rb')
wf3 = wave.open(fileName, 'rb')

# instantiate PyAudio (1)
p = pyaudio.PyAudio()

# define callback (2)
def callback(in_data, frame_count, time_info, status):
    data = wf.readframes(frame_count)
    if data == '': wf.rewind(); data = wf.readframes(frame_count)
    return (data, pyaudio.paContinue)

def callback2(in_data, frame_count, time_info, status):
    data = wf2.readframes(frame_count)
    if data == '': wf2.rewind(); data = wf2.readframes(frame_count)
    return (data, pyaudio.paContinue)

def callback3(in_data, frame_count, time_info, status):
    data = wf3.readframes(frame_count)
    if data == '': wf3.rewind(); data = wf3.readframes(frame_count)
    return (data, pyaudio.paContinue)

def gunshot():
    global slot
    if slot == 0:
        stream = p.open(format=p.get_format_from_width(wf.getsampwidth()),
                        channels=wf.getnchannels(),
                        rate=wf.getframerate(),
                        output=True,
                        stream_callback=callback)
        slot = 1
    elif slot == 1:
        stream = p.open(format=p.get_format_from_width(wf2.getsampwidth()),
                        channels=wf2.getnchannels(),
                        rate=wf2.getframerate(),
                        output=True,
                        stream_callback=callback2)
        slot = 2
    elif slot == 2:
        stream = p.open(format=p.get_format_from_width(wf3.getsampwidth()),
                        channels=wf3.getnchannels(),
                        rate=wf3.getframerate(),
                        output=True,
                        stream_callback=callback3)
        slot = 0
    stream.start_stream()
    return stream

# start the stream (4)
#stream.start_stream()

# wait for stream to finish (5)
#while stream.is_active():
#time.sleep(0.1)

# stop stream (6)
#stream.stop_stream()

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
