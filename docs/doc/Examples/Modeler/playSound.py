# - playSound -
try: import Modeler.Sound as Sound
except: import Modeler.SoundLess as Sound
import time

sound = Sound.registerSound("Images/380.wav")

c = 0
while c < 5:
    print("Playing sound %d"%c)
    Sound.playSound(sound)
    time.sleep(0.5)
    Sound.closeAllSounds()
    c += 1
