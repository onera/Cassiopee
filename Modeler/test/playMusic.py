# - playMusic -
import Modeler.Sound as Sound
import time

music = Sound.playMusic("Images/song.wav")
while music.is_active():
    time.sleep(0.1)
Sound.endMusic(music)
