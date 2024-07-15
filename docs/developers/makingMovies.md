# Making movies

## Creating mp4 from images

      ffmpeg -framerate 20 -pattern_type glob -i "*.png" -b:v 50000k video.mp4

or

      ffmpeg -framerate 20 -i image%02d.png -b:v 50000k video.mp4

## Adding metadata for 360 mp4

      exiftool -XMP-GSpherical:Spherical="true" video.mp4