ffmpeg -framerate 25 -i %05d.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 10 -pix_fmt yuv420p 3D_Heaving-Pitching.mp4
