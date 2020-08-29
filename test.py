import os

os.system("ffmpeg -y -f concat -safe 0 -i segments.txt -c copy output.mp4")