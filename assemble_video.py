import imageio
import numpy as np
import time
import os
import math
import glob

from PIL import Image

frames_between_keyframes = 60
maximum_keyframe_number = 728
zoom_scale = 2.0

# log1.1 of 2 is 7.27
# 60 / that is 8.25, so lets take 8 frames per keyframe

# Arguments for ffmpeg
kargs = {
    'macro_block_size': 8,
    'fps': 60,
    'format': 'FFMPEG',
    'quality': 8
}

t0 = time.time()

segment = 0
segment_names = []

framebuffer = []

writer = imageio.get_writer(f"segment_{segment:08}.mp4", **kargs)
segment_names.append(f"segment_{segment:08}.mp4")

files = glob.glob("output/*.png")
files.sort()

# We start with the previous image, this is modified
previous_keyframe = Image.open(files[0])

(width, height) = (previous_keyframe.width, previous_keyframe.height)
(scaled_width, scaled_height) = (int(zoom_scale * width), int(zoom_scale * height))

val1 = (scaled_width - width) // 2
val2 = (scaled_height - height) // 2

placement_box = (val1, val2, val1 + width, val2 + height)

for i in range(0, maximum_keyframe_number):
    # Creates 10 second segments
    if (i * frames_between_keyframes) % 600 == 0 and i != 0:
        # Add current_framebuffer frames
        for frame in reversed(framebuffer):
            writer.append_data(frame)
        
        segment += 1

        framebuffer = []

        # Writer for adding frames to the video
        writer = imageio.get_writer(f"segment_{segment:08}.mp4", **kargs)
        segment_names.append(f"segment_{segment:08}.mp4")

    next_keyframe = Image.open(files[i + 1])

    # This has to be nearest
    scaled_keyframe = next_keyframe.resize((scaled_width, scaled_height), resample=Image.NEAREST)
    
    scaled_keyframe.paste(previous_keyframe, placement_box)

    current_scale = zoom_scale
    factor = 10**(math.log10(zoom_scale) / frames_between_keyframes)

    for j in range(0, frames_between_keyframes):
        # zoom = 1.0 + (2**(1.0 - j / frames_between_keyframes) - 1.0) * (zoom_scale - 1.0)
        # zoom = zoom_scale**(1.0 - j / frames_between_keyframes)
        # zoom = zoom_scale**((1.0 - j / frames_between_keyframes))
        current_scale /= factor

        temp1 = (1.0 - 1.0 / current_scale) * (scaled_width // 2)
        temp2 = (1.0 - 1.0 / current_scale) * (scaled_height // 2)
        temp3 = scaled_width - temp1
        temp4 = scaled_height - temp2

        next_frame = scaled_keyframe.crop((temp1, temp2, temp3, temp4)).resize((width, height), resample=Image.LANCZOS)

        framebuffer.append(np.array(next_frame))

        elapsed = time.time() - t0
        print(f"{i * frames_between_keyframes + j + 1}/{maximum_keyframe_number * frames_between_keyframes + 1} {elapsed:.2f}s {((i * frames_between_keyframes + j + 1)/elapsed):.2f} fps")

    previous_keyframe = scaled_keyframe.resize((width, height), resample=Image.LANCZOS)


framebuffer.append(np.array(previous_keyframe))

for frame in reversed(framebuffer):
    writer.append_data(frame)

with open('segments.txt', 'w') as f:
    for segment in reversed(segment_names):
        f.write(f"file '{segment}'\n")

writer.close()

elapsed = time.time() - t0
print(f"{maximum_keyframe_number * frames_between_keyframes + 1}/{maximum_keyframe_number * frames_between_keyframes + 1} {elapsed:.2f}s {((maximum_keyframe_number * frames_between_keyframes + 1)/elapsed):.2f} fps")

os.system("ffmpeg -y -f concat -safe 0 -i segments.txt -c copy output.mp4")