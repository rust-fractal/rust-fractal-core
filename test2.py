import imageio
import numpy as np
import time
import os

from PIL import Image

frames_between_keyframes = 60
maximum_keyframe_number = 783
zoom_scale = 2

# Arguments for ffmpeg
kargs = {
    'macro_block_size': 8,
    'fps': 60,
    'format': 'FFMPEG',
    'quality': 5
}

t0 = time.time()

segment = 0
segment_names = []

framebuffer = []

writer = imageio.get_writer(f"segment_{segment:08}.mp4", **kargs)
segment_names.append(f"segment_{segment:08}.mp4")

# We start with the previous image, this is modified
previous_keyframe = Image.open(f"output/output_00000000.png")

(width, height) = (previous_keyframe.width, previous_keyframe.height)
(scaled_width, scaled_height) = (zoom_scale * width, zoom_scale * height)
placement_box = (scaled_width//4, scaled_height//4, scaled_width//4 * 3, scaled_height//4 * 3)

for i in range(0, maximum_keyframe_number):
    # Creates 10 second segments
    if (i * frames_between_keyframes) % 600 == 0 and i != 0:
        for frame in reversed(framebuffer):
            writer.append_data(frame)
        segment += 1
        framebuffer = []
        # Writer for adding frames to the video
        writer = imageio.get_writer(f"segment_{segment:08}.mp4", **kargs)
        segment_names.append(f"segment_{segment:08}.mp4")

    next_keyframe = Image.open(f"output/output_{i + 1:08}.png")

    scaled_keyframe = next_keyframe.resize((scaled_width, scaled_height), resample=Image.LANCZOS)
    
    scaled_keyframe.paste(previous_keyframe, placement_box)

    for j in range(0, frames_between_keyframes):
        zoom = zoom_scale**(1.0 - j / frames_between_keyframes)

        temp1 = (1.0 - 1.0 / zoom) * width
        temp2 = (1.0 - 1.0 / zoom) * height
        temp3 = scaled_width - temp1
        temp4 = scaled_height - temp2

        next_frame = scaled_keyframe.crop((temp1, temp2, temp3, temp4)).resize((width, height), resample=Image.HAMMING)

        framebuffer.append(np.array(next_frame))
        print(f"{i * frames_between_keyframes + j + 1}/{maximum_keyframe_number * frames_between_keyframes + 1} {time.time() - t0:.2f}s")

    previous_keyframe = scaled_keyframe.resize((width, height), resample=Image.LANCZOS)


framebuffer.append(np.array(previous_keyframe))

print(f"{maximum_keyframe_number * frames_between_keyframes + 1}/{maximum_keyframe_number * frames_between_keyframes + 1} {time.time() - t0:.2f}s")

for frame in reversed(framebuffer):
    writer.append_data(frame)

# Command to split up into 
# ffmpeg -i invid.mp4 -threads 3 -vcodec copy -f segment -segment_time 2 cam_out_h264%04d.mp4

with open('segments.txt', 'w') as f:
    for segment in reversed(segment_names):
        f.write(f"file '{segment}'\n")

writer.close()

os.system("ffmpeg -y -f concat -safe 0 -i segments.txt -c copy output.mp4")