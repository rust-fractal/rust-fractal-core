import glob

files = glob.glob("./convert/fractalalex/*.kfp")

for file in files:
    file_read = open(file, 'r')

    colors = []
    palette_iteration_span = 1.0

    for line in file_read.readlines():
        if line.startswith('Colors:'):
            colors = line[7:].replace(" ", "").split(",")[:-1]
        if line.startswith('IterDiv:'):
            palette_iteration_span = 0.1 * 1024 / float(line[8:].replace(" ", ""))

    temp = colors.copy()

    colors[0::3] = temp[2::3]
    colors[2::3] = temp[0::3]

    file_name = file.split('\\')[1].split(".")[0]

    file_new = open(f'./palettes/fractalalex/{file_name}.toml', 'w')
    file_new.write(f"palette = {[int(i) for i in colors]}\n\npalette_iteration_span = {palette_iteration_span}")