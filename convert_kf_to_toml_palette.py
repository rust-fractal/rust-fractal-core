

file_read = open('test.kfp', 'r')

colors = []
palette_iteration_span = 1.0

for line in file_read.readlines():
    if line.startswith('Colors:'):
        colors = line[7:].replace(" ", "").split(",")[:-1]
    if line.startswith('IterDiv:'):
        palette_iteration_span = float(line[8:].replace(" ", ""))

temp = colors.copy()

print(colors)

colors[0::3] = temp[2::3]
colors[2::3] = temp[0::3]

file_new = open('test.toml', 'w')
file_new.write(f"palette = {[int(i) for i in colors]}\n\palette_iteration_span = {1 / palette_iteration_span}")