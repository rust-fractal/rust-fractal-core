import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_theme(style="darkgrid")

f, ax = plt.subplots(figsize=(7, 7))
ax.set(yscale="log")

f = open("test.csv", "r")

for line in f:
    values = np.fromstring(line, sep = ",")

    probe_values = values[0:-1]
    probe_mask = np.logical_or(np.logical_or(probe_values == 0.0, np.isnan(probe_values)), probe_values > 1e200)
    plt.plot(np.flatnonzero(~probe_mask), probe_values[~probe_mask], alpha=0.2)


plt.show()