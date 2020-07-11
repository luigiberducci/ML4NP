import numpy as np
import matplotlib.pyplot as plt

filename = "histo_neutron_danila.txt"

with open(filename, 'r') as f:
    lines = f.readlines()

x, y = [], []
for line in lines:
    line_split = line.split()
    x.append(float(line_split[0]))
    y.append(float(line_split[1]))

plt.step(x, y)
plt.xticks(range(0, 20, 2))
plt.ylim(0, 0.35)
plt.xlim(0, 20)
plt.show()
