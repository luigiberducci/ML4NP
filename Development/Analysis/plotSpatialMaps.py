import numpy as np
import matplotlib.pyplot as plt

def make_circle(r):
    t = np.arange(0, np.pi * 2.0, 0.01)
    t = t.reshape((len(t), 1))
    x = r * np.cos(t)
    y = r * np.sin(t)
    return x, y

def get_ge_center(angle, r=235):
    return 235 * np.cos(angle), 235 * np.sin(angle)

line_x = np.linspace(0, 710, 711)
line_y = np.zeros(711)

inner_x, inner_y = make_circle(175)
outer_x, outer_y = make_circle(295)
ge_x, ge_y = make_circle(235)

# compute center of ge crystals
n_ge_crystals = 14
ge_angle_centers = np.linspace(0, np.pi * 2.0, n_ge_crystals + 1)

plt.plot(line_x, line_y, color='blue')
plt.plot(inner_x, inner_y, color='blue')
plt.plot(outer_x, outer_y, color='blue')
plt.plot(ge_x, ge_y, color='orange')

ge_cry_x, ge_cry_y = make_circle(40)
rotation_angle = 0.22
for ge_ang in ge_angle_centers:
    x0, y0 = get_ge_center(ge_ang)
    x = np.cos(rotation_angle) * x0 - np.sin(rotation_angle) * y0
    y = np.sin(rotation_angle) * x0 + np.cos(rotation_angle) * y0
    plt.scatter(x, y, color='red')
    plt.plot(x + ge_cry_x, y + ge_cry_y, color='red')

x0, y0 = 350, 0
for a in np.linspace(0, 2*np.pi, 100):
    x1, y1 = x0 + np.cos(a), y0 + np.sin(a)
    m = (y1 - y0) / (x1 - x0)
    x = np.linspace(-300, +700, 300)
    y = m * (x - x0) + y0
    plt.plot(x, y, color='green')

plt.xlim(-300, +750)
plt.ylim(-300, +750)
plt.show()
