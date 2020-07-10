import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

x1, y1, z1 = 1, 1, 1

# angles
theta = random.random() * np.pi
phi = random.random() * 2 * np.pi

# cartesian coordinates
x2 = x1 + 1 * np.sin(theta) * np.cos(phi)
y2 = y1 + 1 * np.sin(theta) * np.sin(phi)
z2 = z1 + 1 * np.cos(theta)

# debug
print("X: {}, Y: {}, Z: {}".format(x1, y1, z1))
print("Theta: {}, Phi: {}".format(theta, phi))
print("Theta > Pi/2? {}".format(theta>np.pi/2))

# plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x1, y1, z1, color='r')
ax.scatter(x2, y2, z2, color='k')
ax.plot([x1, x2], [y1, y2], [z1, z2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(-2, 2)
plt.show()

print("***********************************")
x1, y1, z1 = 1, 3, 2
x2, y2, z2 = -5, 0, 4
print("X1: {}, Y1: {}, Z1: {}".format(x1, y1, z1))
print("X2: {}, Y2: {}, Z2: {}".format(x2, y2, z2))
dx = x2-x1
dy = y2-y1
dz = z2-z1
print("Dx: {}, Dy: {}, Dz: {}".format(dx, dy, dz))

linex, liney, linez = [], [], []
for t in range(10):
    x = x1 + t*dx
    y = y1 + t*dy
    z = z1 + t*dz
    linex.append(x)
    liney.append(y)
    linez.append(z)
int_point_x = linex[-2]
int_point_y = liney[-2]
print("Intersection X: {}, Y: {}".format(int_point_x, int_point_y))
int_point_z = linez[-2]
print((int_point_x - x1)/dx)
print((int_point_y - y1)/dy)
der_point_z = z1 + (int_point_x-x1)/dx * dz
print(int_point_z)
print(der_point_z)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x1, y1, z1, color='r')
ax.scatter(x2, y2, z2, color='k')
ax.scatter(int_point_x, int_point_y, int_point_z, color='b')
ax.plot(linex, liney, linez)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

