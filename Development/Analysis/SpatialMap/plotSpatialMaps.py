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
ge_angle_centers = np.linspace(0, np.pi * 2.0, n_ge_crystals+1)

plt.plot(line_x, line_y, color='blue')
plt.plot(inner_x, inner_y, color='blue')
plt.plot(outer_x, outer_y, color='blue')
plt.plot(ge_x, ge_y, color='orange')

ge_cry_x, ge_cry_y = make_circle(40)
rotation_angle = 0.0
for ge_ang in ge_angle_centers[:-1]:
    x0, y0 = get_ge_center(ge_ang)
    x = np.cos(rotation_angle) * x0 - np.sin(rotation_angle) * y0
    y = np.sin(rotation_angle) * x0 + np.cos(rotation_angle) * y0
    print("{} {}".format(x, y))
    plt.scatter(x, y, color='red')
    plt.plot(x + ge_cry_x, y + ge_cry_y, color='red')

x0, y0 = 700, 0
for a in np.linspace(0, 2*np.pi, 1):
    x1, y1 = x0 + np.cos(a), y0 + np.sin(a)
    m = (y1 - y0) / (x1 - x0)
    x = np.linspace(-300, +700, 300)
    y = m * (x - x0) + y0
    plt.plot(x, y, color='green')

minGeRot = -.22
maxGeRot = +.22
minrot_x, minrot_y = 300*np.cos(minGeRot), 300*np.sin(minGeRot)
maxrot_x, maxrot_y = 300*np.cos(maxGeRot), 300*np.sin(maxGeRot)
plt.plot([0, minrot_x], [0, minrot_y])
plt.plot([0, maxrot_x], [0, maxrot_y])
plt.xlim(-300, +750)
plt.ylim(-300, +750)
plt.show()

plt.cla()
ax, ay = 100, -10
bx, by = 0, 845
plt.hlines(y=845, xmin=-300, xmax=+300)
plt.vlines(x=100, ymin=-10, ymax=+845)
plt.scatter(ax, ay)
plt.scatter(bx, by)
plt.plot([0, 100], [845, -10])
alpha = np.arctan2(100, 855)
print(np.rad2deg(alpha))
len_ipo = (845 - (-10)) / np.cos(alpha)
print("Len Ipotenusa: {}".format(len_ipo))
print("Double Check: {}".format(((bx - ax)**2 + (by - ay)**2)**.5))
plt.xlim(-845, 845)
plt.ylim(-845, 845)
plt.show()

exit(0)

def get_angle_three_points(center, a, b):
    angle_a = np.arctan2(a[1] - center[1], a[0] - center[0])
    if angle_a < 0:
        angle_a += 2 * np.pi
    angle_b = np.arctan2(b[1] - center[1], b[0] - center[0])
    if angle_b < 0:
        angle_b += 2 * np.pi
    return angle_b - angle_a

def compute_pr(radius, z):
    # params
    top_z_axis = +845
    bottom_z_axis = -845
    top_z_ge = +425
    bottom_z_ge = -425
    left_ge_margin = -(235+40)
    right_ge_margin = 235+40
    # cases
    if(z > top_z_ge):
        a = (right_ge_margin, top_z_axis)
        b = (left_ge_margin, top_z_ge)
        c = (right_ge_margin, bottom_z_ge)
        d = (right_ge_margin, bottom_z_axis)
    elif(z > bottom_z_ge):
        a = (right_ge_margin, top_z_axis)
        b = (right_ge_margin, top_z_ge)
        c = (right_ge_margin, bottom_z_ge)
        d = (right_ge_margin, bottom_z_axis)
    else:
        a = (right_ge_margin, top_z_axis)
        b = (right_ge_margin, top_z_ge)
        c = (left_ge_margin, bottom_z_ge)
        d = (right_ge_margin, bottom_z_axis)
    # compute angles
    # idea: AD angle is the difference between angle D w.r.t. point (radius,z) and angle A w.r.t point
    ad_angle = get_angle_three_points((radius, z), a, d) 
    bc_angle = get_angle_three_points((radius, z), b, c)
    # debug
    plt.cla()
    plt.scatter(radius, z, label="Point")
    plt.scatter(a[0], a[1], label="A", color='b')
    plt.scatter(b[0], b[1], label="B", color='r')
    plt.scatter(c[0], c[1], label="C", color='r')
    plt.scatter(d[0], d[1], label="D", color='b')
    plt.plot([radius, a[0]], [z, a[1]], color='g')
    plt.plot([radius, b[0]], [z, b[1]], color='g')
    plt.plot([radius, c[0]], [z, c[1]], color='g')
    plt.plot([radius, d[0]], [z, d[1]], color='g')
    plt.hlines(y=z, xmin=left_ge_margin, xmax=radius)
    plt.vlines(x=radius, ymin=bottom_z_axis, ymax=top_z_axis)
    plt.legend()
    #print("AD Angle: {}, BC Angle: {}".format(np.rad2deg(ad_angle), np.rad2deg(bc_angle)))
    plt.pause(.0001)
    return bc_angle / ad_angle

radius = 500
probs = []
for z in range(-845, 850, 5):
    pr = compute_pr(radius, z)
    print("Radius: {}, Z: {} => {}".format(radius, z, pr))
    probs.append(pr)

plt.figure()
plt.plot(range(len(probs)), probs, label="Probs. from -845 to 845")
plt.legend()
plt.show()

print("Done.")
