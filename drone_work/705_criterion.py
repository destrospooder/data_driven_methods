import numpy as np

def within_usable_area(vertex1, vertex2, minor_diameter, point):
    center = np.array([np.mean([vertex1[0], vertex2[0]]), 
                       np.mean([vertex1[1], vertex2[1]])])

    a = np.sqrt((vertex2[0] - vertex1[0])**2 + (vertex2[1] - vertex1[1])**2) / 2
    b = minor_diameter / 2

    theta = np.arctan2(vertex2[1] - vertex1[1], vertex2[0] - vertex1[0])

    # switching coord sys
    x_translated = point[0] - center[0]
    y_translated = point[1] - center[1]
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    x_rotated = x_translated * cos_theta + y_translated * sin_theta
    y_rotated = -x_translated * sin_theta + y_translated * cos_theta

    # check if the point is inside the usable area
    lhs = (x_rotated**2 / a**2) + (y_rotated**2 / b**2)
    return lhs <= 1

# 705 usable area conservative estimate
vertex1 = (-1.15, 2)
vertex2 = (1.7, -1.49)
minor_diameter = 2.25

# replace with whatever point you wanna test
point = (0.5, 0.5)

print(within_usable_area(vertex1, vertex2, minor_diameter, point))