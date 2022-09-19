import math

import numpy as np


def vector_angle(length_b, vector_a: np.array, theta):
    length_a = np.sum(vector_a ** 2) ** 0.5
    return length_a * length_b * math.cos(theta) / vector_a


if __name__ == '__main__':
    a = np.array([0, 2, 0])
    length_b = 1
    theta = math.pi / 2
    print(vector_angle(length_b, a, theta))
