import numpy as np


def rotate_angle(va: np.array, vb: np.array):
    la = np.sum(va ** 2) ** 0.5
    lb = np.sum(vb ** 2) ** 0.5
    return np.arccos(np.dot(va, vb) / la * lb)  # rad


def rotate_axis(va: np.array, vb: np.array):
    return np.cross(va, vb)


def rotate_matrix(va, vb):
    vc = rotate_axis(va, vb)
    theta = rotate_angle(va, vb)
    lc = np.sum(vc ** 2) ** 0.5
    vc = vc / lc
    w = np.array([
        [0, -vc[2], vc[1]],
        [vc[2], 0, -vc[0]],
        [-vc[1], vc[0], 0],
    ])
    return np.eye(3) + w * np.sin(theta) + np.dot(w, w) * (1 - np.cos(theta))


if __name__ == '__main__':
    vector_a = np.array([1, 0, 0])
    vector_b = np.array([0, 1, 0])
    # length_b = 1
    # theta = math.pi / 2
    angle = rotate_angle(vector_a, vector_b)
    axis = rotate_axis(vector_a, vector_b)
    matrix = rotate_matrix(vector_a, vector_b)
    print(angle)
    print(axis)
    print(matrix)
    print(np.dot(matrix, vector_a))
