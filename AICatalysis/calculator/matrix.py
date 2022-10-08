import math
from functools import singledispatch

import numpy as np


def r_angle(va: np.array, vb: np.array):
    la = np.sum(va ** 2) ** 0.5
    lb = np.sum(vb ** 2) ** 0.5
    return np.arccos(np.around(np.dot(va, vb) / (la * lb), 4))  # rad


def r_axis(va: np.array, vb: np.array):
    return np.cross(va, vb)


@singledispatch
def r_matrix(*args):
    pass


@r_matrix.register(np.ndarray)
def _(va: np.ndarray, vb: np.ndarray):
    """

    Args:
        va: vector before rotation
        vb: vector after rotation

    Returns: rotation matrix

    """
    v_axis = r_axis(vb, va)
    angle = r_angle(vb, va)
    l_axis = np.sum(v_axis ** 2) ** 0.5
    v_axis = v_axis / l_axis
    w = np.array([
        [0, -v_axis[2], v_axis[1]],
        [v_axis[2], 0, -v_axis[0]],
        [-v_axis[1], v_axis[0], 0],
    ])
    return np.eye(3) + w * np.sin(angle) + np.dot(w, w) * (1 - np.cos(angle))


@r_matrix.register(float)
def _(angle: float, v_axis: np.ndarray):
    l_axis = np.sum(v_axis ** 2) ** 0.5
    v_axis = v_axis / l_axis
    w = np.array([
        [0, -v_axis[2], v_axis[1]],
        [v_axis[2], 0, -v_axis[0]],
        [-v_axis[1], v_axis[0], 0],
    ])
    return np.eye(3) + w * np.sin(angle) + np.dot(w, w) * (1 - np.cos(angle))


if __name__ == '__main__':
    vector_a = np.array([1, 0, 0])
    vector_b = np.array([0, 1, 0])
    # length_b = 1
    theta = math.pi / 2
    # angle = r_angle(vector_a, vector_b)
    # axis = r_axis(vector_a, vector_b)
    matrix = r_matrix(vector_a, theta)
    # print(angle)
    # print(axis)
    print(matrix)
    # print(np.dot(vector_a, matrix))
