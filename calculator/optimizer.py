import copy
import math
import random

import numpy as np

from common.logger import logger


class Optimizer(object):
    def __init__(self, center, ligand, steps=500, criterion=0.1, lr=0.1):
        self.center = center
        self.ligand = ligand
        self.steps = steps
        self.criterion = criterion
        self.lr = lr

    @staticmethod
    def _get_force(center, ligand):
        distance = np.max((ligand - center) ** 2, axis=0)
        distance = np.where(distance == 0, 0.1, distance)

        direction = np.max((ligand - center), axis=0)
        if np.all(direction) == 0:
            direction = np.where(direction == 0, 0.1, direction)

        return direction / distance

    @staticmethod
    def _rotate(alpha, beta, gamma, ligand: np.array):
        cosa, cosb, cosr = math.cos(alpha), math.cos(beta), math.cos(gamma)
        sina, sinb, sinr = math.sin(alpha), math.sin(beta), math.sin(gamma)
        rotate_matrix = np.array([
            [cosa * cosr - cosb * sina * sinr, -cosb * cosr * sina - cosa * sinr, sina * sinb],
            [cosr * sina + cosa * cosb * sinr, cosa * cosb * cosr - sina * sinr, -cosa * sinb],
            [sinb * sinr, cosr * sinb, cosb],
        ])
        return np.dot(ligand, rotate_matrix)

    @staticmethod
    def _translate(target: np.array, index, ligand):
        delta = target - ligand[index]
        return ligand + delta

    def optimize(self, target, index):
        center = copy.deepcopy(self.center)
        ligand = copy.deepcopy(self.ligand)

        min_force = 10000.
        min_ligand = ligand
        for step in range(self.steps):
            ligand = Optimizer._translate(np.array(target), index, ligand)
            force = Optimizer._get_force(center, ligand)
            min_force = min(np.sum(force), min_force)
            if np.sum(force) <= self.criterion:
                break
            alpha, beta, gamma = force * self.lr * (1 + random.random())
            ligand = Optimizer._rotate(alpha, beta, gamma, ligand)
            ligand = Optimizer._translate(np.array(target), index, ligand)
            if min_force == np.sum(force):
                min_ligand = ligand
        logger.info(f"step: {step + 1} min_force: {min_force} force: {np.sum(force)}")
        return center, min_ligand


if __name__ == '__main__':
    center = np.array([0, 0, 0])
    ligand = np.array([[0.1, 0.2, 0.3], [1, 1.5, 2], [0.2, 1.3, 0]])

    optimizer = Optimizer(center, ligand)
    opt_center, opt_ligand = optimizer.optimize()
    print()
