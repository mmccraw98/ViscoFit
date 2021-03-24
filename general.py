import numpy as np
from numpy import array, sum, sqrt, convolve, exp, ones, cos, dot, pi, arccos, kron, tile, abs, log10, insert, linalg, indices
from scipy.optimize import minimize
from matplotlib import pyplot as plt


def row2mat(row, n):
    '''
    stacks a row vector (numpy (m, )) n times to create a matrix (numpy (m, n)) NOTE: CAN SLOW DOWN COMPUTATION IF DONE MANY TIMES
    :param row: numpy array row vector
    :param n: int number of replications to perform
    :return: numpy matrix (m, n) replicated row vector
    '''
    # do once at the beginning of any calculation to improve performance
    return tile(row, (n, 1)).T


def LR_Maxwell(model_params, time, indentation, radius):
    time_matrix = row2mat(time, model_params[1::2].size)
    relaxance = - sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
    relaxance[0] = (model_params[0] + sum(model_params[1::2])) / (time[1] - time[0])
    return 16 * sqrt(radius) / 3 * convolve(indentation ** (3 / 2), relaxance, mode='full')[:time.size] * (time[1] - time[0])


def LR_Voigt(model_params, time, force, radius):
    time_matrix = row2mat(time, model_params[1::2].size)
    retardance = sum(model_params[1::2] / model_params[2::2] * exp( - time_matrix / model_params[2::2]), axis=1)
    retardance[0] = (model_params[0] + sum(model_params[1::2])) / (time[1] - time[0])
    return 3 / (8 * sqrt(radius)) * convolve(force ** (3 / 2), retardance, mode='full')[:time.size] * (time[1] - time[0])


class SSE_simultaneous_gen_maxwell():
    def __init__(self, forces, times, indentations, radii):
        self.forces = forces
        self.times = times
        self.indentations = indentations
        self.radii = radii

    def test_function(self, model_params):
        # calculate test force data
        test_forces = np.array([LR_Maxwell(model_params, t, h, R) for t, h, R in
                                zip(self.times, self.indentations, self.radii)])
        # calculate global residual
        residual_global = self.forces - test_forces
        # return sse
        sse = np.sum(residual_global**2)
        if (np.where(np.logical_and(model_params[0] <= 1e1, model_params[0] >= 1e9))[0].size > 0) or (np.where(np.logical_and(model_params[1::2] <= 1e1, model_params[1::2] >= 1e9))[0].size > 0) or (np.where(np.logical_and(model_params[2::2] <= 1e-6, model_params[2::2] >= 1e0))[0].size > 0):
            return sse * 1e20
        return sse


t = np.arange(0, 1, 1e-4)
Q_array = np.array([1e6, 1e4, 1e-3, 1e7, 1e-1, 1e6, 1e-2])
J_array = np.array([1e-1, 1e-1, 1e-3])
h = t / max(t) * 50e-9
R = 100e-9
t_matrix = row2mat(t, Q_array[2::2].size)


plt.plot(t, LR_Maxwell(Q_array, t, h / 10, R))
# plt.xscale('log')
# plt.yscale('log')
plt.show()
