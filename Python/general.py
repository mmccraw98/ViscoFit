from numpy import array, sum, sqrt, convolve, exp, ones, cos, dot, pi, arccos, kron, tile, abs, log10, insert, linalg, indices, ndarray, concatenate, argmin
from numpy.random import uniform
from scipy.optimize import minimize, dual_annealing
from time import time


def row2mat(row, n):
    '''
    stacks a row vector (numpy (m, )) n times to create a matrix (numpy (m, n)) NOTE: CAN SLOW DOWN COMPUTATION IF DONE MANY TIMES
    :param row: numpy array row vector
    :param n: int number of replications to perform
    :return: numpy matrix (m, n) replicated row vector
    '''
    # do once at the beginning of any calculation to improve performance
    return tile(row, (n, 1)).T


def tic():
    '''
    start the clock
    :return: a global var for the initial time
    '''
    global current_time
    current_time = time()


def toc(return_numeric=False):
    '''
    stop the clock and return the elapsed time
    :param return_numeric: bool whether or not to return a numeric value or a report string
    :return: numeric value of elapsed time or a report string for elapsed time
    '''
    if return_numeric:
        return time() - current_time
    else:
        print('process completed in {:.2f}s'.format(time() - current_time))


def LR_Maxwell_Force(model_params, time, indentation, radius):
    '''#@TODO incorporate a fluidity check here
    calculates the response force for a generalized maxwell model according to the lee and radok contact mechanics formulation
    :param model_params: numpy array contains the model stiffnesses and time constants in either of the two following forms:
    to incorporate steady state fluidity: array([Elastic Stiffness, Fluidity Time Constant, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
    normal model: array([Elastic Stiffness, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
    :param time: (n,) numpy array time signal of the loading cycle
    :param indentation: (n,) numpy array indentation signal of the loading cycle
    :param radius: float radius of the indenter tip
    :return: (n,) numpy array force signal of the loading cycle
    '''
    time_matrix = row2mat(time, model_params[1::2].size)
    relaxance = - sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
    relaxance[0] = (model_params[0] + sum(model_params[1::2])) / (time[1] - time[0])  # delta function
    return 16 * sqrt(radius) / 3 * convolve(indentation ** (3 / 2), relaxance, mode='full')[:time.size] * (time[1] - time[0])


def LR_Voigt_Indentation(model_params, time, force, radius):
    ''' #@TODO change the names of the stiffnesses and time constants here
    calculates the response indentation for a kelvin-voigt model according to the lee and radok contact mechanics formulation
    :param model_params: numpy array contains the model stiffnesses and time constants in either of the two following forms:
    to incorporate steady state fluidity: array([Elastic Stiffness, Fluidity Time Constant, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
    normal model: array([Elastic Stiffness, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
    :param time: (n,) numpy array time signal of the loading cycle
    :param force: (n,) numpy array force signal of the loading cycle
    :param radius: float radius of the indenter tip
    :return: (n,) numpy array indentation signal of the loading cycle
    '''
    time_matrix = row2mat(time, model_params[1::2].size)
    retardance = sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
    retardance[0] = (model_params[0] + sum(model_params[1::2])) / (time[1] - time[0])  # delta function
    return (3 / (8 * sqrt(radius)) * convolve(force, retardance, mode='full')[:time.size] * (time[1] - time[0])) ** (2 / 3)


class LR_Maxwell():
    def __init__(self, forces, times, indentations, radii):
        if type(forces) not in (ndarray, list):
            forces = [forces]
            times = [times]
            indentations = [indentations]
            radii = [radii]
        if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
            exit('Error: Size Mismatch in Experimental Observables')
        self.force = concatenate(forces)
        self.time = concatenate(times)
        # create a 'mask' of dt to properly integrate each experiment
        self.dts = concatenate([dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)])
        self.indentation = concatenate(indentations)
        # create a 'mask' of radii to scale each experiment
        self.radii = concatenate([radius * ones(arr.shape) for radius, arr in zip(radii, forces)])

    def LR_force(self, relaxance):
        time_matrix = row2mat(self.time, relaxance[1::2].size)
        relaxance = - sum(relaxance[1::2] / relaxance[2::2] * exp(- time_matrix / relaxance[2::2]), axis=1)
        relaxance[0] = (relaxance[0] + sum(relaxance[1::2])) / (time[1] - time[0])  # delta function
        return 16 * sqrt(self.radii) / 3 * convolve(self.indentation ** (3 / 2), relaxance, mode='full')[:self.time.size] * self.dts

    def LR_force_fluidity(self, relaxance): #@TODO implement fluidity here as well as fit and fit_slow
        raise NotImplementedError

    def SSE_fluidity(self, relaxance): #@TODO implement fluidity here as well as fit and fit_slow
        raise NotImplementedError

    def SSE(self, relaxance):
        return sum((self.LR_force(relaxance=relaxance) - self.force) ** 2, axis=0)

    def fit(self, maxiter=1000, max_size=4, E_logbounds=(1, 9), T_logbounds=(-5, 0), fit_sequential=True): #@TODO implement sequential fitting and fluidity
        if fit_sequential:
            pass
        else:
            data = []
            for model_size in range(1, max_size + 1):
                tic()
                initial_guess = concatenate(([uniform(low=E_logbounds[0], high=E_logbounds[1])],
                                             concatenate([[uniform(low=E_logbounds[0], high=E_logbounds[1]),
                                                           uniform(low=T_logbounds[0], high=T_logbounds[1])] for n in range(model_size)])))
                results = minimize(self.SSE, x0=initial_guess, method='Nelder-Mead', options={'maxiter': maxiter,
                                                                                              'maxfev': maxiter,
                                                                                              'xatol': 1e-60,
                                                                                              'fatol': 1e-60})
                data.append([results.x, results.fun, toc(True)])
            data = array(data)
            best_fit = data[argmin(data[:, 1])]
            return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}

    def fit_slow(self, maxiter=1000, max_size=4, E_logbounds=(1, 9), T_logbounds=(-5, 0), fit_sequential=True):  # @TODO implement sequential fitting and fluidity
        if fit_sequential:
            pass
        else:
            data = []
            for model_size in range(1, max_size + 1):
                tic()
                bound = 10 ** concatenate(([E_logbounds[0], E_logbounds[1]],
                                           concatenate([[E_logbounds[0], E_logbounds[1],
                                                         T_logbounds[0], T_logbounds[1]] for n in range(model_size)]))).reshape(-1, 2).astype(float)
                results = dual_annealing(self.SSE, bound, maxiter=maxiter, local_search_options={'method': 'nelder-mead'})
                data.append([results.x, results.fun, toc(True)])
            data = array(data)
            best_fit = data[argmin(data[:, 1])]
            return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}

#@TODO Add Voigt

#@TODO Add Power Law
