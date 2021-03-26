from numpy import array, sum, sqrt, convolve, exp, ones, zeros, tile, ndarray, concatenate, argmin, diff
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
    '''  #@TODO incorporate a fluidity check here
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
    '''  #@TODO change the names of the stiffnesses and time constants here
    calculates the response indentation for a kelvin-voigt model according to the lee and radok contact mechanics formulation
    :param model_params: numpy array contains the model compliances and time constants in either of the two following forms:
    to incorporate steady state fluidity: array([Elastic Compliance, Fluidity Time Constant, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
    normal model: array([Elastic Compliance, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
    :param time: (n,) numpy array time signal of the loading cycle
    :param force: (n,) numpy array force signal of the loading cycle
    :param radius: float radius of the indenter tip
    :return: (n,) numpy array indentation signal of the loading cycle
    '''
    time_matrix = row2mat(time, model_params[1::2].size)
    retardance = sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
    retardance[0] = model_params[0] / (time[1] - time[0])  # delta function
    return (3 / (8 * sqrt(radius)) * convolve(force, retardance, mode='full')[:time.size] * (time[1] - time[0])) ** (2 / 3)


class LR_Maxwell():
    def __init__(self, forces, times, indentations, radii, E_logbounds=(1, 9), T_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the LR_Maxwell class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a maxwell model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param E_logbounds: tuple (float, float) high and low log bound for the elastic elements in the model
        :param T_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.force = concatenate(forces)
            self.time = concatenate(times)
            self.indentation = concatenate(indentations)
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = concatenate([dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)])
            # create a 'mask' of radii to scale each experiment
            self.radii = concatenate([radius * ones(arr.shape) for radius, arr in zip(radii, forces)])
            # create a train of dirac delta functions with magnitude 1 at time 0 for each time signal
            diracs = []
            for t in times:
                temp = zeros(t.shape)
                temp[0] = 1
                diracs.append(temp)
            self.diracs = concatenate(diracs)
        # if there are single inputs
        else:
            self.force = forces
            self.time = times
            self.indentation = indentations
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = self.time[1] - self.time[0]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = radii
            # create a single dirac delta function with magnitude 1 at time 0
            temp = zeros(self.time.shape)
            temp[0] = 1
            self.diracs = temp
        # define the boundaries
        self.E_logbounds = E_logbounds
        self.T_logbounds = T_logbounds

    def LR_force(self, model_params):
        '''
        calculates the response force for a generalized maxwell model according to the lee and radok contact mechanics formulation
        :param model_params: numpy array contains the model stiffnesses and time constants in either of the two following forms:
        to incorporate steady state fluidity: array([Elastic Stiffness, Fluidity Time Constant, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
        normal model: array([Elastic Stiffness, Arm 1 Stiffness, Arm 1 Time Constant, ... Arm N Stiffness, Arm N Time Constant])
        :return: numpy array 'predicted' force signals for all real (experimentally obtained) indentations
        '''
        if model_params.size % 2 == 0:  # fluidity case:
            time_matrix = row2mat(self.time, model_params[0::2].size)
            relaxance = - sum(model_params[0::2] / model_params[1::2] * exp(- time_matrix / model_params[1::2]), axis=1)
            relaxance += (sum(model_params[0::2]) * self.diracs) / self.dts  # add the delta function to the relaxances
            # must divide by dt since the integral of dirac delta MUST be 1 by definiton
        else:  # no fluidity case
            time_matrix = row2mat(self.time, model_params[1::2].size)
            relaxance = - sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
            relaxance += ((model_params[0] + sum(model_params[1::2])) * self.diracs) / self.dts  # add the delta function of the relaxances
            # must divide by dt since the integral of dirac delta MUST be 1 by definiton
        return 16 * sqrt(self.radii) / 3 * convolve(self.indentation ** (3 / 2), relaxance, mode='full')[:self.time.size] * self.dts

    def SSE(self, model_params):
        '''
        gives the sum of squared errors between the 'predicted' force and real (experimentally obtained) force signals
        :param model_params: numpy array of relaxance parameters (refer to LR_force)
        :return: float sum of squared errors between the 'predicted' and real force signals
        '''
        sse = sum((self.LR_force(model_params=model_params) - self.force) ** 2, axis=0)
        if model_params.size % 2 == 0:  # fluidity case
            if (any(10 ** self.E_logbounds[0] < model_params[0::2] < 10 ** self.E_logbounds[1])
                    or any(10 ** self.T_logbounds[0] < model_params[1::2] < 10 ** self.T_logbounds[1])):
                sse *= 1e20
        else:
            if (any(10 ** self.E_logbounds[0] < model_params[1::2] < 10 ** self.E_logbounds[1])
                    or any(10 ** self.T_logbounds[0] < model_params[2::2] < 10 ** self.T_logbounds[1])
                    or (10 ** self.E_logbounds[0] < model_params[0] < 10 ** self.E_logbounds[1])):
                sse *= 1e20
        return sse

    def fit(self, maxiter=1000, max_size=4, fit_sequential=True):
        '''
        fit experimental force distance curve(s) to maxwell model of arbitrary size using a nelder-mead simplex which
        typically gives good fits rather quickly
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        for model_size in range(1, max_size + 1):  # fit without fluidity
            tic()
            initial_guess = 10 ** concatenate(([uniform(low=self.E_logbounds[0], high=self.E_logbounds[1])],
                                               concatenate([[uniform(low=self.E_logbounds[0], high=self.E_logbounds[1]),
                                                             uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size)])))
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * (model_size - 1) + 1] = data[model_size - 2][0][: 2 * (model_size - 1) + 1]
            results = minimize(self.SSE, x0=initial_guess, method='Nelder-Mead', options={'maxiter': maxiter,
                                                                                          'maxfev': maxiter,
                                                                                          'xatol': 1e-60,
                                                                                          'fatol': 1e-60})
            data.append([results.x, results.fun, toc(True)])

        for model_size in range(1, max_size + 1):  # fit with fluidity
            tic()
            initial_guess = 10 ** concatenate([[uniform(low=self.E_logbounds[0], high=self.E_logbounds[1]),
                                                uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size + 1)])
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * model_size] = data[max_size + model_size - 2][0][: 2 * model_size]
            results = minimize(self.SSE, x0=initial_guess, method='Nelder-Mead', options={'maxiter': maxiter,
                                                                                          'maxfev': maxiter,
                                                                                          'xatol': 1e-60,
                                                                                          'fatol': 1e-60})
            data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}

    def fit_slow(self, maxiter=1000, max_size=4, fit_sequential=True):
        '''
        fit experimental force distance curve(s) to maxwell model of arbitrary size using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        for model_size in range(1, max_size + 1):
            tic()
            bound = 10 ** concatenate(([self.E_logbounds[0], self.E_logbounds[1]],
                                           concatenate([[self.E_logbounds[0], self.E_logbounds[1],
                                                         self.T_logbounds[0], self.T_logbounds[1]] for n in range(model_size)]))).reshape(-1, 2).astype(float)
            initial_guess = 10 ** concatenate(([uniform(low=self.E_logbounds[0], high=self.E_logbounds[1])],
                                               concatenate([[uniform(low=self.E_logbounds[0], high=self.E_logbounds[1]),
                                                             uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size)])))
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * (model_size - 1) + 1] = data[model_size - 2][0][: 2 * (model_size - 1) + 1]
            results = dual_annealing(self.SSE, bound, maxiter=maxiter, local_search_options={'method': 'nelder-mead'}, x0=initial_guess)
            data.append([results.x, results.fun, toc(True)])

        for model_size in range(1, max_size + 1):  # fit with fluidity
            tic()
            bound = 10 ** concatenate([[self.E_logbounds[0], self.E_logbounds[1],
                                        self.T_logbounds[0], self.T_logbounds[1]] for n in range(model_size)]).reshape(-1, 2).astype(float)
            initial_guess = 10 ** concatenate([[uniform(low=self.E_logbounds[0], high=self.E_logbounds[1]),
                                                uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size + 1)])
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * model_size] = data[max_size + model_size - 2][0][: 2 * model_size]
            results = dual_annealing(self.SSE, bound, maxiter=maxiter, local_search_options={'method': 'nelder-mead'}, x0=initial_guess)
            data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}


class LR_Voigt():
    def __init__(self, forces, times, indentations, radii, J_logbounds=(-9, -1), T_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the LR_Voigt class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a kelvin-voigt model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param J_logbounds: tuple (float, float) high and low log bound for the compliance elements in the model
        :param T_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.force = concatenate(forces)
            self.time = concatenate(times)
            self.scaled_indentation = concatenate(indentations) ** (3 / 2)
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = concatenate([dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)])
            # create a 'mask' of radii to scale each experiment
            self.radii = concatenate([radius * ones(arr.shape) for radius, arr in zip(radii, forces)])
            # create a train of dirac delta functions with magnitude 1 at time 0 for each time signal
            diracs = []
            for t in times:
                temp = zeros(t.shape)
                temp[0] = 1
                diracs.append(temp)
            self.diracs = concatenate(diracs)
        # if there are single inputs
        else:
            self.force = forces
            self.time = times
            self.scaled_indentation = indentations ** (3 / 2)
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = self.time[1] - self.time[0]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = radii
            # create a single dirac delta function with magnitude 1 at time 0
            temp = zeros(self.time.shape)
            temp[0] = 1
            self.diracs = temp
        # define the boundaries
        self.J_logbounds = J_logbounds
        self.T_logbounds = T_logbounds

    def LR_scaled_indentation(self, model_params):
        '''
        calculates the scaled response indentation (h^3/2) for a generalized maxwell model according to the lee and radok contact mechanics formulation
        :param model_params: numpy array contains the model compliances and time constants in either of the two following forms:
        to incorporate steady state fluidity: array([Elastic Compliance, Fluidity Time Constant, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
        normal model: array([Elastic Compliance, Arm 1 Compliance, Arm 1 Time Constant, ... Arm N Compliance, Arm N Time Constant])
        :return: numpy array scaled 'predicted' indentation signals (h^3/2) for all real (experimentally obtained) forces
        '''
        if model_params.size % 2 == 0:  # fluidity case
            time_matrix = row2mat(self.time, model_params[0::2].size)
            retardance = sum(model_params[2::2] / model_params[3::2] * exp(- time_matrix / model_params[3::2]), axis=1) + model_params[1]
            retardance += (model_params[0] * self.diracs) / self.dts  # add the delta function to the relaxances
            # must divide by dt since the integral of dirac delta MUST be 1 by definiton
        else:  # no fluidity case
            time_matrix = row2mat(self.time, model_params[1::2].size)
            retardance = sum(model_params[1::2] / model_params[2::2] * exp(- time_matrix / model_params[2::2]), axis=1)
            retardance += (model_params[0] * self.diracs) / self.dts  # add the delta function to the relaxances
            # must divide by dt since the integral of dirac delta MUST be 1 by definiton
        return 3 / (8 * sqrt(self.radii)) * convolve(self.force, retardance, mode='full')[:self.time.size] * self.dts

    def SSE(self, model_params):
        '''
        gives the sum of squared errors between the scaled 'predicted' indentation and real scaled (experimentally obtained) indentation signals (h^3/2)
        :param model_params: numpy array of retardance parameters (refer to LR_force)
        :return: float sum of squared errors between the scaled 'predicted' and real indentation signals (h^3/2)
        '''
        sse = sum((self.LR_scaled_indentation(model_params=model_params) - self.scaled_indentation) ** 2, axis=0)
        if model_params.size % 2 == 0:  # fluidity case
            if (any(10 ** self.J_logbounds[0] < model_params[0::2] < 10 ** self.J_logbounds[1])
                    or any(10 ** self.T_logbounds[0] < model_params[1::2] < 10 ** self.T_logbounds[1])):
                sse *= 1e20
        else:
            if (any(10 ** self.J_logbounds[0] < model_params[1::2] < 10 ** self.J_logbounds[1])
                    or any(10 ** self.T_logbounds[0] < model_params[2::2] < 10 ** self.T_logbounds[1])
                    or (10 ** self.J_logbounds[0] < model_params[0] < 10 ** self.J_logbounds[1])):
                sse *= 1e20
        return sse

    def fit(self, maxiter=1000, max_size=4, fit_sequential=True):
        '''
        fit experimental force distance curve(s) to kelvin-voigt model of arbitrary size using a nelder-mead simplex which
        typically gives good fits rather quickly
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        for model_size in range(1, max_size + 1):  # fit without fluidity
            tic()
            initial_guess = 10 ** concatenate(([uniform(low=self.J_logbounds[0], high=self.J_logbounds[1])],
                                               concatenate([[uniform(low=self.J_logbounds[0], high=self.J_logbounds[1]),
                                                             uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size)])))
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * (model_size - 1) + 1] = data[model_size - 2][0][: 2 * (model_size - 1) + 1]
            results = minimize(self.SSE, x0=initial_guess, method='Nelder-Mead', options={'maxiter': maxiter,
                                                                                          'maxfev': maxiter,
                                                                                          'xatol': 1e-60,
                                                                                          'fatol': 1e-60})
            data.append([results.x, results.fun, toc(True)])

        for model_size in range(1, max_size + 1):  # fit with fluidity
            tic()
            initial_guess = 10 ** concatenate([[uniform(low=self.J_logbounds[0], high=self.J_logbounds[1]),
                                                uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size + 1)])
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * model_size] = data[max_size + model_size - 2][0][: 2 * model_size]
            results = minimize(self.SSE, x0=initial_guess, method='Nelder-Mead', options={'maxiter': maxiter,
                                                                                          'maxfev': maxiter,
                                                                                          'xatol': 1e-60,
                                                                                          'fatol': 1e-60})
            data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}

    def fit_slow(self, maxiter=1000, max_size=4, fit_sequential=True):
        '''
        fit experimental force distance curve(s) to kelvin-voigt model of arbitrary size using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :param max_size: int largest number of arms per maxwell model to test (going larger tends to give poor and unphysical fits)
        :param fit_sequential: bool whether or not to fit sequentially (cascade fit from previous model as the initial guess of the next) (RECOMMENDED)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        for model_size in range(1, max_size + 1):
            tic()
            bound = 10 ** concatenate(([self.J_logbounds[0], self.J_logbounds[1]],
                                       concatenate([[self.J_logbounds[0], self.J_logbounds[1],
                                                     self.T_logbounds[0], self.T_logbounds[1]] for n in range(model_size)]))).reshape(-1, 2).astype(float)
            initial_guess = 10 ** concatenate(([uniform(low=self.J_logbounds[0], high=self.J_logbounds[1])],
                                               concatenate([[uniform(low=self.J_logbounds[0], high=self.J_logbounds[1]),
                                                             uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size)])))
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * (model_size - 1) + 1] = data[model_size - 2][0][: 2 * (model_size - 1) + 1]
            results = dual_annealing(self.SSE, bound, maxiter=maxiter, local_search_options={'method': 'nelder-mead'}, x0=initial_guess)
            data.append([results.x, results.fun, toc(True)])

        for model_size in range(1, max_size + 1):  # fit with fluidity
            tic()
            bound = 10 ** concatenate([[self.J_logbounds[0], self.J_logbounds[1],
                                        self.T_logbounds[0], self.T_logbounds[1]] for n in range(model_size)]).reshape(-1, 2).astype(float)
            initial_guess = 10 ** concatenate([[uniform(low=self.J_logbounds[0], high=self.J_logbounds[1]),
                                                uniform(low=self.T_logbounds[0], high=self.T_logbounds[1])] for n in range(model_size + 1)])
            if model_size != 1 and fit_sequential:  # for all guesses past the first guess, use the results from the previous fit as the initial guess
                initial_guess[: 2 * model_size] = data[max_size + model_size - 2][0][: 2 * model_size]
            results = dual_annealing(self.SSE, bound, maxiter=maxiter, local_search_options={'method': 'nelder-mead'}, x0=initial_guess)
            data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}


class LR_PowerLaw():
    def __init__(self, forces, times, indentations, radii, E0_logbounds=(-9, -1), a_logbounds=(-5, 0)):  #@TODO add conical and flat punch indenter options
        '''
        initializes an instance of the Custom_Model class
        used for generating fits, of experimentally obtained force-distance data all belonging to the same sample,
        to a power law rheology model which corresponds to the sample's viscoelastic behavior
        :param forces: either list of numpy arrays or single numpy array corresponding to the force signals from an AFM
        :param times: either list of numpy arrays or single numpy array corresponding to the time signals from an AFM
        :param indentations: either list of numpy arrays or single numpy array corresponding to the indentation signals from an AFM
        :param radii: either list of floats or single float corresponding to the tip radii of an AFM
        :param E0_logbounds: tuple (float, float) high and low log bound for the compliance elements in the model
        :param a_logbounds: tuple (float, float) high and low log bound for the time constants in the model
        '''
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.force = concatenate(forces)
            self.time = concatenate(times)
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = concatenate([dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)])
            # create a 'mask' of radii to scale each experiment
            self.radii = concatenate([radius * ones(arr.shape) for radius, arr in zip(radii, forces)])
            # create a train of dirac delta functions with magnitude 1 at time 0 for each time signal
            diracs = []
            for t in times:
                temp = zeros(t.shape)
                temp[0] = 1
                diracs.append(temp)
            self.diracs = concatenate(diracs)
            # calculate the numerical derivatives of the scaled indentations (h^3/2)
            self.scaled_indentations_deriv = concatenate([concatenate(([0], diff(indentation ** (3 / 2)))) for indentation in indentations]) / self.dts
        # if there are single inputs
        else:
            self.force = forces
            self.time = times
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = self.time[1] - self.time[0]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = radii
            # create a single dirac delta function with magnitude 1 at time 0
            temp = zeros(self.time.shape)
            temp[0] = 1
            self.diracs = temp
            # calculate the numerical derivatives of the scaled indentations (h^3/2)
            self.scaled_indentations_deriv = concatenate(([0], diff(indentations ** (3 / 2)))) / self.dts
        # define the boundaries
        self.E0_logbounds = E0_logbounds
        self.a_logbounds = a_logbounds

    def LR_force(self, model_params):
        '''
        calculates the force respones for a generalized maxwell model according to the lee and radok contact mechanics formulation
        :param model_params: numpy array contains the model's instantaneous relaxation (E0) and power law exponent (a) in the following form
        array([E0, a])
        :return: numpy array scaled 'predicted' force signals for all real (experimentally obtained) forces
        '''
        relaxation = model_params[0] * (1 + self.time / self.dts) ** (- model_params[1])
        return 16 * sqrt(self.radii) / 3 * convolve(self.scaled_indentations_deriv, relaxation, mode='full')[:self.time.size] * self.dts

    def SSE(self, model_params):
        '''
        gives the sum of squared errors between the scaled 'predicted' indentation and real scaled (experimentally obtained) force signals
        :param model_params: numpy array of model parameters (refer to LR_force)
        :return: float sum of squared errors between the scaled 'predicted' and real indentation signals (h^3/2)
        '''
        sse = sum((self.LR_force(model_params=model_params) - self.force) ** 2, axis=0)
        if (any(self.E0_logbounds[0] < model_params[0] < self.E0_logbounds[1])
                or any(self.a_logbounds[0] < model_params[1] < self.a_logbounds[1])):
            sse *= 1e20
        return sse

    def fit(self, maxiter=1000):
        '''
        fit experimental force distance curve(s) to power law rheology model using a nelder-mead simplex which typically gives good fits rather quickly
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        tic()
        initial_guess = 10 ** concatenate([[uniform(low=self.E0_logbounds[0], high=self.E0_logbounds[1]),
                                            uniform(low=self.a_logbounds[0], high=self.a_logbounds[1])]])
        results = minimize(self.SSE, x0=initial_guess, method='Nelder-Mead', options={'maxiter': maxiter,
                                                                                      'maxfev': maxiter,
                                                                                      'xatol': 1e-60,
                                                                                      'fatol': 1e-60})
        data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}

    def fit_slow(self, maxiter=1000):
        '''
        fit experimental force distance curve(s) to power law rheology model using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        tic()
        bound = 10 ** concatenate(([self.E0_logbounds[0], self.E0_logbounds[1]],
                                   [self.a_logbounds[0], self.a_logbounds[1]])).reshape(-1, 2).astype(float)
        initial_guess = 10 ** concatenate([[uniform(low=self.E0_logbounds[0], high=self.E0_logbounds[1]),
                                            uniform(low=self.a_logbounds[0], high=self.a_logbounds[1])]])
        results = dual_annealing(self.SSE, bound, maxiter=maxiter, local_search_options={'method': 'nelder-mead'}, x0=initial_guess)
        data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}


#@TODO add initial guesses and bounds to both nelder mead and dual annealing
class Custom_Model():
    def __init__(self, forces, times, indentations, radii):
        # if there are multiple inputs
        if type(forces) is list:
            # check for any size mismatches
            if any([len(arr) != len(forces) for arr in (times, indentations, radii)]):
                exit('Error: Size Mismatch in Experimental Observables!  All experimental observables must be the same size!')
            # concatenate the lists of experimental observables to put them into a single row vector form
            self.force = concatenate(forces)
            self.time = concatenate(times)
            self.indentation = concatenate(indentations)
            # create a 'mask' of dt to properly integrate each experiment
            self.dts = concatenate([dt * ones(arr.shape) for dt, arr in zip([t[1] - t[0] for t in times], times)])
            # create a 'mask' of radii to scale each experiment
            self.radii = concatenate([radius * ones(arr.shape) for radius, arr in zip(radii, forces)])
        # if there are single inputs
        else:
            self.force = forces
            self.time = times
            self.indentation = indentations
            # dt is a single value rather than a 'mask' array as seen above
            self.dts = self.time[1] - self.time[0]
            # radius is a single value rather than a 'mask' array as seen above
            self.radii = radii
        # defining the currently empty target observable and the observable function
        self.target_observable = None
        self.observable_function = None

    def SSE(self, params):
        return sum((self.observable_function(params) - self.target_observable) ** 2, axis=0)

    def fit(self, maxiter=1000):
        '''
        fit experimental observable of your choice to a custom model for the observable using a nelder-mead simplex which typically gives good fits rather quickly
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        tic()
        results = minimize(self.SSE, method='Nelder-Mead', options={'maxiter': maxiter,
                                                                    'maxfev': maxiter,
                                                                    'xatol': 1e-60,
                                                                    'fatol': 1e-60})
        data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}

    def fit_slow(self, maxiter=1000):
        '''
        fit experimental observable of your choice to a custom model for the observable using simulated annealing with
        a nelder-mead simplex local search, this is very computationally costly and will take a very long time
        though typically results in much better fits
        :param maxiter: int maximum iterations to perform for each fitting attempt (larger number gives longer run time)
        :return: dict {best_fit, (numpy array of final best fit params),
                       final_cost, (float of final cost for the best fit params),
                       time, (float of time taken to generate best fit)}
        '''
        data = []
        tic()
        results = dual_annealing(self.SSE, maxiter=maxiter, local_search_options={'method': 'nelder-mead'})
        data.append([results.x, results.fun, toc(True)])

        data = array(data)
        best_fit = data[argmin(data[:, 1])]
        return {'final_params': best_fit[0], 'final_cost': best_fit[1], 'time': best_fit[2]}


#@TODO test bounds
#@TODO fix bounds to be programmatically better
#@TODO add many fits with statistics
#@TODO add conical and flat punch indenter options
#@TODO add the ibw reader
