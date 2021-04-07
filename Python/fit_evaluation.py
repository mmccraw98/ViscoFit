import viscofit as vf
import general as gmp
import numpy as np
import os
from time import time

# set the root for the data
root = os.path.join('data', 'ViscoVerification-MultiLoadLevel-ExcelConvert')
test_condition_dirs = gmp.get_folders(root)

# loop over all test conditions
for i, test_cond in enumerate(test_condition_dirs):
    print('Fitting on test condition data {} of {}'.format(i + 1, len(test_condition_dirs)))

    test_cond_data = gmp.get_files(test_cond, req_ext='csv')
    test_cond_settings = gmp.get_files(test_cond, req_ext='txt')

    # store the experimental observables for each experiment in the test condition
    fs, hs, ts, rs = [], [], [], []

    for tc_data, tc_settings in zip(test_cond_data, test_cond_settings):
        data_file = gmp.load(tc_data)
        settings_file = gmp.load(tc_settings)
        Ee = float(settings_file.split(sep='E0: ')[1].split(sep=' ')[0])

        # calculate relaxance and retardance parameters
        Es, Ts, Js = [], [], []
        for arm_data in settings_file.split(sep='tau (s)')[1].split(sep='\n')[1:-1]:
            row_data = [float(value) for value in arm_data.split(sep=' ') if value != '']
            Es.append([row_data[1], 0])
            Js.append([row_data[2], 0])
            Ts.append([0, row_data[3]])
        relaxance_params = np.concatenate(([Ee], np.array(Es).ravel() + np.array(Ts).ravel()))
        # omit retardance for now, it is not calculated in this way
        # will test the retardance param error as a post-process
        # retardance_params = np.concatenate(([np.sum(Js) + 1 / Ee], np.array(Js).ravel() + np.array(Ts).ravel()))

        # get the experimental data from the files
        f, z, d, t = data_file['F (N)'], data_file['z (m)'], data_file['d (m)'], data_file['time (s)']
        mask = np.logical_and(f > 0, np.indices(f.shape) < np.argmax(f))[0]
        f, z, d, t = f[mask].values, z[mask].values, d[mask].values, t[mask].values
        h = d - z  # calculating the indentation as the difference between the deflection and the z-sensor
        R = float(settings_file.split(sep='Radius: ')[1].split(sep=' ')[0])  # load the tip radius
        fs.append(f), hs.append(h), ts.append(t), rs.append(R)

    # start the timer
    start = time()

    # initialize the fit for the single test condition
    print('---Maxwell')
    maxwell = vf.maxwellModel(forces=fs, indentations=hs, times=ts, radii=rs)
    print('--Voigt')
    voigt = vf.kelvinVoigtModel(forces=fs, indentations=hs, times=ts, radii=rs)
    print('--Power Law')
    power = vf.powerLawModel(forces=fs, indentations=hs, times=ts, radii=rs)

    # perform the fits
    relaxance_fit = maxwell.fit(maxiter=1000, max_model_size=5, fit_sequential=True, num_attempts=100)
    retardance_fit = voigt.fit(maxiter=1000, max_model_size=5, fit_sequential=True, num_attempts=100)
    power_fit = power.fit(maxiter=1000, num_attempts=100)

    import matplotlib.pyplot as plt
    for f_, t_, h_, r_ in zip(fs, ts, hs, rs):
        plt.plot(vf.forceMaxwell_LeeRadok(relaxance_fit['model_data']['final_params'], t_, h_, r_), label='pred', linestyle='--')
        plt.plot(vf.forceMaxwell_LeeRadok(relaxance_fit['fluidity_data']['final_params'], t_, h_, r_), label='pred_fluid', linestyle=':')
        plt.plot(f, label='real', linestyle='-')
    plt.title('maxwell force comparison')
    plt.legend()
    plt.grid()
    plt.show()

    for f_, t_, h_, r_ in zip(fs, ts, hs, rs):
        plt.plot(vf.indentationKelvinVoigt_LeeRadok(retardance_fit['model_data']['final_params'], t_, f_, r_), label='pred', linestyle='--')
        plt.plot(vf.indentationKelvinVoigt_LeeRadok(retardance_fit['fluidity_data']['final_params'], t_, f_, r_), label='pred_fluid', linestyle=':')
        plt.plot(h, label='real', linestyle='-')
    plt.title('voigt indentation comparison')
    plt.legend()
    plt.grid()
    plt.show()

    for f_, t_, h_, r_ in zip(fs, ts, hs, rs):
        plt.plot(vf.forcePowerLaw_LeeRadok(power_fit['final_params'], t_, h_, r_), label='pred', linestyle='--')
        plt.plot(f, label='real', linestyle='-')
    plt.title('power force comparison')
    plt.legend()
    plt.grid()
    plt.show()

    quit()

    # compare the relaxance_params and relaxance_fit
    # compare the retardance_params and retardance_fit
    gmp.safesave(relaxance_params, os.path.join(test_cond, 'relaxance_real.pkl'), overwrite=True)

    gmp.safesave(relaxance_fit, os.path.join(test_cond, 'simultaneous_relaxance_fit.pkl'), overwrite=True)
    gmp.safesave(retardance_fit, os.path.join(test_cond, 'simultaneous_retardance_fit.pkl'), overwrite=True)
    gmp.safesave(power_fit, os.path.join(test_cond, 'simultaneous_power_fit.pkl'), overwrite=True)

    # stop the timers
    print('Time Elapsed: {}s'.format(time() - start))

    # get the error in the harmonic quantities @TODO

    # get the error between the parameters @TODO
