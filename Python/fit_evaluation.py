import viscofit as vf
import general as gmp
import numpy as np

root = r'data\ViscoVerification-MultiLoadLevel-ExcelConvert'
test_condition_dirs = gmp.get_folders(root)

for test_cond in test_condition_dirs:
    test_cond_data = gmp.get_files(test_cond, req_ext='csv')
    test_cond_settings = gmp.get_files(test_cond, req_ext='txt')

    # store the experimental observables for each experiment in the test condition
    fs, hs, ts, rs = [], [], [], []

    for tc_data, tc_settings in zip(test_cond_data, test_cond_settings):
        data_file = gmp.load(tc_data)
        settings_file = gmp.load(tc_settings)

        f, z, d, t = data_file['F (N)'], data_file['z (m)'], data_file['d (m)'], data_file['time (s)']
        mask = np.logical_and(f > 0, np.indices(f.shape) < np.argmax(f))[0]
        f, z, d, t = f[mask].values, z[mask].values, d[mask].values, t[mask].values
        h = d - z  # calculating the indentation as the difference between the deflection and the z-sensor
        R = float(settings_file.split(sep='Radius: ')[1].split(sep=' ')[0])  # load the tip radius
        # @TODO LOAD THE PARAMETERS
        fs.append(f), hs.append(h), ts.append(t), rs.append(R)

    # initialize the fit for the single test condition
    maxwell = vf.maxwellModel(forces=fs, indentations=hs, times=ts, radii=rs)
    voigt = vf.kelvinVoigtModel(forces=fs, indentations=hs, times=ts, radii=rs)

    # perform the fits
    maxwell.fit(maxiter=100, max_model_size=5, fit_sequential=True, num_attempts=10)
    voigt.fit(maxiter=100, max_model_size=5, fit_sequential=True, num_attempts=10)
