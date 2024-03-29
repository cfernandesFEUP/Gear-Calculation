## GEAR TYPE ##
import numpy as np
def gtype(gear):
    if gear == 'C14':
        alpha = 20.0
        beta = 0.0
        m = 4.5
        z = np.array([16., 24.])
        x = np.array([0.1817, 0.1715])
        b = 14.0
        dsh = 30.0
    if gear == '2020':
        alpha = 20.0
        beta = 0.0
        m = 4.5
        z = np.array([20., 20.])
        x = np.array([0.1766, 0.1766])
        b = 14.0
        dsh = 30.0
    elif gear == 'C40':
        alpha = 20.0
        beta = 0.0
        m = 4.5
        z = np.array([16., 24.])
        x = np.array([0.1817, 0.1715])
        b = 40.0
        dsh = 30.0
    elif gear == '501':
        alpha = 20.0
        beta = 15.0
        m = 3.5
        z = np.array([20., 30.])
        x = np.array([0.1809, 0.0891])
        b = 23.0
        dsh = 30.0
    elif gear == '701':
        alpha = 20.0
        beta = 15.0
        m = 2.5
        z = np.array([28., 42.])
        x = np.array([0.2290, 0.1489])
        b = 17.0
        dsh = 30.0
    elif gear == '951':
        alpha = 20.0
        beta = 15.0
        m = 1.75
        z = np.array([38., 57.])
        x = np.array([1.6915, 2.0003])
        b = 21.2418
        dsh = 30.0
    elif gear == 'M951':
        alpha = 20.0
        beta = 15.0
        m = 1.75
        z = np.array([38., 57.])
        x = np.array([1.5409, 2.1509])
        b = 21.2418
        dsh = 30.0
    elif gear == 'TPA':
        alpha = 20.0 
        beta = 0.0 
        m = 4.5
        z = np.array([16., 24.])
        x = np.array([0.8532, -0.5])
        b = 10.0
        dsh = 30.0
    elif gear == 'EEE':
        alpha = 20.0
        beta = 0.0
        m = 1.75
        z = np.array([38., 57.])
        x = np.array([0., 0.])
        b = 20
        dsh = 20.0
    elif gear == 'OM6':
        alpha = 20.0
        beta = 30.0
        m = 2.
        z = np.array([20., 41.])
        x = np.array([0., 0.])
        b = 20
        dsh = 20.0
    elif gear == 'COM':
        alpha = 20.0
        beta = 0.0
        m = 2.
        z = np.array([25., 40.])
        x = np.array([0.15539, -0.15539])
        b = 20
        dsh = 20.0
    return alpha, beta, m, z, x, b, dsh
