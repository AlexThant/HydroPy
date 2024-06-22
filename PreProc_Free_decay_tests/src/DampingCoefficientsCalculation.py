# -*- coding: utf-8 -*-
# ------------------------------
# Project: 1116-EDF US shared anchors
# Created by: Maxime Chemineau - INNOSEA
# Date: 14/09/2022
# Modified by:
# ------------------------------

"""
Use this script to process OrcaFlex decay tests simulations in order to calculate decay test parameters :
- natural period,
- linear and quadratic damping,
And plot decay figures.

Method used taken from:
2010, Malta et al., DAMPING COEFFICIENT ANALYSES FOR FLOATING OFFSHORE STRUCTURES, OMAE2010-20093

Notations
---------
(M+A) d2x/dt2 + B1       dx/dt + B2       |dx/dt|*dx/dt + C    x = F
      d2x/dt2 + 2*ksi*wn dx/dt + B2/(M+A) |dx/dt|*dx/dt + wn^2 x = F
where:
    M+A   = Structural + Added mass
    B1,B2 = Damping coefficient (Linear and Quadratic)
    C     = Stifness
    ------------------------------
    wn    = sqrt(C/(M+A))    natural frequency of motion
    ksi   = B1/Bcrit         linear damping coefficient (percentage of critical damping)
    Bcrit = 2*wn*(M+A)       critical damping (corresponds to ksi = 1)
          = 2*sqrt((M+A)*C))
"""

import pandas as pd
import OrcFxAPI
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def get_peaks(sigma, decay_var):
    """Function to identify the peaks in a decay
        ---------
        Inputs:
        sigma : [DataFrame] dataframe of the decay
        decay_var : [str] name of the decay variable
        ---------
        Outputs:
        mask : [array] boolean array to get peaks indexes
        peak_values : [list] list of the peak values
        """
    # Identify peaks
    # (values where slope is positive before and negative after)
    sigma[decay_var] = sigma[decay_var] - np.mean(sigma[decay_var])
    d_sigma = np.diff(sigma[decay_var], axis=0)
    mask = (d_sigma[0:-1]*d_sigma[1:] <= 0)  # array of True/False, True for indexes where condition is met
    mask = np.concatenate([[False], mask, [False]], axis=0)  # +1 offset of d_sigma to fall back on sigma indexes
    index = 0
    peak_values = []
    for test in mask:
        if test == True and sigma.loc[index, decay_var] > 0.05:  # avoid very small amplitude peaks
            peak_values.append(sigma.loc[index, decay_var])
        elif test == True and sigma.loc[index, decay_var] < 0.05:  # keep positive peaks only
            mask[index] = False
        index += 1
    return mask, peak_values  # mask = peak indexes (as mask ie boolean array)


def calculate_decay_parameters(df, decay_var, unit):
    """ Function that calculates the linear and quadratic damping coefficients, and the natural period
    ----------
    Inputs:
    df : [DataFrame] dataframe containing the decay time series
    decay_var : [str] variable name of decay
    unit : [str] OrcaFlex result variable unit
    ----------
    Outputs:
    d_decay_results : [dict] dictionnary of the output values
    df[decay_var] : [DataFrame] dataframe of the decay
    """

    # calculate peaks
    df_mean = np.mean(df[decay_var])  # for correction after
    X_p_ind, X_p_vals = get_peaks(df, decay_var)
    df[decay_var] = df[decay_var] + df_mean  # to go back to the original values

    T0k = np.diff(df["Time"][X_p_ind])  # duration of each cycle
    T0 = np.mean(T0k)  # mean duration of cycles

    X0 = max(df[decay_var])  # max value (should be initial decay test value)


    # small condition to help detect problem in post-treatment (if coupling and peak detection issue most of the time)
    cov = np.std(T0k) / T0  # coefficient of variation = std/mean

    if cov > 0.025:  # for a regular decreasing sin wave (decay) expecting small cov, i.e similar periods between each cycle
        print
        "! Warning, check below if the variation of period between each cycle seems valid: (cov = ", '{:.2f}'.format(
            np.std(T0k) / T0 * 100), " %)"
        print
        T0k

    wn = 2 * np.pi / T0  # approximation of natural period: wn = wd/racine(1-ksi^2) FIX : we calcul here the damped natrual period
    # Estimation of the damped natural period, for undamped, wn = wd/racine(1-ksi^2)

    Tk = df["Time"][X_p_ind]  # position of each peak in time
    Xk = df[decay_var][X_p_ind]  # peaks values
    X_mean = np.mean(df[decay_var])  # mean values


    Xk = np.array(Xk) - X_mean  # suppress mean value from future treatments
    X_damp = 4 / (3 * np.pi) * Xk[1:-1]   # 4/3pi*x[k]
    Y_damp = 1 / (2 * np.pi) * np.log(Xk[0:-2] / Xk[2:])  # 1/2pi*ln(x[k-1]/x[k+1])
    if np.isnan(T0):
        m, b = 0, 0
        ksi_lin_equiv = np.zeros(len(df[decay_var]))
        ksi_lin_jordi = np.zeros(len(df[decay_var]))
    else:
        # b is ksi (linear)
        # m is Bquad/(M+Ma)
        m, b = np.polyfit(X_damp, Y_damp, 1)  # linear regression y = m*x+b
        # equivalent linear damping (considering damping is linear only)
        ksi_lin_equiv = 1 / (2 * np.pi) * np.log(Xk[0:-1] / Xk[1:])  # ksi = logarithmic_decrement/2pi
        ksi_lin_jordi = np.log(Xk[0:-1] / Xk[1:]) / np.sqrt(4 * np.pi ** 2 + (
            np.log(Xk[0:-1] / Xk[1:])) ** 2)  # ksi = logarithmic_decrement/racine(2pi^2+logaritmic decrement^2)

    print('m= ', m, '   ', 'b= ', b, '   ', 'peaks= ', X_p_vals)

    # B1/Bcrit has no unit
    # B2/Bcrit is in 1/(m/s)
    d_decay_results = {'variable': decay_var,
                       'unit': unit,
                       'wn': wn,
                       'Period': T0,
                       'Xmean': X_mean,
                       'X0': X0,
                       'B1/Bcrit': b,
                       'B2/Bcrit': m / (2 * wn),
                       'B1_equiv_min/Bcrit': np.min(ksi_lin_equiv),
                       'B1_equiv_mean/Bcrit': np.mean(ksi_lin_equiv),
                       'B1_equiv_max/Bcrit': np.max(ksi_lin_equiv),
                       'B1/Bcrit - Jordi': np.mean(ksi_lin_jordi)}

    return d_decay_results, df[decay_var]

def OrcaVesselName(model):
    """Function to get the name of a vessel in an Orcaflex model
    ----------
    Inputs:
    model = API reference to the Orcaflex model
    ----------
    Outputs:
    floater_name : [str] name of the floater in the OrcaFlex model
    """

    # Listing all objects in the model
    A = model.objects
    # Find the vessel name in the list of objects
    cont = 0
    go = True
    while go:
        # In the list A, the vessel name would come out as a tuple shaped Orcaflex object class like <Vessel: 'Name'>,
        # so we need to convert it to string and do some cleaning in the name.
        if str(A[cont]).split(':')[0].strip().replace('<', '') == 'Vessel':
            floater_name = str(A[cont]).split(':')[-1].replace('>', ' ')[2:-2]
            break
        else:
            cont += 1

    return floater_name

def runstatic(model):
    """ Function to run the static simulation of an OrcaFlex model
    ----------
    Inputs:
    model : API reference to the Orcaflex model
    ---------
    Outputs:
    [x, y, z, rx, ry, rz] : [list of float] list of the static position of the floater
    """
    floater_name = OrcaVesselName(model)
    obj = model[floater_name]
    x_initial = obj.InitialX
    y_initial = obj.InitialY
    z_initial = obj.InitialZ
    rx_initial = obj.InitialHeel
    ry_initial = obj.InitialTrim
    rz_initial = obj.InitialHeading
    obj.InitialX = 0
    obj.InitialY = 0
    obj.InitialZ = 0
    obj.InitialHeel = 0
    obj.InitialTrim = 0
    obj.InitialHeading = 0
    model.CalculateStatics()
    x = obj.TimeHistory('X', OrcFxAPI.pnStaticState, None)
    y = obj.TimeHistory('Y', OrcFxAPI.pnStaticState, None)
    z = obj.TimeHistory('Z', OrcFxAPI.pnStaticState, None)
    rx = obj.TimeHistory('Rotation 1', OrcFxAPI.pnStaticState, None)
    ry = obj.TimeHistory('Rotation 2', OrcFxAPI.pnStaticState, None)
    rz = obj.TimeHistory('Rotation 3', OrcFxAPI.pnStaticState, None)
    obj.InitialX = x_initial
    obj.InitialY = y_initial
    obj.InitialZ = z_initial
    obj.InitialHeel = rx_initial
    obj.InitialTrim = ry_initial
    obj.InitialHeading = rz_initial
    return [x, y, z, rx, ry, rz]

def launch_decays(model_path, decay, duration):
    """This function prepares and launches a decay test
    --------
    Inputs:
    model : [OrcFxAPI model] base model for the decays
    decay : [str] name of the decay to run
    --------
    Outputs:
    None
    """
    print('Preparing the ' + decay + ' decay file...')
    model = OrcFxAPI.Model(model_path)
    [x, y, z, rx, ry, rz] = runstatic(model)  # Static position of the floater
    xd, yd, zd, rxd, ryd, rzd = 0, 0, 0, 0, 0, 0  # Offsets in decays directions
    if decay == 'Surge':
        xd = 5
    elif decay == 'Heave':
        zd = 5
    elif decay == 'Pitch':
        ryd = 5
    elif decay == 'Yaw':
        rzd = 5
    floater_name = OrcaVesselName(model)
    model[floater_name].InitialX = x + xd
    model[floater_name].InitialY = y + yd
    model[floater_name].InitialZ = z + zd
    model[floater_name].InitialHeel = rx + rxd
    model[floater_name].InitialTrim = ry + ryd
    model[floater_name].InitialHeading = rz + rzd
    model[floater_name].IncludedInStatics = 'None'
    model.general.StageDuration[0] = 1
    model.general.StageDuration[1] = duration

    decay_path = model_path.replace('.dat', '_' + decay + '.dat')
    model.SaveData(decay_path)

    print('Launching the ' + decay + ' decay simulation...')
    print('Running...')
    model.RunSimulation()

    print('Saving the ' + decay + ' decay simulation...')
    decay_path = decay_path.replace('.dat', '.sim')
    model.SaveSimulation(decay_path)

def exp_decay(x, a, b, c, omega, phi, epsilon):
    """
    This is an exponential decay function that is fitted into the decay test signal to calculate the natural period
    """
    return a + b * np.cos(omega * x + phi) * np.exp(c * (x - epsilon))

def envelope(b, c, x, epsilon):
    """
    Graphic envelope of he exponential decay function, this is just for the plots
    """
    return b * np.exp(c * (x - epsilon))


if __name__ == "__main__":

    # -------------- USER INPUTS ---------------- #

    base_name = 'EDF Shared Anchor_v5_MoorDesign_staticOK_800m_decays'  # Name of the OrcaFlex .dat base model (w/o '.dat')
    sim_folder = r'C:\Users\MChemineau\Documents\1116-EDF US shared mooring\test_decays'  # Folder where the model is
    decays_duration = 1000  # [s] duration of the decay tests
    decays = ['Surge', 'Heave', 'Pitch', 'Yaw']  # decays directions to run (Possible: 'Surge', 'Heave', 'Pitch', 'Yaw')

    # ------------ END OF USER INPUTS ----------- #

    df_out = pd.DataFrame(columns=['variable',
                                   'unit',
                                   'wn',
                                   'Period',
                                   'Xmean',
                                   'X0',
                                   'B1/Bcrit',
                                   'B2/Bcrit',
                                   'B1_equiv_min/Bcrit',
                                   'B1_equiv_mean/Bcrit',
                                   'B1_equiv_max/Bcrit',
                                   'B1/Bcrit - Jordi'])

    for decay in decays:
        launch_decays(os.path.join(sim_folder, base_name + '.dat'), decay, decays_duration)

    print('Post-processing the decays...')
    simfiles = [file for file in os.listdir(sim_folder) if file.endswith('.sim')]

    timeseries = []
    index = 0
    for file in simfiles:
        print('Post-processing file : ' + file)
        model = OrcFxAPI.Model()
        model.LoadSimulation(sim_folder + '\\' + file)

        floater_name = OrcaVesselName(model)

        time_step = model.general.ImplicitConstantTimeStep
        duration_stage_1 = model.general.StageDuration[1]

        ar = np.zeros((int(duration_stage_1 / time_step) + 1, len(decays) + 1))
        df = pd.DataFrame(ar, columns=['Time']+decays)

        df.loc[:, 'Time'] = np.arange(0, duration_stage_1 + time_step, time_step)

        for obj in model.objects:
            for i in range(len(decays)):
                if obj == model[floater_name] and decays[i] in file:
                    if decays[i] == 'Surge':
                        df.loc[:, decays[i]] = obj.TimeHistory('X', 1, None)
                        unit = 'm'
                    elif decays[i] == 'Heave':
                        df.loc[:, decays[i]] = obj.TimeHistory('Z', 1, None)
                        unit = 'm'
                    elif decays[i] == 'Pitch':
                        df.loc[:, decays[i]] = obj.TimeHistory('Rotation 2', 1, None)
                        unit = 'deg'
                    elif decays[i] == 'Yaw':
                        df.loc[:, decays[i]] = obj.TimeHistory('Rotation 3', 1, None)
                        unit = 'deg'

                    # Calculate natural period and damping coefficients
                    d_decay_results, timeserie = calculate_decay_parameters(df, decays[i], unit)

                    # y0 = Initial guess for the variables of exponential decay, to help the curve_fit algorithm,
                    # just change the index of the second element to match the initial position corresponding
                    # to the degree of freedom
                    y0 = [0, 5, 0, np.round((2 * np.pi / d_decay_results['Period']), decimals=3), 0, 0]
                    # Function curve_fit uses nonlinear least squares method to fit a function f into a data set (x,y)
                    # such that popt, pcov = curve_fit(f, x, y, p0) and popt = values of the other variables of f
                    # excluding x to fit in the y data pcov = covariance of popt. p0 is an initial guess of the values
                    # of exp_decay() to help the algorithm
                    popt1, pcov1 = curve_fit(exp_decay, df['Time'], df[decays[i]], p0=y0)

                    # Plot of the decays and saving figures as .png
                    plt.plot(df['Time'], df[decays[i]], color='red')
                    plt.plot(df['Time'], exp_decay(df['Time'], *popt1), ls="dotted", color="blue", label="Curve Fit")
                    plt.plot(df['Time'], envelope(popt1[1], popt1[2], df['Time'], popt1[-1]), color="green", label="Envelope")
                    plt.plot(df['Time'], envelope(-popt1[1], popt1[2], df['Time'], popt1[-1]), color="green")
                    plt.ylabel(decays[i] + ' [' + unit + ']')
                    plt.xlabel('Time [s]')
                    plt.grid(True)
                    plt.title('Natural Period (s) = ' + str(np.round(d_decay_results['Period'], decimals=2)) + 's')
                    output_name = "decay_figure_" + decays[i]
                    print('Saving figure as ' + output_name + '.png')
                    plt.savefig(sim_folder + '\\' + output_name + '.png')
                    plt.clf()

        for k in d_decay_results:
            df_out.loc[index, k] = d_decay_results[k]

        timeseries.append(timeserie)  # Store Decay variable timeseries

        index += 1

    decay_file_out = os.path.join(sim_folder, 'decays_out.csv')
    df_out.to_csv(decay_file_out, sep='\t', index=True)