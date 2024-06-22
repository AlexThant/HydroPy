"""
author@ Thant Zin Htun
Date@ 10/20/2023
Location@  Tokyo, Japan
Python_version tested@ 3.11.6 & 3.12.0
\\\
This module is to automate freedecay runs using multiprocessing through OrcaFlex API.
This is intended not to be called from another program.
For bug report: thant.htun@owcltd.com
\\\
"""
import os, sys, glob
from scipy.optimize import curve_fit
#import matplotlib; matplotlib.use('TkAgg')  # for compatibility on Mac OS
import matplotlib.pyplot as plt              # for generating plots & graphs
import numpy as np
import OrcFxAPI as OC


plt.rcParams['figure.figsize']=[15,10]
plt.rcParams.update({'font.size': 10})
SpectralDensityPlots = True
PlotLabelList = [r'Surge (m)', r'Sway (m)', r'Heave (m)', r'Roll (deg)', r'Pitch (deg)', r'Yaw (deg)']

#______________________________________________________________________________

class PostProcessing:

    def __init__(self, ResultsFolder, TransientTime, TEnd_Stage0, TimeStep, InitialDisplacements, FloaterName, ObjectType, Position, CalculateDoFs, ResultsReportStyle= 'default'):
        self.ResultsFolder = ResultsFolder
        self.TransientTime, self.TEnd_Stage0, self.TimeStep  = TransientTime, TEnd_Stage0, TimeStep
        self.InitialDisplacements = InitialDisplacements
        self.FloaterName, self.ObjectType, self.Position = FloaterName, ObjectType, Position

        self.RefPeriod = 20            # Estimate of natural periods to be found, in seconds

        self.SimFileList = [r'Surge.sim', r'Sway.sim', r'Heave.sim', r'Roll.sim', r'Pitch.sim', r'Yaw.sim']
        self.DoFList = ['X', 'Y', 'Z', 'Rotation 1', 'Rotation 2', 'Rotation 3']


        #for sim in glob.glob(ResultsFolder + r'\*.sim'):
            #if sim not in SimFileList:
            #    raise Exception('ERROR: File name error (or) .sim file does not exist!')    # This is sloppy. I should modify this algorithm.
        #else:
        #    pass

        if (ResultsReportStyle == 'default'and len(CalculateDoFs) == 6):
                self.PlotResultAllDoFs()
                if (SpectralDensityPlots == True):
                    self.PlotFreqSpectralDensity()

        else:      
            for Index in CalculateDoFs:
                SimFileIndex = Index-1
                SimFile = self.LoadSimFile(self.SimFileList[SimFileIndex])    
                self.PlotResultsEachDoF(SimFileIndex, self.DoFList[SimFileIndex], SimFile)


        return


    def LoadSimFile(self, SimFileName):
        return os.path.join(self.ResultsFolder, SimFileName)


    def ExponentialDecay(self, x, a, b, c, omega, phi, epsilon):
        """ This is an exponential decay function fitting onto the decay test signal to calculate the natural
        period """
        return a + b * np.cos(omega * x + phi) * np.exp(c * (x - epsilon))


    def InitialGuessExpoDecay(self, iD):
        """ Set the values of initial guesses for exponential decay for the curve fit.
        """
        y0 = [0, self.InitialDisplacements[iD], 0, np.round((2.0 * np.pi / self.RefPeriod), decimals= 3),0, 10e-6]
        return y0


    @staticmethod
    def Envelope(b, c, x, epsilon):
        """ Graphic envelope of the exponential decay function, this is just for the plots"""
        return b * np.exp(c * (x - epsilon))


    def PeriodBySpectralDensity(self, SelectedDOF, SimFile):
        """
        To calculate the natural period of the system by using the spectral density plots of OrcaFlex.
        """
        
        SimModel = OC.Model(SimFile)
        Floater  = SimModel[self.FloaterName]
        Freq, Spectra = Floater.SpectralDensity(SelectedDOF)

        # Find the index position of the frequency corresponding to the maximum spectral density.
        MaxSpectralIndex = np.argmax(Spectra)

        return 1.0/Freq[MaxSpectralIndex], Freq, Spectra

        #___*___*___*___*___*___*___*___*___*___*___*___*

    def NaturalPeriodCurveFit(self, Xdata, Ydata, iD):
        """
        To calculate the natural period of the system by using the CurveFit.
        """  
        popt, pcov = curve_fit(self.ExponentialDecay, Xdata, Ydata, p0=(self.InitialGuessExpoDecay(iD)))
        Omega = popt[3];   #Check the element location of returned frequency value in y0 inside InitialGuessExpoDecay function.
       
        return np.round((2.0 * np.pi / Omega), decimals=3), Xdata, Ydata, popt
    

    def TimeHistoryResults(self, iD, SelectedDOF, SimFile):
        """
        To calculate the natural period of the system by using the CurveFit.
        """
        SimModel = OC.Model(SimFile)
        Floater  = SimModel[self.FloaterName]
        Time = SimModel.general.TimeHistory('Time')

        try:
            if (self.ObjectType in ['Buoy', 'buoy']): Displacement = Floater.TimeHistory(SelectedDOF, OC.pnWholeSimulation, OC.oeBuoy(self.Position[0], self.Position[1], self.Position[2]))
            if (self.ObjectType in ['Vessel', 'vessel']): Displacement = Floater.TimeHistory(SelectedDOF, OC.pnWholeSimulation, OC.oeVessel(self.Position[0], self.Position[1], self.Position[2]))
            
        except:
            raise Exception('Error:: incorrect vesel/buoy type.')

        # Set the time to start at zero and reverse the time values of OrcaFlex outputs on plots.
        tmp_t = abs(Time[0])
        for i,t in enumerate(Time):
            Time[i] += tmp_t

        # To make it start from the offset point.
        ''' <<NumPy arrays are faster and more compact than Python lists. An array consumes less memory and is convenient to use. >>'''
        ArraySizeForTransientTime = int(len(Time)/ self.TEnd_Stage0 * self.TransientTime)
        ArraySizeTimeExtract = len(Time) - ArraySizeForTransientTime

        TruncatedDisplacement, TruncatedTime = np.zeros(ArraySizeTimeExtract), np.zeros(ArraySizeTimeExtract)
        for i, time in enumerate(Time):
            if time >= self.TransientTime:

                TruncatedTime[i-ArraySizeForTransientTime], TruncatedDisplacement[i-ArraySizeForTransientTime] = time, Displacement[i]

        return self.NaturalPeriodCurveFit(TruncatedTime, TruncatedDisplacement, iD)

    
    def PlotResultsEachDoF(self, iD, SelectedDOF, SimFile):
        SpectralPeriod, Freq, Spec = self.PeriodBySpectralDensity(SelectedDOF, SimFile)
        NaturalPeriod, Time, Disp, popt = self.TimeHistoryResults(iD, SelectedDOF, SimFile)

        plt.figure(1)
        plt.plot(Time, self.ExponentialDecay(Time, *popt), ls = 'dotted', color = 'r', label = 'Curve_fit')
        plt.plot(Time, Disp, ls='solid' , color = 'g', label = 'Simulation results')

        #plt.legend()
        plt.xlabel('Time (s)'); plt.ylabel(PlotLabelList[iD])
        plt.title('Natural Period (s)=' + str(NaturalPeriod))
        plt.tight_layout()
        plt.show()

        print('Natural Period=' + str(NaturalPeriod))
        print('Spectral Period={:5.2f} DoF:'.format(SpectralPeriod), SelectedDOF)

        return

    def PlotResultAllDoFs(self):

        #SpectralPeriod = self.PeriodBySpectralDensity(DoF)

        #-------------------------- multiple plots ---------------------------#
        FigRow = 3; FigCol = 2
        fig,axs = plt.subplots(3,2) 
        index = -1

        for Col in np.arange(FigCol):
            for Row in np.arange(FigRow):
               index +=1 

               NP, Time, Displacement, popt  = self.TimeHistoryResults(index, self.DoFList[index], os.path.join(self.ResultsFolder, self.SimFileList[index]))

               plt.sca(axs[Row,Col])                                                                 # Plot allocator
               plt.plot(Time, Displacement, ls='solid', color='g', label= PlotLabelList[index])
               plt.plot(Time, self.ExponentialDecay(Time, *popt), ls = 'dotted', color = 'r', label = 'Curve_fit')
               plt.title(r'$T_{N}$=%2.2f' %NP + "s")
               plt.xlim(Time[0], Time[-1])
               #plt.yticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14])
               plt.xlabel('Time (s)'); plt.ylabel(PlotLabelList[index])
               plt.legend(loc='best',prop={'size':8})
        
        plt.tight_layout()
        plt.show()

        return


    def PlotFreqSpectralDensity(self):

        FigRow = 3; FigCol = 2
        fig,axs = plt.subplots(3,2) 
        index = -1


        for Col in np.arange(FigCol):
            for Row in np.arange(FigRow):
               index +=1 

               NP, Freq, Spec  = self.PeriodBySpectralDensity(self.DoFList[index], os.path.join(self.ResultsFolder, self.SimFileList[index]))

               plt.sca(axs[Row,Col])                                                                 # Plot allocator
               plt.plot(Freq, Spec, ls='solid', color='r', label= PlotLabelList[index])
               #plt.title(r'$T_{N}$=%2.2f' %NP + "s")
               plt.xlim(Freq[0], 0.2)
               plt.xticks([0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2])
               plt.xlabel('Freq (Hs)'); plt.ylabel(PlotLabelList[index])
               plt.legend(loc='best',prop={'size':8})
        
        plt.tight_layout()
        plt.show()
        plt.close()
        sys.stdout.flush()

        return


if __name__ == '__main__':
    raise Exception('This is intended to run as a submodule.')