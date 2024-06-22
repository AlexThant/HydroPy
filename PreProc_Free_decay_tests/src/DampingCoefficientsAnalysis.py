
import scipy as sp
import numpy as np
from scipy import signal as SciSig
from matplotlib import pyplot as plt

#*****************************************************************************#
#                                                                             #
#                       Analysis of Wave Statistics                           #
#                                                                             #
#*****************************************************************************#

def get_peaks_using_ZeroCrossings(t, timeseries_data, dt):
    
    data = SciSig.detrend(timeseries_data, type= 'constant')
    mean_val = np.mean(timeseries_data, axis=-1, keepdims=True)	


	#--------------------------------------------------------------------------
    #           Check the starting point of the 1st first wave
    #                       Upgoing or Downgoing
    
 
    if ((data[0] == 0) and (data[1] > 0)):                          # Upgoing /
        first=1

        
    if (((data[0] == 0) and (data[1] < 0)) or (data[0] != 0)):      # Downgoing \
        for i in range(0,  len(t)-1, 1):
            if ((data[i] < 0) and (data[i+1] > 0)):
                first=i
                break
    #--------------------------------------------------------------------------
    #           Check the end point of the last wave   


    if ((data[-1] == 0) and (data[-2] < 0)):                        # Upgoing /
        last= len(t)


    if (((data[-1] == 0) and (data[-2] > 0)) or (data[-1] != 0)):   # Downgoing \
        for i in range( len(t) -2, 0, -1):
            if ((data[i] < 0) and (data[i+1] > 0)):
                last=i
                break
  
    #print('first', first, 'last', last)
    #exit()
   #--------------------------------------------------------------------------
   #            Detecting zero crossing points from data

    m=-1
    n=-1

    xupcross = np.array([]);       yupcross = np.array([]);
    positionxupcross = np.array([]);
   
    xdowncross = np.array([]);     ydowncross = np.array([])
    positionxdowncross = np.array([]);

    for i in range(first, last+1, 1):

       # Detecting up-crossing zero-crossing points
      
       if i == 1:
           m += 1
           xupcross = np.append(xupcross,dt)
           yupcross = np.append(yupcross,0)
           positionxupcross = np.append(positionxupcross,i)

       if ((i > 1) and (i <  len(t) )):
           if ((data[i] < 0) and (data[i+1] > 0)):
               m = m+1
               xupcross = np.append(xupcross,t[i]-(t[i+1]-t[i])/(data[i+1]-data[i])*data[i])
               yupcross = np.append(yupcross,0)
               positionxupcross = np.append(positionxupcross,i)

       if i ==  len(t) :
           m=m+1
           xupcross = np.append(xupcross,t[-1])
           yupcross = np.append(yupcross,0)
           positionxupcross = np.append(positionxupcross,i)

       # Detecting down-crossing zero-crossing points
       
       if ( 1 < i< len(t) ):
           if ((data[i] > 0) and (data[i+1] < 0)):
               n=n+1
               xdowncross = np.append(xdowncross,t[i]-(t[i+1]-t[i])/(data[i+1]-data[i])*data[i])
               ydowncross = np.append(ydowncross,0)
               positionxdowncross = np.append(positionxdowncross,i)
               if positionxdowncross[n] >= len(t) :
                   positionxdowncross[n] = len(t) -1

    #--------------------------------------------------------------------------
    #           Finding the zero-crossing points of time-series data   

    m = -1
    n = -1
    
    xmax = np.array([])           ;       xmin = np.array([])
    ymax = np.array([])           ;       ymin = np.array([])
    positionxmax = np.array([])   ;       positionxmin = np.array([])
    

    for i in range(np.int64(positionxupcross[0]), np.int64(positionxupcross[-1])+1,1):

        if (1 < i <  len(t)):
			
			# Detecting crest
            if (data[i] > 0 and  data[i-1] and data[i+1]):

                if (m == n): # check if after last crest, the program detects the trough or not (m==n mean it detected and the new y(i,1) is the crest in next wave)

                    m += 1
                    xmax = np.append(xmax,t[i])
                    ymax = np.append(ymax,data[i])
                    positionxmax = np.append(positionxmax,i)

                else:

                    if (m != -1 and data[i] > ymax[m]): #replacing the old crest location with new one if the new one is larger

                        xmax[m] = t[i]
                        ymax[m] = data[i]
                        positionxmax[m] = i


        	# Detecting trough
            if (data[i]< 0 and data[i-1] and data[i+1]):

                if n == m-1: # Check if after last trough, the program detects the crest or not (n==m-1 mean it detected and the new y(i,1) is the trough in next wave)

                    n +=1
                    xmin = np.append(xmin,t[i])
                    ymin = np.append(ymin,data[i])
                    positionxmin = np.append(positionxmin,i)

                else:

                    if (n!= -1 and data[i] < ymin[n]): # Replacing the old crest location with new one if the new one is smaller

                        xmin[n] = t[i]
                        ymin[n] = data[i]
                        positionxmin[n] = i

    
    Crests = np.array([]); Troughs = np.array([])
    T_crests = np.array([]); T_troughs = np.array([])
    Threshold = 1.0e-03


    for i in range (0,len(ymax),1):
        
        if (abs(ymax[i]) > Threshold): 
            Crests = np.append(Crests, ymax[i] + mean_val); 
            T_crests = np.append(T_crests, xmax[i])

        if (abs(ymin[i]) > Threshold): 
            Troughs  = np.append(Troughs, ymin[i] + mean_val)
            T_troughs = np.append(T_troughs, xmin[i])


    return T_crests, T_troughs, Crests, Troughs, data


#*****************************************************************************#
#                         Calculation of Decay coefficients

class DataProcessing:

    def __init__(self):
        self.RefPeriod = 20          # Estimate of natural periods to be found, in seconds
        return


    def ExponentialDecay(self, t, a, b, c, omega, phi, epsilon):
        """ This is an exponential decay function fitting onto the decay test signal to calculate the natural
        period """
        return a + b * np.exp(c * t) * np.cos(omega * t - epsilon - phi)


    def InitialGuessExpoDecay(self, x0_value =0.0):
        """ Set the values of initial guesses for exponential decay for the curve fit.
        """
        y0 = [0, x0_value, 0, np.round((2.0 * np.pi / self.RefPeriod), decimals= 3),0, 1e-10]
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
    def NaturalPeriodCurveFit(self, Xdata, Ydata):
        """
        To calculate the natural period of the system by using the CurveFit.
        """  
        popt, pcov = sp.optimize.curve_fit(self.ExponentialDecay, Xdata, Ydata, p0=(self.InitialGuessExpoDecay(Ydata[0])))
        Omega = popt[3];   #Check the element location of returned frequency value in y0 inside InitialGuessExpoDecay function.
        Y_fitted = self.ExponentialDecay(Xdata, *popt)
        return np.round((2.0 * np.pi / Omega), decimals=3), Xdata, Y_fitted
    

#*****************************************************************************#
#                                   Plotting


def generate_figs(SP_JON, eta, freq, nf, display):
    
    if (display=='on'):
        
        #-------------------------- multiple plots ---------------------------#
        
        L1 = np.arange(1,np.floor(nf/2), dtype= 'int')              # Plot the first half of freq range

        fig,axs = plt.subplots(2,1) 
        plt.sca(axs[0])                                             # First plot identifier
        plt.plot(tt,eta,color='r',label="Free-surface elevation")
        plt.title(r'$H_{s}$=%2.2f' %Hs + "m, " + r'$T_{p}$=%2.2f' %Tp +"s,   CASE: " + str(case_name) + "_s" + str(seed))
        plt.xlim(tt[0],tt[-1])
        plt.yticks([-14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14])
        plt.xlabel('Time (s)'); plt.ylabel(r'$\eta (m)$')
        plt.legend(loc='best',prop={'size':10})

	#

        plt.sca(axs[1])                                             # Second plot identifier
        plt.plot(freq,SP_JON,color='b',label="JONSWAP");            # plt.plot(freq,SP_PMS,color='c',label="Pierson-Moskowitz")
        plt.xlim(freq[L1[0]],freq[L1[-1]])
        plt.xlabel('$f$ (Hz)'); plt.ylabel(r'$S_{f} \ (m^{2}/{Hz})$')
        plt.legend(loc='best',prop={'size':12})
	#
        plt.tight_layout()
        plt.show()

        #---------------------------- save plots -----------------------------#

        #plt.savefig("images/Figure__.png", dpi=1200)

    else:
        pass
    
    