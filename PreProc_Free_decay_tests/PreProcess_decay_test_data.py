""" ====================================================================================================================
author@ Thant Zin Htun
Date@ 27/12/2023 \\\modified
Location@  Tokyo, Japan
Python_version tested@ 3.11.6 & 3.12.0
OracaFLex_version tested@11.4a
For bug report: reach out to me at thant.htun@owcltd.com
\\\
///
                         Module to shift time series data and group by
                         ---------------------------------------------

///
\\\
==================================================================================================================== """
from timeit import default_timer as timer
import os, sys
import numpy as np
import pandas as pd
from src import DampingCoefficientsAnalysis as DampCoeff
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize']=[12,6]
plt.rcParams.update({'font.size': 12})
line_color_list= ['w', 'orange', 'magenta', 'teal', 'aquamarine', 'r', 'mediumseagreen','blue', 'cyan']

#*****************************************************************************#
#_______________________________Preparatory___________________________________#

Tstart = timer()
Current_Working_Dir = os.path.dirname(os.path.abspath(__file__))  #os.getcwd() be ware of using the latter syntax when your terminal cwd is different from the file directory 

#*****************************************************************************#
#                                                                             #
#                            User Input data                                  #
#                                                                             #
#*****************************************************************************#
ExcelFileName = 'Opti_FreeDecay_Heave_results.xlsx'  #'Opti_Chain_FreeDecay_Surge_results.xlsx' # #
no_of_sheets_to_read = 3
rows_to_skip = 6
#*****************************************************************************#
#                                                                             #
#                            Reading Input data                               #
#                                                                             #
#*****************************************************************************#

class ReadExcelFile:

    def __init__(self) -> None:
        pass


    def Parse_data(self, ExcelSheetName):
        
        try:
            
            time = self.time_stamps(ExcelSheetName)

            if ("surge" in ExcelFileName.lower()):
                X, DoF = self.surge_displacements(ExcelSheetName)

            elif ("sway" in ExcelFileName.lower()):
                X, DoF = self.sway_displacements(ExcelSheetName)
            
            elif ("heave" in ExcelFileName.lower()):
                X, DoF = self.heave_displacements(ExcelSheetName)
            
            else:
                print('....')

        except FileNotFoundError:
            print('Given excel sheet name not does not match.')

        #print(X)

        return time, X, DoF



    def read_input_dataframe(self, ExcelSheetName):

        # Read Excel data
        try:

            path_to_file = os.path.join(Current_Working_Dir, ExcelFileName)

            if os.path.isfile(path_to_file):

                Excel_File = path_to_file     
                df = pd.read_excel(Excel_File, sheet_name= ExcelSheetName, skiprows= rows_to_skip)
                df = df.replace(r'^\s*$', np.nan, regex= True)                                    # Converts blank cells to NaN
                df.dropna(subset = ['Z'], inplace= True)                             # Removing rows where a cell in a given column has a NaN
                df.replace([np.inf, -np.inf], np.nan, inplace=True)
       
        except  FileExistsError:
            print('ERROR: File name error (or)'+ str(ExcelFileName) + 'file does not exist in the same folder!')
                

        return df

    def time_stamps(self, Sheet):
        df = self.read_input_dataframe(Sheet)
        return df['Time (Seconds)'].values
    
    def surge_displacements(self, Sheet):
        df = self.read_input_dataframe(Sheet)
        DoF = 'Surge'
        print('Surge: I am called.')
        return  df['X'].values, DoF
    
    def sway_displacements(self, Sheet):
        df = self.read_input_dataframe(Sheet)
        DoF = 'Sway'
        print('Sway: I am called.')
        return df['Z'].values, DoF
    
    def heave_displacements(self, Sheet):
        df = self.read_input_dataframe(Sheet)
        DoF = 'Heave'
        print('Heave: I am called.')
        return df['Y'].values, DoF


class DataSeries_TimeShift:

    def __init__(self) -> None:
        pass

    def find_min(self, data_series):
        return min(data_series)
    
    def find_max(self, data_series):
        return max(data_series)
    
    def time_shift(self, time, data_series):

        val_max = self.find_min(data_series[1:]); max_ix, = np.where(data_series == val_max); max_index = max_ix[0]
        val_min = self.find_max(data_series[1:]); min_ix, = np.where(data_series == val_min); min_index = min_ix[0]

        if (val_max and  val_min >= 0.0):

            ini_val = [ min(val_max, val_min) if (max_index > min_index)else
                        max(val_max, val_min)
                      ]
        
        elif (val_max and  val_min < 0.0):
              
              ini_val = [ min(val_max, val_min) if (max_index > min_index)else
                          max(val_max, val_min)
                        ]
                               
        else:
            max(abs(val_max), abs(val_min)) * np.sign(val_min + val_min) 

        
        dt = time[1]-time[0]

        "Comma is important here to extract only an integer value (not as a tuplet)"
        ini_val_index, = np.where(data_series == ini_val)
        #print('initial value index', ini_val_index[0], len(data_series), ini_val)

        slicing_index = ini_val_index[0] 
        truncated_dataseries = data_series[slicing_index:]
        len_dataseries = len(truncated_dataseries)

        t_end = len_dataseries*dt
        truncated_time = np.linspace(start=0, stop= t_end, num=len_dataseries, endpoint = True, dtype = float)
        truncated_time = np.round(truncated_time, decimals=2)

        for idx, val in enumerate(truncated_dataseries):

            # No extrapolation is to be done.
            if not idx in [0, len(truncated_dataseries)-1]:

                try:
                    jump = 1
                    if (np.isnan(val) or np.isinf(val) or np.isneginf(val)):                
                        for jj in range(1,101,1):
                            if not (np.isnan(truncated_dataseries[idx+jj]) and np.isinf(truncated_dataseries[idx+jj])):
                                jump = jj
                                break

                        slope = (truncated_dataseries[idx+jump] - truncated_dataseries[idx-1])/jump

                        # Replace Nan with an interpolated value
                        truncated_dataseries[idx] = truncated_dataseries[idx-1] + slope

                except IndexError:
                    print(IndexError)

        return truncated_time, truncated_dataseries, dt


        #print (truncated_dataseries[0:], truncated_time[0:], len(truncated_time), len(truncated_dataseries))


if __name__ == '__main__':

    Mean_Period = 0.0
    fig, ax = plt.subplots(1,1)

    for id in np.arange(start=1, stop=no_of_sheets_to_read+1, dtype= int):
        time, timeseries_data, DoF_name = ReadExcelFile().Parse_data(ExcelSheetName= 'Sheet' + str(id))
        t_new, data_new, dt = DataSeries_TimeShift().time_shift(time, timeseries_data)
        DampData = DampCoeff.DataProcessing()
        Period_curvefit, t_curvefit, data_curvefit = DampData.NaturalPeriodCurveFit(t_new,data_new)
        Mean_Period += Period_curvefit/no_of_sheets_to_read
        
        print('****************************')
        
        xmax, xmin, Etac, Etat, detrend_data = DampCoeff.get_peaks_using_ZeroCrossings(t_new, data_new, dt)

        # plotting the data
        plt.scatter(xmax, Etac, marker='*',  c='k', facecolors= 'white', s = 38)
        plt.scatter(xmin, Etat, marker= 'H', c='k', facecolors= 'white', s= 38)
        ax.annotate('(%s)' % data_new[0], xy=(t_new[0], data_new[0]), textcoords='data', fontsize = 10)
        plt.plot(time, timeseries_data, color= line_color_list[id], linestyle = 'dotted',linewidth= 1.25)
        plt.plot(t_new, data_new, color= line_color_list[id], linewidth= 1.75, label= DoF_name + str(id) + " (m)" )
        plt.plot(t_curvefit, data_curvefit, color= line_color_list[id], linestyle='dashed', linewidth= 1, label= 'Curve_fit'+ DoF_name + str(id) + " (m)" )
        plt.title('Natural Period (s)=' + str(np.round(Mean_Period, decimals= 3))+'s')

         
    plt.legend(loc='best',prop={'size':12})
    plt.tight_layout()
    plt.show()
    plt.close()
    sys.stdout.flush()
    