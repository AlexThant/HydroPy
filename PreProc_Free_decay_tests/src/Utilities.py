""" ====================================================================================================================
author@ Thant Zin Htun
Date@ 6/12/2023 \\\modified
Location@  Tokyo, Japan
Python_version tested@ 3.11.6 & 3.12.0
OracaFLex_version tested@11.4a
For bug report: reach out to me at thant.htun@owcltd.com
\\\
///
                         Module to set up frquently used functions
                         -----------------------------------------

///
\\\
==================================================================================================================== """

import os
#import src.BatchRunsMultiProcessing as BatchRuns


class Tools:

    def __init__(self, WorkingDir, FileDir, FileName, ExcelSheetName) -> None:
        self.WorkingDir = WorkingDir
        self.FileDir = FileDir
        self.FileName = FileName
        self.ExcelSheetName = ExcelSheetName
        return None
    

    def LoadModelFile(self):
        return os.path.join(self.FileDir, self.FileName)
    

    def RunsFolder(self):
        """Create a folder for FreeDecay runs.
     """ 
        FolderName = self.ExcelSheetName+'RunsFolder'
        path = os.path.join(self.WorkingDir, FolderName)
        if not (os.path.exists(path)):
            os.mkdir(path)
        return path
    
    
    def ResultsFolder(self):
        """Create a folder for FreeDecay runs.
     """ 
        path = os.path.join(self.WorkingDir, r'ResultsFolder')
        if not (os.path.exists(path)):
            os.mkdir(path)
        return path
    
    
    def SaveModel(self, model, caseiD):
            
        caseName = f"Case{caseiD:02d}.dat"
        model.SaveData(os.path.join(self.RunsFolder(), caseName))

        print('Saving '+ caseName +' file .......')

        return 
    
    @staticmethod
    def CheckObjectList(model,OC):
        """Check the available objects in the Orcaflex .dat file.
        """
        #model = OC.Model(self.FileDir)
        for obj in model.objects: 
            print(obj, obj.type)
        #inspect.getmembers(model.general, predicate= inspect.ismethod)
        constraints = [obj for obj in model.objects if obj.type ==
        OC.otConstraint]
        print(obj.type, constraints)
        exit()


    def CallBatchRuns(self, NoOfCores):
        """ Runs in the batch mode in parallel.
        """
        print("Importing jobs for batch runs.")
        JobListFolder = self.RunsFolder()

        return BatchRuns.JobRuns(JobListFolder, NoOfCores)


    @staticmethod
    def ElapsedTime(Tstart, Tend):
        hours, rem = divmod(Tend-Tstart,3600)
        minutes, seconds = divmod(rem, 60)
        return print('Elapsed time (hr:min:sec)= {:0>2}:{:0>2}:{:0>2}'.format(int(hours), int(minutes), int(seconds)))


#__*__*__*__*__*__*__*__*__*
if __name__ == '__main__':
    raise Exception('This is a submodule to be imported into another code.')